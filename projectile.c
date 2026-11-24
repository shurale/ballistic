#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_interp.h>
#include <libconfig.h>

#define DIM 6
#define NHIST 3
#define Y_TARGET 0.0
#define DT_OUT 0.01
#define T_MAX 3600.0
#define NODE_EPS 1e-12

/* --- Физические константы воздуха --- */
#define GAMMA_AIR 1.4
#define R_DRY     287.052874      /* Дж/(кг·К) */
#define MMHG_PA   133.322387415   /* Па в 1 мм рт.ст. */
#define MU_REF    1.716e-5        /* Па·с при 273.15 К */
#define T_MU_REF  273.15
#define SUTH_C    110.4
#define M_LO 0.50                 /* ниже — несжимаемая ветвь */
#define M_HI 0.75                 /* выше — около/сверхзвуковая ветвь */

/* ================= Атмосфера ================= */
typedef struct { double T, p, rho, mu, nu, a; } atmosphere_t;

static void atmosphere_init(atmosphere_t *A, double t_c, double p_mmhg)
{
    A->T   = t_c + 273.15;
    A->p   = p_mmhg * MMHG_PA;
    A->a   = sqrt(GAMMA_AIR * R_DRY * A->T);
    A->rho = A->p / (R_DRY * A->T);
    A->mu  = MU_REF * pow(A->T / T_MU_REF, 1.5) * (T_MU_REF + SUTH_C) / (A->T + SUTH_C);
    A->nu  = A->mu / A->rho;
}

typedef struct { double d, m, area; } body_t;

static void body_init(body_t *b, double d_mm, double mass_g)
{
    b->d = d_mm * 1e-3;
    b->m = mass_g * 1e-3;
    b->area = M_PI * b->d * b->d / 4.0;
}

/* ================= Cd(Re, M) ================= */
#define CD_TAB_MAX 64

typedef struct {
    size_t n;
    double M[CD_TAB_MAX], Cd[CD_TAB_MAX];
    gsl_interp *ip;
    gsl_interp_accel *acc;
} cd_table_t;

static const double cd_m_def[] = {0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,
                                  1.35,1.50,1.75,2.00,2.50,3.00,4.00,5.00,6.00};
static const double cd_c_def[] = {0.480,0.510,0.550,0.630,0.760,0.890,0.950,0.985,
                                  1.000,0.980,0.960,0.930,0.900,0.885,0.875,0.890,0.905};

static int cd_table_init(cd_table_t *t, const double *xm, const double *yc, size_t n)
{
    if (n < 5 || n > CD_TAB_MAX) return -1;
    for (size_t i = 1; i < n; ++i) if (xm[i] <= xm[i-1]) return -1;
    t->n = n;
    memcpy(t->M, xm, n * sizeof(double));
    memcpy(t->Cd, yc, n * sizeof(double));
    t->ip  = gsl_interp_alloc(gsl_interp_akima, n);
    t->acc = gsl_interp_accel_alloc();
    if (!t->ip || !t->acc) return -1;
    gsl_interp_init(t->ip, t->M, t->Cd, n);
    return 0;
}

static void cd_table_free(cd_table_t *t)
{
    if (t->ip)  gsl_interp_free(t->ip);
    if (t->acc) gsl_interp_accel_free(t->acc);
    t->ip = NULL; t->acc = NULL;
}

/* Несжимаемая корреляция Моррисона (Re до 1e6, включает кризис сопротивления) */
static double cd_incompressible(double Re)
{
    if (Re < 1e-6) Re = 1e-6;
    double r5 = Re / 5.0, r263 = Re / 263000.0, r6 = Re / 1e6;
    return 24.0 / Re
         + 2.6 * r5 / (1.0 + pow(r5, 1.52))
         + 0.411 * pow(r263, -7.94) / (1.0 + pow(r263, -8.0))
         + 0.25 * r6 / (1.0 + r6);
}

/* Модифицированная ньютоновская теория: Cd = Cp_max/2, Cp_max — формула Рэлея */
static double cd_newtonian(double M)
{
    const double g = GAMMA_AIR;
    double M2 = M * M;
    double t1 = pow((g + 1) * (g + 1) * M2 / (4 * g * M2 - 2 * (g - 1)), g / (g - 1));
    double t2 = (1 - g + 2 * g * M2) / (g + 1);
    double cp_max = (2.0 / (g * M2)) * (t1 * t2 - 1.0);
    return 0.5 * cp_max;
}

static double cd_sphere(const cd_table_t *t, double Re, double M)
{
    double cd0 = cd_incompressible(Re);
    if (M <= M_LO) return cd0;

    double cd_hi;
    if (M > t->M[t->n - 1]) {
        cd_hi = cd_newtonian(M);
    } else {
        double Me = (M < t->M[0]) ? t->M[0] : M;
        cd_hi = gsl_interp_eval(t->ip, t->M, t->Cd, Me, t->acc);
    }
    if (M >= M_HI) return cd_hi;

    double s = (M - M_LO) / (M_HI - M_LO);
    double w = s * s * (3.0 - 2.0 * s);       /* smoothstep */
    return (1.0 - w) * cd0 + w * cd_hi;
}

static const char *regime_name(double M)
{
    if (M <  M_LO) return "дозвуковой (несжимаемый, Cd(Re))";
    if (M <  M_HI) return "околозвуковой переходный";
    if (M <= 1.40) return "трансзвуковой (скачок-индуцированный отрыв)";
    if (M <= 6.00) return "сверхзвуковой (головная УВ)";
    return "гиперзвуковой (ньютоновский предел)";
}

/* ================= Параметры и правая часть ================= */
struct params {
    double g;
    double vwx, vwy, vwz;
    bool   sphere_mode;
    double k_const;
    atmosphere_t atm;
    body_t body;
    const cd_table_t *tab;
    size_t calls;
    double M_min, M_max, Re_min, Re_max, Cd_min, Cd_max;
};

/* Аэродинамическое состояние в произвольной точке траектории */
static void aero_eval(const struct params *p, const double s[DIM],
                      double *M, double *Re, double *Cd, double *k)
{
    double vrx = s[3] - p->vwx, vry = s[4] - p->vwy, vrz = s[5] - p->vwz;
    double vr = sqrt(vrx*vrx + vry*vry + vrz*vrz);

    *M  = vr / p->atm.a;
    *Re = vr * p->body.d / p->atm.nu;

    if (!p->sphere_mode) { *Cd = NAN; *k = p->k_const; return; }

    if (vr < 1e-12) { *Cd = 0.0; *k = 0.0; *Re = 0.0; return; }  /* стоксов предел */
    *Cd = cd_sphere(p->tab, *Re, *M);
    *k  = *Cd * p->atm.rho * p->body.area / (2.0 * p->body.m);
}

static int func(double t, const double y[], double f[], void *vp)
{
    (void)t;
    struct params *p = vp;
    p->calls++;

    double M, Re, Cd, k;
    aero_eval(p, y, &M, &Re, &Cd, &k);

    double vrx = y[3] - p->vwx, vry = y[4] - p->vwy, vrz = y[5] - p->vwz;
    double vr = sqrt(vrx*vrx + vry*vry + vrz*vrz);

    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];
    f[3] = -k * vr * vrx;
    f[4] = -p->g - k * vr * vry;
    f[5] = -k * vr * vrz;

    if (M  < p->M_min)  p->M_min  = M;   if (M  > p->M_max)  p->M_max  = M;
    if (Re < p->Re_min) p->Re_min = Re;  if (Re > p->Re_max) p->Re_max = Re;
    if (isfinite(Cd)) {
        if (Cd < p->Cd_min) p->Cd_min = Cd;
        if (Cd > p->Cd_max) p->Cd_max = Cd;
    }
    return GSL_SUCCESS;
}

/* ================= Конфигурация ================= */
typedef struct {
    double v0, angle, azimuth, x0, y0, z0;
    double windspeed, winddir;
    double temperature_c, pressure_mmhg;
    double diameter_mm, mass_g, k_over_m;
    char   drag_model[32];
    double tab_m[CD_TAB_MAX], tab_cd[CD_TAB_MAX];
    size_t tab_n;
} settings_t;

/* libconfig различает int и float — принимаем оба */
static bool setting_num(const config_setting_t *s, double *out)
{
    if (!s) return false;
    switch (config_setting_type(s)) {
    case CONFIG_TYPE_FLOAT: *out = config_setting_get_float(s); return true;
    case CONFIG_TYPE_INT:   *out = (double)config_setting_get_int(s); return true;
    case CONFIG_TYPE_INT64: *out = (double)config_setting_get_int64(s); return true;
    default: return false;
    }
}

static bool cfg_num(const config_t *c, const char *path, double *out)
{ return setting_num(config_lookup(c, path), out); }

static void cfg_def(const config_t *c, const char *path, double *dst, double def)
{ if (!cfg_num(c, path, dst)) *dst = def; }

static int load_config(const char *fname, settings_t *S)
{
    config_t cfg;
    config_init(&cfg);
    memset(S, 0, sizeof *S);

    if (!config_read_file(&cfg, fname)) {
        fprintf(stderr, "Ошибка чтения конфигурации: %s:%d - %s\n",
                config_error_file(&cfg), config_error_line(&cfg), config_error_text(&cfg));
        config_destroy(&cfg);
        return -1;
    }

    double dummy;
    if (!cfg_num(&cfg, "v0", &S->v0)) {
        fprintf(stderr, "Ошибка: в %s не задан обязательный параметр v0\n", fname);
        config_destroy(&cfg);
        return -1;
    }
    cfg_def(&cfg, "angle",   &S->angle,   45.0);
    cfg_def(&cfg, "azimuth", &S->azimuth,  0.0);
    cfg_def(&cfg, "x0", &S->x0, 0.0);
    cfg_def(&cfg, "y0", &S->y0, 1.5);
    cfg_def(&cfg, "z0", &S->z0, 0.0);
    cfg_def(&cfg, "windspeed", &S->windspeed, 0.0);
    cfg_def(&cfg, "winddir",   &S->winddir,   0.0);
    cfg_def(&cfg, "temperature_c", &S->temperature_c, 15.0);
    cfg_def(&cfg, "pressure_mmhg", &S->pressure_mmhg, 760.0);
    cfg_def(&cfg, "diameter_mm",   &S->diameter_mm,   6.0);
    cfg_def(&cfg, "mass_g",        &S->mass_g,        0.2);
    cfg_def(&cfg, "k_over_m",      &S->k_over_m,      0.0407);

    const char *dm = NULL;
    if (config_lookup_string(&cfg, "drag_model", &dm) && dm) {
        snprintf(S->drag_model, sizeof S->drag_model, "%s", dm);
    } else {
        /* обратная совместимость: есть k_over_m -> постоянный коэффициент */
        snprintf(S->drag_model, sizeof S->drag_model, "%s",
                 cfg_num(&cfg, "k_over_m", &dummy) ? "constant" : "sphere");
    }

    /* необязательная пользовательская таблица Cd(M) */
    const config_setting_t *sm = config_lookup(&cfg, "cd_table_m");
    const config_setting_t *sc = config_lookup(&cfg, "cd_table_cd");
    if (sm && sc) {
        int n = config_setting_length(sm);
        if (n != config_setting_length(sc) || n < 5 || n > CD_TAB_MAX) {
            fprintf(stderr, "Ошибка: cd_table_m/cd_table_cd разной длины или n<5\n");
            config_destroy(&cfg);
            return -1;
        }
        for (int i = 0; i < n; ++i) {
            if (!setting_num(config_setting_get_elem(sm, i), &S->tab_m[i]) ||
                !setting_num(config_setting_get_elem(sc, i), &S->tab_cd[i])) {
                fprintf(stderr, "Ошибка: элементы cd_table_* должны быть числами\n");
                config_destroy(&cfg);
                return -1;
            }
        }
        S->tab_n = (size_t)n;
    }

    config_destroy(&cfg);
    return 0;
}

static bool parse_double(const char *str, double *out)
{
    char *end = NULL;
    errno = 0;
    double v = strtod(str, &end);
    if (errno != 0) { fprintf(stderr, "Число '%s': %s\n", str, strerror(errno)); return false; }
    if (end == str || *end != '\0') { fprintf(stderr, "'%s' не число\n", str); return false; }
    *out = v;
    return true;
}

static void print_usage(const char *prog)
{
    fprintf(stderr,
        "Использование: %s [ОПЦИИ] [ФАЙЛ_КОНФИГА]\n"
        "  -v  v0 (м/с)          -a  angle (град)      -w  azimuth (град)\n"
        "  -x  x0 (м)            -y  y0 (м)            -z  z0 (м)\n"
        "  -t  температура (°C)  -p  давление (мм рт.ст.)\n"
        "  -d  диаметр (мм)      -m  масса (г)\n"
        "  --wspeed V ветра (м/с)   --wdir направление ветра (град)\n"
        "  --drag constant|sphere   модель сопротивления\n"
        "  -h  справка\n", prog);
}

/* ================= История и обратная интерполяция ================= */
typedef struct { double t, s[DIM]; } sample_t;
typedef struct { double x0, y0, z0, v0, az_rad; } launch_t;

static void hist_push(sample_t *h, size_t *n, double t, const double s[DIM])
{
    for (size_t i = 0; i + 1 < NHIST; ++i) h[i] = h[i+1];
    h[NHIST-1].t = t;
    memcpy(h[NHIST-1].s, s, sizeof(double) * DIM);
    if (*n < NHIST) (*n)++;
}

static size_t impact_interp(const sample_t *h, size_t n, size_t want,
                            double *t_imp, double s_imp[DIM])
{
    double nodes[NHIST]; size_t idx[NHIST]; size_t m = 0;
    for (size_t k = n; k-- > 0 && m < want; ) {
        double yk = h[k].s[1];
        bool dup = false;
        for (size_t j = 0; j < m; ++j)
            if (fabs(yk - nodes[j]) < NODE_EPS) { dup = true; break; }
        if (!dup) { nodes[m] = yk; idx[m] = k; ++m; }
    }
    if (m == 0) return 0;
    if (m == 1) {
        *t_imp = h[idx[0]].t;
        memcpy(s_imp, h[idx[0]].s, sizeof(double) * DIM);
        s_imp[1] = Y_TARGET;
        return 1;
    }
    for (size_t i = 1; i < m; ++i) {
        double yn = nodes[i]; size_t in = idx[i], j = i;
        while (j > 0 && nodes[j-1] > yn) { nodes[j]=nodes[j-1]; idx[j]=idx[j-1]; --j; }
        nodes[j] = yn; idx[j] = in;
    }
    double dd[NHIST], vals[NHIST];
    for (size_t i = 0; i < m; ++i) vals[i] = h[idx[i]].t;
    gsl_poly_dd_init(dd, nodes, vals, m);
    *t_imp = gsl_poly_dd_eval(dd, nodes, m, Y_TARGET);
    s_imp[1] = Y_TARGET;
    for (size_t c = 0; c < DIM; ++c) {
        if (c == 1) continue;
        for (size_t i = 0; i < m; ++i) vals[i] = h[idx[i]].s[c];
        gsl_poly_dd_init(dd, nodes, vals, m);
        s_imp[c] = gsl_poly_dd_eval(dd, nodes, m, Y_TARGET);
    }
    return m;
}

static void report_impact(const launch_t *L, const struct params *p,
                          double t_imp, const double s[DIM], size_t m, bool bracketed,
                          double t_lin, const double sl[DIM], size_t ml)
{
    double dx = s[0] - L->x0, dz = s[2] - L->z0;
    double Rh = hypot(dx, dz), R3 = sqrt(Rh*Rh + L->y0*L->y0);
    double vx = s[3], vy = s[4], vz = s[5];
    double V = sqrt(vx*vx + vy*vy + vz*vz), vh = hypot(vx, vz);
    double elev = atan2(vy, vh) * 180.0 / M_PI;
    double azv  = atan2(vz, vx) * 180.0 / M_PI;
    double ca = cos(L->az_rad), sa = sin(L->az_rad);
    double drift = atan2(ca*vz - sa*vx, ca*vx + sa*vz) * 180.0 / M_PI;

    double M, Re, Cd, k;
    aero_eval(p, s, &M, &Re, &Cd, &k);

    printf("#\n# === Точка падения (y = 0): обратная полиномиальная интерполяция ===\n");
    printf("# Узлов: %zu (степень %zu)\n", m, m - 1);
    if (!bracketed) printf("# ВНИМАНИЕ: y=0 вне интервала узлов (экстраполяция)\n");
    printf("# Время полёта t_imp            = %.6f с\n", t_imp);
    printf("# Координаты x,y,z              = %.6f, %.6f, %.6f м\n", s[0], s[1], s[2]);
    printf("# Дальность горизонтальная R_h   = %.6f м\n", Rh);
    printf("# Расстояние от старта (3D) R    = %.6f м\n", R3);
    printf("# Боковое смещение (снос точки)  = %.6f м\n", -sa*dx + ca*dz);
    printf("# Скорость |V|                   = %.6f м/с", V);
    if (L->v0 > 0) printf("   (%.3f %% от v0)", 100.0 * V / L->v0);
    printf("\n");
    printf("# Угол V к горизонту             = %+.6f град\n", elev);
    printf("# Азимут V                       = %+.6f град\n", azv);
    printf("# Угол сноса V                   = %+.6f град\n", drift);
    printf("# Число Маха в точке падения     = %.6f  (%s)\n", M, regime_name(M));
    printf("# Число Рейнольдса               = %.4e\n", Re);
    if (isfinite(Cd))
        printf("# Cd / k_over_m в падении        = %.6f / %.6e 1/м\n", Cd, k);
    if (m >= 3 && ml >= 2) {
        double dR = fabs(Rh - hypot(sl[0]-L->x0, sl[2]-L->z0));
        double dV = fabs(V - sqrt(sl[3]*sl[3]+sl[4]*sl[4]+sl[5]*sl[5]));
        printf("# Погрешность (квадр-лин): dR=%.3e м, dV=%.3e м/с, dt=%.3e с\n",
               dR, dV, fabs(t_imp - t_lin));
    }
    printf("# IMPACT\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
           s[0], s[1], s[2], t_imp, V, elev, drift, Rh, M);
    printf("#\n");
}

/* ================= main ================= */
int main(int argc, char *argv[])
{
    const char *cfgfile = "projectile.conf";
    double o_v0=0, o_ang=0, o_az=0, o_x=0, o_y=0, o_z=0, o_ws=0, o_wd=0;
    double o_T=0, o_p=0, o_d=0, o_m=0;
    bool h_v0=false,h_ang=false,h_az=false,h_x=false,h_y=false,h_z=false;
    bool h_ws=false,h_wd=false,h_T=false,h_p=false,h_d=false,h_m=false,h_drag=false;
    char o_drag[32] = "";

    static struct option lopts[] = {
        {"wspeed", required_argument, 0, 'S'},
        {"wdir",   required_argument, 0, 'D'},
        {"temp",   required_argument, 0, 't'},
        {"press",  required_argument, 0, 'p'},
        {"diam",   required_argument, 0, 'd'},
        {"mass",   required_argument, 0, 'm'},
        {"drag",   required_argument, 0, 'G'},
        {"help",   no_argument,       0, 'h'},
        {0,0,0,0}
    };

    int opt, oi = 0;
    while ((opt = getopt_long(argc, argv, "v:a:w:x:y:z:t:p:d:m:hS:D:G:", lopts, &oi)) != -1) {
        switch (opt) {
        case 'v': if(!parse_double(optarg,&o_v0)) return EXIT_FAILURE; h_v0=true; break;
        case 'a': if(!parse_double(optarg,&o_ang)) return EXIT_FAILURE; h_ang=true; break;
        case 'w': if(!parse_double(optarg,&o_az)) return EXIT_FAILURE; h_az=true; break;
        case 'x': if(!parse_double(optarg,&o_x)) return EXIT_FAILURE; h_x=true; break;
        case 'y': if(!parse_double(optarg,&o_y)) return EXIT_FAILURE; h_y=true; break;
        case 'z': if(!parse_double(optarg,&o_z)) return EXIT_FAILURE; h_z=true; break;
        case 't': if(!parse_double(optarg,&o_T)) return EXIT_FAILURE; h_T=true; break;
        case 'p': if(!parse_double(optarg,&o_p)) return EXIT_FAILURE; h_p=true; break;
        case 'd': if(!parse_double(optarg,&o_d)) return EXIT_FAILURE; h_d=true; break;
        case 'm': if(!parse_double(optarg,&o_m)) return EXIT_FAILURE; h_m=true; break;
        case 'S': if(!parse_double(optarg,&o_ws)) return EXIT_FAILURE; h_ws=true; break;
        case 'D': if(!parse_double(optarg,&o_wd)) return EXIT_FAILURE; h_wd=true; break;
        case 'G': snprintf(o_drag,sizeof o_drag,"%s",optarg); h_drag=true; break;
        case 'h': print_usage(argv[0]); return EXIT_SUCCESS;
        default:  print_usage(argv[0]); return EXIT_FAILURE;
        }
    }
    if (optind < argc) cfgfile = argv[optind];

    settings_t S;
    fprintf(stderr, "Конфигурация: %s\n", cfgfile);
    if (load_config(cfgfile, &S) != 0) return EXIT_FAILURE;

    if (h_v0) S.v0=o_v0;      if (h_ang) S.angle=o_ang;   if (h_az) S.azimuth=o_az;
    if (h_x)  S.x0=o_x;       if (h_y)   S.y0=o_y;        if (h_z)  S.z0=o_z;
    if (h_ws) S.windspeed=o_ws; if (h_wd) S.winddir=o_wd;
    if (h_T)  S.temperature_c=o_T; if (h_p) S.pressure_mmhg=o_p;
    if (h_d)  S.diameter_mm=o_d;   if (h_m) S.mass_g=o_m;
    if (h_drag) snprintf(S.drag_model,sizeof S.drag_model,"%s",o_drag);

    bool sphere = (strcmp(S.drag_model, "sphere") == 0);
    if (!sphere && strcmp(S.drag_model, "constant") != 0) {
        fprintf(stderr, "drag_model должен быть 'sphere' или 'constant'\n");
        return EXIT_FAILURE;
    }
    if (S.temperature_c <= -273.15 || S.pressure_mmhg <= 0.0) {
        fprintf(stderr, "Некорректные T или p\n"); return EXIT_FAILURE;
    }
    if (sphere && (S.diameter_mm <= 0.0 || S.mass_g <= 0.0)) {
        fprintf(stderr, "Для drag_model=sphere нужны diameter_mm>0 и mass_g>0\n");
        return EXIT_FAILURE;
    }

    cd_table_t tab = {0};
    if (cd_table_init(&tab,
                      S.tab_n ? S.tab_m  : cd_m_def,
                      S.tab_n ? S.tab_cd : cd_c_def,
                      S.tab_n ? S.tab_n  : sizeof cd_m_def / sizeof cd_m_def[0]) != 0) {
        fprintf(stderr, "Ошибка инициализации таблицы Cd(M)\n");
        return EXIT_FAILURE;
    }

    atmosphere_t atm; atmosphere_init(&atm, S.temperature_c, S.pressure_mmhg);
    body_t body;      body_init(&body, S.diameter_mm, S.mass_g);

    double az_rad = S.azimuth * M_PI / 180.0;
    double el_rad = S.angle   * M_PI / 180.0;
    double wd_rad = S.winddir * M_PI / 180.0;

    struct params p;
    memset(&p, 0, sizeof p);
    p.g = 9.81;
    p.vwx = S.windspeed * cos(wd_rad);
    p.vwy = 0.0;
    p.vwz = S.windspeed * sin(wd_rad);
    p.sphere_mode = sphere;
    p.k_const = S.k_over_m;
    p.atm = atm; p.body = body; p.tab = &tab;
    p.M_min = p.Re_min = p.Cd_min = INFINITY;
    p.M_max = p.Re_max = p.Cd_max = 0.0;

    double y[DIM] = { S.x0, S.y0, S.z0,
                      S.v0*cos(el_rad)*cos(az_rad),
                      S.v0*sin(el_rad),
                      S.v0*cos(el_rad)*sin(az_rad) };
    launch_t L = { S.x0, S.y0, S.z0, S.v0, az_rad };

    double M0, Re0, Cd0, k0;
    aero_eval(&p, y, &M0, &Re0, &Cd0, &k0);

    fprintf(stderr, "Атмосфера: T=%.2f K (%.1f °C), p=%.1f Па (%.1f мм рт.ст.)\n",
            atm.T, S.temperature_c, atm.p, S.pressure_mmhg);
    fprintf(stderr, "  a=%.3f м/с, rho=%.5f кг/м^3, mu=%.4e Па·с, nu=%.4e м^2/с\n",
            atm.a, atm.rho, atm.mu, atm.nu);
    fprintf(stderr, "Тело: d=%.2f мм, m=%.3f г, A=%.4e м^2, модель=%s\n",
            S.diameter_mm, S.mass_g, body.area, S.drag_model);
    fprintf(stderr, "Старт: M=%.4f (%s), Re=%.4e, Cd=%.4f, k=%.6e 1/м\n",
            M0, regime_name(M0), Re0, sphere ? Cd0 : NAN, k0);

    gsl_odeiv2_system sys = { func, NULL, DIM, &p };
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
                                                         1e-6, 1e-6, 0.0);
    if (!d) { cd_table_free(&tab); return EXIT_FAILURE; }

    printf("# Атмосфера: T=%.2f K, p=%.1f Pa, a=%.4f m/s, rho=%.6f kg/m^3, nu=%.6e m^2/s\n",
           atm.T, atm.p, atm.a, atm.rho, atm.nu);
    printf("# Тело: d=%.4f m, m=%.6f kg, drag_model=%s\n", body.d, body.m, S.drag_model);
    printf("# t\tx\ty\tz\tvx\tvy\tvz\tMach\tRe\tCd\n");

    sample_t hist[NHIST]; size_t nh = 0;
    memset(hist, 0, sizeof hist);
    double t = 0.0, t1 = DT_OUT;
    size_t steps = 0;
    bool found = false;

    double M, Re, Cd, k;
    aero_eval(&p, y, &M, &Re, &Cd, &k);
    printf("%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.5f\t%.4e\t%.5f\n",
           t, y[0], y[1], y[2], y[3], y[4], y[5], M, Re, isfinite(Cd) ? Cd : 0.0);
    hist_push(hist, &nh, t, y);

    while (t < T_MAX) {
        int st = gsl_odeiv2_driver_apply(d, &t, t1, y);
        steps++;
        if (st != GSL_SUCCESS) {
            fprintf(stderr, "Ошибка интегрирования: %s\n", gsl_strerror(st));
            break;
        }
        aero_eval(&p, y, &M, &Re, &Cd, &k);
        printf("%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.5f\t%.4e\t%.5f\n",
               t, y[0], y[1], y[2], y[3], y[4], y[5], M, Re, isfinite(Cd) ? Cd : 0.0);
        hist_push(hist, &nh, t, y);
        if (y[1] <= Y_TARGET) { found = true; break; }
        t1 += DT_OUT;
    }

    if (found) {
        double ti=0, si[DIM]={0}, tl=0, sl[DIM]={0};
        size_t m  = impact_interp(hist, nh, 3, &ti, si);
        size_t ml = impact_interp(hist, nh, 2, &tl, sl);
        bool br = (nh >= 2) && (hist[nh-1].s[1] <= Y_TARGET) && (hist[nh-2].s[1] >= Y_TARGET);
        if (m) report_impact(&L, &p, ti, si, m, br, tl, sl, ml);
    } else {
        printf("#\n# Точка падения не достигнута\n#\n");
    }

    printf("# === Статистика ===\n");
    printf("# Метод: rk8pd; шагов вывода: %zu; вызовов func: %zu (%.2f на шаг)\n",
           steps, p.calls, steps ? (double)p.calls/steps : 0.0);
    printf("# M:  min=%.4f max=%.4f\n", p.M_min, p.M_max);
    printf("# Re: min=%.4e max=%.4e\n", p.Re_min, p.Re_max);
    if (sphere) printf("# Cd: min=%.4f max=%.4f\n", p.Cd_min, p.Cd_max);

    gsl_odeiv2_driver_free(d);
    cd_table_free(&tab);
    return EXIT_SUCCESS;
}
