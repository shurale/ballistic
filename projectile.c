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
#include <libconfig.h>

#define DIM      6        /* размерность системы: x, y, z, vx, vy, vz */
#define NHIST    3        /* число узлов интерполяции */
#define Y_TARGET 0.0      /* целевая высота (земля) */
#define DT_OUT   0.01     /* шаг вывода данных, с */
#define T_MAX    3600.0   /* предохранительный предел времени, с */
#define NODE_EPS 1e-12    /* порог различимости узлов по высоте, м */

/* --- Статистика ------------------------------------------- */
static size_t g_func_calls = 0;
static size_t g_steps      = 0;

/* --- Параметры системы ОДУ -------------------------------- */
struct params {
    double g;
    double k_over_m;
    double vwx, vwy, vwz;
};

static int func(double t, const double y[], double f[], void *params)
{
    (void)t;                       /* система автономна */
    g_func_calls++;

    const struct params *p = (const struct params *)params;
    double vx = y[3], vy = y[4], vz = y[5];

    double vrx = vx - p->vwx;
    double vry = vy - p->vwy;
    double vrz = vz - p->vwz;
    double vr  = sqrt(vrx*vrx + vry*vry + vrz*vrz);

    f[0] = vx;
    f[1] = vy;
    f[2] = vz;
    f[3] = -p->k_over_m * vr * vrx;
    f[4] = -p->g - p->k_over_m * vr * vry;
    f[5] = -p->k_over_m * vr * vrz;
    return GSL_SUCCESS;
}

/* --- Конфигурация ----------------------------------------- */
static int load_config(const char *filename, double *v0, double *angle, double *azimuth,
                       double *k_over_m, double *x0, double *y0, double *z0,
                       double *windspeed, double *winddir)
{
    config_t cfg;
    config_init(&cfg);

    if (!config_read_file(&cfg, filename)) {
        fprintf(stderr, "Ошибка чтения конфигурации: %s:%d - %s\n",
                config_error_file(&cfg), config_error_line(&cfg),
                config_error_text(&cfg));
        config_destroy(&cfg);
        return -1;
    }

    if (!config_lookup_float(&cfg, "v0", v0) ||
        !config_lookup_float(&cfg, "angle", angle) ||
        !config_lookup_float(&cfg, "azimuth", azimuth) ||
        !config_lookup_float(&cfg, "k_over_m", k_over_m) ||
        !config_lookup_float(&cfg, "x0", x0) ||
        !config_lookup_float(&cfg, "y0", y0) ||
        !config_lookup_float(&cfg, "z0", z0) ||
        !config_lookup_float(&cfg, "windspeed", windspeed) ||
        !config_lookup_float(&cfg, "winddir", winddir)) {
        fprintf(stderr, "Ошибка: в %s отсутствуют обязательные параметры\n", filename);
        config_destroy(&cfg);
        return -1;
    }

    config_destroy(&cfg);
    return 0;
}

static bool parse_double(const char *str, double *out)
{
    char *endptr = NULL;
    errno = 0;
    double val = strtod(str, &endptr);

    if (errno != 0) {
        fprintf(stderr, "Ошибка преобразования числа '%s': %s\n", str, strerror(errno));
        return false;
    }
    if (endptr == str || *endptr != '\0') {
        fprintf(stderr, "Ошибка: '%s' не является корректным числом\n", str);
        return false;
    }
    *out = val;
    return true;
}

static void print_usage(const char *progname)
{
    fprintf(stderr,
        "Использование: %s [ОПЦИИ] [ФАЙЛ_КОНФИГА]\n"
        "\n"
        "Короткие опции:\n"
        "  -v <значение>   Начальная скорость v0 (м/с)\n"
        "  -a <значение>   Угол возвышения angle (градусы)\n"
        "  -w <значение>   Азимут запуска azimuth (градусы)\n"
        "  -x <значение>   Начальная координата x0 (м)\n"
        "  -y <значение>   Начальная координата y0 (м)\n"
        "  -z <значение>   Начальная координата z0 (м)\n"
        "  -h              Показать эту справку\n"
        "\n"
        "Длинные опции:\n"
        "  --wspeed <значение>   Скорость ветра windspeed (м/с)\n"
        "  --wdir <значение>     Направление ветра winddir (градусы)\n"
        "\n"
        "Данные и комментарии -> stdout, диагностика -> stderr.\n",
        progname);
}

/* ==========================================================
   История последних точек и обратная интерполяция
   ========================================================== */
typedef struct { double t; double s[DIM]; } sample_t;

typedef struct { double x0, y0, z0, v0, az_rad; } launch_t;

/* Сдвиг истории: [0] — самая старая точка, [NHIST-1] — самая свежая */
static void hist_push(sample_t *h, size_t *n, double t, const double s[DIM])
{
    for (size_t i = 0; i + 1 < NHIST; ++i) h[i] = h[i + 1];
    h[NHIST - 1].t = t;
    memcpy(h[NHIST - 1].s, s, sizeof(double) * DIM);
    if (*n < NHIST) (*n)++;
}

/* Обратная интерполяция: высота — независимая переменная.
   want = 3 -> квадратичная, want = 2 -> линейная (для оценки погрешности).
   Возвращает число реально использованных узлов. */
static size_t impact_interp(const sample_t *h, size_t n, size_t want,
                            double *t_imp, double s_imp[DIM])
{
    double nodes[NHIST];
    size_t idx[NHIST];
    size_t m = 0;

    /* Отбираем самые свежие точки с различимыми высотами */
    for (size_t k = n; k-- > 0 && m < want; ) {
        double yk = h[k].s[1];
        bool dup = false;
        for (size_t j = 0; j < m; ++j)
            if (fabs(yk - nodes[j]) < NODE_EPS) { dup = true; break; }
        if (!dup) { nodes[m] = yk; idx[m] = k; ++m; }
    }
    if (m == 0) return 0;

    if (m == 1) {                       /* вырожденный случай */
        *t_imp = h[idx[0]].t;
        memcpy(s_imp, h[idx[0]].s, sizeof(double) * DIM);
        s_imp[1] = Y_TARGET;
        return 1;
    }

    /* Сортируем узлы по возрастанию высоты (для устойчивости схемы) */
    for (size_t i = 1; i < m; ++i) {
        double yn = nodes[i];
        size_t in = idx[i], j = i;
        while (j > 0 && nodes[j - 1] > yn) {
            nodes[j] = nodes[j - 1];
            idx[j]   = idx[j - 1];
            --j;
        }
        nodes[j] = yn;
        idx[j]   = in;
    }

    double dd[NHIST], vals[NHIST];

    /* t(y) */
    for (size_t i = 0; i < m; ++i) vals[i] = h[idx[i]].t;
    gsl_poly_dd_init(dd, nodes, vals, m);
    *t_imp = gsl_poly_dd_eval(dd, nodes, m, Y_TARGET);

    /* x(y), z(y), vx(y), vy(y), vz(y); y(0) == 0 по построению */
    s_imp[1] = Y_TARGET;
    for (size_t c = 0; c < DIM; ++c) {
        if (c == 1) continue;
        for (size_t i = 0; i < m; ++i) vals[i] = h[idx[i]].s[c];
        gsl_poly_dd_init(dd, nodes, vals, m);
        s_imp[c] = gsl_poly_dd_eval(dd, nodes, m, Y_TARGET);
    }
    return m;
}

/* --- Отчёт о точке падения (только комментарии '#') -------- */
static void report_impact(const launch_t *L,
                          double t_imp, const double s[DIM], size_t m, bool bracketed,
                          double t_lin, const double sl[DIM], size_t ml)
{
    double dx = s[0] - L->x0;
    double dz = s[2] - L->z0;
    double range_h = hypot(dx, dz);                        /* горизонтальная дальность */
    double dist3d  = sqrt(range_h * range_h + L->y0 * L->y0); /* полное 3D-расстояние */

    double vx = s[3], vy = s[4], vz = s[5];
    double vabs = sqrt(vx * vx + vy * vy + vz * vz);
    double vh   = hypot(vx, vz);
    double elev = atan2(vy, vh) * 180.0 / M_PI;            /* угол к горизонту */
    double az_v = atan2(vz, vx) * 180.0 / M_PI;            /* азимут скорости */

    double ca = cos(L->az_rad), sa = sin(L->az_rad);
    /* знаковый угол между (ca,sa) и (vx,vz) в плоскости XZ */
    double drift   = atan2(ca * vz - sa * vx, ca * vx + sa * vz) * 180.0 / M_PI;
    double along   =  ca * dx + sa * dz;                   /* продольное смещение */
    double lateral = -sa * dx + ca * dz;                   /* боковое смещение */

    printf("#\n");
    printf("# === Точка падения (y = %.1f): обратная полиномиальная интерполяция ===\n",
           Y_TARGET);
    printf("# Узлов интерполяции: %zu (степень полинома: %zu)\n", m, m - 1);
    if (!bracketed)
        printf("# ВНИМАНИЕ: y=0 вне интервала узлов — результат является экстраполяцией\n");
    printf("# Время полёта t_imp             = %.6f с\n", t_imp);
    printf("# Координаты падения x,y,z       = %.6f, %.6f, %.6f м\n", s[0], s[1], s[2]);
    printf("# Дальность горизонтальная R_h    = %.6f м\n", range_h);
    printf("# Расстояние от старта (3D) R     = %.6f м\n", dist3d);
    printf("# Смещение вдоль азимута          = %.6f м\n", along);
    printf("# Боковое смещение (снос точки)   = %.6f м  (>0 — в сторону азимута+90)\n", lateral);
    printf("# Скорость |V|                    = %.6f м/с", vabs);
    if (L->v0 > 0.0)
        printf("   (%.3f %% от v0 = %.3f м/с)", 100.0 * vabs / L->v0, L->v0);
    printf("\n");
    printf("# Угол вектора V к горизонту      = %+.6f град  (<0 — вниз)\n", elev);
    printf("# Азимут вектора V                = %+.6f град\n", az_v);
    printf("# Угол сноса вектора V            = %+.6f град  (>0 — к росту азимута)\n", drift);

    if (m >= 3 && ml >= 2) {
        double dR = fabs(range_h - hypot(sl[0] - L->x0, sl[2] - L->z0));
        double dV = fabs(vabs - sqrt(sl[3]*sl[3] + sl[4]*sl[4] + sl[5]*sl[5]));
        printf("# Оценка погрешности (квадр - лин): dR = %.3e м, dV = %.3e м/с, dt = %.3e с\n",
               dR, dV, fabs(t_imp - t_lin));
    }

    /* Машиночитаемая строка для awk/gnuplot */
    printf("# IMPACT\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
           s[0], s[1], s[2], t_imp, vabs, elev, drift, range_h);
    printf("#\n");
}

/* ========================================================== */
int main(int argc, char *argv[])
{
    const char *config_file = "projectile.conf";

    double opt_v0 = 0.0, opt_angle = 0.0, opt_azimuth = 0.0;
    double opt_x0 = 0.0, opt_y0 = 0.0, opt_z0 = 0.0;
    double opt_windspeed = 0.0, opt_winddir = 0.0;
    bool has_v0 = false, has_angle = false, has_azimuth = false;
    bool has_x0 = false, has_y0 = false, has_z0 = false;
    bool has_windspeed = false, has_winddir = false;

    static struct option long_options[] = {
        {"wspeed", required_argument, 0, 'S'},
        {"wdir",   required_argument, 0, 'D'},
        {"help",   no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int opt, option_index = 0;
    while ((opt = getopt_long(argc, argv, "v:a:w:x:y:z:hS:D:",
                              long_options, &option_index)) != -1) {
        switch (opt) {
        case 'v': if (!parse_double(optarg, &opt_v0)) return EXIT_FAILURE; has_v0 = true; break;
        case 'a': if (!parse_double(optarg, &opt_angle)) return EXIT_FAILURE; has_angle = true; break;
        case 'w': if (!parse_double(optarg, &opt_azimuth)) return EXIT_FAILURE; has_azimuth = true; break;
        case 'x': if (!parse_double(optarg, &opt_x0)) return EXIT_FAILURE; has_x0 = true; break;
        case 'y': if (!parse_double(optarg, &opt_y0)) return EXIT_FAILURE; has_y0 = true; break;
        case 'z': if (!parse_double(optarg, &opt_z0)) return EXIT_FAILURE; has_z0 = true; break;
        case 'S': if (!parse_double(optarg, &opt_windspeed)) return EXIT_FAILURE; has_windspeed = true; break;
        case 'D': if (!parse_double(optarg, &opt_winddir)) return EXIT_FAILURE; has_winddir = true; break;
        case 'h': print_usage(argv[0]); return EXIT_SUCCESS;
        default:  print_usage(argv[0]); return EXIT_FAILURE;
        }
    }
    if (optind < argc) config_file = argv[optind];

    double v0, angle, azimuth, k_over_m, x0, y0, z0, windspeed, winddir;

    /* Вся диагностика — в stderr, чтобы stdout оставался чистым файлом данных */
    fprintf(stderr, "Загрузка конфигурации из файла: %s\n", config_file);
    if (load_config(config_file, &v0, &angle, &azimuth, &k_over_m, &x0, &y0, &z0,
                    &windspeed, &winddir) != 0) {
        fprintf(stderr, "Не удалось загрузить конфигурацию\n");
        return EXIT_FAILURE;
    }

    if (has_v0)        v0        = opt_v0;
    if (has_angle)     angle     = opt_angle;
    if (has_azimuth)   azimuth   = opt_azimuth;
    if (has_x0)        x0        = opt_x0;
    if (has_y0)        y0        = opt_y0;
    if (has_z0)        z0        = opt_z0;
    if (has_windspeed) windspeed = opt_windspeed;
    if (has_winddir)   winddir   = opt_winddir;

    fprintf(stderr, "Итоговые параметры:\n");
    fprintf(stderr, "  v0         = %.2f м/с%s\n", v0, has_v0 ? " (переопределено)" : "");
    fprintf(stderr, "  angle      = %.2f град%s\n", angle, has_angle ? " (переопределено)" : "");
    fprintf(stderr, "  azimuth    = %.2f град%s\n", azimuth, has_azimuth ? " (переопределено)" : "");
    fprintf(stderr, "  k_over_m   = %.6f 1/м\n", k_over_m);
    fprintf(stderr, "  x0,y0,z0   = %.2f, %.2f, %.2f м\n", x0, y0, z0);
    fprintf(stderr, "  windspeed  = %.2f м/с%s\n", windspeed, has_windspeed ? " (переопределено)" : "");
    fprintf(stderr, "  winddir    = %.2f град%s\n", winddir, has_winddir ? " (переопределено)" : "");

    const double g = 9.81;
    double angle_rad   = angle   * M_PI / 180.0;
    double azimuth_rad = azimuth * M_PI / 180.0;
    double winddir_rad = winddir * M_PI / 180.0;

    double vwx = windspeed * cos(winddir_rad);
    double vwy = 0.0;
    double vwz = windspeed * sin(winddir_rad);

    struct params p = {g, k_over_m, vwx, vwy, vwz};

    double vx0 = v0 * cos(angle_rad) * cos(azimuth_rad);
    double vy0 = v0 * sin(angle_rad);
    double vz0 = v0 * cos(angle_rad) * sin(azimuth_rad);
    double y[DIM] = {x0, y0, z0, vx0, vy0, vz0};

    launch_t L = {x0, y0, z0, v0, azimuth_rad};

    gsl_odeiv2_system sys = {func, NULL, DIM, &p};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
                                                        1e-6, 1e-6, 0.0);
    if (!d) { fprintf(stderr, "Не удалось создать драйвер GSL\n"); return EXIT_FAILURE; }

    sample_t hist[NHIST];
    size_t nh = 0;
    memset(hist, 0, sizeof hist);

    double t = 0.0, t1 = DT_OUT;
    bool impact_found = false;

    printf("# t\tx\ty\tz\tvx\tvy\tvz\n");
    printf("%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
           t, y[0], y[1], y[2], y[3], y[4], y[5]);
    hist_push(hist, &nh, t, y);

    while (t < T_MAX) {
        int status = gsl_odeiv2_driver_apply(d, &t, t1, y);
        g_steps++;
        if (status != GSL_SUCCESS) {
            fprintf(stderr, "Ошибка интегрирования: %s\n", gsl_strerror(status));
            break;
        }

        hist_push(hist, &nh, t, y);
        printf("%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
               t, y[0], y[1], y[2], y[3], y[4], y[5]);

        if (y[1] <= Y_TARGET) { impact_found = true; break; }
        t1 += DT_OUT;
    }

    /* === Блок точки падения: после данных, до отладки === */
    if (impact_found) {
        double ti = 0.0, si[DIM] = {0};
        double tl = 0.0, sl[DIM] = {0};
        size_t m  = impact_interp(hist, nh, 3, &ti, si);   /* квадратичная */
        size_t ml = impact_interp(hist, nh, 2, &tl, sl);   /* линейная     */
        bool bracketed = (nh >= 2) &&
                         (hist[nh - 1].s[1] <= Y_TARGET) &&
                         (hist[nh - 2].s[1] >= Y_TARGET);
        if (m > 0)
            report_impact(&L, ti, si, m, bracketed, tl, sl, ml);
        else
            printf("#\n# Не удалось построить интерполяцию точки падения\n#\n");
    } else {
        printf("#\n# Точка падения не достигнута (ошибка интегрирования или t > %.0f с)\n#\n",
               T_MAX);
    }

    /* === Отладочная информация === */
    printf("# === Статистика интегрирования ===\n");
    printf("# Метод              : rk8pd (Prince-Dormand 8(9))\n");
    printf("# Шагов вывода       : %zu\n", g_steps);
    printf("# Вызовов func       : %zu\n", g_func_calls);
    printf("# Вызовов на шаг     : %.2f\n",
           g_steps ? (double)g_func_calls / g_steps : 0.0);
    printf("# Шаг вывода DT_OUT  : %.4f с\n", DT_OUT);

    gsl_odeiv2_driver_free(d);
    return EXIT_SUCCESS;
}
