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
#include <libconfig.h>

// --- Счётчик вызовов функции ОДУ ----------------------------
static size_t g_func_calls = 0;
static size_t g_steps      = 0;

// Структура для передачи параметров в функцию ОДУ
struct params {
    double g;
    double k_over_m;
    double vwx, vwy, vwz;
};

// Функция, описывающая систему ОДУ: dy/dt = f(t, y)
int func(double t, const double y[], double f[], void *params) {
    g_func_calls++;
    struct params *p = (struct params *)params;
    double vx = y[3];
    double vy = y[4];
    double vz = y[5];
    double vrx = vx - p->vwx;
    double vry = vy - p->vwy;
    double vrz = vz - p->vwz;
    double vr = sqrt(vrx * vrx + vry * vry + vrz * vrz);

    f[0] = vx;
    f[1] = vy;
    f[2] = vz;
    f[3] = -p->k_over_m * vr * vrx;
    f[4] = -p->g - p->k_over_m * vr * vry;
    f[5] = -p->k_over_m * vr * vrz;
    return GSL_SUCCESS;
}

// --- Функции интерполяции -----------------------------------

// Точка траектории
typedef struct {
    double t, x, y, z, vx, vy, vz;
} traj_point_t;

// Прямая квадратичная интерполяция Лагранжа: q(t) по 3 точкам
static double lagrange_q(double t,
                         double t0, double q0,
                         double t1, double q1,
                         double t2, double q2) {
    double L0 = (t - t1) * (t - t2) / ((t0 - t1) * (t0 - t2));
    double L1 = (t - t0) * (t - t2) / ((t1 - t0) * (t1 - t2));
    double L2 = (t - t0) * (t - t1) / ((t2 - t0) * (t2 - t1));
    return q0 * L0 + q1 * L1 + q2 * L2;
}

// Обратная квадратичная интерполяция: найти t*, при котором y(t*) = 0
static double inverse_quad_root(double t0, double y0,
                                double t1, double y1,
                                double t2, double y2) {
    double L0 = (0.0 - y1) * (0.0 - y2) / ((y0 - y1) * (y0 - y2));
    double L1 = (0.0 - y0) * (0.0 - y2) / ((y1 - y0) * (y1 - y2));
    double L2 = (0.0 - y0) * (0.0 - y1) / ((y2 - y0) * (y2 - y1));
    return t0 * L0 + t1 * L1 + t2 * L2;
}

// --- Загрузка конфигурации ----------------------------------
int load_config(const char *filename, double *v0, double *angle, double *azimuth,
                double *k_over_m, double *x0, double *y0, double *z0,
                double *windspeed, double *winddir) {
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

bool parse_double(const char *str, double *out) {
    char *endptr = NULL;
    errno = 0;
    double val = strtod(str, &endptr);
    if (errno != 0) {
        fprintf(stderr, "Ошибка преобразования числа '%s': %s\n",
                str, strerror(errno));
        return false;
    }
    if (endptr == str || *endptr != '\0') {
        fprintf(stderr, "Ошибка: '%s' не является корректным числом\n", str);
        return false;
    }
    *out = val;
    return true;
}

void print_usage(const char *progname) {
    fprintf(stderr,
        "Использование: %s [ОПЦИИ] [ФАЙЛ_КОНФИГА]\n"
        "\n"
        "Короткие опции:\n"
        "  -v <значение>   Переопределить начальную скорость v0 (м/с)\n"
        "  -a <значение>   Переопределить угол возвышения angle (градусы)\n"
        "  -w <значение>   Переопределить азимут azimuth (градусы)\n"
        "  -x <значение>   Переопределить начальную координату x0 (м)\n"
        "  -y <значение>   Переопределить начальную координату y0 (м)\n"
        "  -z <значение>   Переопределить начальную координату z0 (м)\n"
        "  -h              Показать эту справку\n"
        "\n"
        "Длинные опции:\n"
        "  --wspeed <значение>   Переопределить скорость ветра windspeed (м/с)\n"
        "  --wdir <значение>     Переопределить направление ветра winddir (градусы)\n"
        "\n"
        "ФАЙЛ_КОНФИГА      Путь к конфигурационному файлу\n"
        "                  (по умолчанию: projectile.conf)\n"
        "\n"
        "Параметры командной строки имеют приоритет над значениями из конфига.\n",
        progname);
}

int main(int argc, char *argv[]) {
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

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "v:a:w:x:y:z:hS:D:",
                              long_options, &option_index)) != -1) {
        switch (opt) {
        case 'v':
            if (!parse_double(optarg, &opt_v0)) return EXIT_FAILURE;
            has_v0 = true;
            break;
        case 'a':
            if (!parse_double(optarg, &opt_angle)) return EXIT_FAILURE;
            has_angle = true;
            break;
        case 'w':
            if (!parse_double(optarg, &opt_azimuth)) return EXIT_FAILURE;
            has_azimuth = true;
            break;
        case 'x':
            if (!parse_double(optarg, &opt_x0)) return EXIT_FAILURE;
            has_x0 = true;
            break;
        case 'y':
            if (!parse_double(optarg, &opt_y0)) return EXIT_FAILURE;
            has_y0 = true;
            break;
        case 'z':
            if (!parse_double(optarg, &opt_z0)) return EXIT_FAILURE;
            has_z0 = true;
            break;
        case 'S':
            if (!parse_double(optarg, &opt_windspeed)) return EXIT_FAILURE;
            has_windspeed = true;
            break;
        case 'D':
            if (!parse_double(optarg, &opt_winddir)) return EXIT_FAILURE;
            has_winddir = true;
            break;
        case 'h':
            print_usage(argv[0]);
            return EXIT_SUCCESS;
        default:
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    if (optind < argc) {
        config_file = argv[optind];
    }

    double v0, angle, azimuth, k_over_m, x0, y0, z0;
    double windspeed, winddir;
    printf("Загрузка конфигурации из файла: %s\n", config_file);
    if (load_config(config_file, &v0, &angle, &azimuth, &k_over_m, &x0, &y0, &z0,
                    &windspeed, &winddir) != 0) {
        fprintf(stderr, "Не удалось загрузить конфигурацию\n");
        return EXIT_FAILURE;
    }

    if (has_v0)       v0        = opt_v0;
    if (has_angle)    angle     = opt_angle;
    if (has_azimuth)  azimuth   = opt_azimuth;
    if (has_x0)       x0        = opt_x0;
    if (has_y0)       y0        = opt_y0;
    if (has_z0)       z0        = opt_z0;
    if (has_windspeed) windspeed = opt_windspeed;
    if (has_winddir)  winddir   = opt_winddir;

    printf("Итоговые параметры:\n");
    printf("  v0         = %.2f м/с%s\n",  v0,        has_v0        ? " (переопределено)" : "");
    printf("  angle      = %.2f°%s\n",     angle,     has_angle     ? " (переопределено)" : "");
    printf("  azimuth    = %.2f°%s\n",     azimuth,   has_azimuth   ? " (переопределено)" : "");
    printf("  k_over_m   = %.4f\n",        k_over_m);
    printf("  x0         = %.2f м%s\n",    x0,        has_x0        ? " (переопределено)" : "");
    printf("  y0         = %.2f м%s\n",    y0,        has_y0        ? " (переопределено)" : "");
    printf("  z0         = %.2f м%s\n",    z0,        has_z0        ? " (переопределено)" : "");
    printf("  windspeed  = %.2f м/с%s\n",  windspeed, has_windspeed ? " (переопределено)" : "");
    printf("  winddir    = %.2f°%s\n",     winddir,   has_winddir   ? " (переопределено)" : "");
    printf("\n");

    double g = 9.81;
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

    double y[6] = {x0, y0, z0, vx0, vy0, vz0};

    gsl_odeiv2_system sys = {func, NULL, 6, &p};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
                                                         1e-6, 1e-6, 0.0);

    double t = 0.0;
    double t1 = 0.01;

    printf("# t\tx\ty\tz\tvx\tvy\tvz\n");
    printf("%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
           t, y[0], y[1], y[2], y[3], y[4], y[5]);

    // Буфер последних 3 точек траектории
    traj_point_t pts[3];
    int n_pts = 0;

    // Сохраняем начальную точку
    pts[0] = (traj_point_t){t, y[0], y[1], y[2], y[3], y[4], y[5]};
    n_pts = 1;

    while (y[1] > 0.0 || t == 0.0) {
        int status = gsl_odeiv2_driver_apply(d, &t, t1, y);
        g_steps++;

        if (status != GSL_SUCCESS) {
            printf("Ошибка интегрирования: %d\n", status);
            break;
        }

        printf("%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
               t, y[0], y[1], y[2], y[3], y[4], y[5]);

        // Сохраняем точку в буфер
        if (n_pts < 3) {
            pts[n_pts++] = (traj_point_t){t, y[0], y[1], y[2], y[3], y[4], y[5]};
        } else {
            pts[0] = pts[1];
            pts[1] = pts[2];
            pts[2] = (traj_point_t){t, y[0], y[1], y[2], y[3], y[4], y[5]};
        }

        t1 += 0.01;
    }

    // --- Точное вычисление точки падения ---------------------
    if (n_pts >= 3 && pts[2].y < 0.0) {
        // Обратная квадратичная интерполяция: t* при y=0
        double t_impact = inverse_quad_root(
            pts[0].t, pts[0].y,
            pts[1].t, pts[1].y,
            pts[2].t, pts[2].y);

        // Прямая интерполяция координат и скоростей в t*
        double x_imp  = lagrange_q(t_impact, pts[0].t, pts[0].x,  pts[1].t, pts[1].x,  pts[2].t, pts[2].x);
        double z_imp  = lagrange_q(t_impact, pts[0].t, pts[0].z,  pts[1].t, pts[1].z,  pts[2].t, pts[2].z);
        double vx_imp = lagrange_q(t_impact, pts[0].t, pts[0].vx, pts[1].t, pts[1].vx, pts[2].t, pts[2].vx);
        double vy_imp = lagrange_q(t_impact, pts[0].t, pts[0].vy, pts[1].t, pts[1].vy, pts[2].t, pts[2].vy);
        double vz_imp = lagrange_q(t_impact, pts[0].t, pts[0].vz, pts[1].t, pts[1].vz, pts[2].t, pts[2].vz);

        // Абсолютная скорость
        double v_imp = sqrt(vx_imp * vx_imp + vy_imp * vy_imp + vz_imp * vz_imp);
        double v_horiz = sqrt(vx_imp * vx_imp + vz_imp * vz_imp);

        // Угол к горизонту (отрицательный при падении)
        double angle_impact = atan2(vy_imp, v_horiz) * 180.0 / M_PI;

        // Азимут вектора скорости в точке падения
        double azimuth_impact = atan2(vz_imp, vx_imp) * 180.0 / M_PI;

        // Угол сноса
        double drift = azimuth_impact - azimuth;
        while (drift > 180.0)  drift -= 360.0;
        while (drift < -180.0) drift += 360.0;

        // Дальность
        double dx = x_imp - x0;
        double dz = z_imp - z0;
        double range_horiz = sqrt(dx * dx + dz * dz);
        double range_3d = sqrt(dx * dx + y0 * y0 + dz * dz);

        // Вывод результатов
        printf("#\n");
        printf("# === Точные характеристики точки падения (квадратичная интерполяция) ===\n");
        printf("# Время падения           : %.4f с\n", t_impact);
        printf("# Координаты падения      : x = %.4f м, y = 0 м, z = %.4f м\n", x_imp, z_imp);
        printf("# Горизонтальная дальность: %.2f м\n", range_horiz);
        printf("# Полная 3D-дальность     : %.2f м\n", range_3d);
        printf("# Скорость в точке падения: %.2f м/с (%.1f%% от начальной)\n",
               v_imp, v_imp / v0 * 100.0);
        printf("# Угол к горизонту        : %.2f°\n", angle_impact);
        printf("# Азимут вектора скорости : %.2f°\n", azimuth_impact);
        printf("# Угол сноса              : %.2f°\n", drift);
    } else {
        printf("#\n");
        printf("# Недостаточно данных для интерполяции точки падения\n");
    }

    // --- Итоговая статистика ---------------------------------
    printf("#\n");
    printf("# === Статистика интегрирования ===\n");
    printf("# Всего шагов        : %zu\n", g_steps);
    printf("# Вызовов func       : %zu\n", g_func_calls);
    printf("# Вызовов на шаг     : %.2f\n",
           g_steps ? (double)g_func_calls / g_steps : 0.0);

    gsl_odeiv2_driver_free(d);
    return EXIT_SUCCESS;
}
