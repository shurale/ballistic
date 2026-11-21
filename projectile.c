#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <libconfig.h>

// Структура для передачи параметров в функцию ОДУ
struct params {
    double g;          // Ускорение свободного падения
    double k_over_m;   // Коэффициент сопротивления, деленный на массу (k)
};

// Функция, описывающая систему ОДУ: dy/dt = f(t, y)
// Вектор состояния: [x, y, z, vx, vy, vz]
int func(double t, const double y[], double f[], void *params) {
    struct params *p = (struct params *)params;
    double vx = y[3];
    double vy = y[4];
    double vz = y[5];
    double v = sqrt(vx * vx + vy * vy + vz * vz); // Модуль скорости в 3D

    f[0] = vx;                            // dx/dt
    f[1] = vy;                            // dy/dt
    f[2] = vz;                            // dz/dt
    f[3] = -p->k_over_m * v * vx;         // dvx/dt
    f[4] = -p->g - p->k_over_m * v * vy;  // dvy/dt
    f[5] = -p->k_over_m * v * vz;         // dvz/dt

    return GSL_SUCCESS;
}

// Функция для чтения конфигурационного файла
int load_config(const char *filename, double *v0, double *angle, double *azimuth,
                double *k_over_m, double *x0, double *y0, double *z0) {
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
        !config_lookup_float(&cfg, "z0", z0)) {
        fprintf(stderr, "Ошибка: в %s отсутствуют обязательные параметры "
                "(v0, angle, azimuth, k_over_m, x0, y0, z0)\n", filename);
        config_destroy(&cfg);
        return -1;
    }

    config_destroy(&cfg);
    return 0;
}

// Безопасное преобразование строки в double с обработкой ошибок
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

// Вывод справки по использованию
void print_usage(const char *progname) {
    fprintf(stderr,
        "Использование: %s [ОПЦИИ] [ФАЙЛ_КОНФИГА]\n"
        "\n"
        "ОПЦИИ:\n"
        "  -v <значение>   Переопределить начальную скорость v0 (м/с)\n"
        "  -a <значение>   Переопределить угол возвышения angle (градусы)\n"
        "  -w <значение>   Переопределить азимут azimuth (градусы)\n"
        "  -x <значение>   Переопределить начальную координату x0 (м)\n"
        "  -y <значение>   Переопределить начальную координату y0 (м)\n"
        "  -z <значение>   Переопределить начальную координату z0 (м)\n"
        "  -h              Показать эту справку\n"
        "\n"
        "ФАЙЛ_КОНФИГА      Путь к конфигурационному файлу\n"
        "                  (по умолчанию: projectile.conf)\n"
        "\n"
        "Параметры командной строки имеют приоритет над значениями из конфига.\n",
        progname);
}

int main(int argc, char *argv[]) {
    const char *config_file = "projectile.conf";

    // Значения из командной строки и флаги их наличия
    double opt_v0 = 0.0, opt_angle = 0.0, opt_azimuth = 0.0;
    double opt_x0 = 0.0, opt_y0 = 0.0, opt_z0 = 0.0;
    bool has_v0 = false, has_angle = false, has_azimuth = false;
    bool has_x0 = false, has_y0 = false, has_z0 = false;

    // Разбор опций командной строки
    int opt;
    while ((opt = getopt(argc, argv, "v:a:w:x:y:z:h")) != -1) {
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
        case 'h':
            print_usage(argv[0]);
            return EXIT_SUCCESS;
        default:
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    // Позиционный аргумент после опций — имя конфигурационного файла
    if (optind < argc) {
        config_file = argv[optind];
    }

    // Загрузка конфигурации
    double v0, angle, azimuth, k_over_m, x0, y0, z0;
    printf("Загрузка конфигурации из файла: %s\n", config_file);
    if (load_config(config_file, &v0, &angle, &azimuth, &k_over_m, &x0, &y0, &z0) != 0) {
        fprintf(stderr, "Не удалось загрузить конфигурацию\n");
        return EXIT_FAILURE;
    }

    // Переопределение значений из командной строки
    if (has_v0)       v0      = opt_v0;
    if (has_angle)    angle   = opt_angle;
    if (has_azimuth)  azimuth = opt_azimuth;
    if (has_x0)       x0      = opt_x0;
    if (has_y0)       y0      = opt_y0;
    if (has_z0)       z0      = opt_z0;

    // Вывод итоговых параметров
    printf("Итоговые параметры:\n");
    printf("  v0       = %.2f м/с%s\n",  v0,      has_v0      ? " (переопределено)" : "");
    printf("  angle    = %.2f°%s\n",     angle,   has_angle   ? " (переопределено)" : "");
    printf("  azimuth  = %.2f°%s\n",     azimuth, has_azimuth ? " (переопределено)" : "");
    printf("  k_over_m = %.4f\n",        k_over_m);
    printf("  x0       = %.2f м%s\n",    x0,      has_x0      ? " (переопределено)" : "");
    printf("  y0       = %.2f м%s\n",    y0,      has_y0      ? " (переопределено)" : "");
    printf("  z0       = %.2f м%s\n",    z0,      has_z0      ? " (переопределено)" : "");
    printf("\n");

    // Константы
    double g = 9.81;
    double angle_rad   = angle   * M_PI / 180.0;
    double azimuth_rad = azimuth * M_PI / 180.0;

    struct params p = {g, k_over_m};

    // Начальные условия: [x0, y0, z0, vx0, vy0, vz0]
    // Компоненты скорости через сферические координаты
    double vx0 = v0 * cos(angle_rad) * cos(azimuth_rad);
    double vy0 = v0 * sin(angle_rad);
    double vz0 = v0 * cos(angle_rad) * sin(azimuth_rad);

    double y[6] = {x0, y0, z0, vx0, vy0, vz0};

    // Инициализация системы ОДУ (размерность 6)
    gsl_odeiv2_system sys = {func, NULL, 6, &p};

    // Создание драйвера с адаптивным шагом
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
                                                         1e-6, 1e-6, 0.0);

    double t = 0.0;
    double t1 = 0.01; // Шаг вывода данных на экран

    // Выводим заголовок для gnuplot
    printf("# t\tx\ty\tz\tvx\tvy\tvz\n");
    printf("%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
           t, y[0], y[1], y[2], y[3], y[4], y[5]);

    // Интегрирование до момента падения на землю (y <= 0)
    while (y[1] > 0.0 || t == 0.0) {
        int status = gsl_odeiv2_driver_apply(d, &t, t1, y);

        if (status != GSL_SUCCESS) {
            printf("Ошибка интегрирования: %d\n", status);
            break;
        }

        printf("%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
               t, y[0], y[1], y[2], y[3], y[4], y[5]);
        t1 += 0.01;
    }

    // Освобождение памяти
    gsl_odeiv2_driver_free(d);
    return EXIT_SUCCESS;
}
