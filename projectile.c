#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
int func(double t, const double y[], double f[], void *params) {
    struct params *p = (struct params *)params;
    double vx = y[2];
    double vy = y[3];
    double v = sqrt(vx * vx + vy * vy); // Модуль скорости

    f[0] = vx;                          // dx/dt
    f[1] = vy;                          // dy/dt
    f[2] = -p->k_over_m * v * vx;       // dvx/dt
    f[3] = -p->g - p->k_over_m * v * vy;// dvy/dt

    return GSL_SUCCESS;
}

// Функция для чтения конфигурационного файла
int load_config(const char *filename, double *v0, double *angle,
                double *k_over_m, double *x0, double *y0) {
    config_t cfg;
    config_init(&cfg);

    // Чтение конфигурационного файла
    if (!config_read_file(&cfg, filename)) {
        fprintf(stderr, "Ошибка чтения конфигурации: %s:%d - %s\n",
                config_error_file(&cfg), config_error_line(&cfg),
                config_error_text(&cfg));
        config_destroy(&cfg);
        return -1;
    }

    // Чтение параметров
    if (!config_lookup_float(&cfg, "v0", v0)) {
        fprintf(stderr, "Ошибка: параметр 'v0' не найден в %s\n", filename);
        config_destroy(&cfg);
        return -1;
    }

    if (!config_lookup_float(&cfg, "angle", angle)) {
        fprintf(stderr, "Ошибка: параметр 'angle' не найден в %s\n", filename);
        config_destroy(&cfg);
        return -1;
    }

    if (!config_lookup_float(&cfg, "k_over_m", k_over_m)) {
        fprintf(stderr, "Ошибка: параметр 'k_over_m' не найден в %s\n", filename);
        config_destroy(&cfg);
        return -1;
    }

    if (!config_lookup_float(&cfg, "x0", x0)) {
        fprintf(stderr, "Ошибка: параметр 'x0' не найден в %s\n", filename);
        config_destroy(&cfg);
        return -1;
    }

    if (!config_lookup_float(&cfg, "y0", y0)) {
        fprintf(stderr, "Ошибка: параметр 'y0' не найден в %s\n", filename);
        config_destroy(&cfg);
        return -1;
    }

    config_destroy(&cfg);
    return 0;
}

int main(int argc, char *argv[]) {
    const char *config_file = "projectile.conf";

    // Если передан аргумент, используем его как имя конфигурационного файла
    if (argc > 1) {
        config_file = argv[1];
    }

    // Параметры из конфигурационного файла
    double v0, angle, k_over_m, x0, y0;

    // Загрузка конфигурации
    printf("Загрузка конфигурации из файла: %s\n", config_file);
    if (load_config(config_file, &v0, &angle, &k_over_m, &x0, &y0) != 0) {
        fprintf(stderr, "Не удалось загрузить конфигурацию\n");
        return EXIT_FAILURE;
    }

    printf("Параметры:\n");
    printf("  v0 = %.2f м/с\n", v0);
    printf("  angle = %.2f°\n", angle);
    printf("  k_over_m = %.4f\n", k_over_m);
    printf("  x0 = %.2f м\n", x0);
    printf("  y0 = %.2f м\n", y0);
    printf("\n");

    // Константы
    double g = 9.81;

    // Преобразование угла в радианы
    double angle_rad = angle * M_PI / 180.0;

    struct params p = {g, k_over_m};

    // Начальные условия: [x0, y0, vx0, vy0]
    double y[4] = {x0, y0, v0 * cos(angle_rad), v0 * sin(angle_rad)};

    // Инициализация системы ОДУ
    gsl_odeiv2_system sys = {func, NULL, 4, &p};

    // Создание драйвера с адаптивным шагом
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
                                                         1e-6, 1e-6, 0.0);

    double t = 0.0;
    double t1 = 0.01; // Шаг вывода данных на экран

    // Выводим заголовок для gnuplot
    printf("# t\tx\ty\tvx\tvy\n");
    printf("%.3f\t%.4f\t%.4f\t%.4f\t%.4f\n", t, y[0], y[1], y[2], y[3]);

    // Интегрирование до момента падения на землю (y <= 0)
    while (y[1] > 0.0 || t == 0.0) {
        int status = gsl_odeiv2_driver_apply(d, &t, t1, y);

        if (status != GSL_SUCCESS) {
            printf("Ошибка интегрирования: %d\n", status);
            break;
        }

        printf("%.3f\t%.4f\t%.4f\t%.4f\t%.4f\n", t, y[0], y[1], y[2], y[3]);
        t1 += 0.01;
    }

    // Освобождение памяти
    gsl_odeiv2_driver_free(d);
    return EXIT_SUCCESS;
}
