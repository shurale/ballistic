#include <stdio.h>
#include <math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

// Структура для передачи параметров в функцию ОДУ
struct params {
    double g;       // Ускорение свободного падения
    double k_over_m; // Коэффициент сопротивления, деленный на массу (k)
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

int main(void) {
    // Физические параметры (пример для страйкбольного шарика 0.2г, 6мм)
    // v0 = 90 м/с, угол 45 градусов
    double v0 = 90.0;
    double angle = 45.0 * M_PI / 180.0;
    
    // k = (Cd * rho * A) / (2 * m)
    // Cd = 0.47, rho = 1.225, A = pi * (0.003)^2, m = 0.0002
    double k_over_m = 0.0407; 
    double g = 9.81;

    struct params p = {g, k_over_m};

    // Начальные условия: [x0, y0, vx0, vy0]
    double y[4] = {0.0, 1.5, v0 * cos(angle), v0 * sin(angle)}; // Старт с высоты 1.5м
    
    // Инициализация системы ОДУ
    // gsl_odeiv2_system: функция, якобиан (для нежёстких можно NULL), размерность, параметры
    gsl_odeiv2_system sys = {func, NULL, 4, &p};

    // Создание драйвера с адаптивным шагом
    // gsl_odeiv2_step_rk8pd - метод Рунге-Кутты 8(9) порядка
    // 1e-6 - абсолютная погрешность, 0.0 - относительная погрешность
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 
                                                         1e-6, 1e-6, 0.0);

    double t = 0.0;
    double t1 = 0.01; // Шаг вывода данных на экран (не шаг интегрирования!)

    // Выводим заголовок для gnuplot
    printf("# t\tx\ty\tvx\tvy\n");
    printf("%.3f\t%.4f\t%.4f\t%.4f\t%.4f\n", t, y[0], y[1], y[2], y[3]);

    // Интегрирование до момента падения на землю (y <= 0)
    while (y[1] > 0.0 || t == 0.0) {
        // gsl_odeiv2_driver_apply интегрирует от t до t1, 
        // автоматически подбирая внутренний шаг метода
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
    return 0;
}
