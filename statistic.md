# Получение статистики работы `gsl_odeiv2_driver`

## 1. Что доступно «из коробки»

Структура `gsl_odeiv2_driver` является непрозрачной (opaque) — её поля не документированы как публичный API, и напрямую к ним обращаться нельзя. Однако GSL предоставляет функции-аксессоры для связанных объектов (`step`, `evolve`, `control`), через которые можно получить нужную информацию.

| Информация | Как получить | Ограничения |
|---|---|---|
| Последний шаг интегрирования `h` | `gsl_odeiv2_step_last_step(step)` | Нужен доступ к объекту `step` |
| Имя численного метода | `gsl_odeiv2_step_name(step)` | ✔ |
| Порядок метода | `gsl_odeiv2_step_order(step)` | ✔ |
| Количество вызовов `func` | Только через обёртку | Нет встроенного счётчика |
| Количество шагов | Только через обёртку | Нет встроенного счётчика |
| Мин/макс допустимый шаг | `gsl_odeiv2_driver_set_minmax()` | Только установка |

## 2. Ключевой приём: переход от `driver` к `evolve`

`gsl_odeiv2_driver` — это высокоуровневая обёртка, которая скрывает от пользователя три компонента: `step`, `control` и `evolve`. Если нам нужен контроль над шагом и статистикой, правильнее собрать эти три объекта вручную:

```
driver  =  step  +  control  +  evolve
(удобно, но непрозрачно)    (полный контроль)
```

## 3. Обновлённый фрагмент `main()` со сбором статистики

```c
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

// --- Счётчик вызовов функции ОДУ ----------------------------
// Глобальная статистика (для простоты; в production — через params)
static size_t g_func_calls = 0;
static size_t g_steps      = 0;
static double g_h_min      =  INFINITY;
static double g_h_max      = 0.0;
static double g_h_sum      = 0.0;

// --- Структура параметров -----------------------------------
struct params {
    double g;
    double k_over_m;
    double vwx, vwy, vwz;
};

// Обёртка над функцией ОДУ: считает вызовы
static int func_wrapper(double t, const double y[], double f[], void *params) {
    g_func_calls++;

    struct params *p = (struct params *)params;
    double vx = y[3];
    double vy = y[4];
    double vz = y[5];
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

// ... (load_config, parse_double, print_usage — без изменений) ...

int main(int argc, char *argv[]) {
    // ... (загрузка конфига, разбор опций — без изменений) ...
    // ... (инициализация y[6], p, sys) ...

    gsl_odeiv2_system sys = {func_wrapper, NULL, 6, &p};

    // --- Ручная сборка драйвера из трёх компонентов ----------
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;

    gsl_odeiv2_step    *step    = gsl_odeiv2_step_alloc(T, 6);
    gsl_odeiv2_control *control = gsl_odeiv2_control_y_new(1e-6, 0.0);
    gsl_odeiv2_evolve  *evolve  = gsl_odeiv2_evolve_alloc(6);

    // Начальный предлагаемый шаг
    double h = 1e-3;
    double t = 0.0;
    double t1 = 0.0;
    const double t_end = 100.0;  // страховочный предел

    printf("# t\tx\ty\tz\th\tfunc_calls\n");

    while (y[1] > 0.0 && t < t_end) {
        t1 += 0.01;  // желаемый момент вывода

        // Интегрируем до t1, делая ровно столько шагов, сколько нужно
        while (t < t1 && y[1] > 0.0) {
            int status = gsl_odeiv2_evolve_apply(
                evolve, control, step,
                &sys, &t, t1, &h, y);

            if (status != GSL_SUCCESS) {
                fprintf(stderr, "Ошибка: %s\n", gsl_strerror(status));
                goto cleanup;
            }

            // --- Собираем статистику ---
            double h_used = gsl_odeiv2_step_last_step(step);
            g_steps++;
            g_h_sum += h_used;
            if (h_used < g_h_min) g_h_min = h_used;
            if (h_used > g_h_max) g_h_max = h_used;
        }

        printf("%.4f\t%.4f\t%.4f\t%.4f\t%.2e\t%zu\n",
               t, y[0], y[1], y[2], gsl_odeiv2_step_last_step(step),
               g_func_calls);
    }

    // --- Итоговая статистика ---------------------------------
    printf("\n# === Статистика интегрирования ===\n");
    printf("# Метод              : %s (порядок %d)\n",
           gsl_odeiv2_step_name(step), gsl_odeiv2_step_order(step));
    printf("# Всего шагов        : %zu\n", g_steps);
    printf("# Вызовов func       : %zu\n", g_func_calls);
    printf("# Вызовов на шаг     : %.2f\n",
           g_steps ? (double)g_func_calls / g_steps : 0.0);
    printf("# Итоговое время t   : %.4f с\n", t);
    printf("# Шаг h: min=%.2e  max=%.2e  avg=%.2e\n",
           g_h_min, g_h_max, g_steps ? g_h_sum / g_steps : 0.0);

cleanup:
    gsl_odeiv2_evolve_free(evolve);
    gsl_odeiv2_control_free(control);
    gsl_odeiv2_step_free(step);
    return EXIT_SUCCESS;
}
```

## 4. Ключевые функции API

### Функции `gsl_odeiv2_step_*`

```c
// Последний реально использованный шаг (после gsl_odeiv2_evolve_apply)
double gsl_odeiv2_step_last_step(const gsl_odeiv2_step *step);

// Имя метода ("rk8pd", "rk4", "rkf45", "msadams", ...)
const char *gsl_odeiv2_step_name(const gsl_odeiv2_step *step);

// Порядок метода (8 для rk8pd)
unsigned int gsl_odeiv2_step_order(const gsl_odeiv2_step *step);
```

### Функции `gsl_odeiv2_evolve_*`

```c
// Один адаптивный шаг: обновляет t, h и y
int gsl_odeiv2_evolve_apply(gsl_odeiv2_evolve *evolve,
                            gsl_odeiv2_control *con,
                            gsl_odeiv2_step *step,
                            const gsl_odeiv2_system *sys,
                            double *t, double t1,
                            double *h, double y[]);

// Сброс состояния (перед новой траекторией)
int gsl_odeiv2_evolve_reset(gsl_odeiv2_evolve *evolve);
```

### Функции `gsl_odeiv2_driver_*` (если всё-таки оставляете driver)

```c
// Сбросить driver и начальный шаг
int gsl_odeiv2_driver_reset_hstart(gsl_odeiv2_driver *d, double hstart);
int gsl_odeiv2_driver_reset(gsl_odeiv2_driver *d);

// Ограничить диапазон шага (полезно для отладки)
int gsl_odeiv2_driver_set_minmax(gsl_odeiv2_driver *d, double hmin, double hmax);
```

## 5. Почему обёртка, а не глобальные переменные — лучше

В примере выше счётчики сделаны глобальными для краткости. В production-коде правильнее хранить их в `struct params` — тогда программа остаётся reentrant-safe и можно запускать несколько интегрирований параллельно:

```c
struct params {
    double g, k_over_m, vwx, vwy, vwz;
    size_t func_calls;   // <-- счётчик внутри параметров
};

static int func_wrapper(double t, const double y[], double f[], void *vp) {
    struct params *p = vp;
    p->func_calls++;
    // ...
}
```

## 6. Что ещё можно извлечь

| Задача | Решение |
|---|---|
| Оценить «жёсткость» задачи | Если `h` резко падает в определённых точках — система приближается к жёсткой; стоит попробовать `gsl_odeiv2_step_msadams` или `gsl_odeiv2_step_rk2imp` |
| Задать минимальный шаг | `gsl_odeiv2_driver_set_minmax(d, 1e-9, 1.0)` |
| Перезапуск с новым начальным шагом | `gsl_odeiv2_driver_reset_hstart(d, 1e-4)` |
| Ловить ошибки GSL | `gsl_set_error_handler(&my_handler)` — свой обработчик вместо `abort()` |
| Получить оценку локальной ошибки | `gsl_odeiv2_evolve_apply` возвращает её неявно через `control`; для явной оценки используйте `gsl_odeiv2_control_errlevel()` |

## 7. Ожидаемая статистика для `rk8pd`

Метод Принца–Дорманда 8(9) — это **embedded method**: на каждом шаге вычисляется 13 стадий (т.е. 13 вызовов `func`), из которых 12 дают решение 8-го порядка, а все 13 — решение 9-го порядка для оценки ошибки. Поэтому отношение `func_calls / steps` для `rk8pd` должно быть **около 13** (плюс несколько вызовов на инициализацию и rejected steps). Если у вас получается сильно больше — значит, контроллер часто отвергает шаги и делает повторные попытки (признак слишком жёстких допусков или сингулярности).

## Источники

1. **GNU Scientific Library Reference Manual — ODE Initial Value Problems**  
   Разделы «Step Functions», «Evolve Functions», «Driver Functions»:  
   [https://www.gnu.org/software/gsl/doc/html/odeiv.html](https://www.gnu.org/software/gsl/doc/html/odeiv.html)

2. **Prince, P. J., Dormand, J. R. (1981). «High order embedded Runge-Kutta formulae», J. Comput. Appl. Math. 7(1), 67–75**  
   Оригинальная статья по методу `rk8pd`:  
   [https://doi.org/10.1016/0771-050X(81)90010-3](https://doi.org/10.1016/0771-050X(81)90010-3)

3. **GSL source code — `ode-initval2/`**  
   Реализация `gsl_odeiv2_step_last_step` и структура `gsl_odeiv2_step`:  
   [https://git.savannah.gnu.org/cgit/gsl.git/tree/ode-initval2](https://git.savannah.gnu.org/cgit/gsl.git/tree/ode-initval2)
