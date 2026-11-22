#!/usr/bin/awk -f
# k_over_m.awk -- расчёт k = Cd*rho_air*A/(2m) для шара
# Пример: awk -f k_over_m.awk -v d=18 -v rho=11340 -v cd=0.47
BEGIN {
    PI  = 4 * atan2(1, 1)
    d   = (d   == "" ? 6      : d)   / 1000   # диаметр, мм -> м
    rho = (rho == "" ? 11340  : rho)         # плотность материала, кг/м^3
    cd  = (cd  == "" ? 0.47   : cd)          # коэффициент сопротивления
    air = (air == "" ? 1.225  : air)         # плотность воздуха, кг/м^3
    g   = 9.81

    A = PI * d^2 / 4
    V = PI * d^3 / 6
    m = rho * V
    k = cd * air * A / (2 * m)               # == 3*cd*air/(4*rho*d)

    printf "d        = %.1f мм\n",  d * 1000
    printf "rho_mat  = %.0f кг/м^3\n", rho
    printf "m        = %.4f кг (%.2f г)\n", m, m * 1000
    printf "A        = %.6e м^2\n", A
    printf "Cd       = %.3f\n", cd
    printf "k_over_m = %.6e 1/м   <-- в projectile.conf\n", k
    printf "1/k      = %.1f м\n", 1 / k
    printf "v_term   = %.1f м/с\n", sqrt(g / k)
}
