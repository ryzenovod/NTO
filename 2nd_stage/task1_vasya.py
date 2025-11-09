import sys
import math

EPS = 1e-9
r_planet = 100.0

def deg2rad (deg):
    return deg * math.pi / 180.0

def sph_to_cart(r, phi_deg, psi_deg):
    """Сферические координаты (r, φ, ψ в градусах) -> декартовы (x, y, z)."""
    phi = deg2rad(phi_deg)      # Перевод широты (φ) из градусов в радианы
    psi = deg2rad(psi_deg)      # Перевод азимута (ψ) из градусов в радианы
    cphi = math.cos(phi)        # cos(φ) используется дважды — считаем один раз
    return [
        r * cphi * math.cos(psi),  # x = r * cosφ * cosψ
        r * cphi * math.sin(psi),  # y = r * cosφ * sinψ
        r * math.sin(phi),         # z = r * sinφ
    ]

def dot (a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] #cкалярное произведение векторов

def norm (a):
    return math.sqrt(dot(a, a)) # длина вектора

def normalize(v):
    l = norm(v)                               # Считаем длину вектора
    if l < EPS:                               # Если длина почти ноль, нормировать нельзя
        return [0.0, 0.0, 0.0]                # Возвращаем нулевой вектор (безопасный возврат)
    return [v[0]/l, v[1]/l, v[2]/l]           # Делим каждую компоненту на длину → единичный вектор

def rodrigues_rotate(axis, angle, v):
    """Поворот вектора v вокруг оси axis на угол angle (радианы) по формуле Родрига."""
    ax = normalize(axis)                      # Ось поворота должна быть единичной
    ux, uy, uz = ax                           # Раскладываем ось по координатам для формул
    c = math.cos(angle)                       # cos(θ)
    s = math.sin(angle)                       # sin(θ)
    x, y, z = v                               # Исходные координаты поворачиваемого вектора
    # Матрица поворота: R = I*c + (1-c)·u·u^T + [u]_×·s
    rx = (c + (1-c)*ux*ux)*x + ((1-c)*ux*uy - s*uz)*y + ((1-c)*ux*uz + s*uy)*z  # Новая компонента x
    ry = ((1-c)*uy*ux + s*uz)*x + (c + (1-c)*uy*uy)*y + ((1-c)*uy*uz - s*ux)*z  # Новая компонента y
    rz = ((1-c)*uz*ux - s*uy)*x + ((1-c)*uz*uy + s*ux)*y + (c + (1-c)*uz*uz)*z  # Новая компонента z
    return [rx, ry, rz]                      # Возвращаем повернутый вектор

def visible_from_point(S, P, R=r_planet):
    """Проверка видимости спутника S из точки P на сфере радиуса R: S·P ≥ R² (с допуском)."""
    return dot(S, P) >= R*R - EPS

def read_floats():
    return [float(x) for x in input().split()]

def main():
    phi_p, psi_p = read_floats()
    P = sph_to_cart(r_planet, phi_p, psi_p)
    t = float(sys.stdin.readline().strip())
    k = int(sys.stdin.readline().strip())
    answers = []
    for i in range(1, k+1):                   # Перебираем орбиты с номерами 1..k
        # Для каждой орбиты нужно прочитать 8 чисел:
        # m_i, r_n, φ_n, ψ_n, R_i, φ_i1, ψ_i1, T_i
        parts = []                             # Временный список для накопления параметров орбиты
        while len(parts) < 8:                  # Некоторые тесты могут разбивать параметры на несколько строк
            parts += read_floats()             # Добавляем числа из следующей строки, пока не наберём 8

        m_i = int(parts[0])                    # m_i — число спутников на этой орбите
        rn, phi_n, psi_n = parts[1], parts[2], parts[3]        # Вектор нормали орбитальной плоскости в сферич. коорд.
        R_i, phi_i1, psi_i1 = parts[4], parts[5], parts[6]     # Радиус орбиты и старт первого спутника (φ, ψ)
        T_i = parts[7]                         # Период обращения по орбите

        axis = sph_to_cart(rn, phi_n, psi_n)   # Переводим нормаль орбиты в декартовы координаты
        s0   = sph_to_cart(R_i, phi_i1, psi_i1)# Начальное положение первого спутника в декартовых координатах

        for j in range(m_i):                   # Перебираем все спутники на орбите: j = 0..m_i-1
            angle = (2.0*math.pi)*t/ T_i + (2.0*math.pi)*j/ m_i  # Угол поворота каждого спутника: ход по времени + сдвиг по номеру
            S = rodrigues_rotate(axis, angle, s0)                # Поворачиваем стартовый вектор вокруг оси орбиты
            if visible_from_point(S, P, r_planet):               # Проверяем критерий видимости S·P ≥ R²
                answers.append((i, j+1))                         # Записываем пару (№ орбиты, № спутника) с 1
    if not answers:                           # Если ни один спутник не виден из точки P
        print(0)                              # Печатаем 0 (по формату задачи)
        return                                

    answers.sort()                            # Сортируем пары по возрастанию (сначала орбита, затем спутник)
    print(len(answers))                       # Печатаем количество видимых спутников
    for o, s in answers:                      # Перебираем все пары по порядку
        print(o, s)

if __name__ == '__main__':
    main()