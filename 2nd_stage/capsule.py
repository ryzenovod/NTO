import sys
import math


def dot(ax, ay, az, bx, by, bz):
    """Скалярное произведение двух 3D-векторов."""
    return ax * bx + ay * by + az * bz


def main():
    data = sys.stdin.read().strip().split()
    if not data:
        print(0)
        return

    it = iter(data)

    # ---------- Чтение вершин ----------
    try:
        n = int(next(it))
    except StopIteration:
        print(0)
        return

    vertices = []
    sum_x = sum_y = sum_z = 0.0

    for _ in range(n):
        x = float(next(it))
        y = float(next(it))
        z = float(next(it))
        vertices.append((x, y, z))
        sum_x += x
        sum_y += y
        sum_z += z

    if n == 0:
        # Нет вершин -> многогранника фактически нет -> цилиндр не впихнуть
        print(0)
        return

    centroid_x = sum_x / n
    centroid_y = sum_y / n
    centroid_z = sum_z / n

    # ---------- Чтение граней и построение плоскостей ----------
    m = int(next(it))
    planes = []

    for _ in range(m):
        k = int(next(it))
        if k < 3:
            # Вырожденная грань, просто пропускаем
            for __ in range(k):
                _ = next(it)
            continue

        idx0 = int(next(it))
        idx1 = int(next(it))
        idx2 = int(next(it))

        # Считываем оставшиеся индексы грани, но для плоскости они не нужны
        for __ in range(k - 3):
            _ = next(it)

        v0 = vertices[idx0]
        v1 = vertices[idx1]
        v2 = vertices[idx2]

        # Векторные произведения для нормали
        ax = v1[0] - v0[0]
        ay = v1[1] - v0[1]
        az = v1[2] - v0[2]
        bx = v2[0] - v0[0]
        by = v2[1] - v0[1]
        bz = v2[2] - v0[2]

        nx = ay * bz - az * by
        ny = az * bx - ax * bz
        nz = ax * by - ay * bx

        len_sq = nx * nx + ny * ny + nz * nz
        if len_sq < 1e-18:
            # Практически нулевая нормаль -> вырожденная грань -> пропускаем
            continue

        # Плоскость: n·x + d = 0 через v0
        d = -(nx * v0[0] + ny * v0[1] + nz * v0[2])

        # Ориентация нормали внутрь многогранника (через центроид)
        s = nx * centroid_x + ny * centroid_y + nz * centroid_z + d
        if s > 0.0:
            nx = -nx
            ny = -ny
            nz = -nz
            d = -d

        length = math.sqrt(len_sq)
        nx /= length
        ny /= length
        nz /= length
        d /= length

        planes.append((nx, ny, nz, d))

    if not planes:
        # Без плоскостей мы не можем определить границы -> консервативно 0
        print(0)
        return

    # ---------- Чтение параметров цилиндра ----------
    try:
        X0 = float(next(it))
        Y0 = float(next(it))
        Z0 = float(next(it))
        U = float(next(it))
        V = float(next(it))
        W = float(next(it))
        ta = float(next(it))
        tb = float(next(it))
    except StopIteration:
        # Некорректный формат входа
        print(0)
        return

    if ta > tb:
        ta, tb = tb, ta

    # Направление оси
    axis_len_sq = U * U + V * V + W * W
    if axis_len_sq < 1e-18:
        # Слишком короткий вектор направления оси -> цилиндр вырожден
        print(0)
        return

    axis_len = math.sqrt(axis_len_sq)
    axis_unit_x = U / axis_len
    axis_unit_y = V / axis_len
    axis_unit_z = W / axis_len

    # Концы оси цилиндра
    Pa_x = X0 + U * ta
    Pa_y = Y0 + V * ta
    Pa_z = Z0 + W * ta

    Pb_x = X0 + U * tb
    Pb_y = Y0 + V * tb
    Pb_z = Z0 + W * tb

    # ---------- Предподсчёт g_max и n_perp для каждой плоскости ----------
    EPS_AXIS = 1e-9       # допуск для проверки оси
    EPS_R_ZERO = 1e-7     # порог, меньше которого радиус считаем нулевым
    SAFETY_EPS = 1e-9     # запас (отступ) от плоскостей для радиуса > 0

    gmax_list = []
    nperp_list = []

    R_hi = float('inf')

    for (nx, ny, nz, d) in planes:
        # Значения плоскости в концах оси
        va = nx * Pa_x + ny * Pa_y + nz * Pa_z + d
        vb = nx * Pb_x + ny * Pb_y + nz * Pb_z + d
        gmax = va if va > vb else vb

        # ----- Проверка: ось целиком внутри многогранника -----
        if gmax > EPS_AXIS:
            # Есть точка оси, лежащая снаружи (по этой грани),
            # значит ни один цилиндр с R > 0 не может целиком поместиться
            print(0)
            return

        gmax_list.append(gmax)

        # Компонента нормали, перпендикулярная оси
        ndotL = nx * axis_unit_x + ny * axis_unit_y + nz * axis_unit_z
        n_perp_sq = 1.0 - ndotL * ndotL

        # Защита от отрицательных нулей из-за округления
        if n_perp_sq < 0.0:
            if n_perp_sq > -1e-12:
                n_perp_sq = 0.0
            else:
                # Сильное отрицательное значение невозможно в корректной геометрии,
                # но консервативно считаем н_perп = 0
                n_perp_sq = 0.0

        n_perp = math.sqrt(n_perp_sq)
        nperp_list.append(n_perp)

        # Аналитическая верхняя граница радиуса по этой плоскости
        if n_perp > 0.0:
            gap = -gmax  # сколько "запаса" до плоскости вдоль нормали на оси
            if gap <= 0.0:
                # Ось уже на границе или вне -> радиуса запаса нет
                R_limit = 0.0
            else:
                R_limit = gap / n_perp
            if R_limit < R_hi:
                R_hi = R_limit

    # ----- Проверка: есть ли вообще положительная аналитическая граница -----
    if not math.isfinite(R_hi) or R_hi <= EPS_R_ZERO:
        # Или радиус вообще ничем не ограничен аналитически (что маловероятно),
        # или максимальный радиус слишком мал -> считаем, что R > 0 недопустим.
        print(0)
        return

    # ---------- Бинарный поиск по радиусу ----------

    def is_radius_ok(R: float) -> bool:
        """Проверка, что цилиндр с данным радиусом R полностью внутри P
        с учетом безопасного отступа SAFETY_EPS.
        Возвращает False, если по какой-либо плоскости цилиндр выходит наружу.
        """
        # Цилиндр R=0 тут не проверяется, используется только для R>0
        for i, (nx, ny, nz, d) in enumerate(planes):
            gmax = gmax_list[i]
            n_perp = nperp_list[i]

            if n_perp == 0.0:
                # Нормаль плоскости параллельна оси цилиндра:
                # расстояние до плоскости не меняется при увеличении R.
                # Достаточно того, что ось внутри (мы уже проверили выше),
                # но на всякий случай проверим небольшим допуском:
                if gmax > EPS_AXIS:
                    # Если ось таки вышла наружу, любой цилиндр тоже выйдет
                    return False
                # Иначе по этой плоскости ограничений на R нет
                continue

            # Максимальное значение плоскости по всем точкам цилиндра:
            # g_max + R * ||n_perp||
            value = gmax + R * n_perp

            # Требуем, чтобы цилиндр был с небольшим отрицательным
            # запасом внутри (value <= -SAFETY_EPS).
            if value > -SAFETY_EPS:
                return False

        return True

    lo = 0.0
    hi = R_hi

    # Бинарный поиск: ищем максимальный R, для которого is_radius_ok(R) == True
    for _ in range(60):  # 2^-60 ~ 1e-18, достаточно для double
        mid = 0.5 * (lo + hi)
        if mid <= 0.0:
            # Механическая защита, но на практике mid > 0 (так как hi > EPS_R_ZERO)
            lo = mid
            continue

        if is_radius_ok(mid):
            lo = mid
        else:
            hi = mid

    best_R = lo

    # ---------- Финальное решение ----------
    if best_R <= EPS_R_ZERO:
        # Радиус либо не существует, либо меньше численного шума:
        # считаем, что нет положительного радиуса
        print(0)
    else:
        print("1 {:.5f}".format(best_R))


if __name__ == "__main__":
    main()
