# Решение задачи "Великая озёрная регата"
# Подробные комментарии к геометрии и шагам алгоритма.

import sys
import math

EPS = 1e-9      # допуск для сравнения с нулём
EPS_T = 1e-7    # минимальный t > 0, чтобы не считать точку наблюдателя B "перекрытием"

# ---------- Вспомогательные функции для чтения ввода ----------

data = sys.stdin.read().strip().split()
if not data:
    # Пустой ввод – на всякий случай
    print(0)
    sys.exit(0)

it = iter(data)

def next_float():
    return float(next(it))

def next_int():
    return int(next(it))

# ---------- Чтение исходных данных ----------

# Точка наблюдателя B
xB = next_float()
yB = next_float()
zB = next_float()

# Направление камеры l
lx = next_float()
ly = next_float()
lz = next_float()

# Кол-во лодок и гор
N = next_int()
M = next_int()

# Лодки: сохраняем вершины пирамидок (флажки)
boats = []
for _ in range(N):
    xb = next_float()
    yb = next_float()
    hb = next_float()
    # Вершина лодки (флажок) находится в (xb, yb, hb)
    boats.append((xb, yb, hb))

# Горы: конусы
mountains = []
for _ in range(M):
    xm = next_float()
    ym = next_float()
    hm = next_float()
    rm = next_float()
    mountains.append((xm, ym, hm, rm))

# ---------- Построение базиса камеры (f, r, u) ----------

# f — направление "вперёд" (нормированный вектор l)
len_l = math.sqrt(lx*lx + ly*ly + lz*lz)
if len_l < EPS:
    # Нулевой вектор направления — камера неопределена.
    # Можно считать, что ничего не видно.
    print(0)
    sys.exit(0)

fx = lx / len_l
fy = ly / len_l
fz = lz / len_l

# r — горизонтальная ось камеры, должна лежать в плоскости z=0
# Обычно берём вектор, перпендикулярный проекции l на плоскость xy:
# r0 = (-ly, lx, 0). Если проекция нулевая (камера строго вверх/вниз), берём (1, 0, 0).
if abs(lx) + abs(ly) > EPS:
    r0x = -ly
    r0y = lx
    r0z = 0.0
else:
    # l вертикален, выбираем произвольное горизонтальное направление
    r0x = 1.0
    r0y = 0.0
    r0z = 0.0

len_r0 = math.sqrt(r0x*r0x + r0y*r0y + r0z*r0z)
rx = r0x / len_r0
ry = r0y / len_r0
rz = r0z / len_r0  # будет 0, т.к. r0z=0

# u — вертикальная ось камеры, получается как векторное произведение f × r
ux = fy*rz - fz*ry
uy = fz*rx - fx*rz
uz = fx*ry - fy*rx

# ---------- Предобработка гор (конусов) относительно точки B ----------

# Для каждой горы заранее считаем:
#  - ex, ey: координаты B относительно центра основания горы
#  - C_d = ex^2 + ey^2
#  - q0 = h_m - zB (для радиусной части)
#  - alpha = (r_m^2)/(h_m^2)
#  - r_m^2 (для быстрого XY-фильтра)
pre_mountains = []
for (xm, ym, hm, rm) in mountains:
    # Дегенеративные горы (нулевая высота или радиус) не могут ничего закрыть
    if hm <= EPS or rm <= EPS:
        pre_mountains.append(None)
        continue
    ex = xB - xm
    ey = yB - ym
    C_d = ex*ex + ey*ey
    q0 = hm - zB
    alpha = (rm*rm) / (hm*hm)
    rm2 = rm * rm
    pre_mountains.append((xm, ym, hm, rm, ex, ey, C_d, q0, alpha, rm2))

# ---------- Проверка: лодка в кадре камеры (FOV 90°, квадратный) ----------

def in_camera_fov(xp, yp, zp):
    """
    Проверка, попадает ли вершина лодки P = (xp, yp, zp) в кадр камеры:
    1) точка впереди камеры (z_c > 0);
    2) |x_c| <= z_c и |y_c| <= z_c (полуугол обзора 45° по X и Y).
    """
    vx = xp - xB
    vy = yp - yB
    vz = zp - zB

    # координаты в базисе камеры
    zc = vx*fx + vy*fy + vz*fz  # вперёд
    if zc <= EPS:
        # сзади или слишком близко в нуле
        return False

    xc = vx*rx + vy*ry + vz*rz  # вправо
    yc = vx*ux + vy*uy + vz*uz  # вверх

    # FOV 90° по обеим осям ⇒ |x_c| <= z_c, |y_c| <= z_c
    if abs(xc) > zc + 1e-12:
        return False
    if abs(yc) > zc + 1e-12:
        return False
    return True

# ---------- Проверка: гора заслоняет лодку? ----------

def mountain_blocks_boat(boat, pm):
    """
    Проверка, закрывает ли одна гора заданную лодку.
    boat = (xb, yb, hb)
    pm = предобработанные параметры горы или None.
    Возвращает True, если гора заслоняет лодку.
    """
    if pm is None:
        return False

    (xm, ym, hm, rm, ex, ey, C_d, q0, alpha, rm2) = pm

    xb, yb, hb = boat
    xp, yp, zp = xb, yb, hb

    # Вектор от наблюдателя B до вершины лодки P
    dx = xp - xB
    dy = yp - yB
    dz = zp - zB

    # ---------- Быстрый фильтр по проекции в плоскость XY ----------
    # Если проекция отрезка B–P в плоскости XY не попадает в круг радиуса r_m вокруг (xm, ym),
    # то и конус не может быть пересечён (горизонтальное расстояние всегда > r_m).
    if abs(dx) < EPS and abs(dy) < EPS:
        # Отрезок вертикальный, ближайшая точка к центру основания горы по XY — это B
        dist2_xy = (xB - xm)**2 + (yB - ym)**2
        if dist2_xy > rm2 + 1e-8:
            return False
    else:
        # общий случай: находим ближайшую к (xm, ym) точку на проекции отрезка B–P в XY
        vx = dx
        vy = dy
        # вектор от B_xy к центру горы
        wx = xm - xB
        wy = ym - yB
        v2 = vx*vx + vy*vy
        t_proj = (vx*wx + vy*wy) / v2
        if t_proj < 0.0:
            t_clamped = 0.0
        elif t_proj > 1.0:
            t_clamped = 1.0
        else:
            t_clamped = t_proj
        closest_x = xB + vx * t_clamped
        closest_y = yB + vy * t_clamped
        dxm = closest_x - xm
        dym = closest_y - ym
        dist2_xy = dxm*dxm + dym*dym
        if dist2_xy > rm2 + 1e-8:
            return False

    # ---------- Интервал по высоте: 0 <= z(t) <= h_m ----------

    # z(t) = zB + dz * t
    # Нужно найти t, при которых z(t) ∈ [0, hm].
    # Затем пересечь этот интервал с [0, 1].

    if abs(dz) < EPS:
        # z(t) постоянна
        z_const = zB
        if z_const < 0.0 - EPS or z_const > hm + EPS:
            # по высоте отрезок вообще не проходит через диапазон горы
            return False
        # по высоте допускается весь [0,1]
        t_z_min = 0.0
        t_z_max = 1.0
    else:
        # решаем неравенства:
        # 1) zB + dz t >= 0
        # 2) zB + dz t <= hm
        t1_min = -1e30
        t1_max =  1e30

        # z >= 0
        t_zero = -zB / dz
        if dz > 0:
            t1_min = max(t1_min, t_zero)
        else:
            t1_max = min(t1_max, t_zero)

        # z <= hm
        t_hm = (hm - zB) / dz
        if dz > 0:
            t1_max = min(t1_max, t_hm)
        else:
            t1_min = max(t1_min, t_hm)

        # Пересечение с [0,1]
        t_z_min = max(0.0, t1_min)
        t_z_max = min(1.0, t1_max)

        if t_z_min > t_z_max + EPS:
            # отрезок B–P не проходит через вертикальный диапазон горы
            return False

    # Теперь по высоте возможные t лежат в [t_z_min, t_z_max].
    # Нам ещё важно, чтобы t > 0 (иначе "перекрытие" только в точке B),
    # поэтому сдвигаем нижнюю границу:
    t0 = max(t_z_min, EPS_T)
    t1 = t_z_max
    if t0 > t1 + EPS:
        return False

    # ---------- Квадратичное неравенство по радиусу: g(t) <= 0 ----------

    # g(t) = X^2 + Y^2 - r_m^2 * (1 - Z/h_m)^2
    # X = ex + dx t, Y = ey + dy t, Z = zB + dz t
    # Мы заранее посчитали:
    #  A1, B1, C1 — для X^2 + Y^2
    #  A2, B2, C2 — для r_m^2 * (1 - Z/h_m)^2
    # g(t) = (A1-A2) t^2 + (B1-B2) t + (C1-C2)

    # Часть X^2 + Y^2
    A1 = dx*dx + dy*dy
    B1 = 2.0 * (ex*dx + ey*dy)
    C1 = C_d  # ex^2 + ey^2

    # Часть радиуса конуса
    # alpha = (r_m^2)/(h_m^2)
    # q0 = h_m - zB
    # q1 = -dz
    # A2 = alpha * q1^2 = alpha * dz^2
    # B2 = 2*alpha*q0*q1 = -2*alpha*q0*dz
    # C2 = alpha*q0^2
    A2 = alpha * dz*dz
    B2 = -2.0 * alpha * q0 * dz
    C2 = alpha * q0 * q0

    A = A1 - A2
    B = B1 - B2
    C = C1 - C2

    # Теперь нужно проверить, существует ли t ∈ [t0, t1], при котором g(t) <= 0.
    # g(t) — квадратичный многочлен.

    def g_val(t):
        return (A*t + B)*t + C  # A*t^2 + B*t + C

    if abs(A) < 1e-12:
        # Линейный (или почти линейный) случай: g(t) = B t + C
        # Минимум на [t0, t1] будет в одном из концов.
        g0 = g_val(t0)
        g1 = g_val(t1)
        g_min = g0 if g0 < g1 else g1
        return g_min <= EPS
    else:
        # Квадратичная парабола. Минимум либо в вершине, либо на границах.
        t_vertex = -B / (2.0 * A)
        # Список точек, в которых будем проверять g(t)
        candidates = [t0, t1]
        if t_vertex > t0 - 1e-12 and t_vertex < t1 + 1e-12:
            candidates.append(t_vertex)

        g_min = float('inf')
        for t in candidates:
            # ограничим t в пределах [t0, t1] на всякий случай
            if t < t0:
                t = t0
            elif t > t1:
                t = t1
            gv = g_val(t)
            if gv < g_min:
                g_min = gv

        return g_min <= EPS

# ---------- Основной цикл по лодкам ----------

visible_indices = []

for i, boat in enumerate(boats, start=1):
    xb, yb, hb = boat
    # Сначала проверяем, попадает ли вершина лодки в кадр по FOV
    if not in_camera_fov(xb, yb, hb):
        continue

    # Проверяем, не закрывает ли лодку хоть одна гора
    blocked = False
    for pm in pre_mountains:
        if pm is None:
            continue
        if mountain_blocks_boat(boat, pm):
            blocked = True
            break

    if not blocked:
        visible_indices.append(i)

# ---------- Вывод результата ----------

print(len(visible_indices))
for idx in visible_indices:
    print(idx)
