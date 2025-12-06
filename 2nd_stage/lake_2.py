import sys
import math

EPS = 1e-9
EPS_T = 1e-7

# ---------- Ввод ----------

data = sys.stdin.read().strip().split()
if not data:
    # Пустой ввод – просто ничего
    print("0.00000")
    print(0)
    sys.exit(0)

it = iter(data)

def next_float():
    return float(next(it))

def next_int():
    return int(next(it))

# Точка A
xA = next_float()
yA = next_float()
zA = next_float()

# Точка B
xB_end = next_float()
yB_end = next_float()
zB_end = next_float()

# Скорость
v = next_float()

# Направление камеры l
lx = next_float()
ly = next_float()
lz = next_float()

# N лодок, M гор
N = next_int()
M = next_int()

boats = []
for _ in range(N):
    xb = next_float()
    yb = next_float()
    hb = next_float()
    boats.append((xb, yb, hb))

mountains = []
for _ in range(M):
    xm = next_float()
    ym = next_float()
    hm = next_float()
    rm = next_float()
    mountains.append((xm, ym, hm, rm))

# ---------- Базис камеры (f, r, u) ----------

len_l = math.sqrt(lx*lx + ly*ly + lz*lz)
if len_l < EPS:
    # камера не определена – считаем, что ничего не видно
    print("0.00000")
    print(0)
    sys.exit(0)

fx = lx / len_l
fy = ly / len_l
fz = lz / len_l

# r – горизонтальная ось (в плоскости z=0)
if abs(lx) + abs(ly) > EPS:
    r0x = -ly
    r0y = lx
    r0z = 0.0
else:
    r0x = 1.0
    r0y = 0.0
    r0z = 0.0

len_r0 = math.sqrt(r0x*r0x + r0y*r0y + r0z*r0z)
rx = r0x / len_r0
ry = r0y / len_r0
rz = r0z / len_r0  # будет 0

# u – вертикальная ось
ux = fy*rz - fz*ry
uy = fz*rx - fx*rz
uz = fx*ry - fy*rx

# ---------- Предобработка гор (часть, не зависящая от камеры) ----------

# Для каждой горы храним: (xm, ym, hm, rm, alpha, rm2)
pre_mountains = []
for (xm, ym, hm, rm) in mountains:
    if hm <= EPS or rm <= EPS:
        pre_mountains.append(None)
        continue
    alpha = (rm * rm) / (hm * hm)
    rm2 = rm * rm
    pre_mountains.append((xm, ym, hm, rm, alpha, rm2))

# ---------- Вспомогательные функции для статической проверки ----------

def in_camera_fov_at(Cx, Cy, Cz, xb, yb, hb):
    """Лодка (xb,yb,hb) в кадре камеры, находящейся в C=(Cx,Cy,Cz)?"""
    dx = xb - Cx
    dy = yb - Cy
    dz = hb - Cz

    zc = dx * fx + dy * fy + dz * fz
    if zc <= EPS:
        return False

    xc = dx * rx + dy * ry + dz * rz
    yc = dx * ux + dy * uy + dz * uz

    if abs(xc) > zc + 1e-12:
        return False
    if abs(yc) > zc + 1e-12:
        return False
    return True


def mountain_blocks_boat_at(Cx, Cy, Cz, boat, pm):
    """Гора (предобработанная как pm) заслоняет лодку из точки C?"""
    if pm is None:
        return False

    xb, yb, hb = boat
    xm, ym, hm, rm, alpha, rm2 = pm

    # Вектор от камеры к лодке
    dx = xb - Cx
    dy = yb - Cy
    dz = hb - Cz

    # Быстрый XY-фильтр
    if abs(dx) < EPS and abs(dy) < EPS:
        # почти вертикальный луч — ближайшая по XY точка это C
        dist2_xy = (Cx - xm) * (Cx - xm) + (Cy - ym) * (Cy - ym)
        if dist2_xy > rm2 + 1e-8:
            return False
    else:
        vx = dx
        vy = dy
        wx = xm - Cx
        wy = ym - Cy
        v2 = vx*vx + vy*vy
        t_proj = (vx*wx + vy*wy) / v2
        if t_proj < 0.0:
            t_clamped = 0.0
        elif t_proj > 1.0:
            t_clamped = 1.0
        else:
            t_clamped = t_proj
        closest_x = Cx + vx * t_clamped
        closest_y = Cy + vy * t_clamped
        dxm = closest_x - xm
        dym = closest_y - ym
        dist2_xy = dxm*dxm + dym*dym
        if dist2_xy > rm2 + 1e-8:
            return False

    # Интервал по высоте: 0 <= z(t) <= h_m
    zC = Cz
    if abs(dz) < EPS:
        # z(t) постоянна
        if zC < -EPS or zC > hm + EPS:
            return False
        t_z_min = 0.0
        t_z_max = 1.0
    else:
        t1_min = -1e30
        t1_max = 1e30

        # z >= 0
        t_zero = -zC / dz
        if dz > 0:
            if t_zero > t1_min:
                t1_min = t_zero
        else:
            if t_zero < t1_max:
                t1_max = t_zero

        # z <= hm
        t_hm = (hm - zC) / dz
        if dz > 0:
            if t_hm < t1_max:
                t1_max = t_hm
        else:
            if t_hm > t1_min:
                t1_min = t_hm

        t_z_min = max(0.0, t1_min)
        t_z_max = min(1.0, t1_max)

        if t_z_min > t_z_max + EPS:
            return False

    t0 = max(t_z_min, EPS_T)
    t1 = t_z_max
    if t0 > t1 + EPS:
        return False

    # Радиусное неравенство g(t) <= 0
    ex = Cx - xm
    ey = Cy - ym
    C_d = ex*ex + ey*ey
    q0 = hm - Cz

    # X^2 + Y^2
    A1 = dx*dx + dy*dy
    B1 = 2.0 * (ex*dx + ey*dy)
    C1 = C_d

    # r_m^2 * (1 - Z/h_m)^2
    A2 = alpha * dz * dz
    B2 = -2.0 * alpha * q0 * dz
    C2 = alpha * q0 * q0

    A = A1 - A2
    B = B1 - B2
    Ccoef = C1 - C2

    def g_val(t):
        return (A * t + B) * t + Ccoef

    if abs(A) < 1e-12:
        g0 = g_val(t0)
        g1 = g_val(t1)
        g_min = g0 if g0 < g1 else g1
        return g_min <= EPS
    else:
        t_vertex = -B / (2.0 * A)
        candidates = [t0, t1]
        if t_vertex > t0 - 1e-12 and t_vertex < t1 + 1e-12:
            candidates.append(t_vertex)

        g_min = float('inf')
        for t in candidates:
            if t < t0:
                tt = t0
            elif t > t1:
                tt = t1
            else:
                tt = t
            gv = g_val(tt)
            if gv < g_min:
                g_min = gv

        return g_min <= EPS


# ---------- Функция видимости для данного s (0..1) ----------

# Вектор AB и длина
ABx = xB_end - xA
ABy = yB_end - yA
ABz = zB_end - zA
AB_len = math.sqrt(ABx*ABx + ABy*ABy + ABz*ABz)

# Если A и B совпадают или скорость нулевая – камера не движется
if AB_len < EPS or v <= 0:
    # Просто решаем статическую задачу из точки A
    Cx0 = xA
    Cy0 = yA
    Cz0 = zA

    visible_indices = []
    for i, boat in enumerate(boats, start=1):
        xb, yb, hb = boat
        if not in_camera_fov_at(Cx0, Cy0, Cz0, xb, yb, hb):
            continue
        blocked = False
        for pm in pre_mountains:
            if pm is None:
                continue
            if mountain_blocks_boat_at(Cx0, Cy0, Cz0, boat, pm):
                blocked = True
                break
        if not blocked:
            visible_indices.append(i)

    print("0.00000")
    print(len(visible_indices))
    for idx in visible_indices:
        print(idx)
    sys.exit(0)

T_total = AB_len / v  # полное время движения

def camera_pos_s(s):
    """Положение камеры для параметра s в [0,1]."""
    return (xA + ABx * s,
            yA + ABy * s,
            zA + ABz * s)

def visible_at_s(s):
    """Возвращает (count, [indices]) для параметра s."""
    Cx, Cy, Cz = camera_pos_s(s)
    visible = []
    for i, boat in enumerate(boats, start=1):
        xb, yb, hb = boat
        if not in_camera_fov_at(Cx, Cy, Cz, xb, yb, hb):
            continue
        blocked = False
        for pm in pre_mountains:
            if pm is None:
                continue
            if mountain_blocks_boat_at(Cx, Cy, Cz, boat, pm):
                blocked = True
                break
        if not blocked:
            visible.append(i)
    return len(visible), visible

# ---------- Адаптивный поиск по s ∈ [0,1] ----------

# Ограничения на рекурсию и число вычислений
MAX_DEPTH = 16         # глубина деления отрезка (2^16 ~ 65k подотрезков в худшем случае)
MAX_EVALS = 1500       # максимум вызовов visible_at_s
EVALS = 0

best_cnt = -1
best_s = 0.0
best_list = []

# Кеш значений, чтобы не пересчитывать одинаковые s
cache = {}

def get_visible_cached(s):
    global EVALS
    if s in cache:
        return cache[s]
    if EVALS >= MAX_EVALS:
        # Дальше не считаем, возвращаем что-нибудь
        return 0, []
    cnt, lst = visible_at_s(s)
    cache[s] = (cnt, lst)
    EVALS += 1
    return cnt, lst

def update_best(s, cnt, lst):
    global best_cnt, best_s, best_list
    if cnt > best_cnt or (cnt == best_cnt and s < best_s):
        best_cnt = cnt
        best_s = s
        best_list = lst

def dfs(l, fl, flist, r, fr, frist, depth):
    if depth >= MAX_DEPTH or EVALS >= MAX_EVALS:
        # обновляем по концам и выходим
        update_best(l, fl, flist)
        update_best(r, fr, frist)
        return

    # оценка: если максимум на концах уже ниже глобального – см. ниже
    local_max = fl if fl > fr else fr
    if local_max < best_cnt:
        # на этом отрезке не улучшим
        update_best(l, fl, flist)
        update_best(r, fr, frist)
        return

    m = 0.5 * (l + r)
    fm, fmlist = get_visible_cached(m)
    update_best(m, fm, fmlist)

    if max(fl, fm, fr) < best_cnt:
        # нет смысла делить дальше
        return

    dfs(l, fl, flist, m, fm, fmlist, depth + 1)
    if EVALS >= MAX_EVALS:
        return
    dfs(m, fm, fmlist, r, fr, frist, depth + 1)

# стартовые точки s=0 и s=1
f0, l0 = get_visible_cached(0.0)
f1, l1 = get_visible_cached(1.0)
update_best(0.0, f0, l0)
update_best(1.0, f1, l1)

dfs(0.0, f0, l0, 1.0, f1, l1, 0)

# ---------- Вывод ----------

t_ans = best_s * T_total
print(f"{t_ans:.5f}")
print(len(best_list))
for idx in best_list:
    print(idx)
