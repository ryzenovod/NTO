import math
import sys

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def norm(v):
    return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def solve_quadratic(A, B, C):
    if abs(A) < 1e-12:
        if abs(B) < 1e-12:
            return [] if abs(C) > 1e-12 else None  # all t or none
        return [-C / B]
    D = B*B - 4*A*C
    if D < 0:
        return []
    sqrtD = math.sqrt(D)
    return [(-B - sqrtD) / (2*A), (-B + sqrtD) / (2*A)]

def intersects_cone(B, P, xm, ym, hm, rm):
    xB, yB, zB = B
    xb, yb, hb = P
    dx = xb - xB
    dy = yb - yB
    dz = hb - zB

    ox = xB - xm
    oy = yB - ym
    oz = zB

    # Коэффициенты f(t) = A t^2 + B t + C
    A = dx*dx + dy*dy - (rm*rm * dz*dz) / (hm*hm)
    B = 2*(ox*dx + oy*dy) + 2*(rm*rm * (hm - oz) * dz) / (hm*hm)
    C = ox*ox + oy*oy - (rm*rm * (hm - oz)*(hm - oz)) / (hm*hm)

    # Решаем f(t) <= 0
    if abs(A) < 1e-12:
        if abs(B) < 1e-12:
            # f(t) = C
            if C > 0:
                return False
            # Проверяем, есть ли t in (0,1) с z in [0, hm]
            z0 = zB
            z1 = zB + dz
            z_min = min(z0, z1)
            z_max = max(z0, z1)
            return not (z_max < 0 or z_min > hm)
        else:
            t0 = -C / B
            if t0 <= 0 or t0 >= 1:
                return False
            z_t = zB + t0 * dz
            return 0 <= z_t <= hm
    else:
        roots = solve_quadratic(A, B, C)
        if not roots:
            # f(t) всегда одного знака
            # Проверим знак в t=0.5
            t_test = 0.5
            f_test = A*t_test*t_test + B*t_test + C
            if f_test > 0:
                return False
            # Всё f(t) <= 0 — проверяем по высоте
            z0 = zB
            z1 = zB + dz
            z_min = min(z0, z1)
            z_max = max(z0, z1)
            return not (z_max < 0 or z_min > hm)
        elif len(roots) == 1:
            t1 = t2 = roots[0]
        else:
            t1, t2 = min(roots), max(roots)

        # f(t) <= 0 на [t1, t2] если A > 0
        if A < 0:
            # f(t) <= 0 вне [t1, t2]
            # Но это невозможно для конуса при движении от внешней точки
            # Поэтому игнорируем
            return False

        t_low = max(t1, 0.0)
        t_high = min(t2, 1.0)
        if t_low >= t_high:
            return False

        z_low = zB + t_low * dz
        z_high = zB + t_high * dz
        z_min = min(z_low, z_high)
        z_max = max(z_low, z_high)

        return not (z_max < 0 or z_min > hm)

def main():
    data = sys.stdin.read().strip().split()
    if not data:
        return

    it = iter(data)
    xB = float(next(it)); yB = float(next(it)); zB = float(next(it))
    xl = float(next(it)); yl = float(next(it)); zl = float(next(it))
    N = int(next(it)); M = int(next(it))

    boats = []
    for _ in range(N):
        xb = float(next(it)); yb = float(next(it)); hb = float(next(it))
        boats.append((xb, yb, hb))

    mountains = []
    for _ in range(M):
        xm = float(next(it)); ym = float(next(it)); hm = float(next(it)); rm = float(next(it))
        mountains.append((xm, ym, hm, rm))

    # Нормируем направление камеры
    l_vec = (xl, yl, zl)
    l_len = norm(l_vec)
    if l_len == 0:
        # Некорректное направление — но по условию такого не будет
        l_norm = (0,0,1)
    else:
        l_norm = (xl/l_len, yl/l_len, zl/l_len)

    B = (xB, yB, zB)
    visible = []

    for idx, boat in enumerate(boats):
        xb, yb, hb = boat
        P = (xb, yb, hb)
        v = (xb - xB, yb - yB, hb - zB)
        d = dot(v, l_norm)
        if d <= 0:
            continue
        v_perp = (v[0] - d*l_norm[0], v[1] - d*l_norm[1], v[2] - d*l_norm[2])
        r = norm(v_perp)
        if r > d:
            continue

        # Проверка на горы
        blocked = False
        for (xm, ym, hm, rm) in mountains:
            if hm <= 0 or rm <= 0:
                continue
            if intersects_cone(B, P, xm, ym, hm, rm):
                blocked = True
                break
    if not blocked:
            visible.append(idx + 1)

    print(len(visible))
    for idx in visible:
        print(idx)

if __name__ == "main":
    main()