import sys
import math


def main():
    data = sys.stdin.read().strip().split()
    if not data:
        print(0)
        return

    it = iter(data)

    try:
        n = int(next(it))
    except StopIteration:
        print(0)
        return

    verts = []
    sx = sy = sz = 0.0
    for _ in range(n):
        x = float(next(it))
        y = float(next(it))
        z = float(next(it))
        verts.append((x, y, z))
        sx += x
        sy += y
        sz += z

    if n == 0:
        print(0)
        return

    cx = sx / n
    cy = sy / n
    cz = sz / n
    centroid = (cx, cy, cz)

    m = int(next(it))

    planes = []

    for _ in range(m):
        k = int(next(it))
        idx0 = int(next(it))
        idx1 = int(next(it))
        idx2 = int(next(it))

        v0 = verts[idx0]
        v1 = verts[idx1]
        v2 = verts[idx2]

        for __ in range(k - 3):
            _ = next(it)

        ax = v1[0] - v0[0]
        ay = v1[1] - v0[1]
        az = v1[2] - v0[2]
        bx = v2[0] - v0[0]
        by = v2[1] - v0[1]
        bz = v2[2] - v0[2]

        nx = ay * bz - az * by
        ny = az * bx - ax * bz
        nz = ax * by - ay * bx

        len2 = nx * nx + ny * ny + nz * nz
        if len2 < 1e-18:
            continue  

        d = -(nx * v0[0] + ny * v0[1] + nz * v0[2])

        s = nx * centroid[0] + ny * centroid[1] + nz * centroid[2] + d
        if s > 0.0:
            nx = -nx
            ny = -ny
            nz = -nz
            d = -d

        length = math.sqrt(nx * nx + ny * ny + nz * nz)
        nx /= length
        ny /= length
        nz /= length
        d /= length

        planes.append((nx, ny, nz, d))

    X0 = float(next(it))
    Y0 = float(next(it))
    Z0 = float(next(it))
    U = float(next(it))
    V = float(next(it))
    W = float(next(it))
    ta = float(next(it))
    tb = float(next(it))

    if ta > tb:
        ta, tb = tb, ta

    L2 = U * U + V * V + W * W
    if L2 < 1e-18:
        print(0)
        return

    L_norm = math.sqrt(L2)
    Lhx = U / L_norm
    Lhy = V / L_norm
    Lhz = W / L_norm

    Pa_x = X0 + U * ta
    Pa_y = Y0 + V * ta
    Pa_z = Z0 + W * ta

    Pb_x = X0 + U * tb
    Pb_y = Y0 + V * tb
    Pb_z = Z0 + W * tb

    EPS_AXIS = 1e-9
    EPS_PARALLEL = 1e-12
    EPS_R = 1e-9

    for nx, ny, nz, d in planes:
        va = nx * Pa_x + ny * Pa_y + nz * Pa_z + d
        vb = nx * Pb_x + ny * Pb_y + nz * Pb_z + d
        worst = va if va > vb else vb
        if worst > EPS_AXIS:
            print(0)
            return

    R_max = float('inf')

    for nx, ny, nz, d in planes:
        va = nx * Pa_x + ny * Pa_y + nz * Pa_z + d
        vb = nx * Pb_x + ny * Pb_y + nz * Pb_z + d
        worst = va if va > vb else vb  

        ndl = nx * Lhx + ny * Lhy + nz * Lhz
        n_perp_sq = 1.0 - ndl * ndl

        if n_perp_sq <= EPS_PARALLEL:
            continue

        if n_perp_sq < 0.0:
            n_perp_sq = 0.0
        n_perp = math.sqrt(n_perp_sq)

        gap = -worst
        if gap < 0.0:
            gap = 0.0

        R_limit = gap / n_perp
        if R_limit < R_max:
            R_max = R_limit

    if not math.isfinite(R_max):
        print(0)
        return

    if R_max <= EPS_R:
        print(0)
    else:
        print("1 {:.5f}".format(R_max))


if __name__ == "__main__":
    main()
