import sys
import math

EPS = 10**-9


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def read_input():
    data = sys.stdin.read().strip().split()
    if not data:
        return None
    it = iter(data)

    n = int(next(it))
    verts = []
    for _ in range(n):
        x = float(next(it))
        y = float(next(it))
        z = float(next(it))
        verts.append((x, y, z))

    m = int(next(it))
    faces = []
    for _ in range(m):
        k = int(next(it))
        idx = [int(next(it)) for __ in range(k)]
        faces.append(idx)

    X0 = float(next(it))
    Y0 = float(next(it))
    Z0 = float(next(it))
    U = float(next(it))
    V = float(next(it))
    W = float(next(it))
    ta = float(next(it))
    tb = float(next(it))

    return verts, faces, (X0, Y0, Z0), (U, V, W), ta, tb


def build_planes(verts, faces):
    n_verts = len(verts)
    cx = sum(p[0] for p in verts) / n_verts
    cy = sum(p[1] for p in verts) / n_verts
    cz = sum(p[2] for p in verts) / n_verts
    centroid = (cx, cy, cz)

    planes = []

    for idx in faces:
        p0 = verts[idx[0]]
        normal = None

        for j in range(1, len(idx) - 1):
            p1 = verts[idx[j]]
            p2 = verts[idx[j + 1]]
            e1 = sub(p1, p0)
            e2 = sub(p2, p0)
            nvec = cross(e1, e2)
            if dot(nvec, nvec) > EPS:
                normal = nvec
                break

        length = math.sqrt(dot(normal, normal))
        nu = (normal[0] / length, normal[1] / length, normal[2] / length)
        h = dot(nu, p0)

        if dot(nu, centroid) - h > 0:
            nu = (-nu[0], -nu[1], -nu[2])
            h = -h

        planes.append((nu, h))

    return planes


def solve():
    parsed = read_input()
    if parsed is None:
        return

    verts, faces, C1, L, ta, tb = parsed
    if ta > tb:
        ta, tb = tb, ta

    planes = build_planes(verts, faces)

    Lx, Ly, Lz = L
    Llen = math.sqrt(Lx * Lx + Ly * Ly + Lz * Lz)
    if Llen < EPS:
        print(0)
        return

    L_hat = (Lx / Llen, Ly / Llen, Lz / Llen)

    def ok(R):
        low = -(10**60)
        high = 10**60

        for n, h in planes:
            a = dot(n, L)
            b = dot(n, C1)
            proj = dot(n, L_hat)
            k_sq = 1.0 - proj * proj
            if k_sq < 0.0:
                k_sq = 0.0
            k = math.sqrt(k_sq)

            rhs = h - b - R * k

            if abs(a) < EPS:
                if rhs < -1e-9:
                    return False
            elif a > 0:
                bound = rhs / a
                if bound < high:
                    high = bound
            else:
                bound = rhs / a
                if bound > low:
                    low = bound

            if low > high + 1e-9:
                return False

        low = max(low, ta)
        high = min(high, tb)
        return low <= high + 1e-9

    if not ok(0.0):
        print(0)
        return

    hi = float("inf")
    for n, h in planes:
        proj = dot(n, L_hat)
        k_sq = 1.0 - proj * proj
        if k_sq < 0.0:
            k_sq = 0.0
        k = math.sqrt(k_sq)
        if k < EPS:
            continue

        a = dot(n, L)
        b = dot(n, C1)

        if abs(a) < EPS:
            cand = (h - b) / k
        elif a > 0:
            cand = (h - b - a * ta) / k
        else:
            cand = (h - b - a * tb) / k

        if cand < hi:
            hi = cand

    if not math.isfinite(hi) or hi <= 0.0:
        print(0)
        return

    lo = 0.0
    for _ in range(70):
        mid = (lo + hi) * 0.5
        if ok(mid):
            lo = mid
        else:
            hi = mid

    if lo <= 1e-7:
        print(0)
        return

    print("1 {:.5f}".format(lo))


if __name__ == "__main__":
    solve()
