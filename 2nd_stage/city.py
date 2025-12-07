import sys
import math


def point_in_view(x, y, z, px, py, pz, dx, dy, dz, len_d, cos_alpha):
    vx = x - px
    vy = y - py
    vz = z - pz

    # Точка в вершине конуса
    if vx == 0 and vy == 0 and vz == 0:
        return True

    dot = vx * dx + vy * dy + vz * dz
    # За спиной
    if dot <= 0:
        return False

    len_v2 = vx * vx + vy * vy + vz * vz
    len_v = math.sqrt(len_v2)
    # cos(theta) = (v·d) / (|v||d|)
    cos_theta = dot / (len_v * len_d)

    # Небольшой запас на погрешности double
    return cos_theta + 1e-12 >= cos_alpha


def main():
    data = sys.stdin.readline().split()
    if not data:
        return

    px, py, pz, dx, dy, dz, alpha = map(int, data)

    n_line = sys.stdin.readline().split()
    while not n_line:
        n_line = sys.stdin.readline().split()
    N = int(n_line[0])

    # Подготовка геометрии
    len_d2 = dx * dx + dy * dy + dz * dz
    # На случай странного ввода с нулевым направлением — тогда ничего не видно.
    if len_d2 == 0:
        # Всё равно надо прочитать вход до конца
        for _ in range(N):
            sys.stdin.readline()
        print(0)
        return

    len_d = math.sqrt(len_d2)
    cos_alpha = math.cos(math.radians(alpha))

    cell_map = {}  # (x,y,z) -> set(building_ids)

    for _ in range(N):
        # Строка модуля: i x1 y1 z1 ... x4 y4 z4
        parts = sys.stdin.readline().split()
        # На случай пустых строк
        while parts and len(parts) < 13:
            more = sys.stdin.readline().split()
            parts += more
        if not parts:
            continue
        vals = list(map(int, parts))
        building_id = vals[0]
        coords = vals[1:]

        # 4 секции
        for i in range(0, 12, 3):
            x = coords[i]
            y = coords[i + 1]
            z = coords[i + 2]

            if point_in_view(x, y, z, px, py, pz, dx, dy, dz, len_d, cos_alpha):
                key = (x, y, z)
                s = cell_map.get(key)
                if s is None:
                    cell_map[key] = {building_id}
                else:
                    s.add(building_id)

    conflicts = []

    for (x, y, z), buildings in cell_map.items():
        if len(buildings) >= 2:
            ids = sorted(buildings)
            conflicts.append((x, y, z, ids))

    conflicts.sort(key=lambda item: (item[0], item[1], item[2]))

    out_lines = []
    out_lines.append(str(len(conflicts)))
    for x, y, z, ids in conflicts:
        line = [str(x), str(y), str(z)] + [str(a) for a in ids]
        out_lines.append(" ".join(line))

    sys.stdout.write("\n".join(out_lines))


if __name__ == "__main__":
    main()
