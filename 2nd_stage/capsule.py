import sys
import math
from typing import List, Tuple, Iterator, Optional


# Константы для численных вычислений
EPS_AXIS = 1e-9       # допуск для проверки оси
EPS_R_ZERO = 1e-7     # порог, меньше которого радиус считаем нулевым
EPS_NORMAL_SQ = 1e-18 # порог для проверки нормали плоскости
EPS_NEGATIVE = 1e-12  # порог для защиты от отрицательных нулей
SAFETY_EPS = 1e-9     # запас (отступ) от плоскостей для радиуса > 0
BINARY_SEARCH_ITERATIONS = 60  # количество итераций бинарного поиска


def read_vertices(iterator: Iterator[str], n: int) -> List[Tuple[float, float, float]]:
    """Читает n вершин многогранника и возвращает список координат."""
    vertices = []
    for _ in range(n):
        x = float(next(iterator))
        y = float(next(iterator))
        z = float(next(iterator))
        vertices.append((x, y, z))
    return vertices


def compute_centroid(vertices: List[Tuple[float, float, float]]) -> Tuple[float, float, float]:
    """Вычисляет центроид (центр масс) многогранника."""
    if not vertices:
        return (0.0, 0.0, 0.0)
    
    n = len(vertices)
    sum_x = sum(v[0] for v in vertices)
    sum_y = sum(v[1] for v in vertices)
    sum_z = sum(v[2] for v in vertices)
    return (sum_x / n, sum_y / n, sum_z / n)


def cross_product(ax: float, ay: float, az: float, 
                  bx: float, by: float, bz: float) -> Tuple[float, float, float]:
    """Вычисляет векторное произведение двух 3D-векторов."""
    nx = ay * bz - az * by
    ny = az * bx - ax * bz
    nz = ax * by - ay * bx
    return (nx, ny, nz)


def build_plane_from_face(vertices: List[Tuple[float, float, float]], 
                         face_indices: List[int],
                         centroid: Tuple[float, float, float]) -> Optional[Tuple[float, float, float, float]]:
    """
    Строит уравнение плоскости из грани многогранника.
    Возвращает (nx, ny, nz, d) для плоскости n·x + d = 0, или None если грань вырождена.
    Нормаль ориентирована внутрь многогранника.
    """
    if len(face_indices) < 3:
        return None
    
    v0 = vertices[face_indices[0]]
    v1 = vertices[face_indices[1]]
    v2 = vertices[face_indices[2]]
    
    # Векторы в плоскости грани
    ax = v1[0] - v0[0]
    ay = v1[1] - v0[1]
    az = v1[2] - v0[2]
    bx = v2[0] - v0[0]
    by = v2[1] - v0[1]
    bz = v2[2] - v0[2]
    
    # Нормаль к плоскости
    nx, ny, nz = cross_product(ax, ay, az, bx, by, bz)
    
    len_sq = nx * nx + ny * ny + nz * nz
    if len_sq < EPS_NORMAL_SQ:
        return None
    
    # Плоскость: n·x + d = 0 через v0
    d = -(nx * v0[0] + ny * v0[1] + nz * v0[2])
    
    # Ориентация нормали внутрь многогранника (через центроид)
    centroid_value = nx * centroid[0] + ny * centroid[1] + nz * centroid[2] + d
    if centroid_value > 0.0:
        nx, ny, nz, d = -nx, -ny, -nz, -d
    
    # Нормализация
    length = math.sqrt(len_sq)
    return (nx / length, ny / length, nz / length, d / length)


def read_planes(iterator: Iterator[str], vertices: List[Tuple[float, float, float]],
                centroid: Tuple[float, float, float]) -> List[Tuple[float, float, float, float]]:
    """Читает грани многогранника и строит уравнения плоскостей."""
    m = int(next(iterator))
    planes = []
    
    for _ in range(m):
        k = int(next(iterator))
        if k < 3:
            # Пропускаем вырожденные грани
            for _ in range(k):
                next(iterator)
            continue
        
        # Читаем индексы вершин грани
        face_indices = [int(next(iterator)) for _ in range(k)]
        
        plane = build_plane_from_face(vertices, face_indices, centroid)
        if plane is not None:
            planes.append(plane)
    
    return planes


def read_cylinder_params(iterator: Iterator[str]) -> Optional[Tuple[float, float, float, 
                                                                    float, float, float, 
                                                                    float, float]]:
    """Читает параметры цилиндра. Возвращает None при ошибке."""
    try:
        X0 = float(next(iterator))
        Y0 = float(next(iterator))
        Z0 = float(next(iterator))
        U = float(next(iterator))
        V = float(next(iterator))
        W = float(next(iterator))
        ta = float(next(iterator))
        tb = float(next(iterator))
        
        if ta > tb:
            ta, tb = tb, ta
        
        return (X0, Y0, Z0, U, V, W, ta, tb)
    except StopIteration:
        return None


def normalize_vector(x: float, y: float, z: float) -> Optional[Tuple[float, float, float]]:
    """Нормализует вектор. Возвращает None если вектор слишком короткий."""
    len_sq = x * x + y * y + z * z
    if len_sq < EPS_NORMAL_SQ:
        return None
    
    length = math.sqrt(len_sq)
    return (x / length, y / length, z / length)


def compute_axis_endpoints(X0: float, Y0: float, Z0: float,
                          U: float, V: float, W: float,
                          ta: float, tb: float) -> Tuple[Tuple[float, float, float],
                                                         Tuple[float, float, float]]:
    """Вычисляет концы оси цилиндра."""
    Pa = (X0 + U * ta, Y0 + V * ta, Z0 + W * ta)
    Pb = (X0 + U * tb, Y0 + V * tb, Z0 + W * tb)
    return (Pa, Pb)


def compute_plane_constraints(planes: List[Tuple[float, float, float, float]],
                             axis_unit: Tuple[float, float, float],
                             axis_endpoints: Tuple[Tuple[float, float, float], Tuple[float, float, float]]
                             ) -> Optional[Tuple[List[float], List[float], float]]:
    """
    Вычисляет ограничения на радиус для каждой плоскости.
    Возвращает (gmax_list, nperp_list, R_hi) или None если ось не внутри многогранника.
    """
    Pa, Pb = axis_endpoints
    gmax_list = []
    nperp_list = []
    R_hi = float('inf')
    
    axis_unit_x, axis_unit_y, axis_unit_z = axis_unit
    
    for nx, ny, nz, d in planes:
        # Значения плоскости в концах оси
        va = nx * Pa[0] + ny * Pa[1] + nz * Pa[2] + d
        vb = nx * Pb[0] + ny * Pb[1] + nz * Pb[2] + d
        gmax = max(va, vb)
        
        # Проверка: ось целиком внутри многогранника
        if gmax > EPS_AXIS:
            return None
        
        gmax_list.append(gmax)
        
        # Компонента нормали, перпендикулярная оси
        ndotL = nx * axis_unit_x + ny * axis_unit_y + nz * axis_unit_z
        n_perp_sq = 1.0 - ndotL * ndotL
        
        # Защита от отрицательных нулей из-за округления
        if n_perp_sq < 0.0:
            if n_perp_sq > -EPS_NEGATIVE:
                n_perp_sq = 0.0
            else:
                n_perp_sq = 0.0
        
        n_perp = math.sqrt(n_perp_sq)
        nperp_list.append(n_perp)
        
        # Аналитическая верхняя граница радиуса по этой плоскости
        if n_perp > 0.0:
            gap = -gmax  # запас до плоскости вдоль нормали на оси
            if gap <= 0.0:
                # Ось уже на границе или вне -> радиуса запаса нет
                R_limit = 0.0
            else:
                R_limit = gap / n_perp
            if R_limit < R_hi:
                R_hi = R_limit
    
    return (gmax_list, nperp_list, R_hi)


def is_radius_valid(R: float, 
                   planes: List[Tuple[float, float, float, float]],
                   gmax_list: List[float],
                   nperp_list: List[float]) -> bool:
    """
    Проверяет, что цилиндр с данным радиусом R полностью внутри многогранника
    с учетом безопасного отступа SAFETY_EPS.
    """
    for i, (nx, ny, nz, d) in enumerate(planes):
        gmax = gmax_list[i]
        n_perp = nperp_list[i]
        
        if n_perp == 0.0:
            # Нормаль плоскости параллельна оси цилиндра
            if gmax > EPS_AXIS:
                return False
            continue
        
        # Максимальное значение плоскости по всем точкам цилиндра
        value = gmax + R * n_perp
        
        # Требуем, чтобы цилиндр был с небольшим отрицательным запасом внутри
        if value > -SAFETY_EPS:
            return False
    
    return True


def binary_search_max_radius(planes: List[Tuple[float, float, float, float]],
                            gmax_list: List[float],
                            nperp_list: List[float],
                            R_hi: float) -> float:
    """Бинарный поиск максимального радиуса цилиндра."""
    lo = 0.0
    hi = R_hi
    
    for _ in range(BINARY_SEARCH_ITERATIONS):
        mid = 0.5 * (lo + hi)
        if mid <= 0.0:
            lo = mid
            continue
        
        if is_radius_valid(mid, planes, gmax_list, nperp_list):
            lo = mid
        else:
            hi = mid
    
    return lo


def main():
    """Основная функция решения задачи."""
    data = sys.stdin.read().strip().split()
    if not data:
        print(0)
        return
    
    iterator = iter(data)
    
    # Чтение вершин
    try:
        n = int(next(iterator))
    except StopIteration:
        print(0)
        return
    
    if n == 0:
        print(0)
        return
    
    vertices = read_vertices(iterator, n)
    centroid = compute_centroid(vertices)
    
    # Чтение граней и построение плоскостей
    planes = read_planes(iterator, vertices, centroid)
    if not planes:
        print(0)
        return
    
    # Чтение параметров цилиндра
    cylinder_params = read_cylinder_params(iterator)
    if cylinder_params is None:
        print(0)
        return
    
    X0, Y0, Z0, U, V, W, ta, tb = cylinder_params
    
    # Нормализация направления оси
    axis_unit = normalize_vector(U, V, W)
    if axis_unit is None:
        print(0)
        return
    
    # Вычисление концов оси
    axis_endpoints = compute_axis_endpoints(X0, Y0, Z0, U, V, W, ta, tb)
    
    # Вычисление ограничений на радиус
    constraints = compute_plane_constraints(planes, axis_unit, axis_endpoints)
    if constraints is None:
        print(0)
        return
    
    gmax_list, nperp_list, R_hi = constraints
    
    # Проверка наличия положительной аналитической границы
    if not math.isfinite(R_hi) or R_hi <= EPS_R_ZERO:
        print(0)
        return
    
    # Бинарный поиск максимального радиуса
    best_R = binary_search_max_radius(planes, gmax_list, nperp_list, R_hi)
    
    # Вывод результата
    if best_R <= EPS_R_ZERO:
        print(0)
    else:
        print(f"1 {best_R:.5f}")


if __name__ == "__main__":
    main()
