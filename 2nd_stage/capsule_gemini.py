import sys
import math

def vec_sub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def vec_add(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])

def vec_dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def vec_cross(a, b):
    return (a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0])

def vec_len(a):
    return math.sqrt(a[0]**2 + a[1]**2 + a[2]**2)

def vec_scale(a, s):
    return (a[0]*s, a[1]*s, a[2]*s)

def vec_normalize(a):
    l = vec_len(a)
    if l < 1e-12: return (0.0, 0.0, 0.0)
    return (a[0]/l, a[1]/l, a[2]/l)

def solve():
    input_data = sys.stdin.read().split()
    if not input_data:
        return
    
    iterator = iter(input_data)
    
    try:
        N = int(next(iterator))
    except StopIteration:
        return

    vertices = []
    for _ in range(N):
        x = float(next(iterator))
        y = float(next(iterator))
        z = float(next(iterator))
        vertices.append((x, y, z))

    M = int(next(iterator))
    faces = []
    for _ in range(M):
        k = int(next(iterator))
        face_idxs = []
        for _ in range(k):
            face_idxs.append(int(next(iterator)))
        faces.append(face_idxs)

    X0 = float(next(iterator))
    Y0 = float(next(iterator))
    Z0 = float(next(iterator))
    U = float(next(iterator))
    V = float(next(iterator))
    W = float(next(iterator))
    ta = float(next(iterator))
    tb = float(next(iterator))

    C1 = (X0, Y0, Z0)
    L_vec = (U, V, W)

    cx, cy, cz = 0.0, 0.0, 0.0
    for v in vertices:
        cx += v[0]
        cy += v[1]
        cz += v[2]
    centroid = (cx/N, cy/N, cz/N)

    planes = []
    
    for f_idxs in faces:
        p0 = vertices[f_idxs[0]]
        
        normal = (0, 0, 0)
        found_normal = False
        
        for i in range(len(f_idxs) - 2):
            p1 = vertices[f_idxs[i+1]]
            p2 = vertices[f_idxs[i+2]]
            edge1 = vec_sub(p1, p0)
            edge2 = vec_sub(p2, p1)
            n = vec_cross(edge1, edge2)
            if vec_len(n) > 1e-9:
                normal = vec_normalize(n)
                found_normal = True
                break
        
        if not found_normal:
            continue 

        vec_to_center = vec_sub(centroid, p0)
        if vec_dot(normal, vec_to_center) < 0:
            normal = vec_scale(normal, -1.0)
            
        d_val = -vec_dot(normal, p0)
        planes.append((normal, d_val))

    P_start = vec_add(C1, vec_scale(L_vec, ta))
    P_end = vec_add(C1, vec_scale(L_vec, tb))
    
    L_unit = vec_normalize(L_vec)
    L_len = vec_len(L_vec)

    R_max = float('inf')

    for normal, d in planes:
        dist_start = vec_dot(normal, P_start) + d
        dist_end = vec_dot(normal, P_end) + d
        
        min_dist = min(dist_start, dist_end)
        
        if min_dist < -1e-9:
            print("0")
            return

        if L_len > 1e-12:
            cos_theta = vec_dot(normal, L_unit)
            if cos_theta > 1.0: cos_theta = 1.0
            if cos_theta < -1.0: cos_theta = -1.0
            
            sin_theta = math.sqrt(1.0 - cos_theta**2)
        else:
            sin_theta = 1.0
        
        if sin_theta > 1e-9:
            r_limit = min_dist / sin_theta
            if r_limit < R_max:
                R_max = r_limit

    if R_max <= 1e-9:
        print("0")
    else:
        print(f"1 {R_max:.5f}")

if __name__ == '__main__':
    solve()