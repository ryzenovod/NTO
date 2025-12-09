import sys
import math

def sub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def add(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])

def mul(a, s):
    return (a[0]*s, a[1]*s, a[2]*s)

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def cross(a, b):
    return (a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0])

def norm_sq(a):
    return a[0]**2 + a[1]**2 + a[2]**2

def normalize(a):
    l = math.sqrt(a[0]**2 + a[1]**2 + a[2]**2)
    if l == 0: return (0,0,0)
    return (a[0]/l, a[1]/l, a[2]/l)

def solve_quadratic(a, b, c):
    if abs(a) < 1e-9:
        if abs(b) < 1e-9:
            return []
        return [-c/b]
    d = b*b - 4*a*c
    if d < 0:
        return []
    sqrt_d = math.sqrt(d)
    return sorted([(-b - sqrt_d)/(2*a), (-b + sqrt_d)/(2*a)])

def solve():
    input_data = sys.stdin.read().split()
    if not input_data:
        return
    iterator = iter(input_data)
    
    try:
        xA, yA, zA = float(next(iterator)), float(next(iterator)), float(next(iterator))
        xB, yB, zB = float(next(iterator)), float(next(iterator)), float(next(iterator))
        v_speed = float(next(iterator))
        xl, yl, zl = float(next(iterator)), float(next(iterator)), float(next(iterator))
        N = int(next(iterator))
        M = int(next(iterator))
        
        boats = []
        for _ in range(N):
            boats.append((float(next(iterator)), float(next(iterator)), float(next(iterator))))
            
        mountains = []
        for _ in range(M):
            mountains.append((float(next(iterator)), float(next(iterator)), float(next(iterator)), float(next(iterator))))
    except StopIteration:
        return

    A = (xA, yA, zA)
    B = (xB, yB, zB)
    vec_AB = sub(B, A)
    dist_AB = math.sqrt(norm_sq(vec_AB))
    dir_AB = normalize(vec_AB) if dist_AB > 1e-9 else (0,0,0)
    
    l_vec = normalize((xl, yl, zl))
    r_vec = cross(l_vec, (0, 0, 1))
    if norm_sq(r_vec) < 1e-9:
        r_vec = (1, 0, 0)
    else:
        r_vec = normalize(r_vec)
    u_vec = cross(r_vec, l_vec)
    
    EPS = 1e-9
    boat_intervals = []

    for i in range(N):
        Q = boats[i]
        QA = sub(Q, A)
        
        c_z0 = dot(QA, l_vec)
        c_z1 = -dot(dir_AB, l_vec)
        
        c_x0 = dot(QA, r_vec)
        c_x1 = -dot(dir_AB, r_vec)
        
        c_y0 = dot(QA, u_vec)
        c_y1 = -dot(dir_AB, u_vec)
        
        def intersect_linear(curr_int, p, q):
            start, end = curr_int
            if abs(q) < EPS:
                if p < -EPS: return (1, 0)
                else: return (start, end)
            root = -p / q
            if q > 0:
                return (max(start, root), end)
            else:
                return (start, min(end, root))
        
        valid_int = (0.0, dist_AB)
        
        valid_int = intersect_linear(valid_int, c_z0 - EPS, c_z1)
        valid_int = intersect_linear(valid_int, c_z0 - c_x0, c_z1 - c_x1)
        valid_int = intersect_linear(valid_int, c_z0 + c_x0, c_z1 + c_x1)
        valid_int = intersect_linear(valid_int, c_z0 - c_y0, c_z1 - c_y1)
        valid_int = intersect_linear(valid_int, c_z0 + c_y0, c_z1 + c_y1)
        
        if valid_int[0] > valid_int[1] + EPS:
            boat_intervals.append([])
            continue
            
        current_intervals = [valid_int]
        
        for mx, my, mh, mr in mountains:
            if not current_intervals:
                break
            
            def is_blocked(d):
                C = add(A, mul(dir_AB, d))
                Cx, Cy, Cz = C[0]-mx, C[1]-my, C[2]
                Qx, Qy, Qz = Q[0]-mx, Q[1]-my, Q[2]
                
                dx, dy, dz = Qx-Cx, Qy-Cy, Qz-Cz
                
                mu = mr / mh
                mu2 = mu * mu
                
                A_c = dx*dx + dy*dy - mu2 * dz*dz
                B_c = 2 * (Cx*dx + Cy*dy - mu2 * (Cz-mh)*dz)
                C_c = Cx*Cx + Cy*Cy - mu2 * (Cz-mh)**2
                
                roots = solve_quadratic(A_c, B_c, C_c)
                
                in_intervals = []
                if not roots:
                    if A_c < 0 or (abs(A_c) < EPS and B_c < 0) or (abs(A_c)<EPS and abs(B_c)<EPS and C_c <= EPS):
                         in_intervals.append((0.0, 1.0))
                else:
                    r1, r2 = roots
                    mid = (r1+r2)/2
                    val = A_c*mid*mid + B_c*mid + C_c
                    if val <= EPS:
                        in_intervals.append((max(0.0, r1), min(1.0, r2)))
                    else:
                        if r1 > 0: in_intervals.append((0.0, min(1.0, r1)))
                        if r2 < 1: in_intervals.append((max(0.0, r2), 1.0))
                
                real_blocked = False
                for s1, s2 in in_intervals:
                    if s1 > s2 + EPS: continue
                    
                    z_start, z_end = 0.0, mh
                    gs = intersect_linear((s1, s2), Cz, dz)
                    gs = intersect_linear(gs, mh - Cz, -dz)
                    
                    if gs[0] <= gs[1] + EPS:
                        real_blocked = True
                        break
                
                if real_blocked: return True
                
                if abs(dz) > EPS:
                    s_base = -Cz / dz
                    if 0 <= s_base <= 1:
                        ix = Cx + s_base*dx
                        iy = Cy + s_base*dy
                        if ix*ix + iy*iy <= mr*mr + EPS:
                            return True
                return False

            Q_s = (Q[0]-mx, Q[1]-my, Q[2]-mh)
            A_s = (A[0]-mx, A[1]-my, A[2]-mh)
            D_vec = dir_AB
            mu = mr / mh
            mu2 = mu * mu
            
            def q_prod(v1, v2): return v1[0]*v2[0] + v1[1]*v2[1] - mu2*v1[2]*v2[2]
            
            QMQ = q_prod(Q_s, Q_s)
            AMA = q_prod(A_s, A_s)
            AMD = q_prod(A_s, D_vec)
            DMD = q_prod(D_vec, D_vec)
            QMA = q_prod(Q_s, A_s)
            QMD = q_prod(Q_s, D_vec)
            
            k2 = QMD**2 - QMQ*DMD
            k1 = 2*(QMA*QMD - QMQ*AMD)
            k0 = QMA**2 - QMQ*AMA
            
            crit_points = solve_quadratic(k2, k1, k0)
            
            Q_r = (Q[0]-mx, Q[1]-my, Q[2])
            A_r = (A[0]-mx, A[1]-my, A[2])
            
            D0 = Q_r[2] - A_r[2]
            D1 = -D_vec[2]
            NX0 = Q_r[2]*A_r[0] - A_r[2]*Q_r[0]
            NX1 = Q_r[2]*D_vec[0] - D_vec[2]*Q_r[0]
            NY0 = Q_r[2]*A_r[1] - A_r[2]*Q_r[1]
            NY1 = Q_r[2]*D_vec[1] - D_vec[2]*Q_r[1]
            mr2 = mr*mr
            
            qk2 = NX1**2 + NY1**2 - mr2 * D1**2
            qk1 = 2*(NX0*NX1 + NY0*NY1 - mr2 * D0 * D1)
            qk0 = NX0**2 + NY0**2 - mr2 * D0**2
            
            crit_points.extend(solve_quadratic(qk2, qk1, qk0))
            
            next_intervals = []
            for start, end in current_intervals:
                local_pts = [start, end]
                for p in crit_points:
                    if start + EPS < p < end - EPS:
                        local_pts.append(p)
                local_pts.sort()
                
                for k in range(len(local_pts)-1):
                    p1, p2 = local_pts[k], local_pts[k+1]
                    mid = (p1+p2)/2
                    if not is_blocked(mid):
                        next_intervals.append((p1, p2))
            
            current_intervals = next_intervals
            
        boat_intervals.append(current_intervals)
        
    events = []
    for i, intervals in enumerate(boat_intervals):
        for start, end in intervals:
            if v_speed > EPS:
                t_start = start / v_speed
                t_end = end / v_speed
            else:
                t_start, t_end = 0.0, 0.0
            
            events.append((t_start, 1, i+1))
            events.append((t_end, -1, i+1))
            
    events.sort(key=lambda x: (x[0], -x[1]))
    
    max_k = 0
    best_t = 0.0
    current_ids = set()
    final_ids = []
    
    if not events:
        print(f"0.00000")
        print(f"0")
        return

    for t, type_, bid in events:
        if type_ == 1:
            current_ids.add(bid)
        else:
            if bid in current_ids:
                current_ids.remove(bid)
        
        if len(current_ids) > max_k:
            max_k = len(current_ids)
            best_t = t
            final_ids = sorted(list(current_ids))
        elif len(current_ids) == max_k and max_k > 0:
            if not final_ids:
                final_ids = sorted(list(current_ids))
                best_t = t
                
    print(f"{best_t:.5f}")
    print(f"{max_k}")
    for bid in final_ids:
        print(bid)

if __name__ == '__main__':
    solve()