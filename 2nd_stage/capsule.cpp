#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
#include <algorithm>

using namespace std;

struct Vec {
    double x, y, z;
    Vec operator+(const Vec& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec operator-(const Vec& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec operator*(double s) const { return {x * s, y * s, z * s}; }
};

double dot(const Vec& a, const Vec& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec cross(const Vec& a, const Vec& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

double len(const Vec& a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

Vec normalize(const Vec& a) {
    double l = len(a);
    if (l < 1e-12) return {0.0, 0.0, 0.0};
    return {a.x / l, a.y / l, a.z / l};
}

struct Plane {
    Vec normal;
    double d;
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    vector<Vec> vertices(N);
    for (int i = 0; i < N; ++i) {
        cin >> vertices[i].x >> vertices[i].y >> vertices[i].z;
    }

    int M;
    cin >> M;
    vector<vector<int>> faces(M);
    for (int i = 0; i < M; ++i) {
        int k;
        cin >> k;
        faces[i].resize(k);
        for (int j = 0; j < k; ++j) {
            cin >> faces[i][j];
        }
    }

    double X0, Y0, Z0, U, V, W, ta, tb;
    cin >> X0 >> Y0 >> Z0 >> U >> V >> W >> ta >> tb;

    Vec C1 = {X0, Y0, Z0};
    Vec L_vec = {U, V, W};

    double cx = 0.0, cy = 0.0, cz = 0.0;
    for (const auto& v : vertices) {
        cx += v.x;
        cy += v.y;
        cz += v.z;
    }
    Vec centroid = {cx / N, cy / N, cz / N};

    vector<Plane> planes;
    for (const auto& f_idxs : faces) {
        Vec p0 = vertices[f_idxs[0]];
        Vec normal = {0.0, 0.0, 0.0};
        bool found_normal = false;

        for (size_t i = 0; i < f_idxs.size() - 2; ++i) {
            Vec p1 = vertices[f_idxs[i+1]];
            Vec p2 = vertices[f_idxs[i+2]];
            Vec edge1 = p1 - p0;
            Vec edge2 = p2 - p1;
            Vec n = cross(edge1, edge2);
            if (len(n) > 1e-9) {
                normal = normalize(n);
                found_normal = true;
                break;
            }
        }

        if (!found_normal) continue;

        Vec vec_to_center = centroid - p0;
        if (dot(normal, vec_to_center) < 0) {
            normal = normal * -1.0;
        }

        double d_val = -dot(normal, p0);
        planes.push_back({normal, d_val});
    }

    Vec P_start = C1 + L_vec * ta;
    Vec P_end = C1 + L_vec * tb;
    Vec L_unit = normalize(L_vec);
    double L_len_val = len(L_vec);

    double R_max = numeric_limits<double>::infinity();

    for (const auto& plane : planes) {
        double dist_start = dot(plane.normal, P_start) + plane.d;
        double dist_end = dot(plane.normal, P_end) + plane.d;
        double min_dist = min(dist_start, dist_end);

        if (min_dist < -1e-9) {
            cout << "0" << endl;
            return 0;
        }

        double sin_theta = 1.0;
        if (L_len_val > 1e-12) {
            double cos_theta = dot(plane.normal, L_unit);
            if (cos_theta > 1.0) cos_theta = 1.0;
            if (cos_theta < -1.0) cos_theta = -1.0;
            sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        }

        if (sin_theta > 1e-9) {
            double r_limit = min_dist / sin_theta;
            if (r_limit < R_max) {
                R_max = r_limit;
            }
        }
    }

    if (R_max <= 1e-9 || R_max == numeric_limits<double>::infinity()) {
        cout << "0" << endl;
    } else {
        cout << "1 " << fixed << setprecision(5) << R_max << endl;
    }

    return 0;
}