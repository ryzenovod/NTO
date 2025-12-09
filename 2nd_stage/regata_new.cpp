#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <set>

using namespace std;

const double EPS = 1e-9;
const double INF = 1e18;

struct Point {
    double x, y, z;
};

Point operator+(const Point& a, const Point& b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
Point operator-(const Point& a, const Point& b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
Point operator*(const Point& a, double s) { return {a.x*s, a.y*s, a.z*s}; }
double dot(const Point& a, const Point& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
Point cross(const Point& a, const Point& b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
double norm_sq(const Point& a) { return dot(a, a); }
Point normalize(const Point& a) {
    double l = sqrt(norm_sq(a));
    if (l < EPS) return {0,0,0};
    return a * (1.0/l);
}

vector<double> solve_quadratic(double a, double b, double c) {
    vector<double> roots;
    if (abs(a) < EPS) {
        if (abs(b) > EPS) roots.push_back(-c/b);
    } else {
        double d = b*b - 4*a*c;
        if (d >= 0) {
            double sqrt_d = sqrt(d);
            roots.push_back((-b - sqrt_d) / (2*a));
            roots.push_back((-b + sqrt_d) / (2*a));
        }
    }
    sort(roots.begin(), roots.end());
    return roots;
}

pair<double, double> intersect(pair<double, double> i1, pair<double, double> i2) {
    return {max(i1.first, i2.first), min(i1.second, i2.second)};
}

pair<double, double> intersect_linear(pair<double, double> interval, double p, double q) {
    if (abs(q) < EPS) {
        if (p < -EPS) return {1.0, 0.0};
        return interval;
    }
    double root = -p/q;
    if (q > 0) return intersect(interval, {root, INF});
    else return intersect(interval, {-INF, root});
}

struct Mountain {
    double x, y, h, r;
};

struct Event {
    double t;
    int type;
    int id;
    
    bool operator<(const Event& other) const {
        if (abs(t - other.t) > EPS) return t < other.t;
        return type > other.type;
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    Point A, B;
    double v;
    Point l_vec_raw;
    int N, M;

    if (!(cin >> A.x >> A.y >> A.z)) return 0;
    cin >> B.x >> B.y >> B.z;
    cin >> v;
    cin >> l_vec_raw.x >> l_vec_raw.y >> l_vec_raw.z;
    cin >> N >> M;

    vector<Point> boats(N);
    for(int i=0; i<N; ++i) cin >> boats[i].x >> boats[i].y >> boats[i].z;

    vector<Mountain> mountains(M);
    for(int i=0; i<M; ++i) cin >> mountains[i].x >> mountains[i].y >> mountains[i].h >> mountains[i].r;

    Point vec_AB = B - A;
    double dist_AB = sqrt(norm_sq(vec_AB));
    Point dir_AB = (dist_AB > EPS) ? normalize(vec_AB) : Point{0,0,0};

    Point l = normalize(l_vec_raw);
    Point r = cross(l, {0, 0, 1});
    if (norm_sq(r) < EPS) r = {1, 0, 0};
    else r = normalize(r);
    Point u = cross(r, l);

    vector<vector<pair<double, double>>> boat_intervals(N);

    for(int i=0; i<N; ++i) {
        Point Q = boats[i];
        Point QA = Q - A;

        double cz0 = dot(QA, l), cz1 = -dot(dir_AB, l);
        double cx0 = dot(QA, r), cx1 = -dot(dir_AB, r);
        double cy0 = dot(QA, u), cy1 = -dot(dir_AB, u);

        pair<double, double> valid = {0.0, dist_AB};
        valid = intersect_linear(valid, cz0 - EPS, cz1);
        valid = intersect_linear(valid, cz0 - cx0, cz1 - cx1);
        valid = intersect_linear(valid, cz0 + cx0, cz1 + cx1);
        valid = intersect_linear(valid, cz0 - cy0, cz1 - cy1);
        valid = intersect_linear(valid, cz0 + cy0, cz1 + cy1);

        if (valid.first > valid.second + EPS) continue;

        vector<pair<double, double>> current_intervals;
        current_intervals.push_back(valid);

        for(const auto& m : mountains) {
            if (current_intervals.empty()) break;

            Point Qs = {Q.x - m.x, Q.y - m.y, Q.z - m.h};
            Point As = {A.x - m.x, A.y - m.y, A.z - m.h};
            Point Qr = {Q.x - m.x, Q.y - m.y, Q.z};
            Point Ar = {A.x - m.x, A.y - m.y, A.z};

            double mu = m.r / m.h;
            double mu2 = mu * mu;

            auto q_dot = [&](Point v1, Point v2) {
                return v1.x*v2.x + v1.y*v2.y - mu2*v1.z*v2.z;
            };

            double QMQ = q_dot(Qs, Qs);
            double AMA = q_dot(As, As);
            double AMD = q_dot(As, dir_AB);
            double DMD = q_dot(dir_AB, dir_AB);
            double QMA = q_dot(Qs, As);
            double QMD = q_dot(Qs, dir_AB);

            double k2 = QMD*QMD - QMQ*DMD;
            double k1 = 2*(QMA*QMD - QMQ*AMD);
            double k0 = QMA*QMA - QMQ*AMA;
            vector<double> crits = solve_quadratic(k2, k1, k0);

            double D0 = Qr.z - Ar.z;
            double D1 = -dir_AB.z;
            double NX0 = Qr.z*Ar.x - Ar.z*Qr.x;
            double NX1 = Qr.z*dir_AB.x - dir_AB.z*Qr.x;
            double NY0 = Qr.z*Ar.y - Ar.z*Qr.y;
            double NY1 = Qr.z*dir_AB.y - dir_AB.z*Qr.y;
            double mr2 = m.r * m.r;

            double qk2 = NX1*NX1 + NY1*NY1 - mr2*D1*D1;
            double qk1 = 2*(NX0*NX1 + NY0*NY1 - mr2*D0*D1);
            double qk0 = NX0*NX0 + NY0*NY0 - mr2*D0*D0;
            
            vector<double> base_crits = solve_quadratic(qk2, qk1, qk0);
            crits.insert(crits.end(), base_crits.begin(), base_crits.end());

            vector<pair<double, double>> next_intervals;
            for(auto interval : current_intervals) {
                vector<double> pts;
                pts.push_back(interval.first);
                pts.push_back(interval.second);
                for(double c : crits) {
                    if (c > interval.first + EPS && c < interval.second - EPS) pts.push_back(c);
                }
                sort(pts.begin(), pts.end());
                pts.erase(unique(pts.begin(), pts.end()), pts.end());

                for(size_t k=0; k < pts.size()-1; ++k) {
                    double t1 = pts[k], t2 = pts[k+1];
                    double mid = (t1+t2)/2.0;

                    Point C = A + dir_AB * mid;
                    Point Cm = {C.x - m.x, C.y - m.y, C.z};
                    Point Qm = {Q.x - m.x, Q.y - m.y, Q.z};
                    Point d_vec = Qm - Cm;

                    double ac = d_vec.x*d_vec.x + d_vec.y*d_vec.y - mu2*d_vec.z*d_vec.z;
                    double bc = 2*(Cm.x*d_vec.x + Cm.y*d_vec.y - mu2*(Cm.z - m.h)*d_vec.z);
                    double cc = Cm.x*Cm.x + Cm.y*Cm.y - mu2*(Cm.z - m.h)*(Cm.z - m.h);
                    
                    vector<double> s_roots = solve_quadratic(ac, bc, cc);
                    vector<pair<double, double>> s_intervals;
                    
                    if (s_roots.empty()) {
                         if (ac < -EPS || (abs(ac) < EPS && bc < -EPS) || (abs(ac)<EPS && abs(bc)<EPS && cc < -EPS)) {
                             s_intervals.push_back({0, 1});
                         }
                    } else {
                        double mv = (s_roots[0] + s_roots[1])/2.0;
                        if (ac*mv*mv + bc*mv + cc < -EPS) {
                            s_intervals.push_back({max(0.0, s_roots[0]), min(1.0, s_roots[1])});
                        } else {
                            if (s_roots[0] > EPS) s_intervals.push_back({0.0, min(1.0, s_roots[0])});
                            if (s_roots[1] < 1.0 - EPS) s_intervals.push_back({max(0.0, s_roots[1]), 1.0});
                        }
                    }

                    bool blocked = false;
                    for(auto si : s_intervals) {
                         if (si.first > si.second + EPS) continue;
                         pair<double, double> z_range = {0.0, m.h};
                         pair<double, double> s_valid = intersect_linear(si, Cm.z, d_vec.z);
                         s_valid = intersect_linear(s_valid, m.h - Cm.z, -d_vec.z);
                         
                         if (s_valid.first <= s_valid.second + EPS) {
                             blocked = true;
                             break;
                         }
                    }

                    if (!blocked && abs(d_vec.z) > EPS) {
                        double s_base = -Cm.z / d_vec.z;
                        if (s_base >= -EPS && s_base <= 1.0 + EPS) {
                            Point inter = Cm + d_vec * s_base;
                            if (inter.x*inter.x + inter.y*inter.y <= m.r*m.r + EPS) blocked = true;
                        }
                    }

                    if (!blocked) next_intervals.push_back({t1, t2});
                }
            }
            current_intervals = next_intervals;
        }
        boat_intervals[i] = current_intervals;
    }

    vector<Event> events;
    for(int i=0; i<N; ++i) {
        for(auto interval : boat_intervals[i]) {
            double t_start = (v > EPS) ? interval.first / v : 0.0;
            double t_end = (v > EPS) ? interval.second / v : 0.0;
            events.push_back({t_start, 1, i+1});
            events.push_back({t_end, -1, i+1});
        }
    }
    sort(events.begin(), events.end());

    int max_k = 0;
    double best_t = 0.0;
    set<int> current_ids;
    vector<int> final_ids;

    if (events.empty()) {
        cout << "0.00000\n0" << endl;
        return 0;
    }

    for(const auto& e : events) {
        if (e.type == 1) current_ids.insert(e.id);
        else current_ids.erase(e.id);

        if ((int)current_ids.size() > max_k) {
            max_k = current_ids.size();
            best_t = e.t;
            final_ids.assign(current_ids.begin(), current_ids.end());
        } else if ((int)current_ids.size() == max_k && max_k > 0) {
            if (final_ids.empty()) {
                best_t = e.t;
                final_ids.assign(current_ids.begin(), current_ids.end());
            }
        }
    }

    cout << fixed << setprecision(5) << best_t << "\n";
    cout << max_k << "\n";
    sort(final_ids.begin(), final_ids.end());
    for(int id : final_ids) cout << id << "\n";

    return 0;
}