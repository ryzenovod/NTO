#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <set>
#include <tuple>
#include <numbers>

using namespace std;

struct Point {
    long long x, y, z;
    auto operator<=>(const Point&) const = default;
};

struct Drone {
    long long px, py, pz;
    long long dx, dy, dz;
    int alpha;
};

bool is_in_cone(const Point& p, const Drone& drone) {
    long long vx = p.x - drone.px;
    long long vy = p.y - drone.py;
    long long vz = p.z - drone.pz;

    long long lenV_sq = vx * vx + vy * vy + vz * vz;
    
    if (lenV_sq == 0) return true;

    long double dot = (long double)vx * drone.dx + (long double)vy * drone.dy + (long double)vz * drone.dz;

    if (drone.alpha < 90 && dot < 0) return false;

    long long lenD_sq = drone.dx * drone.dx + drone.dy * drone.dy + drone.dz * drone.dz;
    
    long double lenV = sqrt((long double)lenV_sq);
    long double lenD = sqrt((long double)lenD_sq);
    
    if (lenD == 0) return false; 

    long double cos_theta = dot / (lenV * lenD);
    
    if (cos_theta > 1.0) cos_theta = 1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;
    
    double theta = acos(cos_theta);
    double alpha_rad = drone.alpha * numbers::pi / 180.0;
    
    return theta <= alpha_rad + 1e-9;
}

struct Entry {
    Point p;
    int building_id;
    
    bool operator<(const Entry& other) const {
        return p < other.p;
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    Drone drone;
    if (!(cin >> drone.px >> drone.py >> drone.pz >> drone.dx >> drone.dy >> drone.dz >> drone.alpha)) return 0;

    int N;
    if (!(cin >> N)) return 0;

    vector<Entry> entries;
    entries.reserve(N * 4);

    for (int i = 0; i < N; ++i) {
        int building_id;
        cin >> building_id;
        for (int j = 0; j < 4; ++j) {
            Point p;
            cin >> p.x >> p.y >> p.z;
            if (is_in_cone(p, drone)) {
                entries.push_back({p, building_id});
            }
        }
    }

    sort(entries.begin(), entries.end());

    vector<Point> conflict_points;
    set<pair<int, int>> conflict_pairs;

    for (size_t i = 0; i < entries.size(); ) {
        size_t j = i;
        vector<int> current_ids;
        
        current_ids.push_back(entries[i].building_id);
        j++;
        while (j < entries.size() && entries[j].p == entries[i].p) {
            current_ids.push_back(entries[j].building_id);
            j++;
        }

        sort(current_ids.begin(), current_ids.end());
        auto last = unique(current_ids.begin(), current_ids.end());
        current_ids.erase(last, current_ids.end());

        if (current_ids.size() > 1) {
            conflict_points.push_back(entries[i].p);
            
            for (size_t a = 0; a < current_ids.size(); ++a) {
                for (size_t b = a + 1; b < current_ids.size(); ++b) {
                    conflict_pairs.insert({current_ids[a], current_ids[b]});
                }
            }
        }
        i = j;
    }

    cout << conflict_points.size() << "\n";
    for (const auto& p : conflict_points) {
        cout << p.x << " " << p.y << " " << p.z << "\n";
    }

    cout << conflict_pairs.size() << "\n";
    for (const auto& p : conflict_pairs) {
        cout << p.first << " " << p.second << "\n";
    }

    return 0;
}