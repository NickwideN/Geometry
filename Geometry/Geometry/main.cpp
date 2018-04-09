#include"Geometry.h"
#include<iostream>
#include<iomanip>
using namespace Geometry;
using namespace std;
// F //wrong
/*
int main() {
    Segment a;
    Segment b;
    cin >> a;
    cin >> b;
    cout << (a.has_intarsection_with(b) ? "YES" : "NO");
    getchar();
    getchar();
    return 0;
}
*/

// G // wrong
/*
int main() {
    Segment a;
    Segment b;
    cin >> a;
    cin >> b;
    cout << setprecision(9) << fixed << distance_between(a, b);
    getchar();
    getchar();
    return 0;
}
*/

// E // wrong
/*
int main() {
    Point C;
    Point A;
    Point B;
    cin >> C >> A >> B;
    cout << setprecision(9) << fixed;
    cout << C.distance_to(Line(A, B)) << '\n';
    cout << C.distance_to(Ray(A, B)) << '\n';
    cout << C.distance_to(Segment(A, B)) << '\n';
    getchar();
    getchar();
    return 0;
}
*/

// C true
/*
int main() {
    int coef[6];
    for (int i = 0; i < 6; ++i) {
        cin >> coef[i];
    }
    Line line_1(coef[0], coef[1], coef[2]);
    Line line_2(coef[3], coef[4], coef[5]);
    cout << std::setprecision(9) << fixed << line_1.get_direction() << '\n' << line_2.get_direction() << '\n';
    if (are_intersecting(line_1, line_2)) {
        cout << intersection(line_1, line_2);
    }
    else if (are_parallel(line_1, line_2)) {
        cout << distance_between(line_1, line_2);
    }
    getchar();
    getchar();
    return 0;
}
*/

// Polygon

int main() {
    Point a(1, 2);
    Point b(1, 4);
    Point c(5, 3);
    Point d(3, 7);
    Point e(6, 2);
    Point p[5] = { a, b, a, d, e };
    Polygon pol1(5, p);
    Polygon pol2(3, a, b, c);
    pol1 = pol2;
    cout << pol1[2];
    getchar();
    getchar();
    return 0;
}
