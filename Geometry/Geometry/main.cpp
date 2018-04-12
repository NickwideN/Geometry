#include"Geometry.h"
#include<iostream>
#include<iomanip>
using namespace Geometry;
using namespace std;
// B
/*
int main() {
Vector radius_vector[4];
for (int i = 0; i < 4; ++i) {
cin >> radius_vector[i];
}
Vector vector_1 = radius_vector[1] - radius_vector[0];
Vector vector_2 = radius_vector[3] - radius_vector[2];
cout.fill('0');
cout.width();
cout << std::setprecision(9) << fixed << abs(vector_1) << ' ' << abs(vector_2) << '\n';
cout << vector_1 + vector_2 << '\n';
cout << scalar_product(vector_1, vector_2) << ' ' << skew_product(vector_1, vector_2) << '\n';
cout << abs(skew_product(vector_1, vector_2)) / 2 << '\n';
getchar();
getchar();
return 0;
}
*/

// C right
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
// D  right
/*
int main() {
    Point C;
    Point A;
    Point B;
    cin >> C >> A >> B;
    Line l(A, B);
    Segment s(A, B);
    Ray r(A, B);
    cout << (l.has_point(C) ? "YES" : "NO") << '\n';
    cout << (r.has_point(C) ? "YES" : "NO") << '\n';
    cout << (s.has_point(C) ? "YES" : "NO") << '\n';
    getchar();
    getchar();
    return 0;
}
*/

// E // right
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

// F //right
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

// G // right
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
// H   wrong
/*
int main() {
    int points_cnt;
    cin >> points_cnt;
    Polygon pol(points_cnt);
    cin >> pol;
    cout << (pol.is_convex() ? "YES" : "NO");
    getchar();
    getchar();
    return 0;
}
*/
// I  //right
/*
int main() {
    int points_cnt;
    cin >> points_cnt;
    Point check_point;
    Polygon pol(points_cnt);
    cin >> check_point;
    cin >> pol;
    cout << (pol.has_point(check_point) ? "YES" : "NO");
    getchar();
    getchar();
    return 0;
}
*/

// J  true
/*
int main() {
    int points_cnt;
    cin >> points_cnt;
    Polygon pol(points_cnt);
    cin >> pol;
    cout << area(pol);
    getchar();
    getchar();
    return 0;
}
*/

// K
int main() {


    getchar();
    getchar();
    return 0;
}

// Checking Polygon
/*
int main() {
Point a(4, 1);
Point b(2, 2);
Point c(0, 6);
Point d(0, 1);
Point e(1, 1);
Segment seg(a, c);
Point p[5] = { a, b, a, d, e };
Polygon pol1(5, p);
Polygon pol2(4, a, b, c, d);
pol1 = pol2;
cout << area(pol2) << '\n';
cout << pol2.has_point(e);
cout << pol2.is_convex();
cout << pol2.has_intarsection_with(seg);
cout << pol2.has_point(a);
cout << pol2.has_point(c);
getchar();
getchar();
return 0;
}*/