#include"Geometry.h"
#include<iostream>
using namespace Geometry;
using namespace std;
int main() {
    Vector a;
    Vector b(4, 3);
    a = Vector(1, 1);
    b = a + b;
    b = a + Geometry::Vector(2, 4);
    std::cout << Geometry::skew_product(Geometry::Vector(6, 3), b) << std::endl;
    Geometry::Vector c(-1, -1);
    std::cout << c << '\n';
    std::cout << agl(c, a) << '\n';
    c[1];
    Geometry::Point p(1, 4);
    std::cout << p << std::endl;
    p.move(c);
    cout << p << endl << endl;
    Point con_point = Vector(2, 3);
    Segment segment(Point(1, 1), Point(3, 5));
    cout << segment.has_point(con_point) << endl;
    Point x(b);
    getchar();
    getchar();
    return 0;
}



