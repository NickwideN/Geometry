#include"Geometry.h"
#include<iostream>
int main() {
    Geometry::Vector a;
    Geometry::Vector b(4, 3);
    a = Geometry::Vector(1, 1);
    b = a + b;
    b = a + Geometry::Vector(2, 4);
    std::cout << Geometry::skew_product(Geometry::Vector(6, 3), b) << std::endl;
    Geometry::Vector c(-1, -1);
    std::cout << c << '\n';
    std::cout << agl(c, a) << '\n';
    c[1];

    Geometry::Point p;
    std::cout << p;
    p.move(c);
    std::cout << p;
    getchar();
    getchar();
    return 0;
}



