#include"Vector.h"
#include<iostream>
int main() {
    Vector a;
    Vector b(4, 3);
    a = Vector(1, 1);
    b = a + b;
    b = a + Vector(2, 4);
    std::cout << area(Vector(6, 3), b) << std::endl;
    Vector c(-1, -1);
    std::cout << c << '\n';
    std::cout << agl(c, a) << '\n';
    c[5];
    getchar();
    getchar();
    return 0;
}



