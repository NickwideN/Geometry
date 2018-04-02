#include"Vector.h"
#include<iostream>
int main() {
    Vector a;
    Vector b(4, 3);
    a = b;
    b = a + b;
    b = a + Vector(2, 4);
    std::cout << area(Vector(6, 3), b) << std::endl;
    Vector c(4, 5);
    std::cout << c << '\n';
    std::cin >> c;
    std::cout << c << '\n';
    c[5];
    getchar();
    getchar();
    return 0;
}



