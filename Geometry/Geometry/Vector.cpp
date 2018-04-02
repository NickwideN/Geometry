#include"Vector.h"
#include<stdarg.h>
Vector::Vector() {
    for (int i = 0; i < DIMENTION; ++i) {
        coordinates[i] = default_value;
    }
}

Vector::Vector(coordinate_t coor_0, coordinate_t coor_1){ //let put a several types
    coordinates[0] = coor_0;
    if (1 < DIMENTION)
        coordinates[1] = coor_1;
    for (int i = 2; i < DIMENTION; ++i) {
        coordinates[i] = default_value;
    }
}

