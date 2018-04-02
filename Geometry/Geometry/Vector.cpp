#include"Vector.h"

Vector::Vector() {
    for (int i = 0; i < DIMENTION; ++i) {
        coordinates[i] = default_value;
    }
}

Vector::Vector(coordinate_t coors, ...) { //doesn't have a cheking if it inputs number of coors that less than DIMENTION
    coordinate_t *pointer_coor = &coors;
    for (int i = 0; i < DIMENTION; ++i) {
        coordinates[i] = *pointer_coor;
        ++pointer_coor;
    }
}
