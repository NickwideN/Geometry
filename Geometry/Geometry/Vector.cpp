#include"Vector.h"
#include<stdarg.h>
Vector::Vector() {
    for (int i = 0; i < DIMENTION; ++i) {
        coordinates[i] = default_value;
    }
}
