#include"Vector.h"
#include<stdarg.h>
Vector::Vector() {
    for (int i = 0; i < DIMENTION; ++i) {
        coordinates[i] = default_value;
    }
}

Vector::Vector(coordinate_t coor_1, ...) { //doesn't have a cheking if it inputs number of coors that less than DIMENTION
        va_list pointer_coor;             //--объ€вление указател€
        va_start(pointer_coor, coor_1);            //--инициализаци€ указател€
        coordinates[0] = coor_1;
        for (int i = 1; i < DIMENTION; ++i)
            coordinates[i] += va_arg(pointer_coor, coordinate_t);        //--перемещение указател€ 
        va_end(pointer_coor);                //--Ђзакрытиеї указател€
}
