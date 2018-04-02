#ifndef _NickwideN_Vector_H
#define _NickwideN_Vector_H

class Vector {
private:
    const static int DIMENTION = 2;
    typedef double coordinate_t;
    coordinate_t coordinates[DIMENTION];
    constexpr static coordinate_t default_value = 0;
public:
    Vector();
    template<typename number_t>
    Vector(number_t coors, ...);
};
#endif // !_NickwideN_Vector_H

template<typename number_t>
inline Vector::Vector(number_t coors, ...){  //doesn't have a cheking if it inputs number of coors that less than DIMENTION
    number_t *pointer_coor = &coors;
    for (int i = 0; i < DIMENTION; ++i) {
        coordinates[i] = *pointer_coor;
        ++pointer_coor;
    }
}
