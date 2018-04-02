#ifndef _NickwideN_Vector_H
#define _NickwideN_Vector_H

class Vector {
private:
    const static int DIMENTION = 4;
    typedef double coordinate_t;
    coordinate_t coordinates[DIMENTION];
    constexpr static coordinate_t default_value = 0;
public:
    Vector();
    template<typename user_t>
    Vector(user_t coor_0, user_t coor_1, user_t coor_2, ...); // only when all types are similar,
                                                                    //doesn't have a cheking if number of coors is less than DIMENTION
    Vector(coordinate_t coor_0, coordinate_t coor_1 = default_value);  //let to pass several types

};
#endif // !_NickwideN_Vector_H

template<typename user_t>
inline Vector::Vector(user_t coor_0, user_t coor_1, user_t coor_2, ...){  //doesn't have a cheking if it inputs number of coors that less than DIMENTION
    user_t *pointer_coor = &coor_0;
    for (int i = 0; i < DIMENTION; ++i) {
        coordinates[i] = *pointer_coor;
        ++pointer_coor;
    }
}


