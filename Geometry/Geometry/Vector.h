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
    Vector(coordinate_t coors, ...);
};
#endif // !_NickwideN_Vector_H
