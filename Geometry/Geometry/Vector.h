#ifndef _NickwideN_Vector_H
#define _NickwideN_Vector_H
#include<iostream>
class Vector {
private:
    const static int DIMENTION = 2;
    typedef double coordinate_t;
    coordinate_t coordinates[DIMENTION];
    constexpr static coordinate_t default_value = 0;
public:
    Vector();
    template<typename user_t>
    Vector(const user_t coor_0, const user_t coor_1, const user_t coor_2, const user_t coor_3, ...); 
                                            // only when all types are similar,
                                            //doesn't have a cheking if number of coors is less than DIMENTION
    Vector(const coordinate_t coor_0, const coordinate_t coor_1 = default_value, const coordinate_t coor_2 = default_value);
                                            //let to pass several types
    Vector & operator += (const Vector & vector);
    friend Vector operator + (const Vector & vector_1, const Vector & vector_2);
    Vector operator + () const;
    Vector operator - () const;
    Vector & operator *= (const coordinate_t & coefficient);
    friend Vector operator * (const Vector & vector, const coordinate_t & coefficient);
    friend coordinate_t operator * (const Vector & vector_1, const Vector & vector_2);
    friend coordinate_t scalar_product (const Vector & vector_1, const Vector & vector_2);
    friend Vector vector_product (const Vector & vector_1, const Vector & vector_2);
    friend coordinate_t area(const Vector & vector_1, const Vector & vector_2);
    coordinate_t operator [] (const int index) const;
    coordinate_t & operator [] (const int index);
    friend double abs (const Vector & vector);
    friend double agl(const Vector & vector_1, const Vector & vector_2);

    
    friend std::ostream & operator << (std::ostream & os, const Vector & vector);
    friend std::istream & operator >> (std::istream & is, Vector & vector);
};
#endif // !_NickwideN_Vector_H

template<typename user_t>
inline Vector::Vector(const user_t coor_0, const user_t coor_1, const user_t coor_2, const user_t coor_3, ...){ 
                                        //doesn't have a cheking if it inputs number of coors that less than DIMENTION
    user_t *pointer_coor = &coor_0;
    for (int i = 0; i < DIMENTION; ++i) {
        coordinates[i] = *pointer_coor;
        ++pointer_coor;
    }
}

