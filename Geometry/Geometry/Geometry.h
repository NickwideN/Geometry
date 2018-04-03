#ifndef _NickwideN_Geometry_H
#define _NickwideN_Geometry_H
#include<iostream>
//-------------------------------------------------------------------------------------------------
//                          README.md here: https://github.com/NickwideN/Geometry
//-------------------------------------------------------------------------------------------------
class Vector {
private:
    const static int DIMENTION = 2;
    typedef double coordinate_t;
    coordinate_t coordinates[DIMENTION];
    constexpr static coordinate_t default_value = 0;

    template <typename T1 = const char *, typename T2 = const char *, typename T3 = const char *, typename T4 = const char *>
    friend void error(T1 p1, T2 p2 = "", T3 p3 = "", T4 p4 = "");
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


class Shape {
public:
    virtual double getArea() {
        return 0.0;
    }
    template<typename shape_t>
    shape_t move(const Vector & vector);

    virtual bool has_point(const Point & point) const = 0;
    virtual bool has_intarsection_with(const Segment & segment) const = 0;
};

class Point {
private:
    Vector radius_vector;
    typedef double coordinate_t;
    class Segment {};
public:
    Point();
    Point(coordinate_t coor_0, coordinate_t coor_1);
    Point & move(const Vector & vector);

    friend std::ostream & operator << (std::ostream & os, const Point & vector);
    friend std::istream & operator >> (std::istream & is, Point & vector);
};

class Segment: public Shape{
private:
    Point point_1;
    Point point_2;
public:
    template<typename Segment>
    Segment move(const Vector & vector);
    virtual bool has_point(const Point & point) const override;
    virtual bool has_intarsection_with(const Segment & segment) const override;

};




#endif // !_NickwideN_Geometry_H

template<typename user_t>
inline Vector::Vector(const user_t coor_0, const user_t coor_1, const user_t coor_2, const user_t coor_3, ...) {
    //doesn't have a cheking if it inputs number of coors that less than DIMENTION
    user_t *pointer_coor = &coor_0;
    for (int i = 0; i < DIMENTION; ++i) {
        coordinates[i] = *pointer_coor;
        ++pointer_coor;
    }
}

template<typename T1, typename T2, typename T3, typename T4>
inline void error(T1 p1, T2 p2, T3 p3, T4 p4) {
    std::cerr << "Error: " << p1 << p2 << p3 << p4 << '\n';
    getchar();
    getchar();
    std::exit(1);
}


template<typename Segment>
Segment Segment::move(const Vector & vector) {
    this->point_1.move(vector);
    this->point_2.move(vector);
    return *this;
}

template<typename shape_t>
inline shape_t Shape::move(const Vector & vector){
    return *this;
}
