#ifndef _NickwideN_Geometry_H
#define _NickwideN_Geometry_H
#include<iostream>
//-------------------------------------------------------------------------------------------------
//                          README.md here: https://github.com/NickwideN/Geometry
//-------------------------------------------------------------------------------------------------
namespace Geometry {
    typedef double coordinate_t;
    const static int DIMENTION = 2;
    constexpr static coordinate_t default_value = 0;
    class Point;
    class Segment;
    class Vector {
    private:
        coordinate_t coordinates[DIMENTION];
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
        Vector & operator -= (const Vector & vector);
        friend Vector operator - (const Vector & vector_1, const Vector & vector_2);
        Vector operator + () const;
        Vector operator - () const;
        Vector & operator *= (const coordinate_t & coefficient);
        friend Vector operator * (const Vector & vector, const coordinate_t & coefficient);
        friend coordinate_t operator * (const Vector & vector_1, const Vector & vector_2);
        friend coordinate_t scalar_product(const Vector & vector_1, const Vector & vector_2);
        friend Vector vector_product(const Vector & vector_1, const Vector & vector_2);
        friend coordinate_t skew_product(const Vector & vector_1, const Vector & vector_2);
        coordinate_t operator [] (const int index) const;
        coordinate_t & operator [] (const int index);
        friend double abs(const Vector & vector);
        friend double sin_agl(const Vector & vector_1, const Vector & vector_2);
        friend double cos_agl(const Vector & vector_1, const Vector & vector_2);
        friend double tan_agl(const Vector & vector_1, const Vector & vector_2);
        friend double agl(const Vector & vector_1, const Vector & vector_2);
        friend bool are_collinear(const Vector & vector_1, const Vector & vector_2);
        friend bool are_coincident(const Vector & vector_1, const Vector & vector_2);
        friend std::ostream & operator << (std::ostream & os, const Vector & vector);
        friend std::istream & operator >> (std::istream & is, Vector & vector);
    };


    class Shape {
    public:
        virtual double getArea() {
            return 0.0;
        }
        virtual Shape & move(const Vector & vector) = 0;
        virtual bool has_point(const Point & point) const = 0;
        virtual bool has_intarsection_with(const Segment & segment) const = 0;

    };

    class Point: public Shape {
    private:
        Vector radius_vector;
    public:
        Point();
        Point(coordinate_t coor_0, coordinate_t coor_1); ///////////////////////////////////////////
        Point & move(const Vector & vector) override;
        bool has_point(const Point & point) const override;
        bool has_intarsection_with(const Segment & segment) const override;
        friend std::ostream & operator << (std::ostream & os, const Point & vector);
        friend std::istream & operator >> (std::istream & is, Point & vector);
        Vector get_radius_vector();
    };

    class Segment: public Shape {
    private:
        Point point_1;
        Point point_2;
    public:
        //Конструкторы
        Segment & move(const Vector & vector) override;
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
}

