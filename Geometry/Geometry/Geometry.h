//-------------------------------------------------------------------------------------------------
//                          README.md is here: https://github.com/NickwideN/Geometry
//-------------------------------------------------------------------------------------------------

#ifndef _NickwideN_Geometry_H
#define _NickwideN_Geometry_H
#include<iostream>
#include<cstring>    // strcmp
#include<cmath>      // sqrt
#include<algorithm>  // qsort

namespace Geometry {
    typedef double coordinate_t;
    typedef double scalar_t;
    const static int DIMENTION = 2;
    constexpr static coordinate_t default_value = 0;
    class Point;
    class Segment;
    class Line;
    class Ray;
    class Polygon;

    class Vector {
    private:
        coordinate_t coordinates[DIMENTION];
    public:
        Vector();
        template<typename user_t>
        Vector(const user_t coor_0, const user_t coor_1, const user_t coor_2, const user_t coor_3, ...);
        // only when all types are similar,
        // doesn't have a cheking if number of coors is less than DIMENTION
        Vector(const coordinate_t coor_0, const coordinate_t coor_1 = default_value, const coordinate_t coor_2 = default_value);
        // let to pass several types
        Vector(const Point & point_0, const Point & point_1);
        Vector & operator += (const Vector & vector);
        friend Vector operator + (const Vector & vector_0, const Vector & vector_1);
        Vector & operator -= (const Vector & vector);
        friend Vector operator - (const Vector & vector_0, const Vector & vector_1);
        Vector operator + () const;
        Vector operator - () const;
        Vector & operator *= (const scalar_t & coefficient);
        Vector & operator /= (const scalar_t & coefficient);
        friend Vector operator * (const Vector & vector, const scalar_t & coefficient);
        friend Vector operator / (const Vector & vector, const scalar_t & coefficient);
        friend Vector operator * (const scalar_t & coefficient, const Vector & vector);
        friend Vector operator / (const scalar_t & coefficient, const Vector & vector);
        friend scalar_t operator * (const Vector & vector_0, const Vector & vector_1);
        friend scalar_t scalar_product(const Vector & vector_0, const Vector & vector_1);
        friend Vector vector_product(const Vector & vector_0, const Vector & vector_1);
        friend coordinate_t skew_product(const Vector & vector_0, const Vector & vector_1);
        friend Vector normal_vector(const Vector & vector);
        coordinate_t operator [] (const int index) const;
        coordinate_t & operator [] (const int index);
        friend scalar_t abs(const Vector & vector);
        friend double sin(const Vector & vector_0, const Vector & vector_1);
        friend double cos(const Vector & vector_0, const Vector & vector_1);
        friend double tan(const Vector & vector_0, const Vector & vector_1);
        friend double agl(const Vector & vector_0, const Vector & vector_1);
        friend bool are_collinear(const Vector & vector_0, const Vector & vector_1);
        friend bool are_coincident(const Vector & vector_0, const Vector & vector_1);
        friend bool are_complanar(const Vector & vector_0, const Vector & vector_1, const Vector & vector_2); 
        friend bool are_co_directed(const Vector & vector_0, const Vector & vector_1);
        bool is_zero();
        friend std::ostream & operator << (std::ostream & os, const Vector & vector);
        friend std::istream & operator >> (std::istream & is, Vector & vector);
    };

    class Shape {
    public:
        virtual Shape & move(const Vector & vector) = 0;
        virtual bool has_point(const Point & point) const = 0;
        virtual bool has_intarsection_with(const Segment & segment) const = 0;
    };

    class Point : public Shape {
    private:
        Vector radius_vector;
    public:
        Point();
        Point(const coordinate_t & coor_0, const coordinate_t & coor_1 = default_value); //no konstruktors for Demention > 2//////////////////
        Point(const Vector radius_vector);

        Point & move(const Vector & vector) override;
        bool has_point(const Point & point) const override;
        bool has_intarsection_with(const Segment & segment) const override;

        friend scalar_t length(const Point & point_0, const Point & point_1);
        scalar_t distance_to(const Line & line) const;
        scalar_t distance_to(const Ray & ray) const;
        scalar_t distance_to(const Segment & segment) const;
        bool operator == (const Point & other) const;

        friend std::ostream & operator << (std::ostream & os, const Point & vector);
        friend std::istream & operator >> (std::istream & is, Point & vector);

        friend Vector;
        friend Segment;
        friend Line;
        friend Ray;
        friend Polygon;
        friend bool are_coincident(const Line & line_0, const Line & line_1);
        friend bool are_intersecting(const Line & line_0, const Line & line_1);
        friend bool are_skew(const Line & line_0, const Line & line_1);
        friend Point intersection(const Line & line_0, const Line & line_1);
        friend scalar_t distance_between(const Line & line_0, const Line & line_1);
        friend scalar_t area(const Polygon & polygon);
        friend Polygon convex_hull(const Polygon & polygon);
    private:
        //function for convex_hull (class Polygon)
        friend int origin_index(const Point * points, const int points_cnt);
    };

    class Segment : public Shape {
    private:
        Point point_0;
        Point point_1;
    public:
        Segment();
        Segment(const Point & point_0, const Point & point_1);
        Segment(const Point & origin, const Vector & direction, const scalar_t & length);

        Segment & move(const Vector & vector) override;
        virtual bool has_point(const Point & point) const override;
        virtual bool has_intarsection_with(const Segment & segment) const override;
        bool has_intarsection_with(const Ray & ray) const;

        friend scalar_t length(const Segment & segment);
        friend scalar_t distance_between(const Segment & segment_0, const Segment & segment_1);
        bool is_point() const;

        friend std::ostream & operator << (std::ostream & os, const Segment & segment);
        friend std::istream & operator >> (std::istream & is, Segment & segmant);

        friend Line;
        friend Point;
        friend Ray;
        friend Polygon;
    };

    class Line : public Shape {
    private:
        Point origin;
        Vector direction;
    public:
        Line(const Point & point_0, const Point & point_1);
        Line(const Point & origin, const Vector & vector, const char * name_vector = "direction"); // vector can be normal or direction
        Line(const coordinate_t & coefficient_of_x, const coordinate_t & coefficient_of_y, const coordinate_t & absolute_term);   // Ax + By + C = 0
        Line(const Segment & segmant);
        Line(const Ray & ray);

        Line & move(const Vector & vector) override;
        bool has_point(const Point & point) const override;
        bool has_intarsection_with(const Segment & segment) const override;

        Vector get_direction() const;
        friend bool are_coincident(const Line & line_0, const Line & line_1);
        friend bool are_parallel(const Line & line_0, const Line & line_1);
        friend bool are_intersecting(const Line & line_0, const Line & line_1);
        friend bool are_skew(const Line & line_0, const Line & line_1);
        friend Point intersection(const Line & line_0, const Line & line_1);
        friend scalar_t distance_between(const Line & line_0, const Line & line_1);
        friend std::ostream & operator << (std::ostream & os, const Line & line);
        friend std::istream & operator >> (std::istream & is, Line & line);

        friend Segment;
        friend Point;
    };

    class Ray : public Shape {
    private:
        Point origin;
        Vector direction;
    public:
        Ray(const Point & origin, const Point & point_of_ray);
        Ray(const Point & origin, const Vector & direction);
        Ray & move(const Vector & vector) override;
        bool has_point(const Point & point) const override;
        bool has_intarsection_with(const Segment & segment) const override; 

        friend std::ostream & operator << (std::ostream & os, const Ray & ray);
        friend std::istream & operator >> (std::istream & is, Ray & ray);

        friend Line;
        friend Point;
        friend Polygon;
    };

    class Polygon : public Shape {
    private:
        Point* points;
        int points_cnt;
        inline Polygon & define_points_cnt(const int & number_of_points);
        inline Polygon & copy_points(const int & number_of_points, const Point * points);
    public:
        Polygon();
        Polygon(const Polygon & other);
        Polygon(Polygon && other);
        Polygon(const int & number_of_points, const Point * points);
        Polygon(const int number_of_points, Point point_0, ...);
        Polygon(const int & number_of_points);
        ~Polygon();
        Polygon & operator = (const Polygon & other);
        Polygon & operator = (Polygon && other);
        Polygon set_point_cnt(const int & point_cnt);
        Polygon & add_point(const Point & point);
        Polygon & remove_point();
        int get_points_cnt() const;
        friend scalar_t area(const Polygon & polygon);
        bool is_convex() const;
        bool has_self_intersection() const;
        bool is_star() const;
    private:
        //functions for convex_hull
        friend int compare_vectors(const void * vector_ptr_0, const void * vector_ptr_1);
        friend int origin_index(const Point * points, const int points_cnt);
        
    public:
        friend Polygon convex_hull(const Polygon & polygon);

        Point operator [] (const int index) const;
        Point & operator [] (const int index);

        Polygon & move(const Vector & vector) override;
        bool has_point(const Point & point) const override; 
        bool has_intarsection_with(const Segment & segment) const override; 

        friend std::ostream & operator << (std::ostream & os, const Polygon & polygon);
        friend std::istream & operator >> (std::istream & is, Polygon & polygon);
    };
}
#endif // !_NickwideN_Geometry_H
    
template<typename user_t>
inline Geometry::Vector::Vector(const user_t coor_0, const user_t coor_1, const user_t coor_2, const user_t coor_3, ...) {
    user_t *pointer_coor = &coor_0;
    for (int i = 0; i < DIMENTION; ++i) {
        coordinates[i] = *pointer_coor;
        ++pointer_coor;
    }
}

