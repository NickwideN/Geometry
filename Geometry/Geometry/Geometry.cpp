#include"Geometry.h"
#include<iostream>
#include<cmath> //sqrt
#include<algorithm> //min, max
#include<cstring>    //strcmp

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Class Vector {};
Geometry::Vector::Vector() {
    for (int i = 0; i < DIMENTION; ++i) {
        this->coordinates[i] = default_value;
    }
}

Geometry::Vector::Vector(const coordinate_t coor_0, const coordinate_t coor_1, const coordinate_t coor_2) {
    coordinates[0] = coor_0;
    if (1 < DIMENTION) {
        this->coordinates[1] = coor_1;
    }
    if (2 < DIMENTION) {
        this->coordinates[2] = coor_2;
    }
    for (int i = 3; i < DIMENTION; ++i) {
        this->coordinates[i] = default_value;
    }
}

Geometry::Vector::Vector(const Point & point_0, const Point & point_1) {
    *this = point_1.radius_vector - point_0.radius_vector;
}

Geometry::Vector & Geometry::Vector::operator+=(const Vector & vector) {
    for (int i = 0; i < DIMENTION; ++i) {
        this->coordinates[i] += vector.coordinates[i];
    }
    return *this;
}

Geometry::Vector Geometry::operator + (const Vector & vector_0, const Vector & vector_1) {
    Vector new_vector(vector_0);
    return new_vector += vector_1;
}

Geometry::Vector & Geometry::Vector::operator -= (const Vector & vector){
    return *this += -vector;
}

Geometry::Vector Geometry::operator - (const Vector & vector_0, const Vector & vector_1) {
    Vector new_vector(vector_0);
    return new_vector -= vector_1;
}

Geometry::Vector Geometry::Vector::operator + () const{
    return *this;
}

Geometry::Vector Geometry::Vector::operator - () const{
    Vector new_vector;
    for (int i = 0; i < DIMENTION; ++i) {
        new_vector[i] = -this->coordinates[i];
    }
    return new_vector;
}

Geometry::Vector & Geometry::Vector::operator *= (const scalar_t & coefficient) {
    for (int i = 0; i < DIMENTION; ++i) {
        this->coordinates[i] *= coefficient;
    }
    return *this;
}

Geometry::Vector Geometry::operator * (const Vector & vector, const scalar_t & coefficient) {
    Geometry::Vector new_vector(vector);
    return new_vector *= coefficient;
}

Geometry::Vector & Geometry::Vector::operator /= (const scalar_t & coefficient) {
    return *this *= 1 / coefficient;
}

Geometry::Vector Geometry::operator / (const Vector & vector, const scalar_t & coefficient) {
    Vector new_vector(vector);
    return new_vector /= coefficient;
}

Geometry::Vector Geometry::operator * (const scalar_t & coefficient, const Vector & vector) {
    return vector * coefficient;
}

Geometry::Vector Geometry::operator / (const scalar_t & coefficient, const Vector & vector) {
    return vector / coefficient;
}

Geometry::scalar_t Geometry::operator * (const Vector & vector_0, const Vector & vector_1) {
    return scalar_product(vector_0, vector_1);
}

Geometry::scalar_t Geometry::scalar_product(const Vector & vector_0, const Vector & vector_1) {
    Geometry::scalar_t product = 0;
    for (int i = 0; i < DIMENTION; ++i) {
        product += vector_0[i] * vector_1[i];
    }
    return product;
}

Geometry::Vector Geometry::vector_product(const Vector & vector_0, const Vector & vector_1) {
    if (DIMENTION != 3) {
        throw "There is no opportunity to take vector product with DIMENTION != 3. For getting an area of parallelogram\
 on the vectors use function \"skew_product(Vector vector1, Vector vector2)\"";
    }
    Vector new_vector;
    new_vector[0] = vector_0[1] * vector_1[2] - vector_0[2] * vector_1[1];
    new_vector[1] = -vector_0[0] * vector_1[2] + vector_0[2] * vector_1[0];
    new_vector[2] = vector_0[0] * vector_1[1] - vector_0[1] * vector_1[0];
    return new_vector;
}

Geometry::scalar_t Geometry::skew_product(const Vector & vector_0, const Vector & vector_1) {
    // TODO: insert the code for DIMENTION > 2
    if (!(DIMENTION == 2 || DIMENTION == 3)) {
        throw "There is no opportunity to take an area with 1 >= DIMENTION >= 4 ";
    }
    if (DIMENTION == 2) {
        return vector_0[0] * vector_1[1] - vector_0[1] * vector_1[0];
    } else { // DIMENTION == 3
        return abs(vector_product(vector_0, vector_1));
    }
}

Geometry::coordinate_t Geometry::Vector::operator[] (const int index) const {
    if (index >= DIMENTION) {
        throw "Attempt to take not exist coordinate of vector";
    }
    return this->coordinates[index];
}

Geometry::coordinate_t & Geometry::Vector::operator[] (const int index) {
    if (index >= DIMENTION) {
        throw "Attempt to take not exist coordinate of vector";
    }
    return this->coordinates[index];
}

Geometry::scalar_t Geometry::abs(const Vector & vector){
    scalar_t abs = 0;
    for (int i = 0; i < DIMENTION; ++i) {
        abs += vector[i] * vector[i];
    }
    abs = sqrt(abs);
    return abs;
}

double Geometry::sin_agl(const Vector & vector_0, const Vector & vector_1) {
    return (std::abs(skew_product(vector_0, vector_1)) / abs(vector_0) / abs(vector_1));
}

double Geometry::cos_agl(const Vector & vector_0, const Vector & vector_1) {
    return ((vector_0 * vector_1) / abs(vector_0) / abs(vector_1));
}

double Geometry::tan_agl(const Vector & vector_0, const Vector & vector_1) {
    return (sin_agl(vector_0, vector_1) / cos_agl(vector_0, vector_1));
}

double Geometry::agl(const Vector & vector_0, const Vector & vector_1) {
    return acos(cos_agl(vector_0, vector_1));
}

bool Geometry::are_collinear(const Vector & vector_0, const Vector & vector_1) {
    return (skew_product(vector_0, vector_1) == 0);
}

bool Geometry::are_coincident(const Vector & vector_0, const Vector & vector_1) {
    bool are_coinc = true;
    for (int i = 0; i < DIMENTION; ++i) {
        if (vector_0[i] != vector_1[i]) {
            are_coinc = false;
            break;
        }
    }
    return are_coinc;
}

bool Geometry::are_complanar(const Vector & vector_0, const Vector & vector_1, const Vector & vector_2) {
    // TODO: insert the code for DIMENTION > 2
    if (DIMENTION > 2) {
        throw "There is no opportunity to do complanarity test with DIMENTION != 2";
    }
   return true;
}

std::ostream & Geometry::operator << (std::ostream & os, const Vector & vector) {
    os << '(';
    for (int i = 0; i < DIMENTION - 1; ++i) {
        os << vector[i] << ", ";
    }
    os << vector[DIMENTION - 1] << ')';
    return os;
}

std::istream & Geometry::operator >> (std::istream & is, Vector & vector) {
    for (int i = 0; i < DIMENTION; ++i) {
        is >> vector[i];
    }
    return is;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Class Point: public Shape {};

Geometry::Point::Point() :
    radius_vector() {
}

Geometry::Point::Point(const coordinate_t & coor_0, const coordinate_t & coor_1) :
    radius_vector(coor_0, coor_1) {
}

Geometry::Point::Point(const Vector ragius_vector) :
    radius_vector(radius_vector) {
}

Geometry::Point & Geometry::Point::move(const Vector & vector) {
    this->radius_vector += vector;
    return *this;
}

bool Geometry::Point::has_point(const Point & point) const{
    return are_coincident(this->radius_vector, point.radius_vector);
}

bool Geometry::Point::has_intarsection_with(const Segment & segment) const {
    return segment.has_point(*this);
}

Geometry::scalar_t Geometry::length(const Point & point_0, const Point & point_1) {
    return abs(Vector(point_0, point_1));
}

std::ostream & Geometry::operator << (std::ostream & os, const Point & point) {
    os << point.radius_vector;
    return os;
}

std::istream & Geometry::operator >> (std::istream & is, Point & point) {
    is >> point.radius_vector;
    return is;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Class Segment: public Shape {};

Geometry::Segment::Segment() :
    point_0(), point_1() {
}

Geometry::Segment::Segment(const Point & point_0, const Point & point_1) :
    point_0(point_0), point_1(point_1) {
}

Geometry::Segment::Segment(const Point & origen, const Vector & direction, const scalar_t & length)
    : point_0(origen) {
    this->point_1 = origen.radius_vector + direction / abs(direction) * length;
}

Geometry::Segment & Geometry::Segment::move(const Vector & vector) {
    this->point_0.move(vector);
    this->point_1.move(vector);
    return *this;
}

bool Geometry::Segment::has_point(const Point & point) const {
    return (std::abs(skew_product(this->point_0.radius_vector, point.radius_vector)) + 
                std::abs(skew_product(point.radius_vector, this->point_1.radius_vector)) == 
                    std::abs(skew_product(this->point_0.radius_vector, this->point_1.radius_vector)));
}

bool Geometry::Segment::has_intarsection_with(const Segment & segment) const {
    return Line(*this).has_intarsection_with(segment) && Line(segment).has_intarsection_with(*this);
}

Geometry::scalar_t Geometry::length(const Segment & segment) {
    return length(segment.point_0, segment.point_1);
}

std::ostream & Geometry::operator << (std::ostream & os, const Segment & segment) {
    os << "point_0: " << segment.point_0 << " point_1: " << segment.point_1;
    return os;
}

std::istream & Geometry::operator >> (std::istream & is, Segment & segmant) {
    is >> segmant.point_0 >> segmant.point_1;
    return is;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Class Line: public Shape {};

Geometry::Line::Line(const Point & point_0, const Point & point_1) :
    origen(point_0), direction(Vector(point_0, point_1)) {
}

Geometry::Line::Line(const Point & origen, const Vector & vector, const char * name_vector) :
    origen(origen) {
    if (!strcmp(name_vector, "direction")) {
        this->direction = vector;
    } else if (!strcmp(name_vector, "normal")) {
        if (DIMENTION == 2) {
            direction = Vector(-vector[1], vector[0]);
        } else {
            throw "There is no opportunity to do define normal vector for line with DIMENTION != 2";
        }
    } else {
        throw "Not correct name of vector. You can choose names \"direction\" (default) or \"vector\"";
        this->direction = vector;
    }
}

Geometry::Line::Line(const coordinate_t & coefficient_of_x, const coordinate_t & coefficient_of_y, const coordinate_t & absolute_term) {
    if (DIMENTION != 2) {
        throw "There is no able to create line use equation Ax + By + C = 0 with DIMENTION != 2";
    }
    direction = Vector(-coefficient_of_y, coefficient_of_x);
    origen = Point(0, -absolute_term/coefficient_of_y);
}

Geometry::Line::Line(const Segment & segment) :
    origen(segment.point_0),
    direction(Vector(segment.point_0, segment.point_0)) {
}

Geometry::Line & Geometry::Line::move(const Vector & vector) {
    this->origen.move(vector);
    return *this;
}

bool Geometry::Line::has_point(const Point & point) const {
    return are_collinear(direction, origen.radius_vector - point.radius_vector);
}

bool Geometry::Line::has_intarsection_with(const Segment & segment) const {
    if (are_skew(*this, Line(segment))) {
        return false;
    }
    scalar_t area_0 = skew_product(segment.point_0.radius_vector, this->direction);
    scalar_t area_1 = skew_product(segment.point_1.radius_vector, this->direction);
    return area_0 * area_1 <= 0;
}

bool Geometry::are_coincident(const Line & line_0, const Line & line_1) {
    return are_parallel(line_0, line_1) && 
        are_collinear(line_0.origen.radius_vector - line_1.origen.radius_vector, line_0.direction);
}

bool Geometry::are_parallel(const Line & line_0, const Line & line_1) {
    return are_collinear(line_0.direction, line_1.direction);;
}

bool Geometry::are_intersecting(const Line & line_0, const Line & line_1){
    return are_complanar(line_0.direction, line_1.direction, line_0.origen.radius_vector - line_1.origen.radius_vector) && 
        !are_parallel(line_0, line_1);
}

bool Geometry::are_skew(const Line & line_0, const Line & line_1) {
    return !are_complanar(line_0.direction, line_1.direction, line_0.origen.radius_vector - line_1.origen.radius_vector);
}

Geometry::Point Geometry::intersection(const Line & line_0, const Line & line_1) {
    //TODO: insert solution for DIMENTION > 2;
    if (DIMENTION > 3) {
        throw "There is no opportunity to find point of intersection for lines with DIMENTION != 2";
    }
    if (!are_intersecting(line_0, line_1)) {
        throw "The lines are not intersecting";
    }
    // r_intersection = r_0 + at  (t is coefficient, r_0 is origen, a is direction)
    scalar_t coefficient = (line_1.direction[1] * (line_0.origen.radius_vector[0] - line_1.origen.radius_vector[0])
        - line_1.direction[1] * (line_0.origen.radius_vector[1] - line_1.origen.radius_vector[1])) /
        (line_0.direction[0] * line_1.direction[1] - line_0.direction[1] * line_1.direction[0]);
    return line_0.origen.radius_vector + line_0.direction * coefficient;
}

Geometry::scalar_t Geometry::distance_between_parallel(const Line & line_0, const Line & line_1) {
    if (!are_parallel(line_0, line_1)) {
        throw "The lines are not parallel";
    }
    return abs(skew_product(line_0.origen.radius_vector - line_1.origen.radius_vector, line_0.direction) / abs(line_0.direction));
}

std::ostream & Geometry::operator << (std::ostream & os, const Line & line) {
    os << "origen: " << line.origen << " direction: " << line.direction;
    return os;
}

std::istream & Geometry::operator >> (std::istream & is, Line & line) {
    is >> line.origen >> line.direction;
    return is;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Class Ray: public Shape {};

Geometry::Ray::Ray(const Point & origen, const Vector & direction) :
    origen(origen), direction(direction) {
}

Geometry::Ray & Geometry::Ray::move(const Vector & vector) {
    this->origen.move(vector);
    return *this;
}

bool Geometry::Ray::has_point(const Point & point) const { ///////////////////////////////////////////////////////////////////
    scalar_t check_product = skew_product(origen.radius_vector, origen.radius_vector + direction);
    return false;
}

bool Geometry::Ray::has_intarsection_with(const Segment & segment) const ////////////////////////////////////////////////////////
{
    return false;
}

std::ostream & Geometry::operator << (std::ostream & os, const Ray & ray) {
    os << "origen: " << ray.origen << " direction :" << ray.direction;
    return os;
}

std::istream & Geometry::operator >> (std::istream & is, Ray & ray) {
    is >> ray.origen >> ray.direction;
    return is;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Class Polygon: public Shape {};

Geometry::Polygon & Geometry::Polygon::move(const Vector & vector) {
    for (int i = 0; i < this->points_cnt; ++i) {
        this->points[i].move(vector);
    }
    return *this;
}

bool Geometry::Polygon::has_point(const Point & point) const {
    bool poligon_contain_point = false;

    return false;
}

bool Geometry::Polygon::has_intarsection_with(const Segment & segment) const ///////////////////////////////////////////////////////////
{
    return false;
}

std::ostream & Geometry::operator << (std::ostream & os, const Polygon & polygon) {
    for (int i = 0; i < polygon.points_cnt; ++i) {
        os << "point_" << i << ": " << polygon.points[i] << " ";
    }
    return os;
}

std::istream & Geometry::operator >> (std::istream & is, Polygon & polygon) {
    for (int i = 0; i < polygon.points_cnt; ++i) {
        is >> polygon.points[i];
    }
    return is;
}