#include"Geometry.h"
#include<iostream>
#include<cmath> //sqrt
#include<algorithm> //min, max
Geometry::Vector::Vector() {
    for (int i = 0; i < DIMENTION; ++i)
        this->coordinates[i] = default_value;
}

Geometry::Vector::Vector(const coordinate_t coor_0, const coordinate_t coor_1, const coordinate_t coor_2){ //let put a several types
    coordinates[0] = coor_0;
    if (1 < DIMENTION)
        this->coordinates[1] = coor_1;
    if (2 < DIMENTION)
        this->coordinates[2] = coor_2;
    for (int i = 3; i < DIMENTION; ++i) {
        this->coordinates[i] = default_value;
    }
}

Geometry::Vector & Geometry::Vector::operator+=(const Vector & vector) {
    for (int i = 0; i < DIMENTION; ++i)
        this->coordinates[i] += vector.coordinates[i];
    return *this;
}

Geometry::Vector Geometry::operator + (const Vector & vector_1, const Vector & vector_2){
    Geometry::Vector new_vector(vector_1);
    return new_vector += vector_2;
}

Geometry::Vector & Geometry::Vector::operator -= (const Vector & vector){
    return *this += -vector;
}

Geometry::Vector Geometry::operator - (const Vector & vector_1, const Vector & vector_2) {
    Geometry::Vector new_vector(vector_1);
    return new_vector -= vector_2;
}

Geometry::Vector Geometry::Vector::operator + () const{
    return *this;
}

Geometry::Vector Geometry::Vector::operator - () const{
    Vector new_vector;
    for (int i = 0; i < DIMENTION; ++i)
        new_vector = -this->coordinates[i];
    return new_vector;
}

Geometry::Vector & Geometry::Vector::operator *= (const coordinate_t & coefficient) {
    for (int i = 0; i < DIMENTION; ++i)
        this->coordinates[i] *= coefficient;
    return *this;
}


Geometry::Vector Geometry::operator * (const Vector & vector, const coordinate_t & coefficient) {
    Geometry::Vector new_vector(vector);
    return new_vector *= coefficient;
}

Geometry::coordinate_t Geometry::operator * (const Vector & vector_1, const Vector & vector_2){
    return scalar_product(vector_1, vector_2);
}

Geometry::coordinate_t Geometry::scalar_product(const Geometry::Vector & vector_1, const Geometry::Vector & vector_2){
    Geometry::coordinate_t result(0);
    for (int i = 0; i < Geometry::DIMENTION; ++i)
        result += vector_1[i] * vector_2[i];
    return result;
}

Geometry::Vector Geometry::vector_product(const Geometry::Vector & vector_1, const Geometry::Vector & vector_2) {
    if (Geometry::DIMENTION != 3)
        Geometry::error("Error: not able to take vector product with DIMENTION != 3. For getting an area of parallelogram on the vectors use function \"skew_product(Vector vector1, Vector vector2)\"");
    Geometry::Vector new_vector;
    new_vector[0] = vector_1[1] * vector_2[2] - vector_1[2] * vector_2[1];
    new_vector[1] = vector_1[0] * vector_2[2] - vector_1[2] * vector_2[0];
    new_vector[2] = vector_1[0] * vector_2[1] - vector_1[1] * vector_2[0];
    return new_vector;
}

Geometry::coordinate_t Geometry::skew_product(const Geometry::Vector & vector_1, const Geometry::Vector & vector_2){
    if (Geometry::DIMENTION != 2)
        Geometry::error("Error : not able to take an area with DIMENTION != 2");
    Geometry::coordinate_t area = 0;
    area = vector_1[0] * vector_2[1] - vector_1[1] * vector_2[0];
    return area;
}

Geometry::coordinate_t Geometry::Vector::operator[] (const int index) const{
    if (index >= Geometry::DIMENTION)
        error("Error: attempt to take not exist coordinate ", index, ". DIMENTION = ", Geometry::DIMENTION);
    return this->coordinates[index];
}

Geometry::coordinate_t & Geometry::Vector::operator[] (const int index) {
    if (index >= Geometry::DIMENTION)
        error("Error: attempt to take not exist coordinate ", index, ". DIMENTION = ", Geometry::DIMENTION);
    return this->coordinates[index];
}

double Geometry::abs(const Geometry::Vector & vector){
    double result = 0;
    for (int i = 0; i < Geometry::DIMENTION; ++i)
        result += vector[i] * vector[i];
    result = sqrt(result);
    return result;
}

double Geometry::sin_agl(const Geometry::Vector & vector_1, const Geometry::Vector & vector_2) {
    return (Geometry::skew_product(vector_1, vector_2) / Geometry::abs(vector_1) / Geometry::abs(vector_2));
}

double Geometry::cos_agl(const Geometry::Vector & vector_1, const Geometry::Vector & vector_2) {
    return ((vector_1 * vector_2) / Geometry::abs(vector_1) / Geometry::abs(vector_2));
}

double Geometry::tan_agl(const Geometry::Vector & vector_1, const Geometry::Vector & vector_2) {
    return (Geometry::sin_agl(vector_1, vector_2) / Geometry::cos_agl(vector_1, vector_2));
}

double Geometry::agl(const Geometry::Vector & vector_1, const Geometry::Vector & vector_2){
    return acos(Geometry::cos_agl(vector_1, vector_2));
}

std::ostream & Geometry::operator<<(std::ostream & os, const Geometry::Vector & vector){
    os << '(';
    for (int i = 0; i < Geometry::DIMENTION - 1; ++i)
        os << vector[i] << ", ";
    os << vector[Geometry::DIMENTION - 1] << ')';
    return os;
}

std::istream & Geometry::operator>>(std::istream & is, Geometry::Vector & vector){
    for (int i = 0; i < Geometry::DIMENTION; ++i)
        is >> vector[i];
    return is;
}

bool Geometry::are_collinear(const Geometry::Vector & vector_1, const Geometry::Vector & vector_2) {
    return (Geometry::skew_product(vector_1, vector_2) == 0);
}

bool Geometry::are_coincident(const Geometry::Vector & vector_1, const Geometry::Vector & vector_2) {
    bool are_coinc = true;
    for (int i = 0; i < Geometry::DIMENTION; ++i) {
        if (vector_1[i] != vector_2[i]) {
            are_coinc = false;
            break;
        }
    }
        return are_coinc;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Class Point: public Shape {};

Geometry::Point::Point():
    radius_vector() {}

Geometry::Point::Point(coordinate_t coor_0, coordinate_t coor_1):
    radius_vector(coor_0, coor_1){}

Geometry::Point & Geometry::Point::move(const Geometry::Vector & vector){
    this->radius_vector += vector;
    return *this;
}

bool Geometry::Point::has_point(const Point & point) const{
    return are_coincident(this->radius_vector, point.radius_vector);
}

bool Geometry::Point::has_intarsection_with(const Segment & segment) const{
    return segment.has_point(*this);
}

std::ostream & Geometry::operator<<(std::ostream & os, const Geometry::Point & point){
    os << point.radius_vector;
    return os;
}

std::istream & Geometry::operator>>(std::istream & is, Geometry::Point & point){
    is >> point.radius_vector;
    return is;
}

Geometry::Vector Geometry::Point::get_radius_vector() {
    return this->radius_vector;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Class Segment: public Shape {};

Geometry::Segment & Geometry::Segment::move(const Geometry::Vector & vector){
    this->point_1.move(vector);
    this->point_2.move(vector);
    return *this;
}

bool Geometry::Segment::has_point(const Point & point) const{
    bool this_segment_has_point = true;
    //Check, has Line(point_1, point_2) the point
    if (!Geometry::are_collinear(this->point_1.get_radius_vector - this->point_2.get_radius_vector, \
                                 this->point_1.get_radius_vector - point.get_radius_vector)) {
        this_segment_has_point = false;
    }
    //Check, has a rectangle with point_1 and point_2 in opposite angles
    for (int i = 0; i < DIMENTION; ++i) {
        if (!((std::min())) {
            this_segment_has_point = false;
        }
    }

    return this_segment_has_point;
}

bool Geometry::Segment::has_intarsection_with(const Segment & segment) const
{
    return false;
}

