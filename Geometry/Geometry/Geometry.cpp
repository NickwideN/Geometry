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

Geometry::Vector::Vector(const Point & point_0, const Point & point_1) {
    *this = point_1.radius_vector - point_0.radius_vector;
}

Geometry::Vector & Geometry::Vector::operator+=(const Vector & vector) {
    for (int i = 0; i < DIMENTION; ++i)
        this->coordinates[i] += vector.coordinates[i];
    return *this;
}

Geometry::Vector Geometry::operator + (const Vector & vector_0, const Vector & vector_1){
    Geometry::Vector new_vector(vector_0);
    return new_vector += vector_1;
}

Geometry::Vector & Geometry::Vector::operator -= (const Vector & vector){
    return *this += -vector;
}

Geometry::Vector Geometry::operator - (const Vector & vector_0, const Vector & vector_1) {
    Geometry::Vector new_vector(vector_0);
    return new_vector -= vector_1;
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

Geometry::Vector & Geometry::Vector::operator *= (const scalar_t & coefficient) {
    for (int i = 0; i < DIMENTION; ++i)
        this->coordinates[i] *= coefficient;
    return *this;
}

Geometry::Vector Geometry::operator * (const Vector & vector, const scalar_t & coefficient) {
    Geometry::Vector new_vector(vector);
    return new_vector *= coefficient;
}

Geometry::Vector & Geometry::Vector::operator /= (const scalar_t & coefficient){
    return *this *= 1 / coefficient;
}

Geometry::Vector Geometry::operator / (const Vector & vector, const scalar_t & coefficient){
    Vector new_vector(vector);
    return new_vector /= coefficient;
}

Geometry::scalar_t Geometry::operator * (const Vector & vector_0, const Vector & vector_1){
    return scalar_product(vector_0, vector_1);
}

Geometry::scalar_t Geometry::scalar_product(const Vector & vector_0, const Vector & vector_1){
    Geometry::scalar_t product = 0;
    for (int i = 0; i < DIMENTION; ++i)
        product += vector_0[i] * vector_1[i];
    return product;
}

Geometry::Vector Geometry::vector_product(const Vector & vector_0, const Vector & vector_1) {
    if (DIMENTION != 3)
        error("Error: not able to take vector product with DIMENTION != 3. For getting an area of parallelogram on the vectors use function \"skew_product(Vector vector1, Vector vector2)\"");
    Vector new_vector;
    new_vector[0] = vector_0[1] * vector_1[2] - vector_0[2] * vector_1[1];
    new_vector[1] = -vector_0[0] * vector_1[2] + vector_0[2] * vector_1[0];
    new_vector[2] = vector_0[0] * vector_1[1] - vector_0[1] * vector_1[0];
    return new_vector;
}

Geometry::scalar_t Geometry::skew_product(const Vector & vector_0, const Vector & vector_1){
    if (DIMENTION != 2)
        error("Error : not able to take an area with DIMENTION != 2");
    return vector_0[0] * vector_1[1] - vector_0[1] * vector_1[0];;
}

Geometry::coordinate_t Geometry::Vector::operator[] (const int index) const{
    if (index >= DIMENTION)
        error("Error: attempt to take not exist coordinate ", index, ". DIMENTION = ", DIMENTION);
    return this->coordinates[index];
}

Geometry::coordinate_t & Geometry::Vector::operator[] (const int index) {
    if (index >= DIMENTION)
        error("Error: attempt to take not exist coordinate ", index, ". DIMENTION = ", DIMENTION);
    return this->coordinates[index];
}

Geometry::scalar_t Geometry::abs(const Vector & vector){
    scalar_t abs = 0;
    for (int i = 0; i < DIMENTION; ++i)
        abs += vector[i] * vector[i];
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

double Geometry::agl(const Vector & vector_0, const Vector & vector_1){
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

bool Geometry::are_complanar(const Vector & vector_0, const Vector & vector_1, const Vector & vector_2){
    if (DIMENTION > 2)
        error("Error : not able to do complanarity test with DIMENTION != 2");
   return true;
}

Geometry::scalar_t Geometry::length(const Point & point_0, const Point & point_1) {
    return abs(Vector(point_0, point_1));
}

std::ostream & Geometry::operator << (std::ostream & os, const Vector & vector) {
    os << '(';
    for (int i = 0; i < DIMENTION - 1; ++i)
        os << vector[i] << ", ";
    os << vector[DIMENTION - 1] << ')';
    return os;
}

std::istream & Geometry::operator >> (std::istream & is, Vector & vector) {
    for (int i = 0; i < DIMENTION; ++i)
        is >> vector[i];
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

Geometry::Point & Geometry::Point::move(const Vector & vector){
    this->radius_vector += vector;
    return *this;
}

bool Geometry::Point::has_point(const Point & point) const{
    return are_coincident(this->radius_vector, point.radius_vector);
}

bool Geometry::Point::has_intarsection_with(const Segment & segment) const{
    return segment.has_point(*this);
}

std::ostream & Geometry::operator<<(std::ostream & os, const Point & point){
    os << point.radius_vector;
    return os;
}

std::istream & Geometry::operator>>(std::istream & is, Point & point){
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

Geometry::Segment & Geometry::Segment::move(const Vector & vector){
    this->point_0.move(vector);
    this->point_1.move(vector);
    return *this;
}

bool Geometry::Segment::has_point(const Point & point) const{
    return (std::abs(skew_product(this->point_0.radius_vector, point.radius_vector)) + 
                std::abs(skew_product(point.radius_vector, this->point_1.radius_vector)) == 
                    std::abs(skew_product(this->point_0.radius_vector, this->point_1.radius_vector)));
}

bool Geometry::Segment::has_intarsection_with(const Segment & segment) const{
    return Line(*this).has_intarsection_with(segment) && Line(segment).has_intarsection_with(*this);
}

Geometry::scalar_t Geometry::length(const Segment & segment){
    return length(segment.point_0, segment.point_1);
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Class Line: public Shape {};

Geometry::Line::Line(const Segment & segment) :
    origen(segment.point_0),
    direction(Vector(segment.point_0, segment.point_0)) {
}

Geometry::Line & Geometry::Line::move(const Vector & vector){
    this->origen.move(vector);
    return *this;
}

bool Geometry::Line::has_point(const Point & point) const {
    return are_collinear(direction, origen.radius_vector - point.radius_vector);
}

bool Geometry::Line::has_intarsection_with(const Segment & segment) const {
    if (are_skew(*this, Line(segment)))
        return false;
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
