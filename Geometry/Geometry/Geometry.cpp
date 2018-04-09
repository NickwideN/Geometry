#include"Geometry.h"
#include<iostream>
#include<cstring>    // strcmp
#include<cmath>      // sqrt
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
    Vector new_vector(vector);
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
    scalar_t product = 0;
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
    // The solution is presented for DIMENTION == (2 & 3);
    // TODO: insert solution for DIMENTION > 3;
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
    return (double(abs(skew_product(vector_0, vector_1))) / abs(vector_0) / abs(vector_1));
}

double Geometry::cos_agl(const Vector & vector_0, const Vector & vector_1) {
    return (double((vector_0 * vector_1)) / abs(vector_0) / abs(vector_1));
}

double Geometry::tan_agl(const Vector & vector_0, const Vector & vector_1) {
    return (double(sin_agl(vector_0, vector_1)) / cos_agl(vector_0, vector_1));
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
    // The solution is presented for DIMENTION == 2;
    // TODO: insert solution for DIMENTION > 2;
    if (DIMENTION > 2) {
        throw "There is no opportunity to do complanarity test with DIMENTION != 2";
    }
   return true;
}

bool Geometry::are_co_directed(const Vector & vector_0, const Vector & vector_1) {
    if (are_collinear(vector_0, vector_1)) {
        for (int i = 0; i < DIMENTION; ++i) {
            if (vector_0[i] * vector_1[i] < 0) {
                return false;
            }
        }
        return true;
    }
    return false;
}

std::ostream & Geometry::operator << (std::ostream & os, const Vector & vector) {
    for (int i = 0; i < DIMENTION; ++i) {
        os << vector[i] << " ";
    }
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

Geometry::Point::Point(const Vector radius_vector) :
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

Geometry::scalar_t Geometry::Point::distance_to(const Line & line) {
    Line parallel_line(*this, line.direction);
    return distance_between(line, parallel_line);
}

Geometry::scalar_t Geometry::Point::distance_to(const Ray & ray) {
    // The solution is presented for DIMENTION == 2;
    // TODO: insert solution for DIMENTION > 2;
    if (DIMENTION > 2) {
        throw "There is no opportunity to take distance from point to ray with DIMENTION > 2";
    }
    Line parallel_line(*this, ray.direction);
    Line normal_line(ray.origen, Vector(-ray.direction[1], ray.direction[0]));
    Point inter_point = intersection(parallel_line, normal_line);
    if (are_co_directed(Vector(inter_point, *this), ray.direction)) {
        return this->distance_to(Line(ray));
    } else {
        return length(*this, ray.origen);
    }
}

Geometry::scalar_t Geometry::Point::distance_to(const Segment & segment) {
    // The solution is presented for DIMENTION == 2;
    // TODO: insert solution for DIMENTION > 2;
    if (DIMENTION > 2) {
        throw "There is no opportunity to take distance from point to segment with DIMENTION > 2";
    }
    if (Line(*this, Vector(-Line(segment).direction[1], Line(segment).direction[0])).has_intarsection_with(segment)) {
        return this->distance_to(Line(segment));
    } else {
        scalar_t distance_to_0 = length(*this, segment.point_0);
        scalar_t distance_to_1 = length(*this, segment.point_1);
        return (distance_to_0 < distance_to_1 ? distance_to_0 : distance_to_1);
    }
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
    Segment tmp_seg = *this;
    Point tmp_point = point;
    if (skew_product(tmp_seg.point_0.radius_vector, tmp_seg.point_1.radius_vector) == 0) {
        Vector move_vector(1, 1);
        if (are_collinear(move_vector, Line(tmp_seg).direction)) {
            move_vector = Vector(1, 2);
        }
        tmp_point.move(move_vector);
        tmp_seg.move(move_vector);
    }
    return (std::abs(skew_product(tmp_seg.point_0.radius_vector, tmp_point.radius_vector)) + 
                std::abs(skew_product(tmp_point.radius_vector, tmp_seg.point_1.radius_vector)) == 
                    std::abs(skew_product(tmp_seg.point_0.radius_vector, tmp_seg.point_1.radius_vector)));
}

bool Geometry::Segment::has_intarsection_with(const Segment & segment) const {
    Line(*this).has_intarsection_with(segment);
    Line(segment).has_intarsection_with(*this);
    return Line(*this).has_intarsection_with(segment) && Line(segment).has_intarsection_with(*this);
}

Geometry::scalar_t Geometry::length(const Segment & segment) {
    return length(segment.point_0, segment.point_1);
}

Geometry::scalar_t Geometry::distance_between(const Segment & segment_0, const Segment & segment_1) {
    Point points[4] = { segment_0.point_0, segment_0.point_1, segment_1.point_0, segment_1.point_1 };
    scalar_t distance[4] = { points[0].distance_to(segment_1), points[1].distance_to(segment_1), 
        points[2].distance_to(segment_0), points[3].distance_to(segment_0) };
    scalar_t min_distance = distance[0];
    for (int i = 1; i < 4; ++i) {
        if (min_distance > distance[i]) {
            min_distance = distance[i];
        }
    }
    return min_distance;
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
    if (are_coincident(direction, Vector(0, 0))) {
        throw "There is no opportunity to creat a line with direction == vector(0)";
    }
}

Geometry::Line::Line(const Point & origen, const Vector & vector, const char * name_vector) :
    origen(origen) {
    if (!strcmp(name_vector, "direction")) {
        this->direction = vector;
        if (are_coincident(vector, Vector(0, 0))) {
            throw "There is no opportunity to creat a line with direction == vector(0)";
        }
    } else if (!strcmp(name_vector, "normal")) {
        if (DIMENTION == 2) {
            this->direction = Vector(-vector[1], vector[0]);
            if (are_coincident(vector, Vector(0, 0))) {
                throw "There is no opportunity to creat a line with direction == vector(0)";
            }
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
    this->direction = Vector(-coefficient_of_y, coefficient_of_x);
    if (are_coincident(this->direction, Vector(0, 0))) {
        throw "There is no opportunity to creat a line with direction == vector(0)";
    }
    if (coefficient_of_y) {
        origen = Point(0, -absolute_term / coefficient_of_y);
    }
    else {
        origen = Point(-absolute_term / coefficient_of_x, 0);
    }
}

Geometry::Line::Line(const Segment & segment) :
    origen(segment.point_0),
    direction(Vector(segment.point_0, segment.point_1)) {
    if (are_coincident(this->direction, Vector(0, 0))) {
        throw "There is no opportunity to creat a line with direction == vector(0)";
    }
}

Geometry::Line::Line(const Ray & ray) :
origen(ray.origen), direction(ray.direction) {
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
    scalar_t area_0 = skew_product(Vector(segment.point_0, this->origen), this->direction);
    scalar_t area_1 = skew_product(Vector(segment.point_1, this->origen), this->direction);
    return area_0 * area_1 <= 0;
}

Geometry::Vector Geometry::Line::get_direction() {
    return this->direction;
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
    // The solution is presented for DIMENTION == 2;
    // TODO: insert solution for DIMENTION > 2;
    if (DIMENTION > 3) {
        throw "There is no opportunity to find point of intersection for lines with DIMENTION != 2";
    }
    if (!are_intersecting(line_0, line_1)) {
        throw "The lines are not intersecting";
    }
    // r_intersection = r_0 + at  (t is coefficient, r_0 is origen, a is direction)
    if (!(line_0.direction[0] * line_1.direction[1] - line_0.direction[1] * line_1.direction[0])) {
        throw "Attempt of division by zero";
    }
    scalar_t coefficient = scalar_t(line_1.direction[0] * (line_0.origen.radius_vector[1] - line_1.origen.radius_vector[1])
        - line_1.direction[1] * (line_0.origen.radius_vector[0] - line_1.origen.radius_vector[0])) /
        (line_0.direction[0] * line_1.direction[1] - line_0.direction[1] * line_1.direction[0]);
    return line_0.origen.radius_vector + line_0.direction * coefficient;
}

Geometry::scalar_t Geometry::distance_between(const Line & line_0, const Line & line_1) {
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

Geometry::Ray::Ray(const Point & origen, const Point & point_of_ray) :
    origen(origen), direction(point_of_ray.radius_vector - origen.radius_vector) {
    if (are_coincident(this->direction, Vector(0, 0))) {
        throw "There is no opportunity to creat a segment with direction == vector(0)";
    }
}

Geometry::Ray::Ray(const Point & origen, const Vector & direction) :
    origen(origen), direction(direction) {
    if (are_coincident(this->direction, Vector(0, 0))) {
        throw "There is no opportunity to creat a segment with direction == vector(0)";
    }
}

Geometry::Ray & Geometry::Ray::move(const Vector & vector) {
    this->origen.move(vector);
    return *this;
}

bool Geometry::Ray::has_point(const Point & point) const { 
    return ((Line(*this).has_point(point)) && 
        (are_co_directed(this->direction, Vector(this->origen, point))));
}

bool Geometry::Ray::has_intarsection_with(const Segment & segment) const {
    Point inter_point = intersection(Line(segment), Line(*this));
    if (segment.has_point(inter_point) && this->has_point(inter_point)) {
        return true;
    }
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

Geometry::Polygon & Geometry::Polygon::define_points_cnt(const int & number_of_points) {
    this->points_cnt = number_of_points;
    this->points = new Point[number_of_points];
    return *this;
}

Geometry::Polygon & Geometry::Polygon::copy_points(const int & number_of_points, const Point * points) {
    for (int i = 0; i < number_of_points; ++i) {
        this->points[i] = points[i];
    }
    return *this;
}

Geometry::Polygon::Polygon() {
    this->define_points_cnt(1);
    this->points[0] = Point();
}

Geometry::Polygon::Polygon(const Polygon & other) {
    this->define_points_cnt(other.points_cnt);
    this->copy_points(this->points_cnt, other.points);
}

Geometry::Polygon::Polygon(Polygon && other) :
    points_cnt(other.points_cnt), points(other.points) {
    other.points_cnt = 0; other.points = 0;
}

Geometry::Polygon::Polygon(const int & number_of_points, const Point * points) {
    this->define_points_cnt(number_of_points);
    this->copy_points(number_of_points, points);
}

Geometry::Polygon::Polygon(const int number_of_points, Point point_0, ...) {
    Point * pointer_point = &point_0;
    this->define_points_cnt(number_of_points);
    for (int i = 0; i < number_of_points; ++i) {
        this->points[i] = *pointer_point;
        ++pointer_point;
    }
}

Geometry::Polygon::Polygon(const int & number_of_points) {
    this->define_points_cnt(number_of_points);
    for (int i = 0; i < this->points_cnt; ++i) {
        this->points[i] = Point();
    }
}

Geometry::Polygon::~Polygon() {
    delete[] this->points;
}

Geometry::Polygon & Geometry::Polygon::operator = (const Polygon & other) {
    Point * tmp_points = new Point[other.points_cnt];
    delete[] this->points;
    this->points = tmp_points;
    this->points_cnt = other.points_cnt;
    this->copy_points(this->points_cnt, other.points);
    return *this;
}

Geometry::Polygon & Geometry::Polygon::operator = (Polygon && other) {
    delete[] this->points;
    this->points = other.points;
    this->points_cnt = other.points_cnt;
    other.points = 0;
    other.points_cnt = 0;
    return *this;
}

Geometry::Polygon Geometry::Polygon::set_point_cnt(const int & point_cnt) {
    Polygon new_polygon(point_cnt);
    int min_point_cnt = (this->points_cnt < new_polygon.points_cnt ? this->points_cnt : new_polygon.points_cnt);
    new_polygon.copy_points(min_point_cnt, this->points);
    for (int i = min_point_cnt; i < new_polygon.points_cnt; ++i) {
        new_polygon.points[i] = this->points[this->points_cnt - 1];
    }
    delete[] this->points;
    return new_polygon;
}

Geometry::Polygon & Geometry::Polygon::add_point(const Point & point) {
    Polygon new_poligon(this->points_cnt + 1);
    new_poligon.copy_points(this->points_cnt, this->points);
    new_poligon.points[new_poligon.points_cnt - 1] = point;
    delete[] this->points;
    return new_poligon;
}

Geometry::Polygon & Geometry::Polygon::remove_point() {
    Polygon new_poligon(this->points_cnt - 1);
    new_poligon.copy_points(this->points_cnt - 1, this->points);
    delete[] this->points;
    return new_poligon;
}

int Geometry::Polygon::get_points_cnt() {
    return this->points_cnt;
}

bool Geometry::Polygon::is_convex() {
    if (this->points_cnt < 3) {
        return true;
    }
    scalar_t check_product = skew_product(Vector(this->points[this->points_cnt - 1], this->points[0]),
        Vector(this->points[0], this->points[1]));
    for (int i = 0; i < this->points_cnt - 2; ++i) {
        if ((check_product * skew_product(Vector(this->points[i], this->points[i + 1]),
            Vector(this->points[i + 1], this->points[i + 2]))) < 0) {
            return false;
        }
    }
    if ((check_product * skew_product(Vector(this->points[this->points_cnt - 2], this->points[this->points_cnt - 1]),
        Vector(this->points[this->points_cnt - 1], this->points[0]))) < 0) {
        return false;
    }
    return true;
}

Geometry::scalar_t Geometry::area(const Polygon & polygon) {
    scalar_t area = 0;
    for (int i = 0; i < polygon.points_cnt - 1; ++i) {
        area += skew_product(polygon.points[i].radius_vector, polygon.points[i + 1].radius_vector);
    }
    area += skew_product(polygon.points[polygon.points_cnt - 1].radius_vector, polygon.points[0].radius_vector);
    area = abs(area / 2);
    return area;
}

Geometry::Point Geometry::Polygon::operator [] (const int index) const {
    return this->points[index];
}

Geometry::Point & Geometry::Polygon::operator [] (const int index) {
    return this->points[index];
}

Geometry::Polygon & Geometry::Polygon::move(const Vector & vector) {
    for (int i = 0; i < this->points_cnt; ++i) {
        this->points[i].move(vector);
    }
    return *this;
}

bool Geometry::Polygon::has_point(const Point & point) const {
    // The solution is presented for cases when polygon is convex;
    // TODO: insert solution for cases when polygon is convex;
    /*if (!this->is_convex) {
        throw "There is no opportunity to study a polygon if it is convex";
    }*/
    scalar_t sum_area = 0;
    for (int i = 0; i < this->points_cnt - 1; ++i) {
        sum_area += abs(skew_product(Vector(point, this->points[i]), Vector(point, this->points[i + 1])));
    }
    sum_area += abs(skew_product(Vector(point, this->points[this->points_cnt - 1]), Vector(point, this->points[0])));
    sum_area /= 2;
    return area(*this) == sum_area;
}

bool Geometry::Polygon::has_intarsection_with(const Segment & segment) const {
    return (this->has_point(segment.point_0) || this->has_point(segment.point_1));
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