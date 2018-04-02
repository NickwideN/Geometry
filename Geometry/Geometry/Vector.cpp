#include"Vector.h"
#include<iostream>
#include<cmath> //sqrt
Vector::Vector() {
    for (int i = 0; i < DIMENTION; ++i)
        this->coordinates[i] = default_value;
}

Vector::Vector(const coordinate_t coor_0, const coordinate_t coor_1, const coordinate_t coor_2){ //let put a several types
    coordinates[0] = coor_0;
    if (1 < DIMENTION)
        this->coordinates[1] = coor_1;
    if (2 < DIMENTION)
        this->coordinates[2] = coor_2;
    for (int i = 3; i < DIMENTION; ++i) {
        this->coordinates[i] = default_value;
    }
}

Vector & Vector::operator+=(const Vector & vector) {
    for (int i = 0; i < DIMENTION; ++i)
        this->coordinates[i] += vector.coordinates[i];
    return *this;
}

Vector operator+(const Vector & vector_1, const Vector & vector_2){
    Vector new_vector(vector_1);
    return new_vector += vector_2;
}

Vector Vector::operator+() const{
    return *this;
}

Vector Vector::operator-() const{
    Vector new_vector;
    for (int i = 0; i < DIMENTION; ++i)
        new_vector = -this->coordinates[i];
    return new_vector;
}

Vector & Vector::operator *= (const coordinate_t & coefficient) {
    for (int i = 0; i < DIMENTION; ++i)
        this->coordinates[i] *= coefficient;
    return *this;
}


Vector operator * (const Vector & vector, const Vector::coordinate_t & coefficient) {
    Vector new_vector(vector);
    return new_vector *= coefficient;
}

Vector::coordinate_t operator*(const Vector & vector_1, const Vector & vector_2){
    return scalar_product(vector_1, vector_2);
}

Vector::coordinate_t scalar_product(const Vector & vector_1, const Vector & vector_2){
    Vector::coordinate_t result(0);
    for (int i = 0; i < Vector::DIMENTION; ++i)
        result += vector_1.coordinates[i] * vector_2.coordinates[i];
    return result;
}

Vector vector_product(const Vector & vector_1, const Vector & vector_2) {
    Vector new_vector;
    try {
        if (Vector::DIMENTION != 3)
            throw 1;
        new_vector[0] = vector_1[1] * vector_2[2] - vector_1[2] * vector_2[1];
        new_vector[1] = vector_1[0] * vector_2[2] - vector_1[2] * vector_2[0];
        new_vector[2] = vector_1[0] * vector_2[1] - vector_1[1] * vector_2[0];
    }
    catch (int error_i) {
        std::cerr << "Error: not able to take vector product with DIMENTION != 3. For getting an area of parallelogram on the vectors use function \"area(Vector vector1, Vector vector2)\"" << std::endl;
        getchar();
        getchar();
        std::exit(1);
    }
    return new_vector;
}

Vector::coordinate_t area(const Vector & vector_1, const Vector & vector_2){
    Vector::coordinate_t area = 0;
    try {
        if (Vector::DIMENTION != 2)
            throw 2;
         area = vector_1[0] * vector_2[1] - vector_1[1] * vector_2[0];
    }
    catch (int error_i) {
        std::cerr << "Error : not able to take an area with DIMENTION != 2" << std::endl;
        getchar();
        getchar();
        std::exit(1);
    }
    return area;
}

Vector::coordinate_t Vector::operator[] (const int index) const{
    try {
        if (index >= Vector::DIMENTION)
            throw 3;
    }
    catch (int error_i) {
        std::cerr << "Error: attempt to take not exist coordinate " << index << ". DIMENTION = " << Vector::DIMENTION << std::endl;
        getchar();
        getchar();
        std::exit(1);
    }
    return this->coordinates[index];
}

Vector::coordinate_t & Vector::operator[] (const int index) {
    try {
        if (index >= Vector::DIMENTION)
            throw 4;
    }
    catch (int error_i) {
        std::cerr << "Error: attempt to change not exist coordinate " << index << ". DIMENTION = " << Vector::DIMENTION << std::endl;
        getchar();
        getchar();
        std::exit(1);
    }
    return this->coordinates[index];
}

double abs(const Vector & vector)
{
    double result = 0;
    for (int i = 0; i < Vector::DIMENTION; ++i)
        result += vector[i] * vector[i];
    result = sqrt(result);
    return result;
}

double agl(const Vector & vector_1, const Vector & vector_2){
    return acos((vector_1 * vector_2) / abs(vector_1) / abs(vector_2));
}


std::ostream & operator<<(std::ostream & os, const Vector & vector){
    std::cout << '(';
    for (int i = 0; i < Vector::DIMENTION - 1; ++i)
        std::cout << vector.coordinates[i] << ", ";
    std::cout << vector.coordinates[Vector::DIMENTION - 1] << ')';
    return os;
}

std::istream & operator>>(std::istream & is, Vector & vector){
    for (int i = 0; i < Vector::DIMENTION; ++i)
        std::cin >> vector.coordinates[i];
    return is;
}
