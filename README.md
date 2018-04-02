# Geometry

Class Vector:
-------------------------------
Класс Vector предназначен для работы с векторами в n-мерном евклидовом пространстве. Размерность такого пространства задается стачическим константным полем класса DIMENTION. Также есть возможность изменить тип координаты вектора.

#### Конструкторы:
```c++
Vector(); 
```    
Конструктор по умолчанию. Создает нулевой вектор.
```c++
template<typename user_t>
Vector(const user_t coor_0, const user_t coor_1, const user_t coor_2, const user_t coor_3, ...);
```
Конструктор принимающий более 3-х координат в качестве аргументов. Для корректной работы конструктора необходимо, чтобы все фактические параметры были одного типа и их количество равнялось размерности данного векторного пространства. 
  
```c++
Vector(const coordinate_t coor_0, const coordinate_t coor_1 = default_value, const coordinate_t coor_2 = default_value);
```
Конструктор, принимающий 1-3 параметров. Фактические параметры могут быть разных типов. При количестве аргументов менее DIMINATION  координаты, которые не инициализировались пользователем, становятся равными нулю.

#### Методы, операции, функции-члены
```c++
    Vector & operator += (const Vector & vector);
```
```c++
    friend Vector operator + (const Vector & vector_1, const Vector & vector_2);
```
```c++
    Vector operator + () const;
```
```c++
    Vector operator - () const;
```
```c++
    Vector & operator *= (const coordinate_t & coefficient);
```
```c++
    friend Vector operator * (const Vector & vector, const coordinate_t & coefficient);
```
```c++
    friend coordinate_t operator * (const Vector & vector_1, const Vector & vector_2);
```
Умножение векторов с помощью операции * возвращает результат скалярного перемножения.
```c++
    friend coordinate_t scalar_product (const Vector & vector_1, const Vector & vector_2);
```
```c++
    friend Vector vector_product (const Vector & vector_1, const Vector & vector_2);
```
Векторное перемножение работает только при размерности 3. В иных случаях вызывается исключение с выходом из программы.
```c++
    friend coordinate_t area(const Vector & vector_1, const Vector & vector_2);
```
area(vector_1, vector_2); возвращает площадь параллелограмма, основанного на векторах 1 и 2. Реализована только для двумерного пространства. В иных случаях вызывается исключение с выходом из программы.
```c++
    coordinate_t operator [] (const int index) const;
    coordinate_t & operator [] (const int index);
```
Доступ по индексу к координатам вектора.
```c++
    friend double abs (const Vector & vector);
```
abs(vector); возвращает модуль вектора vector.
```c++
    friend double agl(const Vector & vector_1, const Vector & vector_2);
```
agl(vector_1, vector_2); возвращает угол между векторами 1 и 2.
```c++
    friend std::ostream & operator << (std::ostream & os, const Vector & vector);
    friend std::istream & operator >> (std::istream & is, Vector & vector);
```
Ввод координат производиться вводом чисел через space-символ.
