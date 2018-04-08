#include"Geometry.h"
#include<iostream>
#include<iomanip>
using namespace Geometry;
using namespace std;
// F
/*
int main() {
    Segment a;
    Segment b;
    cin >> a;
    cin >> b;
    cout << (a.has_intarsection_with(b) ? "YES" : "NO");
    getchar();
    getchar();
    return 0;
}
*/

// G
/*
int main() {
    Segment a;
    Segment b;
    cin >> a;
    cin >> b;
    cout << setprecision(9) << fixed << distance_between(a, b);
    getchar();
    getchar();
    return 0;
}
*/
int main() {
    Point C;
    Point A;
    Point B;
    cin >> C >> A >> B;
    cout << setprecision(9) << fixed;
    cout << C.distance_to(Line(A, B)) << '\n';
    cout << C.distance_to(Ray(A, B)) << '\n';
    cout << C.distance_to(Segment(A, B)) << '\n';
    getchar();
    getchar();
    return 0;
}

/*
int main() {
    Polygon a;
    cin >> a;
    getchar();
    getchar();
    return 0;
}
*/