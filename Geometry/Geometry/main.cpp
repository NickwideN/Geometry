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