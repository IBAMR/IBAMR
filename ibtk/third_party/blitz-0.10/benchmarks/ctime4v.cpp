
#if defined(__GNUC__) && (__GNUC__ < 3)
#ifdef BZ_HAVE_STD
#include <cmath>
#else
#include <math.h>
#endif
inline float cos(float x) { return static_cast<float>(cos(x)); }
inline float sin(float x) { return static_cast<float>(sin(x)); }
inline float tan(float x) { return static_cast<float>(tan(x)); }
inline float log(float x) { return static_cast<float>(log(x)); }
inline float exp(float x) { return static_cast<float>(exp(x)); }
#endif
#include <valarray>


using namespace std;

int main()
{
}


void foo()
{
    const int N = 10;
    valarray<float> A(N), B(N), C(N), D(N), E(N);
    A = B*C + D*E;
    A = B+C + D*cos(E);
    A = B*sin(C) + D*sin(E);
    B = C + D;
    A = A+B+C+D;
    A = B*C*C + D*E;
    A = B-C + log(D)*cos(E);
    A = B*sin(C) + tan(D)/E;
    B = C - D;
    A = A+B*C+D;
}

