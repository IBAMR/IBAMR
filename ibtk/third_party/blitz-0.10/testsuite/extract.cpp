#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

class RGB24 {
public:
    static const int red=0, green=1, blue=2;

    RGB24() { }

    RGB24(char r, char g, char b)
    { 
        colours[0] = r;
        colours[1] = g;
        colours[2] = b;
    }

    unsigned char operator[](int component) const
    { return colours[component]; }

    bool operator==(const RGB24& x)
    { 
        return colours[0] == x.colours[0]
            && colours[1] == x.colours[1]
            && colours[2] == x.colours[2];
    }

private:
    unsigned char colours[3];
};

// declare static data members of RGB24 so we can use them below
const int RGB24::red;
const int RGB24::green;
const int RGB24::blue;

BZ_DECLARE_MULTICOMPONENT_TYPE(RGB24, unsigned char, 3)

BZ_DECLARE_ARRAY_ET_SCALAR_OPS(RGB24)

int main()
{
    Array<TinyVector<int,3>,2> A(8,8);
    A = TinyVector<int,3>(0);
    Array<int,2> B = A(Range(2,5),Range(2,5)).extractComponent(int(),1,3);
    B = 1;
    Array<int,2> C = A.extractComponent(int(),2,3);
    C = 2;

    BZTEST(A(1,1)[0] == 0 && A(1,1)[1] == 0 && A(1,1)[2] == 2);
    BZTEST(A(3,3)[0] == 0 && A(3,3)[1] == 1 && A(3,3)[2] == 2);

    Array<int,2> comp1 = A.extractComponent(int(),0,3);
    Array<int,2> comp2 = A.extractComponent(int(),1,3);
    Array<int,2> comp3 = A.extractComponent(int(),2,3);

    BZTEST(count(comp1 == 0) == 64);
    BZTEST(count(comp2 == 1) == 16);
    BZTEST(count(comp2 == 0) == 64-16);
    BZTEST(count(comp3 == 2) == 64);

#ifdef BZ_HAVE_COMPLEX
    Array<complex<float>,2> D(4,4);
    real(D) = 5;
    imag(D) = 6;

    BZTEST(count(D == complex<float>(5,6)) == 16);

    BZ_USING_NAMESPACE(blitz::tensor)

    D = zip(i, j, complex<float>());
    BZTEST(D(1,2) == complex<float>(1,2));
    BZTEST(D(3,0) == complex<float>(3,0));

    D = zip(real(D), -imag(D), complex<float>());
    BZTEST(D(1,2) == complex<float>(1,-2));
    BZTEST(D(3,0) == complex<float>(3,0));
#endif

    Array<RGB24,2> E(5,5);
    E[RGB24::red] = 30;
    E[RGB24::green] = 50;
    E[RGB24::blue] = 15;
    E[RGB24::green](Range(1,2),Range(1,2)) = 0;

    RGB24 z(30,50,15);
    BZTEST(count(E == z) == 21);
    RGB24 y(30,0,15);
    BZTEST(count(E == y) == 4);

    return 0;
}

