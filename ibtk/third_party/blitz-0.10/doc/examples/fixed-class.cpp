#include <blitz/array.h>
#include <blitz/numinquire.h> // for huge()

using namespace blitz;

// A simple fixed point arithmetic class which represents a point
// in the interval [0,1].
class FixedPoint {

public:
    typedef unsigned int T_mantissa;

    FixedPoint() { }

    FixedPoint(T_mantissa mantissa)
    {  
        mantissa_ = mantissa;
    }

    FixedPoint(double value)
    {
        assert((value >= 0.0) && (value <= 1.0));
        mantissa_ = static_cast<T_mantissa>(value * huge(T_mantissa()));
    }
   
    FixedPoint operator+(FixedPoint x)
    { return FixedPoint(mantissa_ + x.mantissa_); }

    double value() const
    { return mantissa_ / double(huge(T_mantissa())); }

private:
    T_mantissa mantissa_;
};

ostream& operator<<(ostream& os, const FixedPoint& a)
{
    os << a.value();
    return os;
}

int main()
{
    // Create an array using the FixedPoint class:

    Array<FixedPoint, 2> A(4,4), B(4,4);

    A = 0.5, 0.3, 0.8, 0.2,
        0.1, 0.3, 0.2, 0.9,
        0.0, 1.0, 0.7, 0.4,
        0.2, 0.3, 0.8, 0.4;

    B = A + 0.05;

    cout << "B = " << B << endl;

    return 0;
}
