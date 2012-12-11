struct Complex {

    Complex(double _re, double _im)
    {
        re = _re;
        im = _im;
    }

    Complex& operator+=(const Complex& a)
    {
        re += a.re;
        im += a.im;
    }

    Complex operator*(const Complex& a) const
    {
        return Complex(re*a.re-im*a.im, re*a.im+im*a.re);
    }

    double real() const { return re; } 
    double imag() const { return im; }

    double re, im;
};

#if 0
inline Complex sqr1(const Complex& a)
{
    return Complex(a.real() * a.real() - a.imag() * a.imag(),
      2 * a.real() * a.imag());
}

inline Complex sqr2(const Complex& a)
{
    return Complex(a.re * a.re - a.im * a.im, 2 * a.re * a.im);
}

inline Complex sqr3(const Complex& a)
{
    double r = a.re;
    double i = a.im;
    return Complex(r*r-i*i, 2*r*i);
}
#endif

void foo(Complex& a, const Complex& b)
{
    a += sqr(b);
}

void foo2(Complex* __restrict__ a, Complex* __restrict__ b, int N)
{
    for (int i=0; i < N; ++i)
        a[i] += sqr(b[i]);
}

