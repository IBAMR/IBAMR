#include "testsuite.h"
#include <blitz/array.h>
#include <blitz/tinyvec2.h>

BZ_USING_NAMESPACE(blitz)

class crap {
public:
  double operator()(double x) const {return x+1;};
  BZ_DECLARE_FUNCTOR(crap);
};

template<typename T>
void test_1d(T& B)
{
  B = -2, 5, 4, 3, 1, 1;

    BZTEST(sum(B) == 12);
    BZTEST(mean(B) == 2);
    BZTEST(min(B) == -2);
    BZTEST(max(B) == 5);
    BZTEST(product(B) == -120);
    BZTEST(first(B>0) == 1);
    BZTEST(last(B>3) == 2);
    BZTEST(count(B>0) == 5);
    BZTEST(any(B>0));
    BZTEST(!any(B>10));
    BZTEST(!all(B>0));
    BZTEST(all(B>-10));

    BZTEST(all(maxIndex(B)==1));
    BZTEST(all(minIndex(B)==0));

    MinMaxValue<int> mm=minmax(B);
    BZTEST(mm.min==-2);
    BZTEST(mm.max==5);
}

// pass in an expression with value
// (-2, 5, 4, 3, 1, 1)
template<typename T>
void test_1dexpr(const T& B)
{
    BZTEST(sum(B) == 12);
    BZTEST(mean(B) == 2);
    BZTEST(min(B) == -2);
    BZTEST(max(B) == 5);
    BZTEST(product(B) == -120);
    BZTEST(first(B>0) == 1);
    BZTEST(last(B>3) == 2);
    BZTEST(count(B>0) == 5);
    BZTEST(any(B>0));
    BZTEST(!any(B>10));
    BZTEST(!all(B>0));
    BZTEST(all(B>-10));

    BZTEST(all(maxIndex(B)==1));
    BZTEST(all(minIndex(B)==0));

    MinMaxValue<typename T::T_numtype> mm=minmax(B);
    BZTEST(mm.min==-2);
    BZTEST(mm.max==5);
}

void grabner()
{
  BZ_USING_NAMESPACE(blitz::tensor)
  Array<float, 2> a1(2, 2), a2(2, 2), a3(2, 2);
  Array<float, 4> a4(2, 2, 2, 2);
  Array<float, 2> a5(2, 2), a6(2, 2);
  a1 = 1, 0, 0, 1;
  a2 = a1;
  a3 = a1;
  a4 = a1(i, k) * a2(l,k) * a3(j, l);
  a5 = sum(sum(a4, l), k);
  a6 = sum(sum(a1(i, k) * a2(l,k) * a3(j, l), l),k);
  BZTEST(all(a5==a6));
}

int main()
{
    Array<int,2> A(4,3);

    A = 0,  1, 2,
        3,  4, 5,
        6,  7, 8,
        9, 10, 11;

    BZTEST(sum(A) == 11*12/2);
    BZTEST(min(A) == 0);
    BZTEST(max(A) == 11);

    MinMaxValue<int> mm = minmax(A);
    //std::cerr << mm.min << ' ' << mm.max << std::endl;
    BZTEST(mm.min == 0);
    BZTEST(mm.max == 11);

    BZTEST(product(A) == 0);
    BZTEST(all(A >= 0));
    BZTEST(any(A == 7));
    BZTEST(count(A > 1 && A < 5) == 3);
    BZTEST(sum(pow2(A)) == 506);

    BZ_USING_NAMESPACE(blitz::tensor)

    BZTEST(sum(min(A,j)) == 0+3+6+9);
    BZTEST(sum(max(A(j,i),j)) == 9+10+11);

    Array<int,1> B(6);
    test_1d(B);
    TinyVector<int,6> Bv;
    test_1d(Bv);

    Array<int,3> C(2,2,2);
    C = 1, 0, 2, 3, 4, 7, 6, 5;

    BZTEST(sum(C) == 7*8/2);
    BZTEST(all(C <= 7));
    BZTEST(any(C == 5));
    BZTEST(!any(C == 8));
    BZTEST(sum(C(k,j,i)) == 7*8/2);
    BZTEST(sum(C(j,k,i)) == 7*8/2);

    // test behavior reported in bug 2058441
    BZTEST(sum(sum(sum(C,k),j)) == 7*8/2);
    Array<int,1> Cred(2);
    Cred=sum(sum(C,k),j);
    BZTEST(sum(Cred) == 7*8/2);
    grabner();

    mm = minmax(C);
    BZTEST(mm.min == 0);
    BZTEST(mm.max == 7);
    /* these don't work
    mm = minmax(C(i,j,k));
    BZTEST(mm.min == 0);
    BZTEST(mm.max == 7);
    mm = minmax(C(k,j,i));
    BZTEST(mm.min == 0);
    BZTEST(mm.max == 7);
    mm = minmax(C(j,k,i));
    BZTEST(mm.min == 0);
    BZTEST(mm.max == 7);
    */

    // test expression reductions

    cout << "Testing unary expression reductions\n";
    B= 2, -5, -4, -3, -1, -1;
    test_1dexpr(-B);
    cout << "Testing binary expression reductions\n";
    B= -3, 4, 3, 2, 0, 0;
    test_1dexpr(B+1);
    Array<float,1> BB(6);
    cout << "Testing function expression reductions\n";
    B=-2, 25, 16, 9, 1, 0;
    BB = 1, .5, .5, .5, 2, 0;
    test_1dexpr(pow(B,BB));
    cout << "Testing where expression reductions\n";
    B = -2, -30, 4, 3, 0, -1;
    BB = -2, 5, 3, -1, 1, 1;
    test_1dexpr(where(B>BB,B,BB));

    cout << "Testing functor expression reductions\n";
    B = -2, 5, 4, 3, 1, 1;
    crap c;
    test_1dexpr(c(B)-1);

    return 0;
}

