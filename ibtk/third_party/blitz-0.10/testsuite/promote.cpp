#include "testsuite.h"

#include <blitz/promote.h>

BZ_USING_NAMESPACE(blitz)

#ifndef BZ_HAVE_PARTIAL_SPECIALIZATION
int main() { return 0; }
#else

template<typename T1, typename T2>
struct sameType {
    enum { same = 0 };
};

template<typename T>
struct sameType<T,T> {
    enum { same = 1 };
};

template<typename T1, typename T2, typename T_promote>
inline void check(T1,T2,T_promote)
{
#if 0
    cout << BZ_DEBUG_TEMPLATE_AS_STRING_LITERAL(T1)
         << " + "
         << BZ_DEBUG_TEMPLATE_AS_STRING_LITERAL(T2)
         << " -> "
         << BZ_DEBUG_TEMPLATE_AS_STRING_LITERAL(BZ_PROMOTE(T1,T2))
         << " (" << BZ_DEBUG_TEMPLATE_AS_STRING_LITERAL(T_promote)
         << ")" << endl;
#endif

    int same = (int)sameType<BZ_PROMOTE(T1,T2),T_promote>::same;
    BZTEST(same);
}

    struct foo {
        int x[20];
    };
    struct bigfoo {
        int x[40];
    };
    struct tinyfoo {
        unsigned char z;
    };

int main()
{
    typedef short int shortint;
    typedef unsigned char uchar;
    typedef unsigned int uint;
    typedef unsigned short int shortuint;

    check(char(), float(), float());

#ifndef BZ_DISABLE_NEW_PROMOTE  
// Old type-promotion doesn't handle bool
    check(bool(), bool(), int());
    check(bool(), shortint(), int());
#endif

    check(float(), double(), double());
    check(double(), float(), double());
    check(char(), char(), int());
    check(uchar(), char(), int());
    check(uint(), int(), uint());
    check(shortuint(), int(), uint());
    check(shortuint(), shortint(), uint());

#ifdef BZ_HAVE_COMPLEX
    check(complex<float>(),int(),complex<float>());
#ifndef BZ_DISABLE_NEW_PROMOTE  
// Old type-promotion doesn't handle bool
    check(bool(),complex<long double>(),complex<long double>());
// Old type-promotion doesn't handle user-defined types
    check(complex<float>(), foo(), foo());
#endif
    check(float(), complex<float>(), complex<float>());
#endif

#ifndef BZ_DISABLE_NEW_PROMOTE  
// Old type-promotion doesn't handle user-defined types
    check(int(), foo(), foo());
    check(bigfoo(), foo(), bigfoo());
    check(foo(), bigfoo(), bigfoo());
    check(int(), tinyfoo(), tinyfoo());
    check(tinyfoo(), int(), tinyfoo());
#endif

    return 0;
}
#endif
