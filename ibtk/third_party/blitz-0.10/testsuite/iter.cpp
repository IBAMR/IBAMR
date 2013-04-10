#include "testsuite.h"
#include <blitz/array.h>

#ifdef BZ_HAVE_STL
#include <iterator>
#include <algorithm>
BZ_USING_NAMESPACE(std)
#endif

BZ_USING_NAMESPACE(blitz)

void check(const Array<int,2>& A, const Array<int,1>& b)
{
    int i = 0;
    Array<int,2>::const_iterator beg = A.begin(), end = A.end();
    Array<int,2>::const_iterator iter = beg;

    while (iter != end) {
        BZTEST((*iter++) == b(i++));
    }
    BZTEST(i == A.numElements());

    while (iter!=beg) {
        BZTEST((*--iter) == b(--i));
    }
    BZTEST(i == 0);
}

#ifdef BZ_HAVE_STL
template <typename _Iter>
void checkInterface(_Iter iter)
{
    typedef typename iterator_traits<_Iter>::value_type value_type;
    typedef typename iterator_traits<_Iter>::reference reference;
    typedef typename iterator_traits<_Iter>::pointer pointer;
    value_type x = *iter;
    BZTEST(*iter == x);
    reference y(x);
    pointer p = &(*iter);
    BZTEST(*p == y);
}
#endif // BZ_HAVE_STL

int main()
{
  {
    Array<int,2> A(2,3);
    A = 0, 1, 2,
        3, 4, 5;
    Array<int,1> b(6);
    b = 0, 1, 2, 3, 4, 5;

    check(A, b);
  }

  {
    Array<int,2> A(2,3,FortranArray<2>());
    A = tensor::i * 3 + tensor::j;
    Array<int,1> b(6);
    b = 4, 7, 5, 8, 6, 9;
    check(A,b);
  }

#if 0
  {
    Array<int,2> B(5,5,FortranArray<2>());
    B = tensor::i + 5 * tensor::j;

    B.reverseSelf(1);
    check(B);
  }
#endif

  {
    Array<int,2> B(6,6,FortranArray<2>());
    B = tensor::i * 6 + tensor::j;
    Array<int,1> b(12);
    b = 7, 13, 19, 8, 14, 20, 9, 15, 21, 10, 16, 22;
    check(B(Range(1,3),Range(1,4)), b);
  }

  {
    Array<int,2> B;
    Array<int,2>::iterator iter = B.begin(), end = B.end();
    BZTEST(iter == end);
  }

#ifdef BZ_HAVE_STL
  {
    Array<int,2> A(3,3);
    A = 1, 2, 3,
        3, 2, 1,
        4, 3, 2;
    replace(A.begin(), A.end(), 3, 0); // replace each 3 with 0
    Array<int,1> b(9);
    b = 1, 2, 0, 0, 2, 1, 4, 0, 2;
    check(A,b);
    Array<int,2>::iterator iter;
    iter = adjacent_find(A.begin(),A.end());
    TinyVector<int,2> pos = iter.position(), ans(0,2);
    BZTEST(pos[0] == ans[0] && pos[1] == ans[1]);
    checkInterface(iter);
    const Array<int,2>& B(A);
    Array<int,2>::const_iterator citer = B.begin();
    checkInterface(citer);
  }
#endif // BZ_HAVE_STL

}

