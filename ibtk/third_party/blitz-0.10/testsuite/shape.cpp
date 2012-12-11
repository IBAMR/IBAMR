/*
 * Intended coverage: blitz/arrayshape.h
 */

#include "testsuite.h"

#include <blitz/array.h>
#include <blitz/tinyvec2.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    TinyVector<int,1> a1(1);
    BZTEST(all(a1 == shape(1)));

    TinyVector<int,2> a2;
    a2 = 5, 10;
    BZTEST(all(a2 == shape(5,10)));

    TinyVector<int,3> a3;
    a3 = 1, 2, 3;
    BZTEST(sum(a3 == shape(1,2,3)) == 3);

    TinyVector<int,4> a4;
    a4 = 1, 2, 3, 4;
    BZTEST(sum(a4 == shape(1,2,3,4)) == 4);

    TinyVector<int,5> a5;
    a5 = 1, 2, 3, 4, 5;
    BZTEST(count(a5 == shape(1,2,3,4,5)) == 5);

    TinyVector<int,6> a6;
    a6 = 1, 2, 3, 4, 5, 6;
    BZTEST(count(a6 == shape(1,2,3,4,5,6)) == 6);

    TinyVector<int,7> a7;
    a7 = 1, 2, 3 ,4, 5, 6, 7;
    BZTEST(count(a7 == shape(1,2,3,4,5,6,7)) == 7);

    TinyVector<int,8> a8;
    a8 = 1, 2, 3, 4, 5, 6, 7, 8;
    BZTEST(count(a8 == shape(1,2,3,4,5,6,7,8)) == 8);

    TinyVector<int,9> a9;
    a9 = 1, 2, 3, 4, 5, 6, 7, 8, 9;
    BZTEST(count(a9 == shape(1,2,3,4,5,6,7,8,9)) == 9);

    TinyVector<int,10> a10;
    a10 = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
    BZTEST(count(a10 == shape(1,2,3,4,5,6,7,8,9,10)) == 10);

    TinyVector<int,11> a11;
    a11 = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11;
    BZTEST(count(a11 == shape(1,2,3,4,5,6,7,8,9,10,11)) == 11);
}

