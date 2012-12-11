#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
  Array<complex<float>,2> test(4,4) ;
 Array<float,2> test2(2,2) ;

 test2 = real(test(Range(0,1),Range(0,1))) ;
 Array<float,2> test3;
 test3.reference(test2(Range(0,1),Range(0,1)));

 return 0;
}

