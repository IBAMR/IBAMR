#include "testsuite.h"
#include <blitz/array.h>


BZ_USING_NAMESPACE(blitz)

int main()
{
   Array<int,2> A(5,5,FortranArray<2>());

   A = 1,2,3,4,5,
       6,7,8,9,10,
       11,12,13,14,15,
       16,17,18,19,20,
       21,22,23,24,25;

   // cout << A << endl;
   // A.dumpStructureInformation(cout);
   Array<int,2> Asub = A(Range(1,4,3), Range(2,4,2));
   // Asub.dumpStructureInformation(cout);
   // cout << A(Range(1,4,3), Range(2,4,2)) << endl;
   BZTEST(Asub(1,1) == 6);
   BZTEST(Asub(1,2) == 16);
   BZTEST(Asub(2,1) == 9);
   BZTEST(Asub(2,2) == 19);
   return 0;
}

