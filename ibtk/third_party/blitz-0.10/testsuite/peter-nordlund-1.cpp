#include "testsuite.h"

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main() {
  typedef float dtype;
  const int sz = 7;
  Array<dtype,2> A(sz,sz);
  A = 0;

  // cout << "A,0 = " << A << endl;
  BZTEST(sum(A) == 0);
  
  Array<dtype,2> B = A(Range(0,sz-1,3), Range(0,sz-3,2));
  // cout << "B=A, B= " << B << endl;

  B = 1;
  BZTEST(sum(B) == 9);
  BZTEST(sum(A) == 9);

  //cout << "B=1, B= " << B << endl;
  Array<dtype,2> C(3,3);
  C = B; // Seems wrong to 

  BZTEST(sum(C) == 9);

  // cout << "C=B, C=" << C << endl;

  Array<dtype,2> D(sz,sz);
  D = A;

  Array<dtype,2> t1(sz,sz);
  t1 = 1, 0, 1, 0, 1, 0, 0,
       0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0,
       1, 0, 1, 0, 1, 0, 0,
       0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 
       1, 0, 1, 0, 1, 0, 0;

  BZTEST(count(D==t1) == sz*sz);

  // cout << "D=A, D= " << D << endl;
  D(0,0) = 2;
  D(2,2) = 4;
  D(3,5) = 6;

  t1 = 2, 0, 1, 0, 1, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 
       0, 0, 4, 0, 0, 0, 0,
       1, 0, 1, 0, 1, 6, 0,
       0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 
       1, 0, 1, 0, 1, 0, 0;

  BZTEST(count(D == t1) == sz*sz);

  // cout << "D(x,y)=z, D= " << D << endl;

  Array<dtype,2> G(5,5);
  Range I(1,5), J(1,5);
  G = D(I,J).copy(); // Seems wrong to me

  Array<dtype,2> t2(5,5);
  t2 = 0, 0, 0, 0, 0,
       0, 4, 0, 0, 0,
       0, 1, 0, 1, 6,
       0, 0, 0, 0, 0,
       0, 0, 0, 0, 0;

  BZTEST(count(G == t2) == 25);

  // cout << "G=D(I,J).copy(), G= " << G << endl;
  // cout << "----------------------"<<endl;
  // cout << "A = " << A << endl;// Ok
  // cout << "B = " << B << endl;// Ok
  // cout << "C = " << C << endl;// Seems wrong
  // cout << "D = " << D << endl;// Ok
  // cout << "G = " << G << endl;// Seems wrong
  // cout << "D(I,J) = " << D(I,J) << endl;// Ok

  return 0;
}

