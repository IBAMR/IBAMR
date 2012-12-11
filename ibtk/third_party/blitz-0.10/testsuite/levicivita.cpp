#include "testsuite.h"
#include<blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
  Array<double,2> a(3);

  TinyVector<double,3> v,w,z;
  Array<double,1> av(3);

  v=.5,0,0;
  w=0,2,0;
  av=v;

  // manually code cross product using LeviCivita symbol
  z=sum(sum(LeviCivita()*v(tensor::j)*w(tensor::k), tensor::k), tensor::j);

  //cout << z<<endl;
  BZTEST(z[0]==0);
  BZTEST(z[1]==0);
  BZTEST(z[2]==1);

  // this uses dedicated tvcross function
  z=cross(v,w);

  //cout << z << endl;
  BZTEST(z[0]==0);
  BZTEST(z[1]==0);
  BZTEST(z[2]==1);

  // this uses the levicivita expression
  z=cross(av,w);

  //cout << z << endl;
  BZTEST(z[0]==0);
  BZTEST(z[1]==0);
  BZTEST(z[2]==1);

  // see if type promotion works
  TinyVector<char,3> cv;
  TinyVector<int,3> cw;
  cv=0,1,2;
  
  Array<int,1> cav(3);
  cav=1,2,1;

  cw=cross(cv,cav);

  //cout << cw << endl;

  BZTEST(cw[0]==-3);
  BZTEST(cw[1]==2);
  BZTEST(cw[2]==-1);

  return 0;
}
