#include "testsuite.h"
#include <blitz/array.h>

// test the update operators

BZ_USING_NAMESPACE(blitz)

int main()
{
  // First with scalars
  {
    Array<int,1> a(5);
    a=1;

    a+=1;
    BZTEST(all(a==2));
    
    a-=1;
    BZTEST(all(a==1));
    
    a*=8;
    BZTEST(all(a==8));
    
    a/=2;
    BZTEST(all(a==4));
    
    a%=3;
    BZTEST(all(a==1));
    
    a^=2;
    BZTEST(all(a==3));
    
    a&=6;
    BZTEST(all(a==2));
    
    a|=8;
    BZTEST(all(a==10));
    
    a<<=1;
    BZTEST(all(a==20));
    
    a>>=2;
    BZTEST(all(a==5));
  }

  // And now with Arrays
  {
    Array<int,1> a(5),b(5);
    a=1;

    b=1;
    a+=b;
    BZTEST(all(a==2));
    
    a-=b;
    BZTEST(all(a==1));
    
    b=8;
    a*=b;
    BZTEST(all(a==8));
    
    b=2;
    a/=b;
    BZTEST(all(a==4));
    
    b=3;
    a%=b;
    BZTEST(all(a==1));
    
    b=2;
    a^=b;
    BZTEST(all(a==3));
    
    b=6;
    a&=b;
    BZTEST(all(a==2));
    
    b=8;
    a|=b;
    BZTEST(all(a==10));
    
    b=1;
    a<<=b;
    BZTEST(all(a==20));
    
    b=2;
    a>>=b;
    BZTEST(all(a==5));
  }

  // And now with expressions
  {
    Array<int,1> a(5),b(5);
    a=1;

    b=1;
    a+=1*b;
    BZTEST(all(a==2));
    
    a-=1*b;
    BZTEST(all(a==1));
    
    b=8;
    a*=1*b;
    BZTEST(all(a==8));
    
    b=2;
    a/=1*b;
    BZTEST(all(a==4));
    
    b=3;
    a%=1*b;
    BZTEST(all(a==1));
    
    b=2;
    a^=1*b;
    BZTEST(all(a==3));
    
    b=6;
    a&=1*b;
    BZTEST(all(a==2));
    
    b=8;
    a|=1*b;
    BZTEST(all(a==10));
    
    b=1;
    a<<=1*b;
    BZTEST(all(a==20));
    
    b=2;
    a>>=1*b;
    BZTEST(all(a==5));
  }
  



  
  return 0;
}

