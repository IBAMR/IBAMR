#include "testsuite.h"

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

Array< double, 1> a(100); 

int main(){

  BZTEST(a.size() == 100);
  return 0;
}

