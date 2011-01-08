#include <cstdlib>
#include <qd/fpu.h>
#include "config.h"

#define f_main FC_FUNC_(f_main, F_MAIN)

extern "C" void f_main();

int main() {
  fpu_fix_start(NULL);
  f_main();
  return 0;
}

