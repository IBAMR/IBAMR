dnl Tries to determine whether clock_gettime is useable.
dnl
AC_DEFUN([AX_CXX_CLOCK_GETTIME], [
AC_MSG_CHECKING([for clock_gettime useability])
AC_LANG_PUSH(C++)
AC_COMPILE_IFELSE([AC_LANG_SOURCE([
#include <time.h>
int main() {
  struct timespec tv;
  return clock_gettime(CLOCK_REALTIME, &tv);
}
])], [ax_cxx_clock_gettime="yes"], [ax_cxx_clock_gettime="no"])
AC_LANG_POP(C++)
AC_MSG_RESULT([$ax_cxx_clock_gettime])
])

