dnl Tries to determine if the copysign function exists in std namespace, 
dnl global namespace, or doesn't exist (in which case an appropriate
dnl macro is generated).  The result is put into ax_cxx_copysign.
dnl
AC_DEFUN([AX_CXX_COPYSIGN], [
AC_MSG_CHECKING([for copysign])
AC_LANG_PUSH(C++)
ax_cxx_copysign=
AC_COMPILE_IFELSE([AC_LANG_SOURCE([#include <cmath>
                   int main() {
                     std::copysign(1.0, 1.0);
                     return 0;
                   }])], 
                   [AC_MSG_RESULT(std::copysign)
                    ax_cxx_copysign="std::copysign(x, y)"])
if test "x$ax_cxx_copysign" = "x"; then
  AC_COMPILE_IFELSE([AC_LANG_SOURCE([#include <cmath>
                     int main() {
                       ::copysign(1.0, 1.0);
                       return 0;
                     }])], 
                     [AC_MSG_RESULT(::copysign)
                      ax_cxx_copysign="::copysign(x, y)"], 
                     [AC_MSG_RESULT(none)
                      ax_cxx_copysign="( ((y) != 0.0) ? ( ((y) > 0.0) ? (x) : -(x) ) : ( ((1.0 / y) > 0.0) ? (x) : -(x) ) )"])
fi
AC_LANG_POP(C++)
])
