dnl Tries to determine if the isnan function exists in std namespace, 
dnl global namespace, or doesn't exist (in which case an appropriate
dnl macro is generated).  The result is put into ax_cxx_isnan.
dnl
AC_DEFUN([AX_CXX_ISNAN], [
AC_MSG_CHECKING([for isnan])
AC_LANG_PUSH(C++)
ax_cxx_isnan=
AC_COMPILE_IFELSE([#include <cmath>
                   int main() {
                     std::isnan(1.0);
                     return 0;
                   }], 
                   [AC_MSG_RESULT(std::isnan)
                    ax_cxx_isnan="std::isnan(x)"])
if test "x$ax_cxx_isnan" = "x"; then
  AC_COMPILE_IFELSE([#include <cmath>
                     int main() {
                       ::isnan(1.0);
                       return 0;
                     }], 
                     [AC_MSG_RESULT(::isnan)
                      ax_cxx_isnan="::isnan(x)"], 
                     [AC_MSG_RESULT(none)
                      ax_cxx_isnan="((x) != (x))"])
fi
AC_LANG_POP(C++)
])
