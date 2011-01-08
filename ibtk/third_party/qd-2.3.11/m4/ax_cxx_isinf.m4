dnl Tries to determine if the isinf function exists in std namespace, 
dnl global namespace, or doesn't exist (in which case an appropriate
dnl macro is generated).  The result is put into ax_cxx_isinf.
dnl
AC_DEFUN([AX_CXX_ISINF], [
AC_MSG_CHECKING([for isinf])
AC_LANG_PUSH(C++)
ax_cxx_isinf=
AC_COMPILE_IFELSE([#include <cmath>
                   int main() {
                     std::isinf(1.0);
                     return 0;
                   }], 
                   [AC_MSG_RESULT(std::isinf)
                    ax_cxx_isinf="std::isinf(x)"])
if test "x$ax_cxx_isinf" = "x"; then
  AC_COMPILE_IFELSE([#include <cmath>
                     int main() {
                       ::isinf(1.0);
                       return 0;
                     }], 
                     [AC_MSG_RESULT(::isinf)
                      ax_cxx_isinf="::isinf(x)"], 
                     [AC_MSG_RESULT(none)
                      ax_cxx_isinf="( (x) != 0.0 && (x) == 2.0 * (x) )"])
fi
AC_LANG_POP(C++)
])
