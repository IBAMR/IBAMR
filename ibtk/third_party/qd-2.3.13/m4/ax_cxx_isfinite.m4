dnl Tries to determine if the isfinite function exists in std namespace, 
dnl global namespace, or doesn't exist (in which case an appropriate
dnl macro is generated).  The result is put into ax_cxx_isfinite.
dnl
AC_DEFUN([AX_CXX_ISFINITE], [
AC_MSG_CHECKING([for isfinite])
AC_LANG_PUSH(C++)
ax_cxx_isfinite=
AC_COMPILE_IFELSE([AC_LANG_SOURCE([#include <cmath>
                   int main() {
                     std::isfinite(1.0);
                     return 0;
                   }])], 
                   [AC_MSG_RESULT(std::isfinite)
                    ax_cxx_isfinite="std::isfinite(x)"])
if test "x$ax_cxx_isfinite" = "x"; then
  AC_COMPILE_IFELSE([AC_LANG_SOURCE([#include <cmath>
                     int main() {
                       ::isfinite(1.0);
                       return 0;
                     }])], 
                     [AC_MSG_RESULT(::isfinite)
                      ax_cxx_isfinite="::isfinite(x)"], 
                     [AC_MSG_RESULT(none)
                      ax_cxx_isfinite="( ((x) == 0.0) || ((x) != (2.0 * (x))) )"])
fi
AC_LANG_POP(C++)
])
