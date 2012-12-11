dnl Tries to determine appropriate function names for 
dnl fused multiply-add (fma) and fused multiply-subtract (fms).
dnl
dnl Usage: ACX_CXX_FMA(list)
dnl   where fma_list can contain any of ibm, c99.  If list is
dnl   blank it will search for a first compatible function.
dnl
AC_DEFUN([AX_CXX_FMA], [
AC_MSG_CHECKING([for fused multiply-add/subtract])
AC_LANG_PUSH(C++)
ax_cxx_fma_list=$1
if test "x$ax_cxx_fma_list" = "x"; then
  ax_cxx_fma_list="ibm gnu c99 compiler"
fi
ax_cxx_fma=
ax_cxx_fms=
for name in $ax_cxx_fma_list; do
  if test "x$ax_cxx_fma" = "x"; then
    case $name in
      ibm)
        # IBM VisualAge C++ __fmadd / __fmsub.
        AC_RUN_IFELSE([AC_LANG_SOURCE([#include <cmath>
                       #include <builtins.h>
                       int main() {
                         double d = std::ldexp(1.0, -52);
                         double x = __fmadd(1.0 + d, 1.0 - d, -1.0);
                         double y = __fmsub(1.0 + d, 1.0 - d, 1.0);
                         return (x == -d*d && y == -d*d) ? 0 : 1;
                       }])], 
                      [ax_cxx_fma="__fmadd(x,y,z)"
                       ax_cxx_fms="__fmsub(x,y,z)"
                       AC_DEFINE([QD_VACPP_BUILTINS_H], [1], 
                                 [Set to 1 if using VisualAge C++ compiler for __fmadd builtin.])])
      ;;
      gnu)
        # Later gcc (3.4 and later) have __builtin_fma that seems to work.
        AC_RUN_IFELSE([AC_LANG_SOURCE([#include <cmath>
                       int main() {
                         double d = std::ldexp(1.0, -52);
                         return (__builtin_fma(1.0 + d, 1.0 - d, -1.0) == -d*d ? 0 : 1);
                       }])], 
                      [ax_cxx_fma="__builtin_fma(x,y,z)"
                       ax_cxx_fms="__builtin_fma(x,y,-z)"])
      ;;
      ia64)
        # Intel and HP compilers for IA 64 architecture seems to have 
        # _Asm_fma/fms macros.  Not much documentation is available for 
        # these...
        AC_RUN_IFELSE([AC_LANG_SOURCE([#include <cmath>
                       int main() {
                         double d = std::ldexp(1.0, -52);
                         return (_Asm_fma(2, 1.0 + d, 1.0 - d, -1.0) == -d*d ? 0 : 1);
                       }])], 
                      [ax_cxx_fma="_Asm_fma(2, x,y,z)"
                       ax_cxx_fms="_Asm_fms(2, x,y,z)"])
      ;;
      c99)
        # Try C99 fma() function.  Some platforms doesn't seem to implement this
        # correctly (Apple gcc-3.3 for example).
        AC_RUN_IFELSE([AC_LANG_SOURCE([#include <cmath>
                       int main() {
                         double d = std::ldexp(1.0, -52);
                         return (fma(1.0 + d, 1.0 - d, -1.0) == -d*d ? 0 : 1);
                       }])], 
                      [ax_cxx_fma="fma(x,y,z)"
                       ax_cxx_fms="fma(x,y,-z)"])
      ;;
      compiler)
        # Try relying on the compiler to optimize x * y + z into an fma.
        # This method is not recommended since if it is inlined it does not
        # always produce the same correct code.
        AC_RUN_IFELSE([AC_LANG_SOURCE([#include <cmath>
                       int main() {
                         double d = std::ldexp(1.0, -52);
                         return ( (1.0 + d) * (1.0 - d) - 1.0 == -d*d ? 0 : 1);
                       }])],
                       [ax_cxx_fma="((x)*(y) + (z))"
                        ax_cxx_fms="((x)*(y) - (z))"])
      ;;
      *) AC_MSG_ERROR([Unknown option $name to --enable-fma.]) ;;
    esac
  fi
done
AC_LANG_POP(C++)
if test "x$ax_cxx_fma" != "x"; then
  AC_MSG_RESULT([$ax_cxx_fma, $ax_cxx_fms])
else
  AC_MSG_RESULT(none)
fi
])
