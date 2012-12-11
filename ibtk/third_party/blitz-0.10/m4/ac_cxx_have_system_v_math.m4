dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_cxx_have_system_v_math.html
dnl
AC_DEFUN([AC_CXX_HAVE_SYSTEM_V_MATH],
[AC_CACHE_CHECK(whether the compiler supports System V math library,
ac_cv_cxx_have_system_v_math,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 ac_save_LIBS="$LIBS"
 LIBS="$LIBS -lm"
 AC_TRY_LINK([
#ifndef _ALL_SOURCE
 #define _ALL_SOURCE
#endif
#ifndef _XOPEN_SOURCE
 #define _XOPEN_SOURCE
#endif
#ifndef _XOPEN_SOURCE_EXTENDED
 #define _XOPEN_SOURCE_EXTENDED 1
#endif
#include <math.h>],[double x = 1.0; double y = 1.0;
_class(x); trunc(x); finite(x); itrunc(x); nearest(x); rsqrt(x); uitrunc(x);
copysign(x,y); drem(x,y); unordered(x,y);
return 0;],
 ac_cv_cxx_have_system_v_math=yes, ac_cv_cxx_have_system_v_math=no)
 LIBS="$ac_save_LIBS"
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_have_system_v_math" = yes; then
  AC_DEFINE(HAVE_SYSTEM_V_MATH,,[define if the compiler supports System V math library])
fi
])
