

AC_DEFUN([AC_CXX_HAVE_CLIMITS],
[AC_CACHE_CHECK(whether the compiler has <climits> header,
ac_cv_cxx_have_climits,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <climits>],[int i = INT_MIN; return 0;],
 ac_cv_cxx_have_climits=yes, ac_cv_cxx_have_climits=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_have_climits" = yes; then
  AC_DEFINE(HAVE_CLIMITS,,[define if the compiler has <climits> header])
fi
])


