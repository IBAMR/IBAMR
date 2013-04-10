


AC_DEFUN([AC_CXX_HAVE_RUSAGE],
[AC_CACHE_CHECK(whether the compiler has getrusage() function,
ac_cv_cxx_have_rusage,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <sys/resource.h>],[
struct rusage resUsage; getrusage(RUSAGE_SELF, &resUsage); return 0;],
 ac_cv_cxx_have_rusage=yes, ac_cv_cxx_have_rusage=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_have_rusage" = yes; then
  AC_DEFINE(HAVE_RUSAGE,,[define if the compiler has getrusage() function])
fi
])


