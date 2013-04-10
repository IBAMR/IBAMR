AC_DEFUN([AC_CXX_ISNAN_IN_NAMESPACE_STD],
[AC_CACHE_CHECK(whether the compiler has isnan function in namespace std,
ac_cv_cxx_isnan_std,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <cmath>
namespace blitz { int isnan(float x){ return std::isnan(x); } };],[
using namespace blitz; float x = 1.0; blitz::isnan(x); return 0;],
 ac_cv_cxx_isnan_std=yes, ac_cv_cxx_isnan_std=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_isnan_std" = yes; then
  AC_DEFINE(ISNAN_IN_NAMESPACE_STD,,[define if the compiler has isnan function in namespace std])
fi
])
