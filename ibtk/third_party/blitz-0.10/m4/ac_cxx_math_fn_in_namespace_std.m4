AC_DEFUN([AC_CXX_MATH_FN_IN_NAMESPACE_STD],
[AC_CACHE_CHECK(whether the compiler has C math functions in namespace std,
ac_cv_cxx_mathfn_std,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <cmath>
namespace blitz { double pow(double x, double y){ return std::pow(x,y); } };],[
using namespace blitz; double x = 1.0, y = 1.0; blitz::pow(x,y); return 0;],
 ac_cv_cxx_mathfn_std=yes, ac_cv_cxx_mathfn_std=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_mathfn_std" = yes; then
  AC_DEFINE(MATH_FN_IN_NAMESPACE_STD,,[define if the compiler has C math functions in namespace std])
fi
])



