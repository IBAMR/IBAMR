AC_DEFUN([AC_CXX_MATH_ABSINT_IN_NAMESPACE_STD],
[AC_CACHE_CHECK(whether the compiler has C math abs(integer type) in namespace std,
ac_cv_cxx_mathabsint_std,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_REQUIRE([AC_CXX_MATH_FN_IN_NAMESPACE_STD])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <cstdlib>],[
 int i = std::abs(1); long j = std::labs(1L); long k = std::abs(1L); return 0;],
 ac_cv_cxx_mathabsint_std=yes, ac_cv_cxx_mathabsint_std=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_mathabsint_std" = yes; then
  AC_DEFINE(MATH_ABSINT_IN_NAMESPACE_STD,,[define if the compiler has C math abs(integer types) in namespace std])
fi
])
