dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_cxx_template_qualified_return_type.html
dnl
AC_DEFUN([AC_CXX_TEMPLATE_QUALIFIED_RETURN_TYPE],
[AC_CACHE_CHECK(whether the compiler supports template-qualified return types,
ac_cv_cxx_template_qualified_return_type,
[AC_REQUIRE([AC_CXX_TYPENAME])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
#ifndef HAVE_TYPENAME
 #define typename
#endif
template<class X, class Y> struct promote_trait             { typedef X T; };
template<>                 struct promote_trait<int, float> { typedef float T; };
template<class T> class A { public : A () {} };
template<class X, class Y>
A<typename promote_trait<X,Y>::T> operator+ (const A<X>&, const A<Y>&)
{ return A<typename promote_trait<X,Y>::T>(); }
],[A<int> x; A<float> y; A<float> z = x + y; return 0;],
 ac_cv_cxx_template_qualified_return_type=yes, ac_cv_cxx_template_qualified_return_type=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_template_qualified_return_type" = yes; then
  AC_DEFINE(HAVE_TEMPLATE_QUALIFIED_RETURN_TYPE,,
            [define if the compiler supports template-qualified return types])
fi
])
