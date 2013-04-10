dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_cxx_explicit_template_function_qualification.html
dnl
AC_DEFUN([AC_CXX_EXPLICIT_TEMPLATE_FUNCTION_QUALIFICATION],
[AC_CACHE_CHECK(whether the compiler supports explicit template function qualification,
ac_cv_cxx_explicit_template_function_qualification,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
template<class Z> class A { public : A() {} };
template<class X, class Y> A<X> to (const A<Y>&) { return A<X>(); }
],[A<float> x; A<double> y = to<double>(x); return 0;],
 ac_cv_cxx_explicit_template_function_qualification=yes, ac_cv_cxx_explicit_template_function_qualification=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_explicit_template_function_qualification" = yes; then
  AC_DEFINE(HAVE_EXPLICIT_TEMPLATE_FUNCTION_QUALIFICATION,,
            [define if the compiler supports explicit template function qualification])
fi
])
