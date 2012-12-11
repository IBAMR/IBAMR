dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_cxx_default_template_parameters.html
dnl
AC_DEFUN([AC_CXX_DEFAULT_TEMPLATE_PARAMETERS],
[AC_CACHE_CHECK(whether the compiler supports default template parameters,
ac_cv_cxx_default_template_parameters,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
template<class T = double, int N = 10> class A {public: int f() {return 0;}};
],[A<float> a; return a.f();],
 ac_cv_cxx_default_template_parameters=yes, ac_cv_cxx_default_template_parameters=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_default_template_parameters" = yes; then
  AC_DEFINE(HAVE_DEFAULT_TEMPLATE_PARAMETERS,,
            [define if the compiler supports default template parameters])
fi
])
