dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_cxx_full_specialization_syntax.html
dnl
AC_DEFUN([AC_CXX_FULL_SPECIALIZATION_SYNTAX],
[AC_CACHE_CHECK(whether the compiler recognizes the full specialization syntax,
ac_cv_cxx_full_specialization_syntax,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
template<class T> class A        { public : int f () const { return 1; } };
template<>        class A<float> { public:  int f () const { return 0; } };],[
A<float> a; return a.f();],
 ac_cv_cxx_full_specialization_syntax=yes, ac_cv_cxx_full_specialization_syntax=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_full_specialization_syntax" = yes; then
  AC_DEFINE(HAVE_FULL_SPECIALIZATION_SYNTAX,,
            [define if the compiler recognizes the full specialization syntax])
fi
])
