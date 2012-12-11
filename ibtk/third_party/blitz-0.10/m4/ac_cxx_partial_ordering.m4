dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_cxx_partial_ordering.html
dnl
AC_DEFUN([AC_CXX_PARTIAL_ORDERING],
[AC_CACHE_CHECK(whether the compiler supports partial ordering,
ac_cv_cxx_partial_ordering,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
template<int N> struct I {};
template<class T> struct A
{  int r;
   template<class T1, class T2> int operator() (T1, T2)       { r = 0; return r; }
   template<int N1, int N2>     int operator() (I<N1>, I<N2>) { r = 1; return r; }
};],[A<float> x, y; I<0> a; I<1> b; return x (a,b) + y (float(), double());],
 ac_cv_cxx_partial_ordering=yes, ac_cv_cxx_partial_ordering=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_partial_ordering" = yes; then
  AC_DEFINE(HAVE_PARTIAL_ORDERING,,
            [define if the compiler supports partial ordering])
fi
])
