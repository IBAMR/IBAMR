dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_cxx_enum_computations.html
dnl
AC_DEFUN([AC_CXX_ENUM_COMPUTATIONS],
[AC_CACHE_CHECK(whether the compiler handle computations inside an enum,
ac_cv_cxx_enum_computations,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
struct A { enum { a = 5, b = 7, c = 2 }; };
struct B { enum { a = 1, b = 6, c = 9 }; };
template<class T1, class T2> struct Z
{ enum { a = (T1::a > T2::a) ? T1::a : T2::b,
         b = T1::b + T2::b,
         c = (T1::c * T2::c + T2::a + T1::a)
       };
};],[
return (((int)Z<A,B>::a == 5)
     && ((int)Z<A,B>::b == 13)
     && ((int)Z<A,B>::c == 24)) ? 0 : 1;],
 ac_cv_cxx_enum_computations=yes, ac_cv_cxx_enum_computations=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_enum_computations" = yes; then
  AC_DEFINE(HAVE_ENUM_COMPUTATIONS,,
            [define if the compiler handle computations inside an enum])
fi
])
