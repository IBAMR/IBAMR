dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_cxx_enum_computations_with_cast.html
dnl
AC_DEFUN([AC_CXX_ENUM_COMPUTATIONS_WITH_CAST],
[AC_CACHE_CHECK(whether the compiler handles (int) casts in enum computations,
ac_cv_cxx_enum_computations_with_cast,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
struct A { enum { a = 5, b = 7, c = 2 }; };
struct B { enum { a = 1, b = 6, c = 9 }; };
template<class T1, class T2> struct Z
{ enum { a = ((int)T1::a > (int)T2::a) ? (int)T1::a : (int)T2::b,
         b = (int)T1::b + (int)T2::b,
         c = ((int)T1::c * (int)T2::c + (int)T2::a + (int)T1::a)
       };
};],[
return (((int)Z<A,B>::a == 5)
     && ((int)Z<A,B>::b == 13)
     && ((int)Z<A,B>::c == 24)) ? 0 : 1;],
 ac_cv_cxx_enum_computations_with_cast=yes, ac_cv_cxx_enum_computations_with_cast=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_enum_computations_with_cast" = yes; then
  AC_DEFINE(HAVE_ENUM_COMPUTATIONS_WITH_CAST,,
            [define if the compiler handles (int) casts in enum computations])
fi
])
