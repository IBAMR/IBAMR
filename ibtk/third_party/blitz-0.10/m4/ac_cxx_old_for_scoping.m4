dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_cxx_old_for_scoping.html
dnl
AC_DEFUN([AC_CXX_OLD_FOR_SCOPING],
[AC_CACHE_CHECK(whether the compiler accepts the old for scoping rules,
ac_cv_cxx_old_for_scoping,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE(,[int z;for (int i=0; i < 10; ++i)z=z+i;z=i;return z;],
 ac_cv_cxx_old_for_scoping=yes, ac_cv_cxx_old_for_scoping=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_old_for_scoping" = yes; then
  AC_DEFINE(HAVE_OLD_FOR_SCOPING,,[define if the compiler accepts the old for scoping rules])
fi
])
