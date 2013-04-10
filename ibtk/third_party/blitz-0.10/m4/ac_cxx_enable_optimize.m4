
AC_DEFUN([AC_CXX_ENABLE_OPTIMIZE],[

AC_MSG_CHECKING([whether to enable $CXX optimization flags])
AC_ARG_ENABLE(optimize,
AS_HELP_STRING([--enable-optimize],[Enable compiler optimization flags]),
[if test "$enableval" = yes; then
	AC_MSG_RESULT([yes])
  CXXFLAGS="$CXXFLAGS $CXX_OPTIMIZE_FLAGS"
fi],[AC_MSG_RESULT([no])])
AC_ARG_VAR([CXX_OPTIMIZE_FLAGS],[C++ compiler optimization flags])
])
