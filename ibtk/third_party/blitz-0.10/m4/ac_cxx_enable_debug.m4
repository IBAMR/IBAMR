
AC_DEFUN([AC_CXX_ENABLE_DEBUG],[

AC_MSG_CHECKING([whether to enable $CXX debug flags])
AC_ARG_ENABLE(debug,
AS_HELP_STRING([--enable-debug],[Enable compiler debugging flags]),
[if test "$enableval" = yes; then
  AC_MSG_RESULT([yes])
  CXXFLAGS="$CXXFLAGS $CXX_DEBUG_FLAGS"
fi],[AC_MSG_RESULT([no])])
AC_ARG_VAR([CXX_DEBUG_FLAGS],[C++ compiler debugging flags])
])
