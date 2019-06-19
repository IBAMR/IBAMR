# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CHECK_PRAGMA_KEYWORD],[
AC_MSG_CHECKING([whether compiler supports the _Pragma keyword])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[
  _Pragma("GCC diagnostic push")
  _Pragma("GCC diagnostic ignored \"-Wunknown-pragmas\"")
  _Pragma("GCC diagnostic ignored \"-Wpragmas\"")
  _Pragma("GCC diagnostic ignored \"-Wextra\"")
  _Pragma("GCC diagnostic pop")
]])],[
have_pragma_keyword=1
AC_MSG_RESULT(yes)],[
have_pragma_keyword=0
AC_MSG_RESULT(no)])
AC_DEFINE_UNQUOTED(HAVE_PRAGMA_KEYWORD,$have_pragma_keyword,[Boolean value indicating whether the C++ compiler supports the C99 _Pragma syntax for disabling warnings])
])
