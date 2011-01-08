AC_DEFUN([CONFIGURE_HYPRE],[

if test -d "${HYPRE_DIR}/lib" ; then
  LDFLAGS="-L${HYPRE_DIR}/lib $LDFLAGS"
fi
if test -d "${HYPRE_DIR}/include" ; then
  CPPFLAGS="-I${HYPRE_DIR}/include $CPPFLAGS"
fi

AC_CHECK_HEADER([HYPRE.h],,AC_MSG_ERROR([could not find header file HYPRE.h]))

AC_LIB_HAVE_LINKFLAGS([HYPRE])
LIBS="$LIBHYPRE $LIBS"

])