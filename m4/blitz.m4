AC_DEFUN([CONFIGURE_BLITZ],[

if test -d "${BLITZ_DIR}/lib" ; then
  LDFLAGS="-L${BLITZ_DIR}/lib $LDFLAGS"
fi
if test -d "${BLITZ_DIR}/include" ; then
 CPPFLAGS="-I${BLITZ_DIR}/include $CPPFLAGS"
fi

AC_CHECK_HEADER([blitz/blitz.h],,AC_MSG_WARN([could not find header file blitz/blitz.h]))
AC_LIB_HAVE_LINKFLAGS([blitz])
LIBS="$LIBBLITZ $LIBS"

])
