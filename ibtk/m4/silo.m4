AC_DEFUN([CONFIGURE_SILO],[

if test -d "${SILO_DIR}/lib" ; then
  LDFLAGS="-L${SILO_DIR}/lib $LDFLAGS"
fi
if test -d "${SILO_DIR}/include" ; then
  CPPFLAGS="-I${SILO_DIR}/include $CPPFLAGS"
fi

AC_CHECK_HEADER([silo.h],,AC_MSG_WARN([could not find header file silo.h]))
AC_LIB_HAVE_LINKFLAGS([silo])
if test "$HAVE_LIBSILO" == "yes" ; then
  LIBS="$LIBSILO $LIBS"
else
  AC_MSG_WARN([could not find working libsilo; checking for libsiloh5])
  AC_LIB_HAVE_LINKFLAGS([siloh5])
  if test "$HAVE_LIBSILOH5" == "yes" ; then
    LIBS="$LIBSILOH5 $LIBS"
  else
    AC_MSG_WARN([Silo support is enabled, but could not find working libsilo or libsiloh5])
  fi
fi

])