# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_GSL],[
echo
echo "==================================================="
echo "Configuring optional package GNU Scientific Library"
echo "==================================================="

PACKAGE_SETUP_ENVIRONMENT

AC_ARG_ENABLE([gsl],
  AS_HELP_STRING(--enable-gsl,enable support for the optional GSL library @<:@default=yes@:>@),
                 [case "$enableval" in
                    yes)  GSL_ENABLED=yes ;;
                    no)   GSL_ENABLED=no ;;
                    *)    AC_MSG_ERROR(--enable-gsl=$enableval is invalid; choices are "yes" and "no") ;;
                  esac],[GSL_ENABLED=no])

AM_CONDITIONAL([GSL_ENABLED],[test "$GSL_ENABLED" = yes])

AC_ARG_WITH([gsl],
  AS_HELP_STRING(--with-gsl=PATH,location of optional GSL installation),
  [if test "$GSL_ENABLED" = no ; then
     AC_MSG_WARN([--with-gsl is specified, but support for gsl is disabled])
   else
     if test ! -d "$withval" ; then
       AC_MSG_ERROR([it is necessary to specify an existing directory when using --with-gsl=PATH])
     fi
     GSL_DIR=$withval
   fi])

if test "$GSL_ENABLED" = yes; then
  if test x$GSL_DIR != x ; then
    if test -d "${GSL_DIR}/include" ; then
      GSL_CPPFLAGS="-I${GSL_DIR}/include"
    fi
    if test -d "${GSL_DIR}/lib" ; then
      GSL_LDFLAGS="-L${GSL_DIR}/lib"
    fi
  fi

  CPPFLAGS_PREPEND($GSL_CPPFLAGS)
  AC_CHECK_HEADER([gsl/gsl_linalg.h],,AC_MSG_ERROR([could not find header file gsl_linalg.h]))
  AC_CHECK_HEADER([gsl/gsl_integration.h],,AC_MSG_ERROR([could not find header file gsl_integration.h]))

  LDFLAGS_PREPEND($GSL_LDFLAGS)
  AC_LIB_HAVE_LINKFLAGS([gslcblas])
    if test "$HAVE_LIBGSLCBLAS" = no ; then
     AC_MSG_ERROR([could not find working libgslcblas])
    fi
  AC_LIB_HAVE_LINKFLAGS([gsl],[m,gslcblas])
  if test "$HAVE_LIBGSL" = no ; then
    AC_MSG_ERROR([could not find working libgsl])
  fi

  PACKAGE_CPPFLAGS_PREPEND($GSL_CPPFLAGS)
  PACKAGE_LDFLAGS_PREPEND($GSL_LDFLAGS)
  PACKAGE_LIBS_PREPEND("$LIBGSLCBLAS")
  PACKAGE_LIBS_PREPEND("$LIBGSL")
 
else
  AC_MSG_NOTICE([Optional package GNU Scientific Library is DISABLED])
fi

PACKAGE_RESTORE_ENVIRONMENT

])
