# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_SILO],[
echo
echo "================================="
echo "Configuring optional package Silo"
echo "================================="

PACKAGE_SETUP_ENVIRONMENT

AC_ARG_ENABLE([silo],
  AS_HELP_STRING(--enable-silo,enable support for the optional Silo library @<:@default=yes@:>@),
                 [case "$enableval" in
                    yes)  SILO_ENABLED=yes ;;
                    no)   SILO_ENABLED=no ;;
                    *)    AC_MSG_ERROR(--enable-silo=$enableval is invalid; choices are "yes" and "no") ;;
                  esac],[SILO_ENABLED=yes])

AM_CONDITIONAL([SILO_ENABLED],[test "$SILO_ENABLED" = yes])

AC_ARG_WITH([silo],
  AS_HELP_STRING(--with-silo=PATH,location of optional Silo installation),
  [if test "$SILO_ENABLED" = no ; then
     AC_MSG_WARN([--with-silo is specified, but support for Silo is disabled])
   else
     if test ! -d "$withval" ; then
       AC_MSG_ERROR([it is necessary to specify an existing directory when using --with-silo=PATH])
     fi
     SILO_DIR=$withval
   fi])

if test "$SILO_ENABLED" = yes; then
  if test x$SILO_DIR != x ; then
    if test -d "${SILO_DIR}/include" ; then
      SILO_CPPFLAGS="-I${SILO_DIR}/include"
    fi
    if test -d "${SILO_DIR}/lib" ; then
      SILO_LDFLAGS="-L${SILO_DIR}/lib"
    fi
  fi

  CPPFLAGS_PREPEND($SILO_CPPFLAGS)
  AC_CHECK_HEADER([silo.h],,AC_MSG_ERROR([Silo enabled but could not find working silo.h]))

  AC_MSG_CHECKING([for Silo version >= 4.9.1])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <silo.h>
  ]], [[
#if SILO_VERSION_GE(4,9,1)
#else
asdf
#endif
  ]])],[SILO_VERSION_VALID=yes],[SILO_VERSION_VALID=no])
  AC_MSG_RESULT([${SILO_VERSION_VALID}])
  if test "$SILO_VERSION_VALID" = no; then
    AC_MSG_WARN([Silo versions prior to 4.9.1 are likely to be usable but are not officially supported])
    AC_MSG_WARN([suggest upgrading to Silo 4.9.1])
  fi

  LDFLAGS_PREPEND($SILO_LDFLAGS)
  AC_LIB_HAVE_LINKFLAGS([silo])
  if test "$HAVE_LIBSILO" = no ; then
    AC_LIB_HAVE_LINKFLAGS([siloh5])
    if test "$HAVE_LIBSILOH5" = no ; then
      AC_MSG_ERROR([Silo enabled but could not find working libsilo or libsiloh5])
    fi
  fi

  PACKAGE_CPPFLAGS_PREPEND("$SILO_CPPFLAGS")
  PACKAGE_LDFLAGS_PREPEND("$SILO_LDFLAGS")
  if test "$HAVE_LIBSILO" = yes ; then
    PACKAGE_LIBS_PREPEND("$LIBSILO")
  fi
  if test "$HAVE_LIBSILOH5" = yes ; then
    PACKAGE_LIBS_PREPEND("$LIBSILOH5")
  fi
  AC_DEFINE([HAVE_SILO],1,[Define if you have the silo library.])
else
  AC_MSG_NOTICE([Optional package Silo is DISABLED])
fi

PACKAGE_RESTORE_ENVIRONMENT

])
