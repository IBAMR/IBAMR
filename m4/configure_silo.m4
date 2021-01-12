## ---------------------------------------------------------------------
##
## Copyright (c) 2014 - 2020 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

AC_DEFUN([CONFIGURE_SILO],[
echo
echo "================================="
echo "Configuring optional package Silo"
echo "================================="

AC_ARG_ENABLE([silo],
  AS_HELP_STRING(--enable-silo,enable support for the optional Silo library @<:@default=yes@:>@),
                 [case "$enableval" in
                    yes)  SILO_ENABLED=yes ;;
                    no)   SILO_ENABLED=no ;;
                    *)    AC_MSG_ERROR(--enable-silo=$enableval is invalid; choices are "yes" and "no") ;;
                  esac],[SILO_ENABLED=yes])

# while not an immediate dependency of IBAMR, if Silo is statically linked then
# users may need to provide paths to its dependencies (just libz)
AC_ARG_WITH([zlib],
  AS_HELP_STRING(--with-zlib=PATH,location of libz (a Silo dependency)),
  [if test "$SILO_ENABLED" = no ; then
     AC_MSG_WARN([--with-zlib (a path to a silo dependency) is specified, but support for Silo is disabled])
   else
     if test ! -d "$withval" ; then
       AC_MSG_ERROR([it is necessary to specify an existing directory when using --with-zlib=PATH])
     fi
     ZLIB_DIR=$withval
   fi])

# Back to silo:
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
    else
      AC_MSG_ERROR([Unable to find lib directory for silo: ${SILO_DIR}/lib does not exist.])
    fi
  fi

  # handle silo's dependencies:
  if test "x$ZLIB_DIR" = "x" ; then
      :
  else
      LDFLAGS_PREPEND(-L$ZLIB_DIR/lib)
      ADD_RPATH_LDFLAG(${ZLIB_DIR}/lib)
  fi
  AC_SEARCH_LIBS([compress2], [z], [],
                 [AC_MSG_ERROR([Unable to locate a valid zlib installation, which is a dependency of Silo. Since Silo is usually statically linked this must be provided to IBAMR. If zlib is not installed in a standard location then a path to it must be provided via --with-zlib=PATH.])])

  # now handle silo:
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
  AC_SEARCH_LIBS([DBSetDir], [silo siloh5], [],
                 [AC_MSG_ERROR([Silo enabled but could not find working libsilo or libsiloh5])])
  # set up rpath
  ADD_RPATH_LDFLAG(${SILO_DIR}/lib)

  AC_DEFINE([HAVE_SILO],1,[Define if you have the silo library.])
else
  AC_MSG_NOTICE([Optional package Silo is DISABLED])
fi

])
