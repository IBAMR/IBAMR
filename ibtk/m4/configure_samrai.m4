## ---------------------------------------------------------------------
##
## Copyright (c) 2014 - 2019 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

AC_DEFUN([CONFIGURE_SAMRAI],[
echo
echo "==================================="
echo "Configuring required package SAMRAI"
echo "==================================="

PACKAGE_SETUP_ENVIRONMENT

AC_ARG_WITH([samrai],
  AS_HELP_STRING(--with-samrai=PATH,location of required SAMRAI installation),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([it is necessary to specify an existing directory when using --with-samrai=PATH])
   fi
   SAMRAI_DIR=$withval])

if test x$SAMRAI_DIR != x ; then
  if test -d "${SAMRAI_DIR}/include" ; then
    SAMRAI_CPPFLAGS="-I${SAMRAI_DIR}/include"
  fi
  if test -d "${SAMRAI_DIR}/lib" ; then
    SAMRAI_LDFLAGS="-L${SAMRAI_DIR}/lib"
  fi
fi

AC_SUBST(SAMRAI_DIR,[${SAMRAI_DIR}])
AC_SUBST(SAMRAI_FORTDIR,[${SAMRAI_DIR}/include])

CPPFLAGS_PREPEND($SAMRAI_CPPFLAGS)
AC_CHECK_HEADER([SAMRAI_config.h],,AC_MSG_ERROR([could not find header file SAMRAI_config.h]))

AC_MSG_CHECKING([whether SAMRAI is configured with debugging enabled])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <SAMRAI_config.h>
]], [[
#if defined(DEBUG_CHECK_ASSERTIONS)
asdf
#endif
]])],[
SAMRAI_DEBUG_CHECK_ASSERTIONS=false
AC_MSG_RESULT(no)],[
SAMRAI_DEBUG_CHECK_ASSERTIONS=true
AC_MSG_RESULT(yes)
])
if test "$SAMRAI_DEBUG_CHECK_ASSERTIONS" = false; then
  SAMRAI_CPPFLAGS="$SAMRAI_CPPFLAGS -DNDEBUG"
else
  SAMRAI_CPPFLAGS="$SAMRAI_CPPFLAGS -DDEBUG"
fi

LDFLAGS_PREPEND($SAMRAI_LDFLAGS)
AC_LIB_HAVE_LINKFLAGS([SAMRAI])
if test "$HAVE_LIBSAMRAI" = no ; then
  AC_MSG_ERROR([could not find working libSAMRAI])
fi
SAMRAI_LIBS="$LIBSAMRAI"

AC_ARG_ENABLE([samrai-2d],
  AS_HELP_STRING(--enable-samrai-2d,enable optional support for two-dimensional SAMRAI objects @<:@default=yes@:>@),
                 [case "$enableval" in
                    yes)  SAMRAI2D_ENABLED=yes ;;
                    no)   SAMRAI2D_ENABLED=no ;;
                    *)    AC_MSG_ERROR(--enable-samrai-2d=$enableval is invalid; choices are "yes" and "no") ;;
                  esac],[SAMRAI2D_ENABLED=yes])
AM_CONDITIONAL([SAMRAI2D_ENABLED],[test "$SAMRAI2D_ENABLED" = yes])

if test "$SAMRAI2D_ENABLED" = yes; then
  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_hier])
  if test "$HAVE_LIBSAMRAI2D_HIER" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI2d_hier])
  fi
  SAMRAI_LIBS="$LIBSAMRAI2D_HIER $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_xfer])
  if test "$HAVE_LIBSAMRAI2D_XFER" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI2d_xfer])
  fi
  SAMRAI_LIBS="$LIBSAMRAI2D_XFER $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_pdat_std])
  if test "$HAVE_LIBSAMRAI2D_PDAT_STD" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI2d_pdat_std])
  fi
  SAMRAI_LIBS="$LIBSAMRAI2D_PDAT_STD $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_math_std])
  if test "$HAVE_LIBSAMRAI2D_MATH_STD" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI2d_math_std])
  fi
  SAMRAI_LIBS="$LIBSAMRAI2D_MATH_STD $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_mesh])
  if test "$HAVE_LIBSAMRAI2D_MESH" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI2d_mesh])
  fi
  SAMRAI_LIBS="$LIBSAMRAI2D_MESH $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_geom])
  if test "$HAVE_LIBSAMRAI2D_GEOM" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI2d_geom])
  fi
  SAMRAI_LIBS="$LIBSAMRAI2D_GEOM $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_solv])
  if test "$HAVE_LIBSAMRAI2D_SOLV" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI2d_solv])
  fi
  SAMRAI_LIBS="$LIBSAMRAI2D_SOLV $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_algs])
  if test "$HAVE_LIBSAMRAI2D_ALGS" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI2d_algs])
  fi
  SAMRAI_LIBS="$LIBSAMRAI2D_ALGS $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_appu])
  if test "$HAVE_LIBSAMRAI2D_APPU" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI2d_appu])
  fi
  SAMRAI_LIBS="$LIBSAMRAI2D_APPU $SAMRAI_LIBS"
else
  AC_MSG_NOTICE([configuring without the two-dimensional SAMRAI library])
fi

AC_ARG_ENABLE([samrai-3d],
  AS_HELP_STRING(--enable-samrai-3d,enable optional support for three-dimensional SAMRAI objects @<:@default=yes@:>@),
  [SAMRAI3D=$enablevar], [SAMRAI3D=yes])
AM_CONDITIONAL([SAMRAI3D_ENABLED],[test "$SAMRAI3D" = yes])

if test "$SAMRAI3D" = yes; then
  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_hier])
  if test "$HAVE_LIBSAMRAI3D_HIER" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI3d_hier])
  fi
  SAMRAI_LIBS="$LIBSAMRAI3D_HIER $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_xfer])
  if test "$HAVE_LIBSAMRAI3D_XFER" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI3d_xfer])
  fi
  SAMRAI_LIBS="$LIBSAMRAI3D_XFER $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_pdat_std])
  if test "$HAVE_LIBSAMRAI3D_PDAT_STD" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI3d_pdat_std])
  fi
  SAMRAI_LIBS="$LIBSAMRAI3D_PDAT_STD $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_math_std])
  if test "$HAVE_LIBSAMRAI3D_MATH_STD" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI3d_math_std])
  fi
  SAMRAI_LIBS="$LIBSAMRAI3D_MATH_STD $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_mesh])
  if test "$HAVE_LIBSAMRAI3D_MESH" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI3d_mesh])
  fi
  SAMRAI_LIBS="$LIBSAMRAI3D_MESH $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_geom])
  if test "$HAVE_LIBSAMRAI3D_GEOM" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI3d_geom])
  fi
  SAMRAI_LIBS="$LIBSAMRAI3D_GEOM $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_solv])
  if test "$HAVE_LIBSAMRAI3D_SOLV" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI3d_solv])
  fi
  SAMRAI_LIBS="$LIBSAMRAI3D_SOLV $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_algs])
  if test "$HAVE_LIBSAMRAI3D_ALGS" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI3d_algs])
  fi
  SAMRAI_LIBS="$LIBSAMRAI3D_ALGS $SAMRAI_LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_appu])
  if test "$HAVE_LIBSAMRAI3D_APPU" = no ; then
    AC_MSG_ERROR([could not find working libSAMRAI3d_appu])
  fi
  SAMRAI_LIBS="$LIBSAMRAI3D_APPU $SAMRAI_LIBS"
else
  AC_MSG_NOTICE([configuring without the three-dimensional SAMRAI library])
fi

PACKAGE_CPPFLAGS_PREPEND($SAMRAI_CPPFLAGS)
PACKAGE_LDFLAGS_PREPEND($SAMRAI_LDFLAGS)
PACKAGE_LIBS_PREPEND($SAMRAI_LIBS)
PACKAGE_RESTORE_ENVIRONMENT
])
