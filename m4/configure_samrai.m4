# -------------------------------------------------------------
# -------------------------------------------------------------
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
AC_CHECK_HEADER([SAMRAI/SAMRAI_config.h],,AC_MSG_ERROR([could not find header file SAMRAI/SAMRAI_config.h]))

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

AC_ARG_ENABLE([samrai-2d],
  AS_HELP_STRING(--enable-samrai-2d,enable optional support for two-dimensional SAMRAI objects @<:@default=yes@:>@),
                 [case "$enableval" in
                    yes)  SAMRAI2D_ENABLED=yes ;;
                    no)   SAMRAI2D_ENABLED=no ;;
                    *)    AC_MSG_ERROR(--enable-samrai-2d=$enableval is invalid; choices are "yes" and "no") ;;
                  esac],[SAMRAI2D_ENABLED=yes])
AM_CONDITIONAL([SAMRAI2D_ENABLED],[test "$SAMRAI2D_ENABLED" = yes])

AC_ARG_ENABLE([samrai-3d],
  AS_HELP_STRING(--enable-samrai-3d,enable optional support for three-dimensional SAMRAI objects @<:@default=yes@:>@),
  [SAMRAI3D=$enablevar], [SAMRAI3D=yes])
AM_CONDITIONAL([SAMRAI3D_ENABLED],[test "$SAMRAI3D" = yes])


AC_LIB_HAVE_LINKFLAGS([SAMRAI_hier])
if test "$HAVE_LIBSAMRAI_HIER" = no ; then
  AC_MSG_ERROR([could not find working libSAMRAI_hier])
fi
SAMRAI_LIBS="$LIBSAMRAI_HIER $SAMRAI_LIBS"

AC_LIB_HAVE_LINKFLAGS([SAMRAI_xfer])
if test "$HAVE_LIBSAMRAI_XFER" = no ; then
  AC_MSG_ERROR([could not find working libSAMRAI_xfer])
fi
SAMRAI_LIBS="$LIBSAMRAI_XFER $SAMRAI_LIBS"

AC_LIB_HAVE_LINKFLAGS([SAMRAI_pdat])
if test "$HAVE_LIBSAMRAI_PDAT" = no ; then
  AC_MSG_ERROR([could not find working libSAMRAI_pdat])
fi
SAMRAI_LIBS="$LIBSAMRAI_PDAT $SAMRAI_LIBS"

AC_LIB_HAVE_LINKFLAGS([SAMRAI_math])
if test "$HAVE_LIBSAMRAI_MATH" = no ; then
  AC_MSG_ERROR([could not find working libSAMRAI_math])
fi
SAMRAI_LIBS="$LIBSAMRAI_MATH $SAMRAI_LIBS"

AC_LIB_HAVE_LINKFLAGS([SAMRAI_mesh])
if test "$HAVE_LIBSAMRAI_MESH" = no ; then
  AC_MSG_ERROR([could not find working libSAMRAI_mesh])
fi
SAMRAI_LIBS="$LIBSAMRAI_MESH $SAMRAI_LIBS"

AC_LIB_HAVE_LINKFLAGS([SAMRAI_geom])
if test "$HAVE_LIBSAMRAI_GEOM" = no ; then
  AC_MSG_ERROR([could not find working libSAMRAI_geom])
fi
SAMRAI_LIBS="$LIBSAMRAI_GEOM $SAMRAI_LIBS"

AC_LIB_HAVE_LINKFLAGS([SAMRAI_solv])
if test "$HAVE_LIBSAMRAI_SOLV" = no ; then
  AC_MSG_ERROR([could not find working libSAMRAI_solv])
fi
SAMRAI_LIBS="$LIBSAMRAI_SOLV $SAMRAI_LIBS"

AC_LIB_HAVE_LINKFLAGS([SAMRAI_algs])
if test "$HAVE_LIBSAMRAI_ALGS" = no ; then
  AC_MSG_ERROR([could not find working libSAMRAI_algs])
fi
SAMRAI_LIBS="$LIBSAMRAI_ALGS $SAMRAI_LIBS"

AC_LIB_HAVE_LINKFLAGS([SAMRAI_appu])
if test "$HAVE_LIBSAMRAI_APPU" = no ; then
  AC_MSG_ERROR([could not find working libSAMRAI_appu])
fi
SAMRAI_LIBS="$LIBSAMRAI_APPU $SAMRAI_LIBS"

PACKAGE_CPPFLAGS_PREPEND($SAMRAI_CPPFLAGS)
PACKAGE_LDFLAGS_PREPEND($SAMRAI_LDFLAGS)
PACKAGE_LIBS_PREPEND($SAMRAI_LIBS)
PACKAGE_RESTORE_ENVIRONMENT
])
