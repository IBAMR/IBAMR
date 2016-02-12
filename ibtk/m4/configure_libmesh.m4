# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_LIBMESH],[
echo
echo "===================================="
echo "Configuring optional package libMesh"
echo "===================================="

PACKAGE_SETUP_ENVIRONMENT

AC_ARG_ENABLE([libmesh],
  AS_HELP_STRING(--enable-libmesh,enable support for the optional libMesh library @<:@default=yes@:>@),
                 [case "$enableval" in
                    yes)  LIBMESH_ENABLED=yes ;;
                    no)   LIBMESH_ENABLED=no ;;
                    *)    AC_MSG_ERROR(--enable-libmesh=$enableval is invalid; choices are "yes" and "no") ;;
                  esac],[LIBMESH_ENABLED=no])

AM_CONDITIONAL([LIBMESH_ENABLED],[test "$LIBMESH_ENABLED" = yes])

AC_ARG_WITH([libmesh],
  AS_HELP_STRING(--with-libmesh=PATH,location of optional libMesh installation),
  [if test "$LIBMESH_ENABLED" = no ; then
     AC_MSG_WARN([--with-libmesh is specified, but support for libMesh is disabled])
   else
     if test ! -d "$withval" ; then
       AC_MSG_ERROR([it is necessary to specify an existing directory when using --with-libmesh=PATH])
     fi
     LIBMESH_DIR=$withval
   fi])

AC_ARG_WITH([libmesh-method],
  AS_HELP_STRING([--with-libmesh-method=METHOD],[libMesh compile mode options include opt, devel, dbg, prof, oprof @<:@default=opt@:>@]),
  [if test "$LIBMESH_ENABLED" = no ; then
     AC_MSG_WARN([--with-libmesh-method is specified, but support for libMesh is disabled])
     LIBMESH_METHOD=NONE
   fi
   LIBMESH_METHOD=$withval],
  [LIBMESH_METHOD=opt])

if test "$LIBMESH_ENABLED" = yes; then
  if test x$LIBMESH_DIR != x ; then
    if test -d "${LIBMESH_DIR}/include" ; then
      LIBMESH_CPPFLAGS="-I${LIBMESH_DIR}/include"
    fi
    if test -d "${LIBMESH_DIR}/bin" ; then
      LIBMESH_BIN="${LIBMESH_DIR}/bin"
    fi
  fi

  METHOD=$LIBMESH_METHOD
  CPPFLAGS_PREPEND($LIBMESH_CPPFLAGS)
  AC_CHECK_HEADER([libmesh/libmesh.h],,AC_MSG_ERROR([libMesh enabled but could not find working libmesh.h]))
  AC_CHECK_HEADER([libmesh/libmesh_config.h],,AC_MSG_ERROR([libMesh enabled but could not find working libmesh_config.h]))
  AC_MSG_CHECKING([for libMesh version 0.9.5 or newer])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <libmesh/libmesh_config.h>

#if ((LIBMESH_MAJOR_VERSION > 0) || ((LIBMESH_MINOR_VERSION >= 9) && (LIBMES_SUBMINOR_VERSION >= 5)))
#else
asdf
#endif
  ]])],[LIBMESH_VERSION_VALID=yes],[LIBMESH_VERSION_VALID=no])
  AC_MSG_RESULT([${LIBMESH_VERSION_VALID}])
  if test "$LIBMESH_VERSION_VALID" = no; then
    AC_MSG_ERROR([invalid libMesh version detected: please use libMesh 0.9.5 or newer])
  fi
  AC_MSG_NOTICE([obtaining libMesh configuration information from libmesh_common.h])
  AC_RUN_IFELSE([AC_LANG_SOURCE([
#include "libmesh/libmesh_config.h"
#include <iostream>
#define STRINGIFY(str) #str
#define EXPAND_AND_STRINGIFY(str) STRINGIFY(str)
int main()
{
   std::cout << EXPAND_AND_STRINGIFY(LIBMESH_CONFIGURE_INFO) << std::endl;
   return 0;
}
  ])],[
  LIBMESH_CONFIGURE_INFO=`./conftest$EXEEXT | sed -e 's/^"//' -e 's/"$//'`
  for elem in $LIBMESH_CONFIGURE_INFO ; do
    if test `echo $elem | grep -c "PETSC_DIR="` != 0 ; then
      temp="${elem%\'}"
      elem_strip="${temp#\'}"
      LIBMESH_PETSC_DIR=`echo $elem_strip | sed -e 's/.*=//' -e 's/^[ \t]*//'`
    elif test `echo $elem | grep -c "PETSC_ARCH="` != 0 ; then
      temp="${elem%\'}"
      elem_strip="${temp#\'}"
      LIBMESH_PETSC_ARCH=`echo $elem_strip | sed -e 's/.*=//' -e 's/^[ \t]*//'`
    fi
  done
  if test "$PETSC_DIR" != "$LIBMESH_PETSC_DIR" ; then
    AC_MSG_NOTICE([using libMesh PETSC_DIR=${LIBMESH_PETSC_DIR} instead of PETSC_DIR=${PETSC_DIR}])
    PETSC_DIR=$LIBMESH_PETSC_DIR
  fi
  if test "$PETSC_ARCH" != "$LIBMESH_PETSC_ARCH" ; then
    AC_MSG_NOTICE([using libMesh PETSC_ARCH=${LIBMESH_PETSC_ARCH} instead of PETSC_ARCH=${PETSC_ARCH}])
    PETSC_ARCH=$LIBMESH_PETSC_ARCH
  fi
  ],[
  AC_MSG_ERROR([could not execute program to examine settings in libmesh_config.h])])

  AC_PATH_PROG(LIBMESH_CONFIG, libmesh-config, , [$PATH$PATH_SEPARATOR$LIBMESH_BIN])
  if test -z "$LIBMESH_CONFIG" ; then
    AC_MSG_ERROR([cannot find libmesh-config])
  fi

  LIBMESH_CONFIG="env METHOD=$METHOD $LIBMESH_CONFIG"

  LIBMESH_CPPFLAGS="`$LIBMESH_CONFIG --include` `$LIBMESH_CONFIG --cppflags`"
  LIBMESH_CXXFLAGS="`$LIBMESH_CONFIG --cxxflags`"
  LIBMESH_CFLAGS="`$LIBMESH_CONFIG --cflags`"
  LIBMESH_FCFLAGS="`$LIBMESH_CONFIG --fflags`"
  LIBMESH_LIBS="`$LIBMESH_CONFIG --libs`"

  AC_LIB_HAVE_LINKFLAGS([netcdf])
  if test "$HAVE_LIBNETCDF" = yes ; then
    LIBMESH_LIBS="$LIBMESH_LIBS $LIBNETCDF"
  fi

  if test -e "$LIBMESH_DIR/include/Eigen" ; then
    AC_MSG_NOTICE([using Eigen library bundled with libMesh])
    with_eigen="$LIBMESH_DIR"
  fi

  if test -e "$LIBMESH_DIR/include/boost" ; then
    AC_MSG_ERROR([libMesh must be configured to use an external boost library])
  fi

  CPPFLAGS_PREPEND($LIBMESH_CPPFLAGS)
  CXXFLAGS_PREPEND($LIBMESH_CXXFLAGS)
  CFLAGS_PREPEND($LIBMESH_CFLAGS)
  FCFLAGS_PREPEND($LIBMESH_FCFLAGS)
  LIBS_PREPEND($LIBMESH_LIBS)

  case "$METHOD" in
    opt)
      AC_LIB_HAVE_LINKFLAGS([mesh_opt])
      if test "$HAVE_LIBMESH_OPT" = "no" ; then
        AC_MSG_ERROR([libMesh enabled with METHOD=$METHOD but could not find working libmesh_opt])
      fi
      ;;
    devel)
      AC_LIB_HAVE_LINKFLAGS([mesh_devel])
      if test "$HAVE_LIBMESH_DEVEL" = "no" ; then
        AC_MSG_ERROR([libMesh enabled with METHOD=$METHOD but could not find working libmesh_devel])
      fi
      ;;
    dbg)
      AC_LIB_HAVE_LINKFLAGS([mesh_dbg])
      if test "$HAVE_LIBMESH_DBG" = "no" ; then
        AC_MSG_ERROR([libMesh enabled with METHOD=$METHOD but could not find working libmesh_dbg])
      fi
      ;;
    prof)
      AC_LIB_HAVE_LINKFLAGS([mesh_prof])
      if test "$HAVE_LIBMESH_PROF" = "no" ; then
        AC_MSG_ERROR([libMesh enabled with METHOD=$METHOD but could not find working libmesh_prof])
      fi
      ;;
    oprof)
      AC_LIB_HAVE_LINKFLAGS([mesh_oprof])
      if test "$HAVE_LIBMESH_OPROF" = "no" ; then
        AC_MSG_ERROR([libMesh enabled with METHOD=$METHOD but could not find working libmesh_oprof])
      fi
      ;;
    *)
      AC_MSG_ERROR("unknown libMesh METHOD=$METHOD; options are: opt, devel, dbg, prof, oprof") ;;
  esac
  AC_DEFINE([HAVE_LIBMESH],1,[Define if you have the libmesh library.])

  PACKAGE_CPPFLAGS_PREPEND($LIBMESH_CPPFLAGS)
  PACKAGE_CXXFLAGS_PREPEND($LIBMESH_CXXFLAGS)
  PACKAGE_CFLAGS_PREPEND($LIBMESH_CFLAGS)
  PACKAGE_FCFLAGS_PREPEND($LIBMESH_FCFLAGS)
  PACKAGE_LIBS_PREPEND($LIBMESH_LIBS)
else
  AC_MSG_NOTICE([Optional package libMesh is DISABLED])
fi

PACKAGE_RESTORE_ENVIRONMENT

])
