## ---------------------------------------------------------------------
##
## Copyright (c) 2014 - 2021 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

AC_DEFUN([CONFIGURE_LIBMESH],[
echo
echo "===================================="
echo "Configuring optional package libMesh"
echo "===================================="

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
dnl
dnl 1. Check headers:
dnl
  AC_CHECK_HEADER([libmesh/libmesh.h],,AC_MSG_ERROR([libMesh enabled but could not find working libmesh.h]))
  AC_CHECK_HEADER([libmesh/libmesh_config.h],,AC_MSG_ERROR([libMesh enabled but could not find working libmesh_config.h]))
dnl
dnl 2. Check version numbers:
dnl
  AC_MSG_CHECKING([for libMesh version 1.1.0 or newer])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <libmesh/libmesh_config.h>

#if LIBMESH_MAJOR_VERSION >= 1 && LIBMESH_MINOR_VERSION >= 1
// OK
#else
#error
#endif
  ]])],[LIBMESH_VERSION_VALID=yes],[LIBMESH_VERSION_VALID=no])
  AC_MSG_RESULT([${LIBMESH_VERSION_VALID}])
  if test "$LIBMESH_VERSION_VALID" = no; then
    AC_MSG_ERROR([invalid libMesh version detected: please use libMesh 1.1.0 or newer])
  fi
  AC_MSG_NOTICE([obtaining libMesh configuration information from libmesh_common.h])
dnl
dnl 3. Check that libMesh was compiled with PETSc:
dnl
  AC_MSG_CHECKING([for libMesh with PETSc])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <libmesh/libmesh_config.h>

#ifdef LIBMESH_HAVE_PETSC
// OK
#else
#error
#endif
  ]])],[LIBMESH_HAVE_PETSC=yes],[LIBMESH_HAVE_PETSC=no])
  AC_MSG_RESULT([${LIBMESH_HAVE_PETSC}])
  if test "$LIBMESH_HAVE_PETSC" = no; then
    AC_MSG_ERROR([invalid libMesh installation detected: please compile libMesh with PETSc])
  fi
dnl
dnl 4. Make sure that libMesh's PETSc matches the one specified. The linker can
dnl    sometimes get confused by this (the linker claims it cannot find
dnl    libpetsc.so.3.10, which is utterly irrelevant to looking at the
dnl    configuration file), so clear LIBS for the moment:
dnl
  _old_libs="${LIBS}"
  LIBS=""
  AC_RUN_IFELSE([AC_LANG_SOURCE([
#include "libmesh/libmesh_config.h"
#include <iostream>
int main()
{
   std::cout << LIBMESH_CONFIGURE_INFO << std::endl;
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
dnl
dnl 4a. Reset LIBS to its previous value:
dnl
  LIBS="${_old_libs}"

  AC_PATH_PROG(LIBMESH_CONFIG, libmesh-config, , [$PATH$PATH_SEPARATOR$LIBMESH_BIN])
  if test -z "$LIBMESH_CONFIG" ; then
    AC_MSG_ERROR([cannot find libmesh-config])
  fi

  LIBMESH_CONFIG="env METHOD=$METHOD $LIBMESH_CONFIG"

  LIBMESH_CPPFLAGS="`$LIBMESH_CONFIG --include` `$LIBMESH_CONFIG --cppflags`"
  LIBMESH_LIBS="`$LIBMESH_CONFIG --libs`"

dnl
dnl 5. libMesh (for pre-C++11 compatibility) can implement libMesh::UniquePtr
dnl    in several different ways. We are only compatible when libMesh::UniquePtr
dnl    really is std::unique_ptr:
dnl
  AC_MSG_CHECKING([for a usable libMesh UniquePtr/unique_ptr configuration])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <libmesh/auto_ptr.h>
#include <libmesh/libmesh_version.h>

#include <type_traits>

// This check is not necessary for newer versions which completely remove
// UniquePtr
#if LIBMESH_VERSION_LESS_THAN(1, 7, 0)
static_assert(std::is_same<libMesh::UniquePtr<int>, std::unique_ptr<int>>::value,
              "libMesh should use std::unique_ptr");
#else
// OK
#endif
  ]])],[LIBMESH_UNIQUE_PTR_OK=yes],[LIBMESH_UNIQUE_PTR_OK=no])
  AC_MSG_RESULT([${LIBMESH_UNIQUE_PTR_OK}])
  if test "$LIBMESH_UNIQUE_PTR_OK" = no; then
    AC_MSG_ERROR([libMesh must be compiled with C++11 support and without --disable-unique-ptr.])
  fi
dnl
dnl 6. Check bundled libMesh components:
dnl
  if test -e "$LIBMESH_DIR/include/Eigen" ; then
    AC_MSG_NOTICE([using Eigen library bundled with libMesh])
    with_eigen="$LIBMESH_DIR"
  fi

dnl try to see if libMesh has boost:
  AC_MSG_CHECKING([for libMesh with boost])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <libmesh/libmesh_config.h>

#ifdef LIBMESH_HAVE_BOOST
// ok
#else
#error
#endif ]])], [LIBMESH_HAVE_BOOST=yes], [LIBMESH_HAVE_BOOST=no])
  AC_MSG_RESULT([${LIBMESH_HAVE_BOOST}])

dnl: "''
dnl Older versions of libMesh do not define LIBMESH_HAVE_EXTERNAL_BOOST so we check manually
dnl to see if a sufficient version of boost is bundled or not
  if test -e "$LIBMESH_DIR/include/boost/"; then
      if test -e "$LIBMESH_DIR/include/boost/multi_array.hpp"; then
          # ok: perhaps libMesh is installed in /usr/
          :
      else
        AC_MSG_ERROR([The directory $LIBMESH_DIR/include/boost exists and contains a boost installation that does not provide all the headers that IBAMR needs. This can happen when either libMesh is installed with its own bundled version of boost (which is not compatible with IBAMR) or when libMesh is installed into a directory which happens to contain a copy of boost (e.g., when libMesh is recompiled without boost and installed into the same location as a previous copy of libMesh). If you want to use libMesh with boost then both libMesh and IBAMR must use the same external copy of boost. The best way to fix this problem is to delete $LIBMESH_DIR and reinstall libMesh with either no boost support or an external boost library.])
      fi
  fi

  CPPFLAGS_PREPEND($LIBMESH_CPPFLAGS)
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
else
  AC_MSG_NOTICE([Optional package libMesh is DISABLED])
fi

])
