AC_DEFUN([CONFIGURE_LIBMESH],[

AC_ARG_ENABLE([libmesh],
  AS_HELP_STRING(--enable-libmesh,enable support for the optional libMesh library @<:@default=no@:>@),
  [LIBMESH_ENABLED=$enableval], [LIBMESH_ENABLED=no])
if test "$LIBMESH_ENABLED" = "" ; then
  LIBMESH_ENABLED=no
fi

AC_ARG_WITH([libmesh],
  AS_HELP_STRING([--with-libmesh=PATH],[location of optional libMesh installation @<:@default=/usr/local/sfw/libmesh/current@:>@]),
  [if test "$LIBMESH_ENABLED" = no ; then
     AC_MSG_WARN([--with-libmesh is specified, but support for libMesh is DISABLED])
     LIBMESH_DIR=NONE
   else
     if test ! -d "$withval" ; then
       AC_MSG_ERROR([you must specify an existing directory when using --with-libmesh=PATH])
     fi
     LIBMESH_DIR=$withval
   fi],
  [LIBMESH_DIR=/usr/local/sfw/libmesh/current])

AC_ARG_WITH([libmesh-compile-mode],
  AS_HELP_STRING([--with-libmesh-compile-mode=METHOD],[libMesh compile mode options include opt, devel, dbg, prof, oprof @<:@default=dbg@:>@]),
  [if test "$LIBMESH_ENABLED" = no ; then
     AC_MSG_WARN([--with-libmesh-compile-mode is specified, but support for libMesh is DISABLED])
     LIBMESH_COMPILE_MODE=NONE
   else
     LIBMESH_COMPILE_MODE=$withval
   fi],
  [LIBMESH_COMPILE_MODE=dbg])


LIBMESH_HOSTTYPE=`grep "hosttype[ \t]*.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
LIBMESH_LIBDIR=${LIBMESH_DIR}/lib/${LIBMESH_HOSTTYPE}_${LIBMESH_COMPILE_MODE}
LIBMESH_CONTRIBDIR=${LIBMESH_DIR}/contrib
LIBMESH_CONTRIB_LIBDIR=${LIBMESH_DIR}/contrib/lib/${LIBMESH_HOSTTYPE}_${LIBMESH_COMPILE_MODE}

if test "${LIBMESH_COMPILE_MODE}" = "opt" ; then
  CPPFLAGS="${CPPFLAGS} `grep -A3 opt-mode $LIBMESH_DIR/Make.common | grep libmesh_CPPFLAGS | sed -e 's/.*libmesh_CPPFLAGS +=//' | sed -e 's/.*libmesh_CPPFLAGS :=//' | sed -e 's/$(libmesh_CPPFLAGS)//'`"
  CXXFLAGS="${CXXFLAGS} `grep -A3 opt-mode $LIBMESH_DIR/Make.common | grep libmesh_CXXFLAGS | sed -e 's/.*libmesh_CXXFLAGS +=//' | sed -e 's/.*libmesh_CXXFLAGS :=//' | sed -e 's/$(libmesh_CXXFLAGS)//'`"
  CFLAGS="${CFLAGS}     `grep -A3 opt-mode $LIBMESH_DIR/Make.common | grep libmesh_CFLAGS   | sed -e 's/.*libmesh_CFLAGS   +=//' | sed -e 's/.*libmesh_CFLAGS   :=//' | sed -e 's/$(libmesh_CFLAGS)//'  `"
fi
if test "${LIBMESH_COMPILE_MODE}" = "devel" ; then
  CPPFLAGS="${CPPFLAGS} `grep -A3 devel-mode $LIBMESH_DIR/Make.common | grep libmesh_CPPFLAGS | sed -e 's/.*libmesh_CPPFLAGS +=//' | sed -e 's/.*libmesh_CPPFLAGS :=//' | sed -e 's/$(libmesh_CPPFLAGS)//'`"
  CXXFLAGS="${CXXFLAGS} `grep -A3 devel-mode $LIBMESH_DIR/Make.common | grep libmesh_CXXFLAGS | sed -e 's/.*libmesh_CXXFLAGS +=//' | sed -e 's/.*libmesh_CXXFLAGS :=//' | sed -e 's/$(libmesh_CXXFLAGS)//'`"
  CFLAGS="${CFLAGS}     `grep -A3 devel-mode $LIBMESH_DIR/Make.common | grep libmesh_CFLAGS   | sed -e 's/.*libmesh_CFLAGS   +=//' | sed -e 's/.*libmesh_CFLAGS   :=//' | sed -e 's/$(libmesh_CFLAGS)//'  `"
fi
if test "${LIBMESH_COMPILE_MODE}" = "dbg" ; then
  CPPFLAGS="${CPPFLAGS} `grep -A3 debug-mode $LIBMESH_DIR/Make.common | grep libmesh_CPPFLAGS | sed -e 's/.*libmesh_CPPFLAGS +=//' | sed -e 's/.*libmesh_CPPFLAGS :=//' | sed -e 's/$(libmesh_CPPFLAGS)//'`"
  CXXFLAGS="${CXXFLAGS} `grep -A3 debug-mode $LIBMESH_DIR/Make.common | grep libmesh_CXXFLAGS | sed -e 's/.*libmesh_CXXFLAGS +=//' | sed -e 's/.*libmesh_CXXFLAGS :=//' | sed -e 's/$(libmesh_CXXFLAGS)//'`"
  CFLAGS="${CFLAGS}     `grep -A3 debug-mode $LIBMESH_DIR/Make.common | grep libmesh_CFLAGS   | sed -e 's/.*libmesh_CFLAGS   +=//' | sed -e 's/.*libmesh_CFLAGS   :=//' | sed -e 's/$(libmesh_CFLAGS)//'  `"
fi
if test "${LIBMESH_COMPILE_MODE}" = "prof" ; then
  CPPFLAGS="${CPPFLAGS} `grep -A3 prof-mode $LIBMESH_DIR/Make.common | grep libmesh_CPPFLAGS | sed -e 's/.*libmesh_CPPFLAGS +=//' | sed -e 's/.*libmesh_CPPFLAGS :=//' | sed -e 's/$(libmesh_CPPFLAGS)//'`"
  CXXFLAGS="${CXXFLAGS} `grep -A3 prof-mode $LIBMESH_DIR/Make.common | grep libmesh_CXXFLAGS | sed -e 's/.*libmesh_CXXFLAGS +=//' | sed -e 's/.*libmesh_CXXFLAGS :=//' | sed -e 's/$(libmesh_CXXFLAGS)//'`"
  CFLAGS="${CFLAGS}     `grep -A3 prof-mode $LIBMESH_DIR/Make.common | grep libmesh_CFLAGS   | sed -e 's/.*libmesh_CFLAGS   +=//' | sed -e 's/.*libmesh_CFLAGS   :=//' | sed -e 's/$(libmesh_CFLAGS)//'  `"
fi
if test "${LIBMESH_COMPILE_MODE}" = "oprof" ; then
  CPPFLAGS="${CPPFLAGS} `grep -A3 oprof-mode $LIBMESH_DIR/Make.common | grep libmesh_CPPFLAGS | sed -e 's/.*libmesh_CPPFLAGS +=//' | sed -e 's/.*libmesh_CPPFLAGS :=//' | sed -e 's/$(libmesh_CPPFLAGS)//'`"
  CXXFLAGS="${CXXFLAGS} `grep -A3 oprof-mode $LIBMESH_DIR/Make.common | grep libmesh_CXXFLAGS | sed -e 's/.*libmesh_CXXFLAGS +=//' | sed -e 's/.*libmesh_CXXFLAGS :=//' | sed -e 's/$(libmesh_CXXFLAGS)//'`"
  CFLAGS="${CFLAGS}     `grep -A3 oprof-mode $LIBMESH_DIR/Make.common | grep libmesh_CFLAGS   | sed -e 's/.*libmesh_CFLAGS   +=//' | sed -e 's/.*libmesh_CFLAGS   :=//' | sed -e 's/$(libmesh_CFLAGS)//'  `"
fi

LIBMESH_TECPLOT_LIBS=`grep "libmesh_LIBS" $LIBMESH_DIR/Make.common | grep "tecplot" | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
LIBS="${LIBMESH_TECPLOT_LIBS} ${LIBS}"

LIBMESH_TECPLOT_INCLUDE=`grep "libmesh_INCLUDE" $LIBMESH_DIR/Make.common | grep "tecplot" | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
CPPFLAGS="${LIBMESH_TECPLOT_INCLUDE} ${CPPFLAGS}"

LIBMESH_PETSC_ARCH=`grep "PETSC_ARCH[ \t]*.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
if test "${PETSC_ARCH}" != "${LIBMESH_PETSC_ARCH}"; then
  AC_MSG_WARN([using libMesh PETSC_ARCH=${LIBMESH_PETSC_ARCH} instead of PETSC_ARCH=${PETSC_ARCH}])
  PETSC_ARCH=${LIBMESH_PETSC_ARCH}
fi

LIBMESH_PETSC_DIR=`grep "PETSC_DIR[ \t]*.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
if test "${PETSC_DIR}" != "${LIBMESH_PETSC_DIR}"; then
  AC_MSG_WARN([using libMesh PETSC_DIR=${LIBMESH_PETSC_DIR} instead of PETSC_DIR=${PETSC_DIR}])
  PETSC_DIR=${LIBMESH_PETSC_DIR}
fi

if test -d "${LIBMESH_LIBDIR}" ; then
  LDFLAGS="-L${LIBMESH_LIBDIR} $LDFLAGS"
else
  AC_MSG_WARN([libMesh support is enabled, but could not find expected libMesh lib directory ${LIBMESH_LIBDIR}])
fi

if test -d "${LIBMESH_CONTRIB_LIBDIR}" ; then
  LDFLAGS="-L${LIBMESH_CONTRIB_LIBDIR} $LDFLAGS"
else
  AC_MSG_WARN([libMesh support is enabled, but could not find expected libMesh contrib lib directory ${LIBMESH_CONTRIB_LIBDIR}])
fi

contrib_lib_enabled=`grep "enable-laspack.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
AC_MSG_NOTICE([checking whether libMesh contrib package laspack is enabled... ${contrib_lib_enabled}])
if test "${contrib_lib_enabled}" = yes ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/laspack ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([laspack])
  if test "$HAVE_LIBLASPACK" = yes ; then
    LIBS="$LIBLASPACK $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib liblaspack is enabled, but could not find working liblaspack])
  fi
fi

contrib_lib_enabled=`grep "enable-metis.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
AC_MSG_NOTICE([checking whether libMesh contrib package metis is enabled... ${contrib_lib_enabled}])
if test "${contrib_lib_enabled}" = yes ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/metis/Lib ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([GK])
  if test "$HAVE_LIBGK" = yes ; then
    LIBS="$LIBGK $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libmetis is enabled, but could not find working libGK])
  fi
  AC_LIB_HAVE_LINKFLAGS([metis])
  if test "$HAVE_LIBMETIS" = yes ; then
    LIBS="$LIBMETIS $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libmetis is enabled, but could not find working libmetis])
  fi
fi

contrib_lib_enabled=`grep "enable-parmetis.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
AC_MSG_NOTICE([checking whether libMesh contrib package parmetis is enabled... ${contrib_lib_enabled}])
if test "${contrib_lib_enabled}" = yes ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/parmetis/Lib ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([parmetis])
  if test "$HAVE_LIBPARMETIS" = yes ; then
    LIBS="$LIBPARMETIS $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libparmetis is enabled, but could not find working libparmetis])
  fi
fi

contrib_lib_enabled=`grep "enable-sfcurves.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
AC_MSG_NOTICE([checking whether libMesh contrib package sfcurves is enabled... ${contrib_lib_enabled}])
if test "${contrib_lib_enabled}" = yes ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/sfcurves ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([sfcurves])
  if test "$HAVE_LIBSFCURVES" = yes ; then
    LIBS="$LIBSFCURVES $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libsfcurves is enabled, but could not find working libsfcurves])
  fi
fi

contrib_lib_enabled=`grep "enable-gzstream.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
AC_MSG_NOTICE([checking whether libMesh contrib package gzstream is enabled... ${contrib_lib_enabled}])
if test "${contrib_lib_enabled}" = yes ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/gzstream ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([gzstream])
  if test "$HAVE_LIBGZSTREAM" = yes ; then
    LIBS="$LIBGZSTREAM $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libgzstream is enabled, but could not find working libgzstream])
  fi
fi

contrib_lib_enabled=`grep "enable-tetgen.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
AC_MSG_NOTICE([checking whether libMesh contrib package tetgen is enabled... ${contrib_lib_enabled}])
if test "${contrib_lib_enabled}" = yes ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/tetgen ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([tetgen])
  if test "$HAVE_LIBTETGEN" = yes ; then
    LIBS="$LIBTETGEN $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libtetgen is enabled, but could not find working libtetgen])
  fi
fi

contrib_lib_enabled=`grep "enable-triangle.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
AC_MSG_NOTICE([checking whether libMesh contrib package triangle is enabled... ${contrib_lib_enabled}])
if test "${contrib_lib_enabled}" = yes ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/triangle ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([triangle])
  if test "$HAVE_LIBTRIANGLE" = yes ; then
    LIBS="$LIBTRIANGLE $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libtriangle is enabled, but could not find working libtriangle])
  fi
fi

contrib_lib_enabled=`grep "enable-gmv.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
AC_MSG_NOTICE([checking whether libMesh contrib package gmv is enabled... ${contrib_lib_enabled}])
if test "${contrib_lib_enabled}" = yes ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/gmv ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([gmv])
  if test "$HAVE_LIBGMV" = yes ; then
    LIBS="$LIBGMV $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libgmv is enabled, but could not find working libgmv])
  fi
fi

contrib_lib_enabled=`grep "enable-vtk.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
AC_MSG_NOTICE([checking whether libMesh contrib package vtk is enabled... ${contrib_lib_enabled}])
if test "${contrib_lib_enabled}" = yes ; then
  AC_LIB_HAVE_LINKFLAGS([vtk])
  if test "$HAVE_LIBVTK" = yes ; then
    LIBS="$LIBVTK $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libvtk is enabled, but could not find working libvtk])
  fi
fi

contrib_lib_enabled=`grep "enable-netcdf.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
AC_MSG_NOTICE([checking whether libMesh contrib package netcdf is enabled... ${contrib_lib_enabled}])
if test "${contrib_lib_enabled}" = yes ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/netcdf/Lib ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([netcdf])
  if test "$HAVE_LIBNETCDF" = yes ; then
    LIBS="$LIBNETCDF $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libnetcdf is enabled, but could not find working libnetcdf])
  fi
fi

contrib_lib_enabled=`grep "enable-exodus.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
AC_MSG_NOTICE([checking whether libMesh contrib package exodus is enabled... ${contrib_lib_enabled}])
if test "${contrib_lib_enabled}" = yes ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/exodusii/Lib/include ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([exodusii])
  if test "$HAVE_LIBEXODUSII" = yes ; then
    LIBS="$LIBEXODUSII $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libexodusii is enabled, but could not find working libexodusii])
  fi
fi

contrib_lib_enabled=`grep "enable-nemesis.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
AC_MSG_NOTICE([checking whether libMesh contrib package nemesis is enabled... ${contrib_lib_enabled}])
if test "${contrib_lib_enabled}" = yes ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/nemesis/Lib ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([nemesis])
  if test "$HAVE_LIBNEMESIS" = yes ; then
    LIBS="$LIBNEMESIS $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libnemesis is enabled, but could not find working libnemesis])
  fi
fi

contrib_lib_enabled=`grep "enable-libhilbert.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
AC_MSG_NOTICE([checking whether libMesh contrib package libHilbert is enabled... ${contrib_lib_enabled}])
if test "${contrib_lib_enabled}" = yes ; then
  AC_LIB_HAVE_LINKFLAGS([Hilbert])
  if test "$HAVE_LIBHILBERT" = yes ; then
    LIBS="$LIBHILBERT $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libHilbert is enabled, but could not find working libHilbert])
  fi
fi

CPPFLAGS="-I${LIBMESH_DIR}/include ${CPPFLAGS}"
AC_CHECK_HEADER([libmesh/libmesh.h],,AC_MSG_WARN([could not find header file libmesh.h]))
AC_LIB_HAVE_LINKFLAGS([mesh])
if test "$HAVE_LIBMESH" = yes ; then
  LIBS="$LIBMESH $LIBS"
  AC_MSG_CHECKING([for libMesh version 0.8.0])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <libmesh/libmesh_version.h>
]], [[
#if ((LIBMESH_MAJOR_VERSION == 0) && (LIBMESH_MINOR_VERSION == 8) && (LIBMESH_MICRO_VERSION == 0))
#else
asdf
#endif
]])],[LIBMESH_VERSION_0_8_0=yes],[LIBMESH_VERSION_0_8_0=no])
  if test "$LIBMESH_VERSION_0_8_0" = yes; then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_ERROR([incorrect libMesh version detected: please use libMesh 0.8.0-release])
  fi
else
  AC_MSG_WARN([libMesh support is enabled, but could not find working libmesh])
fi

])
