AC_DEFUN([CONFIGURE_LIBMESH],[

LIBMESH_HOSTTYPE=`grep "hosttype[ \t]*.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
LIBMESH_LIBDIR=${LIBMESH_DIR}/lib/${LIBMESH_HOSTTYPE}_${LIBMESH_COMPILE_MODE}
LIBMESH_CONTRIBDIR=${LIBMESH_DIR}/contrib
LIBMESH_CONTRIB_LIBDIR=${LIBMESH_DIR}/contrib/lib/${LIBMESH_HOSTTYPE}_${LIBMESH_COMPILE_MODE}

if test "${LIBMESH_COMPILE_MODE}" == "opt" ; then
  CXXFLAGS="${CXXFLAGS} `grep -A1 opt-mode $LIBMESH_DIR/Make.common | grep libmesh_CXXFLAGS | sed -e 's/.*libmesh_CXXFLAGS +=//'`"
fi
if test "${LIBMESH_COMPILE_MODE}" == "devel" ; then
  CXXFLAGS="${CXXFLAGS} `grep -A1 devel-mode $LIBMESH_DIR/Make.common | grep libmesh_CXXFLAGS | sed -e 's/.*libmesh_CXXFLAGS +=//'`"
fi
if test "${LIBMESH_COMPILE_MODE}" == "dbg" ; then
  CXXFLAGS="${CXXFLAGS} `grep -A1 debug-mode $LIBMESH_DIR/Make.common | grep libmesh_CXXFLAGS | sed -e 's/.*libmesh_CXXFLAGS +=//'`"
fi
if test "${LIBMESH_COMPILE_MODE}" == "prof" ; then
  CXXFLAGS="${CXXFLAGS} `grep -A1 prof-mode $LIBMESH_DIR/Make.common | grep libmesh_CXXFLAGS | sed -e 's/.*libmesh_CXXFLAGS +=//'`"
fi
if test "${LIBMESH_COMPILE_MODE}" == "oprof" ; then
  CXXFLAGS="${CXXFLAGS} `grep -A1 oprof-mode $LIBMESH_DIR/Make.common | grep libmesh_CXXFLAGS | sed -e 's/.*libmesh_CXXFLAGS +=//'`"
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

for subdir in base enums error_estimation fe geom mesh numerics parallel partitioning quadrature solvers systems utils; do
  if test -d "${LIBMESH_DIR}/include/${subdir}" ; then
    CPPFLAGS="-I${LIBMESH_DIR}/include/${subdir} ${CPPFLAGS}"
  else
    AC_MSG_WARN([libMesh support is enabled, but could not find expected libMesh include directory ${LIBMESH_DIR}/include/${subdir}])
  fi
done

contrib_lib_enabled=`grep "enable-laspack.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
echo "checking whether libMesh contrib package laspack is enabled... ${contrib_lib_enabled}"
if test "${contrib_lib_enabled}" == "yes" ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/laspack ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([laspack])
  if test "$HAVE_LIBLASPACK" == "yes" ; then
    LIBS="$LIBLASPACK $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib liblaspack is enabled, but could not find working liblaspack])
  fi
fi

contrib_lib_enabled=`grep "enable-metis.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
echo "checking whether libMesh contrib package metis is enabled... ${contrib_lib_enabled}"
if test "${contrib_lib_enabled}" == "yes" ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/metis/Lib ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([metis])
  if test "$HAVE_LIBMETIS" == "yes" ; then
    LIBS="$LIBMETIS $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libmetis is enabled, but could not find working libmetis])
  fi
fi

contrib_lib_enabled=`grep "enable-parmetis.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
echo "checking whether libMesh contrib package parmetis is enabled... ${contrib_lib_enabled}"
if test "${contrib_lib_enabled}" == "yes" ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/parmetis/Lib ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([parmetis])
  if test "$HAVE_LIBPARMETIS" == "yes" ; then
    LIBS="$LIBPARMETIS $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libparmetis is enabled, but could not find working libparmetis])
  fi
fi

contrib_lib_enabled=`grep "enable-sfcurves.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
echo "checking whether libMesh contrib package sfcurves is enabled... ${contrib_lib_enabled}"
if test "${contrib_lib_enabled}" == "yes" ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/sfcurves ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([sfcurves])
  if test "$HAVE_LIBSFCURVES" == "yes" ; then
    LIBS="$LIBSFCURVES $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libsfcurves is enabled, but could not find working libsfcurves])
  fi
fi

contrib_lib_enabled=`grep "enable-gzstream.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
echo "checking whether libMesh contrib package gzstream is enabled... ${contrib_lib_enabled}"
if test "${contrib_lib_enabled}" == "yes" ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/gzstream ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([gzstream])
  if test "$HAVE_LIBGZSTREAM" == "yes" ; then
    LIBS="$LIBGZSTREAM $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libgzstream is enabled, but could not find working libgzstream])
  fi
fi

contrib_lib_enabled=`grep "enable-tetgen.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
echo "checking whether libMesh contrib package tetgen is enabled... ${contrib_lib_enabled}"
if test "${contrib_lib_enabled}" == "yes" ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/tetgen ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([tetgen])
  if test "$HAVE_LIBTETGEN" == "yes" ; then
    LIBS="$LIBTETGEN $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libtetgen is enabled, but could not find working libtetgen])
  fi
fi

contrib_lib_enabled=`grep "enable-triangle.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
echo "checking whether libMesh contrib package triangle is enabled... ${contrib_lib_enabled}"
if test "${contrib_lib_enabled}" == "yes" ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/triangle ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([triangle])
  if test "$HAVE_LIBTRIANGLE" == "yes" ; then
    LIBS="$LIBTRIANGLE $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libtriangle is enabled, but could not find working libtriangle])
  fi
fi

contrib_lib_enabled=`grep "enable-gmv.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
echo "checking whether libMesh contrib package gmv is enabled... ${contrib_lib_enabled}"
if test "${contrib_lib_enabled}" == "yes" ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/gmv ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([gmv])
  if test "$HAVE_LIBGMV" == "yes" ; then
    LIBS="$LIBGMV $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libgmv is enabled, but could not find working libgmv])
  fi
fi

contrib_lib_enabled=`grep "enable-vtk.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
echo "checking whether libMesh contrib package vtk is enabled... ${contrib_lib_enabled}"
if test "${contrib_lib_enabled}" == "yes" ; then
  AC_LIB_HAVE_LINKFLAGS([vtk])
  if test "$HAVE_LIBVTK" == "yes" ; then
    LIBS="$LIBVTK $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libvtk is enabled, but could not find working libvtk])
  fi
fi

contrib_lib_enabled=`grep "enable-netcdf.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
echo "checking whether libMesh contrib package netcdf is enabled... ${contrib_lib_enabled}"
if test "${contrib_lib_enabled}" == "yes" ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/netcdf/Lib ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([netcdf])
  if test "$HAVE_LIBNETCDF" == "yes" ; then
    LIBS="$LIBNETCDF $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libnetcdf is enabled, but could not find working libnetcdf])
  fi
fi

contrib_lib_enabled=`grep "enable-exodus.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
echo "checking whether libMesh contrib package exodus is enabled... ${contrib_lib_enabled}"
if test "${contrib_lib_enabled}" == "yes" ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/exodusii/Lib/include ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([exodusii])
  if test "$HAVE_LIBEXODUSII" == "yes" ; then
    LIBS="$LIBEXODUSII $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libexodusii is enabled, but could not find working libexodusii])
  fi
fi

contrib_lib_enabled=`grep "enable-nemesis.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
echo "checking whether libMesh contrib package nemesis is enabled... ${contrib_lib_enabled}"
if test "${contrib_lib_enabled}" == "yes" ; then
  CPPFLAGS="-I${LIBMESH_CONTRIBDIR}/nemesis/Lib ${CPPFLAGS}"
  AC_LIB_HAVE_LINKFLAGS([nemesis])
  if test "$HAVE_LIBNEMESIS" == "yes" ; then
    LIBS="$LIBNEMESIS $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libnemesis is enabled, but could not find working libnemesis])
  fi
fi

contrib_lib_enabled=`grep "enable-libhilbert.*=" $LIBMESH_DIR/Make.common | sed -e 's/.*=//' | sed -e 's/[ ]*//' | sed -e 's/[\t]*//'`
echo "checking whether libMesh contrib package libhilbert is enabled... ${contrib_lib_enabled}"
if test "${contrib_lib_enabled}" == "yes" ; then
  AC_LIB_HAVE_LINKFLAGS([hilbert])
  if test "$HAVE_LIBHILBERT" == "yes" ; then
    LIBS="$LIBHILBERT $LIBS"
  else
    AC_MSG_WARN([libMesh contributed lib libhilbert is enabled, but could not find working libhilbert])
  fi
fi

AC_CHECK_HEADER([libmesh.h],,AC_MSG_WARN([could not find header file libmesh.h]))
AC_LIB_HAVE_LINKFLAGS([mesh])
if test "$HAVE_LIBMESH" == "yes" ; then
  LIBS="$LIBMESH $LIBS"
else
  AC_MSG_WARN([libMesh support is enabled, but could not find working libmesh])
fi

])