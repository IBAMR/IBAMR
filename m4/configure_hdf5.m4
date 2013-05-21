# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_HDF5],[
echo
echo "================================="
echo "Configuring required package HDF5"
echo "================================="

PACKAGE_SETUP_ENVIRONMENT

if test `grep -c HDF5 "${PETSC_DIR}/${PETSC_ARCH}/conf/petscvariables"` != 0 ; then
  AC_MSG_NOTICE([PETSc appears to provide hdf5; using PETSc hdf5 library])
  PETSC_BUNDLES_HDF5=yes
else
  PETSC_BUNDLES_HDF5=no

  AC_ARG_WITH([hdf5],
    AS_HELP_STRING(--with-hdf5=PATH,location of required hdf5 installation),
    [if test ! -d "$withval" ; then
       AC_MSG_ERROR([it is necessary to specify an existing directory when using --with-hdf5=PATH])
     fi
     HDF5_DIR=$withval])

  if test x$HDF5_DIR != x ; then
    if test -d "${HDF5_DIR}/include" ; then
      HDF5_CPPFLAGS="-I${HDF5_DIR}/include"
    fi
    if test -d "${HDF5_DIR}/lib" ; then
      HDF5_LDFLAGS="-L${HDF5_DIR}/lib"
    fi
  fi

  CPPFLAGS_PREPEND($HDF5_CPPFLAGS)
  AC_CHECK_HEADER([hdf5.h],,AC_MSG_ERROR([could not find header file hdf5.h]))

  LDFLAGS_PREPEND($HDF5_LDFLAGS)
  AC_LIB_HAVE_LINKFLAGS([hdf5],[z,sz])
  if test "$HAVE_LIBHDF5" = no ; then
    AC_MSG_ERROR([could not find working libhdf5])
  fi

  LIBS_PREPEND($LIBHDF5)
  AC_LIB_HAVE_LINKFLAGS([hdf5_hl],[z,sz])

  PACKAGE_CPPFLAGS_PREPEND($HDF5_CPPFLAGS)
  PACKAGE_LDFLAGS_PREPEND($HDF5_LDFLAGS)
  PACKAGE_LIBS_PREPEND("$LIBHDF5")
  if test "$HAVE_LIBHDF5_HL" = yes ; then
    PACKAGE_LIBS_PREPEND("$LIBHDF5_HL")
  fi
fi

PACKAGE_RESTORE_ENVIRONMENT

])
