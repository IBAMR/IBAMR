AC_DEFUN([CONFIGURE_HDF5],[

if test -d "${HDF5_DIR}/lib" ; then
  LDFLAGS="-L${HDF5_DIR}/lib $LDFLAGS"
fi
if test -d "${HDF5_DIR}/include" ; then
  CPPFLAGS="-I${HDF5_DIR}/include $CPPFLAGS"
fi

AC_CHECK_HEADER([hdf5.h],,AC_MSG_ERROR([could not find header file hdf5.h]))

AC_LIB_HAVE_LINKFLAGS([hdf5],[z,sz])
LIBS="$LIBHDF5 $LIBS"

AC_LIB_HAVE_LINKFLAGS([hdf5_hl])
LIBS="$LIBHDF5_HL $LIBS"

])
