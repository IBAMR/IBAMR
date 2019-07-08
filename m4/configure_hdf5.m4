# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_HDF5],[
echo
echo "================================="
echo "Configuring required package HDF5"
echo "================================="

if test `grep -c HDF5 "${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscvariables"` != 0 ; then
  AC_MSG_NOTICE([PETSc appears to provide HDF5; using PETSc HDF5 library])
  PETSC_BUNDLES_HDF5=yes
else
  PETSC_BUNDLES_HDF5=no

  AC_ARG_WITH([hdf5],
    AS_HELP_STRING(--with-hdf5=PATH,location of required HDF5 installation),
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
  AC_CHECK_LIB([hdf5], H5open, [],
               [AC_MSG_ERROR([could not find working libhdf5])])
  AC_CHECK_LIB([hdf5_hl], H5LTfind_dataset, [],
               [AC_MSG_ERROR([could not find working libhdf5_hl])])
  # set up rpath
  ADD_RPATH_LDFLAG($HDF5_LDFLAGS)
fi

AC_MSG_CHECKING([for HDF5 version >= 1.8.7])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <hdf5.h>
]], [[
#if H5_VERSION_GE(1,8,7)
#else
asdf
#endif
]])],[HDF5_VERSION_VALID=yes],[HDF5_VERSION_VALID=no])
AC_MSG_RESULT([${HDF5_VERSION_VALID}])
if test "$HDF5_VERSION_VALID" = no; then
  AC_MSG_WARN([HDF5 versions prior to 1.8.7 are likely to be usable but are not officially supported])
fi

])
