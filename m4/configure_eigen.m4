# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_EIGEN],[
echo
echo "=================================="
echo "Configuring required package Eigen"
echo "=================================="
PACKAGE_SETUP_ENVIRONMENT
CONTRIB_SRCDIR=$1
AC_ARG_WITH([eigen],
  AS_HELP_STRING(--with-eigen=PATH,location of required Eigen installation),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([it is necessary to specify an existing directory when using --with-eigen=PATH])
   else
     EIGEN_DIR=$withval
   fi])
USING_BUNDLED_EIGEN=no
echo "first check for a system Eigen library; if not found, revert to bundled Eigen library"
CHECK_EIGEN($EIGEN_DIR)
if test "$HAVE_EIGEN" = no ; then
  AC_MSG_NOTICE([could not find system Eigen library, using bundled Eigen library])
  EIGEN_DIR=$CONTRIB_SRCDIR
  CHECK_EIGEN($EIGEN_DIR)
  if test "$HAVE_EIGEN" = no ; then
    AC_MSG_ERROR([could not find bundled Eigen library])
  fi
fi
if test x$EIGEN_DIR != x ; then
  if test -d "${EIGEN_DIR}" ; then
    EIGEN_CPPFLAGS="-I${EIGEN_DIR}"
  fi
fi
AM_CONDITIONAL([USING_BUNDLED_EIGEN],[test "$USING_BUNDLED_EIGEN" = yes])
PACKAGE_CPPFLAGS_PREPEND($EIGEN_CPPFLAGS)
PACKAGE_RESTORE_ENVIRONMENT
])

# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CHECK_EIGEN],[
PACKAGE_SETUP_ENVIRONMENT
CHECK_EIGEN_DIR=$1
if test x$CHECK_EIGEN_DIR != x ; then
  if test -d "${CHECK_EIGEN_DIR}" ; then
    CHECK_EIGEN_CPPFLAGS="-I${CHECK_EIGEN_DIR}"
  fi
fi
CPPFLAGS_PREPEND($CHECK_EIGEN_CPPFLAGS)
$as_unset ac_cv_header_Eigen_Eigen
AC_CHECK_HEADER([Eigen/Eigen],HAVE_EIGEN=yes,HAVE_EIGEN=no)
PACKAGE_RESTORE_ENVIRONMENT
])
