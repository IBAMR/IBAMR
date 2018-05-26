# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_BOOST],[
#echo
#echo "=================================="
#echo "Configuring required package Boost"
#echo "=================================="
PACKAGE_SETUP_ENVIRONMENT
CONTRIB_SRCDIR=$1
AC_ARG_VAR(BOOST_ROOT,[the location of the Boost installation that is to be used.])
#AC_MSG_NOTICE([first check for a system Boost library; if not found, revert to bundled Boost library])
#$as_unset boost_cv_inc_path
#$as_unset boost_cv_lib_version
#$as_unset boost_cv_version
#BOOST_REQUIRE([1.57.0],[AC_MSG_WARN([could not find system Boost library, using bundled Boost library])])
if test x"$USING_BUNDLED_BOOST" = xyes ; then
  if test "$LIBMESH_ENABLED" = yes ; then
    AC_MSG_ERROR([IBAMR is configured to use a bundled Boost library, but when IBAMR is also configured to use libMesh, IBAMR and libMesh both must use the same (external) Boost library])
  fi
  AC_MSG_NOTICE([configuring bundled Boost library])
  BOOST_CPPFLAGS="-I$CONTRIB_SRCDIR"
  AC_DEFINE([HAVE_BOOST], [1], [Defined if the requested minimum BOOST version is satisfied])
else
  USING_BUNDLED_BOOST=no
fi
CPPFLAGS_PREPEND($BOOST_CPPFLAGS)
AC_CHECK_HEADER([boost/array.hpp],,AC_MSG_ERROR([cannot find working boost/array.hpp]))
AC_CHECK_HEADER([boost/multi_array.hpp],,AC_MSG_ERROR([cannot find working boost/multi_array.hpp]))
PACKAGE_CPPFLAGS_PREPEND("$BOOST_CPPFLAGS")
PACKAGE_RESTORE_ENVIRONMENT
])
