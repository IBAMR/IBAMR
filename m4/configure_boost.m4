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

dnl There is some additional boost-checking logic in the libMesh configuration
dnl file since boost is an optional dependency of libMesh.

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
  AC_MSG_NOTICE([configuring bundled Boost library])
  BOOST_CPPFLAGS="-I$CONTRIB_SRCDIR"
  AC_DEFINE([HAVE_BOOST], [1], [Defined if the requested minimum BOOST version is satisfied])
else
  USING_BUNDLED_BOOST=no
fi
CPPFLAGS_PREPEND($BOOST_CPPFLAGS)
# boost will try to find nonexistent headers if we don't specify the correct
# language version
CPPFLAGS_PREPEND($CXX_VERSION_FLAG)
AC_CHECK_HEADER([boost/multi_array.hpp],,AC_MSG_ERROR([cannot find working boost/multi_array.hpp]))
AC_CHECK_HEADER([boost/math/special_functions/round.hpp],,AC_MSG_ERROR([cannot find working boost/math/special_functions/round.hpp]))
AC_CHECK_HEADER([boost/math/tools/roots.hpp],,AC_MSG_ERROR([cannot find working boost/math/tools/roots.hpp]))
AC_CHECK_HEADER([boost/cstdint.hpp],,AC_MSG_ERROR([cannot find working boost/cstdint.h]))
PACKAGE_CPPFLAGS_PREPEND("$BOOST_CPPFLAGS")
PACKAGE_RESTORE_ENVIRONMENT
])
