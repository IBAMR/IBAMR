# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_BOOST],[
echo
echo "=================================="
echo "Configuring required package Boost"
echo "=================================="
PACKAGE_SETUP_ENVIRONMENT
CONTRIB_SRCDIR=$1
AC_ARG_VAR(BOOST_ROOT,[the location of the Boost installation that is to be used.])
USING_BUNDLED_BOOST=no
echo "first check for a system Boost library; if not found, revert to bundled Boost library"
$as_unset boost_cv_inc_path
BOOST_REQUIRE([1.57.0],[
  USING_BUNDLED_BOOST=yes
  AC_MSG_NOTICE([could not find system Boost library, using bundled Boost library])
  $as_unset boost_cv_inc_path
  with_boost="$CONTRIB_SRCDIR"
  BOOST_REQUIRE([1.57.0],[
    AC_MSG_ERROR([could not find bundled Boost library])
  ])
])
BOOST_ARRAY
BOOST_MULTIARRAY
AM_CONDITIONAL([USING_BUNDLED_BOOST],[test "$USING_BUNDLED_BOOST" = yes])
PACKAGE_CPPFLAGS_PREPEND("$BOOST_CPPFLAGS")
PACKAGE_RESTORE_ENVIRONMENT
])
