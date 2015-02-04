# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_MUPARSER],[
echo
echo "====================================="
echo "Configuring required package muParser"
echo "====================================="
PACKAGE_SETUP_ENVIRONMENT
CONTRIB_SRCDIR=$1
CONTRIB_BUILDDIR=$2
AC_ARG_WITH([muparser],
  AS_HELP_STRING(--with-muparser=PATH,location of required muParser installation),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([it is necessary to specify an existing directory when using --with-muparser=PATH])
   fi
   MUPARSER_DIR=$withval])
if test x$MUPARSER_DIR != x ; then
  if test -d "${MUPARSER_DIR}/include" ; then
    MUPARSER_CPPFLAGS="-I${MUPARSER_DIR}/include"
  fi
  if test -d "${MUPARSER_DIR}/lib" ; then
    MUPARSER_LDFLAGS="-L${MUPARSER_DIR}/lib"
  fi
fi
USING_BUNDLED_MUPARSER=no
echo "first check for a system muParser library; if not found, revert to bundled muParser library"
CPPFLAGS_PREPEND($MUPARSER_CPPFLAGS)
$as_unset ac_cv_header_muParser_h
AC_CHECK_HEADER([muParser.h],[
  HAVE_MUPARSER=yes
  LDFLAGS_PREPEND($MUPARSER_LDFLAGS)
  AC_LIB_HAVE_LINKFLAGS([muparser])
  if test "$HAVE_LIBMUPARSER" = no ; then
    AC_LIB_HAVE_LINKFLAGS([muParser])
    if test "$HAVE_LIBMUPARSER" = no ; then
      AC_MSG_ERROR([could not find working libmuParser])
    fi
  fi
  ],[
  PACKAGE_RESTORE_ENVIRONMENT
  PACKAGE_SETUP_ENVIRONMENT
  USING_BUNDLED_MUPARSER=yes
  AC_MSG_NOTICE([could not find system muParser library, using bundled muParser library])
  MUPARSER_SRCDIR=$CONTRIB_SRCDIR
  MUPARSER_BUILDDIR=$CONTRIB_BUILDDIR
  MUPARSER_CPPFLAGS="-I${MUPARSER_SRCDIR}/include"
  MUPARSER_LDFLAGS=""
  CPPFLAGS_PREPEND($MUPARSER_CPPFLAGS)
  $as_unset ac_cv_header_muParser_h
  AC_CHECK_HEADER([muParser.h],[
    HAVE_MUPARSER=yes
    LIBMUPARSER="$MUPARSER_BUILDDIR/lib/libmuparser.a"
    AC_DEFINE(HAVE_LIBMUPARSER, 1, [Define if you have the libmuParser library.])
    ],[
    AC_MSG_ERROR([could not find bundled muParer library])
  ])
])
AM_CONDITIONAL([USING_BUNDLED_MUPARSER],[test "$USING_BUNDLED_MUPARSER" = yes])
PACKAGE_CPPFLAGS_PREPEND($MUPARSER_CPPFLAGS)
PACKAGE_LDFLAGS_PREPEND($MUPARSER_LDFLAGS)
PACKAGE_CONTRIB_LIBS_PREPEND($LIBMUPARSER)
PACKAGE_RESTORE_ENVIRONMENT
])
