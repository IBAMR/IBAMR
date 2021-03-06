## ---------------------------------------------------------------------
##
## Copyright (c) 2014 - 2019 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

AC_DEFUN([CONFIGURE_DOXYGEN],[
DOXYGEN_PATH=$PATH
AC_ARG_WITH([doxygen],
  AS_HELP_STRING([--with-doxygen=PATH],[manually set directory where the doxygen executable resides to PATH]),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([it is necessary to specify an existing directory when using --with-doxygen=PATH])
   fi
   DOXYGEN_DIR=$withval
   DOXYGEN_PATH=$DOXYGEN_DIR$PATH_SEPARATOR$DOXYGEN_PATH])
AC_PATH_PROG(DOXYGEN, doxygen, [], $DOXYGEN_PATH)
if test x$DOXYGEN != x ; then
  HAVE_DOXYGEN=yes
else
  AC_MSG_WARN([doxygen not found])
  HAVE_DOXYGEN=no
  echo "if doxygen is installed, specify its location via --with-doxygen=PATH"
fi
AC_SUBST(HAVE_DOXYGEN,$HAVE_DOXYGEN)
AC_SUBST(DOXYGEN,$DOXYGEN)
AC_SUBST(DOXYGEN_DIR,$DOXYGEN_DIR)
])

# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_DOT],[
DOT_PATH=$PATH
AC_ARG_WITH([dot],
  AS_HELP_STRING([--with-dot=PATH],[manually set directory where the Graphviz dot executable resides to PATH]),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([it is necessary to specify an existing directory when using --with-dot=PATH])
   fi
   DOT_DIR=$withval
   DOT_PATH=$DOT_DIR$PATH_SEPARATOR$DOT_PATH])
AC_PATH_PROG(DOT, dot, [], $DOT_PATH)
if test x$DOT != x ; then
  HAVE_DOT=yes
else
  AC_MSG_WARN([dot not found])
  HAVE_DOT=no
  echo "if dot is installed, specify its location via --with-dot=PATH"
fi
AC_SUBST(HAVE_DOT,$HAVE_DOT)
AC_SUBST(DOT,$DOT)
AC_SUBST(DOT_DIR,$DOT_DIR)
])
