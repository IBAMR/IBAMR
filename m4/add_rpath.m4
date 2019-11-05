## ---------------------------------------------------------------------
##
## Copyright (c) 2019 - 2019 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

# ADD_RPATH_LDFLAG(VALUE)
# ----------------
# Prepend the given library directory with the correct rpath syntax to LDFLAGS.
#
#
AC_DEFUN([ADD_RPATH_LDFLAG],
[
  AC_REQUIRE([AC_LIB_RPATH])
  if test "$enable_rpath" = yes; then
    libdir=$1
    rpath_path=$(eval echo "$acl_cv_hardcode_libdir_flag_spec")
    LDFLAGS_PREPEND("$rpath_path")
  fi
])
