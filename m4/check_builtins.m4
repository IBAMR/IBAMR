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

AC_DEFUN([CHECK_BUILTIN_EXPECT],[
AC_MSG_CHECKING([whether compiler supports __builtin_expect keyword])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[
    if (__builtin_expect(true,1))
    {
        return 0;
    }
    else
    {
        return -1;
    }
]])],[
have_builtin_expect=1
AC_MSG_RESULT(yes)],[
have_builtin_expect=0
AC_MSG_RESULT(no)])
AC_DEFINE_UNQUOTED(HAVE_BUILTIN_EXPECT,$have_builtin_expect,[Boolean value indicating whether the C++ compiler supports the __builtin_expect keyword])
])

# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CHECK_BUILTIN_PREFETCH],[
AC_MSG_CHECKING([whether compiler supports __builtin_prefetch keyword])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[
    __builtin_prefetch(0,0,0);
    return 0;
]])],[
have_builtin_prefetch=1
AC_MSG_RESULT(yes)],[
have_builtin_prefetch=0
AC_MSG_RESULT(no)])
AC_DEFINE_UNQUOTED(HAVE_BUILTIN_PREFETCH,$have_builtin_prefetch,[Boolean value indicating whether the C++ compiler supports the __builtin_prefetch keyword])
])
