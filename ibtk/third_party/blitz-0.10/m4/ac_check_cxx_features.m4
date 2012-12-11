
AC_DEFUN([AC_CHECK_CXX_FEATURES],[

AC_MSG_NOTICE([

C++ compiler ($CXX $CXXFLAGS $LDFLAGS) characteristics

])

OS=`uname -a`
AC_SUBST(OS)
DATE=`date`
AC_SUBST(DATE)

AH_TOP([
/******************************************************************************
 * config.h           Compiler language support flags
 *
 * This file was generated automatically when running the configure script.
 * You should rerun configure each time you switch compilers, install new
 * standard libraries, or change compiler versions.
 *
 */

])

AC_DEFINE_UNQUOTED([_compiler_name],["$CXX"],[CXX])
AC_DEFINE_UNQUOTED([_compiler_options],["$CXXFLAGS"],[CXXFLAGS])
AC_DEFINE_UNQUOTED([_os_name],["$OS"],[uname -a])
AC_DEFINE_UNQUOTED([_config_date],["$DATE"],[date])
AC_DEFINE_UNQUOTED([_platform],["$target"],[target])

AC_CXX_GENERAL
AC_CXX_KEYWORDS
AC_CXX_TYPE_CASTS
AC_CXX_TEMPLATES_FEATURES
AC_CXX_STANDARD_LIBRARY

]))
