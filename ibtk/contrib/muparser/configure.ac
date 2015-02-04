###########################################################################
# Prelude.
###########################################################################
AC_INIT([muParser],[2.2.3],[muparser@beltoforion.de],[muParser],[http://muparser.beltoforion.de])
AC_PREREQ([2.69])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER([config/muparser_config.h.tmp])
AX_PREFIX_CONFIG_H([config/muparser_config.h],[muParser],[config/muparser_config.h.tmp])
AM_INIT_AUTOMAKE([1.12 -Wall -Werror dist-bzip2 foreign -Wno-extra-portability subdir-objects])
AM_SILENT_RULES([yes])
AM_MAINTAINER_MODE([disable])
AC_LANG(C++)

###########################################################################
# Checks for programs.
###########################################################################
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_RANLIB
AM_PROG_AR

###########################################################################
# Output files.
###########################################################################
AC_CONFIG_FILES([Makefile lib/Makefile])
AC_OUTPUT
