
AC_DEFUN([AC_SET_COMPILER_SPECIFIC_HEADER],[

AC_ARG_VAR([COMPILER_SPECIFIC_HEADER],[
automatically set if supported compiler, otherwise
should be set to <some-compiler-specific-path>/bzconfig.h.
See also file blitz/bzconfig.h
])

AC_MSG_CHECKING(compiler specific header)

case "$CXX" in
  *xlc++*)      AX_PREFIX_CONFIG_H([blitz/apple/bzconfig.h],[BZ])
                COMPILER_SPECIFIC_HEADER="apple/bzconfig.h" ;;
  *icpc*|*icc*) AX_PREFIX_CONFIG_H([blitz/intel/bzconfig.h],[BZ])
                COMPILER_SPECIFIC_HEADER="intel/bzconfig.h" ;;
  *xlC*)        AX_PREFIX_CONFIG_H([blitz/ibm/bzconfig.h],[BZ])
                COMPILER_SPECIFIC_HEADER="ibm/bzconfig.h" ;;
  *cxx*)        AX_PREFIX_CONFIG_H([blitz/compaq/bzconfig.h],[BZ])
                COMPILER_SPECIFIC_HEADER="compaq/bzconfig.h" ;;
  *aCC*)        AX_PREFIX_CONFIG_H([blitz/hp/bzconfig.h],[BZ])
                COMPILER_SPECIFIC_HEADER="hp/bzconfig.h" ;;
  *g++*|*c++*)  AX_PREFIX_CONFIG_H([blitz/gnu/bzconfig.h],[BZ])
                COMPILER_SPECIFIC_HEADER="gnu/bzconfig.h" ;;
  *KCC*)        AX_PREFIX_CONFIG_H([blitz/kai/bzconfig.h],[BZ])
                COMPILER_SPECIFIC_HEADER="kai/bzconfig.h" ;;
  *pgCC*)       AX_PREFIX_CONFIG_H([blitz/pgi/bzconfig.h],[BZ])
                COMPILER_SPECIFIC_HEADER="pgi/bzconfig.h" ;;
dnl  *FCC*)        AX_PREFIX_CONFIG_H([blitz/fujitsu/bzconfig.h],[BZ])
dnl                COMPILER_SPECIFIC_HEADER="fujitsu/bzconfig.h" ;;
  *pathCC*)     AX_PREFIX_CONFIG_H([blitz/pathscale/bzconfig.h],[BZ])
                COMPILER_SPECIFIC_HEADER="pathscale/bzconfig.h" ;;
  *CC*)
    case "$target" in
    *sgi*)      AX_PREFIX_CONFIG_H([blitz/sgi/bzconfig.h],[BZ])
                COMPILER_SPECIFIC_HEADER="sgi/bzconfig.h" ;;
    *solaris*)  AX_PREFIX_CONFIG_H([blitz/sun/bzconfig.h],[BZ])
                COMPILER_SPECIFIC_HEADER="sun/bzconfig.h" ;;
    *cray*)     AX_PREFIX_CONFIG_H([blitz/cray/bzconfig.h],[BZ])
                COMPILER_SPECIFIC_HEADER="cray/bzconfig.h" ;;
    *fujitsu*)  AX_PREFIX_CONFIG_H([blitz/fujitsu/bzconfig.h],[BZ])
                COMPILER_SPECIFIC_HEADER="fujitsu/bzconfig.h" ;;
    esac
  ;;
esac

AC_MSG_RESULT($COMPILER_SPECIFIC_HEADER)

AC_SUBST(COMPILER_SPECIFIC_HEADER,$COMPILER_SPECIFIC_HEADER)


])
