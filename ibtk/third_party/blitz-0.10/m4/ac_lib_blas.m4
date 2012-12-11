
AC_DEFUN([AC_LIB_BLAS],[
dnl See if there is a libblas.a which can be used for benchmark comparisons.
dnl blas = Basic Linear Algebra Subroutines

AC_SUBST(BLASLIB)

AC_ARG_WITH(blas,
AS_HELP_STRING([--with-blas@<:@=DIR@:>@],[set the path for the blas library]),
[BLASPATH=-L"$withval"], [BLASPATH=""])

AC_MSG_NOTICE([

Checking to see if there is a usable BLAS library kicking around...
(If so, Blitz++ will use BLAS library to run benchmarks against
if you build the benchmarks.  Blitz++ itself does NOT use BLAS).

])

LDFLAGSsave=$LDFLAGS
LDFLAGS="$BLASPATH $LDFLAGS"

AC_CHECK_LIB(blas, daxpy_,
   [ CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -DUSE_LIBBLAS"
     CXXFFLAGS="-DBZ_FORTRAN_SYMBOLS_WITH_TRAILING_UNDERSCORES"
		 if test -z "$BLASPATH" ; then
       BLASLIB="-lblas"
     else
       BLASLIB="$BLASPATH -lblas"
     fi],[],[$FLIBS])

AC_CHECK_LIB(blas, daxpy,
   [ CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -DUSE_LIBBLAS"
     CXXFFLAGS=""
		 if test -z "$BLASPATH" ; then
       BLASLIB="-lblas"
     else
       BLASLIB="$BLASPATH -lblas"
     fi],[],[$FLIBS])

AC_CHECK_LIB(blas, DAXPY,
   [ CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -DUSE_LIBBLAS"
     CXXFFLAGS="-DBZ_FORTRAN_SYMBOLS_CAPS"
		 if test -z "$BLASPATH" ; then
       BLASLIB="-lblas"
     else
       BLASLIB="$BLASPATH -lblas"
     fi],[],[$FLIBS])

AC_CHECK_LIB(blas, daxpy__,
   [ CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -DUSE_LIBBLAS"
     CXXFFLAGS="-DBZ_FORTRAN_SYMBOLS_WITH_DOUBLE_TRAILING_UNDERSCORES"
		 if test -z "$BLASPATH" ; then
       BLASLIB="-lblas"
     else
       BLASLIB="$BLASPATH -lblas"
     fi],[],[$FLIBS])

LDFLAGS=$LDFLAGSsave

])
