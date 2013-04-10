
AC_DEFUN([AC_COMPILERS_ENABLE_64BITS],[

bz_64bit=false
AC_MSG_CHECKING([whether to enable C++/Fortran 64-bit compilation flags])
AC_ARG_ENABLE(64bit,
AS_HELP_STRING([--enable-64bit],[enable C++/Fortran 64-bit compilation flags]),[
if test "$enableval" = yes; then
	AC_MSG_RESULT([yes])
	bz_64bit=true
fi
],[AC_MSG_RESULT([no])])

if test "$bz_64bit" = true; then
	case "$target" in
	*aix*)
		case "$CXX" in
		*xlC*)
			CXXFLAGS="-q64 $CXXFLAGS"
			AR_FLAGS="-cruX64"
			LDFLAGS="-b64 $LDFLAGS"
		;;
		*KCC)
			CXXFLAGS="-q64 $CXXFLAGS"
			AR_FLAGS="-q64 $AR_FLAGS"
			LDFLAGS="-b64 $LDFLAGS"
		;;
		*g++)
			CXXFLAGS="-maix64 $CXXFLAGS"
			AR_FLAGS="-cruX64"
			LDFLAGS="-Wl,-b64 $LDFLAGS"
		;;
		esac

		case "$F77" in
		*xlf*)
			FFLAGS="-q64 $FFLAGS"
		;;
		esac

		case "$FC" in
		*xlf90*)
			FCFLAGS="-q64 $FCFLAGS"
		;;
		esac
	;;

	*irix*)
		case "$CXX" in
		*KCC)
			CXXFLAGS="-64 $CXXFLAGS"
			LDFLAGS="-64 $LDFLAGS"
			AR_FLAGS="-64 $AR_FLAGS"
		;;
		*CC*)
			CXXFLAGS="-64 $CXXFLAGS"
			LDFLAGS="-64 $LDFLAGS"
			AR_FLAGS="-64 $AR_FLAGS"
		;;
		*g++)
			CXXFLAGS="-mabi=64 $CXXFLAGS"
			LDFLAGS="-Wl,-64 $LDFLAGS"
		;;
		esac

		case "$F77" in
		*f77*)
			FFLAGS="-64 $FFLAGS"
		;;
		esac

		case "$FC" in
		*f90*)
			FCFLAGS="-64 $FCFLAGS"
		;;
		esac 
	;;
	esac
fi

])



