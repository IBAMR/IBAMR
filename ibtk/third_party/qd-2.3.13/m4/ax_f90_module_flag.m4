dnl @synopsis AX_F90_MODULE_FLAG
dnl
dnl @summary Find Fortran 90 modules inclusion flag.
dnl
dnl Find Fortran 90 modules inclusion flag. The module inclusion flag
dnl is stored in the cached variable ax_f90_modflag. An error is
dnl triggered if the flag cannot be found. Supported are the -I GNU
dnl compilers flag, the -M SUN compilers flag, and the -p Absoft Pro
dnl Fortran compiler flag.
dnl
dnl @category Fortran
dnl @author Luc Maisonobe <luc@spaceroots.org>
dnl @author Julian C. Cummings <cummings@cacr.caltech.edu>
dnl @version 2006-01-28
dnl @license AllPermissive

AC_DEFUN([AX_F90_MODULE_FLAG],[
AC_CACHE_CHECK([fortran 90 modules inclusion flag],
ax_cv_f90_modflag,
[AC_LANG_PUSH(Fortran)
i=0
while test \( -f tmpdir_$i \) -o \( -d tmpdir_$i \) ; do
  i=`expr $i + 1`
done
mkdir tmpdir_$i
cd tmpdir_$i
AC_COMPILE_IFELSE([AC_LANG_SOURCE([module conftest_module
   contains
   subroutine conftest_routine
   write(*,'(a)') 'gotcha!'
   end subroutine conftest_routine
   end module conftest_module
  ])],[],[])
cd ..
ax_f90_modflag="not found"
for ax_flag in "-I " "-M" "-p"; do
  if test "$ax_f90_modflag" = "not found" ; then
    ax_save_FCFLAGS="$FCFLAGS"
    FCFLAGS="$ax_save_FCFLAGS ${ax_flag}tmpdir_$i"
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([program conftest_program
       use conftest_module
       call conftest_routine
       end program conftest_program
      ])],[ax_f90_modflag="$ax_flag"],[])
    FCFLAGS="$ax_save_FCFLAGS"
  fi
done
rm -fr tmpdir_$i
if test "$ax_modflag" = "not found" ; then
  AC_MSG_ERROR([unable to find compiler flag for modules inclusion])
fi
AC_LANG_POP(Fortran)
])])
