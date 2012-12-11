dnl ax_f90_module_style
dnl 
dnl Figures out the F90 module naming style:
dnl   - whether the base name is lower or upper case, and
dnl   - the suffix (usually .mod or .MOD).
dnl
dnl Y. Hida  (2006-02-17)
dnl
dnl The code is based on the ax_f90_module_extension code:
dnl
dnl @synopsis AX_F90_MODULE_EXTENSION
dnl
dnl Find Fortran 90 modules file extension. The module extension is
dnl stored in the cached variable ax_f90_modext, or "unknown" if the
dnl extension cannot be found.
dnl
dnl @category Fortran
dnl @author Luc Maisonobe <luc@spaceroots.org>
dnl @version 2005-06-17
dnl @license AllPermissive

AC_DEFUN([AX_F90_MODULE_STYLE],[
AC_CACHE_CHECK([fortran 90 modules naming style],
ax_cv_f90_module_style, 
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
  ])],
  [ax_f90_modext=`ls | sed -n 's,conftest_module\.,,p'`
   if test x$ax_f90_modext = x ; then
dnl Some F90 compilers put module filename in uppercase letters
     ax_f90_modext=`ls | sed -n 's,CONFTEST_MODULE\.,,p'`
     if test x$ax_f90_modext = x ; then
       ax_f90_modext=unknown
     else
       ax_f90_module_style="uppercase, $ax_f90_modext"
     fi
   else
     ax_f90_module_style="lowercase, $ax_f90_modext"
   fi
  ],
  [ax_f90_module_style=unknown])
cd ..
rm -fr tmpdir_$i
AC_LANG_POP(Fortran)
])])
