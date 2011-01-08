dnl Still not quite there yet, as FC_DUMMY_MAIN may be needed
dnl but that complicates matters since FC_DUMMY_MAIN needs the 
dnl correct libraries to link.
dnl
AC_DEFUN([AX_CXX_FC_LIB], [
AC_MSG_CHECKING([whether extra library is needed to link Fortran and C++])
AC_LANG_PUSH(C++)
ax_cxx_fc_lib=
ax_cxx_fc_lib_names="none -lompstubs -lmtsk"
for name in $ax_cxx_fc_lib_names; do
  if test "x$name" = xnone; then
    lib=""
  else
    lib=$name
  fi
  save_LIBS="$LIBS"
  LIBS="$LIBS $FCLIBS $lib"
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], [])], [ax_cxx_fc_lib=$name], [])
  LIBS="$save_LIBS"
  if test "x$ax_cxx_fc_lib" != "x"; then
    break
  fi
done
AC_LANG_POP(C++)
if test "x$ax_cxx_fc_lib" != x; then
  AC_MSG_RESULT($ax_cxx_fc_lib)
  if test "x$ax_cxx_fc_lib" = xnone; then
    ax_cxx_fc_lib=
  fi
  ifelse([$1],,FCLIBS="$FCLIBS $ax_cxx_fc_lib", [$2])
else
  AC_MSG_RESULT(unknown)
  ifelse([$2],,AC_MSG_ERROR(Cannot link C++ and Fortran.), [$2])
fi
])
