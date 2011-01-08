AC_DEFUN([AX_FC_ETIME], [
AC_MSG_CHECKING([for etime])
AC_LANG_PUSH(Fortran)
ax_fc_etime=
ax_fc_etime_names="etime etime_"
for name in $ax_fc_etime_names; do
  AC_LINK_IFELSE([AC_LANG_PROGRAM(, [[
      real*4 t(2), tot
      tot = $name(t)]])], 
    [ax_fc_etime=$name], [])
  if test "x$ax_fc_etime" != "x"; then
    break;
  fi
done
AC_LANG_POP(Fortran)
if test "x$ax_fc_etime" != "x"; then
  AC_MSG_RESULT($ax_fc_etime)
  $1
else
  AC_MSG_RESULT(none)
  ifelse([$2],,AC_MSG_ERROR([Cannot find etime.]), [$2])
fi
])
