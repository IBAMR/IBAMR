# ADD_RPATH_LDFLAG(VALUE)
# ----------------
# Prepend the given library directory with the correct rpath syntax to LDFLAGS.
#
#
AC_DEFUN([ADD_RPATH_LDFLAG],
[
  AC_REQUIRE([AC_LIB_RPATH])
  if test "$enable_rpath" = yes; then
    libdir=$1
    rpath_path=$(eval echo "$acl_cv_hardcode_libdir_flag_spec")
    LDFLAGS_PREPEND("$rpath_path")
  fi
])
