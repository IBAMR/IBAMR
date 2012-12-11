dnl @* AX_CREATE_PKGCONFIG_INFO [(outputfile, [requires [,libs [,summary]]])]
dnl defaults:
dnl   $1 = $PACKAGE_NAME.pc
dnl   $2 = (empty)
dnl   $3 = $PACKAGE_LIBS $LIBS (as set at that point in configure.ac)
dnl   $4 = $PACKAGE_SUMMARY (or $1 Library)
dnl   $5 = $CPPFLAGS $PACKAGE_CFLAGS (as set at the point in configure.ac)
dnl
dnl   PACKAGE_NAME defaults to $PACKAGE if not set.
dnl   PACKAGE_LIBS defaults to -l$PACKAGE_NAME if not set.
dnl
dnl the resulting file is called $PACKAGE.pc.in / $PACKAGE.pc
dnl
dnl You will find this macro most useful in conjunction with ax_spec_defaults
dnl that can read good initializers from the .spec file. In consequencd, most
dnl of the generatable installable stuff can be made from information being
dnl updated in a single place for the whole project.
dnl
dnl @$Id$

AC_DEFUN([AX_CREATE_PKGCONFIG_INFO],[dnl
AS_VAR_PUSHDEF([PKGCONFIG_libdir],[ax_create_pkgconfig_libdir])dnl
AS_VAR_PUSHDEF([PKGCONFIG_libfile],[ax_create_pkgconfig_libfile])dnl
AS_VAR_PUSHDEF([PKGCONFIG_libname],[ax_create_pkgconfig_libname])dnl
AS_VAR_PUSHDEF([PKGCONFIG_version],[ax_create_pkgconfig_version])dnl
AS_VAR_PUSHDEF([PKGCONFIG_description],[ax_create_pkgconfig_description])dnl
AS_VAR_PUSHDEF([PKGCONFIG_requires],[ax_create_pkgconfig_requires])dnl
AS_VAR_PUSHDEF([PKGCONFIG_pkglibs],[ax_create_pkgconfig_pkglibs])dnl
AS_VAR_PUSHDEF([PKGCONFIG_libs],[ax_create_pkgconfig_libs])dnl
AS_VAR_PUSHDEF([PKGCONFIG_cppflags],[ax_create_pkgconfig_cppflags])dnl
AS_VAR_PUSHDEF([PKGCONFIG_generate],[ax_create_pkgconfig_generate])dnl
AS_VAR_PUSHDEF([PKGCONFIG_src_libdir],[ax_create_pkgconfig_src_libdir])dnl
AS_VAR_PUSHDEF([PKGCONFIG_src_headers],[ax_create_pkgconfig_src_headers])dnl

# we need the expanded forms...
test "x$prefix" = xNONE && prefix=$ac_default_prefix
test "x$exec_prefix" = xNONE && exec_prefix='${prefix}'

AC_MSG_CHECKING(our pkgconfig libname)
test ".$PKGCONFIG_libname" != "." || \
PKGCONFIG_libname="ifelse($1,,${PACKAGE_NAME},`basename $1 .pc`)"
test ".$PKGCONFIG_libname" != "." || \
PKGCONFIG_libname="$PACKAGE"
PKGCONFIG_libname=`eval echo "$PKGCONFIG_libname"`
PKGCONFIG_libname=`eval echo "$PKGCONFIG_libname"`
AC_MSG_RESULT($PKGCONFIG_libname)

AC_MSG_CHECKING(our pkgconfig version)
test ".$PKGCONFIG_version" != "." || \
PKGCONFIG_version="${PACKAGE_VERSION}"
test ".$PKGCONFIG_version" != "." || \
PKGCONFIG_version="$VERSION"
PKGCONFIG_version=`eval echo "$PKGCONFIG_version"`
PKGCONFIG_version=`eval echo "$PKGCONFIG_version"`
AC_MSG_RESULT($PKGCONFIG_version)

AC_MSG_CHECKING(our pkgconfig_libdir)
test ".$pkgconfig_libdir" = "." && \
pkgconfig_libdir='${libdir}/pkgconfig'
PKGCONFIG_libdir=`eval echo "$pkgconfig_libdir"`
PKGCONFIG_libdir=`eval echo "$PKGCONFIG_libdir"`
PKGCONFIG_libdir=`eval echo "$PKGCONFIG_libdir"`
AC_MSG_RESULT($pkgconfig_libdir)
test "$pkgconfig_libdir" != "$PKGCONFIG_libdir" && (
AC_MSG_RESULT(expanded our pkgconfig_libdir... $PKGCONFIG_libdir))
AC_SUBST([pkgconfig_libdir])

AC_MSG_CHECKING(our pkgconfig_libfile)
test ".$pkgconfig_libfile" != "." || \
pkgconfig_libfile="ifelse($1,,$PKGCONFIG_libname.pc,`basename $1`)"
PKGCONFIG_libfile=`eval echo "$pkgconfig_libfile"`
PKGCONFIG_libfile=`eval echo "$PKGCONFIG_libfile"`
AC_MSG_RESULT($pkgconfig_libfile)
test "$pkgconfig_libfile" != "$PKGCONFIG_libfile" && (
AC_MSG_RESULT(expanded our pkgconfig_libfile... $PKGCONFIG_libfile))
AC_SUBST([pkgconfig_libfile])

AC_MSG_CHECKING(our pkgconfig description)
PKGCONFIG_description="ifelse($4,,$PACKAGE_SUMMARY,$4)"
test ".$PKGCONFIG_description" != "." || \
PKGCONFIG_description="$PKGCONFIG_libname Library"
PKGCONFIG_description=`eval echo "$PKGCONFIG_description"`
PKGCONFIG_description=`eval echo "$PKGCONFIG_description"`
AC_MSG_RESULT($PKGCONFIG_description)

AC_MSG_CHECKING(our pkgconfig requires)
PKGCONFIG_requires="ifelse($2,,$PACKAGE_REQUIRES,$2)"
PKGCONFIG_requires=`eval echo "$PKGCONFIG_requires"`
PKGCONFIG_requires=`eval echo "$PKGCONFIG_requires"`
AC_MSG_RESULT($PKGCONFIG_requires)

AC_MSG_CHECKING(our pkgconfig ext libs)
PKGCONFIG_pkglibs="$PACKAGE_LIBS"
test ".$PKGCONFIG_pkglibs" != "." || PKGCONFIG_pkglibs="-l$PKGCONFIG_libname"
PKGCONFIG_libs="ifelse($3,,$PKGCONFIG_pkglibs $LIBS,$3)"
PKGCONFIG_libs=`eval echo "$PKGCONFIG_libs"`
PKGCONFIG_libs=`eval echo "$PKGCONFIG_libs"`
AC_MSG_RESULT($PKGCONFIG_libs)

AC_MSG_CHECKING(our pkgconfig cppflags)
PKGCONFIG_cppflags="ifelse($5,,$CPPFLAGS $PACKAGE_CFLAGS,$5)"
PKGCONFIG_cppflags=`eval echo "$PKGCONFIG_cppflags"`
PKGCONFIG_cppflags=`eval echo "$PKGCONFIG_cppflags"`
AC_MSG_RESULT($PKGCONFIG_cppflags)

test ".$PKGCONFIG_generate" != "." || \
PKGCONFIG_generate="ifelse($1,,$PKGCONFIG_libname.pc,$1)"
PKGCONFIG_generate=`eval echo "$PKGCONFIG_generate"`
PKGCONFIG_generate=`eval echo "$PKGCONFIG_generate"`
test "$pkgconfig_libfile" != "$PKGCONFIG_generate" && (
AC_MSG_RESULT(generate the pkgconfig later... $PKGCONFIG_generate))

if test ".$PKGCONFIG_src_libdir" = "." ; then
PKGCONFIG_src_libdir=`pwd`
PKGCONFIG_src_libdir=`AX_DIRNAME("$PKGCONFIG_src_libdir/$PKGCONFIG_generate")`
test ! -d $PKGCONFIG_src_libdir/src || \
PKGCONFIG_src_libdir="$PKGCONFIG_src_libdir/src"
case ".$objdir" in 
*libs) PKGCONFIG_src_libdir="$PKGCONFIG_src_libdir/$objdir" ;; esac
AC_MSG_RESULT(noninstalled pkgconfig -L $PKGCONFIG_src_libdir)
fi

if test ".$PKGCONFIG_src_headers" = "." ; then
PKGCONFIG_src_headers=`pwd`
v="$ac_top_srcdir" ; 
test ".$v" != "." || v="$ax_spec_dir"
test ".$v" != "." || v="$srcdir"
case "$v" in /*) PKG_CONFIG_src_headers="" ;; esac
PKGCONFIG_src_headers=`AX_DIRNAME("$PKGCONFIG_src_headers/$v/x")`
test ! -d $PKGCONFIG_src_headers/incl[]ude || \
PKGCONFIG_src_headers="$PKGCONFIG_src_headers/incl[]ude"
AC_MSG_RESULT(noninstalled pkgconfig -I $PKGCONFIG_src_headers)
fi

dnl AC_CONFIG_COMMANDS crap disallows to use $PKGCONFIG_libfile here...
AC_CONFIG_COMMANDS([$ax_create_pkgconfig_generate],[ 
pkgconfig_generate="$ax_create_pkgconfig_generate"
if test ! -f "$pkgconfig_generate.in" 
then generate="true"
elif grep ' generated by configure ' $pkgconfig_generate.in >/dev/null 
then generate="true"
else generate="false"; 
fi
if $generate ; then
AC_MSG_NOTICE(creating $pkgconfig_generate.in)
cat > $pkgconfig_generate.in <<_EOF
# generated by configure / remove this line to disable regeneration
prefix=@prefix@
exec_prefix=@exec_prefix@
libdir=@libdir@
includedir=@includedir@

Name: @PACKAGE_NAME@
Description: @PACKAGE_DESCRIPTION@
Version: @PACKAGE_VERSION@
Requires: @PACKAGE_REQUIRES@
Libs: -L\${libdir} @LIBS@
Cflags: -I\${includedir} @CPPFLAGS@
_EOF
fi # DONE generate $pkgconfig_generate.in
AC_MSG_NOTICE(creating $pkgconfig_generate)
sed -e "s|@prefix@|${pkgconfig_prefix}|" \
    -e "s|@exec_prefix@|${pkgconfig_execprefix}|" \
    -e "s|@libdir@|${pkgconfig_libdir}|" \
    -e "s|@includedir@|${pkgconfig_includedir}|" \
    -e "s|@PACKAGE_NAME@|${pkgconfig_libname}|" \
    -e "s|@PACKAGE_DESCRIPTION@|${pkgconfig_description}|" \
    -e "s|@PACKAGE_VERSION@|${pkgconfig_version}|" \
    -e "s|@PACKAGE_REQUIRES@|${pkgconfig_requires}|" \
    -e "s|@LIBS@|${pkgconfig_libs}|" \
    -e "s|@CPPFLAGS@|${pkgconfig_cppflags}|" \
    $pkgconfig_generate.in > $pkgconfig_generate
if test ! -s $pkgconfig_generate ; then 
    AC_MSG_ERROR([$pkgconfig_generate is empty])
fi # DONE generate $pkgconfig_generate
pkgconfig_uninstalled=`echo $pkgconfig_generate |sed 's/.pc$/-uninstalled.pc/'`
AC_MSG_NOTICE(creating $pkgconfig_uninstalled)
sed -e "s|@prefix@|${pkgconfig_prefix}|" \
    -e "s|@exec_prefix@|${pkgconfig_execprefix}|" \
    -e "s|@libdir@|${pkgconfig_src_libdir}|" \
    -e "s|@includedir@|${pkgconfig_src_headers}|" \
    -e "s|@PACKAGE_NAME@|${pkgconfig_libname}|" \
    -e "s|@PACKAGE_DESCRIPTION@|${pkgconfig_description}|" \
    -e "s|@PACKAGE_VERSION@|${pkgconfig_version}|" \
    -e "s|@PACKAGE_REQUIRES@|${pkgconfig_requires}|" \
    -e "s|@LIBS@|${pkgconfig_libs}|" \
    -e "s|@CPPFLAGS@|${pkgconfig_cppflags}|" \
    $pkgconfig_generate.in > $pkgconfig_uninstalled
if test ! -s $pkgconfig_uninstalled ; then 
    AC_MSG_ERROR([$pkgconfig_uninstalled is empty])
fi # DONE generate $pkgconfig_uninstalled
],[
dnl AC_CONFIG_COMMANDS crap, the AS_PUSHVAR defines are invalid here...
ax_create_pkgconfig_generate="$ax_create_pkgconfig_generate"
pkgconfig_prefix='$prefix'
pkgconfig_execprefix='$exec_prefix'
pkgconfig_libdir='$libdir'
pkgconfig_includedir='$includedir'
pkgconfig_libname='$ax_create_pkgconfig_libname'
pkgconfig_description='$ax_create_pkgconfig_description'
pkgconfig_version='$ax_create_pkgconfig_version'
pkgconfig_requires='$ax_create_pkgconfig_requires'
pkgconfig_libs='$ax_create_pkgconfig_libs'
pkgconfig_cppflags='$ax_create_pkgconfig_cppflags'
pkgconfig_src_libdir='$ax_create_pkgconfig_src_libdir'
pkgconfig_src_headers='$ax_create_pkgconfig_src_headers'
])dnl
AS_VAR_POPDEF([PKGCONFIG_libdir])dnl
AS_VAR_POPDEF([PKGCONFIG_libfile])dnl
AS_VAR_POPDEF([PKGCONFIG_libname])dnl
AS_VAR_POPDEF([PKGCONFIG_version])dnl
AS_VAR_POPDEF([PKGCONFIG_description])dnl
AS_VAR_POPDEF([PKGCONFIG_requires])dnl
AS_VAR_POPDEF([PKGCONFIG_pkglibs])dnl
AS_VAR_POPDEF([PKGCONFIG_libs])dnl
AS_VAR_POPDEF([PKGCONFIG_cppflags])dnl
AS_VAR_POPDEF([PKGCONFIG_generate])dnl
AS_VAR_POPDEF([PKGCONFIG_src_libdir])dnl
AS_VAR_POPDEF([PKGCONFIG_src_headers])dnl
])
