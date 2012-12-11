dnl Check for doxygen to create API docs
dnl
AC_DEFUN([AC_PROG_DOXYGEN],
[
AC_ARG_ENABLE(doxygen,
AS_HELP_STRING([--enable-doxygen],[enable documentation generation with doxygen (no)]),[],[enable_doxygen=no])

AC_ARG_ENABLE(dot,
AS_HELP_STRING([--enable-dot],[use 'dot' to generate graphs in doxygen (no)]),[],[enable_dot=no])              
AC_ARG_ENABLE(html-docs,
AS_HELP_STRING([--enable-html-docs],[enable HTML generation with doxygen (no)]),
[],[ enable_html_docs=no])              
AC_ARG_ENABLE(latex-docs,
AS_HELP_STRING([--enable-latex-docs],
[enable LaTeX documentation generation with doxygen (no)]),[],[enable_latex_docs=no])              

if test "x$enable_doxygen" = xno; then
 enable_doc=no	
else 
 AC_CHECK_PROG(DOXYGEN, doxygen, doxygen)
 if test x$DOXYGEN = x; then
  if test "x$enable_doxygen" = xyes; then
   AC_MSG_ERROR([could not find doxygen])
  fi
  enable_doc=no
 else
  doxy_ver=`doxygen --version`
  doxy_major=`expr "$doxy_ver" : '\(@<:@0-9@:>@\)\..*'`
  doxy_minor=`expr "$doxy_ver" : '@<:@0-9@:>@\.\(@<:@0-9@:>@\).*'`
  if test $doxy_major -eq "1" -a $doxy_minor -ge "5" ; then
   enable_doc=yes
   AC_CHECK_PROG(DOT, dot, dot)
  else
   AC_MSG_WARN([doxygen version $doxy_ver too old, doxygen will not be used.])
   enable_doc=no
  fi
 fi
fi

AM_CONDITIONAL(DOXYGEN_DOC, test x$enable_doc = xyes)

if test x$DOT = x; then
 if test "x$enable_dot" = xyes; then
  AC_MSG_ERROR([could not find dot])
 fi
 enable_dot=no
else
 enable_dot=yes
fi

AC_SUBST(enable_dot)
AC_SUBST(enable_html_docs)
AC_SUBST(enable_latex_docs)
])

