dnl @* AX_DIRNAME(PATHNAME)
dnl
dnl Parts of the implementation have been taken from AS_DIRNAME from the
dnl main autoconf package in generation 2.5x. However, we do only use 
dnl "sed" to cut out the dirname, and we do additionally clean up some 
dnl dir/.. parts in the resulting pattern.
dnl
dnl this macro may be used in autoconf 2.13 scripts as well.
dnl
dnl @$Id$

AC_DEFUN([AX_DIRNAME],
[echo X[]$1 |
    sed ['s/\/[^\/:][^\/:]*\/..\//\//g
          s/\/[^\/:][^\/:]*\/..\//\//g
          s/\/[^\/:][^\/:]*\/..\//\//g
          s/\/[^\/:][^\/:]*\/..\//\//g
          /^X\(.*[^/]\)\/\/*[^/][^/]*\/*$/{ s//\1/; q; }
          /^X\(\/\/\)[^/].*/{ s//\1/; q; }
          /^X\(\/\/\)$/{ s//\1/; q; }
          /^X\(\/\).*/{ s//\1/; q; }
          s/.*/./; q']])
