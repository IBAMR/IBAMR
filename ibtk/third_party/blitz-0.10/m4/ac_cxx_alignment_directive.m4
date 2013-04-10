dnl By Patrik Jonsson 2011
AC_DEFUN([AC_ALIGNMENT_DIRECTIVE],
[AC_CACHE_CHECK([whether the compiler supports structure alignment hints],
ac_cv_alignment_directive,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
dnl this is the default definition

dnl first test if intel-style directive is recognized
AC_TRY_COMPILE([],[
__declspec(align(16))
int var; var=0;],
   ac_cv_alignment_directive=["__declspec(align(alignment)) vartype varname;"],
   ac_cv_alignment_directive=["no"])
#echo intel test got $ac_cv_alignment_directive

if test "$ac_cv_alignment_directive" == "no"; then
  AC_TRY_COMPILE([],[int __attribute__ ((aligned (16))) var;var=0;],
   ac_cv_alignment_directive=["vartype __attribute__ ((aligned (alignment))) varname;"], )
 #echo gcc test got $ac_cv_alignment_directive
fi
 AC_LANG_RESTORE
]
)

if test "$ac_cv_alignment_directive" == "no"; then
   align_variable_macro="vartype varname;"
else
   align_variable_macro="$ac_cv_alignment_directive"
fi

AH_TEMPLATE([ALIGN_VARIABLE], [Macro for declaring aligned variables])
AC_DEFINE_UNQUOTED([ALIGN_VARIABLE(vartype,varname,alignment)],[$align_variable_macro],[])
])

