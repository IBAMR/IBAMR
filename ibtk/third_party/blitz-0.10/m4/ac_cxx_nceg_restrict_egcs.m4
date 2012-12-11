

AC_DEFUN([AC_CXX_NCEG_RESTRICT_EGCS],
[AC_CACHE_CHECK(whether the compiler recognizes the '__restrict__' keyword,
ac_cv_cxx_nceg_restrict_egcs,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
void add(int length, double * __restrict__ a, 
         const double * __restrict__ b, const double * __restrict__ c)
{ for (int i=0; i < length; ++i) a[i] = b[i] + c[i]; }
],[double a[10], b[10], c[10];
for (int i=0; i < 10; ++i) { a[i] = 0.0; b[i] = 0.0; c[i] = 0.0; }
add(10,a,b,c);
return 0;],
 ac_cv_cxx_nceg_restrict_egcs=yes, ac_cv_cxx_nceg_restrict_egcs=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_nceg_restrict_egcs" = yes; then
  AC_DEFINE(HAVE_NCEG_RESTRICT_EGCS,,
            [define if the compiler supports the __restrict__ keyword])
fi
])


