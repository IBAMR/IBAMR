


AC_DEFUN([AC_CXX_TYPE_PROMOTION],
[AC_CACHE_CHECK(whether the compiler will support the vector type promotion mechanism,
ac_cv_cxx_type_promotion,
[AC_REQUIRE([AC_CXX_TYPENAME])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
#ifndef HAVE_TYPENAME
 #define typename
#endif
template <class T> struct vec3 { T data_[3]; };
template <class T1, class T2> struct promote_trait { typedef T1 T_promote; };
template <> struct promote_trait<int,double> { typedef double T_promote; };
template <class T1, class T2> vec3<typename promote_trait<T1,T2>::T_promote>
operator+(const vec3<T1>& a, const vec3<T2>& b) 
{ vec3<typename promote_trait<T1,T2>::T_promote> c;
  c.data_[0] = a.data_[0] + b.data_[0];
  c.data_[1] = a.data_[1] + b.data_[1];
  c.data_[2] = a.data_[2] + b.data_[2]; return c; }],[
vec3<int> a,b; vec3<double> c,d,e; b=a+a; d=c+c; e=b+d; return 0;],
 ac_cv_cxx_type_promotion=yes, ac_cv_cxx_type_promotion=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_type_promotion" = yes; then
  AC_DEFINE(HAVE_TYPE_PROMOTION,,
            [define if the compiler supports the vector type promotion mechanism])
fi
])


