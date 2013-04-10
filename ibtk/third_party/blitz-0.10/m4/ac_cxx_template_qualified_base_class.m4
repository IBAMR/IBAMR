dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_cxx_template_qualified_base_class.html
dnl
AC_DEFUN([AC_CXX_TEMPLATE_QUALIFIED_BASE_CLASS],
[AC_CACHE_CHECK(whether the compiler supports template-qualified base class specifiers,
ac_cv_cxx_template_qualified_base_class,
[AC_REQUIRE([AC_CXX_TYPENAME])
 AC_REQUIRE([AC_CXX_FULL_SPECIALIZATION_SYNTAX])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
#ifndef HAVE_TYPENAME
 #define typename
#endif
class Base1 { public : int f () const { return 1; } };
class Base2 { public : int f () const { return 0; } };
template<class X> struct base_trait        { typedef Base1 base; };
#ifdef HAVE_FULL_SPECIALIZATION_SYNTAX
template<>        struct base_trait<float> { typedef Base2 base; };
#else
                  struct base_trait<float> { typedef Base2 base; };
#endif
template<class T> class Weird : public base_trait<T>::base
{ public :
  typedef typename base_trait<T>::base base;
  int g () const { return base::f (); }
};],[ Weird<float> z; return z.g ();],
 ac_cv_cxx_template_qualified_base_class=yes, ac_cv_cxx_template_qualified_base_class=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_template_qualified_base_class" = yes; then
  AC_DEFINE(HAVE_TEMPLATE_QUALIFIED_BASE_CLASS,,
            [define if the compiler supports template-qualified base class specifiers])
fi
])
