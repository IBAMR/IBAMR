#ifndef _QD_QD_CONFIG_H
#define _QD_QD_CONFIG_H  1

#ifndef QD_API
#undef QD_API
#endif

/* Set to 1 if using VisualAge C++ compiler for __fmadd builtin. */
#ifndef QD_VACPP_BUILTINS_H
#undef QD_VACPP_BUILTINS_H
#endif

/* If fused multiply-add is available, define to correct macro for
   using it.  It is invoked as QD_FMA(a, b, c) to compute fl(a * b + c). 
   If correctly rounded multiply-add is not available (or if unsure), 
   keep it undefined.*/
#ifndef QD_FMA
#undef QD_FMA
#endif

/* If fused multiply-subtract is available, define to correct macro for
   using it.  It is invoked as QD_FMS(a, b, c) to compute fl(a * b - c). 
   If correctly rounded multiply-add is not available (or if unsure), 
   keep it undefined.*/
#ifndef QD_FMS
#undef QD_FMS
#endif

/* Set the following to 1 to define commonly used function
   to be inlined.  This should be set to 1 unless the compiler 
   does not support the "inline" keyword, or if building for 
   debugging purposes. */
#ifndef QD_INLINE
#undef QD_INLINE
#endif

/* Set the following to 1 to use ANSI C++ standard header files
   such as cmath, iostream, etc.  If set to zero, it will try to 
   include math.h, iostream.h, etc, instead. */
#ifndef QD_HAVE_STD
#undef QD_HAVE_STD 
#endif

/* Set the following to 1 to make the addition and subtraction
   to satisfy the IEEE-style error bound 

      fl(a + b) = (1 + d) * (a + b)

   where |d| <= eps.  If set to 0, the addition and subtraction
   will satisfy the weaker Cray-style error bound

      fl(a + b) = (1 + d1) * a + (1 + d2) * b

   where |d1| <= eps and |d2| eps.           */
#ifndef QD_IEEE_ADD
#undef QD_IEEE_ADD
#endif

/* Set the following to 1 to use slightly inaccurate but faster
   version of multiplication. */
#ifndef QD_SLOPPY_MUL
#undef QD_SLOPPY_MUL
#endif

/* Set the following to 1 to use slightly inaccurate but faster
   version of division. */
#ifndef QD_SLOPPY_DIV
#undef QD_SLOPPY_DIV
#endif

/* Define this macro to be the isfinite(x) function. */
#ifndef QD_ISFINITE
#undef QD_ISFINITE
#endif

/* Define this macro to be the isinf(x) function. */
#ifndef QD_ISINF
#undef QD_ISINF
#endif

/* Define this macro to be the isnan(x) function. */
#ifndef QD_ISNAN
#undef QD_ISNAN
#endif


#endif /* _QD_QD_CONFIG_H */
