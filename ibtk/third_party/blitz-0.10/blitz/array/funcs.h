// -*- C++ -*-
/***************************************************************************
 * blitz/array/funcs.h   Math functions on arrays
 *
 * $Id$
 *
 * Copyright (C) 1997-2011 Todd Veldhuizen <tveldhui@acm.org>
 *
 * This file is a part of Blitz.
 *
 * Blitz is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version.
 *
 * Blitz is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Blitz.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Suggestions:          blitz-devel@lists.sourceforge.net
 * Bugs:                 blitz-support@lists.sourceforge.net
 *
 * For more information, please see the Blitz++ Home Page:
 *    https://sourceforge.net/projects/blitz/
 *
 ****************************************************************************/
#ifndef BZ_ARRAY_FUNCS_H
#define BZ_ARRAY_FUNCS_H

#include <blitz/funcs.h>
#include <blitz/array/newet-macros.h>
#include <blitz/array/reduce.h>
#include <blitz/levicivita.h>

BZ_NAMESPACE(blitz)

/* This file contains definitions of the math functions for ETBase
   types. Note that these definitions are in the blitz namespace and
   will hide the definitions of the builtin math functions in the
   global namespace. For this reason, the calls to the builtin math
   functions in the functors' apply() methods must explicitly qualify
   the namespace.
*/


// unary functions

BZ_DECLARE_ARRAY_ET_UNARY(abs,   Fn_abs)
BZ_DECLARE_ARRAY_ET_UNARY(acos,  Fn_acos)
BZ_DECLARE_ARRAY_ET_UNARY(asin,  Fn_asin)
BZ_DECLARE_ARRAY_ET_UNARY(atan,  Fn_atan)
BZ_DECLARE_ARRAY_ET_UNARY(ceil,  Fn_ceil)
BZ_DECLARE_ARRAY_ET_UNARY(cexp,  Fn_exp)
BZ_DECLARE_ARRAY_ET_UNARY(cos,   Fn_cos)
BZ_DECLARE_ARRAY_ET_UNARY(cosh,  Fn_cosh)
BZ_DECLARE_ARRAY_ET_UNARY(csqrt, Fn_sqrt)
BZ_DECLARE_ARRAY_ET_UNARY(cube,  Fn_cube)
BZ_DECLARE_ARRAY_ET_UNARY(exp,   Fn_exp)
BZ_DECLARE_ARRAY_ET_UNARY(fabs,  Fn_fabs)
BZ_DECLARE_ARRAY_ET_UNARY(floor, Fn_floor)
BZ_DECLARE_ARRAY_ET_UNARY(log,   Fn_log)
BZ_DECLARE_ARRAY_ET_UNARY(log10, Fn_log10)
BZ_DECLARE_ARRAY_ET_UNARY(pow2,  Fn_sqr)
BZ_DECLARE_ARRAY_ET_UNARY(pow3,  Fn_cube)
BZ_DECLARE_ARRAY_ET_UNARY(pow4,  Fn_pow4)
BZ_DECLARE_ARRAY_ET_UNARY(pow5,  Fn_pow5)
BZ_DECLARE_ARRAY_ET_UNARY(pow6,  Fn_pow6)
BZ_DECLARE_ARRAY_ET_UNARY(pow7,  Fn_pow7)
BZ_DECLARE_ARRAY_ET_UNARY(pow8,  Fn_pow8)
BZ_DECLARE_ARRAY_ET_UNARY(sin,   Fn_sin)
BZ_DECLARE_ARRAY_ET_UNARY(sinh,  Fn_sinh)
BZ_DECLARE_ARRAY_ET_UNARY(sqr,   Fn_sqr)
BZ_DECLARE_ARRAY_ET_UNARY(sqrt,  Fn_sqrt)
BZ_DECLARE_ARRAY_ET_UNARY(tan,   Fn_tan)
BZ_DECLARE_ARRAY_ET_UNARY(tanh,  Fn_tanh)

#ifdef BZ_HAVE_COMPLEX_FCNS
BZ_DECLARE_ARRAY_ET_UNARY(arg,   Fn_arg)
BZ_DECLARE_ARRAY_ET_UNARY(conj,  Fn_conj)
BZ_DECLARE_ARRAY_ET_UNARY(imag,  Fn_imag)
BZ_DECLARE_ARRAY_ET_UNARY(norm,  Fn_norm)
BZ_DECLARE_ARRAY_ET_UNARY(real,  Fn_real)
#endif

#ifdef BZ_HAVE_IEEE_MATH
// finite and trunc omitted: blitz-bugs/archive/0189.html
BZ_DECLARE_ARRAY_ET_UNARY(acosh,  Fn_acosh)
BZ_DECLARE_ARRAY_ET_UNARY(asinh,  Fn_asinh)
BZ_DECLARE_ARRAY_ET_UNARY(atanh,  Fn_atanh)
BZ_DECLARE_ARRAY_ET_UNARY(cbrt,   Fn_cbrt)
BZ_DECLARE_ARRAY_ET_UNARY(erf,    Fn_erf)
BZ_DECLARE_ARRAY_ET_UNARY(erfc,   Fn_erfc)
BZ_DECLARE_ARRAY_ET_UNARY(expm1,  Fn_expm1)
// BZ_DECLARE_ARRAY_ET_UNARY(finite, Fn_finite)
BZ_DECLARE_ARRAY_ET_UNARY(ilogb,   Fn_ilogb)
BZ_DECLARE_ARRAY_ET_UNARY(blitz_isnan,  Fn_isnan)
BZ_DECLARE_ARRAY_ET_UNARY(j0,     Fn_j0)
BZ_DECLARE_ARRAY_ET_UNARY(j1,     Fn_j1)
BZ_DECLARE_ARRAY_ET_UNARY(lgamma, Fn_lgamma)
BZ_DECLARE_ARRAY_ET_UNARY(logb,   Fn_logb)
BZ_DECLARE_ARRAY_ET_UNARY(log1p,  Fn_log1p)
BZ_DECLARE_ARRAY_ET_UNARY(rint,   Fn_rint)
// BZ_DECLARE_ARRAY_ET_UNARY(trunc,  Fn_trunc)
BZ_DECLARE_ARRAY_ET_UNARY(y0,     Fn_y0)
BZ_DECLARE_ARRAY_ET_UNARY(y1,     Fn_y1)
#endif

#ifdef BZ_HAVE_SYSTEM_V_MATH
BZ_DECLARE_ARRAY_ET_UNARY(_class,  Fn__class)
BZ_DECLARE_ARRAY_ET_UNARY(itrunc,  Fn_itrunc)
BZ_DECLARE_ARRAY_ET_UNARY(nearest, Fn_nearest)
BZ_DECLARE_ARRAY_ET_UNARY(rsqrt,   Fn_rsqrt)
BZ_DECLARE_ARRAY_ET_UNARY(uitrunc, Fn_uitrunc)
#endif

// cast() function

template<typename T_cast, typename T1>
_bz_inline_et
_bz_ArrayExpr<_bz_ArrayExprUnaryOp<_bz_typename asExpr<T1>::T_expr,
    Cast<_bz_typename asExpr<T1>::T_expr::T_numtype, T_cast> > >
cast(const ETBase<T1>& expr)
{
    return _bz_ArrayExpr<_bz_ArrayExprUnaryOp<
        _bz_typename asExpr<T1>::T_expr,
        Cast<_bz_typename asExpr<T1>::T_expr::T_numtype,T_cast> > >
        (expr.unwrap());
}

// binary functions

BZ_DECLARE_ARRAY_ET_BINARY(atan2,     Fn_atan2)
BZ_DECLARE_ARRAY_ET_BINARY(fmod,      Fn_fmod)
BZ_DECLARE_ARRAY_ET_BINARY(pow,       Fn_pow)

#ifdef BZ_HAVE_COMPLEX_FCNS
BZ_DECLARE_ARRAY_ET_BINARY(polar,     Fn_polar)
#endif

#ifdef BZ_HAVE_SYSTEM_V_MATH
BZ_DECLARE_ARRAY_ET_BINARY(copysign,  Fn_copysign)
BZ_DECLARE_ARRAY_ET_BINARY(drem,      Fn_drem)
BZ_DECLARE_ARRAY_ET_BINARY(hypot,     Fn_hypot)
BZ_DECLARE_ARRAY_ET_BINARY(nextafter, Fn_nextafter)
BZ_DECLARE_ARRAY_ET_BINARY(remainder, Fn_remainder)
BZ_DECLARE_ARRAY_ET_BINARY(scalb,     Fn_scalb)
BZ_DECLARE_ARRAY_ET_BINARY(unordered, Fn_unordered)
#endif

BZ_DECLARE_ARRAY_ET_BINARY((min), Min)
BZ_DECLARE_ARRAY_ET_BINARY((max), Max)

#ifdef BZ_HAVE_SYSTEM_V_MATH

#define BZ_DECLARE_ARRAY_ET_SCALAR_FUNCS(sca)                   \
                                                                \
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(atan2,     Fn_atan2, sca)     \
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(fmod,      Fn_fmod, sca)      \
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(pow,       Fn_pow, sca)       \
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(copysign,  Fn_copysign, sca)  \
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(drem,      Fn_drem, sca)      \
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(hypot,     Fn_hypot, sca)     \
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(nextafter, Fn_nextafter, sca) \
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(remainder, Fn_remainder, sca) \
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(scalb,     Fn_scalb, sca)     \
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(unordered, Fn_unordered, sca) \

#else

#define BZ_DECLARE_ARRAY_ET_SCALAR_FUNCS(sca)                   \
                                                                \
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(atan2,     Fn_atan2, sca)     \
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(fmod,      Fn_fmod, sca)      \
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(pow,       Fn_pow, sca)       \

#endif

BZ_DECLARE_ARRAY_ET_SCALAR_FUNCS(int)
BZ_DECLARE_ARRAY_ET_SCALAR_FUNCS(float)
BZ_DECLARE_ARRAY_ET_SCALAR_FUNCS(double)
BZ_DECLARE_ARRAY_ET_SCALAR_FUNCS(long double)

#ifdef BZ_HAVE_COMPLEX_FCNS
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(polar,     Fn_polar, int)
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(polar,     Fn_polar, float)
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(polar,     Fn_polar, double)
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(polar,     Fn_polar, long double)

template<typename T1, typename T2>
inline _bz_ArrayExprBinaryOp<
    typename asExpr<complex<T1> >::T_expr,
    typename asExpr<T2>::T_expr,
    Fn_pow<complex<T1>,typename asExpr<T2>::T_expr::T_numtype> >
pow(const complex<T1> d1, const ETBase<T2>& d2)
{
    return _bz_ArrayExprBinaryOp<
        typename asExpr<complex<T1> >::T_expr,
        typename asExpr<T2>::T_expr,
        Fn_pow<complex<T1>,typename asExpr<T2>::T_expr::T_numtype> >
        (asExpr<complex<T1> >::getExpr(d1),
         asExpr<T2>::getExpr(d2.unwrap()));
}


// global functions that don't fit anywhere else

// we define a generalized dot product for all classes as sum(a*b)
template<typename T1, typename T2>
inline
_bz_typename ReduceSum<_bz_typename BZ_BLITZ_SCOPE(BzBinaryExprResult)<Multiply,T1,T2>::T_result::T_numtype
>::T_resulttype
dot(const ETBase<T1>& d1, const ETBase<T2>& d2)
{
  return sum(d1 * d2);
}

// we define a generalized cross product for all classes using the
// Levi-Civita symbol. Return type is nice (ever heard of "write-once
// code")... it took 10 times longer to figure out how to write the
// return type than to do everything else...

template<typename T1, typename T2>
inline

_bz_ArrayExpr<
  _bz_ArrayExprReduce<
    _bz_ArrayExpr<
      _bz_ArrayExprReduce<
	_bz_typename BzBinaryExprResult<
	  Multiply,
	  _bz_typename BzBinaryExprResult<
	    Multiply,
	    _bz_ArrayExpr<LeviCivita>,
	    _bz_ArrayExpr<
	      ArrayIndexMapping<
		_bz_typename asExpr<T1>::T_expr, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0>
	      >
	    >::T_result,
	  _bz_ArrayExpr<
	    ArrayIndexMapping<
	      _bz_typename asExpr<T2>::T_expr, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0>
	    >
	  >::T_result,
	2,
	ReduceSum<
	  _bz_typename BzBinaryExprResult<Multiply,_bz_typename BzBinaryExprResult<Multiply,_bz_ArrayExpr<LeviCivita>,_bz_ArrayExpr<ArrayIndexMapping<_bz_typename asExpr<T1>::T_expr, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> > >::T_result,_bz_ArrayExpr<ArrayIndexMapping<_bz_typename asExpr<T2>::T_expr, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> > > ::T_result::T_numtype,
	  BZ_SUMTYPE(bzCC(_bz_typename BzBinaryExprResult<Multiply,_bz_typename BzBinaryExprResult<Multiply,_bz_ArrayExpr<LeviCivita>,_bz_ArrayExpr<ArrayIndexMapping<_bz_typename asExpr<T1>::T_expr, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> > >::T_result,_bz_ArrayExpr<ArrayIndexMapping<_bz_typename asExpr<T2>::T_expr, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> > > ::T_result::T_numtype))>
	>
      >,
    1,
    ReduceSum<BZ_SUMTYPE(bzCC(_bz_typename BzBinaryExprResult<Multiply,_bz_typename BzBinaryExprResult<Multiply,_bz_ArrayExpr<LeviCivita>,_bz_ArrayExpr<ArrayIndexMapping<_bz_typename asExpr<T1>::T_expr, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> > >::T_result,_bz_ArrayExpr<ArrayIndexMapping<_bz_typename asExpr<T2>::T_expr, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> > > ::T_result::T_numtype))>
    >
  >

//int
cross(const ETBase<T1>& d1, const ETBase<T2>& d2)
{
_bz_ArrayExpr<
  _bz_ArrayExprReduce<
    _bz_ArrayExpr<
      _bz_ArrayExprReduce<
	_bz_typename BzBinaryExprResult<
	  Multiply,
	  _bz_typename BzBinaryExprResult<
	    Multiply,
	    _bz_ArrayExpr<LeviCivita>,
	    _bz_ArrayExpr<
	      ArrayIndexMapping<
		_bz_typename asExpr<T1>::T_expr, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0>
	      >
	    >::T_result,
	  _bz_ArrayExpr<
	    ArrayIndexMapping<
	      _bz_typename asExpr<T2>::T_expr, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0>
	    >
	  >::T_result,
	2,
	ReduceSum<
	  _bz_typename BzBinaryExprResult<Multiply,_bz_typename BzBinaryExprResult<Multiply,_bz_ArrayExpr<LeviCivita>,_bz_ArrayExpr<ArrayIndexMapping<_bz_typename asExpr<T1>::T_expr, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> > >::T_result,_bz_ArrayExpr<ArrayIndexMapping<_bz_typename asExpr<T2>::T_expr, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> > > ::T_result::T_numtype,
	  BZ_SUMTYPE(bzCC(_bz_typename BzBinaryExprResult<Multiply,_bz_typename BzBinaryExprResult<Multiply,_bz_ArrayExpr<LeviCivita>,_bz_ArrayExpr<ArrayIndexMapping<_bz_typename asExpr<T1>::T_expr, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> > >::T_result,_bz_ArrayExpr<ArrayIndexMapping<_bz_typename asExpr<T2>::T_expr, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> > > ::T_result::T_numtype))>
	>
      >,
    1,
    ReduceSum<BZ_SUMTYPE(bzCC(_bz_typename BzBinaryExprResult<Multiply,_bz_typename BzBinaryExprResult<Multiply,_bz_ArrayExpr<LeviCivita>,_bz_ArrayExpr<ArrayIndexMapping<_bz_typename asExpr<T1>::T_expr, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> > >::T_result,_bz_ArrayExpr<ArrayIndexMapping<_bz_typename asExpr<T2>::T_expr, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> > > ::T_result::T_numtype))>
    >
  >* x;
//int a=*x;
  return sum(sum(LeviCivita()*d1.unwrap()(tensor::j)*d2.unwrap()(tensor::k),
		 tensor::k),
	     tensor::j);
}


// "scalar" function converts something into a scalar expression. this
// is used to prevent parsing of multicomponent expressions as bitwise
// expressions between dissimilar containers.

template <typename T>
inline
_bz_ArrayExpr<_bz_ArrayExprConstant<T> >
scalar(const T& x) {
  return _bz_ArrayExpr<_bz_ArrayExprConstant<T> >(x);
}

#endif

BZ_NAMESPACE_END

#endif // BZ_ARRAY_FUNCS_H
