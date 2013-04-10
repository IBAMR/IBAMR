/***************************************************************************
 * blitz/array/ops.cc  Basic operators for arrays.
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
#ifndef BZ_ARRAYOPS_CC
#define BZ_ARRAYOPS_CC

#ifndef BZ_ARRAY_H
 #error <blitz/array/ops.cc> must be included via <blitz/array.h>
#endif

#include <blitz/update.h>
#include <blitz/globeval.cc>

BZ_NAMESPACE(blitz)

/*
 * Constant operands
 */

template<typename P_numtype, int N_rank>
_bz_forceinline
Array<P_numtype, N_rank>& Array<P_numtype,N_rank>::initialize(T_numtype x)
{
  // we can't use asExpr here, because if we are initializing an array
  // whose components are also ETBase, it would parse as an array
  // expression, not as an initialization with a scalar.
  (*this) = _bz_ArrayExpr<_bz_ArrayExprConstant<T_numtype> >(x);
  return *this;
}

#ifdef BZ_NEW_EXPRESSION_TEMPLATES

// need to do operator= separately from the generic templates below so
// that operator=(T_numtype) will be the best match for list
// initializations.
template<typename P_numtype, int N_rank> template<typename T_expr>
_bz_forceinline
Array<P_numtype,N_rank>&
Array<P_numtype,N_rank>::operator=(const ETBase<T_expr>& expr)
{
  _bz_evaluate(*this, _bz_typename asExpr<T_expr>::T_expr(expr.unwrap()),
	       _bz_update<T_numtype,
	       _bz_typename asExpr<T_expr>::T_expr::T_result>());
    return *this;
}

// do NOT remove this operator. it won't work without it, trust me...
template<typename P_numtype, int N_rank>
_bz_forceinline
Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator=(const Array<T_numtype,N_rank>& x)
{
  typedef typename asExpr<Array<T_numtype,N_rank> >::T_expr T_expr;
  _bz_evaluate(*this, asExpr<Array<T_numtype,N_rank> >::getExpr(x),
	       _bz_update<T_numtype,
	       _bz_typename T_expr::T_result>());
    return *this;
}

#define BZ_ARRAY_UPDATE(op,name)					\
  template<typename P_numtype, int N_rank>				\
  template<typename T_expr>						\
  _bz_forceinline							\
  Array<P_numtype,N_rank>&						\
  Array<P_numtype,N_rank>::operator op(const ETBase<T_expr>& expr)	\
  {									\
    _bz_evaluate(*this, _bz_typename asExpr<T_expr>::T_expr(expr.unwrap()), \
		 name<T_numtype, _bz_typename asExpr<T_expr>::T_expr::T_result>()); \
    return *this;							\
  }									\
  template<typename P_numtype, int N_rank>				\
  _bz_forceinline							\
  Array<P_numtype,N_rank>&						\
  Array<P_numtype,N_rank>::operator op(const Array<T_numtype, N_rank>& x) \
  {									\
    typedef typename asExpr<Array<T_numtype,N_rank> >::T_expr T_expr;	\
    _bz_evaluate(*this, asExpr<Array<T_numtype, N_rank> >::getExpr(x),	\
		 name<T_numtype, _bz_typename T_expr::T_result>());	\
    return *this;							\
  }									\
  template<typename P_numtype, int N_rank>				\
  _bz_forceinline							\
  Array<P_numtype,N_rank>&						\
  Array<P_numtype,N_rank>::operator op(const T_numtype& x)		\
  {									\
    typedef typename asExpr<T_numtype>::T_expr T_expr;			\
    _bz_evaluate(*this, asExpr<T_numtype>::getExpr(x),			\
		 name<T_numtype, _bz_typename T_expr::T_result>());	\
    return *this;							\
  }

BZ_ARRAY_UPDATE(+=, _bz_plus_update)
BZ_ARRAY_UPDATE(-=, _bz_minus_update)
BZ_ARRAY_UPDATE(*=, _bz_multiply_update)
BZ_ARRAY_UPDATE(/=, _bz_divide_update)
BZ_ARRAY_UPDATE(%=, _bz_mod_update)
BZ_ARRAY_UPDATE(^=, _bz_xor_update)
BZ_ARRAY_UPDATE(&=, _bz_bitand_update)
BZ_ARRAY_UPDATE(|=, _bz_bitor_update)
BZ_ARRAY_UPDATE(<<=, _bz_shiftl_update)
BZ_ARRAY_UPDATE(>>=, _bz_shiftr_update)

#else

template<typename P_numtype, int N_rank>
inline Array<P_numtype, N_rank>&
Array<P_numtype,N_rank>::operator+=(T_numtype x)
{
    (*this) += _bz_ArrayExpr<_bz_ArrayExprConstant<T_numtype> >(x);
    return *this;
}

template<typename P_numtype, int N_rank>
inline Array<P_numtype, N_rank>&
Array<P_numtype,N_rank>::operator-=(T_numtype x)
{
    (*this) -= _bz_ArrayExpr<_bz_ArrayExprConstant<T_numtype> >(x);
    return *this;
}

template<typename P_numtype, int N_rank>
inline Array<P_numtype, N_rank>&
Array<P_numtype,N_rank>::operator*=(T_numtype x)
{
    (*this) *= _bz_ArrayExpr<_bz_ArrayExprConstant<T_numtype> >(x);
    return *this;
}

template<typename P_numtype, int N_rank>
inline Array<P_numtype, N_rank>&
Array<P_numtype,N_rank>::operator/=(T_numtype x)
{
    (*this) /= _bz_ArrayExpr<_bz_ArrayExprConstant<T_numtype> >(x);
    return *this;
}

template<typename P_numtype, int N_rank>
inline Array<P_numtype, N_rank>&
Array<P_numtype,N_rank>::operator%=(T_numtype x)
{
    (*this) %= _bz_ArrayExpr<_bz_ArrayExprConstant<T_numtype> >(x);
    return *this;
}

template<typename P_numtype, int N_rank>
inline Array<P_numtype, N_rank>&
Array<P_numtype,N_rank>::operator^=(T_numtype x)
{
    (*this) ^= _bz_ArrayExpr<_bz_ArrayExprConstant<T_numtype> >(x);
    return *this;
}

template<typename P_numtype, int N_rank>
inline Array<P_numtype, N_rank>&
Array<P_numtype,N_rank>::operator&=(T_numtype x)
{
    (*this) &= _bz_ArrayExpr<_bz_ArrayExprConstant<T_numtype> >(x);
    return *this;
}

template<typename P_numtype, int N_rank>
inline Array<P_numtype, N_rank>&
Array<P_numtype,N_rank>::operator|=(T_numtype x)
{
    (*this) |= _bz_ArrayExpr<_bz_ArrayExprConstant<T_numtype> >(x);
    return *this;
}

template<typename P_numtype, int N_rank>
inline Array<P_numtype, N_rank>&
Array<P_numtype,N_rank>::operator>>=(T_numtype x)
{
    (*this) <<= _bz_ArrayExpr<_bz_ArrayExprConstant<T_numtype> >(x);
    return *this;
}

template<typename P_numtype, int N_rank>
inline Array<P_numtype, N_rank>&
Array<P_numtype,N_rank>::operator<<=(T_numtype x)
{
    (*this) <<= _bz_ArrayExpr<_bz_ArrayExprConstant<T_numtype> >(x);
    return *this;
}

/*
 * Array operands
 */

template<typename P_numtype, int N_rank>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator=(const Array<T_numtype,N_rank>& x)
{
    (*this) = _bz_ArrayExpr<FastArrayIterator<T_numtype, N_rank> >(x.beginFast());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename P_numtype2>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator=(const Array<P_numtype2,N_rank>& x)
{
    (*this) = _bz_ArrayExpr<FastArrayIterator<P_numtype2, N_rank> >(x.beginFast());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename P_numtype2>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator+=(const Array<P_numtype2,N_rank>& x)
{
    (*this) += _bz_ArrayExpr<FastArrayIterator<P_numtype2, N_rank> >(x.beginFast());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename P_numtype2>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator-=(const Array<P_numtype2,N_rank>& x)
{
    (*this) -= _bz_ArrayExpr<FastArrayIterator<P_numtype2, N_rank> >(x.beginFast());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename P_numtype2>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator*=(const Array<P_numtype2,N_rank>& x)
{
    (*this) *= _bz_ArrayExpr<FastArrayIterator<P_numtype2, N_rank> >(x.beginFast());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename P_numtype2>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator/=(const Array<P_numtype2,N_rank>& x)
{
    (*this) /= _bz_ArrayExpr<FastArrayIterator<P_numtype2, N_rank> >(x.beginFast());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename P_numtype2>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator%=(const Array<P_numtype2,N_rank>& x)
{
    (*this) %= _bz_ArrayExpr<FastArrayIterator<P_numtype2, N_rank> >(x.beginFast());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename P_numtype2>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator^=(const Array<P_numtype2,N_rank>& x)
{
    (*this) ^= _bz_ArrayExpr<FastArrayIterator<P_numtype2, N_rank> >(x.beginFast());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename P_numtype2>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator&=(const Array<P_numtype2,N_rank>& x)
{
    (*this) &= _bz_ArrayExpr<FastArrayIterator<P_numtype2, N_rank> >(x.beginFast());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename P_numtype2>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator|=(const Array<P_numtype2,N_rank>& x)
{
    (*this) |= _bz_ArrayExpr<FastArrayIterator<P_numtype2, N_rank> >(x.beginFast());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename P_numtype2>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator>>=(const Array<P_numtype2,N_rank>& x)
{
    (*this) >>= _bz_ArrayExpr<FastArrayIterator<P_numtype2, N_rank> >(x.beginFast());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename P_numtype2>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator<<=(const Array<P_numtype2,N_rank>& x)
{
    (*this) <<= _bz_ArrayExpr<FastArrayIterator<P_numtype2, N_rank> >(x.beginFast());
    return *this;
}

/*
 * Array expression operands
 */

template<typename P_numtype, int N_rank> template<typename T_expr>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator=(BZ_ETPARM(_bz_ArrayExpr<T_expr>) expr)
{
    evaluate(expr, _bz_update<T_numtype, _bz_typename T_expr::T_numtype>());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename T_expr>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator+=(BZ_ETPARM(_bz_ArrayExpr<T_expr>) expr)
{
    evaluate(expr, _bz_plus_update<T_numtype, _bz_typename T_expr::T_numtype>());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename T_expr>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator-=(BZ_ETPARM(_bz_ArrayExpr<T_expr>) expr)
{
    evaluate(expr, _bz_minus_update<T_numtype,
        _bz_typename T_expr::T_numtype>());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename T_expr>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator*=(BZ_ETPARM(_bz_ArrayExpr<T_expr>) expr)
{
    evaluate(expr, _bz_multiply_update<T_numtype,
        _bz_typename T_expr::T_numtype>());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename T_expr>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator/=(BZ_ETPARM(_bz_ArrayExpr<T_expr>) expr)
{
    evaluate(expr, _bz_divide_update<T_numtype,
        _bz_typename T_expr::T_numtype>());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename T_expr>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator%=(BZ_ETPARM(_bz_ArrayExpr<T_expr>) expr)
{
    evaluate(expr, _bz_mod_update<T_numtype,
        _bz_typename T_expr::T_numtype>());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename T_expr>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator^=(BZ_ETPARM(_bz_ArrayExpr<T_expr>) expr)
{
    evaluate(expr, _bz_xor_update<T_numtype,
        _bz_typename T_expr::T_numtype>());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename T_expr>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator&=(BZ_ETPARM(_bz_ArrayExpr<T_expr>) expr)
{
    evaluate(expr, _bz_bitand_update<T_numtype,
        _bz_typename T_expr::T_numtype>());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename T_expr>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator|=(BZ_ETPARM(_bz_ArrayExpr<T_expr>) expr)
{
    evaluate(expr, _bz_bitor_update<T_numtype,
        _bz_typename T_expr::T_numtype>());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename T_expr>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator>>=(BZ_ETPARM(_bz_ArrayExpr<T_expr>) expr)
{
    evaluate(expr, _bz_shiftr_update<T_numtype,
        _bz_typename T_expr::T_numtype>());
    return *this;
}

template<typename P_numtype, int N_rank> template<typename T_expr>
inline Array<P_numtype, N_rank>&
Array<P_numtype, N_rank>::operator<<=(BZ_ETPARM(_bz_ArrayExpr<T_expr>) expr)
{
    evaluate(expr, _bz_shiftl_update<T_numtype,
        _bz_typename T_expr::T_numtype>());
    return *this;
}

#endif // BZ_NEW_EXPRESSION_TEMPLATES

BZ_NAMESPACE_END

#endif // BZ_ARRAYOPS_CC
