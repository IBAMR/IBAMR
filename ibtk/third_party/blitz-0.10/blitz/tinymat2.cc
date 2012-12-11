/***************************************************************************
 * blitz/tinyvec.cc  Declaration of TinyVector methods
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
 ***************************************************************************/

#ifndef BZ_TINYMAT2_CC
#define BZ_TINYMAT2_CC

#include <blitz/tmevaluate.h>
#include <blitz/tinymat2io.cc>
#include <blitz/array/domain.h>

BZ_NAMESPACE(blitz)

template<typename P_numtype, int N_rows, int N_columns>
inline TinyMatrix<P_numtype, N_rows, N_columns>::TinyMatrix(T_numtype initValue)
{
  for (sizeType i=0; i < numElements(); ++i)
    data_[i] = initValue;
}

template<typename P_numtype, int N_rows, int N_columns>
inline TinyMatrix<P_numtype, N_rows, N_columns>::TinyMatrix(const T_matrix& x)
{
  for (sizeType i=0; i < numElements(); ++i)
    data_[i] = x.data_[i];
}

template<typename P_numtype, int N_rows, int N_columns>
template<typename P_numtype2>
inline
TinyMatrix<P_numtype, N_rows, N_columns>::TinyMatrix(const TinyMatrix<P_numtype2, N_rows, N_columns>& x)
{
  for (sizeType i=0; i < numElements(); ++i)
    data_[i] = static_cast<P_numtype>(x.data_[i]);
}

/*
 * Assignment-type operators
 */

template<typename P_numtype, int N_rows, int N_columns>
TinyMatrix<P_numtype, N_rows, N_columns>&
TinyMatrix<P_numtype, N_rows, N_columns>::initialize(T_numtype x)
{
    (*this) = _bz_ArrayExpr<_bz_ArrayExprConstant<T_numtype> >(x);
    return *this;
}

template<typename P_numtype, int N_rows, int N_columns>
template<typename T_expr>
inline
TinyMatrix<P_numtype, N_rows, N_columns>&
TinyMatrix<P_numtype, N_rows, N_columns>::operator=(const ETBase<T_expr>& expr)
{
  _tm_evaluate(_bz_typename asExpr<T_expr>::T_expr(expr.unwrap()),
	       _bz_update<T_numtype, _bz_typename asExpr<T_expr>::T_expr::T_result>());
    return *this;
}

#define BZ_TM2_UPDATE(op,name)						\
  template<typename P_numtype, int N_rows, int N_columns>		\
  template<typename T>							\
  inline TinyMatrix<P_numtype, N_rows, N_columns>&			\
  TinyMatrix<P_numtype, N_rows, N_columns>::operator op(const T& expr)	\
  {									\
    _tm_evaluate(_bz_typename asExpr<T>::T_expr(expr),			\
		 name<T_numtype, _bz_typename asExpr<T>::T_expr::T_result>()); \
  return *this;								\
  }

BZ_TM2_UPDATE(+=, _bz_plus_update)
BZ_TM2_UPDATE(-=, _bz_minus_update)
BZ_TM2_UPDATE(*=, _bz_multiply_update)
BZ_TM2_UPDATE(/=, _bz_divide_update)
BZ_TM2_UPDATE(%=, _bz_mod_update)
BZ_TM2_UPDATE(^=, _bz_xor_update)
BZ_TM2_UPDATE(&=, _bz_bitand_update)
BZ_TM2_UPDATE(|=, _bz_bitor_update)
BZ_TM2_UPDATE(<<=, _bz_shiftl_update)
BZ_TM2_UPDATE(>>=, _bz_shiftr_update)

/*
 * Other member functions
 */


template<typename P_numtype, int N_rows, int N_columns>
inline RectDomain<2>
TinyMatrix<P_numtype, N_rows, N_columns>::domain()
{
  return RectDomain<2>(lbound(), ubound());
}


template<typename P_numtype, int N_rows, int N_columns>
template<int N0, int N1>
inline _bz_ArrayExpr<ArrayIndexMapping<typename asExpr<TinyMatrix<P_numtype, N_rows, N_columns> >::T_expr, N0, N1> >
TinyMatrix<P_numtype, N_rows, N_columns>::operator()(IndexPlaceholder<N0>, IndexPlaceholder<N1>) const
{
  return _bz_ArrayExpr<ArrayIndexMapping<typename asExpr<T_matrix>::T_expr,
					 N0, N1> >(noConst());
}


BZ_NAMESPACE_END

#endif // BZ_TINYVEC_CC
