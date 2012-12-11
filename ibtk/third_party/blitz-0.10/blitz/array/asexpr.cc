// -*- C++ -*-
/***************************************************************************
 * blitz/array/asexpr.h  Declaration of the asExpr helper functions
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
#ifndef BZ_ASEXPR_CC
#define BZ_ASEXPR_cc

#include <blitz/array/asexpr.h>
#include <blitz/array-impl.h>
#include <blitz/levicivita.h>
#include <blitz/range.h>
#include <blitz/array/expr.h>

// implementations of asexpr functions.

BZ_NAMESPACE(blitz)

// default to scalar
template <typename T>
_bz_typename asExpr<T>::T_expr asExpr<T>::getExpr(const T& x)
{ return T_expr(x); }

// expression template term just returns itself
template <typename T>
const _bz_typename asExpr<_bz_ArrayExpr<T> >::T_expr&
asExpr<_bz_ArrayExpr<T> >::getExpr(const T_expr& x)
{ return x; }

// array operand returns iterator
template <typename T,int N>
_bz_typename asExpr<Array<T,N> >::T_expr
asExpr<Array<T,N> >::getExpr(const Array<T,N>& x)
{ return x.beginFast(); }

//  tinyvector operand returns iterator
template <typename T,int N>
_bz_typename asExpr<TinyVector<T,N> >::T_expr
asExpr<TinyVector<T,N> >::getExpr(const TinyVector<T,N>& x)
{ return x.beginFast(); }

//  tinymatrix operands returns iterator
template <typename T,int Nr, int Nc>
_bz_typename asExpr<TinyMatrix<T,Nr, Nc> >::T_expr
asExpr<TinyMatrix<T,Nr, Nc> >::getExpr(const TinyMatrix<T,Nr,Nc>& x)
{ return x.beginFast(); }

//  Index placeholder returns itself
template <int N>
_bz_typename asExpr<IndexPlaceholder<N> >::T_expr
asExpr<IndexPlaceholder<N> >::getExpr(const T_expr& x)
 { return x; }

//  the levi-civita symbol
inline asExpr<LeviCivita>::T_expr
asExpr<LeviCivita>::getExpr(T_expr x)
{ return T_expr(x); }

//  Range
inline asExpr<Range>::T_expr
asExpr<Range>::getExpr(T_expr x)
{ return T_expr(x); }

BZ_NAMESPACE_END

#endif
