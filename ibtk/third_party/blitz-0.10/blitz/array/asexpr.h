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
#ifndef BZ_ASEXPR_H
#define BZ_ASEXPR_H

#include <blitz/et-forward.h>
#include <blitz/numtrait.h>

BZ_NAMESPACE(blitz)

// The traits class asExpr converts arbitrary things to
// expression templatable operands.

//  Default to scalar.
template <typename T>
struct asExpr {
  typedef _bz_ArrayExpr<_bz_ArrayExprConstant<T> > T_expr;
  static T_expr getExpr(const T& x);
};

//  Already an expression template term
template <typename T>
struct asExpr<_bz_ArrayExpr<T> > {
    typedef _bz_ArrayExpr<T> T_expr;
  static const T_expr& getExpr(const T_expr& x);
};

//  Specialization of asExpr for array operands
// why doesn't it wrap iterators in an ArrayExpr?
template <typename T,int N>
struct asExpr<Array<T,N> > {
  //typedef FastArrayIterator<T,N> T_expr;
  typedef _bz_ArrayExpr<FastArrayIterator<T,N> > T_expr;
  static T_expr getExpr(const Array<T,N>& x);
};

//  Specialization of asExpr for tinyvector operands
template <typename T,int N>
struct asExpr<TinyVector<T,N> > {
  typedef _bz_ArrayExpr<FastTV2Iterator<T,N> > T_expr;
  static T_expr getExpr(const TinyVector<T,N>& x);
};

//  Specialization of asExpr for tinymatrix operands
template <typename T,int Nr, int Nc>
struct asExpr<TinyMatrix<T,Nr, Nc> > {
  typedef _bz_ArrayExpr<FastTM2Iterator<T,Nr, Nc> > T_expr;
  static T_expr getExpr(const TinyMatrix<T,Nr,Nc>& x);
};

//  Index placeholder
template <int N>
struct asExpr<IndexPlaceholder<N> > {
  //  typedef _bz_ArrayExpr<IndexPlaceholder<N> > T_expr;
  typedef _bz_ArrayExpr<IndexPlaceholder<N> > T_expr;
  static T_expr getExpr(const T_expr& x);
};

//  the levi-civita symbol
template <>
struct asExpr<LeviCivita> {
  typedef _bz_ArrayExpr<LeviCivita> T_expr;
  static T_expr getExpr(T_expr x);
};

//  Range
template <>
struct asExpr<Range> {
  typedef _bz_ArrayExpr<Range> T_expr;
  static T_expr getExpr(T_expr x);
};


// traits class that transforms ETBase subclasses into the
// ET<>-wrapped superclass and corresponding expression, but unlike
// the asExpr class it leaves POD types alone. This is necessary so
// operators on multicomponent arrays can resolve properly.
template<typename T>
struct asET {
  typedef T T_wrapped;
  typedef T T_expr;
};
template<typename T>
struct asET<ETBase<T> > {
  typedef ETBase<typename asExpr<T>::T_expr> T_wrapped;
  typedef typename asExpr<T>::T_expr T_expr;
};
template<typename T, int N>
struct asET<Array<T,N> > {
  typedef ETBase<typename asExpr<Array<T,N> >::T_expr> T_wrapped;
  typedef typename asExpr<Array<T,N> >::T_expr T_expr;
};
template<typename T, int N>
struct asET<TinyVector<T,N> > {
  typedef ETBase<typename asExpr<TinyVector<T,N> >::T_expr> T_wrapped;
  typedef typename asExpr<TinyVector<T,N> >::T_expr T_expr;
};
template<typename T, int Nr, int Nc>
struct asET<TinyMatrix<T,Nr,Nc> > {
  typedef ETBase<typename asExpr<TinyMatrix<T,Nr,Nc> >::T_expr> T_wrapped;
  typedef typename asExpr<TinyMatrix<T,Nr,Nc> >::T_expr T_expr;
};

// traits class that unwraps an ETBase type, otherwise leaves it untouched.
template<typename T>
struct unwrapET {
  typedef T T_unwrapped;
};
template<typename T>
struct unwrapET<ETBase<T> > {
  typedef T T_unwrapped;
};

// traits classes that are used to switch between an ET type or an
// unknown type. If the supplied type T is an ET type, T_selected will
// be T_ifET, otherwise T.
template<typename T, typename T_ifnotET, typename T_ifET>
struct selectET {
  typedef T_ifnotET T_selected;
};
template<typename T, typename T_ifnotET, typename T_ifET>
struct selectET<ETBase<T>, T_ifnotET, T_ifET> {
  typedef ETBase<T_ifET> T_selected;
};

// for binary exprs, it is more complicated. if T1 or T2 are an ET,
// T_ifET is selected, otherwise T_ifnotET.
template<typename T1, typename T2, typename T_ifnotET, typename T_ifET>
struct selectET2 {
  typedef T_ifnotET T_selected;
};
template<typename T1, typename T2, typename T_ifnotET, typename T_ifET>
struct selectET2<ETBase<T1>, T2, T_ifnotET, T_ifET> {
  typedef ETBase<T_ifET> T_selected;
};
template<typename T1, typename T2, typename T_ifnotET, typename T_ifET>
struct selectET2<T1, ETBase<T2>, T_ifnotET, T_ifET> {
  typedef ETBase<T_ifET> T_selected;
};
template<typename T1, typename T2, typename T_ifnotET, typename T_ifET>
struct selectET2<ETBase<T1>, ETBase<T2>, T_ifnotET, T_ifET> {
  typedef ETBase<T_ifET> T_selected;
};


// traits class that resolves to the ultimate numeric datatype used
// for operations on the container. This is necessary because for
// multicomponent containers we need to determine what the ultimate
// POD data type is.
template<typename T>
struct opType {
  typedef T T_optype;
};
template<typename T>
struct opType<ETBase<T> > {
  typedef typename opType<T>::T_optype T_optype;
};
template<typename T, int N>
struct opType<Array<T,N> > {
  typedef typename opType<T>::T_optype T_optype;
};
template<typename T, int N>
struct opType<TinyVector<T,N> > {
  typedef typename opType<T>::T_optype T_optype;
};
template<typename T, int Nr, int Nc>
struct opType<TinyMatrix<T,Nr,Nc> > {
  typedef typename opType<T>::T_optype T_optype;
};


#ifdef BZ_HAVE_TEMPLATES_AS_TEMPLATE_ARGUMENTS

//  traits classes that provide the return type of operations

template <template <typename T1> class OP, typename O1>
struct BzUnaryExprResult {
    typedef _bz_ArrayExpr<
      _bz_ArrayExprUnaryOp<
	typename asExpr<O1>::T_expr,
	OP<
	  typename asExpr<O1>::T_expr::T_optype
	  >
	> > T_result;
};

template <template <typename T1, typename T2> class OP,
          typename O1, typename O2>
struct BzBinaryExprResult {
    typedef _bz_ArrayExpr<_bz_ArrayExprBinaryOp<
        typename asExpr<O1>::T_expr,
        typename asExpr<O2>::T_expr,
        OP<
	  typename asExpr<O1>::T_expr::T_optype,
	  typename asExpr<O2>::T_expr::T_optype
	  > > > T_result;
};

template <template <typename T1, typename T2, typename T3> class OP,
          typename O1, typename O2, typename O3>
struct BzTernaryExprResult {
    typedef _bz_ArrayExpr<_bz_ArrayExprTernaryOp<
        typename asExpr<O1>::T_expr,
        typename asExpr<O2>::T_expr,
        typename asExpr<O3>::T_expr,
        OP<
	  typename asExpr<O1>::T_expr::T_optype,
	  typename asExpr<O2>::T_expr::T_optype,
	  typename asExpr<O3>::T_expr::T_optype
	  > > > T_result;
};

template <template <typename T1, typename T2, typename T3, typename T4> class OP,
          typename O1, typename O2, typename O3, typename O4>
struct BzQuaternaryExprResult {
  typedef _bz_ArrayExpr<_bz_ArrayExprQuaternaryOp<
			  typename asExpr<O1>::T_expr,
			  typename asExpr<O2>::T_expr,
			  typename asExpr<O3>::T_expr,
			  typename asExpr<O4>::T_expr,
			  OP<
			    typename asExpr<O1>::T_expr::T_optype,
			    typename asExpr<O2>::T_expr::T_optype,
			    typename asExpr<O3>::T_expr::T_optype,
			    typename asExpr<O4>::T_expr::T_optype
			    > > > T_result;
};

template <template <typename T1, typename T2> class RED, int N, typename O1,
	  typename P_result = BZ_SUMTYPE(typename asExpr<O1>::T_expr::T_optype)>
struct BzReductionResult {
  typedef _bz_ArrayExpr<
    _bz_ArrayExprReduce<
      typename asExpr<O1>::T_expr,
      N,
      RED<typename asExpr<O1>::T_expr::T_optype, P_result>
      > > T_result;
};

template<typename O1, int N0, int N1=0, int N2=0, int N3=0, int N4=0,
	 int N5=0, int N6=0, int N7=0, int N8=0, int N9=0, int N10=0>
struct BzIndexmapResult {
  typedef _bz_ArrayExpr<
    ArrayIndexMapping<
      typename asExpr<O1>::T_expr,
      N0, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10
      >
    > T_result;
};

template<template <typename T> class STENCIL, typename O1>
struct BzStencilResult {
  typedef _bz_ArrayExpr<
    STENCIL<
      typename asExpr<O1>::T_expr::T_range_result
      >
    > T_result;
};

template<template <typename T1, typename T2, typename T3> class STENCIL,
	 typename O1, typename O2, typename P_result>
struct BzBinaryStencilResult {
  typedef _bz_ArrayExpr<
    STENCIL<
      typename asExpr<O1>::T_expr::T_range_result,
      typename asExpr<O2>::T_expr::T_range_result,
      P_result
      > > T_result;
};


#endif /* BZ_HAVE_TEMPLATES_AS_TEMPLATE_ARGUMENTS */

BZ_NAMESPACE_END

#endif
