// -*- C++ -*-
/***************************************************************************
 * blitz/array/expr.h     Array<T,N> expression templates
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
#ifndef BZ_EXPR_H
#define BZ_EXPR_H

#include <blitz/prettyprint.h>
#include <blitz/shapecheck.h>
#include <blitz/numinquire.h>
#include <blitz/array/domain.h>
#include <blitz/array/slice.h>
#include <blitz/bounds.h>

/*
 * The array expression templates iterator interface is followed by
 * these classes:
 *
 * FastArrayIterator          <blitz/array/fastiter.h>
 * _bz_ArrayExpr              <blitz/array/expr.h>
 * _bz_ArrayExprUnaryOp               "
 * _bz_ArrayExprBinaryOp              "
 * _bz_ArrayExprTernaryOp             "
 * _bz_ArrayExprConstant              "
 * ArrayIndexMapping          <blitz/array/map.h>
 * _bz_ArrayExprReduce        <blitz/array/reduce.h>
 * _bz_StencilExpr            <blitz/array/stencil-et.h>
 *   ... and derived types            "
 * IndexPlaceholder           <blitz/indexexpr.h>
 * _bz_ArrayWhere             <blitz/array/where.h>
 * _bz_FunctorExpr            <blitz/array/functorExpr.h>
 * _bz_FunctorExpr2                   "
 * _bz_FunctorExpr3                   "
 */

BZ_NAMESPACE(blitz)

#define BZ_MAX(a,b) (a)>(b) ? (a) : (b)
#define BZ_MIN(a,b) (a)<(b) ? (a) : (b)

template<typename T1, typename T2>
class _bz_ExprPair {
public:
    _bz_ExprPair(const T1& a, const T2& b)
        : first_(a), second_(b)
    { }

    const T1& first() const
    { return first_; }

    const T2& second() const
    { return second_; }

protected:
    T1 first_;
    T2 second_;
};

template<typename T1, typename T2>
inline _bz_ExprPair<T1,T2> makeExprPair(const T1& a, const T2& b)
{
    return _bz_ExprPair<T1,T2>(a,b);
}


/** The main expression template class which is used to wrap all other
    elements. This makes an expression easy to recognize.  A container
    that wants to be usable through the expression template machinery
    must must implement the same interface as this class. Once it
    does, it will be able to interoperate with the other containers.
    Note that the expression components (_bz_ArrayExprUnaryOp,
    _bz_ArrayExprBinaryOp, _bz_ArrayExprReduce, etc.) do not inherit
    from ETBase. Since all the functions are duplicated in all ET
    classes, they are documented only in this class. \todo Actually
    document the interface.  */
template<typename P_expr>
class _bz_ArrayExpr
#ifdef BZ_NEW_EXPRESSION_TEMPLATES
    : public ETBase<_bz_ArrayExpr<P_expr> >
#endif
{

public:
    typedef P_expr T_expr;
    typedef _bz_typename T_expr::T_numtype T_numtype;
  // select return type
  typedef typename unwrapET<typename T_expr::T_result>::T_unwrapped test;
  typedef typename selectET<typename T_expr::T_typeprop,
			    T_numtype,
			    _bz_ArrayExpr<test> >::T_selected T_typeprop;
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;

  typedef typename T_expr::T_optype T_optype;
    typedef T_expr T_ctorArg1;
    typedef int    T_ctorArg2;    // dummy
  typedef _bz_ArrayExpr<_bz_typename P_expr::T_range_result> T_range_result;

    static const int
        numArrayOperands = T_expr::numArrayOperands,
        numTVOperands = T_expr::numTVOperands,
        numTMOperands = T_expr::numTMOperands,
        numIndexPlaceholders = T_expr::numIndexPlaceholders,
      minWidth = T_expr::minWidth,
      maxWidth = T_expr::maxWidth,
        rank_ = T_expr::rank_;

  // traits class that computes the vectorized return type for vector
  // width N.
  template<int N> struct tvresult {
    typedef _bz_ArrayExpr<typename T_expr::template tvresult<N>::Type> Type;
  };

    _bz_ArrayExpr(const _bz_ArrayExpr<T_expr>& a)
#ifdef BZ_NEW_EXPRESSION_TEMPLATES
        : ETBase< _bz_ArrayExpr<T_expr> >(a), iter_(a.iter_)
#else
        : iter_(a.iter_)
#endif
    { }

#if defined(BZ_NEW_EXPRESSION_TEMPLATES) && ! defined(__MWERKS__)
    template<typename T>
    _bz_ArrayExpr(const T& a)
        : iter_(a)
    { }
#else

    _bz_ArrayExpr(BZ_ETPARM(T_expr) a)
        : iter_(a)
    { }
#if !defined(__MWERKS__)
    _bz_ArrayExpr(BZ_ETPARM(_bz_typename T_expr::T_ctorArg1) a)
        : iter_(a)
    { }
#endif
#endif

    template<typename T1, typename T2>
    _bz_ArrayExpr(BZ_ETPARM(T1) a, BZ_ETPARM(T2) b)
        : iter_(a, b)
    { }

    template<typename T1, typename T2, typename T3>
    _bz_ArrayExpr(BZ_ETPARM(T1) a, BZ_ETPARM(T2) b, BZ_ETPARM(T3) c)
        : iter_(a, b, c)
    { }

    template<typename T1, typename T2, typename T3, typename T4>
    _bz_ArrayExpr(BZ_ETPARM(T1) a, BZ_ETPARM(T2) b, BZ_ETPARM(T3) c,
        BZ_ETPARM(T4) d) : iter_(a, b, c, d)
    { }

    template<typename T1, typename T2>
    _bz_ArrayExpr(const _bz_ExprPair<T1,T2>& exprpair)
        : iter_(exprpair.first(), exprpair.second())
    { }

  T_result operator*() const { return *iter_; }

  T_result first_value() const { return iter_.first_value(); }

#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
    template<int N_rank>
    T_result operator()(const TinyVector<int, N_rank> i) const { return iter_(i); }
#else
    template<int N_rank>
    T_result operator()(const TinyVector<int, N_rank>& i) const { return iter_(i); }
#endif

  template<int N>
  T_range_result operator()(const RectDomain<N>& d) const
  {
    return T_range_result(iter_(d));
  }

    int ascending(const int rank) const { return iter_.ascending(rank); }
    int ordering(const int rank)  const { return iter_.ordering(rank);  }
    int lbound(const int rank)    const { return iter_.lbound(rank);    }
    int ubound(const int rank)    const { return iter_.ubound(rank);    }
    RectDomain<rank_> domain() const { return iter_.domain(); }

    void push(int position)
    { iter_.push(position); }

    void pop(int position)
    { iter_.pop(position); }

    void advance()
    { iter_.advance(); }

    void advance(int n)
    { iter_.advance(n); }

    void loadStride(int rank)
    { iter_.loadStride(rank); }

    bool isUnitStride(int rank) const
    { return iter_.isUnitStride(rank); }

    bool isUnitStride() const
    { return iter_.isUnitStride(); }

    void advanceUnitStride()
    { iter_.advanceUnitStride(); }

    bool canCollapse(int outerLoopRank, int innerLoopRank) const
    {
        // BZ_DEBUG_MESSAGE("_bz_ArrayExpr<>::canCollapse()");
        return iter_.canCollapse(outerLoopRank, innerLoopRank);
    }

    T_result operator[](int i) const
    { return iter_[i]; }

    T_result fastRead(diffType i) const
    { return iter_.fastRead(i); }

  template<int N>
  typename tvresult<N>::Type fastRead_tv(diffType i) const
  { return iter_.template fastRead_tv<N>(i); }

  bool isVectorAligned(diffType offset) const
  { return iter_.isVectorAligned(offset); }

    // this is needed for the stencil expression fastRead to work
    void _bz_offsetData(sizeType i)
    { iter_._bz_offsetData(i); }

    // and these are needed for stencil expression shift to work
    void _bz_offsetData(sizeType offset, int dim)
    { iter_._bz_offsetData(offset, dim);}

    void _bz_offsetData(sizeType offset1, int dim1, sizeType offset2, int dim2)
    { iter_._bz_offsetData(offset1, dim1, offset2, dim2);}

    diffType suggestStride(int rank) const
    { return iter_.suggestStride(rank); }

    bool isStride(int rank, diffType stride) const
    { return iter_.isStride(rank,stride); }

    void prettyPrint(BZ_STD_SCOPE(string) &str) const
    {
        prettyPrintFormat format(true);  // Terse formatting by default
        iter_.prettyPrint(str, format);
    }

    void prettyPrint(BZ_STD_SCOPE(string) &str,
        prettyPrintFormat& format) const
    { iter_.prettyPrint(str, format); }

    template<typename T_shape>
    bool shapeCheck(const T_shape& shape) const
    { return iter_.shapeCheck(shape); }

    template<int N>
    void moveTo(const TinyVector<int, N>& i)
    {
        iter_.moveTo(i);
    }

    T_result shift(int offset, int dim) const
    {
      return iter_.shift(offset, dim);
    }

    T_result shift(int offset1, int dim1,int offset2, int dim2) const
    {
      return iter_.shift(offset1, dim1, offset2, dim2);
    }

  // sliceinfo for expressions
  template<typename T1, typename T2 = nilArraySection,
	   class T3 = nilArraySection, typename T4 = nilArraySection,
	   class T5 = nilArraySection, typename T6 = nilArraySection,
	   class T7 = nilArraySection, typename T8 = nilArraySection,
	   class T9 = nilArraySection, typename T10 = nilArraySection,
	   class T11 = nilArraySection>
  class SliceInfo {
  public:
    typedef typename T_expr::template SliceInfo<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>::T_slice T_subexpr;
    typedef _bz_ArrayExpr<T_subexpr> T_slice;
};

  // slicing (experimental)

  // because _bz_ArrayExpr is the "top-level" expression class, it has
  // the burden of supplying enough variants of the operator to make
  // it user-friendly. Hence the large numbers that follow
  template<typename T1>
  typename SliceInfo<T1>::T_slice
  operator()(T1 r1) const
  {
    return typename SliceInfo<T1>::T_slice
      (iter_
       (r1,
	nilArraySection(), nilArraySection(), nilArraySection(), nilArraySection(),
	nilArraySection(), nilArraySection(), nilArraySection(),
	nilArraySection(), nilArraySection(), nilArraySection()));
  }

  template<typename T1, typename T2>
  typename SliceInfo<T1,T2>::T_slice
  operator()(T1 r1, T2 r2) const
  {
    typedef typename SliceInfo<T1,T2>::T_slice slice;
    return slice(iter_
		 (r1, r2, nilArraySection(), nilArraySection(), nilArraySection(),
		  nilArraySection(), nilArraySection(), nilArraySection(),
		  nilArraySection(), nilArraySection(), nilArraySection()));
  }


    template<typename T1, typename T2, typename T3>
    typename SliceInfo<T1,T2,T3>::T_slice
    operator()(T1 r1, T2 r2, T3 r3) const
    {
        typedef typename SliceInfo<T1,T2,T3>::T_slice slice;
        return slice(iter_(r1, r2, r3, nilArraySection(), nilArraySection(),
			   nilArraySection(), nilArraySection(), nilArraySection(),
			   nilArraySection(), nilArraySection(), nilArraySection()));
    }

    template<typename T1, typename T2, typename T3, typename T4>
    typename SliceInfo<T1,T2,T3,T4>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4) const
    {
        typedef typename SliceInfo<T1,T2,T3,T4>::T_slice slice;
        return slice(iter_(r1, r2, r3, r4, nilArraySection(), nilArraySection(),
			   nilArraySection(), nilArraySection(), nilArraySection(),
			   nilArraySection(), nilArraySection()));
    }

    template<typename T1, typename T2, typename T3, typename T4, typename T5>
    typename SliceInfo<T1,T2,T3,T4,T5>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5) const
    {
        typedef typename SliceInfo<T1,T2,T3,T4,T5>::T_slice slice;
        return slice(iter_(r1, r2, r3, r4, r5, nilArraySection(),
			   nilArraySection(), nilArraySection(), nilArraySection(),
			   nilArraySection(), nilArraySection()));
    }

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
    typename SliceInfo<T1,T2,T3,T4,T5,T6>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6) const
    {
        typedef typename SliceInfo<T1,T2,T3,T4,T5,T6>::T_slice slice;
        return slice(iter_(r1, r2, r3, r4, r5, r6, nilArraySection(),
			   nilArraySection(), nilArraySection(),
			   nilArraySection(), nilArraySection()));
    }

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
        typename T7>
    typename SliceInfo<T1,T2,T3,T4,T5,T6,T7>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7) const
    {
        typedef typename SliceInfo<T1,T2,T3,T4,T5,T6,T7>::T_slice slice;
        return slice(iter_(r1, r2, r3, r4, r5, r6, r7,
			   nilArraySection(), nilArraySection(),
			   nilArraySection(), nilArraySection()));
    }

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
        typename T7, typename T8>
    typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8) const
    {
        typedef typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8>::T_slice slice;
        return slice(iter_(r1, r2, r3, r4, r5, r6, r7, r8,
			   nilArraySection(), nilArraySection(), nilArraySection()));
    }

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
        typename T7, typename T8, typename T9>
    typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9) const
    {
        typedef typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9>::T_slice slice;
        return slice(iter_(r1, r2, r3, r4, r5, r6, r7, r8, r9,
			   nilArraySection(), nilArraySection()));
    }

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
        typename T7, typename T8, typename T9, typename T10>
    typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10) const
    {
        typedef typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice slice;
        return slice(iter_(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10,
			   nilArraySection()));
    }

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
        typename T7, typename T8, typename T9, typename T10, typename T11>
    typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const
    {
        typedef typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::T_slice slice;
        return slice(iter_(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11));
    }

  // now complete slicings to a scalar, which can't be expressed with the above
  // unlike for arrays, these expressions are rvalues so we return by value
    T_numtype operator()(int i0) const
    {
	return iter_(TinyVector<int, 1>(i0));
    }

  T_numtype operator()(int i0, int i1) const
    {
	return iter_(TinyVector<int, 2>(i0, i1));
    }

  T_numtype operator()(int i0, int i1, int i2) const
    {
      return iter_(TinyVector<int, 3>(i0, i1, i2));
    }

  T_numtype operator()(int i0, int i1, int i2, int i3) const
    {
      return iter_(TinyVector<int, 4>(i0, i1, i2, i3));
    }

  T_numtype operator()(int i0, int i1, int i2, int i3, int i4) const
    {
      return iter_(TinyVector<int, 5>(i0, i1, i2, i3, i4));
    }

  T_numtype operator()(int i0, int i1, int i2, int i3, int i4, int i5) const
    {
      return iter_(TinyVector<int, 6>(i0, i1, i2, i3, i4, i5));
    }

  T_numtype operator()(int i0, int i1, int i2, int i3, int i4, int i5,
		       int i6) const
    {
      return iter_(TinyVector<int, 7>(i0, i1, i2, i3, i4, i5, i6));
    }

  T_numtype operator()(int i0, int i1, int i2, int i3, int i4, int i5,
		       int i6, int i7) const
    {
      return iter_(TinyVector<int, 8>(i0, i1, i2, i3, i4, i5, i6, i7));
    }

  T_numtype operator()(int i0, int i1, int i2, int i3, int i4, int i5,
		       int i6, int i7, int i8) const
    {
      return iter_(TinyVector<int, 9>(i0, i1, i2, i3, i4, i5, i6, i7, i8));
    }

  T_numtype operator()(int i0, int i1, int i2, int i3, int i4, int i5,
		       int i6, int i7, int i8, int i9) const
    {
      return iter_(TinyVector<int, 10>(i0, i1, i2, i3, i4, i5, i6, i7, i8, i9));
    }

  T_numtype operator()(int i0, int i1, int i2, int i3, int i4, int i5,
		       int i6, int i7, int i8, int i9, int i10) const
    {
      return iter_(TinyVector<int, 11>(i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10));
    }

protected:
    _bz_ArrayExpr() { }

    T_expr iter_;
};


template<typename P_expr, typename P_op>
class _bz_ArrayExprUnaryOp {
public:
    typedef P_expr T_expr;
    typedef P_op T_op;
    typedef _bz_typename T_expr::T_numtype T_numtype1;
    typedef _bz_typename T_op::T_numtype T_numtype;

  // select return type
  typedef typename unwrapET<typename T_expr::T_result>::T_unwrapped test;
  typedef typename selectET<typename T_expr::T_typeprop,
			    T_numtype,
			    _bz_ArrayExprUnaryOp<test, T_op> >::T_selected T_typeprop;
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;

  typedef T_numtype T_optype;

    typedef T_expr T_ctorArg1;
    typedef int    T_ctorArg2;    // dummy
  typedef _bz_ArrayExprUnaryOp<_bz_typename P_expr::T_range_result,
			       P_op> T_range_result;

    static const int
        numArrayOperands = T_expr::numArrayOperands,
        numTVOperands = T_expr::numTVOperands,
        numTMOperands = T_expr::numTMOperands,
        numIndexPlaceholders = T_expr::numIndexPlaceholders,
      minWidth = T_expr::minWidth,
      maxWidth = T_expr::maxWidth,
        rank_ = T_expr::rank_;

  template<int N> struct tvresult {
    typedef _bz_ArrayExprUnaryOp<
      typename T_expr::template tvresult<N>::Type,
      T_op> Type;
  };

    _bz_ArrayExprUnaryOp(const _bz_ArrayExprUnaryOp<T_expr, T_op>& a)
        : iter_(a.iter_)
    { }

    _bz_ArrayExprUnaryOp(BZ_ETPARM(T_expr) a)
        : iter_(a)
    { }

  /*
    _bz_ArrayExprUnaryOp(_bz_typename T_expr::T_ctorArg1 a)
        : iter_(a)
    { }
  */

#if BZ_TEMPLATE_CTOR_DOESNT_CAUSE_HAVOC
    template<typename T1>
    explicit _bz_ArrayExprUnaryOp(BZ_ETPARM(T1) a)
        : iter_(a)
    { }
#endif

    int ascending(const int rank) const { return iter_.ascending(rank); }
    int ordering(const int rank)  const { return iter_.ordering(rank);  }
    int lbound(const int rank)    const { return iter_.lbound(rank);    }
    int ubound(const int rank)    const { return iter_.ubound(rank);    }
    RectDomain<rank_> domain() const { return iter_.domain(); }


  /* Functions for reading data. Because they must depend on the
   * result type, they utilize a helper class.
   */

  // For numtypes, apply operator
  template<typename T> struct readHelper {
    static T_result fastRead(const T_expr& iter, diffType i) {
      return (T_op::apply(iter.fastRead(i))); };
    static T_result indexop(const T_expr& iter, int i) {
      return (T_op::apply(iter[i])); };
    static T_result deref(const T_expr& iter) {
      return T_op::apply(*iter); }
    static T_result first_value(const T_expr& iter)  {
      return T_op::apply(iter.first_value()); }
    static T_result shift(const T_expr& iter,
			  int offset, int dim) {
      return T_op::apply(iter.shift(offset, dim)); }
    static T_result shift(const T_expr& iter,
			  int offset1, int dim1, int offset2, int dim2) {
      return T_op::apply(iter.shift(offset1, dim1, offset2, dim2)); }
    template<int N_rank>
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
    static T_result indexop(const T_expr& iter,
			    const TinyVector<int, N_rank> i) {
#else
      static T_result indexop(const T_expr& iter,
			      const TinyVector<int, N_rank>& i) {
#endif
      return T_op::apply(iter(i)); }
  };

  // For ET types, bypass operator and create expression
    template<typename T> struct readHelper<ETBase<T> > {
      static T_result fastRead(const T_expr& iter, diffType i) {
	return iter.fastRead(i); };
      static T_result indexop(const T_expr& iter, int i) {
	return iter[i]; };
      static T_result deref(const T_expr& iter) {
	return *iter; }
    static T_result first_value(const T_expr& iter)  {
      return iter.first_value(); }
    static T_result shift(const T_expr& iter,
			  int offset, int dim) {
      return T_result(iter.shift(offset, dim)); }
    static T_result shift(const T_expr& iter,
			  int offset1, int dim1, int offset2, int dim2) {
      return T_result(iter.shift(offset1, dim1, offset2, dim2)); }
      template<int N_rank>
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
      static T_result indexop(const T_expr& iter,
			      const TinyVector<int, N_rank> i) {
#else
    static T_result indexop(const T_expr& iter,
			    const TinyVector<int, N_rank>& i) {
#endif
      return iter(i); }
    };

    T_result fastRead(diffType i) const {
      return readHelper<T_typeprop>::fastRead(iter_, i); }

  template<int N>
  typename tvresult<N>::Type fastRead_tv(diffType i) const
  { return iter_.template fastRead_tv<N>(i); }


    T_result operator[](int i) const {
      return readHelper<T_typeprop>::indexop(iter_, i); }

    template<int N_rank>
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
    T_result operator()(const TinyVector<int, N_rank> i) const {
#else
      T_result operator()(const TinyVector<int, N_rank>& i) const {
#endif
      return readHelper<T_typeprop>::indexop(iter_,i); }

      T_result operator*() const {
	return readHelper<T_typeprop>::deref(iter_); }

      T_result first_value() const {
	return readHelper<T_typeprop>::first_value(iter_); }

    T_result shift(int offset, int dim) const
    {
      return readHelper<T_typeprop>::shift(iter_, offset, dim);
    }

    T_result shift(int offset1, int dim1,int offset2, int dim2) const
    {
      return readHelper<T_typeprop>::shift(iter_,
					   offset1, dim1, offset2, dim2);
    }

      // ****** end reading

  bool isVectorAligned(diffType offset) const
  { return iter_.isVectorAligned(offset); }

  template<int N>
  T_range_result operator()(const RectDomain<N>& d) const
  {
    return T_range_result(iter_(d));
  }

    void push(int position)
    {
        iter_.push(position);
    }

    void pop(int position)
    {
        iter_.pop(position);
    }

    void advance()
    {
        iter_.advance();
    }

    void advance(int n)
    {
        iter_.advance(n);
    }

    void loadStride(int rank)
    {
        iter_.loadStride(rank);
    }

    bool isUnitStride(int rank) const
    { return iter_.isUnitStride(rank); }

    bool isUnitStride() const
    { return iter_.isUnitStride(); }

    void advanceUnitStride()
    {
        iter_.advanceUnitStride();
    }

  template<int N>
  void moveTo(const TinyVector<int,N>& i)
    {
        iter_.moveTo(i);
    }

    bool canCollapse(int outerLoopRank, int innerLoopRank) const
    {
        // BZ_DEBUG_MESSAGE("_bz_ArrayExprUnaryOp<>::canCollapse");
        return iter_.canCollapse(outerLoopRank, innerLoopRank);
    }

  // this is needed for the stencil expression fastRead to work
  void _bz_offsetData(sizeType i)
  {
    iter_._bz_offsetData(i);
  }

    // and these are needed for stencil expression shift to work
    void _bz_offsetData(sizeType offset, int dim)
    { iter_._bz_offsetData(offset, dim);}

    void _bz_offsetData(sizeType offset1, int dim1, sizeType offset2, int dim2)
    { iter_._bz_offsetData(offset1, dim1, offset2, dim2);}

    diffType suggestStride(int rank) const
    { return iter_.suggestStride(rank); }

    bool isStride(int rank, diffType stride) const
    { return iter_.isStride(rank,stride); }

    void prettyPrint(BZ_STD_SCOPE(string) &str,
        prettyPrintFormat& format) const
    { T_op::prettyPrint(str, format, iter_); }

    template<typename T_shape>
    bool shapeCheck(const T_shape& shape) const
    { return iter_.shapeCheck(shape); }


  // sliceinfo for expressions
  template<typename T1, typename T2 = nilArraySection,
	   class T3 = nilArraySection, typename T4 = nilArraySection,
	   class T5 = nilArraySection, typename T6 = nilArraySection,
	   class T7 = nilArraySection, typename T8 = nilArraySection,
	   class T9 = nilArraySection, typename T10 = nilArraySection,
	   class T11 = nilArraySection>
  class SliceInfo {
  public:
    typedef typename T_expr::template SliceInfo<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>::T_slice T_slice1;
    typedef _bz_ArrayExprUnaryOp<T_slice1, T_op> T_slice;
};

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
        typename T7, typename T8, typename T9, typename T10, typename T11>
    typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const
    {
      return typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::T_slice
	(iter_(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11));
    }

protected:
    _bz_ArrayExprUnaryOp() { }

    T_expr iter_;
};


template<typename P_expr1, typename P_expr2, typename P_op>
class _bz_ArrayExprBinaryOp {
public:
    typedef P_expr1 T_expr1;
    typedef P_expr2 T_expr2;
    typedef P_op T_op;
    typedef _bz_typename T_expr1::T_numtype T_numtype1;
    typedef _bz_typename T_expr2::T_numtype T_numtype2;
    typedef _bz_typename T_op::T_numtype T_numtype;

  // select return type
  typedef typename unwrapET<typename T_expr1::T_result>::T_unwrapped T_unwrapped1;
  typedef typename unwrapET<typename T_expr2::T_result>::T_unwrapped T_unwrapped2;
  typedef typename selectET2<typename T_expr1::T_typeprop,
			     typename T_expr2::T_typeprop,
			     T_numtype,
			     _bz_ArrayExprBinaryOp<typename asExpr<T_unwrapped1>::T_expr,
						   typename asExpr<T_unwrapped2>::T_expr,
						   T_op> >::T_selected T_typeprop;
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;
  typedef typename T_op::T_numtype T_optype;

    typedef T_expr1 T_ctorArg1;
    typedef T_expr2 T_ctorArg2;
  typedef _bz_ArrayExprBinaryOp<_bz_typename P_expr1::T_range_result,
				_bz_typename P_expr2::T_range_result,
				P_op> T_range_result;

    static const int
        numArrayOperands = T_expr1::numArrayOperands
                         + T_expr2::numArrayOperands,
        numTVOperands = T_expr1::numTVOperands +
      T_expr2::numTVOperands,
        numTMOperands = T_expr1::numTMOperands +
      T_expr2::numTMOperands,
        numIndexPlaceholders = T_expr1::numIndexPlaceholders
                             + T_expr2::numIndexPlaceholders,
      minWidth = BZ_MIN(T_expr1::minWidth, T_expr2::minWidth),
      maxWidth = BZ_MAX(T_expr1::maxWidth, T_expr2::maxWidth),
      rank_ = BZ_MAX(T_expr1::rank_, T_expr2::rank_);

  template<int N> struct tvresult {
    typedef _bz_ArrayExprBinaryOp<
      typename T_expr1::template tvresult<N>::Type,
      typename T_expr2::template tvresult<N>::Type,
      T_op> Type;
  };

    _bz_ArrayExprBinaryOp(
        const _bz_ArrayExprBinaryOp<T_expr1, T_expr2, T_op>& a)
        : iter1_(a.iter1_), iter2_(a.iter2_)
    { }

    template<typename T1, typename T2>
    _bz_ArrayExprBinaryOp(BZ_ETPARM(T1) a, BZ_ETPARM(T2) b)
        : iter1_(a), iter2_(b)
    { }


  /* Functions for reading. Because they must depend on the result
   * type, they utilize a helper class.
   */

  // For numtypes, apply operator
  template<typename T> struct readHelper {
    static T_result fastRead(const T_expr1& iter1, const T_expr2& iter2,
			     diffType i) {
      return T_op::apply(iter1.fastRead(i), iter2.fastRead(i)); }
    static T_result indexop(const T_expr1& iter1, const T_expr2& iter2, int i) {
      return T_op::apply(iter1[i], iter2[i]); };
    static T_result deref(const T_expr1& iter1, const T_expr2& iter2) {
      return T_op::apply(*iter1, *iter2); }
    static T_result first_value(const T_expr1& iter1, const T_expr2& iter2)  {
      return T_op::apply(iter1.first_value(), iter2.first_value()); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  int offset, int dim) {
      return T_op::apply(iter1.shift(offset, dim),iter2.shift(offset, dim)); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  int offset1, int dim1, int offset2, int dim2) {
      return T_op::apply(iter1.shift(offset1, dim1, offset2, dim2),
			 iter2.shift(offset1, dim1, offset2, dim2)); }
    template<int N_rank>
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
    static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			    const TinyVector<int, N_rank> i) {
#else
      static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			      const TinyVector<int, N_rank>& i) {
#endif
	return T_op::apply(iter1(i), iter2(i)); }
    };

  // For ET types, bypass operator and create expression
    template<typename T> struct readHelper<ETBase<T> > {
      static T_result fastRead(const T_expr1& iter1, const T_expr2& iter2,
			       diffType i) {
	return T_result(iter1.fastRead(i), iter2.fastRead(i)); }
    static T_result indexop(const T_expr1& iter1, const T_expr2& iter2, int i) {
      return T_result(iter1[i], iter2[i]); };
    static T_result deref(const T_expr1& iter1, const T_expr2& iter2) {
      return T_result(*iter1, *iter2); }
    static T_result first_value(const T_expr1& iter1, const T_expr2& iter2)  {
      return T_result(iter1.first_value(), iter2.first_value()); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  int offset, int dim) {
      return T_result(iter1.shift(offset, dim),iter2.shift(offset, dim)); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  int offset1, int dim1, int offset2, int dim2) {
      return T_result(iter1.shift(offset1, dim1, offset2, dim2),
		      iter2.shift(offset1, dim1, offset2, dim2)); }
      template<int N_rank>
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
      static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			      const TinyVector<int, N_rank> i) {
#else
	static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
				const TinyVector<int, N_rank>& i) {
#endif
	  return T_result(iter1(i), iter2(i)); }
      };

    T_result fastRead(diffType i) const {
      return readHelper<T_typeprop>::fastRead(iter1_, iter2_, i); }

  template<int N>
  typename tvresult<N>::Type fastRead_tv(diffType i) const
      { return typename tvresult<N>::Type(iter1_.template fastRead_tv<N>(i),
					  iter2_.template fastRead_tv<N>(i)); }

    T_result operator[](int i) const {
      return readHelper<T_typeprop>::indexop(iter1_, iter2_, i); }

    template<int N_rank>
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
    T_result operator()(const TinyVector<int, N_rank> i) const {
#else
      T_result operator()(const TinyVector<int, N_rank>& i) const {
#endif
	return readHelper<T_typeprop>::indexop(iter1_, iter2_, i); }

      T_result operator*() const {
	return readHelper<T_typeprop>::deref(iter1_, iter2_); }

      T_result first_value() const {
	return readHelper<T_typeprop>::first_value(iter1_, iter2_); }

    T_result shift(int offset, int dim) const {
      return readHelper<T_typeprop>::shift(iter1_, iter2_, offset, dim); }

    T_result shift(int offset1, int dim1,int offset2, int dim2) const {
      return readHelper<T_typeprop>::shift(iter1_, iter2_,
					   offset1, dim1, offset2, dim2); }


      // ****** end reading

  bool isVectorAligned(diffType offset) const
  { return iter1_.isVectorAligned(offset) &&
      iter2_.isVectorAligned(offset); }

  template<int N>
  T_range_result operator()(const RectDomain<N>& d) const
  {
    return T_range_result(iter1_(d), iter2_(d));
  }

    int ascending(const int rank) const {
        return bounds::compute_ascending(rank, iter1_.ascending(rank), iter2_.ascending(rank));
    }

    int ordering(const int rank) const {
        return bounds::compute_ordering(rank, iter1_.ordering(rank), iter2_.ordering(rank));
    }

    int lbound(const int rank) const {
        return bounds::compute_lbound(rank, iter1_.lbound(rank), iter2_.lbound(rank));
    }

    int ubound(const int rank) const {
        return bounds::compute_ubound(rank, iter1_.ubound(rank), iter2_.ubound(rank));
    }

  // defer calculation to lbound/ubound
  RectDomain<rank_> domain() const
  {
    TinyVector<int, rank_> lb, ub;
    for(int r=0; r<rank_; ++r) {
      lb[r]=lbound(r); ub[r]=ubound(r);
    }
    return RectDomain<rank_>(lb,ub);
  }

    void push(int position)
    {
        iter1_.push(position);
        iter2_.push(position);
    }

    void pop(int position)
    {
        iter1_.pop(position);
        iter2_.pop(position);
    }

    void advance()
    {
        iter1_.advance();
        iter2_.advance();
    }

    void advance(int n)
    {
        iter1_.advance(n);
        iter2_.advance(n);
    }

    void loadStride(int rank)
    {
        iter1_.loadStride(rank);
        iter2_.loadStride(rank);
    }

    bool isUnitStride(int rank) const
    { return iter1_.isUnitStride(rank) && iter2_.isUnitStride(rank); }

    bool isUnitStride() const
    { return iter1_.isUnitStride() && iter2_.isUnitStride(); }

    void advanceUnitStride()
    {
        iter1_.advanceUnitStride();
        iter2_.advanceUnitStride();
    }

    bool canCollapse(int outerLoopRank, int innerLoopRank) const
    {
        // BZ_DEBUG_MESSAGE("_bz_ArrayExprBinaryOp<>::canCollapse");
        return iter1_.canCollapse(outerLoopRank, innerLoopRank)
            && iter2_.canCollapse(outerLoopRank, innerLoopRank);
    }

    // this is needed for the stencil expression fastRead to work
    void _bz_offsetData(sizeType i)
  {
    iter1_._bz_offsetData(i);
    iter2_._bz_offsetData(i);
  }

    // and these are needed for stencil expression shift to work
    void _bz_offsetData(sizeType offset, int dim)
    {
      iter1_._bz_offsetData(offset, dim);
      iter2_._bz_offsetData(offset, dim);
    }

    void _bz_offsetData(sizeType offset1, int dim1, sizeType offset2, int dim2)
    {
      iter1_._bz_offsetData(offset1, dim1, offset2, dim2);
      iter2_._bz_offsetData(offset1, dim1, offset2, dim2);
    }

    diffType suggestStride(int rank) const
    {
        diffType stride1 = iter1_.suggestStride(rank);
        diffType stride2 = iter2_.suggestStride(rank);
        return (stride1 > stride2) ? stride1 : stride2;
    }

    bool isStride(int rank, diffType stride) const
    {
        return iter1_.isStride(rank,stride) && iter2_.isStride(rank,stride);
    }

  template<int N>
  void moveTo(const TinyVector<int,N>& i)
    {
        iter1_.moveTo(i);
        iter2_.moveTo(i);
    }

    void prettyPrint(BZ_STD_SCOPE(string) &str,
        prettyPrintFormat& format) const
    {
        T_op::prettyPrint(str, format, iter1_, iter2_);
    }

    template<typename T_shape>
    bool shapeCheck(const T_shape& shape) const
    { return iter1_.shapeCheck(shape) && iter2_.shapeCheck(shape); }


  // sliceinfo for expressions
  template<typename T1, typename T2 = nilArraySection,
	   class T3 = nilArraySection, typename T4 = nilArraySection,
	   class T5 = nilArraySection, typename T6 = nilArraySection,
	   class T7 = nilArraySection, typename T8 = nilArraySection,
	   class T9 = nilArraySection, typename T10 = nilArraySection,
	   class T11 = nilArraySection>
  class SliceInfo {
  public:
    typedef typename T_expr1::template SliceInfo<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>::T_slice T_slice1;
    typedef typename T_expr2::template SliceInfo<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>::T_slice T_slice2;
    typedef _bz_ArrayExprBinaryOp<T_slice1, T_slice2, T_op> T_slice;
};

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
        typename T7, typename T8, typename T9, typename T10, typename T11>
    typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const
    {
      return typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::T_slice
	(iter1_(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11),
	 iter2_(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11));
    }

protected:
    _bz_ArrayExprBinaryOp() { }

    T_expr1 iter1_;
    T_expr2 iter2_;
};

template<typename P_expr1, typename P_expr2, typename P_expr3, typename P_op>
class _bz_ArrayExprTernaryOp {
public:
    typedef P_expr1 T_expr1;
    typedef P_expr2 T_expr2;
    typedef P_expr3 T_expr3;
    typedef P_op T_op;
    typedef _bz_typename T_expr1::T_numtype T_numtype1;
    typedef _bz_typename T_expr2::T_numtype T_numtype2;
    typedef _bz_typename T_expr3::T_numtype T_numtype3;
    typedef _bz_typename T_op::T_numtype T_numtype;

  // select return type
  typedef typename unwrapET<
    typename T_expr1::T_result>::T_unwrapped T_unwrapped1;
  typedef typename unwrapET<
    typename T_expr2::T_result>::T_unwrapped T_unwrapped2;
  typedef typename unwrapET<
    typename T_expr3::T_result>::T_unwrapped T_unwrapped3;
  typedef typename selectET2<typename T_expr1::T_typeprop,
			     typename T_expr2::T_typeprop,
			     T_numtype,
			     char>::T_selected T_intermediary;

  typedef typename selectET2<
    T_intermediary,
    typename T_expr3::T_typeprop,
    T_numtype,
    _bz_ArrayExprTernaryOp<typename asExpr<T_unwrapped1>::T_expr,
			   typename asExpr<T_unwrapped2>::T_expr,
			   typename asExpr<T_unwrapped3>::T_expr,
			   T_op> >::T_selected T_typeprop;
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;
  typedef typename T_op::T_numtype T_optype;

    typedef T_expr1 T_ctorArg1;
    typedef T_expr2 T_ctorArg2;
    typedef T_expr3 T_ctorArg3;
  typedef _bz_ArrayExprTernaryOp<_bz_typename P_expr1::T_range_result,
				 _bz_typename P_expr2::T_range_result,
				 _bz_typename P_expr3::T_range_result, P_op> T_range_result;

    static const int
        numArrayOperands = T_expr1::numArrayOperands
                         + T_expr2::numArrayOperands
                         + T_expr3::numArrayOperands,
        numTVOperands = T_expr1::numTVOperands +
      T_expr2::numTVOperands +
      T_expr3::numTVOperands,
        numTMOperands = T_expr1::numTMOperands +
      T_expr2::numTMOperands +
      T_expr3::numTMOperands,
        numIndexPlaceholders = T_expr1::numIndexPlaceholders
                             + T_expr2::numIndexPlaceholders
                             + T_expr3::numIndexPlaceholders,
      minWidth = BZ_MIN(BZ_MIN(T_expr1::minWidth, T_expr2::minWidth),
			T_expr3::minWidth),
      maxWidth = BZ_MAX(BZ_MAX(T_expr1::maxWidth, T_expr2::maxWidth),
			T_expr3::maxWidth),
      rank_ = BZ_MAX(BZ_MAX(T_expr1::rank_, T_expr2::rank_),
		     T_expr3::rank_);

  template<int N> struct tvresult {
    typedef _bz_ArrayExprTernaryOp<
      typename T_expr1::template tvresult<N>::Type,
      typename T_expr2::template tvresult<N>::Type,
      typename T_expr3::template tvresult<N>::Type,
      T_op> Type;
  };

    _bz_ArrayExprTernaryOp(
        const _bz_ArrayExprTernaryOp<T_expr1, T_expr2, T_expr3, T_op>& a)
        : iter1_(a.iter1_), iter2_(a.iter2_), iter3_(a.iter3_)
    { }

    template<typename T1, typename T2, typename T3>
    _bz_ArrayExprTernaryOp(BZ_ETPARM(T1) a, BZ_ETPARM(T2) b, BZ_ETPARM(T3) c)
        : iter1_(a), iter2_(b), iter3_(c)
    { }

  /* Functions for reading. Because they must depend on the result
   * type, they utilize a helper class.
   */

  // For numtypes, apply operator
  template<typename T> struct readHelper {
    static T_result fastRead(const T_expr1& iter1, const T_expr2& iter2,
			     const T_expr3& iter3, diffType i) {
      return T_op::apply(iter1.fastRead(i), iter2.fastRead(i),
			 iter3.fastRead(i)); }
    static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			    const T_expr3& iter3, int i) {
      return T_op::apply(iter1[i], iter2[i], iter3[i]); };
    static T_result deref(const T_expr1& iter1, const T_expr2& iter2,
			  const T_expr3& iter3) {
      return T_op::apply(*iter1, *iter2, *iter3); };
    static T_result first_value(const T_expr1& iter1, const T_expr2& iter2,
				const T_expr3& iter3)  {
      return T_op::apply(iter1.first_value(), iter2.first_value(),
			 iter3.first_value()); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  const T_expr3& iter3, int offset, int dim) {
      return T_op::apply(iter1.shift(offset, dim),iter2.shift(offset, dim),
			 iter3.shift(offset, dim)); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  const T_expr3& iter3, int offset1,
			  int dim1, int offset2, int dim2) {
      return T_op::apply(iter1.shift(offset1, dim1, offset2, dim2),
			 iter2.shift(offset1, dim1, offset2, dim2),
			 iter3.shift(offset1, dim1, offset2, dim2)); }
    template<int N_rank>
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
      static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			      const T_expr3& iter3,
			      const TinyVector<int, N_rank> i) {
#else
      static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			      const T_expr3& iter3,
			      const TinyVector<int, N_rank>& i) {
#endif
	return T_op::apply(iter1(i), iter2(i), iter3(i) ); }
    };

    // For ET types, bypass operator and create expression
    template<typename T> struct readHelper<ETBase<T> > {
    static T_result fastRead(const T_expr1& iter1, const T_expr2& iter2,
			     const T_expr3& iter3, diffType i) {
      return T_result(iter1.fastRead(i), iter2.fastRead(i),
		      iter3.fastRead(i)); }
    static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			    const T_expr3& iter3, int i) {
      return T_result(iter1[i], iter2[i], iter3[i]); };
    static T_result deref(const T_expr1& iter1, const T_expr2& iter2,
			    const T_expr3& iter3) {
      return T_result(*iter1, *iter2, *iter3); };
    static T_result first_value(const T_expr1& iter1, const T_expr2& iter2,
				const T_expr3& iter3)  {
      return T_result(iter1.first_value(), iter2.first_value(),
		      iter3.first_value()); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  const T_expr3& iter3, int offset, int dim) {
      return T_result(iter1.shift(offset, dim),iter2.shift(offset, dim),
		      iter3.shift(offset, dim)); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  const T_expr3& iter3, int offset1,
			  int dim1, int offset2, int dim2) {
      return T_result(iter1.shift(offset1, dim1, offset2, dim2),
		      iter2.shift(offset1, dim1, offset2, dim2),
		      iter3.shift(offset1, dim1, offset2, dim2)); }
      template<int N_rank>
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
      static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			      const T_expr3& iter3,
			      const TinyVector<int, N_rank> i) {
#else
      static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			      const T_expr3& iter3,
			      const TinyVector<int, N_rank>& i) {
#endif
	return T_result(iter1(i), iter2(i), iter3(i) ); }
      };

    T_result fastRead(diffType i) const {
      return readHelper<T_typeprop>::fastRead(iter1_, iter2_, iter3_, i); }

      template<int N>
      typename tvresult<N>::Type fastRead_tv(diffType i) const
      { return typename tvresult<N>::Type(iter1_.template fastRead_tv<N>(i),
					  iter2_.template fastRead_tv<N>(i),
					  iter3_.template fastRead_tv<N>(i)); }

    T_result operator[](int i) const {
      return readHelper<T_typeprop>::indexop(iter1_, iter2_, iter3_, i); }

    template<int N_rank>
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
    T_result operator()(const TinyVector<int, N_rank> i) const {
#else
      T_result operator()(const TinyVector<int, N_rank>& i) const {
#endif
	return readHelper<T_typeprop>::indexop(iter1_, iter2_, iter3_, i); }

      T_result operator*() const {
	return readHelper<T_typeprop>::deref(iter1_, iter2_, iter3_); }

      T_result first_value() const {
	return readHelper<T_typeprop>::first_value(iter1_, iter2_, iter3_); }

    T_result shift(int offset, int dim) const {
      return readHelper<T_typeprop>::shift(iter1_, iter2_, iter3_,
					   offset, dim); }

    T_result shift(int offset1, int dim1,int offset2, int dim2) const {
      return readHelper<T_typeprop>::shift(iter1_, iter2_, iter3_,
					   offset1, dim1, offset2, dim2); }


      // ****** end reading

  bool isVectorAligned(diffType offset) const
  { return iter1_.isVectorAligned(offset) &&
      iter2_.isVectorAligned(offset) &&
      iter3_.isVectorAligned(offset); }

  template<int N>
  T_range_result operator()(const RectDomain<N>& d) const
  {
    return T_range_result(iter1_(d), iter2_(d), iter3_(d));
  }

    int ascending(const int rank) const {
        return bounds::compute_ascending(rank, bounds::compute_ascending(
            rank, iter1_.ascending(rank), iter2_.ascending(rank)),
            iter3_.ascending(rank));
    }

    int ordering(const int rank) const {
        return bounds::compute_ordering(rank, bounds::compute_ordering(
            rank, iter1_.ordering(rank), iter2_.ordering(rank)),
            iter3_.ordering(rank));
    }

    int lbound(const int rank) const {
        return bounds::compute_lbound(rank, bounds::compute_lbound(
            rank, iter1_.lbound(rank), iter2_.lbound(rank)),
            iter3_.lbound(rank));
    }

    int ubound(const int rank) const {
        return bounds::compute_ubound(rank, bounds::compute_ubound(
            rank, iter1_.ubound(rank), iter2_.ubound(rank)),
            iter3_.ubound(rank));
    }

  // defer calculation to lbound/ubound
  RectDomain<rank_> domain() const
  {
    TinyVector<int, rank_> lb, ub;
    for(int r=0; r<rank_; ++r) {
      lb[r]=lbound(r); ub[r]=ubound(r);
    }
    return RectDomain<rank_>(lb,ub);
  }

    void push(int position)
    {
        iter1_.push(position);
        iter2_.push(position);
        iter3_.push(position);
    }

    void pop(int position)
    {
        iter1_.pop(position);
        iter2_.pop(position);
        iter3_.pop(position);
    }

    void advance()
    {
        iter1_.advance();
        iter2_.advance();
        iter3_.advance();
    }

    void advance(int n)
    {
        iter1_.advance(n);
        iter2_.advance(n);
        iter3_.advance(n);
    }

    void loadStride(int rank)
    {
        iter1_.loadStride(rank);
        iter2_.loadStride(rank);
        iter3_.loadStride(rank);
    }

    bool isUnitStride(int rank) const
    {
        return iter1_.isUnitStride(rank)
            && iter2_.isUnitStride(rank)
            && iter3_.isUnitStride(rank);
    }

    bool isUnitStride() const
    {
        return iter1_.isUnitStride()
            && iter2_.isUnitStride()
            && iter3_.isUnitStride();
    }

    void advanceUnitStride()
    {
        iter1_.advanceUnitStride();
        iter2_.advanceUnitStride();
        iter3_.advanceUnitStride();
    }

    bool canCollapse(int outerLoopRank, int innerLoopRank) const
    {
        // BZ_DEBUG_MESSAGE("_bz_ArrayExprTernaryOp<>::canCollapse");
        return iter1_.canCollapse(outerLoopRank, innerLoopRank)
            && iter2_.canCollapse(outerLoopRank, innerLoopRank)
            && iter3_.canCollapse(outerLoopRank, innerLoopRank);
    }

    // this is needed for the stencil expression fastRead to work
    void _bz_offsetData(sizeType i)
    {
      iter1_._bz_offsetData(i);
      iter2_._bz_offsetData(i);
      iter3_._bz_offsetData(i);
    }

    // and these are needed for stencil expression shift to work
    void _bz_offsetData(sizeType offset, int dim)
    {
      iter1_._bz_offsetData(offset, dim);
      iter2_._bz_offsetData(offset, dim);
      iter3_._bz_offsetData(offset, dim);
    }

    void _bz_offsetData(sizeType offset1, int dim1, sizeType offset2, int dim2)
    {
      iter1_._bz_offsetData(offset1, dim1, offset2, dim2);
      iter2_._bz_offsetData(offset1, dim1, offset2, dim2);
      iter3_._bz_offsetData(offset1, dim1, offset2, dim2);
    }

    diffType suggestStride(int rank) const
    {
        diffType stride1 = iter1_.suggestStride(rank);
        diffType stride2 = iter2_.suggestStride(rank);
        diffType stride3 = iter3_.suggestStride(rank);
        return stride1 > ( stride2 = (stride2>stride3 ? stride2 : stride3) ) ?
            stride1 : stride2;
    }

    bool isStride(int rank, diffType stride) const
    {
        return iter1_.isStride(rank,stride)
            && iter2_.isStride(rank,stride)
            && iter3_.isStride(rank,stride);
    }

    template<int N>
    void moveTo(const TinyVector<int,N>& i)
    {
        iter1_.moveTo(i);
        iter2_.moveTo(i);
        iter3_.moveTo(i);
    }

    void prettyPrint(BZ_STD_SCOPE(string) &str,
        prettyPrintFormat& format) const
    {
        T_op::prettyPrint(str, format, iter1_, iter2_, iter3_);
    }

    template<typename T_shape>
    bool shapeCheck(const T_shape& shape) const
    {
        return iter1_.shapeCheck(shape)
            && iter2_.shapeCheck(shape)
            && iter3_.shapeCheck(shape);
    }


  // sliceinfo for expressions
  template<typename T1, typename T2 = nilArraySection,
	   class T3 = nilArraySection, typename T4 = nilArraySection,
	   class T5 = nilArraySection, typename T6 = nilArraySection,
	   class T7 = nilArraySection, typename T8 = nilArraySection,
	   class T9 = nilArraySection, typename T10 = nilArraySection,
	   class T11 = nilArraySection>
  class SliceInfo {
  public:
    typedef typename T_expr1::template SliceInfo<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>::T_slice T_slice1;
    typedef typename T_expr2::template SliceInfo<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>::T_slice T_slice2;
    typedef typename T_expr3::template SliceInfo<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>::T_slice T_slice3;
    typedef _bz_ArrayExprTernaryOp<T_slice1, T_slice2, T_slice3, T_op> T_slice;
};

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
        typename T7, typename T8, typename T9, typename T10, typename T11>
    typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const
    {
      return typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::T_slice
	(iter1_(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11),
	 iter2_(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11),
	 iter3_(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11));
    }

protected:
    _bz_ArrayExprTernaryOp() { }

    T_expr1 iter1_;
    T_expr2 iter2_;
    T_expr3 iter3_;
};

template<typename P_expr1, typename P_expr2, typename P_expr3,
	 typename P_expr4, typename P_op>
class _bz_ArrayExprQuaternaryOp {
public:
    typedef P_expr1 T_expr1;
    typedef P_expr2 T_expr2;
    typedef P_expr3 T_expr3;
    typedef P_expr4 T_expr4;
    typedef P_op T_op;
    typedef _bz_typename T_expr1::T_numtype T_numtype1;
    typedef _bz_typename T_expr2::T_numtype T_numtype2;
    typedef _bz_typename T_expr3::T_numtype T_numtype3;
    typedef _bz_typename T_expr4::T_numtype T_numtype4;
    typedef _bz_typename T_op::T_numtype T_numtype;

  // select return type
  typedef typename unwrapET<typename T_expr1::T_result>::T_unwrapped T_unwrapped1;
  typedef typename unwrapET<typename T_expr2::T_result>::T_unwrapped T_unwrapped2;
  typedef typename unwrapET<typename T_expr3::T_result>::T_unwrapped T_unwrapped3;
  typedef typename unwrapET<typename T_expr4::T_result>::T_unwrapped T_unwrapped4;
  typedef typename selectET2<typename T_expr1::T_typeprop,
			     typename T_expr2::T_typeprop,
			     T_numtype,
			     char>::T_selected T_intermediary1;
  typedef typename selectET2<T_intermediary1,
			     typename T_expr3::T_typeprop,
			     T_numtype,
			     char>::T_selected T_intermediary2;
  typedef typename selectET2<
    T_intermediary2,
    typename T_expr4::T_typeprop,
    T_numtype,
    _bz_ArrayExprQuaternaryOp<typename asExpr<T_unwrapped1>::T_expr,
			      typename asExpr<T_unwrapped2>::T_expr,
			      typename asExpr<T_unwrapped3>::T_expr,
			      typename asExpr<T_unwrapped4>::T_expr,
			      T_op> >::T_selected T_typeprop;
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;
  typedef typename T_op::T_numtype T_optype;

    typedef T_expr1 T_ctorArg1;
    typedef T_expr2 T_ctorArg2;
    typedef T_expr3 T_ctorArg3;
    typedef T_expr4 T_ctorArg4;
  typedef _bz_ArrayExprQuaternaryOp<_bz_typename P_expr1::T_range_result,
				    _bz_typename P_expr2::T_range_result,
				    _bz_typename P_expr3::T_range_result,
				    _bz_typename P_expr4::T_range_result,
				    P_op> T_range_result;

  static const int
  numArrayOperands = T_expr1::numArrayOperands
    + T_expr2::numArrayOperands
    + T_expr3::numArrayOperands
    + T_expr4::numArrayOperands,

    numTVOperands = T_expr1::numTVOperands +
    T_expr2::numTVOperands +
    T_expr3::numTVOperands +
    T_expr4::numTVOperands,

    numTMOperands = T_expr1::numTMOperands +
    T_expr2::numTMOperands +
    T_expr3::numTMOperands +
    T_expr4::numTMOperands,

    numIndexPlaceholders = T_expr1::numIndexPlaceholders
    + T_expr2::numIndexPlaceholders
    + T_expr3::numIndexPlaceholders
    + T_expr4::numIndexPlaceholders,

    minWidth = BZ_MIN(BZ_MIN(T_expr1::minWidth, T_expr2::minWidth),
		      BZ_MIN(T_expr3::minWidth, T_expr4::minWidth)),
    maxWidth = BZ_MAX(BZ_MAX(T_expr1::maxWidth, T_expr2::maxWidth),
		      BZ_MAX(T_expr3::maxWidth, T_expr4::maxWidth)),

    rank_ = BZ_MAX(BZ_MAX(T_expr1::rank_, T_expr2::rank_),
		  BZ_MAX(T_expr3::rank_, T_expr4::rank_));

  template<int N> struct tvresult {
    typedef _bz_ArrayExprQuaternaryOp<
      typename T_expr1::template tvresult<N>::Type,
      typename T_expr2::template tvresult<N>::Type,
      typename T_expr3::template tvresult<N>::Type,
      typename T_expr4::template tvresult<N>::Type,
      T_op> Type;
  };

    _bz_ArrayExprQuaternaryOp(
        const _bz_ArrayExprQuaternaryOp<T_expr1, T_expr2, T_expr3, T_expr4, T_op>& a)
        : iter1_(a.iter1_), iter2_(a.iter2_), iter3_(a.iter3_), iter4_(a.iter4_)
    { }

    template<typename T1, typename T2, typename T3, typename T4>
    _bz_ArrayExprQuaternaryOp(BZ_ETPARM(T1) a, BZ_ETPARM(T2) b,
			      BZ_ETPARM(T3) c, BZ_ETPARM(T4) d)
        : iter1_(a), iter2_(b), iter3_(c), iter4_(d)
    { }

  /* Functions for reading. Because they must depend on the result
   * type, they utilize a helper class.
   */

  // For numtypes, apply operator
  template<typename T> struct readHelper {
    static T_result fastRead(const T_expr1& iter1, const T_expr2& iter2,
			     const T_expr3& iter3, const T_expr4& iter4,
			     diffType i) {
      return T_op::apply(iter1.fastRead(i), iter2.fastRead(i),
			 iter3.fastRead(i), iter4.fastRead(i)); }
    static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			    const T_expr3& iter3, const T_expr4& iter4,
			    int i) {
      return T_op::apply(iter1[i], iter2[i], iter3[i], iter4[i]); };
    static T_result deref(const T_expr1& iter1, const T_expr2& iter2,
			  const T_expr3& iter3, const T_expr4& iter4) {
      return T_op::apply(*iter1, *iter2, *iter3, *iter4); };
    static T_result first_value(const T_expr1& iter1, const T_expr2& iter2,
				const T_expr3& iter3, const T_expr4& iter4)  {
      return T_op::apply(iter1.first_value(), iter2.first_value(),
			 iter3.first_value(), iter4.first_value()); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  const T_expr3& iter3, const T_expr4& iter4,
			  int offset, int dim) {
      return T_op::apply(iter1.shift(offset, dim),iter2.shift(offset, dim),
			 iter3.shift(offset, dim), iter4.shift(offset, dim)); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  const T_expr3& iter3, const T_expr4& iter4,
			  int offset1, int dim1, int offset2, int dim2) {
      return T_op::apply(iter1.shift(offset1, dim1, offset2, dim2),
			 iter2.shift(offset1, dim1, offset2, dim2),
			 iter3.shift(offset1, dim1, offset2, dim2),
			 iter4.shift(offset1, dim1, offset2, dim2)); }
    template<int N_rank>
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
      static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			      const T_expr3& iter3, const T_expr4& iter4,
			      const TinyVector<int, N_rank> i) {
#else
      static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			      const T_expr3& iter3, const T_expr4& iter4,
			      const TinyVector<int, N_rank>& i) {
#endif
	return T_op::apply(iter1(i), iter2(i), iter3(i), iter4(i) ); }
    };

    // For ET types, bypass operator and create expression
    template<typename T> struct readHelper<ETBase<T> > {
    static T_result fastRead(const T_expr1& iter1, const T_expr2& iter2,
			     const T_expr3& iter3, const T_expr4& iter4,
			     diffType i) {
      return T_result(iter1.fastRead(i), iter2.fastRead(i),
		      iter3.fastRead(i), iter4.fastRead(i)); }
    static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			    const T_expr3& iter3, const T_expr4& iter4,
			    int i) {
      return T_result(iter1[i], iter2[i], iter3[i], iter4[i]); };
    static T_result deref(const T_expr1& iter1, const T_expr2& iter2,
			  const T_expr3& iter3, const T_expr4& iter4) {
      return T_result(*iter1, *iter2, *iter3, *iter4); };
    static T_result first_value(const T_expr1& iter1, const T_expr2& iter2,
				const T_expr3& iter3, const T_expr3& iter4)  {
      return T_result(iter1.first_value(), iter2.first_value(),
		      iter3.first_value(), iter4.first_value()); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  const T_expr3& iter3, const T_expr4& iter4,
			  int offset, int dim) {
      return T_result(iter1.shift(offset, dim),iter2.shift(offset, dim),
		      iter3.shift(offset, dim), iter4.shift(offset, dim)); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  const T_expr3& iter3, const T_expr4& iter4,
			  int offset1, int dim1, int offset2, int dim2) {
      return T_result(iter1.shift(offset1, dim1, offset2, dim2),
		      iter2.shift(offset1, dim1, offset2, dim2),
		      iter3.shift(offset1, dim1, offset2, dim2),
		      iter4.shift(offset1, dim1, offset2, dim2)); }
      template<int N_rank>
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
      static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			      const T_expr3& iter3, const T_expr4& iter4,
			      const TinyVector<int, N_rank> i) {
#else
      static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			      const T_expr3& iter3, const T_expr4& iter4,
			      const TinyVector<int, N_rank>& i) {
#endif
	return T_result(iter1(i), iter2(i), iter3(i), iter4(i) ); }
      };

    T_result fastRead(diffType i) const {
      return readHelper<T_typeprop>::fastRead(iter1_, iter2_,
					      iter3_, iter4_, i); }

      template<int N>
      typename tvresult<N>::Type fastRead_tv(diffType i) const
      { return typename tvresult<N>::Type(iter1_.template fastRead_tv<N>(i),
					  iter2_.template fastRead_tv<N>(i),
					  iter3_.template fastRead_tv<N>(i),
					  iter4_.template fastRead_tv<N>(i)); }

    T_result operator[](int i) const {
      return readHelper<T_typeprop>::indexop(iter1_, iter2_,
					     iter3_, iter4_, i); }

    template<int N_rank>
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
    T_result operator()(const TinyVector<int, N_rank> i) const {
#else
      T_result operator()(const TinyVector<int, N_rank>& i) const {
#endif
	return readHelper<T_typeprop>::indexop(iter1_, iter2_,
					       iter3_, iter4_, i); }

      T_result operator*() const
      { return readHelper<T_typeprop>::deref(*iter1_, *iter2_,
					     *iter3_, *iter4_); }

      T_result first_value() const {
	return readHelper<T_typeprop>::first_value(iter1_, iter2_,
						   iter3_, iter4_); }

    T_result shift(int offset, int dim) const {
      return readHelper<T_typeprop>::shift(iter1_, iter2_, iter3_,
					   iter4_, offset, dim); }

    T_result shift(int offset1, int dim1,int offset2, int dim2) const {
      return readHelper<T_typeprop>::shift(iter1_, iter2_, iter3_, iter4_,
					   offset1, dim1, offset2, dim2); }

      // ****** end reading

  bool isVectorAligned(diffType offset) const
  { return iter1_.isVectorAligned(offset) &&
      iter2_.isVectorAligned(offset) &&
      iter3_.isVectorAligned(offset) &&
      iter3_.isVectorAligned(offset); }

  template<int N>
  T_range_result operator()(const RectDomain<rank_>& d) const
  {
    return T_range_result(iter1_(d), iter2_(d), iter3_(d), iter4_(d));
  }

    int ascending(const int rank) const {
        return bounds::compute_ascending(rank,
					 bounds::compute_ascending(rank,
								   iter1_.ascending(rank),
								   iter2_.ascending(rank)),
					 bounds::compute_ascending(rank,
								   iter3_.ascending(rank),
								   iter4_.ascending(rank)));
    }

    int ordering(const int rank) const {
        return bounds::compute_ordering(rank,
					bounds::compute_ordering(rank,
								 iter1_.ordering(rank),
								 iter2_.ordering(rank)),
					bounds::compute_ordering(rank,
								 iter3_.ordering(rank),
								 iter4_.ordering(rank)));
    }

    int lbound(const int rank) const {
        return bounds::compute_lbound(rank,
				      bounds::compute_lbound(rank,
							     iter1_.lbound(rank),
							     iter2_.lbound(rank)),
				      bounds::compute_lbound(rank,
							     iter3_.lbound(rank),
							     iter4_.lbound(rank)));
    }

    int ubound(const int rank) const {
        return bounds::compute_ubound(rank,
				      bounds::compute_ubound(rank,
							     iter1_.ubound(rank),
							     iter2_.ubound(rank)),
				      bounds::compute_ubound(rank,
							     iter3_.ubound(rank),
							     iter4_.ubound(rank)));
    }

  // defer calculation to lbound/ubound
  RectDomain<rank_> domain() const
  {
    TinyVector<int, rank_> lb, ub;
    for(int r=0; r<rank_; ++r) {
      lb[r]=lbound(r); ub[r]=ubound(r);
    }
    return RectDomain<rank_>(lb,ub);
  }

    void push(int position)
    {
        iter1_.push(position);
        iter2_.push(position);
        iter3_.push(position);
        iter4_.push(position);
    }

    void pop(int position)
    {
        iter1_.pop(position);
        iter2_.pop(position);
        iter3_.pop(position);
        iter4_.pop(position);
    }

    void advance()
    {
        iter1_.advance();
        iter2_.advance();
        iter3_.advance();
        iter4_.advance();
    }

    void advance(int n)
    {
        iter1_.advance(n);
        iter2_.advance(n);
        iter3_.advance(n);
        iter4_.advance(n);
    }

    void loadStride(int rank)
    {
        iter1_.loadStride(rank);
        iter2_.loadStride(rank);
        iter3_.loadStride(rank);
        iter4_.loadStride(rank);
    }

    bool isUnitStride(int rank) const
    {
        return iter1_.isUnitStride(rank)
            && iter2_.isUnitStride(rank)
            && iter3_.isUnitStride(rank)
            && iter4_.isUnitStride(rank);
    }

    bool isUnitStride() const
    {
        return iter1_.isUnitStride()
            && iter2_.isUnitStride()
            && iter3_.isUnitStride()
            && iter4_.isUnitStride();
    }

    void advanceUnitStride()
    {
        iter1_.advanceUnitStride();
        iter2_.advanceUnitStride();
        iter3_.advanceUnitStride();
        iter4_.advanceUnitStride();
    }

    bool canCollapse(int outerLoopRank, int innerLoopRank) const
    {
        // BZ_DEBUG_MESSAGE("_bz_ArrayExprQuaternaryOp<>::canCollapse");
        return iter1_.canCollapse(outerLoopRank, innerLoopRank)
            && iter2_.canCollapse(outerLoopRank, innerLoopRank)
            && iter3_.canCollapse(outerLoopRank, innerLoopRank)
            && iter4_.canCollapse(outerLoopRank, innerLoopRank);
    }

  // this is needed for the stencil expression fastRead to work
  void _bz_offsetData(sizeType i)
  {
    iter1_._bz_offsetData(i);
    iter2_._bz_offsetData(i);
    iter3_._bz_offsetData(i);
    iter4_._bz_offsetData(i);
  }

    // and these are needed for stencil expression shift to work
    void _bz_offsetData(sizeType offset, int dim)
    {
      iter1_._bz_offsetData(offset, dim);
      iter2_._bz_offsetData(offset, dim);
      iter3_._bz_offsetData(offset, dim);
      iter4_._bz_offsetData(offset, dim);
    }

    void _bz_offsetData(sizeType offset1, int dim1, sizeType offset2, int dim2)
    {
      iter1_._bz_offsetData(offset1, dim1, offset2, dim2);
      iter2_._bz_offsetData(offset1, dim1, offset2, dim2);
      iter3_._bz_offsetData(offset1, dim1, offset2, dim2);
      iter4_._bz_offsetData(offset1, dim1, offset2, dim2);
    }

    diffType suggestStride(int rank) const
    {
        diffType stride1 = iter1_.suggestStride(rank);
        diffType stride2 = iter2_.suggestStride(rank);
        diffType stride3 = iter3_.suggestStride(rank);
        diffType stride4 = iter4_.suggestStride(rank);

        //return stride1 > ( stride2 = (stride2>stride3 ? stride2 : stride3) ) ?
	//  stride1 : stride2;
	return std::max(std::max(stride1, stride2),
			std::max(stride3, stride4));
    }

    bool isStride(int rank, diffType stride) const
    {
        return iter1_.isStride(rank,stride)
            && iter2_.isStride(rank,stride)
            && iter3_.isStride(rank,stride)
            && iter4_.isStride(rank,stride);
    }

    template<int N>
    void moveTo(const TinyVector<int,N>& i)
    {
        iter1_.moveTo(i);
        iter2_.moveTo(i);
        iter3_.moveTo(i);
        iter4_.moveTo(i);
    }

    void prettyPrint(BZ_STD_SCOPE(string) &str,
        prettyPrintFormat& format) const
    {
        T_op::prettyPrint(str, format, iter1_, iter2_, iter3_, iter4_);
    }

    template<typename T_shape>
    bool shapeCheck(const T_shape& shape) const
    {
        return iter1_.shapeCheck(shape)
            && iter2_.shapeCheck(shape)
            && iter3_.shapeCheck(shape)
            && iter4_.shapeCheck(shape);
    }


  // sliceinfo for expressions
  template<typename T1, typename T2 = nilArraySection,
	   class T3 = nilArraySection, typename T4 = nilArraySection,
	   class T5 = nilArraySection, typename T6 = nilArraySection,
	   class T7 = nilArraySection, typename T8 = nilArraySection,
	   class T9 = nilArraySection, typename T10 = nilArraySection,
	   class T11 = nilArraySection>
  class SliceInfo {
  public:
    typedef typename T_expr1::template SliceInfo<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>::T_slice T_slice1;
    typedef typename T_expr2::template SliceInfo<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>::T_slice T_slice2;
    typedef typename T_expr3::template SliceInfo<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>::T_slice T_slice3;
    typedef typename T_expr4::template SliceInfo<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>::T_slice T_slice4;
    typedef _bz_ArrayExprQuaternaryOp<T_slice1, T_slice2, T_slice3, T_slice4, T_op> T_slice;
};

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
        typename T7, typename T8, typename T9, typename T10, typename T11>
    typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const
    {
      return typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::T_slice
	(iter1_(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11),
	 iter2_(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11),
	 iter3_(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11),
	 iter4_(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11));
    }

protected:
    _bz_ArrayExprQuaternaryOp() { }

    T_expr1 iter1_;
    T_expr2 iter2_;
    T_expr3 iter3_;
    T_expr4 iter4_;
};

template<typename P_numtype>
class _bz_ArrayExprConstant {
public:
    typedef P_numtype T_numtype;
  typedef typename opType<T_numtype>::T_optype T_optype;
  typedef typename asET<T_numtype>::T_wrapped T_typeprop;
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;

    typedef T_numtype T_ctorArg1;
    typedef int       T_ctorArg2;    // dummy
  typedef _bz_ArrayExprConstant<P_numtype> T_range_result;
    static const int
        numArrayOperands = 0,
        numTVOperands = 0,
        numTMOperands = 0,
        numIndexPlaceholders = 0,
      minWidth = simdTypes<T_numtype>::vecWidth,
      maxWidth = simdTypes<T_numtype>::vecWidth,
        rank_ = 0;

  /** For the purpose of vectorizing across the container (as opposed
      to for operating on multicomponent types), a constant is always
      a constant. */
  template<int N> struct tvresult {
    typedef _bz_ArrayExprConstant<T_numtype> Type;
  };

    _bz_ArrayExprConstant(const _bz_ArrayExprConstant<T_numtype>& a)
        : value_(a.value_)
    { }

    _bz_ArrayExprConstant(T_numtype value)
        : value_(BZ_NO_PROPAGATE(value))
    {
    }

    // tiny() and huge() return the smallest and largest representable
    // integer values.  See <blitz/numinquire.h>
    // NEEDS_WORK: use tiny(int()) once numeric_limits<T> available on
    // all platforms

    int ascending(const int) const { return INT_MIN; }
    int ordering(const int)  const { return INT_MIN; }
    int lbound(const int)    const { return INT_MIN; }
    int ubound(const int)    const { return INT_MAX; }

  // there is no rank so we use highest possible
  RectDomain<12> domain() const;

    // NEEDS_WORK: use huge(int()) once numeric_limits<T> available on
    // all platforms

    T_result operator*()   const { return value_; }
    T_result first_value() const { return value_; }

#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
    template<int N_rank>
    T_result operator()(const TinyVector<int,N_rank>) const
    { return value_; }
#else
    template<int N_rank>
    T_result operator()(const TinyVector<int,N_rank>&) const
    { return value_; }
#endif

  template<int N_rank>
  const _bz_ArrayExprConstant& operator()(const RectDomain<N_rank>& d) const
  {
    return *this;
  }

    void push(int) { }
    void pop(int) { }
    void advance() { }
    void advance(int) { }
    void loadStride(int) { }

    bool isUnitStride(int) const
    { return true; }

    bool isUnitStride() const
    { return true; }

    void advanceUnitStride()
    { }

    bool canCollapse(int,int) const
    { return true; }

    T_numtype operator[](int) const
    { return value_; }

  // this must return by reference, because for multicomponent arrays
  // it gets turned into an iterator by the containing expression.
  const T_numtype& fastRead(diffType) const
    { return value_; }

  template<int N>
  typename tvresult<N>::Type fastRead_tv(diffType i) const
  { (void)i; return value_; }

  bool isVectorAligned(diffType offset) const
  { (void)offset; return true; }

  // this is needed for the stencil expression fastRead to work
  void _bz_offsetData(sizeType i) const{};

    // and these are needed for stencil expression shift to work
  void _bz_offsetData(sizeType offset, int dim) const {};

  void _bz_offsetData(sizeType offset1, int dim1, sizeType offset2, int dim2) const {};

    diffType suggestStride(int) const
    { return 1; }

    bool isStride(int,diffType) const
    { return true; }

    void moveTo(int) const { }

    T_result shift(int offset, int dim) const {return value_;}

    T_result shift(int offset1, int dim1,int offset2, int dim2) const
    { return value_;}

    template<int N_rank>
    void moveTo(const TinyVector<int,N_rank>&) const { }

    void prettyPrint(BZ_STD_SCOPE(string) &str,
        prettyPrintFormat& format) const
    {
        if (format.tersePrintingSelected())
            str += format.nextScalarOperandSymbol();
        else
            str += BZ_DEBUG_TEMPLATE_AS_STRING_LITERAL(T_numtype);
    }

    template<typename T_shape>
    bool shapeCheck(const T_shape&) const
    { return true; }


  // sliceinfo for expressions
  template<typename T1, typename T2 = nilArraySection,
	   class T3 = nilArraySection, typename T4 = nilArraySection,
	   class T5 = nilArraySection, typename T6 = nilArraySection,
	   class T7 = nilArraySection, typename T8 = nilArraySection,
	   class T9 = nilArraySection, typename T10 = nilArraySection,
	   class T11 = nilArraySection>
  class SliceInfo {
  public:
    typedef _bz_ArrayExprConstant<T_numtype> T_slice;
};

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
        typename T7, typename T8, typename T9, typename T10, typename T11>
    typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const
    {
      return *this;
    }

protected:
    _bz_ArrayExprConstant() { }

    T_numtype value_;
};


BZ_NAMESPACE_END

#include <blitz/array/asexpr.h>

#endif // BZ_ARRAYEXPR_H
