// -*- C++ -*-
/***************************************************************************
 * blitz/array/stencil-et.h  Expression-template-capable stencils
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
#ifndef BZ_ARRAY_STENCIL_ET_H
#define BZ_ARRAY_STENCIL_ET_H

#include <blitz/blitz.h>
#include <blitz/tinyvec2.h>
#include <blitz/prettyprint.h>
#include <blitz/array/domain.h>

BZ_NAMESPACE(blitz)

/* Stencils as currently implemented rely on being able to give an
   iterator type to the operator. Some methods have been implemented
   to make index traversal work so that operands can include index
   remappings, but it's probably not very efficient. Operands
   containing reductions or index placeholders don't work. */

// Utility function for shrinking the domain of the expression the
// stencil operates on. The minb and maxb may be defined in fewer
// dimensions than the expression, in which the shape is the same in
// higher dims.
template<int N_rank, int N_rank2>
RectDomain<N_rank> _bz_shrinkDomain(const RectDomain<N_rank>& d,
				    const TinyVector<int,N_rank2>& minb,
				    const TinyVector<int,N_rank2>& maxb)
{
  TinyVector<int, N_rank> lb(d.lbound()), ub(d.ubound());
  // no slicing w tinyvector
  for(int i=0; i<N_rank2; ++i) {
    lb[i] = lb[i] - ((minb[i]<0) ? minb[i] : 0);
    ub[i] = ub[i] - ((maxb[i]>0) ? maxb[i] : 0);
  }
  return RectDomain<N_rank>(lb, ub);
}


// Helper function that converts slices to unit ranges
template<typename T>
Range _bz_makeRange(const T& r)
{
  return Range(r);
}

inline nilArraySection _bz_makeRange(const nilArraySection& r)
{
  return r;
}


/** necessary because we can't have an #ifdef in the macros
    \todo convert to _bz_Indexpar. */
template<typename T> struct _bz_IndexParameter {
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
  typedef T type;
#else
  typedef const T& type;
#endif
};

/* There are resolution problems with the functions in this file that
   generate ET objects and the functions in stencilops.h that define
   the stencil operators. Because the operators are defined as
   all-encompassing templates, they will match to arrays with higher
   priority than the functions here will match to ETBase<..>. This
   means we must have explicit functions that match Arrays even though
   the general ETBase function would work. */

/** ET base class for applying a stencil to an expression. */
template<typename P_expr, _bz_typename P_result>
class _bz_StencilExpr {
 public:
  typedef P_expr T_expr;
  typedef P_result T_numtype;
  typedef T_expr T_ctorArg1;
  typedef int    T_ctorArg2;    // dummy

  static const int
    numArrayOperands = T_expr::numArrayOperands,
    numTVOperands = T_expr::numTVOperands,
    numTMOperands = T_expr::numTMOperands,
    numIndexPlaceholders = T_expr::numIndexPlaceholders,
  /// vectorization of stencils is disabled by setting this to 0.
    minWidth = 0, //T_expr::minWidth,
    maxWidth = 0, //T_expr::maxWidth,
    rank_ = T_expr::rank_;

 _bz_StencilExpr(const _bz_StencilExpr<T_expr, T_numtype>& a)
   : iter_(a.iter_)
    { }

 _bz_StencilExpr(BZ_ETPARM(T_expr) a)
   : iter_(a)
  { }

 _bz_StencilExpr(_bz_typename T_expr::T_ctorArg1 a)
   : iter_(a)
  { }

#if BZ_TEMPLATE_CTOR_DOESNT_CAUSE_HAVOC
  template<typename T1>
    explicit _bz_StencilExpr(BZ_ETPARM(T1) a)
    : iter_(a)
  { }
#endif

  int ascending(const int rank) const { return iter_.ascending(rank); }
  int ordering(const int rank)  const { return iter_.ordering(rank);  }
  int lbound(const int rank) const
  {
    return iter_.lbound(rank);
  }
  int ubound(const int rank) const
  {
    return iter_.ubound(rank);
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

  //T_numtype first_value() const { return iter_(iter_.lbound()); }

  /** Vectorization doesn't make sense for stencils, so we say so. */
  bool isVectorAligned(diffType offset) const {
    return false; }

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

  void moveTo(const TinyVector<int,_bz_StencilExpr::rank_>& i)
  {
    iter_.moveTo(i);
  }

  bool canCollapse(int outerLoopRank, int innerLoopRank) const
  {
    // BZ_DEBUG_MESSAGE("_bz_StencilExpr<>::canCollapse");
    return iter_.canCollapse(outerLoopRank, innerLoopRank);
  }

  diffType suggestStride(int rank) const
  { return iter_.suggestStride(rank); }

  bool isStride(int rank, diffType stride) const
  { return iter_.isStride(rank,stride); }

  void _bz_offsetData(sizeType i) { iter_._bz_offsetData(i); }

  void prettyPrint(BZ_STD_SCOPE(string) &str) const
  {
    prettyPrintFormat format(true);  /* Terse formatting by default */
    iter_.prettyPrint(str, format);
  }

  template<typename T_shape>
    bool shapeCheck(const T_shape& shape) const
    { return iter_.shapeCheck(shape); }

 protected:
  _bz_StencilExpr() { }

  // it needs to be mutable because in the "conceptually const"
  // methods shift and fastread we offset the iterator but undo it
  // later
  mutable P_expr iter_;
};

/** ET base class for applying a stencil to an expression. */
template<typename P_expr1, typename P_expr2, _bz_typename P_result>
class _bz_StencilExpr2 {
 public:
  typedef P_expr1 T_expr1;
  typedef P_expr2 T_expr2;
  typedef _bz_typename T_expr1::T_numtype T_numtype1;
  typedef _bz_typename T_expr2::T_numtype T_numtype2;
  typedef P_result T_numtype;
  typedef T_expr1 T_ctorArg1;
  typedef T_expr2 T_ctorArg2;

  static const int
  numArrayOperands = T_expr1::numArrayOperands
    + T_expr2::numArrayOperands,
    numTVOperands = T_expr1::numTVOperands + T_expr2::numTVOperands,
    numTMOperands = T_expr1::numTMOperands + T_expr2::numTMOperands,
    numIndexPlaceholders = T_expr1::numIndexPlaceholders
    + T_expr2::numIndexPlaceholders,
  /// vectorization of stencils is disabled by setting this to 0.
    minWidth = 0, //BZ_MIN(T_expr1::minWidth, T_expr2::minWidth),
    maxWidth = 0, //BZ_MAX(T_expr1::maxWidth, T_expr2::maxWidth),
      rank_ = BZ_MAX(T_expr1::rank_, T_expr2::rank_);

  _bz_StencilExpr2(const _bz_StencilExpr2<T_expr1, T_expr2, T_numtype>& a)
    : iter1_(a.iter1_), iter2_(a.iter2_)
    { }

  _bz_StencilExpr2(BZ_ETPARM(T_expr1) a, BZ_ETPARM(T_expr2) b)
    : iter1_(a), iter2_(b)
  { }

  /*
  // what is this for?
 _bz_StencilExpr2(_bz_typename T_expr1::T_ctorArg1 a)
   : iter_(a)
  { }
  */

#if BZ_TEMPLATE_CTOR_DOESNT_CAUSE_HAVOC
  template<typename T1, typename T2>
  explicit _bz_StencilExpr2(BZ_ETPARM(T1) a, BZ_ETPARM(T2) b)
    : iter1_(a), iter2_(b)
  { }
#endif

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

  /** Vectorization doesn't make sense for stencils, so we say so. */
  bool isVectorAligned(diffType offset) const {
    return false; }

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
        return iter1_.canCollapse(outerLoopRank, innerLoopRank)
            && iter2_.canCollapse(outerLoopRank, innerLoopRank);
    }

    void moveTo(const TinyVector<int,rank_>& i)
    {
        iter1_.moveTo(i);
        iter2_.moveTo(i);
    }

  void _bz_offsetData(sizeType i)
  {
    iter1_._bz_offsetData(i);
    iter2_._bz_offsetData(i);
  }

  template<typename T_shape>
    bool shapeCheck(const T_shape& shape) const
  { return iter1_.shapeCheck(shape) && iter2_.shapeCheck(shape); }

 protected:
  _bz_StencilExpr2() { }

  mutable P_expr1 iter1_;
  mutable P_expr2 iter2_;
};

BZ_NAMESPACE_END

// now include the generated file containing all the implementation classes.
#include <blitz/array/stencil-classes.cc>

#endif // BZ_ARRAY_STENCIL_ET_H
