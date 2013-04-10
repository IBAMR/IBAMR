// -*- C++ -*-
/***************************************************************************
 * blitz/array/where.h  where(X,Y,Z) operator for array expressions
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
#ifndef BZ_ARRAYWHERE_H
#define BZ_ARRAYWHERE_H

#include <blitz/blitz.h>
#include <blitz/promote.h>
#include <blitz/prettyprint.h>
#include <blitz/bounds.h>
#include <blitz/meta/metaprog.h>
#include <blitz/tinyvec2.h>
#include <blitz/array/domain.h>
#include <blitz/array/asexpr.h>

BZ_NAMESPACE(blitz)

template<typename P_expr1, typename P_expr2, typename P_expr3>
class _bz_ArrayWhere {

public:
    typedef P_expr1 T_expr1;
    typedef P_expr2 T_expr2;
    typedef P_expr3 T_expr3;
    typedef _bz_typename T_expr2::T_numtype T_numtype2;
    typedef _bz_typename T_expr3::T_numtype T_numtype3;
    typedef BZ_PROMOTE(T_numtype2, T_numtype3) T_numtype;
    typedef T_expr1 T_ctorArg1;
    typedef T_expr2 T_ctorArg2;
    typedef T_expr3 T_ctorArg3;
  typedef _bz_ArrayWhere<_bz_typename P_expr1::T_range_result,
			 _bz_typename P_expr2::T_range_result,
			 _bz_typename P_expr3::T_range_result> T_range_result;

  // select return type
  typedef typename unwrapET<typename T_expr1::T_result>::T_unwrapped T_unwrapped1;
  typedef typename unwrapET<typename T_expr2::T_result>::T_unwrapped T_unwrapped2;
  typedef typename unwrapET<typename T_expr3::T_result>::T_unwrapped T_unwrapped3;
  typedef typename selectET2<typename T_expr1::T_typeprop,
			     typename T_expr2::T_typeprop,
			     T_numtype,
			     char>::T_selected T_intermediary;

  typedef typename selectET2<T_intermediary,
			     typename T_expr3::T_typeprop,
			     T_numtype,
			     _bz_ArrayWhere<typename asExpr<T_unwrapped1>::T_expr,
					    typename asExpr<T_unwrapped2>::T_expr,
					    typename asExpr<T_unwrapped3>::T_expr
					    > >::T_selected T_typeprop;
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;
  typedef T_numtype T_optype;

    static const int
        numArrayOperands = P_expr1::numArrayOperands
                         + P_expr2::numArrayOperands
                         + P_expr3::numArrayOperands,
        numTVOperands = T_expr1::numTVOperands +
      T_expr2::numTVOperands +
      T_expr3::numTVOperands,
        numTMOperands = T_expr1::numTMOperands +
      T_expr2::numTMOperands +
      T_expr3::numTMOperands,
        numIndexPlaceholders = P_expr1::numIndexPlaceholders
                             + P_expr2::numIndexPlaceholders
                             + P_expr3::numIndexPlaceholders,
      minWidth = BZ_MIN(BZ_MIN(T_expr1::minWidth, T_expr2::minWidth),
			T_expr3::minWidth),
      maxWidth = BZ_MAX(BZ_MAX(T_expr1::maxWidth, T_expr2::maxWidth),
			T_expr3::maxWidth),
      rank_ = BZ_MAX(BZ_MAX(T_expr1::rank_, T_expr2::rank_),
		     T_expr3::rank_);

  template<int N> struct tvresult {
    typedef _bz_ArrayWhere<
      typename T_expr1::template tvresult<N>::Type,
      typename T_expr2::template tvresult<N>::Type,
      typename T_expr3::template tvresult<N>::Type> Type;
  };

    _bz_ArrayWhere(const _bz_ArrayWhere<T_expr1,T_expr2,T_expr3>& a)
      : iter1_(a.iter1_), iter2_(a.iter2_), iter3_(a.iter3_)
    { }

    template<typename T1, typename T2, typename T3>
    _bz_ArrayWhere(BZ_ETPARM(T1) a, BZ_ETPARM(T2) b, BZ_ETPARM(T3) c)
      : iter1_(a), iter2_(b), iter3_(c)
    { }

    T_numtype operator*() const
    { return (*iter1_) ? (*iter2_) : (*iter3_); }

  /* Functions for reading. Because they must depend on the result
   * type, they utilize a helper class.
   */

  // For numtypes, apply operator
  template<typename T> struct readHelper {
    static T_result fastRead(const T_expr1& iter1, const T_expr2& iter2,
			     const T_expr3& iter3, diffType i) {
      return iter1.fastRead(i) ? iter2.fastRead(i) : iter3.fastRead(i); }
    static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			    const T_expr3& iter3, int i) {
      return iter1[i] ? iter2[i] : iter3[i]; }
    static T_result first_value(const T_expr1& iter1, const T_expr2& iter2,
				const T_expr3& iter3)  {
      return iter1.first_value() ?
	iter2.first_value() : iter3.first_value(); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  const T_expr3& iter3, int offset, int dim) {
      return iter1.shift(offset, dim) ? iter2.shift(offset, dim) :
	iter3.shift(offset, dim); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  const T_expr3& iter3, int offset1, int dim1,
			  int offset2, int dim2) {
      return iter1.shift(offset1, dim1, offset2, dim2) ?
	iter2.shift(offset1, dim1, offset2, dim2) :
	iter3.shift(offset1, dim1, offset2, dim2); }
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
	return iter1(i) ? iter2(i) : iter3(i); }
      };

    // For ET types, bypass operator and create expression
    template<typename T> struct readHelper<ETBase<T> > {
    static T_result fastRead(const T_expr1& iter1, const T_expr2& iter2,
			     const T_expr3& iter3, diffType i) {
      return T_result(iter1.fastRead(i), iter2.fastRead(i), iter3.fastRead(i)); }
    static T_result indexop(const T_expr1& iter1, const T_expr2& iter2,
			    const T_expr3& iter3, int i) {
      return T_result(iter1[i], iter2[i], iter3[i]); };
    static T_result first_value(const T_expr1& iter1, const T_expr2& iter2,
				const T_expr3& iter3)  {
      return T_result(iter1.first_value(), iter2.first_value(),
		      iter3.first_value()); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  const T_expr3& iter3, int offset, int dim) {
      return T_result(iter1.shift(offset, dim), iter2.shift(offset, dim),
		      iter3.shift(offset, dim)); }
    static T_result shift(const T_expr1& iter1, const T_expr2& iter2,
			  const T_expr3& iter3, int offset1, int dim1,
			  int offset2, int dim2) {
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

    T_range_result operator()(const RectDomain<rank_>& d) const
  { return T_range_result(iter1_(d), iter2_(d), iter3_(d)); }

    int ascending(const int rank) const
    {
        return bounds::compute_ascending(rank, bounds::compute_ascending(
          rank, iter1_.ascending(rank), iter2_.ascending(rank)),
          iter3_.ascending(rank));
    }

    int ordering(const int rank) const
    {
        return bounds::compute_ordering(rank, bounds::compute_ordering(
          rank, iter1_.ordering(rank), iter2_.ordering(rank)),
          iter3_.ordering(rank));
    }

    int lbound(const int rank) const
    {
        return bounds::compute_lbound(rank, bounds::compute_lbound(
          rank, iter1_.lbound(rank), iter2_.lbound(rank)),
          iter3_.lbound(rank));
    }

    int ubound(const int rank) const
    {
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
        return iter1_.canCollapse(outerLoopRank, innerLoopRank)
            && iter2_.canCollapse(outerLoopRank, innerLoopRank)
            && iter3_.canCollapse(outerLoopRank, innerLoopRank);
    }

    template<int N_rank>
    void moveTo(const TinyVector<int,N_rank>& i)
    {
        iter1_.moveTo(i);
        iter2_.moveTo(i);
        iter3_.moveTo(i);
    }

  // this is needed for the stencil expression fastRead to work
  void _bz_offsetData(sizeType i)
  {
    iter1_._bz_offsetData(i);
    iter2_._bz_offsetData(i);
    iter3_._bz_offsetData(i);
  }

    diffType suggestStride(int rank) const
    {
        diffType stride1 = iter1_.suggestStride(rank);
        diffType stride2 = iter2_.suggestStride(rank);
        diffType stride3 = iter3_.suggestStride(rank);
        return stride1>(stride2=(stride2>stride3?stride2:stride3))?stride1:stride2;
    }

    bool isStride(int rank, diffType stride) const
    {
        return iter1_.isStride(rank,stride)
            && iter2_.isStride(rank,stride)
            && iter3_.isStride(rank,stride);
    }

    void prettyPrint(BZ_STD_SCOPE(string) &str,
        prettyPrintFormat& format) const
    {
        str += "where(";
        iter1_.prettyPrint(str,format);
        str += ",";
        iter2_.prettyPrint(str,format);
        str += ",";
        iter3_.prettyPrint(str,format);
        str += ")";
    }

    template<typename T_shape>
    bool shapeCheck(const T_shape& shape) const
    {
        int t1 = iter1_.shapeCheck(shape);
        int t2 = iter2_.shapeCheck(shape);
        int t3 = iter3_.shapeCheck(shape);

        return t1 && t2 && t3;
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
    typedef _bz_ArrayWhere<T_slice1, T_slice2, T_slice3> T_slice;
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

private:
    _bz_ArrayWhere() { }

    T_expr1 iter1_;
    T_expr2 iter2_;
    T_expr3 iter3_;
};

template<typename T1, typename T2, typename T3>
inline
_bz_ArrayExpr<_bz_ArrayWhere<_bz_typename asExpr<T1>::T_expr,
    _bz_typename asExpr<T2>::T_expr, _bz_typename asExpr<T3>::T_expr> >
where(const T1& a, const T2& b, const T3& c)
{
    return _bz_ArrayExpr<_bz_ArrayWhere<_bz_typename asExpr<T1>::T_expr,
       _bz_typename asExpr<T2>::T_expr,
       _bz_typename asExpr<T3>::T_expr> >(a,b,c);
}

BZ_NAMESPACE_END

#endif // BZ_ARRAYWHERE_H

