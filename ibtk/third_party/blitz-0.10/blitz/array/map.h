// -*- C++ -*-
/***************************************************************************
 * blitz/array/map.h      Declaration of the ArrayIndexMapping class
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

/*
 * ArrayIndexMapping is used to implement tensor array notation.  For
 * example:
 *
 * Array<float, 2> A, B;
 * firstIndex i;
 * secondIndex j;
 * thirdIndex k;
 * Array<float, 3> C = A(i,j) * B(j,k);
 *
 * For expression templates purposes, something like B(j,k) is represented
 * by an instance of class ArrayIndexMapping.  This class maps an array onto
 * the destination array coordinate system, e.g. B(j,k) -> C(i,j,k)
 */

#ifndef BZ_ARRAYMAP_H
#define BZ_ARRAYMAP_H

#include <blitz/blitz.h>
#include <blitz/prettyprint.h>
#include <blitz/et-forward.h>
#include <blitz/tinyvec2.h>
#include <blitz/array/domain.h>

BZ_NAMESPACE(blitz)

/*
 * _bz_doArrayIndexMapping is a helper class that does the index
 * remapping.  It is specialized for ranks 1, 2, 3, ..., 11.
 */

template<int N_rank>
struct _bz_doArrayIndexMapping {
  static const int rank=N_rank;
  template<typename T_expr, int N_inputRank>
  static typename T_expr::T_numtype map(const ETBase<T_expr>& expr,
        const TinyVector<int,N_inputRank>&, int, int, int, int, int, int,
        int, int, int, int, int)
    {
        // If you try to use an array index mapping on an array with
        // rank greater than 11, then you'll get a precondition failure
        // here.
      BZPRECHECK(0,"Index mappings for containers of rank>11 not implemented");
      return T_expr::T_numtype();
    }
};

template<>
struct _bz_doArrayIndexMapping<1> {
  static const int rank=1;
  template<typename T_expr, int N_inputRank>
  static typename T_expr::T_numtype map(const ETBase<T_expr>& expr,
        const TinyVector<int,N_inputRank>& index, int i0, int, int, int, int,
        int, int, int, int, int, int)
    {
      // this was the case when it took an array. is it necessary?
      BZPRECHECK(T_expr::rank_==rank,
		 "Rank confusion in _bz_doArrayIndexMapping");
      return expr.unwrap()(index[i0]);
    }

  template<int N_inputRank>
  static TinyVector<int,rank> map_dims(const TinyVector<int,N_inputRank>& index,
			 int i0, int, int, int, int,
			 int, int, int, int, int, int)
    {
      // this might be slower but unlike Array, FAI doesn't have 11
      // overloaded moveTo methods.
      const TinyVector<int,rank> newindex(index[i0]);
      return newindex;
    }
};


template<>
struct _bz_doArrayIndexMapping<2> {
  static const int rank=2;

  template<typename T_expr, int N_inputRank>
  static typename T_expr::T_numtype map(const ETBase<T_expr>& expr,
        const TinyVector<int,N_inputRank>& index, int i0, int i1, int,
        int, int, int, int, int, int, int, int)
    {
      // this was the case when it took an array. is it necessary?
      BZPRECHECK(T_expr::rank_==rank,
		 "Rank confusion in _bz_doArrayIndexMapping");
      return expr.unwrap()(index[i0], index[i1]);
    }

  template<int N_inputRank>
  static TinyVector<int,rank> map_dims(const TinyVector<int,N_inputRank>& index,
			 int i0, int i1, int, int, int,
			 int, int, int, int, int, int)
    {
      return TinyVector<int,rank>(index[i0], index[i1]);
    }
};

template<>
struct _bz_doArrayIndexMapping<3> {
  static const int rank=3;

  template<typename T_expr, int N_inputRank>
  static typename T_expr::T_numtype map(const ETBase<T_expr>& expr,
        const TinyVector<int,N_inputRank>& index, int i0, int i1, int i2,
        int, int, int, int, int, int, int, int)
    {
      // this was the case when it took an array. is it necessary?
      BZPRECHECK(T_expr::rank_==rank,
		 "Rank confusion in _bz_doArrayIndexMapping");
      return expr.unwrap()(index[i0], index[i1], index[i2]);
    }

  template<int N_inputRank>
  static TinyVector<int,rank> map_dims(const TinyVector<int,N_inputRank>& index,
			 int i0, int i1, int i2, int, int,
			 int, int, int, int, int, int)
    {
      return TinyVector<int,rank>(index[i0], index[i1], index[i2]);
    }
};

template<>
struct _bz_doArrayIndexMapping<4> {
  static const int rank=4;

  template<typename T_expr, int N_inputRank>
  static typename T_expr::T_numtype map(const ETBase<T_expr>& expr,
        const TinyVector<int,N_inputRank>& index, int i0, int i1, int i2,
        int i3, int, int, int, int, int, int, int)
    {
      // this was the case when it took an array. is it necessary?
      BZPRECHECK(T_expr::rank_==rank,
		 "Rank confusion in _bz_doArrayIndexMapping");
        return expr.unwrap()(index[i0], index[i1], index[i2], index[i3]);
    }

  template<int N_inputRank>
  static TinyVector<int,rank> map_dims(const TinyVector<int,N_inputRank>& index,
			 int i0, int i1, int i2, int i3, int,
			 int, int, int, int, int, int)
    {
      return TinyVector<int,rank>(index[i0], index[i1], index[i2],
				  index[i3]);
    }
};

template<>
struct _bz_doArrayIndexMapping<5> {
  static const int rank=5;

  template<typename T_expr, int N_inputRank>
  static typename T_expr::T_numtype map(const ETBase<T_expr>& expr,
        const TinyVector<int,N_inputRank>& index, int i0, int i1, int i2,
        int i3, int i4, int, int, int, int, int, int)
    {
      // this was the case when it took an array. is it necessary?
      BZPRECHECK(T_expr::rank_==rank,
		 "Rank confusion in _bz_doArrayIndexMapping");
        return expr.unwrap()(index[i0], index[i1], index[i2], index[i3], index[i4]);
    }

  template<int N_inputRank>
  static TinyVector<int,rank> map_dims(const TinyVector<int,N_inputRank>& index,
			 int i0, int i1, int i2, int i3, int i4,
			 int, int, int, int, int, int)
    {
      return TinyVector<int,rank>(index[i0], index[i1], index[i2],
				  index[i3], index[i4]);
    }
};

template<>
struct _bz_doArrayIndexMapping<6> {
  static const int rank=6;

  template<typename T_expr, int N_inputRank>
  static typename T_expr::T_numtype map(const ETBase<T_expr>& expr,
        const TinyVector<int,N_inputRank>& index, int i0, int i1, int i2,
        int i3, int i4, int i5, int, int, int, int, int)
    {
      // this was the case when it took an array. is it necessary?
      BZPRECHECK(T_expr::rank_==rank,
		 "Rank confusion in _bz_doArrayIndexMapping");
        return expr.unwrap()(index[i0], index[i1], index[i2], index[i3], index[i4],
            index[i5]);
    }

  template<int N_inputRank>
  static TinyVector<int,rank> map_dims(const TinyVector<int,N_inputRank>& index,
			 int i0, int i1, int i2, int i3, int i4,
			 int i5, int, int, int, int, int)
    {
      return TinyVector<int,rank>(index[i0], index[i1], index[i2],
				  index[i3], index[i4], index[i5]);
    }
};

template<>
struct _bz_doArrayIndexMapping<7> {
  static const int rank=7;

  template<typename T_expr, int N_inputRank>
  static typename T_expr::T_numtype map(const ETBase<T_expr>& expr,
        const TinyVector<int,N_inputRank>& index, int i0, int i1, int i2,
        int i3, int i4, int i5, int i6, int, int, int, int)
    {
      // this was the case when it took an array. is it necessary?
      BZPRECHECK(T_expr::rank_==rank,
		 "Rank confusion in _bz_doArrayIndexMapping");
        return expr.unwrap()(index[i0], index[i1], index[i2], index[i3],
			     index[i4], index[i5], index[i6]);
    }

  template<int N_inputRank>
  static TinyVector<int,rank> map_dims(const TinyVector<int,N_inputRank>& index,
			 int i0, int i1, int i2, int i3, int i4,
			 int i5, int i6, int, int, int, int)
    {
      return TinyVector<int,rank>(index[i0], index[i1], index[i2],
				  index[i3], index[i4], index[i5],
				  index[i6]);
    }
};

template<>
struct _bz_doArrayIndexMapping<8> {
  static const int rank=8;

  template<typename T_expr, int N_inputRank>
  static typename T_expr::T_numtype map(const ETBase<T_expr>& expr,
        const TinyVector<int,N_inputRank>& index, int i0, int i1, int i2,
        int i3, int i4, int i5, int i6, int i7, int, int, int)
    {
      // this was the case when it took an array. is it necessary?
      BZPRECHECK(T_expr::rank_==rank,
		 "Rank confusion in _bz_doArrayIndexMapping");
        return expr.unwrap()(index[i0], index[i1], index[i2], index[i3],
			     index[i4], index[i5], index[i6], index[i7]);
    }

  template<int N_inputRank>
  static TinyVector<int,rank> map_dims(const TinyVector<int,N_inputRank>& index,
			 int i0, int i1, int i2, int i3, int i4,
			 int i5, int i6, int i7, int, int, int)
    {
      return TinyVector<int,rank>(index[i0], index[i1], index[i2],
				  index[i3], index[i4], index[i5],
				  index[i6], index[i7]);
    }
};

template<>
struct _bz_doArrayIndexMapping<9> {
  static const int rank=9;

  template<typename T_expr, int N_inputRank>
  static typename T_expr::T_numtype map(const ETBase<T_expr>& expr,
        const TinyVector<int,N_inputRank>& index, int i0, int i1, int i2,
        int i3, int i4, int i5, int i6, int i7, int i8, int, int)
    {
      // this was the case when it took an array. is it necessary?
      BZPRECHECK(T_expr::rank_==rank,
		 "Rank confusion in _bz_doArrayIndexMapping");
        return expr.unwrap()(index[i0], index[i1], index[i2], index[i3],
			     index[i4], index[i5], index[i6], index[i7],
			     index[i8]);
    }

  template<int N_inputRank>
  static TinyVector<int,rank> map_dims(const TinyVector<int,N_inputRank>& index,
			 int i0, int i1, int i2, int i3, int i4,
			 int i5, int i6, int i7, int i8, int, int)
    {
      return TinyVector<int,rank>(index[i0], index[i1], index[i2],
				  index[i3], index[i4], index[i5],
				  index[i6], index[i7], index[i8]);
    }
};

template<>
struct _bz_doArrayIndexMapping<10> {
  static const int rank=10;

  template<typename T_expr, int N_inputRank>
  static typename T_expr::T_numtype map(const ETBase<T_expr>& expr,
        const TinyVector<int,N_inputRank>& index, int i0, int i1, int i2,
        int i3, int i4, int i5, int i6, int i7, int i8, int i9, int)
    {
      // this was the case when it took an array. is it necessary?
      BZPRECHECK(T_expr::rank_==rank,
		 "Rank confusion in _bz_doArrayIndexMapping");
        return expr.unwrap()(index[i0], index[i1], index[i2], index[i3],
			     index[i4], index[i5], index[i6], index[i7],
			     index[i8], index[i9]);
    }

  template<int N_inputRank>
  static TinyVector<int,rank> map_dims(const TinyVector<int,N_inputRank>& index,
			 int i0, int i1, int i2, int i3, int i4,
			 int i5, int i6, int i7, int i8, int i9, int)
    {
      return TinyVector<int,rank>(index[i0], index[i1], index[i2],
				  index[i3], index[i4], index[i5],
				  index[i6], index[i7], index[i8],
				  index[i9]);
    }
};

template<>
struct _bz_doArrayIndexMapping<11> {
  static const int rank=11;

  template<typename T_expr, int N_inputRank>
  static typename T_expr::T_numtype map(const ETBase<T_expr>& expr,
        const TinyVector<int,N_inputRank>& index, int i0, int i1, int i2,
        int i3, int i4, int i5, int i6, int i7, int i8, int i9, int i10)
    {
      // this was the case when it took an array. is it necessary?
      BZPRECHECK(T_expr::rank_==rank,
		 "Rank confusion in _bz_doArrayIndexMapping");
        return expr.unwrap()(index[i0], index[i1], index[i2], index[i3],
			     index[i4], index[i5], index[i6], index[i7],
			     index[i8], index[i9], index[i10]);
    }

  template<int N_inputRank>
  static TinyVector<int,rank> map_dims(const TinyVector<int,N_inputRank>& index,
			 int i0, int i1, int i2, int i3, int i4,
			 int i5, int i6, int i7, int i8, int i9, int i10)
    {
      return TinyVector<int,rank>(index[i0], index[i1], index[i2],
				  index[i3], index[i4], index[i5],
				  index[i6], index[i7], index[i8],
				  index[i9], index[i10]);
    }
};


// default arguments are defined in the fwd header
template<typename P_expr, int N_map0, int N_map1, int N_map2,
    int N_map3, int N_map4, int N_map5, int N_map6, int N_map7,
    int N_map8, int N_map9, int N_map10>
class ArrayIndexMapping {
public:
  typedef P_expr T_expr;
  typedef typename T_expr::T_numtype T_numtype;

  typedef T_numtype T_optype;
  typedef typename asET<T_numtype>::T_wrapped T_typeprop;
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;

  typedef T_expr T_ctorArg1;
    typedef int                            T_ctorArg2;    // dummy
  typedef ArrayIndexMapping<typename T_expr::T_range_result,N_map0,N_map1,N_map2,N_map3,N_map4,N_map5,N_map6,N_map7,N_map8,N_map9,N_map10> T_range_result;

    /*
     * This enum block finds the maximum of the N_map0, N_map1, ..., N_map10
     * parameters and stores it in maxRank10.  The rank of the expression is
     * then maxRank10 + 1, since the IndexPlaceholders start at 0 rather than
     * 1.
     */
    static const int
        maxRank1 = (N_map0 > N_map1) ? N_map0 : N_map1,
        maxRank2 = (N_map2 > maxRank1) ? N_map2 : maxRank1,
        maxRank3 = (N_map3 > maxRank2) ? N_map3 : maxRank2,
        maxRank4 = (N_map4 > maxRank3) ? N_map4 : maxRank3,
        maxRank5 = (N_map5 > maxRank4) ? N_map5 : maxRank4,
        maxRank6 = (N_map6 > maxRank5) ? N_map6 : maxRank5,
        maxRank7 = (N_map7 > maxRank6) ? N_map7 : maxRank6,
        maxRank8 = (N_map8 > maxRank7) ? N_map8 : maxRank7,
        maxRank9 = (N_map9 > maxRank8) ? N_map9 : maxRank8,
        maxRank10 = (N_map10 > maxRank9) ? N_map10 : maxRank9;

    static const int
    numArrayOperands = T_expr::numArrayOperands,
      numTVOperands = T_expr::numTVOperands,
      numTMOperands = T_expr::numTMOperands,
        numIndexPlaceholders = 1,
      minWidth = simdTypes<T_numtype>::vecWidth,
      maxWidth = simdTypes<T_numtype>::vecWidth,
        rank_ = maxRank10 + 1,
        exprRank = T_expr::rank_;

  template<int N> struct tvresult {
    typedef FastTV2Iterator<T_numtype, N> Type;
  };

  /*
    ArrayIndexMapping(const Array<T_numtype, rank>& array)
        : iter_(array)
    {
    }
  */

    ArrayIndexMapping(const ArrayIndexMapping<T_expr,N_map0,
        N_map1,N_map2,N_map3,N_map4,N_map5,N_map6,N_map7,N_map8,N_map9,
        N_map10>& z)
        : iter_(z.iter_)
    {
    }

  ArrayIndexMapping(BZ_ETPARM(T_expr) a)
        : iter_(a)
    { }

  // this is ambiguous with the above
  // ArrayIndexMapping(_bz_typename T_expr::T_ctorArg1 a)
  //       : iter_(a)
  //   { }

  // these bypass the FAI and go directly to the array. That should
  // prevent any performance impact of using the FAI instead of an
  // array directly.

  /* Functions for reading. Because they must depend on the result
   * type, they utilize a helper class.
   */

  // For numtypes, apply operator
  template<typename T> struct readHelper {
    static T_result first_value(const T_expr& iter)  {
      // is the correct thing to do here to return the index zero value?
      return indexop(iter, TinyVector<int,1>(0)); }
    template<int N_rank>
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
    static T_result indexop(const T_expr& iter,
			    const TinyVector<int, N_rank> i) {
#else
      static T_result indexop(const T_expr& iter,
			      const TinyVector<int, N_rank>& i) {
#endif
	return _bz_doArrayIndexMapping<exprRank>::map(iter/*.array()*/, i,
            N_map0, N_map1, N_map2, N_map3, N_map4, N_map5, N_map6,
            N_map7, N_map8, N_map9, N_map10);
  }
    };

  // For ET types, bypass operator and create expression
    template<typename T> struct readHelper<ETBase<T> > {
      template<int N_rank>
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
      static T_result indexop(const T_expr& iter,
			      const TinyVector<int, N_rank> i) {
#else
    static T_result indexop(const T_expr& iter,
			    const TinyVector<int, N_rank>& i) {
#endif
      return iter(i); }
    static T_result first_value(const T_expr& iter)  {
      // is the correct thing to do here to return the index zero value?
      return indexop(iter, TinyVector<int,1>(0)); }
    };

    template<int N_rank>
#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
    T_result operator()(const TinyVector<int, N_rank> i) const {
#else
      T_result operator()(const TinyVector<int, N_rank>& i) const {
#endif
      return readHelper<T_typeprop>::indexop(iter_,i); }

    T_result first_value() const {
      // unclear how to define "first" value for index expressions.
      BZPRECHECK(0,"Minmax reductions of index expressions not implemented");
      return readHelper<T_typeprop>::first_value(iter_); }

  // find which dimension in mapped expression that corresponds to
  // dimension dim. This works such that dimension dim in this
  // expression corresponds to dimension map_dim(dim) in iter_.
    int map_dim(const int dim) const
    {
        if (N_map0 == dim)
	  return 0;
        else if ((N_map1 == dim) && (exprRank > 1))
	  return 1;
        else if ((N_map2 == dim) && (exprRank > 2))
	  return 2;
        else if ((N_map3 == dim) && (exprRank > 3))
	  return 3;
        else if ((N_map4 == dim) && (exprRank > 4))
	  return 4;
        else if ((N_map5 == dim) && (exprRank > 5))
	  return 5;
        else if ((N_map6 == dim) && (exprRank > 6))
	  return 6;
        else if ((N_map7 == dim) && (exprRank > 7))
	  return 7;
        else if ((N_map8 == dim) && (exprRank > 8))
	  return 8;
        else if ((N_map9 == dim) && (exprRank > 9))
	  return 9;
        else if ((N_map10 == dim) && (exprRank > 10))
	  return 10;
        else
	  // means dimension is not in this operand
	  return -1;
    }

  // remaps the dimensions of an index vector so it can be applied to
  // iter_, using the _bz_doArrayIndexMapping helper class.
  template<int N>
  TinyVector<int, exprRank> map_dims(const TinyVector<int, N>& i) const {
    return _bz_doArrayIndexMapping<exprRank>::map_dims
      (i, N_map0, N_map1, N_map2, N_map3, N_map4, N_map5,
       N_map6, N_map7, N_map8, N_map9, N_map10);
  }

  int ascending(const int dim) const
  {
    const int d=map_dim(dim);
    const int o = d>=0 ? iter_.ascending(d) : INT_MIN ;
    return o;
  }

  int ordering(const int dim) const
  {
    // JCC: ignore ordering result from 1d Array
    if (exprRank == 1)
      return INT_MIN;  // tiny(int());

    const int d=map_dim(dim);
    const int o = d>=0 ? iter_.ordering(d) : INT_MIN ;
    return o;
  }

  int lbound(const int dim) const
  {
    const int d=map_dim(dim);
    const int o = d>=0 ? iter_.lbound(d) : INT_MIN ;
    return o;
  }

  int ubound(const int dim) const
  {
    const int d=map_dim(dim);
    const int o = d>=0 ? iter_.ubound(d) : INT_MAX ;
    return o;
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

    // If you have a precondition failure on this routine, it means
    // you are trying to use stack iteration mode on an expression
    // which contains an index placeholder.  You must use index
    // iteration mode instead.
  // (no -- added to support stencils /PJ)
    T_result operator*() const
    {
      return *iter_;
    }

    // See operator*() note
    void push(int)
    {
        BZPRECHECK(0,"Can't use stack iteration on an index mapping.");
    }

    // See operator*() note
    void pop(int)
    {
        BZPRECHECK(0,"Can't use stack iteration on an index mapping.");
    }

    // See operator*() note
    void advance()
    {
        BZPRECHECK(0,"Can't use stack iteration on an index mapping.");
    }

    // See operator*() note
    void advance(int)
    {
        BZPRECHECK(0,"Can't use stack iteration on an index mapping.");
    }

    // See operator*() note
    void loadStride(int)
    {
        BZPRECHECK(0,"Can't use stack iteration on an index mapping.");
    }

    bool isUnitStride(int) const
    {
        BZPRECHECK(0,"Can't use stack iteration on an index mapping.");
        return false;
    }

    bool isUnitStride() const
    {
        BZPRECHECK(0,"Can't use stack iteration on an index mapping.");
        return false;
    }

    void advanceUnitStride()
    {
        BZPRECHECK(0,"Can't use stack iteration on an index mapping.");
    }

    bool canCollapse(int,int) const
    {   BZPRECHECK(0,"Can't use stack iteration on an index mapping.");  return false; }

    T_result operator[](int)
    {
        BZPRECHECK(0,"Can't use stack iteration on an index mapping.");
        return T_result();
    }

    T_result fastRead(int) const
    {
        BZPRECHECK(0,"Can't use stack iteration on an index mapping.");
        return T_result();
    }

  template<int N>
  typename tvresult<N>::Type fastRead_tv(int) const {
    BZPRECHECK(0,"Can't use stack iteration on an index mapping.");
    return TinyVector<T_numtype, N>();
    }

    /** Determining whether the resulting expression is aligned is
	difficult, so to be safe we say no. It shouldn't be attempted
	anyway, though. */
    bool isVectorAligned(diffType offset) const {
      return false; }

    int suggestStride(int) const
    {
        BZPRECHECK(0,"Can't use stack iteration on an index mapping.");
        return 0;
    }

    bool isStride(int,int) const
    {
        BZPRECHECK(0,"Can't use stack iteration on an index mapping.");
        return true;
    }

#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
    template<int N_destrank>
    void moveTo(const TinyVector<int,N_destrank> i)
    {
      iter_.moveTo(map_dims(i));
    }
#else
    template<int N_destrank>
    void moveTo(const TinyVector<int,N_destrank>& i)
    {
      iter_.moveTo(map_dims(i));
    }
#endif

  T_result shift(int offset, int dim) const {
    // need to check if dim is mapped into this expression
    const int d=map_dim(dim);
    if (d<0)
      return *iter_;
    else
      return iter_.shift(offset, d);
  }

  T_result shift(int offset1, int dim1,int offset2, int dim2) const {
    // need to check if dims are mapped into this expression
    int d1=map_dim(dim1);
    int d2=map_dim(dim2);
    if (d1<0) //disable offset
      {d1=0;offset1=0;}
    if (d2<0) //disable offset
      {d2=0;offset2=0;}
    return iter_.shift(offset1, d1, offset2, d2);
  }

  void _bz_offsetData(sizeType i) {
        BZPRECHECK(0,"Can't use stack iteration on an index mapping.");
  }

  template<int N>
  T_range_result operator()(RectDomain<N> d) const
  { // need to reorder dimensions here
    TinyVector<int, exprRank> lb(map_dims(d.lbound())), ub(map_dims(d.ubound()));
    RectDomain<exprRank> newd(lb,ub);
    return T_range_result(iter_(newd));
  }

    void prettyPrint(BZ_STD_SCOPE(string) &str, prettyPrintFormat&) const
    {
        // NEEDS_WORK-- do real formatting for reductions
        str += "map[NEEDS_WORK]";
    }

    template<typename T_shape>
    bool shapeCheck(const T_shape&) const
    {
        // NEEDS_WORK-- do a real shape check (tricky)
        return true;
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
    typedef typename T_expr::template SliceInfo<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>::T_slice T_slice1;
    typedef ArrayIndexMapping<T_slice1, N_map0, N_map1, N_map2,
			      N_map3, N_map4, N_map5, N_map6, N_map7,
			      N_map8, N_map9, N_map10> T_slice;
};

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
        typename T7, typename T8, typename T9, typename T10, typename T11>
    typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const
    {
      /* Slicing for remapped expressions doesn't work. Because of the
	 potential different types (Range vs int) in the expression,
	 it would be very awkward to implement. As far as I can see,
	 it would require manual coding of the 3^11 calling
	 possibilities. /PJ */
      BZPRECONDITION(0);
    }

private:
    ArrayIndexMapping() : iter_( Array<T_numtype, exprRank>() ) { }

    T_expr iter_;
};

BZ_NAMESPACE_END

#endif // BZ_ARRAYMAP_H
