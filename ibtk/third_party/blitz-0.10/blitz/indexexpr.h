// -*- C++ -*-
/***************************************************************************
 * blitz/indexexpr.h     Declaration of the IndexPlaceholder<N> class
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

#ifndef BZ_INDEXEXPR_H
#define BZ_INDEXEXPR_H

#include <blitz/tinyvec2.h>
#include <blitz/prettyprint.h>
#include <blitz/etbase.h>
#include <blitz/array/domain.h>
#include <blitz/array/slice.h>

BZ_NAMESPACE(blitz)

template<int N>
class IndexPlaceholder
#ifdef BZ_NEW_EXPRESSION_TEMPLATES
  : public ETBase<IndexPlaceholder<N> >
#endif
{
public:
    IndexPlaceholder()
    { }

#ifdef BZ_NEW_EXPRESSION_TEMPLATES
    IndexPlaceholder(const IndexPlaceholder<N>& x)
        : ETBase< IndexPlaceholder<N> >(x)
    { }
#else
    IndexPlaceholder(const IndexPlaceholder<N>&)
    { }
#endif

    ~IndexPlaceholder()
    { }

    void operator=(const IndexPlaceholder<N>&)
    { }

    typedef int T_numtype;
    typedef int T_ctorArg1;     // Dummy; not used
    typedef int T_ctorArg2;     // Ditto
  typedef int T_range_result; // dummy

  typedef typename opType<T_numtype>::T_optype T_optype;
  typedef typename asET<T_numtype>::T_wrapped T_typeprop;
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;

    static const int
        numArrayOperands = 0,
        numTVOperands = 0,
        numTMOperands = 0,
        numIndexPlaceholders = 1,
      minWidth = simdTypes<T_numtype>::vecWidth,
      maxWidth = simdTypes<T_numtype>::vecWidth,
        rank_ = N+1;

  /** The vectorized return type for an IndexPlaceholder should be
      some form of range, but that's not useful since a vectorized
      TinyVector assignment can not contain index placeholders. In
      fact, since vectorization doesn't work for index expressions
      anyway, we can just set this to a dummy. */
  template<int M> struct tvresult {
    typedef FastTV2Iterator<T_numtype, M> Type;
  };

    // If you have a precondition failure on this routine, it means
    // you are trying to use stack iteration mode on an expression
    // which contains an index placeholder.  You must use index
    // iteration mode instead.
    int operator*() const {
        BZPRECONDITION(0);
        return 0;
    }

#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
    template<int N_rank>
    T_result operator()(TinyVector<int, N_rank> i) const { return i[N]; }
#else
    template<int N_rank>
    T_result operator()(const TinyVector<int, N_rank>& i) const { return i[N]; }
#endif

    int ascending(int) const { return INT_MIN; }
    int ordering(int)  const { return INT_MIN; }
    int lbound(int)    const { return INT_MIN; }  // tiny(int());
    int ubound(int)    const { return INT_MAX; }  // huge(int());

  RectDomain<rank_> domain() const
  {
    const TinyVector<int, rank_> lb(lbound(0)), ub(ubound(0));
    return RectDomain<rank_>(lb,ub);
  }

    // See operator*() note

    void push(int)       { BZPRECONDITION(0); }
    void pop(int)        { BZPRECONDITION(0); }
    void advance()       { BZPRECONDITION(0); }
    void advance(int)    { BZPRECONDITION(0); }
    void loadStride(int) { BZPRECONDITION(0); }
    template<int N_rank>
    void moveTo(const TinyVector<int,N_rank>& i) { BZPRECONDITION(0); }

    bool isUnitStride(int) const {
        BZPRECONDITION(0);
        return false;
    }

    bool isUnitStride() const {
        BZPRECONDITION(0);
        return false;
    }

    void advanceUnitStride() { BZPRECONDITION(0); }

    bool canCollapse(int,int) const {
        BZPRECONDITION(0);
        return false;
    }

    T_result operator[](int) const {
        BZPRECONDITION(0);
        return T_numtype();
    }

    T_result fastRead(diffType) const {
        BZPRECONDITION(0);
        return T_numtype();
    }

  template<int M>
  typename tvresult<M>::Type fastRead_tv(diffType) const {
    BZPRECONDITION(0);
    return TinyVector<T_numtype, M>();
    }

  /** There are no alignment issues here, so just return true. */
  bool isVectorAligned(diffType offset) const {
    (void)offset;
    return true; }

    diffType suggestStride(int) const {
        BZPRECONDITION(0);
        return 0;
    }

    bool isStride(int,diffType) const {
        BZPRECONDITION(0);
        return true;
    }

  // don't know how to define shift, as it relies on having an
  // implicit position. thus stencils won't work
  T_result shift(int offset, int dim) const
  { BZPRECHECK(0,"Stencils of index expressions are not implemented");
    return T_numtype(); }
  T_result shift(int offset1, int dim1,int offset2, int dim2) const {
    BZPRECHECK(0,"Stencils of index expressions are not implemented");
    return T_numtype(); }
  void _bz_offsetData(sizeType i) { BZPRECONDITION(0); }

  // Unclear how to define this, and stencils don't work anyway
  T_range_result operator()(RectDomain<rank_> d) const
  { BZPRECHECK(0,"Stencils of index expressions are not implemented");
    return T_range_result(); }

    void prettyPrint(BZ_STD_SCOPE(string) &str, prettyPrintFormat&) const {
        // NEEDS_WORK-- do real formatting for reductions
        str += "index-expr[NEEDS_WORK]";
    }

    template<typename T_shape>
    bool shapeCheck(const T_shape&) const { return true; }

  // sliceinfo for index placeholder does nothing
  template<typename T1, typename T2 = nilArraySection,
	   class T3 = nilArraySection, typename T4 = nilArraySection,
	   class T5 = nilArraySection, typename T6 = nilArraySection,
	   class T7 = nilArraySection, typename T8 = nilArraySection,
	   class T9 = nilArraySection, typename T10 = nilArraySection,
	   class T11 = nilArraySection>
  class SliceInfo {
  public:
    static const int new_rank = (N>0) ? ArraySectionInfo<T1>::rank : 0
      + (N>1) ? ArraySectionInfo<T2>::rank : 0
      + (N>2) ? ArraySectionInfo<T3>::rank : 0
      + (N>4) ? ArraySectionInfo<T4>::rank : 0
      + (N>5) ? ArraySectionInfo<T5>::rank : 0
      + (N>6) ? ArraySectionInfo<T6>::rank : 0
      + (N>7) ? ArraySectionInfo<T7>::rank : 0
      + (N>8) ? ArraySectionInfo<T8>::rank : 0
      + (N>9) ? ArraySectionInfo<T9>::rank : 0
      + (N>10) ? ArraySectionInfo<T10>::rank : 0
      + (N>11) ? ArraySectionInfo<T11>::rank : 0;
    typedef IndexPlaceholder<new_rank> T_slice;
  };

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
        typename T7, typename T8, typename T9, typename T10, typename T11>
    typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const
    {
      // slicing of index placeholders is not implemented. there are
      // two problems: First, if you slice the dimension of the index
      // placeholder, it should be replaced with a constant. Second,
      // if you restrict the range of that dimension, the index
      // placeholder should start from a nonzero value.
      BZPRECONDITION(0);
    }
};

typedef IndexPlaceholder<0> firstIndex;
typedef IndexPlaceholder<1> secondIndex;
typedef IndexPlaceholder<2> thirdIndex;
typedef IndexPlaceholder<3> fourthIndex;
typedef IndexPlaceholder<4> fifthIndex;
typedef IndexPlaceholder<5> sixthIndex;
typedef IndexPlaceholder<6> seventhIndex;
typedef IndexPlaceholder<7> eighthIndex;
typedef IndexPlaceholder<8> ninthIndex;
typedef IndexPlaceholder<9> tenthIndex;
typedef IndexPlaceholder<10> eleventhIndex;

#ifndef BZ_NO_TENSOR_INDEX_OBJECTS

BZ_NAMESPACE(tensor)
    _bz_global blitz::IndexPlaceholder<0> i;
    _bz_global blitz::IndexPlaceholder<1> j;
    _bz_global blitz::IndexPlaceholder<2> k;
    _bz_global blitz::IndexPlaceholder<3> l;
    _bz_global blitz::IndexPlaceholder<4> m;
    _bz_global blitz::IndexPlaceholder<5> n;
    _bz_global blitz::IndexPlaceholder<6> o;
    _bz_global blitz::IndexPlaceholder<7> p;
    _bz_global blitz::IndexPlaceholder<8> q;
    _bz_global blitz::IndexPlaceholder<9> r;
    _bz_global blitz::IndexPlaceholder<10> s;
    _bz_global blitz::IndexPlaceholder<11> t;
BZ_NAMESPACE_END // tensor

#endif

BZ_NAMESPACE_END

#endif // BZ_INDEXEXPR_H
