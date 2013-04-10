// -*- C++ -*-
/***************************************************************************
 * blitz/range.h      Declaration of the Range class
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

#ifndef BZ_RANGE_H
#define BZ_RANGE_H


#include <blitz/blitz.h>
#include <blitz/etbase.h>
#include <blitz/array/asexpr.h>
#include <blitz/prettyprint.h>
#include <blitz/tinyvec2.h>
#include <climits>                  // for INT_MIN, INT_MAX

BZ_NAMESPACE(blitz)

// A Range object is an ET that generates the specified sequence.

// Examples:
// Array<int,1> A(7);
// A = 0,1,2,3,4,5,6;
// A(Range::all());                   [0,1,2,3,4,5,6]
// A(Range(3,5));                     [3,4,5]
// A(Range(3,toEnd));                 [3,4,5,6]
// A(Range(fromStart,3));             [0,1,2,3]
// A(Range(1,5,2));                   [1,3,5]
// A(Range(5,1,-2));                  [5,3,1]
// A(Range(fromStart,toEnd,2));       [0,2,4,6]

template<int N_rank> class RectDomain;
class nilArraySection;

const int fromStart = INT_MIN;
const int toEnd = INT_MAX;

// Class Range
class Range : public ETBase<Range> {
public:
    typedef int T_numtype;
  typedef opType<T_numtype>::T_optype T_optype;

  typedef asET<T_numtype>::T_wrapped T_typeprop;
  typedef unwrapET<T_typeprop>::T_unwrapped T_result;

    typedef void T_ctorArg1;
    typedef char       T_ctorArg2;    // dummy
  typedef TinyVector<int, 1> T_index;
  typedef Range T_range_result;
    static const int
    numArrayOperands = 0,
        numTVOperands = 0,
        numTMOperands = 0,
      numIndexPlaceholders = 1,
      minWidth = simdTypes<T_numtype>::vecWidth,
      maxWidth = simdTypes<T_numtype>::vecWidth,
      rank_ = 1;

  /** The vectorized return type for a Range should be another range,
      but that's not useful since a vectorized TinyVector assignment
      can not contain index placeholders. In fact, since vectorization
      doesn't work for index expressions anyway, we can just set this
      to a dummy. */
  template<int N> struct tvresult {
    typedef FastTV2Iterator<T_numtype, N> Type;
  };

    Range()
    {
        first_ = fromStart;
        last_ = toEnd;
        stride_ = 1;
    }

    Range(const Range& r)
    {
        first_ = r.first_;
        last_ = r.last_;
        stride_ = r.stride_;
    }

    explicit Range(T_numtype slicePosition)
    {
        first_ = slicePosition;
        last_ = slicePosition;
        stride_ = 1;
    }

    Range(T_numtype first, T_numtype last, diffType stride=1)
        : first_(first), last_(last), stride_(stride)
    {
        BZPRECHECK((first == fromStart) || (last == toEnd) ||
            ((first < last) && (stride > 0)) ||
            ((first > last) && (stride < 0)) ||
            (first == last), (*this) << " is an invalid range.");
        BZPRECHECK((first == fromStart) || (last == toEnd) ||
            (last-first) % stride == 0,
            (*this) << ": the stride must evenly divide the range");
    }

  T_numtype operator*()   const { BZPRECONDITION(0); return 0; }
  T_numtype first_value() const { BZPRECONDITION(0); return 0; }

    int ascending(const int) const { return true; }
    int ordering(const int)  const { return 0; }
    int lbound(const int)    const { return 0; }
  int ubound(const int)    const { return length()-1; }

  RectDomain<rank_> domain() const;

  bool assertInRange(const T_index& BZ_DEBUG_PARAM(index)) const;

#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
    template<int N_rank>
    T_numtype operator()(const TinyVector<int,N_rank> i) const
  { assertInRange(i); return operator[](i[0]); }
#else
    template<int N_rank>
    T_numtype operator()(const TinyVector<int,N_rank>& i) const
  { assertInRange(i); return operator[](i[0]); }
#endif

  T_numtype operator[](int i) const;

    T_numtype operator()(int i) const
    {
      return operator[](i);
    }

  // we could work out how this should work.
  template<int N_rank>
  const Range operator()(const RectDomain<N_rank>& d) const
  {
    BZPRECONDITION(0); return *this;
  }

    void push(int) { }
    void pop(int) { }
    void advance() { }
    void advance(int) { }
    void loadStride(int) { }

    bool isUnitStride(int) const
  { BZPRECONDITION(0); return 0; }

    void advanceUnitStride()
    { }

    bool canCollapse(int,int) const
    { return true; }

    T_numtype fastRead(diffType) const
  { BZPRECONDITION(0); return 0; }

  template<int N>
  typename tvresult<N>::Type fastRead_tv(diffType) const
  { BZPRECONDITION(0); return 0; }

  // this is needed for the stencil expression fastRead to work
  void _bz_offsetData(sizeType i) const{BZPRECONDITION(0);};

    // and these are needed for stencil expression shift to work
  void _bz_offsetData(sizeType offset, int dim) const {BZPRECONDITION(0);};

  void _bz_offsetData(sizeType offset1, int dim1, sizeType offset2, int dim2) const {BZPRECONDITION(0);};

    diffType suggestStride(int) const
    { return 1; }

    bool isStride(int,diffType) const
    { return true; }

  void moveTo(int) const { BZPRECONDITION(0); }

  T_numtype shift(int offset, int dim) const {
    BZPRECONDITION(0); return T_numtype(); }

    T_numtype shift(int offset1, int dim1,int offset2, int dim2) const
  { BZPRECONDITION(0); return T_numtype(); }

    template<int N_rank>
    void moveTo(const TinyVector<int,N_rank>&) const { BZPRECONDITION(0); }

    void prettyPrint(BZ_STD_SCOPE(string) &str,
        prettyPrintFormat& format) const
    {
#ifdef BZ_HAVE_STD
	    BZ_STD_SCOPE(ostringstream) ostr;
#else
            ostrstream ostr;
#endif
	    ostr << "Range(" << first_ << ", " << last_
		 << ", " << stride_ << ")";
            str += ostr.str();
    }


  // old range stuff below. what is needed for the range interface and
  // what is old vecexpr stuff?

    T_numtype first(T_numtype lowRange = 0) const
    {
        if (first_ == fromStart)
            return lowRange;
        return first_;
    }

    T_numtype last(T_numtype highRange = 0) const
    {
        if (last_ == toEnd)
            return highRange;
        return last_;
    }

  int length(int =0) const
    {
        BZPRECONDITION(first_ != fromStart);
        BZPRECONDITION(last_ != toEnd);
        BZPRECONDITION((last_ - first_) % stride_ == 0);
        return (last_ - first_) / stride_ + 1;
    }

    diffType stride() const
    { return stride_; }

    bool isAscendingContiguous() const
    {
        return (((first_ < last_) && (stride_ == 1)) || (first_ == last_));
    }

    void setRange(T_numtype first, T_numtype last, diffType stride=1)
    {
        BZPRECONDITION(((first < last) && (stride > 0)) ||
            ((first > last) && (stride < 0)) ||
            (first == last));
        BZPRECONDITION((last-first) % stride == 0);
        first_ = first;
        last_ = last;
        stride_ = stride;
    }

    static Range all()
    { return Range(fromStart,toEnd,1); }

  /// \todo this talks about the stride of the RANGE, not the expression stride.
    bool isUnitStride() const
    { return stride_ == 1; }

    // Operators
    Range operator-(T_numtype shift) const
    {
        BZPRECONDITION(first_ != fromStart);
        BZPRECONDITION(last_ != toEnd);
        return Range(first_ - shift, last_ - shift, stride_);
    }

    Range operator+(T_numtype shift) const
    {
        BZPRECONDITION(first_ != fromStart);
        BZPRECONDITION(last_ != toEnd);
        return Range(first_ + shift, last_ + shift, stride_);
    }

    friend inline ostream& operator<<(ostream& os, const Range& range)
    {
        os << "Range(" << range.first() << "," << range.last() << ","
           << range.stride() << ")";

        return os;
    }

  // we can't reduce the rank of a range, so we can't slice it
  template<typename T1, typename T2 = nilArraySection,
	   class T3 = nilArraySection, typename T4 = nilArraySection,
	   class T5 = nilArraySection, typename T6 = nilArraySection,
	   class T7 = nilArraySection, typename T8 = nilArraySection,
	   class T9 = nilArraySection, typename T10 = nilArraySection,
	   class T11 = nilArraySection>
  class SliceInfo {
  public:
    typedef void T_slice;
  };

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
        typename T7, typename T8, typename T9, typename T10, typename T11>
    typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const
    {
      return *this;
    }

private:
  T_numtype first_, last_;
  diffType stride_;
};

BZ_NAMESPACE_END

#endif // BZ_RANGE_H
