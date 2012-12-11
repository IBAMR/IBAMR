// -*- C++ -*-
/***************************************************************************
 * blitz/tinymat.h       Declaration of TinyMatrix<T, N, M>
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

#ifndef BZ_TINYMAT_H
#define BZ_TINYMAT_H

#ifdef __GNUC__
#if (__GNUC_MINOR__ >= 6)
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#endif

#include <blitz/blitz.h>
#include <blitz/tinyvec2.h>
#include <blitz/tinyvec2io.cc> // fix this - for some reason it doesn't get included through tinyvec2.h
#include <blitz/listinit.h>
#include <blitz/etbase.h>
#include <blitz/array/slice.h>
#include <blitz/indexmap-forward.h>
//#include <blitz/meta/vecassign.h>
//#include <blitz/update.h>

#ifdef BZ_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/serialization.hpp>
#endif


BZ_NAMESPACE(blitz)

/*****************************************************************************
 * Forward declarations
 */

template<typename P_numtype, int N_rows, int N_columns>
class FastTM2Iterator;
template<typename P_numtype, int N_rows, int N_columns>
class FastTM2CopyIterator;

template<typename P_expr>
class _bz_ArrayExpr;

template<int N0>
class IndexPlaceholder;

/*****************************************************************************
 * Declaration of class TinyMatrix
 */

template<typename P_numtype, int N_rows, int N_columns>
class TinyMatrix : public ETBase<TinyMatrix<P_numtype, N_rows, N_columns> >
{
public:
    typedef P_numtype T_numtype;
    // typedef _bz_tinyMatrixRef<T_numtype, N_rows, N_columns, N_columns, 1>
    //     T_reference;
    typedef TinyVector<int, 2>  T_index;
    typedef TinyMatrix<T_numtype, N_rows, N_columns> T_matrix;
  typedef FastTM2Iterator<T_numtype,N_rows, N_columns> T_iterator;
    typedef T_numtype*                                   iterator;
    typedef const T_numtype*                             const_iterator;
  typedef FastTM2CopyIterator<P_numtype, N_rows, N_columns> T_range_result;

    static const int
    //numArrayOperands = 1,
    //numIndexPlaceholders = 0,
        rank_ = 2;

    TinyMatrix() { }

  inline TinyMatrix(const T_matrix&);

  template <typename T_numtype2>
  inline TinyMatrix(const TinyMatrix<T_numtype2, N_rows, N_columns>& x);

  inline TinyMatrix(T_numtype initValue);

  static TinyVector<int, 2>    base()
  { return 0; }

  static int                               base(int rank)
  { return 0; }

    iterator                          begin()
    { return iterator(*this); }

    const_iterator                    begin() const
    { return const_iterator(*this); }

    T_iterator                        beginFast() const
    { return T_iterator(*this); }

  static int                               cols()
  { return N_columns; }

  static int                               columns()
  { return N_columns; }

    const T_numtype* restrict         data() const
    { return data_; }

    T_numtype* restrict               data()
    { return data_; }

    const T_numtype* restrict         dataFirst() const
    { return data_; }

    T_numtype* restrict               dataFirst()
    { return data_; }

  static int                               dimensions()
    { return 2; }

  static RectDomain<2>                domain();

    iterator                          end()
  { return data_+size(); }

    const_iterator                    end() const
  { return data_+size(); }

  static int                               extent(int rank)
  { BZPRECONDITION(rank<2); return rank==0 ? N_rows : N_columns; }

  static const TinyVector<int,2>     extent()
  { return TinyVector<int,2>(N_rows, N_columns); }

  static bool                              isMajorRank(int rank)
  { BZPRECONDITION(rank<2); return rank==0; }
  static bool                              isMinorRank(int rank)
  { BZPRECONDITION(rank<2); return rank==1; }
  static bool                              isRankStoredAscending(int rank)
  { BZPRECONDITION(rank<2); return true; }

  static int                               lbound(int rank)
  { BZPRECONDITION(rank<2); return 0; }
  static TinyVector<int,2>            lbound()
  { return 0; }

  static int                               length(int rank)
  { BZPRECONDITION(rank<2); return rank==0 ? N_rows : N_columns; }
  static TinyVector<int, 2>    length()
  { return TinyVector<int,2>(N_rows, N_columns); }

  static sizeType                               numElements()
  { return N_rows*N_columns; }

  static int                               ordering(int rank)
  { BZPRECONDITION(rank<2); return 1-rank; }

  static TinyVector<int, 2>    ordering()
  { return TinyVector<int,2>(0,1); }

  static int                               rank()
    { return rank_; }

  static int                               rows()
    { return N_rows; }

  static TinyVector<int, 2>    shape()
  { return length(); }

  static sizeType                               size()
    { return numElements(); }

  static TinyVector<diffType, 2>    stride()
  { return TinyVector<int,2>(N_columns,1); }

  static diffType                               stride(int rank)
  { BZPRECONDITION(rank<2); return rank==0 ? N_columns : 1; }

  static int                               ubound(int rank)
    { return length(rank) - 1; }

  static TinyVector<int, 2>           ubound()
  { return length()-1; }

  static int                               zeroOffset()
    { return 0; }

    //////////////////////////////////////////////
    // Debugging routines
    //////////////////////////////////////////////

    bool isInRange(const T_index& index) const {
        for (int i=0; i < rank_; ++i)
	  if (index[i] < base(i) || (index[i] - base(i)) >= length(i))
                return false;

        return true;
    }

    bool isInRange(int i0, int i1) const {
      return i0 >= base(0) && (i0 - base(0)) < length(0)
	&& i1 >= base(1) && (i1 - base(1)) < length(1);
    }

    bool assertInRange(const T_index& BZ_DEBUG_PARAM(index)) const {
        BZPRECHECK(isInRange(index), "TinyMatrix index out of range: " << index
            << endl << "Lower bounds: " << base() << endl
		   <<         "Length:       " << length() << endl);
        return true;
    }

    bool assertInRange(int BZ_DEBUG_PARAM(i0), int BZ_DEBUG_PARAM(i1)) const {
        BZPRECHECK(isInRange(i0,i1), "TinyMatrix index out of range: ("
            << i0 << ", " << i1 << ")"
            << endl << "Lower bounds: " << base() << endl
		   <<         "Length:       " << length() << endl);
        return true;
    }

    //////////////////////////////////////////////
    // Subscripting operators
    //////////////////////////////////////////////


    T_matrix& noConst() const
    { return const_cast<T_matrix&>(*this); }

  const T_numtype& restrict operator()(unsigned i, unsigned j) const
    {
      assertInRange(i,j);
      return data_[i*N_columns + j];
    }

  T_numtype& restrict operator()(unsigned i, unsigned j)
    {
      assertInRange(i,j);
      return data_[i*N_columns + j];
    }

  const T_numtype& restrict operator()(T_index i) const
    {
      assertInRange(i);
      return data_[i[0]*N_columns + i[1]];
    }

  T_numtype& restrict operator()(T_index i)
    {
      assertInRange(i);
      return data_[i[0]*N_columns + i[1]];
    }

    template<int N0, int N1>
    inline
    _bz_ArrayExpr<ArrayIndexMapping<typename asExpr<T_matrix>::T_expr, N0, N1> >
    operator()(IndexPlaceholder<N0>, IndexPlaceholder<N1>) const;

  /** Fastread must return reference so the iterator can turn it into an
      iterator for the contained in case we have a multicomponent. */
  const T_numtype& fastRead(diffType i) const
  { return data_[i]; }

  /** Since data_ is simd aligned by construction, we just have
      to check the offest. */
  bool isVectorAligned(diffType offset) const
  { return (offset%simdTypes<T_numtype>::vecWidth)==0; }

    // T_reference getRef()
    // { return T_reference((T_numtype*)data_); }

    // const T_reference getRef() const
    // { return T_reference((T_numtype*)data_); }

    //////////////////////////////////////////////
    // Assignment Operators
    //////////////////////////////////////////////

    // Scalar operand
    ListInitializationSwitch<T_matrix,T_numtype*>
    operator=(T_numtype x)
    {
        return ListInitializationSwitch<T_matrix,T_numtype*>(*this, x);
    }

  T_matrix& initialize(T_numtype x);

    T_numtype* restrict getInitializationIterator()
    { return dataFirst(); }

    template<typename T_expr>
    T_matrix& operator=(const ETBase<T_expr>&);

    template<typename T> T_matrix& operator+=(const T&);
    template<typename T> T_matrix& operator-=(const T&);
    template<typename T> T_matrix& operator*=(const T&);
    template<typename T> T_matrix& operator/=(const T&);
    template<typename T> T_matrix& operator%=(const T&);
    template<typename T> T_matrix& operator^=(const T&);
    template<typename T> T_matrix& operator&=(const T&);
    template<typename T> T_matrix& operator|=(const T&);
    template<typename T> T_matrix& operator>>=(const T&);
    template<typename T> T_matrix& operator<<=(const T&);

  static bool canCollapse(int outerRank, int innerRank) {
#ifdef BZ_DEBUG_TRAVERSE
        BZ_DEBUG_MESSAGE("stride(" << innerRank << ")=" << stride(innerRank)
          << ", extent()=" << extent(innerRank) << ", stride(outerRank)="
          << stride(outerRank));
#endif
        return (stride(innerRank) * extent(innerRank) == stride(outerRank));
    }


private:
  template<typename T_expr, typename T_update>
  void _tm_evaluate(const T_expr& expr, T_update);

#ifdef BZ_HAVE_BOOST_SERIALIZATION
  friend class boost::serialization::access;

  template<class T_arch>
  void serialize(T_arch& ar, const unsigned int version) {
    ar & data_;
  };
#endif

  BZ_ALIGN_VARIABLE(T_numtype, data_[N_rows * N_columns], BZ_SIMD_WIDTH)
};

BZ_NAMESPACE_END

#ifdef BZ_HAVE_BOOST_SERIALIZATION
namespace boost {
  namespace mpi {
    template<typename T> struct is_mpi_datatype;
    template <typename T, int N, int M>
    struct is_mpi_datatype<blitz::TinyMatrix<T, N, M> >
      : public is_mpi_datatype<T> { };
  } };
#endif

#include <blitz/tm2fastiter.h>  // Iterators
#include <blitz/tinymat2.cc>

#if (__GNUC__ && (__GNUC_MINOR__ >= 6))
#pragma GCC diagnostic pop
#endif

#endif // BZ_TINYMAT_H
