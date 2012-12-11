// -*- C++ -*-
/***************************************************************************
 * blitz/array/storage.h  Memory layout of Arrays.
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
#ifndef BZ_ARRAY_STORAGE_H
#define BZ_ARRAY_STORAGE_H

#ifdef BZ_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/serialization.hpp>
#endif
#ifdef BZ_HAVE_BOOST_MPI
#include <boost/mpi/datatype.hpp>
#endif

BZ_NAMESPACE(blitz)


/// Enum for specifying whether padding for alignment should be done.
enum paddingPolicy {
  /// The container should have contiguous data.
  contiguousData,
  /** The container data should be padded so minor rank dimensions
      always are SIMD aligned. */
  paddedData
};


/** The default padding policy, set by the configure script
    --enable-array-length-padding argument. */
const paddingPolicy defaultPadding = BZ_PADDING_DEFAULT;


/**
  Declaration of class GeneralArrayStorage<N_rank>

  This class describes a storage format for an N-dimensional array.
  The dimensions can be stored in an arbitrary order (for example, as
  a C-style row major array or Fortran-style column major array, or
  something else entirely).  Each dimension can be stored in either
  ascending (the most common) or descending order.  Each dimension can
  have its own base (starting index value: e.g. 0 for C-style arrays,
  1 for Fortran arrays).  The GeneralArrayStorage class also
  determines the padding policy when SIMD width is set.

  GeneralArrayStorage<N> defaults to C-style arrays.  To implement
  other storage formats, subclass and modify the constructor.  The
  class FortranArray, below, is an example.

  Objects inheriting from GeneralArrayStorage<N> can be passed as
  an optional constructor argument to Array objects.
  e.g. Array<int,3> A(16,16,16, FortranArray<3>());
  will create a 3-dimensional 16x16x16 Fortran-style array.
 */

template<int N_rank>
class GeneralArrayStorage {
public:
    class noInitializeFlag { };

    GeneralArrayStorage(noInitializeFlag)
    { }

    GeneralArrayStorage(paddingPolicy pp = defaultPadding)
    {
        for (int i=0; i < N_rank; ++i)
          ordering_(i) = N_rank - 1 - i;
        ascendingFlag_ = true;
        base_ = 0;
	paddingPolicy_ = pp;
    }

    GeneralArrayStorage(const GeneralArrayStorage<N_rank>& x)
      : paddingPolicy_(x.paddingPolicy_), ascendingFlag_(x.ascendingFlag_),
	ordering_(x.ordering_), base_(x.base_)
    {
    }

    GeneralArrayStorage(TinyVector<int,N_rank> ordering,
			TinyVector<bool,N_rank> ascendingFlag,
			paddingPolicy pp = defaultPadding)
      : paddingPolicy_(pp), ascendingFlag_(ascendingFlag),
	ordering_(ordering)
    {
        base_ = 0;
    }

    ~GeneralArrayStorage()
    { }

    GeneralArrayStorage<N_rank>& operator=(
        const GeneralArrayStorage<N_rank>& rhs)
    {
        ordering_ = rhs.ordering();
        ascendingFlag_ = rhs.ascendingFlag();
        base_ = rhs.base();
	paddingPolicy_ = rhs.paddingPolicy_;
        return *this;
    }

    TinyVector<int, N_rank>& ordering()
    { return ordering_; }

    const TinyVector<int, N_rank>& ordering() const
    { return ordering_; }

    int ordering(int i) const
    { return ordering_[i]; }

    void setOrdering(int i, int order)
    { ordering_[i] = order; }

    bool allRanksStoredAscending() const
    {
        for (int i=0; i < N_rank; ++i)
            if (!ascendingFlag_[i]) return false;
        return true;
    }

    bool isRankStoredAscending(int i) const
    { return ascendingFlag_[i]; }

    TinyVector<bool, N_rank>& ascendingFlag()
    { return ascendingFlag_; }

    const TinyVector<bool, N_rank>& ascendingFlag() const
    { return ascendingFlag_; }

    void setAscendingFlag(int i, bool ascendingFlag)
    { ascendingFlag_[i] = ascendingFlag; }

    TinyVector<int, N_rank>& base()
    { return base_; }

    const TinyVector<int, N_rank>& base() const
    { return base_; }

    int base(int i) const
    { return base_[i]; }

    void setBase(int i, int base)
    { base_[i] = base; }

    void setBase(const TinyVector<int, N_rank>& base)
    { base_ = base; }

  const paddingPolicy& padding() const { return paddingPolicy_; }

private:
#ifdef BZ_HAVE_BOOST_SERIALIZATION
  friend class boost::serialization::access;

    template<class T_arch>
    void serialize(T_arch& ar, const unsigned int version) {
      ar & paddingPolicy_;
      ar & ascendingFlag_;
      ar & ordering_;
      ar & base_;
    };
#endif

protected:
    /*
     * ordering_[] specifies the order in which the array is stored in
     * memory.  For a newly allocated array, ordering_(0) will give the
     * rank with unit stride, and ordering_(N_rank-1) will be the rank
     * with largest stride.  An order like [2, 1, 0] corresponds to
     * C-style array storage; an order like [0, 1, 2] corresponds to
     * Fortran array storage.
     *
     * ascendingFlag_[] indicates whether the data in a rank is stored
     * in ascending or descending order.  Most of the time these values
     * will all be true (indicating ascending order).  Some peculiar
     * formats (e.g. MS-Windows BMP image format) store the data in
     * descending order.
     *
     * base_[] gives the first valid index for each rank.  For a C-style
     * array, all the base_ elements will be zero; for a Fortran-style
     * array, they will be one.  base_[] can be set arbitrarily using
     * the Array constructor which takes a Range argument, e.g.
     * Array<float,2> A(Range(30,40),Range(23,33));
     * will create an array with base_[] = { 30, 23 }.
     */

    // declare the tinyvector<bool> first, because it doesn't have alignment
    paddingPolicy paddingPolicy_;
    TinyVector<bool, N_rank> ascendingFlag_;
    TinyVector<int,  N_rank> ordering_;
    TinyVector<int,  N_rank> base_;
};


/** This tag class can be used to provide a nicer notation for
    constructing padded arrays: instead of
      Array<int,2> A(3, 3, GeneralArrayStorage<2>(paddedData));
    one can simply write:
      Array<int,2> A(3, 3, paddedArray);
    where paddedArray is an object of type _bz_paddedTag.
*/
class _bz_paddedTag {
public:
    operator GeneralArrayStorage<1>()
    { return GeneralArrayStorage<1>(paddedData); }

    operator GeneralArrayStorage<2>()
    { return GeneralArrayStorage<2>(paddedData); }

    operator GeneralArrayStorage<3>()
    { return GeneralArrayStorage<3>(paddedData); }

    operator GeneralArrayStorage<4>()
    { return GeneralArrayStorage<4>(paddedData); }

    operator GeneralArrayStorage<5>()
    { return GeneralArrayStorage<5>(paddedData); }

    operator GeneralArrayStorage<6>()
    { return GeneralArrayStorage<6>(paddedData); }

    operator GeneralArrayStorage<7>()
    { return GeneralArrayStorage<7>(paddedData); }

    operator GeneralArrayStorage<8>()
    { return GeneralArrayStorage<8>(paddedData); }

    operator GeneralArrayStorage<9>()
    { return GeneralArrayStorage<9>(paddedData); }

    operator GeneralArrayStorage<10>()
    { return GeneralArrayStorage<10>(paddedData); }

    operator GeneralArrayStorage<11>()
    { return GeneralArrayStorage<11>(paddedData); }
};

// A global instance of this class will be placed in
// the blitz library (libblitz.a on unix machines).

_bz_global _bz_paddedTag paddedArray;


/** This tag class can be used to provide a nicer notation for
    constructing unpadded arrays: instead of
      Array<int,2> A(3, 3, GeneralArrayStorage<2>(contiguousData));
    one can simply write:
      Array<int,2> A(3, 3, contiguousArray);
    where contiguousArray is an object of type _bz_contiguousTag.
*/
class _bz_contiguousTag {
public:
    operator GeneralArrayStorage<1>()
    { return GeneralArrayStorage<1>(contiguousData); }

    operator GeneralArrayStorage<2>()
    { return GeneralArrayStorage<2>(contiguousData); }

    operator GeneralArrayStorage<3>()
    { return GeneralArrayStorage<3>(contiguousData); }

    operator GeneralArrayStorage<4>()
    { return GeneralArrayStorage<4>(contiguousData); }

    operator GeneralArrayStorage<5>()
    { return GeneralArrayStorage<5>(contiguousData); }

    operator GeneralArrayStorage<6>()
    { return GeneralArrayStorage<6>(contiguousData); }

    operator GeneralArrayStorage<7>()
    { return GeneralArrayStorage<7>(contiguousData); }

    operator GeneralArrayStorage<8>()
    { return GeneralArrayStorage<8>(contiguousData); }

    operator GeneralArrayStorage<9>()
    { return GeneralArrayStorage<9>(contiguousData); }

    operator GeneralArrayStorage<10>()
    { return GeneralArrayStorage<10>(contiguousData); }

    operator GeneralArrayStorage<11>()
    { return GeneralArrayStorage<11>(contiguousData); }
};

// A global instance of this class will be placed in
// the blitz library (libblitz.a on unix machines).

_bz_global _bz_contiguousTag contiguousArray;


/**
   Class FortranArray specializes GeneralArrayStorage to provide
   Fortran style arrays (column major ordering, base of 1).  The
   noInitializeFlag() passed to the base constructor indicates that
   the subclass will take care of initializing the ordering_,
   ascendingFlag_ and base_ members.
 */
template<int N_rank>
class FortranArray : public GeneralArrayStorage<N_rank> {
private:
    typedef GeneralArrayStorage<N_rank> T_base;
    typedef _bz_typename T_base::noInitializeFlag noInitializeFlag;
    using T_base::ordering_;
    using T_base::ascendingFlag_;
    using T_base::base_;
    using T_base::paddingPolicy_;
public:
    FortranArray(paddingPolicy pp = defaultPadding)
        : GeneralArrayStorage<N_rank>(noInitializeFlag())
    {
        for (int i=0; i < N_rank; ++i)
          ordering_(i) = i;
        ascendingFlag_ = true;
        base_ = 1;
	paddingPolicy_ = pp;
    }
};


// This tag class can be used to provide a nicer notation for
// constructing Fortran-style arrays: instead of
//     Array<int,2> A(3, 3, FortranArray<2>());
// one can simply write:
//     Array<int,2> A(3, 3, fortranArray);
// where fortranArray is an object of type _bz_fortranTag.

class _bz_fortranTag {
public:
    operator GeneralArrayStorage<1>()
    { return FortranArray<1>(); }

    operator GeneralArrayStorage<2>()
    { return FortranArray<2>(); }

    operator GeneralArrayStorage<3>()
    { return FortranArray<3>(); }

    operator GeneralArrayStorage<4>()
    { return FortranArray<4>(); }

    operator GeneralArrayStorage<5>()
    { return FortranArray<5>(); }

    operator GeneralArrayStorage<6>()
    { return FortranArray<6>(); }

    operator GeneralArrayStorage<7>()
    { return FortranArray<7>(); }

    operator GeneralArrayStorage<8>()
    { return FortranArray<8>(); }

    operator GeneralArrayStorage<9>()
    { return FortranArray<9>(); }

    operator GeneralArrayStorage<10>()
    { return FortranArray<10>(); }

    operator GeneralArrayStorage<11>()
    { return FortranArray<11>(); }
};

// A global instance of this class will be placed in
// the blitz library (libblitz.a on unix machines).

_bz_global _bz_fortranTag fortranArray;


/**
   Class ColumnMajorArray specializes GeneralArrayStorage to provide
   column major arrays (column major ordering, base of 0).
 */

template<int N_rank>
class ColumnMajorArray : public GeneralArrayStorage<N_rank> {
private:
    typedef GeneralArrayStorage<N_rank> T_base;
    typedef _bz_typename T_base::noInitializeFlag noInitializeFlag;
    using T_base::ordering_;
    using T_base::ascendingFlag_;
    using T_base::base_;
    using T_base::paddingPolicy_;
public:
  ColumnMajorArray(paddingPolicy pp = defaultPadding)
    : GeneralArrayStorage<N_rank>(noInitializeFlag())
    {
      ordering_ = tensor::i;//Range(0, N_rank - 1);
        ascendingFlag_ = true;
        base_ = 0;
	paddingPolicy_ = pp;
    }
};

// This tag class can be used to provide a nicer notation for
// constructing column major arrays: instead of
//     Array<int,2> A(3, 3, ColumnMajorArray<2>());
// one can simply write:
//     Array<int,2> A(3, 3, columnMajorArray);
// where columnMajorArray is an object of type _bz_columnMajorTag.

class _bz_columnMajorTag {

public:
    operator GeneralArrayStorage<1>()
    { return ColumnMajorArray<1>(); }

    operator GeneralArrayStorage<2>()
    { return ColumnMajorArray<2>(); }

    operator GeneralArrayStorage<3>()
    { return ColumnMajorArray<3>(); }

    operator GeneralArrayStorage<4>()
    { return ColumnMajorArray<4>(); }

    operator GeneralArrayStorage<5>()
    { return ColumnMajorArray<5>(); }

    operator GeneralArrayStorage<6>()
    { return ColumnMajorArray<6>(); }

    operator GeneralArrayStorage<7>()
    { return ColumnMajorArray<7>(); }

    operator GeneralArrayStorage<8>()
    { return ColumnMajorArray<8>(); }

    operator GeneralArrayStorage<9>()
    { return ColumnMajorArray<9>(); }

    operator GeneralArrayStorage<10>()
    { return ColumnMajorArray<10>(); }

    operator GeneralArrayStorage<11>()
    { return ColumnMajorArray<11>(); }
};

// A global instance of this class will be placed in
// the blitz library (libblitz.a on unix machines).

_bz_global _bz_columnMajorTag columnMajorArray;


BZ_NAMESPACE_END

#ifdef BZ_HAVE_BOOST_MPI
namespace boost { namespace mpi {
  template <int N>
  struct is_mpi_datatype<blitz::GeneralArrayStorage<N> > : boost::mpl::true_ { };
} }
#endif

#endif // BZ_ARRAY_STORAGE_H

