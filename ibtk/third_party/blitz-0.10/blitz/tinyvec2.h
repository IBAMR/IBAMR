// -*- C++ -*-
/***************************************************************************
 * blitz/tinyvec.h      Declaration of the TinyVector<T, N> class
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

#ifndef BZ_TINYVEC_H
#define BZ_TINYVEC_H

#ifndef __INTEL_COMPILER
#ifdef __GNUC__
#if (__GNUC_MINOR__ >= 6)
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#endif
#endif

#include <blitz/blitz.h>
#include <blitz/listinit.h>
#include <blitz/et-forward.h>
#include <blitz/etbase.h>
#include <blitz/simdtypes.h>
#include <blitz/array/slice.h>

#ifdef BZ_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/serialization.hpp>
#endif
#ifdef BZ_HAVE_BOOST_MPI
#include <boost/mpi/datatype.hpp>
#endif

#ifdef BZ_HAVE_CSTRING
#include <cstring> // For memcpy
#endif

BZ_NAMESPACE(blitz)

/*****************************************************************************
 * Forward declarations
 */

template<typename P_numtype, int N_length>
class FastTV2Iterator;
template<typename P_numtype, int N_length>
class FastTV2CopyIterator;


/** The TinyVector class is a one-dimensional, fixed length vector
    that implements the blitz expression template
    machinery. TinyVector-only expressions are very fast because they
    usually get reduced to just the unrolled (and vectorized, if
    enabled) assembly instructions. TinyVectors can also be used in
    mixed expressions with other ET classes. */
template<typename P_numtype, int N_length>
class TinyVector : public ETBase<TinyVector<P_numtype, N_length> >
{
public:

    //////////////////////////////////////////////
    // Public Types
    //////////////////////////////////////////////

    typedef P_numtype                                    T_numtype;
    typedef TinyVector<T_numtype,N_length>               T_vector;
    typedef FastTV2Iterator<T_numtype,N_length>         T_iterator;
    typedef T_numtype*                                   iterator;
    typedef const T_numtype*                             const_iterator;
  typedef FastTV2CopyIterator<P_numtype, N_length> T_range_result;

    static const int
    //numArrayOperands = 1,
    //numIndexPlaceholders = 0,
        rank_ = 1;

    TinyVector()  { }
    ~TinyVector() { }

  TinyVector(const TinyVector<T_numtype,N_length>& x);

    template <typename T_numtype2>
    TinyVector(const TinyVector<T_numtype2,N_length>& x);

  /** This constructor creates a TinyVector from another ETBase
      object. It needs to be explicit to avoid all kinds of
      ambiguities. */
    template <typename T_expr>
    inline explicit TinyVector(const ETBase<T_expr>& expr) {
      *this = expr; }

  /** This constructor creates a TinyVector specifically from an
      expression. This one we do NOT want to be explicit because that
      breaks simple construction assignments like "TinyVector<double,
      1> v = a+b;", forcing the user to explicitly write it like a
      construction. */
    template <typename T_expr>
    inline TinyVector(const _bz_ArrayExpr<T_expr>& expr) {
      *this = expr; }

    inline TinyVector(const T_numtype initValue);

    inline TinyVector(const T_numtype x[]) {
        memcpy(data_,x,N_length*sizeof(T_numtype));
    }


    TinyVector(T_numtype x0, T_numtype x1)
    {
        data_[0] = x0;
        data_[1] = x1;
    }

    TinyVector(T_numtype x0, T_numtype x1, T_numtype x2)
    {
        data_[0] = x0;
        data_[1] = x1;
        data_[2] = x2;
    }

    TinyVector(T_numtype x0, T_numtype x1, T_numtype x2,
        T_numtype x3)
    {
        data_[0] = x0;
        data_[1] = x1;
        data_[2] = x2;
        data_[3] = x3;
    }

    TinyVector(T_numtype x0, T_numtype x1, T_numtype x2,
        T_numtype x3, T_numtype x4)
    {
        data_[0] = x0;
        data_[1] = x1;
        data_[2] = x2;
        data_[3] = x3;
        data_[4] = x4;
    }

    TinyVector(T_numtype x0, T_numtype x1, T_numtype x2,
        T_numtype x3, T_numtype x4, T_numtype x5)
    {
        data_[0] = x0;
        data_[1] = x1;
        data_[2] = x2;
        data_[3] = x3;
        data_[4] = x4;
        data_[5] = x5;
    }

    TinyVector(T_numtype x0, T_numtype x1, T_numtype x2,
        T_numtype x3, T_numtype x4, T_numtype x5, T_numtype x6)
    {
        data_[0] = x0;
        data_[1] = x1;
        data_[2] = x2;
        data_[3] = x3;
        data_[4] = x4;
        data_[5] = x5;
        data_[6] = x6;
    }

    TinyVector(T_numtype x0, T_numtype x1, T_numtype x2,
        T_numtype x3, T_numtype x4, T_numtype x5, T_numtype x6,
        T_numtype x7)
    {
        data_[0] = x0;
        data_[1] = x1;
        data_[2] = x2;
        data_[3] = x3;
        data_[4] = x4;
        data_[5] = x5;
        data_[6] = x6;
        data_[7] = x7;
    }

    TinyVector(T_numtype x0, T_numtype x1, T_numtype x2,
        T_numtype x3, T_numtype x4, T_numtype x5, T_numtype x6,
        T_numtype x7, T_numtype x8)
    {
        data_[0] = x0;
        data_[1] = x1;
        data_[2] = x2;
        data_[3] = x3;
        data_[4] = x4;
        data_[5] = x5;
        data_[6] = x6;
        data_[7] = x7;
        data_[8] = x8;
    }

    TinyVector(T_numtype x0, T_numtype x1, T_numtype x2,
        T_numtype x3, T_numtype x4, T_numtype x5, T_numtype x6,
        T_numtype x7, T_numtype x8, T_numtype x9)
    {
        data_[0] = x0;
        data_[1] = x1;
        data_[2] = x2;
        data_[3] = x3;
        data_[4] = x4;
        data_[5] = x5;
        data_[6] = x6;
        data_[7] = x7;
        data_[8] = x8;
        data_[9] = x9;
    }

    TinyVector(T_numtype x0, T_numtype x1, T_numtype x2,
        T_numtype x3, T_numtype x4, T_numtype x5, T_numtype x6,
        T_numtype x7, T_numtype x8, T_numtype x9, T_numtype x10)
    {
        data_[0] = x0;
        data_[1] = x1;
        data_[2] = x2;
        data_[3] = x3;
        data_[4] = x4;
        data_[5] = x5;
        data_[6] = x6;
        data_[7] = x7;
        data_[8] = x8;
        data_[9] = x9;
        data_[10] = x10;
    }

  static int base()
  { return 0; }

  static int                               base(int rank)
  { BZPRECONDITION(rank==0); return 0; }


    T_iterator      beginFast() const       { return T_iterator(*this);      }

    iterator       begin()       { return data_; }
    const_iterator begin() const { return data_; }

  static int                               dimensions()
    { return 1; }

    iterator       end()       { return data_ + N_length; }
    const_iterator end() const { return data_ + N_length; }

    T_numtype * restrict data()
    { return data_; }

    const T_numtype * restrict data() const
    { return data_; }

    T_numtype * restrict dataFirst()
    { return data_; }

    const T_numtype * restrict dataFirst() const
    { return data_; }

    const TinyVector<int, rank_>    shape() const
    { return N_length; }

  static int                               lbound(int rank)
  { BZPRECONDITION(rank==0); return 0; }
  static int            lbound()
  { return 0; }

  static int                               length(int rank)
  { BZPRECONDITION(rank==0); return N_length; }
  static int    length()
  { return N_length; }

  static int                               extent(int rank)
  { BZPRECONDITION(rank==0); return N_length; }

  static int                               ordering(int storageRankIndex)
  { return 0; }

  static int    ordering()
  { return 0; }

  static  int                               rank()
    { return rank_; }

    static sizeType                               numElements()
  { return length(); }

    static diffType    stride()
    { return 1; }

  static diffType                               stride(int rank)
    { BZPRECONDITION(rank==0); return 1; }

  static int                               ubound(int rank)
  { BZPRECONDITION(rank==0); return length()-1; }

  static int           ubound()
  { return length()-1; }

     template<typename P_expr, typename P_updater>
     void _bz_assign(P_expr, P_updater);

    T_numtype operator*() const
    { return *data_; }

    //////////////////////////////////////////////
    // Subscripting operators
    //////////////////////////////////////////////

    T_vector& noConst() const
    { return const_cast<T_vector&>(*this); }

  static bool lengthCheck(unsigned i)
    {
        BZPRECHECK(i < N_length,
            "TinyVector<" << BZ_DEBUG_TEMPLATE_AS_STRING_LITERAL(T_numtype)
            << "," << N_length << "> index out of bounds: " << i);
        return true;
    }

    const T_numtype& operator()(unsigned i) const
    {
        BZPRECONDITION(lengthCheck(i));
        return data_[i];
    }

    T_numtype& restrict operator()(unsigned i)
    {
        BZPRECONDITION(lengthCheck(i));
        return data_[i];
    }

  T_numtype operator()(TinyVector<int,1> i) const
    {
        BZPRECONDITION(lengthCheck(i[0]));
        return data_[i[0]];
    }

    template<int N0>
    _bz_ArrayExpr<ArrayIndexMapping<typename asExpr<T_vector>::T_expr, N0> >
    operator()(IndexPlaceholder<N0>) const;

    const T_numtype& operator[](unsigned i) const
    {
        BZPRECONDITION(lengthCheck(i));
        return data_[i];
    }

    T_numtype& restrict operator[](unsigned i)
    {
        BZPRECONDITION(lengthCheck(i));
        return data_[i];
    }

  // must return reference so the iterator can turn it into an
  // iterator for the contained in case we have a multicomponent.
  const T_numtype& fastRead(diffType i) const
    { return data_[i]; }

  /** Since data_ is simd aligned by construction, we just have
      to check the offest. */
  bool isVectorAligned(diffType offset) const
  { return (offset%simdTypes<T_numtype>::vecWidth)==0; }

  bool canCollapse(int outerLoopRank, int innerLoopRank) const
  {
    BZPRECONDITION(outerLoopRank==0);
    BZPRECONDITION(innerLoopRank==0);
    return true;
  }

    //////////////////////////////////////////////
    // Assignment operators
    //////////////////////////////////////////////

    // Scalar operand
    ListInitializationSwitch<T_vector,T_numtype*> operator=(T_numtype x)
    {
        return ListInitializationSwitch<T_vector,T_numtype*>(*this, x);
    }

  T_vector& initialize(T_numtype);

    template<typename T_expr>
    T_vector& operator=(const ETBase<T_expr>&);

    template<typename T> T_vector& operator+=(const T&);
    template<typename T> T_vector& operator-=(const T&);
    template<typename T> T_vector& operator*=(const T&);
    template<typename T> T_vector& operator/=(const T&);
    template<typename T> T_vector& operator%=(const T&);
    template<typename T> T_vector& operator^=(const T&);
    template<typename T> T_vector& operator&=(const T&);
    template<typename T> T_vector& operator|=(const T&);
    template<typename T> T_vector& operator>>=(const T&);
    template<typename T> T_vector& operator<<=(const T&);

    T_numtype* restrict getInitializationIterator()
    { return dataFirst(); }

  // // vectors can't be sliced
  // template<typename T1, typename T2 = nilArraySection,
  // 	   class T3 = nilArraySection, typename T4 = nilArraySection,
  // 	   class T5 = nilArraySection, typename T6 = nilArraySection,
  // 	   class T7 = nilArraySection, typename T8 = nilArraySection,
  // 	   class T9 = nilArraySection, typename T10 = nilArraySection,
  // 	   class T11 = nilArraySection>
  // class SliceInfo {
  // public:
  //   typedef void T_slice;
  // };

private:
  template<typename T_expr, typename T_update>
  void _tv_evaluate(const T_expr& expr, T_update);

#ifdef BZ_HAVE_BOOST_SERIALIZATION
  friend class boost::serialization::access;

  template<class T_arch>
  void serialize(T_arch& ar, const unsigned int version) {
    ar & data_;
  };
#endif


  BZ_ALIGN_VARIABLE(T_numtype, data_[N_length], BZ_SIMD_WIDTH)
};

// Specialization for N = 0: KCC is giving some
// peculiar errors, perhaps this will fix.

template<typename T>
class TinyVector<T,0> {
};

BZ_NAMESPACE_END

#ifdef BZ_HAVE_BOOST_SERIALIZATION
namespace boost {
  namespace mpi {
    template<typename T> struct is_mpi_datatype;
    template <typename T, int N>
    struct is_mpi_datatype<blitz::TinyVector<T, N> >
      : public is_mpi_datatype<T> { };
  } };
#endif

#include <blitz/tinyvec2.cc>

#ifndef __INTEL_COMPILER
#if (__GNUC__ && (__GNUC_MINOR__ >= 6))
#pragma GCC diagnostic pop
#endif
#endif

#endif // BZ_TINYVEC_H
