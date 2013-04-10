// -*- C++ -*-
/***************************************************************************
 * blitz/vector.h      Declaration of the Vector<P_numtype> class
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
 * KNOWN BUGS
 *
 * 1. operator[](Vector<int>) won't match; compiler complains of no
 *       suitable copy constructor for VectorPick<T>
 * 2. Vector<T>(_bz_VecExpr<E>) constructor generates warning
 * 3. operator+=,-=,..etc.(Random<D>) won't match; compiler complains of
 *       no suitable copy constructor for _bz_VecExpr(...).
 */

#ifndef BZ_VECTOR_H
#define BZ_VECTOR_H

#include <blitz/blitz.h>
#include <blitz/array.h>

BZ_NAMESPACE(blitz)

// Forward declarations
// template<typename P_numtype> class VectorIter;
// template<typename P_numtype> class VectorIterConst;
// template<typename P_expr>    class _bz_VecExpr;
// template<typename P_numtype> class VectorPick;
// template<typename P_numtype> class Random;

// Declaration of class Vector<P_numtype>, which just is a rank-1 array.

template<typename P_numtype>
class Vector : public Array<P_numtype,1>
{
 public:
    typedef Array<P_numtype,1> T_base;
    typedef typename T_base::T_numtype T_numtype;
    typedef typename T_base::T_index T_index;
    typedef typename T_base::T_array T_array;
    typedef typename T_base::T_iterator T_iterator;
    typedef typename T_base::iterator iterator;
    typedef typename T_base::const_iterator const_iterator;

    using T_base::rank_;

    template<typename T_expr>
      explicit Vector(_bz_ArrayExpr<T_expr> expr) :
      T_base(expr) {};

    // forward constructors

    explicit Vector(int length0,
		    GeneralArrayStorage<rank_> storage = GeneralArrayStorage<rank_>()) :
      T_base(length0, storage) {};

  Vector(GeneralArrayStorage<rank_> storage = GeneralArrayStorage<rank_>()) :
      T_base(storage) {};

    Vector(T_numtype* restrict dataFirst, TinyVector<int, rank_> shape,
	   GeneralArrayStorage<rank_> storage = GeneralArrayStorage<rank_>()) :
      T_base(dataFirst, shape, storage) {};

    Vector(Range r0,
	   GeneralArrayStorage<rank_> storage = GeneralArrayStorage<rank_>()) :
      T_base(r0, storage) {};

    Vector(const Array<T_numtype, rank_>& array) : T_base(array) {};

    template<int N_rank2, typename R0, typename R1, typename R2, typename R3, typename R4,
	     typename R5, typename R6, typename R7, typename R8, typename R9, typename R10>
      Vector(Array<T_numtype,N_rank2>& array, R0 r0, R1 r1, R2 r2,
	     R3 r3, R4 r4, R5 r5, R6 r6, R7 r7, R8 r8, R9 r9, R10 r10) :
      T_base(array, r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10) {};

  using T_base::operator=;

};


BZ_NAMESPACE_END

#endif // BZ_VECTOR_H
