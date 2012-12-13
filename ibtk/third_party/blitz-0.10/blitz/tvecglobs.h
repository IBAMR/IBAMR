// -*- C++ -*-
/***************************************************************************
 * blitz/tvecglobs.h     TinyVector global functions
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

#ifndef BZ_TVECGLOBS_H
#define BZ_TVECGLOBS_H

#include <blitz/meta/metaprog.h>
#include <blitz/numtrait.h>

#include <blitz/tvcross.h>       // Cross products
#include <blitz/meta/dot.h>
#include <blitz/meta/product.h>
#include <blitz/meta/sum.h>

BZ_NAMESPACE(blitz)

template<typename T_numtype1, typename T_numtype2, int N_length>
inline BZ_PROMOTE(T_numtype1, T_numtype2)
dot(const TinyVector<T_numtype1, N_length>& a,
    const TinyVector<T_numtype2, N_length>& b)
{
    return _bz_meta_vectorDot<N_length, 0>::f(a,b);
}

template<typename T_numtype1, int N_length>
inline BZ_SUMTYPE(T_numtype1)
product(const TinyVector<T_numtype1, N_length>& a)
{
    return _bz_meta_vectorProduct<N_length, 0>::f(a);
}

template<typename T_numtype, int N_length>
inline BZ_SUMTYPE(T_numtype)
sum(const TinyVector<T_numtype, N_length>& a)
{
    return _bz_meta_vectorSum<N_length, 0>::f(a);
}

template<typename T_numtype, int N_length>
inline BZ_SUMTYPE(T_numtype)
norm2(const TinyVector<T_numtype, N_length>& a)
{
    return std::sqrt(sum(pow2(a)));
}

// explicit returntype functions follow

template<typename T_ret>
class _bz_returntype {
public:

  template<typename T_numtype1, typename T_numtype2, int N_length>
static inline T_ret
dot(const TinyVector<T_numtype1, N_length>& a,
    const TinyVector<T_numtype2, N_length>& b)
{
  return _bz_meta_vectorDotRet<N_length, 0, T_ret>::f(a,b);
}

template<typename T_numtype1, int N_length>
static inline T_ret
product(const TinyVector<T_numtype1, N_length>& a)
{
    return _bz_meta_vectorProductRet<N_length, 0, T_ret>::f(a);
}

template<typename T_numtype, int N_length>
static inline T_ret
sum(const TinyVector<T_numtype, N_length>& a)
{
    return _bz_meta_vectorSumRet<N_length, 0, T_ret>::f(a);
}

};

BZ_NAMESPACE_END

#endif // BZ_TVECGLOBS_H
