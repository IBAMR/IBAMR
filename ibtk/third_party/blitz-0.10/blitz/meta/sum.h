// -*- C++ -*-
/***************************************************************************
 * blitz/meta/sum.h      TinyVector sum metaprogram
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

#ifndef BZ_META_SUM_H
#define BZ_META_SUM_H

#ifndef BZ_METAPROG_H
 #include <blitz/meta/metaprog.h>
#endif

BZ_NAMESPACE(blitz)

template<int N, int I>
class _bz_meta_vectorSum {
public:
    static const int loopFlag = (I < N-1) ? 1 : 0;

    template<typename T_expr1>
    static inline BZ_SUMTYPE(_bz_typename T_expr1::T_numtype)
    f(const T_expr1& a)
    {
        return a[I] +
            _bz_meta_vectorSum<loopFlag * N, loopFlag * (I+1)>::f(a);
    }
};

template<>
class _bz_meta_vectorSum<0,0> {
public:
    template<typename T_expr1>
    static inline _bz_meta_nullOperand f(const T_expr1&)
    { return _bz_meta_nullOperand(); }

};

template<int N, int I, typename T_ret>
class _bz_meta_vectorSumRet {
public:
    static const int loopFlag = (I < N-1) ? 1 : 0;

    template<typename T_expr1>
    static inline T_ret
    f(const T_expr1& a)
    {
      return static_cast<T_ret>(a[I]) +
	_bz_meta_vectorSumRet<loopFlag * N, loopFlag * (I+1), T_ret>::f(a);
    }
};

template<typename T_ret>
class _bz_meta_vectorSumRet<0,0, T_ret> {
public:
    template<typename T_expr1>
    static inline _bz_meta_nullOperand f(const T_expr1&)
    { return _bz_meta_nullOperand(); }

};

BZ_NAMESPACE_END

#endif // BZ_META_SUM_H
