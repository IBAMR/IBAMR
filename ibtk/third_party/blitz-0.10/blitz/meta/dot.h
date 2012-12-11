// -*- C++ -*-
/***************************************************************************
 * blitz/meta/dot.h      Tiny vector dot product metaprogram
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

#ifndef BZ_META_DOT_H
#define BZ_META_DOT_H

#ifndef BZ_PROMOTE_H
 #include <blitz/promote.h>
#endif

#ifndef BZ_METAPROG_H
 #include <blitz/meta/metaprog.h>
#endif

BZ_NAMESPACE(blitz)

template<int N, int I>
class _bz_meta_vectorDot {
public:
    static const int loopFlag = (I < N-1) ? 1 : 0;

    template<typename T_expr1, typename T_expr2>
    static inline BZ_PROMOTE(_bz_typename T_expr1::T_numtype,
        _bz_typename T_expr2::T_numtype)
    f(const T_expr1& a, const T_expr2& b)
    {
        return a[I] * b[I]
            + _bz_meta_vectorDot<loopFlag * N, loopFlag * (I+1)>::f(a,b);
    }

    template<typename T_expr1, typename T_expr2>
    static inline BZ_PROMOTE(_bz_typename T_expr1::T_numtype,
        _bz_typename T_expr2::T_numtype)
    f_value_ref(T_expr1 a, const T_expr2& b)
    {
        return a[I] * b[I]
            + _bz_meta_vectorDot<loopFlag * N, loopFlag * (I+1)>::f(a,b);
    }

    template<typename T_expr1, typename T_expr2>
    static inline BZ_PROMOTE(_bz_typename T_expr1::T_numtype,
        _bz_typename T_expr2::T_numtype)
    f_ref_value(const T_expr1& a, T_expr2 b)
    {
        return a[I] * b[I]
            + _bz_meta_vectorDot<loopFlag * N, loopFlag * (I+1)>::f(a,b);
    }

    template<typename T_expr1, typename P_numtype2>
    static inline BZ_PROMOTE(_bz_typename T_expr1::T_numtype,
        P_numtype2)
    dotWithArgs(const T_expr1& a, P_numtype2 i1, P_numtype2 i2=0,
        P_numtype2 i3=0, P_numtype2 i4=0, P_numtype2 i5=0, P_numtype2 i6=0,
        P_numtype2 i7=0, P_numtype2 i8=0, P_numtype2 i9=0, P_numtype2 i10=0)
    {
        return a[I] * i1
            + _bz_meta_vectorDot<loopFlag * N, loopFlag * (I+1)>::dotWithArgs
                 (a, i2, i3, i4, i5, i6, i7, i8, i9);
    }
};

template<>
class _bz_meta_vectorDot<0,0> {
public:
    template<typename T_expr1, typename T_expr2>
    static inline _bz_meta_nullOperand f(const T_expr1&, const T_expr2&)
    { return _bz_meta_nullOperand(); }

    template<typename T_expr1, typename P_numtype2>
    static inline _bz_meta_nullOperand
    dotWithArgs(const T_expr1& a, P_numtype2 i1, P_numtype2 i2=0,
        P_numtype2 i3=0, P_numtype2 i4=0, P_numtype2 i5=0, P_numtype2 i6=0,
        P_numtype2 i7=0, P_numtype2 i8=0, P_numtype2 i9=0, P_numtype2 i10=0)
    {
        return _bz_meta_nullOperand();
    }

};


template<int N, int I, typename T_ret>
class _bz_meta_vectorDotRet {
public:
    static const int loopFlag = (I < N-1) ? 1 : 0;

    template<typename T_expr1, typename T_expr2>
    static inline T_ret
    f(const T_expr1& a, const T_expr2& b)
    {
      return static_cast<T_ret>(a[I]) * static_cast<T_ret>(b[I])
	+ _bz_meta_vectorDotRet<loopFlag * N, loopFlag * (I+1), T_ret>::f(a,b);
    }

    template<typename T_expr1, typename T_expr2>
    static inline T_ret
    f_value_ref(T_expr1 a, const T_expr2& b)
    {
      return static_cast<T_ret>(a[I]) * static_cast<T_ret>(b[I])
	+ _bz_meta_vectorDotRet<loopFlag * N, loopFlag * (I+1), T_ret>::f(a,b);
    }

    template<typename T_expr1, typename T_expr2>
    static inline T_ret
    f_ref_value(const T_expr1& a, T_expr2 b)
    {
      return static_cast<T_ret>(a[I]) * static_cast<T_ret>(b[I])
	+ _bz_meta_vectorDotRet<loopFlag * N, loopFlag * (I+1), T_ret>::f(a,b);
    }

    template<typename T_expr1, typename P_numtype2>
    static inline T_ret
    dotWithArgs(const T_expr1& a, P_numtype2 i1, P_numtype2 i2=0,
        P_numtype2 i3=0, P_numtype2 i4=0, P_numtype2 i5=0, P_numtype2 i6=0,
        P_numtype2 i7=0, P_numtype2 i8=0, P_numtype2 i9=0, P_numtype2 i10=0)
    {
      return static_cast<T_ret>(a[I]) * static_cast<T_ret>(i1)
	+ _bz_meta_vectorDotRet<loopFlag * N, loopFlag * (I+1), T_ret>::dotWithArgs
                 (a, i2, i3, i4, i5, i6, i7, i8, i9);
    }
};

template<typename T_ret>
class _bz_meta_vectorDotRet<0,0, T_ret> {
public:
    template<typename T_expr1, typename T_expr2>
    static inline _bz_meta_nullOperand f(const T_expr1&, const T_expr2&)
    { return _bz_meta_nullOperand(); }

    template<typename T_expr1, typename P_numtype2>
    static inline _bz_meta_nullOperand
    dotWithArgs(const T_expr1& a, P_numtype2 i1, P_numtype2 i2=0,
        P_numtype2 i3=0, P_numtype2 i4=0, P_numtype2 i5=0, P_numtype2 i6=0,
        P_numtype2 i7=0, P_numtype2 i8=0, P_numtype2 i9=0, P_numtype2 i10=0)
    {
        return _bz_meta_nullOperand();
    }

};

BZ_NAMESPACE_END

#endif // BZ_META_DOT_H
