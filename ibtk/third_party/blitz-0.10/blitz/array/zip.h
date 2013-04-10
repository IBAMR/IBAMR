// -*- C++ -*-
/***************************************************************************
 * blitz/array/zip.h  "zip" scalar arrays into a multicomponent array expr
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
#ifndef BZ_ARRAYZIP_H
#define BZ_ARRAYZIP_H

#ifndef BZ_ARRAY_H
 #error <blitz/array/zip.h> must be included via <blitz/array.h>
#endif

BZ_NAMESPACE(blitz)

template<typename P_component, typename T1, typename T2>
struct Zip2 {
    typedef P_component T_numtype;

    static inline T_numtype apply(T1 a, T2 b)
    { return T_numtype(a,b); }

    template<typename T_left, typename T_right>
    static inline void prettyPrint(BZ_STD_SCOPE(string) &str,
        prettyPrintFormat& format, const T_left& t1,
        const T_right& t2)
    {
        str += "zip(";
        t1.prettyPrint(str, format);
        str += ",";
        t2.prettyPrint(str, format);
        str += ")";
    }
};

template<typename T_component, typename T1, typename T2>
inline _bz_ArrayExpr<_bz_ArrayExprBinaryOp<
    _bz_typename asExpr<T1>::T_expr,
    _bz_typename asExpr<T2>::T_expr,
    Zip2<T_component,
         _bz_typename asExpr<T1>::T_expr::T_numtype,
         _bz_typename asExpr<T2>::T_expr::T_numtype> > >
zip(const T1& a, const T2& b, T_component)
{
    return _bz_ArrayExpr<_bz_ArrayExprBinaryOp<
        _bz_typename asExpr<T1>::T_expr,
        _bz_typename asExpr<T2>::T_expr,
        Zip2<T_component,
             _bz_typename asExpr<T1>::T_expr::T_numtype,
             _bz_typename asExpr<T2>::T_expr::T_numtype> > >(a,b);
}

BZ_NAMESPACE_END

#endif // BZ_ARRAYZIP_H

