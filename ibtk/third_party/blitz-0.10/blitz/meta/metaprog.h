// -*- C++ -*-
/***************************************************************************
 * blitz/meta/metaprog.h   Useful metaprogram declarations
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

#ifndef BZ_META_METAPROG_H
#define BZ_META_METAPROG_H

BZ_NAMESPACE(blitz)

// Null Operand

class _bz_meta_nullOperand {
public:
    _bz_meta_nullOperand() { }
};

template<typename T> inline T operator+(const T& a, _bz_meta_nullOperand)
{ return a; }
template<typename T> inline T operator*(const T& a, _bz_meta_nullOperand)
{ return a; }

// MetaMax

template<int N1, int N2>
class _bz_meta_max {
public:
    static const int max = (N1 > N2) ? N1 : N2;
};

// MetaMin

template<int N1, int N2>
class _bz_meta_min {
public:
    static const int min = (N1 < N2) ? N1 : N2;
};

BZ_NAMESPACE_END

#endif // BZ_META_METAPROG_H
