// -*- C++ -*-
/***************************************************************************
 * blitz/tvcross.h      Cross product of TinyVector<N,3>'s
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

#ifndef BZ_TVCROSS_H
#define BZ_TVCROSS_H

#ifndef BZ_TINYVEC_H
 #error <blitz/tvcross.h> must be included via <blitz/tinyvec.h>
#endif

BZ_NAMESPACE(blitz)

/*
 * cross product.
 *
 * NEEDS_WORK: - cross product of two different vector types
 *             - cross product involving expressions
 */

template<typename T_numtype>
TinyVector<T_numtype,3> cross(const TinyVector<T_numtype,3>& x,
    const TinyVector<T_numtype,3>& y)
{
    return TinyVector<T_numtype,3>(x[1]*y[2] - y[1]*x[2],
        y[0]*x[2] - x[0]*y[2], x[0]*y[1] - y[0]*x[1]);
}


BZ_NAMESPACE_END

#endif // BZ_TVCROSS_H
