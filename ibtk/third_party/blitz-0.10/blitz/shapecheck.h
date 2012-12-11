// -*- C++ -*-
/***************************************************************************
 * blitz/shapecheck.h    Functions for checking conformability of arrays
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

#ifndef BZ_SHAPECHECK_H
#define BZ_SHAPECHECK_H

#include <iostream>
#include <blitz/blitz.h>

BZ_NAMESPACE(blitz)

/*
 * The function areShapesConformable(A,B) checks that the shapes
 * A and B are conformable (i.e. the same size/geometry).  Typically
 * the A and B parameters are of type TinyVector<int,N_rank> and represent
 * the extent of the arrays.  It's possible that in the future jagged-edged
 * arrays will be supported, in which case shapes may be lists
 * of subdomains.
 */

template<typename T_shape1, typename T_shape2>
inline bool areShapesConformable(const T_shape1&, const T_shape2&)
{
    // If the shape objects are different types, this means
    // that the arrays are different ranks, or one is jagged
    // edged, etc.  In this case the two arrays are not
    // conformable.
    return false;
}

template<typename T_shape>
inline bool areShapesConformable(const T_shape& a, const T_shape& b)
{
    // The shape objects are the same type, so compare them.

    // NEEDS_WORK-- once the "all" reduction is implemented, should
    // use it.
    // return all(a == b);

    for (int i=0; i < a.length(); ++i)
    {
        if (a[i] != b[i])
        {
	  BZ_DEBUG_MESSAGE("Incompatible shapes detected: " << std::endl
			   << a << std::endl << b << std::endl);
            return false;
        }
    }

    return true;
}

BZ_NAMESPACE_END

#endif
