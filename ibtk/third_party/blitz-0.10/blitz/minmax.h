// -*- C++ -*-
/***************************************************************************
 * blitz/minmax.h  Declaration of min and max functions
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

#ifndef BZ_MINMAX_H
#define BZ_MINMAX_H

#include <blitz/promote.h>

BZ_NAMESPACE(blitz)

/*
 * These functions are in their own namespace (blitz::minmax) to avoid
 * conflicts with the array reduction operations min and max.
 */

BZ_NAMESPACE(extrema)

template<typename T1, typename T2>
BZ_PROMOTE(T1,T2) (min)(const T1& a, const T2& b)
{
    typedef BZ_PROMOTE(T1,T2) T_promote;

    if (a <= b)
        return T_promote(a);
    else
        return T_promote(b);
}

template<typename T1, typename T2>
BZ_PROMOTE(T1,T2) (max)(const T1& a, const T2& b)
{
    typedef BZ_PROMOTE(T1,T2) T_promote;

    if (a >= b)
        return T_promote(a);
    else
        return T_promote(b);
}

BZ_NAMESPACE_END

BZ_NAMESPACE_END

#endif
