/***************************************************************************
 * blitz/tinyvecio.cc      TinyVector I/O methods
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
#ifndef BZ_TINYVEC2IO_CC
#define BZ_TINYVEC2IO_CC

#include <blitz/tinyvec2.h>

BZ_NAMESPACE(blitz)

// NEEDS_WORK???
// This version of operator<< is updated on August 2005
// by Sergei Mingaleev <mingaleev@gmail.com>. The output
// format for 2D TinyVector is the same as for complex
// numbers - so that one can read Array of complex numbers
// into Array of TinyVectors for speed-up of calculations.
// Also, the corresponding operator>> is updated.

template<typename P_numtype, int N_length>
ostream& operator<<(ostream& os, const TinyVector<P_numtype, N_length>& x)
{
    os << "(" << x[0];
    for (int i=1; i < N_length; ++i)
    {
        os << "," << x[i];
    }
    os << ")";
    return os;
}

// Input of tinyvec contribute by Adam Levar <adaml@mcneilhouse.com>
// and updated by Sergei Mingaleev <mingaleev@gmail.com>
template <typename T_numtype, int N_length>
istream& operator>>(istream& is, TinyVector<T_numtype, N_length>& x)
{
    char sep;

    is >> sep;
    BZPRECHECK(sep == '(', "Format error while scanning input TinyVector"
        << endl << " (expected '(' opening TinyVector)");

    is >> x(0);
    for (int i = 1; i < N_length; ++i)
    {
        is >> sep;
        BZPRECHECK(sep == ',', "Format error while scanning input TinyVector"
             << endl << " (expected ',' between TinyVector components)");
        BZPRECHECK(!is.bad(), "Premature end of input while scanning TinyVector");
        is >> x(i);
    }
    is >> sep;
    BZPRECHECK(sep == ')', "Format error while scanning input TinyVector"
       << endl << " (expected ')' closing TinyVector)");

    return is;
}

BZ_NAMESPACE_END

#endif // BZ_TINYVECIO_CC
