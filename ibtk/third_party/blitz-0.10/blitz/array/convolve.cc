/***************************************************************************
 * blitz/array/convolve.cc  One-dimensional convolution
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

#ifndef BZ_ARRAY_CONVOLVE_CC
#define BZ_ARRAY_CONVOLVE_CC

BZ_NAMESPACE(blitz)

template<typename T>
Array<T,1> convolve(const Array<T,1>& B, const Array<T,1>& C)
{
    int Bl = B.lbound(0), Bh = B.ubound(0);
    int Cl = C.lbound(0), Ch = C.ubound(0);

    int lbound = Bl + Cl;
    int ubound = Bh + Ch;

    Array<T,1> A(Range(lbound,ubound));

    T result;
    for (int i=lbound; i <= ubound; ++i)
    {
        int jl = i - Ch;
        if (jl < Bl)
            jl = Bl;

        int jh = i - Cl;
        if (jh > Bh)
            jh = Bh;

        result = 0;
        for (int j=jl; j <= jh; ++j)
            result += B(j) * C(i-j);

        A(i) = result;
    }

    return A;
}

BZ_NAMESPACE_END

#endif // BZ_ARRAY_CONVOLVE_CC

