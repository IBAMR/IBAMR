/***************************************************************************
 * blitz/tinymatio.cc      TinyMatrix I/O methods
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
#ifndef BZ_TINYMAT2IO_CC
#define BZ_TINYMAT2IO_CC

#include <blitz/tinymat2.h>

BZ_NAMESPACE(blitz)

// NEEDS_WORK???
// This version of operator<< is updated on August 2005
// by Sergei Mingaleev <mingaleev@gmail.com>.
// Also, the corresponding operator>> is updated.

template <typename P_numtype, int N_rows, int N_columns>
ostream& operator<<(ostream& os,
    const TinyMatrix<P_numtype, N_rows, N_columns>& x)
{
    os << "(";
    for (int i=0; i < N_rows-1; ++i)
    {
        os << x(i,0);
        for (int j=1; j < N_columns; ++j)
        {
            os << "," << x(i,j);
        }
        os << "; ";
    }
    os << x(N_rows-1,0);
    for (int j=1; j < N_columns; ++j)
    {
       os << "," << x(N_rows-1,j);
    }
    os << ")";
    return os;
}

template <typename P_numtype, int N_rows, int N_columns>
istream& operator>>(istream& is,
    TinyMatrix<P_numtype, N_rows, N_columns>& x)
{
    char sep;
    is >> sep;
    BZPRECHECK(sep == '(', "Format error while scanning input TinyMatrix"
      << endl << " (expected '(' opening TinyMatrix)");

    for (int i=0; i<N_rows-1; i++) {
      is >> x(i,0);
      for (int j=1; j<N_columns; j++) {
        is >> sep;
        BZPRECHECK(sep == ',', "Format error while scanning input TinyMatrix"
	  << endl << " (expected ',' between TinyMatrix components)");
        is >> x(i,j);
      }
      is >> sep;
      BZPRECHECK(sep == ';', "Format error while scanning input TinyMatrix"
        << endl << " (expected ';' between TinyMatrix rows)");
    }
    is >> x(N_rows-1,0);
    for (int j=1; j<N_columns; j++) {
      is >> sep;
      BZPRECHECK(sep == ',', "Format error while scanning input TinyMatrix"
        << endl << " (expected ',' between TinyMatrix components)");
      is >> x(N_rows-1,j);
    }
    is >> sep;
    BZPRECHECK(sep == ')', "Format error while scanning input TinyMatrix"
      << endl << " (expected ')' closing TinyMatrix)");

    return is;
}

BZ_NAMESPACE_END

#endif // BZ_TINYMATIO_CC

