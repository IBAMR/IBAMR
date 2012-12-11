/***************************************************************************
 * blitz/array/io.cc  Input/output of arrays.
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
#ifndef BZ_ARRAYIO_CC
#define BZ_ARRAYIO_CC

#ifndef BZ_ARRAY_H
 #error <blitz/array/io.cc> must be included via <blitz/array.h>
#endif

BZ_NAMESPACE(blitz)

// NEEDS_WORK???
// This version of operator<< is updated on August 2005
// by Sergei Mingaleev <mingaleev@gmail.com>.
// Also, the corresponding operator>> is updated.

template<typename T_numtype>
ostream& operator<<(ostream& os, const Array<T_numtype,1>& x)
{
  // Write the extent vector: e.g., (-4, 4)

  os << "(" << x.lbound(0) << "," << x.ubound(0) << ")";
  os << endl << "[ ";

  for (int i1=x.lbound(0); i1<=x.ubound(0); i1++) {
    os << x(i1) << " ";
  }
  os << "]" << endl;
  return os;
}

template<typename T_numtype, int N_rank>
ostream& operator<<(ostream& os, const Array<T_numtype,N_rank>& x)
{
  // Write the extent vector: this is separated by 'x's, e.g.
  // (1, 10) x (-4, 4) x (-5, 5)

  for (int i=0; i < N_rank; ++i) {
    os << "(";
    os << x.lbound(i);
    os << ",";
    os << x.ubound(i);
    os << ")";
    if (i != N_rank-1) os << " x ";
  }
  os << endl << "[ ";

  switch (N_rank) {
    case 2:
      for (int i1=x.lbound(0); i1<=x.ubound(0); i1++) {
        for (int i2=x.lbound(1); i2<=x.ubound(1); i2++) {
          os << x(i1,i2) << " ";
        }
        if (i1 != x.ubound(0)) os << endl << "  ";
      }
      break;
    case 3:
      for (int i1=x.lbound(0); i1<=x.ubound(0); i1++) {
	for (int i2=x.lbound(1); i2<=x.ubound(1); i2++) {
          for (int i3=x.lbound(2); i3<=x.ubound(2); i3++) {
              os << x(i1,i2,i3) << " ";
	  }
	  if (i1 != x.ubound(0) || i2 != x.ubound(1)) os << endl << "  ";
	}
      }
      break;
    case 4:
      for (int i1=x.lbound(0); i1<=x.ubound(0); i1++) {
	for (int i2=x.lbound(1); i2<=x.ubound(1); i2++) {
          for (int i3=x.lbound(2); i3<=x.ubound(2); i3++) {
            for (int i4=x.lbound(3); i4<=x.ubound(3); i4++) {
                os << x(i1,i2,i3,i4) << " ";
	    }
	    if (i1 != x.ubound(0) || i2 != x.ubound(1) || i3 != x.ubound(2))
              os << endl << "  ";
	  }
	}
      }
      break;
    default:
      cout << "Error: operator<< for " << N_rank
           << "D Array is not supported!" << endl;
      BZASSERT("ERROR!");
      break;
  };

  os << "]" << endl;
  return os;
}

/*
 *  Input
 */

template<typename T_numtype, int N_rank>
istream& operator>>(istream& is, Array<T_numtype,N_rank>& x)
{
  TinyVector<int,N_rank> lower_bounds, upper_bounds, extent;
  char sep;

  // Read the extent vector: this is separated by 'x's, e.g.
  // (1, 10) x (-4, 4) x (-5, 5)

  for (int i=0; i < N_rank; ++i) {
    is >> sep;
    BZPRECHECK(!is.bad(), "Premature end of input while scanning Array");
    BZPRECHECK(sep == '(', "Format error while scanning input \
Array \n -- expected '(' opening Array extents");

    is >> lower_bounds(i);
    is >> sep;
    BZPRECHECK(sep == ',', "Format error while scanning input \
Array \n -- expected ',' between Array extents");
    is >> upper_bounds(i);

    is >> sep;
    BZPRECHECK(sep == ')', "Format error while scanning input \
Array \n -- expected ',' closing Array extents");

    if (i != N_rank-1) {
      is >> sep;
      BZPRECHECK(sep == 'x', "Format error while scanning input \
Array \n (expected 'x' between Array extents)");
    }
  }

  is >> sep;
  BZPRECHECK(sep == '[', "Format error while scanning input \
Array \n (expected '[' before beginning of Array data)");

  for (int i=0; i < N_rank; ++i)
      extent(i) = upper_bounds(i) - lower_bounds(i) + 1;
  x.resize(extent);
  x.reindexSelf(lower_bounds);

  switch (N_rank) {
    case 1:
      for (int i1=x.lbound(0); i1<=x.ubound(0); i1++) {
        BZPRECHECK(!is.bad(), "Premature end of input while scanning Array");
	is >> x(i1);
      }
      break;
    case 2:
      for (int i1=x.lbound(0); i1<=x.ubound(0); i1++) {
        for (int i2=x.lbound(1); i2<=x.ubound(1); i2++) {
          BZPRECHECK(!is.bad(), "Premature end of input while scanning Array");
	  is >> x(i1,i2);
	}
      }
      break;
    case 3:
      for (int i1=x.lbound(0); i1<=x.ubound(0); i1++) {
	for (int i2=x.lbound(1); i2<=x.ubound(1); i2++) {
          for (int i3=x.lbound(2); i3<=x.ubound(2); i3++) {
            BZPRECHECK(!is.bad(), "Premature end of input while scanning Array");
	    is >> x(i1,i2,i3);
	  }
	}
      }
      break;
    case 4:
      for (int i1=x.lbound(0); i1<=x.ubound(0); i1++) {
	for (int i2=x.lbound(1); i2<=x.ubound(1); i2++) {
          for (int i3=x.lbound(2); i3<=x.ubound(2); i3++) {
            for (int i4=x.lbound(3); i4<=x.ubound(3); i4++) {
              BZPRECHECK(!is.bad(), "Premature end of input while scanning Array");
  	      is >> x(i1,i2,i3,i4);
	    }
	  }
	}
      }
      break;
    default:
      cout << "Error: read() for " << N_rank
           << "D Array is not supported!" << endl;
      BZASSERT("ERROR!");
      break;
  };

  is >> sep;
  BZPRECHECK(sep == ']', "Format error while scanning input \
Array \n (expected ']' after end of Array data)");
  return is;
}

BZ_NAMESPACE_END

#endif // BZ_ARRAYIO_CC
