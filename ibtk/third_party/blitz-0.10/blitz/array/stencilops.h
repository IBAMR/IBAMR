// -*- C++ -*-
/***************************************************************************
 * blitz/array/stencilops.h  Stencil operators
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
#ifndef BZ_ARRAYSTENCILOPS_H
#define BZ_ARRAYSTENCILOPS_H

// NEEDS_WORK: need to factor many of the stencils in terms of the
// integer constants, e.g. 16*(A(-1,0)+A(0,-1)+A(0,1)+A(1,0))

#ifndef BZ_ARRAYSTENCILS_H
 #error <blitz/array/stencilops.h> must be included via <blitz/array/stencils.h>
#endif

#include <blitz/array/geometry.h>
#include <blitz/tinymat2.h>

BZ_NAMESPACE(blitz)

/* These operators aren't really meant to be applied by the user and
   they being named the same as the stencil operators applied to
   arrays cause difficulties with the resolution. For this reason, the
   operators are defined with _stencilop appended to their name. */
#define BZ_DECLARE_STENCIL_OPERATOR1(name,A)				\
  template<typename T>							\
  inline _bz_typename T::T_numtype name ## _stencilop(const T& A)		\
  {

#define BZ_END_STENCIL_OPERATOR   }

#define BZ_DECLARE_STENCIL_OPERATOR2(name,A,B)				\
  template<typename T1, typename T2>					\
  inline BZ_PROMOTE(_bz_typename T1::T_numtype,				\
		    _bz_typename T2::T_numtype) name ## _stencilop(const T1& A, \
								   const T2& B) \
  {

#define BZ_DECLARE_STENCIL_OPERATOR3(name,A,B,C)			\
  template<typename T1, typename T2, typename T3>			\
  inline BZ_PROMOTE(BZ_PROMOTE(_bz_typename T1::T_numtype,		\
			       _bz_typename T2::T_numtype),		\
		    _bz_typename T3::T_numtype) name ## _stencilop(const T1& A,	\
								   const T2& B,	\
								   const T3& C) \
  {

// These constants are accurate to 45 decimal places = 149 bits of mantissa
const double recip_2 = .500000000000000000000000000000000000000000000;
const double recip_4 = .250000000000000000000000000000000000000000000;
const double recip_6 = .166666666666666666666666666666666666666666667;
const double recip_8 = .125000000000000000000000000000000000000000000;
const double recip_12 = .0833333333333333333333333333333333333333333333;
const double recip_144 = .00694444444444444444444444444444444444444444444;

/****************************************************************************
 * Laplacian Operators
 ****************************************************************************/

BZ_DECLARE_STENCIL_OPERATOR1(Laplacian2D, A)
  return -4.0 * (*A)
    + A.shift(-1,0) + A.shift(1,0)
    + A.shift(-1,1) + A.shift(1,1);
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR1(Laplacian3D, A)
  return -6.0 * (*A)
    + A.shift(-1,0) + A.shift(1,0)
    + A.shift(-1,1) + A.shift(1,1)
    + A.shift(-1,2) + A.shift(1,2);
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR1(Laplacian2D4, A)
  return -60.0 * (*A)
    + 16.0 * (A.shift(-1,0) + A.shift(1,0) + A.shift(-1,1) + A.shift(1,1))
    -        (A.shift(-2,0) + A.shift(2,0) + A.shift(-2,1) + A.shift(2,1));
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR1(Laplacian2D4n, A)
  return Laplacian2D4(A) * recip_12;
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR1(Laplacian3D4, A)
  return -90.0 * (*A)
    + 16.0 * (A.shift(-1,0) + A.shift(1,0) + A.shift(-1,1) + A.shift(1,1) +
              A.shift(-1,2) + A.shift(1,2))
    -        (A.shift(-2,0) + A.shift(2,0) + A.shift(-2,1) + A.shift(2,1) +
              A.shift(-2,2) + A.shift(2,2));
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR1(Laplacian3D4n, A)
  return Laplacian3D4(A) * recip_12;
BZ_END_STENCIL_OPERATOR

/****************************************************************************
 * Derivatives
 ****************************************************************************/

#define BZ_DECLARE_DIFF(name)                                                 \
  template<typename T>                                                        \
  inline _bz_typename T::T_numtype name ## _stencilop(const T& A, int dim = firstDim)

#define BZ_DECLARE_MULTIDIFF(name)                                            \
  template<typename T>                                                        \
  inline _bz_typename multicomponent_traits<_bz_typename                      \
     T::T_numtype>::T_element name ## _stencilop(const T& A, int comp, int dim)

/****************************************************************************
 * Central differences with accuracy O(h^2)
 ****************************************************************************/

BZ_DECLARE_DIFF(central12) {
  return A.shift(1,dim) - A.shift(-1,dim);
}

BZ_DECLARE_DIFF(central22) {
  return -2.0 * (*A) + A.shift(1,dim) + A.shift(-1,dim);
}

BZ_DECLARE_DIFF(central32) {
  return -2.0 * (A.shift(1,dim) - A.shift(-1,dim))
    +           (A.shift(2,dim) - A.shift(-2,dim));
}

BZ_DECLARE_DIFF(central42) {
  return 6.0 * (*A)
    - 4.0 * (A.shift(1,dim) + A.shift(-1,dim))
    +       (A.shift(2,dim) + A.shift(-2,dim));
}

/****************************************************************************
 * Central differences with accuracy O(h^2)  (multicomponent versions)
 ****************************************************************************/

BZ_DECLARE_MULTIDIFF(central12) {
 return A.shift(1,dim)[comp] - A.shift(-1,dim)[comp];
}

BZ_DECLARE_MULTIDIFF(central22) {
  return -2.0 * (*A)[comp]
    + A.shift(1,dim)[comp] + A.shift(-1,dim)[comp];
}

BZ_DECLARE_MULTIDIFF(central32) {
  return -2.0 * (A.shift(1,dim)[comp] - A.shift(-1,dim)[comp])
    +           (A.shift(2,dim)[comp] - A.shift(-2,dim)[comp]);
}

BZ_DECLARE_MULTIDIFF(central42) {
  return 6.0 * (*A)[comp]
    -4.0 * (A.shift(1,dim)[comp] + A.shift(-1,dim)[comp])
    +      (A.shift(2,dim)[comp] + A.shift(-2,dim)[comp]);
}

/****************************************************************************
 * Central differences with accuracy O(h^2)  (normalized versions)
 ****************************************************************************/

BZ_DECLARE_DIFF(central12n) {
  return central12_stencilop(A,dim) * recip_2;
}

BZ_DECLARE_DIFF(central22n) {
  return central22_stencilop(A,dim);
}

BZ_DECLARE_DIFF(central32n) {
  return central32_stencilop(A,dim) * recip_2;
}

BZ_DECLARE_DIFF(central42n) {
  return central42_stencilop(A,dim);
}

/****************************************************************************
 * Central differences with accuracy O(h^2)  (normalized multicomponent)
 ****************************************************************************/

BZ_DECLARE_MULTIDIFF(central12n) {
  return central12_stencilop(A,comp,dim) * recip_2;
}

BZ_DECLARE_MULTIDIFF(central22n) {
  return central22_stencilop(A,comp,dim);
}

BZ_DECLARE_MULTIDIFF(central32n) {
  return central32_stencilop(A,comp,dim) * recip_2;
}

BZ_DECLARE_MULTIDIFF(central42n) {
  return central42_stencilop(A,comp,dim);
}

/****************************************************************************
 * Central differences with accuracy O(h^4)
 ****************************************************************************/

BZ_DECLARE_DIFF(central14) {
  return 8.0 * (A.shift(1,dim) - A.shift(-1,dim))
    -          (A.shift(2,dim) - A.shift(-2,dim));
}

BZ_DECLARE_DIFF(central24) {
  return -30.0 * (*A)
    + 16.0 * (A.shift(1,dim) + A.shift(-1,dim))
    -        (A.shift(2,dim) + A.shift(-2,dim));
}

BZ_DECLARE_DIFF(central34) {
  return -13.0 * (A.shift(1,dim) - A.shift(-1,dim))
    +      8.0 * (A.shift(2,dim) - A.shift(-2,dim))
    -            (A.shift(3,dim) - A.shift(-3,dim));
}

BZ_DECLARE_DIFF(central44) {
  return 56.0 * (*A)
    - 39.0 * (A.shift(1,dim) + A.shift(-1,dim))
    + 12.0 * (A.shift(2,dim) + A.shift(-2,dim))
    -        (A.shift(3,dim) + A.shift(-3,dim));
}

/****************************************************************************
 * Central differences with accuracy O(h^4)  (multicomponent versions)
 ****************************************************************************/

BZ_DECLARE_MULTIDIFF(central14) {
  return 8.0 * (A.shift(1,dim)[comp] - A.shift(-1,dim)[comp])
    -          (A.shift(2,dim)[comp] - A.shift(-2,dim)[comp]);
}

BZ_DECLARE_MULTIDIFF(central24) {
  return - 30.0 * (*A)[comp]
    + 16.0 * (A.shift(1,dim)[comp] + A.shift(-1,dim)[comp])
    -        (A.shift(2,dim)[comp] + A.shift(-2,dim)[comp]);
}

BZ_DECLARE_MULTIDIFF(central34) {
  return -13.0 * (A.shift(1,dim)[comp] - A.shift(-1,dim)[comp])
    +      8.0 * (A.shift(2,dim)[comp] - A.shift(-2,dim)[comp])
    -            (A.shift(3,dim)[comp] - A.shift(-3,dim)[comp]);
}

BZ_DECLARE_MULTIDIFF(central44) {
  return 56.0 * (*A)[comp]
    - 39.0 * (A.shift(1,dim)[comp] + A.shift(-1,dim)[comp])
    + 12.0 * (A.shift(2,dim)[comp] + A.shift(-2,dim)[comp])
    -        (A.shift(3,dim)[comp] + A.shift(-3,dim)[comp]);
}

/****************************************************************************
 * Central differences with accuracy O(h^4)  (normalized)
 ****************************************************************************/

BZ_DECLARE_DIFF(central14n) {
  return central14_stencilop(A,dim) * recip_12;
}

BZ_DECLARE_DIFF(central24n) {
  return central24_stencilop(A,dim) * recip_12;
}

BZ_DECLARE_DIFF(central34n) {
  return central34_stencilop(A,dim) * recip_8;
}

BZ_DECLARE_DIFF(central44n) {
  return central44_stencilop(A,dim) * recip_6;
}

/****************************************************************************
 * Central differences with accuracy O(h^4)  (normalized, multicomponent)
 ****************************************************************************/

BZ_DECLARE_MULTIDIFF(central14n) {
  return central14_stencilop(A,comp,dim) * recip_12;
}

BZ_DECLARE_MULTIDIFF(central24n) {
  return central24_stencilop(A,comp,dim) * recip_12;
}

BZ_DECLARE_MULTIDIFF(central34n) {
  return central34_stencilop(A,comp,dim) * recip_8;
}

BZ_DECLARE_MULTIDIFF(central44n) {
  return central44_stencilop(A,comp,dim) * recip_6;
}

/****************************************************************************
 * Backward differences with accuracy O(h)
 ****************************************************************************/

BZ_DECLARE_DIFF(backward11) {
  return (*A) - A.shift(-1,dim);
}

BZ_DECLARE_DIFF(backward21) {
  return (*A) - 2.0 * A.shift(-1,dim) + A.shift(-2,dim);
}

BZ_DECLARE_DIFF(backward31) {
  return (*A) - 3.0 * A.shift(-1,dim) + 3.0 * A.shift(-2,dim)
    - A.shift(-3,dim);
}

BZ_DECLARE_DIFF(backward41) {
  return (*A) - 4.0 * A.shift(-1,dim) + 6.0 * A.shift(-2,dim)
    - 4.0 * A.shift(-3,dim) + A.shift(-4,dim);
}

/****************************************************************************
 * Backward differences with accuracy O(h) (multicomponent versions)
 ****************************************************************************/

BZ_DECLARE_MULTIDIFF(backward11) {
  return (*A)[comp] - A.shift(-1,dim)[comp];
}

BZ_DECLARE_MULTIDIFF(backward21) {
  return (*A)[comp] - 2.0 * A.shift(-1,dim)[comp] + A.shift(-2,dim)[comp];
}

BZ_DECLARE_MULTIDIFF(backward31) {
  return (*A)[comp] - 3.0 * A.shift(-1,dim)[comp] + 3.0 * A.shift(-2,dim)[comp]
    - A.shift(-3,dim)[comp];
}

BZ_DECLARE_MULTIDIFF(backward41) {
  return (*A)[comp] - 4.0 * A.shift(-1,dim)[comp] + 6.0 * A.shift(-2,dim)[comp]
    - 4.0 * A.shift(-3,dim)[comp] + A.shift(-4,dim)[comp];
}

/****************************************************************************
 * Backward differences with accuracy O(h)  (normalized)
 ****************************************************************************/

BZ_DECLARE_DIFF(backward11n) { return backward11_stencilop(A,dim); }
BZ_DECLARE_DIFF(backward21n) { return backward21_stencilop(A,dim); }
BZ_DECLARE_DIFF(backward31n) { return backward31_stencilop(A,dim); }
BZ_DECLARE_DIFF(backward41n) { return backward41_stencilop(A,dim); }

/****************************************************************************
 * Backward differences with accuracy O(h)  (normalized, multicomponent)
 ****************************************************************************/

BZ_DECLARE_MULTIDIFF(backward11n) { return backward11_stencilop(A,comp,dim); }
BZ_DECLARE_MULTIDIFF(backward21n) { return backward21_stencilop(A,comp,dim); }
BZ_DECLARE_MULTIDIFF(backward31n) { return backward31_stencilop(A,comp,dim); }
BZ_DECLARE_MULTIDIFF(backward41n) { return backward41_stencilop(A,comp,dim); }

/****************************************************************************
 * Backward differences with accuracy O(h^2)
 ****************************************************************************/

BZ_DECLARE_DIFF(backward12) {
  return 3.0 * (*A) - 4.0 * A.shift(-1,dim) + A.shift(-2,dim);
}

BZ_DECLARE_DIFF(backward22) {
  return 2.0 * (*A) - 5.0 * A.shift(-1,dim) + 4.0 * A.shift(-2,dim)
    - A.shift(-3,dim);
}

BZ_DECLARE_DIFF(backward32) {
  return 5.0 * (*A) - 18.0 * A.shift(-1,dim) + 24.0 * A.shift(-2,dim)
    - 14.0 * A.shift(-3,dim) + 3.0 * A.shift(-4,dim);
}

BZ_DECLARE_DIFF(backward42) {
  return 3.0 * (*A) - 14.0 * A.shift(-1,dim) + 26.0 * A.shift(-2,dim)
    - 24.0 * A.shift(-3,dim) + 11.0 * A.shift(-4,dim) - 2.0 * A.shift(-5,dim);
}

/****************************************************************************
 * Backward differences with accuracy O(h^2) (multicomponent versions)
 ****************************************************************************/

BZ_DECLARE_MULTIDIFF(backward12) {
  return 3.0 * (*A)[comp] - 4.0 * A.shift(-1,dim)[comp]
    + A.shift(-2,dim)[comp];
}

BZ_DECLARE_MULTIDIFF(backward22) {
  return 2.0 * (*A)[comp] - 5.0 * A.shift(-1,dim)[comp]
    + 4.0 * A.shift(-2,dim)[comp] - A.shift(-3,dim)[comp];
}

BZ_DECLARE_MULTIDIFF(backward32) {
  return 5.0 * (*A)[comp] - 18.0 * A.shift(-1,dim)[comp]
    + 24.0 * A.shift(-2,dim)[comp] - 14.0 * A.shift(-3,dim)[comp]
    + 3.0 * A.shift(-4,dim)[comp];
}

BZ_DECLARE_MULTIDIFF(backward42) {
  return 3.0 * (*A)[comp] - 14.0 * A.shift(-1,dim)[comp]
    + 26.0 * A.shift(-2,dim)[comp] - 24.0 * A.shift(-3,dim)[comp]
    + 11.0 * A.shift(-4,dim)[comp] - 2.0 * A.shift(-5,dim)[comp];
}

/****************************************************************************
 * Backward differences with accuracy O(h^2)  (normalized)
 ****************************************************************************/

BZ_DECLARE_DIFF(backward12n) { return backward12_stencilop(A,dim) * recip_2; }
BZ_DECLARE_DIFF(backward22n) { return backward22_stencilop(A,dim); }
BZ_DECLARE_DIFF(backward32n) { return backward32_stencilop(A,dim) * recip_2; }
BZ_DECLARE_DIFF(backward42n) { return backward42_stencilop(A,dim); }

/****************************************************************************
 * Backward differences with accuracy O(h^2)  (normalized, multicomponent)
 ****************************************************************************/

BZ_DECLARE_MULTIDIFF(backward12n) { return backward12_stencilop(A,comp,dim) * recip_2; }
BZ_DECLARE_MULTIDIFF(backward22n) { return backward22_stencilop(A,comp,dim); }
BZ_DECLARE_MULTIDIFF(backward32n) { return backward32_stencilop(A,comp,dim) * recip_2; }
BZ_DECLARE_MULTIDIFF(backward42n) { return backward42_stencilop(A,comp,dim); }

/****************************************************************************
 * Forward differences with accuracy O(h)
 ****************************************************************************/

BZ_DECLARE_DIFF(forward11) {
  return -(*A) + A.shift(1,dim);
}

BZ_DECLARE_DIFF(forward21) {
  return (*A) - 2.0 * A.shift(1,dim) + A.shift(2,dim);
}

BZ_DECLARE_DIFF(forward31) {
  return -(*A) + 3.0 * A.shift(1,dim) - 3.0 * A.shift(2,dim) + A.shift(3,dim);
}

BZ_DECLARE_DIFF(forward41) {
  return (*A) - 4.0 * A.shift(1,dim) + 6.0 * A.shift(2,dim)
    - 4.0 * A.shift(3,dim) + A.shift(4,dim);
}

/****************************************************************************
 * Forward differences with accuracy O(h)  (multicomponent versions)
 ****************************************************************************/

BZ_DECLARE_MULTIDIFF(forward11) {
  return  -(*A)[comp] + A.shift(1,dim)[comp];
}

BZ_DECLARE_MULTIDIFF(forward21) {
  return (*A)[comp] - 2.0 * A.shift(1,dim)[comp] + A.shift(2,dim)[comp];
}

BZ_DECLARE_MULTIDIFF(forward31) {
  return -(*A)[comp] + 3.0 * A.shift(1,dim)[comp] - 3.0 * A.shift(2,dim)[comp]
    + A.shift(3,dim)[comp];
}

BZ_DECLARE_MULTIDIFF(forward41) {
  return (*A)[comp] - 4.0 * A.shift(1,dim)[comp] + 6.0 * A.shift(2,dim)[comp]
    - 4.0 * A.shift(3,dim)[comp] + A.shift(4,dim)[comp];
}

/****************************************************************************
 * Forward differences with accuracy O(h)     (normalized)
 ****************************************************************************/

BZ_DECLARE_DIFF(forward11n) { return forward11_stencilop(A,dim); }
BZ_DECLARE_DIFF(forward21n) { return forward21_stencilop(A,dim); }
BZ_DECLARE_DIFF(forward31n) { return forward31_stencilop(A,dim); }
BZ_DECLARE_DIFF(forward41n) { return forward41_stencilop(A,dim); }

/****************************************************************************
 * Forward differences with accuracy O(h)     (multicomponent,normalized)
 ****************************************************************************/

BZ_DECLARE_MULTIDIFF(forward11n) { return forward11_stencilop(A,comp,dim); }
BZ_DECLARE_MULTIDIFF(forward21n) { return forward21_stencilop(A,comp,dim); }
BZ_DECLARE_MULTIDIFF(forward31n) { return forward31_stencilop(A,comp,dim); }
BZ_DECLARE_MULTIDIFF(forward41n) { return forward41_stencilop(A,comp,dim); }

/****************************************************************************
 * Forward differences with accuracy O(h^2)
 ****************************************************************************/

BZ_DECLARE_DIFF(forward12) {
  return -3.0 * (*A) + 4.0 * A.shift(1,dim) - A.shift(2,dim);
}

BZ_DECLARE_DIFF(forward22) {
  return 2.0 * (*A) - 5.0 * A.shift(1,dim) + 4.0 * A.shift(2,dim)
    - A.shift(3,dim);
}

BZ_DECLARE_DIFF(forward32) {
  return -5.0 * (*A) + 18.0 * A.shift(1,dim) - 24.0 * A.shift(2,dim)
    + 14.0 * A.shift(3,dim) - 3.0 * A.shift(4,dim);
}

BZ_DECLARE_DIFF(forward42) {
  return 3.0 * (*A) - 14.0 * A.shift(1,dim) + 26.0 * A.shift(2,dim)
    - 24.0 * A.shift(3,dim) + 11.0 * A.shift(4,dim) - 2.0 * A.shift(5,dim);
}

/****************************************************************************
 * Forward differences with accuracy O(h^2)   (multicomponent versions)
 ****************************************************************************/

BZ_DECLARE_MULTIDIFF(forward12) {
  return -3.0 * (*A)[comp] + 4.0 * A.shift(1,dim)[comp] - A.shift(2,dim)[comp];
}

BZ_DECLARE_MULTIDIFF(forward22) {
  return 2.0 * (*A)[comp] - 5.0 * A.shift(1,dim)[comp]
    + 4.0 * A.shift(2,dim)[comp] - A.shift(3,dim)[comp];
}

BZ_DECLARE_MULTIDIFF(forward32) {
  return -5.0 * (*A)[comp] + 18.0 * A.shift(1,dim)[comp]
    - 24.0 * A.shift(2,dim)[comp] + 14.0 * A.shift(3,dim)[comp]
    - 3.0 * A.shift(4,dim)[comp];
}

BZ_DECLARE_MULTIDIFF(forward42) {
  return 3.0 * (*A)[comp] - 14.0 * A.shift(1,dim)[comp]
    + 26.0 * A.shift(2,dim)[comp] - 24.0 * A.shift(3,dim)[comp]
    + 11.0 * A.shift(4,dim)[comp] - 2.0 * A.shift(5,dim)[comp];
}


/****************************************************************************
 * Forward differences with accuracy O(h^2)     (normalized)
 ****************************************************************************/

BZ_DECLARE_DIFF(forward12n) { return forward12_stencilop(A,dim) * recip_2; }
BZ_DECLARE_DIFF(forward22n) { return forward22_stencilop(A,dim); }
BZ_DECLARE_DIFF(forward32n) { return forward32_stencilop(A,dim) * recip_2; }
BZ_DECLARE_DIFF(forward42n) { return forward42_stencilop(A,dim); }

/****************************************************************************
 * Forward differences with accuracy O(h^2)     (normalized)
 ****************************************************************************/

BZ_DECLARE_MULTIDIFF(forward12n) { return forward12_stencilop(A,comp,dim) * recip_2; }
BZ_DECLARE_MULTIDIFF(forward22n) { return forward22_stencilop(A,comp,dim); }
BZ_DECLARE_MULTIDIFF(forward32n) { return forward32_stencilop(A,comp,dim) * recip_2; }
BZ_DECLARE_MULTIDIFF(forward42n) { return forward42_stencilop(A,comp,dim); }

/****************************************************************************
 * Gradient operators
 ****************************************************************************/

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,2> grad2D_stencilop(const T& A) {
  return TinyVector<_bz_typename T::T_numtype,2>(
    central12_stencilop(A,firstDim),
    central12_stencilop(A,secondDim));
}

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,2> grad2D4_stencilop(const T& A) {
  return TinyVector<_bz_typename T::T_numtype,2>(
    central14_stencilop(A,firstDim),
    central14_stencilop(A,secondDim));
}

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,3> grad3D_stencilop(const T& A) {
  return TinyVector<_bz_typename T::T_numtype,3>(
    central12_stencilop(A,firstDim),
    central12_stencilop(A,secondDim),
    central12_stencilop(A,thirdDim));
}

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,3> grad3D4_stencilop(const T& A) {
  return TinyVector<_bz_typename T::T_numtype,3>(
    central14_stencilop(A,firstDim),
    central14_stencilop(A,secondDim),
    central14_stencilop(A,thirdDim));
}

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,2> grad2Dn_stencilop(const T& A) {
  return TinyVector<_bz_typename T::T_numtype,2>(
    central12n_stencilop(A,firstDim),
    central12n_stencilop(A,secondDim));
}

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,2> grad2D4n_stencilop(const T& A) {
  return TinyVector<_bz_typename T::T_numtype,2>(
    central14n_stencilop(A,firstDim),
    central14n_stencilop(A,secondDim));
}

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,3> grad3Dn_stencilop(const T& A) {
  return TinyVector<_bz_typename T::T_numtype,3>(
    central12n_stencilop(A,firstDim),
    central12n_stencilop(A,secondDim),
    central12n_stencilop(A,thirdDim));
}

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,3> grad3D4n_stencilop(const T& A) {
  return TinyVector<_bz_typename T::T_numtype,3>(
    central14n_stencilop(A,firstDim),
    central14n_stencilop(A,secondDim),
    central14n_stencilop(A,thirdDim));
}

/****************************************************************************
 * Grad-squared operators
 ****************************************************************************/

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,2> gradSqr2D_stencilop(const T& A) {
  return TinyVector<_bz_typename T::T_numtype,2>(
    central22_stencilop(A,firstDim),
    central22_stencilop(A,secondDim));
}

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,2> gradSqr2D4_stencilop(const T& A) {
  return TinyVector<_bz_typename T::T_numtype,2>(
    central24_stencilop(A,firstDim),
    central24_stencilop(A,secondDim));
}

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,3> gradSqr3D_stencilop(const T& A) {
  return TinyVector<_bz_typename T::T_numtype,3>(
    central22_stencilop(A,firstDim),
    central22_stencilop(A,secondDim),
    central22_stencilop(A,thirdDim));
}

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,3> gradSqr3D4_stencilop(const T& A) {
  return TinyVector<_bz_typename T::T_numtype,3>(
    central24_stencilop(A,firstDim),
    central24_stencilop(A,secondDim),
    central24_stencilop(A,thirdDim));
}

/****************************************************************************
 * Grad-squared operators (normalized)
 ****************************************************************************/

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,2> gradSqr2Dn_stencilop(const T& A) {
  return gradSqr2D_stencilop(A);
}

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,2> gradSqr2D4n_stencilop(const T& A) {
  return TinyVector<_bz_typename T::T_numtype,2>(
    central24_stencilop(A,firstDim) * recip_12,
    central24_stencilop(A,secondDim) * recip_12);
}

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,3> gradSqr3Dn_stencilop(const T& A) {
  return gradSqr3D_stencilop(A);(A);
}

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,3> gradSqr3D4n_stencilop(const T& A) {
  return TinyVector<_bz_typename T::T_numtype,3>(
    central24_stencilop(A,firstDim) * recip_12,
    central24_stencilop(A,secondDim) * recip_12,
    central24_stencilop(A,thirdDim) * recip_12);
}

/****************************************************************************
 * Gradient operators on a vector field
 ****************************************************************************/

template<typename T>
inline TinyMatrix<_bz_typename multicomponent_traits<_bz_typename
    T::T_numtype>::T_element, 3, 3>
Jacobian3D_stencilop(const T& A)
{
    const int x=0, y=1, z=2;
    const int u=0, v=1, w=2;

    TinyMatrix<_bz_typename multicomponent_traits<_bz_typename
        T::T_numtype>::T_element, 3, 3> grad;

    grad(u,x) = central12_stencilop(A,u,x);
    grad(u,y) = central12_stencilop(A,u,y);
    grad(u,z) = central12_stencilop(A,u,z);
    grad(v,x) = central12_stencilop(A,v,x);
    grad(v,y) = central12_stencilop(A,v,y);
    grad(v,z) = central12_stencilop(A,v,z);
    grad(w,x) = central12_stencilop(A,w,x);
    grad(w,y) = central12_stencilop(A,w,y);
    grad(w,z) = central12_stencilop(A,w,z);

    return grad;
}

template<typename T>
inline TinyMatrix<_bz_typename multicomponent_traits<_bz_typename
    T::T_numtype>::T_element, 3, 3>
Jacobian3Dn_stencilop(const T& A)
{
    const int x=0, y=1, z=2;
    const int u=0, v=1, w=2;

    TinyMatrix<_bz_typename multicomponent_traits<_bz_typename
        T::T_numtype>::T_element, 3, 3> grad;

    grad(u,x) = central12n_stencilop(A,u,x);
    grad(u,y) = central12n_stencilop(A,u,y);
    grad(u,z) = central12n_stencilop(A,u,z);
    grad(v,x) = central12n_stencilop(A,v,x);
    grad(v,y) = central12n_stencilop(A,v,y);
    grad(v,z) = central12n_stencilop(A,v,z);
    grad(w,x) = central12n_stencilop(A,w,x);
    grad(w,y) = central12n_stencilop(A,w,y);
    grad(w,z) = central12n_stencilop(A,w,z);

    return grad;
}

template<typename T>
inline TinyMatrix<_bz_typename multicomponent_traits<_bz_typename
    T::T_numtype>::T_element, 3, 3>
Jacobian3D4_stencilop(const T& A)
{
    const int x=0, y=1, z=2;
    const int u=0, v=1, w=2;

    TinyMatrix<_bz_typename multicomponent_traits<_bz_typename
        T::T_numtype>::T_element, 3, 3> grad;

    grad(u,x) = central14_stencilop(A,u,x);
    grad(u,y) = central14_stencilop(A,u,y);
    grad(u,z) = central14_stencilop(A,u,z);
    grad(v,x) = central14_stencilop(A,v,x);
    grad(v,y) = central14_stencilop(A,v,y);
    grad(v,z) = central14_stencilop(A,v,z);
    grad(w,x) = central14_stencilop(A,w,x);
    grad(w,y) = central14_stencilop(A,w,y);
    grad(w,z) = central14_stencilop(A,w,z);

    return grad;
}

template<typename T>
inline TinyMatrix<_bz_typename multicomponent_traits<_bz_typename
    T::T_numtype>::T_element, 3, 3>
Jacobian3D4n_stencilop(const T& A)
{
    const int x=0, y=1, z=2;
    const int u=0, v=1, w=2;

    TinyMatrix<_bz_typename multicomponent_traits<_bz_typename
        T::T_numtype>::T_element, 3, 3> grad;

    grad(u,x) = central14n_stencilop(A,u,x);
    grad(u,y) = central14n_stencilop(A,u,y);
    grad(u,z) = central14n_stencilop(A,u,z);
    grad(v,x) = central14n_stencilop(A,v,x);
    grad(v,y) = central14n_stencilop(A,v,y);
    grad(v,z) = central14n_stencilop(A,v,z);
    grad(w,x) = central14n_stencilop(A,w,x);
    grad(w,y) = central14n_stencilop(A,w,y);
    grad(w,z) = central14n_stencilop(A,w,z);

    return grad;
}

/****************************************************************************
 * Curl operators
 ****************************************************************************/

// O(h^2) curl, using central difference

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,3>
curl_stencilop(const T& vx, const T& vy, const T& vz) {
  const int x = firstDim, y = secondDim, z = thirdDim;

  return TinyVector<_bz_typename T::T_numtype,3>(
    central12_stencilop(vz,y)-central12_stencilop(vy,z),
    central12_stencilop(vx,z)-central12_stencilop(vz,x),
    central12_stencilop(vy,x)-central12_stencilop(vx,y));
}

// Normalized O(h^2) curl, using central difference
template<typename T>
inline TinyVector<_bz_typename T::T_numtype,3>
curln_stencilop(const T& vx, const T& vy, const T& vz) {
  const int x = firstDim, y = secondDim, z = thirdDim;

  return TinyVector<_bz_typename T::T_numtype,3>(
    (central12_stencilop(vz,y)-central12_stencilop(vy,z)) * recip_2,
    (central12_stencilop(vx,z)-central12_stencilop(vz,x)) * recip_2,
    (central12_stencilop(vy,x)-central12_stencilop(vx,y)) * recip_2);
}

// Multicomponent curl
template<typename T>
inline _bz_typename T::T_numtype curl3D_stencilop(const T& A) {
  const int x = firstDim, y = secondDim, z = thirdDim;

  return _bz_typename T::T_numtype(
    central12_stencilop(A,z,y)-central12_stencilop(A,y,z),
    central12_stencilop(A,x,z)-central12_stencilop(A,z,x),
    central12_stencilop(A,y,x)-central12_stencilop(A,x,y));
}

// Normalized multicomponent curl
template<typename T>
inline _bz_typename T::T_numtype curl3Dn_stencilop(const T& A) {
  const int x = firstDim, y = secondDim, z = thirdDim;

  return _bz_typename T::T_numtype(
    (central12_stencilop(A,z,y)-central12_stencilop(A,y,z)) * recip_2,
    (central12_stencilop(A,x,z)-central12_stencilop(A,z,x)) * recip_2,
    (central12_stencilop(A,y,x)-central12_stencilop(A,x,y)) * recip_2);
}

// O(h^4) curl, using 4th order central difference
template<typename T>
inline TinyVector<_bz_typename T::T_numtype,3>
curl4_stencilop(const T& vx, const T& vy, const T& vz) {
  const int x = firstDim, y = secondDim, z = thirdDim;

  return TinyVector<_bz_typename T::T_numtype,3>(
    central14_stencilop(vz,y)-central14_stencilop(vy,z),
    central14_stencilop(vx,z)-central14_stencilop(vz,x),
    central14_stencilop(vy,x)-central14_stencilop(vx,y));
}

// Normalized O(h^4) curl, using 4th order central difference
template<typename T>
inline TinyVector<_bz_typename T::T_numtype,3>
curl4n_stencilop(const T& vx, const T& vy, const T& vz) {
  const int x = firstDim, y = secondDim, z = thirdDim;

  return TinyVector<_bz_typename T::T_numtype,3>(
    (central14_stencilop(vz,y)-central14_stencilop(vy,z)) * recip_2,
    (central14_stencilop(vx,z)-central14_stencilop(vz,x)) * recip_2,
    (central14_stencilop(vy,x)-central14_stencilop(vx,y)) * recip_2);
}

// O(h^4) curl, using 4th order central difference (multicomponent version)
template<typename T>
inline _bz_typename T::T_numtype
curl3D4_stencilop(const T& A) {
  const int x = firstDim, y = secondDim, z = thirdDim;

  return _bz_typename T::T_numtype(
    central14_stencilop(A,z,y)-central14_stencilop(A,y,z),
    central14_stencilop(A,x,z)-central14_stencilop(A,z,x),
    central14_stencilop(A,y,x)-central14_stencilop(A,x,y));
}

// O(h^4) curl, using 4th order central difference (normalized multicomponent)
template<typename T>
inline _bz_typename T::T_numtype
curl3D4n_stencilop(const T& A) {
  const int x = firstDim, y = secondDim, z = thirdDim;

  return _bz_typename T::T_numtype(
    (central14_stencilop(A,z,y)-central14_stencilop(A,y,z)) * recip_2,
    (central14_stencilop(A,x,z)-central14_stencilop(A,z,x)) * recip_2,
    (central14_stencilop(A,y,x)-central14_stencilop(A,x,y)) * recip_2);
}



// Two-dimensional curl

template<typename T>
inline _bz_typename T::T_numtype
curl_stencilop(const T& vx, const T& vy) {
  const int x = firstDim, y = secondDim;

  return central12_stencilop(vy,x)-central12_stencilop(vx,y);
}

template<typename T>
inline _bz_typename T::T_numtype
curln_stencilop(const T& vx, const T& vy) {
  const int x = firstDim, y = secondDim;

  return (central12_stencilop(vy,x)-central12_stencilop(vx,y)) * recip_2;
}

// Multicomponent curl
template<typename T>
inline _bz_typename T::T_numtype::T_numtype curl2D_stencilop(const T& A) {
  const int x = firstDim, y = secondDim;
  return central12_stencilop(A,y,x)-central12_stencilop(A,x,y);
}

template<typename T>
inline _bz_typename T::T_numtype::T_numtype curl2Dn_stencilop(const T& A) {
  const int x = firstDim, y = secondDim;
  return (central12_stencilop(A,y,x)-central12_stencilop(A,x,y)) * recip_2;
}


// 4th order versions

template<typename T>
inline _bz_typename T::T_numtype
curl4_stencilop(const T& vx, const T& vy) {
  const int x = firstDim, y = secondDim;

  return central14_stencilop(vy,x)-central14_stencilop(vx,y);
}

template<typename T>
inline _bz_typename T::T_numtype
curl4n_stencilop(const T& vx, const T& vy) {
  const int x = firstDim, y = secondDim;

  return (central14_stencilop(vy,x)-central14_stencilop(vx,y)) * recip_12;
}

// Multicomponent curl
template<typename T>
inline _bz_typename T::T_numtype::T_numtype curl2D4_stencilop(const T& A) {
  const int x = firstDim, y = secondDim;
  return central14_stencilop(A,y,x)-central14_stencilop(A,x,y);
}

template<typename T>
inline _bz_typename T::T_numtype::T_numtype curl2D4n_stencilop(const T& A) {
  const int x = firstDim, y = secondDim;
  return (central14_stencilop(A,y,x)-central14_stencilop(A,x,y)) * recip_12;
}

/****************************************************************************
 * Divergence
 ****************************************************************************/


BZ_DECLARE_STENCIL_OPERATOR2(div,vx,vy)
  return central12_stencilop(vx,firstDim) + central12_stencilop(vy,secondDim);
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR2(divn,vx,vy)
  return (central12_stencilop(vx,firstDim) + central12_stencilop(vy,secondDim))
     * recip_2;
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR2(div4,vx,vy)
  return central14_stencilop(vx,firstDim) + central14_stencilop(vy,secondDim);
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR2(div4n,vx,vy)
  return (central14_stencilop(vx,firstDim) + central14_stencilop(vy,secondDim))
    * recip_12;
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR3(div,vx,vy,vz)
  return central12_stencilop(vx,firstDim) + central12_stencilop(vy,secondDim)
    + central12_stencilop(vz,thirdDim);
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR3(divn,vx,vy,vz)
  return (central12_stencilop(vx,firstDim) + central12_stencilop(vy,secondDim)
    + central12_stencilop(vz,thirdDim)) * recip_2;
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR3(div4,vx,vy,vz)
  return central14_stencilop(vx,firstDim) + central14_stencilop(vy,secondDim)
    + central14_stencilop(vz,thirdDim);
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR3(div4n,vx,vy,vz)
  return (central14_stencilop(vx,firstDim) + central14_stencilop(vy,secondDim)
    + central14_stencilop(vz,thirdDim)) * recip_12;
BZ_END_STENCIL_OPERATOR

// these return a scalar, which is T_numtype of T::T_result (which may be an ET)
template<typename T>
inline _bz_typename T::T_result::T_numtype
div2D_stencilop(const T& A)
{
    const int x = firstDim, y = secondDim;
    return central12_stencilop(A,x,x) + central12_stencilop(A,y,y);
}

template<typename T>
inline _bz_typename T::T_result::T_numtype
div2Dn_stencilop(const T& A)
{
    const int x = firstDim, y = secondDim;
    return (central12_stencilop(A,x,x) + central12_stencilop(A,y,y)) * recip_2;
}

template<typename T>
inline _bz_typename T::T_result::T_numtype
div2D4_stencilop(const T& A)
{
    const int x = firstDim, y = secondDim;
    return central14_stencilop(A,x,x) + central14_stencilop(A,y,y);
}

template<typename T>
inline _bz_typename T::T_result::T_numtype
div2D4n_stencilop(const T& A)
{
    const int x = firstDim, y = secondDim;
    return (central14_stencilop(A,x,x) + central14_stencilop(A,y,y)) * recip_12;
}

template<typename T>
inline _bz_typename T::T_result::T_numtype
div3D_stencilop(const T& A)
{
    const int x = firstDim, y = secondDim, z = thirdDim;
    return central12_stencilop(A,x,x) + central12_stencilop(A,y,y) + central12_stencilop(A,z,z);
}

template<typename T>
inline _bz_typename T::T_result::T_numtype
div3Dn_stencilop(const T& A)
{
    const int x = firstDim, y = secondDim, z = thirdDim;
    return (central12_stencilop(A,x,x) + central12_stencilop(A,y,y) + central12_stencilop(A,z,z)) * recip_2;
}

template<typename T>
inline _bz_typename T::T_result::T_numtype
div3D4_stencilop(const T& A)
{
    const int x = firstDim, y = secondDim, z = thirdDim;
    return central14_stencilop(A,x,x) + central14_stencilop(A,y,y) + central14_stencilop(A,z,z);
}

template<typename T>
inline _bz_typename T::T_result::T_numtype
div3D4n_stencilop(const T& A)
{
    const int x = firstDim, y = secondDim, z = thirdDim;
    return (central14_stencilop(A,x,x) + central14_stencilop(A,y,y) + central14_stencilop(A,z,z)) * recip_12;
}

/****************************************************************************
 * Mixed Partial derivatives
 ****************************************************************************/

template<typename T>
inline _bz_typename T::T_numtype
mixed22_stencilop(const T& A, int x, int y)
{
    return A.shift(-1,x,-1,y) - A.shift(-1,x,1,y)
        -  A.shift(1,x,-1,y) + A.shift(1,x,1,y);
}

template<typename T>
inline _bz_typename T::T_numtype
mixed22n_stencilop(const T& A, int x, int y)
{
    return mixed22_stencilop(A,x,y) * recip_4;
}

template<typename T>
inline _bz_typename T::T_numtype
mixed24_stencilop(const T& A, int x, int y)
{
    return 64.0 * (A.shift(-1,x,-1,y) - A.shift(-1,x,1,y) -
                   A.shift(1,x,-1,y) + A.shift(1,x,1,y))
        +         (A.shift(-2,x,1,y) - A.shift(-1,x,2,y) -
                   A.shift(1,x,2,y) - A.shift(2,x,1,y) +
                   A.shift(2,x,-1,y) + A.shift(1,x,-2,y) -
                   A.shift(-1,x,-2,y) + A.shift(-2,x,-1,y))
        +   8.0 * (A.shift(-1,x,1,y) + A.shift(-1,x,2,y) -
                   A.shift(2,x,-2,y) + A.shift(2,x,2,y));
}

template<typename T>
inline _bz_typename T::T_numtype
mixed24n_stencilop(const T& A, int x, int y)
{
    return mixed24_stencilop(A,x,y) * recip_144;
}

/****************************************************************************
 * Smoothers
 ****************************************************************************/

// NEEDS_WORK-- put other stencil operators here:
//   Average5pt2D
//   Average7pt3D
// etc.

/****************************************************************************
 * Stencil operators with geometry (experimental)
 ****************************************************************************/

template<typename T>
inline _bz_typename multicomponent_traits<_bz_typename
    T::T_numtype>::T_element div3DVec4_stencilop(const T& A,
    const UniformCubicGeometry<3>& geom)
{
    const int x = 0, y = 1, z = 2;

    return (central14_stencilop(A, x, firstDim) + central14_stencilop(A, y, secondDim)
        + central14_stencilop(A, z, thirdDim)) * recip_12 * geom.recipSpatialStep();
}

template<typename T>
inline _bz_typename T::T_numtype Laplacian3D4_stencilop(const T& A,
    const UniformCubicGeometry<3>& geom)
{
    return Laplacian3D4n_stencilop(A) * geom.recipSpatialStepPow2();
}

template<typename T>
inline _bz_typename T::T_numtype Laplacian3DVec4_stencilop(const T& A,
    const UniformCubicGeometry<3>& geom)
{
    typedef _bz_typename T::T_numtype vector3d;
    typedef _bz_typename multicomponent_traits<vector3d>::T_element
        T_element;
    const int u = 0, v = 1, w = 2;
    const int x = 0, y = 1, z = 2;

    // central24 is a 5-point stencil
    // This is a 9*5 = 45 point stencil

    T_element t1 = (central24_stencilop(A,u,x) + central24_stencilop(A,u,y) + central24_stencilop(A,u,z))
        * recip_12 * geom.recipSpatialStepPow2();

    T_element t2 = (central24_stencilop(A,v,x) + central24_stencilop(A,v,y) + central24_stencilop(A,v,z))
        * recip_12 * geom.recipSpatialStepPow2();

    T_element t3 = (central24_stencilop(A,w,x) + central24_stencilop(A,w,y) + central24_stencilop(A,w,z))
        * recip_12 * geom.recipSpatialStepPow2();

    return vector3d(t1,t2,t3);
}

template<typename T>
inline TinyMatrix<_bz_typename multicomponent_traits<_bz_typename
    T::T_numtype>::T_element, 3, 3>
grad3DVec4_stencilop(const T& A, const UniformCubicGeometry<3>& geom)
{
    const int x=0, y=1, z=2;
    const int u=0, v=1, w=2;

    TinyMatrix<_bz_typename multicomponent_traits<_bz_typename
        T::T_numtype>::T_element, 3, 3> grad;

    // This is a 9*4 = 36 point stencil
    grad(u,x) = central14n_stencilop(A,u,x) * geom.recipSpatialStep();
    grad(u,y) = central14n_stencilop(A,u,y) * geom.recipSpatialStep();
    grad(u,z) = central14n_stencilop(A,u,z) * geom.recipSpatialStep();
    grad(v,x) = central14n_stencilop(A,v,x) * geom.recipSpatialStep();
    grad(v,y) = central14n_stencilop(A,v,y) * geom.recipSpatialStep();
    grad(v,z) = central14n_stencilop(A,v,z) * geom.recipSpatialStep();
    grad(w,x) = central14n_stencilop(A,w,x) * geom.recipSpatialStep();
    grad(w,y) = central14n_stencilop(A,w,y) * geom.recipSpatialStep();
    grad(w,z) = central14n_stencilop(A,w,z) * geom.recipSpatialStep();

    return grad;
}

template<typename T>
inline TinyVector<_bz_typename T::T_numtype,3> grad3D4_stencilop(const T& A,
    const UniformCubicGeometry<3>& geom) {
  return TinyVector<_bz_typename T::T_numtype,3>(
    central14_stencilop(A,firstDim) * recip_12 * geom.recipSpatialStep(),
    central14_stencilop(A,secondDim) * recip_12 * geom.recipSpatialStep(),
    central14_stencilop(A,thirdDim) * recip_12 * geom.recipSpatialStep());
}

BZ_NAMESPACE_END

#endif // BZ_ARRAYSTENCILOPS_H

