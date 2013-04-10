// -*- C++ -*-
/***************************************************************************
 * blitz/tuning.h      Platform-specific code tuning
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

#ifndef BZ_TUNING_H
#define BZ_TUNING_H

// These estimates should be conservative (i.e. underestimate the
// cache sizes). \todo these can be const ints instead of macros.
#define BZ_L1_CACHE_ESTIMATED_SIZE    32768
#define BZ_L2_CACHE_ESTIMATED_SIZE    6291456
// This will work for 32, 16 also
#define BZ_L1_CACHE_LINE_SIZE         64
#define BZ_CACHE_LINES_TO_ALIGN       16

#undef  BZ_PARTIAL_LOOP_UNROLL
#define BZ_PASS_EXPR_BY_VALUE
#undef  BZ_PTR_INC_FASTER_THAN_INDIRECTION
#define BZ_MANUAL_VECEXPR_COPY_CONSTRUCTOR
#undef  BZ_KCC_COPY_PROPAGATION_KLUDGE
#undef  BZ_ALTERNATE_FORWARD_BACKWARD_TRAVERSALS
#undef  BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
#define BZ_INLINE_GROUP1
#define BZ_INLINE_GROUP2
#define BZ_COLLAPSE_LOOPS
#define BZ_USE_FAST_READ_ARRAY_EXPR
#define BZ_ARRAY_EXPR_USE_COMMON_STRIDE
#undef  BZ_ARRAY_SPACE_FILLING_TRAVERSAL
#undef  BZ_ARRAY_FAST_TRAVERSAL_UNROLL
#undef  BZ_ARRAY_STACK_TRAVERSAL_CSE_AND_ANTIALIAS
#undef  BZ_ARRAY_STACK_TRAVERSAL_UNROLL
#define BZ_ARRAY_2D_STENCIL_TILING
#define BZ_ARRAY_2D_STENCIL_TILE_SIZE       128
#undef  BZ_INTERLACE_ARRAYS
#define  BZ_ALIGN_BLOCKS_ON_CACHELINE_BOUNDARY
#define BZ_FAST_COMPILE
#define BZ_TV_EVALUATE_UNROLL_LENGTH 0
#define BZ_MAX_BITS_FOR_BINARY_UNROLL 8
#define BZ_VECTORIZED_LOOP_WIDTH 32


#ifndef BZ_DISABLE_NEW_ET
 #define BZ_NEW_EXPRESSION_TEMPLATES
#endif

#ifdef BZ_FAST_COMPILE
#define BZ_ETPARMS_CONSTREF
#define BZ_NO_INLINE_ET
#endif

// possibly overridden by specific compilers below
#define _bz_forceinline inline
#define _bz_inline_et inline


/*
 * Platform-specific tuning
 */

#ifdef _CRAYT3E
 // The backend compiler on the T3E does a better job of
 // loop unrolling.
 #undef BZ_PARTIAL_LOOP_UNROLL
 #undef BZ_ARRAY_FAST_TRAVERSAL_UNROLL
 #undef BZ_ARRAY_STACK_TRAVERSAL_UNROLL
#endif

#ifdef __INTEL_COMPILER
 // icpc does not vectorize the unrolled loop so this is def. bad
 #define BZ_TV_EVALUATE_UNROLL_LENGTH 0

 // defines for inlining
 #undef _bz_forceinline
 #undef _bz_inline_et
 #define _bz_forceinline __forceinline
 #define _bz_inline_et __forceinline

#else // need this since icpc also defines __GNUC__
#ifdef __GNUC__
 // The egcs compiler does a good job of loop unrolling, if
 // -funroll-loops is used.
 #undef BZ_PARTIAL_LOOP_UNROLL
 #undef BZ_ARRAY_FAST_TRAVERSAL_UNROLL
 #undef BZ_ARRAY_STACK_TRAVERSAL_UNROLL
#endif
#endif

#ifdef  BZ_DISABLE_KCC_COPY_PROPAGATION_KLUDGE
 #undef BZ_KCC_COPY_PROPAGATION_KLUDGE
#endif

#ifdef  BZ_INLINE_GROUP1
 #define _bz_inline1 inline
#else
 #define _bz_inline1
#endif

#ifdef  BZ_INLINE_GROUP2
 #define _bz_inline2 inline
#else
 #define _bz_inline2
#endif

// override definitions above
#ifdef  BZ_NO_INLINE_ET
 #undef _bz_inline_et
 #define _bz_inline_et
#endif

#ifdef  BZ_ETPARMS_CONSTREF
 #define BZ_ETPARM(X) const X&
#else
 #define BZ_ETPARM(X) X
#endif

#ifdef __DECCXX
 // The DEC cxx compiler has problems with loop unrolling
 // because of aliasing.  Loop unrolling and anti-aliasing
 // is done by Blitz++.

  #define  BZ_PARTIAL_LOOP_UNROLL
  #define  BZ_ARRAY_STACK_TRAVERSAL_CSE_AND_ANTIALIAS
  #define  BZ_ARRAY_STACK_TRAVERSAL_UNROLL
#endif

/*
 * BZ_NO_PROPAGATE(X) prevents the compiler from performing
 * copy propagation on a variable.  This is used for loop
 * unrolling to prevent KAI C++ from rearranging the
 * ordering of memory accesses.
 */

#define BZ_NO_PROPAGATE(X)   X

#ifdef __KCC
#ifdef BZ_USE_NO_PROPAGATE
    extern "C" int __kai_apply(const char*, ...);

    #undef  BZ_NO_PROPAGATE(X)
    #define BZ_NO_PROPAGATE(X)  __kai_apply("(%a)",&X)
#endif
#endif

#endif // BZ_TUNING_H
