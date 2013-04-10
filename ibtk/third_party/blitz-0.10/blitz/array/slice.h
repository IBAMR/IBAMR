// -*- C++ -*-
/***************************************************************************
 * blitz/array/slice.h    Helper classes for slicing arrays
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
#ifndef BZ_ARRAYSLICE_H
#define BZ_ARRAYSLICE_H

// #ifndef BZ_ARRAY_H
//  #error <blitz/array/slice.h> must be included via <blitz/array.h>
// #endif

#include <blitz/range.h>

BZ_NAMESPACE(blitz)

// Forward declarations
template<typename T, int N>
class Array;

class nilArraySection { };

template<typename T>
class ArraySectionInfo {
public:
    static const int isValidType = 0, rank = 0, isPick = 0;
};

template<>
class ArraySectionInfo<Range> {
public:
    static const int isValidType = 1, rank = 1, isPick = 0;
};

template<>
class ArraySectionInfo<int> {
public:
    static const int isValidType = 1, rank = 0, isPick = 0;
};

template<>
class ArraySectionInfo<nilArraySection> {
public:
    static const int isValidType = 1, rank = 0, isPick = 0;
};

template<typename T_numtype, typename T1, typename T2 = nilArraySection,
    class T3 = nilArraySection, typename T4 = nilArraySection,
    class T5 = nilArraySection, typename T6 = nilArraySection,
    class T7 = nilArraySection, typename T8 = nilArraySection,
    class T9 = nilArraySection, typename T10 = nilArraySection,
    class T11 = nilArraySection>
class SliceInfo {
public:
    static const int
        numValidTypes = ArraySectionInfo<T1>::isValidType
                      + ArraySectionInfo<T2>::isValidType
                      + ArraySectionInfo<T3>::isValidType
                      + ArraySectionInfo<T4>::isValidType
                      + ArraySectionInfo<T5>::isValidType
                      + ArraySectionInfo<T6>::isValidType
                      + ArraySectionInfo<T7>::isValidType
                      + ArraySectionInfo<T8>::isValidType
                      + ArraySectionInfo<T9>::isValidType
                      + ArraySectionInfo<T10>::isValidType
                      + ArraySectionInfo<T11>::isValidType,

        rank          = ArraySectionInfo<T1>::rank
                      + ArraySectionInfo<T2>::rank
                      + ArraySectionInfo<T3>::rank
                      + ArraySectionInfo<T4>::rank
                      + ArraySectionInfo<T5>::rank
                      + ArraySectionInfo<T6>::rank
                      + ArraySectionInfo<T7>::rank
                      + ArraySectionInfo<T8>::rank
                      + ArraySectionInfo<T9>::rank
                      + ArraySectionInfo<T10>::rank
                      + ArraySectionInfo<T11>::rank,

        isPick        = ArraySectionInfo<T1>::isPick
                      + ArraySectionInfo<T2>::isPick
                      + ArraySectionInfo<T3>::isPick
                      + ArraySectionInfo<T4>::isPick
                      + ArraySectionInfo<T5>::isPick
                      + ArraySectionInfo<T6>::isPick
                      + ArraySectionInfo<T7>::isPick
                      + ArraySectionInfo<T8>::isPick
                      + ArraySectionInfo<T9>::isPick
                      + ArraySectionInfo<T10>::isPick
                      + ArraySectionInfo<T11>::isPick;

    typedef Array<T_numtype,numValidTypes> T_array;
    typedef Array<T_numtype,rank> T_slice;
};


BZ_NAMESPACE_END

#endif // BZ_ARRAYSLICE_H
