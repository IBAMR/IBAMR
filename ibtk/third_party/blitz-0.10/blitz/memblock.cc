/***************************************************************************
 * blitz/memblock.cc      MemoryBlock<T> and MemoryBlockReference<T> methods
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
#ifndef BZ_MEMBLOCK_CC
#define BZ_MEMBLOCK_CC

#include <blitz/numtrait.h>

BZ_NAMESPACE(blitz)

template<typename P_type>
void MemoryBlock<P_type>::deallocate()
{
  const sizeType byteWidth= simdTypes<P_type>::byteWidth;
  const int cacheLineSize = BZ_L1_CACHE_LINE_SIZE;

#ifdef BZ_ALIGN_BLOCKS_ON_CACHELINE_BOUNDARY
  const sizeType minLengthToAlign =
    BZ_CACHE_LINES_TO_ALIGN*cacheLineSize;
#else
  const sizeType minLengthToAlign = blitz::huge(size_t());
#endif

  const sizeType numBytes = length_ * sizeof(T_type);

  // depending on the situation, the block was allocated as any of
  // three different types, so we need to take care to get it right
  // when freeing

  if (!allocatedByUs_)
    // the block was allocated by someone else, so it should just be freed
    delete [] dataBlockAddress_;

  else {
    // we allocated the block, so it may have been shifted
    if ((numBytes < minLengthToAlign) || (cacheLineSize%byteWidth !=0)) {
      // block was not cache line shifted
      delete [] dBA_tv_;
    }
    else {
      // block was cache line shifted, manually call destructors
      if (!NumericTypeTraits<T_type>::hasTrivialCtor)
	for (sizeType i=0; i < length_; ++i)
	  data_[i].~T_type();
      delete [] dBA_char_;
    }
  }
}


/** Allocate a memory block. If we're allocating a large array it may
    be desireable for performance reasons to force the array to start
    on a cache line boundary.  We do this by allocating a little more
    memory than necessary, then shifting the pointer to the next cache
    line boundary. Note that this only works if the SIMD byte width of
    P_type fits evenly into the cache line size, otherwise we may have
    to allocate a lot more memory to be sure to get to a common
    multiple of both the cache line size and the size of the
    object. This can easily be the case for multicomponent containers,
    and in this case, nothing is done.x  */
template<typename P_type>
inline void MemoryBlock<P_type>::allocate(sizeType length)
{
    TAU_TYPE_STRING(p1, "MemoryBlock<T>::allocate() [T="
        + CT(P_type) + "]");
    TAU_PROFILE(p1, "void ()", TAU_BLITZ);

    const sizeType vecWidth= simdTypes<P_type>::vecWidth;
    const sizeType byteWidth= simdTypes<P_type>::byteWidth;
    const int cacheLineSize = BZ_L1_CACHE_LINE_SIZE;

    BZASSERT(length%vecWidth==0);

#ifdef BZ_ALIGN_BLOCKS_ON_CACHELINE_BOUNDARY
    const sizeType minLengthToAlign =
      BZ_CACHE_LINES_TO_ALIGN*cacheLineSize;
#else
    const sizeType minLengthToAlign = blitz::huge(size_t());
#endif

    const sizeType numBytes = length * sizeof(T_type);

    //cout << "Size: " << sizeof(P_type) << ", alignment " << __alignof__(P_type) << endl;

    if ((numBytes < minLengthToAlign) || (cacheLineSize%byteWidth !=0)) {
      // no shifting to cache line
      dBA_tv_ =
	new typename simdTypes<P_type>::vecType[length/vecWidth];
      data_= dataBlockAddress_;
    }
    else {
      // shift to cache line

      dBA_char_ = new char[numBytes + cacheLineSize + 1];

      // Shift to the next cache line boundary.
      diffType offsetToCacheLine =
	diffType(dBA_char_) % cacheLineSize;
      diffType shift = (offsetToCacheLine == 0) ? 0 :
	(cacheLineSize - offsetToCacheLine);
      data_char_ = dBA_char_ + shift;

      // Use placement new to construct types with nontrival ctors
      if (!NumericTypeTraits<T_type>::hasTrivialCtor) {
	for (sizeType i=0; i < length; ++i)
	  new(&data_[i]) T_type;
      }
    }

    allocatedByUs_ = true;

#ifdef BZ_DEBUG_LOG_ALLOCATIONS
    cout << "MemoryBlock: allocated " << setw(8) << length
         << " at " << ((void *)dataBlockAddress_) << endl;
#endif

    BZASSERT(isVectorAligned(data_));
}


BZ_NAMESPACE_END

#endif // BZ_MEMBLOCK_CC
