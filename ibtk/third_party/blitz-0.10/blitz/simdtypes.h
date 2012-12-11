#ifndef BZ_SIMDTYPES_H
#define BZ_SIMDTYPES_H

#include <blitz/blitz.h>
#include <stdint.h>

BZ_NAMESPACE(blitz)

template<class T,int N> class TinyVector;

/** This metaprogram returns the number of bits necessary to fit the
    specified number. The metaprogram result for I is obtained as
    _bz_meta_bitwidth<I,0>::width. */
template<size_t I, size_t L> struct _bz_meta_bitwidth {
  static const size_t width = _bz_meta_bitwidth<(I>>1), L+1>::width;
};

template<size_t L> struct _bz_meta_bitwidth<0,L> {
  static const size_t width = L;
};


/** Helper class that defines the width of the simd instructions for a
    given type. It also defines a type that is a TinyVector of the
    appropriate length to fill the simd width. This is used as the
    "evaluation type" to facilitate automatic vectorization. Because
    TinyVectors are aligned on BZ_SIMD_WIDTH, an array of vecType will
    be dense in memory, which is required to reinterpret an array of T
    to an array of vecType. \todo How do we check that the type fits
    evenly into the simd width? */
template<typename T> class simdTypes {
public:
  /// SIMD width of type in bytes (sizeof(T) if simd width does not
  /// fit a T)
  static const size_t byteWidth =
    BZ_SIMD_WIDTH>sizeof(T) ? BZ_SIMD_WIDTH : sizeof(T);

  /// SIMD width of types in number of elements.
  static const size_t vecWidth =
    BZ_SIMD_WIDTH>sizeof(T) ? BZ_SIMD_WIDTH/sizeof(T) : 1;

  /// TinyVector type of T that fills the simd width.
  typedef TinyVector<T, vecWidth> vecType;

  /** Test if a pointer to T is simd aligned. */
  static inline bool isVectorAligned(const T* restrict pointer)
  { return (uintptr_t)((void*)pointer) % BZ_SIMD_WIDTH == 0; }

  /** Return number of elements from pointer to next simd width
      boundary. This is used to figure out how many scalar operations
      need to be done before beginning vectorized operations. */
  static inline diffType offsetToAlignment(const T* restrict pointer) {
    const uintptr_t m = (uintptr_t)((void*)pointer) & (byteWidth-1);
    return m ? (byteWidth - m)/sizeof(T) : 0;
}

  /** Return a length which has been padded to next larger even SIMD
      width. */
  static inline size_t paddedLength(size_t length) {
    return (length & (vecWidth-1)) ?
      (length & ~(vecWidth-1)) + vecWidth : length;
  };
};

/** General function just forwards to the simdTypes class. */
template <typename T>
inline bool isVectorAligned(const T* restrict pointer) {
  return simdTypes<T>::isVectorAligned(pointer); }


BZ_NAMESPACE_END

#endif
