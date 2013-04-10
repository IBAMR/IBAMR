#ifndef BZ_RANGE_CC
#define BZ_RANGE_CC


#include <blitz/range.h>
#include <blitz/array/ops.h>

BZ_NAMESPACE(blitz)


inline Range::T_numtype Range::operator[](int i) const
{
  return first_ + i * stride_;
}

inline RectDomain<Range::rank_> Range::domain() const
{
  return RectDomain<rank_>(TinyVector<int,1>(0),
			   TinyVector<int,1>(ubound(0)));
}

inline bool Range::assertInRange(const T_index& BZ_DEBUG_PARAM(index)) const
{
  BZPRECHECK((index[0]<=ubound(0)) && (index[0]>=0),
	     "Range index out of range: " << index[0]
	     << endl << "Lower bounds: 0" << endl
	     <<         "Length:      " << length() << endl);
  return true;
}

BZ_NAMESPACE_END

#endif
