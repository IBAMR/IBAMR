#ifndef BZ_EXPR_CC
#define BZ_EXPR_CC

#include <blitz/array/expr.h>

BZ_NAMESPACE(blitz)

template<typename P_numtype>
inline RectDomain<12>
_bz_ArrayExprConstant<P_numtype>::domain() const
{
  TinyVector<int, 12> lb(lbound(0)), ub(ubound(0));
  return RectDomain<12>(lb,ub);
}

BZ_NAMESPACE_END

#endif
