#ifndef BZ_ET_FORWARD_H
#define BZ_ET_FORWARD_H

#include <blitz/blitz.h>

// forward declarations of ET classes

BZ_NAMESPACE(blitz)

template<typename T> class ETBase;

// the actual array/vector classes
template<typename P_numtype, int N_rank> class Array;
template<typename P_numtype, int N_length> class TinyVector;
template<typename P_numtype, int N_rows, int N_columns> class TinyMatrix;

// and their fast iterators
template<typename P_numtype, int N_rank> class FastArrayIterator;
template<typename P_numtype, int N_rank> class FastArrayCopyIterator;
template<typename P_numtype, int N_length> class FastTV2Iterator;
template<typename P_numtype, int N_length> class FastTV2CopyIterator;
template<typename P_numtype, int N_rows,
	 int N_columns> class FastTM2Iterator;
template<typename P_numtype, int N_rows,
	 int N_columns> class FastTM2CopyIterator;

// expression container classes
template<typename P_expr> class _bz_ArrayExpr;
template<typename P_expr, typename P_op> class _bz_ArrayExprUnaryOp;
template<typename P_expr1, typename P_expr2,
	 typename P_op> class _bz_ArrayExprBinaryOp;
template<typename P_expr1, typename P_expr2,
	 typename P_expr3, typename P_op> class _bz_ArrayExprTernaryOp;
template<typename P_expr1, typename P_expr2, typename P_expr3,
	 typename P_expr4, typename P_op> class _bz_ArrayExprQuaternaryOp;
template<typename P_numtype> class _bz_ArrayExprConstant;

// more specialized expression types
template<typename P_expr, int N_map0, int N_map1=0, int N_map2=0,
	 int N_map3=0, int N_map4=0, int N_map5=0, int N_map6=0, int N_map7=0,
	 int N_map8=0, int N_map9=0, int N_map10=0> class ArrayIndexMapping;

template<typename T_expr, int N_index,
	 typename T_reduction> class _bz_ArrayExprReduce;

template<typename P_expr, _bz_typename P_result> class _bz_StencilExpr;

template<int N> class IndexPlaceholder;

template<typename P_expr1, typename P_expr2,
	 typename P_expr3> class _bz_ArrayWhere;

// functors
template<typename P_functor, typename P_expr,
	 typename P_result> class _bz_FunctorExpr;

template<typename P_functor, typename P_expr1,
	 typename P_expr2, typename P_result> class _bz_FunctorExpr2;

template<typename P_functor, typename P_expr1,
	 typename P_expr2, typename P_expr3,
	 class P_result> class _bz_FunctorExpr3;

class LeviCivita;

class Range;

BZ_NAMESPACE_END


#endif
