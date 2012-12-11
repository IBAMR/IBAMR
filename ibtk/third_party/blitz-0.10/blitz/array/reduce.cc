/***************************************************************************
 * blitz/array/reduce.cc  Array reductions.
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
#ifndef BZ_ARRAYREDUCE_H
 #error <blitz/array/reduce.cc> must be included via <blitz/array/reduce.h>
#endif

BZ_NAMESPACE(blitz)

template<typename T_expr, typename T_reduction>
_bz_typename T_reduction::T_resulttype
_bz_ArrayExprFullReduce(T_expr expr, T_reduction reduction)
{
#ifdef BZ_TAU_PROFILING
    // Tau profiling code.  Provide Tau with a pretty-printed version of
    // the expression.
    static BZ_STD_SCOPE(string) exprDescription;
    if (!exprDescription.length())      // faked static initializer
    {
        exprDescription = T_reduction::name();
        exprDescription += "(";
        prettyPrintFormat format(true);   // Terse mode on
        expr.prettyPrint(exprDescription, format);
        exprDescription += ")";
    }
    TAU_PROFILE(" ", exprDescription, TAU_BLITZ);
#endif // BZ_TAU_PROFILING

    return _bz_reduceWithIndexTraversal(expr, reduction);

#ifdef BZ_NOT_IMPLEMENTED_FLAG
    if ((T_expr::numIndexPlaceholders > 0) || (T_reduction::needIndex))
    {
        // The expression involves index placeholders, so have to
        // use index traversal rather than stack traversal.
        return reduceWithIndexTraversal(expr, reduction);
    }
    else {
        // Use a stack traversal
        return reduceWithStackTraversal(expr, reduction);
    }
#endif
}

template <typename T_index> struct _bz_IndexingVariant;

template <>
struct _bz_IndexingVariant<int> {
    template <int N>
    static int index(const TinyVector<int,N>& ind,const int i) {
        return ind[i];
    }
};

template <int N>
struct _bz_IndexingVariant<TinyVector<int,N> > {
    static const TinyVector<int,N>& index(const TinyVector<int,N>& ind,const int) {
        return ind;
    }
};

template<typename T_index, typename T_expr, typename T_reduction>
_bz_typename T_reduction::T_resulttype
_bz_reduceWithIndexTraversalGeneric(T_expr expr, T_reduction reduction)
{
    // This is optimized assuming C-style arrays.

    const int rank = T_expr::rank_;

    TinyVector<int,T_expr::rank_> index, first, last;

    unsigned long count = 1;

    for (int i=0; i < rank; ++i) {
        first(i) = expr.lbound(i);
        last(i) = expr.ubound(i) + 1;
        index(i) = first(i);
        count *= last(i) - first(i);
    }

    const int maxRank = rank - 1;
    const int lastlbound = expr.lbound(maxRank);
    const int lastubound = expr.ubound(maxRank);

    const int lastIndex = lastubound + 1;

    typedef _bz_IndexingVariant<T_index> adapter;

    _bz_ReduceReset<T_reduction::needIndex,T_reduction::needInit> reset;
    reset(reduction,first,expr);

    while(true) {
        for (index[maxRank]=lastlbound;index[maxRank]<lastIndex;++index[maxRank])
            if (!reduction(expr(index),adapter::index(index,maxRank)))
                return reduction.result(count);

        int j = rank-2;
        for (;j>=0;--j) {
            index(j+1) = first(j+1);
            ++index(j);
            if (index(j) < last(j))
                break;
        }

        if (j<0)
            return reduction.result(count);
    }
}

template<typename T_expr, typename T_reduction>
_bz_typename T_reduction::T_resulttype
_bz_reduceWithIndexTraversal(T_expr expr, T_reduction reduction)
{
    return _bz_reduceWithIndexTraversalGeneric<int>(expr,reduction);
}

// This version is for reductions that require a vector of index positions.

template<typename T_expr, typename T_reduction>
_bz_typename T_reduction::T_resulttype
_bz_reduceWithIndexVectorTraversal(T_expr expr, T_reduction reduction)
{
    // We are doing minIndex/maxIndex, so initialize with lower bound
    return _bz_reduceWithIndexTraversalGeneric<TinyVector<int,T_expr::rank_> >(expr,reduction);
}

BZ_NAMESPACE_END

