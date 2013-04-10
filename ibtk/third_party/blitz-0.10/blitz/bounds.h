#ifndef BZ_BOUNDS_H
#define BZ_BOUNDS_H

#include <blitz/blitz.h>

BZ_NAMESPACE(blitz)

struct bounds {
    static int compute_ascending(const int BZ_DEBUG_PARAM(rank),
                                 const int ascending1, const int ascending2)
    {
        // The value INT_MIN indicates that there are no arrays
        // in a subtree of the expression.  This logic returns
        // whichever ascending is available.  If there are two
        // conflicting ascending values, this is an error.

        if (ascending1 == ascending2)
            return ascending1;
        else if (ascending1 == INT_MIN)
            return ascending2;
        else if (ascending2 == INT_MIN)
            return ascending1;

        BZ_DEBUG_MESSAGE("Two array operands have different"
            << endl << "ascending flags: for rank " << rank
            << ", the flags are " << ascending1 << " and "
            << ascending2 << endl);
        BZ_PRE_FAIL;
        return 0;
    }

    static int compute_ordering(const int BZ_DEBUG_PARAM(rank),
                                const int order1, const int order2)
    {
        // The value INT_MIN indicates that there are no arrays
        // in a subtree of the expression.  This logic returns
        // whichever ordering is available.  If there are two
        // conflicting ordering values, this is an error.

        if (order1 == order2)
            return order1;
        else if (order1 == INT_MIN)
            return order2;
        else if (order2 == INT_MIN)
            return order1;

        BZ_DEBUG_MESSAGE("Two array operands have different"
            << endl << "orders: for rank " << rank << ", the orders are "
            << order1 << " and " << order2 << endl);
        BZ_PRE_FAIL;
        return 0;
    }

    static int compute_lbound(const int BZ_DEBUG_PARAM(rank),
                              const int lbound1, const int lbound2)
    {
        // The value INT_MIN indicates that there are no arrays
        // in a subtree of the expression.  This logic returns
        // whichever lbound is available.  If there are two
        // conflicting lbound values, this is an error.

        if (lbound1 == lbound2)
            return lbound1;
        else if (lbound1 == INT_MIN)
            return lbound2;
        else if (lbound2 == INT_MIN)
            return lbound1;

        BZ_DEBUG_MESSAGE("Two array operands have different"
            << endl << "lower bounds: in rank " << rank << ", the bounds are "
            << lbound1 << " and " << lbound2 << endl);
        BZ_PRE_FAIL;
        return 0;
    }

    static int compute_ubound(const int BZ_DEBUG_PARAM(rank),
                              const int ubound1, const int ubound2)
    {
        // The value INT_MAX indicates that there are no arrays
        // in a subtree of the expression.  This logic returns
        // whichever ubound is available.  If there are two
        // conflicting ubound values, this is an error.

        if (ubound1 == ubound2)
            return ubound1;
        else if (ubound1 == INT_MAX)
            return ubound2;
        else if (ubound2 == INT_MAX)
            return ubound1;

        BZ_DEBUG_MESSAGE("Two array operands have different"
            << endl << "upper bounds: in rank " << rank << ", the bounds are "
            << ubound1 << " and " << ubound2 << endl);
        BZ_PRE_FAIL;
        return 0;
    }
};

BZ_NAMESPACE_END

#endif
