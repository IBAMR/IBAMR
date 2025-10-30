// ---------------------------------------------------------------------
//
// Copyright (c) 2025 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Verify that nested tbox::Array<tbox::Array<bool>> works correctly

#include <ibtk/IBTKInit.h>

#include <tbox/Array.h>

#include <fstream>
#include <string>
#include <vector>

#include <ibtk/app_namespaces.h>

int
main()
{
    std::ofstream output("output");

    for (int i = 0; i < 100000; ++i)
    {
        Array<Array<bool> > touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            touches_regular_bdry[axis].resizeArray(2);
            touches_periodic_bdry[axis].resizeArray(2);
            for (int upperlower = 0; upperlower < 2; ++upperlower)
            {
                touches_regular_bdry[axis][upperlower] = (i + 1) % (upperlower + 1) == NDIM;
                touches_periodic_bdry[axis][upperlower] = (i + 1) % (upperlower + 1) == NDIM;
            }

            Array<Array<bool> > touches_regular_bdry_2 = touches_regular_bdry;
            Array<Array<bool> > touches_regular_bdry_3 = touches_regular_bdry;
            Array<Array<bool> > touches_regular_bdry_4 = touches_regular_bdry;
            Array<Array<bool> > touches_regular_bdry_5 = touches_regular_bdry;
        }
    }

    output << "total number of Array<int> allocations = " << Array<int>::getNumberOfAllocations() << std::endl;
    output << "total number of Array<bool> allocations = " << Array<bool>::getNumberOfAllocations() << std::endl;
    output << "total number of Array<Array<bool>> allocations = " << Array<Array<bool> >::getNumberOfAllocations()
           << std::endl;
}
