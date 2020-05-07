// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "USourceFunction.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <IBAMR_config.h>

#include <SAMRAI_config.h>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

USourceFunction::USourceFunction(const string& object_name,
                                 Pointer<GridGeometry<NDIM> > grid_geom,
                                 Pointer<Database> input_db)
    : CartGridFunction(object_name), d_object_name(object_name), d_grid_geom(grid_geom)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(grid_geom);
#endif

    d_analytical_pr_gradient = input_db->getDouble("AnalyticalPressureGradient");
    return;
} // USourceFunction

USourceFunction::~USourceFunction()
{
    // intentionally blank
    return;
} // ~USourceFunction

void
USourceFunction::setDataOnPatch(const int data_idx,
                                Pointer<Variable<NDIM> > /*var*/,
                                Pointer<Patch<NDIM> > patch,
                                const double /*data_time*/,
                                const bool /*initial_time*/,
                                Pointer<PatchLevel<NDIM> > /*level*/)
{
    Pointer<SideData<NDIM, double> > F_data = patch->getPatchData(data_idx);

    const Box<NDIM>& patch_box = patch->getBox();
    for (unsigned int d = 0; d < NDIM; d++)
    {
        for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, d)); it; it++)
        {
            SideIndex<NDIM> s_i(it(), d, SideIndex<NDIM>::Lower);

            if (d == 0)
            {
                (*F_data)(s_i) = d_analytical_pr_gradient;
            }
            else
            {
                (*F_data)(s_i) = 0.0;
            }
        }
    }
    return;
} // setDataOnPatch
