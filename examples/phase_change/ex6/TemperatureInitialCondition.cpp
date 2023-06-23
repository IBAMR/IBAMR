// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "TemperatureInitialCondition.h"

#include <SAMRAI_config.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

TemperatureInitialCondition::TemperatureInitialCondition(const std::string& object_name,
                                                         const double initial_liquid_solid_interface_position,
                                                         const double initial_liquid_temperature,
                                                         const double initial_solid_temperature)
    : d_object_name(object_name),
      d_initial_liquid_solid_interface_position(initial_liquid_solid_interface_position),
      d_initial_liquid_temperature(initial_liquid_temperature),
      d_initial_solid_temperature(initial_solid_temperature)
{
    // intentionally blank
    return;
} // TemperatureInitialCondition

bool
TemperatureInitialCondition::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
TemperatureInitialCondition::setDataOnPatch(const int data_idx,
                                            Pointer<Variable<NDIM> > /*var*/,
                                            Pointer<Patch<NDIM> > patch,
                                            const double /*data_time*/,
                                            const bool initial_time,
                                            Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    // Set the temperature function throughout the domain
    if (initial_time)
    {
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(data_idx);

        Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* patch_X_lower = patch_geom->getXLower();
        const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
        const double* const patch_dx = patch_geom->getDx();

        for (Box<NDIM>::Iterator it(patch_box); it; it++)
        {
            CellIndex<NDIM> ci(it());

            // Get physical coordinates
            IBTK::Vector coord = IBTK::Vector::Zero();
            for (int d = 0; d < NDIM; ++d)
            {
                coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
            }

            if (coord[1] >= d_initial_liquid_solid_interface_position)
            {
                (*D_data)(ci) = d_initial_liquid_temperature;
            }
            else
            {
                (*D_data)(ci) = d_initial_solid_temperature;
            }
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
