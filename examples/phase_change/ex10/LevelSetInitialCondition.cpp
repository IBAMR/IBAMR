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
#include "LevelSetInitialCondition.h"

#include <SAMRAI_config.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetInitialCondition::LevelSetInitialCondition(const std::string& object_name,
                                                   const double initial_horizontal_interface_position,
                                                   std::vector<std::pair<double, IBTK::VectorNd> > bubbles_position,
                                                   const bool center_bubble_required)
    : d_object_name(object_name),
      d_initial_horizontal_interface_position(initial_horizontal_interface_position),
      d_bubbles_position(bubbles_position),
      d_center_bubble_required(center_bubble_required)

{
    // intentionally blank
    return;
} // LevelSetInitialCondition

bool
LevelSetInitialCondition::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
LevelSetInitialCondition::setDataOnPatch(const int data_idx,
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

            // left bubble.
            std::pair<double, IBTK::VectorNd> left_bubble = d_bubbles_position[0];
            const double left_bubble_radius = left_bubble.first;
            IBTK::Vector left_bubble_coord = left_bubble.second;
            double left_bubble_signed_distance = std::sqrt(std::pow(coord[0] - left_bubble_coord[0], 2.0) +
                                                           std::pow(coord[1] - left_bubble_coord[1], 2.0));
            left_bubble_signed_distance = -(left_bubble_radius - left_bubble_signed_distance);

            // center bubble.
            std::pair<double, IBTK::VectorNd> center_bubble = d_bubbles_position[1];
            const double center_bubble_radius = center_bubble.first;
            IBTK::Vector center_bubble_coord = center_bubble.second;
            double center_bubble_signed_distance = std::sqrt(std::pow(coord[0] - center_bubble_coord[0], 2.0) +
                                                             std::pow(coord[1] - center_bubble_coord[1], 2.0));
            center_bubble_signed_distance = -(center_bubble_radius - center_bubble_signed_distance);

            // right bubble.
            std::pair<double, IBTK::VectorNd> right_bubble = d_bubbles_position[2];
            const double right_bubble_radius = right_bubble.first;
            IBTK::Vector right_bubble_coord = right_bubble.second;
            double right_bubble_signed_distance = std::sqrt(std::pow(coord[0] - right_bubble_coord[0], 2.0) +
                                                            std::pow(coord[1] - right_bubble_coord[1], 2.0));
            right_bubble_signed_distance = -(right_bubble_radius - right_bubble_signed_distance);

            // Initial horizontal position
            const double horizontal_interface_signed_distance = d_initial_horizontal_interface_position - coord[1];

            if (d_center_bubble_required)
            {
                (*D_data)(ci) = std::min({ left_bubble_signed_distance,
                                           center_bubble_signed_distance,
                                           right_bubble_signed_distance,
                                           horizontal_interface_signed_distance });
            }
            else
            {
                (*D_data)(ci) = std::min({ left_bubble_signed_distance,
                                           right_bubble_signed_distance,
                                           horizontal_interface_signed_distance });
            }
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
