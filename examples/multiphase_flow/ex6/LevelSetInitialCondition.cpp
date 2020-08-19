// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "LevelSetInitialCondition.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetInitialCondition::LevelSetInitialCondition(const std::string& object_name,
                                                   FilmInterface init_film,
                                                   CircularInterface init_circle)
    : d_object_name(object_name), d_init_film(init_film), d_init_circle(init_circle)
{
    // intentionally blank
    return;
} // LevelSetInitialCondition

LevelSetInitialCondition::~LevelSetInitialCondition()
{
    // intentionally blank
    return;
} // ~LevelSetInitialCondition

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
    // Set the level set function throughout the domain
    if (initial_time)
    {
        // Get the parameters for the interface
        const double& R = d_init_circle.R;
        const IBTK::Vector& X0 = d_init_circle.X0;
        const double& film_height = d_init_film.height;
        const int height_dim = NDIM - 1;

        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(data_idx);
        for (Box<NDIM>::Iterator it(patch_box); it; it++)
        {
            CellIndex<NDIM> ci(it());

            // Get physical coordinates
            IBTK::Vector coord = IBTK::Vector::Zero();
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();
            for (int d = 0; d < NDIM; ++d)
            {
                coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
            }

            // Distance from the bubble
            const double distance_bubble =
                std::sqrt(std::pow((coord[0] - X0(0)), 2.0) + std::pow((coord[1] - X0(1)), 2.0)
#if (NDIM == 3)
                          + std::pow((coord[2] - X0(2)), 2.0)
#endif
                              ) -
                R;

            // Distance from the film
            const double distance_film = coord[height_dim] - film_height;

            if (distance_film <= 0)
            {
                // If within the film, set LS as the negative distance away
                (*D_data)(ci) = distance_film;
            }
            else if (distance_bubble <= 0)
            {
                // If within the bubble, again set the LS as the negative distance away
                (*D_data)(ci) = distance_bubble;
            }
            else
            {
                // Otherwise, set the distance as the minimum between the two
                (*D_data)(ci) = std::min(distance_bubble, distance_film);
            }
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
