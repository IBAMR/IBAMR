// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_falling_sphere_inactivate_structures
#define included_falling_sphere_inactivate_structures

#include <ibamr/ConstraintIBMethod.h>

#include <ibtk/LDataManager.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <CartesianGridGeometry.h>
#include <IntVector.h>

#include <array>
#include <vector>

inline std::array<bool, 2 * NDIM>
get_open_boundaries(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                    const SAMRAI::hier::IntVector<NDIM>& periodic_shift)
{
    std::array<bool, 2 * NDIM> open_boundaries{};
    if (input_db->keyExists("open_boundaries"))
    {
        input_db->getBoolArray("open_boundaries", open_boundaries.data(), 2 * NDIM);
    }
    for (int d = 0; d < NDIM; ++d)
    {
        if (periodic_shift[d])
        {
            open_boundaries[2 * d] = open_boundaries[2 * d + 1] = false;
        }
    }
    return open_boundaries;
}

inline void
inactivate_structures_at_outlets(IBAMR::ConstraintIBMethod& ib_method,
                                 const SAMRAI::geom::CartesianGridGeometry<NDIM>& grid_geometry,
                                 const std::array<bool, 2 * NDIM>& open_boundaries,
                                 const double radius,
                                 const int level_number)
{
    auto* l_data_manager = ib_method.getLDataManager();
    const double* const x_lower = grid_geometry.getXLower();
    const double* const x_upper = grid_geometry.getXUpper();
    std::vector<int> inactive_structures;
    for (const int structure_id : l_data_manager->getLagrangianStructureIDs(level_number))
    {
        if (!l_data_manager->getLagrangianStructureIsActivated(structure_id, level_number))
        {
            continue;
        }
        // Read live positions since ConstraintIBMethod's cached center can lag by a time step.
        const auto structure_com = l_data_manager->computeLagrangianStructureCenterOfMass(structure_id, level_number);
        for (int d = 0; d < NDIM; ++d)
        {
            if ((open_boundaries[2 * d] && structure_com[d] - radius < x_lower[d]) ||
                (open_boundaries[2 * d + 1] && structure_com[d] + radius > x_upper[d]))
            {
                inactive_structures.push_back(structure_id);
                break;
            }
        }
    }
    if (!inactive_structures.empty())
    {
        l_data_manager->inactivateLagrangianStructures(inactive_structures, level_number);
    }
}

#endif
