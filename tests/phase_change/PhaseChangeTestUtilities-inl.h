// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_PhaseChangeTestUtilities_inl
#define included_PhaseChangeTestUtilities_inl

#include <ibamr/config.h>

template <class Integrator>
int
RegridCountingIntegrator<Integrator>::getConfigurationResetCount() const
{
    return d_configuration_reset_count;
}

template <class Integrator>
int
RegridCountingIntegrator<Integrator>::getMeshChangeCount() const
{
    return d_mesh_change_count;
}

template <class Integrator>
void
RegridCountingIntegrator<Integrator>::resetHierarchyConfigurationSpecialized(
    SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
    int coarsest_level,
    int finest_level)
{
    Integrator::resetHierarchyConfigurationSpecialized(hierarchy, coarsest_level, finest_level);
    ++d_configuration_reset_count;
    std::ostringstream boxes;
    if (hierarchy->getFinestLevelNumber() > 0)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(1);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            boxes << level->getPatch(p())->getBox();
        }
    }
    if (!d_fine_boxes.empty() && boxes.str() != d_fine_boxes)
    {
        ++d_mesh_change_count;
    }
    d_fine_boxes = boxes.str();
}

#endif
