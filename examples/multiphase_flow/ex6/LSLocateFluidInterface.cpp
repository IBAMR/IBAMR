// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// SAMRAI INCLUDES
#include "ibtk/samrai_compatibility_names.h"
#include <ibtk/HierarchyMathOps.h>

#include "LSLocateFluidInterface.h"
#include "SAMRAIBox.h"
#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAICartesianPatchGeometry.h"
#include "SAMRAICellData.h"
#include "SAMRAICellIndex.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIHierarchyCellDataOpsReal.h"
#include "SAMRAIIndex.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIVariableDatabase.h"

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

void
callLSLocateFluidInterfaceCallbackFunction(int D_idx,
                                           SAMRAIPointer<HierarchyMathOps> hier_math_ops,
                                           double time,
                                           bool initial_time,
                                           void* ctx)
{
    // Set the level set information
    static LSLocateFluidInterface* ptr_LSLocateFluidInterface = static_cast<LSLocateFluidInterface*>(ctx);
    ptr_LSLocateFluidInterface->setLevelSetPatchData(D_idx, hier_math_ops, time, initial_time);

    return;
} // callLSLocateFluidInterfaceCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////
LSLocateFluidInterface::LSLocateFluidInterface(const std::string& object_name,
                                               SAMRAIPointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                               SAMRAIPointer<SAMRAICellVariable<double>> ls_var,
                                               FilmInterface init_film,
                                               CircularInterface init_circle)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_ls_var(ls_var),
      d_init_film(init_film),
      d_init_circle(init_circle)
{
    // intentionally left blank
    return;
} // LSLocateFluidInterface

LSLocateFluidInterface::~LSLocateFluidInterface()
{
    // intentionally left blank
    return;
}

void
LSLocateFluidInterface::setLevelSetPatchData(int D_idx,
                                             SAMRAIPointer<HierarchyMathOps> hier_math_ops,
                                             double /*time*/,
                                             bool initial_time)
{
    SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // If not the intial time, set the level set to the current value maintained by the integrator
    if (!initial_time)
    {
        SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
        const int ls_current_idx =
            var_db->mapVariableAndContextToIndex(d_ls_var, d_adv_diff_solver->getCurrentContext());
        SAMRAIHierarchyCellDataOpsReal<double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

        hier_cc_data_ops.copyData(D_idx, ls_current_idx);

        return;
    }

    // Set the initial condition for locating the interface
    const double& R = d_init_circle.R;
    const IBTK::Vector& X0 = d_init_circle.X0;
    const double& film_height = d_init_film.height;
    const int height_dim = NDIM - 1;

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
        for (SAMRAIPatchLevel::Iterator p(level); p; p++)
        {
            SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
            const SAMRAIBox& patch_box = patch->getBox();
            SAMRAIPointer<SAMRAICellData<double>> D_data = patch->getPatchData(D_idx);
            for (SAMRAIBox::Iterator it(patch_box); it; it++)
            {
                SAMRAICellIndex ci(it());

                // Get physical coordinates
                IBTK::Vector coord = IBTK::Vector::Zero();
                SAMRAIPointer<SAMRAICartesianPatchGeometry> patch_geom = patch->getPatchGeometry();
                const double* patch_X_lower = patch_geom->getXLower();
                const SAMRAIIndex& patch_lower_idx = patch_box.lower();
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
    }
    return;
} // setLevelSetPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
