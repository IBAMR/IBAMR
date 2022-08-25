// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/HierarchyMathOps.h>

#include "LSLocateInterface.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

void
callLSLocateInterfaceCallbackFunction(int D_idx,
                                      Pointer<HierarchyMathOps> hier_math_ops,
                                      double time,
                                      bool initial_time,
                                      void* ctx)
{
    // Set the level set information
    static LSLocateInterface* ptr_LSLocateInterface = static_cast<LSLocateInterface*>(ctx);
    ptr_LSLocateInterface->setLevelSetPatchData(D_idx, hier_math_ops, time, initial_time);

    return;
} // callLSLocateInterfaceCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////
LSLocateInterface::LSLocateInterface(const std::string& object_name,
                                     Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                     Pointer<CellVariable<NDIM, double> > ls_var,
                                     const double initial_horizontal_interface_position,
                                     std::vector<std::pair<double, IBTK::VectorNd> > bubbles_position)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_ls_var(ls_var),
      d_initial_horizontal_interface_position(initial_horizontal_interface_position),
      d_bubbles_position(bubbles_position)
{
    // intentionally left blank
    return;
} // LSLocateInterface

LSLocateInterface::~LSLocateInterface()
{
    // intentionally left blank
    return;
}

void
LSLocateInterface::setLevelSetPatchData(int D_idx,
                                        Pointer<HierarchyMathOps> hier_math_ops,
                                        double /*time*/,
                                        bool initial_time)
{
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // If not the intial time, set the level set to the current value maintained
    // by the integrator
    if (!initial_time)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int ls_current_idx =
            var_db->mapVariableAndContextToIndex(d_ls_var, d_adv_diff_solver->getCurrentContext());
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

        hier_cc_data_ops.copyData(D_idx, ls_current_idx);

        return;
    }

    // Set the initial condition for locating the interface

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(D_idx);
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

                (*D_data)(ci) = std::min({ left_bubble_signed_distance,
                                           center_bubble_signed_distance,
                                           right_bubble_signed_distance,
                                           horizontal_interface_signed_distance });
            }
        }
    }
    return;
} // setLevelSetPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
