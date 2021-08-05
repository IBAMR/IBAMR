// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2019 by the IBAMR developers
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

#include "LSLocateStructureInterface.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

void
callLSLocateStructureInterfaceCallbackFunction(int D_idx,
                                               Pointer<HierarchyMathOps> hier_math_ops,
                                               double time,
                                               bool initial_time,
                                               void* ctx)
{
    // Set the level set information
    static LSLocateStructureInterface* ptr_LSLocateStructureInterface = static_cast<LSLocateStructureInterface*>(ctx);
    ptr_LSLocateStructureInterface->setLevelSetPatchData(D_idx, hier_math_ops, time, initial_time);

    return;
} // callLSLocateStructureInterfaceCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////
LSLocateStructureInterface::LSLocateStructureInterface(const std::string& object_name,
                                                       Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                                       Pointer<CellVariable<NDIM, double> > ls_var,
                                                       LDataManager* lag_data_manager,
                                                       double vol_elem,
                                                       RectangleInterface* rectangle)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_ls_var(ls_var),
      d_lag_data_manager(lag_data_manager),
      d_vol_elem(vol_elem),
      d_rectangle(rectangle)
{
    // intentionally left blank
    return;
} // LSLocateStructureInterface

LSLocateStructureInterface::~LSLocateStructureInterface()
{
    // intentionally left blank
    return;
}

void
LSLocateStructureInterface::setLevelSetPatchData(int D_idx,
                                                 Pointer<HierarchyMathOps> hier_math_ops,
                                                 double time,
                                                 bool initial_time)
{
    setLevelSetPatchDataByGeometry(D_idx, hier_math_ops, time, initial_time);
    return;
} // setLevelSetPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
LSLocateStructureInterface::setLevelSetPatchDataByGeometry(int D_idx,
                                                           Pointer<HierarchyMathOps> hier_math_ops,
                                                           double /*time*/,
                                                           bool /*initial_time*/)
{
    // In this version of this class, the initial level set location is set to be
    // exact. Also we are assuming the rectangle does not move.

    if (NDIM == 2)
    {
        TBOX_ERROR("Presently not implemented for NDIM != 3");
    }

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    const double x0 = d_rectangle->X0(0);
    const double y0 = d_rectangle->X0(1);
    const double z0 = d_rectangle->X0(2);
    const double l = d_rectangle->length;
    const double w = d_rectangle->width;
    const double h = d_rectangle->height;

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
                IBTK::Vector X = IBTK::Vector::Zero();
                Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* patch_X_lower = patch_geom->getXLower();
                const SAMRAI::hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
                const double* const patch_dx = patch_geom->getDx();
                for (int d = 0; d < NDIM; ++d)
                {
                    X[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }

                const double relx = X[0] - x0;
                const double rely = X[1] - y0;
                const double relz = X[2] - z0;

                bool inside = (-0.5 * l < relx && relx < 0.5 * l && -0.5 * w < rely && rely < 0.5 * w &&
                               -0.5 * h < relz && relz < 0.5 * h);

                if (!inside)
                {
                    const double dx = std::max(abs(relx) - 0.5 * l, 0.0);
                    const double dy = std::max(abs(rely) - 0.5 * w, 0.0);
                    const double dz = std::max(abs(relz) - 0.5 * h, 0.0);
                    (*D_data)(ci) = std::sqrt(dx * dx + dy * dy + dz * dz);
                }
                else
                {
                    double dx_min = std::min(abs(relx + 0.5 * l), abs(0.5 * l - relx));
                    double dy_min = std::min(abs(rely + 0.5 * w), abs(0.5 * w - rely));
                    double dz_min = std::min(abs(relz + 0.5 * h), abs(0.5 * h - relz));
                    (*D_data)(ci) = -std::min(dx_min, std::min(dy_min, dz_min));
                }
            }
        }
    }
    return;
} // setLevelSetPatchDataByGeometry
