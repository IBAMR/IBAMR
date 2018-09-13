// Filename LSLocateLayerInterface.cpp
// Created on Jul 5, 2018 by Nishant Nangia

#include "LSLocateLayerInterface.h"

#include <CartesianGridGeometry.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/HierarchyMathOps.h>

/////////////////////////////// STATIC ///////////////////////////////////////

void
callLSLocateLayerInterfaceCallbackFunction(int D_idx,
                                           Pointer<HierarchyMathOps> hier_math_ops,
                                           double time,
                                           bool initial_time,
                                           void* ctx)
{
    // Set the level set information
    static LSLocateLayerInterface* ptr_LSLocateLayerInterface = static_cast<LSLocateLayerInterface*>(ctx);
    ptr_LSLocateLayerInterface->setLevelSetPatchData(D_idx, hier_math_ops, time, initial_time);

    return;
} // callLSLocateLayerInterfaceCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////
LSLocateLayerInterface::LSLocateLayerInterface(const std::string& object_name,
                                               Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                               Pointer<CellVariable<NDIM, double> > ls_var,
                                               LayerInterface init_layer)
    : d_object_name(object_name), d_adv_diff_solver(adv_diff_solver), d_ls_var(ls_var), d_init_layer(init_layer)
{
    // intentionally left blank
    return;
} // LSLocateLayerInterface

LSLocateLayerInterface::~LSLocateLayerInterface()
{
    // intentionally left blank
    return;
}

void
LSLocateLayerInterface::setLevelSetPatchData(int D_idx,
                                             Pointer<HierarchyMathOps> hier_math_ops,
                                             double /*time*/,
                                             bool initial_time)
{
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // If not the intial time, set the level set to the current value maintained by the integrator
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
    const double height = d_init_layer.height;

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
                Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* patch_X_lower = patch_geom->getXLower();
                const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
                const double* const patch_dx = patch_geom->getDx();
                const int d = NDIM - 1;
                double coord = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);

                (*D_data)(ci) = coord - height;
            }
        }
    }
    return;
} // setLevelSetPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
