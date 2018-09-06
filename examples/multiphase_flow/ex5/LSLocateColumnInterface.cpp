// Filename LSLocateColumnInterface.cpp
// Created on Jul 5, 2018 by Nishant Nangia

#include "LSLocateColumnInterface.h"

#include <CartesianGridGeometry.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/HierarchyMathOps.h>

/////////////////////////////// STATIC ///////////////////////////////////////

void
callLSLocateColumnInterfaceCallbackFunction(int D_idx,
                                            Pointer<HierarchyMathOps> hier_math_ops,
                                            double time,
                                            bool initial_time,
                                            void* ctx)
{
    // Set the level set information
    static LSLocateColumnInterface* ptr_LSLocateColumnInterface = static_cast<LSLocateColumnInterface*>(ctx);
    ptr_LSLocateColumnInterface->setLevelSetPatchData(D_idx, hier_math_ops, time, initial_time);

    return;
} // callLSLocateColumnInterfaceCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////
LSLocateColumnInterface::LSLocateColumnInterface(const std::string& object_name,
                                                 Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                                 Pointer<CellVariable<NDIM, double> > ls_var,
                                                 ColumnInterface init_column)
    : d_object_name(object_name), d_adv_diff_solver(adv_diff_solver), d_ls_var(ls_var), d_init_column(init_column)
{
    // intentionally left blank
    return;
} // LSLocateColumnInterface

LSLocateColumnInterface::~LSLocateColumnInterface()
{
    // intentionally left blank
    return;
}

void
LSLocateColumnInterface::setLevelSetPatchData(int D_idx,
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
    const IBTK::Vector& X_UR = d_init_column.X_UR;

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

                // Check if the coordinate is inside the interface
                const bool inside_interface = (coord[0] <= X_UR[0]) && (coord[1] <= X_UR[1])
#if (NDIM == 3)
                                              && (coord[2] <= X_UR[2])
#endif
                    ;
                if (inside_interface)
                {
                    // If inside the interface, simply set the distance to be the minimum distance from all faces of the
                    // column
                    double abs_dist[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        abs_dist[d] = std::abs(coord[d] - X_UR[d]);
                    }
                    (*D_data)(ci) = -(*std::min_element(abs_dist, abs_dist + NDIM));
                }
                else
                {
                    // If outside the interface, figure out the closest face and figure out the distance from that.
                    // Note that this will make a slight error in distances near the corner, but likely does not matter.
                    if (coord[0] >= X_UR[0])
                        (*D_data)(ci) = std::abs(coord[0] - X_UR[0]);
                    else if (coord[1] >= X_UR[1])
                        (*D_data)(ci) = std::abs(coord[1] - X_UR[1]);
#if (NDIM == 3)
                    else if (coord[2] >= X_UR[2])
                        (*D_data)(ci) = std::abs(coord[2] - X_UR[2]);
#endif
                    else
                        TBOX_ERROR("This statement should not be reached");
                }
            }
        }
    }
    return;
} // setLevelSetPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
