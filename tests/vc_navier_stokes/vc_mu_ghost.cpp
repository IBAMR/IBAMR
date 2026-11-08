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

#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <EdgeData.h>
#include <EdgeGeometry.h>
#include <LoadBalancer.h>
#include <NodeData.h>
#include <NodeGeometry.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
// Set the viscosity at the requested time. The integrator sets the viscosity at the new time only with such functions.
void
set_viscosity(int mu_idx,
              Pointer<hier::Variable<NDIM>> mu_var,
              Pointer<HierarchyMathOps> hier_math_ops,
              const int /*cycle_num*/,
              const double time,
              const double /*current_time*/,
              const double /*new_time*/,
              void* ctx)
{
    static_cast<CartGridFunction*>(ctx)->setDataOnPatchHierarchy(
        mu_idx, mu_var, hier_math_ops->getPatchHierarchy(), time);
}
} // namespace

// Check that the viscosity set at the start of a time step is positive in the ghost cells and at the ghost nodes (2D)
// or edges (3D), if the boundary data are far below the viscosity or if no viscosity boundary conditions are
// registered. Only values that the integrator sets are checked.
int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<INSVCStaggeredHierarchyIntegrator> ins_integrator =
            new INSVCStaggeredNonConservativeHierarchyIntegrator(
                "INSVCStaggeredNonConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredNonConservativeHierarchyIntegrator"));
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               ins_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            u_bc_coefs[d] =
                new muParserRobinBcCoefs("u_bc_coefs_" + std::to_string(d),
                                         app_initializer->getComponentDatabase("VelocityBcCoefs_" + std::to_string(d)),
                                         grid_geometry);
        }
        ins_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);

        Pointer<CellVariable<NDIM, double>> mu_var = new CellVariable<NDIM, double>("mu");
        Pointer<CartGridFunction> mu_fcn = new muParserCartGridFunction(
            "mu_fcn", app_initializer->getComponentDatabase("ViscosityFunction"), grid_geometry);
        ins_integrator->registerViscosityVariable(mu_var);
        ins_integrator->registerViscosityInitialConditions(mu_fcn);
        ins_integrator->registerResetFluidViscosityFcn(&set_viscosity, static_cast<void*>(mu_fcn.getPointer()));
        std::unique_ptr<muParserRobinBcCoefs> mu_bc_coef;
        if (app_initializer->getInputDatabase()->getBoolWithDefault("register_viscosity_boundary_conditions", true))
        {
            mu_bc_coef = std::make_unique<muParserRobinBcCoefs>(
                "mu_bc_coef", app_initializer->getComponentDatabase("ViscosityBoundaryConditions"), grid_geometry);
            ins_integrator->registerViscosityBoundaryConditions(mu_bc_coef.get());
        }

        ins_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        const double current_time = ins_integrator->getIntegratorTime();
        ins_integrator->preprocessIntegrateHierarchy(
            current_time, current_time + ins_integrator->getMaximumTimeStepSize(), 1);

        // Find the smallest value, or NaN if any value is NaN. A value that is not positive is a failure.
        const double nan = std::numeric_limits<double>::quiet_NaN();
        double min_cell = std::numeric_limits<double>::max(), min_node = std::numeric_limits<double>::max();
        const auto update_min = [&](double& minimum, const double value)
        { minimum = std::isnan(value) ? nan : (std::isnan(minimum) ? minimum : std::min(minimum, value)); };
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);
        const Box<NDIM> domain_box = patch_hierarchy->getGridGeometry()->getPhysicalDomain()[0];
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> cell_data =
                patch->getPatchData(ins_integrator->getLinearOperatorMuPatchDataIndex());
            for (Box<NDIM>::Iterator b(cell_data->getGhostBox()); b; b++)
            {
                update_min(min_cell, (*cell_data)(CellIndex<NDIM>(b())));
            }
            const Box<NDIM> patch_ghost_box = Box<NDIM>::grow(patch->getBox(), 1);
#if (NDIM == 2)
            Pointer<NodeData<NDIM, double>> mu_data =
                patch->getPatchData(ins_integrator->getInterpolatedLinearOperatorMuPatchDataIndex());
            const Box<NDIM> set_box = NodeGeometry<NDIM>::toNodeBox(patch_ghost_box);
            const Box<NDIM> domain_set_box = NodeGeometry<NDIM>::toNodeBox(domain_box);
            for (Box<NDIM>::Iterator b(mu_data->getGhostBox()); b; b++)
            {
                if (set_box.contains(b()) || domain_set_box.contains(b()))
                {
                    update_min(min_node, (*mu_data)(NodeIndex<NDIM>(b(), NodeIndex<NDIM>::LowerLeft)));
                }
            }
#elif (NDIM == 3)
            Pointer<EdgeData<NDIM, double>> mu_data =
                patch->getPatchData(ins_integrator->getInterpolatedLinearOperatorMuPatchDataIndex());
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> set_box = EdgeGeometry<NDIM>::toEdgeBox(patch_ghost_box, axis);
                const Box<NDIM> domain_set_box = EdgeGeometry<NDIM>::toEdgeBox(domain_box, axis);
                const ArrayData<NDIM, double>& array = mu_data->getArrayData(axis);
                for (Box<NDIM>::Iterator b(array.getBox()); b; b++)
                {
                    if (set_box.contains(b()) || domain_set_box.contains(b()))
                    {
                        update_min(min_node, array(b(), 0));
                    }
                }
            }
#endif
        }
#if (NDIM == 2)
        const char* const node_name = "nodes";
#elif (NDIM == 3)
        const char* const node_name = "edges";
#endif
        pout << std::setprecision(6) << "smallest viscosity in the cells, including ghost cells = " << min_cell << '\n'
             << "smallest viscosity at the " << node_name << ", including ghost " << node_name << " = " << min_node
             << '\n';
        if (!(min_cell > 0.0) || !(min_node > 0.0))
        {
            TBOX_ERROR("the viscosity is not positive everywhere\n");
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }
    }
} // main
