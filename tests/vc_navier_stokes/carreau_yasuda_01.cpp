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

#include <ibamr/INSVCStaggeredCarreauYasudaHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

// Advance a Carreau-Yasuda fluid and report the final velocity norms and viscosity range.
int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Logger::getInstance()->setAbortAppender(new TestAppender());

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    // The integrator borrows the boundary condition objects.
    std::vector<std::unique_ptr<RobinBcCoefStrategy<NDIM>>> u_bc_coefs_owner(NDIM);
    Pointer<INSVCStaggeredCarreauYasudaHierarchyIntegrator> time_integrator =
        new INSVCStaggeredCarreauYasudaHierarchyIntegrator(
            "INSVCStaggeredCarreauYasudaHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSVCStaggeredCarreauYasudaHierarchyIntegrator"));
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", time_integrator, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> load_balancer =
        new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
        new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

    Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
        "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
    time_integrator->registerVelocityInitialConditions(u_init);

    // Directions without periodicity use the velocity boundary conditions in the input file.
    std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
    if (grid_geometry->getPeriodicShift().min() == 0)
    {
        for (int d = 0; d < NDIM; ++d)
        {
            u_bc_coefs_owner[d] = std::make_unique<muParserRobinBcCoefs>(
                "u_bc_coefs_" + std::to_string(d),
                app_initializer->getComponentDatabase("VelocityBcCoefs_" + std::to_string(d)),
                grid_geometry);
            u_bc_coefs[d] = u_bc_coefs_owner[d].get();
        }
        time_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
    }

    time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
    app_initializer.setNull();

    while (!IBTK::rel_equal_eps(time_integrator->getIntegratorTime(), time_integrator->getEndTime()) &&
           time_integrator->stepsRemaining())
    {
        time_integrator->advanceHierarchy(time_integrator->getMaximumTimeStepSize());
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_idx = var_db->mapVariableAndContextToIndex(time_integrator->getVelocityVariable(),
                                                           time_integrator->getCurrentContext());
    const int mu_idx = var_db->mapVariableAndContextToIndex(time_integrator->getViscosityVariable(),
                                                            time_integrator->getCurrentContext());
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, 0, finest_ln);
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, 0, finest_ln);
    const int wgt_sc_idx = time_integrator->getHierarchyMathOps()->getSideWeightPatchDescriptorIndex();

    plog << std::setprecision(10) << "time = " << time_integrator->getIntegratorTime() << '\n'
         << "u L2-norm = " << hier_sc_data_ops.L2Norm(u_idx, wgt_sc_idx) << '\n'
         << "u max-norm = " << hier_sc_data_ops.maxNorm(u_idx, wgt_sc_idx) << '\n'
         << "mu min = " << hier_cc_data_ops.min(mu_idx) << '\n'
         << "mu max = " << hier_cc_data_ops.max(mu_idx) << '\n';
    return 0;
}
