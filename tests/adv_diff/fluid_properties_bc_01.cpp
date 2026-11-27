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

// Config files

#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <HierarchySideDataOpsReal.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/vc_ins_utilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <memory>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "fluid_properties_bc.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create the grid and the integrator that owns the level set and its boundary condition.
        Pointer<AdvDiffHierarchyIntegrator> time_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
            "AdvDiffSemiImplicitHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
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

        // The level set is constant, so only its ghost cells differ between the two boundary conditions.
        Pointer<CellVariable<NDIM, double>> phi_var = new CellVariable<NDIM, double>("phi");
        time_integrator->registerTransportedQuantity(phi_var);
        time_integrator->setInitialConditions(
            phi_var,
            new muParserCartGridFunction(
                "phi_init", app_initializer->getComponentDatabase("LevelSetInitialCondition"), grid_geometry));
        auto bc_coef_a = std::make_unique<muParserRobinBcCoefs>(
            "bc_coef_a", app_initializer->getComponentDatabase("LevelSetBcCoefsA"), grid_geometry);
        auto bc_coef_b = std::make_unique<muParserRobinBcCoefs>(
            "bc_coef_b", app_initializer->getComponentDatabase("LevelSetBcCoefsB"), grid_geometry);
        time_integrator->setPhysicalBcCoef(phi_var, bc_coef_a.get());

        Pointer<FaceVariable<NDIM, double>> u_adv_var = new FaceVariable<NDIM, double>("u_adv");
        time_integrator->registerAdvectionVelocity(u_adv_var);
        time_integrator->setAdvectionVelocityFunction(
            u_adv_var,
            new muParserCartGridFunction(
                "u_fcn", app_initializer->getComponentDatabase("AdvectionVelocityFunction"), grid_geometry));
        time_integrator->setAdvectionVelocity(phi_var, u_adv_var);

        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Set the side-centered density, which fills the ghost cells of the level set with its boundary condition.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<SideVariable<NDIM, double>> rho_var = new SideVariable<NDIM, double>("rho");
        const int rho_idx = var_db->registerVariableAndContext(rho_var, var_db->getContext("rho"), IntVector<NDIM>(0));
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        const double time = time_integrator->getIntegratorTime();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(rho_idx, time);
        }
        Pointer<HierarchyMathOps> hier_math_ops = new HierarchyMathOps("HierarchyMathOps", patch_hierarchy);
        VCINSUtilities::SetFluidProperties fluid_properties("SetFluidProperties",
                                                            time_integrator,
                                                            phi_var,
                                                            input_db->getDouble("RHO_LIQUID"),
                                                            input_db->getDouble("RHO_GAS"),
                                                            input_db->getDouble("MU_LIQUID"),
                                                            input_db->getDouble("MU_GAS"),
                                                            input_db->getDouble("NUM_INTERFACE_CELLS"));
        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

        plog << std::setprecision(10);
        fluid_properties.setDensityPatchData(rho_idx, rho_var, hier_math_ops, 0, time, time, time);
        plog << "rho with the first boundary condition:\n"
             << "  min = " << hier_sc_data_ops.min(rho_idx) << "\n"
             << "  max = " << hier_sc_data_ops.max(rho_idx) << "\n";

        // Replace the boundary condition on the unchanged hierarchy, and destroy the first one.
        time_integrator->setPhysicalBcCoef(phi_var, bc_coef_b.get());
        bc_coef_a.reset();
        fluid_properties.setDensityPatchData(rho_idx, rho_var, hier_math_ops, 0, time, time, time);
        plog << "rho with the replacement boundary condition:\n"
             << "  min = " << hier_sc_data_ops.min(rho_idx) << "\n"
             << "  max = " << hier_sc_data_ops.max(rho_idx) << "\n";

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(rho_idx);
        }
        var_db->removePatchDataIndex(rho_idx);

    } // cleanup dynamically allocated objects prior to shutdown
} // main
