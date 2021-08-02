// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include "ibamr/CFINSForcing.h"
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <fstream>
#include <iostream>

#include <ibamr/app_namespaces.h>

#include <ibamr/namespaces.h>

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown
        TimerManager::createManager(nullptr);
        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> time_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator;
        adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
            "AdvDiffSemiImplicitHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        time_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
        // Create initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        time_integrator->registerVelocityInitialConditions(u_init);
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        time_integrator->registerPressureInitialConditions(p_init);

        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        // Create body force function specification objects (when necessary).
        Pointer<CartGridFunctionSet> external_functions = new CartGridFunctionSet("ExternalFunctions");
        Pointer<CFINSForcing> polymericStressForcing;
        Pointer<CartGridFunction> f_fcn;
        if (input_db->keyExists("ForcingFunction"))
        {
            f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            external_functions->addFunction(f_fcn);
        }
        if (input_db->keyExists("ComplexFluid"))
        {
            polymericStressForcing = new CFINSForcing("PolymericStressForcing",
                                                      app_initializer->getComponentDatabase("ComplexFluid"),
                                                      time_integrator,
                                                      grid_geometry,
                                                      adv_diff_integrator,
                                                      visit_data_writer);

            external_functions->addFunction(polymericStressForcing);
        }
        time_integrator->registerBodyForceFunction(external_functions);
        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            iteration_num += 1;
        }

        Pointer<CellVariable<NDIM, double> > s_var = polymericStressForcing->getVariable();
        Pointer<CellVariable<NDIM, double> > sxx_var = new CellVariable<NDIM, double>("Sxx");
        Pointer<CellVariable<NDIM, double> > syy_var = new CellVariable<NDIM, double>("Syy");
        Pointer<CellVariable<NDIM, double> > sxy_var = new CellVariable<NDIM, double>("Sxy");
        const Pointer<VariableContext> s_ctx = adv_diff_integrator->getCurrentContext();

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int s_idx = var_db->mapVariableAndContextToIndex(s_var, s_ctx);
        const int sxx_idx = var_db->registerVariableAndContext(sxx_var, s_ctx);
        const int syy_idx = var_db->registerVariableAndContext(syy_var, s_ctx);
        const int sxy_idx = var_db->registerVariableAndContext(sxy_var, s_ctx);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(sxx_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(syy_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(sxy_idx, loop_time);
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM, double> > s_data = patch->getPatchData(s_idx);
                Pointer<CellData<NDIM, double> > sxx_data = patch->getPatchData(sxx_idx);
                Pointer<CellData<NDIM, double> > syy_data = patch->getPatchData(syy_idx);
                Pointer<CellData<NDIM, double> > sxy_data = patch->getPatchData(sxy_idx);

                sxx_data->copyDepth(0, *s_data, 0);
                syy_data->copyDepth(0, *s_data, 1);
                sxy_data->copyDepth(0, *s_data, 2);
            }
        }

        HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

        pout << "Soln norms in sxx at time " << loop_time << ":\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(sxx_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << std::setprecision(10) << hier_cc_data_ops.L2Norm(sxx_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << std::setprecision(10) << hier_cc_data_ops.maxNorm(sxx_idx, wgt_cc_idx) << "\n";

        pout << "Soln norms in syy at time " << loop_time << ":\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(syy_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << std::setprecision(10) << hier_cc_data_ops.L2Norm(syy_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << std::setprecision(10) << hier_cc_data_ops.maxNorm(syy_idx, wgt_cc_idx) << "\n";

        pout << "Soln norms in sxy at time " << loop_time << ":\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(sxy_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << std::setprecision(10) << hier_cc_data_ops.L2Norm(sxy_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << std::setprecision(10) << hier_cc_data_ops.maxNorm(sxy_idx, wgt_cc_idx) << "\n";

    } // cleanup dynamically allocated objects prior to shutdown
} // main
