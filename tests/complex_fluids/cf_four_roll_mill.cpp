// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
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
        auto app_initializer = make_samrai_shared<AppInitializer>(argc, argv, "INS.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        SAMRAIPointer<INSHierarchyIntegrator> time_integrator = make_samrai_shared<INSStaggeredHierarchyIntegrator>(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        SAMRAIPointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator;
        adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
            "AdvDiffSemiImplicitHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        auto grid_geometry = make_samrai_shared<CartesianGridGeometryNd>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        auto patch_hierarchy = make_samrai_shared<PatchHierarchyNd>("PatchHierarchy", grid_geometry);
        auto error_detector = make_samrai_shared<StandardTagAndInitializeNd>(
            "StandardTagAndInitialize",
            time_integrator,
            app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        auto box_generator = make_samrai_shared<BergerRigoutsosNd>();
        auto load_balancer =
            make_samrai_shared<LoadBalancerNd>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        auto gridding_algorithm =
            make_samrai_shared<GriddingAlgorithmNd>("GriddingAlgorithm",
                                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                                    error_detector,
                                                    box_generator,
                                                    load_balancer);
        time_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
        // Create initial condition specification objects.
        SAMRAIPointer<CartGridFunction> u_init = make_samrai_shared<muParserCartGridFunction>(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        time_integrator->registerVelocityInitialConditions(u_init);
        SAMRAIPointer<CartGridFunction> p_init = make_samrai_shared<muParserCartGridFunction>(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        time_integrator->registerPressureInitialConditions(p_init);

        SAMRAIPointer<VisItDataWriterNd> visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        // Create body force function specification objects (when necessary).
        auto external_functions = make_samrai_shared<CartGridFunctionSet>("ExternalFunctions");
        SAMRAIPointer<CFINSForcing> polymericStressForcing;
        SAMRAIPointer<CartGridFunction> f_fcn;
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
        double loop_time = time_integrator->getIntegratorTime();

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            loop_time = time_integrator->getIntegratorTime();

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;
        }

        SAMRAIPointer<CellVariableNd<double> > s_var = polymericStressForcing->getVariable();
        SAMRAIPointer<CellVariableNd<double> > sxx_var = make_samrai_shared<CellVariableNd<double> >("Sxx");
        SAMRAIPointer<CellVariableNd<double> > syy_var = make_samrai_shared<CellVariableNd<double> >("Syy");
        SAMRAIPointer<CellVariableNd<double> > sxy_var = make_samrai_shared<CellVariableNd<double> >("Sxy");
        const SAMRAIPointer<VariableContext> s_ctx = adv_diff_integrator->getCurrentContext();

        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
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
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                SAMRAIPointer<PatchNd> patch = level->getPatch(p());
                SAMRAIPointer<CellDataNd<double> > s_data = patch->getPatchData(s_idx);
                SAMRAIPointer<CellDataNd<double> > sxx_data = patch->getPatchData(sxx_idx);
                SAMRAIPointer<CellDataNd<double> > syy_data = patch->getPatchData(syy_idx);
                SAMRAIPointer<CellDataNd<double> > sxy_data = patch->getPatchData(sxy_idx);

                sxx_data->copyDepth(0, *s_data, 0);
                syy_data->copyDepth(0, *s_data, 1);
                sxy_data->copyDepth(0, *s_data, 2);
            }
        }

        HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        HierarchyCellDataOpsRealNd<double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

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
