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
#include <ibamr/AdvDiffPredictorCorrectorHierarchyIntegrator.h>
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include <LocationIndexRobinBcCoefs.h>

// Set up application namespace declarations
#include "ibamr/AdvDiffConvectiveOperatorManager.h"

#include "ibtk/muParserCartGridFunction.h"
#include "ibtk/muParserRobinBcCoefs.h"

#include <SAMRAIVectorReal.h>

#include <ibamr/app_namespaces.h>

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        // prevent a warning about timer initialization
        TimerManager::createManager(nullptr);

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "adv_diff.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");

        // Create major algorithm and data objects that comprise the
        // application.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        Pointer<VisItDataWriter<NDIM> > visit_writer = app_initializer->getVisItDataWriter();

        auto var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> var_ctx = var_db->getContext("Context");
        Pointer<CellVariable<NDIM, double> > q_var = new CellVariable<NDIM, double>("CC_var");
        Pointer<CellVariable<NDIM, double> > convec_var = new CellVariable<NDIM, double>("Convec var");
        Pointer<CellVariable<NDIM, double> > exact_var = new CellVariable<NDIM, double>("Exact var");
        Pointer<FaceVariable<NDIM, double> > u_var = new FaceVariable<NDIM, double>("U");

        const int q_idx = var_db->registerVariableAndContext(q_var, var_ctx);
        const int convec_idx = var_db->registerVariableAndContext(convec_var, var_ctx);
        const int exact_idx = var_db->registerVariableAndContext(exact_var, var_ctx);
        const int u_idx = var_db->registerVariableAndContext(u_var, var_ctx);
//#define OUTPUT_VIZ_FILES // Comment out if you want to draw things
#ifdef OUTPUT_VIZ_FILES
        visit_writer->registerPlotQuantity("Q", "SCALAR", q_idx);
        visit_writer->registerPlotQuantity("Convec", "SCALAR", convec_idx);
        visit_writer->registerPlotQuantity("Error", "SCALAR", exact_idx);
#endif

        Pointer<muParserCartGridFunction> u_fcn =
            new muParserCartGridFunction("U", app_initializer->getComponentDatabase("U"), grid_geometry);
        Pointer<muParserCartGridFunction> q_fcn =
            new muParserCartGridFunction("Q", app_initializer->getComponentDatabase("Q"), grid_geometry);
        Pointer<muParserCartGridFunction> exact_fcn =
            new muParserCartGridFunction("Exact", app_initializer->getComponentDatabase("Exact"), grid_geometry);

        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        std::vector<RobinBcCoefStrategy<NDIM>*> q_bc_coefs(1);
        if (periodic_shift.min() > 0)
            q_bc_coefs[0] = nullptr;
        else
            q_bc_coefs[0] =
                new muParserRobinBcCoefs("Q_bcs", app_initializer->getComponentDatabase("Q_bcs"), grid_geometry);

        std::vector<std::string> convec_oper_types = { "CENTERED", "CUI", "PPM", "WAVE_PROP" };
        std::vector<Pointer<ConvectiveOperator> > convec_opers(convec_oper_types.size());
        auto oper_manager = AdvDiffConvectiveOperatorManager::getManager();
        int i = 0;
        for (const auto& convec_oper_type : convec_oper_types)
        {
            auto differencing_type =
                IBAMR::string_to_enum<ConvectiveDifferencingType>(input_db->getString("CONVECTIVE_DIFFERENCING_TYPE"));
            convec_opers[i++] = oper_manager->allocateOperator(convec_oper_type,
                                                               convec_oper_type,
                                                               q_var,
                                                               app_initializer->getComponentDatabase("ConvecOper"),
                                                               differencing_type,
                                                               q_bc_coefs);
        }

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        const int tag_buffer = std::numeric_limits<int>::max();
        int level_number = 0;
        while (gridding_algorithm->levelCanBeRefined(level_number))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            ++level_number;
        }

        const int finest_level = patch_hierarchy->getFinestLevelNumber();
        for (int ln = 0; ln <= finest_level; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(q_idx, 0.0);
            level->allocatePatchData(convec_idx, 0.0);
            level->allocatePatchData(exact_idx, 0.0);
            level->allocatePatchData(u_idx, 0.0);
        }

        solv::SAMRAIVectorReal<NDIM, double> q_vec("Q_vec", patch_hierarchy, 0, finest_level);
        q_vec.addComponent(q_var, q_idx);

        u_fcn->setDataOnPatchHierarchy(u_idx, u_var, patch_hierarchy, 0.0, false, 0, finest_level);
        q_fcn->setDataOnPatchHierarchy(q_idx, q_var, patch_hierarchy, 0.0, false, 0, finest_level);

#ifdef OUTPUT_VIZ_FILES
        int step = 0;
#endif
        auto do_test = [&](Pointer<ConvectiveOperator> convec_oper)
        {
            convec_oper->initializeOperatorState(q_vec, q_vec);
            convec_oper->setAdvectionVelocity(u_idx);
            convec_oper->applyConvectiveOperator(q_idx, convec_idx);
            exact_fcn->setDataOnPatchHierarchy(exact_idx, exact_var, patch_hierarchy, 0.0, false, 0, finest_level);

            HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
            hier_math_ops.resetLevels(0, finest_level);
            const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

            HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, 0, finest_level);
            hier_cc_data_ops.subtract(exact_idx, convec_idx, exact_idx);

            pout << "Error using: " << convec_oper->getName() << ":\n"
                 << "  L1-norm :" << hier_cc_data_ops.L1Norm(exact_idx, wgt_cc_idx) << "\n"
                 << "  L2-norm :" << hier_cc_data_ops.L2Norm(exact_idx, wgt_cc_idx) << "\n"
                 << "  max-norm:" << hier_cc_data_ops.maxNorm(exact_idx, wgt_cc_idx) << "\n";
#ifdef OUTPUT_VIZ_FILES
            visit_writer->writePlotData(patch_hierarchy, step++, 0.0);
#endif
        };

        for (const auto& convec_oper : convec_opers) do_test(convec_oper);

        for (int ln = 0; ln <= finest_level; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(q_idx);
            level->deallocatePatchData(convec_idx);
            level->deallocatePatchData(exact_idx);
            level->deallocatePatchData(u_idx);
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main
