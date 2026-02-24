// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// GENERAL CONFIGURATION

// SAMRAI INCLUDES
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIBox.h"
#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAICartesianPatchGeometry.h"
#include "SAMRAICellData.h"
#include "SAMRAICellIndex.h"
#include "SAMRAICellIterator.h"
#include "SAMRAICellVariable.h"
#include "SAMRAICoarseFineBoundary.h"
#include "SAMRAICoarsenAlgorithm.h"
#include "SAMRAICoarsenOperator.h"
#include "SAMRAIHierarchyCellDataOpsReal.h"
#include "SAMRAIHierarchySideDataOpsReal.h"
#include "SAMRAIIndex.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchGeometry.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIRefineAlgorithm.h"
#include "SAMRAIRefineOperator.h"
#include "SAMRAISideData.h"
#include "SAMRAISideIndex.h"
#include "SAMRAISideIterator.h"
#include "SAMRAISideVariable.h"
#include "SAMRAIVariableDatabase.h"
#include "SAMRAIVisItDataWriter.h"

#include <petscsys.h>

#include <SAMRAI_config.h>

// IBAMR INCLUDES
#include <ibamr/app_namespaces.h>

// IBTK INCLUDES
#include <ibtk/AppInitializer.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include "BoxArray.h"

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name> [PETSc options]                             *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Parse command line options, set some standard options from the input
    // file, and enable file logging.
    SAMRAIPointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
    SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

    // Retrieve "Main" section of the input database.
    SAMRAIPointer<Database> main_db = app_initializer->getComponentDatabase("Main");

    int coarse_hier_dump_interval = 0;
    int fine_hier_dump_interval = 0;
    if (main_db->keyExists("hier_dump_interval"))
    {
        coarse_hier_dump_interval = main_db->getInteger("hier_dump_interval");
        fine_hier_dump_interval = main_db->getInteger("hier_dump_interval");
    }
    else if (main_db->keyExists("coarse_hier_dump_interval") && main_db->keyExists("fine_hier_dump_interval"))
    {
        coarse_hier_dump_interval = main_db->getInteger("coarse_hier_dump_interval");
        fine_hier_dump_interval = main_db->getInteger("fine_hier_dump_interval");
    }
    else
    {
        TBOX_ERROR("hierarchy dump intervals not specified in input file. . .\n");
    }

    string coarse_hier_dump_dirname;
    if (main_db->keyExists("coarse_hier_dump_dirname"))
    {
        coarse_hier_dump_dirname = main_db->getString("coarse_hier_dump_dirname");
    }
    else
    {
        TBOX_ERROR("key `coarse_hier_dump_dirname' not specified in input file");
    }

    string fine_hier_dump_dirname;
    if (main_db->keyExists("fine_hier_dump_dirname"))
    {
        fine_hier_dump_dirname = main_db->getString("fine_hier_dump_dirname");
    }
    else
    {
        TBOX_ERROR("key `fine_hier_dump_dirname' not specified in input file");
    }

    // Create major algorithm and data objects that comprise application.
    SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geom = new SAMRAICartesianGridGeometry(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

    // Initialize variables.
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();

    SAMRAIPointer<VariableContext> current_ctx = var_db->getContext("INSStaggeredHierarchyIntegrator::CURRENT");
    SAMRAIPointer<VariableContext> scratch_ctx = var_db->getContext("INSStaggeredHierarchyIntegrator::SCRATCH");

    SAMRAIPointer<SAMRAISideVariable<double>> U_var =
        new SAMRAISideVariable<double>("INSStaggeredHierarchyIntegrator::U");
    const int U_idx = var_db->registerVariableAndContext(U_var, current_ctx);
    const int U_interp_idx = var_db->registerClonedPatchDataIndex(U_var, U_idx);
    const int U_scratch_idx = var_db->registerVariableAndContext(U_var, scratch_ctx, 2);

    SAMRAIPointer<SAMRAICellVariable<double>> IDX_var = new SAMRAICellVariable<double>("Indicator");
    const int I_idx = var_db->registerVariableAndContext(IDX_var, current_ctx);

    SAMRAIPointer<SAMRAICellVariable<double>> P_var =
        new SAMRAICellVariable<double>("INSStaggeredHierarchyIntegrator::P");
    const int P_idx = var_db->registerVariableAndContext(P_var, current_ctx);
    const int P_interp_idx = var_db->registerClonedPatchDataIndex(P_var, P_idx);
    const int P_scratch_idx = var_db->registerVariableAndContext(P_var, scratch_ctx, 2);

    SAMRAIPointer<VariableContext> s_c_ctx = var_db->getContext("AdvDiffSemiImplicitHierarchyIntegrator::CURRENT");
    SAMRAIPointer<VariableContext> s_s_ctx = var_db->getContext("AdvDiffSemiImplicitHierarchyIntegrator::SCRATCH");

    SAMRAIPointer<SAMRAICellVariable<double>> S_var =
        new SAMRAICellVariable<double>("ComplexFluidForcing::W_cc", NDIM * (NDIM + 1) / 2);
    const int S_idx = var_db->registerVariableAndContext(S_var, s_c_ctx);
    const int S_interp_idx = var_db->registerClonedPatchDataIndex(S_var, S_idx);
    const int S_scratch_idx = var_db->registerVariableAndContext(S_var, scratch_ctx, 2);

    SAMRAIPointer<SAMRAICellVariable<double>> Sxx_var = new SAMRAICellVariable<double>("Sxx");
    SAMRAIPointer<SAMRAICellVariable<double>> Syy_var = new SAMRAICellVariable<double>("Syy");
    SAMRAIPointer<SAMRAICellVariable<double>> Sxy_var = new SAMRAICellVariable<double>("Sxy");
    SAMRAIPointer<SAMRAICellVariable<double>> U_draw_var = new SAMRAICellVariable<double>("U", NDIM);
    const int Sxx_idx = var_db->registerVariableAndContext(Sxx_var, s_c_ctx);
    const int Syy_idx = var_db->registerVariableAndContext(Syy_var, s_c_ctx);
    const int Sxy_idx = var_db->registerVariableAndContext(Sxy_var, s_c_ctx);
    const int U_draw_idx = var_db->registerVariableAndContext(U_draw_var, s_c_ctx);
    // Set up visualization plot file writer.
    SAMRAIPointer<SAMRAIVisItDataWriter> visit_data_writer =
        new SAMRAIVisItDataWriter("VisIt Writer", main_db->getString("viz_dump_dirname"), 1);
    visit_data_writer->registerPlotQuantity("Indicator", "SCALAR", I_idx);
    visit_data_writer->registerPlotQuantity("P", "SCALAR", P_idx);
    visit_data_writer->registerPlotQuantity("P interp", "SCALAR", P_interp_idx);
    visit_data_writer->registerPlotQuantity("Sxx", "SCALAR", S_idx);
    visit_data_writer->registerPlotQuantity("Sxx interp", "SCALAR", S_interp_idx);
    visit_data_writer->registerPlotQuantity("Syy", "SCALAR", S_idx, 1);
    visit_data_writer->registerPlotQuantity("Syy interp", "SCALAR", S_interp_idx, 1);
    visit_data_writer->registerPlotQuantity("Sxy", "SCALAR", S_idx, 2);
    visit_data_writer->registerPlotQuantity("Sxy interp", "SCALAR", S_interp_idx, 2);
    visit_data_writer->registerPlotQuantity("U", "VECTOR", U_draw_idx);

    // Time step loop.
    double loop_time = 0.0;
    int coarse_iteration_num = coarse_hier_dump_interval;
    int fine_iteration_num = fine_hier_dump_interval;

    bool files_exist = true;
    for (; files_exist;
         coarse_iteration_num += coarse_hier_dump_interval, fine_iteration_num += fine_hier_dump_interval)
    {
        char temp_buf[128];

        std::snprintf(temp_buf, sizeof(temp_buf), "%05d.samrai.%05d", coarse_iteration_num, IBTK_MPI::getRank());
        string coarse_file_name = coarse_hier_dump_dirname + "/" + "hier_data.";
        coarse_file_name += temp_buf;

        std::snprintf(temp_buf, sizeof(temp_buf), "%05d.samrai.%05d", fine_iteration_num, IBTK_MPI::getRank());
        string fine_file_name = fine_hier_dump_dirname + "/" + "hier_data.";
        fine_file_name += temp_buf;

        for (int rank = 0; rank < IBTK_MPI::getNodes(); ++rank)
        {
            if (rank == IBTK_MPI::getRank())
            {
                fstream coarse_fin, fine_fin;
                coarse_fin.open(coarse_file_name.c_str(), ios::in);
                fine_fin.open(fine_file_name.c_str(), ios::in);
                if (!coarse_fin.is_open() || !fine_fin.is_open())
                {
                    std::cout << "couldnt find file " << rank << "\n";
                    std::cout << "File name was: " << coarse_file_name << "\n and " << fine_file_name << "\n";
                    files_exist = false;
                }
                else
                {
                    std::cout << "found file " << rank << "\n";
                }
                coarse_fin.close();
                fine_fin.close();
            }
            IBTK_MPI::barrier();
        }

        if (!files_exist) break;

        pout << endl;
        pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        pout << "processing data" << endl;
        pout << "     coarse iteration number = " << coarse_iteration_num << endl;
        pout << "     fine iteration number = " << fine_iteration_num << endl;
        pout << "     coarse file name = " << coarse_file_name << endl;
        pout << "     fine file name = " << fine_file_name << endl;

        // Read in data to post-process.
        ComponentSelector hier_data;
        hier_data.setFlag(U_idx);
        hier_data.setFlag(P_idx);
        hier_data.setFlag(S_idx);

        SAMRAIPointer<HDFDatabase> coarse_hier_db = new HDFDatabase("coarse_hier_db");
        coarse_hier_db->open(coarse_file_name);

        SAMRAIPointer<SAMRAIPatchHierarchy> coarse_patch_hierarchy =
            new SAMRAIPatchHierarchy("CoarsePatchHierarchy", grid_geom, false);
        coarse_patch_hierarchy->getFromDatabase(coarse_hier_db->getDatabase("PatchHierarchy"), hier_data);

        const double coarse_loop_time = coarse_hier_db->getDouble("loop_time");

        coarse_hier_db->close();

        SAMRAIPointer<HDFDatabase> fine_hier_db = new HDFDatabase("fine_hier_db");
        fine_hier_db->open(fine_file_name);

        SAMRAIPointer<SAMRAIPatchHierarchy> fine_patch_hierarchy = new SAMRAIPatchHierarchy(
            "FinePatchHierarchy", grid_geom->makeRefinedGridGeometry("FineGridGeometry", 2, false), false);
        fine_patch_hierarchy->getFromDatabase(fine_hier_db->getDatabase("PatchHierarchy"), hier_data);

        const double fine_loop_time = fine_hier_db->getDouble("loop_time");

        fine_hier_db->close();

        TBOX_ASSERT(IBTK::rel_equal_eps(coarse_loop_time, fine_loop_time));
        loop_time = fine_loop_time;
        pout << "     loop time = " << loop_time << endl;

        SAMRAIPointer<SAMRAIPatchHierarchy> coarsened_fine_patch_hierarchy =
            fine_patch_hierarchy->makeCoarsenedPatchHierarchy("CoarsenedFinePatchHierarchy", 2, false);

        // Setup hierarchy operations objects.
        SAMRAIHierarchyCellDataOpsReal<double> coarse_hier_cc_data_ops(
            coarse_patch_hierarchy, 0, coarse_patch_hierarchy->getFinestLevelNumber());
        SAMRAIHierarchySideDataOpsReal<double> coarse_hier_sc_data_ops(
            coarse_patch_hierarchy, 0, coarse_patch_hierarchy->getFinestLevelNumber());
        HierarchyMathOps hier_math_ops("hier_math_ops", coarse_patch_hierarchy);
        hier_math_ops.setPatchHierarchy(coarse_patch_hierarchy);
        hier_math_ops.resetLevels(0, coarse_patch_hierarchy->getFinestLevelNumber());
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        // Allocate patch data.
        for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = coarse_patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(U_interp_idx, loop_time);
            level->allocatePatchData(P_interp_idx, loop_time);
            level->allocatePatchData(S_interp_idx, loop_time);
            level->allocatePatchData(U_scratch_idx, loop_time);
            level->allocatePatchData(P_scratch_idx, loop_time);
            level->allocatePatchData(S_scratch_idx, loop_time);
            level->allocatePatchData(Sxx_idx, loop_time);
            level->allocatePatchData(Syy_idx, loop_time);
            level->allocatePatchData(Sxy_idx, loop_time);
            level->allocatePatchData(I_idx, loop_time);
            level->allocatePatchData(U_draw_idx, loop_time);
        }

        for (int ln = 0; ln <= fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = fine_patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(U_interp_idx, loop_time);
            level->allocatePatchData(P_interp_idx, loop_time);
            level->allocatePatchData(S_interp_idx, loop_time);
            level->allocatePatchData(U_scratch_idx, loop_time);
            level->allocatePatchData(P_scratch_idx, loop_time);
            level->allocatePatchData(S_scratch_idx, loop_time);
        }

        for (int ln = 0; ln <= coarsened_fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(U_idx, loop_time);
            level->allocatePatchData(P_idx, loop_time);
            level->allocatePatchData(S_idx, loop_time);
            level->allocatePatchData(U_interp_idx, loop_time);
            level->allocatePatchData(P_interp_idx, loop_time);
            level->allocatePatchData(S_interp_idx, loop_time);
            level->allocatePatchData(U_scratch_idx, loop_time);
            level->allocatePatchData(P_scratch_idx, loop_time);
            level->allocatePatchData(S_scratch_idx, loop_time);
        }

        // Synchronize the coarse hierarchy data.
        for (int ln = coarse_patch_hierarchy->getFinestLevelNumber(); ln > 0; --ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> coarser_level = coarse_patch_hierarchy->getPatchLevel(ln - 1);
            SAMRAIPointer<SAMRAIPatchLevel> finer_level = coarse_patch_hierarchy->getPatchLevel(ln);

            SAMRAICoarsenAlgorithm coarsen_alg;
            SAMRAIPointer<SAMRAICoarsenOperator> coarsen_op;

            coarsen_op = grid_geom->lookupCoarsenOperator(U_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(U_idx, U_idx, coarsen_op);

            coarsen_op = grid_geom->lookupCoarsenOperator(P_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(P_idx, P_idx, coarsen_op);

            coarsen_op = grid_geom->lookupCoarsenOperator(S_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(S_idx, S_idx, coarsen_op);

            coarsen_alg.createSchedule(coarser_level, finer_level)->coarsenData();
        }

        // Synchronize the fine hierarchy data.
        for (int ln = fine_patch_hierarchy->getFinestLevelNumber(); ln > 0; --ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> coarser_level = fine_patch_hierarchy->getPatchLevel(ln - 1);
            SAMRAIPointer<SAMRAIPatchLevel> finer_level = fine_patch_hierarchy->getPatchLevel(ln);

            SAMRAICoarsenAlgorithm coarsen_alg;
            SAMRAIPointer<SAMRAICoarsenOperator> coarsen_op;

            coarsen_op = grid_geom->lookupCoarsenOperator(U_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(U_idx, U_idx, coarsen_op);

            coarsen_op = grid_geom->lookupCoarsenOperator(P_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(P_idx, P_idx, coarsen_op);

            coarsen_op = grid_geom->lookupCoarsenOperator(S_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(S_idx, S_idx, coarsen_op);

            coarsen_alg.createSchedule(coarser_level, finer_level)->coarsenData();
        }

        // Coarsen data from the fine hierarchy to the coarsened fine hierarchy.
        for (int ln = 0; ln <= fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> dst_level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);
            SAMRAIPointer<SAMRAIPatchLevel> src_level = fine_patch_hierarchy->getPatchLevel(ln);

            SAMRAIPointer<SAMRAICoarsenOperator> coarsen_op;
            for (SAMRAIPatchLevel::Iterator p(dst_level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> dst_patch = dst_level->getPatch(p());
                SAMRAIPointer<SAMRAIPatch> src_patch = src_level->getPatch(p());
                const SAMRAIBox& coarse_box = dst_patch->getBox();
                TBOX_ASSERT(SAMRAIBox::coarsen(src_patch->getBox(), 2) == coarse_box);

                coarsen_op = grid_geom->lookupCoarsenOperator(U_var, "CONSERVATIVE_COARSEN");
                coarsen_op->coarsen(*dst_patch, *src_patch, U_interp_idx, U_idx, coarse_box, 2);

                coarsen_op = grid_geom->lookupCoarsenOperator(P_var, "CONSERVATIVE_COARSEN");
                coarsen_op->coarsen(*dst_patch, *src_patch, P_interp_idx, P_idx, coarse_box, 2);

                coarsen_op = grid_geom->lookupCoarsenOperator(S_var, "CONSERVATIVE_COARSEN");
                coarsen_op->coarsen(*dst_patch, *src_patch, S_interp_idx, S_idx, coarse_box, 2);
            }
        }

        // Interpolate and copy data from the coarsened fine patch hierarchy to
        // the coarse patch hierarchy.
        for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            pout << "Interpolating on level " << ln << "\n";
            SAMRAIPointer<SAMRAIPatchLevel> dst_level = coarse_patch_hierarchy->getPatchLevel(ln);
            SAMRAIPointer<SAMRAIPatchLevel> src_level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);

            SAMRAIRefineAlgorithm refine_alg;
            SAMRAIPointer<SAMRAIRefineOperator> refine_op;

            refine_op = grid_geom->lookupRefineOperator(U_var, "CONSERVATIVE_LINEAR_REFINE");
            refine_alg.registerRefine(U_interp_idx, U_interp_idx, U_scratch_idx, refine_op);

            refine_op = grid_geom->lookupRefineOperator(P_var, "LINEAR_REFINE");
            refine_alg.registerRefine(P_interp_idx, P_interp_idx, P_scratch_idx, refine_op);

            refine_op = grid_geom->lookupRefineOperator(S_var, "LINEAR_REFINE");
            refine_alg.registerRefine(S_interp_idx, S_interp_idx, S_scratch_idx, refine_op);

            ComponentSelector data_indices;
            data_indices.setFlag(U_scratch_idx);
            data_indices.setFlag(P_scratch_idx);
            data_indices.setFlag(S_scratch_idx);
            CartExtrapPhysBdryOp bc_helper(data_indices, "LINEAR");

            refine_alg.createSchedule(dst_level, src_level, ln - 1, coarse_patch_hierarchy, &bc_helper)
                ->fillData(loop_time);
        }

        for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = coarse_patch_hierarchy->getPatchLevel(ln);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                const SAMRAIBox& box = patch->getBox();
                SAMRAIPointer<SAMRAICellData<double>> U_draw_data = patch->getPatchData(U_draw_idx);
                SAMRAIPointer<SAMRAISideData<double>> U_data = patch->getPatchData(U_interp_idx);
                for (SAMRAICellIterator i(box); i; i++)
                {
                    SAMRAICellIndex idx = i();
                    SAMRAISideIndex idx_b(idx, 1, 0);
                    SAMRAISideIndex idx_u(idx, 1, 1);
                    SAMRAISideIndex idx_l(idx, 0, 0);
                    SAMRAISideIndex idx_r(idx, 0, 1);
                    (*U_draw_data)(idx, 0) = 0.5 * ((*U_data)(idx_l) + (*U_data)(idx_r));
                    (*U_draw_data)(idx, 1) = 0.5 * ((*U_data)(idx_b) + (*U_data)(idx_u));
                }
            }
        }

        // Output plot data before taking norms of differences.
        visit_data_writer->writePlotData(coarse_patch_hierarchy, coarse_iteration_num, loop_time);

        // Compute norms of differences.
        coarse_hier_sc_data_ops.subtract(U_interp_idx, U_idx, U_interp_idx);
        coarse_hier_cc_data_ops.subtract(P_interp_idx, P_idx, P_interp_idx);
        coarse_hier_cc_data_ops.subtract(S_interp_idx, S_idx, S_interp_idx);

        for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = coarse_patch_hierarchy->getPatchLevel(ln);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                SAMRAIPointer<SAMRAICellData<double>> S_data = patch->getPatchData(S_interp_idx);
                SAMRAIPointer<SAMRAICellData<double>> Sxx_data = patch->getPatchData(Sxx_idx);
                SAMRAIPointer<SAMRAICellData<double>> Syy_data = patch->getPatchData(Syy_idx);
                SAMRAIPointer<SAMRAICellData<double>> Sxy_data = patch->getPatchData(Sxy_idx);
                Sxx_data->copyDepth(0, *S_data, 0);
                Syy_data->copyDepth(0, *S_data, 1);
                Sxy_data->copyDepth(0, *S_data, 2);
            }
        }

        for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = coarse_patch_hierarchy->getPatchLevel(ln);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                const SAMRAIBox& box = patch->getBox();
                const SAMRAIIndex& idx_low = box.lower();
                SAMRAIPointer<SAMRAICartesianPatchGeometry> pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double* const xlow = pgeom->getXLower();
                SAMRAIPointer<SAMRAICellData<double>> wgt_cc_data = patch->getPatchData(wgt_cc_idx);
                SAMRAIPointer<SAMRAICellData<double>> I_data = patch->getPatchData(I_idx);
                SAMRAIPointer<SAMRAISideData<double>> wgt_sc_data = patch->getPatchData(wgt_sc_idx);
                I_data->fillAll(1.0);
                for (SAMRAICellIterator i(box); i; i++)
                {
                    SAMRAICellIndex idx = *i;
                    std::vector<double> x;
                    x.resize(NDIM);
                    for (int d = 0; d < NDIM; ++d) x[d] = xlow[d] + dx[d] * (idx(d) - idx_low(d) + 0.5);
                    double dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1]);
                    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
                    if ((r > (1.0 - 4 * dr)) && (r < (1.0 + 4 * dr)))
                    {
                        pout << "Found a center cell\n";
                        (*wgt_cc_data)(idx) = 0.0;
                        (*I_data)(idx) = 0.0;
                    }
                }
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (SAMRAISideIterator i(box, axis); i; i++)
                    {
                        const SAMRAISideIndex& idx = *i;
                        std::vector<double> x;
                        x.resize(NDIM);
                        for (int d = 0; d < NDIM; ++d)
                            x[d] = xlow[d] + dx[d] * (idx(d) - idx_low(d) + (axis == 0 ? 0.0 : 0.5));
                        double dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1]);
                        double r = sqrt(x[0] + x[0] + x[1] * x[1]);
                        if ((r > (1.0 - 4 * dr)) && (r < (1.0 + 4 * dr)))
                        {
                            pout << "Found a side cell\n";
                            (*wgt_sc_data)(idx) = 0.0;
                        }
                    }
                }
            }
        }

        for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = coarse_patch_hierarchy->getPatchLevel(ln);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                const SAMRAIBox& box = patch->getBox();
                SAMRAIPointer<SAMRAICellData<double>> U_draw_data = patch->getPatchData(U_draw_idx);
                SAMRAIPointer<SAMRAISideData<double>> U_data = patch->getPatchData(U_interp_idx);
                for (SAMRAICellIterator i(box); i; i++)
                {
                    SAMRAICellIndex idx = i();
                    SAMRAISideIndex idx_b(idx, 1, 0);
                    SAMRAISideIndex idx_u(idx, 1, 1);
                    SAMRAISideIndex idx_l(idx, 0, 0);
                    SAMRAISideIndex idx_r(idx, 0, 1);
                    (*U_draw_data)(idx, 0) = 0.5 * ((*U_data)(idx_l) + (*U_data)(idx_r));
                    (*U_draw_data)(idx, 1) = 0.5 * ((*U_data)(idx_b) + (*U_data)(idx_u));
                }
            }
        }

        pout << "\n"
             << "Error in " << Sxx_var->getName() << " at time " << loop_time << ":\n"
             << "  L1-norm:  " << coarse_hier_cc_data_ops.L1Norm(Sxx_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << coarse_hier_cc_data_ops.L2Norm(Sxx_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << coarse_hier_cc_data_ops.maxNorm(Sxx_idx, wgt_cc_idx) << "\n";

        pout << "\n"
             << "Error in " << Syy_var->getName() << " at time " << loop_time << ":\n"
             << "  L1-norm:  " << coarse_hier_cc_data_ops.L1Norm(Syy_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << coarse_hier_cc_data_ops.L2Norm(Syy_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << coarse_hier_cc_data_ops.maxNorm(Syy_idx, wgt_cc_idx) << "\n";

        pout << "\n"
             << "Error in " << Sxy_var->getName() << " at time " << loop_time << ":\n"
             << "  L1-norm:  " << coarse_hier_cc_data_ops.L1Norm(Sxy_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << coarse_hier_cc_data_ops.L2Norm(Sxy_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << coarse_hier_cc_data_ops.maxNorm(Sxy_idx, wgt_cc_idx) << "\n";

        pout << "\n"
             << "Error in " << U_var->getName() << " at time " << loop_time << ":\n"
             << "  L1-norm:  " << coarse_hier_sc_data_ops.L1Norm(U_interp_idx, wgt_sc_idx) << "\n"
             << "  L2-norm:  " << coarse_hier_sc_data_ops.L2Norm(U_interp_idx, wgt_sc_idx) << "\n"
             << "  max-norm: " << coarse_hier_sc_data_ops.maxNorm(U_interp_idx, wgt_sc_idx) << "\n";

        pout << "\n"
             << "Error in " << P_var->getName() << " at time " << loop_time << ":\n"
             << "  L1-norm:  " << coarse_hier_cc_data_ops.L1Norm(P_interp_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << coarse_hier_cc_data_ops.L2Norm(P_interp_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << coarse_hier_cc_data_ops.maxNorm(P_interp_idx, wgt_cc_idx) << "\n";

        // Output plot data after taking norms of differences.
        visit_data_writer->writePlotData(coarse_patch_hierarchy, coarse_iteration_num + 1, loop_time);

        pout << endl;
        pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        pout << endl;
    }
    return 0;
} // main
