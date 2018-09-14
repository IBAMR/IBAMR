// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// GENERAL CONFIGURATION
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>
#include <petscsys.h>

// IBAMR INCLUDES
#include <ibamr/app_namespaces.h>

// IBTK INCLUDES
#include <ibtk/AppInitializer.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/HierarchyMathOps.h>

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
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    // Parse command line options, set some standard options from the input
    // file, and enable file logging.
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    // Retrieve "Main" section of the input database.
    Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");

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
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

    // Initialize variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    Pointer<VariableContext> current_ctx = var_db->getContext("INSStaggeredHierarchyIntegrator::CURRENT");
    Pointer<VariableContext> scratch_ctx = var_db->getContext("INSStaggeredHierarchyIntegrator::SCRATCH");

    Pointer<SideVariable<NDIM, double> > U_var = new SideVariable<NDIM, double>("INSStaggeredHierarchyIntegrator::U");
    const int U_idx = var_db->registerVariableAndContext(U_var, current_ctx);
    const int U_interp_idx = var_db->registerClonedPatchDataIndex(U_var, U_idx);
    const int U_scratch_idx = var_db->registerVariableAndContext(U_var, scratch_ctx, 2);

    Pointer<CellVariable<NDIM, double> > P_var = new CellVariable<NDIM, double>("INSStaggeredHierarchyIntegrator::P");
    const int P_idx = var_db->registerVariableAndContext(P_var, current_ctx);
    const int P_interp_idx = var_db->registerClonedPatchDataIndex(P_var, P_idx);
    const int P_scratch_idx = var_db->registerVariableAndContext(P_var, scratch_ctx, 2);

    // Set up visualization plot file writer.
    Pointer<VisItDataWriter<NDIM> > visit_data_writer =
        new VisItDataWriter<NDIM>("VisIt Writer", main_db->getString("viz_dump_dirname"), 1);
    visit_data_writer->registerPlotQuantity("P", "SCALAR", P_idx);
    visit_data_writer->registerPlotQuantity("P interp", "SCALAR", P_interp_idx);

    // Time step loop.
    double loop_time = 0.0;
    int coarse_iteration_num = coarse_hier_dump_interval;
    int fine_iteration_num = fine_hier_dump_interval;

    bool files_exist = true;
    for (; files_exist;
         coarse_iteration_num += coarse_hier_dump_interval, fine_iteration_num += fine_hier_dump_interval)
    {
        char temp_buf[128];

        sprintf(temp_buf, "%05d.samrai.%05d", coarse_iteration_num, SAMRAI_MPI::getRank());
        string coarse_file_name = coarse_hier_dump_dirname + "/" + "hier_data.";
        coarse_file_name += temp_buf;

        sprintf(temp_buf, "%05d.samrai.%05d", fine_iteration_num, SAMRAI_MPI::getRank());
        string fine_file_name = fine_hier_dump_dirname + "/" + "hier_data.";
        fine_file_name += temp_buf;

        for (int rank = 0; rank < SAMRAI_MPI::getNodes(); ++rank)
        {
            if (rank == SAMRAI_MPI::getRank())
            {
                fstream coarse_fin, fine_fin;
                coarse_fin.open(coarse_file_name.c_str(), ios::in);
                fine_fin.open(fine_file_name.c_str(), ios::in);
                if (!coarse_fin.is_open() || !fine_fin.is_open())
                {
                    files_exist = false;
                }
                coarse_fin.close();
                fine_fin.close();
            }
            SAMRAI_MPI::barrier();
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

        Pointer<HDFDatabase> coarse_hier_db = new HDFDatabase("coarse_hier_db");
        coarse_hier_db->open(coarse_file_name);

        Pointer<PatchHierarchy<NDIM> > coarse_patch_hierarchy =
            new PatchHierarchy<NDIM>("CoarsePatchHierarchy", grid_geom, false);
        coarse_patch_hierarchy->getFromDatabase(coarse_hier_db->getDatabase("PatchHierarchy"), hier_data);

        const double coarse_loop_time = coarse_hier_db->getDouble("loop_time");

        coarse_hier_db->close();

        Pointer<HDFDatabase> fine_hier_db = new HDFDatabase("fine_hier_db");
        fine_hier_db->open(fine_file_name);

        Pointer<PatchHierarchy<NDIM> > fine_patch_hierarchy = new PatchHierarchy<NDIM>(
            "FinePatchHierarchy", grid_geom->makeRefinedGridGeometry("FineGridGeometry", 2, false), false);
        fine_patch_hierarchy->getFromDatabase(fine_hier_db->getDatabase("PatchHierarchy"), hier_data);

        const double fine_loop_time = fine_hier_db->getDouble("loop_time");

        fine_hier_db->close();

        TBOX_ASSERT(MathUtilities<double>::equalEps(coarse_loop_time, fine_loop_time));
        loop_time = fine_loop_time;
        pout << "     loop time = " << loop_time << endl;

        Pointer<PatchHierarchy<NDIM> > coarsened_fine_patch_hierarchy =
            fine_patch_hierarchy->makeCoarsenedPatchHierarchy("CoarsenedFinePatchHierarchy", 2, false);

        // Setup hierarchy operations objects.
        HierarchyCellDataOpsReal<NDIM, double> coarse_hier_cc_data_ops(
            coarse_patch_hierarchy, 0, coarse_patch_hierarchy->getFinestLevelNumber());
        HierarchySideDataOpsReal<NDIM, double> coarse_hier_sc_data_ops(
            coarse_patch_hierarchy, 0, coarse_patch_hierarchy->getFinestLevelNumber());
        HierarchyMathOps hier_math_ops("hier_math_ops", coarse_patch_hierarchy);
        hier_math_ops.setPatchHierarchy(coarse_patch_hierarchy);
        hier_math_ops.resetLevels(0, coarse_patch_hierarchy->getFinestLevelNumber());
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        // Allocate patch data.
        for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = coarse_patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(U_interp_idx, loop_time);
            level->allocatePatchData(P_interp_idx, loop_time);
            level->allocatePatchData(U_scratch_idx, loop_time);
            level->allocatePatchData(P_scratch_idx, loop_time);
        }

        for (int ln = 0; ln <= fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = fine_patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(U_interp_idx, loop_time);
            level->allocatePatchData(P_interp_idx, loop_time);
            level->allocatePatchData(U_scratch_idx, loop_time);
            level->allocatePatchData(P_scratch_idx, loop_time);
        }

        for (int ln = 0; ln <= coarsened_fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(U_idx, loop_time);
            level->allocatePatchData(P_idx, loop_time);
            level->allocatePatchData(U_interp_idx, loop_time);
            level->allocatePatchData(P_interp_idx, loop_time);
            level->allocatePatchData(U_scratch_idx, loop_time);
            level->allocatePatchData(P_scratch_idx, loop_time);
        }

        // Synchronize the coarse hierarchy data.
        for (int ln = coarse_patch_hierarchy->getFinestLevelNumber(); ln > 0; --ln)
        {
            Pointer<PatchLevel<NDIM> > coarser_level = coarse_patch_hierarchy->getPatchLevel(ln - 1);
            Pointer<PatchLevel<NDIM> > finer_level = coarse_patch_hierarchy->getPatchLevel(ln);

            CoarsenAlgorithm<NDIM> coarsen_alg;
            Pointer<CoarsenOperator<NDIM> > coarsen_op;

            coarsen_op = grid_geom->lookupCoarsenOperator(U_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(U_idx, U_idx, coarsen_op);

            coarsen_op = grid_geom->lookupCoarsenOperator(P_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(P_idx, P_idx, coarsen_op);

            coarsen_alg.createSchedule(coarser_level, finer_level)->coarsenData();
        }

        // Synchronize the fine hierarchy data.
        for (int ln = fine_patch_hierarchy->getFinestLevelNumber(); ln > 0; --ln)
        {
            Pointer<PatchLevel<NDIM> > coarser_level = fine_patch_hierarchy->getPatchLevel(ln - 1);
            Pointer<PatchLevel<NDIM> > finer_level = fine_patch_hierarchy->getPatchLevel(ln);

            CoarsenAlgorithm<NDIM> coarsen_alg;
            Pointer<CoarsenOperator<NDIM> > coarsen_op;

            coarsen_op = grid_geom->lookupCoarsenOperator(U_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(U_idx, U_idx, coarsen_op);

            coarsen_op = grid_geom->lookupCoarsenOperator(P_var, "CONSERVATIVE_COARSEN");
            coarsen_alg.registerCoarsen(P_idx, P_idx, coarsen_op);

            coarsen_alg.createSchedule(coarser_level, finer_level)->coarsenData();
        }

        // Coarsen data from the fine hierarchy to the coarsened fine hierarchy.
        for (int ln = 0; ln <= fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > dst_level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);
            Pointer<PatchLevel<NDIM> > src_level = fine_patch_hierarchy->getPatchLevel(ln);

            Pointer<CoarsenOperator<NDIM> > coarsen_op;
            for (PatchLevel<NDIM>::Iterator p(dst_level); p; p++)
            {
                Pointer<Patch<NDIM> > dst_patch = dst_level->getPatch(p());
                Pointer<Patch<NDIM> > src_patch = src_level->getPatch(p());
                const Box<NDIM>& coarse_box = dst_patch->getBox();
                TBOX_ASSERT(Box<NDIM>::coarsen(src_patch->getBox(), 2) == coarse_box);

                coarsen_op = grid_geom->lookupCoarsenOperator(U_var, "CONSERVATIVE_COARSEN");
                coarsen_op->coarsen(*dst_patch, *src_patch, U_interp_idx, U_idx, coarse_box, 2);

                coarsen_op = grid_geom->lookupCoarsenOperator(P_var, "CONSERVATIVE_COARSEN");
                coarsen_op->coarsen(*dst_patch, *src_patch, P_interp_idx, P_idx, coarse_box, 2);
            }
        }

        // Interpolate and copy data from the coarsened fine patch hierarchy to
        // the coarse patch hierarchy.
        for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > dst_level = coarse_patch_hierarchy->getPatchLevel(ln);
            Pointer<PatchLevel<NDIM> > src_level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);

            RefineAlgorithm<NDIM> refine_alg;
            Pointer<RefineOperator<NDIM> > refine_op;

            refine_op = grid_geom->lookupRefineOperator(U_var, "CONSERVATIVE_LINEAR_REFINE");
            refine_alg.registerRefine(U_interp_idx, U_interp_idx, U_scratch_idx, refine_op);

            refine_op = grid_geom->lookupRefineOperator(P_var, "LINEAR_REFINE");
            refine_alg.registerRefine(P_interp_idx, P_interp_idx, P_scratch_idx, refine_op);

            ComponentSelector data_indices;
            data_indices.setFlag(U_scratch_idx);
            data_indices.setFlag(P_scratch_idx);
            CartExtrapPhysBdryOp bc_helper(data_indices, "LINEAR");

            refine_alg.createSchedule(dst_level, src_level, ln - 1, coarse_patch_hierarchy, &bc_helper)
                ->fillData(loop_time);
        }

        // Output plot data before taking norms of differences.
        visit_data_writer->writePlotData(coarse_patch_hierarchy, coarse_iteration_num, loop_time);

        // Compute norms of differences.
        coarse_hier_sc_data_ops.subtract(U_interp_idx, U_idx, U_interp_idx);
        coarse_hier_cc_data_ops.subtract(P_interp_idx, P_idx, P_interp_idx);

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

    SAMRAIManager::shutdown();
    SAMRAI_MPI::finalize();
    return 0;
} // main
