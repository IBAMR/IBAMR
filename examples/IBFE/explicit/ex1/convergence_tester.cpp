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
#include <SAMRAI_config.h>

// PETSC INCLUDES
#include <petsc.h>

// IBTK INCLUDES
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/HierarchyMathOps.h>

// LIBMESH INCLUDES
#include <libmesh/equation_systems.h>
#include <libmesh/exact_solution.h>
#include <libmesh/mesh.h>
using namespace libMesh;

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <ComponentSelector.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <PatchHierarchy.h>
#include <SideVariable.h>
#include <VariableDatabase.h>
#include <VisItDataWriter.h>
#include <tbox/Database.h>
#include <tbox/HDFDatabase.h>
#include <tbox/InputDatabase.h>
#include <tbox/InputManager.h>
#include <tbox/MathUtilities.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>
#include <tbox/SAMRAIManager.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/Utilities.h>

using namespace IBTK;
using namespace SAMRAI;
using namespace std;

int
main(int argc, char* argv[])
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    {
        tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
        tbox::SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
        tbox::SAMRAIManager::startup();

        if (argc != 2)
        {
            tbox::pout << "USAGE:  " << argv[0] << " <input filename>\n"
                       << "  options:\n"
                       << "  none at this time" << endl;
            tbox::SAMRAI_MPI::abort();
            return (-1);
        }

        string input_filename = argv[1];
        tbox::plog << "input_filename = " << input_filename << endl;

        // Create input database and parse all data in input file.
        tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
        tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

        // Retrieve "Main" section of the input database.
        tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

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

        // Create major algorithm and data objects which comprise application.
        tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geom =
            new geom::CartesianGridGeometry<NDIM>("CartesianGeometry", input_db->getDatabase("CartesianGeometry"));

        // Initialize variables.
        hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();

        tbox::Pointer<hier::VariableContext> current_ctx =
            var_db->getContext("INSStaggeredHierarchyIntegrator::CURRENT");
        tbox::Pointer<hier::VariableContext> scratch_ctx =
            var_db->getContext("INSStaggeredHierarchyIntegrator::SCRATCH");

        tbox::Pointer<pdat::SideVariable<NDIM, double> > U_var =
            new pdat::SideVariable<NDIM, double>("INSStaggeredHierarchyIntegrator::U");
        const int U_idx = var_db->registerVariableAndContext(U_var, current_ctx);
        const int U_interp_idx = var_db->registerClonedPatchDataIndex(U_var, U_idx);
        const int U_scratch_idx = var_db->registerVariableAndContext(U_var, scratch_ctx, 2);

        tbox::Pointer<pdat::CellVariable<NDIM, double> > P_var =
            new pdat::CellVariable<NDIM, double>("INSStaggeredHierarchyIntegrator::P");
        //     tbox::Pointer<pdat::CellVariable<NDIM,double> > P_var = new
        //     pdat::CellVariable<NDIM,double>("INSStaggeredHierarchyIntegrator::P_extrap");
        const int P_idx = var_db->registerVariableAndContext(P_var, current_ctx);
        const int P_interp_idx = var_db->registerClonedPatchDataIndex(P_var, P_idx);
        const int P_scratch_idx = var_db->registerVariableAndContext(P_var, scratch_ctx, 2);

        // Set up visualization plot file writer.
        tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer =
            new appu::VisItDataWriter<NDIM>("VisIt Writer", main_db->getString("viz_dump_dirname"), 1);
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

            sprintf(temp_buf, "%05d.samrai.%05d", coarse_iteration_num, tbox::SAMRAI_MPI::getRank());
            string coarse_file_name = coarse_hier_dump_dirname + "/" + "hier_data.";
            coarse_file_name += temp_buf;

            sprintf(temp_buf, "%05d.samrai.%05d", fine_iteration_num, tbox::SAMRAI_MPI::getRank());
            string fine_file_name = fine_hier_dump_dirname + "/" + "hier_data.";
            fine_file_name += temp_buf;

            for (int rank = 0; rank < tbox::SAMRAI_MPI::getNodes(); ++rank)
            {
                if (rank == tbox::SAMRAI_MPI::getRank())
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
                tbox::SAMRAI_MPI::barrier();
            }

            if (!files_exist) break;

            tbox::pout << endl;
            tbox::pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            tbox::pout << "processing data" << endl;
            tbox::pout << "     coarse iteration number = " << coarse_iteration_num << endl;
            tbox::pout << "     fine iteration number = " << fine_iteration_num << endl;
            tbox::pout << "     coarse file name = " << coarse_file_name << endl;
            tbox::pout << "     fine file name = " << fine_file_name << endl;

            // Read in data to post-process.
            hier::ComponentSelector hier_data;
            hier_data.setFlag(U_idx);
            hier_data.setFlag(P_idx);

            tbox::Pointer<tbox::HDFDatabase> coarse_hier_db = new tbox::HDFDatabase("coarse_hier_db");
            coarse_hier_db->open(coarse_file_name);

            tbox::Pointer<hier::PatchHierarchy<NDIM> > coarse_patch_hierarchy =
                new hier::PatchHierarchy<NDIM>("CoarsePatchHierarchy", grid_geom, false);
            coarse_patch_hierarchy->getFromDatabase(coarse_hier_db->getDatabase("PatchHierarchy"), hier_data);

            const double coarse_loop_time = coarse_hier_db->getDouble("loop_time");

            coarse_hier_db->close();

            tbox::Pointer<tbox::HDFDatabase> fine_hier_db = new tbox::HDFDatabase("fine_hier_db");
            fine_hier_db->open(fine_file_name);

            tbox::Pointer<hier::PatchHierarchy<NDIM> > fine_patch_hierarchy = new hier::PatchHierarchy<NDIM>(
                "FinePatchHierarchy", grid_geom->makeRefinedGridGeometry("FineGridGeometry", 2, false), false);
            fine_patch_hierarchy->getFromDatabase(fine_hier_db->getDatabase("PatchHierarchy"), hier_data);

            const double fine_loop_time = fine_hier_db->getDouble("loop_time");

            fine_hier_db->close();

            TBOX_ASSERT(tbox::MathUtilities<double>::equalEps(coarse_loop_time, fine_loop_time));
            loop_time = fine_loop_time;
            tbox::pout << "     loop time = " << loop_time << endl;

            tbox::Pointer<hier::PatchHierarchy<NDIM> > coarsened_fine_patch_hierarchy =
                fine_patch_hierarchy->makeCoarsenedPatchHierarchy("CoarsenedFinePatchHierarchy", 2, false);

            // Setup hierarchy operations objects.
            math::HierarchyCellDataOpsReal<NDIM, double> coarse_hier_cc_data_ops(
                coarse_patch_hierarchy, 0, coarse_patch_hierarchy->getFinestLevelNumber());
            math::HierarchySideDataOpsReal<NDIM, double> coarse_hier_sc_data_ops(
                coarse_patch_hierarchy, 0, coarse_patch_hierarchy->getFinestLevelNumber());
            HierarchyMathOps hier_math_ops("hier_math_ops", coarse_patch_hierarchy);
            hier_math_ops.setPatchHierarchy(coarse_patch_hierarchy);
            hier_math_ops.resetLevels(0, coarse_patch_hierarchy->getFinestLevelNumber());
            const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
            const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

            // Allocate patch data.
            for (int ln = 0; ln <= coarse_patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM> > level = coarse_patch_hierarchy->getPatchLevel(ln);
                level->allocatePatchData(U_interp_idx, loop_time);
                level->allocatePatchData(P_interp_idx, loop_time);
                level->allocatePatchData(U_scratch_idx, loop_time);
                level->allocatePatchData(P_scratch_idx, loop_time);
            }

            for (int ln = 0; ln <= fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM> > level = fine_patch_hierarchy->getPatchLevel(ln);
                level->allocatePatchData(U_interp_idx, loop_time);
                level->allocatePatchData(P_interp_idx, loop_time);
                level->allocatePatchData(U_scratch_idx, loop_time);
                level->allocatePatchData(P_scratch_idx, loop_time);
            }

            for (int ln = 0; ln <= coarsened_fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM> > level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);
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
                tbox::Pointer<hier::PatchLevel<NDIM> > coarser_level = coarse_patch_hierarchy->getPatchLevel(ln - 1);
                tbox::Pointer<hier::PatchLevel<NDIM> > finer_level = coarse_patch_hierarchy->getPatchLevel(ln);

                xfer::CoarsenAlgorithm<NDIM> coarsen_alg;
                tbox::Pointer<xfer::CoarsenOperator<NDIM> > coarsen_op;

                coarsen_op = grid_geom->lookupCoarsenOperator(U_var, "CONSERVATIVE_COARSEN");
                coarsen_alg.registerCoarsen(U_idx, U_idx, coarsen_op);

                coarsen_op = grid_geom->lookupCoarsenOperator(P_var, "CONSERVATIVE_COARSEN");
                coarsen_alg.registerCoarsen(P_idx, P_idx, coarsen_op);

                coarsen_alg.createSchedule(coarser_level, finer_level)->coarsenData();
            }

            // Synchronize the fine hierarchy data.
            for (int ln = fine_patch_hierarchy->getFinestLevelNumber(); ln > 0; --ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM> > coarser_level = fine_patch_hierarchy->getPatchLevel(ln - 1);
                tbox::Pointer<hier::PatchLevel<NDIM> > finer_level = fine_patch_hierarchy->getPatchLevel(ln);

                xfer::CoarsenAlgorithm<NDIM> coarsen_alg;
                tbox::Pointer<xfer::CoarsenOperator<NDIM> > coarsen_op;

                coarsen_op = grid_geom->lookupCoarsenOperator(U_var, "CONSERVATIVE_COARSEN");
                coarsen_alg.registerCoarsen(U_idx, U_idx, coarsen_op);

                coarsen_op = grid_geom->lookupCoarsenOperator(P_var, "CONSERVATIVE_COARSEN");
                coarsen_alg.registerCoarsen(P_idx, P_idx, coarsen_op);

                coarsen_alg.createSchedule(coarser_level, finer_level)->coarsenData();
            }

            // Coarsen data from the fine hierarchy to the coarsened fine hierarchy.
            for (int ln = 0; ln <= fine_patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM> > dst_level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);
                tbox::Pointer<hier::PatchLevel<NDIM> > src_level = fine_patch_hierarchy->getPatchLevel(ln);

                tbox::Pointer<xfer::CoarsenOperator<NDIM> > coarsen_op;
                for (hier::PatchLevel<NDIM>::Iterator p(dst_level); p; p++)
                {
                    tbox::Pointer<hier::Patch<NDIM> > dst_patch = dst_level->getPatch(p());
                    tbox::Pointer<hier::Patch<NDIM> > src_patch = src_level->getPatch(p());
                    const hier::Box<NDIM>& coarse_box = dst_patch->getBox();
                    TBOX_ASSERT(hier::Box<NDIM>::coarsen(src_patch->getBox(), 2) == coarse_box);

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
                tbox::Pointer<hier::PatchLevel<NDIM> > dst_level = coarse_patch_hierarchy->getPatchLevel(ln);
                tbox::Pointer<hier::PatchLevel<NDIM> > src_level = coarsened_fine_patch_hierarchy->getPatchLevel(ln);

                xfer::RefineAlgorithm<NDIM> refine_alg;
                tbox::Pointer<xfer::RefineOperator<NDIM> > refine_op;

                refine_op = grid_geom->lookupRefineOperator(U_var, "CONSERVATIVE_LINEAR_REFINE");
                refine_alg.registerRefine(U_interp_idx, U_interp_idx, U_scratch_idx, refine_op);

                refine_op = grid_geom->lookupRefineOperator(P_var, "LINEAR_REFINE");
                refine_alg.registerRefine(P_interp_idx, P_interp_idx, P_scratch_idx, refine_op);

                hier::ComponentSelector data_indices;
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

            tbox::pout << "\n"
                       << "Error in " << U_var->getName() << " at time " << loop_time << ":\n"
                       << "  L1-norm:  " << coarse_hier_sc_data_ops.L1Norm(U_interp_idx, wgt_sc_idx) << "\n"
                       << "  L2-norm:  " << coarse_hier_sc_data_ops.L2Norm(U_interp_idx, wgt_sc_idx) << "\n"
                       << "  max-norm: " << coarse_hier_sc_data_ops.maxNorm(U_interp_idx, wgt_sc_idx) << "\n";

            tbox::pout << "\n"
                       << "Error in " << P_var->getName() << " at time " << loop_time << ":\n"
                       << "  L1-norm:  " << coarse_hier_cc_data_ops.L1Norm(P_interp_idx, wgt_cc_idx) << "\n"
                       << "  L2-norm:  " << coarse_hier_cc_data_ops.L2Norm(P_interp_idx, wgt_cc_idx) << "\n"
                       << "  max-norm: " << coarse_hier_cc_data_ops.maxNorm(P_interp_idx, wgt_cc_idx) << "\n";

            // Output plot data after taking norms of differences.
            visit_data_writer->writePlotData(coarse_patch_hierarchy, coarse_iteration_num + 1, loop_time);

            // Do the same thing for the FE data.
            string file_name;

            Mesh mesh_coarse(init.comm(), NDIM);
            file_name = coarse_hier_dump_dirname + "/" + "fe_mesh.";
            sprintf(temp_buf, "%05d", coarse_iteration_num);
            file_name += temp_buf;
            file_name += ".xda";
            mesh_coarse.read(file_name);

            Mesh mesh_fine(init.comm(), NDIM);
            file_name = fine_hier_dump_dirname + "/" + "fe_mesh.";
            sprintf(temp_buf, "%05d", fine_iteration_num);
            file_name += temp_buf;
            file_name += ".xda";
            mesh_fine.read(file_name);

            EquationSystems equation_systems_coarse(mesh_coarse);
            file_name = coarse_hier_dump_dirname + "/" + "fe_equation_systems.";
            sprintf(temp_buf, "%05d", coarse_iteration_num);
            file_name += temp_buf;
            equation_systems_coarse.read(
                file_name,
                (EquationSystems::READ_HEADER | EquationSystems::READ_DATA | EquationSystems::READ_ADDITIONAL_DATA));

            EquationSystems equation_systems_fine(mesh_fine);
            file_name = fine_hier_dump_dirname + "/" + "fe_equation_systems.";
            sprintf(temp_buf, "%05d", fine_iteration_num);
            file_name += temp_buf;
            equation_systems_fine.read(
                file_name,
                (EquationSystems::READ_HEADER | EquationSystems::READ_DATA | EquationSystems::READ_ADDITIONAL_DATA));

            ExactSolution error_estimator(equation_systems_coarse);
            error_estimator.attach_reference_solution(&equation_systems_fine);

            error_estimator.compute_error("IB coordinates system", "X_0");
            double X0_error[3];
            X0_error[0] = error_estimator.l1_error("IB coordinates system", "X_0");
            X0_error[1] = error_estimator.l2_error("IB coordinates system", "X_0");
            X0_error[2] = error_estimator.l_inf_error("IB coordinates system", "X_0");

            error_estimator.compute_error("IB coordinates system", "X_1");
            double X1_error[3];
            X1_error[0] = error_estimator.l1_error("IB coordinates system", "X_1");
            X1_error[1] = error_estimator.l2_error("IB coordinates system", "X_1");
            X1_error[2] = error_estimator.l_inf_error("IB coordinates system", "X_1");

            double X_error[3];
            X_error[0] = X0_error[0] + X1_error[0];
            X_error[1] = sqrt(X0_error[1] * X0_error[1] + X1_error[1] * X1_error[1]);
            X_error[2] = max(X0_error[2], X1_error[2]);

            tbox::pout << "\n"
                       << "Error in X at time " << loop_time << ":\n"
                       << "  L1-norm:  " << X_error[0] << "\n"
                       << "  L2-norm:  " << X_error[1] << "\n"
                       << "  max-norm: " << X_error[2] << "\n";

            tbox::pout << endl;
            tbox::pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            tbox::pout << endl;
        }

        tbox::SAMRAIManager::shutdown();
    }
    return 0;
} // main
