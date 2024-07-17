// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2022 by the IBAMR developers
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
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// C++ includes
#include <numeric>

// Headers for basic PETSc objects
#include "petscao.h"
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <SAMRAIVectorReal.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/muParserCartGridFunction.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "sc_prolongation.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometryNd> grid_geometry = new CartesianGridGeometryNd(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchyNd> patch_hierarchy = new PatchHierarchyNd("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitializeNd> error_detector = new StandardTagAndInitializeNd(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsosNd> box_generator = new BergerRigoutsosNd();
        Pointer<LoadBalancerNd> load_balancer =
            new LoadBalancerNd("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithmNd> gridding_algorithm =
            new GriddingAlgorithmNd("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

        // Create variables and register them with the variable database.
        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        Pointer<SideVariableNd<double> > u_sc_var = new SideVariableNd<double>("u_sc");
        Pointer<SideVariableNd<double> > f_sc_var = new SideVariableNd<double>("f_sc");
        Pointer<SideVariableNd<double> > e_sc_var = new SideVariableNd<double>("e_sc");

        const int u_sc_idx = var_db->registerVariableAndContext(u_sc_var, ctx, IntVectorNd(1));
        const int f_sc_idx = var_db->registerVariableAndContext(f_sc_var, ctx, IntVectorNd(1));
        const int e_sc_idx = var_db->registerVariableAndContext(e_sc_var, ctx, IntVectorNd(1));

        Pointer<CellVariableNd<double> > u_cc_var = new CellVariableNd<double>("u_cc", NDIM);
        Pointer<CellVariableNd<double> > f_cc_var = new CellVariableNd<double>("f_cc", NDIM);
        Pointer<CellVariableNd<double> > e_cc_var = new CellVariableNd<double>("e_cc", NDIM);

        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, IntVectorNd(0));
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, IntVectorNd(0));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVectorNd(0));

        // Register variables for plotting.
        Pointer<VisItDataWriterNd> visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity(u_cc_var->getName(), "VECTOR", u_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(u_cc_var->getName() + std::to_string(d), "SCALAR", u_cc_idx, d);
        }

        visit_data_writer->registerPlotQuantity(f_cc_var->getName(), "VECTOR", f_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(f_cc_var->getName() + std::to_string(d), "SCALAR", f_cc_idx, d);
        }

        visit_data_writer->registerPlotQuantity(e_cc_var->getName(), "VECTOR", e_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(e_cc_var->getName() + std::to_string(d), "SCALAR", e_cc_idx, d);
        }

        // Initialize the AMR patch hierarchy.
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        int tag_buffer = 1;
        int level_number = 0;
        bool done = false;
        while (!done && (gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_sc_idx, 0.0);
            level->allocatePatchData(f_sc_idx, 0.0);
            level->allocatePatchData(e_sc_idx, 0.0);
            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
        }

        // Setup exact solutions.
        muParserCartGridFunction fcn("f", app_initializer->getComponentDatabase("f"), grid_geometry);
        fcn.setDataOnPatchHierarchy(u_sc_idx, u_sc_var, patch_hierarchy, 0.0);
        fcn.setDataOnPatchHierarchy(f_sc_idx, f_sc_var, patch_hierarchy, 0.0);
        fcn.setDataOnPatchHierarchy(e_sc_idx, e_sc_var, patch_hierarchy, 0.0);

        // Compute u DOFs per processor.
        std::vector<std::vector<int> > num_dofs_per_proc;
        Pointer<SideVariableNd<int> > u_dof_index_var = new SideVariableNd<int>("u_dof_index");
        ;
        const IntVectorNd no_ghosts = 0;
        const int u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, ctx, no_ghosts);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        if (finest_ln != 1)
        {
            TBOX_ERROR("This is a 2 level example \n. Please set max_levels = 2 in the input file \n");
        }
        num_dofs_per_proc.resize(finest_ln + 1);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_dof_index_idx, 0.0);
            PETScVecUtilities::constructPatchLevelDOFIndices(num_dofs_per_proc[ln], u_dof_index_idx, level);
        }

        // Construct the coarse and fine level PETSc Vecs
        const int mpi_rank = IBTK_MPI::getRank();
        const int n_local_coarsest = num_dofs_per_proc[coarsest_ln][mpi_rank];
        const int n_total_coarsest =
            std::accumulate(num_dofs_per_proc[coarsest_ln].begin(), num_dofs_per_proc[coarsest_ln].end(), 0);
        const int n_local_finest = num_dofs_per_proc[finest_ln][mpi_rank];
        const int n_total_finest =
            std::accumulate(num_dofs_per_proc[finest_ln].begin(), num_dofs_per_proc[finest_ln].end(), 0);

        Vec x;
        VecCreateMPI(PETSC_COMM_WORLD, n_local_coarsest, n_total_coarsest, &x);
        Vec X, Y, E;
        VecCreateMPI(PETSC_COMM_WORLD, n_local_finest, n_total_finest, &X);
        PETScVecUtilities::copyToPatchLevelVec(
            x, u_sc_idx, u_dof_index_idx, patch_hierarchy->getPatchLevel(coarsest_ln));
        PETScVecUtilities::copyToPatchLevelVec(X, u_sc_idx, u_dof_index_idx, patch_hierarchy->getPatchLevel(finest_ln));
        VecDuplicate(X, &Y);
        VecDuplicate(X, &E);

        // Construct PETSc prolongation mat.
        std::string prolongation_op_type = input_db->getString("prolongation_op_type");
        AO coarse_level_ao = NULL;
        PETScVecUtilities::constructPatchLevelAO(coarse_level_ao,
                                                 num_dofs_per_proc[coarsest_ln],
                                                 u_dof_index_idx,
                                                 patch_hierarchy->getPatchLevel(coarsest_ln),
                                                 /*ao_offset*/ 0);
        Mat prolongation_mat = NULL;
        PETScMatUtilities::constructProlongationOp(prolongation_mat,
                                                   prolongation_op_type,
                                                   u_dof_index_idx,
                                                   num_dofs_per_proc[finest_ln],
                                                   num_dofs_per_proc[coarsest_ln],
                                                   patch_hierarchy->getPatchLevel(finest_ln),
                                                   patch_hierarchy->getPatchLevel(coarsest_ln),
                                                   coarse_level_ao,
                                                   /*coarse_ao_offset*/ 0);

        // Obtain finest DOFs from prolongation matrix.
        // P * x = Y
        MatMult(prolongation_mat, x, Y);

        PetscViewer matlab_viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, "PROLONG.MAT", FILE_MODE_WRITE, &matlab_viewer);
        PetscViewerPushFormat(matlab_viewer, PETSC_VIEWER_ASCII_MATLAB);
        MatView(prolongation_mat, matlab_viewer);
        PetscViewerDestroy(&matlab_viewer);

        // Compute error of prolongation E = Y - X, and find its norm.
        double norm_1, norm_2, norm_oo;
        VecWAXPY(E, -1.0, X, Y);
        VecNorm(E, NORM_1, &norm_1);
        VecNorm(E, NORM_2, &norm_2);
        VecNorm(E, NORM_INFINITY, &norm_oo);
        pout << "|| E ||_1 of the error vector  = " << norm_1 << "\n"
             << "|| E ||_2 of the error vector  = " << norm_2 << "\n"
             << "|| E ||_oo of the error vector = " << norm_oo << "\n";

        // Copy PETSc Vec to SAMRAI vec
        Pointer<RefineScheduleNd> data_synch_sched =
            PETScVecUtilities::constructDataSynchSchedule(f_sc_idx, patch_hierarchy->getPatchLevel(finest_ln));
        PETScVecUtilities::copyFromPatchLevelVec(Y,
                                                 f_sc_idx,
                                                 u_dof_index_idx,
                                                 patch_hierarchy->getPatchLevel(finest_ln),
                                                 data_synch_sched,
                                                 Pointer<RefineScheduleNd>(NULL));

        // Setup SAMRAI vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int h_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        SAMRAIVectorRealNd<double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorRealNd<double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        f_vec.addComponent(f_sc_var, f_sc_idx, h_sc_idx);
        e_vec.addComponent(e_sc_var, e_sc_idx, h_sc_idx);

        // Compute error and print error norms.
        e_vec.subtract(Pointer<SAMRAIVectorRealNd<double> >(&e_vec, false),
                       Pointer<SAMRAIVectorRealNd<double> >(&f_vec, false));
        pout << "|e|_oo = " << e_vec.maxNorm() << "\n";
        pout << "|e|_2  = " << e_vec.L2Norm() << "\n";
        pout << "|e|_1  = " << e_vec.L1Norm() << "\n";

        // Interpolate the side-centered data to cell centers for output.
        static const bool synch_cf_interface = false;
        hier_math_ops.interp(u_cc_idx,
                             u_cc_var,
                             u_sc_idx,
                             u_sc_var,
                             Pointer<HierarchyGhostCellInterpolation>(NULL),
                             0.0,
                             synch_cf_interface);
        hier_math_ops.interp(f_cc_idx,
                             f_cc_var,
                             f_sc_idx,
                             f_sc_var,
                             Pointer<HierarchyGhostCellInterpolation>(NULL),
                             0.0,
                             synch_cf_interface);
        hier_math_ops.interp(e_cc_idx,
                             e_cc_var,
                             e_sc_idx,
                             e_sc_var,
                             Pointer<HierarchyGhostCellInterpolation>(NULL),
                             0.0,
                             synch_cf_interface);

        // Set invalid values on coarse levels (i.e., coarse-grid values that
        // are covered by finer grid patches) to equal zero.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber() - 1; ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            BoxArrayNd refined_region_boxes;
            Pointer<PatchLevelNd> next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                Pointer<PatchNd> patch = level->getPatch(p());
                const BoxNd& patch_box = patch->getBox();
                Pointer<CellDataNd<double> > e_cc_data = patch->getPatchData(e_cc_idx);
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const BoxNd refined_box = refined_region_boxes[i];
                    const BoxNd intersection = BoxNd::grow(patch_box, 1) * refined_box;
                    if (!intersection.empty())
                    {
                        e_cc_data->fillAll(0.0, intersection);
                    }
                }
            }
        }

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

    } // cleanup dynamically allocated objects prior to shutdown
} // main
