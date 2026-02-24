// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2024 by the IBAMR developers
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
#include "ibtk/samrai_compatibility_names.h"

#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include "SAMRAIBergerRigoutsos.h"
#include "SAMRAIBox.h"
#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAICellData.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIGriddingAlgorithm.h"
#include "SAMRAIIntVector.h"
#include "SAMRAILoadBalancer.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIRobinBcCoefStrategy.h"
#include "SAMRAISAMRAIVectorReal.h"
#include "SAMRAIStandardTagAndInitialize.h"
#include "SAMRAIVariableDatabase.h"
#include "SAMRAIVisItDataWriter.h"

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CCPoissonSolverManager.h>
#include <ibtk/IBTKInit.h>
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
        SAMRAIPointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_poisson.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geometry = new SAMRAICartesianGridGeometry(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy = new SAMRAIPatchHierarchy("PatchHierarchy", grid_geometry);
        SAMRAIPointer<SAMRAIStandardTagAndInitialize> error_detector = new SAMRAIStandardTagAndInitialize(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        SAMRAIPointer<SAMRAIBergerRigoutsos> box_generator = new SAMRAIBergerRigoutsos();
        SAMRAIPointer<SAMRAILoadBalancer> load_balancer =
            new SAMRAILoadBalancer("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        SAMRAIPointer<SAMRAIGriddingAlgorithm> gridding_algorithm =
            new SAMRAIGriddingAlgorithm("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create variables and register them with the variable database.
        SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
        SAMRAIPointer<VariableContext> ctx = var_db->getContext("context");

        SAMRAIPointer<SAMRAICellVariable<double>> u_cc_var = new SAMRAICellVariable<double>("u_cc");
        SAMRAIPointer<SAMRAICellVariable<double>> f_cc_var = new SAMRAICellVariable<double>("f_cc");
        SAMRAIPointer<SAMRAICellVariable<double>> e_cc_var = new SAMRAICellVariable<double>("e_cc");
        SAMRAIPointer<SAMRAICellVariable<double>> r_cc_var = new SAMRAICellVariable<double>("r_cc");

        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, SAMRAIIntVector(1));
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, SAMRAIIntVector(1));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, SAMRAIIntVector(1));
        const int r_cc_idx = var_db->registerVariableAndContext(r_cc_var, ctx, SAMRAIIntVector(1));

        // Register variables for plotting.
        SAMRAIPointer<SAMRAIVisItDataWriter> visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity(u_cc_var->getName(), "SCALAR", u_cc_idx);
        visit_data_writer->registerPlotQuantity(f_cc_var->getName(), "SCALAR", f_cc_idx);
        visit_data_writer->registerPlotQuantity(e_cc_var->getName(), "SCALAR", e_cc_idx);
        visit_data_writer->registerPlotQuantity(r_cc_var->getName(), "SCALAR", r_cc_idx);

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
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
            level->allocatePatchData(r_cc_idx, 0.0);
        }

        // Setup vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int h_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAISAMRAIVectorReal<double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAISAMRAIVectorReal<double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAISAMRAIVectorReal<double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAISAMRAIVectorReal<double> r_vec("r", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        u_vec.addComponent(u_cc_var, u_cc_idx, h_cc_idx);
        f_vec.addComponent(f_cc_var, f_cc_idx, h_cc_idx);
        e_vec.addComponent(e_cc_var, e_cc_idx, h_cc_idx);
        r_vec.addComponent(r_cc_var, r_cc_idx, h_cc_idx);

        u_vec.setToScalar(0.0);
        f_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);
        r_vec.setToScalar(1.0);

        // Setup exact solutions.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        muParserCartGridFunction f_fcn("f", app_initializer->getComponentDatabase("f"), grid_geometry);

        u_fcn.setDataOnPatchHierarchy(e_cc_idx, e_cc_var, patch_hierarchy, 0.0);
        f_fcn.setDataOnPatchHierarchy(f_cc_idx, f_cc_var, patch_hierarchy, 0.0);

        // Ensure that the right-hand-side vector has no components in the
        // nullspace of the operator.
        f_vec.addScalar(SAMRAIPointer<SAMRAISAMRAIVectorReal<double>>(&f_vec, false),
                        -f_vec.dot(SAMRAIPointer<SAMRAISAMRAIVectorReal<double>>(&r_vec, false)) /
                            r_vec.dot(SAMRAIPointer<SAMRAISAMRAIVectorReal<double>>(&r_vec, false)));

        // Setup the Poisson solver.
        PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setCZero();
        poisson_spec.setDConstant(-1.0);
        SAMRAIRobinBcCoefStrategy* bc_coef = nullptr;
        CCLaplaceOperator laplace_op("laplace_op");
        laplace_op.setPoissonSpecifications(poisson_spec);
        laplace_op.setPhysicalBcCoef(bc_coef);
        laplace_op.initializeOperatorState(u_vec, f_vec);

        string solver_type = input_db->getString("solver_type");
        SAMRAIPointer<Database> solver_db = input_db->getDatabase("solver_db");
        string precond_type = input_db->getString("precond_type");
        SAMRAIPointer<Database> precond_db = input_db->getDatabase("precond_db");
        SAMRAIPointer<PoissonSolver> poisson_solver = CCPoissonSolverManager::getManager()->allocateSolver(
            solver_type, "poisson_solver", solver_db, "", precond_type, "poisson_precond", precond_db, "");
        poisson_solver->setPoissonSpecifications(poisson_spec);
        poisson_solver->setPhysicalBcCoef(bc_coef);
        poisson_solver->initializeSolverState(u_vec, f_vec);

        // Solve -L*u = f.
        u_vec.setToScalar(0.0);
        poisson_solver->solveSystem(u_vec, f_vec);

        // Compute error and print error norms.
        e_vec.subtract(SAMRAIPointer<SAMRAISAMRAIVectorReal<double>>(&e_vec, false),
                       SAMRAIPointer<SAMRAISAMRAIVectorReal<double>>(&u_vec, false));
        pout << "|e|_oo = " << e_vec.maxNorm() << "\n";
        pout << "|e|_2  = " << e_vec.L2Norm() << "\n";
        pout << "|e|_1  = " << e_vec.L1Norm() << "\n";

        // Compute the residual and print residual norms.
        laplace_op.apply(u_vec, r_vec);
        r_vec.subtract(SAMRAIPointer<SAMRAISAMRAIVectorReal<double>>(&f_vec, false),
                       SAMRAIPointer<SAMRAISAMRAIVectorReal<double>>(&r_vec, false));
        pout << "|r|_oo = " << r_vec.maxNorm() << "\n";
        pout << "|r|_2  = " << r_vec.L2Norm() << "\n";
        pout << "|r|_1  = " << r_vec.L1Norm() << "\n";

        // Set invalid values on coarse levels (i.e., coarse-grid values that
        // are covered by finer grid patches) to equal zero.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber() - 1; ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            BoxArray<NDIM> refined_region_boxes;
            SAMRAIPointer<SAMRAIPatchLevel> next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                const SAMRAIBox& patch_box = patch->getBox();
                SAMRAIPointer<SAMRAICellData<double>> e_cc_data = patch->getPatchData(e_cc_idx);
                SAMRAIPointer<SAMRAICellData<double>> r_cc_data = patch->getPatchData(r_cc_idx);
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const SAMRAIBox refined_box = refined_region_boxes[i];
                    const SAMRAIBox intersection = SAMRAIBox::grow(patch_box, 1) * refined_box;
                    if (!intersection.empty())
                    {
                        e_cc_data->fillAll(0.0, intersection);
                        r_cc_data->fillAll(0.0, intersection);
                    }
                }
            }
        }

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

    } // cleanup dynamically allocated objects prior to shutdown
} // main
