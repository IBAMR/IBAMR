// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2026 by the IBAMR developers
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

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellIterator.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CCPoissonSolverManager.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LinearSolver.h>
#include <ibtk/muParserCartGridFunction.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

// Verify output for the Poisson example program in IBTK.

namespace
{
// Print the norms of each component of a cell-centered quantity.
void
print_component_norms(const std::string& name,
                      const int data_idx,
                      const int wgt_idx,
                      Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                      const int depth)
{
    for (int d = 0; d < depth; ++d)
    {
        double norm_oo = 0.0, norm_2 = 0.0, norm_1 = 0.0;
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> data = patch->getPatchData(data_idx);
                Pointer<CellData<NDIM, double>> wgt = patch->getPatchData(wgt_idx);
                for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                {
                    const double w = (*wgt)(ci());
                    if (w > 0.0) norm_oo = std::max(norm_oo, std::abs((*data)(ci(), d)));
                    norm_2 += w * std::pow((*data)(ci(), d), 2);
                    norm_1 += w * std::abs((*data)(ci(), d));
                }
            }
        }
        norm_oo = IBTK_MPI::maxReduction(norm_oo);
        norm_2 = std::sqrt(IBTK_MPI::sumReduction(norm_2));
        norm_1 = IBTK_MPI::sumReduction(norm_1);
        plog << "|" << name << "[" << d << "]|_oo = " << norm_oo << "\n";
        plog << "|" << name << "[" << d << "]|_2  = " << norm_2 << "\n";
        plog << "|" << name << "[" << d << "]|_1  = " << norm_1 << "\n";
    }
}
} // namespace

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_poisson.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        const int depth = input_db->getIntegerWithDefault("depth", 1);
        Pointer<CellVariable<NDIM, double>> u_cc_var = new CellVariable<NDIM, double>("u_cc", depth);
        Pointer<CellVariable<NDIM, double>> f_cc_var = new CellVariable<NDIM, double>("f_cc", depth);
        Pointer<CellVariable<NDIM, double>> b_cc_var = new CellVariable<NDIM, double>("b_cc", depth);
        Pointer<CellVariable<NDIM, double>> e_cc_var = new CellVariable<NDIM, double>("e_cc", depth);
        Pointer<CellVariable<NDIM, double>> r_cc_var = new CellVariable<NDIM, double>("r_cc", depth);

        const IntVector<NDIM> gcw(input_db->getIntegerWithDefault("vector_ghost_cell_width", 1));
        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, gcw);
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, gcw);
        const int b_cc_idx = var_db->registerVariableAndContext(b_cc_var, ctx, gcw);
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVector<NDIM>(1));
        const int r_cc_idx = var_db->registerVariableAndContext(r_cc_var, ctx, IntVector<NDIM>(1));

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
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(b_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
            level->allocatePatchData(r_cc_idx, 0.0);
        }

        // Setup vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int h_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> b_vec("b", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> r_vec("r", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        u_vec.addComponent(u_cc_var, u_cc_idx, h_cc_idx);
        f_vec.addComponent(f_cc_var, f_cc_idx, h_cc_idx);
        b_vec.addComponent(b_cc_var, b_cc_idx, h_cc_idx);
        e_vec.addComponent(e_cc_var, e_cc_idx, h_cc_idx);
        r_vec.addComponent(r_cc_var, r_cc_idx, h_cc_idx);

        u_vec.setToScalar(0.0);
        f_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);
        r_vec.setToScalar(1.0);

        // Setup exact solutions.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        muParserCartGridFunction f_fcn("f", app_initializer->getComponentDatabase("f"), grid_geometry);

        f_fcn.setDataOnPatchHierarchy(f_cc_idx, f_cc_var, patch_hierarchy, 0.0);

        // Ensure that the right-hand-side vector has no components in the
        // nullspace of the operator. Keep the right-hand side with its mean in b
        // for a solver that is told the null space.
        b_vec.copyVector(Pointer<SAMRAIVectorReal<NDIM, double>>(&f_vec, false));
        f_vec.addScalar(Pointer<SAMRAIVectorReal<NDIM, double>>(&f_vec, false),
                        -f_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double>>(&r_vec, false)) /
                            r_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double>>(&r_vec, false)));

        // Setup the Poisson solver.
        PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setCZero();
        poisson_spec.setDConstant(-1.0);
        const std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs(depth, nullptr);
        CCLaplaceOperator laplace_op("laplace_op");
        laplace_op.setPoissonSpecifications(poisson_spec);
        laplace_op.setPhysicalBcCoefs(bc_coefs);
        laplace_op.initializeOperatorState(u_vec, f_vec);

        string solver_type = input_db->getString("solver_type");
        Pointer<Database> solver_db = input_db->getDatabase("solver_db");
        string precond_type = input_db->getString("precond_type");
        Pointer<Database> precond_db = input_db->getDatabase("precond_db");
        Pointer<PoissonSolver> poisson_solver = CCPoissonSolverManager::getManager()->allocateSolver(
            solver_type, "poisson_solver", solver_db, "", precond_type, "poisson_precond", precond_db, "");
        poisson_solver->setPoissonSpecifications(poisson_spec);
        poisson_solver->setPhysicalBcCoefs(bc_coefs);

        // With null_space = "AFTER_FIRST_SOLVE", the solver is told after a first solve that the null space of the
        // operator contains the constant vector, and solves again. With null_space = "REINITIALIZE", it is told
        // before it is initialized, and is deallocated and initialized again between two solves.
        auto& linear_solver = dynamic_cast<LinearSolver&>(*poisson_solver);
        const string null_space = input_db->getStringWithDefault("null_space", "");
        if (null_space == "REINITIALIZE")
        {
            linear_solver.setNullSpace(true);
        }
        poisson_solver->initializeSolverState(u_vec, f_vec);

        const auto solve_and_print_norms = [&]()
        {
            // Solve -L*u = f. A solver that is told the null space gets the right-hand side with its mean, b. The
            // residual below is that of f in either case.
            u_vec.setToScalar(0.0);
            if (linear_solver.getNullSpaceContainsConstantVector())
            {
                poisson_solver->solveSystem(u_vec, b_vec);
            }
            else
            {
                poisson_solver->solveSystem(u_vec, f_vec);
            }
            if (!null_space.empty())
            {
                plog << "iterations = " << linear_solver.getNumIterations() << "\n";
            }

            // Compute error and print error norms.
            u_fcn.setDataOnPatchHierarchy(e_cc_idx, e_cc_var, patch_hierarchy, 0.0);
            e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&e_vec, false),
                           Pointer<SAMRAIVectorReal<NDIM, double>>(&u_vec, false));
            plog << "|e|_oo = " << e_vec.maxNorm() << "\n";
            plog << "|e|_2  = " << e_vec.L2Norm() << "\n";
            plog << "|e|_1  = " << e_vec.L1Norm() << "\n";
            if (depth > 1) print_component_norms("e", e_cc_idx, h_cc_idx, patch_hierarchy, depth);

            // Compute the residual and print residual norms.
            laplace_op.apply(u_vec, r_vec);
            r_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&f_vec, false),
                           Pointer<SAMRAIVectorReal<NDIM, double>>(&r_vec, false));
            plog << "|r|_oo = " << r_vec.maxNorm() << "\n";
            plog << "|r|_2  = " << r_vec.L2Norm() << "\n";
            plog << "|r|_1  = " << r_vec.L1Norm() << "\n";
            if (depth > 1) print_component_norms("r", r_cc_idx, h_cc_idx, patch_hierarchy, depth);
        };
        solve_and_print_norms();
        if (null_space == "AFTER_FIRST_SOLVE")
        {
            linear_solver.setNullSpace(true);
            solve_and_print_norms();
        }
        else if (null_space == "REINITIALIZE")
        {
            poisson_solver->deallocateSolverState();
            poisson_solver->initializeSolverState(u_vec, f_vec);
            solve_and_print_norms();
        }

        // Set invalid values on coarse levels (i.e., coarse-grid values that
        // are covered by finer grid patches) to equal zero.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber() - 1; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            BoxArray<NDIM> refined_region_boxes;
            Pointer<PatchLevel<NDIM>> next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, double>> e_cc_data = patch->getPatchData(e_cc_idx);
                Pointer<CellData<NDIM, double>> r_cc_data = patch->getPatchData(r_cc_idx);
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const Box<NDIM> refined_box = refined_region_boxes[i];
                    const Box<NDIM> intersection = Box<NDIM>::grow(patch_box, 1) * refined_box;
                    if (!intersection.empty())
                    {
                        e_cc_data->fillAll(0.0, intersection);
                        r_cc_data->fillAll(0.0, intersection);
                    }
                }
            }
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main
