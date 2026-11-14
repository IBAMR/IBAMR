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

#include <ibtk/AppInitializer.h>
#include <ibtk/CCScaledLaplaceOperator.h>
#include <ibtk/CCScaledPoissonPETScLevelSolver.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SAMRAI_config.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>
#include <VisItDataWriter.h>

#include <iomanip>
#include <string>
#include <vector>

#include <ibtk/app_namespaces.h>

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_scaled_poisson.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        const double data_time = 0.0;
        const bool run_discrete_test = input_db->getBool("run_discrete_test");

        // Construct the hierarchy infrastructure.
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

        Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

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

        // Construct only level zero.
        gridding_algorithm->makeCoarsestLevel(hierarchy, data_time);

        if (hierarchy->getFinestLevelNumber() != 0)
        {
            TBOX_ERROR("This test requires exactly one hierarchy level.\n");
        }

        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);

        // Register solution, forcing, diagnostic, and coefficient data.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = var_db->getContext("TEST");

        Pointer<CellVariable<NDIM, double>> u_var = new CellVariable<NDIM, double>("u");
        Pointer<CellVariable<NDIM, double>> f_var = new CellVariable<NDIM, double>("f");
        Pointer<CellVariable<NDIM, double>> exact_var = new CellVariable<NDIM, double>("u_exact");
        Pointer<CellVariable<NDIM, double>> error_var = new CellVariable<NDIM, double>("error");
        Pointer<CellVariable<NDIM, double>> residual_var = new CellVariable<NDIM, double>("residual");
        Pointer<CellVariable<NDIM, double>> C_var = new CellVariable<NDIM, double>("C");
        Pointer<CellVariable<NDIM, double>> S_var = new CellVariable<NDIM, double>("S");
        Pointer<SideVariable<NDIM, double>> D_var = new SideVariable<NDIM, double>("D");

        const IntVector<NDIM> one_ghost(1);
        const IntVector<NDIM> no_ghosts(0);

        const int u_idx = var_db->registerVariableAndContext(u_var, context, one_ghost);
        const int f_idx = var_db->registerVariableAndContext(f_var, context, one_ghost);
        const int exact_idx = var_db->registerVariableAndContext(exact_var, context, one_ghost);
        const int error_idx = var_db->registerVariableAndContext(error_var, context, one_ghost);
        const int residual_idx = var_db->registerVariableAndContext(residual_var, context, one_ghost);
        const int C_idx = var_db->registerVariableAndContext(C_var, context, no_ghosts);
        const int S_idx = var_db->registerVariableAndContext(S_var, context, no_ghosts);
        const int D_idx = var_db->registerVariableAndContext(D_var, context, no_ghosts);

        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();

        if (visit_data_writer)
        {
            visit_data_writer->registerPlotQuantity("u_solution", "SCALAR", u_idx);
            visit_data_writer->registerPlotQuantity("u_error", "SCALAR", error_idx);
            visit_data_writer->registerPlotQuantity("u_exact", "SCALAR", exact_idx);
        }

        const std::vector<int> data_indices = { u_idx, f_idx, exact_idx, error_idx, residual_idx, C_idx, S_idx, D_idx };

        for (const int idx : data_indices) level->allocatePatchData(idx, data_time);

        // Construct weighted vectors restricted to level zero.
        HierarchyMathOps hier_math_ops("test_math_ops", hierarchy, 0, 0);
        const int weight_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> u_vec("u", hierarchy, 0, 0);
        SAMRAIVectorReal<NDIM, double> f_vec("f", hierarchy, 0, 0);
        SAMRAIVectorReal<NDIM, double> exact_vec("exact", hierarchy, 0, 0);
        SAMRAIVectorReal<NDIM, double> error_vec("error", hierarchy, 0, 0);
        SAMRAIVectorReal<NDIM, double> residual_vec("residual", hierarchy, 0, 0);

        u_vec.addComponent(u_var, u_idx, weight_idx);
        f_vec.addComponent(f_var, f_idx, weight_idx);
        exact_vec.addComponent(exact_var, exact_idx, weight_idx);
        error_vec.addComponent(error_var, error_idx, weight_idx);
        residual_vec.addComponent(residual_var, residual_idx, weight_idx);

        u_vec.setToScalar(0.0);
        f_vec.setToScalar(0.0);
        exact_vec.setToScalar(0.0);
        error_vec.setToScalar(0.0);
        residual_vec.setToScalar(0.0);

        HierarchyCellDataOpsReal<NDIM, double> data_ops(hierarchy, 0, 0);

        // Read all analytic fields from the input database.
        muParserCartGridFunction C_fcn("C_fcn", app_initializer->getComponentDatabase("C"), grid_geometry);
        muParserCartGridFunction S_fcn("S_fcn", app_initializer->getComponentDatabase("S"), grid_geometry);
        muParserCartGridFunction D_fcn("D_fcn", app_initializer->getComponentDatabase("D"), grid_geometry);
        muParserCartGridFunction exact_fcn(
            "exact_fcn", app_initializer->getComponentDatabase("u_exact"), grid_geometry);
        muParserCartGridFunction f_fcn("f_fcn", app_initializer->getComponentDatabase("f"), grid_geometry);

        C_fcn.setDataOnPatchHierarchy(C_idx, C_var, hierarchy, data_time);
        S_fcn.setDataOnPatchHierarchy(S_idx, S_var, hierarchy, data_time);
        D_fcn.setDataOnPatchHierarchy(D_idx, D_var, hierarchy, data_time);
        exact_fcn.setDataOnPatchHierarchy(exact_idx, exact_var, hierarchy, data_time);
        f_fcn.setDataOnPatchHierarchy(f_idx, f_var, hierarchy, data_time);

        // Physical Robin boundary conditions are also specified in the input.
        muParserRobinBcCoefs bc_coefs(
            "bc_coefs", app_initializer->getComponentDatabase("BoundaryConditions"), grid_geometry);

        PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setCPatchDataId(C_idx);
        poisson_spec.setDPatchDataId(D_idx);

        // Independent hierarchy operator for residual and discrete-RHS tests.
        CCScaledLaplaceOperator scaled_op("scaled_op", nullptr, false);
        scaled_op.setPoissonSpecifications(poisson_spec);
        scaled_op.setScalePatchDataIndex(S_idx);
        scaled_op.setPhysicalBcCoef(&bc_coefs);
        scaled_op.setHomogeneousBc(false);
        scaled_op.setSolutionTime(data_time);

        CCScaledPoissonPETScLevelSolver solver(
            "scaled_solver", app_initializer->getComponentDatabase("Solver"), "scaled_poisson_");
        solver.setPoissonSpecifications(poisson_spec);
        solver.setScalePatchDataIndex(S_idx);
        solver.setPhysicalBcCoef(&bc_coefs);
        solver.setHomogeneousBc(false);
        solver.setSolutionTime(data_time);
        solver.setInitialGuessNonzero(false);

        scaled_op.initializeOperatorState(u_vec, f_vec);
        solver.initializeSolverState(u_vec, f_vec);

        pout << std::scientific << std::setprecision(12);
        pout << "\nSingle-level scaled Poisson test\n";
        pout << "number of levels = " << hierarchy->getFinestLevelNumber() + 1 << "\n";

        auto solve_and_report = [&](const std::string& label)
        {
            u_vec.setToScalar(0.0);
            solver.solveSystem(u_vec, f_vec);

            data_ops.subtract(error_idx, exact_idx, u_idx);

            // Include the nonhomogeneous physical boundary conditions.
            scaled_op.apply(u_vec, residual_vec);
            data_ops.subtract(residual_idx, f_idx, residual_idx);

            pout << "\n" << label << "\n";
            pout << "|error|_oo = " << error_vec.maxNorm() << "\n";
            pout << "|error|_2   = " << error_vec.L2Norm() << "\n";
            pout << "|error|_1   = " << error_vec.L1Norm() << "\n";

            const double exact_l2 = exact_vec.L2Norm();
            if (exact_l2 > 0.0)
            {
                pout << "relative L2 error = " << error_vec.L2Norm() / exact_l2 << "\n";
            }

            pout << "|residual|_oo = " << residual_vec.maxNorm() << "\n";
            pout << "|residual|_2   = " << residual_vec.L2Norm() << "\n";
        };

        // Test 1: continuum forcing read from the input file.
        solve_and_report("CONTINUOUS MANUFACTURED SOLUTION");

        if (visit_data_writer)
        {
            visit_data_writer->writePlotData(hierarchy, 0, data_time);
        }

        if (run_discrete_test)
        {
            // Test 2: replace the analytic forcing with the discrete operator action.
            scaled_op.apply(exact_vec, f_vec);
            solve_and_report("DISCRETE MANUFACTURED SOLUTION");
        }

        solver.deallocateSolverState();
        scaled_op.deallocateOperatorState();
    }

    return 0;
}
