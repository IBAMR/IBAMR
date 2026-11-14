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
#include <ibtk/CCPoissonSolverManager.h>
#include <ibtk/CCScaledLaplaceOperator.h>
#include <ibtk/CCScaledPoissonPETScLevelSolver.h>
#include <ibtk/CCScaledPoissonPointRelaxationFACOperator.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/KrylovLinearSolver.h>
#include <ibtk/PETScKrylovPoissonSolver.h>
#include <ibtk/PoissonFACPreconditioner.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <IntVector.h>
#include <LoadBalancer.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SAMRAI_config.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>
#include <VisItDataWriter.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

#include <ibtk/app_namespaces.h>

#if (NDIM != 2)
#error "This test and its input file require NDIM == 2."
#endif

namespace
{
struct CompositeNorms
{
    double max_norm = 0.0;
    double L2_norm = 0.0;
    double L1_norm = 0.0;
};

Pointer<PoissonSolver>
allocate_scaled_poisson_petsc_krylov_solver(const std::string& object_name,
                                            Pointer<Database> input_db,
                                            const std::string& default_options_prefix)
{
    Pointer<PETScKrylovPoissonSolver> solver =
        new PETScKrylovPoissonSolver(object_name, input_db, default_options_prefix);

    Pointer<Database> operator_db = nullptr;
    if (input_db && input_db->isDatabase("laplace_op"))
    {
        operator_db = input_db->getDatabase("laplace_op");
    }

    solver->setOperator(new CCScaledLaplaceOperator(object_name + "::CCScaledLaplaceOperator", operator_db));

    return solver;
} // allocate_scaled_poisson_petsc_krylov_solver

void
register_scaled_poisson_solvers()
{
    CCPoissonSolverManager* manager = CCPoissonSolverManager::getManager();

    manager->registerSolverFactoryFunction("SCALED_PETSC_LEVEL_SOLVER",
                                           CCScaledPoissonPETScLevelSolver::allocate_solver);

    manager->registerSolverFactoryFunction("SCALED_POINT_RELAXATION_FAC_PRECONDITIONER",
                                           CCScaledPoissonPointRelaxationFACOperator::allocate_solver);

    manager->registerSolverFactoryFunction("SCALED_PETSC_KRYLOV_SOLVER", allocate_scaled_poisson_petsc_krylov_solver);

    return;
} // register_scaled_poisson_solvers

CompositeNorms
compute_composite_norms(Pointer<PatchHierarchy<NDIM>> hierarchy, const int data_idx, const int weight_idx)
{
    CompositeNorms norms;

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> data = patch->getPatchData(data_idx);
            Pointer<CellData<NDIM, double>> weights = patch->getPatchData(weight_idx);

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const double weight = (*weights)(ci());

                // Do not count coarse cells covered by finer levels.
                if (weight <= 0.0) continue;

                const double value = (*data)(ci());

                norms.max_norm = std::max(norms.max_norm, std::abs(value));
                norms.L2_norm += weight * value * value;
                norms.L1_norm += weight * std::abs(value);
            }
        }
    }

    norms.max_norm = IBTK_MPI::maxReduction(norms.max_norm);
    norms.L2_norm = std::sqrt(IBTK_MPI::sumReduction(norms.L2_norm));
    norms.L1_norm = IBTK_MPI::sumReduction(norms.L1_norm);

    return norms;
} // compute_composite_norms
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_scaled_poisson_fac.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        const double data_time = 0.0;
        const bool run_discrete_test = input_db->getBool("run_discrete_test");

        // Register all factories before allocating the outer solver or FAC strategy.
        register_scaled_poisson_solvers();

        // Construct hierarchy infrastructure.
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

        // Construct all requested hierarchy levels.
        gridding_algorithm->makeCoarsestLevel(hierarchy, data_time);

        const int tag_buffer = 1;
        int level_number = 0;

        while (gridding_algorithm->levelCanBeRefined(level_number))
        {
            gridding_algorithm->makeFinerLevel(hierarchy, data_time, data_time, tag_buffer);

            if (!hierarchy->finerLevelExists(level_number)) break;

            ++level_number;
        }

        const int finest_ln = hierarchy->getFinestLevelNumber();

        // if (finest_ln < 1)
        // {
        //     TBOX_ERROR("This test requires at least two hierarchy levels.\n");
        // }

        // Register solution, forcing, diagnostics, and coefficient variables.
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

        const std::vector<int> data_indices = { u_idx, f_idx, exact_idx, error_idx, residual_idx, C_idx, S_idx, D_idx };

        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);

            for (const int idx : data_indices) level->allocatePatchData(idx, data_time);
        }

        // Composite-grid weights are zero on covered coarse cells.
        HierarchyMathOps hier_math_ops("test_math_ops", hierarchy, 0, finest_ln);
        const int weight_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> u_vec("u", hierarchy, 0, finest_ln);
        SAMRAIVectorReal<NDIM, double> f_vec("f", hierarchy, 0, finest_ln);
        SAMRAIVectorReal<NDIM, double> exact_vec("exact", hierarchy, 0, finest_ln);
        SAMRAIVectorReal<NDIM, double> error_vec("error", hierarchy, 0, finest_ln);
        SAMRAIVectorReal<NDIM, double> residual_vec("residual", hierarchy, 0, finest_ln);

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

        HierarchyCellDataOpsReal<NDIM, double> data_ops(hierarchy, 0, finest_ln);

        // Read all analytic fields from the input file.
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

        muParserRobinBcCoefs bc_coefs(
            "bc_coefs", app_initializer->getComponentDatabase("BoundaryConditions"), grid_geometry);

        PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setCPatchDataId(C_idx);
        poisson_spec.setDPatchDataId(D_idx);

        // Allocate the outer solver and its preconditioner through the manager.
        Pointer<PoissonSolver> poisson_solver =
            CCPoissonSolverManager::getManager()->allocateSolver(input_db->getString("solver_type"),
                                                                 "scaled_krylov",
                                                                 app_initializer->getComponentDatabase("KrylovSolver"),
                                                                 "scaled_poisson_",
                                                                 input_db->getString("precond_type"),
                                                                 "scaled_fac",
                                                                 app_initializer->getComponentDatabase("FAC"),
                                                                 "fac");

        // Use non-owning raw casts to inspect the manager-created components.
        auto* krylov_solver = dynamic_cast<KrylovLinearSolver*>(poisson_solver.getPointer());

        if (!krylov_solver)
        {
            TBOX_ERROR("The selected outer solver is not a KrylovLinearSolver.\n");
        }

        auto* scaled_op = dynamic_cast<CCScaledLaplaceOperator*>(krylov_solver->getOperator().getPointer());

        if (!scaled_op)
        {
            TBOX_ERROR("The outer solver does not contain CCScaledLaplaceOperator.\n");
        }

        auto* fac_preconditioner =
            dynamic_cast<PoissonFACPreconditioner*>(krylov_solver->getPreconditioner().getPointer());

        if (!fac_preconditioner)
        {
            TBOX_ERROR("The selected preconditioner is not a PoissonFACPreconditioner.\n");
        }

        auto* fac_strategy = dynamic_cast<CCScaledPoissonPointRelaxationFACOperator*>(
            fac_preconditioner->getFACPreconditionerStrategy().getPointer());

        if (!fac_strategy)
        {
            TBOX_ERROR("The FAC preconditioner does not contain the scaled point-relaxation strategy.\n");
        }

        // Standard configuration propagates through the PoissonSolver interface.
        poisson_solver->setPoissonSpecifications(poisson_spec);
        poisson_solver->setPhysicalBcCoef(&bc_coefs);
        poisson_solver->setHomogeneousBc(false);
        poisson_solver->setSolutionTime(data_time);
        poisson_solver->setTimeInterval(data_time, data_time + 1.0);

        krylov_solver->setInitialGuessNonzero(false);

        // S is additional configuration not represented by PoissonSpecifications.
        scaled_op->setScalePatchDataIndex(S_idx);
        fac_strategy->setScalePatchDataIndex(S_idx);

        // Initializes the operator, FAC strategy, and scaled coarse solver.
        poisson_solver->initializeSolverState(u_vec, f_vec);

        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();

        if (visit_data_writer)
        {
            visit_data_writer->registerPlotQuantity("u_solution", "SCALAR", u_idx);
            visit_data_writer->registerPlotQuantity("u_exact", "SCALAR", exact_idx);
            visit_data_writer->registerPlotQuantity("u_error", "SCALAR", error_idx);
            visit_data_writer->registerPlotQuantity("residual", "SCALAR", residual_idx);
            visit_data_writer->registerPlotQuantity("C", "SCALAR", C_idx);
            visit_data_writer->registerPlotQuantity("S", "SCALAR", S_idx);
        }

        pout << std::scientific << std::setprecision(12);
        pout << "\nManager-allocated scaled Poisson Krylov/FAC test\n";
        pout << "number of hierarchy levels = " << finest_ln + 1 << "\n";

        auto solve_and_report = [&](const std::string& label)
        {
            u_vec.setToScalar(0.0);

            poisson_solver->solveSystem(u_vec, f_vec);

            data_ops.subtract(error_idx, exact_idx, u_idx);

            // Evaluate the physical operator with its nonhomogeneous Robin data.
            scaled_op->setHomogeneousBc(false);
            scaled_op->apply(u_vec, residual_vec);
            data_ops.subtract(residual_idx, f_idx, residual_idx);

            const CompositeNorms error_norms = compute_composite_norms(hierarchy, error_idx, weight_idx);
            const CompositeNorms residual_norms = compute_composite_norms(hierarchy, residual_idx, weight_idx);
            const CompositeNorms exact_norms = compute_composite_norms(hierarchy, exact_idx, weight_idx);

            pout << "\n" << label << "\n";
            pout << "|error|_oo = " << error_norms.max_norm << "\n";
            pout << "|error|_2   = " << error_norms.L2_norm << "\n";
            pout << "|error|_1   = " << error_norms.L1_norm << "\n";

            if (exact_norms.L2_norm > 0.0)
            {
                pout << "relative L2 error = " << error_norms.L2_norm / exact_norms.L2_norm << "\n";
            }

            pout << "|residual|_oo = " << residual_norms.max_norm << "\n";
            pout << "|residual|_2   = " << residual_norms.L2_norm << "\n";
            pout << "|residual|_1   = " << residual_norms.L1_norm << "\n";
        };

        // Continuous manufactured forcing: use this error for refinement studies.
        solve_and_report("CONTINUOUS MANUFACTURED SOLUTION");

        // Output the continuous-test fields before the second test overwrites them.
        if (visit_data_writer)
        {
            visit_data_writer->writePlotData(hierarchy, 0, data_time);
        }

        if (run_discrete_test)
        {
            // Manufacture a discrete composite-grid RHS using the exact grid data.
            scaled_op->setHomogeneousBc(false);
            scaled_op->apply(exact_vec, f_vec);

            solve_and_report("DISCRETE MANUFACTURED SOLUTION");
        }

        // Also deallocates the attached operator, FAC state, and coarse solver.
        poisson_solver->deallocateSolverState();
    }

    return 0;
}
