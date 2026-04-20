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

#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/StaggeredStokesIBJacobianFACPreconditioner.h>
#include <ibamr/StaggeredStokesIBJacobianOperator.h>
#include <ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesIBOperator.h>
#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/LData.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <IntVector.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SAMRAI_config.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
struct StructureSpec
{
    int num_curve_points = 64;
    double ds = 1.0 / 64.0;
    double x_center = 0.5;
    double y_center = 0.5;
    double x_radius = 0.2;
    double y_radius = 0.2;
    double spring_stiffness = 2.0e2;
    int finest_ln = 0;
};

void
generate_structure(const unsigned int& strct_num,
                   const int& ln,
                   int& num_vertices,
                   std::vector<IBTK::Point>& vertex_posn,
                   void* ctx)
{
    auto* spec = static_cast<StructureSpec*>(ctx);
    if (!spec)
    {
        TBOX_ERROR("generate_structure(): missing structure specification context\n");
    }

    if (ln != spec->finest_ln || strct_num != 0)
    {
        num_vertices = 0;
        vertex_posn.resize(0);
        return;
    }

    num_vertices = spec->num_curve_points;
    vertex_posn.resize(num_vertices);
    for (int k = 0; k < num_vertices; ++k)
    {
        const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_vertices);
        vertex_posn[k](0) = spec->x_center + spec->x_radius * std::cos(theta);
        vertex_posn[k](1) = spec->y_center + spec->y_radius * std::sin(theta);
    }
}

void
generate_springs(
    const unsigned int& strct_num,
    const int& ln,
    std::multimap<int, IBRedundantInitializer::Edge>& spring_map,
    std::map<IBRedundantInitializer::Edge, IBRedundantInitializer::SpringSpec, IBRedundantInitializer::EdgeComp>&
        spring_spec,
    void* ctx)
{
    auto* spec = static_cast<StructureSpec*>(ctx);
    if (!spec)
    {
        TBOX_ERROR("generate_springs(): missing structure specification context\n");
    }
    if (ln != spec->finest_ln || strct_num != 0) return;

    for (int k = 0; k < spec->num_curve_points; ++k)
    {
        IBRedundantInitializer::Edge edge = { k, (k + 1) % spec->num_curve_points };
        if (edge.first > edge.second) std::swap(edge.first, edge.second);
        spring_map.insert(std::make_pair(edge.first, edge));

        IBRedundantInitializer::SpringSpec spec_data;
        spec_data.force_fcn_idx = 0;
        spec_data.parameters.resize(2);
        spec_data.parameters[0] = spec->spring_stiffness;
        spec_data.parameters[1] = 0.0;
        spring_spec.insert(std::make_pair(edge, spec_data));
    }
}

} // namespace

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line. For non-restarted case, command line is:      *
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
    // Initialize IBAMR and dependent libraries.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

#ifndef IBTK_HAVE_SILO
    // Suppress warnings caused by running without Silo.
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
#endif

    // Cleanup dynamically allocated objects prior to shutdown.
    {
        // Parse command-line options and the input database.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        const double current_time = 0.0;
        const double dt = input_db->getDoubleWithDefault("DT", 0.005);
        const double rho = input_db->getDoubleWithDefault("RHO", 1.0);
        const double mu = input_db->getDoubleWithDefault("MU", 1.0);
        const double new_time = current_time + dt;
        const int num_linear_solves = input_db->getIntegerWithDefault("NUM_LINEAR_SOLVES", 1);
        const bool use_matrix_based_saj = input_db->getBoolWithDefault("USE_MATRIX_BASED_SAJ", true);

        StructureSpec structure_spec;
        structure_spec.ds = input_db->getDoubleWithDefault("DS", 1.0 / 64.0);
        structure_spec.x_center = input_db->getDoubleWithDefault("X_CENTER", 0.5);
        structure_spec.y_center = input_db->getDoubleWithDefault("Y_CENTER", 0.5);
        structure_spec.x_radius = input_db->getDoubleWithDefault("X_RADIUS", 0.2);
        structure_spec.y_radius = input_db->getDoubleWithDefault("Y_RADIUS", 0.2);
        structure_spec.spring_stiffness = input_db->getDoubleWithDefault("SPRING_STIFFNESS", 2.0e2);
        if (!(structure_spec.ds > 0.0)) TBOX_ERROR("DS must be positive\n");
        if (!(structure_spec.x_radius > 0.0) || !(structure_spec.y_radius > 0.0))
        {
            TBOX_ERROR("X_RADIUS and Y_RADIUS must be positive\n");
        }

        const double a = structure_spec.x_radius;
        const double b = structure_spec.y_radius;
        const double h = std::pow((a - b) / (a + b), 2.0);
        const double circumference =
            M_PI * (a + b) * (1.0 + 3.0 * h / (10.0 + std::sqrt(std::max(0.0, 4.0 - 3.0 * h))));
        structure_spec.num_curve_points = std::max(3, static_cast<int>(circumference / structure_spec.ds));

        // Construct the core IBAMR algorithm objects.
        Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
        ib_method_ops->setUseFixedLEOperators(true);

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               ib_method_ops,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        structure_spec.finest_ln = input_db->getIntegerWithDefault("MAX_LEVELS", 1) - 1;

        Pointer<IBRedundantInitializer> ib_initializer = new IBRedundantInitializer(
            "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
        ib_initializer->setStructureNamesOnLevel(structure_spec.finest_ln, { "curve2d" });
        ib_initializer->registerInitStructureFunction(generate_structure, &structure_spec);
        ib_initializer->registerInitSpringDataFunction(generate_springs, &structure_spec);
        ib_method_ops->registerLInitStrategy(ib_initializer);

        Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
        ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);

        // Build the patch hierarchy.
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, current_time);
        int tag_buffer = input_db->getIntegerWithDefault("TAG_BUFFER", 1);
        int level_number = 0;
        bool done = false;
        while (!done && gridding_algorithm->levelCanBeRefined(level_number))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, current_time, true, tag_buffer);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        // Register Eulerian variables and contexts.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> current_ctx = var_db->getContext("current_ctx");
        Pointer<VariableContext> scratch_ctx = var_db->getContext("scratch_ctx");
        Pointer<VariableContext> solver_ctx = var_db->getContext("solver_ctx");

        Pointer<SideVariable<NDIM, double>> f_var = new SideVariable<NDIM, double>("f_var");
        Pointer<CellVariable<NDIM, double>> g_var = new CellVariable<NDIM, double>("g_var");
        Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p_var");
        Pointer<CellVariable<NDIM, int>> p_dof_index_var = new CellVariable<NDIM, int>("p_dof_index");
        Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u_var");
        Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("u_dof_index");

        const IntVector<NDIM> ib_ghosts = ib_method_ops->getMinimumGhostCellWidth();
        const IntVector<NDIM> one_ghost = IntVector<NDIM>(1);
        const IntVector<NDIM> no_ghosts = IntVector<NDIM>(0);

        const int u_current_idx = var_db->registerVariableAndContext(u_var, current_ctx, ib_ghosts);
        const int u_sol_idx = var_db->registerVariableAndContext(u_var, solver_ctx, one_ghost);
        const int f_rhs_idx = var_db->registerVariableAndContext(f_var, solver_ctx, one_ghost);
        const int p_sol_idx = var_db->registerVariableAndContext(p_var, solver_ctx, one_ghost);
        const int g_rhs_idx = var_db->registerVariableAndContext(g_var, solver_ctx, one_ghost);
        const int u_scratch_idx = var_db->registerVariableAndContext(u_var, scratch_ctx, ib_ghosts);
        const int f_scratch_idx = var_db->registerVariableAndContext(f_var, scratch_ctx, ib_ghosts);
        const int u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, scratch_ctx, ib_ghosts);
        const int p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, scratch_ctx, no_ghosts);

        const std::vector<int> allocated_patch_data_indices = {
            u_current_idx, u_scratch_idx, f_scratch_idx, u_dof_index_idx, p_dof_index_idx
        };
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (const int data_idx : allocated_patch_data_indices) level->allocatePatchData(data_idx, current_time);
        }

        Pointer<HierarchySideDataOpsReal<NDIM, double>> hier_velocity_data_ops =
            new HierarchySideDataOpsReal<NDIM, double>(patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        Pointer<HierarchyCellDataOpsReal<NDIM, double>> hier_pressure_data_ops =
            new HierarchyCellDataOpsReal<NDIM, double>(patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        // Initialize velocity and force data.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            u_init->setDataOnPatchHierarchy(u_current_idx, u_var, patch_hierarchy, current_time);
        }
        else
        {
            hier_velocity_data_ops->setToScalar(u_current_idx, 0.0, false);
        }
        hier_velocity_data_ops->setToScalar(u_scratch_idx, 0.0, false);
        hier_velocity_data_ops->setToScalar(f_scratch_idx, 0.0, false);

        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        std::vector<Pointer<CoarsenSchedule<NDIM>>> u_synch_scheds(finest_ln + 1);
        std::vector<Pointer<RefineSchedule<NDIM>>> u_ghost_fill_scheds(finest_ln + 1);
        std::vector<Pointer<RefineSchedule<NDIM>>> f_prolongation_scheds(finest_ln + 1);

        ib_method_ops->initializePatchHierarchy(patch_hierarchy,
                                                gridding_algorithm,
                                                u_current_idx,
                                                u_synch_scheds,
                                                u_ghost_fill_scheds,
                                                0,
                                                current_time,
                                                true);
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();

        // Assemble fixed Lagrangian-Eulerian operators for the implicit solve.
        ib_method_ops->preprocessIntegrateData(current_time, new_time, 1);
        ib_method_ops->updateFixedLEOperators();

        std::vector<std::vector<int>> num_dofs_per_proc(finest_ln + 1);
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
                num_dofs_per_proc[ln], u_dof_index_idx, p_dof_index_idx, level);
        }

        Mat A = nullptr;
        ib_method_ops->constructLagrangianForceJacobian(A, MATAIJ, new_time);
        Mat J = nullptr;
        ib_method_ops->constructInterpOp(J,
                                         PETScMatUtilities::ib_4_delta_fcn,
                                         PETScMatUtilities::ib_4_delta_stencil,
                                         num_dofs_per_proc[finest_ln],
                                         u_dof_index_idx,
                                         new_time);

        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        Pointer<SAMRAIVectorReal<NDIM, double>> eul_sol_vec =
            new SAMRAIVectorReal<NDIM, double>("eul_sol_vec", patch_hierarchy, 0, finest_ln);
        eul_sol_vec->addComponent(u_var, u_sol_idx, wgt_sc_idx, hier_velocity_data_ops);
        eul_sol_vec->addComponent(p_var, p_sol_idx, wgt_cc_idx, hier_pressure_data_ops);
        eul_sol_vec->allocateVectorData();

        Pointer<SAMRAIVectorReal<NDIM, double>> eul_rhs_vec =
            new SAMRAIVectorReal<NDIM, double>("eul_rhs_vec", patch_hierarchy, 0, finest_ln);
        eul_rhs_vec->addComponent(f_var, f_rhs_idx, wgt_sc_idx, hier_velocity_data_ops);
        eul_rhs_vec->addComponent(g_var, g_rhs_idx, wgt_cc_idx, hier_pressure_data_ops);
        eul_rhs_vec->allocateVectorData();

        hier_velocity_data_ops->copyData(u_sol_idx, u_current_idx);
        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            p_init->setDataOnPatchHierarchy(p_sol_idx, p_var, patch_hierarchy, current_time);
        }
        else
        {
            hier_pressure_data_ops->setToScalar(p_sol_idx, 0.0, false);
        }
        hier_velocity_data_ops->setToScalar(f_rhs_idx, 0.0, false);
        hier_pressure_data_ops->setToScalar(g_rhs_idx, 0.0, false);

        // Configure the staggered Stokes operator and boundary conditions.
        const double lambda = 0.0;
        PoissonSpecifications U_problem_coefs("stokes_ib_implicit_example::U_problem_coefs");
        U_problem_coefs.setCConstant(rho / dt + lambda);
        U_problem_coefs.setDConstant(-mu);

        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
        if (periodic_shift.min() <= 0)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                const std::string bc_name = "u_bc_coefs_" + std::to_string(d);
                const std::string db_name = "VelocityBcCoefs_" + std::to_string(d);
                u_bc_coefs[d] =
                    new muParserRobinBcCoefs(bc_name, app_initializer->getComponentDatabase(db_name), grid_geometry);
            }
        }

        Pointer<StaggeredStokesOperator> stokes_op =
            new StaggeredStokesOperator("stokes_ib_implicit_example::stokes_op", false);
        stokes_op->setVelocityPoissonSpecifications(U_problem_coefs);
        stokes_op->setPhysicalBcCoefs(u_bc_coefs, nullptr);
        stokes_op->setTimeInterval(current_time, new_time);
        stokes_op->setSolutionTime(new_time);

        StaggeredStokesIBOperator::Context ctx;
        ctx.ib_implicit_ops = ib_method_ops;
        ctx.stokes_op = stokes_op;
        ctx.u_phys_bdry_op = nullptr;
        ctx.hier_velocity_data_ops = hier_velocity_data_ops;
        ctx.u_synch_scheds = u_synch_scheds;
        ctx.u_ghost_fill_scheds = u_ghost_fill_scheds;
        ctx.f_prolongation_scheds = f_prolongation_scheds;
        ctx.patch_level = patch_hierarchy->getPatchLevel(finest_ln);
        ctx.u_idx = u_scratch_idx;
        ctx.f_idx = f_scratch_idx;
        ctx.u_current_idx = u_current_idx;
        ctx.u_dof_index_idx = u_dof_index_idx;
        ctx.p_dof_index_idx = p_dof_index_idx;
        ctx.use_fixed_le_operators = true;
        const std::string ib_time_stepping_type =
            input_db->getStringWithDefault("IB_TIME_STEPPING_TYPE", "MIDPOINT_RULE");
        ctx.time_stepping_type = IBAMR::string_to_enum<TimeSteppingType>(ib_time_stepping_type);

        Pointer<StaggeredStokesIBJacobianOperator> jac_op =
            new StaggeredStokesIBJacobianOperator("stokes_ib_implicit_example::jacobian_op");
        jac_op->setOperatorContext(ctx);
        jac_op->setTimeInterval(current_time, new_time);
        jac_op->setSolutionTime(new_time);
        jac_op->initializeOperatorState(*eul_sol_vec, *eul_rhs_vec);
        jac_op->formJacobian(*eul_sol_vec);

        Pointer<Database> stokes_ib_precond_db =
            input_db->isDatabase("stokes_ib_precond_db") ? input_db->getDatabase("stokes_ib_precond_db") : nullptr;
        if (!stokes_ib_precond_db)
        {
            TBOX_ERROR("missing stokes_ib_precond_db in input file\n");
        }
        const bool has_pressure_nullspace = stokes_ib_precond_db->getBoolWithDefault("has_pressure_nullspace", true);

        Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_op = new StaggeredStokesIBLevelRelaxationFACOperator(
            "stokes_ib_implicit_example::fac_op", stokes_ib_precond_db, "stokes_ib_pc_");
        Pointer<StaggeredStokesIBJacobianFACPreconditioner> fac_pc = new StaggeredStokesIBJacobianFACPreconditioner(
            "stokes_ib_implicit_example::fac_pc", fac_op, stokes_ib_precond_db, "stokes_ib_pc_");
        Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();

        fac_pc->setVelocityPoissonSpecifications(U_problem_coefs);
        fac_pc->setPhysicalBcCoefs(u_bc_coefs, nullptr);
        fac_pc->setPhysicalBoundaryHelper(bc_helper);
        fac_pc->setTimeInterval(current_time, new_time);
        fac_pc->setSolutionTime(new_time);
        fac_pc->setHomogeneousBc(true);
        fac_pc->setComponentsHaveNullSpace(false, has_pressure_nullspace);
        fac_pc->setIBTimeSteppingType(ctx.time_stepping_type);
        fac_pc->setIBForceJacobian(A);
        fac_pc->setIBInterpOp(J);
        fac_pc->setIBImplicitStrategy(ib_method_ops);
        fac_pc->initializeSolverState(*eul_sol_vec, *eul_rhs_vec);

        if (use_matrix_based_saj)
        {
            Mat SAJ = fac_op->getEulerianElasticityLevelOp(finest_ln);
            jac_op->setIBCouplingJacobian(SAJ);
            pout << "Using matrix-based SAJ Jacobian apply path" << std::endl;
        }
        else
        {
            pout << "Using matrix-free Jacobian apply path" << std::endl;
        }

        Pointer<SAMRAIVectorReal<NDIM, double>> linear_rhs_seed_vec = eul_sol_vec->cloneVector("linear_rhs_seed_vec");
        linear_rhs_seed_vec->allocateVectorData();
        linear_rhs_seed_vec->setToScalar(0.0);

        if (input_db->isDatabase("LinearSolveRhsSeedVelocity"))
        {
            Pointer<CartGridFunction> rhs_seed_u_fcn = new muParserCartGridFunction(
                "rhs_seed_u_fcn", app_initializer->getComponentDatabase("LinearSolveRhsSeedVelocity"), grid_geometry);
            rhs_seed_u_fcn->setDataOnPatchHierarchy(
                linear_rhs_seed_vec->getComponentDescriptorIndex(0), u_var, patch_hierarchy, current_time);
        }
        else
        {
            hier_velocity_data_ops->setToScalar(linear_rhs_seed_vec->getComponentDescriptorIndex(0), 1.0, false);
        }

        if (input_db->isDatabase("LinearSolveRhsSeedPressure"))
        {
            Pointer<CartGridFunction> rhs_seed_p_fcn = new muParserCartGridFunction(
                "rhs_seed_p_fcn", app_initializer->getComponentDatabase("LinearSolveRhsSeedPressure"), grid_geometry);
            rhs_seed_p_fcn->setDataOnPatchHierarchy(
                linear_rhs_seed_vec->getComponentDescriptorIndex(1), p_var, patch_hierarchy, current_time);
        }

        Pointer<SAMRAIVectorReal<NDIM, double>> linear_rhs = eul_rhs_vec->cloneVector("linear_rhs");
        linear_rhs->allocateVectorData();
        linear_rhs->setToScalar(0.0);
        jac_op->apply(*linear_rhs_seed_vec, *linear_rhs);

        Pointer<PETScKrylovLinearSolver> linear_solver =
            new PETScKrylovLinearSolver("stokes_ib_implicit_example::linear_solver", nullptr, "ib_");
        linear_solver->setOperator(jac_op);
        linear_solver->setPreconditioner(fac_pc);
        linear_solver->setTimeInterval(current_time, new_time);
        linear_solver->setSolutionTime(new_time);
        linear_solver->setInitialGuessNonzero(false);
        Pointer<SAMRAIVectorReal<NDIM, double>> pressure_nullspace_vec;
        if (has_pressure_nullspace)
        {
            pressure_nullspace_vec = eul_sol_vec->cloneVector("pressure_nullspace_vec");
            pressure_nullspace_vec->allocateVectorData();
            hier_velocity_data_ops->setToScalar(pressure_nullspace_vec->getComponentDescriptorIndex(0), 0.0, false);
            hier_pressure_data_ops->setToScalar(pressure_nullspace_vec->getComponentDescriptorIndex(1), 1.0, false);
            linear_solver->setNullSpace(
                false, std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>>(1, pressure_nullspace_vec));
        }

        Pointer<SAMRAIVectorReal<NDIM, double>> linear_sol = eul_sol_vec->cloneVector("linear_sol");
        linear_sol->allocateVectorData();

        bool all_solves_success = true;
        // Repeatedly solve the linearized implicit Stokes-IB system.
        for (int solve_n = 0; solve_n < num_linear_solves; ++solve_n)
        {
            linear_sol->setToScalar(0.0);
            const bool linear_success = linear_solver->solveSystem(*linear_sol, *linear_rhs);
            all_solves_success = all_solves_success && linear_success;
            pout << "linear_solve_" << solve_n << "_success = " << static_cast<int>(linear_success) << std::endl;
            pout << "linear_solve_" << solve_n << "_iterations = " << linear_solver->getNumIterations() << std::endl;
            pout << "linear_solve_" << solve_n << "_residual_norm = " << linear_solver->getResidualNorm() << std::endl;
        }

        double rhs_u_norm = std::numeric_limits<double>::quiet_NaN();
        double rhs_p_norm = std::numeric_limits<double>::quiet_NaN();
        double sol_u_norm = std::numeric_limits<double>::quiet_NaN();
        double sol_p_norm = std::numeric_limits<double>::quiet_NaN();

        rhs_u_norm = hier_velocity_data_ops->L2Norm(linear_rhs->getComponentDescriptorIndex(0), wgt_sc_idx);
        rhs_p_norm = hier_pressure_data_ops->L2Norm(linear_rhs->getComponentDescriptorIndex(1), wgt_cc_idx);
        sol_u_norm = hier_velocity_data_ops->L2Norm(linear_sol->getComponentDescriptorIndex(0), wgt_sc_idx);
        sol_p_norm = hier_pressure_data_ops->L2Norm(linear_sol->getComponentDescriptorIndex(1), wgt_cc_idx);

        pout << "rhs_velocity_l2_norm = " << rhs_u_norm << std::endl;
        pout << "rhs_pressure_l2_norm = " << rhs_p_norm << std::endl;
        pout << "solution_velocity_l2_norm = " << sol_u_norm << std::endl;
        pout << "solution_pressure_l2_norm = " << sol_p_norm << std::endl;

        linear_sol->deallocateVectorData();
        linear_rhs->deallocateVectorData();
        linear_rhs_seed_vec->deallocateVectorData();
        if (pressure_nullspace_vec) pressure_nullspace_vec->deallocateVectorData();

        fac_pc->deallocateSolverState();
        jac_op->deallocateOperatorState();

        ib_method_ops->postprocessIntegrateData(current_time, new_time, 1);

        eul_sol_vec->deallocateVectorData();
        eul_rhs_vec->deallocateVectorData();

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (const int data_idx : allocated_patch_data_indices)
            {
                if (level->checkAllocated(data_idx)) level->deallocatePatchData(data_idx);
            }
        }

        PetscErrorCode ierr = MatDestroy(&A);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&J);
        IBTK_CHKERRQ(ierr);

        // Deallocate dynamically allocated boundary-condition objects.
        for (int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

        plog << "Input database:\n";
        input_db->printClassData(plog);

        if (!all_solves_success)
        {
            TBOX_ERROR("one or more linear solves failed\n");
        }
        if (!std::isfinite(rhs_u_norm) || !std::isfinite(rhs_p_norm) || !std::isfinite(sol_u_norm) ||
            !std::isfinite(sol_p_norm))
        {
            TBOX_ERROR("linear solve produced non-finite norms\n");
        }
    }

    return 0;
}
