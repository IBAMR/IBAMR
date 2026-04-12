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
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
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
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
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
    double x_radius = 0.23;
    double y_radius = 0.25 * 0.25 / 0.23;
    double spring_stiffness = 1.0e2;
    int finest_ln = 0;
};

struct CaseResult
{
    int depth = 1;
    double spring_stiffness = 0.0;
    double op_rel_residual = std::numeric_limits<double>::quiet_NaN();
    double prec_rel_residual = std::numeric_limits<double>::quiet_NaN();
    double prec_op_rel_residual = std::numeric_limits<double>::quiet_NaN();
    bool preconditioner_success = false;
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

void
load_petsc_options_file(Pointer<Database> input_db, const int argc, char* argv[])
{
    if (!input_db->keyExists("PETSC_OPTIONS_FILE")) return;

    const std::string petsc_options_file = input_db->getString("PETSC_OPTIONS_FILE");
    std::string resolved_options_file = petsc_options_file;
    auto path_exists = [](const std::string& path) -> bool
    {
        std::ifstream stream(path.c_str());
        return stream.good();
    };
    if (!path_exists(resolved_options_file) && argc > 1)
    {
        const std::string input_filename = argv[1];
        const std::size_t last_sep = input_filename.find_last_of("/\\");
        if (last_sep != std::string::npos)
        {
            const std::string candidate = input_filename.substr(0, last_sep + 1) + petsc_options_file;
            if (path_exists(candidate)) resolved_options_file = candidate;
        }
    }
    if (!path_exists(resolved_options_file)) resolved_options_file = std::string(SOURCE_DIR) + petsc_options_file;

    std::ifstream options_stream(resolved_options_file.c_str());
    if (!options_stream.good())
    {
        TBOX_ERROR("could not open PETSc options file: " << petsc_options_file << "\n");
    }
    options_stream.close();

    const PetscErrorCode ierr =
        PetscOptionsInsertFile(PETSC_COMM_WORLD, nullptr, resolved_options_file.c_str(), PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
}

void
set_nontrivial_probe_vector(const int u_idx, const int p_idx, Pointer<PatchHierarchy<NDIM>> patch_hierarchy)
{
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();
            const double* const x_lower = patch_geom->getXLower();
            const SAMRAI::hier::Index<NDIM>& lower = patch_box.lower();

            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                for (Box<NDIM>::Iterator b(side_box); b; b++)
                {
                    const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
                    double x[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        const double offset = (d == axis) ? 0.0 : 0.5;
                        x[d] = x_lower[d] + (static_cast<double>(i_s(d) - lower(d)) + offset) * dx[d];
                    }
#if (NDIM == 2)
                    const double val = (axis == 0) ?
                                           (std::sin(2.0 * M_PI * x[1]) + 0.125 * std::cos(4.0 * M_PI * x[0])) :
                                           (-std::sin(2.0 * M_PI * x[0]) + 0.125 * std::cos(4.0 * M_PI * x[1]));
#else
                    const double val = (axis == 0) ? std::sin(2.0 * M_PI * x[1]) :
                                       ((axis == 1) ? std::sin(2.0 * M_PI * x[2]) : std::sin(2.0 * M_PI * x[0]));
#endif
                    (*u_data)(i_s) = val;
                }
            }

            for (Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                const CellIndex<NDIM> i_c(b());
                const double x = x_lower[0] + (static_cast<double>(i_c(0) - lower(0)) + 0.5) * dx[0];
                const double y = x_lower[1] + (static_cast<double>(i_c(1) - lower(1)) + 0.5) * dx[1];
                (*p_data)(i_c) = std::cos(2.0 * M_PI * x) + std::sin(2.0 * M_PI * y) +
                                 0.1 * std::sin(2.0 * M_PI * (x + y));
            }
        }
    }
}

double
vector_relative_residual(const Pointer<SAMRAIVectorReal<NDIM, double>>& lhs,
                         const Pointer<SAMRAIVectorReal<NDIM, double>>& rhs,
                         const Pointer<SAMRAIVectorReal<NDIM, double>>& scratch)
{
    scratch->subtract(lhs, rhs);
    const double rhs_norm = rhs->L2Norm();
    return scratch->L2Norm() / std::max(rhs_norm, 1.0e-14);
}

CaseResult
run_case(Pointer<AppInitializer> app_initializer,
         Pointer<Database> input_db,
         const int depth,
         const double spring_stiffness,
         const int case_id,
         int& test_failures)
{
    CaseResult result;
    result.depth = depth;
    result.spring_stiffness = spring_stiffness;
    const std::string case_suffix = "_case_" + std::to_string(case_id);

    const double current_time = 0.0;
    const double dt = input_db->getDoubleWithDefault("DT", 0.005);
    const double rho = input_db->getDoubleWithDefault("RHO", 1.0);
    const double mu = input_db->getDoubleWithDefault("MU", 1.0);
    const double new_time = current_time + dt;
    const bool use_fixed_le_operators = input_db->getBoolWithDefault("USE_FIXED_LE_OPERATORS", true);

    const int N = input_db->getIntegerWithDefault("N", 8);
    const int ref_ratio = input_db->getIntegerWithDefault("REF_RATIO", 2);
    const double L = input_db->getDoubleWithDefault("L", 1.0);

    int finest_grid_cells = N;
    for (int ln = 1; ln < depth; ++ln)
    {
        finest_grid_cells *= ref_ratio;
    }
    const double finest_dx = L / static_cast<double>(finest_grid_cells);

    StructureSpec structure_spec;
    structure_spec.ds = finest_dx;
    structure_spec.x_center = input_db->getDoubleWithDefault("X_CENTER", 0.5);
    structure_spec.y_center = input_db->getDoubleWithDefault("Y_CENTER", 0.5);
    structure_spec.x_radius = input_db->getDoubleWithDefault("X_RADIUS", 0.23);
    structure_spec.y_radius = input_db->getDoubleWithDefault("Y_RADIUS", 0.25 * 0.25 / 0.23);
    structure_spec.spring_stiffness = spring_stiffness;
    structure_spec.finest_ln = depth - 1;

    if (!(structure_spec.ds > 0.0))
    {
        TBOX_ERROR("DS must be positive\n");
    }
    if (!(structure_spec.x_radius > 0.0) || !(structure_spec.y_radius > 0.0))
    {
        TBOX_ERROR("X_RADIUS and Y_RADIUS must be positive\n");
    }

    const double a = structure_spec.x_radius;
    const double b = structure_spec.y_radius;
    const double h = std::pow((a - b) / (a + b), 2.0);
    const double circumference =
        M_PI * (a + b) * (1.0 + 3.0 * h / (10.0 + std::sqrt(std::max(0.0, 4.0 - 3.0 * h))));
    structure_spec.num_curve_points = std::max(3, static_cast<int>(std::llround(circumference / structure_spec.ds)));

    Pointer<IBMethod> ib_method_ops =
        new IBMethod("IBMethod" + case_suffix, app_initializer->getComponentDatabase("IBMethod"));
    ib_method_ops->setUseFixedLEOperators(use_fixed_le_operators);

    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry" + case_suffix, app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy =
        new PatchHierarchy<NDIM>("PatchHierarchy" + case_suffix, grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector =
        new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize" + case_suffix,
                                           ib_method_ops,
                                           app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> load_balancer =
        new LoadBalancer<NDIM>("LoadBalancer" + case_suffix, app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
        new GriddingAlgorithm<NDIM>("GriddingAlgorithm" + case_suffix,
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

    Pointer<IBRedundantInitializer> ib_initializer =
        new IBRedundantInitializer("IBRedundantInitializer" + case_suffix,
                                   app_initializer->getComponentDatabase("IBRedundantInitializer"));
    ib_initializer->setStructureNamesOnLevel(structure_spec.finest_ln, { "curve2d" });
    ib_initializer->registerInitStructureFunction(generate_structure, &structure_spec);
    ib_initializer->registerInitSpringDataFunction(generate_springs, &structure_spec);
    ib_method_ops->registerLInitStrategy(ib_initializer);

    Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
    ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);

    gridding_algorithm->makeCoarsestLevel(patch_hierarchy, current_time);
    const int target_finest_ln = depth - 1;
    const int tag_buffer = input_db->getIntegerWithDefault("TAG_BUFFER", 2);
    int level_number = 0;
    while (level_number < target_finest_ln && gridding_algorithm->levelCanBeRefined(level_number))
    {
        gridding_algorithm->makeFinerLevel(patch_hierarchy, current_time, true, tag_buffer);
        ++level_number;
    }
    if (patch_hierarchy->getFinestLevelNumber() != target_finest_ln)
    {
        TBOX_ERROR("Failed to build hierarchy with requested depth = " << depth << "\n");
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    static Pointer<VariableContext> current_ctx = nullptr;
    static Pointer<VariableContext> scratch_ctx = nullptr;
    static Pointer<VariableContext> solver_ctx = nullptr;
    static Pointer<SideVariable<NDIM, double>> f_var = nullptr;
    static Pointer<CellVariable<NDIM, double>> g_var = nullptr;
    static Pointer<CellVariable<NDIM, double>> p_var = nullptr;
    static Pointer<CellVariable<NDIM, int>> p_dof_index_var = nullptr;
    static Pointer<SideVariable<NDIM, double>> u_var = nullptr;
    static Pointer<SideVariable<NDIM, int>> u_dof_index_var = nullptr;
    static int u_current_idx = -1;
    static int u_sol_idx = -1;
    static int f_rhs_idx = -1;
    static int p_sol_idx = -1;
    static int g_rhs_idx = -1;
    static int u_scratch_idx = -1;
    static int f_scratch_idx = -1;
    static int u_dof_index_idx = -1;
    static int p_dof_index_idx = -1;
    static bool vars_initialized = false;

    if (!vars_initialized)
    {
        current_ctx = var_db->getContext("current_ctx");
        scratch_ctx = var_db->getContext("scratch_ctx");
        solver_ctx = var_db->getContext("solver_ctx");
        f_var = new SideVariable<NDIM, double>("f_var");
        g_var = new CellVariable<NDIM, double>("g_var");
        p_var = new CellVariable<NDIM, double>("p_var");
        p_dof_index_var = new CellVariable<NDIM, int>("p_dof_index");
        u_var = new SideVariable<NDIM, double>("u_var");
        u_dof_index_var = new SideVariable<NDIM, int>("u_dof_index");

        const IntVector<NDIM> ib_ghosts = ib_method_ops->getMinimumGhostCellWidth();
        const IntVector<NDIM> one_ghost = IntVector<NDIM>(1);
        const IntVector<NDIM> no_ghosts = IntVector<NDIM>(0);

        u_current_idx = var_db->registerVariableAndContext(u_var, current_ctx, ib_ghosts);
        u_sol_idx = var_db->registerVariableAndContext(u_var, solver_ctx, one_ghost);
        f_rhs_idx = var_db->registerVariableAndContext(f_var, solver_ctx, one_ghost);
        p_sol_idx = var_db->registerVariableAndContext(p_var, solver_ctx, one_ghost);
        g_rhs_idx = var_db->registerVariableAndContext(g_var, solver_ctx, one_ghost);
        u_scratch_idx = var_db->registerVariableAndContext(u_var, scratch_ctx, ib_ghosts);
        f_scratch_idx = var_db->registerVariableAndContext(f_var, scratch_ctx, ib_ghosts);
        u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, scratch_ctx, ib_ghosts);
        p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, scratch_ctx, no_ghosts);
        vars_initialized = true;
    }

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

    hier_velocity_data_ops->setToScalar(u_current_idx, 0.0, false);
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
                                     PETScMatUtilities::ib_4_interp_fcn,
                                     PETScMatUtilities::ib_4_interp_stencil,
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

    eul_sol_vec->setToScalar(0.0);
    eul_rhs_vec->setToScalar(0.0);

    const double lambda = 0.0;
    PoissonSpecifications U_problem_coefs("stokes_ib_operator_preconditioner_residual_sweep_01::U_problem_coefs");
    U_problem_coefs.setCConstant(rho / dt + lambda);
    U_problem_coefs.setDConstant(-mu);

    const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
    std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
    if (periodic_shift.min() <= 0)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);
            const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);
            u_bc_coefs[d] = new muParserRobinBcCoefs(
                bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
        }
    }

    Pointer<StaggeredStokesOperator> stokes_op =
        new StaggeredStokesOperator("stokes_ib_operator_preconditioner_residual_sweep_01::stokes_op", false);
    stokes_op->setVelocityPoissonSpecifications(U_problem_coefs);
    stokes_op->setPhysicalBcCoefs(u_bc_coefs, nullptr);
    stokes_op->setTimeInterval(current_time, new_time);
    stokes_op->setSolutionTime(new_time);

    StaggeredStokesIBOperatorContext ctx;
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
    ctx.use_fixed_le_operators = use_fixed_le_operators;
    ctx.time_stepping_type = IBAMR::string_to_enum<TimeSteppingType>(
        input_db->getStringWithDefault("IB_TIME_STEPPING", "BACKWARD_EULER"));

    Pointer<StaggeredStokesIBJacobianOperator> jac_op =
        new StaggeredStokesIBJacobianOperator("stokes_ib_operator_preconditioner_residual_sweep_01::jacobian_op");
    jac_op->setOperatorContext(ctx);
    jac_op->setTimeInterval(current_time, new_time);
    jac_op->setSolutionTime(new_time);
    jac_op->initializeOperatorState(*eul_sol_vec, *eul_rhs_vec);
    jac_op->formJacobian(*eul_sol_vec);

    Pointer<Database> stokes_ib_precond_db =
        input_db->isDatabase("stokes_ib_precond_db") ? input_db->getDatabase("stokes_ib_precond_db") : nullptr;

    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_op = new StaggeredStokesIBLevelRelaxationFACOperator(
        "stokes_ib_operator_preconditioner_residual_sweep_01::fac_op", stokes_ib_precond_db, "stokes_ib_pc_");
    Pointer<StaggeredStokesIBJacobianFACPreconditioner> fac_pc = new StaggeredStokesIBJacobianFACPreconditioner(
        "stokes_ib_operator_preconditioner_residual_sweep_01::fac_pc", fac_op, stokes_ib_precond_db, "stokes_ib_pc_");
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();

    fac_pc->setVelocityPoissonSpecifications(U_problem_coefs);
    fac_pc->setPhysicalBcCoefs(u_bc_coefs, nullptr);
    fac_pc->setPhysicalBoundaryHelper(bc_helper);
    fac_pc->setTimeInterval(current_time, new_time);
    fac_pc->setSolutionTime(new_time);
    fac_pc->setHomogeneousBc(true);
    fac_pc->setComponentsHaveNullSpace(false, true);
    fac_pc->setIBTimeSteppingType(ctx.time_stepping_type);
    fac_pc->setIBForceJacobian(A);
    fac_pc->setIBInterpOp(J);
    fac_pc->setIBImplicitStrategy(ib_method_ops);
    fac_pc->initializeSolverState(*eul_sol_vec, *eul_rhs_vec);

    Pointer<SAMRAIVectorReal<NDIM, double>> x_sol = eul_sol_vec->cloneVector("x_sol");
    x_sol->allocateVectorData();
    x_sol->setToScalar(0.0);
    set_nontrivial_probe_vector(x_sol->getComponentDescriptorIndex(0), x_sol->getComponentDescriptorIndex(1), patch_hierarchy);

    Pointer<SAMRAIVectorReal<NDIM, double>> x_rhs = eul_rhs_vec->cloneVector("x_rhs");
    x_rhs->allocateVectorData();
    x_rhs->setToScalar(0.0);
    hier_velocity_data_ops->copyData(x_rhs->getComponentDescriptorIndex(0), x_sol->getComponentDescriptorIndex(0));
    hier_pressure_data_ops->copyData(x_rhs->getComponentDescriptorIndex(1), x_sol->getComponentDescriptorIndex(1));

    Pointer<SAMRAIVectorReal<NDIM, double>> Ax_rhs = eul_rhs_vec->cloneVector("Ax_rhs");
    Ax_rhs->allocateVectorData();
    Ax_rhs->setToScalar(0.0);
    jac_op->apply(*x_sol, *Ax_rhs);

    Pointer<SAMRAIVectorReal<NDIM, double>> rhs_diff = eul_rhs_vec->cloneVector("rhs_diff");
    rhs_diff->allocateVectorData();
    rhs_diff->setToScalar(0.0);
    result.op_rel_residual = vector_relative_residual(Ax_rhs, x_rhs, rhs_diff);

    Pointer<SAMRAIVectorReal<NDIM, double>> Mx_sol = eul_sol_vec->cloneVector("Mx_sol");
    Mx_sol->allocateVectorData();
    Mx_sol->setToScalar(0.0);
    result.preconditioner_success = fac_pc->solveSystem(*Mx_sol, *x_rhs);
    if (!result.preconditioner_success)
    {
        ++test_failures;
        pout << "preconditioner solve failed for depth=" << depth << " K=" << spring_stiffness << std::endl;
    }

    Pointer<SAMRAIVectorReal<NDIM, double>> AMx_rhs = eul_rhs_vec->cloneVector("AMx_rhs");
    AMx_rhs->allocateVectorData();
    AMx_rhs->setToScalar(0.0);
    jac_op->apply(*Mx_sol, *AMx_rhs);
    result.prec_rel_residual = vector_relative_residual(AMx_rhs, x_rhs, rhs_diff);

    Pointer<SAMRAIVectorReal<NDIM, double>> MAx_sol = eul_sol_vec->cloneVector("MAx_sol");
    MAx_sol->allocateVectorData();
    MAx_sol->setToScalar(0.0);
    result.preconditioner_success = fac_pc->solveSystem(*MAx_sol, *Ax_rhs) && result.preconditioner_success;

    Pointer<SAMRAIVectorReal<NDIM, double>> sol_diff = eul_sol_vec->cloneVector("sol_diff");
    sol_diff->allocateVectorData();
    sol_diff->setToScalar(0.0);
    result.prec_op_rel_residual = vector_relative_residual(MAx_sol, x_sol, sol_diff);

    const bool finite_residuals = std::isfinite(result.op_rel_residual) && std::isfinite(result.prec_rel_residual) &&
                                  std::isfinite(result.prec_op_rel_residual);
    if (!finite_residuals)
    {
        ++test_failures;
        pout << "non-finite residual detected for depth=" << depth << " K=" << spring_stiffness << std::endl;
    }

    fac_pc->deallocateSolverState();
    jac_op->deallocateOperatorState();

    ib_method_ops->postprocessIntegrateData(current_time, new_time, 1);

    deallocate_vector_data(*eul_sol_vec);
    deallocate_vector_data(*eul_rhs_vec);
    free_vector_components(*eul_sol_vec);
    free_vector_components(*eul_rhs_vec);

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

    for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    return result;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

#ifndef IBTK_HAVE_SILO
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
#endif

    int test_failures = 0;

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        load_petsc_options_file(input_db, argc, argv);

        const int max_levels = input_db->getIntegerWithDefault("MAX_LEVELS", 1);
        const int depth = max_levels;
        const double k_val = input_db->getDoubleWithDefault("K", input_db->getDoubleWithDefault("SPRING_STIFFNESS", 1.0e2));
        const bool valid_depth = 0 < depth && depth <= max_levels;
        if (!valid_depth)
        {
            ++test_failures;
            pout << "invalid MAX_LEVELS: " << max_levels << std::endl;
        }
        else
        {
            const CaseResult result = run_case(app_initializer, input_db, depth, k_val, /*case_id*/ 0, test_failures);
            pout << std::scientific << std::setprecision(15)
                 << "depth=" << result.depth << " K=" << result.spring_stiffness
                 << " op_rel_residual=" << result.op_rel_residual
                 << " prec_rel_residual=" << result.prec_rel_residual
                 << " prec_op_rel_residual=" << result.prec_op_rel_residual
                 << " preconditioner_success=" << static_cast<int>(result.preconditioner_success) << std::endl;
        }

        pout << "test_failures = " << test_failures << std::endl;
    }

    return test_failures;
}
