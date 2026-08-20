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

/*
 * This test targets the implicit Stokes-IB solver components directly rather
 * than the full hierarchy integrator. It builds the nonlinear residual,
 * compares the Jacobian against a matrix-free finite-difference action, and
 * exercises the FAC-preconditioned linear solve on controlled single-level and
 * multilevel configurations.
 */

#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/StaggeredStokesIBJacobianFACPreconditioner.h>
#include <ibamr/StaggeredStokesIBJacobianOperator.h>
#include <ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesIBOperator.h>
#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/LData.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PETScMFFDJacobianOperator.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
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

bool
side_l2_norm_is_finite(const Pointer<HierarchySideDataOpsReal<NDIM, double>>& side_data_ops,
                       const int data_idx,
                       const int weight_idx,
                       double& l2_norm)
{
    l2_norm = side_data_ops->L2Norm(data_idx, weight_idx);
    return std::isfinite(l2_norm);
}

bool
cell_l2_norm_is_finite(const Pointer<HierarchyCellDataOpsReal<NDIM, double>>& cell_data_ops,
                       const int data_idx,
                       const int weight_idx,
                       double& l2_norm)
{
    l2_norm = cell_data_ops->L2Norm(data_idx, weight_idx);
    return std::isfinite(l2_norm);
}

void
set_divergence_free_probe_velocity(const int u_idx, Pointer<PatchHierarchy<NDIM>> patch_hierarchy)
{
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
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
                    const double val = (axis == 0) ? std::sin(2.0 * M_PI * x[1]) : -std::sin(2.0 * M_PI * x[0]);
#else
                    const double val = (axis == 0) ?
                                           std::sin(2.0 * M_PI * x[1]) :
                                           ((axis == 1) ? std::sin(2.0 * M_PI * x[2]) : std::sin(2.0 * M_PI * x[0]));
#endif
                    (*u_data)(i_s) = val;
                }
            }
        }
    }
}

} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Suppress warnings so the test output does not depend on build-path
    // details embedded in warning diagnostics.
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);

    int test_failures = 0;

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        const double current_time = 0.0;
        const double dt = input_db->getDoubleWithDefault("DT", 0.005);
        const double rho = input_db->getDoubleWithDefault("RHO", 1.0);
        const double mu = input_db->getDoubleWithDefault("MU", 1.0);
        const double new_time = current_time + dt;
        const bool use_fixed_le_operators = input_db->getBoolWithDefault("USE_FIXED_LE_OPERATORS", true);
        const bool verify_fgmres_physical_residuals = input_db->keyExists("VERIFY_FGMRES_PHYSICAL_RESIDUALS") ?
                                                          input_db->getBool("VERIFY_FGMRES_PHYSICAL_RESIDUALS") :
                                                          false;
        StructureSpec structure_spec;
        structure_spec.ds = input_db->getDoubleWithDefault("DS", 1.0 / 64.0);
        structure_spec.x_center = input_db->getDoubleWithDefault("X_CENTER", 0.5);
        structure_spec.y_center = input_db->getDoubleWithDefault("Y_CENTER", 0.5);
        structure_spec.x_radius = input_db->getDoubleWithDefault("X_RADIUS", 0.2);
        structure_spec.y_radius = input_db->getDoubleWithDefault("Y_RADIUS", 0.2);
        structure_spec.spring_stiffness = input_db->getDoubleWithDefault("SPRING_STIFFNESS", 2.0e2);
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
        structure_spec.num_curve_points = std::max(3, static_cast<int>(circumference / structure_spec.ds));
        if (structure_spec.num_curve_points < 3)
        {
            TBOX_ERROR("computed num_curve_points must be >= 3\n");
        }

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

        Pointer<IBRedundantInitializer> ib_initializer = new IBRedundantInitializer(
            "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
        structure_spec.finest_ln = input_db->getIntegerWithDefault("MAX_LEVELS", 1) - 1;
        ib_initializer->setStructureNamesOnLevel(structure_spec.finest_ln, { "curve2d" });
        ib_initializer->registerInitStructureFunction(generate_structure, &structure_spec);
        ib_initializer->registerInitSpringDataFunction(generate_springs, &structure_spec);
        ib_method_ops->registerLInitStrategy(ib_initializer);
        Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
        ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);

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

        ib_method_ops->preprocessIntegrateData(current_time, new_time, /*num_cycles*/ 1);
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

        const double lambda = 0.0;
        PoissonSpecifications U_problem_coefs("stokes_ib_solver_components::U_problem_coefs");
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
            new StaggeredStokesOperator("stokes_ib_solver_components::stokes_op", false);
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
        ctx.use_fixed_le_operators = use_fixed_le_operators;
        ctx.time_stepping_type = IBAMR::string_to_enum<TimeSteppingType>(
            input_db->getStringWithDefault("IB_TIME_STEPPING", "MIDPOINT_RULE"));

        StaggeredStokesIBOperator nonlinear_op("stokes_ib_solver_components::nonlinear_op", false);
        nonlinear_op.setOperatorContext(ctx);
        nonlinear_op.setTimeInterval(current_time, new_time);
        nonlinear_op.setSolutionTime(new_time);
        nonlinear_op.initializeOperatorState(*eul_sol_vec, *eul_rhs_vec);

        Pointer<SAMRAIVectorReal<NDIM, double>> nonlinear_probe = eul_sol_vec->cloneVector("nonlinear_probe");
        nonlinear_probe->allocateVectorData();
        nonlinear_probe->setToScalar(0.0);
        set_divergence_free_probe_velocity(nonlinear_probe->getComponentDescriptorIndex(0), patch_hierarchy);

        Pointer<SAMRAIVectorReal<NDIM, double>> f_probe = eul_rhs_vec->cloneVector("f_probe");
        f_probe->allocateVectorData();
        f_probe->setToScalar(0.0);
        nonlinear_op.apply(*nonlinear_probe, *f_probe);

        double nonlinear_side_norm = std::numeric_limits<double>::quiet_NaN();
        double nonlinear_cell_norm = std::numeric_limits<double>::quiet_NaN();
        const bool expect_trivial_nonlinear = (std::abs(structure_spec.spring_stiffness) <= 1.0e-14) &&
                                              (std::abs(rho) <= 1.0e-14) && (std::abs(mu) <= 1.0e-14);
        if (!side_l2_norm_is_finite(
                hier_velocity_data_ops, f_probe->getComponentDescriptorIndex(0), wgt_sc_idx, nonlinear_side_norm) ||
            !cell_l2_norm_is_finite(
                hier_pressure_data_ops, f_probe->getComponentDescriptorIndex(1), wgt_cc_idx, nonlinear_cell_norm))
        {
            ++test_failures;
            pout << "nonlinear operator produced non-finite norm" << std::endl;
        }
        else if (nonlinear_side_norm <= 1.0e-14 && nonlinear_cell_norm <= 1.0e-14)
        {
            if (!expect_trivial_nonlinear)
            {
                ++test_failures;
                pout << "nonlinear operator action is trivial" << std::endl;
            }
        }
        else if (expect_trivial_nonlinear)
        {
            ++test_failures;
            pout << "nonlinear operator action is nontrivial when SPRING_STIFFNESS, RHO, and MU are zero" << std::endl;
        }

        Pointer<StaggeredStokesIBJacobianOperator> jac_op =
            new StaggeredStokesIBJacobianOperator("stokes_ib_solver_components::jacobian_op");
        jac_op->setOperatorContext(ctx);
        jac_op->setTimeInterval(current_time, new_time);
        jac_op->setSolutionTime(new_time);
        jac_op->initializeOperatorState(*eul_sol_vec, *eul_rhs_vec);
        jac_op->formJacobian(*eul_sol_vec);

        Pointer<SAMRAIVectorReal<NDIM, double>> v = eul_sol_vec->cloneVector("v");
        v->allocateVectorData();
        v->setToScalar(0.0);
        if (verify_fgmres_physical_residuals)
        {
            // A rigid translation lies in the elastic coupling nullspace and
            // cannot provide a well-scaled end-to-end IB residual check.
            set_divergence_free_probe_velocity(v->getComponentDescriptorIndex(0), patch_hierarchy);
        }
        else
        {
            hier_velocity_data_ops->setToScalar(v->getComponentDescriptorIndex(0), 1.0, false);
        }
        hier_pressure_data_ops->setToScalar(v->getComponentDescriptorIndex(1), -0.25, false);

        Pointer<SAMRAIVectorReal<NDIM, double>> jv = eul_rhs_vec->cloneVector("jv");
        jv->allocateVectorData();
        jv->setToScalar(0.0);
        jac_op->apply(*v, *jv);

        const double fd_rel_tol = input_db->getDoubleWithDefault("FD_REL_TOL", 5.0e-2);
        Pointer<PETScMFFDJacobianOperator> mffd_jac_op =
            new PETScMFFDJacobianOperator("stokes_ib_solver_components::mffd_jacobian_op", "ib_jac_mffd_");
        mffd_jac_op->setOperator(Pointer<GeneralOperator>(&nonlinear_op, false));
        mffd_jac_op->setTimeInterval(current_time, new_time);
        mffd_jac_op->setSolutionTime(new_time);
        mffd_jac_op->initializeOperatorState(*eul_sol_vec, *eul_rhs_vec);
        mffd_jac_op->formJacobian(*eul_sol_vec);

        Pointer<SAMRAIVectorReal<NDIM, double>> fd_jv = eul_rhs_vec->cloneVector("fd_jv");
        fd_jv->allocateVectorData();
        fd_jv->setToScalar(0.0);
        mffd_jac_op->apply(*v, *fd_jv);

        Pointer<SAMRAIVectorReal<NDIM, double>> diff = eul_rhs_vec->cloneVector("diff");
        diff->allocateVectorData();
        diff->subtract(fd_jv, jv);

        double jv_side_norm = std::numeric_limits<double>::quiet_NaN();
        double jv_cell_norm = std::numeric_limits<double>::quiet_NaN();
        double diff_side_norm = std::numeric_limits<double>::quiet_NaN();
        double diff_cell_norm = std::numeric_limits<double>::quiet_NaN();
        const bool jv_finite =
            side_l2_norm_is_finite(
                hier_velocity_data_ops, jv->getComponentDescriptorIndex(0), wgt_sc_idx, jv_side_norm) &&
            cell_l2_norm_is_finite(
                hier_pressure_data_ops, jv->getComponentDescriptorIndex(1), wgt_cc_idx, jv_cell_norm);
        const bool diff_finite =
            side_l2_norm_is_finite(
                hier_velocity_data_ops, diff->getComponentDescriptorIndex(0), wgt_sc_idx, diff_side_norm) &&
            cell_l2_norm_is_finite(
                hier_pressure_data_ops, diff->getComponentDescriptorIndex(1), wgt_cc_idx, diff_cell_norm);
        if (!jv_finite || !diff_finite)
        {
            ++test_failures;
            pout << "jacobian norms are non-finite" << std::endl;
        }
        else if (jv_side_norm <= 1.0e-14 && jv_cell_norm <= 1.0e-14)
        {
            ++test_failures;
            pout << "jacobian action is trivial" << std::endl;
        }
        else
        {
            const double rel_error =
                std::sqrt(diff_side_norm * diff_side_norm + diff_cell_norm * diff_cell_norm) /
                std::max(std::sqrt(jv_side_norm * jv_side_norm + jv_cell_norm * jv_cell_norm), 1.0e-14);
            pout << "fd_relative_error = " << rel_error << std::endl;
            if (!(rel_error <= fd_rel_tol))
            {
                ++test_failures;
                pout << "fd_relative_error exceeds tolerance: " << fd_rel_tol << std::endl;
            }
        }
        mffd_jac_op->deallocateOperatorState();

        Pointer<Database> stokes_ib_precond_db =
            input_db->isDatabase("stokes_ib_precond_db") ? input_db->getDatabase("stokes_ib_precond_db") : nullptr;

        Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_op = new StaggeredStokesIBLevelRelaxationFACOperator(
            "stokes_ib_solver_components::fac_op", stokes_ib_precond_db, "stokes_ib_pc_");
        Pointer<StaggeredStokesIBJacobianFACPreconditioner> fac_pc = new StaggeredStokesIBJacobianFACPreconditioner(
            "stokes_ib_solver_components::fac_pc", fac_op, stokes_ib_precond_db, "stokes_ib_pc_");
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
        const bool verify_pressure_cav =
            input_db->keyExists("VERIFY_PRESSURE_CAV") ? input_db->getBool("VERIFY_PRESSURE_CAV") : false;
        const bool verify_fac_cycle_observer =
            input_db->keyExists("VERIFY_FAC_CYCLE_OBSERVER") ? input_db->getBool("VERIFY_FAC_CYCLE_OBSERVER") : false;
        Pointer<StaggeredStokesPETScLevelSolver> pressure_cav_level_solver;
        fac_pc->initializeSolverState(*eul_sol_vec, *eul_rhs_vec);
        if (verify_pressure_cav)
        {
            pressure_cav_level_solver = fac_op->getStaggeredStokesPETScLevelSolver(finest_ln);
            const auto& overlap_dofs = pressure_cav_level_solver->getASMSubdomains();
            const auto& nonoverlap_dofs = pressure_cav_level_solver->getASMNonoverlapSubdomains();
            const auto& pressure_seeds = pressure_cav_level_solver->getCouplingAwareASMPressureSeedDOFs();
            const int expected_patch_count = input_db->getInteger("EXPECTED_PRESSURE_CAV_PATCH_COUNT");

            bool seed_order_valid = static_cast<int>(pressure_seeds.size()) == expected_patch_count &&
                                    pressure_seeds.size() == overlap_dofs.size();
            int enlarged_patch_count = 0;
            for (std::size_t k = 0; k < std::min(pressure_seeds.size(), overlap_dofs.size()); ++k)
            {
                seed_order_valid =
                    seed_order_valid &&
                    std::binary_search(overlap_dofs[k].begin(), overlap_dofs[k].end(), pressure_seeds[k]);
                if (overlap_dofs[k].size() > 2 * NDIM + 1) ++enlarged_patch_count;
            }
            const bool partition_absent = nonoverlap_dofs.size() == overlap_dofs.size() &&
                                          std::all_of(nonoverlap_dofs.begin(),
                                                      nonoverlap_dofs.end(),
                                                      [](const auto& dofs) { return dofs.empty(); });
            if (!seed_order_valid || enlarged_patch_count == 0 || !partition_absent) ++test_failures;

            pout << "pressure_cav_patch_count = " << pressure_seeds.size() << std::endl;
            pout << "pressure_cav_enlarged_patch_count = " << enlarged_patch_count << std::endl;
            pout << "pressure_cav_seed_order_valid = " << (seed_order_valid ? "true" : "false") << std::endl;
            pout << "pressure_cav_partition_absent = " << (partition_absent ? "true" : "false") << std::endl;
        }
        Mat SAJ = fac_op->getEulerianElasticityLevelOp(finest_ln);
        jac_op->setIBCouplingJacobian(SAJ);

        using CycleStage = IBTK::FACPreconditioner::CycleStage;
        std::vector<CycleStage> observed_fac_stages;
        std::vector<int> observed_fac_levels;
        bool fac_observer_views_valid = true;
        IBTK::FACPreconditioner::CycleObserver fac_cycle_observer = [&](const CycleStage stage,
                                                                        const int level_num,
                                                                        const SAMRAIVectorReal<NDIM, double>& solution,
                                                                        const SAMRAIVectorReal<NDIM, double>& rhs)
        {
            observed_fac_stages.push_back(stage);
            observed_fac_levels.push_back(level_num);
            fac_observer_views_valid = fac_observer_views_valid && solution.getPatchHierarchy() == patch_hierarchy &&
                                       rhs.getPatchHierarchy() == patch_hierarchy &&
                                       solution.getCoarsestLevelNumber() == 0 && rhs.getCoarsestLevelNumber() == 0 &&
                                       solution.getFinestLevelNumber() == finest_ln &&
                                       rhs.getFinestLevelNumber() == finest_ln && std::isfinite(solution.L2Norm()) &&
                                       std::isfinite(rhs.L2Norm());
        };

        const std::vector<CycleStage> expected_fac_stages = {
            CycleStage::PRE_SMOOTH_INPUT,  CycleStage::PRE_SMOOTH_OUTPUT, CycleStage::COARSE_RHS,
            CycleStage::COARSE_CORRECTION, CycleStage::POST_SMOOTH_INPUT, CycleStage::POST_SMOOTH_OUTPUT
        };
        const std::vector<int> expected_fac_levels = { finest_ln, finest_ln, 0, 0, finest_ln, finest_ln };
        const std::vector<CycleStage> expected_fac_no_pre_stages = { CycleStage::COARSE_RHS,
                                                                     CycleStage::COARSE_CORRECTION,
                                                                     CycleStage::POST_SMOOTH_INPUT,
                                                                     CycleStage::POST_SMOOTH_OUTPUT };
        const std::vector<int> expected_fac_no_pre_levels = { 0, 0, finest_ln, finest_ln };
        Pointer<SAMRAIVectorReal<NDIM, double>> fac_observer_sol;
        bool fac_observer_first_cycle_valid = true;
        bool fac_observer_disabled_silent = true;
        bool fac_observer_no_pre_cycle_valid = true;
        if (verify_fac_cycle_observer)
        {
            const int configured_num_pre_sweeps = fac_pc->getNumPreSmoothingSweeps();
            fac_observer_sol = eul_sol_vec->cloneVector("fac_observer_sol");
            fac_observer_sol->allocateVectorData();
            fac_observer_sol->setToScalar(0.0);
            fac_pc->setCycleObserver(fac_cycle_observer);
            const bool observer_solve_success = fac_pc->solveSystem(*fac_observer_sol, *jv);
            fac_observer_first_cycle_valid = observer_solve_success && observed_fac_stages == expected_fac_stages &&
                                             observed_fac_levels == expected_fac_levels && fac_observer_views_valid &&
                                             std::isfinite(fac_observer_sol->L2Norm()) &&
                                             fac_observer_sol->L2Norm() > 1.0e-14;
            if (!fac_observer_first_cycle_valid) ++test_failures;

            const std::size_t observed_stage_count = observed_fac_stages.size();
            fac_pc->setCycleObserver({});
            fac_observer_sol->setToScalar(0.0);
            const bool disabled_solve_success = fac_pc->solveSystem(*fac_observer_sol, *jv);
            fac_observer_disabled_silent = disabled_solve_success && observed_fac_stages.size() == observed_stage_count;
            if (!fac_observer_disabled_silent) ++test_failures;

            observed_fac_stages.clear();
            observed_fac_levels.clear();
            fac_observer_views_valid = true;
            fac_pc->setNumPreSmoothingSweeps(0);
            fac_pc->setCycleObserver(fac_cycle_observer);
            fac_observer_sol->setToScalar(0.0);
            const bool no_pre_solve_success = fac_pc->solveSystem(*fac_observer_sol, *jv);
            fac_observer_no_pre_cycle_valid =
                no_pre_solve_success && observed_fac_stages == expected_fac_no_pre_stages &&
                observed_fac_levels == expected_fac_no_pre_levels && fac_observer_views_valid &&
                std::isfinite(fac_observer_sol->L2Norm()) && fac_observer_sol->L2Norm() > 1.0e-14;
            if (!fac_observer_no_pre_cycle_valid) ++test_failures;
            fac_pc->setCycleObserver({});
            fac_pc->setNumPreSmoothingSweeps(configured_num_pre_sweeps);
        }

        {
            PetscErrorCode ierr = 0;
            Mat saj_unscaled = nullptr;
            ierr = MatPtAP(A, J, MAT_INITIAL_MATRIX, 1.0, &saj_unscaled);
            IBTK_CHKERRQ(ierr);
            Pointer<PatchLevel<NDIM>> finest_level = patch_hierarchy->getPatchLevel(finest_ln);
            Pointer<CartesianGridGeometry<NDIM>> finest_grid_geom = patch_hierarchy->getGridGeometry();
            const double* const dx0 = finest_grid_geom->getDx();
            const IntVector<NDIM>& ratio = finest_level->getRatio();
            double cell_volume = 1.0;
            for (unsigned d = 0; d < NDIM; ++d)
            {
                cell_volume *= dx0[d] / static_cast<double>(ratio(d));
            }
            const double theta_ds = 2.0 * M_PI / static_cast<double>(structure_spec.num_curve_points);

            Mat saj_cell_scaled = nullptr;
            ierr = MatDuplicate(saj_unscaled, MAT_COPY_VALUES, &saj_cell_scaled);
            IBTK_CHKERRQ(ierr);
            ierr = MatScale(saj_cell_scaled, -dt / cell_volume);
            IBTK_CHKERRQ(ierr);

            Mat saj_theta_scaled = nullptr;
            ierr = MatDuplicate(saj_unscaled, MAT_COPY_VALUES, &saj_theta_scaled);
            IBTK_CHKERRQ(ierr);
            ierr = MatScale(saj_theta_scaled, -dt * theta_ds / cell_volume);
            IBTK_CHKERRQ(ierr);

            auto compute_rel_mat_error = [](Mat lhs, Mat rhs) -> double
            {
                Mat diff = nullptr;
                PetscErrorCode ierr_local = MatDuplicate(lhs, MAT_COPY_VALUES, &diff);
                IBTK_CHKERRQ(ierr_local);
                ierr_local = MatAXPY(diff, -1.0, rhs, DIFFERENT_NONZERO_PATTERN);
                IBTK_CHKERRQ(ierr_local);
                PetscReal diff_norm = 0.0;
                PetscReal rhs_norm = 0.0;
                ierr_local = MatNorm(diff, NORM_FROBENIUS, &diff_norm);
                IBTK_CHKERRQ(ierr_local);
                ierr_local = MatNorm(lhs, NORM_FROBENIUS, &rhs_norm);
                IBTK_CHKERRQ(ierr_local);
                ierr_local = MatDestroy(&diff);
                IBTK_CHKERRQ(ierr_local);
                return static_cast<double>(diff_norm) / std::max(static_cast<double>(rhs_norm), 1.0e-14);
            };

            const double saj_cell_scaled_rel_error = compute_rel_mat_error(SAJ, saj_cell_scaled);
            const double saj_theta_scaled_rel_error = compute_rel_mat_error(SAJ, saj_theta_scaled);
            pout << "saj_cell_scaled_relative_error = " << saj_cell_scaled_rel_error << std::endl;
            pout << "saj_theta_scaled_relative_error = " << saj_theta_scaled_rel_error << std::endl;

            ierr = MatDestroy(&saj_unscaled);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&saj_cell_scaled);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&saj_theta_scaled);
            IBTK_CHKERRQ(ierr);
        }

        const bool run_saj_vector_compare = input_db->getBoolWithDefault("RUN_SAJ_VECTOR_COMPARE", false);
        if (run_saj_vector_compare)
        {
            Pointer<SAMRAIVectorReal<NDIM, double>> saj_jv = eul_rhs_vec->cloneVector("saj_jv");
            saj_jv->allocateVectorData();
            saj_jv->setToScalar(0.0);
            jac_op->apply(*v, *saj_jv);

            Pointer<SAMRAIVectorReal<NDIM, double>> saj_diff = eul_rhs_vec->cloneVector("saj_diff");
            saj_diff->allocateVectorData();
            saj_diff->subtract(saj_jv, jv);

            double saj_jv_side_norm = std::numeric_limits<double>::quiet_NaN();
            double saj_jv_cell_norm = std::numeric_limits<double>::quiet_NaN();
            double saj_diff_side_norm = std::numeric_limits<double>::quiet_NaN();
            double saj_diff_cell_norm = std::numeric_limits<double>::quiet_NaN();
            const bool saj_jv_finite =
                side_l2_norm_is_finite(
                    hier_velocity_data_ops, saj_jv->getComponentDescriptorIndex(0), wgt_sc_idx, saj_jv_side_norm) &&
                cell_l2_norm_is_finite(
                    hier_pressure_data_ops, saj_jv->getComponentDescriptorIndex(1), wgt_cc_idx, saj_jv_cell_norm);
            const bool saj_diff_finite =
                side_l2_norm_is_finite(
                    hier_velocity_data_ops, saj_diff->getComponentDescriptorIndex(0), wgt_sc_idx, saj_diff_side_norm) &&
                cell_l2_norm_is_finite(
                    hier_pressure_data_ops, saj_diff->getComponentDescriptorIndex(1), wgt_cc_idx, saj_diff_cell_norm);
            if (!saj_jv_finite || !saj_diff_finite)
            {
                ++test_failures;
                pout << "saj jacobian comparison norms are non-finite" << std::endl;
            }
            else
            {
                const double saj_rel_error =
                    std::sqrt(saj_diff_side_norm * saj_diff_side_norm + saj_diff_cell_norm * saj_diff_cell_norm) /
                    std::max(std::sqrt(saj_jv_side_norm * saj_jv_side_norm + saj_jv_cell_norm * saj_jv_cell_norm),
                             1.0e-14);
                pout << "saj_relative_error = " << saj_rel_error << std::endl;
                const double saj_rel_tol = input_db->getDoubleWithDefault("SAJ_REL_TOL", 5.0e-12);
                if (!(saj_rel_error <= saj_rel_tol))
                {
                    ++test_failures;
                    pout << "saj_relative_error exceeds tolerance: " << saj_rel_tol << std::endl;
                }
            }
        }

        Pointer<PETScKrylovLinearSolver> linear_solver =
            new PETScKrylovLinearSolver("stokes_ib_solver_components::linear_solver", nullptr, "ib_");
        linear_solver->setOperator(jac_op);
        linear_solver->setPreconditioner(fac_pc);
        linear_solver->setTimeInterval(current_time, new_time);
        linear_solver->setSolutionTime(new_time);
        linear_solver->setInitialGuessNonzero(false);

        Pointer<SAMRAIVectorReal<NDIM, double>> linear_sol = eul_sol_vec->cloneVector("linear_sol");
        linear_sol->allocateVectorData();
        linear_sol->setToScalar(0.0);
        if (verify_fgmres_physical_residuals)
        {
            linear_solver->initializeSolverState(*linear_sol, *jv);
            SAJ = fac_op->getEulerianElasticityLevelOp(finest_ln);
            PetscErrorCode ierr = KSPSetResidualHistory(linear_solver->getPETScKSP(), nullptr, 201, PETSC_TRUE);
            IBTK_CHKERRQ(ierr);
        }
        const bool linear_success = linear_solver->solveSystem(*linear_sol, *jv);
        if (!linear_success)
        {
            ++test_failures;
            pout << "krylov linear solve failed" << std::endl;
        }
        pout << "krylov_linear_iterations = " << linear_solver->getNumIterations() << std::endl;
        pout << "krylov_linear_residual_norm = " << linear_solver->getResidualNorm() << std::endl;

        double linear_side_norm = std::numeric_limits<double>::quiet_NaN();
        double linear_cell_norm = std::numeric_limits<double>::quiet_NaN();
        if (!side_l2_norm_is_finite(
                hier_velocity_data_ops, linear_sol->getComponentDescriptorIndex(0), wgt_sc_idx, linear_side_norm) ||
            !cell_l2_norm_is_finite(
                hier_pressure_data_ops, linear_sol->getComponentDescriptorIndex(1), wgt_cc_idx, linear_cell_norm))
        {
            ++test_failures;
            pout << "krylov linear solve produced non-finite norm" << std::endl;
        }
        else if (linear_side_norm <= 1.0e-14 && linear_cell_norm <= 1.0e-14)
        {
            ++test_failures;
            pout << "krylov linear solve action is trivial" << std::endl;
        }

        if (verify_fgmres_physical_residuals)
        {
            const double original_relative_residual_tol = input_db->getDouble("FGMRES_ORIGINAL_RELATIVE_RESIDUAL_TOL");
            const double ib_coupling_relative_residual_tol = input_db->getDouble("IB_COUPLING_RELATIVE_RESIDUAL_TOL");
            const char* ksp_type = nullptr;
            PetscErrorCode ierr = KSPGetType(linear_solver->getPETScKSP(), &ksp_type);
            IBTK_CHKERRQ(ierr);
            const bool fgmres_type_valid = ksp_type && std::string(ksp_type) == KSPFGMRES;

            const PetscReal* residual_history = nullptr;
            PetscInt residual_history_size = 0;
            ierr = KSPGetResidualHistory(linear_solver->getPETScKSP(), &residual_history, &residual_history_size);
            IBTK_CHKERRQ(ierr);
            bool fgmres_history_valid = residual_history_size >= 2;
            for (PetscInt k = 0; k < residual_history_size; ++k)
            {
                fgmres_history_valid = fgmres_history_valid && std::isfinite(residual_history[k]);
            }
            const double fgmres_history_final_relative =
                residual_history_size >= 2 ?
                    residual_history[residual_history_size - 1] / std::max(residual_history[0], 1.0e-30) :
                    std::numeric_limits<double>::quiet_NaN();
            fgmres_history_valid = fgmres_history_valid && std::isfinite(fgmres_history_final_relative) &&
                                   fgmres_history_final_relative <= original_relative_residual_tol;

            Pointer<SAMRAIVectorReal<NDIM, double>> linear_action = eul_rhs_vec->cloneVector("linear_action");
            Pointer<SAMRAIVectorReal<NDIM, double>> original_residual = eul_rhs_vec->cloneVector("original_residual");
            Pointer<SAMRAIVectorReal<NDIM, double>> ib_action_matrix_free =
                eul_rhs_vec->cloneVector("ib_action_matrix_free");
            for (const Pointer<SAMRAIVectorReal<NDIM, double>>& diagnostic_vec :
                 { linear_action, original_residual, ib_action_matrix_free })
            {
                diagnostic_vec->allocateVectorData();
                diagnostic_vec->setToScalar(0.0);
            }

            // Re-evaluate the original matrix-free Jacobian after the solve;
            // the KSP-reported residual is not used as a physical substitute.
            jac_op->apply(*linear_sol, *linear_action);
            original_residual->subtract(jv, linear_action);

            // Evaluate the matrix-free IB action directly. With zero velocity
            // coefficients and zero pressure, the auxiliary Jacobian's
            // momentum component contains only the live IB coupling action;
            // its divergence component is intentionally discarded.
            PoissonSpecifications ib_only_coefs("stokes_ib_solver_components::ib_only_coefs");
            ib_only_coefs.setCConstant(0.0);
            ib_only_coefs.setDConstant(0.0);
            Pointer<StaggeredStokesOperator> ib_only_stokes_op =
                new StaggeredStokesOperator("stokes_ib_solver_components::ib_only_stokes_op", false);
            ib_only_stokes_op->setVelocityPoissonSpecifications(ib_only_coefs);
            ib_only_stokes_op->setPhysicalBcCoefs(u_bc_coefs, nullptr);
            ib_only_stokes_op->setTimeInterval(current_time, new_time);
            ib_only_stokes_op->setSolutionTime(new_time);
            StaggeredStokesIBOperator::Context ib_only_ctx = ctx;
            ib_only_ctx.stokes_op = ib_only_stokes_op;
            Pointer<StaggeredStokesIBJacobianOperator> ib_only_jac_op =
                new StaggeredStokesIBJacobianOperator("stokes_ib_solver_components::ib_only_jacobian_op");
            ib_only_jac_op->setOperatorContext(ib_only_ctx);
            ib_only_jac_op->setTimeInterval(current_time, new_time);
            ib_only_jac_op->setSolutionTime(new_time);
            ib_only_jac_op->initializeOperatorState(*eul_sol_vec, *eul_rhs_vec);
            ib_only_jac_op->formJacobian(*eul_sol_vec);
            Pointer<SAMRAIVectorReal<NDIM, double>> ib_input = linear_sol->cloneVector("ib_input");
            ib_input->allocateVectorData();
            ib_input->copyVector(linear_sol);
            hier_pressure_data_ops->setToScalar(ib_input->getComponentDescriptorIndex(1), 0.0, false);
            ib_only_jac_op->apply(*ib_input, *ib_action_matrix_free);
            hier_pressure_data_ops->setToScalar(ib_action_matrix_free->getComponentDescriptorIndex(1), 0.0, false);

            const double rhs_momentum_norm =
                hier_velocity_data_ops->L2Norm(jv->getComponentDescriptorIndex(0), wgt_sc_idx);
            const double rhs_divergence_norm =
                hier_pressure_data_ops->L2Norm(jv->getComponentDescriptorIndex(1), wgt_cc_idx);
            const double momentum_residual_norm =
                hier_velocity_data_ops->L2Norm(original_residual->getComponentDescriptorIndex(0), wgt_sc_idx);
            const double divergence_residual_norm =
                hier_pressure_data_ops->L2Norm(original_residual->getComponentDescriptorIndex(1), wgt_cc_idx);
            const double rhs_norm = std::hypot(rhs_momentum_norm, rhs_divergence_norm);
            const double original_residual_norm = std::hypot(momentum_residual_norm, divergence_residual_norm);
            const double denominator_floor = std::sqrt(std::numeric_limits<double>::epsilon()) * rhs_norm;
            const double original_relative_residual = original_residual_norm / std::max(rhs_norm, 1.0e-30);
            const double momentum_relative_residual =
                momentum_residual_norm / std::max({ rhs_momentum_norm, denominator_floor, 1.0e-30 });
            const double divergence_relative_residual =
                divergence_residual_norm / std::max({ rhs_divergence_norm, denominator_floor, 1.0e-30 });

            Vec saj_input = nullptr;
            Vec saj_output = nullptr;
            Vec rhs_output = nullptr;
            Vec ib_matrix_free_output = nullptr;
            ierr = MatCreateVecs(SAJ, &saj_input, &saj_output);
            IBTK_CHKERRQ(ierr);
            ierr = VecDuplicate(saj_output, &rhs_output);
            IBTK_CHKERRQ(ierr);
            ierr = VecDuplicate(saj_output, &ib_matrix_free_output);
            IBTK_CHKERRQ(ierr);
            StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(saj_input,
                                                                  linear_sol->getComponentDescriptorIndex(0),
                                                                  u_dof_index_idx,
                                                                  linear_sol->getComponentDescriptorIndex(1),
                                                                  p_dof_index_idx,
                                                                  patch_hierarchy->getPatchLevel(finest_ln));
            ierr = MatMult(SAJ, saj_input, saj_output);
            IBTK_CHKERRQ(ierr);
            StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(rhs_output,
                                                                  jv->getComponentDescriptorIndex(0),
                                                                  u_dof_index_idx,
                                                                  jv->getComponentDescriptorIndex(1),
                                                                  p_dof_index_idx,
                                                                  patch_hierarchy->getPatchLevel(finest_ln));
            StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(ib_matrix_free_output,
                                                                  ib_action_matrix_free->getComponentDescriptorIndex(0),
                                                                  u_dof_index_idx,
                                                                  ib_action_matrix_free->getComponentDescriptorIndex(1),
                                                                  p_dof_index_idx,
                                                                  patch_hierarchy->getPatchLevel(finest_ln));
            double ib_matrix_free_action_norm = std::numeric_limits<double>::quiet_NaN();
            double ib_saj_action_norm = std::numeric_limits<double>::quiet_NaN();
            double rhs_algebraic_norm = std::numeric_limits<double>::quiet_NaN();
            ierr = VecNorm(rhs_output, NORM_2, &rhs_algebraic_norm);
            IBTK_CHKERRQ(ierr);
            ierr = VecNorm(ib_matrix_free_output, NORM_2, &ib_matrix_free_action_norm);
            IBTK_CHKERRQ(ierr);
            ierr = VecNorm(saj_output, NORM_2, &ib_saj_action_norm);
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPY(ib_matrix_free_output, -1.0, saj_output);
            IBTK_CHKERRQ(ierr);
            double ib_coupling_residual_norm = std::numeric_limits<double>::quiet_NaN();
            ierr = VecNorm(ib_matrix_free_output, NORM_2, &ib_coupling_residual_norm);
            IBTK_CHKERRQ(ierr);
            const double ib_coupling_relative_residual =
                ib_coupling_residual_norm /
                std::max({ ib_matrix_free_action_norm,
                           ib_saj_action_norm,
                           std::sqrt(std::numeric_limits<double>::epsilon()) * rhs_algebraic_norm,
                           1.0e-30 });

            const bool original_residual_valid = std::isfinite(original_relative_residual) &&
                                                 original_relative_residual <= original_relative_residual_tol;
            const bool ib_coupling_residual_valid = std::isfinite(ib_coupling_relative_residual) &&
                                                    ib_coupling_relative_residual <= ib_coupling_relative_residual_tol;
            if (!fgmres_type_valid || !fgmres_history_valid || !original_residual_valid || !ib_coupling_residual_valid)
            {
                ++test_failures;
            }

            pout << "fgmres_type_valid = " << (fgmres_type_valid ? "true" : "false") << std::endl;
            pout << "fgmres_residual_history_size = " << residual_history_size << std::endl;
            pout << "fgmres_history_final_relative = " << fgmres_history_final_relative << std::endl;
            pout << "fgmres_history_valid = " << (fgmres_history_valid ? "true" : "false") << std::endl;
            pout << "original_momentum_residual_l2 = " << momentum_residual_norm << std::endl;
            pout << "original_momentum_relative_residual = " << momentum_relative_residual << std::endl;
            pout << "original_divergence_residual_l2 = " << divergence_residual_norm << std::endl;
            pout << "original_divergence_relative_residual = " << divergence_relative_residual << std::endl;
            pout << "original_relative_residual = " << original_relative_residual << std::endl;
            pout << "ib_coupling_residual_l2 = " << ib_coupling_residual_norm << std::endl;
            pout << "ib_coupling_relative_residual = " << ib_coupling_relative_residual << std::endl;
            pout << "fgmres_original_residual_valid = " << (original_residual_valid ? "true" : "false") << std::endl;
            pout << "ib_coupling_residual_valid = " << (ib_coupling_residual_valid ? "true" : "false") << std::endl;

            ierr = VecDestroy(&ib_matrix_free_output);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&rhs_output);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&saj_output);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&saj_input);
            IBTK_CHKERRQ(ierr);
            ib_input->deallocateVectorData();
            ib_only_jac_op->deallocateOperatorState();
            for (const Pointer<SAMRAIVectorReal<NDIM, double>>& diagnostic_vec :
                 { linear_action, original_residual, ib_action_matrix_free })
            {
                diagnostic_vec->deallocateVectorData();
            }
            linear_solver->deallocateSolverState();
        }

        bool fac_observer_reinitialize_valid = true;
        if (verify_fac_cycle_observer)
        {
            observed_fac_stages.clear();
            observed_fac_levels.clear();
            fac_observer_views_valid = true;
            fac_pc->setCycleObserver(fac_cycle_observer);
        }
        fac_pc->deallocateSolverState();
        if (verify_fac_cycle_observer)
        {
            fac_pc->initializeSolverState(*eul_sol_vec, *eul_rhs_vec);
            fac_observer_sol->setToScalar(0.0);
            const bool reinitialized_solve_success = fac_pc->solveSystem(*fac_observer_sol, *jv);
            fac_observer_reinitialize_valid =
                reinitialized_solve_success && observed_fac_stages == expected_fac_stages &&
                observed_fac_levels == expected_fac_levels && fac_observer_views_valid &&
                std::isfinite(fac_observer_sol->L2Norm()) && fac_observer_sol->L2Norm() > 1.0e-14;
            if (!fac_observer_reinitialize_valid) ++test_failures;
            fac_pc->setCycleObserver({});
            fac_pc->deallocateSolverState();

            pout << "fac_cycle_observer_first_cycle_valid = " << (fac_observer_first_cycle_valid ? "true" : "false")
                 << std::endl;
            pout << "fac_cycle_observer_disabled_silent = " << (fac_observer_disabled_silent ? "true" : "false")
                 << std::endl;
            pout << "fac_cycle_observer_no_pre_cycle_valid = " << (fac_observer_no_pre_cycle_valid ? "true" : "false")
                 << std::endl;
            pout << "fac_cycle_observer_reinitialize_valid = " << (fac_observer_reinitialize_valid ? "true" : "false")
                 << std::endl;
        }
        if (verify_pressure_cav)
        {
            const bool state_cleared = pressure_cav_level_solver->getCouplingAwareASMPressureSeedDOFs().empty();
            if (!state_cleared) ++test_failures;
            pout << "pressure_cav_state_cleared = " << (state_cleared ? "true" : "false") << std::endl;
        }

        jac_op->deallocateOperatorState();
        nonlinear_op.deallocateOperatorState();

        ib_method_ops->postprocessIntegrateData(current_time, new_time, /*num_cycles*/ 1);

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

        plog << "Input database:\n";
        input_db->printClassData(plog);
        pout << "test_failures = " << test_failures << std::endl;
    }

    return test_failures;
}
