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
#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartSideDoubleRT0Coarsen.h>
#include <ibtk/CartSideDoubleRT0Refine.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LDataManager.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/ibtk_utilities.h>

#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenOperator.h>
#include <CoarsenSchedule.h>
#include <GriddingAlgorithm.h>
#include <HierarchySideDataOpsReal.h>
#include <IntVector.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>
#include <RefineSchedule.h>
#include <SAMRAI_config.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <string>
#include <utility>
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
ib4_interp_fcn(const double r, double* const w)
{
    const double q = std::sqrt(-7.0 + 12.0 * r - 4.0 * r * r);
    w[0] = 0.125 * (5.0 - 2.0 * r - q);
    w[1] = 0.125 * (5.0 - 2.0 * r + q);
    w[2] = 0.125 * (-1.0 + 2.0 * r + q);
    w[3] = 0.125 * (-1.0 + 2.0 * r - q);
}

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
        vertex_posn.clear();
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
build_coarse_op_from_samrai_transfers(Mat coarse_op,
                                      Mat fine_op,
                                      const std::vector<int>& coarse_num_dofs_per_proc,
                                      const std::vector<int>& fine_num_dofs_per_proc,
                                      Pointer<SideVariable<NDIM, double>> u_var,
                                      const int u_idx,
                                      const int u_dof_index_idx,
                                      Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                      IntVector<NDIM> gcw)
{
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    Pointer<PatchLevel<NDIM>> coarse_level = patch_hierarchy->getPatchLevel(coarsest_ln);
    Pointer<PatchLevel<NDIM>> fine_level = patch_hierarchy->getPatchLevel(finest_ln);

    Pointer<CartesianGridGeometry<NDIM>> geometry = patch_hierarchy->getGridGeometry();
    geometry->addSpatialCoarsenOperator(new CartSideDoubleRT0Coarsen(gcw));
    geometry->addSpatialRefineOperator(new CartSideDoubleRT0Refine());

    Pointer<RefineOperator<NDIM>> prolongation_op = geometry->lookupRefineOperator(u_var, "RT0_REFINE");
    Pointer<CoarsenOperator<NDIM>> restriction_op = geometry->lookupCoarsenOperator(u_var, "RT0_COARSEN");

    Pointer<RefineAlgorithm<NDIM>> prolongation_refine_algorithm = new RefineAlgorithm<NDIM>();
    Pointer<CoarsenAlgorithm<NDIM>> restriction_coarsen_algorithm = new CoarsenAlgorithm<NDIM>();
    prolongation_refine_algorithm->registerRefine(u_idx, u_idx, u_idx, prolongation_op, nullptr);
    restriction_coarsen_algorithm->registerCoarsen(u_idx, u_idx, restriction_op, IntVector<NDIM>(0));

    Pointer<RefineSchedule<NDIM>> prolongation_schedule = prolongation_refine_algorithm->createSchedule(
        fine_level, Pointer<PatchLevel<NDIM>>(), coarsest_ln, patch_hierarchy, nullptr);
    Pointer<CoarsenSchedule<NDIM>> restriction_schedule =
        restriction_coarsen_algorithm->createSchedule(coarse_level, fine_level);

    const int mpi_rank = IBTK_MPI::getRank();
    const int n_local_coarse = coarse_num_dofs_per_proc[mpi_rank];
    const int i_lower_coarse =
        std::accumulate(coarse_num_dofs_per_proc.begin(), coarse_num_dofs_per_proc.begin() + mpi_rank, 0);
    const int n_total_coarse = std::accumulate(coarse_num_dofs_per_proc.begin(), coarse_num_dofs_per_proc.end(), 0);

    const int n_local_fine = fine_num_dofs_per_proc[mpi_rank];
    const int n_total_fine = std::accumulate(fine_num_dofs_per_proc.begin(), fine_num_dofs_per_proc.end(), 0);

    Vec x = nullptr, y = nullptr, X = nullptr, Y = nullptr;
    int ierr = VecCreateMPI(PETSC_COMM_WORLD, n_local_coarse, n_total_coarse, &x);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, n_local_coarse, n_total_coarse, &y);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, n_local_fine, n_total_fine, &X);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, n_local_fine, n_total_fine, &Y);
    IBTK_CHKERRQ(ierr);

    Pointer<RefineSchedule<NDIM>> ghost_fill_sched_coarse =
        IBTK::PETScVecUtilities::constructGhostFillSchedule(u_idx, coarse_level);
    Pointer<RefineSchedule<NDIM>> ghost_fill_sched_fine =
        IBTK::PETScVecUtilities::constructGhostFillSchedule(u_idx, fine_level);
    Pointer<RefineSchedule<NDIM>> data_synch_sched_coarse =
        IBTK::PETScVecUtilities::constructDataSynchSchedule(u_idx, coarse_level);
    Pointer<RefineSchedule<NDIM>> data_synch_sched_fine =
        IBTK::PETScVecUtilities::constructDataSynchSchedule(u_idx, fine_level);

    Pointer<HierarchySideDataOpsReal<NDIM, double>> hier_velocity_data_ops =
        new HierarchySideDataOpsReal<NDIM, double>(patch_hierarchy, coarsest_ln, finest_ln);

    for (int col = 0; col < n_total_coarse; ++col)
    {
        ierr = VecZeroEntries(x);
        IBTK_CHKERRQ(ierr);
        ierr = VecSetValue(x, col, 1.0, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyBegin(x);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(x);
        IBTK_CHKERRQ(ierr);

        IBTK::PETScVecUtilities::copyFromPatchLevelVec(
            x, u_idx, u_dof_index_idx, coarse_level, data_synch_sched_coarse, ghost_fill_sched_coarse);

        hier_velocity_data_ops->resetLevels(finest_ln, finest_ln);
        hier_velocity_data_ops->setToScalar(u_idx, 0.0, false);
        hier_velocity_data_ops->resetLevels(coarsest_ln, finest_ln);

        prolongation_schedule->fillData(0.0);

        IBTK::PETScVecUtilities::copyToPatchLevelVec(X, u_idx, u_dof_index_idx, fine_level);
        ierr = MatMult(fine_op, X, Y);
        IBTK_CHKERRQ(ierr);

        hier_velocity_data_ops->resetLevels(finest_ln, finest_ln);
        hier_velocity_data_ops->setToScalar(u_idx, 0.0, false);
        hier_velocity_data_ops->resetLevels(coarsest_ln, finest_ln);

        IBTK::PETScVecUtilities::copyFromPatchLevelVec(
            Y, u_idx, u_dof_index_idx, fine_level, data_synch_sched_fine, ghost_fill_sched_fine);

        hier_velocity_data_ops->resetLevels(coarsest_ln, coarsest_ln);
        hier_velocity_data_ops->setToScalar(u_idx, 0.0, false);
        hier_velocity_data_ops->resetLevels(coarsest_ln, finest_ln);

        restriction_schedule->coarsenData();
        IBTK::PETScVecUtilities::copyToPatchLevelVec(y, u_idx, u_dof_index_idx, coarse_level);

        const PetscScalar* y_array = nullptr;
        ierr = VecGetArrayRead(y, &y_array);
        IBTK_CHKERRQ(ierr);

        std::vector<PetscInt> nnz_indices;
        std::vector<PetscScalar> nnz_values;
        nnz_indices.reserve(static_cast<std::size_t>(n_local_coarse));
        nnz_values.reserve(static_cast<std::size_t>(n_local_coarse));

        for (int j = 0; j < n_local_coarse; ++j)
        {
            const PetscScalar val = y_array[j];
            if (IBTK::abs_equal_eps(val, 0.0)) continue;
            const PetscInt global_idx = i_lower_coarse + j;
            nnz_indices.push_back(global_idx);
            nnz_values.push_back(val);
        }

        if (!nnz_indices.empty())
        {
            const PetscInt nnz_size = static_cast<PetscInt>(nnz_indices.size());
            const PetscInt col_idx = static_cast<PetscInt>(col);
            ierr = MatSetValues(coarse_op, nnz_size, nnz_indices.data(), 1, &col_idx, nnz_values.data(), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }

        ierr = VecRestoreArrayRead(y, &y_array);
        IBTK_CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(coarse_op, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(coarse_op, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    ierr = VecDestroy(&x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&y);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&X);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&Y);
    IBTK_CHKERRQ(ierr);
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR("implicit_stokes_ib_coarse_saj_transfer_parity_01 requires serial execution (np=1).\n");
    }

#ifndef IBTK_HAVE_SILO
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
#endif

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    load_petsc_options_file(input_db, argc, argv);

    const double current_time = input_db->getDoubleWithDefault("START_TIME", 0.0);
    const double dt = input_db->getDoubleWithDefault("DT", 0.03125);
    const double new_time = current_time + dt;
    const int num_cycles = input_db->getIntegerWithDefault("NUM_CYCLES", 1);
    const int tag_buffer = input_db->getIntegerWithDefault("TAG_BUFFER", 1);
    StructureSpec structure_spec;
    const double l_domain = input_db->getDoubleWithDefault("L_DOMAIN", 1.0);
    const int n_grid = input_db->getIntegerWithDefault("N_GRID", 32);
    const double r_cyl = input_db->getDoubleWithDefault("R_CYL", 0.25);
    const double alpha = input_db->getDoubleWithDefault("ALPHA", 0.23);
    const double beta = input_db->getDoubleWithDefault("BETA", (r_cyl * r_cyl) / alpha);
    const double k_spring = input_db->getDoubleWithDefault("K_SPRING", 1.0e2);

    if (!(l_domain > 0.0) || !(n_grid > 0) || !(r_cyl > 0.0) || !(alpha > 0.0) || !(beta > 0.0) || !(k_spring > 0.0))
    {
        TBOX_ERROR("L_DOMAIN, N_GRID, R_CYL, ALPHA, BETA, and K_SPRING must be positive\n");
    }

    const double ds_default = (l_domain / static_cast<double>(n_grid)) / r_cyl;
    structure_spec.ds = input_db->getDoubleWithDefault("DS", ds_default);
    structure_spec.x_center = input_db->getDoubleWithDefault("X_CENTER", 0.5);
    structure_spec.y_center = input_db->getDoubleWithDefault("Y_CENTER", 0.5);
    structure_spec.x_radius = input_db->getDoubleWithDefault("X_RADIUS", alpha);
    structure_spec.y_radius = input_db->getDoubleWithDefault("Y_RADIUS", beta);
    structure_spec.spring_stiffness = input_db->getDoubleWithDefault("SPRING_STIFFNESS", k_spring);
    if (!(structure_spec.ds > 0.0))
    {
        TBOX_ERROR("DS must be positive\n");
    }
    if (!(structure_spec.x_radius > 0.0) || !(structure_spec.y_radius > 0.0))
    {
        TBOX_ERROR("X_RADIUS and Y_RADIUS must be positive\n");
    }
    if (!(structure_spec.spring_stiffness > 0.0))
    {
        TBOX_ERROR("SPRING_STIFFNESS must be positive\n");
    }

    const double a = structure_spec.x_radius;
    const double b = structure_spec.y_radius;
    const double h = std::pow((a - b) / (a + b), 2.0);
    const double circumference = M_PI * (a + b) * (1.0 + 3.0 * h / (10.0 + std::sqrt(std::max(0.0, 4.0 - 3.0 * h))));
    structure_spec.num_curve_points = std::max(3, static_cast<int>(circumference / structure_spec.ds));

    Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
    ib_method_ops->setUseFixedLEOperators(true);

    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", ib_method_ops, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
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
    ib_initializer->setStructureNamesOnLevel(structure_spec.finest_ln, { "parity_curve2d" });
    ib_initializer->registerInitStructureFunction(generate_structure, &structure_spec);
    ib_initializer->registerInitSpringDataFunction(generate_springs, &structure_spec);
    ib_method_ops->registerLInitStrategy(ib_initializer);
    Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
    ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);

    gridding_algorithm->makeCoarsestLevel(patch_hierarchy, current_time);
    int level_number = 0;
    bool done = false;
    while (!done && gridding_algorithm->levelCanBeRefined(level_number))
    {
        gridding_algorithm->makeFinerLevel(patch_hierarchy, current_time, true, tag_buffer);
        done = !patch_hierarchy->finerLevelExists(level_number);
        ++level_number;
    }

    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    if (finest_ln < 1)
    {
        TBOX_ERROR("implicit_stokes_ib_coarse_saj_transfer_parity_01 requires at least two levels.\n");
    }

    const int coarsest_ln = 0;
    const double data_time = new_time;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> current_ctx = var_db->getContext("coarse_saj_transfer_parity_current_ctx");
    Pointer<VariableContext> scratch_ctx = var_db->getContext("coarse_saj_transfer_parity_scratch_ctx");
    IntVector<NDIM> gcw = ib_method_ops->getMinimumGhostCellWidth();

    Pointer<SideVariable<NDIM, int>> u_full_dof_var = new SideVariable<NDIM, int>("coarse_saj_u_full_dof");
    Pointer<CellVariable<NDIM, int>> p_full_dof_var = new CellVariable<NDIM, int>("coarse_saj_p_full_dof");
    Pointer<SideVariable<NDIM, int>> u_vel_dof_var = new SideVariable<NDIM, int>("coarse_saj_u_vel_dof");
    Pointer<SideVariable<NDIM, double>> u_data_var = new SideVariable<NDIM, double>("coarse_saj_u_data");

    const int u_current_idx = var_db->registerVariableAndContext(u_data_var, current_ctx, gcw);
    const int u_data_idx = var_db->registerVariableAndContext(u_data_var, scratch_ctx, gcw);
    const int u_full_dof_index_idx =
        var_db->registerVariableAndContext(u_full_dof_var, scratch_ctx, IntVector<NDIM>(1));
    const int p_full_dof_index_idx =
        var_db->registerVariableAndContext(p_full_dof_var, scratch_ctx, IntVector<NDIM>(1));
    const int u_vel_dof_index_idx = var_db->registerVariableAndContext(u_vel_dof_var, scratch_ctx, IntVector<NDIM>(1));

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(u_current_idx, current_time);
        level->allocatePatchData(u_data_idx, current_time);
        level->allocatePatchData(u_full_dof_index_idx, data_time);
        level->allocatePatchData(p_full_dof_index_idx, data_time);
        level->allocatePatchData(u_vel_dof_index_idx, data_time);
    }

    Pointer<HierarchySideDataOpsReal<NDIM, double>> hier_velocity_data_ops =
        new HierarchySideDataOpsReal<NDIM, double>(patch_hierarchy, 0, finest_ln);
    hier_velocity_data_ops->setToScalar(u_current_idx, 0.0, false);
    hier_velocity_data_ops->setToScalar(u_data_idx, 0.0, false);

    std::vector<Pointer<CoarsenSchedule<NDIM>>> u_synch_scheds(finest_ln + 1);
    std::vector<Pointer<RefineSchedule<NDIM>>> u_ghost_fill_scheds(finest_ln + 1);
    ib_method_ops->initializePatchHierarchy(
        patch_hierarchy, gridding_algorithm, u_current_idx, u_synch_scheds, u_ghost_fill_scheds, 0, current_time, true);
    ib_method_ops->freeLInitStrategy();
    ib_initializer.setNull();
    ib_method_ops->preprocessIntegrateData(current_time, new_time, num_cycles);

    LDataManager* l_data_manager = ib_method_ops->getLDataManager();
    gcw = l_data_manager->getGhostCellWidth();

    Mat A = nullptr;
    ib_method_ops->updateFixedLEOperators();
    ib_method_ops->constructLagrangianForceJacobian(A, MATAIJ, data_time);

    Pointer<PatchLevel<NDIM>> fine_level = patch_hierarchy->getPatchLevel(finest_ln);
    std::vector<int> full_num_dofs_fine;
    IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        full_num_dofs_fine, u_full_dof_index_idx, p_full_dof_index_idx, fine_level);

    Mat J = nullptr;
    ib_method_ops->constructInterpOp(J, ib4_interp_fcn, 4, full_num_dofs_fine, u_full_dof_index_idx, data_time);

    Mat SAJ_full = nullptr;
    int ierr = MatPtAP(A, J, MAT_INITIAL_MATRIX, 1.0, &SAJ_full);
    IBTK_CHKERRQ(ierr);

    std::vector<std::set<int>> full_field_is;
    std::vector<std::string> full_field_names;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
        full_field_is, full_field_names, full_num_dofs_fine, u_full_dof_index_idx, p_full_dof_index_idx, fine_level);
    const auto velocity_name_it = std::find(full_field_names.begin(), full_field_names.end(), "velocity");
    if (velocity_name_it == full_field_names.end())
    {
        TBOX_ERROR("implicit_stokes_ib_coarse_saj_transfer_parity_01: velocity field not found.\n");
    }
    const std::size_t velocity_idx =
        static_cast<std::size_t>(std::distance(full_field_names.begin(), velocity_name_it));
    std::vector<PetscInt> velocity_dofs(full_field_is[velocity_idx].begin(), full_field_is[velocity_idx].end());
    IS velocity_is = nullptr;
    ierr = ISCreateGeneral(PETSC_COMM_WORLD,
                           static_cast<PetscInt>(velocity_dofs.size()),
                           velocity_dofs.empty() ? nullptr : velocity_dofs.data(),
                           PETSC_COPY_VALUES,
                           &velocity_is);
    IBTK_CHKERRQ(ierr);

    Mat SAJ_u_fine = nullptr;
    ierr = MatCreateSubMatrix(SAJ_full, velocity_is, velocity_is, MAT_INITIAL_MATRIX, &SAJ_u_fine);
    IBTK_CHKERRQ(ierr);
    ierr = ISDestroy(&velocity_is);
    IBTK_CHKERRQ(ierr);
    Vec saj_row_max_abs = nullptr;
    ierr = MatCreateVecs(SAJ_u_fine, &saj_row_max_abs, nullptr);
    IBTK_CHKERRQ(ierr);
    ierr = MatGetRowMaxAbs(SAJ_u_fine, saj_row_max_abs, nullptr);
    IBTK_CHKERRQ(ierr);
    PetscReal saj_inf_norm = 0.0;
    ierr = VecNorm(saj_row_max_abs, NORM_INFINITY, &saj_inf_norm);
    IBTK_CHKERRQ(ierr);
    const double saj_nontrivial_tol = input_db->getDoubleWithDefault("SAJ_NONTRIVIAL_TOL", 1.0e-14);
    int test_failures = 0;
    if (!(saj_inf_norm > saj_nontrivial_tol)) ++test_failures;

    std::vector<int> num_u_dofs_coarse;
    std::vector<int> num_u_dofs_fine;
    Pointer<PatchLevel<NDIM>> coarse_level = patch_hierarchy->getPatchLevel(coarsest_ln);
    IBTK::PETScVecUtilities::constructPatchLevelDOFIndices(num_u_dofs_coarse, u_vel_dof_index_idx, coarse_level);
    IBTK::PETScVecUtilities::constructPatchLevelDOFIndices(num_u_dofs_fine, u_vel_dof_index_idx, fine_level);

    AO coarse_u_ao = nullptr;
    IBTK::PETScVecUtilities::constructPatchLevelAO(
        coarse_u_ao, num_u_dofs_coarse, u_vel_dof_index_idx, coarse_level, 0);

    Mat u_prolong = nullptr;
    IBTK::PETScMatUtilities::constructProlongationOp(u_prolong,
                                                     "RT0",
                                                     u_vel_dof_index_idx,
                                                     num_u_dofs_fine,
                                                     num_u_dofs_coarse,
                                                     fine_level,
                                                     coarse_level,
                                                     coarse_u_ao,
                                                     0);

    Vec u_restriction_scale = nullptr;
    IBTK::PETScMatUtilities::constructRestrictionScalingOp(u_prolong, u_restriction_scale);

    Mat SAJ_coarse_matrix = nullptr;
    ierr = MatPtAP(SAJ_u_fine, u_prolong, MAT_INITIAL_MATRIX, 1.0, &SAJ_coarse_matrix);
    IBTK_CHKERRQ(ierr);
    ierr = MatDiagonalScale(SAJ_coarse_matrix, u_restriction_scale, nullptr);
    IBTK_CHKERRQ(ierr);

    Mat SAJ_coarse_samrai = nullptr;
    ierr = MatDuplicate(SAJ_coarse_matrix, MAT_DO_NOT_COPY_VALUES, &SAJ_coarse_samrai);
    IBTK_CHKERRQ(ierr);
    ierr = MatZeroEntries(SAJ_coarse_samrai);
    IBTK_CHKERRQ(ierr);

    build_coarse_op_from_samrai_transfers(SAJ_coarse_samrai,
                                          SAJ_u_fine,
                                          num_u_dofs_coarse,
                                          num_u_dofs_fine,
                                          u_data_var,
                                          u_data_idx,
                                          u_vel_dof_index_idx,
                                          patch_hierarchy,
                                          gcw);

    Mat diff = nullptr;
    ierr = MatDuplicate(SAJ_coarse_matrix, MAT_COPY_VALUES, &diff);
    IBTK_CHKERRQ(ierr);
    ierr = MatAXPY(diff, -1.0, SAJ_coarse_samrai, DIFFERENT_NONZERO_PATTERN);
    IBTK_CHKERRQ(ierr);
    Vec row_max_abs = nullptr;
    ierr = MatCreateVecs(diff, &row_max_abs, nullptr);
    IBTK_CHKERRQ(ierr);
    ierr = MatGetRowMaxAbs(diff, row_max_abs, nullptr);
    IBTK_CHKERRQ(ierr);
    PetscReal inf_norm = 0.0;
    ierr = VecNorm(row_max_abs, NORM_INFINITY, &inf_norm);
    IBTK_CHKERRQ(ierr);

    const double parity_tol = input_db->getDoubleWithDefault("PARITY_TOL", 1.0e-12);
    if (!(inf_norm <= parity_tol)) ++test_failures;

    pout << "saj_fine_inf_norm = " << saj_inf_norm << std::endl;
    pout << "saj_nontrivial_tol = " << saj_nontrivial_tol << std::endl;
    pout << "coarse_saj_matrix_vs_samrai_inf_norm = " << inf_norm << std::endl;
    pout << "coarse_saj_parity_tol = " << parity_tol << std::endl;
    pout << "test_failures = " << test_failures << std::endl;

    ierr = MatDestroy(&diff);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&row_max_abs);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&SAJ_coarse_samrai);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&SAJ_coarse_matrix);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&u_restriction_scale);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&u_prolong);
    IBTK_CHKERRQ(ierr);
    ierr = AODestroy(&coarse_u_ao);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&SAJ_u_fine);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&saj_row_max_abs);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&SAJ_full);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&J);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    IBTK_CHKERRQ(ierr);

    ib_method_ops->postprocessIntegrateData(current_time, new_time, num_cycles);

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(u_current_idx);
        level->deallocatePatchData(u_full_dof_index_idx);
        level->deallocatePatchData(p_full_dof_index_idx);
        level->deallocatePatchData(u_vel_dof_index_idx);
        level->deallocatePatchData(u_data_idx);
    }

    plog << "Input database:\n";
    input_db->printClassData(plog);
    return test_failures;
}
