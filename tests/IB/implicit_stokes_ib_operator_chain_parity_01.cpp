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

#include <ibamr/IBImplicitStaggeredHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
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
#include <IntVector.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <PoissonSpecifications.h>
#include <SAMRAI_config.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
// Match the implicit_ib_coupling_aware_vanka test.m defaults.
constexpr double L_DOMAIN = 1.0;
constexpr int N_GRID = 32;
constexpr double R_CYL = 0.25;
constexpr double ALPHA = 0.23;
constexpr double BETA = (R_CYL * R_CYL) / ALPHA;
constexpr double RHO = 1.0;
constexpr double MU = 1.0e-2;
constexpr double K_SPRING = 1.0e2;
constexpr double DX = L_DOMAIN / static_cast<double>(N_GRID);
constexpr double DT = 0.5 * DX;

int s_finest_ln = 0;

struct MatrixEntry
{
    int row = 0;
    int col = 0;
    double val = 0.0;
};

struct MatrixData
{
    int nrows = 0;
    int ncols = 0;
    std::vector<MatrixEntry> entries;
};

struct EulerianDofMapEntry
{
    char type = '?'; // 'u', 'v', or 'p'
    int axis = -1;   // 0/1 for velocity, -1 for pressure
    int i = 0;
    int j = 0;
    bool set = false;
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

int
compute_num_lag_nodes()
{
    const double ds0 = DX / R_CYL;
    const int approx = static_cast<int>(std::llround((2.0 * M_PI) / ds0));
    return approx;
}

void
generate_structure(const unsigned int& strct_num,
                   const int& ln,
                   int& num_vertices,
                   std::vector<IBTK::Point>& vertex_posn,
                   void* /*ctx*/)
{
    if (ln != s_finest_ln || strct_num != 0)
    {
        num_vertices = 0;
        vertex_posn.clear();
        return;
    }

    const int n_nodes = compute_num_lag_nodes();
    num_vertices = n_nodes;
    vertex_posn.resize(num_vertices);
    for (int k = 0; k < num_vertices; ++k)
    {
        const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_vertices);
        vertex_posn[k](0) = 0.5 + ALPHA * std::cos(theta);
        vertex_posn[k](1) = 0.5 + BETA * std::sin(theta);
    }
}

void
generate_springs(
    const unsigned int& strct_num,
    const int& ln,
    std::multimap<int, IBRedundantInitializer::Edge>& spring_map,
    std::map<IBRedundantInitializer::Edge, IBRedundantInitializer::SpringSpec, IBRedundantInitializer::EdgeComp>&
        spring_spec,
    void* /*ctx*/)
{
    if (ln != s_finest_ln || strct_num != 0) return;

    const int n_nodes = compute_num_lag_nodes();
    const double ds = (2.0 * M_PI) / static_cast<double>(n_nodes);
    const double spring_k = K_SPRING / (ds * ds);
    for (int k = 0; k < n_nodes; ++k)
    {
        IBRedundantInitializer::Edge edge = { k, (k + 1) % n_nodes };
        if (edge.first > edge.second) std::swap(edge.first, edge.second);
        spring_map.insert(std::make_pair(edge.first, edge));

        IBRedundantInitializer::SpringSpec spec_data;
        spec_data.force_fcn_idx = 0;
        spec_data.parameters.resize(2);
        spec_data.parameters[0] = spring_k;
        spec_data.parameters[1] = 0.0;
        spring_spec.insert(std::make_pair(edge, spec_data));
    }
}

MatrixData
extract_matrix_data(Mat mat)
{
    MatrixData data;
    PetscInt nrows = 0, ncols = 0;
    int ierr = MatGetSize(mat, &nrows, &ncols);
    IBTK_CHKERRQ(ierr);
    data.nrows = static_cast<int>(nrows);
    data.ncols = static_cast<int>(ncols);

    PetscInt row_start = 0, row_end = 0;
    ierr = MatGetOwnershipRange(mat, &row_start, &row_end);
    IBTK_CHKERRQ(ierr);

    for (PetscInt row = row_start; row < row_end; ++row)
    {
        PetscInt ncols_row = 0;
        const PetscInt* cols = nullptr;
        const PetscScalar* vals = nullptr;
        ierr = MatGetRow(mat, row, &ncols_row, &cols, &vals);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < ncols_row; ++k)
        {
            if (std::abs(static_cast<double>(vals[k])) <= 1.0e-13) continue;
            MatrixEntry e;
            e.row = static_cast<int>(row);
            e.col = static_cast<int>(cols[k]);
            e.val = static_cast<double>(vals[k]);
            data.entries.push_back(e);
        }
        ierr = MatRestoreRow(mat, row, &ncols_row, &cols, &vals);
        IBTK_CHKERRQ(ierr);
    }

    std::sort(data.entries.begin(),
              data.entries.end(),
              [](const MatrixEntry& a, const MatrixEntry& b)
              {
                  if (a.row != b.row) return a.row < b.row;
                  return a.col < b.col;
              });

    return data;
}

std::vector<double>
extract_vector_data(Vec vec)
{
    PetscInt n_global = 0;
    PetscInt n_local = 0;
    int ierr = VecGetSize(vec, &n_global);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetLocalSize(vec, &n_local);
    IBTK_CHKERRQ(ierr);

    PetscInt row_start = 0, row_end = 0;
    ierr = VecGetOwnershipRange(vec, &row_start, &row_end);
    IBTK_CHKERRQ(ierr);

    std::vector<double> vals(static_cast<std::size_t>(n_global), 0.0);
    const PetscScalar* arr = nullptr;
    ierr = VecGetArrayRead(vec, &arr);
    IBTK_CHKERRQ(ierr);
    for (PetscInt i = 0; i < n_local; ++i)
    {
        vals[static_cast<std::size_t>(row_start + i)] = static_cast<double>(arr[i]);
    }
    ierr = VecRestoreArrayRead(vec, &arr);
    IBTK_CHKERRQ(ierr);

    if (IBTK_MPI::getNodes() > 1)
    {
        TBOX_ERROR("This test currently supports serial execution only.\n");
    }

    if (row_start != 0 || row_end != n_global)
    {
        TBOX_ERROR("Unexpected non-contiguous local ownership in serial vector extraction.\n");
    }

    return vals;
}

void
write_matrix_ref(const std::string& filename, const MatrixData& data)
{
    std::ofstream os(filename);
    os << data.nrows << " " << data.ncols << " " << data.entries.size() << "\n";
    os << std::setprecision(17);
    for (const auto& e : data.entries)
    {
        os << e.row << " " << e.col << " " << e.val << "\n";
    }
}

void
write_vector_ref(const std::string& filename, const std::vector<double>& data)
{
    std::ofstream os(filename);
    os << data.size() << "\n";
    os << std::setprecision(17);
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        os << i << " " << data[i] << "\n";
    }
}

bool
read_matrix_ref(const std::string& filename, MatrixData& data)
{
    std::ifstream is(filename);
    if (!is) return false;

    std::size_t nnz = 0;
    is >> data.nrows >> data.ncols >> nnz;
    data.entries.resize(nnz);
    for (std::size_t k = 0; k < nnz; ++k)
    {
        is >> data.entries[k].row >> data.entries[k].col >> data.entries[k].val;
    }
    return true;
}

bool
read_vector_ref(const std::string& filename, std::vector<double>& data)
{
    std::ifstream is(filename);
    if (!is) return false;

    std::size_t n = 0;
    is >> n;
    data.resize(n);
    for (std::size_t k = 0; k < n; ++k)
    {
        std::size_t idx = 0;
        is >> idx >> data[k];
    }
    return true;
}

bool
compare_matrix(const std::string& name,
               const MatrixData& got,
               const MatrixData& ref,
               const double tol,
               std::ostream& os,
               int& test_failures)
{
    if (got.nrows != ref.nrows || got.ncols != ref.ncols)
    {
        ++test_failures;
        os << name << " dimension mismatch: got (" << got.nrows << "," << got.ncols << ") ref (" << ref.nrows << ","
           << ref.ncols << ")\n";
        return false;
    }

    std::vector<MatrixEntry> got_entries, ref_entries;
    got_entries.reserve(got.entries.size());
    ref_entries.reserve(ref.entries.size());
    for (const auto& e : got.entries)
    {
        if (std::abs(e.val) > tol) got_entries.push_back(e);
    }
    for (const auto& e : ref.entries)
    {
        if (std::abs(e.val) > tol) ref_entries.push_back(e);
    }

    if (got_entries.size() != ref_entries.size())
    {
        ++test_failures;
        os << name << " nnz mismatch after tol-filtering: got " << got_entries.size() << " ref " << ref_entries.size()
           << " (tol=" << tol << ")\n";
        return false;
    }

    double max_abs_err = 0.0;
    for (std::size_t k = 0; k < got_entries.size(); ++k)
    {
        if (got_entries[k].row != ref_entries[k].row || got_entries[k].col != ref_entries[k].col)
        {
            ++test_failures;
            os << name << " sparsity mismatch at entry " << k << ": got (" << got_entries[k].row << ","
               << got_entries[k].col << ") ref (" << ref_entries[k].row << "," << ref_entries[k].col << ")\n";
            return false;
        }
        const double err = std::abs(got_entries[k].val - ref_entries[k].val);
        const double scale = std::max(std::abs(ref_entries[k].val), 1.0);
        max_abs_err = std::max(max_abs_err, err);
        if (err > tol * (scale + 1.0))
        {
            ++test_failures;
            os << name << " value mismatch at entry " << k << ": err=" << err << " ref=" << ref_entries[k].val
               << " scale=" << scale << " tol=" << tol << " (criterion err <= tol*(scale+1))\n";
            return false;
        }
    }

    return true;
}

bool
compare_vector(const std::string& name,
               const std::vector<double>& got,
               const std::vector<double>& ref,
               const double tol,
               std::ostream& os,
               int& test_failures)
{
    if (got.size() != ref.size())
    {
        ++test_failures;
        os << name << " size mismatch: got " << got.size() << " ref " << ref.size() << "\n";
        return false;
    }

    double max_abs_err = 0.0;
    for (std::size_t k = 0; k < got.size(); ++k)
    {
        const double err = std::abs(got[k] - ref[k]);
        const double scale = std::max(std::abs(ref[k]), 1.0);
        if (err > tol * (scale + 1.0))
        {
            ++test_failures;
            os << name << " value mismatch: err=" << err << " at idx=" << k << " got=" << got[k] << " ref=" << ref[k]
               << " scale=" << scale << " tol=" << tol << " (criterion err <= tol*(scale+1))\n";
            return false;
        }
        if (err > max_abs_err)
        {
            max_abs_err = err;
        }
    }

    return true;
}

bool
check_lagrangian_a_row_nnz(const MatrixData& A_data, const double nz_tol, std::ostream& os, int& test_failures)
{
    if (A_data.nrows != A_data.ncols)
    {
        ++test_failures;
        os << "A structural check failed: matrix is not square (" << A_data.nrows << "x" << A_data.ncols << ")\n";
        return false;
    }

    std::vector<int> row_nnz(static_cast<std::size_t>(A_data.nrows), 0);
    for (const auto& e : A_data.entries)
    {
        if (std::abs(e.val) > nz_tol) ++row_nnz[static_cast<std::size_t>(e.row)];
    }

    for (int row = 0; row < A_data.nrows; ++row)
    {
        if (row_nnz[static_cast<std::size_t>(row)] != 3)
        {
            ++test_failures;
            os << "A structural check failed: row " << row << " has " << row_nnz[static_cast<std::size_t>(row)]
               << " clearly nonzero entries; expected 3 (nz_tol=" << nz_tol << ")\n";
            return false;
        }
    }

    return true;
}

std::string
ref_path(const std::string& ref_dir, const std::string& base)
{
    return ref_dir + "/" + base + ".ref";
}

std::string
ref_path_level(const std::string& ref_dir, const std::string& base, const int ln)
{
    return ref_dir + "/" + base + ".level" + std::to_string(ln) + ".ref";
}

std::vector<EulerianDofMapEntry>
build_eulerian_dof_map(const int total_dofs,
                       const int u_dof_index_idx,
                       const int p_dof_index_idx,
                       Pointer<PatchLevel<NDIM>> level)
{
    std::vector<EulerianDofMapEntry> dof_map(static_cast<std::size_t>(total_dofs));

    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);

        for (int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM>& side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
            for (Box<NDIM>::Iterator b(side_box); b; b++)
            {
                const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
                const int dof_idx = (*u_dof_data)(i_s);
                if (dof_idx < 0) continue;
                if (dof_idx >= total_dofs)
                {
                    TBOX_ERROR("Eulerian DOF index out of range while building map.\n");
                }

                EulerianDofMapEntry& e = dof_map[static_cast<std::size_t>(dof_idx)];
                if (e.set)
                {
                    continue;
                }
                e.type = 'u';
                e.axis = axis;
                e.i = i_s(0);
                e.j = i_s(1);
                e.set = true;
            }
        }

        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& ic = b();
            const int dof_idx = (*p_dof_data)(ic);
            if (dof_idx < 0) continue;
            if (dof_idx >= total_dofs)
            {
                TBOX_ERROR("Eulerian DOF index out of range while building map.\n");
            }

            EulerianDofMapEntry& e = dof_map[static_cast<std::size_t>(dof_idx)];
            if (e.set)
            {
                continue;
            }
            e.type = 'p';
            e.axis = -1;
            e.i = ic(0);
            e.j = ic(1);
            e.set = true;
        }
    }

    for (int dof_idx = 0; dof_idx < total_dofs; ++dof_idx)
    {
        if (!dof_map[static_cast<std::size_t>(dof_idx)].set)
        {
            TBOX_ERROR("Missing Eulerian DOF map entry while exporting reference data.\n");
        }
    }

    return dof_map;
}

void
write_eulerian_dof_map_ref(const std::string& filename, const std::vector<EulerianDofMapEntry>& dof_map)
{
    std::ofstream os(filename);
    if (!os)
    {
        TBOX_ERROR("Unable to open Eulerian DOF map reference file for writing: " << filename << "\n");
    }
    os << dof_map.size() << "\n";
    for (std::size_t k = 0; k < dof_map.size(); ++k)
    {
        const EulerianDofMapEntry& e = dof_map[k];
        os << k << " " << e.type << " " << e.axis << " " << e.i << " " << e.j << "\n";
    }
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR("implicit_stokes_ib_operator_chain_parity_01 requires serial execution (np=1).\n");
    }

#ifndef IBTK_HAVE_SILO
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
#endif

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    Pointer<Database> test_db =
        input_db->isDatabase("test") ? input_db->getDatabase("test") : Pointer<Database>(input_db, false);

    const bool write_reference = test_db->getBoolWithDefault("write_reference", false);
    const bool check_all_levels = test_db->getBoolWithDefault("check_all_levels", false);
    const bool compare_base_chain = test_db->getBoolWithDefault("compare_base_chain", true);
    const double parity_tol = test_db->getDoubleWithDefault("parity_tol", 1.0e-12);
    const std::string ref_dir_cfg =
        test_db->getStringWithDefault("reference_dir", "reference-generated/single_level_case_01");
    const std::filesystem::path ref_dir_path = std::filesystem::path(ref_dir_cfg).is_absolute() ?
                                                   std::filesystem::path(ref_dir_cfg) :
                                                   std::filesystem::path(SOURCE_DIR) / ref_dir_cfg;
    const std::string ref_dir = ref_dir_path.string();
    if (write_reference)
    {
        std::filesystem::create_directories(ref_dir);
    }

    if (input_db->keyExists("petsc_options_file"))
    {
        const std::string petsc_options_file = input_db->getString("petsc_options_file");
        PetscOptionsInsertFile(PETSC_COMM_WORLD, nullptr, petsc_options_file.c_str(), PETSC_TRUE);
    }

    Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
        "INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
    Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
    ib_method_ops->setUseFixedLEOperators(true);

    Pointer<IBHierarchyIntegrator> time_integrator =
        new IBImplicitStaggeredHierarchyIntegrator("IBHierarchyIntegrator",
                                                   app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                                   ib_method_ops,
                                                   navier_stokes_integrator);

    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", time_integrator, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
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
    s_finest_ln = input_db->getInteger("MAX_LEVELS") - 1;
    ib_initializer->setStructureNamesOnLevel(s_finest_ln, { "parity_curve2d" });
    ib_initializer->registerInitStructureFunction(generate_structure);
    ib_initializer->registerInitSpringDataFunction(generate_springs);
    ib_method_ops->registerLInitStrategy(ib_initializer);
    Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
    ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);

    time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
    const double current_time = time_integrator->getIntegratorTime();
    const double new_time = current_time + time_integrator->getMaximumTimeStepSize();
    const int num_cycles = time_integrator->getNumberOfCycles();
    time_integrator->preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(finest_ln);
    const double data_time = new_time;

    LDataManager* l_data_manager = ib_method_ops->getLDataManager();
    Pointer<IBTK::LData> X_data = l_data_manager->getLData(IBTK::LDataManager::POSN_DATA_NAME, finest_ln);
    Vec X_vec = X_data->getVec();

    Mat A = nullptr;
    ib_method_ops->updateFixedLEOperators();
    ib_method_ops->constructLagrangianForceJacobian(A, MATAIJ, data_time);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("operator_chain_parity_ctx");
    Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("operator_chain_u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof_index_var = new CellVariable<NDIM, int>("operator_chain_p_dof");
    const int u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, ctx, IntVector<NDIM>(1));
    const int p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, ctx, IntVector<NDIM>(1));
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> ln_level = patch_hierarchy->getPatchLevel(ln);
        ln_level->allocatePatchData(u_dof_index_idx, data_time);
        ln_level->allocatePatchData(p_dof_index_idx, data_time);
    }

    std::vector<int> num_dofs_per_proc;
    IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);
    const int total_dofs = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);

    Mat J = nullptr;
    ib_method_ops->constructInterpOp(J, ib4_interp_fcn, 4, num_dofs_per_proc, u_dof_index_idx, data_time);

    Mat S = nullptr;
    int ierr = MatTranspose(J, MAT_INITIAL_MATRIX, &S);
    IBTK_CHKERRQ(ierr);

    Mat SAJ = nullptr;
    ierr = MatPtAP(A, J, MAT_INITIAL_MATRIX, 1.0, &SAJ);
    IBTK_CHKERRQ(ierr);

    SAMRAI::solv::PoissonSpecifications U_problem_coefs("U_problem_coefs");
    U_problem_coefs.setCConstant(RHO / DT);
    U_problem_coefs.setDConstant(-MU);
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);

    Mat stokes_mat = nullptr;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(
        stokes_mat, U_problem_coefs, u_bc_coefs, data_time, num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);

    Mat A00_full = nullptr;
    ierr = MatDuplicate(stokes_mat, MAT_COPY_VALUES, &A00_full);
    IBTK_CHKERRQ(ierr);
    ierr = MatAXPY(A00_full, 1.0, SAJ, DIFFERENT_NONZERO_PATTERN);
    IBTK_CHKERRQ(ierr);

    Mat A00 = nullptr;
    IBAMR::StaggeredStokesPETScMatUtilities::constructA00VelocitySubmatrix(
        A00, A00_full, num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);

    Vec F_vec = nullptr;
    ierr = VecDuplicate(X_vec, &F_vec);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(A, X_vec, F_vec);
    IBTK_CHKERRQ(ierr);

    const std::vector<double> X_vals = extract_vector_data(X_vec);
    const std::vector<double> F_vals = extract_vector_data(F_vec);
    const MatrixData A_data = extract_matrix_data(A);
    const MatrixData S_data = extract_matrix_data(S);
    const MatrixData J_data = extract_matrix_data(J);
    const MatrixData SAJ_data = extract_matrix_data(SAJ);
    const MatrixData A00_data = extract_matrix_data(A00);
    std::vector<MatrixData> SAJ_level_data(static_cast<std::size_t>(finest_ln + 1));
    std::vector<MatrixData> A00_level_data(static_cast<std::size_t>(finest_ln + 1));

    if (check_all_levels)
    {
        std::vector<std::vector<int>> level_num_dofs_per_proc(static_cast<std::size_t>(finest_ln + 1));
        std::vector<std::vector<int>> level_num_u_dofs_per_proc(static_cast<std::size_t>(finest_ln + 1));
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> ln_level = patch_hierarchy->getPatchLevel(ln);
            IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
                level_num_dofs_per_proc[static_cast<std::size_t>(ln)], u_dof_index_idx, p_dof_index_idx, ln_level);
            IBTK::PETScVecUtilities::constructPatchLevelDOFIndices(
                level_num_u_dofs_per_proc[static_cast<std::size_t>(ln)], u_dof_index_idx, ln_level);
        }

        std::vector<Mat> SAJ_u_level(static_cast<std::size_t>(finest_ln + 1), nullptr);
        IBAMR::StaggeredStokesPETScMatUtilities::constructA00VelocitySubmatrix(
            SAJ_u_level[static_cast<std::size_t>(finest_ln)],
            SAJ,
            level_num_dofs_per_proc[static_cast<std::size_t>(finest_ln)],
            u_dof_index_idx,
            p_dof_index_idx,
            patch_hierarchy->getPatchLevel(finest_ln));
        for (int ln = finest_ln - 1; ln >= 0; --ln)
        {
            Pointer<PatchLevel<NDIM>> fine_level = patch_hierarchy->getPatchLevel(ln + 1);
            Pointer<PatchLevel<NDIM>> coarse_level = patch_hierarchy->getPatchLevel(ln);

            AO coarse_u_ao = nullptr;
            IBTK::PETScVecUtilities::constructPatchLevelAO(
                coarse_u_ao, level_num_u_dofs_per_proc[static_cast<std::size_t>(ln)], u_dof_index_idx, coarse_level, 0);
            Mat u_prolong = nullptr;
            IBTK::PETScMatUtilities::constructProlongationOp(
                u_prolong,
                "RT0",
                u_dof_index_idx,
                level_num_u_dofs_per_proc[static_cast<std::size_t>(ln + 1)],
                level_num_u_dofs_per_proc[static_cast<std::size_t>(ln)],
                fine_level,
                coarse_level,
                coarse_u_ao,
                0);
            Vec u_restriction_scale = nullptr;
            IBTK::PETScMatUtilities::constructRestrictionScalingOp(u_prolong, u_restriction_scale);

            ierr = MatPtAP(SAJ_u_level[static_cast<std::size_t>(ln + 1)],
                           u_prolong,
                           MAT_INITIAL_MATRIX,
                           1.0,
                           &SAJ_u_level[static_cast<std::size_t>(ln)]);
            IBTK_CHKERRQ(ierr);
            ierr = MatDiagonalScale(SAJ_u_level[static_cast<std::size_t>(ln)], u_restriction_scale, nullptr);
            IBTK_CHKERRQ(ierr);

            ierr = VecDestroy(&u_restriction_scale);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&u_prolong);
            IBTK_CHKERRQ(ierr);
            ierr = AODestroy(&coarse_u_ao);
            IBTK_CHKERRQ(ierr);
        }

        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> ln_level = patch_hierarchy->getPatchLevel(ln);
            const std::vector<int>& ln_num_dofs_per_proc = level_num_dofs_per_proc[static_cast<std::size_t>(ln)];

            Mat ln_stokes_mat = nullptr;
            IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(ln_stokes_mat,
                                                                                    U_problem_coefs,
                                                                                    u_bc_coefs,
                                                                                    data_time,
                                                                                    ln_num_dofs_per_proc,
                                                                                    u_dof_index_idx,
                                                                                    p_dof_index_idx,
                                                                                    ln_level);
            Mat ln_stokes_A00 = nullptr;
            IBAMR::StaggeredStokesPETScMatUtilities::constructA00VelocitySubmatrix(
                ln_stokes_A00, ln_stokes_mat, ln_num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, ln_level);
            Mat ln_A00 = nullptr;
            ierr = MatDuplicate(ln_stokes_A00, MAT_COPY_VALUES, &ln_A00);
            IBTK_CHKERRQ(ierr);
            ierr = MatAXPY(ln_A00, 1.0, SAJ_u_level[static_cast<std::size_t>(ln)], DIFFERENT_NONZERO_PATTERN);
            IBTK_CHKERRQ(ierr);

            SAJ_level_data[static_cast<std::size_t>(ln)] =
                extract_matrix_data(SAJ_u_level[static_cast<std::size_t>(ln)]);
            A00_level_data[static_cast<std::size_t>(ln)] = extract_matrix_data(ln_A00);

            if (write_reference)
            {
                const int ln_total_dofs = std::accumulate(ln_num_dofs_per_proc.begin(), ln_num_dofs_per_proc.end(), 0);
                const std::vector<EulerianDofMapEntry> ln_dof_map =
                    build_eulerian_dof_map(ln_total_dofs, u_dof_index_idx, p_dof_index_idx, ln_level);
                write_eulerian_dof_map_ref(ref_path_level(ref_dir, "eulerian_dof_map", ln), ln_dof_map);
            }

            ierr = MatDestroy(&ln_A00);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&ln_stokes_A00);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&ln_stokes_mat);
            IBTK_CHKERRQ(ierr);
        }

        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            ierr = MatDestroy(&SAJ_u_level[static_cast<std::size_t>(ln)]);
            IBTK_CHKERRQ(ierr);
        }
    }

    int test_failures = 0;
    std::ostringstream diag;
    const double a_nz_tol = std::max(1.0e-12, 100.0 * parity_tol);
    check_lagrangian_a_row_nnz(A_data, a_nz_tol, diag, test_failures);
    if (write_reference)
    {
        const std::vector<EulerianDofMapEntry> eulerian_dof_map =
            build_eulerian_dof_map(total_dofs, u_dof_index_idx, p_dof_index_idx, level);
        write_eulerian_dof_map_ref(ref_path(ref_dir, "eulerian_dof_map"), eulerian_dof_map);
        if (compare_base_chain)
        {
            write_vector_ref(ref_path(ref_dir, "X"), X_vals);
            write_matrix_ref(ref_path(ref_dir, "A"), A_data);
            write_vector_ref(ref_path(ref_dir, "F"), F_vals);
            write_matrix_ref(ref_path(ref_dir, "S"), S_data);
            write_matrix_ref(ref_path(ref_dir, "J"), J_data);
            write_matrix_ref(ref_path(ref_dir, "SAJ"), SAJ_data);
            write_matrix_ref(ref_path(ref_dir, "A00"), A00_data);
        }
        if (check_all_levels)
        {
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                write_matrix_ref(ref_path_level(ref_dir, "SAJ", ln), SAJ_level_data[static_cast<std::size_t>(ln)]);
                write_matrix_ref(ref_path_level(ref_dir, "A00", ln), A00_level_data[static_cast<std::size_t>(ln)]);
            }
        }
    }
    else
    {
        if (compare_base_chain)
        {
            MatrixData A_ref, S_ref, J_ref, SAJ_ref, A00_ref;
            std::vector<double> X_ref, F_ref;

            if (!read_vector_ref(ref_path(ref_dir, "X"), X_ref))
            {
                ++test_failures;
                diag << "missing reference file: " << ref_path(ref_dir, "X") << "\n";
            }
            if (!read_matrix_ref(ref_path(ref_dir, "A"), A_ref))
            {
                ++test_failures;
                diag << "missing reference file: " << ref_path(ref_dir, "A") << "\n";
            }
            if (!read_vector_ref(ref_path(ref_dir, "F"), F_ref))
            {
                ++test_failures;
                diag << "missing reference file: " << ref_path(ref_dir, "F") << "\n";
            }
            if (!read_matrix_ref(ref_path(ref_dir, "S"), S_ref))
            {
                ++test_failures;
                diag << "missing reference file: " << ref_path(ref_dir, "S") << "\n";
            }
            if (!read_matrix_ref(ref_path(ref_dir, "J"), J_ref))
            {
                ++test_failures;
                diag << "missing reference file: " << ref_path(ref_dir, "J") << "\n";
            }
            if (!read_matrix_ref(ref_path(ref_dir, "SAJ"), SAJ_ref))
            {
                ++test_failures;
                diag << "missing reference file: " << ref_path(ref_dir, "SAJ") << "\n";
            }
            if (!read_matrix_ref(ref_path(ref_dir, "A00"), A00_ref))
            {
                ++test_failures;
                diag << "missing reference file: " << ref_path(ref_dir, "A00") << "\n";
            }

            if (test_failures == 0)
            {
                compare_vector("X", X_vals, X_ref, parity_tol, diag, test_failures);
                compare_matrix("A", A_data, A_ref, parity_tol, diag, test_failures);
                compare_vector("F", F_vals, F_ref, parity_tol, diag, test_failures);
                compare_matrix("S", S_data, S_ref, parity_tol, diag, test_failures);
                compare_matrix("J", J_data, J_ref, parity_tol, diag, test_failures);
                compare_matrix("SAJ", SAJ_data, SAJ_ref, parity_tol, diag, test_failures);
                compare_matrix("A00", A00_data, A00_ref, parity_tol, diag, test_failures);
            }
        }
        if (test_failures == 0 && check_all_levels)
        {
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                MatrixData SAJ_ln_ref, A00_ln_ref;
                if (!read_matrix_ref(ref_path_level(ref_dir, "SAJ", ln), SAJ_ln_ref))
                {
                    ++test_failures;
                    diag << "missing reference file: " << ref_path_level(ref_dir, "SAJ", ln) << "\n";
                    continue;
                }
                if (!read_matrix_ref(ref_path_level(ref_dir, "A00", ln), A00_ln_ref))
                {
                    ++test_failures;
                    diag << "missing reference file: " << ref_path_level(ref_dir, "A00", ln) << "\n";
                    continue;
                }
                compare_matrix("SAJ.level" + std::to_string(ln),
                               SAJ_level_data[static_cast<std::size_t>(ln)],
                               SAJ_ln_ref,
                               parity_tol,
                               diag,
                               test_failures);
                compare_matrix("A00.level" + std::to_string(ln),
                               A00_level_data[static_cast<std::size_t>(ln)],
                               A00_ln_ref,
                               parity_tol,
                               diag,
                               test_failures);
            }
        }
    }

    if (test_failures > 0)
    {
        pout << diag.str();
    }

    ierr = VecDestroy(&F_vec);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&A00);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&A00_full);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&stokes_mat);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&SAJ);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&S);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&J);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    IBTK_CHKERRQ(ierr);

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> ln_level = patch_hierarchy->getPatchLevel(ln);
        ln_level->deallocatePatchData(u_dof_index_idx);
        ln_level->deallocatePatchData(p_dof_index_idx);
    }

    std::ofstream out("output");
    out << "test_failures = " << test_failures << "\n";
    return test_failures > 0 ? 1 : 0;
}
