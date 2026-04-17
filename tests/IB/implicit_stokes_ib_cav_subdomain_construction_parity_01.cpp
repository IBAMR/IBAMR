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

#include <petscao.h>
#include <petscmat.h>

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

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
constexpr double L_DOMAIN = 1.0;
constexpr int N_GRID = 32;
constexpr double R_CYL = 0.25;
constexpr double ALPHA = 0.23;
constexpr double BETA = (R_CYL * R_CYL) / ALPHA;
constexpr double RHO = 1.0;
constexpr double MU = 1.0e-2;
constexpr double K_SPRING_DEFAULT = 1.0e2;
constexpr double DX = L_DOMAIN / static_cast<double>(N_GRID);
constexpr double DT = 0.5 * DX;

int s_finest_ln = 0;
double s_spring_stiffness = K_SPRING_DEFAULT;

std::string
level_suffix(const int ln)
{
    return ".level" + std::to_string(ln);
}

void
write_matrix_market(const std::string& path, Mat mat)
{
    PetscInt nrows = 0, ncols = 0;
    int ierr = MatGetSize(mat, &nrows, &ncols);
    IBTK_CHKERRQ(ierr);

    std::vector<std::tuple<int, int, double>> entries;
    for (PetscInt i = 0; i < nrows; ++i)
    {
        PetscInt ncols_row = 0;
        const PetscInt* cols = nullptr;
        const PetscScalar* vals = nullptr;
        ierr = MatGetRow(mat, i, &ncols_row, &cols, &vals);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < ncols_row; ++k)
        {
            entries.emplace_back(static_cast<int>(i), static_cast<int>(cols[k]), static_cast<double>(vals[k]));
        }
        ierr = MatRestoreRow(mat, i, &ncols_row, &cols, &vals);
        IBTK_CHKERRQ(ierr);
    }
    std::sort(entries.begin(), entries.end());

    std::ofstream os(path);
    os << "%%MatrixMarket matrix coordinate real general\n";
    os << nrows << " " << ncols << " " << entries.size() << "\n";
    os << std::setprecision(17);
    for (const auto& [i, j, val] : entries)
    {
        os << (i + 1) << " " << (j + 1) << " " << val << "\n";
    }
}

void
write_matrix_row_pattern(const std::string& path, Mat mat)
{
    PetscInt nrows = 0, ncols = 0;
    int ierr = MatGetSize(mat, &nrows, &ncols);
    IBTK_CHKERRQ(ierr);

    std::ofstream os(path);
    os << nrows << "\n";
    for (PetscInt i = 0; i < nrows; ++i)
    {
        PetscInt ncols_row = 0;
        const PetscInt* cols = nullptr;
        ierr = MatGetRow(mat, i, &ncols_row, &cols, nullptr);
        IBTK_CHKERRQ(ierr);
        os << i << " " << ncols_row;
        for (PetscInt k = 0; k < ncols_row; ++k) os << " " << cols[k];
        os << "\n";
        ierr = MatRestoreRow(mat, i, &ncols_row, &cols, nullptr);
        IBTK_CHKERRQ(ierr);
    }
}

void
write_int_vector(const std::string& path, const std::vector<int>& vals)
{
    std::ofstream os(path);
    os << vals.size() << "\n";
    for (const int v : vals) os << v << "\n";
}

void
write_set_vector(const std::string& path, const std::vector<std::set<int>>& sets)
{
    std::ofstream os(path);
    os << sets.size() << "\n";
    for (std::size_t k = 0; k < sets.size(); ++k)
    {
        os << k << " " << sets[k].size();
        for (const int dof : sets[k]) os << " " << dof;
        os << "\n";
    }
}

void
write_axis_map(const std::string& path, const std::unordered_map<int, int>& map_data)
{
    std::ofstream os(path);
    os << map_data.size() << "\n";
    std::vector<std::pair<int, int>> ordered_map(map_data.begin(), map_data.end());
    std::sort(ordered_map.begin(), ordered_map.end());
    for (const auto& kv : ordered_map) os << kv.first << " " << kv.second << "\n";
}

void
write_set_map(const std::string& path, const std::unordered_map<int, std::vector<int>>& map_data)
{
    std::ofstream os(path);
    os << map_data.size() << "\n";
    std::vector<std::pair<int, std::vector<int>>> ordered_map(map_data.begin(), map_data.end());
    std::sort(
        ordered_map.begin(), ordered_map.end(), [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });
    for (const auto& kv : ordered_map)
    {
        os << kv.first << " " << kv.second.size();
        for (const int v : kv.second) os << " " << v;
        os << "\n";
    }
}

std::string
set_signature(const std::set<int>& values)
{
    std::ostringstream os;
    bool first = true;
    for (const int value : values)
    {
        if (!first) os << ",";
        os << value;
        first = false;
    }
    return os.str();
}

std::vector<std::string>
subdomain_signatures(const std::vector<std::set<int>>& subdomains)
{
    std::vector<std::string> signatures;
    signatures.reserve(subdomains.size());
    for (const auto& subdomain : subdomains) signatures.push_back(set_signature(subdomain));
    std::sort(signatures.begin(), signatures.end());
    return signatures;
}

std::vector<int>
compute_ordered_seed_velocity_dofs(Pointer<PatchLevel<NDIM>> patch_level,
                                   const int seed_velocity_axis,
                                   const int u_dof_index_idx,
                                   const std::unordered_map<int, int>& velocity_dof_to_component_axis)
{
    std::set<int> axis_velocity_dofs;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Box<NDIM> side_patch_box = SideGeometry<NDIM>::toSideBox(patch_box, seed_velocity_axis);
        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        for (Box<NDIM>::Iterator b(side_patch_box); b; b++)
        {
            const SideIndex<NDIM> i_s(b(), seed_velocity_axis, SideIndex<NDIM>::Lower);
            const int dof = (*u_dof_data)(i_s);
            if (dof < 0) continue;
            const auto axis_it = velocity_dof_to_component_axis.find(dof);
            if (axis_it == velocity_dof_to_component_axis.end()) continue;
            if (axis_it->second != seed_velocity_axis) continue;
            axis_velocity_dofs.insert(dof);
        }
    }

    std::vector<int> ordered_seeds;
    ordered_seeds.reserve(axis_velocity_dofs.size());
    for (const int dof : axis_velocity_dofs) ordered_seeds.push_back(dof);
    return ordered_seeds;
}

using SubmatrixTriplet = std::tuple<int, int, double>;

std::vector<SubmatrixTriplet>
compute_submatrix_triplets(Mat A00_mat, const std::set<int>& dofs)
{
    PetscInt nrows = 0, ncols_mat = 0;
    int ierr = MatGetSize(A00_mat, &nrows, &ncols_mat);
    IBTK_CHKERRQ(ierr);

    std::vector<int> sorted_dofs;
    sorted_dofs.reserve(dofs.size());
    for (const int dof : dofs)
    {
        if (0 <= dof && dof < static_cast<int>(nrows) && dof < static_cast<int>(ncols_mat))
        {
            sorted_dofs.push_back(dof);
        }
    }
    std::unordered_set<int> dof_lookup(sorted_dofs.begin(), sorted_dofs.end());
    std::vector<SubmatrixTriplet> entries;
    for (const int row : sorted_dofs)
    {
        PetscInt ncols = 0;
        const PetscInt* cols = nullptr;
        const PetscScalar* vals = nullptr;
        ierr = MatGetRow(A00_mat, static_cast<PetscInt>(row), &ncols, &cols, &vals);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < ncols; ++k)
        {
            const int col = static_cast<int>(cols[k]);
            if (dof_lookup.find(col) == dof_lookup.end()) continue;
            entries.emplace_back(row, col, static_cast<double>(vals[k]));
        }
        ierr = MatRestoreRow(A00_mat, static_cast<PetscInt>(row), &ncols, &cols, &vals);
        IBTK_CHKERRQ(ierr);
    }
    std::sort(entries.begin(), entries.end());
    return entries;
}

void
write_overlap_submatrix_blocks(const std::string& path, Mat A00_mat, const std::vector<std::set<int>>& overlap_sets)
{
    PetscInt nrows = 0, ncols = 0;
    int ierr = MatGetSize(A00_mat, &nrows, &ncols);
    IBTK_CHKERRQ(ierr);

    std::ofstream os(path);
    os << overlap_sets.size() << "\n";
    os << std::setprecision(17);
    for (std::size_t k = 0; k < overlap_sets.size(); ++k)
    {
        std::set<int> dof_set;
        for (const int dof : overlap_sets[k])
        {
            if (0 <= dof && dof < static_cast<int>(nrows) && dof < static_cast<int>(ncols))
            {
                dof_set.insert(dof);
            }
        }
        const std::vector<SubmatrixTriplet> entries = compute_submatrix_triplets(A00_mat, dof_set);
        os << "subdomain " << k << " " << dof_set.size() << " " << entries.size() << "\n";
        os << "dofs " << dof_set.size();
        for (const int dof : dof_set) os << " " << dof;
        os << "\n";
        for (const auto& e : entries)
        {
            os << std::get<0>(e) << " " << std::get<1>(e) << " " << std::get<2>(e) << "\n";
        }
    }
}

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
    return static_cast<int>(std::llround((2.0 * M_PI) / ds0));
}

void
generate_structure(const unsigned int& strct_num,
                   const int& ln,
                   int& num_vertices,
                   std::vector<IBTK::Point>& vertex_posn,
                   void*)
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
    void*)
{
    if (ln != s_finest_ln || strct_num != 0) return;

    const int n_nodes = compute_num_lag_nodes();
    const double ds = (2.0 * M_PI) / static_cast<double>(n_nodes);
    const double spring_k = s_spring_stiffness / (ds * ds);
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

std::set<int>
expand_velocity_dofs_from_seed_components(const std::set<int>& seed_velocity_components,
                                          Mat A00_mat,
                                          const bool use_structural_coupling = false)
{
    std::set<int> initial_velocity_dofs = seed_velocity_components;
    for (const int velocity_dof : seed_velocity_components)
    {
        PetscInt ncols = 0;
        const PetscInt* cols = nullptr;
        const PetscScalar* vals = nullptr;
        const PetscInt row = static_cast<PetscInt>(velocity_dof);
        int ierr = MatGetRow(A00_mat, row, &ncols, &cols, &vals);
        IBTK_CHKERRQ(ierr);

        double row_max_abs = 0.0;
        for (PetscInt k = 0; k < ncols; ++k)
        {
            row_max_abs = std::max(row_max_abs, static_cast<double>(PetscAbsScalar(vals[k])));
        }
        const double numerical_zero_tol =
            std::max(static_cast<double>(ncols) * std::numeric_limits<double>::epsilon() * row_max_abs,
                     IBTK_RELATIVE_NUMERICAL_ZERO_TOL * row_max_abs);

        for (PetscInt k = 0; k < ncols; ++k)
        {
            const double value = PetscRealPart(vals[k]);
            if (use_structural_coupling || std::abs(value) > numerical_zero_tol)
            {
                initial_velocity_dofs.insert(static_cast<int>(cols[k]));
            }
        }

        ierr = MatRestoreRow(A00_mat, row, &ncols, &cols, &vals);
        IBTK_CHKERRQ(ierr);
    }
    return initial_velocity_dofs;
}

std::set<int>
matlab_extract_coupled_dofs_relaxed(const int seed_velocity_dof,
                                    Mat A00_mat,
                                    const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
                                    const std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs,
                                    const std::unordered_map<int, int>& velocity_dof_to_component_axis)
{
    const std::set<int> initial_velocity_dofs =
        expand_velocity_dofs_from_seed_components({ seed_velocity_dof }, A00_mat);

    std::set<int> involved_cells;
    for (const int vel_dof : initial_velocity_dofs)
    {
        const auto it = velocity_dof_to_adjacent_cell_dofs.find(vel_dof);
        if (it == velocity_dof_to_adjacent_cell_dofs.end()) continue;
        involved_cells.insert(it->second.begin(), it->second.end());
    }

    std::set<int> closure;
    for (const int cell_dof : involved_cells)
    {
        const auto it = cell_dof_to_closure_dofs.find(cell_dof);
        if (it == cell_dof_to_closure_dofs.end()) continue;
        closure.insert(it->second.begin(), it->second.end());
    }

    std::set<int> out = closure;
    for (const int dof : initial_velocity_dofs)
    {
        if (velocity_dof_to_component_axis.find(dof) != velocity_dof_to_component_axis.end()) out.insert(dof);
    }
    return out;
}

std::set<int>
matlab_extract_coupled_dofs_strict(const int seed_velocity_dof,
                                   Mat A00_mat,
                                   const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
                                   const std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs,
                                   const std::unordered_map<int, int>& velocity_dof_to_component_axis,
                                   const std::unordered_map<int, std::vector<int>>& velocity_dof_to_paired_seed_velocity_dofs)
{
    std::set<int> initial_seed_components = { seed_velocity_dof };
    const auto pair_it = velocity_dof_to_paired_seed_velocity_dofs.find(seed_velocity_dof);
    if (pair_it != velocity_dof_to_paired_seed_velocity_dofs.end())
    {
        initial_seed_components.insert(pair_it->second.begin(), pair_it->second.end());
    }

    const std::set<int> initial_velocity_dofs =
        expand_velocity_dofs_from_seed_components(initial_seed_components, A00_mat);

    std::set<int> candidate_cells;
    for (const int vel_dof : initial_velocity_dofs)
    {
        const auto it = velocity_dof_to_adjacent_cell_dofs.find(vel_dof);
        if (it == velocity_dof_to_adjacent_cell_dofs.end()) continue;
        candidate_cells.insert(it->second.begin(), it->second.end());
    }

    std::set<int> accepted_cells;
    for (const int cell_dof : candidate_cells)
    {
        const auto closure_it = cell_dof_to_closure_dofs.find(cell_dof);
        if (closure_it == cell_dof_to_closure_dofs.end()) continue;
        bool valid = true;
        for (const int dof : closure_it->second)
        {
            if (velocity_dof_to_component_axis.find(dof) == velocity_dof_to_component_axis.end()) continue;
            if (initial_velocity_dofs.find(dof) == initial_velocity_dofs.end())
            {
                valid = false;
                break;
            }
        }
        if (valid) accepted_cells.insert(cell_dof);
    }

    std::set<int> out;
    for (const int cell_dof : accepted_cells)
    {
        const auto it = cell_dof_to_closure_dofs.find(cell_dof);
        if (it == cell_dof_to_closure_dofs.end()) continue;
        out.insert(it->second.begin(), it->second.end());
    }
    return out;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR("implicit_stokes_ib_cav_subdomain_construction_parity_01 requires serial execution (np=1).\n");
    }

#ifndef IBTK_HAVE_SILO
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
#endif

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    Pointer<Database> test_db =
        input_db->isDatabase("test") ? input_db->getDatabase("test") : Pointer<Database>(input_db, false);
    bool export_bridge_data = false;
    std::string export_bridge_dir = "bridge-generated";
    if (test_db->keyExists("export_bridge_data")) export_bridge_data = test_db->getBool("export_bridge_data");
    if (test_db->keyExists("export_bridge_dir")) export_bridge_dir = test_db->getString("export_bridge_dir");
    if (input_db->keyExists("export_bridge_data")) export_bridge_data = input_db->getBool("export_bridge_data");
    if (input_db->keyExists("export_bridge_dir")) export_bridge_dir = input_db->getString("export_bridge_dir");
    s_spring_stiffness = input_db->getDoubleWithDefault("SPRING_STIFFNESS", K_SPRING_DEFAULT);
    if (!(s_spring_stiffness > 0.0))
    {
        TBOX_ERROR("implicit_stokes_ib_cav_subdomain_construction_parity_01 requires SPRING_STIFFNESS > 0.0.\n");
    }
    if (export_bridge_data) std::filesystem::create_directories(export_bridge_dir);

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

    Mat A = nullptr;
    ib_method_ops->updateFixedLEOperators();
    ib_method_ops->constructLagrangianForceJacobian(A, MATAIJ, data_time);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("cav_subdomain_parity_ctx");
    Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("cav_parity_u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof_index_var = new CellVariable<NDIM, int>("cav_parity_p_dof");
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

    Mat J = nullptr;
    ib_method_ops->constructInterpOp(J, ib4_interp_fcn, 4, num_dofs_per_proc, u_dof_index_idx, data_time);

    Mat SAJ = nullptr;
    int ierr = MatPtAP(A, J, MAT_INITIAL_MATRIX, 1.0, &SAJ);
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

    IBAMR::StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData map_data;
    IBAMR::StaggeredStokesPETScMatUtilities::buildPatchLevelCellClosureMaps(
        map_data, u_dof_index_idx, p_dof_index_idx, level);

    std::vector<std::set<int>> overlap_relaxed, nonoverlap_relaxed;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains(
        overlap_relaxed,
        nonoverlap_relaxed,
        num_dofs_per_proc,
        u_dof_index_idx,
        level,
        Pointer<CoarseFineBoundary<NDIM>>(nullptr),
        A00_full,
        map_data,
        0,
        1,
        IBAMR::CouplingAwareASMSeedTraversalOrder::I_J,
        IBAMR::CouplingAwareASMClosurePolicy::RELAXED);

    std::vector<std::set<int>> overlap_strict, nonoverlap_strict;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains(
        overlap_strict,
        nonoverlap_strict,
        num_dofs_per_proc,
        u_dof_index_idx,
        level,
        Pointer<CoarseFineBoundary<NDIM>>(nullptr),
        A00_full,
        map_data,
        0,
        1,
        IBAMR::CouplingAwareASMSeedTraversalOrder::I_J,
        IBAMR::CouplingAwareASMClosurePolicy::STRICT);

    int test_failures = 0;
    if (overlap_relaxed.size() != overlap_strict.size())
    {
        ++test_failures;
    }

    const std::vector<int> ordered_seeds =
        compute_ordered_seed_velocity_dofs(level, 0, u_dof_index_idx, map_data.velocity_dof_to_component_axis);

    if (overlap_relaxed.size() != ordered_seeds.size() || overlap_strict.size() != ordered_seeds.size())
    {
        ++test_failures;
    }

    std::vector<std::set<int>> matlab_relaxed_subdomains;
    matlab_relaxed_subdomains.reserve(ordered_seeds.size());
    for (const int seed : ordered_seeds)
    {
        matlab_relaxed_subdomains.push_back(matlab_extract_coupled_dofs_relaxed(seed,
                                                                                 A00_full,
                                                                                 map_data.velocity_dof_to_adjacent_cell_dofs,
                                                                                 map_data.cell_dof_to_closure_dofs,
                                                                                 map_data.velocity_dof_to_component_axis));
    }
    const bool relaxed_matches = subdomain_signatures(overlap_relaxed) == subdomain_signatures(matlab_relaxed_subdomains);
    if (!relaxed_matches) ++test_failures;

    std::vector<std::set<int>> matlab_strict_subdomains;
    matlab_strict_subdomains.reserve(ordered_seeds.size());
    for (const int seed : ordered_seeds)
    {
        matlab_strict_subdomains.push_back(matlab_extract_coupled_dofs_strict(
            seed,
            A00_full,
            map_data.velocity_dof_to_adjacent_cell_dofs,
            map_data.cell_dof_to_closure_dofs,
            map_data.velocity_dof_to_component_axis,
            map_data.velocity_dof_to_paired_seed_velocity_dofs));
    }
    const bool strict_matches = subdomain_signatures(overlap_strict) == subdomain_signatures(matlab_strict_subdomains);
    if (!strict_matches) ++test_failures;

    if (export_bridge_data)
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

            IBAMR::StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData ln_map_data;
            IBAMR::StaggeredStokesPETScMatUtilities::buildPatchLevelCellClosureMaps(
                ln_map_data, u_dof_index_idx, p_dof_index_idx, ln_level);

            std::vector<std::set<int>> ln_overlap_relaxed, ln_nonoverlap_relaxed;
            IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains(
                ln_overlap_relaxed,
                ln_nonoverlap_relaxed,
                ln_num_dofs_per_proc,
                u_dof_index_idx,
                ln_level,
                Pointer<CoarseFineBoundary<NDIM>>(nullptr),
                ln_A00,
                ln_map_data,
                0,
                1,
                IBAMR::CouplingAwareASMSeedTraversalOrder::I_J,
                IBAMR::CouplingAwareASMClosurePolicy::RELAXED);

            std::vector<std::set<int>> ln_overlap_strict, ln_nonoverlap_strict;
            IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains(
                ln_overlap_strict,
                ln_nonoverlap_strict,
                ln_num_dofs_per_proc,
                u_dof_index_idx,
                ln_level,
                Pointer<CoarseFineBoundary<NDIM>>(nullptr),
                ln_A00,
                ln_map_data,
                0,
                1,
                IBAMR::CouplingAwareASMSeedTraversalOrder::I_J,
                IBAMR::CouplingAwareASMClosurePolicy::STRICT);

            const std::vector<int> ln_ordered_seeds = compute_ordered_seed_velocity_dofs(
                ln_level, 0, u_dof_index_idx, ln_map_data.velocity_dof_to_component_axis);
            const std::string suffix = level_suffix(ln);
            write_matrix_market(export_bridge_dir + "/A00" + suffix + ".mtx", ln_A00);
            write_matrix_row_pattern(export_bridge_dir + "/A00" + suffix + ".row_pattern.txt", ln_A00);
            write_int_vector(export_bridge_dir + "/seed_velocity_dofs_axis0" + suffix + ".txt", ln_ordered_seeds);
            write_int_vector(export_bridge_dir + "/num_dofs_per_proc" + suffix + ".txt", ln_num_dofs_per_proc);
            write_set_vector(export_bridge_dir + "/ibamr_overlap_relaxed" + suffix + ".txt", ln_overlap_relaxed);
            write_set_vector(export_bridge_dir + "/ibamr_overlap_strict" + suffix + ".txt", ln_overlap_strict);
            write_set_vector(export_bridge_dir + "/ibamr_nonoverlap_relaxed" + suffix + ".txt", ln_nonoverlap_relaxed);
            write_set_vector(export_bridge_dir + "/ibamr_nonoverlap_strict" + suffix + ".txt", ln_nonoverlap_strict);
            write_axis_map(export_bridge_dir + "/velocity_dof_to_component_axis" + suffix + ".txt",
                           ln_map_data.velocity_dof_to_component_axis);
            write_set_map(export_bridge_dir + "/velocity_dof_to_adjacent_cell_dofs" + suffix + ".txt",
                          ln_map_data.velocity_dof_to_adjacent_cell_dofs);
            write_set_map(export_bridge_dir + "/cell_dof_to_closure_dofs" + suffix + ".txt",
                          ln_map_data.cell_dof_to_closure_dofs);
            write_set_map(export_bridge_dir + "/velocity_dof_to_paired_seed_velocity_dofs" + suffix + ".txt",
                          ln_map_data.velocity_dof_to_paired_seed_velocity_dofs);
            write_overlap_submatrix_blocks(
                export_bridge_dir + "/ibamr_overlap_submat_relaxed" + suffix + ".txt", ln_A00, ln_overlap_relaxed);
            write_overlap_submatrix_blocks(
                export_bridge_dir + "/ibamr_overlap_submat_strict" + suffix + ".txt", ln_A00, ln_overlap_strict);

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
        pout << "wrote_bridge_data_dir = " << export_bridge_dir << "\n";
    }

    ierr = MatDestroy(&A00_full);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&stokes_mat);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&SAJ);
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

    plog << "Input database:\n";
    input_db->printClassData(plog);
    pout << "test_failures = " << test_failures << "\n";
    return test_failures > 0 ? 1 : 0;
}
