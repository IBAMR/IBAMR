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
#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/LMesh.h>
#include <ibtk/LNode.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <petscviewer.h>

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
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
struct StructureSpec
{
    int num_curve_points = 64;
    double ds = 1.0 / 64.0;
    double elasticity_k = 1.0e8;
    double x_center = 0.5;
    double y_center = 0.5;
    double x_radius = 0.2;
    double y_radius = 0.2;
    double spring_stiffness = 0.0;
    int finest_ln = 0;
    bool use_theta_spacing = false;
};

struct ParityAuditConfig
{
    bool enabled = false;
    std::string output_dir = "/tmp/ibamr_cav_live_parity";
    std::string case_id = "default_case";
    std::string closure_policy = "RELAXED";
    bool export_preconditioned_operator = true;
    std::string shell_pc_type = "multiplicative-eigen-reference";
    bool has_eigen_subdomain_solver_type = false;
    std::string eigen_subdomain_solver_type;
    bool has_eigen_subdomain_solver_threshold = false;
    double eigen_subdomain_solver_threshold = 0.0;
    bool has_a00_solver_type = false;
    std::string a00_solver_type;
    bool has_a00_solver_threshold = false;
    double a00_solver_threshold = 0.0;
    bool has_schur_solver_type = false;
    std::string schur_solver_type;
    bool has_schur_solver_threshold = false;
    double schur_solver_threshold = 0.0;
    bool has_blas_lapack_subdomain_solver_type = false;
    std::string blas_lapack_subdomain_solver_type;
    bool has_blas_lapack_subdomain_solver_rcond = false;
    double blas_lapack_subdomain_solver_rcond = 0.0;
};

struct MarkerRecord
{
    int marker_index = -1;
    double x = 0.0;
    double y = 0.0;
    double fx = 0.0;
    double fy = 0.0;
};

struct DofRecord
{
    int dof = -1;
    std::string kind;
    int axis = -1;
    int i = 0;
    int j = 0;
    double x = 0.0;
    double y = 0.0;
};

std::string
join_path(const std::string& a, const std::string& b)
{
    if (a.empty()) return b;
    if (b.empty()) return a;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
}

std::string
json_escape(const std::string& s)
{
    std::ostringstream os;
    for (char c : s)
    {
        switch (c)
        {
        case '\\':
            os << "\\\\";
            break;
        case '"':
            os << "\\\"";
            break;
        case '\n':
            os << "\\n";
            break;
        case '\r':
            os << "\\r";
            break;
        case '\t':
            os << "\\t";
            break;
        default:
            os << c;
            break;
        }
    }
    return os.str();
}

std::ofstream
open_output_file(const std::string& filename)
{
    std::ofstream out(filename.c_str());
    if (!out)
    {
        TBOX_ERROR("unable to open output file: " << filename << "\n");
    }
    out << std::scientific << std::setprecision(17);
    return out;
}

void
write_matrix_market(const std::string& filename, Mat mat)
{
    PetscInt nrows = 0;
    PetscInt ncols = 0;
    PetscErrorCode ierr = MatGetSize(mat, &nrows, &ncols);
    IBTK_CHKERRQ(ierr);

    PetscInt nnz = 0;
    for (PetscInt i = 0; i < nrows; ++i)
    {
        PetscInt row_nnz = 0;
        const PetscInt* row_cols = nullptr;
        const PetscScalar* row_vals = nullptr;
        ierr = MatGetRow(mat, i, &row_nnz, &row_cols, &row_vals);
        IBTK_CHKERRQ(ierr);
        nnz += row_nnz;
        ierr = MatRestoreRow(mat, i, &row_nnz, &row_cols, &row_vals);
        IBTK_CHKERRQ(ierr);
    }

    std::ofstream out = open_output_file(filename);
    out << "%%MatrixMarket matrix coordinate real general\n";
    out << nrows << " " << ncols << " " << nnz << "\n";
    for (PetscInt i = 0; i < nrows; ++i)
    {
        PetscInt row_nnz = 0;
        const PetscInt* row_cols = nullptr;
        const PetscScalar* row_vals = nullptr;
        ierr = MatGetRow(mat, i, &row_nnz, &row_cols, &row_vals);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < row_nnz; ++k)
        {
            out << (i + 1) << " " << (row_cols[k] + 1) << " " << PetscRealPart(row_vals[k]) << "\n";
        }
        ierr = MatRestoreRow(mat, i, &row_nnz, &row_cols, &row_vals);
        IBTK_CHKERRQ(ierr);
    }
}

double
stage_d_test_vector_value(const DofRecord& rec)
{
    const double i = static_cast<double>(rec.i + 1);
    const double j = static_cast<double>(rec.j + 1);
    if (rec.kind == "velocity")
    {
        if (rec.axis == 0)
        {
            return std::sin(0.31 * i) + std::cos(0.17 * j);
        }
        if (rec.axis == 1)
        {
            return -0.5 * std::cos(0.29 * i) + 0.75 * std::sin(0.13 * j);
        }
    }
    return 0.25 * std::sin(0.11 * i + 0.07 * j) + 0.1;
}

void
write_vector_matrix_market(const std::string& filename, Vec vec)
{
    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR("write_vector_matrix_market() requires serial execution.\n");
    }

    PetscInt n = 0;
    PetscErrorCode ierr = VecGetSize(vec, &n);
    IBTK_CHKERRQ(ierr);
    const PetscScalar* vals = nullptr;
    ierr = VecGetArrayRead(vec, &vals);
    IBTK_CHKERRQ(ierr);

    PetscInt nnz = 0;
    for (PetscInt i = 0; i < n; ++i)
    {
        if (PetscRealPart(vals[i]) != 0.0) ++nnz;
    }

    std::ofstream out = open_output_file(filename);
    out << "%%MatrixMarket matrix coordinate real general\n";
    out << n << " 1 " << nnz << "\n";
    for (PetscInt i = 0; i < n; ++i)
    {
        const double value = PetscRealPart(vals[i]);
        if (value == 0.0) continue;
        out << (i + 1) << " 1 " << value << "\n";
    }

    ierr = VecRestoreArrayRead(vec, &vals);
    IBTK_CHKERRQ(ierr);
}

std::string
uppercase(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return s;
}

ParityAuditConfig
parse_parity_audit_config(Pointer<Database> input_db)
{
    ParityAuditConfig cfg;
    if (!input_db->isDatabase("parity_audit_db")) return cfg;

    Pointer<Database> parity_db = input_db->getDatabase("parity_audit_db");
    cfg.enabled = parity_db->getBoolWithDefault("enabled", false);
    cfg.output_dir = parity_db->getStringWithDefault("output_dir", cfg.output_dir);
    cfg.case_id = parity_db->getStringWithDefault("case_id", cfg.case_id);
    cfg.closure_policy = uppercase(parity_db->getStringWithDefault("closure_policy", cfg.closure_policy));
    cfg.export_preconditioned_operator =
        parity_db->getBoolWithDefault("export_preconditioned_operator", cfg.export_preconditioned_operator);
    cfg.shell_pc_type = parity_db->getStringWithDefault("shell_pc_type", cfg.shell_pc_type);
    if (parity_db->keyExists("eigen_subdomain_solver_type"))
    {
        cfg.has_eigen_subdomain_solver_type = true;
        cfg.eigen_subdomain_solver_type = parity_db->getString("eigen_subdomain_solver_type");
    }
    if (parity_db->keyExists("eigen_subdomain_solver_threshold"))
    {
        cfg.has_eigen_subdomain_solver_threshold = true;
        cfg.eigen_subdomain_solver_threshold = parity_db->getDouble("eigen_subdomain_solver_threshold");
    }
    if (parity_db->keyExists("a00_solver_type"))
    {
        cfg.has_a00_solver_type = true;
        cfg.a00_solver_type = parity_db->getString("a00_solver_type");
    }
    if (parity_db->keyExists("a00_solver_threshold"))
    {
        cfg.has_a00_solver_threshold = true;
        cfg.a00_solver_threshold = parity_db->getDouble("a00_solver_threshold");
    }
    if (parity_db->keyExists("schur_solver_type"))
    {
        cfg.has_schur_solver_type = true;
        cfg.schur_solver_type = parity_db->getString("schur_solver_type");
    }
    if (parity_db->keyExists("schur_solver_threshold"))
    {
        cfg.has_schur_solver_threshold = true;
        cfg.schur_solver_threshold = parity_db->getDouble("schur_solver_threshold");
    }
    if (parity_db->keyExists("blas_lapack_subdomain_solver_type"))
    {
        cfg.has_blas_lapack_subdomain_solver_type = true;
        cfg.blas_lapack_subdomain_solver_type = parity_db->getString("blas_lapack_subdomain_solver_type");
    }
    if (parity_db->keyExists("blas_lapack_subdomain_solver_rcond"))
    {
        cfg.has_blas_lapack_subdomain_solver_rcond = true;
        cfg.blas_lapack_subdomain_solver_rcond = parity_db->getDouble("blas_lapack_subdomain_solver_rcond");
    }
    if (cfg.closure_policy != "RELAXED" && cfg.closure_policy != "STRICT")
    {
        TBOX_ERROR("parity_audit_db.closure_policy must be RELAXED or STRICT\n");
    }
    return cfg;
}

void
enforce_parity_preconditioner_setup(Pointer<Database> stokes_ib_precond_db, const ParityAuditConfig& parity_cfg)
{
    Pointer<Database> level_solver_db = stokes_ib_precond_db->isDatabase("level_solver_db") ?
                                            stokes_ib_precond_db->getDatabase("level_solver_db") :
                                            nullptr;
    Pointer<Database> coarse_solver_db = stokes_ib_precond_db->isDatabase("coarse_solver_db") ?
                                             stokes_ib_precond_db->getDatabase("coarse_solver_db") :
                                             nullptr;
    if (!level_solver_db || !coarse_solver_db)
    {
        TBOX_ERROR("parity audit requires level_solver_db and coarse_solver_db in stokes_ib_precond_db\n");
    }

    stokes_ib_precond_db->putBool("has_pressure_nullspace", true);
    stokes_ib_precond_db->putInteger("num_pre_sweeps", 1);
    stokes_ib_precond_db->putInteger("num_post_sweeps", 1);
    stokes_ib_precond_db->putBool("rediscretize_stokes", true);
    stokes_ib_precond_db->putBool("res_rediscretized_stokes", true);
    stokes_ib_precond_db->putString("level_solver_type", "PETSC_LEVEL_SOLVER");
    stokes_ib_precond_db->putDouble("level_solver_rel_residual_tol", 1.0e-10);
    stokes_ib_precond_db->putDouble("level_solver_abs_residual_tol", 1.0e-50);
    stokes_ib_precond_db->putInteger("level_solver_max_iterations", 1);
    stokes_ib_precond_db->putString("coarse_solver_type", "PETSC_LEVEL_SOLVER");
    stokes_ib_precond_db->putDouble("coarse_solver_rel_residual_tol", 1.0e-10);
    stokes_ib_precond_db->putDouble("coarse_solver_abs_residual_tol", 1.0e-50);
    stokes_ib_precond_db->putInteger("coarse_solver_max_iterations", 1);

    level_solver_db->putBool("initial_guess_nonzero", true);
    level_solver_db->putString("ksp_type", "richardson");
    level_solver_db->putString("pc_type", "shell");
    level_solver_db->putString("shell_pc_type", parity_cfg.shell_pc_type);
    level_solver_db->putString("asm_subdomain_construction_mode", "COUPLING_AWARE");
    level_solver_db->putString("coupling_aware_asm_closure_policy", parity_cfg.closure_policy);
    level_solver_db->putInteger("coupling_aware_asm_seed_axis", 0);
    level_solver_db->putInteger("coupling_aware_asm_seed_stride", 1);
    level_solver_db->putString("coupling_aware_asm_seed_traversal_order", "I_J");
    if (parity_cfg.has_eigen_subdomain_solver_type)
        level_solver_db->putString("eigen_subdomain_solver_type", parity_cfg.eigen_subdomain_solver_type);
    if (parity_cfg.has_eigen_subdomain_solver_threshold)
        level_solver_db->putDouble("eigen_subdomain_solver_threshold", parity_cfg.eigen_subdomain_solver_threshold);
    if (parity_cfg.has_a00_solver_type) level_solver_db->putString("a00_solver_type", parity_cfg.a00_solver_type);
    if (parity_cfg.has_a00_solver_threshold)
        level_solver_db->putDouble("a00_solver_threshold", parity_cfg.a00_solver_threshold);
    if (parity_cfg.has_schur_solver_type)
        level_solver_db->putString("schur_solver_type", parity_cfg.schur_solver_type);
    if (parity_cfg.has_schur_solver_threshold)
        level_solver_db->putDouble("schur_solver_threshold", parity_cfg.schur_solver_threshold);
    if (parity_cfg.has_blas_lapack_subdomain_solver_type)
    {
        level_solver_db->putString("blas_lapack_subdomain_solver_type", parity_cfg.blas_lapack_subdomain_solver_type);
    }
    if (parity_cfg.has_blas_lapack_subdomain_solver_rcond)
    {
        level_solver_db->putDouble("blas_lapack_subdomain_solver_rcond", parity_cfg.blas_lapack_subdomain_solver_rcond);
    }

    coarse_solver_db->putBool("initial_guess_nonzero", true);
    coarse_solver_db->putString("ksp_type", "preonly");
    coarse_solver_db->putString("pc_type", "svd");
}

std::vector<MarkerRecord>
collect_marker_records(Pointer<IBMethod> ib_method_ops, const int level_num, const double data_time)
{
    (void)data_time;
    std::vector<Pointer<IBTK::LData>>* X_data = nullptr;
    std::vector<Pointer<IBTK::LData>>* F_data = nullptr;
    bool* X_needs_ghost_fill = nullptr;
    bool* F_needs_ghost_fill = nullptr;
    ib_method_ops->getPositionData(&X_data, &X_needs_ghost_fill, IBTK::TimePoint::NEW_TIME);
    ib_method_ops->getForceData(&F_data, &F_needs_ghost_fill, IBTK::TimePoint::NEW_TIME);

    if (!X_data || !F_data || level_num >= static_cast<int>(X_data->size()) ||
        level_num >= static_cast<int>(F_data->size()) || !(*X_data)[level_num] || !(*F_data)[level_num])
    {
        ib_method_ops->getPositionData(&X_data, &X_needs_ghost_fill, IBTK::TimePoint::CURRENT_TIME);
        ib_method_ops->getForceData(&F_data, &F_needs_ghost_fill, IBTK::TimePoint::CURRENT_TIME);
    }

    if (!X_data || !F_data || level_num >= static_cast<int>(X_data->size()) ||
        level_num >= static_cast<int>(F_data->size()))
    {
        return {};
    }
    Pointer<IBTK::LData> X_level_data = (*X_data)[level_num];
    Pointer<IBTK::LData> F_level_data = (*F_data)[level_num];
    if (!X_level_data || !F_level_data) return {};

    IBTK::LDataManager* l_data_manager = ib_method_ops->getLDataManager();
    if (!l_data_manager) return {};
    auto* X_arr = X_level_data->getLocalFormVecArray();
    auto* F_arr = F_level_data->getLocalFormVecArray();
    if (!X_arr || !F_arr) return {};
    Pointer<IBTK::LMesh> l_mesh = l_data_manager ? l_data_manager->getLMesh(level_num) : nullptr;
    if (!l_mesh) return {};

    std::vector<MarkerRecord> marker_records;
    const std::vector<IBTK::LNode*>& local_nodes = l_mesh->getLocalNodes();
    marker_records.reserve(local_nodes.size());
    for (const IBTK::LNode* node : local_nodes)
    {
        if (!node) continue;
        const int local_petsc_idx = node->getLocalPETScIndex();
        MarkerRecord rec;
        rec.marker_index = node->getLagrangianIndex();
        rec.x = (*X_arr)[local_petsc_idx][0];
        rec.y = (*X_arr)[local_petsc_idx][1];
        rec.fx = (*F_arr)[local_petsc_idx][0];
        rec.fy = (*F_arr)[local_petsc_idx][1];
        marker_records.push_back(rec);
    }
    std::sort(marker_records.begin(),
              marker_records.end(),
              [](const MarkerRecord& a, const MarkerRecord& b) { return a.marker_index < b.marker_index; });
    return marker_records;
}

void
write_markers_json(const std::string& filename, const int level_num, const std::vector<MarkerRecord>& marker_records)
{
    std::ofstream out = open_output_file(filename);
    out << "{\n";
    out << "  \"level\": " << level_num << ",\n";
    out << "  \"count\": " << marker_records.size() << ",\n";
    out << "  \"markers\": [\n";
    for (std::size_t k = 0; k < marker_records.size(); ++k)
    {
        const MarkerRecord& rec = marker_records[k];
        out << "    {\"marker_index\": " << rec.marker_index << ", \"x\": " << rec.x << ", \"y\": " << rec.y
            << ", \"fx\": " << rec.fx << ", \"fy\": " << rec.fy << "}";
        if (k + 1 < marker_records.size()) out << ",";
        out << "\n";
    }
    out << "  ]\n";
    out << "}\n";
}

std::vector<DofRecord>
collect_level_dof_records(Pointer<PatchLevel<NDIM>> level, const int u_dof_index_idx, const int p_dof_index_idx)
{
    std::map<int, DofRecord> dof_map;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);
        Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
        const Box<NDIM>& patch_box = patch->getBox();
        const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
        const double* x_lower = patch_geom->getXLower();
        const double* dx = patch_geom->getDx();

        for (int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
            for (Box<NDIM>::Iterator b(side_box); b; b++)
            {
                const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
                const int dof = (*u_dof_data)(i_s);
                if (dof < 0 || dof_map.count(dof) != 0) continue;

                DofRecord rec;
                rec.dof = dof;
                rec.kind = "velocity";
                rec.axis = axis;
                rec.i = i_s(0);
                rec.j = i_s(1);
                rec.x = x_lower[0] + (static_cast<double>(i_s(0) - patch_lower(0)) + (axis == 0 ? 0.0 : 0.5)) * dx[0];
                rec.y = x_lower[1] + (static_cast<double>(i_s(1) - patch_lower(1)) + (axis == 1 ? 0.0 : 0.5)) * dx[1];
                dof_map.insert(std::make_pair(dof, rec));
            }
        }

        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            const CellIndex<NDIM>& i = b();
            const int dof = (*p_dof_data)(i);
            if (dof < 0 || dof_map.count(dof) != 0) continue;

            DofRecord rec;
            rec.dof = dof;
            rec.kind = "pressure";
            rec.axis = -1;
            rec.i = i(0);
            rec.j = i(1);
            rec.x = x_lower[0] + (static_cast<double>(i(0) - patch_lower(0)) + 0.5) * dx[0];
            rec.y = x_lower[1] + (static_cast<double>(i(1) - patch_lower(1)) + 0.5) * dx[1];
            dof_map.insert(std::make_pair(dof, rec));
        }
    }

    std::vector<DofRecord> records;
    records.reserve(dof_map.size());
    for (const auto& kv : dof_map)
    {
        records.push_back(kv.second);
    }
    return records;
}

void
write_dof_map_json(const std::string& filename, const int level_num, const std::vector<DofRecord>& dof_records)
{
    std::ofstream out = open_output_file(filename);
    out << "{\n";
    out << "  \"level\": " << level_num << ",\n";
    out << "  \"count\": " << dof_records.size() << ",\n";
    out << "  \"dofs\": [\n";
    for (std::size_t k = 0; k < dof_records.size(); ++k)
    {
        const DofRecord& rec = dof_records[k];
        out << "    {\"dof\": " << rec.dof << ", \"kind\": \"" << rec.kind << "\", \"axis\": " << rec.axis
            << ", \"i\": " << rec.i << ", \"j\": " << rec.j << ", \"x\": " << rec.x << ", \"y\": " << rec.y << "}";
        if (k + 1 < dof_records.size()) out << ",";
        out << "\n";
    }
    out << "  ]\n";
    out << "}\n";
}

void
write_subdomains_json(const std::string& filename,
                      const int level_num,
                      const std::vector<std::vector<int>>& overlap_subdomains,
                      const std::vector<std::vector<int>>& nonoverlap_subdomains,
                      const std::vector<int>& seed_velocity_dofs)
{
    std::ofstream out = open_output_file(filename);
    out << "{\n";
    out << "  \"level\": " << level_num << ",\n";
    out << "  \"num_subdomains\": " << overlap_subdomains.size() << ",\n";
    out << "  \"seed_velocity_dofs\": [";
    for (std::size_t k = 0; k < seed_velocity_dofs.size(); ++k)
    {
        out << seed_velocity_dofs[k];
        if (k + 1 < seed_velocity_dofs.size()) out << ", ";
    }
    out << "],\n";
    out << "  \"overlap\": [\n";
    for (std::size_t k = 0; k < overlap_subdomains.size(); ++k)
    {
        out << "    [";
        for (std::size_t q = 0; q < overlap_subdomains[k].size(); ++q)
        {
            out << overlap_subdomains[k][q];
            if (q + 1 < overlap_subdomains[k].size()) out << ", ";
        }
        out << "]";
        if (k + 1 < overlap_subdomains.size()) out << ",";
        out << "\n";
    }
    out << "  ],\n";
    out << "  \"nonoverlap\": [\n";
    for (std::size_t k = 0; k < nonoverlap_subdomains.size(); ++k)
    {
        out << "    [";
        for (std::size_t q = 0; q < nonoverlap_subdomains[k].size(); ++q)
        {
            out << nonoverlap_subdomains[k][q];
            if (q + 1 < nonoverlap_subdomains[k].size()) out << ", ";
        }
        out << "]";
        if (k + 1 < nonoverlap_subdomains.size()) out << ",";
        out << "\n";
    }
    out << "  ]\n";
    out << "}\n";
}

void
write_metadata_json(const std::string& filename,
                    const ParityAuditConfig& parity_cfg,
                    const StructureSpec& structure_spec,
                    const int finest_ln,
                    const double dt,
                    const double rho,
                    const double mu,
                    const std::vector<std::vector<int>>& num_dofs_per_proc)
{
    std::ofstream out = open_output_file(filename);
    out << "{\n";
    out << "  \"case_id\": \"" << json_escape(parity_cfg.case_id) << "\",\n";
    out << "  \"closure_policy\": \"" << json_escape(parity_cfg.closure_policy) << "\",\n";
    out << "  \"shell_pc_type\": \"" << json_escape(parity_cfg.shell_pc_type) << "\",\n";
    if (parity_cfg.has_eigen_subdomain_solver_type)
    {
        out << "  \"eigen_subdomain_solver_type\": \"" << json_escape(parity_cfg.eigen_subdomain_solver_type)
            << "\",\n";
    }
    if (parity_cfg.has_eigen_subdomain_solver_threshold)
    {
        out << "  \"eigen_subdomain_solver_threshold\": " << parity_cfg.eigen_subdomain_solver_threshold << ",\n";
    }
    if (parity_cfg.has_a00_solver_type)
    {
        out << "  \"a00_solver_type\": \"" << json_escape(parity_cfg.a00_solver_type) << "\",\n";
    }
    if (parity_cfg.has_a00_solver_threshold)
    {
        out << "  \"a00_solver_threshold\": " << parity_cfg.a00_solver_threshold << ",\n";
    }
    if (parity_cfg.has_schur_solver_type)
    {
        out << "  \"schur_solver_type\": \"" << json_escape(parity_cfg.schur_solver_type) << "\",\n";
    }
    if (parity_cfg.has_schur_solver_threshold)
    {
        out << "  \"schur_solver_threshold\": " << parity_cfg.schur_solver_threshold << ",\n";
    }
    if (parity_cfg.has_blas_lapack_subdomain_solver_type)
    {
        out << "  \"blas_lapack_subdomain_solver_type\": \""
            << json_escape(parity_cfg.blas_lapack_subdomain_solver_type) << "\",\n";
    }
    if (parity_cfg.has_blas_lapack_subdomain_solver_rcond)
    {
        out << "  \"blas_lapack_subdomain_solver_rcond\": " << parity_cfg.blas_lapack_subdomain_solver_rcond
            << ",\n";
    }
    out << "  \"finest_level\": " << finest_ln << ",\n";
    out << "  \"dt\": " << dt << ",\n";
    out << "  \"rho\": " << rho << ",\n";
    out << "  \"mu\": " << mu << ",\n";
    out << "  \"marker_spacing_ds\": " << structure_spec.ds << ",\n";
    out << "  \"num_curve_points\": " << structure_spec.num_curve_points << ",\n";
    out << "  \"level_global_dof_counts\": [";
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        const int n_global = std::accumulate(num_dofs_per_proc[ln].begin(), num_dofs_per_proc[ln].end(), 0);
        out << n_global;
        if (ln < finest_ln) out << ", ";
    }
    out << "]\n";
    out << "}\n";
}

void
write_level_subdomain_matrices(const std::string& output_dir,
                               const int level_num,
                               Mat level_mat,
                               const std::vector<std::vector<int>>& overlap_subdomains)
{
    for (std::size_t k = 0; k < overlap_subdomains.size(); ++k)
    {
        const std::vector<int>& dofs = overlap_subdomains[k];
        if (dofs.empty()) continue;
        std::vector<PetscInt> petsc_dofs(dofs.begin(), dofs.end());
        IS overlap_is = nullptr;
        Mat subdomain_mat = nullptr;
        PetscErrorCode ierr = ISCreateGeneral(PETSC_COMM_SELF,
                                              static_cast<PetscInt>(petsc_dofs.size()),
                                              petsc_dofs.data(),
                                              PETSC_COPY_VALUES,
                                              &overlap_is);
        IBTK_CHKERRQ(ierr);
        ierr = MatCreateSubMatrix(level_mat, overlap_is, overlap_is, MAT_INITIAL_MATRIX, &subdomain_mat);
        IBTK_CHKERRQ(ierr);
        const std::string filename =
            join_path(output_dir, "A_subdomain_level" + std::to_string(level_num) + "_k" + std::to_string(k) + ".mtx");
        write_matrix_market(filename, subdomain_mat);
        ierr = MatDestroy(&subdomain_mat);
        IBTK_CHKERRQ(ierr);
        ierr = ISDestroy(&overlap_is);
        IBTK_CHKERRQ(ierr);
    }
}

void
write_level_transfer_operators(const std::string& output_dir,
                               Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_op,
                               Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                               const int coarsest_ln,
                               const int finest_ln,
                               Pointer<HierarchySideDataOpsReal<NDIM, double>> hier_velocity_data_ops,
                               Pointer<HierarchyCellDataOpsReal<NDIM, double>> hier_pressure_data_ops,
                               const int u_dof_index_idx,
                               const int p_dof_index_idx,
                               Pointer<SAMRAIVectorReal<NDIM, double>> eul_template,
                               const std::vector<std::vector<int>>& num_dofs_per_proc)
{
    auto collect_pressure_dof_set = [p_dof_index_idx](Pointer<PatchLevel<NDIM>> level)
    {
        std::unordered_set<int> pressure_dofs;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);
            const Box<NDIM>& patch_box = patch->getBox();
            for (Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                const int dof = (*p_dof_data)(b());
                if (dof >= 0) pressure_dofs.insert(dof);
            }
        }
        return pressure_dofs;
    };

    Pointer<SAMRAIVectorReal<NDIM, double>> transfer_vec = eul_template->cloneVector("parity_transfer");
    transfer_vec->allocateVectorData();
    transfer_vec->setToScalar(0.0);
    const int transfer_u_idx = transfer_vec->getComponentDescriptorIndex(0);
    const int transfer_p_idx = transfer_vec->getComponentDescriptorIndex(1);

    for (int ln = coarsest_ln; ln < finest_ln; ++ln)
    {
        PetscErrorCode ierr = 0;
        const int n_coarse = std::accumulate(num_dofs_per_proc[ln].begin(), num_dofs_per_proc[ln].end(), 0);
        const int n_fine = std::accumulate(num_dofs_per_proc[ln + 1].begin(), num_dofs_per_proc[ln + 1].end(), 0);
        Pointer<PatchLevel<NDIM>> coarse_level = patch_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM>> fine_level = patch_hierarchy->getPatchLevel(ln + 1);
        const std::unordered_set<int> coarse_pressure_dofs = collect_pressure_dof_set(coarse_level);

        Mat P = fac_op->getProlongationOp(ln, StaggeredStokesFACPreconditionerStrategy::OperatorComponent::FULL);
        Vec L = fac_op->getCoarseningScalingOp(ln, StaggeredStokesFACPreconditionerStrategy::OperatorComponent::FULL);
        const std::string p_filename = join_path(output_dir, "P_level" + std::to_string(ln) + ".mtx");
        write_matrix_market(p_filename, P);

        Mat R_adjoint = nullptr;
        Mat Pt = nullptr;
        ierr = MatTranspose(P, MAT_INITIAL_MATRIX, &Pt);
        IBTK_CHKERRQ(ierr);
        ierr = MatDuplicate(Pt, MAT_COPY_VALUES, &R_adjoint);
        IBTK_CHKERRQ(ierr);
        ierr = MatDiagonalScale(R_adjoint, L, nullptr);
        IBTK_CHKERRQ(ierr);

        Vec fine_basis_col = nullptr;
        Vec coarse_output_col = nullptr;
        Mat R_live = nullptr;
        ierr = VecCreateMPI(PETSC_COMM_WORLD, n_fine, n_fine, &fine_basis_col);
        IBTK_CHKERRQ(ierr);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, n_coarse, n_coarse, &coarse_output_col);
        IBTK_CHKERRQ(ierr);
        ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, n_coarse, n_fine, 16, nullptr, &R_live);
        IBTK_CHKERRQ(ierr);

        std::vector<PetscInt> row_ids;
        std::vector<PetscScalar> row_vals;
        row_ids.reserve(static_cast<std::size_t>(n_coarse));
        row_vals.reserve(static_cast<std::size_t>(n_coarse));

        for (int j = 0; j < n_fine; ++j)
        {
            const PetscInt col_id = static_cast<PetscInt>(j);
            ierr = VecSet(fine_basis_col, 0.0);
            IBTK_CHKERRQ(ierr);
            ierr = VecSetValue(fine_basis_col, j, 1.0, INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyBegin(fine_basis_col);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyEnd(fine_basis_col);
            IBTK_CHKERRQ(ierr);

            hier_velocity_data_ops->setToScalar(transfer_u_idx, 0.0, false);
            hier_pressure_data_ops->setToScalar(transfer_p_idx, 0.0, false);
            StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(fine_basis_col,
                                                                    transfer_u_idx,
                                                                    u_dof_index_idx,
                                                                    transfer_p_idx,
                                                                    p_dof_index_idx,
                                                                    fine_level,
                                                                    nullptr,
                                                                    nullptr);

            fac_op->restrictResidual(*transfer_vec, *transfer_vec, ln);

            ierr = VecSet(coarse_output_col, 0.0);
            IBTK_CHKERRQ(ierr);
            StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
                coarse_output_col, transfer_u_idx, u_dof_index_idx, transfer_p_idx, p_dof_index_idx, coarse_level);
            ierr = VecAssemblyBegin(coarse_output_col);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyEnd(coarse_output_col);
            IBTK_CHKERRQ(ierr);

            const PetscScalar* out_arr = nullptr;
            ierr = VecGetArrayRead(coarse_output_col, &out_arr);
            IBTK_CHKERRQ(ierr);
            row_ids.clear();
            row_vals.clear();
            for (int i = 0; i < n_coarse; ++i)
            {
                const double val = PetscRealPart(out_arr[i]);
                if (val == 0.0) continue;
                row_ids.push_back(static_cast<PetscInt>(i));
                row_vals.push_back(static_cast<PetscScalar>(val));
            }
            if (!row_ids.empty())
            {
                ierr = MatSetValues(R_live,
                                    static_cast<PetscInt>(row_ids.size()),
                                    row_ids.data(),
                                    1,
                                    &col_id,
                                    row_vals.data(),
                                    INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
            ierr = VecRestoreArrayRead(coarse_output_col, &out_arr);
            IBTK_CHKERRQ(ierr);
        }

        ierr = MatAssemblyBegin(R_live, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(R_live, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);

        Mat R = nullptr;
        ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, n_coarse, n_fine, 16, nullptr, &R);
        IBTK_CHKERRQ(ierr);
        for (PetscInt i = 0; i < n_coarse; ++i)
        {
            const Mat row_source = coarse_pressure_dofs.count(static_cast<int>(i)) ? R_live : R_adjoint;
            PetscInt row_nnz = 0;
            const PetscInt* row_cols = nullptr;
            const PetscScalar* row_vals_ptr = nullptr;
            ierr = MatGetRow(row_source, i, &row_nnz, &row_cols, &row_vals_ptr);
            IBTK_CHKERRQ(ierr);
            if (row_nnz > 0)
            {
                ierr = MatSetValues(R, 1, &i, row_nnz, row_cols, row_vals_ptr, INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
            ierr = MatRestoreRow(row_source, i, &row_nnz, &row_cols, &row_vals_ptr);
            IBTK_CHKERRQ(ierr);
        }
        ierr = MatAssemblyBegin(R, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);

        const std::string r_filename = join_path(output_dir, "R_level" + std::to_string(ln) + ".mtx");
        write_matrix_market(r_filename, R);
        ierr = MatDestroy(&R);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&R_live);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&R_adjoint);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&Pt);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&coarse_output_col);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&fine_basis_col);
        IBTK_CHKERRQ(ierr);
    }
    transfer_vec->deallocateVectorData();
}

void
write_preconditioned_operator_apply_vectors(const std::string& preconditioned_input_filename,
                                            const std::string& preconditioned_output_filename,
                                            const std::string& preconditioned_coarse_rhs_filename,
                                            const std::string& preconditioned_coarse_correction_filename,
                                            const std::string& preconditioner_input_filename,
                                            const std::string& preconditioner_output_filename,
                                            const std::string& parity_case_dir,
                                            Pointer<StaggeredStokesIBJacobianOperator> /*jac_op*/,
                                            Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_op,
                                            Pointer<StaggeredStokesIBJacobianFACPreconditioner> fac_pc,
                                            Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                            const int finest_ln,
                                            Pointer<HierarchySideDataOpsReal<NDIM, double>> hier_velocity_data_ops,
                                            Pointer<HierarchyCellDataOpsReal<NDIM, double>> hier_pressure_data_ops,
                                            const int u_dof_index_idx,
                                            const int p_dof_index_idx,
                                            Pointer<SAMRAIVectorReal<NDIM, double>> eul_template,
                                            const std::vector<DofRecord>& dof_records)
{
    const int n_dofs = static_cast<int>(dof_records.size());
    const char* prior_first_sweep_dump_dir_env = std::getenv("IBAMR_CAV_FIRST_SWEEP_DUMP_DIR");
    const std::string prior_first_sweep_dump_dir =
        prior_first_sweep_dump_dir_env ? std::string(prior_first_sweep_dump_dir_env) : std::string();
    const int set_first_sweep_dump_dir_result =
        setenv("IBAMR_CAV_FIRST_SWEEP_DUMP_DIR", parity_case_dir.c_str(), 1);
    if (set_first_sweep_dump_dir_result != 0)
    {
        TBOX_ERROR("failed to set IBAMR_CAV_FIRST_SWEEP_DUMP_DIR for parity diagnostic export\n");
    }

    Pointer<SAMRAIVectorReal<NDIM, double>> seed = eul_template->cloneVector("parity_seed");
    Pointer<SAMRAIVectorReal<NDIM, double>> rhs = eul_template->cloneVector("parity_rhs");
    Pointer<SAMRAIVectorReal<NDIM, double>> sol = eul_template->cloneVector("parity_sol");
    seed->allocateVectorData();
    rhs->allocateVectorData();
    sol->allocateVectorData();
    seed->setToScalar(0.0);
    rhs->setToScalar(0.0);
    sol->setToScalar(0.0);

    const int seed_u_idx = seed->getComponentDescriptorIndex(0);
    const int seed_p_idx = seed->getComponentDescriptorIndex(1);
    const int rhs_u_idx = rhs->getComponentDescriptorIndex(0);
    const int rhs_p_idx = rhs->getComponentDescriptorIndex(1);
    const int sol_u_idx = sol->getComponentDescriptorIndex(0);
    const int sol_p_idx = sol->getComponentDescriptorIndex(1);
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(finest_ln);

    const int coarsest_ln = eul_template->getCoarsestLevelNumber();
    const int level_min = std::min(coarsest_ln, finest_ln);
    const int level_max = std::max(coarsest_ln, finest_ln);
    Pointer<PatchLevel<NDIM>> coarse_level = patch_hierarchy->getPatchLevel(coarsest_ln);
    const std::vector<DofRecord> coarse_dof_records =
        collect_level_dof_records(coarse_level, u_dof_index_idx, p_dof_index_idx);
    const int n_coarse_dofs = static_cast<int>(coarse_dof_records.size());
    std::vector<Pointer<PatchLevel<NDIM>>> patch_level_by_ln(level_max + 1, nullptr);
    std::vector<std::vector<DofRecord>> dof_records_by_ln(level_max + 1);
    for (int ln = level_min; ln <= level_max; ++ln)
    {
        patch_level_by_ln[ln] = patch_hierarchy->getPatchLevel(ln);
        dof_records_by_ln[ln] = collect_level_dof_records(patch_level_by_ln[ln], u_dof_index_idx, p_dof_index_idx);
    }

    Pointer<SAMRAIVectorReal<NDIM, double>> coarse_rhs_diag = eul_template->cloneVector("parity_coarse_rhs_diag");
    Pointer<SAMRAIVectorReal<NDIM, double>> coarse_sol_diag = eul_template->cloneVector("parity_coarse_sol_diag");
    coarse_rhs_diag->allocateVectorData();
    coarse_sol_diag->allocateVectorData();
    coarse_rhs_diag->setToScalar(0.0);
    coarse_sol_diag->setToScalar(0.0);
    const int coarse_rhs_u_idx = coarse_rhs_diag->getComponentDescriptorIndex(0);
    const int coarse_rhs_p_idx = coarse_rhs_diag->getComponentDescriptorIndex(1);
    const int coarse_sol_u_idx = coarse_sol_diag->getComponentDescriptorIndex(0);
    const int coarse_sol_p_idx = coarse_sol_diag->getComponentDescriptorIndex(1);
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>> pre_smooth_input_diag_by_ln(level_max + 1, nullptr);
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>> pre_smooth_output_diag_by_ln(level_max + 1, nullptr);
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>> post_smooth_input_diag_by_ln(level_max + 1, nullptr);
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>> post_smooth_output_diag_by_ln(level_max + 1, nullptr);
    std::vector<SAMRAIVectorReal<NDIM, double>*> pre_smooth_input_diag_ptr_by_ln(level_max + 1, nullptr);
    std::vector<SAMRAIVectorReal<NDIM, double>*> pre_smooth_output_diag_ptr_by_ln(level_max + 1, nullptr);
    std::vector<SAMRAIVectorReal<NDIM, double>*> post_smooth_input_diag_ptr_by_ln(level_max + 1, nullptr);
    std::vector<SAMRAIVectorReal<NDIM, double>*> post_smooth_output_diag_ptr_by_ln(level_max + 1, nullptr);
    for (int ln = level_min; ln <= level_max; ++ln)
    {
        if (ln == coarsest_ln) continue;
        pre_smooth_input_diag_by_ln[ln] =
            eul_template->cloneVector("parity_pre_smooth_input_diag_level_" + std::to_string(ln));
        pre_smooth_output_diag_by_ln[ln] =
            eul_template->cloneVector("parity_pre_smooth_output_diag_level_" + std::to_string(ln));
        post_smooth_input_diag_by_ln[ln] =
            eul_template->cloneVector("parity_post_smooth_input_diag_level_" + std::to_string(ln));
        post_smooth_output_diag_by_ln[ln] =
            eul_template->cloneVector("parity_post_smooth_output_diag_level_" + std::to_string(ln));
        pre_smooth_input_diag_by_ln[ln]->allocateVectorData();
        pre_smooth_output_diag_by_ln[ln]->allocateVectorData();
        post_smooth_input_diag_by_ln[ln]->allocateVectorData();
        post_smooth_output_diag_by_ln[ln]->allocateVectorData();
        pre_smooth_input_diag_by_ln[ln]->setToScalar(0.0);
        pre_smooth_output_diag_by_ln[ln]->setToScalar(0.0);
        post_smooth_input_diag_by_ln[ln]->setToScalar(0.0);
        post_smooth_output_diag_by_ln[ln]->setToScalar(0.0);
        pre_smooth_input_diag_ptr_by_ln[ln] = pre_smooth_input_diag_by_ln[ln].getPointer();
        pre_smooth_output_diag_ptr_by_ln[ln] = pre_smooth_output_diag_by_ln[ln].getPointer();
        post_smooth_input_diag_ptr_by_ln[ln] = post_smooth_input_diag_by_ln[ln].getPointer();
        post_smooth_output_diag_ptr_by_ln[ln] = post_smooth_output_diag_by_ln[ln].getPointer();
    }

    Vec apply_input = nullptr;
    Vec apply_rhs = nullptr;
    Vec apply_preconditioned_output = nullptr;
    Vec apply_preconditioner_output = nullptr;
    Vec apply_preconditioned_coarse_rhs = nullptr;
    Vec apply_preconditioned_coarse_correction = nullptr;
    PetscErrorCode ierr = VecCreateMPI(PETSC_COMM_WORLD, n_dofs, n_dofs, &apply_input);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, n_dofs, n_dofs, &apply_rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, n_dofs, n_dofs, &apply_preconditioned_output);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, n_dofs, n_dofs, &apply_preconditioner_output);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, n_coarse_dofs, n_coarse_dofs, &apply_preconditioned_coarse_rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, n_coarse_dofs, n_coarse_dofs, &apply_preconditioned_coarse_correction);
    IBTK_CHKERRQ(ierr);

    ierr = VecSet(apply_input, 0.0);
    IBTK_CHKERRQ(ierr);
    for (const auto& rec : dof_records)
    {
        const PetscInt dof = static_cast<PetscInt>(rec.dof);
        const double value = stage_d_test_vector_value(rec);
        ierr = VecSetValue(apply_input, dof, value, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(apply_input);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(apply_input);
    IBTK_CHKERRQ(ierr);

    hier_velocity_data_ops->setToScalar(seed_u_idx, 0.0, false);
    hier_pressure_data_ops->setToScalar(seed_p_idx, 0.0, false);
    StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(
        apply_input, seed_u_idx, u_dof_index_idx, seed_p_idx, p_dof_index_idx, level, nullptr, nullptr);

    Pointer<StaggeredStokesPETScLevelSolver> finest_level_solver =
        fac_op->getStaggeredStokesPETScLevelSolver(finest_ln);
    if (!finest_level_solver)
    {
        TBOX_ERROR("failed to access finest-level PETSc solver for parity Stage D matrix apply\n");
    }
    const KSP& finest_level_ksp = finest_level_solver->getPETScKSP();
    Mat finest_level_A = nullptr;
    ierr = KSPGetOperators(finest_level_ksp, &finest_level_A, nullptr);
    IBTK_CHKERRQ(ierr);

    ierr = VecSet(apply_rhs, 0.0);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(finest_level_A, apply_input, apply_rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(apply_rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(apply_rhs);
    IBTK_CHKERRQ(ierr);

    rhs->setToScalar(0.0);
    StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(
        apply_rhs, rhs_u_idx, u_dof_index_idx, rhs_p_idx, p_dof_index_idx, level, nullptr, nullptr);
    sol->setToScalar(0.0);
    const bool solve_success = fac_pc->solveSystemWithCoarseLevelDiagnostics(
        *sol,
        *rhs,
        coarse_sol_diag.getPointer(),
        coarse_rhs_diag.getPointer(),
        coarsest_ln,
        &pre_smooth_input_diag_ptr_by_ln,
        &pre_smooth_output_diag_ptr_by_ln,
        &post_smooth_input_diag_ptr_by_ln,
        &post_smooth_output_diag_ptr_by_ln);
    if (!solve_success)
    {
        TBOX_ERROR("failed to apply FAC preconditioner while exporting vector apply for M^{-1}A\n");
    }

    if (prior_first_sweep_dump_dir.empty())
    {
        const int unset_first_sweep_dump_dir_result = unsetenv("IBAMR_CAV_FIRST_SWEEP_DUMP_DIR");
        if (unset_first_sweep_dump_dir_result != 0)
        {
            TBOX_ERROR("failed to unset IBAMR_CAV_FIRST_SWEEP_DUMP_DIR after parity diagnostic export\n");
        }
    }
    else
    {
        const int restore_first_sweep_dump_dir_result =
            setenv("IBAMR_CAV_FIRST_SWEEP_DUMP_DIR", prior_first_sweep_dump_dir.c_str(), 1);
        if (restore_first_sweep_dump_dir_result != 0)
        {
            TBOX_ERROR("failed to restore IBAMR_CAV_FIRST_SWEEP_DUMP_DIR after parity diagnostic export\n");
        }
    }

    ierr = VecSet(apply_preconditioned_output, 0.0);
    IBTK_CHKERRQ(ierr);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        apply_preconditioned_output, sol_u_idx, u_dof_index_idx, sol_p_idx, p_dof_index_idx, level);
    ierr = VecAssemblyBegin(apply_preconditioned_output);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(apply_preconditioned_output);
    IBTK_CHKERRQ(ierr);

    ierr = VecSet(apply_preconditioned_coarse_rhs, 0.0);
    IBTK_CHKERRQ(ierr);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(apply_preconditioned_coarse_rhs,
                                                          coarse_rhs_u_idx,
                                                          u_dof_index_idx,
                                                          coarse_rhs_p_idx,
                                                          p_dof_index_idx,
                                                          coarse_level);
    ierr = VecAssemblyBegin(apply_preconditioned_coarse_rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(apply_preconditioned_coarse_rhs);
    IBTK_CHKERRQ(ierr);

    ierr = VecSet(apply_preconditioned_coarse_correction, 0.0);
    IBTK_CHKERRQ(ierr);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(apply_preconditioned_coarse_correction,
                                                          coarse_sol_u_idx,
                                                          u_dof_index_idx,
                                                          coarse_sol_p_idx,
                                                          p_dof_index_idx,
                                                          coarse_level);
    ierr = VecAssemblyBegin(apply_preconditioned_coarse_correction);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(apply_preconditioned_coarse_correction);
    IBTK_CHKERRQ(ierr);

    for (int ln = level_min; ln <= level_max; ++ln)
    {
        if (ln == coarsest_ln) continue;
        const int n_level_dofs = static_cast<int>(dof_records_by_ln[ln].size());
        Pointer<PatchLevel<NDIM>> level_ln = patch_level_by_ln[ln];
        Pointer<SAMRAIVectorReal<NDIM, double>> pre_smooth_input_diag = pre_smooth_input_diag_by_ln[ln];
        Pointer<SAMRAIVectorReal<NDIM, double>> pre_smooth_output_diag = pre_smooth_output_diag_by_ln[ln];
        Pointer<SAMRAIVectorReal<NDIM, double>> post_smooth_input_diag = post_smooth_input_diag_by_ln[ln];
        Pointer<SAMRAIVectorReal<NDIM, double>> post_smooth_output_diag = post_smooth_output_diag_by_ln[ln];
        if (!level_ln || !pre_smooth_input_diag || !pre_smooth_output_diag || !post_smooth_input_diag ||
            !post_smooth_output_diag)
        {
            TBOX_ERROR("missing smoother diagnostics state for level " << ln << "\n");
        }

        const int pre_input_u_idx = pre_smooth_input_diag->getComponentDescriptorIndex(0);
        const int pre_input_p_idx = pre_smooth_input_diag->getComponentDescriptorIndex(1);
        const int pre_output_u_idx = pre_smooth_output_diag->getComponentDescriptorIndex(0);
        const int pre_output_p_idx = pre_smooth_output_diag->getComponentDescriptorIndex(1);
        const int post_input_u_idx = post_smooth_input_diag->getComponentDescriptorIndex(0);
        const int post_input_p_idx = post_smooth_input_diag->getComponentDescriptorIndex(1);
        const int post_output_u_idx = post_smooth_output_diag->getComponentDescriptorIndex(0);
        const int post_output_p_idx = post_smooth_output_diag->getComponentDescriptorIndex(1);

        Vec apply_pre_smooth_input_level = nullptr;
        Vec apply_pre_smooth_output_level = nullptr;
        Vec apply_post_smooth_input_level = nullptr;
        Vec apply_post_smooth_output_level = nullptr;
        ierr = VecCreateMPI(PETSC_COMM_WORLD, n_level_dofs, n_level_dofs, &apply_pre_smooth_input_level);
        IBTK_CHKERRQ(ierr);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, n_level_dofs, n_level_dofs, &apply_pre_smooth_output_level);
        IBTK_CHKERRQ(ierr);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, n_level_dofs, n_level_dofs, &apply_post_smooth_input_level);
        IBTK_CHKERRQ(ierr);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, n_level_dofs, n_level_dofs, &apply_post_smooth_output_level);
        IBTK_CHKERRQ(ierr);

        ierr = VecSet(apply_pre_smooth_input_level, 0.0);
        IBTK_CHKERRQ(ierr);
        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            apply_pre_smooth_input_level, pre_input_u_idx, u_dof_index_idx, pre_input_p_idx, p_dof_index_idx, level_ln);
        ierr = VecAssemblyBegin(apply_pre_smooth_input_level);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(apply_pre_smooth_input_level);
        IBTK_CHKERRQ(ierr);

        ierr = VecSet(apply_pre_smooth_output_level, 0.0);
        IBTK_CHKERRQ(ierr);
        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            apply_pre_smooth_output_level,
            pre_output_u_idx,
            u_dof_index_idx,
            pre_output_p_idx,
            p_dof_index_idx,
            level_ln);
        ierr = VecAssemblyBegin(apply_pre_smooth_output_level);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(apply_pre_smooth_output_level);
        IBTK_CHKERRQ(ierr);

        ierr = VecSet(apply_post_smooth_input_level, 0.0);
        IBTK_CHKERRQ(ierr);
        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            apply_post_smooth_input_level,
            post_input_u_idx,
            u_dof_index_idx,
            post_input_p_idx,
            p_dof_index_idx,
            level_ln);
        ierr = VecAssemblyBegin(apply_post_smooth_input_level);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(apply_post_smooth_input_level);
        IBTK_CHKERRQ(ierr);

        ierr = VecSet(apply_post_smooth_output_level, 0.0);
        IBTK_CHKERRQ(ierr);
        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            apply_post_smooth_output_level,
            post_output_u_idx,
            u_dof_index_idx,
            post_output_p_idx,
            p_dof_index_idx,
            level_ln);
        ierr = VecAssemblyBegin(apply_post_smooth_output_level);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(apply_post_smooth_output_level);
        IBTK_CHKERRQ(ierr);

        write_vector_matrix_market(
            join_path(parity_case_dir, "preconditioned_apply_pre_smooth_input_level" + std::to_string(ln) + ".mtx"),
            apply_pre_smooth_input_level);
        write_vector_matrix_market(
            join_path(parity_case_dir, "preconditioned_apply_pre_smooth_output_level" + std::to_string(ln) + ".mtx"),
            apply_pre_smooth_output_level);
        write_vector_matrix_market(
            join_path(parity_case_dir, "preconditioned_apply_post_smooth_input_level" + std::to_string(ln) + ".mtx"),
            apply_post_smooth_input_level);
        write_vector_matrix_market(
            join_path(parity_case_dir, "preconditioned_apply_post_smooth_output_level" + std::to_string(ln) + ".mtx"),
            apply_post_smooth_output_level);

        // Backward-compatible file names currently interpreted as smoother outputs.
        write_vector_matrix_market(
            join_path(parity_case_dir, "preconditioned_apply_pre_smooth_level" + std::to_string(ln) + ".mtx"),
            apply_pre_smooth_output_level);
        write_vector_matrix_market(
            join_path(parity_case_dir, "preconditioned_apply_post_smooth_level" + std::to_string(ln) + ".mtx"),
            apply_post_smooth_output_level);

        ierr = VecDestroy(&apply_post_smooth_output_level);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&apply_post_smooth_input_level);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&apply_pre_smooth_output_level);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&apply_pre_smooth_input_level);
        IBTK_CHKERRQ(ierr);
    }

    rhs->setToScalar(0.0);
    StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(
        apply_input, rhs_u_idx, u_dof_index_idx, rhs_p_idx, p_dof_index_idx, level, nullptr, nullptr);
    sol->setToScalar(0.0);
    const bool preconditioner_only_solve_success = fac_pc->solveSystem(*sol, *rhs);
    if (!preconditioner_only_solve_success)
    {
        TBOX_ERROR("failed to apply FAC preconditioner while exporting vector apply for M^{-1}\n");
    }

    ierr = VecSet(apply_preconditioner_output, 0.0);
    IBTK_CHKERRQ(ierr);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        apply_preconditioner_output, sol_u_idx, u_dof_index_idx, sol_p_idx, p_dof_index_idx, level);
    ierr = VecAssemblyBegin(apply_preconditioner_output);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(apply_preconditioner_output);
    IBTK_CHKERRQ(ierr);

    write_vector_matrix_market(preconditioned_input_filename, apply_input);
    write_vector_matrix_market(preconditioned_output_filename, apply_preconditioned_output);
    write_vector_matrix_market(preconditioned_coarse_rhs_filename, apply_preconditioned_coarse_rhs);
    write_vector_matrix_market(preconditioned_coarse_correction_filename, apply_preconditioned_coarse_correction);
    write_vector_matrix_market(preconditioner_input_filename, apply_input);
    write_vector_matrix_market(preconditioner_output_filename, apply_preconditioner_output);

    ierr = VecDestroy(&apply_preconditioned_coarse_correction);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&apply_preconditioned_coarse_rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&apply_rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&apply_preconditioner_output);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&apply_preconditioned_output);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&apply_input);
    IBTK_CHKERRQ(ierr);
    seed->deallocateVectorData();
    for (int ln = level_min; ln <= level_max; ++ln)
    {
        if (ln == coarsest_ln) continue;
        pre_smooth_input_diag_by_ln[ln]->deallocateVectorData();
        pre_smooth_output_diag_by_ln[ln]->deallocateVectorData();
        post_smooth_input_diag_by_ln[ln]->deallocateVectorData();
        post_smooth_output_diag_by_ln[ln]->deallocateVectorData();
    }
    coarse_sol_diag->deallocateVectorData();
    coarse_rhs_diag->deallocateVectorData();
    rhs->deallocateVectorData();
    sol->deallocateVectorData();
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

    if (!path_exists(resolved_options_file))
    {
        TBOX_ERROR("could not open PETSc options file: " << petsc_options_file << "\n");
    }

    const PetscErrorCode ierr =
        PetscOptionsInsertFile(PETSC_COMM_WORLD, nullptr, resolved_options_file.c_str(), PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
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

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

#ifndef IBTK_HAVE_SILO
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
#endif

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        load_petsc_options_file(input_db, argc, argv);

        const double current_time = 0.0;
        const double dt = input_db->getDoubleWithDefault("DT", 0.005);
        const double rho = input_db->getDoubleWithDefault("RHO", 1.0);
        const double mu = input_db->getDoubleWithDefault("MU", 1.0);
        const double new_time = current_time + dt;
        const int num_linear_solves = input_db->getIntegerWithDefault("NUM_LINEAR_SOLVES", 1);
        const bool use_matrix_based_saj = input_db->getBoolWithDefault("USE_MATRIX_BASED_SAJ", true);
        const ParityAuditConfig parity_cfg = parse_parity_audit_config(input_db);
        if (parity_cfg.enabled && IBTK_MPI::getNodes() != 1)
        {
            TBOX_ERROR("parity_audit_db requires serial execution (np=1)\n");
        }

        StructureSpec structure_spec;
        const int n_coarse = input_db->getIntegerWithDefault("N", 16);
        const int ref_ratio = input_db->getIntegerWithDefault("REF_RATIO", 2);
        const int max_levels = input_db->getIntegerWithDefault("MAX_LEVELS", 1);
        const double domain_length = input_db->getDoubleWithDefault("L", 1.0);
        const double dx_fine =
            domain_length / (static_cast<double>(n_coarse) * std::pow(static_cast<double>(ref_ratio), max_levels - 1));
        structure_spec.ds = input_db->getDoubleWithDefault("DS", dx_fine);
        structure_spec.x_center = input_db->getDoubleWithDefault("X_CENTER", 0.5);
        structure_spec.y_center = input_db->getDoubleWithDefault("Y_CENTER", 0.5);
        structure_spec.x_radius = input_db->getDoubleWithDefault("X_RADIUS", 0.2);
        structure_spec.y_radius = input_db->getDoubleWithDefault("Y_RADIUS", 0.2);
        if (!input_db->keyExists("K"))
        {
            TBOX_ERROR("K must be specified in the input file\n");
        }
        structure_spec.elasticity_k = input_db->getDouble("K");
        structure_spec.use_theta_spacing = input_db->getBoolWithDefault("USE_THETA_SPACING_FOR_STRUCTURE", false);
        if (!(structure_spec.ds > 0.0)) TBOX_ERROR("DS must be positive\n");
        if (structure_spec.elasticity_k < 0.0) TBOX_ERROR("K must be nonnegative\n");
        structure_spec.spring_stiffness = structure_spec.elasticity_k / structure_spec.ds;
        if (!(structure_spec.x_radius > 0.0) || !(structure_spec.y_radius > 0.0))
        {
            TBOX_ERROR("X_RADIUS and Y_RADIUS must be positive\n");
        }

        const double a = structure_spec.x_radius;
        const double b = structure_spec.y_radius;
        const double h = std::pow((a - b) / (a + b), 2.0);
        const double circumference =
            M_PI * (a + b) * (1.0 + 3.0 * h / (10.0 + std::sqrt(std::max(0.0, 4.0 - 3.0 * h))));
        if (structure_spec.use_theta_spacing)
        {
            structure_spec.num_curve_points =
                std::max(3, static_cast<int>(std::llround((2.0 * M_PI) / structure_spec.ds)));
        }
        else
        {
            structure_spec.num_curve_points = std::max(3, static_cast<int>(circumference / structure_spec.ds));
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

        structure_spec.finest_ln = input_db->getIntegerWithDefault("MAX_LEVELS", 1) - 1;

        Pointer<IBRedundantInitializer> ib_initializer = new IBRedundantInitializer(
            "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
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
        const int p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, scratch_ctx, one_ghost);

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

        ib_method_ops->preprocessIntegrateData(current_time, new_time, 1);
        ib_method_ops->updateFixedLEOperators();
        ib_method_ops->computeLagrangianForce(new_time);

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
        if (parity_cfg.enabled) enforce_parity_preconditioner_setup(stokes_ib_precond_db, parity_cfg);
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

        std::vector<DofRecord> finest_level_dof_records;
        if (parity_cfg.enabled)
        {
            const std::string parity_case_dir =
                join_path(join_path(parity_cfg.output_dir, parity_cfg.case_id), "policy=" + parity_cfg.closure_policy);
            Utilities::recursiveMkdir(parity_case_dir);

            write_metadata_json(join_path(parity_case_dir, "metadata.json"),
                                parity_cfg,
                                structure_spec,
                                finest_ln,
                                dt,
                                rho,
                                mu,
                                num_dofs_per_proc);

            const std::vector<MarkerRecord> marker_records = collect_marker_records(ib_method_ops, finest_ln, new_time);
            write_markers_json(join_path(parity_case_dir, "markers_level" + std::to_string(finest_ln) + ".json"),
                               finest_ln,
                               marker_records);

            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
                const std::vector<DofRecord> level_dof_records =
                    collect_level_dof_records(level, u_dof_index_idx, p_dof_index_idx);
                write_dof_map_json(
                    join_path(parity_case_dir, "dof_map_level" + std::to_string(ln) + ".json"), ln, level_dof_records);

                Pointer<StaggeredStokesPETScLevelSolver> level_solver = fac_op->getStaggeredStokesPETScLevelSolver(ln);
                std::vector<std::vector<int>>* nonoverlap_subdomains = nullptr;
                std::vector<std::vector<int>>* overlap_subdomains = nullptr;
                level_solver->getASMSubdomains(&nonoverlap_subdomains, &overlap_subdomains);
                if (!nonoverlap_subdomains || !overlap_subdomains)
                {
                    TBOX_ERROR("failed to collect ASM subdomains for level " << ln << "\n");
                }

                const std::vector<int>& seed_velocity_dofs = level_solver->getCouplingAwareASMSeedVelocityDOFs();
                write_subdomains_json(join_path(parity_case_dir, "subdomains_level" + std::to_string(ln) + ".json"),
                                      ln,
                                      *overlap_subdomains,
                                      *nonoverlap_subdomains,
                                      seed_velocity_dofs);

                const KSP& level_ksp = level_solver->getPETScKSP();
                Mat level_A = nullptr;
                PetscErrorCode ierr = KSPGetOperators(level_ksp, &level_A, nullptr);
                IBTK_CHKERRQ(ierr);
                write_matrix_market(join_path(parity_case_dir, "A_level" + std::to_string(ln) + ".mtx"), level_A);
                write_level_subdomain_matrices(parity_case_dir, ln, level_A, *overlap_subdomains);

                if (ln == finest_ln)
                {
                    finest_level_dof_records = level_dof_records;
                }
            }

            write_level_transfer_operators(parity_case_dir,
                                           fac_op,
                                           patch_hierarchy,
                                           0,
                                           finest_ln,
                                           hier_velocity_data_ops,
                                           hier_pressure_data_ops,
                                           u_dof_index_idx,
                                           p_dof_index_idx,
                                           eul_sol_vec,
                                           num_dofs_per_proc);

            if (parity_cfg.export_preconditioned_operator)
            {
                write_preconditioned_operator_apply_vectors(
                    join_path(parity_case_dir, "preconditioned_apply_input_level_fine.mtx"),
                    join_path(parity_case_dir, "preconditioned_apply_output_level_fine.mtx"),
                    join_path(parity_case_dir, "preconditioned_apply_coarse_rhs_level0.mtx"),
                    join_path(parity_case_dir, "preconditioned_apply_coarse_correction_level0.mtx"),
                    join_path(parity_case_dir, "preconditioner_apply_input_level_fine.mtx"),
                    join_path(parity_case_dir, "preconditioner_apply_output_level_fine.mtx"),
                    parity_case_dir,
                    jac_op,
                    fac_op,
                    fac_pc,
                    patch_hierarchy,
                    finest_ln,
                    hier_velocity_data_ops,
                    hier_pressure_data_ops,
                    u_dof_index_idx,
                    p_dof_index_idx,
                    eul_sol_vec,
                    finest_level_dof_records);
            }

            pout << "parity_audit_export_dir = " << parity_case_dir << std::endl;
        }

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
