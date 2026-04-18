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
#include <ibamr/StaggeredStokesIBJacobianFACPreconditioner.h>
#include <ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/ibtk_utilities.h>

#include <petscblaslapack.h>
#include <petscis.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
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
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
// Match the implicit_ib_coupling_aware_vanka test.m defaults.
constexpr double L_DOMAIN = 1.0;
constexpr double R_CYL = 0.25;
constexpr double ALPHA = 0.23;
constexpr double BETA = (R_CYL * R_CYL) / ALPHA;
constexpr double RHO = 1.0;
constexpr double MU = 1.0e-2;
constexpr double K_SPRING_DEFAULT = 1.0e2;

int s_finest_ln = 0;
double s_spring_stiffness = K_SPRING_DEFAULT;
double s_domain_length = L_DOMAIN;
int s_finest_grid_cells = 32;
double s_dx = L_DOMAIN / static_cast<double>(s_finest_grid_cells);
double s_dt = 0.5 * s_dx;

enum class SpringNormalization
{
    JuliaDensity,
    IBAMRNodal
};

SpringNormalization s_spring_normalization = SpringNormalization::JuliaDensity;

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

struct SweepTraceFrame
{
    int sweep = -1;
    int subdomain = -1;
    int stage = -1; // 0 = after subdomain update, 1 = after pressure projection
    std::vector<double> y;
};

struct SweepResidualTraceFrame
{
    int sweep = -1;
    int subdomain = -1;
    int stage = -1; // currently only stage 0 is emitted
    std::vector<double> residual;
};

struct SweepDeltaTraceFrame
{
    int sweep = -1;
    int subdomain = -1;
    int stage = -1; // currently only stage 0 is emitted
    std::vector<double> delta;
};

struct VCycleTrace
{
    std::vector<double> y_pre;
    std::vector<double> residual_fine;
    std::vector<double> coarse_rhs;
    std::vector<double> coarse_sol;
    std::vector<double> fine_correction;
    std::vector<double> y_after_correction;
    std::vector<double> y_post;
};

struct LocalSolveDebugDump
{
    bool enabled = false;
    int target_sweep = -1;
    int target_subdomain = -1;
    int target_global_row = -1;
    std::string output_dir;
    bool dumped = false;
};

struct MatrixVcycleShellContext
{
    const std::vector<Mat>* level_l = nullptr;
    const std::vector<Mat>* level_restrict = nullptr;
    const std::vector<Mat>* level_prolong = nullptr;
    const std::vector<std::vector<std::set<int>>>* level_overlap = nullptr;
    const std::vector<std::vector<int>>* level_pressure_dofs = nullptr;
    int num_pre_sweeps = 1;
    int num_post_sweeps = 1;
    double alpha = 1.0;
};

struct MatrixKrylovSummary
{
    bool converged = false;
    int iterations = 0;
    int converged_reason = 0;
    double residual_norm = std::numeric_limits<double>::quiet_NaN();
    double relative_residual_norm = std::numeric_limits<double>::quiet_NaN();
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
    const double ds0 = s_dx / R_CYL;
    const int approx = static_cast<int>(std::llround((2.0 * M_PI) / ds0));
    return approx;
}

double
compute_lag_ds()
{
    return (2.0 * M_PI) / static_cast<double>(compute_num_lag_nodes());
}

double
compute_matrix_spread_scale()
{
    double cell_volume = 1.0;
    for (int d = 0; d < NDIM; ++d)
    {
        cell_volume *= s_dx;
    }
    switch (s_spring_normalization)
    {
    case SpringNormalization::JuliaDensity:
        return compute_lag_ds() / cell_volume;
    case SpringNormalization::IBAMRNodal:
        return 1.0 / cell_volume;
    default:
        TBOX_ERROR("Unsupported spring normalization mode.\n");
    }
}

double
compute_pair_spring_stiffness()
{
    const double lag_ds = compute_lag_ds();
    switch (s_spring_normalization)
    {
    case SpringNormalization::JuliaDensity:
        return s_spring_stiffness / (lag_ds * lag_ds);
    case SpringNormalization::IBAMRNodal:
        return s_spring_stiffness / lag_ds;
    default:
        TBOX_ERROR("Unsupported spring normalization mode.\n");
    }
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
        vertex_posn[k](0) = s_domain_length * (0.5 + ALPHA * std::cos(theta));
        vertex_posn[k](1) = s_domain_length * (0.5 + BETA * std::sin(theta));
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
    const double spring_k = compute_pair_spring_stiffness();
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

int
compute_finest_grid_cells(Pointer<Database> input_db)
{
    const int base_n = input_db->getInteger("N");
    const int max_levels = input_db->getInteger("MAX_LEVELS");
    const int ref_ratio = input_db->getIntegerWithDefault("REF_RATIO", 1);
    int finest_grid_cells = base_n;
    for (int ln = 1; ln < max_levels; ++ln)
    {
        finest_grid_cells *= ref_ratio;
    }
    return finest_grid_cells;
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

std::vector<double>
solve_pseudoinverse(Mat mat, Vec rhs)
{
    PetscInt nrows = 0, ncols = 0;
    int ierr = MatGetSize(mat, &nrows, &ncols);
    IBTK_CHKERRQ(ierr);
    if (nrows <= 0 || ncols <= 0)
    {
        return {};
    }

    PetscBLASInt m = 0, n = 0;
    ierr = PetscBLASIntCast(nrows, &m);
    IBTK_CHKERRQ(ierr);
    ierr = PetscBLASIntCast(ncols, &n);
    IBTK_CHKERRQ(ierr);
    const PetscBLASInt minmn = std::min(m, n);
    const PetscBLASInt lda = std::max<PetscBLASInt>(1, m);
    const PetscBLASInt ldu = std::max<PetscBLASInt>(1, m);
    const PetscBLASInt ldvt = std::max<PetscBLASInt>(1, n);

    std::vector<PetscScalar> A(static_cast<std::size_t>(lda) * static_cast<std::size_t>(n), 0.0);
    for (PetscInt i = 0; i < nrows; ++i)
    {
        PetscInt ncols_row = 0;
        const PetscInt* cols = nullptr;
        const PetscScalar* vals = nullptr;
        ierr = MatGetRow(mat, i, &ncols_row, &cols, &vals);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < ncols_row; ++k)
        {
            const PetscInt j = cols[k];
            if (0 <= j && j < ncols)
            {
                A[static_cast<std::size_t>(i) + static_cast<std::size_t>(j) * static_cast<std::size_t>(lda)] = vals[k];
            }
        }
        ierr = MatRestoreRow(mat, i, &ncols_row, &cols, &vals);
        IBTK_CHKERRQ(ierr);
    }

    PetscInt rhs_size = 0;
    ierr = VecGetSize(rhs, &rhs_size);
    IBTK_CHKERRQ(ierr);
    if (rhs_size != nrows)
    {
        TBOX_ERROR("solve_pseudoinverse(): rhs size mismatch.\n");
    }

    const PetscScalar* rhs_arr = nullptr;
    ierr = VecGetArrayRead(rhs, &rhs_arr);
    IBTK_CHKERRQ(ierr);
    std::vector<PetscScalar> b(static_cast<std::size_t>(m), 0.0);
    for (PetscInt i = 0; i < nrows; ++i)
    {
        b[static_cast<std::size_t>(i)] = rhs_arr[i];
    }
    ierr = VecRestoreArrayRead(rhs, &rhs_arr);
    IBTK_CHKERRQ(ierr);

    std::vector<PetscReal> sigma(static_cast<std::size_t>(minmn), 0.0);
    std::vector<PetscScalar> U(static_cast<std::size_t>(ldu) * static_cast<std::size_t>(m), 0.0);
    std::vector<PetscScalar> VT(static_cast<std::size_t>(ldvt) * static_cast<std::size_t>(n), 0.0);
#if defined(PETSC_USE_COMPLEX)
    std::vector<PetscReal> rwork(static_cast<std::size_t>(5 * std::max<PetscBLASInt>(1, minmn)), 0.0);
#endif

    PetscBLASInt info = 0;
    PetscBLASInt lwork = -1;
    PetscScalar work_query = 0.0;
    const char jobu = 'A';
    const char jobvt = 'A';
    LAPACKgesvd_(&jobu,
                 &jobvt,
                 &m,
                 &n,
                 A.data(),
                 &lda,
                 sigma.data(),
                 U.data(),
                 &ldu,
                 VT.data(),
                 &ldvt,
                 &work_query,
                 &lwork,
#if defined(PETSC_USE_COMPLEX)
                 rwork.data(),
#endif
                 &info);
    if (info != 0)
    {
        TBOX_ERROR("solve_pseudoinverse(): LAPACKgesvd workspace query failed with info = " << info << ".\n");
    }

    lwork = std::max<PetscBLASInt>(1, static_cast<PetscBLASInt>(PetscRealPart(work_query)));
    std::vector<PetscScalar> work(static_cast<std::size_t>(lwork), 0.0);
    LAPACKgesvd_(&jobu,
                 &jobvt,
                 &m,
                 &n,
                 A.data(),
                 &lda,
                 sigma.data(),
                 U.data(),
                 &ldu,
                 VT.data(),
                 &ldvt,
                 work.data(),
                 &lwork,
#if defined(PETSC_USE_COMPLEX)
                 rwork.data(),
#endif
                 &info);
    if (info != 0)
    {
        TBOX_ERROR("solve_pseudoinverse(): LAPACKgesvd failed with info = " << info << ".\n");
    }

    PetscReal sigma_max = 0.0;
    for (PetscBLASInt i = 0; i < minmn; ++i)
    {
        sigma_max = std::max(sigma_max, PetscAbsReal(sigma[static_cast<std::size_t>(i)]));
    }
    const double tol =
        std::max(static_cast<double>(std::max(nrows, ncols)) * std::numeric_limits<double>::epsilon() * sigma_max,
                 1.0e-14 * sigma_max);

    std::vector<PetscScalar> ub(static_cast<std::size_t>(m), 0.0);
    for (PetscBLASInt i = 0; i < m; ++i)
    {
        PetscScalar sum = 0.0;
        for (PetscBLASInt j = 0; j < m; ++j)
        {
            sum += PetscConj(
                       U[static_cast<std::size_t>(j) + static_cast<std::size_t>(i) * static_cast<std::size_t>(ldu)]) *
                   b[static_cast<std::size_t>(j)];
        }
        ub[static_cast<std::size_t>(i)] = sum;
    }

    std::vector<PetscScalar> z(static_cast<std::size_t>(minmn), 0.0);
    for (PetscBLASInt i = 0; i < minmn; ++i)
    {
        const double s = sigma[static_cast<std::size_t>(i)];
        if (std::abs(s) > tol)
        {
            z[static_cast<std::size_t>(i)] = ub[static_cast<std::size_t>(i)] / s;
        }
    }

    std::vector<PetscScalar> x(static_cast<std::size_t>(n), 0.0);
    for (PetscBLASInt i = 0; i < n; ++i)
    {
        PetscScalar sum = 0.0;
        for (PetscBLASInt k = 0; k < minmn; ++k)
        {
            sum += PetscConj(
                       VT[static_cast<std::size_t>(k) + static_cast<std::size_t>(i) * static_cast<std::size_t>(ldvt)]) *
                   z[static_cast<std::size_t>(k)];
        }
        x[static_cast<std::size_t>(i)] = sum;
    }

    std::vector<double> out(static_cast<std::size_t>(n), 0.0);
    for (PetscBLASInt i = 0; i < n; ++i)
    {
        out[static_cast<std::size_t>(i)] = PetscRealPart(x[static_cast<std::size_t>(i)]);
    }
    return out;
}

void
configure_double_stream(std::ostream& os)
{
    os << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
}

void
configure_long_double_stream(std::ostream& os)
{
    os << std::scientific << std::setprecision(std::numeric_limits<long double>::max_digits10);
}

std::uint64_t
double_to_bits(const double value)
{
    std::uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value), "unexpected double size");
    std::memcpy(&bits, &value, sizeof(value));
    return bits;
}

void
write_matrix_ref(const std::string& filename, const MatrixData& data)
{
    std::ofstream os(filename);
    os << data.nrows << " " << data.ncols << " " << data.entries.size() << "\n";
    configure_double_stream(os);
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
    configure_double_stream(os);
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
check_lagrangian_a_row_nnz(const MatrixData& A_data,
                           const double abs_nz_tol,
                           const double rel_nz_tol,
                           std::ostream& os,
                           int& test_failures)
{
    if (A_data.nrows != A_data.ncols)
    {
        ++test_failures;
        os << "A structural check failed: matrix is not square (" << A_data.nrows << "x" << A_data.ncols << ")\n";
        return false;
    }

    std::vector<double> row_max_abs(static_cast<std::size_t>(A_data.nrows), 0.0);
    for (const auto& e : A_data.entries)
    {
        auto& max_abs = row_max_abs[static_cast<std::size_t>(e.row)];
        max_abs = std::max(max_abs, std::abs(e.val));
    }

    std::vector<int> row_nnz(static_cast<std::size_t>(A_data.nrows), 0);
    for (const auto& e : A_data.entries)
    {
        const std::size_t row = static_cast<std::size_t>(e.row);
        const double row_tol = std::max(abs_nz_tol, rel_nz_tol * row_max_abs[row]);
        if (std::abs(e.val) > row_tol) ++row_nnz[row];
    }

    for (int row = 0; row < A_data.nrows; ++row)
    {
        if (row_nnz[static_cast<std::size_t>(row)] != 3)
        {
            ++test_failures;
            os << "A structural check failed: row " << row << " has " << row_nnz[static_cast<std::size_t>(row)]
               << " clearly nonzero entries; expected 3 (abs_nz_tol=" << abs_nz_tol << ", rel_nz_tol=" << rel_nz_tol
               << ", row_max_abs=" << row_max_abs[static_cast<std::size_t>(row)] << ")\n";
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
    configure_double_stream(os);
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
write_double_vector(const std::string& path, const std::vector<double>& vals)
{
    std::ofstream os(path);
    os << vals.size() << "\n";
    configure_double_stream(os);
    for (const double v : vals) os << v << "\n";
}

void
write_sweep_trace(const std::string& path, const std::vector<SweepTraceFrame>& trace)
{
    std::ofstream os(path);
    int n_dofs = 0;
    if (!trace.empty()) n_dofs = static_cast<int>(trace.front().y.size());
    os << trace.size() << " " << n_dofs << "\n";
    configure_double_stream(os);
    for (const auto& frame : trace)
    {
        os << frame.sweep << " " << frame.subdomain << " " << frame.stage << "\n";
        for (int i = 0; i < static_cast<int>(frame.y.size()); ++i)
        {
            os << frame.y[static_cast<std::size_t>(i)];
            if (i + 1 < static_cast<int>(frame.y.size())) os << " ";
        }
        os << "\n";
    }
}

void
write_sweep_residual_trace(const std::string& path, const std::vector<SweepResidualTraceFrame>& trace)
{
    std::ofstream os(path);
    int n_dofs = 0;
    if (!trace.empty()) n_dofs = static_cast<int>(trace.front().residual.size());
    os << trace.size() << " " << n_dofs << "\n";
    configure_double_stream(os);
    for (const auto& frame : trace)
    {
        os << frame.sweep << " " << frame.subdomain << " " << frame.stage << "\n";
        for (int i = 0; i < static_cast<int>(frame.residual.size()); ++i)
        {
            os << frame.residual[static_cast<std::size_t>(i)];
            if (i + 1 < static_cast<int>(frame.residual.size())) os << " ";
        }
        os << "\n";
    }
}

void
write_sweep_delta_trace(const std::string& path, const std::vector<SweepDeltaTraceFrame>& trace)
{
    std::ofstream os(path);
    int n_dofs = 0;
    if (!trace.empty()) n_dofs = static_cast<int>(trace.front().delta.size());
    os << trace.size() << " " << n_dofs << "\n";
    configure_double_stream(os);
    for (const auto& frame : trace)
    {
        os << frame.sweep << " " << frame.subdomain << " " << frame.stage << "\n";
        for (int i = 0; i < static_cast<int>(frame.delta.size()); ++i)
        {
            os << frame.delta[static_cast<std::size_t>(i)];
            if (i + 1 < static_cast<int>(frame.delta.size())) os << " ";
        }
        os << "\n";
    }
}

void
write_vec(const std::string& path, Vec vec)
{
    write_double_vector(path, extract_vector_data(vec));
}

std::vector<double>
expand_subvector_to_global(const std::vector<PetscInt>& is_vals, Vec sub_vec, const PetscInt n_global)
{
    PetscInt n_sub = 0;
    int ierr = VecGetSize(sub_vec, &n_sub);
    IBTK_CHKERRQ(ierr);
    if (n_sub != static_cast<PetscInt>(is_vals.size()))
    {
        TBOX_ERROR("Subvector/global-index size mismatch in expand_subvector_to_global.\n");
    }

    std::vector<double> out(static_cast<std::size_t>(n_global), 0.0);
    const PetscScalar* arr = nullptr;
    ierr = VecGetArrayRead(sub_vec, &arr);
    IBTK_CHKERRQ(ierr);
    for (PetscInt i = 0; i < n_sub; ++i)
    {
        out[static_cast<std::size_t>(is_vals[static_cast<std::size_t>(i)])] = static_cast<double>(arr[i]);
    }
    ierr = VecRestoreArrayRead(sub_vec, &arr);
    IBTK_CHKERRQ(ierr);
    return out;
}

void
write_debug_row_dump(const std::string& path, Mat L, PetscInt row, Vec rhs_full, Vec y_before, Vec residual_full)
{
    std::ofstream os(path);
    configure_long_double_stream(os);

    PetscInt nrows = 0, ncols_mat = 0;
    int ierr = MatGetSize(L, &nrows, &ncols_mat);
    IBTK_CHKERRQ(ierr);
    if (row < 0 || row >= nrows)
    {
        os << "status invalid_row\n";
        os << "row " << static_cast<long long>(row) << "\n";
        os << "matrix_nrows " << static_cast<long long>(nrows) << "\n";
        os << "matrix_ncols " << static_cast<long long>(ncols_mat) << "\n";
        return;
    }

    std::vector<double> y_vals = extract_vector_data(y_before);

    PetscScalar rhs_val = 0.0;
    ierr = VecGetValues(rhs_full, 1, &row, &rhs_val);
    IBTK_CHKERRQ(ierr);
    PetscScalar residual_val = 0.0;
    ierr = VecGetValues(residual_full, 1, &row, &residual_val);
    IBTK_CHKERRQ(ierr);

    PetscInt ncols = 0;
    const PetscInt* cols = nullptr;
    const PetscScalar* vals = nullptr;
    ierr = MatGetRow(L, row, &ncols, &cols, &vals);
    IBTK_CHKERRQ(ierr);

    double dot_double = 0.0;
    long double dot_long_double = 0.0L;
    for (PetscInt k = 0; k < ncols; ++k)
    {
        const PetscInt col = cols[k];
        const double a = static_cast<double>(vals[k]);
        const double y = y_vals[static_cast<std::size_t>(col)];
        dot_double += a * y;
        dot_long_double += static_cast<long double>(a) * static_cast<long double>(y);
    }

    const double residual_manual_double = static_cast<double>(rhs_val) - dot_double;
    const long double residual_manual_long_double = static_cast<long double>(rhs_val) - dot_long_double;

    os << "status ok\n";
    os << "row " << static_cast<long long>(row) << "\n";
    os << "rhs " << static_cast<long double>(rhs_val) << "\n";
    os << "petsc_residual " << static_cast<long double>(residual_val) << "\n";
    os << "manual_dot_double " << static_cast<long double>(dot_double) << "\n";
    os << "manual_dot_long_double " << dot_long_double << "\n";
    os << "manual_residual_double " << static_cast<long double>(residual_manual_double) << "\n";
    os << "manual_residual_long_double " << residual_manual_long_double << "\n";
    os << "manual_minus_petsc_double "
       << (static_cast<long double>(residual_manual_double) - static_cast<long double>(residual_val)) << "\n";
    os << "manual_minus_petsc_long_double " << (residual_manual_long_double - static_cast<long double>(residual_val))
       << "\n";
    os << "ncols " << static_cast<long long>(ncols) << "\n";
    os << "entries_begin\n";
    os << "index col value value_bits y y_bits product\n";
    for (PetscInt k = 0; k < ncols; ++k)
    {
        const PetscInt col = cols[k];
        const double a = static_cast<double>(vals[k]);
        const double y = y_vals[static_cast<std::size_t>(col)];
        os << static_cast<long long>(k) << " " << static_cast<long long>(col) << " " << a << " " << double_to_bits(a)
           << " " << y << " " << double_to_bits(y) << " " << (a * y) << "\n";
    }
    os << "entries_end\n";

    ierr = MatRestoreRow(L, row, &ncols, &cols, &vals);
    IBTK_CHKERRQ(ierr);
}

void
write_local_solve_debug_dump(const LocalSolveDebugDump& dump_cfg,
                             const std::vector<PetscInt>& is_vals,
                             Vec rhs_full,
                             Vec y_before,
                             Vec residual_full,
                             Vec r_sub,
                             Mat L_sub,
                             Mat L_full,
                             const std::vector<double>& delta_vals,
                             Vec y_after)
{
    std::filesystem::create_directories(dump_cfg.output_dir);
    std::vector<int> dofs;
    dofs.reserve(is_vals.size());
    for (const PetscInt idx : is_vals) dofs.push_back(static_cast<int>(idx));
    write_int_vector(dump_cfg.output_dir + "/subdomain_dofs.txt", dofs);
    write_vec(dump_cfg.output_dir + "/rhs_full.txt", rhs_full);
    write_vec(dump_cfg.output_dir + "/y_before_full.txt", y_before);
    write_vec(dump_cfg.output_dir + "/residual_full.txt", residual_full);
    write_vec(dump_cfg.output_dir + "/local_rhs.txt", r_sub);
    write_matrix_market(dump_cfg.output_dir + "/local_A.mtx", L_sub);
    write_double_vector(dump_cfg.output_dir + "/local_delta_cpp.txt", delta_vals);
    write_vec(dump_cfg.output_dir + "/y_after_full.txt", y_after);
    if (dump_cfg.target_global_row >= 0)
    {
        write_debug_row_dump(dump_cfg.output_dir + "/debug_row.txt",
                             L_full,
                             dump_cfg.target_global_row,
                             rhs_full,
                             y_before,
                             residual_full);
    }
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

template <class MapType>
void
write_axis_map(const std::string& path, const MapType& map_data)
{
    std::ofstream os(path);
    os << map_data.size() << "\n";
    for (const auto& kv : map_data) os << kv.first << " " << kv.second << "\n";
}

template <class Values>
void
write_dof_values(std::ofstream& os, const Values& values)
{
    for (const int v : values) os << " " << v;
}

void
write_dof_values(std::ofstream& os, const std::set<int>& values)
{
    for (const int v : values) os << " " << v;
}

template <class MapType>
void
write_set_map(const std::string& path, const MapType& map_data)
{
    std::ofstream os(path);
    os << map_data.size() << "\n";
    for (const auto& kv : map_data)
    {
        os << kv.first << " " << kv.second.size();
        write_dof_values(os, kv.second);
        os << "\n";
    }
}

using SubmatrixTriplet = std::tuple<int, int, double>;

std::vector<SubmatrixTriplet>
compute_submatrix_triplets(Mat mat, const std::set<int>& dofs)
{
    PetscInt nrows = 0, ncols_mat = 0;
    int ierr = MatGetSize(mat, &nrows, &ncols_mat);
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
        ierr = MatGetRow(mat, static_cast<PetscInt>(row), &ncols, &cols, &vals);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < ncols; ++k)
        {
            const int col = static_cast<int>(cols[k]);
            if (dof_lookup.find(col) == dof_lookup.end()) continue;
            entries.emplace_back(row, col, static_cast<double>(vals[k]));
        }
        ierr = MatRestoreRow(mat, static_cast<PetscInt>(row), &ncols, &cols, &vals);
        IBTK_CHKERRQ(ierr);
    }
    std::sort(entries.begin(), entries.end());
    return entries;
}

void
write_overlap_submatrix_blocks(const std::string& path, Mat mat, const std::vector<std::set<int>>& overlap_sets)
{
    std::ofstream os(path);
    os << overlap_sets.size() << "\n";
    configure_double_stream(os);
    for (std::size_t k = 0; k < overlap_sets.size(); ++k)
    {
        std::set<int> dof_set = overlap_sets[k];
        const std::vector<SubmatrixTriplet> entries = compute_submatrix_triplets(mat, dof_set);
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
    return std::vector<int>(axis_velocity_dofs.begin(), axis_velocity_dofs.end());
}

std::vector<int>
extract_field_dofs(const std::vector<std::set<int>>& field_is,
                   const std::vector<std::string>& field_names,
                   const std::string& name)
{
    const auto it = std::find(field_names.begin(), field_names.end(), name);
    if (it == field_names.end())
    {
        TBOX_ERROR("Missing field '" << name << "' in field index sets.\n");
    }
    const std::size_t idx = static_cast<std::size_t>(std::distance(field_names.begin(), it));
    return std::vector<int>(field_is[idx].begin(), field_is[idx].end());
}

void
assert_subdomain_structure(const std::vector<std::set<int>>& overlap_sets,
                           const std::set<int>& velocity_dofs,
                           const std::set<int>& pressure_dofs,
                           const int ln,
                           const std::string& policy_name)
{
    if (overlap_sets.size() <= 1)
    {
        TBOX_ERROR("Expected more than one CAV subdomain at level " << ln << " for policy " << policy_name << ".\n");
    }
    for (std::size_t k = 0; k < overlap_sets.size(); ++k)
    {
        bool has_velocity = false;
        bool has_pressure = false;
        for (const int dof : overlap_sets[k])
        {
            if (velocity_dofs.find(dof) != velocity_dofs.end()) has_velocity = true;
            if (pressure_dofs.find(dof) != pressure_dofs.end()) has_pressure = true;
            if (has_velocity && has_pressure) break;
        }
        if (!has_velocity || !has_pressure)
        {
            TBOX_ERROR("Invalid CAV subdomain structure at level " << ln << " subdomain " << k << " for policy "
                                                                   << policy_name << ": has_velocity=" << has_velocity
                                                                   << ", has_pressure=" << has_pressure << ".\n");
        }
    }
}

void
fill_nonzero_vector(Vec vec, const double scale)
{
    PetscInt row_start = 0, row_end = 0;
    int ierr = VecGetOwnershipRange(vec, &row_start, &row_end);
    IBTK_CHKERRQ(ierr);
    PetscScalar* arr = nullptr;
    ierr = VecGetArray(vec, &arr);
    IBTK_CHKERRQ(ierr);
    for (PetscInt i = row_start; i < row_end; ++i)
    {
        const double val = scale * (0.25 + std::sin(0.137 * static_cast<double>(i + 1)) +
                                    0.5 * std::cos(0.173 * static_cast<double>(i + 2)));
        arr[i - row_start] = static_cast<PetscScalar>(val);
    }
    ierr = VecRestoreArray(vec, &arr);
    IBTK_CHKERRQ(ierr);
}

std::pair<int, int>
infer_eulerian_grid_shape(const std::vector<EulerianDofMapEntry>& dof_map)
{
    int nx = -1;
    int ny = -1;
    for (const EulerianDofMapEntry& e : dof_map)
    {
        if (!e.set || e.type != 'p') continue;
        nx = std::max(nx, e.i + 1);
        ny = std::max(ny, e.j + 1);
    }
    if (nx <= 0 || ny <= 0)
    {
        for (const EulerianDofMapEntry& e : dof_map)
        {
            if (!e.set) continue;
            nx = std::max(nx, e.i + 1);
            ny = std::max(ny, e.j + 1);
        }
    }
    if (nx <= 0 || ny <= 0)
    {
        TBOX_ERROR("Unable to infer Eulerian grid shape for smooth V-cycle probe.\n");
    }
    return std::make_pair(nx, ny);
}

double
periodic_coordinate(const int idx, const int n, const double offset)
{
    const int wrapped = ((idx % n) + n) % n;
    return (static_cast<double>(wrapped) + offset) / static_cast<double>(n);
}

void
fill_smooth_eulerian_mode_vector(Vec vec,
                                 const std::vector<EulerianDofMapEntry>& dof_map,
                                 const double scale,
                                 const int mode_x,
                                 const int mode_y)
{
    const std::pair<int, int> grid_shape = infer_eulerian_grid_shape(dof_map);
    const int nx = grid_shape.first;
    const int ny = grid_shape.second;

    PetscInt row_start = 0, row_end = 0;
    int ierr = VecGetOwnershipRange(vec, &row_start, &row_end);
    IBTK_CHKERRQ(ierr);
    PetscScalar* arr = nullptr;
    ierr = VecGetArray(vec, &arr);
    IBTK_CHKERRQ(ierr);
    for (PetscInt i = row_start; i < row_end; ++i)
    {
        const EulerianDofMapEntry& e = dof_map[static_cast<std::size_t>(i)];
        double x = 0.0;
        double y = 0.0;
        if (e.type == 'u' && e.axis == 0)
        {
            x = periodic_coordinate(e.i, nx, 0.0);
            y = periodic_coordinate(e.j, ny, 0.5);
        }
        else if (e.type == 'u' && e.axis == 1)
        {
            x = periodic_coordinate(e.i, nx, 0.5);
            y = periodic_coordinate(e.j, ny, 0.0);
        }
        else
        {
            x = periodic_coordinate(e.i, nx, 0.5);
            y = periodic_coordinate(e.j, ny, 0.5);
        }

        const double theta_x = 2.0 * M_PI * static_cast<double>(mode_x) * x;
        const double theta_y = 2.0 * M_PI * static_cast<double>(mode_y) * y;
        double val = 0.0;
        if (e.type == 'u' && e.axis == 0)
        {
            val = std::sin(theta_x) * std::cos(theta_y);
        }
        else if (e.type == 'u' && e.axis == 1)
        {
            val = -std::cos(theta_x) * std::sin(theta_y);
        }
        else if (e.type == 'p')
        {
            val = 0.35 * std::cos(theta_x + 0.13) * std::cos(theta_y - 0.11);
        }
        arr[i - row_start] = static_cast<PetscScalar>(scale * val);
    }
    ierr = VecRestoreArray(vec, &arr);
    IBTK_CHKERRQ(ierr);
}

void
fill_vcycle_probe_state(Vec vec,
                        const std::string& probe_type,
                        const std::vector<EulerianDofMapEntry>& dof_map,
                        const double scale,
                        const int mode_x,
                        const int mode_y)
{
    if (probe_type == "index_sinusoid")
    {
        fill_nonzero_vector(vec, scale);
    }
    else if (probe_type == "smooth_eulerian_mode")
    {
        fill_smooth_eulerian_mode_vector(vec, dof_map, scale, mode_x, mode_y);
    }
    else
    {
        TBOX_ERROR("Unknown vcycle_probe_type '" << probe_type << "'.\n");
    }
}

Vec
create_serial_vec(const PetscInt size)
{
    Vec vec = nullptr;
    const int ierr = VecCreateSeq(PETSC_COMM_SELF, size, &vec);
    IBTK_CHKERRQ(ierr);
    return vec;
}

double
compute_vec_norm_inf(Vec vec)
{
    PetscReal norm = 0.0;
    const int ierr = VecNorm(vec, NORM_INFINITY, &norm);
    IBTK_CHKERRQ(ierr);
    return static_cast<double>(norm);
}

void
assert_vec_finite_nontrivial(Vec vec, const std::string& label)
{
    const double norm = compute_vec_norm_inf(vec);
    if (!std::isfinite(norm))
    {
        TBOX_ERROR(label << " has non-finite infinity norm.\n");
    }
    if (!(norm > 0.0))
    {
        TBOX_ERROR(label << " is trivial (zero infinity norm).\n");
    }
}

void
subtract_pressure_mean(Vec y, const std::vector<int>& pressure_dofs)
{
    if (pressure_dofs.empty()) return;
    std::vector<PetscInt> p_is_vals(pressure_dofs.begin(), pressure_dofs.end());
    IS pressure_is = nullptr;
    int ierr = ISCreateGeneral(
        PETSC_COMM_WORLD, static_cast<PetscInt>(p_is_vals.size()), p_is_vals.data(), PETSC_COPY_VALUES, &pressure_is);
    IBTK_CHKERRQ(ierr);
    Vec p_sub = nullptr;
    ierr = VecGetSubVector(y, pressure_is, &p_sub);
    IBTK_CHKERRQ(ierr);
    PetscScalar p_sum = 0.0;
    ierr = VecSum(p_sub, &p_sum);
    IBTK_CHKERRQ(ierr);
    PetscInt n_p = 0;
    ierr = VecGetSize(p_sub, &n_p);
    IBTK_CHKERRQ(ierr);
    if (n_p > 0)
    {
        const PetscScalar p_mean = p_sum / static_cast<PetscScalar>(n_p);
        ierr = VecShift(p_sub, -p_mean);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecRestoreSubVector(y, pressure_is, &p_sub);
    IBTK_CHKERRQ(ierr);
    ierr = ISDestroy(&pressure_is);
    IBTK_CHKERRQ(ierr);
}

void
apply_coupling_aware_sweep(Vec y,
                           Vec b,
                           Mat L,
                           const std::vector<std::set<int>>& overlap_sets,
                           const std::vector<int>& pressure_dofs,
                           const int n_sweeps,
                           const double alpha,
                           std::vector<SweepTraceFrame>* trace_frames = nullptr,
                           std::vector<SweepResidualTraceFrame>* residual_trace_frames = nullptr,
                           std::vector<SweepDeltaTraceFrame>* delta_trace_frames = nullptr,
                           LocalSolveDebugDump* debug_dump = nullptr)
{
    if (n_sweeps <= 0) return;
    Vec residual = nullptr;
    int ierr = VecDuplicate(y, &residual);
    IBTK_CHKERRQ(ierr);
    PetscInt n_global = 0;
    ierr = VecGetSize(y, &n_global);
    IBTK_CHKERRQ(ierr);

    for (int sweep = 0; sweep < n_sweeps; ++sweep)
    {
        for (int k = 0; k < static_cast<int>(overlap_sets.size()); ++k)
        {
            const std::set<int>& dof_set = overlap_sets[static_cast<std::size_t>(k)];
            std::vector<PetscInt> is_vals(dof_set.begin(), dof_set.end());
            if (is_vals.empty()) continue;

            ierr = MatMult(L, y, residual);
            IBTK_CHKERRQ(ierr);
            ierr = VecAYPX(residual, -1.0, b);
            IBTK_CHKERRQ(ierr);
            if (residual_trace_frames)
            {
                SweepResidualTraceFrame frame;
                frame.sweep = sweep;
                frame.subdomain = k;
                frame.stage = 0;
                frame.residual = extract_vector_data(residual);
                residual_trace_frames->push_back(std::move(frame));
            }

            IS sub_is = nullptr;
            ierr = ISCreateGeneral(
                PETSC_COMM_WORLD, static_cast<PetscInt>(is_vals.size()), is_vals.data(), PETSC_COPY_VALUES, &sub_is);
            IBTK_CHKERRQ(ierr);

            Vec r_sub_view = nullptr;
            ierr = VecGetSubVector(residual, sub_is, &r_sub_view);
            IBTK_CHKERRQ(ierr);
            Vec r_sub = nullptr;
            ierr = VecDuplicate(r_sub_view, &r_sub);
            IBTK_CHKERRQ(ierr);
            ierr = VecCopy(r_sub_view, r_sub);
            IBTK_CHKERRQ(ierr);
            ierr = VecRestoreSubVector(residual, sub_is, &r_sub_view);
            IBTK_CHKERRQ(ierr);

            Mat L_sub = nullptr;
            ierr = MatCreateSubMatrix(L, sub_is, sub_is, MAT_INITIAL_MATRIX, &L_sub);
            IBTK_CHKERRQ(ierr);

            Vec delta_sub = nullptr;
            ierr = VecDuplicate(r_sub, &delta_sub);
            IBTK_CHKERRQ(ierr);

            Vec y_before = nullptr;
            if (debug_dump && debug_dump->enabled && !debug_dump->dumped && debug_dump->target_sweep == sweep &&
                debug_dump->target_subdomain == k)
            {
                ierr = VecDuplicate(y, &y_before);
                IBTK_CHKERRQ(ierr);
                ierr = VecCopy(y, y_before);
                IBTK_CHKERRQ(ierr);
            }

            const std::vector<double> delta_vals = solve_pseudoinverse(L_sub, r_sub);
            if (static_cast<PetscInt>(delta_vals.size()) != static_cast<PetscInt>(is_vals.size()))
            {
                TBOX_ERROR("Local CAV subdomain solve produced unexpected size.\n");
            }
            PetscScalar* delta_arr = nullptr;
            ierr = VecGetArray(delta_sub, &delta_arr);
            IBTK_CHKERRQ(ierr);
            for (std::size_t d = 0; d < delta_vals.size(); ++d)
            {
                delta_arr[d] = static_cast<PetscScalar>(delta_vals[d]);
            }
            ierr = VecRestoreArray(delta_sub, &delta_arr);
            IBTK_CHKERRQ(ierr);
            if (alpha != 1.0)
            {
                ierr = VecScale(delta_sub, alpha);
                IBTK_CHKERRQ(ierr);
            }
            if (delta_trace_frames)
            {
                SweepDeltaTraceFrame frame;
                frame.sweep = sweep;
                frame.subdomain = k;
                frame.stage = 0;
                frame.delta = expand_subvector_to_global(is_vals, delta_sub, n_global);
                delta_trace_frames->push_back(std::move(frame));
            }

            Vec y_sub = nullptr;
            ierr = VecGetSubVector(y, sub_is, &y_sub);
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPY(y_sub, 1.0, delta_sub);
            IBTK_CHKERRQ(ierr);
            ierr = VecRestoreSubVector(y, sub_is, &y_sub);
            IBTK_CHKERRQ(ierr);

            if (y_before)
            {
                write_local_solve_debug_dump(
                    *debug_dump, is_vals, b, y_before, residual, r_sub, L_sub, L, delta_vals, y);
                debug_dump->dumped = true;
                ierr = VecDestroy(&y_before);
                IBTK_CHKERRQ(ierr);
            }

            if (trace_frames)
            {
                SweepTraceFrame frame;
                frame.sweep = sweep;
                frame.subdomain = k;
                frame.stage = 0;
                frame.y = extract_vector_data(y);
                trace_frames->push_back(std::move(frame));
            }

            ierr = VecDestroy(&delta_sub);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&L_sub);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&r_sub);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&sub_is);
            IBTK_CHKERRQ(ierr);
        }
        subtract_pressure_mean(y, pressure_dofs);
        if (trace_frames)
        {
            SweepTraceFrame frame;
            frame.sweep = sweep;
            frame.subdomain = -1;
            frame.stage = 1;
            frame.y = extract_vector_data(y);
            trace_frames->push_back(std::move(frame));
        }
    }

    ierr = VecDestroy(&residual);
    IBTK_CHKERRQ(ierr);
}

void
apply_two_level_vcycle(Vec y_out,
                       Vec rhs_fine,
                       Mat L_fine,
                       Mat L_coarse,
                       Mat R_fine_to_coarse,
                       Mat P_coarse_to_fine,
                       const std::vector<std::set<int>>& overlap_fine,
                       const std::vector<int>& pressure_dofs_fine,
                       const int num_pre_sweeps,
                       const int num_post_sweeps,
                       const double alpha,
                       VCycleTrace* trace = nullptr)
{
    int ierr;
    Vec y = nullptr;
    ierr = VecDuplicate(rhs_fine, &y);
    IBTK_CHKERRQ(ierr);
    ierr = VecZeroEntries(y);
    IBTK_CHKERRQ(ierr);

    apply_coupling_aware_sweep(y, rhs_fine, L_fine, overlap_fine, pressure_dofs_fine, num_pre_sweeps, alpha);
    if (trace) trace->y_pre = extract_vector_data(y);

    Vec residual_fine = nullptr;
    ierr = VecDuplicate(rhs_fine, &residual_fine);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(L_fine, y, residual_fine);
    IBTK_CHKERRQ(ierr);
    ierr = VecAYPX(residual_fine, -1.0, rhs_fine);
    IBTK_CHKERRQ(ierr);
    if (trace) trace->residual_fine = extract_vector_data(residual_fine);

    Vec coarse_sol = nullptr;
    Vec coarse_rhs = nullptr;
    ierr = MatCreateVecs(L_coarse, &coarse_sol, &coarse_rhs);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(R_fine_to_coarse, residual_fine, coarse_rhs);
    IBTK_CHKERRQ(ierr);
    if (trace) trace->coarse_rhs = extract_vector_data(coarse_rhs);

    const std::vector<double> coarse_vals = solve_pseudoinverse(L_coarse, coarse_rhs);
    PetscInt coarse_size = 0;
    ierr = VecGetSize(coarse_sol, &coarse_size);
    IBTK_CHKERRQ(ierr);
    if (static_cast<PetscInt>(coarse_vals.size()) != coarse_size)
    {
        TBOX_ERROR("Coarse solve produced unexpected size in apply_two_level_vcycle.\n");
    }
    PetscScalar* coarse_arr = nullptr;
    ierr = VecGetArray(coarse_sol, &coarse_arr);
    IBTK_CHKERRQ(ierr);
    for (PetscInt i = 0; i < coarse_size; ++i)
    {
        coarse_arr[i] = static_cast<PetscScalar>(coarse_vals[static_cast<std::size_t>(i)]);
    }
    ierr = VecRestoreArray(coarse_sol, &coarse_arr);
    IBTK_CHKERRQ(ierr);
    if (trace) trace->coarse_sol = extract_vector_data(coarse_sol);

    Vec fine_correction = nullptr;
    ierr = VecDuplicate(rhs_fine, &fine_correction);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(P_coarse_to_fine, coarse_sol, fine_correction);
    IBTK_CHKERRQ(ierr);
    if (trace) trace->fine_correction = extract_vector_data(fine_correction);
    ierr = VecAXPY(y, 1.0, fine_correction);
    IBTK_CHKERRQ(ierr);
    if (trace) trace->y_after_correction = extract_vector_data(y);

    apply_coupling_aware_sweep(y, rhs_fine, L_fine, overlap_fine, pressure_dofs_fine, num_post_sweeps, alpha);
    if (trace) trace->y_post = extract_vector_data(y);

    ierr = VecCopy(y, y_out);
    IBTK_CHKERRQ(ierr);

    ierr = VecDestroy(&fine_correction);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&coarse_rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&coarse_sol);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&residual_fine);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&y);
    IBTK_CHKERRQ(ierr);
}

void
apply_multilevel_vcycle(Vec y_out,
                        Vec rhs,
                        const std::vector<Mat>& level_l,
                        const std::vector<Mat>& level_restrict,
                        const std::vector<Mat>& level_prolong,
                        const std::vector<std::vector<std::set<int>>>& level_overlap,
                        const std::vector<std::vector<int>>& level_pressure_dofs,
                        const int level_num,
                        const int num_pre_sweeps,
                        const int num_post_sweeps,
                        const double alpha)
{
    int ierr;
    if (level_num < 0 || level_num >= static_cast<int>(level_l.size()))
    {
        TBOX_ERROR("apply_multilevel_vcycle(): invalid level index " << level_num << "\n");
    }

    if (level_num == 0)
    {
        const std::vector<double> coarse_vals = solve_pseudoinverse(level_l[0], rhs);
        PetscScalar* coarse_arr = nullptr;
        ierr = VecGetArray(y_out, &coarse_arr);
        IBTK_CHKERRQ(ierr);
        for (std::size_t i = 0; i < coarse_vals.size(); ++i)
        {
            coarse_arr[i] = static_cast<PetscScalar>(coarse_vals[i]);
        }
        ierr = VecRestoreArray(y_out, &coarse_arr);
        IBTK_CHKERRQ(ierr);
        subtract_pressure_mean(y_out, level_pressure_dofs[0]);
        return;
    }

    Vec y = nullptr;
    ierr = VecDuplicate(rhs, &y);
    IBTK_CHKERRQ(ierr);
    ierr = VecZeroEntries(y);
    IBTK_CHKERRQ(ierr);

    apply_coupling_aware_sweep(y,
                               rhs,
                               level_l[static_cast<std::size_t>(level_num)],
                               level_overlap[static_cast<std::size_t>(level_num)],
                               level_pressure_dofs[static_cast<std::size_t>(level_num)],
                               num_pre_sweeps,
                               alpha);

    Vec residual = nullptr;
    ierr = VecDuplicate(rhs, &residual);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(level_l[static_cast<std::size_t>(level_num)], y, residual);
    IBTK_CHKERRQ(ierr);
    ierr = VecAYPX(residual, -1.0, rhs);
    IBTK_CHKERRQ(ierr);

    Vec coarse_rhs = nullptr;
    ierr = MatCreateVecs(level_l[static_cast<std::size_t>(level_num - 1)], &coarse_rhs, nullptr);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(level_restrict[static_cast<std::size_t>(level_num)], residual, coarse_rhs);
    IBTK_CHKERRQ(ierr);

    Vec coarse_sol = nullptr;
    ierr = VecDuplicate(coarse_rhs, &coarse_sol);
    IBTK_CHKERRQ(ierr);
    apply_multilevel_vcycle(coarse_sol,
                            coarse_rhs,
                            level_l,
                            level_restrict,
                            level_prolong,
                            level_overlap,
                            level_pressure_dofs,
                            level_num - 1,
                            num_pre_sweeps,
                            num_post_sweeps,
                            alpha);

    Vec fine_correction = nullptr;
    ierr = VecDuplicate(y, &fine_correction);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(level_prolong[static_cast<std::size_t>(level_num)], coarse_sol, fine_correction);
    IBTK_CHKERRQ(ierr);
    ierr = VecAXPY(y, 1.0, fine_correction);
    IBTK_CHKERRQ(ierr);

    apply_coupling_aware_sweep(y,
                               rhs,
                               level_l[static_cast<std::size_t>(level_num)],
                               level_overlap[static_cast<std::size_t>(level_num)],
                               level_pressure_dofs[static_cast<std::size_t>(level_num)],
                               num_post_sweeps,
                               alpha);

    ierr = VecCopy(y, y_out);
    IBTK_CHKERRQ(ierr);

    ierr = VecDestroy(&fine_correction);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&coarse_sol);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&coarse_rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&residual);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&y);
    IBTK_CHKERRQ(ierr);
}

PetscErrorCode
apply_matrix_vcycle_shell_pc(PC pc, Vec x, Vec y)
{
    void* ctx_void = nullptr;
    PetscErrorCode ierr = PCShellGetContext(pc, &ctx_void);
    CHKERRQ(ierr);
    auto* ctx = static_cast<MatrixVcycleShellContext*>(ctx_void);
    if (!ctx || !ctx->level_l || !ctx->level_overlap || !ctx->level_pressure_dofs)
    {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL, "matrix V-cycle shell context is not initialized");
    }

    apply_multilevel_vcycle(y,
                            x,
                            *ctx->level_l,
                            *ctx->level_restrict,
                            *ctx->level_prolong,
                            *ctx->level_overlap,
                            *ctx->level_pressure_dofs,
                            static_cast<int>(ctx->level_l->size()) - 1,
                            ctx->num_pre_sweeps,
                            ctx->num_post_sweeps,
                            ctx->alpha);
    return 0;
}

MatrixKrylovSummary
run_matrix_gmres_smoke(Vec rhs,
                       const std::vector<Mat>& level_l,
                       const std::vector<Mat>& level_restrict,
                       const std::vector<Mat>& level_prolong,
                       const std::vector<std::vector<std::set<int>>>& level_overlap,
                       const std::vector<std::vector<int>>& level_pressure_dofs,
                       const int num_pre_sweeps,
                       const int num_post_sweeps,
                       const double alpha,
                       const double tol,
                       const int max_it,
                       const int restart)
{
    if (level_l.empty())
    {
        TBOX_ERROR("run_matrix_gmres_smoke(): no level operators are available.\n");
    }

    PetscErrorCode ierr;
    MatrixKrylovSummary summary;
    MatrixVcycleShellContext shell_ctx;
    shell_ctx.level_l = &level_l;
    shell_ctx.level_restrict = &level_restrict;
    shell_ctx.level_prolong = &level_prolong;
    shell_ctx.level_overlap = &level_overlap;
    shell_ctx.level_pressure_dofs = &level_pressure_dofs;
    shell_ctx.num_pre_sweeps = num_pre_sweeps;
    shell_ctx.num_post_sweeps = num_post_sweeps;
    shell_ctx.alpha = alpha;

    KSP ksp = nullptr;
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, level_l.back(), level_l.back());
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPGMRES);
    IBTK_CHKERRQ(ierr);
    ierr = KSPGMRESSetRestart(ksp, restart);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetPCSide(ksp, PC_LEFT);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetNormType(ksp, KSP_NORM_PRECONDITIONED);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, max_it);
    IBTK_CHKERRQ(ierr);

    PC pc = nullptr;
    ierr = KSPGetPC(ksp, &pc);
    IBTK_CHKERRQ(ierr);
    ierr = PCSetType(pc, PCSHELL);
    IBTK_CHKERRQ(ierr);
    ierr = PCShellSetName(pc, "matrix_cav_vcycle");
    IBTK_CHKERRQ(ierr);
    ierr = PCShellSetContext(pc, &shell_ctx);
    IBTK_CHKERRQ(ierr);
    ierr = PCShellSetApply(pc, apply_matrix_vcycle_shell_pc);
    IBTK_CHKERRQ(ierr);

    ierr = KSPSetUp(ksp);
    IBTK_CHKERRQ(ierr);

    Vec x = nullptr;
    ierr = VecDuplicate(rhs, &x);
    IBTK_CHKERRQ(ierr);
    ierr = VecZeroEntries(x);
    IBTK_CHKERRQ(ierr);

    Vec mb = nullptr;
    ierr = VecDuplicate(rhs, &mb);
    IBTK_CHKERRQ(ierr);
    ierr = PCApply(pc, rhs, mb);
    IBTK_CHKERRQ(ierr);
    double preconditioned_rhs_norm = 0.0;
    ierr = VecNorm(mb, NORM_2, &preconditioned_rhs_norm);
    IBTK_CHKERRQ(ierr);
    if (!(preconditioned_rhs_norm > 0.0))
    {
        preconditioned_rhs_norm = 1.0;
    }

    ierr = KSPSolve(ksp, rhs, x);
    IBTK_CHKERRQ(ierr);

    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(ksp, &reason);
    IBTK_CHKERRQ(ierr);
    summary.converged = static_cast<int>(reason) > 0;
    summary.converged_reason = static_cast<int>(reason);

    ierr = KSPGetIterationNumber(ksp, &summary.iterations);
    IBTK_CHKERRQ(ierr);
    ierr = KSPGetResidualNorm(ksp, &summary.residual_norm);
    IBTK_CHKERRQ(ierr);
    summary.relative_residual_norm = summary.residual_norm / preconditioned_rhs_norm;

    ierr = VecDestroy(&mb);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    IBTK_CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp);
    IBTK_CHKERRQ(ierr);

    return summary;
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
    const bool export_bridge_data = test_db->getBoolWithDefault("export_bridge_data", false);
    const std::string export_bridge_dir = test_db->getStringWithDefault("export_bridge_dir", "bridge-generated");
    const int seed_velocity_axis = test_db->getIntegerWithDefault("coupling_aware_asm_seed_axis", 0);
    const int seed_velocity_stride = test_db->getIntegerWithDefault("coupling_aware_asm_seed_stride", 1);
    const double shell_pc_relaxation_factor = test_db->getDoubleWithDefault("shell_pc_relaxation_factor", 1.0);
    const int num_pre_sweeps = test_db->getIntegerWithDefault("num_pre_sweeps", 1);
    const int num_post_sweeps = test_db->getIntegerWithDefault("num_post_sweeps", 1);
    const bool run_matrix_krylov_smoke = test_db->getBoolWithDefault("run_matrix_krylov_smoke", false);
    const std::string matrix_krylov_policy = test_db->getStringWithDefault("matrix_krylov_policy", "relaxed");
    const double matrix_krylov_tol = test_db->getDoubleWithDefault("matrix_krylov_tol", 1.0e-8);
    const int matrix_krylov_max_it = test_db->getIntegerWithDefault("matrix_krylov_max_it", 150);
    const int matrix_krylov_restart = test_db->getIntegerWithDefault("matrix_krylov_restart", 150);
    const std::string matrix_krylov_summary_file =
        test_db->getStringWithDefault("matrix_krylov_summary_file", "matrix_gmres_summary.txt");
    const std::string matrix_krylov_export_dir = test_db->getStringWithDefault("matrix_krylov_export_dir", "");
    const bool run_live_fac_probe = test_db->getBoolWithDefault("run_live_fac_probe", false);
    const double live_fac_probe_scale = test_db->getDoubleWithDefault("live_fac_probe_scale", 1.0);
    const std::string live_fac_probe_export_dir = test_db->getStringWithDefault("live_fac_probe_export_dir", "");
    const std::string vcycle_probe_type = test_db->getStringWithDefault("vcycle_probe_type", "index_sinusoid");
    const int vcycle_probe_mode_x = test_db->getIntegerWithDefault("vcycle_probe_mode_x", 1);
    const int vcycle_probe_mode_y = test_db->getIntegerWithDefault("vcycle_probe_mode_y", 1);
    const bool debug_local_solve_dump = test_db->getBoolWithDefault("debug_local_solve_dump", false);
    const int debug_local_solve_dump_level = test_db->getIntegerWithDefault("debug_local_solve_dump_level", -1);
    const int debug_local_solve_dump_sweep = test_db->getIntegerWithDefault("debug_local_solve_dump_sweep", -1);
    const int debug_local_solve_dump_subdomain = test_db->getIntegerWithDefault("debug_local_solve_dump_subdomain", -1);
    const int debug_local_solve_dump_global_row =
        test_db->getIntegerWithDefault("debug_local_solve_dump_global_row", -1);
    const std::string debug_local_solve_dump_policy =
        test_db->getStringWithDefault("debug_local_solve_dump_policy", "relaxed");
    const std::string debug_local_solve_dump_dir =
        test_db->getStringWithDefault("debug_local_solve_dump_dir", "local-solve-debug");
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
    if (export_bridge_data)
    {
        std::filesystem::create_directories(export_bridge_dir);
    }

    s_domain_length = input_db->getDoubleWithDefault("L", L_DOMAIN);
    s_finest_grid_cells = compute_finest_grid_cells(input_db);
    if (!(s_domain_length > 0.0) || s_finest_grid_cells <= 0)
    {
        TBOX_ERROR(
            "implicit_stokes_ib_operator_chain_parity_01 requires positive domain length and finest grid size.\n");
    }
    s_dx = s_domain_length / static_cast<double>(s_finest_grid_cells);
    s_dt = 0.5 * s_dx;

    s_spring_stiffness = input_db->getDoubleWithDefault("SPRING_STIFFNESS", K_SPRING_DEFAULT);
    if (!(s_spring_stiffness > 0.0))
    {
        TBOX_ERROR("implicit_stokes_ib_operator_chain_parity_01 requires SPRING_STIFFNESS > 0.\n");
    }
    const std::string spring_normalization = input_db->getStringWithDefault("SPRING_NORMALIZATION", "JULIA_DENSITY");
    if (spring_normalization == "JULIA_DENSITY")
    {
        s_spring_normalization = SpringNormalization::JuliaDensity;
    }
    else if (spring_normalization == "IBAMR_NODAL")
    {
        s_spring_normalization = SpringNormalization::IBAMRNodal;
    }
    else
    {
        TBOX_ERROR("SPRING_NORMALIZATION must be JULIA_DENSITY or IBAMR_NODAL.\n");
    }
    if ((vcycle_probe_type == "smooth_eulerian_mode") &&
        !(vcycle_probe_mode_x >= 0 && vcycle_probe_mode_y >= 0 && (vcycle_probe_mode_x > 0 || vcycle_probe_mode_y > 0)))
    {
        TBOX_ERROR("smooth_eulerian_mode requires nonnegative probe modes with at least one positive mode.\n");
    }
    if (matrix_krylov_policy != "relaxed" && matrix_krylov_policy != "strict")
    {
        TBOX_ERROR("matrix_krylov_policy must be either 'relaxed' or 'strict'.\n");
    }
    if (run_matrix_krylov_smoke &&
        (!(matrix_krylov_tol > 0.0) || matrix_krylov_max_it <= 0 || matrix_krylov_restart <= 0))
    {
        TBOX_ERROR("matrix GMRES smoke mode requires positive tolerance, max_it, and restart.\n");
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
    const double matrix_spread_scale = compute_matrix_spread_scale();
    ierr = MatScale(SAJ, -s_dt * matrix_spread_scale);
    IBTK_CHKERRQ(ierr);

    SAMRAI::solv::PoissonSpecifications U_problem_coefs("U_problem_coefs");
    U_problem_coefs.setCConstant(RHO / s_dt);
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
    const MatrixData L_data = extract_matrix_data(A00_full);
    std::vector<MatrixData> SAJ_level_data(static_cast<std::size_t>(finest_ln + 1));
    std::vector<MatrixData> A00_level_data(static_cast<std::size_t>(finest_ln + 1));
    std::vector<MatrixData> L_level_data(static_cast<std::size_t>(finest_ln + 1));

    const bool build_all_levels =
        check_all_levels || export_bridge_data || run_matrix_krylov_smoke || run_live_fac_probe;
    if (build_all_levels)
    {
        const std::size_t n_levels = static_cast<std::size_t>(finest_ln + 1);
        std::vector<std::vector<int>> level_num_dofs_per_proc(n_levels);
        std::vector<Mat> level_saj(n_levels, nullptr);
        std::vector<Mat> level_stokes(n_levels, nullptr);
        std::vector<Mat> level_l(n_levels, nullptr);
        std::vector<Mat> level_a00(n_levels, nullptr);
        std::vector<Mat> level_restrict(n_levels, nullptr); // maps level ln -> ln-1 for ln >= 1
        std::vector<Mat> level_prolong(n_levels, nullptr);  // maps level ln-1 -> ln for ln >= 1
        std::vector<StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData> level_map_data(n_levels);
        std::vector<std::vector<std::set<int>>> level_overlap_relaxed(n_levels), level_nonoverlap_relaxed(n_levels);
        std::vector<std::vector<std::set<int>>> level_overlap_strict(n_levels), level_nonoverlap_strict(n_levels);
        std::vector<std::vector<int>> level_seed_velocity_dofs_axis0(n_levels);
        std::vector<std::vector<int>> level_velocity_dofs(n_levels), level_pressure_dofs(n_levels);
        std::vector<std::vector<EulerianDofMapEntry>> level_eulerian_dof_maps(n_levels);

        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> ln_level = patch_hierarchy->getPatchLevel(ln);
            IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
                level_num_dofs_per_proc[static_cast<std::size_t>(ln)], u_dof_index_idx, p_dof_index_idx, ln_level);
        }

        ierr = MatDuplicate(SAJ, MAT_COPY_VALUES, &level_saj[static_cast<std::size_t>(finest_ln)]);
        IBTK_CHKERRQ(ierr);
        for (int ln = finest_ln - 1; ln >= 0; --ln)
        {
            Pointer<PatchLevel<NDIM>> fine_level = patch_hierarchy->getPatchLevel(ln + 1);
            Pointer<PatchLevel<NDIM>> coarse_level = patch_hierarchy->getPatchLevel(ln);

            AO coarse_ao = nullptr;
            int u_coarse_ao_offset = 0, p_coarse_ao_offset = 0;
            IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelAO(
                coarse_ao,
                level_num_dofs_per_proc[static_cast<std::size_t>(ln)],
                u_dof_index_idx,
                p_dof_index_idx,
                coarse_level,
                u_coarse_ao_offset,
                p_coarse_ao_offset);

            IBAMR::StaggeredStokesPETScMatUtilities::constructProlongationOp(
                level_prolong[static_cast<std::size_t>(ln + 1)],
                "RT0",
                "CONSERVATIVE",
                u_dof_index_idx,
                p_dof_index_idx,
                level_num_dofs_per_proc[static_cast<std::size_t>(ln + 1)],
                level_num_dofs_per_proc[static_cast<std::size_t>(ln)],
                fine_level,
                coarse_level,
                coarse_ao,
                u_coarse_ao_offset,
                p_coarse_ao_offset);

            Vec restriction_scale = nullptr;
            IBTK::PETScMatUtilities::constructRestrictionScalingOp(level_prolong[static_cast<std::size_t>(ln + 1)],
                                                                   restriction_scale);
            ierr = MatTranspose(
                level_prolong[static_cast<std::size_t>(ln + 1)], MAT_INITIAL_MATRIX, &level_restrict[ln + 1]);
            IBTK_CHKERRQ(ierr);
            ierr = MatDiagonalScale(level_restrict[static_cast<std::size_t>(ln + 1)], restriction_scale, nullptr);
            IBTK_CHKERRQ(ierr);

            ierr = MatPtAP(level_saj[static_cast<std::size_t>(ln + 1)],
                           level_prolong[static_cast<std::size_t>(ln + 1)],
                           MAT_INITIAL_MATRIX,
                           1.0,
                           &level_saj[static_cast<std::size_t>(ln)]);
            IBTK_CHKERRQ(ierr);
            ierr = MatDiagonalScale(level_saj[static_cast<std::size_t>(ln)], restriction_scale, nullptr);
            IBTK_CHKERRQ(ierr);

            ierr = VecDestroy(&restriction_scale);
            IBTK_CHKERRQ(ierr);
            ierr = AODestroy(&coarse_ao);
            IBTK_CHKERRQ(ierr);
        }

        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            const std::size_t ln_idx = static_cast<std::size_t>(ln);
            Pointer<PatchLevel<NDIM>> ln_level = patch_hierarchy->getPatchLevel(ln);
            const std::vector<int>& ln_num_dofs_per_proc = level_num_dofs_per_proc[ln_idx];

            IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(level_stokes[ln_idx],
                                                                                    U_problem_coefs,
                                                                                    u_bc_coefs,
                                                                                    data_time,
                                                                                    ln_num_dofs_per_proc,
                                                                                    u_dof_index_idx,
                                                                                    p_dof_index_idx,
                                                                                    ln_level);
            ierr = MatDuplicate(level_stokes[ln_idx], MAT_COPY_VALUES, &level_l[ln_idx]);
            IBTK_CHKERRQ(ierr);
            ierr = MatAXPY(level_l[ln_idx], 1.0, level_saj[ln_idx], DIFFERENT_NONZERO_PATTERN);
            IBTK_CHKERRQ(ierr);
            IBAMR::StaggeredStokesPETScMatUtilities::constructA00VelocitySubmatrix(
                level_a00[ln_idx], level_l[ln_idx], ln_num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, ln_level);

            Mat ln_saj_u = nullptr;
            IBAMR::StaggeredStokesPETScMatUtilities::constructA00VelocitySubmatrix(
                ln_saj_u, level_saj[ln_idx], ln_num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, ln_level);
            SAJ_level_data[ln_idx] = extract_matrix_data(ln_saj_u);
            A00_level_data[ln_idx] = extract_matrix_data(level_a00[ln_idx]);
            L_level_data[ln_idx] = extract_matrix_data(level_l[ln_idx]);
            ierr = MatDestroy(&ln_saj_u);
            IBTK_CHKERRQ(ierr);

            StaggeredStokesPETScMatUtilities::buildPatchLevelCellClosureMaps(
                level_map_data[ln_idx], u_dof_index_idx, p_dof_index_idx, ln_level);

            StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains(
                level_overlap_relaxed[ln_idx],
                level_nonoverlap_relaxed[ln_idx],
                ln_num_dofs_per_proc,
                u_dof_index_idx,
                ln_level,
                Pointer<CoarseFineBoundary<NDIM>>(nullptr),
                level_a00[ln_idx],
                level_map_data[ln_idx],
                seed_velocity_axis,
                seed_velocity_stride,
                CouplingAwareASMSeedTraversalOrder::I_J,
                CouplingAwareASMClosurePolicy::RELAXED);

            StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains(
                level_overlap_strict[ln_idx],
                level_nonoverlap_strict[ln_idx],
                ln_num_dofs_per_proc,
                u_dof_index_idx,
                ln_level,
                Pointer<CoarseFineBoundary<NDIM>>(nullptr),
                level_a00[ln_idx],
                level_map_data[ln_idx],
                seed_velocity_axis,
                seed_velocity_stride,
                CouplingAwareASMSeedTraversalOrder::I_J,
                CouplingAwareASMClosurePolicy::STRICT);

            std::vector<std::set<int>> field_is;
            std::vector<std::string> field_names;
            StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
                field_is, field_names, ln_num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, ln_level);
            level_velocity_dofs[ln_idx] = extract_field_dofs(field_is, field_names, "velocity");
            level_pressure_dofs[ln_idx] = extract_field_dofs(field_is, field_names, "pressure");
            const std::set<int> velocity_dof_set(level_velocity_dofs[ln_idx].begin(),
                                                 level_velocity_dofs[ln_idx].end());
            const std::set<int> pressure_dof_set(level_pressure_dofs[ln_idx].begin(),
                                                 level_pressure_dofs[ln_idx].end());

            assert_subdomain_structure(
                level_overlap_relaxed[ln_idx], velocity_dof_set, pressure_dof_set, ln, "RELAXED");
            assert_subdomain_structure(level_overlap_strict[ln_idx], velocity_dof_set, pressure_dof_set, ln, "STRICT");

            level_seed_velocity_dofs_axis0[ln_idx] = compute_ordered_seed_velocity_dofs(
                ln_level, seed_velocity_axis, u_dof_index_idx, level_map_data[ln_idx].velocity_dof_to_component_axis);
            if (level_seed_velocity_dofs_axis0[ln_idx].size() != level_overlap_relaxed[ln_idx].size() ||
                level_seed_velocity_dofs_axis0[ln_idx].size() != level_overlap_strict[ln_idx].size())
            {
                TBOX_ERROR("Level " << ln << " seed/subdomain mismatch in coupling-aware export path.\n");
            }

            const int ln_total_dofs = std::accumulate(ln_num_dofs_per_proc.begin(), ln_num_dofs_per_proc.end(), 0);
            level_eulerian_dof_maps[ln_idx] =
                build_eulerian_dof_map(ln_total_dofs, u_dof_index_idx, p_dof_index_idx, ln_level);

            if (write_reference)
            {
                write_eulerian_dof_map_ref(ref_path_level(ref_dir, "eulerian_dof_map", ln),
                                           level_eulerian_dof_maps[ln_idx]);
            }
        }

        if (run_live_fac_probe)
        {
            if (!live_fac_probe_export_dir.empty())
            {
                std::filesystem::create_directories(live_fac_probe_export_dir);
                const int num_lag_nodes = compute_num_lag_nodes();
                const double lag_ds = compute_lag_ds();
                std::ofstream os(std::filesystem::path(live_fac_probe_export_dir) / "problem_setup.txt");
                configure_double_stream(os);
                os << "domain_length " << s_domain_length << "\n";
                os << "finest_grid_cells " << s_finest_grid_cells << "\n";
                os << "finest_level_number " << finest_ln << "\n";
                os << "depth " << (finest_ln + 1) << "\n";
                os << "dx " << s_dx << "\n";
                os << "dt " << s_dt << "\n";
                os << "num_lag_nodes " << num_lag_nodes << "\n";
                os << "lag_ds " << lag_ds << "\n";
                os << "spring_stiffness_input " << s_spring_stiffness << "\n";
                os << "spring_pair_k " << compute_pair_spring_stiffness() << "\n";
                os << "spring_normalization "
                   << (s_spring_normalization == SpringNormalization::JuliaDensity ? "JULIA_DENSITY" : "IBAMR_NODAL")
                   << "\n";
                os << "rho " << RHO << "\n";
                os << "mu " << MU << "\n";
            }

            VariableDatabase<NDIM>* live_var_db = VariableDatabase<NDIM>::getDatabase();
            Pointer<VariableContext> scratch_ctx = navier_stokes_integrator->getScratchContext();
            Pointer<Variable<NDIM>> u_var = navier_stokes_integrator->getVelocityVariable();
            Pointer<Variable<NDIM>> p_var = navier_stokes_integrator->getPressureVariable();
            const int u_scratch_idx = live_var_db->mapVariableAndContextToIndex(u_var, scratch_ctx);
            const int p_scratch_idx = live_var_db->mapVariableAndContextToIndex(p_var, scratch_ctx);

            Pointer<VariableContext> live_rhs_ctx = live_var_db->getContext("operator_chain_live_fac_rhs_ctx");
            Pointer<SideVariable<NDIM, double>> f_var = new SideVariable<NDIM, double>("operator_chain_live_fac_f");
            Pointer<CellVariable<NDIM, double>> g_var = new CellVariable<NDIM, double>("operator_chain_live_fac_g");
            const int f_rhs_idx = live_var_db->registerVariableAndContext(f_var, live_rhs_ctx, IntVector<NDIM>(1));
            const int g_rhs_idx = live_var_db->registerVariableAndContext(g_var, live_rhs_ctx, IntVector<NDIM>(1));

            Pointer<HierarchySideDataOpsReal<NDIM, double>> hier_velocity_data_ops =
                new HierarchySideDataOpsReal<NDIM, double>(patch_hierarchy, 0, finest_ln);
            Pointer<HierarchyCellDataOpsReal<NDIM, double>> hier_pressure_data_ops =
                new HierarchyCellDataOpsReal<NDIM, double>(patch_hierarchy, 0, finest_ln);
            HierarchyMathOps hier_math_ops("operator_chain_live_fac_math_ops", patch_hierarchy);
            const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
            const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

            Pointer<SAMRAIVectorReal<NDIM, double>> live_sol_vec =
                new SAMRAIVectorReal<NDIM, double>("operator_chain_live_fac_sol", patch_hierarchy, 0, finest_ln);
            live_sol_vec->addComponent(u_var, u_scratch_idx, wgt_sc_idx, hier_velocity_data_ops);
            live_sol_vec->addComponent(p_var, p_scratch_idx, wgt_cc_idx, hier_pressure_data_ops);
            live_sol_vec->allocateVectorData(data_time);
            live_sol_vec->setToScalar(0.0);

            Pointer<SAMRAIVectorReal<NDIM, double>> live_rhs_vec =
                new SAMRAIVectorReal<NDIM, double>("operator_chain_live_fac_rhs", patch_hierarchy, 0, finest_ln);
            live_rhs_vec->addComponent(f_var, f_rhs_idx, wgt_sc_idx, hier_velocity_data_ops);
            live_rhs_vec->addComponent(g_var, g_rhs_idx, wgt_cc_idx, hier_pressure_data_ops);
            live_rhs_vec->allocateVectorData(data_time);
            live_rhs_vec->setToScalar(0.0);

            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM>> ln_level = patch_hierarchy->getPatchLevel(ln);
                if (!ln_level->checkAllocated(u_scratch_idx)) ln_level->allocatePatchData(u_scratch_idx, data_time);
                if (!ln_level->checkAllocated(p_scratch_idx)) ln_level->allocatePatchData(p_scratch_idx, data_time);
            }

            auto export_level_patch_vec =
                [&](const std::string& basename, const int u_data_idx, const int p_data_idx, const int ln)
            {
                const std::size_t ln_idx = static_cast<std::size_t>(ln);
                const int total_dofs =
                    std::accumulate(level_num_dofs_per_proc[ln_idx].begin(), level_num_dofs_per_proc[ln_idx].end(), 0);
                Vec level_vec = create_serial_vec(total_dofs);
                StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(level_vec,
                                                                      u_data_idx,
                                                                      u_dof_index_idx,
                                                                      p_data_idx,
                                                                      p_dof_index_idx,
                                                                      patch_hierarchy->getPatchLevel(ln));
                if (!live_fac_probe_export_dir.empty())
                {
                    const std::string suffix = level_suffix(ln);
                    write_vec(
                        (std::filesystem::path(live_fac_probe_export_dir) / (basename + suffix + ".txt")).string(),
                        level_vec);
                }
                PetscReal norm = 0.0;
                const int ierr = VecNorm(level_vec, NORM_2, &norm);
                IBTK_CHKERRQ(ierr);
                const std::string label = basename + "_l2_norm" + level_suffix(ln);
                pout << label << " = " << std::setprecision(17) << static_cast<double>(norm) << "\n";
                const int destroy_ierr = VecDestroy(&level_vec);
                IBTK_CHKERRQ(destroy_ierr);
            };

            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                const std::size_t ln_idx = static_cast<std::size_t>(ln);
                const int total_dofs =
                    std::accumulate(level_num_dofs_per_proc[ln_idx].begin(), level_num_dofs_per_proc[ln_idx].end(), 0);
                if (static_cast<int>(level_eulerian_dof_maps[ln_idx].size()) != total_dofs)
                {
                    TBOX_ERROR("Live FAC probe level " << ln << " DOF-map size mismatch.\n");
                }
                if (!live_fac_probe_export_dir.empty())
                {
                    write_eulerian_dof_map_ref((std::filesystem::path(live_fac_probe_export_dir) /
                                                ("eulerian_dof_map" + level_suffix(ln) + ".txt"))
                                                   .string(),
                                               level_eulerian_dof_maps[ln_idx]);
                }

                Vec level_probe = create_serial_vec(total_dofs);
                fill_vcycle_probe_state(level_probe,
                                        vcycle_probe_type,
                                        level_eulerian_dof_maps[ln_idx],
                                        live_fac_probe_scale,
                                        vcycle_probe_mode_x,
                                        vcycle_probe_mode_y);
                Pointer<RefineSchedule<NDIM>> data_synch_sched =
                    StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(
                        f_rhs_idx, g_rhs_idx, patch_hierarchy->getPatchLevel(ln));
                Pointer<RefineSchedule<NDIM>> ghost_fill_sched =
                    StaggeredStokesPETScVecUtilities::constructGhostFillSchedule(
                        f_rhs_idx, g_rhs_idx, patch_hierarchy->getPatchLevel(ln));
                StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(level_probe,
                                                                        f_rhs_idx,
                                                                        u_dof_index_idx,
                                                                        g_rhs_idx,
                                                                        p_dof_index_idx,
                                                                        patch_hierarchy->getPatchLevel(ln),
                                                                        data_synch_sched,
                                                                        ghost_fill_sched);
                const int ierr = VecDestroy(&level_probe);
                IBTK_CHKERRQ(ierr);
            }

            Pointer<Database> stokes_ib_precond_db =
                input_db->isDatabase("stokes_ib_precond_db") ? input_db->getDatabase("stokes_ib_precond_db") : nullptr;
            Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_op =
                new StaggeredStokesIBLevelRelaxationFACOperator(
                    "operator_chain_live_fac::fac_op", stokes_ib_precond_db, "stokes_ib_pc_");
            Pointer<StaggeredStokesIBJacobianFACPreconditioner> fac_pc = new StaggeredStokesIBJacobianFACPreconditioner(
                "operator_chain_live_fac::fac_pc", fac_op, stokes_ib_precond_db, "stokes_ib_pc_");
            Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
            const std::vector<RobinBcCoefStrategy<NDIM>*> live_u_bc_coefs =
                navier_stokes_integrator->getIntermediateVelocityBoundaryConditions();
            RobinBcCoefStrategy<NDIM>* live_p_bc_coef = navier_stokes_integrator->getProjectionBoundaryConditions();
            const TimeSteppingType ib_time_stepping_type = IBAMR::string_to_enum<TimeSteppingType>(
                input_db->getStringWithDefault("IB_TIME_STEPPING", "MIDPOINT_RULE"));

            fac_pc->setVelocityPoissonSpecifications(U_problem_coefs);
            fac_pc->setPhysicalBcCoefs(live_u_bc_coefs, live_p_bc_coef);
            fac_pc->setPhysicalBoundaryHelper(bc_helper);
            fac_pc->setTimeInterval(current_time, new_time);
            fac_pc->setSolutionTime(new_time);
            fac_pc->setHomogeneousBc(true);
            fac_pc->setComponentsHaveNullSpace(false, true);
            fac_pc->setIBTimeSteppingType(ib_time_stepping_type);
            fac_pc->setIBForceJacobian(A);
            fac_pc->setIBInterpOp(J);
            fac_pc->setIBImplicitStrategy(ib_method_ops);
            fac_pc->initializeSolverState(*live_sol_vec, *live_rhs_vec);
            const bool live_fac_success = fac_pc->solveSystem(*live_sol_vec, *live_rhs_vec);

            pout << "live_fac_probe_success = " << static_cast<int>(live_fac_success) << "\n";
            pout << "live_fac_probe_type = " << vcycle_probe_type << "\n";
            pout << "live_fac_probe_scale = " << std::setprecision(17) << live_fac_probe_scale << "\n";
            export_level_patch_vec("live_fac_probe_rhs", f_rhs_idx, g_rhs_idx, 0);
            export_level_patch_vec("live_fac_probe_sol", u_scratch_idx, p_scratch_idx, 0);
            for (int ln = 1; ln <= finest_ln; ++ln)
            {
                export_level_patch_vec("live_fac_probe_rhs", f_rhs_idx, g_rhs_idx, ln);
                export_level_patch_vec("live_fac_probe_sol", u_scratch_idx, p_scratch_idx, ln);
            }

            fac_pc->deallocateSolverState();
        }

        if (run_matrix_krylov_smoke)
        {
            Vec matrix_krylov_rhs = nullptr;
            ierr = MatCreateVecs(level_l[static_cast<std::size_t>(finest_ln)], &matrix_krylov_rhs, nullptr);
            IBTK_CHKERRQ(ierr);
            ierr = MatMult(S, F_vec, matrix_krylov_rhs);
            IBTK_CHKERRQ(ierr);
            assert_vec_finite_nontrivial(matrix_krylov_rhs, "matrix GMRES rhs");

            if (!matrix_krylov_export_dir.empty())
            {
                std::filesystem::create_directories(matrix_krylov_export_dir);
                {
                    const int num_lag_nodes = compute_num_lag_nodes();
                    const double lag_ds = compute_lag_ds();
                    std::ofstream os(std::filesystem::path(matrix_krylov_export_dir) / "problem_setup.txt");
                    configure_double_stream(os);
                    os << "domain_length " << s_domain_length << "\n";
                    os << "finest_grid_cells " << s_finest_grid_cells << "\n";
                    os << "dx " << s_dx << "\n";
                    os << "dt " << s_dt << "\n";
                    os << "num_lag_nodes " << num_lag_nodes << "\n";
                    os << "lag_ds " << lag_ds << "\n";
                    os << "spring_stiffness_input " << s_spring_stiffness << "\n";
                    os << "spring_pair_k " << compute_pair_spring_stiffness() << "\n";
                    os << "spring_normalization "
                       << (s_spring_normalization == SpringNormalization::JuliaDensity ? "JULIA_DENSITY" :
                                                                                         "IBAMR_NODAL")
                       << "\n";
                    os << "matrix_spread_scale " << matrix_spread_scale << "\n";
                    os << "max_levels " << (finest_ln + 1) << "\n";
                    os << "base_grid_cells " << input_db->getInteger("N") << "\n";
                    os << "ref_ratio " << input_db->getIntegerWithDefault("REF_RATIO", 1) << "\n";
                }
                {
                    Vec matrix_krylov_rhs_export = nullptr;
                    ierr = VecDuplicate(matrix_krylov_rhs, &matrix_krylov_rhs_export);
                    IBTK_CHKERRQ(ierr);
                    ierr = VecCopy(matrix_krylov_rhs, matrix_krylov_rhs_export);
                    IBTK_CHKERRQ(ierr);
                    ierr = VecScale(matrix_krylov_rhs_export, matrix_spread_scale);
                    IBTK_CHKERRQ(ierr);
                    write_vec((std::filesystem::path(matrix_krylov_export_dir) / "matrix_krylov_rhs.txt").string(),
                              matrix_krylov_rhs_export);
                    ierr = VecDestroy(&matrix_krylov_rhs_export);
                    IBTK_CHKERRQ(ierr);
                }
                for (int ln = 0; ln <= finest_ln; ++ln)
                {
                    const std::size_t ln_idx = static_cast<std::size_t>(ln);
                    const std::string suffix = level_suffix(ln);
                    write_matrix_market(
                        (std::filesystem::path(matrix_krylov_export_dir) / ("L" + suffix + ".mtx")).string(),
                        level_l[ln_idx]);
                    write_eulerian_dof_map_ref(
                        (std::filesystem::path(matrix_krylov_export_dir) / ("eulerian_dof_map" + suffix + ".txt"))
                            .string(),
                        level_eulerian_dof_maps[ln_idx]);
                    write_int_vector(
                        (std::filesystem::path(matrix_krylov_export_dir) / ("pressure_dofs" + suffix + ".txt"))
                            .string(),
                        level_pressure_dofs[ln_idx]);
                    if (ln > 0)
                    {
                        write_matrix_market(
                            (std::filesystem::path(matrix_krylov_export_dir) / ("R" + suffix + ".mtx")).string(),
                            level_restrict[ln_idx]);
                        write_matrix_market(
                            (std::filesystem::path(matrix_krylov_export_dir) / ("P" + suffix + ".mtx")).string(),
                            level_prolong[ln_idx]);
                    }
                }
            }

            ierr = VecScale(matrix_krylov_rhs, matrix_spread_scale);
            IBTK_CHKERRQ(ierr);

            const auto& matrix_krylov_overlap =
                (matrix_krylov_policy == "relaxed") ? level_overlap_relaxed : level_overlap_strict;
            const MatrixKrylovSummary matrix_summary = run_matrix_gmres_smoke(matrix_krylov_rhs,
                                                                              level_l,
                                                                              level_restrict,
                                                                              level_prolong,
                                                                              matrix_krylov_overlap,
                                                                              level_pressure_dofs,
                                                                              num_pre_sweeps,
                                                                              num_post_sweeps,
                                                                              shell_pc_relaxation_factor,
                                                                              matrix_krylov_tol,
                                                                              matrix_krylov_max_it,
                                                                              matrix_krylov_restart);

            pout << "matrix_gmres_policy = " << matrix_krylov_policy << "\n";
            pout << "matrix_gmres_converged = " << static_cast<int>(matrix_summary.converged) << "\n";
            pout << "matrix_gmres_iterations = " << matrix_summary.iterations << "\n";
            pout << "matrix_gmres_converged_reason = " << matrix_summary.converged_reason << "\n";
            pout << "matrix_gmres_residual_norm = " << std::setprecision(17) << matrix_summary.residual_norm << "\n";
            pout << "matrix_gmres_relative_residual = " << std::setprecision(17)
                 << matrix_summary.relative_residual_norm << "\n";

            std::filesystem::path summary_path = matrix_krylov_summary_file;
            if (!summary_path.is_absolute() && export_bridge_data)
            {
                summary_path = std::filesystem::path(export_bridge_dir) / summary_path;
            }
            if (summary_path.has_parent_path())
            {
                std::filesystem::create_directories(summary_path.parent_path());
            }
            std::ofstream summary_os(summary_path);
            configure_double_stream(summary_os);
            summary_os << "policy " << matrix_krylov_policy << "\n";
            summary_os << "converged " << static_cast<int>(matrix_summary.converged) << "\n";
            summary_os << "iterations " << matrix_summary.iterations << "\n";
            summary_os << "converged_reason " << matrix_summary.converged_reason << "\n";
            summary_os << "residual_norm " << matrix_summary.residual_norm << "\n";
            summary_os << "relative_residual_norm " << matrix_summary.relative_residual_norm << "\n";
            summary_os << "tol " << matrix_krylov_tol << "\n";
            summary_os << "max_it " << matrix_krylov_max_it << "\n";
            summary_os << "restart " << matrix_krylov_restart << "\n";
            summary_os << "num_pre_sweeps " << num_pre_sweeps << "\n";
            summary_os << "num_post_sweeps " << num_post_sweeps << "\n";
            summary_os << "alpha " << shell_pc_relaxation_factor << "\n";
            summary_os << "finest_level " << finest_ln << "\n";
            summary_os << "finest_grid_cells " << s_finest_grid_cells << "\n";
            summary_os << "spring_stiffness_input " << s_spring_stiffness << "\n";
            summary_os.close();

            ierr = VecDestroy(&matrix_krylov_rhs);
            IBTK_CHKERRQ(ierr);
        }

        if (export_bridge_data)
        {
            {
                std::ofstream os(export_bridge_dir + "/smoother_params.txt");
                os << num_pre_sweeps << " " << num_post_sweeps << " ";
                configure_double_stream(os);
                os << shell_pc_relaxation_factor << "\n";
            }
            {
                const int num_lag_nodes = compute_num_lag_nodes();
                const double lag_ds = (2.0 * M_PI) / static_cast<double>(num_lag_nodes);
                std::ofstream os(export_bridge_dir + "/problem_setup.txt");
                configure_double_stream(os);
                os << "domain_length " << s_domain_length << "\n";
                os << "finest_grid_cells " << s_finest_grid_cells << "\n";
                os << "dx " << s_dx << "\n";
                os << "dt " << s_dt << "\n";
                os << "num_lag_nodes " << num_lag_nodes << "\n";
                os << "lag_ds " << lag_ds << "\n";
                os << "spring_stiffness_input " << s_spring_stiffness << "\n";
                os << "spring_pair_k " << s_spring_stiffness / (lag_ds * lag_ds) << "\n";
            }
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                const std::size_t ln_idx = static_cast<std::size_t>(ln);
                const std::string suffix = level_suffix(ln);
                write_matrix_market(export_bridge_dir + "/L" + suffix + ".mtx", level_l[ln_idx]);
                write_matrix_row_pattern(export_bridge_dir + "/L" + suffix + ".row_pattern.txt", level_l[ln_idx]);
                write_matrix_market(export_bridge_dir + "/A00" + suffix + ".mtx", level_a00[ln_idx]);
                write_matrix_row_pattern(export_bridge_dir + "/A00" + suffix + ".row_pattern.txt", level_a00[ln_idx]);
                write_matrix_market(export_bridge_dir + "/SAJ" + suffix + ".mtx", level_saj[ln_idx]);

                write_int_vector(export_bridge_dir + "/seed_velocity_dofs_axis0" + suffix + ".txt",
                                 level_seed_velocity_dofs_axis0[ln_idx]);
                write_int_vector(export_bridge_dir + "/num_dofs_per_proc" + suffix + ".txt",
                                 level_num_dofs_per_proc[ln_idx]);
                write_int_vector(export_bridge_dir + "/velocity_dofs" + suffix + ".txt", level_velocity_dofs[ln_idx]);
                write_int_vector(export_bridge_dir + "/pressure_dofs" + suffix + ".txt", level_pressure_dofs[ln_idx]);

                write_set_vector(export_bridge_dir + "/ibamr_overlap_relaxed" + suffix + ".txt",
                                 level_overlap_relaxed[ln_idx]);
                write_set_vector(export_bridge_dir + "/ibamr_overlap_strict" + suffix + ".txt",
                                 level_overlap_strict[ln_idx]);
                write_set_vector(export_bridge_dir + "/ibamr_nonoverlap_relaxed" + suffix + ".txt",
                                 level_nonoverlap_relaxed[ln_idx]);
                write_set_vector(export_bridge_dir + "/ibamr_nonoverlap_strict" + suffix + ".txt",
                                 level_nonoverlap_strict[ln_idx]);
                write_axis_map(export_bridge_dir + "/velocity_dof_to_component_axis" + suffix + ".txt",
                               level_map_data[ln_idx].velocity_dof_to_component_axis);
                write_set_map(export_bridge_dir + "/velocity_dof_to_adjacent_cell_dofs" + suffix + ".txt",
                              level_map_data[ln_idx].velocity_dof_to_adjacent_cell_dofs);
                write_set_map(export_bridge_dir + "/cell_dof_to_closure_dofs" + suffix + ".txt",
                              level_map_data[ln_idx].cell_dof_to_closure_dofs);
                write_set_map(export_bridge_dir + "/velocity_dof_to_paired_seed_velocity_dofs" + suffix + ".txt",
                              level_map_data[ln_idx].velocity_dof_to_paired_seed_velocity_dofs);

                write_overlap_submatrix_blocks(export_bridge_dir + "/ibamr_overlap_submat_relaxed" + suffix + ".txt",
                                               level_a00[ln_idx],
                                               level_overlap_relaxed[ln_idx]);
                write_overlap_submatrix_blocks(export_bridge_dir + "/ibamr_overlap_submat_strict" + suffix + ".txt",
                                               level_a00[ln_idx],
                                               level_overlap_strict[ln_idx]);
                write_overlap_submatrix_blocks(export_bridge_dir + "/ibamr_overlap_submat_full_saddle_relaxed" +
                                                   suffix + ".txt",
                                               level_l[ln_idx],
                                               level_overlap_relaxed[ln_idx]);
                write_overlap_submatrix_blocks(export_bridge_dir + "/ibamr_overlap_submat_full_saddle_strict" + suffix +
                                                   ".txt",
                                               level_l[ln_idx],
                                               level_overlap_strict[ln_idx]);

                Vec smoother_input = nullptr;
                ierr = MatCreateVecs(level_l[ln_idx], &smoother_input, nullptr);
                IBTK_CHKERRQ(ierr);
                Vec smoother_state = nullptr;
                ierr = VecDuplicate(smoother_input, &smoother_state);
                IBTK_CHKERRQ(ierr);
                fill_nonzero_vector(smoother_state, 1.0 + 0.1 * static_cast<double>(ln));
                subtract_pressure_mean(smoother_state, level_pressure_dofs[ln_idx]);
                ierr = MatMult(level_l[ln_idx], smoother_state, smoother_input);
                IBTK_CHKERRQ(ierr);
                assert_vec_finite_nontrivial(smoother_input, "smoother input");
                Vec smoother_relaxed = nullptr;
                Vec smoother_strict = nullptr;
                ierr = VecDuplicate(smoother_input, &smoother_relaxed);
                IBTK_CHKERRQ(ierr);
                ierr = VecDuplicate(smoother_input, &smoother_strict);
                IBTK_CHKERRQ(ierr);
                ierr = VecZeroEntries(smoother_relaxed);
                IBTK_CHKERRQ(ierr);
                ierr = VecZeroEntries(smoother_strict);
                IBTK_CHKERRQ(ierr);
                std::vector<SweepTraceFrame> smoother_trace_relaxed;
                std::vector<SweepTraceFrame> smoother_trace_strict;
                std::vector<SweepResidualTraceFrame> smoother_residual_trace_relaxed;
                std::vector<SweepResidualTraceFrame> smoother_residual_trace_strict;
                std::vector<SweepDeltaTraceFrame> smoother_delta_trace_relaxed;
                std::vector<SweepDeltaTraceFrame> smoother_delta_trace_strict;
                LocalSolveDebugDump relaxed_debug_dump;
                relaxed_debug_dump.enabled = debug_local_solve_dump && ln == debug_local_solve_dump_level &&
                                             debug_local_solve_dump_policy == "relaxed";
                relaxed_debug_dump.target_sweep = debug_local_solve_dump_sweep;
                relaxed_debug_dump.target_subdomain = debug_local_solve_dump_subdomain;
                relaxed_debug_dump.target_global_row = debug_local_solve_dump_global_row;
                relaxed_debug_dump.output_dir = debug_local_solve_dump_dir;
                LocalSolveDebugDump strict_debug_dump;
                strict_debug_dump.enabled = debug_local_solve_dump && ln == debug_local_solve_dump_level &&
                                            debug_local_solve_dump_policy == "strict";
                strict_debug_dump.target_sweep = debug_local_solve_dump_sweep;
                strict_debug_dump.target_subdomain = debug_local_solve_dump_subdomain;
                strict_debug_dump.target_global_row = debug_local_solve_dump_global_row;
                strict_debug_dump.output_dir = debug_local_solve_dump_dir;
                apply_coupling_aware_sweep(smoother_relaxed,
                                           smoother_input,
                                           level_l[ln_idx],
                                           level_overlap_relaxed[ln_idx],
                                           level_pressure_dofs[ln_idx],
                                           1,
                                           shell_pc_relaxation_factor,
                                           &smoother_trace_relaxed,
                                           &smoother_residual_trace_relaxed,
                                           &smoother_delta_trace_relaxed,
                                           &relaxed_debug_dump);
                apply_coupling_aware_sweep(smoother_strict,
                                           smoother_input,
                                           level_l[ln_idx],
                                           level_overlap_strict[ln_idx],
                                           level_pressure_dofs[ln_idx],
                                           1,
                                           shell_pc_relaxation_factor,
                                           &smoother_trace_strict,
                                           &smoother_residual_trace_strict,
                                           &smoother_delta_trace_strict,
                                           &strict_debug_dump);
                assert_vec_finite_nontrivial(smoother_relaxed, "smoother relaxed output");
                assert_vec_finite_nontrivial(smoother_strict, "smoother strict output");

                write_vec(export_bridge_dir + "/smoother_input" + suffix + ".txt", smoother_input);
                write_vec(export_bridge_dir + "/smoother_output_relaxed" + suffix + ".txt", smoother_relaxed);
                write_vec(export_bridge_dir + "/smoother_output_strict" + suffix + ".txt", smoother_strict);
                write_sweep_trace(export_bridge_dir + "/smoother_trace_relaxed" + suffix + ".txt",
                                  smoother_trace_relaxed);
                write_sweep_trace(export_bridge_dir + "/smoother_trace_strict" + suffix + ".txt",
                                  smoother_trace_strict);
                write_sweep_residual_trace(export_bridge_dir + "/smoother_residual_trace_relaxed" + suffix + ".txt",
                                           smoother_residual_trace_relaxed);
                write_sweep_residual_trace(export_bridge_dir + "/smoother_residual_trace_strict" + suffix + ".txt",
                                           smoother_residual_trace_strict);
                write_sweep_delta_trace(export_bridge_dir + "/smoother_delta_trace_relaxed" + suffix + ".txt",
                                        smoother_delta_trace_relaxed);
                write_sweep_delta_trace(export_bridge_dir + "/smoother_delta_trace_strict" + suffix + ".txt",
                                        smoother_delta_trace_strict);

                ierr = VecDestroy(&smoother_strict);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&smoother_relaxed);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&smoother_input);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&smoother_state);
                IBTK_CHKERRQ(ierr);
            }

            for (int ln = 1; ln <= finest_ln; ++ln)
            {
                const std::size_t fine_idx = static_cast<std::size_t>(ln);
                const std::size_t coarse_idx = static_cast<std::size_t>(ln - 1);
                const std::string fine_suffix = level_suffix(ln);

                write_matrix_market(export_bridge_dir + "/P" + fine_suffix + ".mtx", level_prolong[fine_idx]);
                write_matrix_market(export_bridge_dir + "/R" + fine_suffix + ".mtx", level_restrict[fine_idx]);

                Vec fine_input = nullptr;
                ierr = MatCreateVecs(level_l[fine_idx], &fine_input, nullptr);
                IBTK_CHKERRQ(ierr);
                fill_nonzero_vector(fine_input, 1.6 + 0.1 * static_cast<double>(ln));
                Vec coarse_restricted = nullptr;
                ierr = MatCreateVecs(level_l[coarse_idx], &coarse_restricted, nullptr);
                IBTK_CHKERRQ(ierr);
                ierr = MatMult(level_restrict[fine_idx], fine_input, coarse_restricted);
                IBTK_CHKERRQ(ierr);
                assert_vec_finite_nontrivial(fine_input, "restriction input fine");
                assert_vec_finite_nontrivial(coarse_restricted, "restriction output coarse");
                write_vec(export_bridge_dir + "/restriction_input_fine" + fine_suffix + ".txt", fine_input);
                write_vec(export_bridge_dir + "/restriction_output_coarse" + fine_suffix + ".txt", coarse_restricted);

                Vec coarse_input = nullptr;
                ierr = VecDuplicate(coarse_restricted, &coarse_input);
                IBTK_CHKERRQ(ierr);
                fill_nonzero_vector(coarse_input, 1.9 + 0.1 * static_cast<double>(ln));
                Vec fine_prolonged = nullptr;
                ierr = VecDuplicate(fine_input, &fine_prolonged);
                IBTK_CHKERRQ(ierr);
                ierr = MatMult(level_prolong[fine_idx], coarse_input, fine_prolonged);
                IBTK_CHKERRQ(ierr);
                assert_vec_finite_nontrivial(coarse_input, "prolongation input coarse");
                assert_vec_finite_nontrivial(fine_prolonged, "prolongation output fine");
                write_vec(export_bridge_dir + "/prolongation_input_coarse" + fine_suffix + ".txt", coarse_input);
                write_vec(export_bridge_dir + "/prolongation_output_fine" + fine_suffix + ".txt", fine_prolonged);

                Vec vcycle_rhs = nullptr;
                ierr = VecDuplicate(fine_input, &vcycle_rhs);
                IBTK_CHKERRQ(ierr);
                Vec vcycle_state = nullptr;
                ierr = VecDuplicate(vcycle_rhs, &vcycle_state);
                IBTK_CHKERRQ(ierr);
                fill_vcycle_probe_state(vcycle_state,
                                        vcycle_probe_type,
                                        level_eulerian_dof_maps[fine_idx],
                                        2.3 + 0.1 * static_cast<double>(ln),
                                        vcycle_probe_mode_x,
                                        vcycle_probe_mode_y);
                subtract_pressure_mean(vcycle_state, level_pressure_dofs[fine_idx]);
                ierr = MatMult(level_l[fine_idx], vcycle_state, vcycle_rhs);
                IBTK_CHKERRQ(ierr);
                Vec vcycle_relaxed = nullptr;
                Vec vcycle_strict = nullptr;
                ierr = VecDuplicate(vcycle_rhs, &vcycle_relaxed);
                IBTK_CHKERRQ(ierr);
                ierr = VecDuplicate(vcycle_rhs, &vcycle_strict);
                IBTK_CHKERRQ(ierr);
                VCycleTrace vcycle_trace_relaxed;
                VCycleTrace vcycle_trace_strict;
                apply_two_level_vcycle(vcycle_relaxed,
                                       vcycle_rhs,
                                       level_l[fine_idx],
                                       level_l[coarse_idx],
                                       level_restrict[fine_idx],
                                       level_prolong[fine_idx],
                                       level_overlap_relaxed[fine_idx],
                                       level_pressure_dofs[fine_idx],
                                       num_pre_sweeps,
                                       num_post_sweeps,
                                       shell_pc_relaxation_factor,
                                       &vcycle_trace_relaxed);
                apply_two_level_vcycle(vcycle_strict,
                                       vcycle_rhs,
                                       level_l[fine_idx],
                                       level_l[coarse_idx],
                                       level_restrict[fine_idx],
                                       level_prolong[fine_idx],
                                       level_overlap_strict[fine_idx],
                                       level_pressure_dofs[fine_idx],
                                       num_pre_sweeps,
                                       num_post_sweeps,
                                       shell_pc_relaxation_factor,
                                       &vcycle_trace_strict);
                assert_vec_finite_nontrivial(vcycle_rhs, "v-cycle input rhs");
                assert_vec_finite_nontrivial(vcycle_relaxed, "v-cycle relaxed output");
                assert_vec_finite_nontrivial(vcycle_strict, "v-cycle strict output");
                write_vec(export_bridge_dir + "/vcycle_input_state" + fine_suffix + ".txt", vcycle_state);
                write_vec(export_bridge_dir + "/vcycle_input_rhs" + fine_suffix + ".txt", vcycle_rhs);
                write_vec(export_bridge_dir + "/vcycle_output_relaxed" + fine_suffix + ".txt", vcycle_relaxed);
                write_vec(export_bridge_dir + "/vcycle_output_strict" + fine_suffix + ".txt", vcycle_strict);
                write_double_vector(export_bridge_dir + "/vcycle_trace_relaxed_y_pre" + fine_suffix + ".txt",
                                    vcycle_trace_relaxed.y_pre);
                write_double_vector(export_bridge_dir + "/vcycle_trace_relaxed_residual_fine" + fine_suffix + ".txt",
                                    vcycle_trace_relaxed.residual_fine);
                write_double_vector(export_bridge_dir + "/vcycle_trace_relaxed_coarse_rhs" + fine_suffix + ".txt",
                                    vcycle_trace_relaxed.coarse_rhs);
                write_double_vector(export_bridge_dir + "/vcycle_trace_relaxed_coarse_sol" + fine_suffix + ".txt",
                                    vcycle_trace_relaxed.coarse_sol);
                write_double_vector(export_bridge_dir + "/vcycle_trace_relaxed_fine_correction" + fine_suffix + ".txt",
                                    vcycle_trace_relaxed.fine_correction);
                write_double_vector(export_bridge_dir + "/vcycle_trace_relaxed_y_after_correction" + fine_suffix +
                                        ".txt",
                                    vcycle_trace_relaxed.y_after_correction);
                write_double_vector(export_bridge_dir + "/vcycle_trace_relaxed_y_post" + fine_suffix + ".txt",
                                    vcycle_trace_relaxed.y_post);
                write_double_vector(export_bridge_dir + "/vcycle_trace_strict_y_pre" + fine_suffix + ".txt",
                                    vcycle_trace_strict.y_pre);
                write_double_vector(export_bridge_dir + "/vcycle_trace_strict_residual_fine" + fine_suffix + ".txt",
                                    vcycle_trace_strict.residual_fine);
                write_double_vector(export_bridge_dir + "/vcycle_trace_strict_coarse_rhs" + fine_suffix + ".txt",
                                    vcycle_trace_strict.coarse_rhs);
                write_double_vector(export_bridge_dir + "/vcycle_trace_strict_coarse_sol" + fine_suffix + ".txt",
                                    vcycle_trace_strict.coarse_sol);
                write_double_vector(export_bridge_dir + "/vcycle_trace_strict_fine_correction" + fine_suffix + ".txt",
                                    vcycle_trace_strict.fine_correction);
                write_double_vector(export_bridge_dir + "/vcycle_trace_strict_y_after_correction" + fine_suffix +
                                        ".txt",
                                    vcycle_trace_strict.y_after_correction);
                write_double_vector(export_bridge_dir + "/vcycle_trace_strict_y_post" + fine_suffix + ".txt",
                                    vcycle_trace_strict.y_post);

                ierr = VecDestroy(&vcycle_strict);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&vcycle_relaxed);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&vcycle_rhs);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&vcycle_state);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&fine_prolonged);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&coarse_input);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&coarse_restricted);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&fine_input);
                IBTK_CHKERRQ(ierr);
            }

            pout << "wrote_bridge_data_dir = " << export_bridge_dir << "\n";
        }

        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            const std::size_t ln_idx = static_cast<std::size_t>(ln);
            ierr = MatDestroy(&level_a00[ln_idx]);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&level_l[ln_idx]);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&level_stokes[ln_idx]);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&level_saj[ln_idx]);
            IBTK_CHKERRQ(ierr);
            if (ln > 0)
            {
                ierr = MatDestroy(&level_restrict[ln_idx]);
                IBTK_CHKERRQ(ierr);
                ierr = MatDestroy(&level_prolong[ln_idx]);
                IBTK_CHKERRQ(ierr);
            }
        }
    }

    int test_failures = 0;
    std::ostringstream diag;
    const double a_abs_nz_tol = std::max(1.0e-12, 100.0 * parity_tol);
    const double a_rel_nz_tol = std::max(1.0e-12, 100.0 * parity_tol);
    check_lagrangian_a_row_nnz(A_data, a_abs_nz_tol, a_rel_nz_tol, diag, test_failures);
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

    plog << "Input database:\n";
    input_db->printClassData(plog);
    pout << "test_failures = " << test_failures << "\n";
    return test_failures > 0 ? 1 : 0;
}
