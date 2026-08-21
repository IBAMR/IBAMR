// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_cav_raw_operator_comparator
#define included_cav_raw_operator_comparator

#include <ibamr/StaggeredStokesPETScLevelSolver.h>

#include <tbox/Pointer.h>

#include <petscmat.h>
#include <petscvec.h>

#include <array>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
} // namespace hier
} // namespace SAMRAI

namespace IBAMR
{
namespace TestSupport
{
enum class CAVRawDofKind
{
    VELOCITY,
    PRESSURE
};

struct CAVRawDofRecord
{
    PetscInt dof = -1;
    CAVRawDofKind kind = CAVRawDofKind::VELOCITY;
    int axis = -1;
    std::array<int, NDIM> index = {};
    std::array<double, NDIM> position = {};
};

struct CAVRawMatrixEntry
{
    PetscInt row = -1;
    PetscInt column = -1;
    double value = 0.0;
};

struct CAVRawMatrixMarketData
{
    PetscInt nrows = 0;
    PetscInt ncols = 0;
    std::vector<CAVRawMatrixEntry> entries;
};

struct CAVRawDofMapData
{
    int dimension = NDIM;
    PetscInt size = 0;
    std::vector<CAVRawDofRecord> dofs;
};

struct CAVRawOperatorBundle
{
    int dimension = NDIM;
    PetscInt nrows = 0;
    PetscInt ncols = 0;
    std::vector<CAVRawDofRecord> dofs;
    std::vector<CAVRawMatrixEntry> matrix_entries;
    std::vector<double> equation_values;
    std::vector<double> state_values;
    std::vector<PetscInt> ordered_dofs;
};

struct CAVLiveExportManifest
{
    std::string candidate_sha;
    bool candidate_dirty = false;
    std::string oracle_sha;
    std::string case_id;
    int dimension = NDIM;
    int mpi_ranks = 1;
    std::string pressure_equation;
    double pressure_equation_row_multiplier_to_oracle = 1.0;
    std::string pressure_gauge;
    std::string patch_seed_type;
    std::string closure_policy;
    int seed_stride = 1;
    std::string traversal_order;
    std::string composition;
    std::string local_solver_backend;
};

struct CAVLocalSolveTraceRecord
{
    int sweep = -1;
    int patch_ordinal = -1;
    std::string artifact_stem;
};

struct CAVRawMappingSpec
{
    double coordinate_tolerance = 1.0e-14;
    double comparison_tolerance = 0.0;
    // This scale applies only to candidate pressure-equation rows and RHS entries. It never changes columns or
    // pressure-state values.
    double candidate_pressure_equation_scale = 1.0;
    // A pressure constant is the only legal state gauge adjustment.
    bool align_pressure_gauge = false;
};

struct CAVRawComparison
{
    bool dof_bijection = false;
    bool matrix_structure = false;
    bool matrix_values = false;
    bool equation_values = false;
    bool state_values = false;
    bool ordered_dofs = false;
    bool matched = false;
    double matrix_max_abs_error = 0.0;
    double matrix_scaled_error = 0.0;
    double equation_max_abs_error = 0.0;
    double equation_scaled_error = 0.0;
    double state_max_abs_error = 0.0;
    double state_scaled_error = 0.0;
    std::string first_mismatch;
};

struct CAVRawComparatorControlResults
{
    bool exact_control = false;
    bool declared_permutation_sign_gauge = false;
    bool undeclared_permutation_detected = false;
    bool velocity_sign_detected = false;
    bool pressure_sign_detected = false;
    bool legal_pressure_gauge = false;
    bool illegal_velocity_shift_detected = false;
    bool omitted_dof_detected = false;
    bool omitted_row_detected = false;
    bool ordered_omission_detected = false;
    bool ordered_reordering_detected = false;
};

CAVRawOperatorBundle captureCAVRawOperatorBundle(const StaggeredStokesPETScLevelSolver::LiveOperatorStateView& view,
                                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level,
                                                 int u_dof_index_idx,
                                                 int p_dof_index_idx,
                                                 Vec equation_values = nullptr,
                                                 Vec state_values = nullptr);

CAVRawDofMapData captureCAVRawDofMap(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level,
                                     int u_dof_index_idx,
                                     int p_dof_index_idx);

void writeCAVRawDofMap(const CAVRawDofMapData& dof_map, const std::string& filename);

CAVRawDofMapData readCAVRawDofMap(const std::string& filename);

bool sameCAVRawDofMap(const CAVRawDofMapData& lhs, const CAVRawDofMapData& rhs);

void writeCAVRawOperatorBundle(const CAVRawOperatorBundle& bundle, const std::string& prefix);

CAVRawOperatorBundle readCAVRawOperatorBundle(const std::string& prefix);

void writeCAVRawMatrixMarket(Mat matrix, const std::string& filename);

CAVRawMatrixMarketData readCAVRawMatrixMarket(const std::string& filename);

bool sameCAVRawMatrixMarket(Mat matrix, const CAVRawMatrixMarketData& data);

void writeCAVRawVectorMarket(Vec vector, const std::string& filename);

std::vector<double> readCAVRawVectorMarket(const std::string& filename);

bool sameCAVRawVectorMarket(Vec vector, const std::vector<double>& values);

void writeCAVRawIndexList(const std::vector<int>& indices, const std::string& filename);

std::vector<int> readCAVRawIndexList(const std::string& filename);

void writeCAVRawIndexSets(const std::vector<std::vector<int>>& index_sets, const std::string& filename);

std::vector<std::vector<int>> readCAVRawIndexSets(const std::string& filename);

void writeCAVLiveExportManifest(const CAVLiveExportManifest& manifest, const std::string& filename);

CAVLiveExportManifest readCAVLiveExportManifest(const std::string& filename);

bool sameCAVLiveExportManifest(const CAVLiveExportManifest& lhs, const CAVLiveExportManifest& rhs);

void writeCAVLocalSolveTraceIndex(const std::vector<CAVLocalSolveTraceRecord>& records, const std::string& filename);

std::vector<CAVLocalSolveTraceRecord> readCAVLocalSolveTraceIndex(const std::string& filename);

bool sameCAVLocalSolveTraceIndex(const std::vector<CAVLocalSolveTraceRecord>& lhs,
                                 const std::vector<CAVLocalSolveTraceRecord>& rhs);

CAVRawComparison compareCAVRawOperatorBundles(const CAVRawOperatorBundle& candidate,
                                              const CAVRawOperatorBundle& reference,
                                              const CAVRawMappingSpec& mapping);

PetscInt
findCAVRawDof(const CAVRawOperatorBundle& bundle, CAVRawDofKind kind, int axis, const std::array<int, NDIM>& index);

double getCAVRawMatrixValue(const CAVRawOperatorBundle& bundle, PetscInt row, PetscInt column);

std::size_t countCAVRawMatrixRowNonzeros(const CAVRawOperatorBundle& bundle, PetscInt row);

CAVRawComparatorControlResults runCAVRawComparatorControls();

} // namespace TestSupport
} // namespace IBAMR

#endif
