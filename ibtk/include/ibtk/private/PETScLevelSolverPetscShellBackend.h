// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_private_PETScLevelSolverPetscShellBackend
#define included_IBTK_private_PETScLevelSolverPetscShellBackend

#include <petscksp.h>

#include <memory>
#include <vector>

namespace IBTK
{
class PETScLevelSolver;

class PETScLevelSolverPetscShellBackend
{
public:
    explicit PETScLevelSolverPetscShellBackend(PETScLevelSolver& solver);

    void initialize();
    void deallocate();
    void beginAccumulateCorrection(int subdomain_num, Vec sub_y, Vec y);
    void endAccumulateCorrection(int subdomain_num, Vec sub_y, Vec y);
    void accumulateCorrection(int subdomain_num, Vec sub_y, Vec y);
    void applyAdditive(Vec x, Vec y);
    void applyMultiplicative(Vec x, Vec y);

private:
    struct Data
    {
        Vec shell_r = nullptr;
        InsertMode prolongation_insert_mode = INSERT_VALUES;
        IS owned_residual_update_rows_is = nullptr;
        std::vector<IS> global_overlap_is, global_nonoverlap_is;
        std::vector<IS> local_overlap_is, local_nonoverlap_is;
        std::vector<VecScatter> restriction, prolongation;
        std::vector<KSP> sub_ksp;
        Mat* sub_mat = nullptr;
        Mat* active_residual_update_mat = nullptr;
        std::vector<Vec> sub_x, sub_y, active_residual_update_x, active_residual_update_y;
        std::vector<std::vector<int>> active_update_local_positions;
    };

    void updateResidual(int subdomain_num, Vec sub_y, Vec residual);

    PETScLevelSolver& d_solver;
    std::unique_ptr<Data> d_data;
};
} // namespace IBTK

#endif
