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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/CCScaledPoissonPETScLevelSolver.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/PoissonUtilities.h>

#include <tbox/Array.h>

#include <petscmat.h>
#include <petscvec.h>

#include <BoundaryBox.h>
#include <CellData.h>
#include <CellDataFactory.h>
#include <CellIndex.h>
#include <CoarseFineBoundary.h>
#include <Patch.h>
#include <PatchDescriptor.h>
#include <PatchGeometry.h>
#include <PatchLevel.h>
#include <PoissonSpecifications.h>
#include <SideDataFactory.h>
#include <VariableDatabase.h>

#include <cmath>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

CCScaledPoissonPETScLevelSolver::CCScaledPoissonPETScLevelSolver(const std::string& object_name,
                                                                 Pointer<Database> input_db,
                                                                 const std::string& default_options_prefix)
    : CCPoissonPETScLevelSolver(object_name, input_db, default_options_prefix)
{
    return;
} // CCScaledPoissonPETScLevelSolver

CCScaledPoissonPETScLevelSolver::~CCScaledPoissonPETScLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~CCScaledPoissonPETScLevelSolver

void
CCScaledPoissonPETScLevelSolver::setScalePatchDataIndex(const int scale_idx)
{
#if !defined(NDEBUG)
    if (scale_idx < 0)
    {
        TBOX_ERROR(d_object_name << "::setScalePatchDataIndex():\n"
                                 << "  scale_idx must be a valid patch data index.\n");
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<CellDataFactory<NDIM, double>> factory = var_db->getPatchDescriptor()->getPatchDataFactory(scale_idx);

    if (!factory || factory->getDefaultDepth() != 1)
    {
        TBOX_ERROR(d_object_name << "::setScalePatchDataIndex():\n"
                                 << "  scale must refer to scalar cell-centered double data.\n");
    }
#endif
    d_scale_idx = scale_idx;

    return;
} // setScalePatchDataIndex

/////////////////////////////// PROTECTED ////////////////////////////////////

void
CCScaledPoissonPETScLevelSolver::createOperatorSpecialized()
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    Pointer<CellDataFactory<NDIM, int>> dof_factory =
        var_db->getPatchDescriptor()->getPatchDataFactory(d_dof_index_idx);

    if (!dof_factory)
    {
        TBOX_ERROR(d_object_name << "::createOperatorSpecialized():\n"
                                 << "  DOF index data have not been configured.\n");
    }

    // Follow the coefficient-depth convention of the existing matrix utilities.
    if (!d_poisson_spec.dIsConstant())
    {
        Pointer<SideDataFactory<NDIM, double>> D_factory =
            var_db->getPatchDescriptor()->getPatchDataFactory(d_poisson_spec.getDPatchDataId());

        if (!D_factory || D_factory->getDefaultDepth() != dof_factory->getDefaultDepth())
        {
            TBOX_ERROR(d_object_name << "::createOperatorSpecialized():\n"
                                     << "  variable D must be side-centered double data with depth "
                                        "matching the solution.\n");
        }
    }

    // Construct only the diffusion matrix, including homogeneous Robin elimination.
    PoissonSpecifications diffusion_spec(d_poisson_spec);
    diffusion_spec.setCZero();

    PETScMatUtilities::constructPatchLevelCCLaplaceOp(
        d_petsc_mat, diffusion_spec, d_bc_coefs, d_solution_time, d_num_dofs_per_proc, d_dof_index_idx, d_level);

    Vec scale_vec = nullptr;
    Vec C_vec = nullptr;

    int ierr = VecDuplicate(d_petsc_x, &scale_vec);
    IBTK_CHKERRQ(ierr);

    ierr = VecDuplicate(d_petsc_x, &C_vec);
    IBTK_CHKERRQ(ierr);

    copyCoefficientToPETScVec(scale_vec, d_scale_idx, 1.0);

    const int C_idx = d_poisson_spec.cIsVariable() ? d_poisson_spec.getCPatchDataId() : invalid_index;
    const double C_constant =
        d_poisson_spec.cIsVariable() || d_poisson_spec.cIsZero() ? 0.0 : d_poisson_spec.getCConstant();

    copyCoefficientToPETScVec(C_vec, C_idx, C_constant);

    // Left scaling multiplies rows: A = diag(S)*L_D.
    ierr = MatDiagonalScale(d_petsc_mat, scale_vec, nullptr);
    IBTK_CHKERRQ(ierr);

    // Add C after scaling, so C is not multiplied by S.
    ierr = MatDiagonalSet(d_petsc_mat, C_vec, ADD_VALUES);
    IBTK_CHKERRQ(ierr);

    ierr = MatAssemblyBegin(d_petsc_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    ierr = MatAssemblyEnd(d_petsc_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    ierr = VecDestroy(&scale_vec);
    IBTK_CHKERRQ(ierr);

    ierr = VecDestroy(&C_vec);
    IBTK_CHKERRQ(ierr);

    // The parent sets d_petsc_pc = d_petsc_mat after this hook returns.
    return;
} // createOperatorSpecialized

void
CCScaledPoissonPETScLevelSolver::setupKSPVecs(Vec& petsc_x,
                                              Vec& petsc_b,
                                              SAMRAIVectorReal<NDIM, double>& x,
                                              SAMRAIVectorReal<NDIM, double>& b)
{
    if (d_initial_guess_nonzero) copyToPETScVec(petsc_x, x);

    const int x_idx = x.getComponentDescriptorIndex(0);
    const int b_idx = b.getComponentDescriptorIndex(0);
    const auto b_adj_idx = d_cached_eulerian_data.getCachedPatchDataIndex(b_idx);

    PoissonSpecifications diffusion_spec(d_poisson_spec);
    diffusion_spec.setCZero();

    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
        Pointer<PatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();

        Pointer<CellData<NDIM, double>> x_data = patch->getPatchData(x_idx);
        Pointer<CellData<NDIM, double>> b_data = patch->getPatchData(b_idx);
        Pointer<CellData<NDIM, double>> b_adj_data = patch->getPatchData(b_adj_idx);

        Pointer<CellData<NDIM, double>> scale_data;
        if (d_scale_idx != invalid_index) scale_data = patch->getPatchData(d_scale_idx);

        // First accumulate only the unscaled boundary contribution.
        b_adj_data->fillAll(0.0);

        if (pgeom->intersectsPhysicalBoundary())
        {
            PoissonUtilities::adjustRHSAtPhysicalBoundary(
                *b_adj_data, patch, diffusion_spec, d_bc_coefs, d_solution_time, d_homogeneous_bc);
        }

        const Array<BoundaryBox<NDIM>>& type_1_cf_bdry =
            d_level_num == 0 ? Array<BoundaryBox<NDIM>>() : d_cf_boundary->getBoundaries(patch->getPatchNumber(), 1);

        if (type_1_cf_bdry.size() > 0)
        {
            PoissonUtilities::adjustRHSAtCoarseFineBoundary(
                *b_adj_data, *x_data, patch, diffusion_spec, type_1_cf_bdry);
        }

        // b_adjusted = b + S*(physical and coarse-fine boundary contributions).
        const Box<NDIM>& patch_box = patch->getBox();

        for (Box<NDIM>::Iterator it(patch_box); it; it++)
        {
            const CellIndex<NDIM> ci(it());
            const double scale = scale_data ? (*scale_data)(ci) : 1.0;

            for (int depth = 0; depth < b_data->getDepth(); ++depth)
            {
                if (scale == 0.0)
                    (*b_adj_data)(ci, depth) = (*b_data)(ci, depth);
                else
                    (*b_adj_data)(ci, depth) = (*b_data)(ci, depth) + scale * (*b_adj_data)(ci, depth);
            }
        }
    }

    PETScVecUtilities::copyToPatchLevelVec(petsc_b, b_adj_idx, d_dof_index_idx, d_level);

    return;
} // setupKSPVecs

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CCScaledPoissonPETScLevelSolver::copyCoefficientToPETScVec(Vec vec,
                                                           const int coefficient_idx,
                                                           const double constant) const
{
    if (coefficient_idx == invalid_index)
    {
        const int ierr = VecSet(vec, constant);
        IBTK_CHKERRQ(ierr);
        return;
    }

    int ierr = VecSet(vec, 0.0);
    IBTK_CHKERRQ(ierr);

    PetscInt first_row = 0;
    PetscInt last_row = 0;

    ierr = VecGetOwnershipRange(vec, &first_row, &last_row);
    IBTK_CHKERRQ(ierr);

    PetscScalar* values = nullptr;
    ierr = VecGetArray(vec, &values);
    IBTK_CHKERRQ(ierr);

    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
        Pointer<CellData<NDIM, double>> coefficient_data = patch->getPatchData(coefficient_idx);
        Pointer<CellData<NDIM, int>> dof_data = patch->getPatchData(d_dof_index_idx);
        const Box<NDIM>& patch_box = patch->getBox();

        for (Box<NDIM>::Iterator it(patch_box); it; it++)
        {
            const CellIndex<NDIM> ci(it());
            const double value = (*coefficient_data)(ci);

            for (int depth = 0; depth < dof_data->getDepth(); ++depth)
            {
                const PetscInt row = (*dof_data)(ci, depth);

                if (first_row <= row && row < last_row)
                {
                    values[row - first_row] = value;
                }
            }
        }
    }

    ierr = VecRestoreArray(vec, &values);
    IBTK_CHKERRQ(ierr);

    return;
} // copyCoefficientToPETScVec

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK
