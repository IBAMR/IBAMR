// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#include <ibamr/StaggeredStokesPETScVecUtilities.h>

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/SideSynchCopyFillPattern.h>
#include <ibtk/samrai_compatibility_names.h>

#include <petscao.h>
#include <petsclog.h>
#include <petscsys.h>
#include <petscvec.h>

#include <BoxArray.h>
#include <MultiblockDataTranslator.h>
#include <SAMRAIBox.h>
#include <SAMRAIBoxList.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAICellData.h>
#include <SAMRAICellGeometry.h>
#include <SAMRAICellIndex.h>
#include <SAMRAICellVariable.h>
#include <SAMRAIIndex.h>
#include <SAMRAIIntVector.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIPointer.h>
#include <SAMRAIRefineAlgorithm.h>
#include <SAMRAIRefineClasses.h>
#include <SAMRAIRefineOperator.h>
#include <SAMRAIRefineSchedule.h>
#include <SAMRAISideData.h>
#include <SAMRAISideGeometry.h>
#include <SAMRAISideIndex.h>
#include <SAMRAISideVariable.h>
#include <SAMRAIUtilities.h>
#include <SAMRAIVariable.h>
#include <SAMRAIVariableDatabase.h>
#include <SAMRAIVariableFillPattern.h>

#include <algorithm>
#include <array>
#include <numeric>
#include <ostream>
#include <string>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

// FORTRAN ROUTINES
#if (NDIM == 2)
#define COPY_TO_PATCHLEVEL_VEC_MAC_FC IBAMR_FC_FUNC_(copy_to_patchlevel_vec_mac2d, COPY_TO_PATCHLEVEL_VEC_MAC2D)

#define COPY_FROM_PATCHLEVEL_VEC_MAC_FC IBAMR_FC_FUNC_(copy_from_patchlevel_vec_mac2d, COPY_FROM_PATCHLEVEL_VEC_MAC2D)
#endif

#if (NDIM == 3)
#define COPY_TO_PATCHLEVEL_VEC_MAC_FC IBAMR_FC_FUNC_(copy_to_patchlevel_vec_mac3d, COPY_TO_PATCHLEVEL_VEC_MAC3D)

#define COPY_FROM_PATCHLEVEL_VEC_MAC_FC IBAMR_FC_FUNC_(copy_from_patchlevel_vec_mac3d, COPY_FROM_PATCHLEVEL_VEC_MAC3D)
#endif

extern "C"
{
    void COPY_TO_PATCHLEVEL_VEC_MAC_FC(const int&,
                                       const int&,
                                       const int&,
                                       const int&,
#if (NDIM == 3)
                                       const int&,
                                       const int&,
#endif
                                       const int&,
                                       const int&,
                                       const double*,
                                       const int&,
                                       const double*,
                                       const double*,
#if (NDIM == 3)
                                       const double*,
#endif
                                       const int&,
                                       const int*,
                                       const int&,
                                       const int*,
                                       const int*,
#if (NDIM == 3)
                                       const int*,
#endif
                                       const int&,
                                       double*);

    void COPY_FROM_PATCHLEVEL_VEC_MAC_FC(const int&,
                                         const int&,
                                         const int&,
                                         const int&,
#if (NDIM == 3)
                                         const int&,
                                         const int&,
#endif
                                         const int&,
                                         const int&,
                                         double*,
                                         const int&,
                                         double*,
                                         double*,
#if (NDIM == 3)
                                         double*,
#endif
                                         const int&,
                                         const int*,
                                         const int&,
                                         const int*,
                                         const int*,
#if (NDIM == 3)
                                         const int*,
#endif
                                         const int&,
                                         const double*);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(Vec& vec,
                                                      const int u_data_idx,
                                                      const int u_dof_index_idx,
                                                      const int p_data_idx,
                                                      const int p_dof_index_idx,
                                                      SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    SAMRAIPointer<SAMRAIVariable> u_data_var;
    var_db->mapIndexToVariable(u_data_idx, u_data_var);
    SAMRAIPointer<SAMRAISideVariable<double>> u_data_sc_var = u_data_var;
    SAMRAIPointer<SAMRAIVariable> p_data_var;
    var_db->mapIndexToVariable(p_data_idx, p_data_var);
    SAMRAIPointer<SAMRAICellVariable<double>> p_data_cc_var = p_data_var;
    if (u_data_sc_var && p_data_cc_var)
    {
#if !defined(NDEBUG)
        SAMRAIPointer<SAMRAIVariable> u_dof_index_var;
        var_db->mapIndexToVariable(u_dof_index_idx, u_dof_index_var);
        SAMRAIPointer<SAMRAISideVariable<int>> u_dof_index_sc_var = u_dof_index_var;
        TBOX_ASSERT(u_dof_index_sc_var);
        SAMRAIPointer<SAMRAIVariable> p_dof_index_var;
        var_db->mapIndexToVariable(p_dof_index_idx, p_dof_index_var);
        SAMRAIPointer<SAMRAICellVariable<int>> p_dof_index_cc_var = p_dof_index_var;
        TBOX_ASSERT(p_dof_index_cc_var);
#endif
        copyToPatchLevelVec_MAC(vec, u_data_idx, u_dof_index_idx, p_data_idx, p_dof_index_idx, patch_level);
    }
    else
    {
        TBOX_ERROR("StaggeredStokesPETScVecUtilities::copyToPatchLevelVec():\n"
                   << "  unsupported data centering types for variables " << u_data_var->getName() << " and "
                   << p_data_var->getName() << "\n");
    }
    return;
} // copyToPatchLevelVec

void
StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(Vec& vec,
                                                        const int u_data_idx,
                                                        const int u_dof_index_idx,
                                                        const int p_data_idx,
                                                        const int p_dof_index_idx,
                                                        SAMRAIPointer<SAMRAIPatchLevel> patch_level,
                                                        SAMRAIPointer<SAMRAIRefineSchedule> data_synch_sched,
                                                        SAMRAIPointer<SAMRAIRefineSchedule> ghost_fill_sched)
{
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    SAMRAIPointer<SAMRAIVariable> u_data_var;
    var_db->mapIndexToVariable(u_data_idx, u_data_var);
    SAMRAIPointer<SAMRAISideVariable<double>> u_data_sc_var = u_data_var;
    SAMRAIPointer<SAMRAIVariable> p_data_var;
    var_db->mapIndexToVariable(p_data_idx, p_data_var);
    SAMRAIPointer<SAMRAICellVariable<double>> p_data_cc_var = p_data_var;
    if (u_data_sc_var && p_data_cc_var)
    {
#if !defined(NDEBUG)
        SAMRAIPointer<SAMRAIVariable> u_dof_index_var;
        var_db->mapIndexToVariable(u_dof_index_idx, u_dof_index_var);
        SAMRAIPointer<SAMRAISideVariable<int>> u_dof_index_sc_var = u_dof_index_var;
        TBOX_ASSERT(u_dof_index_sc_var);
        SAMRAIPointer<SAMRAIVariable> p_dof_index_var;
        var_db->mapIndexToVariable(p_dof_index_idx, p_dof_index_var);
        SAMRAIPointer<SAMRAICellVariable<int>> p_dof_index_cc_var = p_dof_index_var;
        TBOX_ASSERT(p_dof_index_cc_var);
#endif
        copyFromPatchLevelVec_MAC(vec, u_data_idx, u_dof_index_idx, p_data_idx, p_dof_index_idx, patch_level);
        if (data_synch_sched)
        {
            SAMRAIPointer<SAMRAIRefineClasses> data_synch_config = data_synch_sched->getEquivalenceClasses();
            SAMRAIRefineAlgorithm data_synch_alg;
            data_synch_alg.registerRefine(u_data_idx, u_data_idx, u_data_idx, nullptr, new SideSynchCopyFillPattern());
            data_synch_alg.resetSchedule(data_synch_sched);
            data_synch_sched->fillData(0.0);
            data_synch_sched->reset(data_synch_config);
        }
    }
    else
    {
        TBOX_ERROR("StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec():\n"
                   << "  unsupported data centering types for variables " << u_data_var->getName() << " and "
                   << p_data_var->getName() << "\n");
    }
    if (ghost_fill_sched)
    {
        SAMRAIPointer<SAMRAIRefineClasses> ghost_fill_config = ghost_fill_sched->getEquivalenceClasses();
        SAMRAIRefineAlgorithm ghost_fill_alg;
        ghost_fill_alg.registerRefine(u_data_idx, u_data_idx, u_data_idx, nullptr);
        ghost_fill_alg.registerRefine(p_data_idx, p_data_idx, p_data_idx, nullptr);
        ghost_fill_alg.resetSchedule(ghost_fill_sched);
        ghost_fill_sched->fillData(0.0);
        ghost_fill_sched->reset(ghost_fill_config);
    }
    return;
} // copyFromPatchLevelVec

SAMRAIPointer<SAMRAIRefineSchedule>
StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(const int u_data_idx,
                                                             const int p_data_idx,
                                                             SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    SAMRAIPointer<SAMRAIVariable> u_data_var;
    var_db->mapIndexToVariable(u_data_idx, u_data_var);
    SAMRAIPointer<SAMRAISideVariable<double>> u_data_sc_var = u_data_var;
    SAMRAIPointer<SAMRAIVariable> p_data_var;
    var_db->mapIndexToVariable(p_data_idx, p_data_var);
    SAMRAIPointer<SAMRAICellVariable<double>> p_data_cc_var = p_data_var;
    SAMRAIPointer<SAMRAIRefineSchedule> data_synch_sched;
    if (u_data_sc_var && p_data_cc_var)
    {
        SAMRAIRefineAlgorithm data_synch_alg;
        data_synch_alg.registerRefine(u_data_idx, u_data_idx, u_data_idx, nullptr, new SideSynchCopyFillPattern());
        data_synch_sched = data_synch_alg.createSchedule(patch_level);
    }
    else
    {
        TBOX_ERROR("StaggeredStokesPETScVecUtilities::constructDataSynchSchedule():\n"
                   << "  unsupported data centering types for variables " << u_data_var->getName() << " and "
                   << p_data_var->getName() << "\n");
    }
    return data_synch_sched;
} // constructDataSynchSchedule

SAMRAIPointer<SAMRAIRefineSchedule>
StaggeredStokesPETScVecUtilities::constructGhostFillSchedule(const int u_data_idx,
                                                             const int p_data_idx,
                                                             SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    SAMRAIRefineAlgorithm ghost_fill_alg;
    ghost_fill_alg.registerRefine(u_data_idx, u_data_idx, u_data_idx, nullptr);
    ghost_fill_alg.registerRefine(p_data_idx, p_data_idx, p_data_idx, nullptr);
    return ghost_fill_alg.createSchedule(patch_level);
} // constructGhostFillSchedule

void
StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(std::vector<int>& num_dofs_per_proc,
                                                                const int u_dof_index_idx,
                                                                const int p_dof_index_idx,
                                                                SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    SAMRAIPointer<SAMRAIVariable> u_dof_index_var;
    var_db->mapIndexToVariable(u_dof_index_idx, u_dof_index_var);
    SAMRAIPointer<SAMRAISideVariable<int>> u_dof_index_sc_var = u_dof_index_var;
    SAMRAIPointer<SAMRAIVariable> p_dof_index_var;
    var_db->mapIndexToVariable(p_dof_index_idx, p_dof_index_var);
    SAMRAIPointer<SAMRAICellVariable<int>> p_dof_index_cc_var = p_dof_index_var;
    if (u_dof_index_sc_var && p_dof_index_cc_var)
    {
        constructPatchLevelDOFIndices_MAC(num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, patch_level);
    }
    else
    {
        TBOX_ERROR("StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices():\n"
                   << "  unsupported data centering types for variables " << u_dof_index_var->getName() << " and "
                   << p_dof_index_var->getName() << "\n");
    }
    return;
} // constructPatchLevelDOFIndices

void
StaggeredStokesPETScVecUtilities::constructPatchLevelAO(AO& ao,
                                                        std::vector<int>& num_dofs_per_proc,
                                                        int u_dof_index_idx,
                                                        int p_dof_index_idx,
                                                        SAMRAIPointer<SAMRAIPatchLevel> patch_level,
                                                        int& u_ao_offset,
                                                        int& p_ao_offset)
{
    if (ao)
    {
        int ierr = AODestroy(&ao);
        IBTK_CHKERRQ(ierr);
    }

    // Determine the grid extents.
    const BoxArray<NDIM>& domain_boxes = patch_level->getPhysicalDomain();
#if !defined(NDEBUG)
    TBOX_ASSERT(domain_boxes.size() == 1);
#endif
    const SAMRAIIndex& domain_lower = domain_boxes[0].lower();
    const SAMRAIIndex& domain_upper = domain_boxes[0].upper();
    SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geom = patch_level->getGridGeometry();
    SAMRAIIntVector periodic_shift = grid_geom->getPeriodicShift(patch_level->getRatio());

    const SAMRAIIndex p_num_cells = domain_upper - domain_lower + 1;
    std::array<SAMRAIIndex, NDIM> u_num_cells;
    for (unsigned d = 0; d < NDIM; ++d)
    {
        SAMRAIIndex offset = 1;
        offset(d) = periodic_shift(d) ? 1 : 2;
        u_num_cells[d] = domain_upper - domain_lower + offset;
    }
    int n_u_samrai_dofs = 0;
    for (unsigned d = 0; d < NDIM; ++d) n_u_samrai_dofs += u_num_cells[d].getProduct();

    // Compute offsets.
    u_ao_offset = 0;
    p_ao_offset = u_ao_offset + n_u_samrai_dofs;

    // Compute PETSc to SAMRAI index mapping
    const int mpi_rank = IBTK_MPI::getRank();
    const int n_local = num_dofs_per_proc[mpi_rank];
    const int i_lower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int i_upper = i_lower + n_local;
    std::vector<int> petsc_idxs(n_local, -1), samrai_idxs(n_local, -1);

    for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
        const SAMRAIBox& patch_box = patch->getBox();
        SAMRAIPointer<SAMRAISideData<int>> u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        SAMRAIPointer<SAMRAICellData<int>> p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        const int depth = u_dof_index_data->getDepth();
#if !defined(NDEBUG)
        TBOX_ASSERT(depth == 1);
        TBOX_ASSERT(depth == p_dof_index_data);
#endif

        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            int data_offset = 0;
            for (unsigned side = 0; side < component_axis; ++side)
            {
                int side_offset = depth;
                for (unsigned d = 0; d < NDIM; ++d) side_offset *= u_num_cells[side](d);
                data_offset += side_offset;
            }

            for (SAMRAIBox::Iterator b(SAMRAISideGeometry::toSideBox(patch_box, component_axis)); b; b++)
            {
                const SAMRAICellIndex& i = b();
                const SAMRAISideIndex is(i, component_axis, SAMRAISideIndex::Lower);

                for (int d = 0; d < depth; ++d)
                {
                    const int u_dof_idx = (*u_dof_index_data)(is, d);

                    if (u_dof_idx < i_lower || u_dof_idx >= i_upper) continue;
                    petsc_idxs[u_dof_idx - i_lower] = u_dof_idx;
                    samrai_idxs[u_dof_idx - i_lower] = IndexUtilities::mapIndexToInteger(
                        i, domain_lower, u_num_cells[component_axis], d, data_offset + u_ao_offset, periodic_shift);
                }
            }
        }

        for (SAMRAIBox::Iterator b(patch_box); b; b++)
        {
            const SAMRAICellIndex& i = b();
            for (int d = 0; d < depth; ++d)
            {
                const int p_dof_idx = (*p_dof_index_data)(i, d);

                if (p_dof_idx < i_lower || p_dof_idx >= i_upper) continue;
                petsc_idxs[p_dof_idx - i_lower] = p_dof_idx;
                samrai_idxs[p_dof_idx - i_lower] =
                    IndexUtilities::mapIndexToInteger(i, domain_lower, p_num_cells, d, p_ao_offset, periodic_shift);
            }
        }
    }

    int counter = 0;
    std::vector<int> petsc_map(n_local, -1), samrai_map(n_local, -1);
    for (int k = 0; k < n_local; ++k)
    {
        if (petsc_idxs[k] < 0) continue;
        petsc_map[counter] = petsc_idxs[k];
        samrai_map[counter] = samrai_idxs[k];
        ++counter;
    }

    AOCreateMapping(PETSC_COMM_WORLD, counter, &samrai_map[0], &petsc_map[0], &ao);

    return;
} // constructPatchLevelAO

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
StaggeredStokesPETScVecUtilities::copyToPatchLevelVec_MAC(Vec& vec,
                                                          const int u_data_idx,
                                                          const int u_dof_index_idx,
                                                          const int p_data_idx,
                                                          const int p_dof_index_idx,
                                                          SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    int ierr;
    int first_local, last_local;
    ierr = VecGetOwnershipRange(vec, &first_local, &last_local);
    IBTK_CHKERRQ(ierr);
    PetscScalar* array;
    ierr = VecGetArray(vec, &array);
    IBTK_CHKERRQ(ierr);

    for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
        const SAMRAIBox& patch_box = patch->getBox();
        const SAMRAIIntVector& ilower = patch_box.lower();
        const SAMRAIIntVector& iupper = patch_box.upper();

        SAMRAIPointer<SAMRAISideData<double>> u_data = patch->getPatchData(u_data_idx);
        SAMRAIPointer<SAMRAISideData<int>> u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        const SAMRAIIntVector u_gcw = u_data->getGhostCellWidth();
        const SAMRAIIntVector u_dof_gcw = u_dof_index_data->getGhostCellWidth();
#if !defined(NDEBUG)
        TBOX_ASSERT(u_gcw.min() == u_gcw.max());
        TBOX_ASSERT(u_dof_gcw.min() == u_dof_gcw.max());
#endif

        SAMRAIPointer<SAMRAICellData<double>> p_data = patch->getPatchData(p_data_idx);
        SAMRAIPointer<SAMRAICellData<int>> p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        const SAMRAIIntVector p_gcw = p_data->getGhostCellWidth();
        const SAMRAIIntVector p_dof_gcw = p_dof_index_data->getGhostCellWidth();
#if !defined(NDEBUG)
        TBOX_ASSERT(p_gcw.min() == p_gcw.max());
        TBOX_ASSERT(p_dof_gcw.min() == p_dof_gcw.max());
#endif

        COPY_TO_PATCHLEVEL_VEC_MAC_FC(ilower(0),
                                      iupper(0),
                                      ilower(1),
                                      iupper(1),
#if (NDIM == 3)
                                      ilower(2),
                                      iupper(2),
#endif
                                      first_local,
                                      last_local,
                                      p_data->getPointer(),
                                      p_gcw.min(),
                                      u_data->getPointer(0),
                                      u_data->getPointer(1),
#if (NDIM == 3)
                                      u_data->getPointer(2),
#endif
                                      u_gcw.min(),
                                      p_dof_index_data->getPointer(),
                                      p_dof_gcw.min(),
                                      u_dof_index_data->getPointer(0),
                                      u_dof_index_data->getPointer(1),
#if (NDIM == 3)
                                      u_dof_index_data->getPointer(2),
#endif
                                      u_dof_gcw.min(),
                                      array);
    }
    ierr = VecRestoreArray(vec, &array);
    IBTK_CHKERRQ(ierr);

    return;
} // copyToPatchLevelVec_MAC

void
StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec_MAC(Vec& vec,
                                                            const int u_data_idx,
                                                            const int u_dof_index_idx,
                                                            const int p_data_idx,
                                                            const int p_dof_index_idx,
                                                            SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    int ierr;
    int first_local, last_local;
    ierr = VecGetOwnershipRange(vec, &first_local, &last_local);
    IBTK_CHKERRQ(ierr);
    const PetscScalar* array;
    ierr = VecGetArrayRead(vec, &array);
    IBTK_CHKERRQ(ierr);

    for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
        const SAMRAIBox& patch_box = patch->getBox();
        const SAMRAIIntVector& ilower = patch_box.lower();
        const SAMRAIIntVector& iupper = patch_box.upper();

        SAMRAIPointer<SAMRAISideData<double>> u_data = patch->getPatchData(u_data_idx);
        SAMRAIPointer<SAMRAISideData<int>> u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        const SAMRAIIntVector u_gcw = u_data->getGhostCellWidth();
        const SAMRAIIntVector u_dof_gcw = u_dof_index_data->getGhostCellWidth();
#if !defined(NDEBUG)
        TBOX_ASSERT(u_gcw.min() == u_gcw.max());
        TBOX_ASSERT(u_dof_gcw.min() == u_dof_gcw.max());
#endif

        SAMRAIPointer<SAMRAICellData<double>> p_data = patch->getPatchData(p_data_idx);
        SAMRAIPointer<SAMRAICellData<int>> p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        const SAMRAIIntVector p_gcw = p_data->getGhostCellWidth();
        const SAMRAIIntVector p_dof_gcw = p_dof_index_data->getGhostCellWidth();
#if !defined(NDEBUG)
        TBOX_ASSERT(p_gcw.min() == p_gcw.max());
        TBOX_ASSERT(p_dof_gcw.min() == p_dof_gcw.max());
#endif

        COPY_FROM_PATCHLEVEL_VEC_MAC_FC(ilower(0),
                                        iupper(0),
                                        ilower(1),
                                        iupper(1),
#if (NDIM == 3)
                                        ilower(2),
                                        iupper(2),
#endif
                                        first_local,
                                        last_local,
                                        p_data->getPointer(),
                                        p_gcw.min(),
                                        u_data->getPointer(0),
                                        u_data->getPointer(1),
#if (NDIM == 3)
                                        u_data->getPointer(2),
#endif
                                        u_gcw.min(),
                                        p_dof_index_data->getPointer(),
                                        p_dof_gcw.min(),
                                        u_dof_index_data->getPointer(0),
                                        u_dof_index_data->getPointer(1),
#if (NDIM == 3)
                                        u_dof_index_data->getPointer(2),
#endif
                                        u_dof_gcw.min(),
                                        array);
    }
    ierr = VecRestoreArrayRead(vec, &array);
    IBTK_CHKERRQ(ierr);

    return;
} // copyFromPatchLevelVec_MAC

void
StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices_MAC(std::vector<int>& num_dofs_per_proc,
                                                                    const int u_dof_index_idx,
                                                                    const int p_dof_index_idx,
                                                                    SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    // Create variables to keep track of whether a particular velocity location
    // is the "master" location.
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    SAMRAIPointer<SAMRAISideVariable<int>> patch_num_var = new SAMRAISideVariable<int>(
        "StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices_side()::"
        "patch_num_var");
    static const int patch_num_idx = var_db->registerPatchDataIndex(patch_num_var);
    patch_level->allocatePatchData(patch_num_idx);
    SAMRAIPointer<SAMRAISideVariable<bool>> u_mastr_loc_var = new SAMRAISideVariable<bool>(
        "StaggeredStokesPETScVecUtilities::"
        "constructPatchLevelDOFIndices_side()::u_"
        "mastr_loc_var");
    static const int u_mastr_loc_idx = var_db->registerPatchDataIndex(u_mastr_loc_var);
    patch_level->allocatePatchData(u_mastr_loc_idx);
    for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
        SAMRAIPointer<SAMRAISideData<int>> patch_num_data = patch->getPatchData(patch_num_idx);
        SAMRAIPointer<SAMRAISideData<bool>> u_mastr_loc_data = patch->getPatchData(u_mastr_loc_idx);
        patch_num_data->fillAll(patch_num);
        u_mastr_loc_data->fillAll(false);
    }

    // Synchronize the patch number at patch boundaries to determine which patch
    // owns a given DOF along patch boundaries.
    SAMRAIRefineAlgorithm bdry_synch_alg;
    bdry_synch_alg.registerRefine(patch_num_idx, patch_num_idx, patch_num_idx, nullptr, new SideSynchCopyFillPattern());
    bdry_synch_alg.registerRefine(
        u_dof_index_idx, u_dof_index_idx, u_dof_index_idx, nullptr, new SideSynchCopyFillPattern());
    bdry_synch_alg.createSchedule(patch_level)->fillData(0.0);

    // For a single patch in a periodic domain, the far side DOFs are not master.
    SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geom = patch_level->getGridGeometry();
    SAMRAIIntVector periodic_shift = grid_geom->getPeriodicShift(patch_level->getRatio());
    const BoxArray<NDIM>& domain_boxes = patch_level->getPhysicalDomain();
#if !defined(NDEBUG)
    TBOX_ASSERT(domain_boxes.size() == 1);
#endif
    const SAMRAIIndex& domain_upper = domain_boxes[0].upper();

    // Determine the number of local DOFs.
    int local_dof_count = 0;
    for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
        const SAMRAIBox& patch_box = patch->getBox();
        const SAMRAIIntVector patch_size = patch_box.numberCells();

        SAMRAIPointer<SAMRAISideData<int>> u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        SAMRAIPointer<SAMRAISideData<int>> patch_num_data = patch->getPatchData(patch_num_idx);
        SAMRAIPointer<SAMRAISideData<bool>> u_mastr_loc_data = patch->getPatchData(u_mastr_loc_idx);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            const int upper_domain_side_idx = domain_upper(component_axis) + 1;
            for (SAMRAIBox::Iterator b(SAMRAISideGeometry::toSideBox(patch_box, component_axis)); b; b++)
            {
                const SAMRAICellIndex& i = b();
                const SAMRAISideIndex is(i, component_axis, SAMRAISideIndex::Lower);
                bool fully_periodic_patch_in_axis =
                    periodic_shift(component_axis) && (patch_size(component_axis) == periodic_shift(component_axis));
                bool periodic_image = fully_periodic_patch_in_axis && (i(component_axis) == upper_domain_side_idx);
                if ((*patch_num_data)(is) == patch_num && !periodic_image)
                {
                    (*u_mastr_loc_data)(is) = true;
                    ++local_dof_count;
                }
            }
        }
        local_dof_count += SAMRAICellGeometry::toCellBox(patch_box).size();
    }

    // Determine the number of DOFs local to each MPI process and compute the
    // local DOF index offset.
    const int mpi_size = IBTK_MPI::getNodes();
    const int mpi_rank = IBTK_MPI::getRank();
    num_dofs_per_proc.resize(mpi_size);
    std::fill(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);
    IBTK_MPI::allGather(local_dof_count, &num_dofs_per_proc[0]);
    const int local_dof_offset = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);

    // Assign local DOF indices.
    int counter = local_dof_offset;
    for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
        const SAMRAIBox& patch_box = patch->getBox();
        SAMRAIPointer<SAMRAISideData<int>> u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        u_dof_index_data->fillAll(-1);
        SAMRAIPointer<SAMRAISideData<bool>> u_mastr_loc_data = patch->getPatchData(u_mastr_loc_idx);
        SAMRAIPointer<SAMRAICellData<int>> p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        p_dof_index_data->fillAll(-1);
        std::array<SAMRAIBox, NDIM> data_boxes;
        SAMRAIBoxList data_box_union(patch_box);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            data_boxes[component_axis] = SAMRAISideGeometry::toSideBox(patch_box, component_axis);
            data_box_union.unionBoxes(data_boxes[component_axis]);
        }
        data_box_union.simplifyBoxes();
        for (SAMRAIBoxList::Iterator bl(data_box_union); bl; bl++)
        {
            for (SAMRAIBox::Iterator b(bl()); b; b++)
            {
                const SAMRAIIndex& ic = b();
                for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
                {
                    if (data_boxes[component_axis].contains(ic))
                    {
                        const SAMRAISideIndex is(ic, component_axis, SAMRAISideIndex::Lower);
                        if ((*u_mastr_loc_data)(is))
                        {
                            (*u_dof_index_data)(is) = counter++;
                        }
                    }
                }
            }
        }
        for (SAMRAIBox::Iterator b(SAMRAICellGeometry::toCellBox(patch_box)); b; b++)
        {
            (*p_dof_index_data)(b()) = counter++;
        }
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(counter - local_dof_offset == num_dofs_per_proc[mpi_rank]);
#endif

    // Deallocate patch_num variable data.
    patch_level->deallocatePatchData(patch_num_idx);
    patch_level->deallocatePatchData(u_mastr_loc_idx);

    // Communicate ghost DOF indices.
    SAMRAIRefineAlgorithm dof_synch_alg;
    dof_synch_alg.registerRefine(
        u_dof_index_idx, u_dof_index_idx, u_dof_index_idx, nullptr, new SideSynchCopyFillPattern());
    dof_synch_alg.createSchedule(patch_level)->fillData(0.0);
    SAMRAIRefineAlgorithm ghost_fill_alg;
    ghost_fill_alg.registerRefine(u_dof_index_idx, u_dof_index_idx, u_dof_index_idx, nullptr);
    ghost_fill_alg.registerRefine(p_dof_index_idx, p_dof_index_idx, p_dof_index_idx, nullptr);
    ghost_fill_alg.createSchedule(patch_level)->fillData(0.0);
    return;
} // constructPatchLevelDOFIndices_MAC

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
