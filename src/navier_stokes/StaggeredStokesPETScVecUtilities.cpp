// Filename: StaggeredStokesPETScVecUtilities.cpp
// Created on 03 Apr 2012 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <algorithm>
#include <array>
#include <numeric>
#include <ostream>
#include <string>
#include <vector>

#include "Box.h"
#include "BoxList.h"
#include "CartesianGridGeometry.h"
#include "CellData.h"
#include "CellGeometry.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "IBAMR_config.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineClasses.h"
#include "RefineOperator.h"
#include "RefineSchedule.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "VariableFillPattern.h"
#include "ibamr/StaggeredStokesPETScVecUtilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/SideSynchCopyFillPattern.h"
#include "ibtk/compiler_hints.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

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

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Index;
} // namespace hier
} // namespace SAMRAI

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
                                                      Pointer<PatchLevel<NDIM> > patch_level)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > u_data_var;
    var_db->mapIndexToVariable(u_data_idx, u_data_var);
    Pointer<SideVariable<NDIM, double> > u_data_sc_var = u_data_var;
    Pointer<Variable<NDIM> > p_data_var;
    var_db->mapIndexToVariable(p_data_idx, p_data_var);
    Pointer<CellVariable<NDIM, double> > p_data_cc_var = p_data_var;
    if (u_data_sc_var && p_data_cc_var)
    {
#if !defined(NDEBUG)
        Pointer<Variable<NDIM> > u_dof_index_var;
        var_db->mapIndexToVariable(u_dof_index_idx, u_dof_index_var);
        Pointer<SideVariable<NDIM, int> > u_dof_index_sc_var = u_dof_index_var;
        TBOX_ASSERT(u_dof_index_sc_var);
        Pointer<Variable<NDIM> > p_dof_index_var;
        var_db->mapIndexToVariable(p_dof_index_idx, p_dof_index_var);
        Pointer<CellVariable<NDIM, int> > p_dof_index_cc_var = p_dof_index_var;
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
                                                        Pointer<PatchLevel<NDIM> > patch_level,
                                                        Pointer<RefineSchedule<NDIM> > data_synch_sched,
                                                        Pointer<RefineSchedule<NDIM> > ghost_fill_sched)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > u_data_var;
    var_db->mapIndexToVariable(u_data_idx, u_data_var);
    Pointer<SideVariable<NDIM, double> > u_data_sc_var = u_data_var;
    Pointer<Variable<NDIM> > p_data_var;
    var_db->mapIndexToVariable(p_data_idx, p_data_var);
    Pointer<CellVariable<NDIM, double> > p_data_cc_var = p_data_var;
    if (u_data_sc_var && p_data_cc_var)
    {
#if !defined(NDEBUG)
        Pointer<Variable<NDIM> > u_dof_index_var;
        var_db->mapIndexToVariable(u_dof_index_idx, u_dof_index_var);
        Pointer<SideVariable<NDIM, int> > u_dof_index_sc_var = u_dof_index_var;
        TBOX_ASSERT(u_dof_index_sc_var);
        Pointer<Variable<NDIM> > p_dof_index_var;
        var_db->mapIndexToVariable(p_dof_index_idx, p_dof_index_var);
        Pointer<CellVariable<NDIM, int> > p_dof_index_cc_var = p_dof_index_var;
        TBOX_ASSERT(p_dof_index_cc_var);
#endif
        copyFromPatchLevelVec_MAC(vec, u_data_idx, u_dof_index_idx, p_data_idx, p_dof_index_idx, patch_level);
        if (data_synch_sched)
        {
            Pointer<RefineClasses<NDIM> > data_synch_config = data_synch_sched->getEquivalenceClasses();
            RefineAlgorithm<NDIM> data_synch_alg;
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
        Pointer<RefineClasses<NDIM> > ghost_fill_config = ghost_fill_sched->getEquivalenceClasses();
        RefineAlgorithm<NDIM> ghost_fill_alg;
        ghost_fill_alg.registerRefine(u_data_idx, u_data_idx, u_data_idx, nullptr);
        ghost_fill_alg.registerRefine(p_data_idx, p_data_idx, p_data_idx, nullptr);
        ghost_fill_alg.resetSchedule(ghost_fill_sched);
        ghost_fill_sched->fillData(0.0);
        ghost_fill_sched->reset(ghost_fill_config);
    }
    return;
} // copyFromPatchLevelVec

Pointer<RefineSchedule<NDIM> >
StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(const int u_data_idx,
                                                             const int p_data_idx,
                                                             Pointer<PatchLevel<NDIM> > patch_level)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > u_data_var;
    var_db->mapIndexToVariable(u_data_idx, u_data_var);
    Pointer<SideVariable<NDIM, double> > u_data_sc_var = u_data_var;
    Pointer<Variable<NDIM> > p_data_var;
    var_db->mapIndexToVariable(p_data_idx, p_data_var);
    Pointer<CellVariable<NDIM, double> > p_data_cc_var = p_data_var;
    Pointer<RefineSchedule<NDIM> > data_synch_sched;
    if (u_data_sc_var && p_data_cc_var)
    {
        RefineAlgorithm<NDIM> data_synch_alg;
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

Pointer<RefineSchedule<NDIM> >
StaggeredStokesPETScVecUtilities::constructGhostFillSchedule(const int u_data_idx,
                                                             const int p_data_idx,
                                                             Pointer<PatchLevel<NDIM> > patch_level)
{
    RefineAlgorithm<NDIM> ghost_fill_alg;
    ghost_fill_alg.registerRefine(u_data_idx, u_data_idx, u_data_idx, nullptr);
    ghost_fill_alg.registerRefine(p_data_idx, p_data_idx, p_data_idx, nullptr);
    return ghost_fill_alg.createSchedule(patch_level);
} // constructGhostFillSchedule

void
StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(std::vector<int>& num_dofs_per_proc,
                                                                const int u_dof_index_idx,
                                                                const int p_dof_index_idx,
                                                                Pointer<PatchLevel<NDIM> > patch_level)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > u_dof_index_var;
    var_db->mapIndexToVariable(u_dof_index_idx, u_dof_index_var);
    Pointer<SideVariable<NDIM, int> > u_dof_index_sc_var = u_dof_index_var;
    Pointer<Variable<NDIM> > p_dof_index_var;
    var_db->mapIndexToVariable(p_dof_index_idx, p_dof_index_var);
    Pointer<CellVariable<NDIM, int> > p_dof_index_cc_var = p_dof_index_var;
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
                                                        Pointer<PatchLevel<NDIM> > patch_level,
                                                        int& u_ao_offset,
                                                        int& p_ao_offset)
{
    int ierr;
    if (ao)
    {
        ierr = AODestroy(&ao);
        IBTK_CHKERRQ(ierr);
    }

    // Determine the grid extents.
    const BoxArray<NDIM>& domain_boxes = patch_level->getPhysicalDomain();
#if !defined(NDEBUG)
    TBOX_ASSERT(domain_boxes.size() == 1);
#endif
    const Index<NDIM>& domain_lower = domain_boxes[0].lower();
    const Index<NDIM>& domain_upper = domain_boxes[0].upper();
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = patch_level->getGridGeometry();
    IntVector<NDIM> periodic_shift = grid_geom->getPeriodicShift(patch_level->getRatio());

    const Index<NDIM> p_num_cells = domain_upper - domain_lower + 1;
    std::array<Index<NDIM>, NDIM> u_num_cells;
    for (unsigned d = 0; d < NDIM; ++d)
    {
        Index<NDIM> offset = 1;
        offset(d) = periodic_shift(d) ? 1 : 2;
        u_num_cells[d] = domain_upper - domain_lower + offset;
    }
    int n_u_samrai_dofs = 0;
    for (unsigned d = 0; d < NDIM; ++d) n_u_samrai_dofs += u_num_cells[d].getProduct();

    // Compute offsets.
    u_ao_offset = 0;
    p_ao_offset = u_ao_offset + n_u_samrai_dofs;

    // Compute PETSc to SAMRAI index mapping
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int n_local = num_dofs_per_proc[mpi_rank];
    const int i_lower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int i_upper = i_lower + n_local;
    std::vector<int> petsc_idxs(n_local, -1), samrai_idxs(n_local, -1);

    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int> > p_dof_index_data = patch->getPatchData(p_dof_index_idx);
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

            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, component_axis)); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const SideIndex<NDIM> is(i, component_axis, SideIndex<NDIM>::Lower);

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

        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            const CellIndex<NDIM>& i = b();
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
                                                          Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    int first_local, last_local;
    ierr = VecGetOwnershipRange(vec, &first_local, &last_local);
    IBTK_CHKERRQ(ierr);
    PetscScalar* array;
    ierr = VecGetArray(vec, &array);
    IBTK_CHKERRQ(ierr);

    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const IntVector<NDIM>& ilower = patch_box.lower();
        const IntVector<NDIM>& iupper = patch_box.upper();

        Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_data_idx);
        Pointer<SideData<NDIM, int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        const IntVector<NDIM> u_gcw = u_data->getGhostCellWidth();
        const IntVector<NDIM> u_dof_gcw = u_dof_index_data->getGhostCellWidth();
#if !defined(NDEBUG)
        TBOX_ASSERT(u_gcw.min() == u_gcw.max());
        TBOX_ASSERT(u_dof_gcw.min() == u_dof_gcw.max());
#endif

        Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(p_data_idx);
        Pointer<CellData<NDIM, int> > p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        const IntVector<NDIM> p_gcw = p_data->getGhostCellWidth();
        const IntVector<NDIM> p_dof_gcw = p_dof_index_data->getGhostCellWidth();
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
                                                            Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    int first_local, last_local;
    ierr = VecGetOwnershipRange(vec, &first_local, &last_local);
    IBTK_CHKERRQ(ierr);
    const PetscScalar* array;
    ierr = VecGetArrayRead(vec, &array);
    IBTK_CHKERRQ(ierr);

    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const IntVector<NDIM>& ilower = patch_box.lower();
        const IntVector<NDIM>& iupper = patch_box.upper();

        Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_data_idx);
        Pointer<SideData<NDIM, int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        const IntVector<NDIM> u_gcw = u_data->getGhostCellWidth();
        const IntVector<NDIM> u_dof_gcw = u_dof_index_data->getGhostCellWidth();
#if !defined(NDEBUG)
        TBOX_ASSERT(u_gcw.min() == u_gcw.max());
        TBOX_ASSERT(u_dof_gcw.min() == u_dof_gcw.max());
#endif

        Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(p_data_idx);
        Pointer<CellData<NDIM, int> > p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        const IntVector<NDIM> p_gcw = p_data->getGhostCellWidth();
        const IntVector<NDIM> p_dof_gcw = p_dof_index_data->getGhostCellWidth();
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
                                                                    Pointer<PatchLevel<NDIM> > patch_level)
{
    // Create variables to keep track of whether a particular velocity location
    // is the "master" location.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, int> > patch_num_var = new SideVariable<NDIM, int>(
        "StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices_side()::patch_num_var");
    static const int patch_num_idx = var_db->registerPatchDataIndex(patch_num_var);
    patch_level->allocatePatchData(patch_num_idx);
    Pointer<SideVariable<NDIM, bool> > u_mastr_loc_var = new SideVariable<NDIM, bool>(
        "StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices_side()::u_mastr_loc_var");
    static const int u_mastr_loc_idx = var_db->registerPatchDataIndex(u_mastr_loc_var);
    patch_level->allocatePatchData(u_mastr_loc_idx);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
        Pointer<SideData<NDIM, int> > patch_num_data = patch->getPatchData(patch_num_idx);
        Pointer<SideData<NDIM, bool> > u_mastr_loc_data = patch->getPatchData(u_mastr_loc_idx);
        patch_num_data->fillAll(patch_num);
        u_mastr_loc_data->fillAll(false);
    }

    // Synchronize the patch number at patch boundaries to determine which patch
    // owns a given DOF along patch boundaries.
    RefineAlgorithm<NDIM> bdry_synch_alg;
    bdry_synch_alg.registerRefine(patch_num_idx, patch_num_idx, patch_num_idx, nullptr, new SideSynchCopyFillPattern());
    bdry_synch_alg.registerRefine(
        u_dof_index_idx, u_dof_index_idx, u_dof_index_idx, nullptr, new SideSynchCopyFillPattern());
    bdry_synch_alg.createSchedule(patch_level)->fillData(0.0);

    // For a single patch in a periodic domain, the far side DOFs are not master.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = patch_level->getGridGeometry();
    IntVector<NDIM> periodic_shift = grid_geom->getPeriodicShift(patch_level->getRatio());
    const BoxArray<NDIM>& domain_boxes = patch_level->getPhysicalDomain();
#if !defined(NDEBUG)
    TBOX_ASSERT(domain_boxes.size() == 1);
#endif
    const Index<NDIM>& domain_upper = domain_boxes[0].upper();

    // Determine the number of local DOFs.
    int local_dof_count = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
        const Box<NDIM>& patch_box = patch->getBox();
        const IntVector<NDIM> patch_size = patch_box.numberCells();

        Pointer<SideData<NDIM, int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<SideData<NDIM, int> > patch_num_data = patch->getPatchData(patch_num_idx);
        Pointer<SideData<NDIM, bool> > u_mastr_loc_data = patch->getPatchData(u_mastr_loc_idx);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            const int upper_domain_side_idx = domain_upper(component_axis) + 1;
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, component_axis)); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const SideIndex<NDIM> is(i, component_axis, SideIndex<NDIM>::Lower);
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
        local_dof_count += CellGeometry<NDIM>::toCellBox(patch_box).size();
    }

    // Determine the number of DOFs local to each MPI process and compute the
    // local DOF index offset.
    const int mpi_size = SAMRAI_MPI::getNodes();
    const int mpi_rank = SAMRAI_MPI::getRank();
    num_dofs_per_proc.resize(mpi_size);
    std::fill(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);
    SAMRAI_MPI::allGather(local_dof_count, &num_dofs_per_proc[0]);
    const int local_dof_offset = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);

    // Assign local DOF indices.
    int counter = local_dof_offset;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        u_dof_index_data->fillAll(-1);
        Pointer<SideData<NDIM, bool> > u_mastr_loc_data = patch->getPatchData(u_mastr_loc_idx);
        Pointer<CellData<NDIM, int> > p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        p_dof_index_data->fillAll(-1);
        std::array<Box<NDIM>, NDIM> data_boxes;
        BoxList<NDIM> data_box_union(patch_box);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            data_boxes[component_axis] = SideGeometry<NDIM>::toSideBox(patch_box, component_axis);
            data_box_union.unionBoxes(data_boxes[component_axis]);
        }
        data_box_union.simplifyBoxes();
        for (BoxList<NDIM>::Iterator bl(data_box_union); bl; bl++)
        {
            for (Box<NDIM>::Iterator b(bl()); b; b++)
            {
                const Index<NDIM>& ic = b();
                for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
                {
                    if (data_boxes[component_axis].contains(ic))
                    {
                        const SideIndex<NDIM> is(ic, component_axis, SideIndex<NDIM>::Lower);
                        if ((*u_mastr_loc_data)(is))
                        {
                            (*u_dof_index_data)(is) = counter++;
                        }
                    }
                }
            }
        }
        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
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
    RefineAlgorithm<NDIM> dof_synch_alg;
    dof_synch_alg.registerRefine(
        u_dof_index_idx, u_dof_index_idx, u_dof_index_idx, nullptr, new SideSynchCopyFillPattern());
    dof_synch_alg.createSchedule(patch_level)->fillData(0.0);
    RefineAlgorithm<NDIM> ghost_fill_alg;
    ghost_fill_alg.registerRefine(u_dof_index_idx, u_dof_index_idx, u_dof_index_idx, nullptr);
    ghost_fill_alg.registerRefine(p_dof_index_idx, p_dof_index_idx, p_dof_index_idx, nullptr);
    ghost_fill_alg.createSchedule(patch_level)->fillData(0.0);
    return;
} // constructPatchLevelDOFIndices_MAC

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
