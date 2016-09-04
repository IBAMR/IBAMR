// Filename: PETScVecUtilities.cpp
// Created on 23 Aug 2010 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

#include <stddef.h>
#include <algorithm>
#include <numeric>
#include <ostream>
#include <string>
#include <vector>

#include "Box.h"
#include "BoxList.h"
#include "CellData.h"
#include "CellGeometry.h"
#include "CellIndex.h"
#include "CellVariable.h"
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
#include "boost/array.hpp"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/PETScVecUtilities.h"
#include "ibtk/SideSynchCopyFillPattern.h"
#include "ibtk/compiler_hints.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "ibtk/IndexUtilities.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Index;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
PETScVecUtilities::copyToPatchLevelVec(Vec& vec,
                                       const int data_idx,
                                       const int dof_index_idx,
                                       Pointer<PatchLevel<NDIM> > patch_level)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > data_var;
    var_db->mapIndexToVariable(data_idx, data_var);
    Pointer<CellVariable<NDIM, double> > data_cc_var = data_var;
    Pointer<SideVariable<NDIM, double> > data_sc_var = data_var;
    if (data_cc_var)
    {
#if !defined(NDEBUG)
        Pointer<Variable<NDIM> > dof_index_var;
        var_db->mapIndexToVariable(dof_index_idx, dof_index_var);
        Pointer<CellVariable<NDIM, int> > dof_index_cc_var = dof_index_var;
        TBOX_ASSERT(dof_index_cc_var);
#endif
        copyToPatchLevelVec_cell(vec, data_idx, dof_index_idx, patch_level);
    }
    else if (data_sc_var)
    {
#if !defined(NDEBUG)
        Pointer<Variable<NDIM> > dof_index_var;
        var_db->mapIndexToVariable(dof_index_idx, dof_index_var);
        Pointer<SideVariable<NDIM, int> > dof_index_sc_var = dof_index_var;
        TBOX_ASSERT(dof_index_sc_var);
#endif
        copyToPatchLevelVec_side(vec, data_idx, dof_index_idx, patch_level);
    }
    else
    {
        TBOX_ERROR("PETScVecUtilities::copyToPatchLevelVec():\n"
                   << "  unsupported data centering type for variable "
                   << data_var->getName()
                   << "\n");
    }
    return;
} // copyToPatchLevelVec

void
PETScVecUtilities::copyFromPatchLevelVec(Vec& vec,
                                         const int data_idx,
                                         const int dof_index_idx,
                                         Pointer<PatchLevel<NDIM> > patch_level,
                                         Pointer<RefineSchedule<NDIM> > data_synch_sched,
                                         Pointer<RefineSchedule<NDIM> > ghost_fill_sched)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > data_var;
    var_db->mapIndexToVariable(data_idx, data_var);
    Pointer<CellVariable<NDIM, double> > data_cc_var = data_var;
    Pointer<SideVariable<NDIM, double> > data_sc_var = data_var;
    if (data_cc_var)
    {
#if !defined(NDEBUG)
        Pointer<Variable<NDIM> > dof_index_var;
        var_db->mapIndexToVariable(dof_index_idx, dof_index_var);
        Pointer<CellVariable<NDIM, int> > dof_index_cc_var = dof_index_var;
        TBOX_ASSERT(dof_index_cc_var);
#endif
        copyFromPatchLevelVec_cell(vec, data_idx, dof_index_idx, patch_level);
    }
    else if (data_sc_var)
    {
#if !defined(NDEBUG)
        Pointer<Variable<NDIM> > dof_index_var;
        var_db->mapIndexToVariable(dof_index_idx, dof_index_var);
        Pointer<SideVariable<NDIM, int> > dof_index_sc_var = dof_index_var;
        TBOX_ASSERT(dof_index_sc_var);
#endif
        copyFromPatchLevelVec_side(vec, data_idx, dof_index_idx, patch_level);
        if (data_synch_sched)
        {
            Pointer<RefineClasses<NDIM> > data_synch_config = data_synch_sched->getEquivalenceClasses();
            RefineAlgorithm<NDIM> data_synch_alg;
            data_synch_alg.registerRefine(data_idx, data_idx, data_idx, NULL, new SideSynchCopyFillPattern());
            data_synch_alg.resetSchedule(data_synch_sched);
            data_synch_sched->fillData(0.0);
            data_synch_sched->reset(data_synch_config);
        }
    }
    else
    {
        TBOX_ERROR("PETScVecUtilities::copyFromPatchLevelVec():\n"
                   << "  unsupported data centering type for variable "
                   << data_var->getName()
                   << "\n");
    }
    if (ghost_fill_sched)
    {
        Pointer<RefineClasses<NDIM> > ghost_fill_config = ghost_fill_sched->getEquivalenceClasses();
        RefineAlgorithm<NDIM> ghost_fill_alg;
        ghost_fill_alg.registerRefine(data_idx, data_idx, data_idx, NULL);
        ghost_fill_alg.resetSchedule(ghost_fill_sched);
        ghost_fill_sched->fillData(0.0);
        ghost_fill_sched->reset(ghost_fill_config);
    }
    return;
} // copyFromPatchLevelVec

Pointer<RefineSchedule<NDIM> >
PETScVecUtilities::constructDataSynchSchedule(const int data_idx, Pointer<PatchLevel<NDIM> > patch_level)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > data_var;
    var_db->mapIndexToVariable(data_idx, data_var);
    Pointer<CellVariable<NDIM, double> > data_cc_var = data_var;
    Pointer<SideVariable<NDIM, double> > data_sc_var = data_var;
    Pointer<RefineSchedule<NDIM> > data_synch_sched;
    if (data_cc_var)
    {
        // intentionally blank
        //
        // NOTE: This is the only standard SAMRAI data centering that does not
        // require synchronization.
    }
    else if (data_sc_var)
    {
        RefineAlgorithm<NDIM> data_synch_alg;
        data_synch_alg.registerRefine(data_idx, data_idx, data_idx, NULL, new SideSynchCopyFillPattern());
        data_synch_sched = data_synch_alg.createSchedule(patch_level);
    }
    else
    {
        TBOX_ERROR("PETScVecUtilities::constructDataSynchSchedule():\n"
                   << "  unsupported data centering type for variable "
                   << data_var->getName()
                   << "\n");
    }
    return data_synch_sched;
} // constructDataSynchSchedule

Pointer<RefineSchedule<NDIM> >
PETScVecUtilities::constructGhostFillSchedule(const int data_idx, Pointer<PatchLevel<NDIM> > patch_level)
{
    RefineAlgorithm<NDIM> ghost_fill_alg;
    ghost_fill_alg.registerRefine(data_idx, data_idx, data_idx, NULL);
    return ghost_fill_alg.createSchedule(patch_level);
} // constructGhostFillSchedule

void
PETScVecUtilities::constructPatchLevelDOFIndices(std::vector<int>& num_dofs_per_proc,
                                                 const int dof_index_idx,
                                                 Pointer<PatchLevel<NDIM> > patch_level)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > dof_index_var;
    var_db->mapIndexToVariable(dof_index_idx, dof_index_var);
    Pointer<CellVariable<NDIM, int> > dof_index_cc_var = dof_index_var;
    Pointer<SideVariable<NDIM, int> > dof_index_sc_var = dof_index_var;
    if (dof_index_cc_var)
    {
        constructPatchLevelDOFIndices_cell(num_dofs_per_proc, dof_index_idx, patch_level);
    }
    else if (dof_index_sc_var)
    {
        constructPatchLevelDOFIndices_side(num_dofs_per_proc, dof_index_idx, patch_level);
    }
    else
    {
        TBOX_ERROR("PETScVecUtilities::constructPatchLevelDOFIndices():\n"
                   << "  unsupported data centering type for variable "
                   << dof_index_var->getName()
                   << "\n");
    }
    return;
} // constructPatchLevelDOFIndices

void
PETScVecUtilities::constructPatchLevelAO(AO& ao,
                                         std::vector<int>& num_dofs_per_proc,
                                         int dof_index_idx,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
                                         const int ao_offset)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > dof_index_var;
    var_db->mapIndexToVariable(dof_index_idx, dof_index_var);
    Pointer<CellVariable<NDIM, int> > dof_index_cc_var = dof_index_var;
    Pointer<SideVariable<NDIM, int> > dof_index_sc_var = dof_index_var;
    if (dof_index_cc_var)
    {
        constructPatchLevelAO_cell(ao, num_dofs_per_proc, dof_index_idx, patch_level, ao_offset);
    }
    else if (dof_index_sc_var)
    {
        constructPatchLevelAO_side(ao, num_dofs_per_proc, dof_index_idx, patch_level, ao_offset);
    }
    else
    {
        TBOX_ERROR("PETScVecUtilities::constructPatchLevelAO():\n"
                   << "  unsupported data centering type for variable "
                   << dof_index_var->getName()
                   << "\n");
    }

} // constructPatchLevelAO

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
PETScVecUtilities::copyToPatchLevelVec_cell(Vec& vec,
                                            const int data_idx,
                                            const int dof_index_idx,
                                            Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    int i_lower, i_upper;
    ierr = VecGetOwnershipRange(vec, &i_lower, &i_upper);
    IBTK_CHKERRQ(ierr);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, double> > data = patch->getPatchData(data_idx);
        const int depth = data->getDepth();
        Pointer<CellData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(depth == dof_index_data->getDepth());
#endif
        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& i = b();
            for (int d = 0; d < depth; ++d)
            {
                const int dof_index = (*dof_index_data)(i, d);
                if (LIKELY(i_lower <= dof_index && dof_index < i_upper))
                {
                    ierr = VecSetValues(vec, 1, &dof_index, &(*data)(i, d), INSERT_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
            }
        }
    }
    ierr = VecAssemblyBegin(vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vec);
    IBTK_CHKERRQ(ierr);
    return;
} // copyToPatchLevelVec_cell

void
PETScVecUtilities::copyToPatchLevelVec_side(Vec& vec,
                                            const int data_idx,
                                            const int dof_index_idx,
                                            Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    int i_lower, i_upper;
    ierr = VecGetOwnershipRange(vec, &i_lower, &i_upper);
    IBTK_CHKERRQ(ierr);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, double> > data = patch->getPatchData(data_idx);
        const int depth = data->getDepth();
        Pointer<SideData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(depth == dof_index_data->getDepth());
#endif
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, component_axis)); b; b++)
            {
                const SideIndex<NDIM> i(b(), component_axis, SideIndex<NDIM>::Lower);
                for (int d = 0; d < depth; ++d)
                {
                    const int dof_index = (*dof_index_data)(i, d);
                    if (LIKELY(i_lower <= dof_index && dof_index < i_upper))
                    {
                        ierr = VecSetValues(vec, 1, &dof_index, &(*data)(i, d), INSERT_VALUES);
                        IBTK_CHKERRQ(ierr);
                    }
                }
            }
        }
    }
    ierr = VecAssemblyBegin(vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vec);
    IBTK_CHKERRQ(ierr);
    return;
} // copyToPatchLevelVec_side

void
PETScVecUtilities::copyFromPatchLevelVec_cell(Vec& vec,
                                              const int data_idx,
                                              const int dof_index_idx,
                                              Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    int i_lower, i_upper;
    ierr = VecGetOwnershipRange(vec, &i_lower, &i_upper);
    IBTK_CHKERRQ(ierr);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, double> > data = patch->getPatchData(data_idx);
        const int depth = data->getDepth();
        Pointer<CellData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(depth == dof_index_data->getDepth());
#endif
        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& i = b();
            for (int d = 0; d < depth; ++d)
            {
                const int dof_index = (*dof_index_data)(i, d);
                if (LIKELY(i_lower <= dof_index && dof_index < i_upper))
                {
                    ierr = VecGetValues(vec, 1, &dof_index, &(*data)(i, d));
                    IBTK_CHKERRQ(ierr);
                }
            }
        }
    }
    return;
} // copyFromPatchLevelVec_cell

void
PETScVecUtilities::copyFromPatchLevelVec_side(Vec& vec,
                                              const int data_idx,
                                              const int dof_index_idx,
                                              Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    int i_lower, i_upper;
    ierr = VecGetOwnershipRange(vec, &i_lower, &i_upper);
    IBTK_CHKERRQ(ierr);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, double> > data = patch->getPatchData(data_idx);
        const int depth = data->getDepth();
        Pointer<SideData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(depth == dof_index_data->getDepth());
#endif
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, component_axis)); b; b++)
            {
                const SideIndex<NDIM> i(b(), component_axis, SideIndex<NDIM>::Lower);
                for (int d = 0; d < depth; ++d)
                {
                    const int dof_index = (*dof_index_data)(i, d);
                    if (LIKELY(i_lower <= dof_index && dof_index < i_upper))
                    {
                        ierr = VecGetValues(vec, 1, &dof_index, &(*data)(i, d));
                        IBTK_CHKERRQ(ierr);
                    }
                }
            }
        }
    }
    return;
} // copyFromPatchLevelVec_side

void
PETScVecUtilities::constructPatchLevelDOFIndices_cell(std::vector<int>& num_dofs_per_proc,
                                                      const int dof_index_idx,
                                                      Pointer<PatchLevel<NDIM> > patch_level)
{
    // Determine the number of local DOFs.
    int local_dof_count = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
        const int depth = dof_index_data->getDepth();
        local_dof_count += depth * CellGeometry<NDIM>::toCellBox(patch_box).size();
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
        Pointer<CellData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
        dof_index_data->fillAll(-1);
        const int depth = dof_index_data->getDepth();
        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& i = b();
            for (int d = 0; d < depth; ++d)
            {
                (*dof_index_data)(i, d) = counter++;
            }
        }
    }

    // Communicate ghost DOF indices.
    RefineAlgorithm<NDIM> ghost_fill_alg;
    ghost_fill_alg.registerRefine(dof_index_idx, dof_index_idx, dof_index_idx, NULL);
    ghost_fill_alg.createSchedule(patch_level)->fillData(0.0);
    return;
} // constructPatchLevelDOFIndices_cell

void
PETScVecUtilities::constructPatchLevelDOFIndices_side(std::vector<int>& num_dofs_per_proc,
                                                      const int dof_index_idx,
                                                      Pointer<PatchLevel<NDIM> > patch_level)
{
    // Create variables to keep track of whether a particular location is the
    // "master" location.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, int> > patch_num_var =
        new SideVariable<NDIM, int>("PETScVecUtilities::constructPatchLevelDOFIndices_side()::patch_num_var");
    static const int patch_num_idx = var_db->registerPatchDataIndex(patch_num_var);
    patch_level->allocatePatchData(patch_num_idx);
    Pointer<SideVariable<NDIM, bool> > mastr_loc_var =
        new SideVariable<NDIM, bool>("PETScVecUtilities::constructPatchLevelDOFIndices_side()::mastr_loc_var");
    static const int mastr_loc_idx = var_db->registerPatchDataIndex(mastr_loc_var);
    patch_level->allocatePatchData(mastr_loc_idx);
    int counter = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
        const int depth = dof_index_data->getDepth();
        Pointer<SideData<NDIM, int> > patch_num_data = patch->getPatchData(patch_num_idx);
        patch_num_data->fillAll(patch_num);
        Pointer<SideData<NDIM, bool> > mastr_loc_data = patch->getPatchData(mastr_loc_idx);
        mastr_loc_data->fillAll(false);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, component_axis)); b; b++)
            {
                const SideIndex<NDIM> i(b(), component_axis, SideIndex<NDIM>::Lower);
                for (int d = 0; d < depth; ++d)
                {
                    (*dof_index_data)(i, d) = counter++;
                }
            }
        }
    }

    // Synchronize the patch number and preliminary DOF index data at patch
    // boundaries to determine which patch owns a given DOF along patch
    // boundaries.
    RefineAlgorithm<NDIM> bdry_synch_alg;
    bdry_synch_alg.registerRefine(patch_num_idx, patch_num_idx, patch_num_idx, NULL, new SideSynchCopyFillPattern());
    bdry_synch_alg.registerRefine(dof_index_idx, dof_index_idx, dof_index_idx, NULL, new SideSynchCopyFillPattern());
    bdry_synch_alg.createSchedule(patch_level)->fillData(0.0);

    // Determine the number of local DOFs.
    int local_dof_count = 0;
    counter = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
        const int depth = dof_index_data->getDepth();
        Pointer<SideData<NDIM, int> > patch_num_data = patch->getPatchData(patch_num_idx);
        Pointer<SideData<NDIM, bool> > mastr_loc_data = patch->getPatchData(mastr_loc_idx);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, component_axis)); b; b++)
            {
                const SideIndex<NDIM> i(b(), component_axis, SideIndex<NDIM>::Lower);
                bool mastr_loc = (*patch_num_data)(i) == patch_num;
                for (int d = 0; d < depth; ++d)
                {
                    mastr_loc = ((*dof_index_data)(i, d) == counter++) && mastr_loc;
                }
                (*mastr_loc_data)(i) = mastr_loc;
                if (LIKELY(mastr_loc)) local_dof_count += depth;
            }
        }
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
    counter = local_dof_offset;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
        const int depth = dof_index_data->getDepth();
        dof_index_data->fillAll(-1);
        Pointer<SideData<NDIM, bool> > mastr_loc_data = patch->getPatchData(mastr_loc_idx);
        boost::array<Box<NDIM>, NDIM> data_boxes;
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
                    if (UNLIKELY(!data_boxes[component_axis].contains(ic))) continue;
                    const SideIndex<NDIM> is(ic, component_axis, SideIndex<NDIM>::Lower);
                    if (UNLIKELY(!(*mastr_loc_data)(is))) continue;
                    for (int d = 0; d < depth; ++d)
                    {
                        (*dof_index_data)(is, d) = counter++;
                    }
                }
            }
        }
    }

    // Deallocate temporary variable data.
    patch_level->deallocatePatchData(patch_num_idx);
    patch_level->deallocatePatchData(mastr_loc_idx);

    // Communicate ghost DOF indices.
    RefineAlgorithm<NDIM> dof_synch_alg;
    dof_synch_alg.registerRefine(dof_index_idx, dof_index_idx, dof_index_idx, NULL, new SideSynchCopyFillPattern());
    dof_synch_alg.createSchedule(patch_level)->fillData(0.0);
    RefineAlgorithm<NDIM> ghost_fill_alg;
    ghost_fill_alg.registerRefine(dof_index_idx, dof_index_idx, dof_index_idx, NULL);
    ghost_fill_alg.createSchedule(patch_level)->fillData(0.0);
    return;
} // constructPatchLevelDOFIndices_side

void
PETScVecUtilities::constructPatchLevelAO_cell(AO& ao,
                                              std::vector<int>& num_dofs_per_proc,
                                              const int dof_index_idx,
                                              Pointer<PatchLevel<NDIM> > patch_level,
                                              const int ao_offset)
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
    Index<NDIM> num_cells = 1;
    num_cells += domain_upper - domain_lower;

    // Compute PETSc to SAMRAI index mapping.
    // Note that num of local dofs can be greater than the local
    // mapping size, i.e., it is possible to map a sub-component of the
    // vector.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int n_local = num_dofs_per_proc[mpi_rank];
    const int i_lower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int i_upper = i_lower + n_local;
    std::vector<int> petsc_idxs(n_local, -1), samrai_idxs(n_local, -1);

    int counter = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
        const int depth = dof_index_data->getDepth();
        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& i = b();

            for (int d = 0; d < depth; ++d, ++counter)
            {
                const int dof_idx = (*dof_index_data)(i, d);
#if !defined(NDEBUG)
                TBOX_ASSERT(dof_idx >= i_lower && dof_idx < i_upper);
#endif
                petsc_idxs[counter] = dof_idx;
                samrai_idxs[counter] = IndexUtilities::mapIndexToInteger(i, domain_lower, num_cells, d, ao_offset);
            }
        }
    }

    AOCreateMapping(PETSC_COMM_WORLD, counter, &samrai_idxs[0], &petsc_idxs[0], &ao);

    return;

} // constructPatchLevelAO_cell

void
PETScVecUtilities::constructPatchLevelAO_side(AO& ao,
                                              std::vector<int>& num_dofs_per_proc,
                                              const int dof_index_idx,
                                              Pointer<PatchLevel<NDIM> > patch_level,
                                              const int ao_offset)
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
    boost::array<Index<NDIM>, NDIM> num_cells;
    for (unsigned d = 0; d < NDIM; ++d)
    {
        Index<NDIM> offset = 1;
        offset(d) = periodic_shift(d) ? 1 : 2;
        num_cells[d] = domain_upper - domain_lower + offset;
    }

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
        Pointer<SideData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
        const int depth = dof_index_data->getDepth();

        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            int data_offset = 0;
            for (unsigned side = 0; side < component_axis; ++side)
            {
                int side_offset = depth;
                for (unsigned d = 0; d < NDIM; ++d) side_offset *= num_cells[side](d);
                data_offset += side_offset;
            }

            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, component_axis)); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const SideIndex<NDIM> is(i, component_axis, SideIndex<NDIM>::Lower);

                for (int d = 0; d < depth; ++d)
                {
                    const int dof_idx = (*dof_index_data)(is, d);

                    if (dof_idx < i_lower || dof_idx >= i_upper) continue;
                    petsc_idxs[dof_idx - i_lower] = dof_idx;
                    samrai_idxs[dof_idx - i_lower] = IndexUtilities::mapIndexToInteger(
                        i, domain_lower, num_cells[component_axis], d, data_offset + ao_offset, periodic_shift);
                }
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

} // constructPatchLevelAO_side

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
