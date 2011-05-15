// Filename: PETScVecUtilities.C
// Created on 23 Aug 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "PETScVecUtilities.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/SideSynchCopyFillPattern.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <ArrayDataBasicOps.h>
#include <RefineAlgorithm.h>

// C++ STDLIB INCLUDES
#include <numeric>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
PETScVecUtilities::constructPatchVecWrapper(
    Vec& vec,
    CellData<NDIM,double>& data)
{
    int ierr;
    if (vec != static_cast<Vec>(NULL))
    {
        ierr = VecDestroy(vec); IBTK_CHKERRQ(ierr);
    }
    const int nvals = data.getDepth()*CellGeometry<NDIM>::toCellBox(data.getGhostBox()).size();
    double* data_ptr = data.getPointer();
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, nvals, data_ptr, &vec); IBTK_CHKERRQ(ierr);
    return;
}// constructPatchVecWrapper

void
PETScVecUtilities::constructPatchVecWrappers(
    blitz::TinyVector<Vec,NDIM>& vecs,
    SideData<NDIM,double>& data)
{
    int ierr;
    for (int component_axis = 0; component_axis < NDIM; ++component_axis)
    {
        if (vecs[component_axis] != static_cast<Vec>(NULL))
        {
            ierr = VecDestroy(vecs[component_axis]); IBTK_CHKERRQ(ierr);
        }
        const int nvals = data.getDepth()*SideGeometry<NDIM>::toSideBox(data.getGhostBox(),component_axis).size();
        double* data_ptr = data.getPointer(component_axis);
        ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, nvals, data_ptr, &vecs[component_axis]); IBTK_CHKERRQ(ierr);
    }
    return;
}// constructPatchVecWrappers

int
PETScVecUtilities::constructPatchDOFIndices(
    CellData<NDIM,int>& dof_index,
    CellData<NDIM,double>& data)
{
    int counter = 0;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(dof_index.getBox() == data.getBox());
    TBOX_ASSERT(dof_index.getGhostBox().contains(data.getGhostBox()));
    TBOX_ASSERT(dof_index.getDepth() == data.getDepth());
#endif
    dof_index.fillAll(-1);
    const Box<NDIM>& data_box = CellGeometry<NDIM>::toCellBox(data.getGhostBox());
    const int depth = data.getDepth();
    for (int d = 0; d < depth; ++d)
    {
        for (Box<NDIM>::Iterator b(data_box); b; b++)
        {
            const CellIndex<NDIM>& i = b();
            dof_index(i,d) = counter++;
        }
    }
    return counter;
}// constructPatchDOFIndices

blitz::TinyVector<int,NDIM>
PETScVecUtilities::constructPatchDOFIndices(
    SideData<NDIM,int>& dof_index,
    SideData<NDIM,double>& data)
{
    blitz::TinyVector<int,NDIM> axis_counter(0);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(dof_index.getBox() == data.getBox());
    TBOX_ASSERT(dof_index.getGhostBox().contains(data.getGhostBox()));
    TBOX_ASSERT(dof_index.getDepth() == data.getDepth());
#endif
    dof_index.fillAll(-1);
    for (int component_axis = 0; component_axis < NDIM; ++component_axis)
    {
        const Box<NDIM>& data_box = SideGeometry<NDIM>::toSideBox(data.getGhostBox(), component_axis);
        const int depth = data.getDepth();
        for (int d = 0; d < depth; ++d)
        {
            for (Box<NDIM>::Iterator b(data_box); b; b++)
            {
                const SideIndex<NDIM> i(b(), component_axis, SideIndex<NDIM>::Lower);
                dof_index(i,d) = axis_counter[component_axis]++;
            }
        }
    }
    return axis_counter;
}// constructPatchDOFIndices

void
PETScVecUtilities::constructPatchLevelVec(
    Vec& vec,
    const int data_idx,
    Pointer<CellVariable<NDIM,double> > data_var,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    if (vec != static_cast<Vec>(NULL))
    {
        ierr = VecDestroy(vec); IBTK_CHKERRQ(ierr);
    }

    int nvals_local = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<CellData<NDIM,double> > data = patch->getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!data.isNull());
#endif
        nvals_local += data->getDepth()*CellGeometry<NDIM>::toCellBox(data->getGhostBox()).size();
    }
    ierr = VecCreateMPI(PETSC_COMM_WORLD, nvals_local, PETSC_DETERMINE, &vec); IBTK_CHKERRQ(ierr);
    return;
}// constructPatchLevelVec

void
PETScVecUtilities::constructPatchLevelVec(
    Vec& vec,
    const int data_idx,
    Pointer<SideVariable<NDIM,double> > data_var,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    if (vec != static_cast<Vec>(NULL))
    {
        ierr = VecDestroy(vec); IBTK_CHKERRQ(ierr);
    }

    int nvals_local = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<SideData<NDIM,double> > data = patch->getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!data.isNull());
#endif
        for (int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            nvals_local += data->getDepth()*SideGeometry<NDIM>::toSideBox(data->getGhostBox(),component_axis).size();
        }
    }
    ierr = VecCreateMPI(PETSC_COMM_WORLD, nvals_local, PETSC_DETERMINE, &vec); IBTK_CHKERRQ(ierr);
    return;
}// constructPatchLevelVec

void
PETScVecUtilities::copyToPatchLevelVec(
    Vec& vec,
    const int data_idx,
    Pointer<CellVariable<NDIM,double> > data_var,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    int local_offset;
    ierr = VecGetOwnershipRange(vec, &local_offset, PETSC_NULL); IBTK_CHKERRQ(ierr);
    int current_offset = local_offset;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<CellData<NDIM,double> > data = patch->getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!data.isNull());
#endif
        const int nvals_patch = data->getDepth()*CellGeometry<NDIM>::toCellBox(data->getGhostBox()).size();
        std::vector<int> idxs(nvals_patch);
        for (int k = 0; k < nvals_patch; ++k)
        {
            idxs[k] = k+current_offset;
        }
        current_offset += nvals_patch;
        ierr = VecSetValues(vec, nvals_patch, &idxs[0], data->getPointer(), INSERT_VALUES); IBTK_CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(vec); IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vec); IBTK_CHKERRQ(ierr);
    return;
}// copyToPatchLevelVec

void
PETScVecUtilities::copyFromPatchLevelVec(
    Vec& vec,
    const int data_idx,
    Pointer<CellVariable<NDIM,double> > data_var,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    int local_offset;
    ierr = VecGetOwnershipRange(vec, &local_offset, PETSC_NULL); IBTK_CHKERRQ(ierr);
    int current_offset = local_offset;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<CellData<NDIM,double> > data = patch->getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!data.isNull());
#endif
        const int nvals_patch = data->getDepth()*CellGeometry<NDIM>::toCellBox(data->getGhostBox()).size();
        std::vector<int> idxs(nvals_patch);
        for (int k = 0; k < nvals_patch; ++k)
        {
            idxs[k] = k+current_offset;
        }
        current_offset += nvals_patch;
        ierr = VecGetValues(vec, nvals_patch, &idxs[0], data->getPointer()); IBTK_CHKERRQ(ierr);
    }
    return;
}// copyFromPatchLevelVec

void
PETScVecUtilities::copyToPatchLevelVec(
    Vec& vec,
    const int data_idx,
    Pointer<SideVariable<NDIM,double> > data_var,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    int local_offset;
    ierr = VecGetOwnershipRange(vec, &local_offset, PETSC_NULL); IBTK_CHKERRQ(ierr);
    int current_offset = local_offset;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<SideData<NDIM,double> > data = patch->getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!data.isNull());
#endif
        for (int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            const int nvals_patch = data->getDepth()*SideGeometry<NDIM>::toSideBox(data->getGhostBox(),component_axis).size();
            std::vector<int> idxs(nvals_patch);
            for (int k = 0; k < nvals_patch; ++k)
            {
                idxs[k] = k+current_offset;
            }
            current_offset += nvals_patch;
            ierr = VecSetValues(vec, nvals_patch, &idxs[0], data->getPointer(component_axis), INSERT_VALUES); IBTK_CHKERRQ(ierr);
        }
    }
    ierr = VecAssemblyBegin(vec); IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vec); IBTK_CHKERRQ(ierr);
    return;
}// copyToPatchLevelVec

void
PETScVecUtilities::copyFromPatchLevelVec(
    Vec& vec,
    const int data_idx,
    Pointer<SideVariable<NDIM,double> > data_var,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    int local_offset;
    ierr = VecGetOwnershipRange(vec, &local_offset, PETSC_NULL); IBTK_CHKERRQ(ierr);
    int current_offset = local_offset;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<SideData<NDIM,double> > data = patch->getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!data.isNull());
#endif
        for (int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            const int nvals_patch = data->getDepth()*SideGeometry<NDIM>::toSideBox(data->getGhostBox(),component_axis).size();
            std::vector<int> idxs(nvals_patch);
            for (int k = 0; k < nvals_patch; ++k)
            {
                idxs[k] = k+current_offset;
            }
            current_offset += nvals_patch;
            ierr = VecGetValues(vec, nvals_patch, &idxs[0], data->getPointer(component_axis)); IBTK_CHKERRQ(ierr);
        }
    }
    return;
}// copyFromPatchLevelVec

void
PETScVecUtilities::constructPatchLevelDOFIndices(
    const int dof_index_idx,
    Pointer<CellVariable<NDIM,int> > dof_index_var,
    const int data_idx,
    Pointer<CellVariable<NDIM,double> > data_var,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    ArrayDataBasicOps<NDIM,int> patch_ops;

    // Initialize the DOF indices and determine the number of local DOF indices.
    int local_dof_count = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<CellData<NDIM,int> > dof_index = patch->getPatchData(dof_index_idx);
        Pointer<CellData<NDIM,double> > data = patch->getPatchData(data_idx);
        const int patch_dof_count = constructPatchDOFIndices(*dof_index, *data);
        patch_ops.addScalar(dof_index->getArrayData(),
                            dof_index->getArrayData(),
                            local_dof_count,
                            data->getArrayData().getBox());
        local_dof_count += patch_dof_count;
    }

    // Determine the DOF index offset.
    const int mpi_size = SAMRAI_MPI::getNodes();
    const int mpi_rank = SAMRAI_MPI::getRank();
    std::vector<int> num_dofs_proc(mpi_size,0);
    SAMRAI_MPI::allGather(local_dof_count, &num_dofs_proc[0]);
    const int local_dof_offset = std::accumulate(num_dofs_proc.begin(), num_dofs_proc.begin()+mpi_rank, 0);

    // Shift the DOF indices.
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<CellData<NDIM,int> > dof_index = patch->getPatchData(dof_index_idx);
        Pointer<CellData<NDIM,double> > data = patch->getPatchData(data_idx);
        patch_ops.addScalar(dof_index->getArrayData(),
                            dof_index->getArrayData(),
                            local_dof_offset,
                            data->getArrayData().getBox());
    }
    return;
}// constructPatchLevelDOFIndices

void
PETScVecUtilities::constructPatchLevelDOFIndices(
    const int dof_index_idx,
    Pointer<SideVariable<NDIM,int> > dof_index_var,
    const int data_idx,
    Pointer<SideVariable<NDIM,double> > data_var,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    ArrayDataBasicOps<NDIM,int> patch_ops;

    // Initialize the DOF indices and determine the number of local DOF indices.
    int local_dof_count = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<SideData<NDIM,int> > dof_index = patch->getPatchData(dof_index_idx);
        Pointer<SideData<NDIM,double> > data = patch->getPatchData(data_idx);
        const blitz::TinyVector<int,NDIM> patch_dof_counts = constructPatchDOFIndices(*dof_index, *data);
        for (int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            patch_ops.addScalar(dof_index->getArrayData(component_axis),
                                dof_index->getArrayData(component_axis),
                                local_dof_count,
                                data->getArrayData(component_axis).getBox());
            local_dof_count += patch_dof_counts[component_axis];
        }
    }

    // Determine the DOF index offset.
    const int mpi_size = SAMRAI_MPI::getNodes();
    const int mpi_rank = SAMRAI_MPI::getRank();
    std::vector<int> num_dofs_proc(mpi_size,0);
    SAMRAI_MPI::allGather(local_dof_count, &num_dofs_proc[0]);
    const int local_dof_offset = std::accumulate(num_dofs_proc.begin(), num_dofs_proc.begin()+mpi_rank, 0);

    // Shift the DOF indices.
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<SideData<NDIM,int> > dof_index = patch->getPatchData(dof_index_idx);
        Pointer<SideData<NDIM,double> > data = patch->getPatchData(data_idx);
        for (int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            patch_ops.addScalar(dof_index->getArrayData(component_axis),
                                dof_index->getArrayData(component_axis),
                                local_dof_offset,
                                data->getArrayData(component_axis).getBox());
        }
    }
    return;
}// constructPatchLevelDOFIndices

void
PETScVecUtilities::constrainPatchLevelVec(
    Vec& vec,
    const int dof_index_idx,
    Pointer<CellVariable<NDIM,int> > dof_index_var,
    Pointer<PatchLevel<NDIM> > patch_level,
    Pointer<RefineSchedule<NDIM> > dof_index_fill)
{
    int ierr;

    // Create a clone of the DOF index data and fill ghost cell values using the
    // original DOF index data.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int dof_master_index_idx = var_db->registerClonedPatchDataIndex(dof_index_var, dof_index_idx);
    patch_level->allocatePatchData(dof_master_index_idx);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<CellData<NDIM,int> > dof_index = patch->getPatchData(dof_index_idx);
        Pointer<CellData<NDIM,int> > dof_master_index = patch->getPatchData(dof_master_index_idx);
        dof_master_index->copy(*dof_index);
    }
    RefineAlgorithm<NDIM> ref_algorithm;
    ref_algorithm.registerRefine(dof_master_index_idx, dof_index_idx, dof_master_index_idx, NULL);
    if (dof_index_fill.isNull())
    {
        dof_index_fill = ref_algorithm.createSchedule(patch_level);
    }
    else
    {
        ref_algorithm.resetSchedule(dof_index_fill);
    }
    dof_index_fill->fillData(0.0);

    // Loop over the patches and constrain the ghost DOFs to their "true"
    // interior values.
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<CellData<NDIM,int> > dof_index = patch->getPatchData(dof_index_idx);
        Pointer<CellData<NDIM,int> > dof_master_index = patch->getPatchData(dof_master_index_idx);
        const Box<NDIM>& patch_box = CellGeometry<NDIM>::toCellBox(patch->getBox());
        const Box<NDIM>& ghost_box = CellGeometry<NDIM>::toCellBox(dof_index->getGhostBox());
        BoxList<NDIM> ghost_boxes(ghost_box);
        ghost_boxes.removeIntersections(patch_box);
        const int data_depth = dof_index->getDepth();

        // Ghost region values.
        for (int d = 0; d < data_depth; ++d)
        {
            for (BoxList<NDIM>::Iterator bl(ghost_boxes); bl; bl++)
            {
                for (Box<NDIM>::Iterator b(bl()); b; b++)
                {
                    const CellIndex<NDIM>& i = b();
                    const int idx = (*dof_index)(i,d);
                    const int mastr_idx = (*dof_master_index)(i,d);
                    if (idx != mastr_idx)
                    {
                        ierr = VecSetValue(vec, idx, 0.0, INSERT_VALUES); IBTK_CHKERRQ(ierr);
                    }
                }
            }
        }
    }

    // Free the temporary data and assemble the vector.
    ierr = VecAssemblyBegin(vec); IBTK_CHKERRQ(ierr);
    patch_level->deallocatePatchData(dof_master_index_idx);
    var_db->removePatchDataIndex(dof_master_index_idx);
    ierr = VecAssemblyEnd(vec); IBTK_CHKERRQ(ierr);
    return;
}// constrainPatchLevelVec

void
PETScVecUtilities::constrainPatchLevelVec(
    Vec& vec,
    const int dof_index_idx,
    Pointer<SideVariable<NDIM,int> > dof_index_var,
    Pointer<PatchLevel<NDIM> > patch_level,
    Pointer<RefineSchedule<NDIM> > dof_index_fill)
{
    int ierr;

    // Create a clone of the DOF index data and fill ghost cell values using the
    // original DOF index data.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int dof_master_index_idx = var_db->registerClonedPatchDataIndex(dof_index_var, dof_index_idx);
    patch_level->allocatePatchData(dof_master_index_idx);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<SideData<NDIM,int> > dof_index = patch->getPatchData(dof_index_idx);
        Pointer<SideData<NDIM,int> > dof_master_index = patch->getPatchData(dof_master_index_idx);
        dof_master_index->copy(*dof_index);
    }
    RefineAlgorithm<NDIM> ref_algorithm;
    ref_algorithm.registerRefine(dof_master_index_idx, dof_index_idx, dof_master_index_idx, NULL);
    ref_algorithm.registerRefine(dof_master_index_idx, dof_master_index_idx, dof_master_index_idx, NULL, new SideSynchCopyFillPattern());
    if (dof_index_fill.isNull())
    {
        dof_index_fill = ref_algorithm.createSchedule(patch_level);
    }
    else
    {
        ref_algorithm.resetSchedule(dof_index_fill);
    }
    dof_index_fill->fillData(0.0);

    // Loop over the patches and constrain the ghost DOFs to their "true"
    // interior values.
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<SideData<NDIM,int> > dof_index = patch->getPatchData(dof_index_idx);
        Pointer<SideData<NDIM,int> > dof_master_index = patch->getPatchData(dof_master_index_idx);
        for (int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            const Box<NDIM>& patch_box = SideGeometry<NDIM>::toSideBox(patch->getBox(),component_axis);
            const Box<NDIM>& ghost_box = SideGeometry<NDIM>::toSideBox(dof_index->getGhostBox(),component_axis);
            BoxList<NDIM> ghost_boxes(ghost_box);
            ghost_boxes.removeIntersections(patch_box);
            const int data_depth = dof_index->getDepth();

            // Ghost region values.
            for (int d = 0; d < data_depth; ++d)
            {
                for (BoxList<NDIM>::Iterator bl(ghost_boxes); bl; bl++)
                {
                    for (Box<NDIM>::Iterator b(bl()); b; b++)
                    {
                        const SideIndex<NDIM> i(b(),component_axis,SideIndex<NDIM>::Lower);
                        const int idx = (*dof_index)(i,d);
                        const int mastr_idx = (*dof_master_index)(i,d);
                        if (idx != mastr_idx)
                        {
                            ierr = VecSetValue(vec, idx, 0.0, INSERT_VALUES); IBTK_CHKERRQ(ierr);
                        }
                    }
                }
            }
        }
    }

    // Free the temporary data and assemble the vector.
    ierr = VecAssemblyBegin(vec); IBTK_CHKERRQ(ierr);
    patch_level->deallocatePatchData(dof_master_index_idx);
    var_db->removePatchDataIndex(dof_master_index_idx);
    ierr = VecAssemblyEnd(vec); IBTK_CHKERRQ(ierr);
    return;
}// constrainPatchLevelVec

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
