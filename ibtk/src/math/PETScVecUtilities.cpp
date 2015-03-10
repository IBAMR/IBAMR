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

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/IntVector.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineClasses.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/xfer/VariableFillPattern.h"
#include "boost/array.hpp"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/PETScVecUtilities.h"
#include "ibtk/SideSynchCopyFillPattern.h"
#include "ibtk/compiler_hints.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "petscsys.h"
#include "petscvec.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{

class Index;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void PETScVecUtilities::copyToPatchLevelVec(Vec& vec,
                                            const int data_idx,
                                            const int dof_index_idx,
                                            boost::shared_ptr<PatchLevel> patch_level)
{
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    boost::shared_ptr<Variable> data_var;
    var_db->mapIndexToVariable(data_idx, data_var);
    auto data_cc_var = boost::dynamic_pointer_cast<CellVariable<double> >(data_var);
    auto data_sc_var = boost::dynamic_pointer_cast<SideVariable<double> >(data_var);
    if (data_cc_var)
    {
        copyToPatchLevelVec_cell(vec, data_idx, dof_index_idx, patch_level);
    }
    else if (data_sc_var)
    {
        copyToPatchLevelVec_side(vec, data_idx, dof_index_idx, patch_level);
    }
    else
    {
        TBOX_ERROR("PETScVecUtilities::copyToPatchLevelVec():\n"
                   << "  unsupported data centering type for variable " << data_var->getName() << "\n");
    }
    return;
} // copyToPatchLevelVec

void PETScVecUtilities::copyFromPatchLevelVec(Vec& vec,
                                              const int data_idx,
                                              const int dof_index_idx,
                                              boost::shared_ptr<PatchLevel> patch_level,
                                              boost::shared_ptr<RefineSchedule> data_synch_sched,
                                              boost::shared_ptr<RefineSchedule> ghost_fill_sched)
{
    boost::shared_ptr<RefineOperator> no_refine_op;
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    boost::shared_ptr<Variable> data_var;
    var_db->mapIndexToVariable(data_idx, data_var);
    auto data_cc_var = boost::dynamic_pointer_cast<CellVariable<double> >(data_var);
    auto data_sc_var = boost::dynamic_pointer_cast<SideVariable<double> >(data_var);
    if (data_cc_var)
    {
        copyFromPatchLevelVec_cell(vec, data_idx, dof_index_idx, patch_level);
    }
    else if (data_sc_var)
    {
        copyFromPatchLevelVec_side(vec, data_idx, dof_index_idx, patch_level);
        if (data_synch_sched)
        {
            auto data_synch_config = data_synch_sched->getEquivalenceClasses();
            auto synch_op = boost::make_shared<SideSynchCopyFillPattern>();
            RefineAlgorithm data_synch_alg;
            data_synch_alg.registerRefine(data_idx, data_idx, data_idx, no_refine_op, synch_op);
            data_synch_alg.resetSchedule(data_synch_sched);
            data_synch_sched->fillData(0.0);
            data_synch_sched->reset(data_synch_config);
        }
    }
    else
    {
        TBOX_ERROR("PETScVecUtilities::copyFromPatchLevelVec():\n"
                   << "  unsupported data centering type for variable " << data_var->getName() << "\n");
    }
    if (ghost_fill_sched)
    {
        auto ghost_fill_config = ghost_fill_sched->getEquivalenceClasses();
        RefineAlgorithm ghost_fill_alg;
        ghost_fill_alg.registerRefine(data_idx, data_idx, data_idx, no_refine_op);
        ghost_fill_alg.resetSchedule(ghost_fill_sched);
        ghost_fill_sched->fillData(0.0);
        ghost_fill_sched->reset(ghost_fill_config);
    }
    return;
} // copyFromPatchLevelVec

boost::shared_ptr<RefineSchedule>
PETScVecUtilities::constructDataSynchSchedule(const int data_idx, boost::shared_ptr<PatchLevel> patch_level)
{
    boost::shared_ptr<RefineOperator> no_refine_op;
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    boost::shared_ptr<Variable> data_var;
    var_db->mapIndexToVariable(data_idx, data_var);
    auto data_cc_var = boost::dynamic_pointer_cast<CellVariable<double> >(data_var);
    auto data_sc_var = boost::dynamic_pointer_cast<SideVariable<double> >(data_var);
    boost::shared_ptr<RefineSchedule> data_synch_sched;
    if (data_cc_var)
    {
        // intentionally blank
        //
        // NOTE: This is the only standard SAMRAI data centering that does not
        // require synchronization.
    }
    else if (data_sc_var)
    {
        auto synch_op = boost::make_shared<SideSynchCopyFillPattern>();
        RefineAlgorithm data_synch_alg;
        data_synch_alg.registerRefine(data_idx, data_idx, data_idx, no_refine_op, synch_op);
        data_synch_sched = data_synch_alg.createSchedule(patch_level);
    }
    else
    {
        TBOX_ERROR("PETScVecUtilities::constructDataSynchSchedule():\n"
                   << "  unsupported data centering type for variable " << data_var->getName() << "\n");
    }
    return data_synch_sched;
} // constructDataSynchSchedule

boost::shared_ptr<RefineSchedule>
PETScVecUtilities::constructGhostFillSchedule(const int data_idx, boost::shared_ptr<PatchLevel> patch_level)
{
    boost::shared_ptr<RefineOperator> no_refine_op;
    RefineAlgorithm ghost_fill_alg;
    ghost_fill_alg.registerRefine(data_idx, data_idx, data_idx, no_refine_op);
    return ghost_fill_alg.createSchedule(patch_level);
} // constructGhostFillSchedule

void PETScVecUtilities::constructPatchLevelDOFIndices(std::vector<int>& num_dofs_per_proc,
                                                      const int dof_index_idx,
                                                      boost::shared_ptr<PatchLevel> patch_level)
{
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    boost::shared_ptr<Variable> dof_index_var;
    var_db->mapIndexToVariable(dof_index_idx, dof_index_var);
    auto dof_index_cc_var = boost::dynamic_pointer_cast<CellVariable<int> >(dof_index_var);
    auto dof_index_sc_var = boost::dynamic_pointer_cast<SideVariable<int> >(dof_index_var);
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
                   << "  unsupported data centering type for variable " << dof_index_var->getName() << "\n");
    }
    return;
} // constructPatchLevelDOFIndices

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void PETScVecUtilities::copyToPatchLevelVec_cell(Vec& vec,
                                                 const int data_idx,
                                                 const int dof_index_idx,
                                                 boost::shared_ptr<PatchLevel> patch_level)
{
    int ierr;
    int i_lower, i_upper;
    ierr = VecGetOwnershipRange(vec, &i_lower, &i_upper);
    IBTK_CHKERRQ(ierr);
    for (auto p = patch_level->begin(); p != patch_level->end(); ++p)
    {
        auto patch =*p;
        const Box& patch_box = patch->getBox();
        auto data = BOOST_CAST<CellData<double> >(patch->getPatchData(data_idx));
        const int depth = data->getDepth();
        auto dof_index_data = BOOST_CAST<CellData<int> >(patch->getPatchData(dof_index_idx));
        for (auto b = CellGeometry::begin(patch_box); b != CellGeometry::end(patch_box); ++b)
        {
            const CellIndex& i = *b;
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

void PETScVecUtilities::copyToPatchLevelVec_side(Vec& vec,
                                                 const int data_idx,
                                                 const int dof_index_idx,
                                                 boost::shared_ptr<PatchLevel> patch_level)
{
    int ierr;
    int i_lower, i_upper;
    ierr = VecGetOwnershipRange(vec, &i_lower, &i_upper);
    IBTK_CHKERRQ(ierr);
    for (auto p = patch_level->begin(); p != patch_level->end(); ++p)
    {
        auto patch =*p;
        const Box& patch_box = patch->getBox();
        auto data = BOOST_CAST<SideData<double> >(patch->getPatchData(data_idx));
        const int depth = data->getDepth();
        auto dof_index_data = BOOST_CAST<SideData<int> >(patch->getPatchData(dof_index_idx));
        TBOX_ASSERT(depth == dof_index_data->getDepth());
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (auto b = SideGeometry::begin(patch_box, component_axis);
                 b != SideGeometry::end(patch_box, component_axis);
                 ++b)
            {
                const SideIndex& i = *b;
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

void PETScVecUtilities::copyFromPatchLevelVec_cell(Vec& vec,
                                                   const int data_idx,
                                                   const int dof_index_idx,
                                                   boost::shared_ptr<PatchLevel> patch_level)
{
    int ierr;
    int i_lower, i_upper;
    ierr = VecGetOwnershipRange(vec, &i_lower, &i_upper);
    IBTK_CHKERRQ(ierr);
    for (auto p = patch_level->begin(); p != patch_level->end(); ++p)
    {
        auto patch =*p;
        const Box& patch_box = patch->getBox();
        auto data = BOOST_CAST<CellData<double> >(patch->getPatchData(data_idx));
        const int depth = data->getDepth();
        auto dof_index_data = BOOST_CAST<CellData<int> >(patch->getPatchData(dof_index_idx));
        TBOX_ASSERT(depth == dof_index_data->getDepth());
        for (auto b = CellGeometry::begin(patch_box); b != CellGeometry::end(patch_box); ++b)
        {
            const CellIndex& i = *b;
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

void PETScVecUtilities::copyFromPatchLevelVec_side(Vec& vec,
                                                   const int data_idx,
                                                   const int dof_index_idx,
                                                   boost::shared_ptr<PatchLevel> patch_level)
{
    int ierr;
    int i_lower, i_upper;
    ierr = VecGetOwnershipRange(vec, &i_lower, &i_upper);
    IBTK_CHKERRQ(ierr);
    for (auto p = patch_level->begin(); p != patch_level->end(); ++p)
    {
        auto patch =*p;
        const Box& patch_box = patch->getBox();
        auto data = BOOST_CAST<SideData<double> >(patch->getPatchData(data_idx));
        const int depth = data->getDepth();
        auto dof_index_data = BOOST_CAST<SideData<int> >(patch->getPatchData(dof_index_idx));
        TBOX_ASSERT(depth == dof_index_data->getDepth());
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (auto b = SideGeometry::begin(patch_box, component_axis);
                 b != SideGeometry::end(patch_box, component_axis);
                 ++b)
            {
                const SideIndex& i = *b;
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

void PETScVecUtilities::constructPatchLevelDOFIndices_cell(std::vector<int>& num_dofs_per_proc,
                                                           const int dof_index_idx,
                                                           boost::shared_ptr<PatchLevel> patch_level)
{
    // Determine the number of local DOFs.
    int local_dof_count = 0;
    for (auto p = patch_level->begin(); p != patch_level->end(); ++p)
    {
        auto patch =*p;
        const Box& patch_box = patch->getBox();
        auto dof_index_data = BOOST_CAST<CellData<int> >(patch->getPatchData(dof_index_idx));
        const int depth = dof_index_data->getDepth();
        local_dof_count += depth * CellGeometry::toCellBox(patch_box).size();
    }

    // Determine the number of DOFs local to each MPI process and compute the
    // local DOF index offset.
    tbox::SAMRAI_MPI comm(PETSC_COMM_WORLD);
    const int mpi_size = comm.getSize();
    const int mpi_rank = comm.getRank();
    num_dofs_per_proc.resize(mpi_size);
    std::fill(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);
    comm.Allgather(&local_dof_count, 1, MPI_INT, &num_dofs_per_proc[0], mpi_size, MPI_INT);
    const int local_dof_offset = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);

    // Assign local DOF indices.
    int counter = local_dof_offset;
    for (auto p = patch_level->begin(); p != patch_level->end(); ++p)
    {
        auto patch =*p;
        const Box& patch_box = patch->getBox();
        auto dof_index_data = BOOST_CAST<CellData<int> >(patch->getPatchData(dof_index_idx));
        dof_index_data->fillAll(-1);
        const int depth = dof_index_data->getDepth();
        for (CellIterator b = CellGeometry::begin(patch_box); b != CellGeometry::end(patch_box); ++b)
        {
            const CellIndex& i = *b;
            for (int d = 0; d < depth; ++d)
            {
                (*dof_index_data)(i, d) = counter++;
            }
        }
    }

    // Communicate ghost DOF indices.
    boost::shared_ptr<RefineOperator> no_refine_op;
    RefineAlgorithm ghost_fill_alg;
    ghost_fill_alg.registerRefine(dof_index_idx, dof_index_idx, dof_index_idx, no_refine_op);
    ghost_fill_alg.createSchedule(patch_level)->fillData(0.0);
    return;
} // constructPatchLevelDOFIndices_cell

void PETScVecUtilities::constructPatchLevelDOFIndices_side(std::vector<int>& num_dofs_per_proc,
                                                           const int dof_index_idx,
                                                           boost::shared_ptr<PatchLevel> patch_level)
{
    static const int ID_OWNER_RANK_DEPTH = 0;
    static const int ID_LOCAL_VALUE_DEPTH = 1;

    // Create variables to keep track of whether a particular location is the
    // "master" location.
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    auto patch_id_var = boost::make_shared<SideVariable<int> >(
        DIM, "PETScVecUtilities::constructPatchLevelDOFIndices_side()::patch_id_var", 2);
    static const int patch_id_idx = var_db->registerPatchDataIndex(patch_id_var);
    patch_level->allocatePatchData(patch_id_idx);
    auto mastr_loc_var = boost::make_shared<SideVariable<int> >(
        DIM, "PETScVecUtilities::constructPatchLevelDOFIndices_side()::mastr_loc_var");
    static const int mastr_loc_idx = var_db->registerPatchDataIndex(mastr_loc_var);
    patch_level->allocatePatchData(mastr_loc_idx);
    int counter = 0;
    for (auto p = patch_level->begin(); p != patch_level->end(); ++p)
    {
        auto patch =*p;
        const GlobalId& patch_id = patch->getGlobalId();
        const Box& patch_box = patch->getBox();
        auto dof_index_data = BOOST_CAST<SideData<int> >(patch->getPatchData(dof_index_idx));
        const int depth = dof_index_data->getDepth();
        auto patch_id_data = BOOST_CAST<SideData<int> >(patch->getPatchData(patch_id_idx));
        patch_id_data->fill(patch_id.getOwnerRank(), patch_id_data->getGhostBox(), ID_OWNER_RANK_DEPTH);
        patch_id_data->fill(patch_id.getLocalId().getValue(), patch_id_data->getGhostBox(), ID_LOCAL_VALUE_DEPTH);
        auto mastr_loc_data = BOOST_CAST<SideData<int> >(patch->getPatchData(mastr_loc_idx));
        mastr_loc_data->fillAll(false);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (auto b = SideGeometry::begin(patch_box, component_axis);
                 b != SideGeometry::end(patch_box, component_axis);
                 ++b)
            {
                const SideIndex& i = *b;
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
    boost::shared_ptr<RefineOperator> no_refine_op;
    auto synch_op = boost::make_shared<SideSynchCopyFillPattern>();
    RefineAlgorithm bdry_synch_alg;
    bdry_synch_alg.registerRefine(patch_id_idx, patch_id_idx, patch_id_idx, no_refine_op, synch_op);
    bdry_synch_alg.registerRefine(dof_index_idx, dof_index_idx, dof_index_idx, no_refine_op, synch_op);
    bdry_synch_alg.createSchedule(patch_level)->fillData(0.0);

    // Determine the number of local DOFs.
    int local_dof_count = 0;
    counter = 0;
    for (auto p = patch_level->begin(); p != patch_level->end(); ++p)
    {
        auto patch =*p;
        const GlobalId& patch_id = patch->getGlobalId();
        const Box& patch_box = patch->getBox();
        auto dof_index_data = BOOST_CAST<SideData<int> >(patch->getPatchData(dof_index_idx));
        const int depth = dof_index_data->getDepth();
        auto patch_id_data = BOOST_CAST<SideData<int> >(patch->getPatchData(patch_id_idx));
        auto mastr_loc_data = BOOST_CAST<SideData<int> >(patch->getPatchData(mastr_loc_idx));
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (auto b = SideGeometry::begin(patch_box, component_axis);
                 b != SideGeometry::end(patch_box, component_axis);
                 ++b)
            {
                const SideIndex& i = *b;
                const int global_id_owner = (*patch_id_data)(i, ID_OWNER_RANK_DEPTH);
                const int local_id_val = (*patch_id_data)(i, ID_LOCAL_VALUE_DEPTH);
                bool mastr_loc = GlobalId(LocalId(local_id_val), global_id_owner) == patch_id;
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
    tbox::SAMRAI_MPI comm(PETSC_COMM_WORLD);
    const int mpi_size = comm.getSize();
    const int mpi_rank = comm.getRank();
    num_dofs_per_proc.resize(mpi_size);
    std::fill(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);
    comm.Allgather(&local_dof_count, 1, MPI_INT, &num_dofs_per_proc[0], mpi_size, MPI_INT);
    const int local_dof_offset = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);

    // Assign local DOF indices.
    counter = local_dof_offset;
    for (auto p = patch_level->begin(); p != patch_level->end(); ++p)
    {
        auto patch =*p;
        const Box& patch_box = patch->getBox();
        auto dof_index_data = BOOST_CAST<SideData<int> >(patch->getPatchData(dof_index_idx));
        const int depth = dof_index_data->getDepth();
        dof_index_data->fillAll(-1);
        auto mastr_loc_data = BOOST_CAST<SideData<int> >(patch->getPatchData(mastr_loc_idx));
        std::vector<Box> data_boxes(NDIM, Box(DIM));
        BoxContainer data_box_union(patch_box);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            data_boxes[component_axis] = SideGeometry::toSideBox(patch_box, component_axis);
            BoxContainer new_boxes(data_boxes[component_axis]);
            new_boxes.removeIntersections(data_box_union);
            for (auto bl = new_boxes.begin(); bl != new_boxes.end(); ++bl)
            {
                data_box_union.push_back(*bl);
            }
        }
        for (auto bl = data_box_union.begin(); bl != data_box_union.end(); ++bl)
        {
            for (auto b = CellGeometry::begin(*bl); b != CellGeometry::end(*bl); ++b)
            {
                const Index& ic = *b;
                for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
                {
                    if (UNLIKELY(!data_boxes[component_axis].contains(ic))) continue;
                    const SideIndex is(ic, component_axis, SideIndex::Lower);
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
    patch_level->deallocatePatchData(patch_id_idx);
    patch_level->deallocatePatchData(mastr_loc_idx);

    // Communicate ghost DOF indices.
    RefineAlgorithm dof_synch_alg;
    dof_synch_alg.registerRefine(dof_index_idx, dof_index_idx, dof_index_idx, no_refine_op, synch_op);
    dof_synch_alg.createSchedule(patch_level)->fillData(0.0);
    RefineAlgorithm ghost_fill_alg;
    ghost_fill_alg.registerRefine(dof_index_idx, dof_index_idx, dof_index_idx, no_refine_op);
    ghost_fill_alg.createSchedule(patch_level)->fillData(0.0);
    return;
} // constructPatchLevelDOFIndices_side

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
