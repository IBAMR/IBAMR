// Filename: StaggeredStokesPETScVecUtilities.cpp
// Created on 03 Apr 2012 by Boyce Griffith
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
#include "ibamr/StaggeredStokesPETScVecUtilities.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/SideSynchCopyFillPattern.h"
#include "ibtk/compiler_hints.h"
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

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(Vec& vec,
                                                           const int u_data_idx,
                                                           const int u_dof_index_idx,
                                                           const int p_data_idx,
                                                           const int p_dof_index_idx,
                                                           boost::shared_ptr<PatchLevel> patch_level)
{
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    boost::shared_ptr<Variable> u_data_var;
    var_db->mapIndexToVariable(u_data_idx, u_data_var);
    auto u_data_sc_var = BOOST_CAST<SideVariable<double> >(u_data_var);
    boost::shared_ptr<Variable> p_data_var;
    var_db->mapIndexToVariable(p_data_idx, p_data_var);
    auto p_data_cc_var = BOOST_CAST<CellVariable<double> >(p_data_var);
    if (u_data_sc_var && p_data_cc_var)
    {
        copyToPatchLevelVec_MAC(vec, u_data_idx, u_dof_index_idx, p_data_idx, p_dof_index_idx, patch_level);
    }
    else
    {
        TBOX_ERROR("StaggeredStokesPETScVecUtilities::copyToPatchLevelVec():\n"
                   << "  unsupported data centering types for variables " << u_data_var->getName() << " and "
                   << p_data_var->getName() << "\n");
    }
    return;
}

void StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(Vec& vec,
                                                             const int u_data_idx,
                                                             const int u_dof_index_idx,
                                                             const int p_data_idx,
                                                             const int p_dof_index_idx,
                                                             boost::shared_ptr<PatchLevel> patch_level,
                                                             boost::shared_ptr<RefineSchedule> data_synch_sched,
                                                             boost::shared_ptr<RefineSchedule> ghost_fill_sched)
{
    boost::shared_ptr<RefineOperator> no_refine_op;
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    boost::shared_ptr<Variable> u_data_var;
    var_db->mapIndexToVariable(u_data_idx, u_data_var);
    auto u_data_sc_var = BOOST_CAST<SideVariable<double> >(u_data_var);
    boost::shared_ptr<Variable> p_data_var;
    var_db->mapIndexToVariable(p_data_idx, p_data_var);
    auto p_data_cc_var = BOOST_CAST<CellVariable<double> >(p_data_var);
    if (u_data_sc_var && p_data_cc_var)
    {
        copyFromPatchLevelVec_MAC(vec, u_data_idx, u_dof_index_idx, p_data_idx, p_dof_index_idx, patch_level);
        if (data_synch_sched)
        {
            auto data_synch_config = data_synch_sched->getEquivalenceClasses();
            RefineAlgorithm data_synch_alg(DIM);
            auto synch_op = boost::make_shared<SideSynchCopyFillPattern>();
            data_synch_alg.registerRefine(u_data_idx, u_data_idx, u_data_idx, no_refine_op, synch_op);
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
        auto ghost_fill_config = ghost_fill_sched->getEquivalenceClasses();
        RefineAlgorithm ghost_fill_alg(DIM);
        ghost_fill_alg.registerRefine(u_data_idx, u_data_idx, u_data_idx, no_refine_op);
        ghost_fill_alg.registerRefine(p_data_idx, p_data_idx, p_data_idx, no_refine_op);
        ghost_fill_alg.resetSchedule(ghost_fill_sched);
        ghost_fill_sched->fillData(0.0);
        ghost_fill_sched->reset(ghost_fill_config);
    }
    return;
}

boost::shared_ptr<RefineSchedule>
StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(const int u_data_idx,
                                                             const int p_data_idx,
                                                             boost::shared_ptr<PatchLevel> patch_level)
{
    boost::shared_ptr<RefineOperator> no_refine_op;
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    boost::shared_ptr<Variable> u_data_var;
    var_db->mapIndexToVariable(u_data_idx, u_data_var);
    auto u_data_sc_var = BOOST_CAST<SideVariable<double> >(u_data_var);
    boost::shared_ptr<Variable> p_data_var;
    var_db->mapIndexToVariable(p_data_idx, p_data_var);
    auto p_data_cc_var = BOOST_CAST<CellVariable<double> >(p_data_var);
    boost::shared_ptr<RefineSchedule> data_synch_sched;
    if (u_data_sc_var && p_data_cc_var)
    {
        auto synch_op = boost::make_shared<SideSynchCopyFillPattern>();
        RefineAlgorithm data_synch_alg(DIM);
        data_synch_alg.registerRefine(u_data_idx, u_data_idx, u_data_idx, no_refine_op, synch_op);
        data_synch_sched = data_synch_alg.createSchedule(patch_level);
    }
    else
    {
        TBOX_ERROR("StaggeredStokesPETScVecUtilities::constructDataSynchSchedule():\n"
                   << "  unsupported data centering types for variables " << u_data_var->getName() << " and "
                   << p_data_var->getName() << "\n");
    }
    return data_synch_sched;
}

boost::shared_ptr<RefineSchedule>
StaggeredStokesPETScVecUtilities::constructGhostFillSchedule(const int u_data_idx,
                                                             const int p_data_idx,
                                                             boost::shared_ptr<PatchLevel> patch_level)
{
    boost::shared_ptr<RefineOperator> no_refine_op;
    RefineAlgorithm ghost_fill_alg(DIM);
    ghost_fill_alg.registerRefine(u_data_idx, u_data_idx, u_data_idx, no_refine_op);
    ghost_fill_alg.registerRefine(p_data_idx, p_data_idx, p_data_idx, no_refine_op);
    return ghost_fill_alg.createSchedule(patch_level);
}

void StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(std::vector<int>& num_dofs_per_proc,
                                                                     const int u_dof_index_idx,
                                                                     const int p_dof_index_idx,
                                                                     boost::shared_ptr<PatchLevel> patch_level)
{
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    boost::shared_ptr<Variable> u_dof_index_var;
    var_db->mapIndexToVariable(u_dof_index_idx, u_dof_index_var);
    auto u_dof_index_sc_var = BOOST_CAST<SideVariable<int> >(u_dof_index_var);
    boost::shared_ptr<Variable> p_dof_index_var;
    var_db->mapIndexToVariable(p_dof_index_idx, p_dof_index_var);
    auto p_dof_index_cc_var = BOOST_CAST<CellVariable<int> >(p_dof_index_var);
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
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void StaggeredStokesPETScVecUtilities::copyToPatchLevelVec_MAC(Vec& vec,
                                                               const int u_data_idx,
                                                               const int u_dof_index_idx,
                                                               const int p_data_idx,
                                                               const int p_dof_index_idx,
                                                               boost::shared_ptr<PatchLevel> patch_level)
{
    int ierr;
    int ilower, iupper;
    ierr = VecGetOwnershipRange(vec, &ilower, &iupper);
    IBTK_CHKERRQ(ierr);
    for (auto p = patch_level->begin(); p != patch_level->end(); ++p)
    {
        auto patch =*p;
        const Box& patch_box = patch->getBox();
        boost::shared_ptr<SideData<double> > u_data = patch->getPatchData(u_data_idx);
        boost::shared_ptr<SideData<int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (SideIterator b(patch_box, component_axis); b; b++)
            {
                const SideIndex& is = b();
                const int dof_index = (*u_dof_index_data)(is);
                if (LIKELY(ilower <= dof_index && dof_index < iupper))
                {
                    ierr = VecSetValues(vec, 1, &dof_index, &(*u_data)(is), INSERT_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
            }
        }
        boost::shared_ptr<CellData<double> > p_data = patch->getPatchData(p_data_idx);
        boost::shared_ptr<CellData<int> > p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        for (CellIterator b = CellGeometry::begin(patch_box); b != CellGeometry::end(patch_box); ++b)
        {
            const CellIndex& ic = b();
            const int dof_index = (*p_dof_index_data)(ic);
            if (LIKELY(ilower <= dof_index && dof_index < iupper))
            {
                ierr = VecSetValues(vec, 1, &dof_index, &(*p_data)(ic), INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
        }
    }
    ierr = VecAssemblyBegin(vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vec);
    IBTK_CHKERRQ(ierr);
    return;
}

void StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec_MAC(Vec& vec,
                                                                 const int u_data_idx,
                                                                 const int u_dof_index_idx,
                                                                 const int p_data_idx,
                                                                 const int p_dof_index_idx,
                                                                 boost::shared_ptr<PatchLevel> patch_level)
{
    int ierr;
    int ilower, iupper;
    ierr = VecGetOwnershipRange(vec, &ilower, &iupper);
    IBTK_CHKERRQ(ierr);
    for (auto p = patch_level->begin(); p != patch_level->end(); ++p)
    {
        auto patch =*p;
        const Box& patch_box = patch->getBox();
        boost::shared_ptr<SideData<double> > u_data = patch->getPatchData(u_data_idx);
        boost::shared_ptr<SideData<int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (SideIterator b(patch_box, component_axis); b; b++)
            {
                const SideIndex& is = b();
                const int dof_index = (*u_dof_index_data)(is);
                if (LIKELY(ilower <= dof_index && dof_index < iupper))
                {
                    ierr = VecGetValues(vec, 1, &dof_index, &(*u_data)(is));
                    IBTK_CHKERRQ(ierr);
                }
            }
        }
        boost::shared_ptr<CellData<double> > p_data = patch->getPatchData(p_data_idx);
        boost::shared_ptr<CellData<int> > p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        for (CellIterator b = CellGeometry::begin(patch_box); b != CellGeometry::end(patch_box); ++b)
        {
            const CellIndex& ic = b();
            const int dof_index = (*p_dof_index_data)(ic);
            if (LIKELY(ilower <= dof_index && dof_index < iupper))
            {
                ierr = VecGetValues(vec, 1, &dof_index, &(*p_data)(ic));
                IBTK_CHKERRQ(ierr);
            }
        }
    }
    ierr = VecAssemblyBegin(vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vec);
    IBTK_CHKERRQ(ierr);
    return;
}

void StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices_MAC(std::vector<int>& num_dofs_per_proc,
                                                                         const int u_dof_index_idx,
                                                                         const int p_dof_index_idx,
                                                                         boost::shared_ptr<PatchLevel> patch_level)
{
    static const int ID_OWNER_RANK_DEPTH = 0;
    static const int ID_LOCAL_VALUE_DEPTH = 1;

    // Create variables to keep track of whether a particular velocity location
    // is the "master" location.
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    auto patch_id_var = boost::make_shared<SideVariable<int> >(
        DIM, "StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices_side()::patch_id_var", 2);
    static const int patch_id_idx = var_db->registerPatchDataIndex(patch_id_var);
    patch_level->allocatePatchData(patch_id_idx);
    auto u_mastr_loc_var = boost::make_shared<SideVariable<bool> >(
        DIM, "StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices_side()::u_mastr_loc_var");
    static const int u_mastr_loc_idx = var_db->registerPatchDataIndex(u_mastr_loc_var);
    patch_level->allocatePatchData(u_mastr_loc_idx);
    int counter = 0;
    for (auto p = patch_level->begin(); p != patch_level->end(); ++p)
    {
        auto patch =*p;
        const GlobalId& patch_id = patch->getGlobalId();
        const Box& patch_box = patch->getBox();
        boost::shared_ptr<SideData<int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        boost::shared_ptr<SideData<int> > patch_id_data = patch->getPatchData(patch_id_idx);
        boost::shared_ptr<SideData<bool> > u_mastr_loc_data = patch->getPatchData(u_mastr_loc_idx);
        patch_id_data->fill(patch_id.getOwnerRank(), patch_id_data->getGhostBox(), ID_OWNER_RANK_DEPTH);
        patch_id_data->fill(patch_id.getLocalId().getValue(), patch_id_data->getGhostBox(), ID_LOCAL_VALUE_DEPTH);
        u_mastr_loc_data->fillAll(false);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (SideIterator b(patch_box, component_axis); b; b++)
            {
                (*u_dof_index_data)(b()) = counter++;
            }
        }
    }

    // Synchronize the patch number and preliminary DOF index data at patch
    // boundaries to determine which patch owns a given DOF along patch
    // boundaries.
    boost::shared_ptr<RefineOperator> no_refine_op;
    auto synch_op = boost::make_shared<SideSynchCopyFillPattern>();
    RefineAlgorithm bdry_synch_alg(DIM);
    bdry_synch_alg.registerRefine(patch_id_idx, patch_id_idx, patch_id_idx, no_refine_op, synch_op);
    bdry_synch_alg.registerRefine(u_dof_index_idx, u_dof_index_idx, u_dof_index_idx, no_refine_op, synch_op);
    bdry_synch_alg.createSchedule(patch_level)->fillData(0.0);

    // Determine the number of local DOFs.
    int local_dof_count = 0;
    counter = 0;
    for (auto p = patch_level->begin(); p != patch_level->end(); ++p)
    {
        auto patch =*p;
        const GlobalId& patch_id = patch->getGlobalId();
        const Box& patch_box = patch->getBox();
        boost::shared_ptr<SideData<int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        boost::shared_ptr<SideData<int> > patch_id_data = patch->getPatchData(patch_id_idx);
        boost::shared_ptr<SideData<bool> > u_mastr_loc_data = patch->getPatchData(u_mastr_loc_idx);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (SideIterator b(patch_box, component_axis); b; b++)
            {
                const SideIndex& is = b();
                const int global_id_owner = (*patch_id_data)(is, ID_OWNER_RANK_DEPTH);
                const int local_id_val = (*patch_id_data)(is, ID_LOCAL_VALUE_DEPTH);
                const bool u_mastr_loc = ((*u_dof_index_data)(is) == counter++) &&
                                         (GlobalId(LocalId(local_id_val), global_id_owner) == patch_id);
                (*u_mastr_loc_data)(is) = u_mastr_loc;
                if (LIKELY(u_mastr_loc)) ++local_dof_count;
            }
        }
        local_dof_count += CellGeometry::toCellBox(patch_box).size();
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
        boost::shared_ptr<SideData<int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        u_dof_index_data->fillAll(-1);
        boost::shared_ptr<SideData<bool> > u_mastr_loc_data = patch->getPatchData(u_mastr_loc_idx);
        boost::shared_ptr<CellData<int> > p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        p_dof_index_data->fillAll(-1);
        std::vector<Box> data_boxes(NDIM, Box(DIM));
        BoxContainer data_box_union(patch_box);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            data_boxes[component_axis] = SideGeometry::toSideBox(patch_box, component_axis);
            data_box_union.unionBoxes(data_boxes[component_axis]);
        }
        data_box_union.simplifyBoxes();
        for (auto bl(data_box_union); bl; bl++)
        {
            for (CellIterator b(bl()); b; b++)
            {
                const CellIndex& ic = b();
                for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
                {
                    if (UNLIKELY(!data_boxes[component_axis].contains(ic))) continue;
                    const SideIndex is(ic, component_axis, SideIndex::Lower);
                    if (UNLIKELY(!(*u_mastr_loc_data)(is))) continue;
                    (*u_dof_index_data)(is) = counter++;
                }
            }
        }
        for (CellIterator b = CellGeometry::begin(patch_box); b != CellGeometry::end(patch_box); ++b)
        {
            (*p_dof_index_data)(b()) = counter++;
        }
    }

    // Deallocate patch_id variable data.
    patch_level->deallocatePatchData(patch_id_idx);
    patch_level->deallocatePatchData(u_mastr_loc_idx);

    // Communicate ghost DOF indices.
    RefineAlgorithm dof_synch_alg(DIM);
    dof_synch_alg.registerRefine(u_dof_index_idx, u_dof_index_idx, u_dof_index_idx, no_refine_op, synch_op);
    dof_synch_alg.createSchedule(patch_level)->fillData(0.0);
    RefineAlgorithm ghost_fill_alg(DIM);
    ghost_fill_alg.registerRefine(u_dof_index_idx, u_dof_index_idx, u_dof_index_idx, no_refine_op);
    ghost_fill_alg.registerRefine(p_dof_index_idx, p_dof_index_idx, p_dof_index_idx, no_refine_op);
    ghost_fill_alg.createSchedule(patch_level)->fillData(0.0);
    return;
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
