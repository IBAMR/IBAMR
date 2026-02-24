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

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/PETScVecUtilities.h"
#include "ibtk/SideSynchCopyFillPattern.h"
#include "ibtk/compiler_hints.h"
#include "ibtk/samrai_compatibility_names.h"

#include "BoxArray.h"
#include "MultiblockDataTranslator.h"
#include "SAMRAIBox.h"
#include "SAMRAIBoxList.h"
#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAICellData.h"
#include "SAMRAICellGeometry.h"
#include "SAMRAICellIndex.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIIndex.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIRefineAlgorithm.h"
#include "SAMRAIRefineClasses.h"
#include "SAMRAIRefineOperator.h"
#include "SAMRAIRefineSchedule.h"
#include "SAMRAISideData.h"
#include "SAMRAISideGeometry.h"
#include "SAMRAISideIndex.h"
#include "SAMRAISideVariable.h"
#include "SAMRAIUtilities.h"
#include "SAMRAIVariable.h"
#include "SAMRAIVariableDatabase.h"
#include "SAMRAIVariableFillPattern.h"

#include "petscao.h"
#include "petscvec.h"
#include <petsclog.h>

#include <algorithm>
#include <array>
#include <numeric>
#include <ostream>
#include <string>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
PETScVecUtilities::copyToPatchLevelVec(Vec& vec,
                                       const int data_idx,
                                       const int dof_index_idx,
                                       SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    SAMRAIPointer<SAMRAIVariable> data_var;
    var_db->mapIndexToVariable(data_idx, data_var);
    SAMRAIPointer<SAMRAICellVariable<double> > data_cc_var = data_var;
    SAMRAIPointer<SAMRAISideVariable<double> > data_sc_var = data_var;
    if (data_cc_var)
    {
#if !defined(NDEBUG)
        SAMRAIPointer<SAMRAIVariable> dof_index_var;
        var_db->mapIndexToVariable(dof_index_idx, dof_index_var);
        SAMRAIPointer<SAMRAICellVariable<int> > dof_index_cc_var = dof_index_var;
        TBOX_ASSERT(dof_index_cc_var);
#endif
        copyToPatchLevelVec_cell(vec, data_idx, dof_index_idx, patch_level);
    }
    else if (data_sc_var)
    {
#if !defined(NDEBUG)
        SAMRAIPointer<SAMRAIVariable> dof_index_var;
        var_db->mapIndexToVariable(dof_index_idx, dof_index_var);
        SAMRAIPointer<SAMRAISideVariable<int> > dof_index_sc_var = dof_index_var;
        TBOX_ASSERT(dof_index_sc_var);
#endif
        copyToPatchLevelVec_side(vec, data_idx, dof_index_idx, patch_level);
    }
    else
    {
        TBOX_ERROR("PETScVecUtilities::copyToPatchLevelVec():\n"
                   << "  unsupported data centering type for variable " << data_var->getName() << "\n");
    }
    return;
} // copyToPatchLevelVec

void
PETScVecUtilities::copyFromPatchLevelVec(Vec& vec,
                                         const int data_idx,
                                         const int dof_index_idx,
                                         SAMRAIPointer<SAMRAIPatchLevel> patch_level,
                                         SAMRAIPointer<SAMRAIRefineSchedule> data_synch_sched,
                                         SAMRAIPointer<SAMRAIRefineSchedule> ghost_fill_sched)
{
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    SAMRAIPointer<SAMRAIVariable> data_var;
    var_db->mapIndexToVariable(data_idx, data_var);
    SAMRAIPointer<SAMRAICellVariable<double> > data_cc_var = data_var;
    SAMRAIPointer<SAMRAISideVariable<double> > data_sc_var = data_var;
    if (data_cc_var)
    {
#if !defined(NDEBUG)
        SAMRAIPointer<SAMRAIVariable> dof_index_var;
        var_db->mapIndexToVariable(dof_index_idx, dof_index_var);
        SAMRAIPointer<SAMRAICellVariable<int> > dof_index_cc_var = dof_index_var;
        TBOX_ASSERT(dof_index_cc_var);
#endif
        copyFromPatchLevelVec_cell(vec, data_idx, dof_index_idx, patch_level);
    }
    else if (data_sc_var)
    {
#if !defined(NDEBUG)
        SAMRAIPointer<SAMRAIVariable> dof_index_var;
        var_db->mapIndexToVariable(dof_index_idx, dof_index_var);
        SAMRAIPointer<SAMRAISideVariable<int> > dof_index_sc_var = dof_index_var;
        TBOX_ASSERT(dof_index_sc_var);
#endif
        copyFromPatchLevelVec_side(vec, data_idx, dof_index_idx, patch_level);
        if (data_synch_sched)
        {
            SAMRAIPointer<SAMRAIRefineClasses> data_synch_config = data_synch_sched->getEquivalenceClasses();
            SAMRAIRefineAlgorithm data_synch_alg;
            data_synch_alg.registerRefine(data_idx, data_idx, data_idx, nullptr, new SideSynchCopyFillPattern());
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
        SAMRAIPointer<SAMRAIRefineClasses> ghost_fill_config = ghost_fill_sched->getEquivalenceClasses();
        SAMRAIRefineAlgorithm ghost_fill_alg;
        ghost_fill_alg.registerRefine(data_idx, data_idx, data_idx, nullptr);
        ghost_fill_alg.resetSchedule(ghost_fill_sched);
        ghost_fill_sched->fillData(0.0);
        ghost_fill_sched->reset(ghost_fill_config);
    }
    return;
} // copyFromPatchLevelVec

SAMRAIPointer<SAMRAIRefineSchedule>
PETScVecUtilities::constructDataSynchSchedule(const int data_idx, SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    SAMRAIPointer<SAMRAIVariable> data_var;
    var_db->mapIndexToVariable(data_idx, data_var);
    SAMRAIPointer<SAMRAICellVariable<double> > data_cc_var = data_var;
    SAMRAIPointer<SAMRAISideVariable<double> > data_sc_var = data_var;
    SAMRAIPointer<SAMRAIRefineSchedule> data_synch_sched;
    if (data_cc_var)
    {
        // intentionally blank
        //
        // NOTE: This is the only standard SAMRAI data centering that does not
        // require synchronization.
    }
    else if (data_sc_var)
    {
        SAMRAIRefineAlgorithm data_synch_alg;
        data_synch_alg.registerRefine(data_idx, data_idx, data_idx, nullptr, new SideSynchCopyFillPattern());
        data_synch_sched = data_synch_alg.createSchedule(patch_level);
    }
    else
    {
        TBOX_ERROR("PETScVecUtilities::constructDataSynchSchedule():\n"
                   << "  unsupported data centering type for variable " << data_var->getName() << "\n");
    }
    return data_synch_sched;
} // constructDataSynchSchedule

SAMRAIPointer<SAMRAIRefineSchedule>
PETScVecUtilities::constructGhostFillSchedule(const int data_idx, SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    SAMRAIRefineAlgorithm ghost_fill_alg;
    ghost_fill_alg.registerRefine(data_idx, data_idx, data_idx, nullptr);
    return ghost_fill_alg.createSchedule(patch_level);
} // constructGhostFillSchedule

void
PETScVecUtilities::constructPatchLevelDOFIndices(std::vector<int>& num_dofs_per_proc,
                                                 const int dof_index_idx,
                                                 SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    SAMRAIPointer<SAMRAIVariable> dof_index_var;
    var_db->mapIndexToVariable(dof_index_idx, dof_index_var);
    SAMRAIPointer<SAMRAICellVariable<int> > dof_index_cc_var = dof_index_var;
    SAMRAIPointer<SAMRAISideVariable<int> > dof_index_sc_var = dof_index_var;
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

void
PETScVecUtilities::constructPatchLevelAO(AO& ao,
                                         std::vector<int>& num_dofs_per_proc,
                                         int dof_index_idx,
                                         SAMRAIPointer<SAMRAIPatchLevel> patch_level,
                                         const int ao_offset)
{
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    SAMRAIPointer<SAMRAIVariable> dof_index_var;
    var_db->mapIndexToVariable(dof_index_idx, dof_index_var);
    SAMRAIPointer<SAMRAICellVariable<int> > dof_index_cc_var = dof_index_var;
    SAMRAIPointer<SAMRAISideVariable<int> > dof_index_sc_var = dof_index_var;
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
                   << "  unsupported data centering type for variable " << dof_index_var->getName() << "\n");
    }

} // constructPatchLevelAO

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
PETScVecUtilities::copyToPatchLevelVec_cell(Vec& vec,
                                            const int data_idx,
                                            const int dof_index_idx,
                                            SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    int ierr;
    int i_lower, i_upper;
    ierr = VecGetOwnershipRange(vec, &i_lower, &i_upper);
    IBTK_CHKERRQ(ierr);
    for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
        const SAMRAIBox& patch_box = patch->getBox();
        SAMRAIPointer<SAMRAICellData<double> > data = patch->getPatchData(data_idx);
        const int depth = data->getDepth();
        SAMRAIPointer<SAMRAICellData<int> > dof_index_data = patch->getPatchData(dof_index_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(depth == dof_index_data->getDepth());
#endif
        for (SAMRAIBox::Iterator b(SAMRAICellGeometry::toCellBox(patch_box)); b; b++)
        {
            const SAMRAICellIndex& i = b();
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
                                            SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    int ierr;
    int i_lower, i_upper;
    ierr = VecGetOwnershipRange(vec, &i_lower, &i_upper);
    IBTK_CHKERRQ(ierr);
    for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
        const SAMRAIBox& patch_box = patch->getBox();
        SAMRAIPointer<SAMRAISideData<double> > data = patch->getPatchData(data_idx);
        const int depth = data->getDepth();
        SAMRAIPointer<SAMRAISideData<int> > dof_index_data = patch->getPatchData(dof_index_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(depth == dof_index_data->getDepth());
#endif
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (SAMRAIBox::Iterator b(SAMRAISideGeometry::toSideBox(patch_box, component_axis)); b; b++)
            {
                const SAMRAISideIndex i(b(), component_axis, SAMRAISideIndex::Lower);
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
                                              SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    int ierr;
    int i_lower, i_upper;
    ierr = VecGetOwnershipRange(vec, &i_lower, &i_upper);
    IBTK_CHKERRQ(ierr);
    for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
        const SAMRAIBox& patch_box = patch->getBox();
        SAMRAIPointer<SAMRAICellData<double> > data = patch->getPatchData(data_idx);
        const int depth = data->getDepth();
        SAMRAIPointer<SAMRAICellData<int> > dof_index_data = patch->getPatchData(dof_index_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(depth == dof_index_data->getDepth());
#endif
        for (SAMRAIBox::Iterator b(SAMRAICellGeometry::toCellBox(patch_box)); b; b++)
        {
            const SAMRAICellIndex& i = b();
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
                                              SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    int ierr;
    int i_lower, i_upper;
    ierr = VecGetOwnershipRange(vec, &i_lower, &i_upper);
    IBTK_CHKERRQ(ierr);
    for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
        const SAMRAIBox& patch_box = patch->getBox();
        SAMRAIPointer<SAMRAISideData<double> > data = patch->getPatchData(data_idx);
        const int depth = data->getDepth();
        SAMRAIPointer<SAMRAISideData<int> > dof_index_data = patch->getPatchData(dof_index_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(depth == dof_index_data->getDepth());
#endif
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (SAMRAIBox::Iterator b(SAMRAISideGeometry::toSideBox(patch_box, component_axis)); b; b++)
            {
                const SAMRAISideIndex i(b(), component_axis, SAMRAISideIndex::Lower);
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
                                                      SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    // Determine the number of local DOFs.
    int local_dof_count = 0;
    for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
        const SAMRAIBox& patch_box = patch->getBox();
        SAMRAIPointer<SAMRAICellData<int> > dof_index_data = patch->getPatchData(dof_index_idx);
        const int depth = dof_index_data->getDepth();
        local_dof_count += depth * SAMRAICellGeometry::toCellBox(patch_box).size();
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
        SAMRAIPointer<SAMRAICellData<int> > dof_index_data = patch->getPatchData(dof_index_idx);
        dof_index_data->fillAll(-1);
        const int depth = dof_index_data->getDepth();
        for (SAMRAIBox::Iterator b(SAMRAICellGeometry::toCellBox(patch_box)); b; b++)
        {
            const SAMRAICellIndex& i = b();
            for (int d = 0; d < depth; ++d)
            {
                (*dof_index_data)(i, d) = counter++;
            }
        }
    }

    // Communicate ghost DOF indices.
    SAMRAIRefineAlgorithm ghost_fill_alg;
    ghost_fill_alg.registerRefine(dof_index_idx, dof_index_idx, dof_index_idx, nullptr);
    ghost_fill_alg.createSchedule(patch_level)->fillData(0.0);
    return;
} // constructPatchLevelDOFIndices_cell

void
PETScVecUtilities::constructPatchLevelDOFIndices_side(std::vector<int>& num_dofs_per_proc,
                                                      const int dof_index_idx,
                                                      SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
    // Create variables to keep track of whether a particular location is the
    // "master" location.
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    SAMRAIPointer<SAMRAISideVariable<int> > patch_num_var =
        new SAMRAISideVariable<int>("PETScVecUtilities::constructPatchLevelDOFIndices_side()::patch_num_var");
    static const int patch_num_idx = var_db->registerPatchDataIndex(patch_num_var);
    patch_level->allocatePatchData(patch_num_idx);
    SAMRAIPointer<SAMRAISideVariable<bool> > mastr_loc_var =
        new SAMRAISideVariable<bool>("PETScVecUtilities::constructPatchLevelDOFIndices_side()::mastr_loc_var");
    static const int mastr_loc_idx = var_db->registerPatchDataIndex(mastr_loc_var);
    patch_level->allocatePatchData(mastr_loc_idx);
    int counter = 0;
    for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
        const SAMRAIBox& patch_box = patch->getBox();
        SAMRAIPointer<SAMRAISideData<int> > dof_index_data = patch->getPatchData(dof_index_idx);
        const int depth = dof_index_data->getDepth();
        SAMRAIPointer<SAMRAISideData<int> > patch_num_data = patch->getPatchData(patch_num_idx);
        patch_num_data->fillAll(patch_num);
        SAMRAIPointer<SAMRAISideData<bool> > mastr_loc_data = patch->getPatchData(mastr_loc_idx);
        mastr_loc_data->fillAll(false);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (SAMRAIBox::Iterator b(SAMRAISideGeometry::toSideBox(patch_box, component_axis)); b; b++)
            {
                const SAMRAISideIndex i(b(), component_axis, SAMRAISideIndex::Lower);
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
    SAMRAIRefineAlgorithm bdry_synch_alg;
    bdry_synch_alg.registerRefine(patch_num_idx, patch_num_idx, patch_num_idx, nullptr, new SideSynchCopyFillPattern());
    bdry_synch_alg.registerRefine(dof_index_idx, dof_index_idx, dof_index_idx, nullptr, new SideSynchCopyFillPattern());
    bdry_synch_alg.createSchedule(patch_level)->fillData(0.0);

    // Determine the number of local DOFs.
    int local_dof_count = 0;
    counter = 0;
    for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
        const SAMRAIBox& patch_box = patch->getBox();
        SAMRAIPointer<SAMRAISideData<int> > dof_index_data = patch->getPatchData(dof_index_idx);
        const int depth = dof_index_data->getDepth();
        SAMRAIPointer<SAMRAISideData<int> > patch_num_data = patch->getPatchData(patch_num_idx);
        SAMRAIPointer<SAMRAISideData<bool> > mastr_loc_data = patch->getPatchData(mastr_loc_idx);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            for (SAMRAIBox::Iterator b(SAMRAISideGeometry::toSideBox(patch_box, component_axis)); b; b++)
            {
                const SAMRAISideIndex i(b(), component_axis, SAMRAISideIndex::Lower);
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
    const int mpi_size = IBTK_MPI::getNodes();
    const int mpi_rank = IBTK_MPI::getRank();
    num_dofs_per_proc.resize(mpi_size);
    std::fill(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);
    IBTK_MPI::allGather(local_dof_count, &num_dofs_per_proc[0]);
    const int local_dof_offset = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);

    // Assign local DOF indices.
    counter = local_dof_offset;
    for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
        const SAMRAIBox& patch_box = patch->getBox();
        SAMRAIPointer<SAMRAISideData<int> > dof_index_data = patch->getPatchData(dof_index_idx);
        const int depth = dof_index_data->getDepth();
        dof_index_data->fillAll(-1);
        SAMRAIPointer<SAMRAISideData<bool> > mastr_loc_data = patch->getPatchData(mastr_loc_idx);
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
                    if (UNLIKELY(!data_boxes[component_axis].contains(ic))) continue;
                    const SAMRAISideIndex is(ic, component_axis, SAMRAISideIndex::Lower);
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
    SAMRAIRefineAlgorithm dof_synch_alg;
    dof_synch_alg.registerRefine(dof_index_idx, dof_index_idx, dof_index_idx, nullptr, new SideSynchCopyFillPattern());
    dof_synch_alg.createSchedule(patch_level)->fillData(0.0);
    SAMRAIRefineAlgorithm ghost_fill_alg;
    ghost_fill_alg.registerRefine(dof_index_idx, dof_index_idx, dof_index_idx, nullptr);
    ghost_fill_alg.createSchedule(patch_level)->fillData(0.0);
    return;
} // constructPatchLevelDOFIndices_side

void
PETScVecUtilities::constructPatchLevelAO_cell(AO& ao,
                                              std::vector<int>& num_dofs_per_proc,
                                              const int dof_index_idx,
                                              SAMRAIPointer<SAMRAIPatchLevel> patch_level,
                                              const int ao_offset)
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
    SAMRAIIndex num_cells = 1;
    num_cells += domain_upper - domain_lower;

    // Compute PETSc to SAMRAI index mapping.
    // Note that num of local dofs can be greater than the local
    // mapping size, i.e., it is possible to map a sub-component of the
    // vector.
    const int mpi_rank = IBTK_MPI::getRank();
    const int n_local = num_dofs_per_proc[mpi_rank];
    const int i_lower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int i_upper = i_lower + n_local;
    std::vector<int> petsc_idxs(n_local, -1), samrai_idxs(n_local, -1);

    int counter = 0;
    for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
        const SAMRAIBox& patch_box = patch->getBox();
        SAMRAIPointer<SAMRAICellData<int> > dof_index_data = patch->getPatchData(dof_index_idx);
        const int depth = dof_index_data->getDepth();
        for (SAMRAIBox::Iterator b(SAMRAICellGeometry::toCellBox(patch_box)); b; b++)
        {
            const SAMRAICellIndex& i = b();

            for (int d = 0; d < depth; ++d, ++counter)
            {
                const int dof_idx = (*dof_index_data)(i, d);
#if !defined(NDEBUG)
                TBOX_ASSERT(dof_idx >= i_lower && dof_idx < i_upper);
#else
                NULL_USE(i_upper);
#endif
                petsc_idxs[counter] = dof_idx;
                samrai_idxs[counter] =
                    IndexUtilities::mapIndexToInteger(i, domain_lower, num_cells, d, ao_offset, periodic_shift);
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
                                              SAMRAIPointer<SAMRAIPatchLevel> patch_level,
                                              const int ao_offset)
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
    std::array<SAMRAIIndex, NDIM> num_cells;
    for (unsigned d = 0; d < NDIM; ++d)
    {
        SAMRAIIndex offset = 1;
        offset(d) = periodic_shift(d) ? 1 : 2;
        num_cells[d] = domain_upper - domain_lower + offset;
    }

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
        SAMRAIPointer<SAMRAISideData<int> > dof_index_data = patch->getPatchData(dof_index_idx);
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

            for (SAMRAIBox::Iterator b(SAMRAISideGeometry::toSideBox(patch_box, component_axis)); b; b++)
            {
                const SAMRAICellIndex& i = b();
                const SAMRAISideIndex is(i, component_axis, SAMRAISideIndex::Lower);

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
