// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
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

#include "ibtk/IBTK_MPI.h"
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/SAMRAIGhostDataAccumulator.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <petscis.h>
#include <petscistypes.h>
#include <petsclog.h>
#include <petscvec.h>

#include <ArrayData.h>
#include <BasePatchHierarchy.h>
#include <Box.h>
#include <CellData.h>
#include <CellDataFactory.h>
#include <CellVariable.h>
#include <IntVector.h>
#include <MultiblockDataTranslator.h>
#include <Patch.h>
#include <PatchData.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <SideDataFactory.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <numeric>
#include <string>
#include <type_traits>
#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////
namespace
{
static Timer *t_constructor;
static Timer *t_accumulate_ghost_data;
static Timer *t_add_or_get;

// The code leading to VecGetValues and VecSetValues is identical, so use one
// function for both. Here the last argument is true if we set and false if we
// get.
void
add_or_get(Pointer<PatchLevel<NDIM> >& level,
           const IntVector<NDIM> gcw,
           const bool cc_data,
           const int local_dof_idx,
           const int value_idx,
           Vec vec,
           const bool add)
{
    IBTK_TIMER_START(t_add_or_get);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        TBOX_ASSERT(gcw == patch->getPatchData(value_idx)->getGhostCellWidth());
        std::vector<ArrayData<NDIM, double>*> values_ptrs;
        std::vector<ArrayData<NDIM, int>*> dofs_ptrs;
        if (cc_data) // semantics are slightly different for side-centered data
        {
            Pointer<CellData<NDIM, double> > value_data = patch->getPatchData(value_idx);
            Pointer<CellData<NDIM, int> > dof_data = patch->getPatchData(local_dof_idx);
            TBOX_ASSERT(value_data);
            TBOX_ASSERT(dof_data);

            dofs_ptrs.push_back(&dof_data->getArrayData());
            values_ptrs.push_back(&value_data->getArrayData());
        }
        else
        {
            Pointer<SideData<NDIM, double> > value_data = patch->getPatchData(value_idx);
            Pointer<SideData<NDIM, int> > dof_data = patch->getPatchData(local_dof_idx);
            TBOX_ASSERT(value_data);
            TBOX_ASSERT(dof_data);

            for (int d = 0; d < NDIM; ++d)
            {
                dofs_ptrs.push_back(&dof_data->getArrayData(d));
                values_ptrs.push_back(&value_data->getArrayData(d));
            }
        }

        static_assert(std::is_same<PetscInt, int>::value, "only implemented for 32-bit PETSc indices");
        Vec local;
        int ierr = VecGhostGetLocalForm(vec, &local);
        IBTK_CHKERRQ(ierr);
        ierr = VecSetOption(local, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
        IBTK_CHKERRQ(ierr);
        for (unsigned int d = 0; d < dofs_ptrs.size(); ++d)
        {
            const PetscInt size = dofs_ptrs[d]->getDepth() * dofs_ptrs[d]->getOffset();
            TBOX_ASSERT(size == values_ptrs[d]->getDepth() * values_ptrs[d]->getOffset());
            if (add)
                ierr = VecSetValues(local, size, dofs_ptrs[d]->getPointer(), values_ptrs[d]->getPointer(), ADD_VALUES);
            else
                ierr = VecGetValues(local, size, dofs_ptrs[d]->getPointer(), values_ptrs[d]->getPointer());

            IBTK_CHKERRQ(ierr);
        }
        ierr = VecGhostRestoreLocalForm(vec, &local);
        IBTK_CHKERRQ(ierr);
    }
    IBTK_TIMER_STOP(t_add_or_get);
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////
SAMRAIGhostDataAccumulator::SAMRAIGhostDataAccumulator(Pointer<BasePatchHierarchy<NDIM> > patch_hierarchy,
                                                       Pointer<Variable<NDIM> > var,
                                                       const IntVector<NDIM> gcw,
                                                       const int coarsest_ln,
                                                       const int finest_ln)
    : d_hierarchy(patch_hierarchy), d_var(var), d_gcw(gcw), d_coarsest_ln(coarsest_ln), d_finest_ln(finest_ln)
{
    auto set_timer = [&](const char *name)
    { return TimerManager::getManager()->getTimer(name); };
    t_constructor = set_timer("IBTK::SAMRAIGhostDataAccumulator::SAMRAIGhostDataAccumulator()");
    t_accumulate_ghost_data = set_timer("IBTK::SAMRAIGhostDataAccumulator::accumulateGhostData()");
    t_add_or_get = set_timer("IBTK::SAMRAIGhostDataAccumulator::accumulateGhostData()[add_or_get]");

    IBTK_TIMER_START(t_constructor);
    // Determine data layout:
    Pointer<CellVariable<NDIM, double> > cc_var = var;
    Pointer<SideVariable<NDIM, double> > sc_var = var;
    const bool cc_data = cc_var;
    const bool sc_data = sc_var;
    TBOX_ASSERT(cc_data || sc_data);
    d_cc_data = cc_data;

    // Create a context into which all indexing variables are grouped for this class:
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("SAMRAIGhostDataAccumulator");

    // Determine depth:
    int depth = 0;
    if (d_cc_data)
    {
        Pointer<CellDataFactory<NDIM, double> > cc_data_factory = cc_var->getPatchDataFactory();
        depth = cc_data_factory->getDefaultDepth();
    }
    else
    {
        Pointer<SideDataFactory<NDIM, double> > sc_data_factory = sc_var->getPatchDataFactory();
        depth = sc_data_factory->getDefaultDepth();
    }
    TBOX_ASSERT(depth != 0);

    const std::string name = "SAMRAIGhostDataAccumulator::dof_" + var->getName();

    // Create the dof indexing variable and its data:
    d_vecs.resize(d_finest_ln + 1); // be lazy and index this array directly by level number
    Pointer<Variable<NDIM> > dof_var;
    if (var_db->checkVariableExists(name))
        dof_var = var_db->getVariable(name);
    else
        dof_var = d_cc_data ? Pointer<Variable<NDIM> >(new CellVariable<NDIM, int>(name, depth)) :
                              Pointer<Variable<NDIM> >(new SideVariable<NDIM, int>(name, depth));

    d_global_dof_idx = var_db->registerVariableAndContext(dof_var, context, d_gcw);
    d_local_dof_idx = var_db->registerClonedPatchDataIndex(dof_var, d_global_dof_idx);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_global_dof_idx);
        level->allocatePatchData(d_local_dof_idx);

        std::vector<int> num_dofs_per_proc;
        PETScVecUtilities::constructPatchLevelDOFIndices(num_dofs_per_proc, d_global_dof_idx, level);
        const int mpi_rank = IBTK_MPI::getRank();

        // half-open range of DoFs on the current processor
        const int local_dofs_begin =
            std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
        const int local_dofs_end = local_dofs_begin + num_dofs_per_proc[mpi_rank];

        std::vector<PetscInt> ghosts;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            std::vector<ArrayData<NDIM, int>*> dofs_ptrs;
            if (d_cc_data)
            {
                Pointer<CellData<NDIM, int> > dof_data = patch->getPatchData(d_global_dof_idx);
                dofs_ptrs.push_back(&dof_data->getArrayData());
            }
            else
            {
                Pointer<SideData<NDIM, int> > dof_data = patch->getPatchData(d_global_dof_idx);
                // SideData stores multiple arrays (one for each dimension)
                for (int d = 0; d < NDIM; ++d) dofs_ptrs.push_back(&dof_data->getArrayData(d));
            }
            for (ArrayData<NDIM, int>* dofs : dofs_ptrs)
            {
                const int* dofs_ptr = dofs->getPointer();
                const int size = dofs->getBox().size() * dofs->getDepth();
                for (int n = 0; n < size; ++n)
                    if (dofs_ptr[n] < local_dofs_begin || local_dofs_end <= dofs_ptr[n]) ghosts.push_back(dofs_ptr[n]);
            }
        }

        // deduplicate
        std::sort(ghosts.begin(), ghosts.end());
        ghosts.erase(std::unique(ghosts.begin(), ghosts.end()), ghosts.end());
        // make sure that ghosts doesn't have any negative entries
        ghosts.erase(ghosts.begin(), std::upper_bound(ghosts.begin(), ghosts.end(), -1));
        int ierr = VecCreateGhost(PETSC_COMM_WORLD,
                                  num_dofs_per_proc[mpi_rank],
                                  PETSC_DECIDE,
                                  static_cast<PetscInt>(ghosts.size()),
                                  ghosts.data(),
                                  &d_vecs[ln]);
        IBTK_CHKERRQ(ierr);
        ierr = VecSetOption(d_vecs[ln], VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
        IBTK_CHKERRQ(ierr);

        // set up local indices
        ISLocalToGlobalMapping mapping;
        ierr = VecGetLocalToGlobalMapping(d_vecs[ln], &mapping);
        IBTK_CHKERRQ(ierr);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            std::vector<ArrayData<NDIM, int>*> dofs_ptrs;
            if (d_cc_data)
            {
                Pointer<CellData<NDIM, int> > global_dofs = patch->getPatchData(d_global_dof_idx);
                Pointer<CellData<NDIM, int> > local_dofs = patch->getPatchData(d_local_dof_idx);
                ArrayData<NDIM, int>& global_dof_data = global_dofs->getArrayData();
                ArrayData<NDIM, int>& local_dof_data = local_dofs->getArrayData();
                const auto size = local_dof_data.getBox().size() * local_dof_data.getDepth();
                // since these arrays contain ghost cells outside the physical domain they will always contain -1s
                ierr = ISGlobalToLocalMappingApply(
                    mapping, IS_GTOLM_MASK, size, global_dof_data.getPointer(), nullptr, local_dof_data.getPointer());
                IBTK_CHKERRQ(ierr);
            }
            else
            {
                Pointer<SideData<NDIM, int> > global_dofs = patch->getPatchData(d_global_dof_idx);
                Pointer<SideData<NDIM, int> > local_dofs = patch->getPatchData(d_local_dof_idx);
                // SideData stores multiple arrays (one for each dimension)
                for (int d = 0; d < NDIM; ++d)
                {
                    ArrayData<NDIM, int>& global_dof_data = global_dofs->getArrayData(d);
                    ArrayData<NDIM, int>& local_dof_data = local_dofs->getArrayData(d);
                    const auto size = local_dof_data.getBox().size() * local_dof_data.getDepth();
                    // since these arrays contain ghost cells outside the physical domain they will always contain -1s
                    ierr = ISGlobalToLocalMappingApply(mapping,
                                                       IS_GTOLM_MASK,
                                                       size,
                                                       global_dof_data.getPointer(),
                                                       nullptr,
                                                       local_dof_data.getPointer());
                    IBTK_CHKERRQ(ierr);
                }
            }
        }
    } // loop over levels
    IBTK_TIMER_STOP(t_constructor);
}

void
SAMRAIGhostDataAccumulator::accumulateGhostData(const int idx)
{
    IBTK_TIMER_START(t_accumulate_ghost_data);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;
    var_db->mapIndexToVariable(idx, var);
    TBOX_ASSERT(var == d_var);

    // the next part looks like PETScVecUtilities::copyToPatchLevelVec()
    // except here we include the ghost box and use ADD_VALUES since we
    // are accumulating.

    // 1. Copy data to the Vec:
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        // VecSet doesn't modify ghost values - use the local form instead
        Vec local;
        int ierr = VecGhostGetLocalForm(d_vecs[ln], &local);
        IBTK_CHKERRQ(ierr);
        ierr = VecSet(local, 0.0);
        IBTK_CHKERRQ(ierr);
        ierr = VecGhostRestoreLocalForm(d_vecs[ln], &local);
        IBTK_CHKERRQ(ierr);

        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        add_or_get(level, d_gcw, d_cc_data, d_local_dof_idx, idx, d_vecs[ln], true);
    }

    // 2. Accumulate:
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        const int ierr = VecGhostUpdateBegin(d_vecs[ln], ADD_VALUES, SCATTER_REVERSE);
        IBTK_CHKERRQ(ierr);
    }
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        const int ierr = VecGhostUpdateEnd(d_vecs[ln], ADD_VALUES, SCATTER_REVERSE);
        IBTK_CHKERRQ(ierr);
    }
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        const int ierr = VecGhostUpdateBegin(d_vecs[ln], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
    }
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        const int ierr = VecGhostUpdateEnd(d_vecs[ln], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
    }

    // 3. copy back:
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        add_or_get(level, d_gcw, d_cc_data, d_local_dof_idx, idx, d_vecs[ln], false);
    }
    IBTK_TIMER_STOP(t_accumulate_ghost_data);
}

SAMRAIGhostDataAccumulator::~SAMRAIGhostDataAccumulator()
{
    for (Vec& vec : d_vecs)
        if (vec != nullptr) VecDestroy(&vec);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (d_global_dof_idx != IBTK::invalid_index && level->checkAllocated(d_global_dof_idx))
            level->deallocatePatchData(d_global_dof_idx);
        if (d_local_dof_idx != IBTK::invalid_index && level->checkAllocated(d_local_dof_idx))
            level->deallocatePatchData(d_local_dof_idx);
    }
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    var_db->removePatchDataIndex(d_local_dof_idx);
    var_db->removePatchDataIndex(d_global_dof_idx);
}
/////////////////////////////// PRIVATE //////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
