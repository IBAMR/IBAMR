// Filename: EdgeDataSynchronization.cpp
// Created on 02 Feb 2011 by Boyce Griffith
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
#include <ostream>
#include <string>
#include <vector>

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/xfer/VariableFillPattern.h"
#include "boost/array.hpp"
#include "ibtk/EdgeDataSynchronization.h"
#include "ibtk/EdgeSynchCopyFillPattern.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI
{
namespace xfer
{

class CoarsenPatchStrategy;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

EdgeDataSynchronization::EdgeDataSynchronization()
    : d_is_initialized(false), d_transaction_comps(), d_coarsest_ln(-1), d_finest_ln(-1), d_coarsen_alg(NULL),
      d_coarsen_scheds(), d_refine_alg(), d_refine_scheds()
{
    // intentionally blank
    return;
}

EdgeDataSynchronization::~EdgeDataSynchronization()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
}

void EdgeDataSynchronization::initializeOperatorState(const SynchronizationTransactionComponent& transaction_comp,
                                                      const boost::shared_ptr<PatchHierarchy>& hierarchy)
{
    initializeOperatorState(std::vector<SynchronizationTransactionComponent>(1, transaction_comp), hierarchy);
    return;
}

void EdgeDataSynchronization::initializeOperatorState(
    const std::vector<SynchronizationTransactionComponent>& transaction_comps,
    const boost::shared_ptr<PatchHierarchy>& hierarchy)
{
    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Reset the transaction components.
    d_transaction_comps = transaction_comps;

    // Cache hierarchy data.
    d_hierarchy = hierarchy;
    d_grid_geom = BOOST_CAST<CartesianGridGeometry>(d_hierarchy->getGridGeometry());
    d_coarsest_ln = 0;
    d_finest_ln = d_hierarchy->getFinestLevelNumber();

    // Setup cached coarsen algorithms and schedules.
    auto var_db = VariableDatabase::getDatabase();
    bool registered_coarsen_op = false;
    d_coarsen_alg = boost::make_shared<CoarsenAlgorithm>(DIM);
    for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
    {
        const std::string& coarsen_op_name = d_transaction_comps[comp_idx].d_coarsen_op_name;
        if (coarsen_op_name != "NONE")
        {
            const int data_idx = d_transaction_comps[comp_idx].d_data_idx;
            boost::shared_ptr<Variable> var;
            var_db->mapIndexToVariable(data_idx, var);
            auto coarsen_op = d_grid_geom->lookupCoarsenOperator(var, coarsen_op_name);
            d_coarsen_alg->registerCoarsen(data_idx, data_idx, coarsen_op);
            registered_coarsen_op = true;
        }
    }

    CoarsenPatchStrategy* coarsen_strategy = NULL;
    d_coarsen_scheds.resize(d_finest_ln + 1);
    if (registered_coarsen_op)
    {
        for (int ln = d_coarsest_ln + 1; ln <= d_finest_ln; ++ln)
        {
            auto level = d_hierarchy->getPatchLevel(ln);
            auto coarser_level = d_hierarchy->getPatchLevel(ln - 1);
            d_coarsen_scheds[ln] = d_coarsen_alg->createSchedule(coarser_level, level, coarsen_strategy);
        }
    }

    // Setup cached refine algorithms and schedules.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        d_refine_alg[axis] = boost::make_shared<RefineAlgorithm>();
        for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
        {
            const int data_idx = d_transaction_comps[comp_idx].d_data_idx;
#ifndef NDEBUG
            boost::shared_ptr<Variable> var;
            var_db->mapIndexToVariable(data_idx, var);
            auto ec_var = boost::dynamic_pointer_cast<EdgeVariable<double> >(var);
            if (!ec_var)
            {
                TBOX_ERROR("EdgeDataSynchronization::initializeOperatorState():\n"
                           << "  only double-precision edge-centered data is supported." << std::endl);
            }
#endif
            const boost::shared_ptr<RefineOperator> no_refine_op;
            auto fill_pattern = boost::make_shared<EdgeSynchCopyFillPattern>(axis);
            d_refine_alg[axis]->registerRefine(data_idx, data_idx, data_idx, no_refine_op, fill_pattern);
        }

        d_refine_scheds[axis].resize(d_finest_ln + 1);
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            auto level = d_hierarchy->getPatchLevel(ln);
            d_refine_scheds[axis][ln] = d_refine_alg[axis]->createSchedule(level);
        }
    }

    // Indicate the operator is initialized.
    d_is_initialized = true;
    return;
}

void EdgeDataSynchronization::resetTransactionComponent(const SynchronizationTransactionComponent& transaction_comp)
{
    TBOX_ASSERT(d_is_initialized);
    if (d_transaction_comps.size() != 1)
    {
        TBOX_ERROR("EdgeDataSynchronization::resetTransactionComponent():"
                   << "  invalid reset operation.  attempting to change the number of registered "
                      "synchronization transaction components.\n");
    }
    resetTransactionComponents(std::vector<SynchronizationTransactionComponent>(1, transaction_comp));
    return;
}

void EdgeDataSynchronization::resetTransactionComponents(
    const std::vector<SynchronizationTransactionComponent>& transaction_comps)
{
    TBOX_ASSERT(d_is_initialized);
    if (d_transaction_comps.size() != transaction_comps.size())
    {
        TBOX_ERROR("EdgeDataSynchronization::resetTransactionComponents():"
                   << "  invalid reset operation.  attempting to change the number of registered "
                      "synchronization transaction components.\n");
    }

    // Reset the transaction components.
    d_transaction_comps = transaction_comps;

    // Reset cached coarsen algorithms and schedules.
    auto var_db = VariableDatabase::getDatabase();
    bool registered_coarsen_op = false;
    d_coarsen_alg = boost::make_shared<CoarsenAlgorithm>(DIM);
    for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
    {
        const std::string& coarsen_op_name = d_transaction_comps[comp_idx].d_coarsen_op_name;
        if (coarsen_op_name != "NONE")
        {
            const int data_idx = d_transaction_comps[comp_idx].d_data_idx;
            boost::shared_ptr<Variable> var;
            var_db->mapIndexToVariable(data_idx, var);
            auto coarsen_op = d_grid_geom->lookupCoarsenOperator(var, coarsen_op_name);
            d_coarsen_alg->registerCoarsen(data_idx, data_idx, coarsen_op);
            registered_coarsen_op = true;
        }
    }

    if (registered_coarsen_op)
    {
        for (int ln = d_coarsest_ln + 1; ln <= d_finest_ln; ++ln)
        {
            d_coarsen_alg->resetSchedule(d_coarsen_scheds[ln]);
        }
    }

    // Reset cached refine algorithms and schedules.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        d_refine_alg[axis] = boost::make_shared<RefineAlgorithm>();
        for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
        {
            const int data_idx = d_transaction_comps[comp_idx].d_data_idx;
#ifndef NDEBUG
            boost::shared_ptr<Variable> var;
            var_db->mapIndexToVariable(data_idx, var);
            auto ec_var = boost::dynamic_pointer_cast<EdgeVariable<double> >(var);
            if (!ec_var)
            {
                TBOX_ERROR("EdgeDataSynchronization::resetTransactionComponents():\n"
                           << "  only double-precision edge-centered data is supported." << std::endl);
            }
#endif
            const boost::shared_ptr<RefineOperator> no_refine_op;
            auto fill_pattern = boost::make_shared<EdgeSynchCopyFillPattern>(axis);
            d_refine_alg[axis]->registerRefine(data_idx, data_idx, data_idx, no_refine_op, fill_pattern);
        }

        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_refine_alg[axis]->resetSchedule(d_refine_scheds[axis][ln]);
        }
    }
    return;
}

void EdgeDataSynchronization::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    // Clear cached communication schedules.
    d_coarsen_alg.reset();
    d_coarsen_scheds.clear();

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        d_refine_alg[axis].reset();
        d_refine_scheds[axis].clear();
    }

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;
    return;
}

void EdgeDataSynchronization::synchronizeData(const double fill_time)
{
    TBOX_ASSERT(d_is_initialized);
    for (int ln = d_finest_ln; ln >= d_coarsest_ln; --ln)
    {
        auto level = d_hierarchy->getPatchLevel(ln);

        // Synchronize data on the current level.
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            d_refine_scheds[axis][ln]->fillData(fill_time);
        }

        // When appropriate, coarsen data from the current level to the next
        // coarser level.
        if (ln > d_coarsest_ln && d_coarsen_scheds[ln]) d_coarsen_scheds[ln]->coarsenData();
    }
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
