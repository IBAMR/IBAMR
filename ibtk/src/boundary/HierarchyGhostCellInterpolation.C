// Filename: HierarchyGhostCellInterpolation.C
// Created on 05 Nov 2007 by Boyce Griffith
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

#include "HierarchyGhostCellInterpolation.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/CartCellDoubleCubicCoarsen.h>
#include <ibtk/CartCellDoubleQuadraticCFInterpolation.h>
#include <ibtk/CartSideDoubleCubicCoarsen.h>
#include <ibtk/CartSideDoubleQuadraticCFInterpolation.h>
#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/RefinePatchStrategySet.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <CellVariable.h>
#include <SideVariable.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

HierarchyGhostCellInterpolation::HierarchyGhostCellInterpolation()
    : d_is_initialized(false),
      d_homogeneous_bc(false),
      d_transaction_comps(),
      d_hierarchy(NULL),
      d_grid_geom(NULL),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_coarsen_alg(NULL),
      d_coarsen_strategy(NULL),
      d_coarsen_scheds(),
      d_refine_alg(NULL),
      d_refine_strategy(NULL),
      d_refine_scheds(),
      d_cf_bdry_ops(),
      d_extrap_bc_ops(),
      d_cc_robin_bc_ops(),
      d_sc_robin_bc_ops()
{
    // intentionally blank
    return;
}// HierarchyGhostCellInterpolation

HierarchyGhostCellInterpolation::~HierarchyGhostCellInterpolation()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
}// ~HierarchyGhostCellInterpolation

void
HierarchyGhostCellInterpolation::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    for (unsigned comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
    {
        if (!d_cc_robin_bc_ops[comp_idx].isNull()) d_cc_robin_bc_ops[comp_idx]->setHomogeneousBc(d_homogeneous_bc);
        for (std::vector<RobinBcCoefStrategy<NDIM>*>::iterator it =
                 d_transaction_comps[comp_idx].d_robin_bc_coefs.begin();
             it != d_transaction_comps[comp_idx].d_robin_bc_coefs.end(); ++it)
        {
            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(*it);
            if (extended_bc_coef != NULL)
            {
                extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
            }
        }
        if (!d_sc_robin_bc_ops[comp_idx].isNull()) d_sc_robin_bc_ops[comp_idx]->setHomogeneousBc(d_homogeneous_bc);
        for (std::vector<RobinBcCoefStrategy<NDIM>*>::iterator it =
                 d_transaction_comps[comp_idx].d_robin_bc_coefs.begin();
             it != d_transaction_comps[comp_idx].d_robin_bc_coefs.end(); ++it)
        {
            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(*it);
            if (extended_bc_coef != NULL)
            {
                extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
            }
        }
    }
    return;
}// setHomogeneousBc

void
HierarchyGhostCellInterpolation::initializeOperatorState(
    const InterpolationTransactionComponent transaction_comp,
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int coarsest_ln,
    const int finest_ln)
{
    initializeOperatorState(std::vector<InterpolationTransactionComponent>(1,transaction_comp), hierarchy, coarsest_ln, finest_ln);
    return;
}// initializeOperatorState

void
HierarchyGhostCellInterpolation::initializeOperatorState(
    const std::vector<InterpolationTransactionComponent> transaction_comps,
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int coarsest_ln,
    const int finest_ln)
{
    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Reset the transaction components.
    d_transaction_comps = transaction_comps;

    // Cache hierarchy data.
    d_hierarchy   = hierarchy;
    d_grid_geom   = d_hierarchy->getGridGeometry();
    d_coarsest_ln = coarsest_ln == -1 ? 0                                   : coarsest_ln;
    d_finest_ln   =   finest_ln == -1 ? d_hierarchy->getFinestLevelNumber() :   finest_ln;

    // Register the cubic coarsen operators with the grid geometry object.
    IBTK_DO_ONCE(
        d_grid_geom->addSpatialCoarsenOperator(new CartCellDoubleCubicCoarsen());
        d_grid_geom->addSpatialCoarsenOperator(new CartSideDoubleCubicCoarsen());
                 );

    // Setup cached coarsen algorithms and schedules.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    bool registered_coarsen_op = false;
    d_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    for (unsigned comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
    {
        const std::string& coarsen_op_name = d_transaction_comps[comp_idx].d_coarsen_op_name;
        if (coarsen_op_name != "NONE")
        {
            const int data_idx = d_transaction_comps[comp_idx].d_data_idx;
            Pointer<Variable<NDIM> > var;
            var_db->mapIndexToVariable(data_idx, var);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!var.isNull());
#endif
            Pointer<CoarsenOperator<NDIM> > coarsen_op =
                d_grid_geom->lookupCoarsenOperator(var, coarsen_op_name);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!coarsen_op.isNull());
#endif
            d_coarsen_alg->registerCoarsen(data_idx,  // destination
                                           data_idx,  // source
                                           coarsen_op);
            registered_coarsen_op = true;
        }
    }

    d_coarsen_strategy = NULL;

    d_coarsen_scheds.resize(d_finest_ln+1);
    if (registered_coarsen_op)
    {
        for (int src_ln = std::max(1,d_coarsest_ln); src_ln <= d_finest_ln; ++src_ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(src_ln);
            Pointer<PatchLevel<NDIM> > coarser_level = d_hierarchy->getPatchLevel(src_ln-1);
            d_coarsen_scheds[src_ln] = d_coarsen_alg->createSchedule(coarser_level, level, d_coarsen_strategy);
        }
    }

    // Setup cached refine algorithms and schedules.
    d_cf_bdry_ops.resize(d_transaction_comps.size());
    d_extrap_bc_ops.resize(d_transaction_comps.size());
    d_cc_robin_bc_ops.resize(d_transaction_comps.size());
    d_sc_robin_bc_ops.resize(d_transaction_comps.size());
    d_refine_alg = new RefineAlgorithm<NDIM>();
    std::vector<RefinePatchStrategy<NDIM>*> refine_patch_strategies;
    for (unsigned comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
    {
        const int data_idx = d_transaction_comps[comp_idx].d_data_idx;
        Pointer<Variable<NDIM> > var;
        var_db->mapIndexToVariable(data_idx, var);
        Pointer<CellVariable<NDIM,double> > cc_var = var;
        Pointer<NodeVariable<NDIM,double> > nc_var = var;
        Pointer<SideVariable<NDIM,double> > sc_var = var;
        Pointer<RefineOperator<NDIM> > refine_op = NULL;
        Pointer<VariableFillPattern<NDIM> > fill_pattern = d_transaction_comps[comp_idx].d_fill_pattern;
        if (!cc_var.isNull())
        {
            refine_op = NULL;
            d_cf_bdry_ops[comp_idx] = new CartCellDoubleQuadraticCFInterpolation();
            d_cf_bdry_ops[comp_idx]->setConsistentInterpolationScheme(d_transaction_comps[comp_idx].d_consistent_type_2_bdry);
            d_cf_bdry_ops[comp_idx]->setPatchDataIndex(data_idx);
            d_cf_bdry_ops[comp_idx]->setPatchHierarchy(d_hierarchy);
            refine_patch_strategies.push_back(d_cf_bdry_ops[comp_idx]);
        }
        else if (!nc_var.isNull())
        {
            refine_op = d_grid_geom->lookupRefineOperator(nc_var, "LINEAR_REFINE");
            d_cf_bdry_ops[comp_idx] = NULL;
        }
        else if (!sc_var.isNull())
        {
            refine_op = NULL;
            d_cf_bdry_ops[comp_idx] = new CartSideDoubleQuadraticCFInterpolation();
            d_cf_bdry_ops[comp_idx]->setConsistentInterpolationScheme(d_transaction_comps[comp_idx].d_consistent_type_2_bdry);
            d_cf_bdry_ops[comp_idx]->setPatchDataIndex(data_idx);
            d_cf_bdry_ops[comp_idx]->setPatchHierarchy(d_hierarchy);
            refine_patch_strategies.push_back(d_cf_bdry_ops[comp_idx]);
        }
        else
        {
            TBOX_ERROR("HierarchyGhostCellInterpolation::initializeOperatorState():\n"
                       << "  only double-precision cell-, node-, or side-centered data is presently supported." << std::endl);
        }

        d_refine_alg->registerRefine(data_idx,  // destination
                                     data_idx,  // source
                                     data_idx,  // temporary work space
                                     refine_op,
                                     fill_pattern);

        const std::string& phys_bdry_extrap_type = d_transaction_comps[comp_idx].d_phys_bdry_extrap_type;
        if (phys_bdry_extrap_type != "NONE")
        {
            d_extrap_bc_ops[comp_idx] = new CartExtrapPhysBdryOp(data_idx, phys_bdry_extrap_type);
            refine_patch_strategies.push_back(d_extrap_bc_ops[comp_idx]);
        }

        const std::vector<RobinBcCoefStrategy<NDIM>*>& robin_bc_coefs = d_transaction_comps[comp_idx].d_robin_bc_coefs;
        bool null_bc_coefs = true;
        for (std::vector<RobinBcCoefStrategy<NDIM>*>::const_iterator cit = robin_bc_coefs.begin();
             cit != robin_bc_coefs.end(); ++cit)
        {
            if (*cit != NULL) null_bc_coefs = false;
        }
        if (!null_bc_coefs && !cc_var.isNull())
        {
            d_cc_robin_bc_ops[comp_idx] = new CartCellRobinPhysBdryOp(data_idx, robin_bc_coefs, d_homogeneous_bc);
        }
        if (!null_bc_coefs && !sc_var.isNull())
        {
            d_sc_robin_bc_ops[comp_idx] = new CartSideRobinPhysBdryOp(data_idx, robin_bc_coefs, d_homogeneous_bc);
        }
    }

    d_refine_strategy = new RefinePatchStrategySet(refine_patch_strategies.begin(), refine_patch_strategies.end(), false);

    d_refine_scheds.resize(d_finest_ln+1);
    for (int dst_ln = d_coarsest_ln; dst_ln <= d_finest_ln; ++dst_ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(dst_ln);
        d_refine_scheds[dst_ln] = d_refine_alg->createSchedule(level, dst_ln-1, d_hierarchy, d_refine_strategy);
    }

    // Setup physical BC type.
    setHomogeneousBc(d_homogeneous_bc);

    // Indicate the operator is initialized.
    d_is_initialized = true;
    return;
}// initializeOperatorState

void
HierarchyGhostCellInterpolation::resetTransactionComponent(
    const InterpolationTransactionComponent& transaction_comp)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_is_initialized);
#endif
    if (d_transaction_comps.size() != 1)
    {
        TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponent():"
                   << "  invalid reset operation.  attempting to change the number of registered interpolation transaction components.\n");
    }
    resetTransactionComponents(std::vector<InterpolationTransactionComponent>(1,transaction_comp));
    return;
}// resetTransactionComponent

void
HierarchyGhostCellInterpolation::resetTransactionComponents(
    const std::vector<InterpolationTransactionComponent>& transaction_comps)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_is_initialized);
#endif
    if (d_transaction_comps.size() != transaction_comps.size())
    {
        TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():"
                   << "  invalid reset operation.  attempting to change the number of registered interpolation transaction components.\n");
    }

    // Reset the transaction components.
    d_transaction_comps = transaction_comps;

    // Reset cached coarsen algorithms and schedules.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    bool registered_coarsen_op = false;
    d_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    for (unsigned comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
    {
        const std::string& coarsen_op_name = d_transaction_comps[comp_idx].d_coarsen_op_name;
        if (coarsen_op_name != "NONE")
        {
            const int data_idx = d_transaction_comps[comp_idx].d_data_idx;
            Pointer<Variable<NDIM> > var;
            var_db->mapIndexToVariable(data_idx, var);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!var.isNull());
#endif
            Pointer<CoarsenOperator<NDIM> > coarsen_op =
                d_grid_geom->lookupCoarsenOperator(var, coarsen_op_name);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!coarsen_op.isNull());
#endif
            d_coarsen_alg->registerCoarsen(data_idx,  // destination
                                           data_idx,  // source
                                           coarsen_op);
            registered_coarsen_op = true;
        }
    }

    if (registered_coarsen_op)
    {
        for (int src_ln = std::max(1,d_coarsest_ln); src_ln <= d_finest_ln; ++src_ln)
        {
            d_coarsen_alg->resetSchedule(d_coarsen_scheds[src_ln]);
        }
    }

    // Reset cached refine algorithms and schedules.
    d_refine_alg = new RefineAlgorithm<NDIM>();
    for (unsigned comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
    {
        const int data_idx = d_transaction_comps[comp_idx].d_data_idx;
        Pointer<Variable<NDIM> > var;
        var_db->mapIndexToVariable(data_idx, var);
        Pointer<CellVariable<NDIM,double> > cc_var = var;
        Pointer<NodeVariable<NDIM,double> > nc_var = var;
        Pointer<SideVariable<NDIM,double> > sc_var = var;
        Pointer<RefineOperator<NDIM> > refine_op = NULL;
        Pointer<VariableFillPattern<NDIM> > fill_pattern = d_transaction_comps[comp_idx].d_fill_pattern;
        if (!d_cf_bdry_ops[comp_idx].isNull()) d_cf_bdry_ops[comp_idx]->setPatchDataIndex(data_idx);
        if (!cc_var.isNull())
        {
            refine_op = NULL;
        }
        else if (!nc_var.isNull())
        {
            refine_op = d_grid_geom->lookupRefineOperator(nc_var, "LINEAR_REFINE");
        }
        else if (!sc_var.isNull())
        {
            refine_op = NULL;
        }
        else
        {
            TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                       << "  only double-precision cell-, node-, or side-centered data is presently supported." << std::endl);
        }

        d_refine_alg->registerRefine(data_idx,  // destination
                                     data_idx,  // source
                                     data_idx,  // temporary work space
                                     refine_op,
                                     fill_pattern);

        const std::string& phys_bdry_extrap_type = d_transaction_comps[comp_idx].d_phys_bdry_extrap_type;
        if (!d_extrap_bc_ops[comp_idx].isNull())
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(phys_bdry_extrap_type != "NONE");
#endif
            d_extrap_bc_ops[comp_idx]->setPatchDataIndex(data_idx);
            d_extrap_bc_ops[comp_idx]->setExtrapolationType(phys_bdry_extrap_type);
        }
#ifdef DEBUG_CHECK_ASSERTIONS
        else
        {
            TBOX_ASSERT(phys_bdry_extrap_type == "NONE");
        }

#endif
        const std::vector<RobinBcCoefStrategy<NDIM>*>& robin_bc_coefs = d_transaction_comps[comp_idx].d_robin_bc_coefs;
        bool null_bc_coefs = true;
        for (std::vector<RobinBcCoefStrategy<NDIM>*>::const_iterator cit = robin_bc_coefs.begin();
             cit != robin_bc_coefs.end(); ++cit)
        {
            if (*cit != NULL) null_bc_coefs = false;
        }
        if (!d_cc_robin_bc_ops[comp_idx].isNull())
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!null_bc_coefs);
            TBOX_ASSERT(!cc_var.isNull());
#endif
            d_cc_robin_bc_ops[comp_idx]->setPhysicalBcCoefs(robin_bc_coefs);
            d_cc_robin_bc_ops[comp_idx]->setPatchDataIndex(data_idx);
        }
        if (!d_sc_robin_bc_ops[comp_idx].isNull())
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!null_bc_coefs);
            TBOX_ASSERT(!sc_var.isNull());
#endif
            d_sc_robin_bc_ops[comp_idx]->setPhysicalBcCoefs(robin_bc_coefs);
            d_sc_robin_bc_ops[comp_idx]->setPatchDataIndex(data_idx);
        }
    }

    for (int dst_ln = d_coarsest_ln; dst_ln <= d_finest_ln; ++dst_ln)
    {
        d_refine_alg->resetSchedule(d_refine_scheds[dst_ln]);
    }
    return;
}// resetTransactionComponents

void
HierarchyGhostCellInterpolation::reinitializeOperatorState(
    Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    if (!d_is_initialized) return;
    initializeOperatorState(d_transaction_comps, hierarchy);
    return;
}// reinitializeOperatorState

void
HierarchyGhostCellInterpolation::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    // Clear cached refinement operators.
    d_cf_bdry_ops.clear();
    d_extrap_bc_ops.clear();
    d_cc_robin_bc_ops.clear();
    d_sc_robin_bc_ops.clear();

    // Clear cached communication schedules.
    d_coarsen_alg.setNull();
    if (d_coarsen_strategy != NULL) delete d_coarsen_strategy;
    d_coarsen_strategy = NULL;
    d_coarsen_scheds.clear();

    d_refine_alg.setNull();
    if (d_refine_strategy != NULL) delete d_refine_strategy;
    d_refine_strategy = NULL;
    d_refine_scheds.clear();

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;
    return;
}// deallocateOperatorState

void
HierarchyGhostCellInterpolation::fillData(
    const double& fill_time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_is_initialized);
#endif
    // Ensure the boundary condition objects are in the correct state.
    for (unsigned comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
    {
        if (!d_cc_robin_bc_ops[comp_idx].isNull()) d_cc_robin_bc_ops[comp_idx]->setHomogeneousBc(d_homogeneous_bc);
        for (std::vector<RobinBcCoefStrategy<NDIM>*>::iterator it =
                 d_transaction_comps[comp_idx].d_robin_bc_coefs.begin();
             it != d_transaction_comps[comp_idx].d_robin_bc_coefs.end(); ++it)
        {
            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(*it);
            if (extended_bc_coef != NULL)
            {
                extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
            }
        }
        if (!d_sc_robin_bc_ops[comp_idx].isNull()) d_sc_robin_bc_ops[comp_idx]->setHomogeneousBc(d_homogeneous_bc);
        for (std::vector<RobinBcCoefStrategy<NDIM>*>::iterator it =
                 d_transaction_comps[comp_idx].d_robin_bc_coefs.begin();
             it != d_transaction_comps[comp_idx].d_robin_bc_coefs.end(); ++it)
        {
            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(*it);
            if (extended_bc_coef != NULL)
            {
                extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
            }
        }
    }

    // Synchronize data on the patch hierarchy prior to filling ghost cell
    // values.
    for (int src_ln = d_finest_ln; src_ln >= std::max(1,d_coarsest_ln); --src_ln)
    {
        if (!d_coarsen_scheds[src_ln].isNull()) d_coarsen_scheds[src_ln]->coarsenData();
    }

    // Perform the initial data fill, using extrapolation to determine ghost
    // cell values at physical boundaries.
    for (int dst_ln = d_coarsest_ln; dst_ln <= d_finest_ln; ++dst_ln)
    {
        if (!d_refine_scheds[dst_ln].isNull())
        {
            d_refine_scheds[dst_ln]->fillData(fill_time);
        }
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(dst_ln);
        const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            for (unsigned comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
            {
                if (!d_cf_bdry_ops[comp_idx].isNull())
                {
                    const int data_idx = d_transaction_comps[comp_idx].d_data_idx;
                    const IntVector<NDIM>& ghost_width_to_fill = patch->getPatchData(data_idx)->getGhostCellWidth();
                    d_cf_bdry_ops[comp_idx]->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
                }
            }
        }
    }

    // Set Robin boundary conditions at physical boundaries.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            if (patch->getPatchGeometry()->getTouchesRegularBoundary())
            {
                for (unsigned comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
                {
                    if (!d_cc_robin_bc_ops[comp_idx].isNull())
                    {
                        const int data_idx = d_transaction_comps[comp_idx].d_data_idx;
                        const IntVector<NDIM>& ghost_width_to_fill = patch->getPatchData(data_idx)->getGhostCellWidth();
                        d_cc_robin_bc_ops[comp_idx]->setPhysicalBoundaryConditions(*patch, fill_time, ghost_width_to_fill);
                    }
                    if (!d_sc_robin_bc_ops[comp_idx].isNull())
                    {
                        const int data_idx = d_transaction_comps[comp_idx].d_data_idx;
                        const IntVector<NDIM>& ghost_width_to_fill = patch->getPatchData(data_idx)->getGhostCellWidth();
                        d_sc_robin_bc_ops[comp_idx]->setPhysicalBoundaryConditions(*patch, fill_time, ghost_width_to_fill);
                    }
                }
            }
        }
    }
    return;
}// fillData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::HierarchyGhostCellInterpolation>;

//////////////////////////////////////////////////////////////////////////////
