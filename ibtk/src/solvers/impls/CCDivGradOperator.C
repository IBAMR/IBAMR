// Filename: CCDivGradOperator.C
// Created on 18 Feb 2010 by Boyce Griffith
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

#include "CCDivGradOperator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <CellDataFactory.h>
#include <CellVariable.h>
#include <HierarchyDataOpsManager.h>
#include <Variable.h>
#include <tbox/Database.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Timers.
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CCDivGradOperator::CCDivGradOperator(
    const std::string& object_name)
    : LinearOperator(true),
      d_object_name(object_name),
      d_is_initialized(false),
      d_scalar_hier_bdry_fill(NULL),
      d_vector_hier_bdry_fill(NULL),
      d_scalar_scratch_idx(),
      d_vector_scratch_idx(),
      d_scalar_scratch_var(),
      d_vector_scratch_var(),
      d_hier_cc_data_ops(),
      d_hier_math_ops(),
      d_hier_math_ops_external(false),
      d_hierarchy(),
      d_coarsest_ln(-1),
      d_finest_ln(-1)
{
    // Setup Timers.
    IBTK_DO_ONCE(
        t_apply                     = TimerManager::getManager()->getTimer("IBTK::CCDivGradOperator::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBTK::CCDivGradOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBTK::CCDivGradOperator::deallocateOperatorState()");
                 );
    return;
}// CCDivGradOperator()

CCDivGradOperator::~CCDivGradOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
}// ~CCDivGradOperator()

void
CCDivGradOperator::setHierarchyMathOps(
    Pointer<HierarchyMathOps> hier_math_ops)
{
    d_hier_math_ops = hier_math_ops;
    d_hier_math_ops_external = !d_hier_math_ops.isNull();
    return;
}// setHierarchyMathOps

void
CCDivGradOperator::modifyRhsForInhomogeneousBc(
    SAMRAIVectorReal<NDIM,double>& y)
{
    // intentionally blank
    return;
}// modifyRhsForInhomogeneousBc

void
CCDivGradOperator::apply(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    SAMRAI_MPI::barrier();
    t_apply->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_is_initialized);
#endif

    static const int comp = 0;

    const Pointer<Variable<NDIM> >& x_var = x.getComponentVariable(comp);
    const Pointer<Variable<NDIM> >& y_var = y.getComponentVariable(comp);

    Pointer<CellVariable<NDIM,double> > x_cc_var = x_var;
    Pointer<CellVariable<NDIM,double> > y_cc_var = y_var;

    if (x_cc_var.isNull() || y_cc_var.isNull())
    {
        TBOX_ERROR(d_object_name << "::apply()\n"
                   << "  encountered non-cell centered vector components" << std::endl);
    }

#ifdef DEBUG_CHECK_ASSERTIONS
    Pointer<CellDataFactory<NDIM,double> > x_factory = x_cc_var->getPatchDataFactory();
    Pointer<CellDataFactory<NDIM,double> > y_factory = y_cc_var->getPatchDataFactory();

    TBOX_ASSERT(!x_factory.isNull());
    TBOX_ASSERT(!y_factory.isNull());

    const unsigned int x_depth = x_factory->getDefaultDepth();
    const unsigned int y_depth = y_factory->getDefaultDepth();

    TBOX_ASSERT(x_depth == 1);
    TBOX_ASSERT(y_depth == 1);
#endif

    const int x_idx = x.getComponentDescriptorIndex(comp);
    const int y_idx = y.getComponentDescriptorIndex(comp);

    // Compute the action of the operator.
    d_hier_cc_data_ops->copyData(d_scalar_scratch_idx, x_idx);
    d_hier_math_ops->grad(
        d_vector_scratch_idx, d_vector_scratch_var,
        +1.0, d_scalar_scratch_idx, d_scalar_scratch_var,
        d_scalar_hier_bdry_fill, 0.0,
        0.0, -1, Pointer<CellVariable<NDIM,double> >(NULL),
        0);
    d_hier_math_ops->div(
        y_idx, y_cc_var,
        -1.0, d_vector_scratch_idx, d_vector_scratch_var,
        d_vector_hier_bdry_fill, 0.0,
        0.0, -1, Pointer<CellVariable<NDIM,double> >(NULL),
        0, 0);

    t_apply->stop();
    return;
}// apply

void
CCDivGradOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAIVectorReal<NDIM,double>& out)
{
    SAMRAI_MPI::barrier();
    t_initialize_operator_state->start();

    static const int comp = 0;

    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup operator state.
    d_hierarchy   = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln   = in.getFinestLevelNumber();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
#endif

    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<CellVariable<NDIM,double> > cc_var = new CellVariable<NDIM,double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, d_hierarchy, true);
    d_hier_cc_data_ops->resetLevels(d_coarsest_ln, d_finest_ln);

    if (!d_hier_math_ops_external)
    {
        d_hier_math_ops = new HierarchyMathOps(
            d_object_name+"::HierarchyMathOps", d_hierarchy, d_coarsest_ln, d_finest_ln);
    }
    else
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!d_hier_math_ops.isNull());
#endif
    }

    // Get variables and patch data descriptors for temporary data.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const Pointer<Variable<NDIM> >& in_var = in.getComponentVariable(comp);
    const int in_idx = in.getComponentDescriptorIndex(comp);
    d_scalar_scratch_idx = var_db->registerClonedPatchDataIndex(in_var, in_idx);
    d_vector_scratch_idx = var_db->registerClonedPatchDataIndex(in_var, in_idx);
    var_db->mapIndexToVariable(d_scalar_scratch_idx, d_scalar_scratch_var);
    var_db->mapIndexToVariable(d_vector_scratch_idx, d_vector_scratch_var);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(in_idx >= 0);
    TBOX_ASSERT(d_scalar_scratch_idx >= 0);
    TBOX_ASSERT(d_vector_scratch_idx >= 0);

    Pointer<CellDataFactory<NDIM,double> > scalar_scratch_factory =
        var_db->getPatchDescriptor()->getPatchDataFactory(d_scalar_scratch_idx);
    TBOX_ASSERT(!scalar_scratch_factory.isNull());
    TBOX_ASSERT(scalar_scratch_factory->getGhostCellWidth() == CELLG);
    TBOX_ASSERT(scalar_scratch_factory->getDefaultDepth() == 1);
#endif

    Pointer<CellDataFactory<NDIM,double> > vector_scratch_factory =
        var_db->getPatchDescriptor()->getPatchDataFactory(d_vector_scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!vector_scratch_factory.isNull());
    TBOX_ASSERT(vector_scratch_factory->getGhostCellWidth() == CELLG);
    TBOX_ASSERT(vector_scratch_factory->getDefaultDepth() == 1);
#endif
    vector_scratch_factory->setDefaultDepth(NDIM);

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scalar_scratch_idx, 0.0);
        level->allocatePatchData(d_vector_scratch_idx, 0.0);
    }

    // Setup the interpolation transaction information.
    Pointer<VariableFillPattern<NDIM> > fill_pattern = new CellNoCornersFillPattern(CELLG, false, false, true);

    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;

    std::vector<InterpolationTransactionComponent> scalar_transaction_comps;
    InterpolationTransactionComponent scalar_component(d_scalar_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL, fill_pattern);
    scalar_transaction_comps.push_back(scalar_component);

    std::vector<InterpolationTransactionComponent> vector_transaction_comps;
    InterpolationTransactionComponent vector_component(d_vector_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL, fill_pattern);
    vector_transaction_comps.push_back(vector_component);

    // Initialize the interpolation operators.
    d_scalar_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_scalar_hier_bdry_fill->initializeOperatorState(scalar_transaction_comps, d_hierarchy);

    d_vector_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_vector_hier_bdry_fill->initializeOperatorState(vector_transaction_comps, d_hierarchy);

    // Indicate the operator is initialized.
    d_is_initialized = true;

    t_initialize_operator_state->stop();
    return;
}// initializeOperatorState

void
CCDivGradOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    SAMRAI_MPI::barrier();
    t_deallocate_operator_state->start();

    // Deallocate the interpolation operators.
    d_scalar_hier_bdry_fill->deallocateOperatorState();
    d_vector_hier_bdry_fill->deallocateOperatorState();
    d_scalar_hier_bdry_fill.setNull();
    d_vector_hier_bdry_fill.setNull();

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        if (ln <= d_hierarchy->getFinestLevelNumber())
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(d_scalar_scratch_idx);
            level->deallocatePatchData(d_vector_scratch_idx);
        }
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    var_db->removePatchDataIndex(d_scalar_scratch_idx);
    var_db->removePatchDataIndex(d_vector_scratch_idx);

    if (!d_hier_math_ops_external) d_hier_math_ops.setNull();

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;

    t_deallocate_operator_state->stop();
    return;
}// deallocateOperatorState

void
CCDivGradOperator::enableLogging(
    bool enabled)
{
    TBOX_WARNING("CCDivGradOperator::enableLogging() not supported" << std::endl);
    return;
}// enableLogging

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::CCDivGradOperator>;

//////////////////////////////////////////////////////////////////////////////
