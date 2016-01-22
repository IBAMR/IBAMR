// Filename: CCLaplaceOperator.cpp
// Created on 24 Oct 2003 by Boyce Griffith
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

#include "CellDataFactory.h"
#include "CellVariable.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "PoissonSpecifications.h"
#include "SAMRAIVectorReal.h"
#include "VariableFillPattern.h"
#include "ibtk/CCLaplaceOperator.h"
#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LaplaceOperator.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
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

CCLaplaceOperator::CCLaplaceOperator(const std::string& object_name, const bool homogeneous_bc)
    : LaplaceOperator(object_name, homogeneous_bc),
      d_ncomp(0),
      d_fill_pattern(NULL),
      d_transaction_comps(),
      d_hier_bdry_fill(NULL),
      d_no_fill(NULL),
      d_x(NULL),
      d_b(NULL),
      d_hierarchy(),
      d_coarsest_ln(-1),
      d_finest_ln(-1)
{
    // Setup the operator to use default scalar-valued boundary conditions.
    setPhysicalBcCoef(NULL);

    // Setup Timers.
    IBTK_DO_ONCE(t_apply = TimerManager::getManager()->getTimer("IBTK::CCLaplaceOperator::apply()");
                 t_initialize_operator_state =
                     TimerManager::getManager()->getTimer("IBTK::CCLaplaceOperator::initializeOperatorState()");
                 t_deallocate_operator_state =
                     TimerManager::getManager()->getTimer("IBTK::CCLaplaceOperator::deallocateOperatorState()"););
    return;
} // CCLaplaceOperator()

CCLaplaceOperator::~CCLaplaceOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
} // ~CCLaplaceOperator()

void
CCLaplaceOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    IBTK_TIMER_START(t_apply);

#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        Pointer<CellVariable<NDIM, double> > x_cc_var = x.getComponentVariable(comp);
        Pointer<CellVariable<NDIM, double> > y_cc_var = y.getComponentVariable(comp);
        if (!x_cc_var || !y_cc_var)
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                                     << "  encountered non-cell centered vector components"
                                     << std::endl);
        }
        Pointer<CellDataFactory<NDIM, double> > x_factory = x_cc_var->getPatchDataFactory();
        Pointer<CellDataFactory<NDIM, double> > y_factory = y_cc_var->getPatchDataFactory();
        TBOX_ASSERT(x_factory);
        TBOX_ASSERT(y_factory);
        const unsigned int x_depth = x_factory->getDefaultDepth();
        const unsigned int y_depth = y_factory->getDefaultDepth();
        TBOX_ASSERT(x_depth == y_depth);
        if (x_depth != d_bc_coefs.size() || y_depth != d_bc_coefs.size())
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                                     << "  each vector component must have data depth == "
                                     << d_bc_coefs.size()
                                     << "\n"
                                     << "  since d_bc_coefs.size() == "
                                     << d_bc_coefs.size()
                                     << std::endl);
        }
    }
#endif

    // Simultaneously fill ghost cell values for all components.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps;
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        InterpolationTransactionComponent x_component(x.getComponentDescriptorIndex(comp),
                                                      DATA_REFINE_TYPE,
                                                      USE_CF_INTERPOLATION,
                                                      DATA_COARSEN_TYPE,
                                                      BDRY_EXTRAP_TYPE,
                                                      CONSISTENT_TYPE_2_BDRY,
                                                      d_bc_coefs,
                                                      d_fill_pattern);
        transaction_comps.push_back(x_component);
    }
    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time);
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);

    // Compute the action of the operator.
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        Pointer<CellVariable<NDIM, double> > x_cc_var = x.getComponentVariable(comp);
        Pointer<CellVariable<NDIM, double> > y_cc_var = y.getComponentVariable(comp);
        const int x_idx = x.getComponentDescriptorIndex(comp);
        const int y_idx = y.getComponentDescriptorIndex(comp);
        for (unsigned int l = 0; l < d_bc_coefs.size(); ++l)
        {
            d_hier_math_ops->laplace(y_idx,
                                     y_cc_var,
                                     d_poisson_spec,
                                     x_idx,
                                     x_cc_var,
                                     d_no_fill,
                                     0.0,
                                     0.0,
                                     -1,
                                     Pointer<CellVariable<NDIM, double> >(NULL),
                                     l,
                                     l);
        }
    }

    IBTK_TIMER_STOP(t_apply);
    return;
} // apply

void
CCLaplaceOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                           const SAMRAIVectorReal<NDIM, double>& out)
{
    IBTK_TIMER_START(t_initialize_operator_state);

    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup solution and rhs vectors.
    d_x = in.cloneVector(in.getName());
    d_b = out.cloneVector(out.getName());

    // Setup operator state.
    d_hierarchy = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln = in.getFinestLevelNumber();

    d_ncomp = in.getNumberOfComponents();

#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
    TBOX_ASSERT(d_ncomp == out.getNumberOfComponents());
#endif

    if (!d_hier_math_ops_external)
    {
        d_hier_math_ops =
            new HierarchyMathOps(d_object_name + "::HierarchyMathOps", d_hierarchy, d_coarsest_ln, d_finest_ln);
    }
    else
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_hier_math_ops);
#endif
    }

    // Setup the interpolation transaction information.
    d_fill_pattern = NULL;
    if (d_poisson_spec.dIsConstant())
    {
        d_fill_pattern = new CellNoCornersFillPattern(CELLG, false, false, true);
    }
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    d_transaction_comps.clear();
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        InterpolationTransactionComponent component(d_x->getComponentDescriptorIndex(comp),
                                                    DATA_REFINE_TYPE,
                                                    USE_CF_INTERPOLATION,
                                                    DATA_COARSEN_TYPE,
                                                    BDRY_EXTRAP_TYPE,
                                                    CONSISTENT_TYPE_2_BDRY,
                                                    d_bc_coefs,
                                                    d_fill_pattern);
        d_transaction_comps.push_back(component);
    }

    // Initialize the interpolation operators.
    d_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(d_transaction_comps, d_hierarchy, d_coarsest_ln, d_finest_ln);

    // Indicate the operator is initialized.
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
CCLaplaceOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_operator_state);

    // Deallocate the interpolation operators.
    d_hier_bdry_fill->deallocateOperatorState();
    d_hier_bdry_fill.setNull();
    d_transaction_comps.clear();
    d_fill_pattern.setNull();

    // Deallocate hierarchy math operations object.
    if (!d_hier_math_ops_external) d_hier_math_ops.setNull();

    // Delete the solution and rhs vectors.
    d_x->freeVectorComponents();
    d_x.setNull();

    d_b->freeVectorComponents();
    d_b.setNull();

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
