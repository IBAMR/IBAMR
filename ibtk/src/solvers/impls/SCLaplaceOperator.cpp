// Filename: SCLaplaceOperator.cpp
// Created on 11 Apr 2008 by Boyce Griffith
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

#include "SAMRAI/hier/IntVector.h"

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/pdat/SideDataFactory.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/xfer/VariableFillPattern.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LaplaceOperator.h"
#include "ibtk/SCLaplaceOperator.h"
#include "ibtk/SideNoCornersFillPattern.h"
#include "ibtk/StaggeredPhysicalBoundaryHelper.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI
{
namespace solv
{

class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int SIDEG = 1;

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

SCLaplaceOperator::SCLaplaceOperator(const std::string& object_name, const bool homogeneous_bc)
    : LaplaceOperator(object_name, homogeneous_bc), d_ncomp(0), d_fill_pattern(NULL), d_transaction_comps(),
      d_hier_bdry_fill(NULL), d_no_fill(NULL), d_x(NULL), d_b(NULL), d_hierarchy(), d_coarsest_ln(-1), d_finest_ln(-1)
{
    // Setup the operator to use default vector-valued boundary conditions.
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy*>(NDIM, static_cast<RobinBcCoefStrategy*>(NULL)));

    // Setup Timers.
    IBTK_DO_ONCE(t_apply = TimerManager::getManager()->getTimer("IBTK::SCLaplaceOperator::apply()");
                 t_initialize_operator_state =
                     TimerManager::getManager()->getTimer("IBTK::SCLaplaceOperator::initializeOperatorState()");
                 t_deallocate_operator_state =
                     TimerManager::getManager()->getTimer("IBTK::SCLaplaceOperator::deallocateOperatorState()"););
    return;
} // SCLaplaceOperator()

SCLaplaceOperator::~SCLaplaceOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
} // ~SCLaplaceOperator()

void SCLaplaceOperator::apply(SAMRAIVectorReal<double>& x, SAMRAIVectorReal<double>& y)
{
    IBTK_TIMER_START(t_apply);

    TBOX_ASSERT(d_is_initialized);
    TBOX_ASSERT(d_bc_coefs.size() == NDIM);
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        boost::shared_ptr<SideVariable<double> > x_sc_var = x.getComponentVariable(comp);
        boost::shared_ptr<SideVariable<double> > y_sc_var = y.getComponentVariable(comp);
        if (!x_sc_var || !y_sc_var)
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                                     << "  encountered non-side centered vector components" << std::endl);
        }
        const unsigned int x_depth = x_sc_var->getDepth();
        const unsigned int y_depth = y_sc_var->getDepth();
        TBOX_ASSERT(x_depth == y_depth);
        if (x_depth != 1 || y_depth != 1)
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                                     << "  each vector component must have data depth == 1" << std::endl);
        }
    }

    // Allocate scratch data.
    if (d_x) d_x->allocateVectorData();

    // Simultaneously fill ghost cell values for all components.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps;
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        InterpolationTransactionComponent x_component(d_x->getComponentDescriptorIndex(comp),
                                                      x.getComponentDescriptorIndex(comp),
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
        boost::shared_ptr<SideVariable<double> > x_sc_var = x.getComponentVariable(comp);
        boost::shared_ptr<SideVariable<double> > y_sc_var = y.getComponentVariable(comp);
        const int x_scratch_idx = d_x->getComponentDescriptorIndex(comp);
        const int y_idx = y.getComponentDescriptorIndex(comp);
        d_hier_math_ops->laplace(y_idx, y_sc_var, d_poisson_spec, x_scratch_idx, x_sc_var, d_no_fill, 0.0);
        const int x_idx = x.getComponentDescriptorIndex(comp);
        d_bc_helpers[comp]->copyDataAtDirichletBoundaries(y_idx, x_idx);
    }

    // Deallocate scratch data.
    if (d_x) d_x->deallocateVectorData();

    IBTK_TIMER_STOP(t_apply);
    return;
} // apply

void SCLaplaceOperator::initializeOperatorState(const SAMRAIVectorReal<double>& in,
                                                const SAMRAIVectorReal<double>& out)
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

    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
    TBOX_ASSERT(d_ncomp == out.getNumberOfComponents());

    if (!d_hier_math_ops_external)
    {
        d_hier_math_ops =
            new HierarchyMathOps(d_object_name + "::HierarchyMathOps", d_hierarchy, d_coarsest_ln, d_finest_ln);
    }
    else
    {
        TBOX_ASSERT(d_hier_math_ops);
    }

    // Setup cached BC data.
    d_bc_helpers.resize(d_ncomp);
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        d_bc_helpers[comp] = new StaggeredPhysicalBoundaryHelper();
        d_bc_helpers[comp]->cacheBcCoefData(d_bc_coefs, d_solution_time, d_hierarchy);
    }

    // Setup the interpolation transaction information.
    d_fill_pattern = NULL;
    if (d_poisson_spec.dIsConstant())
    {
        d_fill_pattern = new SideNoCornersFillPattern(SIDEG, false, false, true);
    }
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    d_transaction_comps.clear();
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        InterpolationTransactionComponent component(d_x->getComponentDescriptorIndex(comp),
                                                    in.getComponentDescriptorIndex(comp),
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

void SCLaplaceOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_operator_state);

    // Deallocate the interpolation operators.
    d_hier_bdry_fill->deallocateOperatorState();
    d_hier_bdry_fill.reset();
    d_transaction_comps.clear();
    d_fill_pattern.reset();

    // Deallocate hierarchy math operations object.
    if (!d_hier_math_ops_external) d_hier_math_ops.reset();

    // Delete the solution and rhs vectors.
    d_x->freeVectorComponents();
    d_x.reset();

    d_b->freeVectorComponents();
    d_b.reset();

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
