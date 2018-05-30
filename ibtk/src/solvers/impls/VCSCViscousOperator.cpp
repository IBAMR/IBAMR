// Filename: VCSCViscousOperator.cpp
// Created on 17 Aug 2017 by Amneet Bhalla
//
// Copyright (c) 2002-2014, Amneet Bhalla and Nishant Nangia
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

#include <ostream>
#include <stddef.h>
#include <string>
#include <vector>

#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "PoissonSpecifications.h"
#include "SAMRAIVectorReal.h"
#include "SideDataFactory.h"
#include "SideVariable.h"
#include "VariableFillPattern.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LaplaceOperator.h"
#include "ibtk/SideNoCornersFillPattern.h"
#include "ibtk/StaggeredPhysicalBoundaryHelper.h"
#include "ibtk/VCSCViscousOperator.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{

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

VCSCViscousOperator::VCSCViscousOperator(const std::string& object_name, const bool homogeneous_bc)
    : SCLaplaceOperator(object_name, homogeneous_bc)
{
    // Setup the operator to use default vector-valued boundary conditions.
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, static_cast<RobinBcCoefStrategy<NDIM>*>(NULL)));

    // Setup Timers.
    IBTK_DO_ONCE(t_apply = TimerManager::getManager()->getTimer("IBTK::VCSCViscousOperator::apply()");
                 t_initialize_operator_state =
                     TimerManager::getManager()->getTimer("IBTK::VCSCViscousOperator::initializeOperatorState()");
                 t_deallocate_operator_state =
                     TimerManager::getManager()->getTimer("IBTK::VCSCViscousOperator::deallocateOperatorState()"););

    // Set a default interpolation type.
    d_D_interp_type = IBTK::VC_HARMONIC_INTERP;
    return;
} // VCSCViscousOperator()

VCSCViscousOperator::~VCSCViscousOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
} // ~VCSCViscousOperator()

void
VCSCViscousOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    IBTK_TIMER_START(t_apply);

#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
    TBOX_ASSERT(d_bc_coefs.size() == NDIM);
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        Pointer<SideVariable<NDIM, double> > x_sc_var = x.getComponentVariable(comp);
        Pointer<SideVariable<NDIM, double> > y_sc_var = y.getComponentVariable(comp);
        if (!x_sc_var || !y_sc_var)
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                                     << "  encountered non-side centered vector components"
                                     << std::endl);
        }
        Pointer<SideDataFactory<NDIM, double> > x_factory = x_sc_var->getPatchDataFactory();
        Pointer<SideDataFactory<NDIM, double> > y_factory = y_sc_var->getPatchDataFactory();
        TBOX_ASSERT(x_factory);
        TBOX_ASSERT(y_factory);
        const unsigned int x_depth = x_factory->getDefaultDepth();
        const unsigned int y_depth = y_factory->getDefaultDepth();
        TBOX_ASSERT(x_depth == y_depth);
        if (x_depth != 1 || y_depth != 1)
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                                     << "  each vector component must have data depth == 1"
                                     << std::endl);
        }
    }
#endif

    // Allocate scratch data.
    d_x->allocateVectorData();

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

    double alpha = 1.0;
    double beta = 1.0;
    if (d_poisson_spec.cIsZero() || d_poisson_spec.cIsConstant())
    {
        beta = d_poisson_spec.cIsZero() ? 0.0 : d_poisson_spec.getCConstant();
    }

    // Compute the action of the operator.
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        Pointer<SideVariable<NDIM, double> > x_sc_var = x.getComponentVariable(comp);
        Pointer<SideVariable<NDIM, double> > y_sc_var = y.getComponentVariable(comp);
        Pointer<SideVariable<NDIM, double> > x_scratch_var = d_x->getComponentVariable(comp);
        const int x_scratch_idx = d_x->getComponentDescriptorIndex(comp);
        const int y_idx = y.getComponentDescriptorIndex(comp);
        d_hier_math_ops->vc_laplace(y_idx,
                                    y_sc_var,
                                    alpha,
                                    beta,
                                    d_poisson_spec.getDPatchDataId(),
#if (NDIM == 2)
                                    Pointer<NodeVariable<NDIM, double> >(NULL),
#elif (NDIM == 3)
                                    Pointer<EdgeVariable<NDIM, double> >(NULL),
#endif
                                    x_scratch_idx,
                                    x_sc_var,
                                    Pointer<HierarchyGhostCellInterpolation>(NULL),
                                    d_solution_time,
                                    d_D_interp_type,
                                    d_poisson_spec.cIsVariable() ? d_poisson_spec.getCPatchDataId() : -1);
        const int x_idx = x.getComponentDescriptorIndex(comp);
        d_bc_helpers[comp]->copyDataAtDirichletBoundaries(y_idx, x_idx);
    }

    // Deallocate scratch data.
    d_x->deallocateVectorData();

    IBTK_TIMER_STOP(t_apply);
    return;
} // apply

void
VCSCViscousOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
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

    // Setup cached BC data.
    d_bc_helpers.resize(d_ncomp);
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        d_bc_helpers[comp] = new StaggeredPhysicalBoundaryHelper();
        d_bc_helpers[comp]->cacheBcCoefData(d_bc_coefs, d_solution_time, d_hierarchy);
    }

    // Setup the interpolation transaction information.
    d_fill_pattern = NULL;
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

void
VCSCViscousOperator::deallocateOperatorState()
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

void
VCSCViscousOperator::setDPatchDataInterpolationType(const IBTK::VCInterpType D_interp_type)
{
    d_D_interp_type = D_interp_type;
    return;
} // setDPatchDataInterpolationType

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
