// Filename: INSStaggeredStokesOperator.C
// Created on 29 Apr 2008 by Boyce Griffith
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

#include "INSStaggeredStokesOperator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/INSStaggeredPressureBcCoef.h>
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/SideNoCornersFillPattern.h>

// C++ STDLIB INCLUDES
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
static const int SIDEG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

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

INSStaggeredStokesOperator::INSStaggeredStokesOperator(
    const INSCoefs& problem_coefs,
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& U_bc_coefs,
    Pointer<INSStaggeredPhysicalBoundaryHelper> U_bc_helper,
    RobinBcCoefStrategy<NDIM>* P_bc_coef,
    Pointer<HierarchyMathOps> hier_math_ops)
    : d_is_initialized(false),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_dt(std::numeric_limits<double>::quiet_NaN()),
      d_problem_coefs(problem_coefs),
      d_helmholtz_spec("INSStaggeredStokesOperator::helmholtz_spec"),
      d_hier_math_ops(hier_math_ops),
      d_homogeneous_bc(false),
      d_correcting_rhs(false),
      d_U_bc_coefs(U_bc_coefs),
      d_U_bc_helper(U_bc_helper),
      d_P_bc_coef(P_bc_coef),
      d_U_P_bdry_fill_op(Pointer<HierarchyGhostCellInterpolation>(NULL)),
      d_no_fill_op(Pointer<HierarchyGhostCellInterpolation>(NULL)),
      d_x_scratch(NULL)
{
    // Setup Timers.
    IBAMR_DO_ONCE(
        t_apply                     = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredStokesOperator::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredStokesOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredStokesOperator::deallocateOperatorState()");
                  );
    return;
}// INSStaggeredStokesOperator

INSStaggeredStokesOperator::~INSStaggeredStokesOperator()
{
    deallocateOperatorState();
    return;
}// ~INSStaggeredStokesOperator

void
INSStaggeredStokesOperator::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
INSStaggeredStokesOperator::setTimeInterval(
    const double current_time,
    const double new_time)
{
    const double rho    = d_problem_coefs.getRho();
    const double mu     = d_problem_coefs.getMu();
    const double lambda = d_problem_coefs.getLambda();
    d_current_time = current_time;
    d_new_time = new_time;
    d_dt = d_new_time-d_current_time;
    d_helmholtz_spec.setCConstant((rho/d_dt)+0.5*lambda);
    d_helmholtz_spec.setDConstant(          -0.5*mu    );
    return;
}// setTimeInterval

void
INSStaggeredStokesOperator::modifyRhsForInhomogeneousBc(
    SAMRAIVectorReal<NDIM,double>& y)
{
    // Set y := y - A*0, i.e., shift the right-hand-side vector to account for
    // inhomogeneous boundary conditions.
    if (!d_homogeneous_bc)
    {
        d_correcting_rhs = true;

        Pointer<SAMRAIVectorReal<NDIM,double> > x = y.cloneVector("");
        x->allocateVectorData();
        x->setToScalar(0.0);

        Pointer<SAMRAIVectorReal<NDIM,double> > b = y.cloneVector("");
        b->allocateVectorData();
        b->setToScalar(0.0);

        apply(*x,*b);
        y.subtract(Pointer<SAMRAIVectorReal<NDIM,double> >(&y, false), b);

        x->freeVectorComponents();
        x.setNull();

        b->freeVectorComponents();
        b.setNull();

        d_correcting_rhs = false;
    }
    return;
}// modifyRhsForInhomogeneousBc

void
INSStaggeredStokesOperator::apply(
    const bool homogeneous_bc,
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    IBAMR_TIMER_START(t_apply);

    // Get the vector components.
//  const int U_in_idx       =            x.getComponentDescriptorIndex(0);
//  const int P_in_idx       =            x.getComponentDescriptorIndex(1);
    const int U_out_idx      =            y.getComponentDescriptorIndex(0);
    const int P_out_idx      =            y.getComponentDescriptorIndex(1);
    const int U_scratch_idx  = d_x_scratch->getComponentDescriptorIndex(0);
    const int P_scratch_idx  = d_x_scratch->getComponentDescriptorIndex(1);

    const Pointer<Variable<NDIM> >& U_out_var = y.getComponentVariable(0);
    const Pointer<Variable<NDIM> >& P_out_var = y.getComponentVariable(1);

    Pointer<SideVariable<NDIM,double> > U_out_sc_var = U_out_var;
    Pointer<CellVariable<NDIM,double> > P_out_cc_var = P_out_var;

    const Pointer<Variable<NDIM> >& U_scratch_var = d_x_scratch->getComponentVariable(0);
    const Pointer<Variable<NDIM> >& P_scratch_var = d_x_scratch->getComponentVariable(1);

    Pointer<SideVariable<NDIM,double> > U_scratch_sc_var = U_scratch_var;
    Pointer<CellVariable<NDIM,double> > P_scratch_cc_var = P_scratch_var;

    d_x_scratch->copyVector(Pointer<SAMRAIVectorReal<NDIM,double> >(&x,false));

    // Reset the interpolation operators and fill the data.
    Pointer<VariableFillPattern<NDIM> > sc_fill_pattern = new SideNoCornersFillPattern(SIDEG, false, false, true);
    Pointer<VariableFillPattern<NDIM> > cc_fill_pattern = new CellNoCornersFillPattern(CELLG, false, false, true);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent U_scratch_component(U_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_U_bc_coefs, sc_fill_pattern);
    InterpolationTransactionComponent P_scratch_component(P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_P_bc_coef , cc_fill_pattern);
    std::vector<InterpolationTransactionComponent> U_P_components(2);
    U_P_components[0] = U_scratch_component;
    U_P_components[1] = P_scratch_component;
    INSStaggeredPressureBcCoef* P_bc_coef = dynamic_cast<INSStaggeredPressureBcCoef*>(d_P_bc_coef);
    if (P_bc_coef != NULL) P_bc_coef->setVelocityNewPatchDataIndex(U_scratch_idx);
    d_U_P_bdry_fill_op->setHomogeneousBc(homogeneous_bc);
    d_U_P_bdry_fill_op->resetTransactionComponents(U_P_components);
    d_U_P_bdry_fill_op->fillData(d_new_time);

    // Compute the action of the operator:
    //      A*[u;p] = [((rho/dt)*I-0.5*mu*L)*u + grad p; -div u].
    static const bool cf_bdry_synch = true;
    d_hier_math_ops->grad(
        U_out_idx, U_out_sc_var,
        cf_bdry_synch,
        1.0, P_scratch_idx, P_scratch_cc_var, d_no_fill_op, d_new_time);
    d_hier_math_ops->laplace(
        U_out_idx, U_out_sc_var,
        d_helmholtz_spec,
        U_scratch_idx, U_scratch_sc_var, d_no_fill_op, d_new_time,
        1.0,
        U_out_idx, U_out_sc_var);
    if (!d_U_bc_helper.isNull()) d_U_bc_helper->zeroValuesAtDirichletBoundaries(U_out_idx);

    d_hier_math_ops->div(
        P_out_idx, P_out_cc_var,
        -1.0, U_scratch_idx, U_scratch_sc_var, d_no_fill_op, d_new_time,
        cf_bdry_synch);

    IBAMR_TIMER_STOP(t_apply);
    return;
}// apply

void
INSStaggeredStokesOperator::apply(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    apply(d_correcting_rhs ? d_homogeneous_bc : true,x,y);
    return;
}// apply

void
INSStaggeredStokesOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAIVectorReal<NDIM,double>& out)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    if (d_is_initialized) deallocateOperatorState();

    d_x_scratch = in.cloneVector("INSStaggeredStokesOperator::x_scratch");
    d_x_scratch->allocateVectorData();

    Pointer<VariableFillPattern<NDIM> > sc_fill_pattern = new SideNoCornersFillPattern(SIDEG, false, false, true);
    Pointer<VariableFillPattern<NDIM> > cc_fill_pattern = new CellNoCornersFillPattern(CELLG, false, false, true);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent U_scratch_component(d_x_scratch->getComponentDescriptorIndex(0), DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_U_bc_coefs, sc_fill_pattern);
    InterpolationTransactionComponent P_scratch_component(d_x_scratch->getComponentDescriptorIndex(1), DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_P_bc_coef , cc_fill_pattern);

    std::vector<InterpolationTransactionComponent> U_P_components(2);
    U_P_components[0] = U_scratch_component;
    U_P_components[1] = P_scratch_component;

    d_U_P_bdry_fill_op = new HierarchyGhostCellInterpolation();
    d_U_P_bdry_fill_op->initializeOperatorState(U_P_components, d_x_scratch->getPatchHierarchy());

    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
}// initializeOperatorState

void
INSStaggeredStokesOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

    d_x_scratch->freeVectorComponents();
    d_x_scratch.setNull();

    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
}// deallocateOperatorState

void
INSStaggeredStokesOperator::enableLogging(
    bool enabled)
{
    // intentionally blank
    return;
}// enableLogging

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::INSStaggeredStokesOperator>;

//////////////////////////////////////////////////////////////////////////////
