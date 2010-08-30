// Filename: INSStaggeredVCStokesOperator.C
// Created on 15 Jun 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "INSStaggeredVCStokesOperator.h"

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
static Pointer<Timer> t_apply;
static Pointer<Timer> t_initialize_operator_state;
static Pointer<Timer> t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredVCStokesOperator::INSStaggeredVCStokesOperator(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
    Pointer<INSStaggeredPhysicalBoundaryHelper> U_bc_helper,
    RobinBcCoefStrategy<NDIM>* P_bc_coef,
    Pointer<HierarchyMathOps> hier_math_ops)
    : d_is_initialized(false),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_dt(std::numeric_limits<double>::quiet_NaN()),
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
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_apply                     = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredVCStokesOperator::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredVCStokesOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredVCStokesOperator::deallocateOperatorState()");
        timers_need_init = false;
    }
    return;
}// INSStaggeredVCStokesOperator

INSStaggeredVCStokesOperator::~INSStaggeredVCStokesOperator()
{
    deallocateOperatorState();
    return;
}// ~INSStaggeredVCStokesOperator

void
INSStaggeredVCStokesOperator::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
INSStaggeredVCStokesOperator::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_dt = d_new_time-d_current_time;
    return;
}// setTimeInterval

void
INSStaggeredVCStokesOperator::registerViscosityVariable(
    Pointer<NodeVariable<NDIM,double> > mu_var,
    const int mu_data_idx)
{
    d_mu_var = mu_var;
    d_mu_data_idx = mu_data_idx;
    return;
}// registerViscosityVariable

void
INSStaggeredVCStokesOperator::modifyRhsForInhomogeneousBc(
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
INSStaggeredVCStokesOperator::apply(
    const bool homogeneous_bc,
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    t_apply->start();

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
    Pointer<VariableFillPattern<NDIM> > sc_fill_pattern = new SideNoCornersFillPattern(SIDEG);
    Pointer<VariableFillPattern<NDIM> > cc_fill_pattern = new CellNoCornersFillPattern(CELLG);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent U_scratch_component(U_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_U_bc_coefs, sc_fill_pattern);
    InterpolationTransactionComponent P_scratch_component(P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_P_bc_coef , cc_fill_pattern);
    std::vector<InterpolationTransactionComponent> U_P_components(2);
    U_P_components[0] = U_scratch_component;
    U_P_components[1] = P_scratch_component;
    INSStaggeredPressureBcCoef* P_bc_coef = dynamic_cast<INSStaggeredPressureBcCoef*>(d_P_bc_coef);
    P_bc_coef->setVelocityNewPatchDataIndex(U_scratch_idx);
    d_U_P_bdry_fill_op->setHomogeneousBc(homogeneous_bc);
    d_U_P_bdry_fill_op->resetTransactionComponents(U_P_components);
    d_U_P_bdry_fill_op->fillData(d_new_time);

    // Compute the action of the operator:
    //      A*[u;p] = [((rho/dt)*I-0.5*mu*L)*u + grad p; -div u].
    //
    // Thomas: Please update the forgoing comment to reflect the form of the
    // variable-density and variable-viscosity operator.
    static const bool cf_bdry_synch = true;
    d_hier_math_ops->grad(
        U_out_idx, U_out_sc_var,
        cf_bdry_synch,
        1.0, P_scratch_idx, P_scratch_cc_var, d_no_fill_op, d_new_time);

//  This was from the uniform-density and uniform-viscosity implementation:
//
//     d_hier_math_ops->laplace(
//         U_out_idx, U_out_sc_var,
//         d_helmholtz_spec,
//         U_scratch_idx, U_scratch_sc_var, d_no_fill_op, d_new_time,
//         1.0,
//         U_out_idx, U_out_sc_var);
//
//  Thomas: the foregoing needs to be replaced by the appropriate
//  variable-density and variable-viscosity code.

    d_U_bc_helper->zeroValuesAtDirichletBoundaries(U_out_idx);

    d_hier_math_ops->div(
        P_out_idx, P_out_cc_var,
        -1.0, U_scratch_idx, U_scratch_sc_var, d_no_fill_op, d_new_time,
        cf_bdry_synch);

    t_apply->stop();
    return;
}// apply

void
INSStaggeredVCStokesOperator::apply(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    apply(d_correcting_rhs ? d_homogeneous_bc : true,x,y);
    return;
}// apply

void
INSStaggeredVCStokesOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAIVectorReal<NDIM,double>& out)
{
    t_initialize_operator_state->start();

    if (d_is_initialized) deallocateOperatorState();

    d_x_scratch = in.cloneVector("INSStaggeredVCStokesOperator::x_scratch");
    d_x_scratch->allocateVectorData();

    Pointer<VariableFillPattern<NDIM> > sc_fill_pattern = new SideNoCornersFillPattern(SIDEG);
    Pointer<VariableFillPattern<NDIM> > cc_fill_pattern = new CellNoCornersFillPattern(CELLG);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent U_scratch_component(d_x_scratch->getComponentDescriptorIndex(0), DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_U_bc_coefs, sc_fill_pattern);
    InterpolationTransactionComponent P_scratch_component(d_x_scratch->getComponentDescriptorIndex(1), DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_P_bc_coef , cc_fill_pattern);

    std::vector<InterpolationTransactionComponent> U_P_components(2);
    U_P_components[0] = U_scratch_component;
    U_P_components[1] = P_scratch_component;

    d_U_P_bdry_fill_op = new HierarchyGhostCellInterpolation();
    d_U_P_bdry_fill_op->initializeOperatorState(U_P_components, d_x_scratch->getPatchHierarchy());

    d_is_initialized = true;

    t_initialize_operator_state->stop();
    return;
}// initializeOperatorState

void
INSStaggeredVCStokesOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    t_deallocate_operator_state->start();

    d_x_scratch->freeVectorComponents();
    d_x_scratch.setNull();

    d_is_initialized = false;

    t_deallocate_operator_state->stop();
    return;
}// deallocateOperatorState

void
INSStaggeredVCStokesOperator::enableLogging(
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
template class Pointer<IBAMR::INSStaggeredVCStokesOperator>;

//////////////////////////////////////////////////////////////////////////////
