// Filename: INSStaggeredVCStokesOperator.C
// Created on 15 Jun 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith, Thomas Fai
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
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// C++ STDLIB INCLUDES
#include <limits>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

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
static const std::string  U_DATA_COARSEN_TYPE =    "CUBIC_COARSEN";
static const std::string  P_DATA_COARSEN_TYPE =    "CUBIC_COARSEN";
static const std::string MU_DATA_COARSEN_TYPE = "CONSTANT_COARSEN";

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

INSStaggeredVCStokesOperator::INSStaggeredVCStokesOperator(
    Pointer<HierarchyMathOps> hier_math_ops)
    : d_is_initialized(false),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_dt(std::numeric_limits<double>::quiet_NaN()),
      d_hier_math_ops(hier_math_ops),
      d_homogeneous_bc(false),
      d_correcting_rhs(false),
      d_U_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM,static_cast<RobinBcCoefStrategy<NDIM>*>(NULL))),
      d_P_bc_coef(NULL),
      d_U_P_MU_bdry_fill_op(Pointer<HierarchyGhostCellInterpolation>(NULL)),
      d_no_fill_op(Pointer<HierarchyGhostCellInterpolation>(NULL)),
      d_x_scratch(NULL)
{
    // Setup Timers.
    IBAMR_DO_ONCE(
        t_apply                     = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredVCStokesOperator::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredVCStokesOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredVCStokesOperator::deallocateOperatorState()");
                  );
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
INSStaggeredVCStokesOperator::registerDensityVariable(
    Pointer<SideVariable<NDIM,double> > rho_var,
    const int rho_data_idx)
{
    d_rho_var = rho_var;
    d_rho_data_idx = rho_data_idx;
    return;
}// registerDensityVariable

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
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent U_scratch_component(U_scratch_idx, U_DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY);
    InterpolationTransactionComponent P_scratch_component(P_scratch_idx, P_DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY);
    InterpolationTransactionComponent mu_component(d_mu_data_idx, MU_DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY);
    std::vector<InterpolationTransactionComponent> U_P_MU_components(3);
    U_P_MU_components[0] = U_scratch_component;
    U_P_MU_components[1] = P_scratch_component;
    U_P_MU_components[2] = mu_component;
    d_U_P_MU_bdry_fill_op->resetTransactionComponents(U_P_MU_components);
    d_U_P_MU_bdry_fill_op->fillData(d_new_time);

    // Setup hierarchy data ops object for U.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<PatchHierarchy<NDIM> > hierarchy = y.getPatchHierarchy();
    Pointer<HierarchySideDataOpsReal<NDIM,double> > hier_sc_data_ops =
        hier_ops_manager->getOperationsDouble(U_out_var, hierarchy, true);

    // Compute the action of the operator:
    //      A*[u;p] = [(rho/dt)*u-0.5*div*(mu*(grad u + (grad u)^T)) + grad p; -div u].
    //
    static const bool cf_bdry_synch = true;

    d_hier_math_ops->pointwiseMultiply(
	U_out_idx, U_out_sc_var,
        d_rho_data_idx, d_rho_var,
        U_scratch_idx, U_scratch_sc_var,
        0.0, -1, NULL);

    d_hier_math_ops->vc_laplace(
        U_out_idx, U_out_sc_var,
        -0.5, 0.0, d_mu_data_idx, d_mu_var,
        U_scratch_idx, U_scratch_sc_var, d_no_fill_op, d_new_time,
        1.0/d_dt, U_out_idx, U_out_sc_var);

    d_hier_math_ops->grad(
        U_out_idx, U_out_sc_var,
        cf_bdry_synch,
        1.0, P_scratch_idx, P_scratch_cc_var, d_no_fill_op, d_new_time,
        1.0, U_out_idx, U_out_sc_var);

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

    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent U_scratch_component(d_x_scratch->getComponentDescriptorIndex(0), U_DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY);
    InterpolationTransactionComponent P_scratch_component(d_x_scratch->getComponentDescriptorIndex(1), P_DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY);
    InterpolationTransactionComponent mu_component(d_mu_data_idx, MU_DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY);

    std::vector<InterpolationTransactionComponent> U_P_MU_components(3);
    U_P_MU_components[0] = U_scratch_component;
    U_P_MU_components[1] = P_scratch_component;
    U_P_MU_components[2] = mu_component;

    d_U_P_MU_bdry_fill_op = new HierarchyGhostCellInterpolation();
    d_U_P_MU_bdry_fill_op->initializeOperatorState(U_P_MU_components, d_x_scratch->getPatchHierarchy());

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
