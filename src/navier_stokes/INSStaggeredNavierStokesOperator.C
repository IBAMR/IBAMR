// Filename: INSStaggeredNavierStokesOperator.C
// Last modified: <17.Jul.2008 19:07:16 griffith@box230.cims.nyu.edu>
// Created on 08 May 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

#include "INSStaggeredNavierStokesOperator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredNavierStokesOperator::INSStaggeredNavierStokesOperator(
    const double rho,
    const double mu,
    const double lambda,
    SAMRAI::tbox::Pointer<INSStaggeredStokesOperator> stokes_op,
    SAMRAI::tbox::Pointer<INSStaggeredConvectiveOperator> convective_op,
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> > hier_sc_data_ops)
    : d_is_initialized(false),
      d_U_current_idx(-1),
      d_rho(rho),
      d_mu(mu),
      d_lambda(lambda),
      d_stokes_op(stokes_op),
      d_convective_op(convective_op),
      d_hier_sc_data_ops(hier_sc_data_ops),
      d_x_scratch(NULL),
      d_y_scratch(NULL)
{
    // intentionally blank
    return;
}// INSStaggeredNavierStokesOperator

INSStaggeredNavierStokesOperator::~INSStaggeredNavierStokesOperator()
{
    deallocateOperatorState();
    return;
}// ~INSStaggeredNavierStokesOperator

void
INSStaggeredNavierStokesOperator::setVelocityCurrentPatchDataIndex(
    const int U_current_idx)
{
    d_U_current_idx = U_current_idx;
    return;
}// setVelocityCurrentPatchDataIndex

void
INSStaggeredNavierStokesOperator::apply(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y)
{
    // Initialize the operator (if necessary).
    const bool deallocate_at_completion = !d_is_initialized;
    if (!d_is_initialized) initializeOperatorState(x,y);

    // Get the vector components.
    const int U_in_idx          =            x.getComponentDescriptorIndex(0);
    const int U_out_idx         =            y.getComponentDescriptorIndex(0);
    const int U_in_scratch_idx  = d_x_scratch->getComponentDescriptorIndex(0);
    const int U_out_scratch_idx = d_y_scratch->getComponentDescriptorIndex(0);

    d_x_scratch->copyVector(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >(&x,false));
    d_y_scratch->copyVector(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >(&y,false));

    // Set U_in_scratch := 0.5*(U(n) + U(n+1)) = U(n+1/2).
    d_hier_sc_data_ops->linearSum(U_in_scratch_idx, 0.5, d_U_current_idx, 0.5, U_in_idx);

    // Compute the action of the Stokes operator (the linear part of the
    // problem).
    d_stokes_op->apply(x,y);

    // Compute the action of the convective operator (the nonlinear part of the
    // problem)
    d_convective_op->apply(*d_x_scratch, *d_y_scratch);

    // Compute the final result.
    d_hier_sc_data_ops->axpy(U_out_idx, d_rho, U_out_scratch_idx, U_out_idx);

    // Deallocate the operator (if necessary).
    if (deallocate_at_completion) deallocateOperatorState();
    return;
}// apply

void
INSStaggeredNavierStokesOperator::initializeOperatorState(
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& out)
{
    if (d_is_initialized) deallocateOperatorState();

    d_x_scratch =  in.cloneVector("INSStaggeredNavierStokesOperator::x_scratch");
    d_y_scratch = out.cloneVector("INSStaggeredNavierStokesOperator::y_scratch");

    d_x_scratch->allocateVectorData();
    d_y_scratch->allocateVectorData();

    d_stokes_op    ->initializeOperatorState(*d_x_scratch,*d_y_scratch);
    d_convective_op->initializeOperatorState(*d_x_scratch,*d_y_scratch);

    d_is_initialized = true;
    return;
}// initializeOperatorState

void
INSStaggeredNavierStokesOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    d_stokes_op    ->deallocateOperatorState();
    d_convective_op->deallocateOperatorState();

    d_x_scratch->freeVectorComponents();
    d_y_scratch->freeVectorComponents();

    d_x_scratch.setNull();
    d_y_scratch.setNull();

    d_is_initialized = false;
    return;
}// deallocateOperatorState

void
INSStaggeredNavierStokesOperator::enableLogging(
    bool enabled)
{
    d_stokes_op    ->enableLogging(enabled);
    d_convective_op->enableLogging(enabled);
    return;
}// enableLogging

void
INSStaggeredNavierStokesOperator::printClassData(
    std::ostream& os) const
{
    os << "INSStaggeredNavierStokesOperator\n";
    os << "d_stokes_op->printClassData():\n";
    d_stokes_op->printClassData(os);
    os << "d_convective_op->printClassData():\n";
    d_convective_op->printClassData(os);
    return;
}// printClassData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::INSStaggeredNavierStokesOperator>;

//////////////////////////////////////////////////////////////////////////////
