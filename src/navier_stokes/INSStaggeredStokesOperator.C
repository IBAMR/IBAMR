// Filename: INSStaggeredLinearOperator.C
// Last modified: <29.Apr.2008 14:59:50 griffith@box230.cims.nyu.edu>
// Created on 29 Apr 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

#include "INSStaggeredLinearOperator.h"

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

namespace
{
// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_COARSEN_TYPE = "CONSERVATIVE_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
INSStaggeredLinearOperator::apply(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y)
{
    // Get the vector components.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > x_scratch = x.cloneVector("INSStaggeredLinearOperator::scratch");
    x_scratch->allocateVectorData(d_new_time);
    x_scratch->copyVector(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >(&x,false));

    const int U_scratch_idx = x_scratch->getComponentDescriptorIndex(0);
    const int P_scratch_idx = x_scratch->getComponentDescriptorIndex(1);

    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_scratch_var = x_scratch->getComponentVariable(0);
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_scratch_var = x_scratch->getComponentVariable(1);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_scratch_sc_var = U_scratch_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_scratch_cc_var = P_scratch_var;

    const int U_out_idx = y.getComponentDescriptorIndex(0);
    const int P_out_idx = y.getComponentDescriptorIndex(1);

    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_out_var = y.getComponentVariable(0);
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_out_var = y.getComponentVariable(1);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_out_sc_var = U_out_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_out_cc_var = P_out_var;

    // Reset the interpolation operators and fill the data.
    typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent U_component(U_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);  // XXXX
    InterpolationTransactionComponent P_component(P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);

    std::vector<InterpolationTransactionComponent> transaction_comps(2);
    transaction_comps[0] = U_component;
    transaction_comps[1] = P_component;

    d_U_P_bdry_fill_op->resetTransactionComponents(transaction_comps);
    d_U_P_bdry_fill_op->fillData(d_new_time);

    // Compute the action of the operator:
    //      A*[u;p] = [((rho/dt)*I-0.5*mu*L)*u + grad p; -div u].
    bool cf_bdry_synch;
    cf_bdry_synch = true;
    d_hier_math_ops->grad(
        U_out_idx, U_out_sc_var,
        cf_bdry_synch,
        1.0, P_scratch_idx, P_scratch_cc_var, d_no_fill_op, d_new_time);
    cf_bdry_synch = false;
    d_hier_math_ops->div(
        P_out_idx, P_out_cc_var,
        -1.0, U_scratch_idx, U_scratch_sc_var, d_no_fill_op, d_new_time,
        cf_bdry_synch);

    SAMRAI::solv::PoissonSpecifications helmholtz_spec("INSStaggeredLinearOperator::helmholtz_spec");
    helmholtz_spec.setCConstant((d_rho/d_dt)+0.5*d_lambda);
    helmholtz_spec.setDConstant(            -0.5*d_mu    );
    d_hier_math_ops->laplace(
        U_out_idx, U_out_sc_var,
        helmholtz_spec,
        U_scratch_idx, U_scratch_sc_var,
        d_no_fill_op, d_new_time,
        1.0,
        U_out_idx, U_out_sc_var);

    // Deallocate scratch data.
    x_scratch->deallocateVectorData();
    return;
}// apply

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::INSStaggeredLinearOperator>;

//////////////////////////////////////////////////////////////////////////////
