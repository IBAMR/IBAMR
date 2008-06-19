// Filename: INSStaggeredStokesOperator.C
// Last modified: <18.Jun.2008 17:24:23 griffith@box230.cims.nyu.edu>
// Created on 29 Apr 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

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
INSStaggeredStokesOperator::modifyRhsForInhomogeneousBc(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y)
{
    // Set y := y - A*0, i.e., shift the right-hand-side vector to account for
    // inhomogeneous boundary conditions.
    if (!d_homogeneous_bc)
    {
        d_correcting_rhs = true;

        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > x = y.cloneVector("");
        x->allocateVectorData();
        x->setToScalar(0.0);

        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > b = y.cloneVector("");
        b->allocateVectorData();
        b->setToScalar(0.0);

        apply(*x,*b);
        y.subtract(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >(&y, false), b);

        x->deallocateVectorData();
        x->freeVectorComponents();
        x.setNull();

        b->deallocateVectorData();
        b->freeVectorComponents();
        b.setNull();

        d_correcting_rhs = false;
    }
    return;
}// modifyRhsForInhomogeneousBc

void
INSStaggeredStokesOperator::apply(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y)
{
    // Initialize the operator (if necessary).
    const bool deallocate_at_completion = !d_is_initialized;
    if (!d_is_initialized) initializeOperatorState(x,y);

    // Get the vector components.
    const int U_out_idx      =            y.getComponentDescriptorIndex(0);
    const int P_out_idx      =            y.getComponentDescriptorIndex(1);
    const int U_scratch_idx  = d_x_scratch->getComponentDescriptorIndex(0);
    const int P_scratch_idx  = d_x_scratch->getComponentDescriptorIndex(1);

    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_out_var = y.getComponentVariable(0);
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_out_var = y.getComponentVariable(1);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_out_sc_var = U_out_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_out_cc_var = P_out_var;

    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_scratch_var = d_x_scratch->getComponentVariable(0);
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_scratch_var = d_x_scratch->getComponentVariable(1);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_scratch_sc_var = U_scratch_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_scratch_cc_var = P_scratch_var;

    d_x_scratch->copyVector(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >(&x,false));

    // Reset the interpolation operators and fill the data.
    typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent U_scratch_component(U_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_U_bc_coefs);
    InterpolationTransactionComponent P_scratch_component(P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);

    std::vector<InterpolationTransactionComponent> U_P_scratch_components(2);
    U_P_scratch_components[0] = U_scratch_component;
    U_P_scratch_components[1] = P_scratch_component;

    const bool homogeneous_bc = d_correcting_rhs ? d_homogeneous_bc : true;
    d_U_P_bdry_fill_op->resetTransactionComponents(U_P_scratch_components);
    d_U_P_bdry_fill_op->setHomogeneousBc(homogeneous_bc);
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
    d_hier_math_ops->laplace(
        U_out_idx, U_out_sc_var,
        d_helmholtz_spec,
        U_scratch_idx, U_scratch_sc_var,
        d_no_fill_op, d_new_time,
        1.0,
        U_out_idx, U_out_sc_var);

    // XXXX
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = y.getPatchHierarchy();
    const int coarsest_ln = y.getCoarsestLevelNumber();
    const int finest_ln = y.getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
            const SAMRAI::hier::Index<NDIM>& patch_upper = patch_box.upper();
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_out_data = patch->getPatchData(U_out_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                static const int lower = 0;
                if (pgeom->getTouchesRegularBoundary(axis,lower))
                {
                    SAMRAI::hier::Box<NDIM> lower_box = patch_box;
                    lower_box.lower()(axis) = patch_lower(axis)-1;
                    lower_box.upper()(axis) = patch_lower(axis)-1;
                    U_out_data->fillAll(0.0,lower_box);
                }
                static const int upper = 1;
                if (pgeom->getTouchesRegularBoundary(axis,upper))
                {
                    SAMRAI::hier::Box<NDIM> upper_box = patch_box;
                    upper_box.lower()(axis) = patch_upper(axis)+1;
                    upper_box.upper()(axis) = patch_upper(axis)+1;
                    U_out_data->fillAll(0.0,upper_box);
                }
            }
        }
    }
    // XXXX

    // Deallocate the operator (if necessary).
    if (deallocate_at_completion) deallocateOperatorState();
    return;
}// apply

void
INSStaggeredStokesOperator::initializeOperatorState(
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& out)
{
    if (d_is_initialized) deallocateOperatorState();

    d_x_scratch = in.cloneVector("INSStaggeredStokesOperator::x_scratch");
    d_x_scratch->allocateVectorData();

    typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent U_scratch_component(d_x_scratch->getComponentDescriptorIndex(0), DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_U_bc_coefs);
    InterpolationTransactionComponent P_scratch_component(d_x_scratch->getComponentDescriptorIndex(1), DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);

    std::vector<InterpolationTransactionComponent> U_P_scratch_components(2);
    U_P_scratch_components[0] = U_scratch_component;
    U_P_scratch_components[1] = P_scratch_component;

    d_U_P_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_U_P_bdry_fill_op->initializeOperatorState(U_P_scratch_components, d_x_scratch->getPatchHierarchy());

    d_is_initialized = true;
    return;
}// initializeOperatorState

void
INSStaggeredStokesOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    d_x_scratch->deallocateVectorData();
    d_x_scratch->freeVectorComponents();
    d_x_scratch.setNull();

    d_is_initialized = false;
    return;
}// deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::INSStaggeredStokesOperator>;

//////////////////////////////////////////////////////////////////////////////
