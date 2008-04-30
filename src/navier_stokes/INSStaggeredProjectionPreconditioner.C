// Filename: INSStaggeredProjectionPreconditioner.C
// Last modified: <29.Apr.2008 15:39:02 griffith@box230.cims.nyu.edu>
// Created on 29 Apr 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

#include "INSStaggeredProjectionPreconditioner.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif


// XXXX
#include <ibtk/SCLaplaceOperator.h>
#include <ibtk/PETScKrylovLinearSolver.h>

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

bool
INSStaggeredProjectionPreconditioner::solveSystem(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x)
{
    // Get the patch hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

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

    const int P_in_idx = x.getComponentDescriptorIndex(1);

    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_in_var = x.getComponentVariable(1);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_in_cc_var = P_in_var;

    const int U_out_idx = y.getComponentDescriptorIndex(0);
    const int P_out_idx = y.getComponentDescriptorIndex(1);

    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_out_var = y.getComponentVariable(0);
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_out_var = y.getComponentVariable(1);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_out_sc_var = U_out_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_out_cc_var = P_out_var;

    // Setup the solver vectors.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > U_scratch;
    U_scratch = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredProjectionPreconditioner::U_scratch", hierarchy, coarsest_ln, finest_ln);
    U_scratch->addComponent(U_scratch_sc_var, U_scratch_idx, d_wgt_sc_idx, d_hier_sc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > U_out;
    U_out = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredProjectionPreconditioner::U_out", hierarchy, coarsest_ln, finest_ln);
    U_out->addComponent(U_out_sc_var, U_out_idx, d_wgt_sc_idx, d_hier_sc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > P_out;
    P_out = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredProjectionPreconditioner::P_out", hierarchy, coarsest_ln, finest_ln);
    P_out->addComponent(P_out_cc_var, P_out_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    // Setup the scratch vectors.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > Grad_P_scratch;
    Grad_P_scratch = U_out->cloneVector("INSStaggeredProjectionPreconditioner::Grad_P_scratch");
    Grad_P_scratch->allocateVectorData(d_new_time);
    const int Grad_P_scratch_idx = Grad_P_scratch->getComponentDescriptorIndex(0);
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > Grad_P_scratch_sc_var = Grad_P_scratch->getComponentVariable(0);

    // Include the pressure gradient in the right-hand-side when employing a
    // pressure-increment (BCG-like) projection.
    if (d_projection_type == "pressure_increment")
    {
        typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        InterpolationTransactionComponent P_scratch_component(P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);  // XXXX
        d_P_bdry_fill_op->resetTransactionComponent(P_scratch_component);
        static const bool Grad_P_scratch_cf_bdry_synch = true;
        d_hier_math_ops->grad(
            Grad_P_scratch_idx, Grad_P_scratch_sc_var,
            Grad_P_scratch_cf_bdry_synch,
            1.0,
            P_scratch_idx, P_scratch_cc_var,
            d_P_bdry_fill_op, d_current_time);
        d_hier_sc_data_ops->subtract(U_scratch_idx, U_scratch_idx, Grad_P_scratch_idx);
    }

    // Solve for u^{*}.
    d_helmholtz_solver->solveSystem(*U_out,*U_scratch);
    SAMRAI::tbox::plog << "INSStaggeredProjectionPreconditioner::solveSystem(): Helmholtz solve number of iterations = " << d_helmholtz_solver->getNumIterations() << "\n";
    SAMRAI::tbox::plog << "INSStaggeredProjectionPreconditioner::solveSystem(): Helmholtz solve residual norm        = " << d_helmholtz_solver->getResidualNorm()  << "\n";
    if (d_helmholtz_solver->getNumIterations() == d_helmholtz_solver->getMaxIterations())
    {
        SAMRAI::tbox::pout << "INSStaggeredProjectionPreconditioner::solveSystem():"
                           <<"  WARNING: Helmholtz solver iterations == max iterations\n";
    }

    // Project the intermediate velocity u^{*} to obtain u^{n+1}.
    d_hier_projector->projectHierarchy(
        d_rho, d_dt, d_current_time+0.5*d_dt,
        d_projection_type,
        U_out_idx, U_out_sc_var,
        P_scratch_idx, P_scratch_cc_var,
        P_scratch_idx, P_scratch_cc_var,
        Grad_P_scratch_idx, Grad_P_scratch_sc_var,
        U_out_idx, U_out_sc_var);

    // Compute p^{n+1/2}.
    SAMRAI::solv::PoissonSpecifications pressure_helmholtz_spec("INSStaggeredProjectionPreconditioner::pressure_helmholtz_spec");
    pressure_helmholtz_spec.setCConstant(1.0+0.5*d_dt*d_lambda/d_rho);
    pressure_helmholtz_spec.setDConstant(   -0.5*d_dt*d_mu    /d_rho);
    d_hier_math_ops->laplace(
        P_out_idx, P_out_cc_var,
        pressure_helmholtz_spec,
        P_scratch_idx, P_scratch_cc_var,
        d_no_fill_op, d_current_time,
        (d_projection_type == "pressure_increment" ? 1.0 : 0.0),
        P_in_idx, P_in_cc_var);
    if (d_normalize_pressure)
    {
        const double P_mean = (1.0/d_volume)*d_hier_cc_data_ops->integral(P_out_idx, d_wgt_cc_idx);
        d_hier_cc_data_ops->addScalar(P_out_idx, P_out_idx, -P_mean);
    }

    // Deallocate scratch data.
    x_scratch->deallocateVectorData();
    Grad_P_scratch->deallocateVectorData();
    return true;
}// solveSystem

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::INSStaggeredProjectionPreconditioner>;

//////////////////////////////////////////////////////////////////////////////
