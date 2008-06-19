// Filename: INSStaggeredProjectionPreconditioner.C
// Last modified: <18.Jun.2008 18:58:51 griffith@box230.cims.nyu.edu>
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
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b)
{
    // Initialize the solver (if necessary).
    const bool deallocate_at_completion = !d_is_initialized;
    if (!d_is_initialized) initializeSolverState(x,b);

    // Get the vector components.
    const int U_scratch_idx      = d_b_scratch->getComponentDescriptorIndex(0);
    const int Grad_P_scratch_idx = d_x_scratch->getComponentDescriptorIndex(0);
    const int P_scratch_idx      = d_b_scratch->getComponentDescriptorIndex(1);

    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_scratch_var      = d_b_scratch->getComponentVariable(0);
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& Grad_P_scratch_var = d_x_scratch->getComponentVariable(0);
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_scratch_var      = d_b_scratch->getComponentVariable(1);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_scratch_sc_var      = U_scratch_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > Grad_P_scratch_sc_var = Grad_P_scratch_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_scratch_cc_var      = P_scratch_var;

    const int U_in_idx = b.getComponentDescriptorIndex(0);
    const int P_in_idx = b.getComponentDescriptorIndex(1);

    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_in_var = b.getComponentVariable(0);
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_in_var = b.getComponentVariable(1);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_in_sc_var = U_in_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_in_cc_var = P_in_var;

    (void) U_in_idx;
    (void) U_in_var;
    (void) U_in_sc_var;

    const int U_out_idx = x.getComponentDescriptorIndex(0);
    const int P_out_idx = x.getComponentDescriptorIndex(1);

    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_out_var = x.getComponentVariable(0);
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_out_var = x.getComponentVariable(1);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_out_sc_var = U_out_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_out_cc_var = P_out_var;

    d_b_scratch->copyVector(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >(&b,false));

    // Setup the component solver vectors.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > U_scratch;
    U_scratch = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredProjectionPreconditioner::U_scratch", d_hierarchy, d_coarsest_ln, d_finest_ln);
    U_scratch->addComponent(U_scratch_sc_var, U_scratch_idx, d_wgt_sc_idx, d_hier_sc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > U_out;
    U_out = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredProjectionPreconditioner::U_out", d_hierarchy, d_coarsest_ln, d_finest_ln);
    U_out->addComponent(U_out_sc_var, U_out_idx, d_wgt_sc_idx, d_hier_sc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > P_out;
    P_out = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredProjectionPreconditioner::P_out", d_hierarchy, d_coarsest_ln, d_finest_ln);
    P_out->addComponent(P_out_cc_var, P_out_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > Grad_P_scratch;
    Grad_P_scratch = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredProjectionPreconditioner::Grad_P_scratch", d_hierarchy, d_coarsest_ln, d_finest_ln);
    Grad_P_scratch->addComponent(Grad_P_scratch_sc_var, Grad_P_scratch_idx, d_wgt_sc_idx, d_hier_sc_data_ops);

    // Include the pressure gradient in the right-hand-side when employing a
    // pressure-increment (BCG-like) projection.
    if (d_projection_type == "pressure_increment")
    {
        typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        InterpolationTransactionComponent P_scratch_component(P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);
        d_P_bdry_fill_op->resetTransactionComponent(P_scratch_component);
        d_P_bdry_fill_op->fillData(d_current_time);

        static const bool Grad_P_scratch_cf_bdry_synch = true;
        d_hier_math_ops->grad(
            Grad_P_scratch_idx, Grad_P_scratch_sc_var,
            Grad_P_scratch_cf_bdry_synch,
            1.0,
            P_scratch_idx, P_scratch_cc_var,
            d_no_fill_op, d_current_time);
        d_hier_sc_data_ops->subtract(U_scratch_idx, U_scratch_idx, Grad_P_scratch_idx);
    }

    // Solve for u^{*}.
    TBOX_ASSERT(U_scratch->getComponentDescriptorIndex(0) == d_b_scratch->getComponentDescriptorIndex(0));
    d_U_bdry_fill_op->setHomogeneousBc(true);
    d_U_bdry_fill_op->fillData(d_new_time);
    U_out->copyVector(U_scratch);
    d_helmholtz_solver->setInitialGuessNonzero(true);
    d_helmholtz_solver->solveSystem(*U_out,*U_scratch);

    if (d_do_log) SAMRAI::tbox::plog << "INSStaggeredProjectionPreconditioner::solveSystem(): Helmholtz solve number of iterations = " << d_helmholtz_solver->getNumIterations() << "\n";
    if (d_do_log) SAMRAI::tbox::plog << "INSStaggeredProjectionPreconditioner::solveSystem(): Helmholtz solve residual norm        = " << d_helmholtz_solver->getResidualNorm()  << "\n";
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
    d_hier_math_ops->laplace(
        P_out_idx, P_out_cc_var,
        d_pressure_helmholtz_spec,
        P_scratch_idx, P_scratch_cc_var,
        d_no_fill_op, d_current_time,
        (d_projection_type == "pressure_increment" ? 1.0 : 0.0),
        P_in_idx, P_in_cc_var);
    if (d_normalize_pressure)
    {
        const double P_mean = (1.0/d_volume)*d_hier_cc_data_ops->integral(P_out_idx, d_wgt_cc_idx);
        d_hier_cc_data_ops->addScalar(P_out_idx, P_out_idx, -P_mean);
    }

    // Deallocate the solver (if necessary).
    if (deallocate_at_completion) deallocateSolverState();
    return true;
}// solveSystem

void
INSStaggeredProjectionPreconditioner::initializeSolverState(
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b)
{
    if (d_is_initialized) deallocateSolverState();

    d_x_scratch = x.cloneVector("INSStaggeredProjectionPreconditioner::x_scratch");
    d_b_scratch = b.cloneVector("INSStaggeredProjectionPreconditioner::b_scratch");

    d_x_scratch->allocateVectorData();
    d_b_scratch->allocateVectorData();

    typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent U_scratch_component(d_b_scratch->getComponentDescriptorIndex(0), DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_U_bc_coefs);
    d_U_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_U_bdry_fill_op->initializeOperatorState(U_scratch_component, d_b_scratch->getPatchHierarchy());

    InterpolationTransactionComponent P_scratch_component(d_b_scratch->getComponentDescriptorIndex(1), DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);
    d_P_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_P_bdry_fill_op->initializeOperatorState(P_scratch_component, d_b_scratch->getPatchHierarchy());

    d_hierarchy = x.getPatchHierarchy();
    d_coarsest_ln = x.getCoarsestLevelNumber();
    d_finest_ln = x.getFinestLevelNumber();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_hierarchy == b.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == b.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == b.getFinestLevelNumber());
#endif

    d_wgt_cc_var = d_hier_math_ops->getCellWeightVariable();
    d_wgt_sc_var = d_hier_math_ops->getSideWeightVariable();
    d_wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    d_wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
    d_volume = d_hier_math_ops->getVolumeOfPhysicalDomain();

    d_is_initialized = true;
    return;
}// initializeSolverState

void
INSStaggeredProjectionPreconditioner::deallocateSolverState()
{
    if (!d_is_initialized) return;

    d_x_scratch->deallocateVectorData();
    d_b_scratch->deallocateVectorData();

    d_x_scratch->freeVectorComponents();
    d_b_scratch->freeVectorComponents();

    d_x_scratch.setNull();
    d_b_scratch.setNull();

    d_is_initialized = false;
    return;
}// deallocateSolverState

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::INSStaggeredProjectionPreconditioner>;

//////////////////////////////////////////////////////////////////////////////
