// Filename: INSStaggeredHierarchyIntegrator.C
// Last modified: <23.Apr.2008 18:16:05 boyce@trasnaform2.local>
// Created on 20 Mar 2008 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "INSStaggeredHierarchyIntegrator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

#define COUPLED_SOLVER 1

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PETScSAMRAIVectorReal.h>
#include <ibtk/SCLaplaceOperator.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <HierarchyDataOpsManager.h>
#include <PatchCellDataOpsReal.h>
#include <tbox/MathUtilities.h>
#include <tbox/RestartManager.h>
#include <tbox/TimerManager.h>

// C++ STDLIB INCLUDES
#include <iterator>
#include <limits>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define ADVECT_DERIVATIVE_F77 F77_FUNC_(advect_derivative2d, ADVECT_DERIVATIVE2D)
#define GODUNOV_PREDICT_WITH_SOURCE_F77 F77_FUNC_(godunov_predict_with_source2d, GODUNOV_PREDICT_WITH_SOURCE2D)
#define NAVIER_STOKES_INTERP_COMPS_F77 F77_FUNC_(navier_stokes_interp_comps2d, NAVIER_STOKES_INTERP_COMPS2D)
#define NAVIER_STOKES_RESET_ADV_VELOCITY_F77 F77_FUNC_(navier_stokes_reset_adv_velocity2d, NAVIER_STOKES_RESET_ADV_VELOCITY2D)
#define NAVIER_STOKES_SC_STABLEDT_F77 F77_FUNC_(navier_stokes_sc_stabledt2d, NAVIER_STOKES_SC_STABLEDT2D)
#endif

#if (NDIM == 3)
#define ADVECT_DERIVATIVE_F77 F77_FUNC_(advect_derivative3d, ADVECT_DERIVATIVE3D)
#define NAVIER_STOKES_INTERP_COMPS_F77 F77_FUNC_(navier_stokes_interp_comps_3d, NAVIER_STOKES_INTERP_COMPS_3D)
#define NAVIER_STOKES_RESET_ADV_VELOCITY_F77 F77_FUNC_(navier_stokes_reset_adv_velocity_3d, NAVIER_STOKES_RESET_ADV_VELOCITY_3D)
#define NAVIER_STOKES_SC_STABLEDT_F77 F77_FUNC_(navier_stokes_sc_stabledt3d, NAVIER_STOKES_SC_STABLEDT3D)
#endif

extern "C"
{
    void
    ADVECT_DERIVATIVE_F77(
        const double*,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const double* , const double* ,
        const int& , const int& ,
        double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
        double*
#endif
                          );

    void
    GODUNOV_PREDICT_WITH_SOURCE_F77(
        const double* , const double& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , double* , double* , double* , double* ,
        const double* , double* ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        double* , double* ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , double* , double* ,
        const double* , double* , double* ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        double* , double* , double* ,
        double* , double* , double*
#endif
                                       );

    void
    NAVIER_STOKES_INTERP_COMPS_F77(
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        double* , double* ,
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        double* , double*
#endif
#if (NDIM == 3)
#endif
                                    );

    void
    NAVIER_STOKES_RESET_ADV_VELOCITY_F77(
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        double* , double* ,
        const int& , const int& ,
        const double* , const double* ,
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        double* , double* ,
        const int& , const int& ,
        const double* , const double*
#endif
#if (NDIM == 3)
#endif
                                          );

    void
    NAVIER_STOKES_SC_STABLEDT_F77(
        const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
#endif
        double&
                                  );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_hierarchy_integrator;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_advance_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_get_stable_timestep;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_regrid_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_integrate_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_synchronize_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_synchronize_new_levels;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_time_dependent_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_data_to_preadvance_state;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_hierarchy_configuration;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_apply_gradient_detector;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_put_to_database;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
static const int SIDEG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// NOTE: The number of ghost cells required by the Godunov advection scheme
// depends on the number of cycles we perform and the order of the
// reconstruction.  These values were chosen to work with xsPPM7 (the modified
// piecewise parabolic method of Rider, Greenough, and Kamm).
static const int NUM_GODUNOV_CYCLES = 2;
static const int GADVECTG = NUM_GODUNOV_CYCLES+3;

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_COARSEN_TYPE = "CONSERVATIVE_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Version of INSStaggeredHierarchyIntegrator restart file data.
static const int INS_STAGGERED_HIERARCHY_INTEGRATOR_VERSION = 1;
}

namespace
{

// Linear operator used with the Krylov solver.
class INSOperator
    : public IBTK::LinearOperator
{
public:
    /*!
     * XXXX
     */
    INSOperator(
        const double rho,
        const double mu,
        const double lambda,
        const double dt,
        const double apply_time,
        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
        SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> U_P_bdry_fill_op)
        : d_rho(rho),
          d_mu(mu),
          d_lambda(lambda),
          d_dt(dt),
          d_apply_time(apply_time),
          d_hier_math_ops(hier_math_ops),
          d_U_P_bdry_fill_op(U_P_bdry_fill_op),
          d_no_fill_op(SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation>(NULL))
        {
            // intentionally blank
            return;
        }// INSOperator

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~INSOperator()
        {
            // intentionally blank
            return;
        }// ~INSOperator

    /*!
     * \name Linear operator functionality.
     */
    //\{

    /*!
     * \brief Compute y=Ax.
     *
     * Before calling this function, the form of the vectors x and y should be
     * set properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in these vectors should be allocated.
     * Thus, the user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must
     * have enough ghost cells for computation of Ax).
     *
     * \note In general, the vectors x and y \em cannot be the same.
     *
     * Upon return from this function, the y vector will contain the result of
     * the application of A to x.
     *
     * initializeOperatorState must be called prior to any calls to
     * applyOperator.
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y output: y=Ax
     */
    virtual void
    apply(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y)
        {
            // Get the vector components.
            const int U_in_idx = x.getComponentDescriptorIndex(0);
            const int P_in_idx = x.getComponentDescriptorIndex(1);

            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_in_var = x.getComponentVariable(0);
            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_in_var = x.getComponentVariable(1);

            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_in_sc_var = U_in_var;
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_in_cc_var = P_in_var;

            const int U_out_idx = y.getComponentDescriptorIndex(0);
            const int P_out_idx = y.getComponentDescriptorIndex(1);

            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_out_var = y.getComponentVariable(0);
            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_out_var = y.getComponentVariable(1);

            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_out_sc_var = U_out_var;
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_out_cc_var = P_out_var;

            // Reset the interpolation operators and fill the data.
            typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
            InterpolationTransactionComponent U_component(U_in_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);  // XXXX
            InterpolationTransactionComponent P_component(P_in_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);

            std::vector<InterpolationTransactionComponent> transaction_comps(2);
            transaction_comps[0] = U_component;
            transaction_comps[1] = P_component;

            d_U_P_bdry_fill_op->resetTransactionComponents(transaction_comps);
            d_U_P_bdry_fill_op->fillData(d_apply_time);

            // Compute the action of the operator:
            //      A*[u;p] = [((rho/dt)*I-0.5*mu*L)*u + grad p; -div u].
            bool cf_bdry_synch;
            cf_bdry_synch = true;
            d_hier_math_ops->grad(
                U_out_idx, U_out_sc_var,
                cf_bdry_synch,
                1.0, P_in_idx, P_in_cc_var, d_no_fill_op, d_apply_time);
            cf_bdry_synch = false;
            d_hier_math_ops->div(
                P_out_idx, P_out_cc_var,
                -1.0, U_in_idx, U_in_sc_var, d_no_fill_op, d_apply_time,
                cf_bdry_synch);

            SAMRAI::solv::PoissonSpecifications helmholtz_spec("INSOperator::helmholtz_spec");
            helmholtz_spec.setCConstant((d_rho/d_dt)+0.5*d_lambda);
            helmholtz_spec.setDConstant(            -0.5*d_mu    );
            d_hier_math_ops->laplace(
                U_out_idx, U_out_sc_var,
                helmholtz_spec,
                U_in_idx, U_in_sc_var,
                d_no_fill_op, d_apply_time,
                1.0,
                U_out_idx, U_out_sc_var);
            return;
        }// apply

    //\}

    /*!
     * \name Logging functions.
     */
    //\{

    /*!
     * \brief Enable or disable logging.
     *
     * \param enabled logging state: true=on, false=off
     */
    virtual void
    enableLogging(
        bool enabled=true)
        {
            // intentionally blank
            return;
        }// enableLogging

    /*!
     * \brief Print out internal class data for debugging.
     */
    virtual void
    printClassData(
        std::ostream& os) const
        {
            // intentionally blank
            return;
        }// printClassData

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSOperator(
        const INSOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSOperator&
    operator=(
        const INSOperator& that);

    // Problem coefficients.
    const double d_rho;
    const double d_mu;
    const double d_lambda;

    // The timestep size.
    const double d_dt;

    // The simulation time.
    const double d_apply_time;

    // Math objects.
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;

    // Boundary condition objects.
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_U_P_bdry_fill_op, d_no_fill_op;
};

// Preconditioner used with the Krylov solver.
class INSProjectionPreconditioner
    : public IBTK::LinearSolver
{
public:
    /*!
     * XXXX
     */
    INSProjectionPreconditioner(
        const double rho,
        const double mu,
        const double lambda,
        const double current_time,
        const double new_time,
        SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver> helmholtz_solver,
        SAMRAI::tbox::Pointer<HierarchyProjector> hier_projector,
        SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > hier_cc_data_ops,
        SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> > hier_sc_data_ops,
        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
        SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> U_bdry_fill_op,
        SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> P_bdry_fill_op)
        : d_rho(rho),
          d_mu(mu),
          d_lambda(lambda),
          d_current_time(current_time),
          d_new_time(new_time),
          d_helmholtz_solver(helmholtz_solver),
          d_hier_projector(hier_projector),
          d_hier_cc_data_ops(hier_cc_data_ops),
          d_hier_sc_data_ops(hier_sc_data_ops),
          d_hier_math_ops(hier_math_ops),
          d_U_bdry_fill_op(U_bdry_fill_op),
          d_P_bdry_fill_op(P_bdry_fill_op),
          d_no_fill_op(SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation>(NULL))
        {
            // Get the control volume weight variables and patch data descriptor
            // indices.
            d_wgt_cc_var = d_hier_math_ops->getCellWeightVariable();
            d_wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

            d_wgt_sc_var = d_hier_math_ops->getSideWeightVariable();
            d_wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();

            // Get the volume of the physical domain.
            d_volume = d_hier_math_ops->getVolumeOfPhysicalDomain();
            return;
        }// INSProjectionPreconditioner

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~INSProjectionPreconditioner()
        {
            // intentionally blank
            return;
        }// ~INSProjectionPreconditioner

    /*!
     * \name Linear solver functionality.
     */
    //\{

    /*!
     * \brief Set y = P[x].
     */
    virtual bool
    solveSystem(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x)
        {
            const double dt = d_new_time-d_current_time;

            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
            const int coarsest_ln = 0;
            const int finest_ln = hierarchy->getFinestLevelNumber();

            // Get the vector components.
            const int U_in_idx = x.getComponentDescriptorIndex(0);
            const int P_in_idx = x.getComponentDescriptorIndex(1);

            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_in_var = x.getComponentVariable(0);
            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_in_var = x.getComponentVariable(1);

            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_in_sc_var = U_in_var;
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_in_cc_var = P_in_var;

            const int U_out_idx = y.getComponentDescriptorIndex(0);
            const int P_out_idx = y.getComponentDescriptorIndex(1);

            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_out_var = y.getComponentVariable(0);
            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_out_var = y.getComponentVariable(1);

            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_out_sc_var = U_out_var;
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_out_cc_var = P_out_var;

            // Setup the solver vectors.
            SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > U_out_vec, U_rhs_vec;
            U_out_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
                "INSProjectionPreconditioner::U_out_vec", hierarchy, coarsest_ln, finest_ln);
            U_out_vec->addComponent(U_out_sc_var, U_out_idx, d_wgt_sc_idx, d_hier_sc_data_ops);
            U_rhs_vec = U_out_vec->cloneVector("INSProjectionPreconditioner::U_rhs_vec");
            U_rhs_vec->allocateVectorData(d_new_time);

            const int U_rhs_idx = U_rhs_vec->getComponentDescriptorIndex(0);
            const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_rhs_sc_var = U_rhs_vec->getComponentVariable(0);

            SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > Grad_P_vec;
            Grad_P_vec = U_out_vec->cloneVector("INSProjectionPreconditioner::Grad_P_vec");
            Grad_P_vec->allocateVectorData(d_new_time);
            const int Grad_P_idx = Grad_P_vec->getComponentDescriptorIndex(0);
            const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > Grad_P_sc_var = Grad_P_vec->getComponentVariable(0);

            SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > P_out_vec, P_scratch_vec;
            P_out_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
                "INSProjectionPreconditioner::P_out_vec", hierarchy, coarsest_ln, finest_ln);
            P_out_vec->addComponent(P_out_cc_var, P_out_idx, d_wgt_cc_idx, d_hier_cc_data_ops);
            P_scratch_vec = P_out_vec->cloneVector("INSProjectionPreconditioner::P_scratch_vec");
            P_scratch_vec->allocateVectorData(d_current_time);

            const int P_scratch_idx = P_scratch_vec->getComponentDescriptorIndex(0);
            const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_scratch_cc_var = P_scratch_vec->getComponentVariable(0);

            // Reset the interpolation operators.
            typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
            InterpolationTransactionComponent U_out_component(U_out_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);  // XXXX
            InterpolationTransactionComponent P_scratch_component(P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);

            d_U_bdry_fill_op->resetTransactionComponent(U_out_component);
            d_P_bdry_fill_op->resetTransactionComponent(P_scratch_component);

            // Initialize the right-hand side terms.
            static const bool Grad_P_cf_bdry_synch = true;
            d_hier_cc_data_ops->copyData(P_scratch_idx, P_in_idx);
            d_hier_math_ops->grad(
                Grad_P_idx, Grad_P_sc_var,
                Grad_P_cf_bdry_synch,
                1.0,
                P_scratch_idx, P_scratch_cc_var,
                d_P_bdry_fill_op, d_current_time);

            SAMRAI::solv::PoissonSpecifications rhs_spec("INSProjectionPreconditioner::rhs_spec");
            rhs_spec.setCConstant((d_rho/dt)-0.5*d_lambda);
            rhs_spec.setDConstant(          +0.5*d_mu    );
            d_hier_sc_data_ops->copyData(U_out_idx, U_in_idx);
            d_hier_math_ops->laplace(
                U_rhs_idx, U_rhs_sc_var,
                rhs_spec,
                U_out_idx, U_out_sc_var,
                d_U_bdry_fill_op, d_current_time,
                -1.0,
                Grad_P_idx, Grad_P_sc_var);

            // Setup the linear solver.
            SAMRAI::solv::PoissonSpecifications helmholtz_spec("INSProjectionPreconditioner::helmholtz_spec");
            helmholtz_spec.setCConstant((d_rho/dt)+0.5*d_lambda);
            helmholtz_spec.setDConstant(          -0.5*d_mu    );
            SAMRAI::tbox::Pointer<IBTK::SCLaplaceOperator> helmholtz_operator = new IBTK::SCLaplaceOperator(
                "INSProjectionPreconditioner::Helmholtz Operator", helmholtz_spec, NULL);  // XXXX

            helmholtz_operator->setPoissonSpecifications(helmholtz_spec);
            //helmholtz_operator->setPhysicalBcCoefs(U_bc_coefs);  // XXXX
            helmholtz_operator->setHomogeneousBc(true);
            helmholtz_operator->setTime(d_new_time);
            helmholtz_operator->setHierarchyMathOps(d_hier_math_ops);

            SAMRAI::tbox::Pointer<IBTK::KrylovLinearSolver> linear_solver = new IBTK::PETScKrylovLinearSolver(
                "INSProjectionPreconditioner::PETSc Krylov solver", "adv_diff_");
            linear_solver->setInitialGuessNonzero(true);
            linear_solver->setOperator(helmholtz_operator);

            // Solve for u^{*}.
            linear_solver->solveSystem(*U_out_vec,*U_rhs_vec);
            SAMRAI::tbox::plog << "INSProjectionPreconditioner::integrateHierarchy(): linear solve number of iterations = " << linear_solver->getNumIterations() << "\n";
            SAMRAI::tbox::plog << "INSProjectionPreconditioner::integrateHierarchy(): linear solve residual norm        = " << linear_solver->getResidualNorm()  << "\n";
            if (linear_solver->getNumIterations() == linear_solver->getMaxIterations())
            {
                SAMRAI::tbox::pout << "INSProjectionPreconditioner::integrateHierarchy():"
                                   <<"  WARNING: linear solver iterations == max iterations\n";
            }

            // Reset the intermediate velocity u^{*} := u^{*} + dt/rho grad p^{n-1/2}.
            //
            // XXXX: Switch back to standard formulation to allow for correct
            // BCs for phi...?
            //
            // NOTE: If this switch is made, be sure to do a full copy of P_scratch
            // into P_scratch (interior + ghost cell region) to get BCs correct!
            d_hier_sc_data_ops->axpy(U_out_idx, dt/d_rho, Grad_P_idx, U_out_idx);

            // Project the intermediate velocity u^{*}.
            d_hier_projector->projectHierarchy(
                d_rho, dt, d_current_time+0.5*dt,
                "pressure_update",
                U_out_idx, U_out_sc_var,
                P_in_idx, P_in_cc_var,
                P_scratch_idx, P_scratch_cc_var,
                Grad_P_idx, Grad_P_sc_var,
                U_out_idx, U_out_sc_var);

            // Update the pressure.
            d_hier_cc_data_ops->copyData(P_out_idx, P_scratch_idx);
#if 0 // XXXX
            SAMRAI::solv::PoissonSpecifications pressure_helmholtz_spec("INSProjectionPreconditioner::helmholtz_spec");
            pressure_helmholtz_spec.setCConstant(+0.5*dt*d_lambda/d_rho);
            pressure_helmholtz_spec.setDConstant(-0.5*dt*d_mu    /d_rho);
            d_hier_cc_data_ops->subtract(P_scratch_idx, P_scratch_idx, P_in_idx);
            d_hier_math_ops->laplace(
                P_out_idx, P_out_cc_var,
                pressure_helmholtz_spec,
                P_scratch_idx, P_scratch_cc_var,
                d_P_bdry_fill_op, d_current_time,
                1.0,
                P_out_idx, P_out_cc_var);
#endif
            // Normalize p^{n+1/2} to have mean (discrete integral) zero.
            if (true) // XXXX (d_normalize_pressure)
            {
                const double P_mean = (1.0/d_volume)*d_hier_cc_data_ops->integral(P_out_idx, d_wgt_cc_idx);
                d_hier_cc_data_ops->addScalar(P_out_idx, P_out_idx, -P_mean);
            }

            // Deallocate scratch data.
            U_rhs_vec->deallocateVectorData();
            Grad_P_vec->deallocateVectorData();
            P_scratch_vec->deallocateVectorData();
            return true;
        }// solveSystem

    //\}

    /*!
     * \name Functions to access solver parameters.
     */
    //\{

    /*!
     * \brief Set whether the initial guess is non-zero.
     */
    virtual void
    setInitialGuessNonzero(
        bool initial_guess_nonzero=true)
        {
            // intentionally blank
            return;
        }// setInitialGuessNonzero


    /*!
     * \brief Get whether the initial guess is non-zero.
     */
    virtual bool
    getInitialGuessNonzero() const
        {
            // intentionally blank
            return true;
        }// getInitialGuessNonzero

    /*!
     * \brief Set the maximum number of iterations to use per solve.
     */
    virtual void
    setMaxIterations(
        int max_iterations)
        {
            // intentionally blank
            return;
        }// setMaxIterations

    /*!
     * \brief Get the maximum number of iterations to use per solve.
     */
    virtual int
    getMaxIterations() const
        {
            // intentionally blank
            return 1;
        }// getMaxIterations

    /*!
     * \brief Set the absolute residual tolerance for convergence.
     */
    virtual void
    setAbsoluteTolerance(
        double abs_residual_tol)
        {
            // intentionally blank
            return;
        }// setAbsoluteTolerance

    /*!
     * \brief Get the absolute residual tolerance for convergence.
     */
    virtual double
    getAbsoluteTolerance() const
        {
            // intentionally blank
            return 0.0;
        }// getAbsoluteTolerance

    /*!
     * \brief Set the relative residual tolerance for convergence.
     */
    virtual void
    setRelativeTolerance(
        double rel_residual_tol)
        {
            // intentionally blank
            return;
        }// setRelativeTolerance

    /*!
     * \brief Get the relative residual tolerance for convergence.
     */
    virtual double
    getRelativeTolerance() const
        {
            // intentionally blank
            return 0.0;
        }// getRelativeTolerance

    //\}

    /*!
     * \name Functions to access data on the most recent solve.
     */
    //\{

    /*!
     * \brief Return the iteration count from the most recent linear solve.
     */
    virtual int
    getNumIterations() const
        {
            // intentionally blank
            return 0;
        }// getNumIterations

    /*!
     * \brief Return the residual norm from the most recent iteration.
     */
    virtual double
    getResidualNorm() const
        {
            return 0.0;
        }// getResidualNorm

    //\}

    /*!
     * \name Logging functions.
     */
    //\{

    /*!
     * \brief Enable or disable logging.
     *
     * \param enabled logging state: true=on, false=off
     */
    virtual void
    enableLogging(
        bool enabled=true)
        {
            // intentionally blank
            return;
        }// enableLogging

    /*!
     * \brief Print out internal class data for debugging.
     */
    virtual void
    printClassData(
        std::ostream& os) const
        {
            // intentionally blank
            return;
        }// printClassData

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSProjectionPreconditioner();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSProjectionPreconditioner(
        const INSProjectionPreconditioner& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSProjectionPreconditioner&
    operator=(
        const INSProjectionPreconditioner& that);

    // Problem coefficients.
    const double d_rho;
    const double d_mu;
    const double d_lambda;

    // The simulation time.
    const double d_current_time, d_new_time;

    // Linear solvers.
    SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver> d_helmholtz_solver;

    // Projection functionality.
    SAMRAI::tbox::Pointer<HierarchyProjector> d_hier_projector;

    // Math objects.
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> > d_hier_sc_data_ops;
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_wgt_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_wgt_sc_var;
    int d_wgt_cc_idx, d_wgt_sc_idx;
    double d_volume;

    // Boundary condition objects.
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_U_bdry_fill_op, d_P_bdry_fill_op, d_no_fill_op;
};

}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredHierarchyIntegrator::INSStaggeredHierarchyIntegrator(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    SAMRAI::tbox::Pointer<HierarchyProjector> hier_projector,
    bool register_for_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT(!hier_projector.isNull());
#endif
    d_object_name = object_name;
    d_registered_for_restart = register_for_restart;

    d_hierarchy = hierarchy;

    d_hier_projector = hier_projector;

    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Set some default values.
    d_U_scale = 1.0;
    d_P_scale = 1.0;
    d_F_scale = 1.0;

    d_start_time = 0.0;
    d_end_time = std::numeric_limits<double>::max();
    d_grow_dt = 2.0;
    d_max_integrator_steps = std::numeric_limits<int>::max();

    d_using_vorticity_tagging = false;
    d_Omega_max = 0.0;

    d_second_order_pressure_update = false;
    d_normalize_pressure = false;

    d_num_init_cycles = 5;

    d_regrid_interval = 1;
    d_old_dt = -1.0;
    d_integrator_time = std::numeric_limits<double>::quiet_NaN();
    d_integrator_step = std::numeric_limits<int>::max();

    d_conservation_form = false;

    d_output_U = false;
    d_output_P = false;
    d_output_F = false;
    d_output_Omega = false;
    d_output_Div_U = false;

    d_rho    = std::numeric_limits<double>::quiet_NaN();
    d_mu     = std::numeric_limits<double>::quiet_NaN();
    d_nu     = std::numeric_limits<double>::quiet_NaN();
    d_lambda = std::numeric_limits<double>::quiet_NaN();

    d_cfl = 0.9;

    d_dt_max = std::numeric_limits<double>::max();
    d_dt_max_time_max = std::numeric_limits<double>::max();
    d_dt_max_time_min = -(d_dt_max_time_max-std::numeric_limits<double>::epsilon());

    d_is_initialized = false;

    d_do_log = false;

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, from_restart);

    // Determine the properties of the projection scheme.
    SAMRAI::tbox::pout << d_object_name << ": using a STAGGERED GRID exact projection method\n"
                       << d_object_name << ": pressure update is " << (d_second_order_pressure_update ? "SECOND-ORDER" : "FIRST-ORDER") << " accurate\n"
                       << d_object_name << ": advection discretization employs " << (d_conservation_form ? "CONSERVATIVE" : "NON-CONSERVATIVE") << " differencing\n";

    // Obtain the Hierarchy data operations objects.
    SAMRAI::math::HierarchyDataOpsManager<NDIM>* hier_ops_manager =
        SAMRAI::math::HierarchyDataOpsManager<NDIM>::getManager();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cc_var =
        new SAMRAI::pdat::CellVariable<NDIM,double>("cc_var");
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > sc_var =
        new SAMRAI::pdat::SideVariable<NDIM,double>("sc_var");

    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, hierarchy);
    d_hier_sc_data_ops = hier_ops_manager->getOperationsDouble(sc_var, hierarchy);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_initialize_hierarchy_integrator = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::advanceHierarchy()");
        t_get_stable_timestep = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::getStableTimestep()");
        t_regrid_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::regridHierarchy()");
        t_integrate_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::integrateHierarchy()");
        t_synchronize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_data_to_preadvance_state = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSStaggeredHierarchyIntegrator::putToDatabase()");
        timers_need_init = false;
    }
    return;
}// INSStaggeredHierarchyIntegrator

INSStaggeredHierarchyIntegrator::~INSStaggeredHierarchyIntegrator()
{
    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }

    for (RefinePatchStrategyMap::iterator it = d_rstrategies.begin();
         it != d_rstrategies.end(); ++it)
    {
        delete (*it).second;
    }

    for (CoarsenPatchStrategyMap::iterator it = d_cstrategies.begin();
         it != d_cstrategies.end(); ++it)
    {
        delete (*it).second;
    }
    return;
}// ~INSStaggeredHierarchyIntegrator

const std::string&
INSStaggeredHierarchyIntegrator::getName() const
{
    return d_object_name;
}// getName

void
INSStaggeredHierarchyIntegrator::registerVelocityInitialConditions(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> U_init)
{
    d_U_init = U_init;
    return;
}// registerVelocityInitialConditions

void
INSStaggeredHierarchyIntegrator::registerPressureInitialConditions(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> P_init)
{
    d_P_init = P_init;
    return;
}// registerPressureInitialConditions

void
INSStaggeredHierarchyIntegrator::registerBodyForceSpecification(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> F_set)
{
    d_F_set = F_set;
    return;
}// registerBodyForceSpecification

void
INSStaggeredHierarchyIntegrator::registerVisItDataWriter(
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!visit_writer.isNull());
#endif
    d_visit_writer = visit_writer;
    return;
}// registerVisItDataWriter

///
///  The following routines:
///
///      getHierarchyMathOps(),
///      setHierarchyMathOps(),
///      isManagingHierarchyMathOps()
///
///  allow for the sharing of a single HierarchyMathOps object between mutiple
///  HierarchyIntegrator objects.
///

SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps>
INSStaggeredHierarchyIntegrator::getHierarchyMathOps() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hier_math_ops.isNull());
#endif
    return d_hier_math_ops;
}// getHierarchyMathOps

void
INSStaggeredHierarchyIntegrator::setHierarchyMathOps(
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
    const bool manage_ops)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hier_math_ops.isNull());
#endif
    d_hier_math_ops = hier_math_ops;
    d_is_managing_hier_math_ops = manage_ops;
    return;
}// setHierarchyMathOps

bool
INSStaggeredHierarchyIntegrator::isManagingHierarchyMathOps() const
{
    return d_is_managing_hier_math_ops;
}// isManagingHierarchyMathOps

///
///  The following routines:
///
///      initializeHierarchyIntegrator(),
///      initializeHierarchy(),
///      advanceHierarchy(),
///      atRegridPoint(),
///      getIntegratorTime(),
///      getStartTime(),
///      getEndTime(),
///      getIntegratorStep(),
///      getMaxIntegratorSteps(),
///      stepsRemaining(),
///      getPatchHierarchy(),
///      getGriddingAlgorithm(),
///      getHierarchyProjector()
///
///  allow the INSStaggeredHierarchyIntegrator to be used as a hierarchy integrator.
///

void
INSStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg)
{
    t_initialize_hierarchy_integrator->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!gridding_alg.isNull());
#endif
    d_gridding_alg = gridding_alg;

    // Setup the tag buffer.
    if (d_using_default_tag_buffer)
    {
        d_tag_buffer.resizeArray(d_gridding_alg->getMaxLevels());
        for (int i = 0; i < d_gridding_alg->getMaxLevels(); ++i)
        {
            d_tag_buffer[i] = d_regrid_interval;
        }
    }
    else
    {
        if (d_tag_buffer.getSize() < d_gridding_alg->getMaxLevels())
        {
            int tsize = d_tag_buffer.getSize();
            d_tag_buffer.resizeArray(d_gridding_alg->getMaxLevels());
            for (int i = tsize; i < d_gridding_alg->getMaxLevels(); ++i)
            {
                d_tag_buffer[i] = d_tag_buffer[tsize-1];
            }
        }
    }

    // Initialize all variable contexts.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    d_current_context = var_db->getContext(d_object_name+"::CURRENT");
    d_new_context     = var_db->getContext(d_object_name+"::NEW"    );
    d_scratch_context = var_db->getContext(d_object_name+"::SCRATCH");

    // Initialize all variables.
    d_U_var          = new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::U"          );
    d_U_cc_var       = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::U_cc",  NDIM);
    d_P_var          = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::P"          );
    d_F_var          = new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::F"          );
    d_F_cc_var       = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::F_cc",  NDIM);
#if ( NDIM == 2)
    d_Omega_var      = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Omega"      );
#endif
#if ( NDIM == 3)
    d_Omega_var      = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Omega", NDIM);
    d_Omega_Norm_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::||Omega||_2");
#endif
    d_Div_U_var      = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Div_U"      );
    d_Grad_P_var     = new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::Grad_P"     );
    d_gadvect_U_var  = new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::gadvect_U"  );
    d_gadvect_F_var  = new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::gadvect_F"  );

    // Create the default communication algorithms.
    d_fill_after_regrid = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    d_calgs["SYNCH_CURRENT_STATE_DATA"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    d_calgs["SYNCH_NEW_STATE_DATA"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();

    // Register state variables that are maintained by the
    // INSStaggeredHierarchyIntegrator.

    const SAMRAI::hier::IntVector<NDIM> cell_ghosts = CELLG;
    const SAMRAI::hier::IntVector<NDIM> side_ghosts = SIDEG;
    const SAMRAI::hier::IntVector<NDIM> gadvect_velocity_ghosts = GADVECTG;
    const SAMRAI::hier::IntVector<NDIM> gadvect_force_ghosts = NUM_GODUNOV_CYCLES;
    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    registerVariable(d_U_current_idx, d_U_new_idx, d_U_scratch_idx,
                     d_U_var, side_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_U_cc_current_idx, d_U_cc_new_idx, d_U_cc_scratch_idx,
                     d_U_cc_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_P_current_idx, d_P_new_idx, d_P_scratch_idx,
                     d_P_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE");

    registerVariable(d_F_current_idx, d_F_new_idx, d_F_scratch_idx,
                     d_F_var, side_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_F_cc_current_idx, d_F_cc_new_idx, d_F_cc_scratch_idx,
                     d_F_cc_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_Omega_current_idx, d_Omega_new_idx, d_Omega_scratch_idx,
                     d_Omega_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE");
#if (NDIM == 3)
    registerVariable(d_Omega_Norm_current_idx, d_Omega_Norm_new_idx, d_Omega_Norm_scratch_idx,
                     d_Omega_Norm_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE");
#endif
    registerVariable(d_Div_U_current_idx, d_Div_U_new_idx, d_Div_U_scratch_idx,
                     d_Div_U_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    // Register scratch variables that are maintained by the
    // INSStaggeredHierarchyIntegrator.
    registerVariable(d_Grad_P_scratch_idx, d_Grad_P_var, side_ghosts);

    // Register other scratch variables.
    d_gadvect_U_scratch_idx = var_db->registerVariableAndContext(
        d_gadvect_U_var, getScratchContext(), gadvect_velocity_ghosts);
    d_gadvect_F_scratch_idx = var_db->registerVariableAndContext(
        d_gadvect_F_var, getScratchContext(), gadvect_force_ghosts);

    // Register variables for plotting.
    if (!d_visit_writer.isNull())
    {
        if (d_output_U)
        {
            d_visit_writer->registerPlotQuantity(
                d_U_var->getName(), "VECTOR", d_U_cc_current_idx, 0, d_U_scale);

            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(
                    d_U_var->getName()+stream.str(), "SCALAR", d_U_cc_current_idx, d, d_U_scale);
            }
        }

        if (d_output_P)
        {
            d_visit_writer->registerPlotQuantity(
                d_P_var->getName(), "SCALAR", d_P_current_idx, 0, d_P_scale);
        }

        if (d_output_F)
        {
            d_visit_writer->registerPlotQuantity(
                d_F_var->getName(), "VECTOR", d_F_cc_current_idx, 0, d_F_scale);

            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(
                    d_F_var->getName()+stream.str(), "SCALAR", d_F_cc_current_idx, d, d_F_scale);
            }
        }

        if (d_output_Omega)
        {
            d_visit_writer->registerPlotQuantity(
                d_Omega_var->getName(), (NDIM == 2) ? "SCALAR" : "VECTOR", d_Omega_current_idx);
#if (NDIM == 3)
            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(
                    d_Omega_var->getName()+stream.str(), "SCALAR", d_Omega_current_idx, d);
            }
#endif
        }

        if (d_output_Div_U)
        {
            d_visit_writer->registerPlotQuantity(
                d_Div_U_var->getName(), "SCALAR", d_Div_U_current_idx);
        }
    }

    // Create several refinement communications algorithms, used in filling
    // ghost cell data.
    std::string ralg_name;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator;

    ralg_name = "gadvect_U_scratch_bdry_fill";
    d_ralgs[ralg_name] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(d_gadvect_U_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs[ralg_name]->registerRefine(d_gadvect_U_scratch_idx, // destination
                                       d_U_current_idx,         // source
                                       d_gadvect_U_scratch_idx, // temporary work space
                                       refine_operator);
    d_rstrategies[ralg_name] = new IBTK::CartExtrapPhysBdryOp(d_gadvect_U_scratch_idx, BDRY_EXTRAP_TYPE);

    ralg_name = "gadvect_F_scratch_bdry_fill";
    d_ralgs[ralg_name] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(d_gadvect_F_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs[ralg_name]->registerRefine(d_gadvect_F_scratch_idx, // destination
                                       d_gadvect_F_scratch_idx, // source
                                       d_gadvect_F_scratch_idx, // temporary work space
                                       refine_operator);
    d_rstrategies[ralg_name] = new IBTK::CartExtrapPhysBdryOp(d_gadvect_F_scratch_idx, BDRY_EXTRAP_TYPE);

    // Setup the Hierarchy math operations object.
    d_hier_math_ops = new IBTK::HierarchyMathOps(d_object_name+"::HierarchyMathOps", d_hierarchy);
    d_is_managing_hier_math_ops = true;
    d_hier_projector->setHierarchyMathOps(d_hier_math_ops);

    // Set the current integration time.
    if (!SAMRAI::tbox::RestartManager::getManager()->isFromRestart())
    {
        d_integrator_time = d_start_time;
        d_integrator_step = 0;
    }

    // Indicate that the integrator has been initialized.
    d_is_initialized = true;

    t_initialize_hierarchy_integrator->stop();
    return;
}// initializeHierarchyIntegrator

double
INSStaggeredHierarchyIntegrator::initializeHierarchy()
{
    t_initialize_hierarchy->start();

    if (!d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::initializeHierarchy():\n" <<
                   "  must call initializeHierarchyIntegrator() prior to calling initializeHierarchy()." << std::endl);
    }

    // Initialize the patch hierarchy.
    const bool initial_time = !SAMRAI::tbox::RestartManager::getManager()->isFromRestart();

    if (!initial_time)
    {
        d_hierarchy->getFromRestart(d_gridding_alg->getMaxLevels());
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        d_gridding_alg->getTagAndInitializeStrategy()->
            resetHierarchyConfiguration(d_hierarchy, coarsest_ln, finest_ln);
    }
    else
    {
        d_gridding_alg->makeCoarsestLevel(d_hierarchy,d_start_time);

        int level_number = 0;
        bool done = false;
        while (!done && (d_gridding_alg->levelCanBeRefined(level_number)))
        {
            d_gridding_alg->
                makeFinerLevel(d_hierarchy,
                               d_integrator_time, initial_time,
                               d_tag_buffer[level_number]);

            done = !d_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        // After data on each level is initialized at simulation start time,
        // coarser levels are synchronized with finer levels that didn't exist
        // when the coarser level initial data was set.
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();

        if (finest_ln > 0)
        {
            synchronizeNewLevels(d_hierarchy, coarsest_ln, finest_ln,
                                 d_start_time, initial_time);
        }
    }

    // The next timestep is given by the minimum allowable timestep over all
    // levels in the patch hierarchy.
    const double dt_next = getStableTimestep(getCurrentContext());

    t_initialize_hierarchy->stop();
    return dt_next;
}// initializeHierarchy

double
INSStaggeredHierarchyIntegrator::advanceHierarchy(
    const double dt)
{
    t_advance_hierarchy->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_end_time >= d_integrator_time+dt);
#endif

    const double current_time = d_integrator_time;
    const double new_time = d_integrator_time+dt;
    const bool initial_time = SAMRAI::tbox::MathUtilities<double>::equalEps(d_integrator_time,d_start_time);

    // Set the guess for the initial pressure to zero.
    if (initial_time)
    {
        d_hier_cc_data_ops->setToScalar(d_P_current_idx, 0.0);
    }

    // Regrid the patch hierarchy.
    const bool do_regrid = ((d_regrid_interval == 0)
                            ? false
                            : (d_integrator_step % d_regrid_interval == 0));
    if (do_regrid)
    {
        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): regridding prior to timestep " << d_integrator_step << "\n";
        regridHierarchy();
    }

    // Integrate all time dependent data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): integrating time dependent data\n";
    integrateHierarchy(current_time, new_time);

    // Synchronize all data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): synchronizing data\n";
    synchronizeHierarchy();

    // Reset all time dependent data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): resetting time dependent data\n";
    resetTimeDependentHierData(new_time);

    // Determine the next stable timestep from u(n+1).
    const double dt_next = getStableTimestep(getCurrentContext());

    t_advance_hierarchy->stop();
    return dt_next;
}// advanceHierarchy

double
INSStaggeredHierarchyIntegrator::getStableTimestep(
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx)
{
    t_get_stable_timestep->start();

    const bool initial_time = SAMRAI::tbox::MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    double dt_next = std::numeric_limits<double>::max();

    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        dt_next = std::min(d_cfl*getLevelDt(level,ctx),dt_next);
    }

    if (d_integrator_time+dt_next >= d_end_time)
    {
        dt_next = d_end_time - d_integrator_time;
    }

    if (d_integrator_time >= d_dt_max_time_min &&
        d_integrator_time <= d_dt_max_time_max)
    {
        dt_next = std::min(d_dt_max,dt_next);
    }

    if (!initial_time)
    {
        dt_next = std::min(dt_next,d_grow_dt*d_old_dt);
    }

    t_get_stable_timestep->stop();
    return dt_next;
}// getStableTimestep()

bool
INSStaggeredHierarchyIntegrator::atRegridPoint() const
{
    const int level_number = 0;
    return ((d_integrator_step > 0)
            && d_gridding_alg->levelCanBeRefined(level_number)
            && (d_regrid_interval == 0
                ? false
                : (d_integrator_step % d_regrid_interval == 0)));
}// atRegridPoint

double
INSStaggeredHierarchyIntegrator::getIntegratorTime() const
{
    return d_integrator_time;
}// getIntegratorTime

double
INSStaggeredHierarchyIntegrator::getStartTime() const
{
    return d_start_time;
}// getStartTime

double
INSStaggeredHierarchyIntegrator::getEndTime() const
{
    return d_end_time;
}// getEndTime

int
INSStaggeredHierarchyIntegrator::getIntegratorStep() const
{
    return d_integrator_step;
}// getIntegratorStep

int
INSStaggeredHierarchyIntegrator::getMaxIntegratorSteps() const
{
    return d_max_integrator_steps;
}// getMaxIntegratorSteps

bool
INSStaggeredHierarchyIntegrator::stepsRemaining() const
{
    return (d_integrator_step < d_max_integrator_steps);
}// stepsRemaining

const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >
INSStaggeredHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> >
INSStaggeredHierarchyIntegrator::getGriddingAlgorithm() const
{
    return d_gridding_alg;
}// getGriddingAlgorithm

SAMRAI::tbox::Pointer<HierarchyProjector>
INSStaggeredHierarchyIntegrator::getHierarchyProjector() const
{
    return d_hier_projector;
}// getHierarchyProjector

///
///  The following routines:
///
///      regridHierarchy(),
///      integrateHierarchy(),
///      synchronizeHierarchy(),
///      synchronizeNewLevels(),
///      resetTimeDependentHierData(),
///      resetHierDataToPreadvanceState()
///
///  allow the INSStaggeredHierarchyIntegrator to provide data management for a time
///  integrator which making use of this class.
///

void
INSStaggeredHierarchyIntegrator::regridHierarchy()
{
    t_regrid_hierarchy->start();

    const int coarsest_ln = 0;
    d_gridding_alg->regridAllFinerLevels(d_hierarchy, coarsest_ln, d_integrator_time, d_tag_buffer);
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Synchronize the state data on the patch hierarchy.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["SYNCH_CURRENT_STATE_DATA"][ln]->coarsenData();
    }

    t_regrid_hierarchy->stop();
    return;
}// regridHierarchy

double
INSStaggeredHierarchyIntegrator::integrateHierarchy(
    const double current_time,
    const double new_time)
{
    t_integrate_hierarchy->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    const bool initial_time = SAMRAI::tbox::MathUtilities<double>::equalEps(d_integrator_time,d_start_time);

    // Synchronize current state data.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["SYNCH_CURRENT_STATE_DATA"][ln]->coarsenData();
    }

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data    ,     new_time);
    }

    const int num_cycles = initial_time ? d_num_init_cycles : 1;
    for (int cycle = 0; cycle < num_cycles; ++cycle)
    {
        if (d_do_log && initial_time)
        {
            SAMRAI::tbox::plog << "\n\n++++++++++++++++++++++++++++++++++++++++++++++++\n";
            SAMRAI::tbox::plog << "+\n";
            SAMRAI::tbox::plog << "+ Performing cycle " << cycle+1 << " of " << d_num_init_cycles << " to initialize the pressure at time t=dt/2\n";
            SAMRAI::tbox::plog << "+\n";
            SAMRAI::tbox::plog << "++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
        }

        // Setup the solver vectors.
        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > U_scratch_vec, U_rhs_vec;
        U_scratch_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
            d_object_name+"::U_scratch_vec", d_hierarchy, 0, finest_ln);
        U_scratch_vec->addComponent(d_U_var, d_U_scratch_idx, d_wgt_sc_idx, d_hier_sc_data_ops);
        U_rhs_vec = U_scratch_vec->cloneVector(d_object_name+"::U_rhs_vec");
        U_rhs_vec->allocateVectorData(current_time);

        const int U_rhs_idx = U_rhs_vec->getComponentDescriptorIndex(0);
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_rhs_var = U_rhs_vec->getComponentVariable(0);
#if COUPLED_SOLVER
        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > P_scratch_vec, P_rhs_vec;
        P_scratch_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
            d_object_name+"::P_scratch_vec", d_hierarchy, 0, finest_ln);
        P_scratch_vec->addComponent(d_P_var, d_P_scratch_idx, d_wgt_cc_idx, d_hier_cc_data_ops);
        P_rhs_vec = P_scratch_vec->cloneVector(d_object_name+"::P_rhs_vec");
        P_rhs_vec->allocateVectorData(current_time);

        const int P_rhs_idx = P_rhs_vec->getComponentDescriptorIndex(0);
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_rhs_var = P_rhs_vec->getComponentVariable(0);
        d_hier_cc_data_ops->setToScalar(P_rhs_idx,0.0);
#endif
        // Initialize the right-hand side terms.
        static const bool Grad_P_scratch_cf_bdry_synch = true;
        d_hier_cc_data_ops->copyData(d_P_scratch_idx, d_P_current_idx);
        d_hier_math_ops->grad(
            d_Grad_P_scratch_idx, d_Grad_P_var,
            Grad_P_scratch_cf_bdry_synch,
            1.0,
            d_P_scratch_idx, d_P_var,
            d_P_scratch_bdry_fill_op, current_time);

        computeConvectiveDerivative(
            U_rhs_idx,
            d_U_current_idx,
            d_F_current_idx,
            d_Grad_P_scratch_idx,
            current_time, new_time, d_conservation_form);

        SAMRAI::solv::PoissonSpecifications rhs_spec(d_object_name+"::rhs_spec");
        rhs_spec.setCConstant((d_rho/dt)-0.5*d_lambda);
        rhs_spec.setDConstant(          +0.5*d_mu    );
        d_hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
        d_hier_math_ops->laplace(
            U_rhs_idx, U_rhs_var,
            rhs_spec,
            d_U_scratch_idx, d_U_var,
            d_U_scratch_bdry_fill_op, current_time,
            -d_rho,
            U_rhs_idx, U_rhs_var);
#if (!COUPLED_SOLVER)
        d_hier_sc_data_ops->subtract(U_rhs_idx, U_rhs_idx, d_Grad_P_scratch_idx);
#endif

        if (!d_F_set.isNull())
        {
            if (d_F_set->isTimeDependent())
            {
                d_F_set->setDataOnPatchHierarchy(
                    d_F_new_idx, d_F_var, d_hierarchy, new_time);
            }
            d_hier_sc_data_ops->axpy(U_rhs_idx, +0.5, d_F_current_idx, U_rhs_idx);
            d_hier_sc_data_ops->axpy(U_rhs_idx, +0.5, d_F_new_idx    , U_rhs_idx);
        }
        else
        {
            d_hier_sc_data_ops->setToScalar(d_F_new_idx, 0.0);
        }

        // Setup the linear solver.
        SAMRAI::solv::PoissonSpecifications helmholtz_spec(d_object_name+"::helmholtz_spec");
        helmholtz_spec.setCConstant((d_rho/dt)+0.5*d_lambda);
        helmholtz_spec.setDConstant(          -0.5*d_mu    );
        SAMRAI::tbox::Pointer<IBTK::SCLaplaceOperator> helmholtz_operator = new IBTK::SCLaplaceOperator(
            d_object_name+"::Helmholtz Operator", helmholtz_spec, NULL);  // XXXX

        helmholtz_operator->setPoissonSpecifications(helmholtz_spec);
        //helmholtz_operator->setPhysicalBcCoefs(U_bc_coefs);  // XXXX
        helmholtz_operator->setHomogeneousBc(false);
        helmholtz_operator->setTime(new_time);
        helmholtz_operator->setHierarchyMathOps(d_hier_math_ops);

        SAMRAI::tbox::Pointer<IBTK::KrylovLinearSolver> helmholtz_solver = new IBTK::PETScKrylovLinearSolver(
            d_object_name+"::PETSc Krylov solver", "adv_diff_");
        helmholtz_solver->setInitialGuessNonzero(true);
        helmholtz_solver->setOperator(helmholtz_operator);

#if (!COUPLED_SOLVER)
        // Solve for u^{*}.
        helmholtz_solver->solveSystem(*U_scratch_vec,*U_rhs_vec);
        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::integrateHierarchy(): linear solve number of iterations = " << helmholtz_solver->getNumIterations() << "\n";
        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::integrateHierarchy(): linear solve residual norm        = " << helmholtz_solver->getResidualNorm()  << "\n";
        if (helmholtz_solver->getNumIterations() == helmholtz_solver->getMaxIterations())
        {
            SAMRAI::tbox::pout << d_object_name << "::integrateHierarchy():"
                               <<"  WARNING: linear solver iterations == max iterations\n";
        }

        // Deallocate scratch data.
        U_rhs_vec->deallocateVectorData();

        // Reset the intermediate velocity u^{*} := u^{*} + dt/rho grad p^{n-1/2}.
        //
        // XXXX: Switch back to original formulation to allow for correct BCs
        // for phi...?
        d_hier_sc_data_ops->axpy(d_U_scratch_idx, dt/d_rho, d_Grad_P_scratch_idx, d_U_scratch_idx);

        // Project the intermediate velocity u^{*}.
        d_hier_projector->projectHierarchy(
            d_rho, dt, current_time+0.5*dt,
            "pressure_update",
            d_U_scratch_idx, d_U_var,
            d_P_current_idx, d_P_var,
            d_P_scratch_idx, d_P_var,
            d_Grad_P_scratch_idx, d_Grad_P_var,
            d_U_scratch_idx, d_U_var);
        d_hier_sc_data_ops->copyData(d_U_new_idx, d_U_scratch_idx); // XXXX

        // Update the pressure.
        d_hier_cc_data_ops->copyData(d_P_new_idx, d_P_scratch_idx);
        if (d_second_order_pressure_update)
        {
            SAMRAI::solv::PoissonSpecifications pressure_helmholtz_spec(d_object_name+"::helmholtz_spec");
            pressure_helmholtz_spec.setCConstant(+0.5*dt*d_lambda/d_rho);
            pressure_helmholtz_spec.setDConstant(-0.5*dt*d_mu    /d_rho);
            d_hier_cc_data_ops->subtract(d_P_scratch_idx, d_P_scratch_idx, d_P_current_idx);
            d_hier_math_ops->laplace(
                d_P_new_idx, d_P_var,
                pressure_helmholtz_spec,
                d_P_scratch_idx, d_P_var,
                d_P_scratch_bdry_fill_op, current_time,
                1.0,
                d_P_new_idx, d_P_var);
        }
#else
        // Reset the solution, rhs, and nullspace vectors.
        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > sol_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
            d_object_name+"::sol_vec", d_hierarchy, 0, finest_ln);
        sol_vec->addComponent(d_U_var,d_U_scratch_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
        sol_vec->addComponent(d_P_var,d_P_scratch_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > rhs_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
            d_object_name+"::rhs_vec", d_hierarchy, 0, finest_ln);
        rhs_vec->addComponent(d_U_var,U_rhs_idx,d_wgt_sc_idx,d_hier_sc_data_ops);
        rhs_vec->addComponent(d_P_var,P_rhs_idx,d_wgt_cc_idx,d_hier_cc_data_ops);

        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > nul_vec = sol_vec->cloneVector(d_object_name+"::nul_vec");
        nul_vec->allocateVectorData(current_time);
        d_hier_sc_data_ops->setToScalar(nul_vec->getComponentDescriptorIndex(0), 0.0);
        d_hier_cc_data_ops->setToScalar(nul_vec->getComponentDescriptorIndex(1), 1.0);

        // Setup the nullspace object.
        PetscErrorCode ierr;
        MatNullSpace petsc_nullsp;
        Vec petsc_nullsp_vec = IBTK::PETScSAMRAIVectorReal<double>::createPETScVector(nul_vec, PETSC_COMM_SELF);
        Vec vecs[] = {petsc_nullsp_vec};
        static const PetscTruth has_cnst = PETSC_FALSE;
        ierr = MatNullSpaceCreate(PETSC_COMM_SELF, has_cnst, 1, vecs, &petsc_nullsp); IBTK_CHKERRQ(ierr);

        // Setup the patch boundary filling objects.
        typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        InterpolationTransactionComponent U_scratch_component(d_U_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL); // XXXX
        InterpolationTransactionComponent P_scratch_component(d_P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);

        std::vector<InterpolationTransactionComponent> U_P_transaction_comps(2);
        U_P_transaction_comps[0] = U_scratch_component;
        U_P_transaction_comps[1] = P_scratch_component;

        SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> U_P_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
        U_P_bdry_fill_op->initializeOperatorState(U_P_transaction_comps, d_hierarchy);

        SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> U_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
        U_bdry_fill_op->initializeOperatorState(U_scratch_component, d_hierarchy);

        SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> P_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
        P_bdry_fill_op->initializeOperatorState(P_scratch_component, d_hierarchy);

        // Setup the linear operator.
        SAMRAI::tbox::Pointer<IBTK::LinearOperator> linear_op = new INSOperator(
            d_rho, d_mu, d_lambda,
            dt,
            new_time,
            d_hier_math_ops, U_P_bdry_fill_op);

        // Setup the linear solver.
        SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver> linear_solver = new IBTK::PETScKrylovLinearSolver(
            d_object_name+"::linear_solver", "ins_");
        linear_solver->setInitialGuessNonzero(true);
        linear_solver->setOperator(linear_op);
        linear_solver->initializeSolverState(*sol_vec,*rhs_vec);
        KSP petsc_ksp = linear_solver->getPETScKSP();
        ierr = KSPSetNullSpace(petsc_ksp, petsc_nullsp); IBTK_CHKERRQ(ierr);

        SAMRAI::tbox::Pointer<IBTK::LinearSolver> pc_op = new INSProjectionPreconditioner(
            d_rho, d_mu, d_lambda,
            current_time, new_time,
            helmholtz_solver, d_hier_projector,
            d_hier_cc_data_ops, d_hier_sc_data_ops, d_hier_math_ops,
            U_bdry_fill_op, P_bdry_fill_op);
        linear_solver->setPreconditioner(pc_op);

        // Solve system.
        linear_solver->solveSystem(*sol_vec,*rhs_vec);

        // Pull out solution components.
        d_hier_sc_data_ops->copyData(d_U_new_idx, sol_vec->getComponentDescriptorIndex(0));
        d_hier_cc_data_ops->copyData(d_P_new_idx, sol_vec->getComponentDescriptorIndex(1));

        // Deallocate the nullspace object.
        if (petsc_nullsp != static_cast<MatNullSpace>(NULL))
        {
            ierr = MatNullSpaceDestroy(petsc_nullsp); IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(petsc_nullsp_vec); IBTK_CHKERRQ(ierr);
        }

        // Deallocate scratch data.
        U_rhs_vec->deallocateVectorData();
        P_rhs_vec->deallocateVectorData();
        nul_vec->deallocateVectorData();
#endif
        // Normalize p^{n+1/2} to have mean (discrete integral) zero.
        if (d_normalize_pressure)
        {
            const double P_mean = (1.0/d_volume)*d_hier_cc_data_ops->integral(d_P_new_idx, d_wgt_cc_idx);
            d_hier_cc_data_ops->addScalar(d_P_new_idx, d_P_new_idx, -P_mean);
        }

        // Reset the value the current estimate of the pressure when performing
        // multiple cycles.
        if (cycle < num_cycles-1)
        {
            d_hier_cc_data_ops->copyData(d_P_current_idx, d_P_new_idx);
        }
    }

    // Compute the cell-centered approximation to u^{n+1} (used for
    // visualization only).
    reinterpolateVelocity(getNewContext());

    // Compute the cell-centered approximation to f^{n+1} (used for
    // visualization only).
    reinterpolateForce(getNewContext());

    // Setup U_scratch to allow for ghost cell filling.
    d_hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_new_idx);
    d_U_scratch_bdry_fill_op->fillData(new_time);

    // Compute Omega = curl U.
    d_hier_math_ops->curl(
        d_Omega_new_idx, d_Omega_var,
        d_U_scratch_idx, d_U_var,
        d_no_fill_op, new_time);
#if (NDIM == 3)
    d_hier_math_ops->pointwise_L2Norm(
        d_Omega_Norm_new_idx, d_Omega_Norm_var,
        d_Omega_new_idx, d_Omega_var);
#endif

    // Compute max ||Omega||_2.
#if (NDIM == 2)
    d_Omega_max = std::max(+d_hier_cc_data_ops->max(d_Omega_new_idx),
                           -d_hier_cc_data_ops->min(d_Omega_new_idx));
#endif
#if (NDIM == 3)
    d_Omega_max = d_hier_cc_data_ops->max(d_Omega_Norm_new_idx);
#endif
    // Compute Div U.
    d_hier_math_ops->div(
        d_Div_U_new_idx, d_Div_U_var,
        1.0, d_U_scratch_idx, d_U_var, d_no_fill_op, new_time, false,
        0.0, -1, SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >(NULL),
        0, 0);
#if 0
    // Compute the kinetic energy of the fluid.
    SAMRAI::tbox::pout << "\nkinetic energy = " << 0.5*d_rho*d_hier_cc_data_ops->dot(d_U_cc_new_idx, d_U_cc_new_idx, d_wgt_cc_idx) << std::endl;
#endif
    t_integrate_hierarchy->stop();
    return getStableTimestep(getNewContext());
}// integrateHierarchy

void
INSStaggeredHierarchyIntegrator::synchronizeHierarchy()
{
    t_synchronize_hierarchy->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["SYNCH_NEW_STATE_DATA"][ln]->coarsenData();
    }

    t_synchronize_hierarchy->stop();
    return;
}// synchronizeHierarchy

void
INSStaggeredHierarchyIntegrator::synchronizeNewLevels(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level,
    const double sync_time,
    const bool initial_time)
{
    t_synchronize_new_levels->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0)
                && (coarsest_level < finest_level)
                && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    // Synchronize initial data on the hierarchy.
    if (finest_level > 0)
    {
        for (int ln = finest_level; ln > coarsest_level; --ln)
        {
            d_cscheds["SYNCH_CURRENT_STATE_DATA"][ln]->coarsenData();
        }
    }

    t_synchronize_new_levels->stop();
    return;
}// synchronizeNewLevels

void
INSStaggeredHierarchyIntegrator::resetTimeDependentHierData(
    const double new_time)
{
    t_reset_time_dependent_data->start();

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Advance the simulation time.
    d_old_dt = new_time - d_integrator_time;
    d_integrator_time = new_time;
    ++d_integrator_step;

    // Swap SAMRAI::hier::PatchData<NDIM> pointers between the current and new contexts.
    for (std::list<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > >::const_iterator sv =
             d_state_variables.begin();
         sv != d_state_variables.end(); ++sv)
    {
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& v = *sv;

        const int src_idx = var_db->mapVariableAndContextToIndex(
            v, getNewContext());
        const int dst_idx = var_db->mapVariableAndContextToIndex(
            v, getCurrentContext());

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > src_data =
                    patch->getPatchData(src_idx);
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > dst_data =
                    patch->getPatchData(dst_idx);

#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(src_data->getBox() == dst_data->getBox());
                TBOX_ASSERT(src_data->getGhostCellWidth() ==
                       dst_data->getGhostCellWidth());
#endif

                patch->setPatchData(dst_idx, src_data);
                patch->setPatchData(src_idx, dst_data);
            }
        }
    }

    // Deallocate the scratch and new data and reset the time of the current
    // data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->setTime(new_time, d_current_data);
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data);
    }

    t_reset_time_dependent_data->stop();
    return;
}// resetTimeDependentHierData

void
INSStaggeredHierarchyIntegrator::resetHierDataToPreadvanceState()
{
    t_reset_data_to_preadvance_state->start();

    // Deallocate the scratch and new data and reset the time of the current
    // data.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->setTime(d_integrator_time, d_current_data);
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data);
    }

    t_reset_data_to_preadvance_state->stop();
    return;
}// resetHierDataToPreadvanceState

///
///  The following routines:
///
///      initializeLevelData(),
///      resetHierarchyConfiguration(),
///      applyGradientDetector()
///
///  are concrete implementations of functions declared in the
///  mesh::StandardTagAndInitStrategy<NDIM> abstract base class.
///

void
INSStaggeredHierarchyIntegrator::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > base_hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > base_old_level,
    const bool allocate_data)
{
    t_initialize_level_data->start();

    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > old_level = base_old_level;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif
    // Allocate storage needed to initialize the level and fill data from
    // coarser levels in AMR hierarchy, if any.
    //
    // Since time gets set when we allocate data, re-stamp it to current time if
    // we don't need to allocate.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    if (allocate_data)
    {
        level->allocatePatchData(d_current_data, init_data_time);
    }
    else
    {
        level->setTime(init_data_time, d_current_data);
    }

    // Fill data from coarser levels in AMR hierarchy.
    if (level_number > 0 || !old_level.isNull())
    {
        level->allocatePatchData(d_scratch_data, init_data_time);

        IBTK::CartExtrapPhysBdryOp fill_after_regrid_bc_op(d_fill_after_regrid_bc_idxs, BDRY_EXTRAP_TYPE);
        d_fill_after_regrid->
            createSchedule(level,
                           old_level,
                           level_number-1,
                           hierarchy,
                           &fill_after_regrid_bc_op)->fillData(init_data_time);

        level->deallocatePatchData(d_scratch_data);
    }

    // Initialize level data at the initial time.
    if (initial_time)
    {
        // If no initialization object is provided, initialize the velocity,
        // divergance, and vorticity to zero.  Otherwise, use the initialization
        // object to set the velocity to some specified value and compute the
        // divergance and vorticity corresponding to the initial velocity.
        if (d_U_init.isNull())
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_current_data =
                    patch->getPatchData(d_U_current_idx);
                U_current_data->fillAll(0.0);

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_cc_current_data =
                    patch->getPatchData(d_U_cc_current_idx);
                U_cc_current_data->fillAll(0.0);

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_current_data =
                    patch->getPatchData(d_Omega_current_idx);
                Omega_current_data->fillAll(0.0);
#if (NDIM == 3)
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_Norm_current_data =
                    patch->getPatchData(d_Omega_Norm_current_idx);
                Omega_Norm_current_data->fillAll(0.0);
#endif
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Div_U_current_data =
                    patch->getPatchData(d_Div_U_current_idx);
                Div_U_current_data->fillAll(0.0);
            }
        }
        else
        {
            level->allocatePatchData(d_U_scratch_idx, init_data_time);

            // Initialize U.
            d_U_init->setDataOnPatchLevel(
                d_U_current_idx, d_U_var, level,
                init_data_time, initial_time);
            IBTK::PatchMathOps patch_math_ops;
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_current_data =
                    patch->getPatchData(d_U_current_idx);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_cc_current_data =
                    patch->getPatchData(d_U_cc_current_idx);

                patch_math_ops.interp(U_cc_current_data, U_current_data, patch);
            }

            // Fill in U boundary data from coarser levels.
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom =
                d_hierarchy->getGridGeometry();
            SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > ralg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
            SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator = grid_geom->lookupRefineOperator(
                d_U_var, "CONSERVATIVE_LINEAR_REFINE");
            ralg->registerRefine(d_U_scratch_idx, // destination
                                 d_U_current_idx, // source
                                 d_U_scratch_idx, // temporary work space
                                 refine_operator);
            IBTK::CartExtrapPhysBdryOp bc_op(d_U_scratch_idx, BDRY_EXTRAP_TYPE);
            ralg->createSchedule(level, level_number-1, hierarchy, &bc_op)->fillData(init_data_time);

            // Initialize quantities derived from the initial value of U.
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_scratch_data =
                    patch->getPatchData(d_U_scratch_idx);

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_current_data =
                    patch->getPatchData(d_Omega_current_idx);
                patch_math_ops.curl(Omega_current_data, U_scratch_data, patch);
#if (NDIM == 3)
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_Norm_current_data =
                    patch->getPatchData(d_Omega_Norm_current_idx);
                patch_math_ops.pointwise_L2Norm(Omega_Norm_current_data, Omega_current_data, patch);
#endif
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Div_U_current_data =
                    patch->getPatchData(d_Div_U_current_idx);
                patch_math_ops.div(Div_U_current_data, 1.0, U_scratch_data, 0.0, SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> >(NULL), patch);
            }

            level->deallocatePatchData(d_U_scratch_idx);
        }

        // Initialize the maximum value of ||Omega||_2 on the grid.
        if (level_number == 0)
        {
            d_Omega_max = 0.0;
        }

        SAMRAI::math::PatchCellDataOpsReal<NDIM,double> patch_cc_data_ops;
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
#if (NDIM == 2)
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_current_data =
                patch->getPatchData(d_Omega_current_idx);
            d_Omega_max = std::max(d_Omega_max, +patch_cc_data_ops.max(Omega_current_data, patch_box));
            d_Omega_max = std::max(d_Omega_max, -patch_cc_data_ops.min(Omega_current_data, patch_box));
#endif
#if (NDIM == 3)
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_Norm_current_data =
                patch->getPatchData(d_Omega_Norm_current_idx);
            d_Omega_max = std::max(d_Omega_max, patch_cc_data_ops.max(Omega_Norm_current_data, patch_box));
#endif
        }

        // If no initialization object is provided, initialize the pressure to
        // zero.  Otherwise, use the initialization object to set the pressure
        // to some specified value.
        //
        // NOTE: This initial value for the pressure IS NOT USED by the time
        // integrator and is only specified for purposes of visualization.
        if (d_P_init.isNull())
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > P_current_data =
                    patch->getPatchData(d_P_current_idx);
                P_current_data->fillAll(0.0);
            }
        }
        else
        {
            // Initialize P.
            d_P_init->setDataOnPatchLevel(
                d_P_current_idx, d_P_var, level,
                init_data_time, initial_time);
        }

        // If no initialization object is provided, initialize the body force to
        // zero.  Otherwise, use the initialization object to set the body force
        // to some specified value.
        if (d_F_set.isNull())
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > F_current_data =
                    patch->getPatchData(d_F_current_idx);
                F_current_data->fillAll(0.0);

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > F_cc_current_data =
                    patch->getPatchData(d_F_cc_current_idx);
                F_cc_current_data->fillAll(0.0);
            }
        }
        else
        {
            // Initialize F.
            d_F_set->setDataOnPatchLevel(
                d_F_current_idx, d_F_var, level,
                init_data_time, initial_time);
            IBTK::PatchMathOps patch_math_ops;
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > F_current_data =
                    patch->getPatchData(d_F_current_idx);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > F_cc_current_data =
                    patch->getPatchData(d_F_cc_current_idx);

                patch_math_ops.interp(F_cc_current_data, F_current_data, patch);
            }

        }
    }

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
INSStaggeredHierarchyIntegrator::resetHierarchyConfiguration(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    t_reset_hierarchy_configuration->start();

    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0)
                && (coarsest_level <= finest_level)
                && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // Reset the HierarchyProjector object.
    d_hier_projector->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // Reset the Hierarchy data operations for the new hierarchy configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);

    d_hier_sc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_sc_data_ops->resetLevels(0, finest_hier_level);

    // Reset the Hierarchy math operations for the new configuration.
    if (d_is_managing_hier_math_ops)
    {
        d_hier_math_ops->setPatchHierarchy(hierarchy);
        d_hier_math_ops->resetLevels(0, finest_hier_level);
    }

    // Get the control volume weight variables and patch data descriptor
    // indices.
    d_wgt_cc_var = d_hier_math_ops->getCellWeightVariable();
    d_wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

    d_wgt_sc_var = d_hier_math_ops->getSideWeightVariable();
    d_wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();

    // Get the volume of the physical domain.
    d_volume = d_hier_math_ops->getVolumeOfPhysicalDomain();

    // Setup the patch boundary filling objects.
    typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;

    InterpolationTransactionComponent U_scratch_component(d_U_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);  // XXXX
    d_U_scratch_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_U_scratch_bdry_fill_op->initializeOperatorState(U_scratch_component, d_hierarchy);

    InterpolationTransactionComponent P_scratch_component(d_P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);  // XXXX
    d_P_scratch_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_P_scratch_bdry_fill_op->initializeOperatorState(P_scratch_component, d_hierarchy);

    InterpolationTransactionComponent rhs_component(d_Grad_P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);  // XXXX
    d_rhs_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_rhs_bdry_fill_op->initializeOperatorState(rhs_component, d_hierarchy);

    // If we have added or removed a level, resize the schedule vectors.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin();
         it!= d_ralgs.end(); ++it)
    {
        d_rscheds[(*it).first].resize(finest_hier_level+1);
    }

    for (CoarsenAlgMap::const_iterator it = d_calgs.begin();
         it!= d_calgs.end(); ++it)
    {
        d_cscheds[(*it).first].resize(finest_hier_level+1);
    }

    // (Re)build refine communication schedules.  These are created for all
    // levels in the hierarchy.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin();
         it!= d_ralgs.end(); ++it)
    {
        for (int ln = coarsest_level; ln <= finest_hier_level; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

            d_rscheds[(*it).first][ln] = (*it).second->
                createSchedule(level, ln-1, hierarchy, d_rstrategies[(*it).first]);
        }
    }

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin();
         it!= d_calgs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level,1); ln <= finest_hier_level; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarser_level =
                hierarchy->getPatchLevel(ln-1);

            d_cscheds[(*it).first][ln] = (*it).second->
                createSchedule(coarser_level, level, d_cstrategies[(*it).first]);
        }
    }

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
INSStaggeredHierarchyIntegrator::applyGradientDetector(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    t_apply_gradient_detector->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // It is necessary to untag all cells prior to tagging.
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
        tags_data->fillAll(0);
    }

    // Tag cells based on the magnatude of the vorticity.
    //
    // Note that if either the relative or absolute threshold is zero for a
    // particular level, no tagging is performed on that level.
    if (d_using_vorticity_tagging)
    {
        const double Omega_rel_thresh =
            (level_number >= 0 && level_number < d_Omega_rel_thresh.getSize()
             ? d_Omega_rel_thresh[level_number]
             : (level_number < 0
                ? d_Omega_rel_thresh[0]
                : d_Omega_rel_thresh[d_Omega_rel_thresh.size()-1]));
        const double Omega_abs_thresh =
            (level_number >= 0 && level_number < d_Omega_abs_thresh.getSize()
             ? d_Omega_abs_thresh[level_number]
             : (level_number < 0
                ? d_Omega_abs_thresh[0]
                : d_Omega_abs_thresh[d_Omega_abs_thresh.size()-1]));
        if (Omega_rel_thresh > 0.0 && Omega_abs_thresh > 0.0)
        {
            const double thresh = sqrt(std::numeric_limits<double>::epsilon()) +
                std::min(Omega_rel_thresh*d_Omega_max, Omega_abs_thresh);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > tags_data = patch->
                    getPatchData(tag_index);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Omega_current_data = patch->
                    getPatchData(d_Omega_current_idx);

                for (SAMRAI::pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
                {
                    const SAMRAI::hier::Index<NDIM>& i = ic();
#if (NDIM == 2)
                    if (std::abs((*Omega_current_data)(i)) > thresh)
                    {
                        (*tags_data)(i) = 1;
                    }
#endif
#if (NDIM == 3)
                    double norm_Omega_sq = 0.0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        norm_Omega_sq += (*Omega_current_data)(i,d)*(*Omega_current_data)(i,d);
                    }
                    const double norm_Omega = sqrt(norm_Omega_sq);
                    if (norm_Omega > thresh)
                    {
                        (*tags_data)(i) = 1;
                    }
#endif
                }
            }
        }
    }

    t_apply_gradient_detector->stop();
    return;
}// applyGradientDetector

///
///  The following routines:
///
///      getVelocityVar(),
///      getPressureVar(),
///      getForceVar()
///
///  allows access to the various state variables maintained by the integrator.
///

SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >
INSStaggeredHierarchyIntegrator::getVelocityVar()
{
    return d_U_var;
}// getVelocityVar

SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >
INSStaggeredHierarchyIntegrator::getPressureVar()
{
    return d_P_var;
}// getPressureVar

SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >
INSStaggeredHierarchyIntegrator::getForceVar()
{
    return d_F_var;
}// getForceVar

///
///  The following routines:
///
///      getCurrentContext(),
///      getNewContext(),
///      getScratchContext()
///
///  allow access to the various variable contexts maintained by the integrator.
///

///
/// We simply reuse the SAMRAI::hier::VariableContext objects defined in the
/// AdvDiffIntegrator object.
///

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
INSStaggeredHierarchyIntegrator::getCurrentContext() const
{
    return d_current_context;
}// getCurrentContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
INSStaggeredHierarchyIntegrator::getNewContext() const
{
    return d_new_context;
}// getNewContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
INSStaggeredHierarchyIntegrator::getScratchContext() const
{
    return d_scratch_context;
}// getScratchContext

///
/// The following routines:
///
///      reinterpolateVelocity(),
///      reinterpolateForce()
///
/// are miscelaneous utility functions.

void
INSStaggeredHierarchyIntegrator::reinterpolateVelocity(
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx)
{
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    const int    U_idx = var_db->mapVariableAndContextToIndex(   d_U_var, ctx);
    const int U_cc_idx = var_db->mapVariableAndContextToIndex(d_U_cc_var, ctx);
    static const bool synch_cf_interface = true;
    d_hier_math_ops->interp(
        U_cc_idx, d_U_cc_var,
        U_idx   , d_U_var   ,
        d_no_fill_op, d_integrator_time, synch_cf_interface);
    return;
}// reinterpolateVelocity

void
INSStaggeredHierarchyIntegrator::reinterpolateForce(
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx)
{
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    const int    F_idx = var_db->mapVariableAndContextToIndex(   d_F_var, ctx);
    const int F_cc_idx = var_db->mapVariableAndContextToIndex(d_F_cc_var, ctx);
    static const bool synch_cf_interface = true;
    d_hier_math_ops->interp(
        F_cc_idx, d_F_cc_var,
        F_idx   , d_F_var   ,
        d_no_fill_op, d_integrator_time, synch_cf_interface);
    return;
}// reinterpolateForce

///
///  The following routines:
///
///      putToDatabase()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::tbox::Serializable abstract base class.
///

void
INSStaggeredHierarchyIntegrator::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("INS_STAGGERED_HIERARCHY_INTEGRATOR_VERSION",
                   INS_STAGGERED_HIERARCHY_INTEGRATOR_VERSION);
    db->putDouble("d_start_time", d_start_time);
    db->putDouble("d_end_time", d_end_time);
    db->putDouble("d_grow_dt", d_grow_dt);
    db->putInteger("d_max_integrator_steps", d_max_integrator_steps);
    db->putInteger("d_num_cycles", d_num_cycles);
    db->putInteger("d_num_init_cycles", d_num_init_cycles);
    db->putInteger("d_regrid_interval", d_regrid_interval);
    db->putBool("d_using_default_tag_buffer", d_using_default_tag_buffer);
    db->putIntegerArray("d_tag_buffer", d_tag_buffer);
    db->putBool("d_conservation_form", d_conservation_form);
    db->putBool("d_using_vorticity_tagging", d_using_vorticity_tagging);
    db->putDoubleArray("d_Omega_rel_thresh", d_Omega_rel_thresh);
    db->putDoubleArray("d_Omega_abs_thresh", d_Omega_abs_thresh);
    db->putDouble("d_Omega_max", d_Omega_max);
    db->putBool("d_second_order_pressure_update", d_second_order_pressure_update);
    db->putBool("d_normalize_pressure", d_normalize_pressure);
    db->putDouble("d_old_dt", d_old_dt);
    db->putDouble("d_integrator_time", d_integrator_time);
    db->putInteger("d_integrator_step", d_integrator_step);
    db->putDouble("d_dt_max", d_dt_max);
    db->putDouble("d_dt_max_time_max", d_dt_max_time_max);
    db->putDouble("d_dt_max_time_min", d_dt_max_time_min);
    db->putDouble("d_rho", d_rho);
    db->putDouble("d_mu", d_mu);
    db->putDouble("d_nu", d_nu);
    db->putDouble("d_lambda", d_lambda);

    t_put_to_database->stop();
    return;
}// putToDatabase

///
///  The following routines:
///
///      printClassData()
///
///  are provided for your viewing pleasure.
///

void
INSStaggeredHierarchyIntegrator::printClassData(
    std::ostream& os) const
{
    os << "\nINSStaggeredHierarchyIntegrator::printClassData..." << std::endl;
    os << "this = " << const_cast<INSStaggeredHierarchyIntegrator*>(this) << std::endl;
    os << "d_object_name = " << d_object_name << "\n"
       << "d_registered_for_restart = " << d_registered_for_restart << std::endl;
    os << "d_hierarchy = " << d_hierarchy.getPointer() << "\n"
       << "d_gridding_alg = " << d_gridding_alg.getPointer() << std::endl;
    os << "d_visit_writer = " << d_visit_writer.getPointer() << "\n"
       << "d_U_scale = " << d_U_scale << "\n"
       << "d_P_scale = " << d_P_scale << "\n"
       << "d_F_scale = " << d_F_scale << std::endl;
    os << "d_hier_projector = " << d_hier_projector.getPointer() << std::endl;
    os << "d_start_time = " << d_start_time << "\n"
       << "d_end_time = " << d_end_time << "\n"
       << "d_grow_dt = " << d_grow_dt << "\n"
       << "d_max_integrator_steps = " << d_max_integrator_steps << std::endl;
    os << "d_num_cycles = " << d_num_cycles << std::endl;
    os << "d_num_init_cycles = " << d_num_init_cycles << std::endl;
    os << "d_regrid_interval = " << d_regrid_interval << std::endl;
    os << "d_using_default_tag_buffer = " << d_using_default_tag_buffer << "\n"
       << "d_tag_buffer = [ ";
    std::copy(d_tag_buffer.getPointer(), d_tag_buffer.getPointer()+d_tag_buffer.size(), std::ostream_iterator<int>(os, " , "));
    os << " ]" << std::endl;
    os << "d_conservation_form = " << d_conservation_form << std::endl;
    os << "d_using_vorticity_tagging = " << d_using_vorticity_tagging << "\n"
       << "d_Omega_rel_thresh = [ ";
    std::copy(d_Omega_rel_thresh.getPointer(), d_Omega_rel_thresh.getPointer()+d_Omega_rel_thresh.size(), std::ostream_iterator<double>(os, " , "));
    os << " ]\n"
       << "d_Omega_abs_thresh = [ ";
    std::copy(d_Omega_abs_thresh.getPointer(), d_Omega_abs_thresh.getPointer()+d_Omega_abs_thresh.size(), std::ostream_iterator<double>(os, " , "));
    os << " ]\n"
       << "d_Omega_max = " << d_Omega_max << std::endl;
    os << "d_second_order_pressure_update = " << d_second_order_pressure_update << std::endl;
    os << "d_normalize_pressure = " << d_normalize_pressure << std::endl;
    os << "d_output_U = " << d_output_U << "\n"
       << "d_output_P = " << d_output_P << "\n"
       << "d_output_F = " << d_output_F << "\n"
       << "d_output_Omega = " << d_output_Omega << "\n"
       << "d_output_Div_U = " << d_output_Div_U << std::endl;
    os << "d_old_dt = " << d_old_dt << "\n"
       << "d_integrator_time = " << d_integrator_time << "\n"
       << "d_integrator_step = " << d_integrator_step << std::endl;
    os << "d_cfl = " << d_cfl << std::endl;
    os << "d_dt_max = " << d_dt_max << "\n"
       << "d_dt_max_time_max = " << d_dt_max_time_max << "\n"
       << "d_dt_max_time_min = " << d_dt_max_time_min << std::endl;
    os << "d_do_log = " << d_do_log << std::endl;
    os << "d_rho = " << d_rho << "\n"
       << "d_mu = " << d_mu << "\n"
       << "d_nu = " << d_nu << "\n"
       << "d_lambda = " << d_lambda << std::endl;
    os << "d_hier_cc_data_ops = " << d_hier_cc_data_ops.getPointer() << "\n"
       << "d_hier_sc_data_ops = " << d_hier_sc_data_ops.getPointer() << "=n"
       << "d_hier_math_ops = " << d_hier_math_ops.getPointer() << "\n"
       << "d_is_managing_hier_math_ops = " << d_is_managing_hier_math_ops << std::endl;
    os << "d_wgt_cc_var = " << d_wgt_cc_var.getPointer() << "\n"
       << "d_wgt_sc_var = " << d_wgt_sc_var.getPointer() << "\n"
       << "d_volume = " << d_volume << std::endl;
    os << "Skipping variables, patch data descriptors, communications algorithms, etc." << std::endl;
    return;
}// printClassData

/////////////////////////////// PROTECTED ////////////////////////////////////

void
INSStaggeredHierarchyIntegrator::registerVariable(
    int& current_idx,
    int& new_idx,
    int& scratch_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > variable,
    const SAMRAI::hier::IntVector<NDIM>& scratch_ghosts,
    const std::string& coarsen_name,
    const std::string& refine_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!variable.isNull());
#endif
    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    current_idx = -1; // insure that uninitialized variable patch data
    new_idx     = -1; // descriptor indices cause errors
    scratch_idx = -1;

    d_state_variables.push_back(variable);

    // Setup the current context.
    current_idx = var_db->registerVariableAndContext(
        variable, getCurrentContext(), no_ghosts);
    d_current_data.setFlag(current_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(current_idx);
    }

    // Setup the new context.
    new_idx = var_db->registerVariableAndContext(
        variable, getNewContext(), no_ghosts);
    d_new_data.setFlag(new_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(new_idx);
    }

    // Setup the scratch context.
    scratch_idx = var_db->registerVariableAndContext(
        variable, getScratchContext(), scratch_ghosts);
    d_scratch_data.setFlag(scratch_idx);

    // Setup the communication algorithms.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator =
        grid_geom->lookupRefineOperator(variable, refine_name);
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > coarsen_operator =
        grid_geom->lookupCoarsenOperator(variable, coarsen_name);

    // Setup the refine algorithm used to fill data in new or modified patch
    // levels following a regrid operation.
    if (!refine_operator.isNull())
    {
        d_fill_after_regrid->
            registerRefine(current_idx, // destination
                           current_idx, // source
                           scratch_idx, // temporary work space
                           refine_operator);

        // Keep track of the cell-centered scratch data indices, for use in the
        // refinement schedule used to fill data in new or modified patch levels
        // following a regrid operation.
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cc_var = variable;
        if (!cc_var.isNull())
        {
            d_fill_after_regrid_bc_idxs.setFlag(scratch_idx);
        }
    }

    // Setup the SYNCH_CURRENT_STATE_DATA and SYNCH_NEW_STATE_DATA algorithms,
    // used to synchronize the data on the hierarchy.
    if (!coarsen_operator.isNull())
    {
        d_calgs["SYNCH_CURRENT_STATE_DATA"]->
            registerCoarsen(current_idx, // destination
                            current_idx, // source
                            coarsen_operator);

        d_calgs["SYNCH_NEW_STATE_DATA"]->
            registerCoarsen(new_idx, // destination
                            new_idx, // source
                            coarsen_operator);
    }
    return;
}// registerVariable

void
INSStaggeredHierarchyIntegrator::registerVariable(
    int& scratch_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > variable,
    const SAMRAI::hier::IntVector<NDIM>& scratch_ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!variable.isNull());
#endif

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    scratch_idx = -1; // insure that uninitialized variable patch data
                      // descriptor indices cause errors

    d_scratch_variables.push_back(variable);

    // Setup the scratch context.
    scratch_idx = var_db->registerVariableAndContext(
        variable, getScratchContext(), scratch_ghosts);
    d_scratch_data.setFlag(scratch_idx);
    return;
}// registerVariable

/////////////////////////////// PRIVATE //////////////////////////////////////

void
INSStaggeredHierarchyIntegrator::computeConvectiveDerivative(
    const int& N_idx,
    const int& U_idx,
    const int& F_idx,
    const int& Grad_P_idx,
    const double current_time,
    const double new_time,
    const bool conservation_form)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Allocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_gadvect_U_scratch_idx, current_time);
        level->allocatePatchData(d_gadvect_F_scratch_idx, current_time);
    }

    // Setup communications schedules and fill velocity boundary data.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator;

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > ralg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    refine_operator = grid_geom->lookupRefineOperator(
        d_gadvect_U_var, "CONSERVATIVE_LINEAR_REFINE");
    ralg->registerRefine(d_gadvect_U_scratch_idx, // destination
                         U_idx,                   // source
                         d_gadvect_U_scratch_idx, // temporary work space
                         refine_operator);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        static const std::string ralg_name = "gadvect_U_scratch_bdry_fill";
        ralg->resetSchedule(d_rscheds[ralg_name][ln]);
        d_rscheds[ralg_name][ln]->fillData(current_time);
        d_ralgs[ralg_name]->resetSchedule(d_rscheds[ralg_name][ln]);
    }

    // Setup right-hand-side terms for velocity extrapolation.
    SAMRAI::solv::PoissonSpecifications gadvect_rhs_spec(d_object_name+"::gadvect_rhs_spec");
    gadvect_rhs_spec.setCConstant(-d_lambda/d_rho);
    gadvect_rhs_spec.setDConstant(+d_mu    /d_rho);
    d_hier_math_ops->laplace(
        d_gadvect_F_scratch_idx, d_gadvect_F_var,
        gadvect_rhs_spec,
        d_gadvect_U_scratch_idx, d_gadvect_U_var,
        d_no_fill_op, current_time,
        -1.0/d_rho,
        Grad_P_idx, d_Grad_P_var);
    d_hier_sc_data_ops->add(d_gadvect_F_scratch_idx, d_gadvect_F_scratch_idx, F_idx);

    // Fill right-hand-side boundary conditions.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        static const std::string ralg_name = "gadvect_F_scratch_bdry_fill";
        d_rscheds[ralg_name][ln]->fillData(current_time);
    }

#if (NDIM == 2)
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

            const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            static const SAMRAI::hier::IntVector<NDIM> growth_factor = NUM_GODUNOV_CYCLES-1;

            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::hier::Box<NDIM>& grown_patch_box = SAMRAI::hier::Box<NDIM>::grow(patch_box,growth_factor);
            const SAMRAI::hier::IntVector<NDIM>& grown_patch_lower = grown_patch_box.lower();
            const SAMRAI::hier::IntVector<NDIM>& grown_patch_upper = grown_patch_box.upper();

            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > N_data = patch->getPatchData(N_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_data = patch->getPatchData(d_gadvect_U_scratch_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > F_data = patch->getPatchData(d_gadvect_F_scratch_idx);

            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > N_grown_data = new SAMRAI::pdat::SideData<NDIM,double>(grown_patch_box,N_data->getDepth(),N_data->getGhostCellWidth()-growth_factor);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_grown_data = new SAMRAI::pdat::SideData<NDIM,double>(grown_patch_box,U_data->getDepth(),U_data->getGhostCellWidth()-growth_factor);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > F_grown_data = new SAMRAI::pdat::SideData<NDIM,double>(grown_patch_box,F_data->getDepth(),F_data->getGhostCellWidth()-growth_factor);
            U_grown_data->copy(*U_data);
            F_grown_data->copy(*F_data);

            const SAMRAI::hier::IntVector<NDIM> ghosts = SAMRAI::hier::IntVector<NDIM>(GADVECTG)-growth_factor;
            SAMRAI::hier::Box<NDIM> side_boxes[NDIM];
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> >  U_adv_data[NDIM];
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > U_half_data[NDIM];
            for (int axis = 0; axis < NDIM; ++axis)
            {
                side_boxes [axis] = SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(grown_patch_box,axis);
                U_adv_data [axis] = new SAMRAI::pdat::FaceData<NDIM,double>(side_boxes[axis],1,ghosts);
                U_half_data[axis] = new SAMRAI::pdat::FaceData<NDIM,double>(side_boxes[axis],1,ghosts);
            }

            NAVIER_STOKES_INTERP_COMPS_F77(
                grown_patch_lower(0), grown_patch_upper(0),
                grown_patch_lower(1), grown_patch_upper(1),
                U_grown_data->getGhostCellWidth()(0),   U_grown_data->getGhostCellWidth()(1),
                U_grown_data->getPointer(0),            U_grown_data->getPointer(1),
                side_boxes[0].lower(0),                 side_boxes[0].upper(0),
                side_boxes[0].lower(1),                 side_boxes[0].upper(1),
                U_adv_data[0]->getGhostCellWidth()(0),  U_adv_data[0]->getGhostCellWidth()(1),
                U_adv_data[0]->getPointer(0),           U_adv_data[0]->getPointer(1),
                side_boxes[1].lower(0),                 side_boxes[1].upper(0),
                side_boxes[1].lower(1),                 side_boxes[1].upper(1),
                U_adv_data[1]->getGhostCellWidth()(0),  U_adv_data[1]->getGhostCellWidth()(1),
                U_adv_data[1]->getPointer(0),           U_adv_data[1]->getPointer(1));

            for (int count = 0; count < NUM_GODUNOV_CYCLES; ++count)
            {
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > dU_grown_data =
                        new SAMRAI::pdat::SideData<NDIM,double>(U_grown_data->getBox(), U_grown_data->getDepth(), U_grown_data->getGhostCellWidth());
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_L_grown_data =
                        new SAMRAI::pdat::SideData<NDIM,double>(U_grown_data->getBox(), U_grown_data->getDepth(), U_grown_data->getGhostCellWidth());
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_R_grown_data =
                        new SAMRAI::pdat::SideData<NDIM,double>(U_grown_data->getBox(), U_grown_data->getDepth(), U_grown_data->getGhostCellWidth());
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_grown_scratch_data =
                        new SAMRAI::pdat::SideData<NDIM,double>(U_grown_data->getBox(), U_grown_data->getDepth(), U_grown_data->getGhostCellWidth());
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > F_grown_scratch_data =
                        new SAMRAI::pdat::SideData<NDIM,double>(F_grown_data->getBox(), F_grown_data->getDepth(), F_grown_data->getGhostCellWidth());
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > U_half_scratch_data =
                        new SAMRAI::pdat::FaceData<NDIM,double>( U_half_data[axis]->getBox(), U_half_data[axis]->getDepth(), U_half_data[axis]->getGhostCellWidth());
                    GODUNOV_PREDICT_WITH_SOURCE_F77(
                        dx, dt,
                        side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                        side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                        U_grown_data->getGhostCellWidth()(0), U_grown_data->getGhostCellWidth()(1),
                        F_grown_data->getGhostCellWidth()(0), F_grown_data->getGhostCellWidth()(1),
                        U_grown_data->getPointer(axis), dU_grown_data->getPointer(axis), U_L_grown_data->getPointer(axis), U_R_grown_data->getPointer(axis), U_grown_scratch_data->getPointer(axis),
                        F_grown_data       ->getPointer(axis),       F_grown_scratch_data->getPointer(axis),
                        U_adv_data [axis]  ->getGhostCellWidth()(0), U_adv_data [axis]   ->getGhostCellWidth()(1),
                        U_half_data[axis]  ->getGhostCellWidth()(0), U_half_data[axis]   ->getGhostCellWidth()(1),
                        U_adv_data [axis]  ->getPointer(0),          U_adv_data [axis]   ->getPointer(1),
                        U_half_scratch_data->getPointer(0),          U_half_scratch_data ->getPointer(1),
                        U_half_data[axis]  ->getPointer(0),          U_half_data[axis]   ->getPointer(1));
                }

                NAVIER_STOKES_RESET_ADV_VELOCITY_F77(
                    side_boxes[0].lower(0), side_boxes[0].upper(0),
                    side_boxes[0].lower(1), side_boxes[0].upper(1),
                    U_adv_data [0]->getGhostCellWidth()(0), U_adv_data [0]->getGhostCellWidth()(1),
                    U_adv_data [0]->getPointer(0),          U_adv_data [0]->getPointer(1),
                    U_half_data[0]->getGhostCellWidth()(0), U_half_data[0]->getGhostCellWidth()(1),
                    U_half_data[0]->getPointer(0),          U_half_data[0]->getPointer(1),
                    side_boxes[1].lower(0), side_boxes[1].upper(0),
                    side_boxes[1].lower(1), side_boxes[1].upper(1),
                    U_adv_data [1]->getGhostCellWidth()(0), U_adv_data [1]->getGhostCellWidth()(1),
                    U_adv_data [1]->getPointer(0),          U_adv_data [1]->getPointer(1),
                    U_half_data[1]->getGhostCellWidth()(0), U_half_data[1]->getGhostCellWidth()(1),
                    U_half_data[1]->getPointer(0),          U_half_data[1]->getPointer(1));
            }

            for (int axis = 0; axis < NDIM; ++axis)
            {
                ADVECT_DERIVATIVE_F77(
                    dx,
                    side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                    side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                    U_adv_data [axis]->getGhostCellWidth()(0), U_adv_data [axis]->getGhostCellWidth()(1),
                    U_half_data[axis]->getGhostCellWidth()(0), U_half_data[axis]->getGhostCellWidth()(1),
                    U_adv_data [axis]->getPointer(0),          U_adv_data [axis]->getPointer(1),
                    U_half_data[axis]->getPointer(0),          U_half_data[axis]->getPointer(1),
                    N_grown_data->getGhostCellWidth()(0), N_grown_data->getGhostCellWidth()(1),
                    N_grown_data->getPointer(axis));
            }
            N_data->copy(*N_grown_data);
        }
    }
#endif

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_gadvect_U_scratch_idx);
        level->deallocatePatchData(d_gadvect_F_scratch_idx);
    }
    return;
}// computeConvectiveDerivative

double
INSStaggeredHierarchyIntegrator::getLevelDt(
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx) const
{
    double stable_dt = std::numeric_limits<double>::max();
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        stable_dt = std::min(stable_dt,getPatchDt(patch,ctx));
    }
    return stable_dt;
}// getLevelDt

double
INSStaggeredHierarchyIntegrator::getPatchDt(
    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx) const
{
    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
        patch->getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const SAMRAI::hier::Index<NDIM>& ilower = patch->getBox().lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch->getBox().upper();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_data =
        patch->getPatchData(d_U_var, ctx);
    const SAMRAI::hier::IntVector<NDIM>& U_ghost_cells = U_data->getGhostCellWidth();

    double stable_dt = std::numeric_limits<double>::max();
#if (NDIM == 2)
    NAVIER_STOKES_SC_STABLEDT_F77(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),
        U_ghost_cells(0),U_ghost_cells(1),
        U_data->getPointer(0),U_data->getPointer(1),
        stable_dt);
#endif
#if (NDIM == 3)
    NAVIER_STOKES_SC_STABLEDT_F77(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
        U_ghost_cells(0),U_ghost_cells(1),U_ghost_cells(2),
        U_data->getPointer(0),U_data->getPointer(1),U_data->getPointer(2),
        stable_dt);
#endif
    return stable_dt;
}// getPatchDt

void
INSStaggeredHierarchyIntegrator::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif

    // Read in data members from input database.
    d_end_time = db->getDoubleWithDefault("end_time", d_end_time);
    d_grow_dt = db->getDoubleWithDefault("grow_dt", d_grow_dt);

    d_max_integrator_steps = db->getIntegerWithDefault(
        "max_integrator_steps", d_max_integrator_steps);

    d_num_cycles = db->getIntegerWithDefault("num_cycles", d_num_cycles);

    d_regrid_interval = db->getIntegerWithDefault(
        "regrid_interval", d_regrid_interval);

    if (db->keyExists("tag_buffer"))
    {
        d_tag_buffer = db->getIntegerArray("tag_buffer");
        d_using_default_tag_buffer = false;
    }
    else
    {
        d_using_default_tag_buffer = true;
        TBOX_WARNING(d_object_name << ":  "
                     << "Key data `tag_buffer' not found in input.  "
                     << "Default values used.  See class header for details.");
    }

    d_conservation_form = db->getBoolWithDefault(
        "conservation_form", d_conservation_form);

    d_using_vorticity_tagging = db->getBoolWithDefault(
        "using_vorticity_tagging", d_using_vorticity_tagging);

    if (d_using_vorticity_tagging)
    {
        if (db->keyExists("vorticity_rel_thresh"))
        {
            d_Omega_rel_thresh = db->getDoubleArray("vorticity_rel_thresh");
        }
        else
        {
            TBOX_WARNING(d_object_name << ":\n"
                         << "  Vorticity tagging is enabled but key data `vorticity_rel_thresh' not found in input.\n"
                         << "  Using default values for all levels in the locally refined grid.\n");
            d_Omega_rel_thresh.resizeArray(1);
            d_Omega_rel_thresh[0] = 1.0;
        }

        for (int i = 0; i < d_Omega_rel_thresh.getSize(); ++i)
        {
            if (d_Omega_rel_thresh[i] < 0.0 || d_Omega_rel_thresh[i] > 1.0)
            {
                TBOX_ERROR(d_object_name << ":  "
                           << "relative vorticity thresholds for each level must lie in the interval [0,1].\n");
            }
        }

        if (db->keyExists("vorticity_abs_thresh"))
        {
            d_Omega_abs_thresh = db->getDoubleArray("vorticity_abs_thresh");
        }
        else
        {
            TBOX_WARNING(d_object_name << ":\n"
                         << "  Vorticity tagging is enabled but key data `vorticity_abs_thresh' not found in input.\n"
                         << "  Using default values for all levels in the locally refined grid.\n");
            d_Omega_abs_thresh.resizeArray(1);
            d_Omega_abs_thresh[0] = std::numeric_limits<double>::max();
        }

        for (int i = 0; i < d_Omega_abs_thresh.getSize(); ++i)
        {
            if (d_Omega_abs_thresh[i] < 0.0)
            {
                TBOX_ERROR(d_object_name << ":  "
                           << "absolute vorticity thresholds for each level must be nonnegative.\n");
            }
        }
    }

    d_output_U      = db->getBoolWithDefault("output_U"     , d_output_U     );
    d_output_P      = db->getBoolWithDefault("output_P"     , d_output_P     );
    d_output_F      = db->getBoolWithDefault("output_F"     , d_output_F     );
    d_output_Omega  = db->getBoolWithDefault("output_Omega" , d_output_Omega );
    d_output_Div_U  = db->getBoolWithDefault("output_Div_U" , d_output_Div_U );

    d_U_scale = db->getDoubleWithDefault("U_scale", d_U_scale);
    d_P_scale = db->getDoubleWithDefault("P_scale", d_P_scale);
    d_F_scale = db->getDoubleWithDefault("F_scale", d_F_scale);

    d_cfl = db->getDoubleWithDefault("cfl",d_cfl);

    d_dt_max = db->getDoubleWithDefault("dt_max",d_dt_max);
    d_dt_max_time_max = db->getDoubleWithDefault(
        "dt_max_time_max", d_dt_max_time_max);
    d_dt_max_time_min = db->getDoubleWithDefault(
        "dt_max_time_min", d_dt_max_time_min);

    d_do_log = db->getBoolWithDefault("enable_logging", d_do_log);

    if (!is_from_restart)
    {
        d_start_time = db->getDoubleWithDefault("start_time", d_start_time);

        d_num_init_cycles = db->getIntegerWithDefault(
            "num_init_cycles", d_num_init_cycles);

        d_second_order_pressure_update = db->getBoolWithDefault(
            "second_order_pressure_update", d_second_order_pressure_update);

        d_normalize_pressure = db->getBoolWithDefault(
            "normalize_pressure", d_normalize_pressure);

        if (db->keyExists("rho"))
        {
            d_rho = db->getDouble("rho");
        }
        else
        {
            TBOX_ERROR(d_object_name << ":  "
                       << "Key data `rho' not found in input.");
        }

        if (db->keyExists("mu"))
        {
            d_mu = db->getDouble("mu");
        }
        else
        {
            TBOX_ERROR(d_object_name << ":  "
                       << "Key data `mu' not found in input.");
        }

        d_nu = d_mu/d_rho;  // the kinematic viscosity

        if (db->keyExists("lambda"))
        {
            d_lambda = db->getDouble("lambda");
        }
        else
        {
            d_lambda = 0.0;
        }

        d_lambda = d_lambda/d_rho;
    }
    return;
}// getFromInput

void
INSStaggeredHierarchyIntegrator::getFromRestart()
{
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> restart_db =
        SAMRAI::tbox::RestartManager::getManager()->getRootDatabase();

    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to "
                   << d_object_name << " not found in restart file.");
    }

    int ver = db->getInteger("INS_STAGGERED_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != INS_STAGGERED_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different than class version.");
    }

    d_start_time = db->getDouble("d_start_time");
    d_end_time = db->getDouble("d_end_time");
    d_grow_dt = db->getDouble("d_grow_dt");
    d_max_integrator_steps = db->getInteger("d_max_integrator_steps");
    d_num_cycles = db->getInteger("d_num_cycles");
    d_num_init_cycles = db->getInteger("d_num_init_cycles");
    d_regrid_interval = db->getInteger("d_regrid_interval");
    d_using_default_tag_buffer = db->getBool("d_using_default_tag_buffer");
    d_tag_buffer = db->getIntegerArray("d_tag_buffer");
    d_conservation_form = db->getBool("d_conservation_form");
    d_using_vorticity_tagging = db->getBool("d_using_vorticity_tagging");
    d_Omega_rel_thresh = db->getDoubleArray("d_Omega_rel_thresh");
    d_Omega_abs_thresh = db->getDoubleArray("d_Omega_abs_thresh");
    d_Omega_max = db->getDouble("d_Omega_max");
    d_second_order_pressure_update = db->getBool("d_second_order_pressure_update");
    d_normalize_pressure = db->getBool("d_normalize_pressure");
    d_old_dt = db->getDouble("d_old_dt");
    d_integrator_time = db->getDouble("d_integrator_time");
    d_integrator_step = db->getInteger("d_integrator_step");
    d_dt_max = db->getDouble("d_dt_max");
    d_dt_max_time_max = db->getDouble("d_dt_max_time_max");
    d_dt_max_time_min = db->getDouble("d_dt_max_time_min");
    d_rho = db->getDouble("d_rho");
    d_mu = db->getDouble("d_mu");
    d_nu = db->getDouble("d_nu");
    d_lambda = db->getDouble("d_lambda");
    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::INSStaggeredHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
