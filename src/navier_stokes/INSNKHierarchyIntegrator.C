// Filename: INSNKHierarchyIntegrator.C
// Last modified: <21.Mar.2008 19:09:32 griffith@box221.cims.nyu.edu>
// Created on 20 Mar 2008 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "INSNKHierarchyIntegrator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/FACPreconditionerLSWrapper.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScSAMRAIVectorReal.h>

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
#define NAVIER_STOKES_STABLEDT_F77 F77_FUNC_(navier_stokes_stabledt2d, NAVIER_STOKES_STABLEDT2D)
#define WENO_CONVECTIVE_FLUXES_F77 F77_FUNC_(weno_convective_fluxes2d, WENO_CONVECTIVE_FLUXES2D)
#endif

#if (NDIM == 3)
#define NAVIER_STOKES_STABLEDT_F77 F77_FUNC_(navier_stokes_stabledt3d, NAVIER_STOKES_STABLEDT3D)
#endif

extern "C"
{
    void
    NAVIER_STOKES_STABLEDT_F77(
        const double*,
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

#if (NDIM == 2)
    void
    WENO_CONVECTIVE_FLUXES_F77(
        double* , double* , const int& ,
        double* , double* , const int& ,
        double* , double* , const int& ,
        double* , double* , const int& ,
        double* , double* , const int& ,
        double* , const int& ,
        const double* , const int& ,
        const double* , const int&,
        const int& , const int& ,
        const int& , const int&
                               );
#endif
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
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_integrate_hierarchy_1st_order;
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
static const int FACEG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

static const int WENO_CELLG = 4;
static const int WENO_FACEG = 4;

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_COARSEN_TYPE = "CONSERVATIVE_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Version of INSNKHierarchyIntegrator restart file data.
static const int INS_NK_HIERARCHY_INTEGRATOR_VERSION = 1;

// Linear operator used with the Krylov solver.
class INSNKOperator
    : public IBTK::LinearOperator
{
public:
    /*!
     * XXXX
     */
    INSNKOperator(
        const std::string& object_name,
        const double rho,
        const double mu,
        const double lambda,
        const double dt,
        const double apply_time,
        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
        SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> U_P_bdry_fill_op)
        : d_object_name(object_name),
          d_rho(rho),
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
        }// INSNKOperator

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~INSNKOperator()
        {
            // intentionally blank
            return;
        }// ~INSNKOperator

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
            // Setup the boundary condition objects.
            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* U_bc_coef = NULL;  // XXXX
            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef = NULL;  // XXXX

            // Get the vector components.
            const int U_in_idx = x.getComponentDescriptorIndex(0);
            const int P_in_idx = x.getComponentDescriptorIndex(1);

            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_in_var = x.getComponentVariable(0);
            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_in_var = x.getComponentVariable(1);

            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > U_in_cc_var = U_in_var;
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_in_cc_var = P_in_var;

            const int U_out_idx = y.getComponentDescriptorIndex(0);
            const int P_out_idx = y.getComponentDescriptorIndex(1);

            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_out_var = y.getComponentVariable(0);
            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_out_var = y.getComponentVariable(1);

            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > U_out_cc_var = U_out_var;
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_out_cc_var = P_out_var;

            // Reset the interpolation operators and fill the data.
            //
            // XXXX: This will not work correctly for inhomogeneous boundary
            // conditions.
            typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
            InterpolationTransactionComponent U_component(U_in_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, U_bc_coef);
            InterpolationTransactionComponent P_component(P_in_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, P_bc_coef);

            std::vector<InterpolationTransactionComponent> transaction_comps(2);
            transaction_comps[0] = U_component;
            transaction_comps[1] = P_component;

            d_U_P_bdry_fill_op->resetTransactionComponents(transaction_comps);
            d_U_P_bdry_fill_op->fillData(d_apply_time);

            // Compute the action of the operator.
            SAMRAI::solv::PoissonSpecifications helmholtz_spec(d_object_name+"::helmholtz_spec");
            helmholtz_spec.setCConstant((d_rho/d_dt)+d_lambda);
            helmholtz_spec.setDConstant(-d_mu);

            for (int d = 0; d < NDIM; ++d)
            {
                d_hier_math_ops->laplace(
                    U_out_idx, U_out_cc_var,
                    helmholtz_spec, U_in_idx, U_in_cc_var, d_no_fill_op, d_apply_time,
                    0.0, -1, SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >(NULL),
                    d, d);
            }

            d_hier_math_ops->grad(
                U_out_idx, U_out_cc_var,
                1.0, P_in_idx, P_in_cc_var, d_no_fill_op, d_apply_time,
                1.0, U_out_idx, U_out_cc_var);

            d_hier_math_ops->div(
                P_out_idx, P_out_cc_var,
                1.0, U_in_idx, U_in_cc_var, d_no_fill_op, d_apply_time,
                0.0, -1, SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >(NULL));
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
    INSNKOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSNKOperator(
        const INSNKOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSNKOperator&
    operator=(
        const INSNKOperator& that);

    // Housekeeping.
    std::string d_object_name;

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
class INSNKPreconditioner
    : public IBTK::LinearSolver
{
public:
    /*!
     * XXXX
     */
    INSNKPreconditioner(
        const std::string& object_name,
        const double rho,
        const double mu,
        const double lambda,
        const double dt,
        const double apply_time,
        SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver> poisson_solver,
        SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver> helmholtz_solver,
        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
        SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> U_bdry_fill_op,
        SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> P_bdry_fill_op)
        : d_object_name(object_name),
          d_rho(rho),
          d_mu(mu),
          d_lambda(lambda),
          d_dt(dt),
          d_apply_time(apply_time),
          d_poisson_solver(poisson_solver),
          d_helmholtz_solver(helmholtz_solver),
          d_hier_math_ops(hier_math_ops),
          d_U_bdry_fill_op(U_bdry_fill_op),
          d_P_bdry_fill_op(P_bdry_fill_op),
          d_no_fill_op(SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation>(NULL))
        {
            // intentionally blank
            return;
        }// INSNKPreconditioner

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~INSNKPreconditioner()
        {
            // intentionally blank
            return;
        }// ~INSNKPreconditioner

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
            // Set the initial guess to equal x.
            y.copyVector(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >(&x,false));

            // Setup the boundary condition objects.
            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* U_bc_coef = NULL;  // XXXX
            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef = NULL;  // XXXX

            // Get the vector patch hierarchy.
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
            const int coarsest_ln = x.getCoarsestLevelNumber();
            const int finest_ln = x.getFinestLevelNumber();

            // Get the vector components.
            const int U_in_idx = x.getComponentDescriptorIndex(0);
            const int P_in_idx = x.getComponentDescriptorIndex(1);

            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_in_var = x.getComponentVariable(0);
            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_in_var = x.getComponentVariable(1);

            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > U_in_cc_var = U_in_var;
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_in_cc_var = P_in_var;

            const int U_out_idx = y.getComponentDescriptorIndex(0);
            const int P_out_idx = y.getComponentDescriptorIndex(1);

            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_out_var = y.getComponentVariable(0);
            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_out_var = y.getComponentVariable(1);

            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > U_out_cc_var = U_out_var;
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_out_cc_var = P_out_var;

            // Compute the action of the operator.
            SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > helmholtz_sol_vec =
                new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
                    d_object_name+"::helmholtz_sol_vec", hierarchy, coarsest_ln, finest_ln);
            helmholtz_sol_vec->addComponent(U_out_var,U_out_idx,y.getControlVolumeIndex(0));

            SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > helmholtz_rhs_vec =
                helmholtz_sol_vec->cloneVector(d_object_name+"::helmholtz_rhs_vec");
            const int helmholtz_rhs_idx = helmholtz_rhs_vec->getComponentDescriptorIndex(0);
            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& helmholtz_rhs_var = helmholtz_rhs_vec->getComponentVariable(0);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > helmholtz_rhs_cc_var = helmholtz_rhs_var;
            helmholtz_rhs_vec->allocateVectorData(d_apply_time);
            helmholtz_rhs_vec->copyVector(helmholtz_sol_vec);
            helmholtz_rhs_vec->scale(d_rho/d_dt, helmholtz_rhs_vec);

            typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
            InterpolationTransactionComponent P_in_component(P_in_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, P_bc_coef);
            d_P_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
            d_P_bdry_fill_op->initializeOperatorState(P_in_component, hierarchy);
            d_P_bdry_fill_op->fillData(d_apply_time);

            d_hier_math_ops->grad(
                helmholtz_rhs_idx, helmholtz_rhs_cc_var,
                -1.0, P_in_idx, P_in_cc_var, d_no_fill_op, d_apply_time,
                1.0, helmholtz_rhs_idx, helmholtz_rhs_cc_var);

            d_helmholtz_solver->setInitialGuessNonzero(true);
            d_helmholtz_solver->initializeSolverState(*helmholtz_sol_vec,*helmholtz_rhs_vec);
            d_helmholtz_solver->solveSystem(*helmholtz_sol_vec,*helmholtz_rhs_vec);
            helmholtz_rhs_vec->deallocateVectorData();

            InterpolationTransactionComponent U_out_component(U_out_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, U_bc_coef);
            d_U_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
            d_U_bdry_fill_op->initializeOperatorState(U_out_component, hierarchy);
            d_U_bdry_fill_op->fillData(d_apply_time);

            SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > poisson_sol_vec =
                new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
                    d_object_name+"::poisson_sol_vec", hierarchy, coarsest_ln, finest_ln);
            poisson_sol_vec->addComponent(P_out_var,P_out_idx,y.getControlVolumeIndex(1));

            SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > poisson_rhs_vec =
                poisson_sol_vec->cloneVector(d_object_name+"::poisson_rhs_vec");
            const int poisson_rhs_idx = poisson_rhs_vec->getComponentDescriptorIndex(0);
            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& poisson_rhs_var = poisson_rhs_vec->getComponentVariable(0);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > poisson_rhs_cc_var = poisson_rhs_var;
            poisson_rhs_vec->allocateVectorData(d_apply_time);

            d_hier_math_ops->div(
                poisson_rhs_idx, poisson_rhs_cc_var,
                -d_rho/d_dt, U_out_idx, U_out_cc_var, d_no_fill_op, d_apply_time,
                0.0, -1, SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >(NULL));

            d_poisson_solver->setInitialGuessNonzero(false);
            d_poisson_solver->initializeSolverState(*poisson_sol_vec,*poisson_rhs_vec);
            d_poisson_solver->solveSystem(*poisson_sol_vec,*poisson_rhs_vec);
            poisson_rhs_vec->deallocateVectorData();

            InterpolationTransactionComponent P_out_component(P_out_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, P_bc_coef);
            d_P_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
            d_P_bdry_fill_op->initializeOperatorState(P_out_component, hierarchy);
            d_P_bdry_fill_op->fillData(d_apply_time);

            d_hier_math_ops->grad(
                U_out_idx, U_out_cc_var,
                -d_dt/d_rho, P_out_idx, P_out_cc_var, d_no_fill_op, d_apply_time,
                1.0, U_out_idx, U_out_cc_var);
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
    INSNKPreconditioner();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSNKPreconditioner(
        const INSNKPreconditioner& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSNKPreconditioner&
    operator=(
        const INSNKPreconditioner& that);

    // Housekeeping.
    std::string d_object_name;

    // Problem coefficients.
    const double d_rho;
    const double d_mu;
    const double d_lambda;

    // The timestep size.
    const double d_dt;

    // The simulation time.
    const double d_apply_time;

    // Linear solvers.
    SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver> d_poisson_solver, d_helmholtz_solver;

    // Math objects.
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;

    // Boundary condition objects.
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_U_bdry_fill_op, d_P_bdry_fill_op, d_no_fill_op;
};

}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSNKHierarchyIntegrator::INSNKHierarchyIntegrator(
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

    d_normalize_pressure = false;

    d_regrid_interval = 1;
    d_old_dt = -1.0;
    d_integrator_time = std::numeric_limits<double>::quiet_NaN();
    d_integrator_step = std::numeric_limits<int>::max();

    d_output_U = false;
    d_output_P = false;
    d_output_F = false;
    d_output_Omega = false;
    d_output_Div_U = false;

    d_rho    = std::numeric_limits<double>::quiet_NaN();
    d_mu     = std::numeric_limits<double>::quiet_NaN();
    d_nu     = std::numeric_limits<double>::quiet_NaN();
    d_lambda = std::numeric_limits<double>::quiet_NaN();

    d_cfl = 0.5;

    d_dt_max = std::numeric_limits<double>::max();
    d_dt_max_time_max = std::numeric_limits<double>::max();
    d_dt_max_time_min = -(d_dt_max_time_max-std::numeric_limits<double>::epsilon());

    d_is_initialized = false;

    d_do_log = false;

    d_petsc_nullsp = static_cast<MatNullSpace>(NULL);
    d_petsc_nullsp_vec = static_cast<Vec>(NULL);

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, from_restart);

    // Obtain the Hierarchy data operations objects.
    SAMRAI::math::HierarchyDataOpsManager<NDIM>* hier_ops_manager =
        SAMRAI::math::HierarchyDataOpsManager<NDIM>::getManager();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cc_var =
        new SAMRAI::pdat::CellVariable<NDIM,double>("cc_var");
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > fc_var =
        new SAMRAI::pdat::FaceVariable<NDIM,double>("fc_var");

    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, hierarchy);
    d_hier_fc_data_ops = hier_ops_manager->getOperationsDouble(fc_var, hierarchy);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_initialize_hierarchy_integrator = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSNKHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSNKHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSNKHierarchyIntegrator::advanceHierarchy()");
        t_get_stable_timestep = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSNKHierarchyIntegrator::getStableTimestep()");
        t_regrid_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSNKHierarchyIntegrator::regridHierarchy()");
        t_integrate_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSNKHierarchyIntegrator::integrateHierarchy()");
        t_integrate_hierarchy_1st_order = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSNKHierarchyIntegrator::integrateHierarchy_1st_order()");
        t_synchronize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSNKHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSNKHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSNKHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_data_to_preadvance_state = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSNKHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSNKHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSNKHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSNKHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSNKHierarchyIntegrator::putToDatabase()");
        timers_need_init = false;
    }
    return;
}// INSNKHierarchyIntegrator

INSNKHierarchyIntegrator::~INSNKHierarchyIntegrator()
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

    if (d_petsc_nullsp != static_cast<MatNullSpace>(NULL))
    {
        PetscErrorCode ierr;
        ierr = MatNullSpaceDestroy(d_petsc_nullsp); IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(d_petsc_nullsp_vec); IBTK_CHKERRQ(ierr);
    }
    return;
}// ~INSNKHierarchyIntegrator

const std::string&
INSNKHierarchyIntegrator::getName() const
{
    return d_object_name;
}// getName

void
INSNKHierarchyIntegrator::registerVelocityInitialConditions(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> U_init)
{
    d_U_init = U_init;
    return;
}// registerVelocityInitialConditions

void
INSNKHierarchyIntegrator::registerPressureInitialConditions(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> P_init)
{
    d_P_init = P_init;
    return;
}// registerPressureInitialConditions

void
INSNKHierarchyIntegrator::registerBodyForceSpecification(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> F_set)
{
    d_F_set = F_set;
    return;
}// registerBodyForceSpecification

void
INSNKHierarchyIntegrator::registerVisItDataWriter(
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
INSNKHierarchyIntegrator::getHierarchyMathOps() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hier_math_ops.isNull());
#endif
    return d_hier_math_ops;
}// getHierarchyMathOps

void
INSNKHierarchyIntegrator::setHierarchyMathOps(
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
INSNKHierarchyIntegrator::isManagingHierarchyMathOps() const
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
///  allow the INSNKHierarchyIntegrator to be used as a hierarchy integrator.
///

void
INSNKHierarchyIntegrator::initializeHierarchyIntegrator(
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
    d_U_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::U",NDIM);
    d_P_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::P");
#if (NDIM == 2)
    d_Omega_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Omega");
#endif
#if ( NDIM == 3)
    d_Omega_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Omega",NDIM);
    d_Omega_Norm_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Omega_Norm");
#endif
    d_Div_U_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Div_U");

    d_V_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::V",NDIM);
    d_convective_flux_var = new SAMRAI::pdat::FaceVariable<NDIM,double>(d_object_name+"::convective_flux");

    // Create the default communication algorithms.
    d_fill_after_regrid = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    d_calgs["SYNCH_CURRENT_STATE_DATA"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    d_calgs["SYNCH_NEW_STATE_DATA"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();

    // Register state variables that are maintained by the
    // INSNKHierarchyIntegrator.
    const SAMRAI::hier::IntVector<NDIM> cell_ghosts = CELLG;
    const SAMRAI::hier::IntVector<NDIM> face_ghosts = FACEG;
    const SAMRAI::hier::IntVector<NDIM> weno_cell_ghosts = WENO_CELLG;
    const SAMRAI::hier::IntVector<NDIM> weno_face_ghosts = WENO_FACEG;
    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    registerVariable(d_U_current_idx, d_U_new_idx, d_U_scratch_idx,
                     d_U_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_P_current_idx, d_P_new_idx, d_P_scratch_idx,
                     d_P_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE");

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
    // INSNKHierarchyIntegrator.
    d_V_idx = var_db->registerVariableAndContext(
        d_V_var, getScratchContext(), weno_cell_ghosts);

    d_convective_flux_idx[0] = var_db->registerVariableAndContext(
        d_convective_flux_var, getScratchContext(), weno_face_ghosts);
    for (int d = 1; d < NDIM; ++d)
    {
        d_convective_flux_idx[d] = var_db->registerClonedPatchDataIndex(
            d_convective_flux_var, d_convective_flux_idx[0]);
    }

    // Register variables for plotting.
    if (!d_visit_writer.isNull())
    {
        if (d_output_U)
        {
            d_visit_writer->registerPlotQuantity(
                d_U_var->getName(), "VECTOR", d_U_current_idx, 0, d_U_scale);

            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(
                    d_U_var->getName()+stream.str(), "SCALAR", d_U_current_idx, d, d_U_scale);
            }
        }

        if (d_output_P)
        {
            d_visit_writer->registerPlotQuantity(
                d_P_var->getName(), "SCALAR", d_P_current_idx, 0, d_P_scale);
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

    // Create refinement communications algorithms.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator;
    refine_operator = grid_geom->lookupRefineOperator(
        d_V_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["SYNCH_V_DATA"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    d_ralgs["SYNCH_V_DATA"]->registerRefine(d_V_idx, // destination
                                            d_V_idx, // source
                                            d_V_idx, // temporary work space
                                            refine_operator);
    d_rstrategies["SYNCH_V_DATA"] = new IBTK::CartExtrapPhysBdryOp(
        d_V_idx, BDRY_EXTRAP_TYPE);

    // Create coarsening communications algorithms.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > coarsen_operator;
    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_V_var, "CONSERVATIVE_COARSEN");
    d_calgs["SYNCH_V_DATA"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    d_calgs["SYNCH_V_DATA"]->registerCoarsen(d_V_idx, // destination
                                             d_V_idx, // source
                                             coarsen_operator);

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
INSNKHierarchyIntegrator::initializeHierarchy()
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
INSNKHierarchyIntegrator::advanceHierarchy(
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

    if (initial_time)
    {
        for (int cycle = 0; cycle < 4; ++cycle)
        {
            d_hier_cc_data_ops->copyData(d_P_current_idx, d_P_new_idx);
            resetHierDataToPreadvanceState();
            integrateHierarchy(current_time, new_time);
        }
    }

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
INSNKHierarchyIntegrator::getStableTimestep(
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
INSNKHierarchyIntegrator::atRegridPoint() const
{
    const int level_number = 0;
    return ((d_integrator_step > 0)
            && d_gridding_alg->levelCanBeRefined(level_number)
            && (d_regrid_interval == 0
                ? false
                : (d_integrator_step % d_regrid_interval == 0)));
}// atRegridPoint

double
INSNKHierarchyIntegrator::getIntegratorTime() const
{
    return d_integrator_time;
}// getIntegratorTime

double
INSNKHierarchyIntegrator::getStartTime() const
{
    return d_start_time;
}// getStartTime

double
INSNKHierarchyIntegrator::getEndTime() const
{
    return d_end_time;
}// getEndTime

int
INSNKHierarchyIntegrator::getIntegratorStep() const
{
    return d_integrator_step;
}// getIntegratorStep

int
INSNKHierarchyIntegrator::getMaxIntegratorSteps() const
{
    return d_max_integrator_steps;
}// getMaxIntegratorSteps

bool
INSNKHierarchyIntegrator::stepsRemaining() const
{
    return (d_integrator_step < d_max_integrator_steps);
}// stepsRemaining

const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >
INSNKHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> >
INSNKHierarchyIntegrator::getGriddingAlgorithm() const
{
    return d_gridding_alg;
}// getGriddingAlgorithm

SAMRAI::tbox::Pointer<HierarchyProjector>
INSNKHierarchyIntegrator::getHierarchyProjector() const
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
///  allow the INSNKHierarchyIntegrator to provide data management for a time
///  integrator which making use of this class.
///

void
INSNKHierarchyIntegrator::regridHierarchy()
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
INSNKHierarchyIntegrator::integrateHierarchy(
    const double current_time,
    const double new_time)
{
    t_integrate_hierarchy->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double dt = new_time - current_time;

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

#if 0
    // 2nd-order SSP Runge-Kutta.
    integrateHierarchy_1st_order(
        current_time, current_time+dt, d_U_current_idx, d_P_current_idx, d_U_new_idx, d_P_new_idx);

    integrateHierarchy_1st_order(
        current_time+dt, current_time+2.0*dt, d_U_new_idx, d_P_new_idx, d_U_new_idx, d_P_new_idx);

    d_hier_cc_data_ops->linearSum(d_U_new_idx, 0.5, d_U_current_idx, 0.5, d_U_new_idx);
    d_hier_cc_data_ops->linearSum(d_P_new_idx, 0.5, d_P_current_idx, 0.5, d_P_new_idx);
#endif
    // 3rd-order SSP Runge-Kutta.
    integrateHierarchy_1st_order(
        current_time, current_time+dt, d_U_current_idx, d_P_current_idx, d_U_new_idx, d_P_new_idx);

    integrateHierarchy_1st_order(
        current_time+dt, current_time+2.0*dt, d_U_new_idx, d_P_new_idx, d_U_new_idx, d_P_new_idx);

    d_hier_cc_data_ops->linearSum(d_U_new_idx, 0.75, d_U_current_idx, 0.25, d_U_new_idx);
    d_hier_cc_data_ops->linearSum(d_P_new_idx, 0.75, d_P_current_idx, 0.25, d_P_new_idx);

    integrateHierarchy_1st_order(
        current_time+0.5*dt, current_time+1.5*dt, d_U_new_idx, d_P_new_idx, d_U_new_idx, d_P_new_idx);

    d_hier_cc_data_ops->linearSum(d_U_new_idx, 1.0/3.0, d_U_current_idx, 2.0/3.0, d_U_new_idx);
    d_hier_cc_data_ops->linearSum(d_P_new_idx, 1.0/3.0, d_P_current_idx, 2.0/3.0, d_P_new_idx);

    // Setup U_scratch to allow for ghost cell filling.
    d_hier_cc_data_ops->copyData(d_U_scratch_idx, d_U_new_idx);
    d_U_bdry_fill_op->fillData(new_time);

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
        1.0, d_U_scratch_idx, d_U_var, d_no_fill_op, new_time,
        0.0, -1, SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >(NULL),
        0, 0);

    t_integrate_hierarchy->stop();
    return getStableTimestep(getNewContext());
}// integrateHierarchy

void
INSNKHierarchyIntegrator::synchronizeHierarchy()
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
INSNKHierarchyIntegrator::synchronizeNewLevels(
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
INSNKHierarchyIntegrator::resetTimeDependentHierData(
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
INSNKHierarchyIntegrator::resetHierDataToPreadvanceState()
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
INSNKHierarchyIntegrator::initializeLevelData(
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

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_current_data =
                    patch->getPatchData(d_U_current_idx);
                U_current_data->fillAll(0.0);

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

            // Fill in U boundary data from coarser levels.
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom =
                d_hierarchy->getGridGeometry();
            SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator =
                grid_geom->lookupRefineOperator(
                    d_U_var, "CONSERVATIVE_LINEAR_REFINE");
            SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > ralg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
            ralg->registerRefine(d_U_scratch_idx, // destination
                                 d_U_current_idx, // source
                                 d_U_scratch_idx, // temporary work space
                                 refine_operator);
            IBTK::CartExtrapPhysBdryOp bc_op(d_U_scratch_idx, BDRY_EXTRAP_TYPE);
            ralg->createSchedule(level, level_number-1, hierarchy, &bc_op)->fillData(init_data_time);

            // Initialize quantities derived from the initial value of U.
            IBTK::PatchMathOps patch_math_ops;
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_scratch_data =
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
    }

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
INSNKHierarchyIntegrator::resetHierarchyConfiguration(
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

    d_hier_fc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_fc_data_ops->resetLevels(0, finest_hier_level);

    // Reset the Hierarchy math operations for the new configuration.
    if (d_is_managing_hier_math_ops)
    {
        d_hier_math_ops->setPatchHierarchy(hierarchy);
        d_hier_math_ops->resetLevels(0, finest_hier_level);
    }

    // Get the cell weight variable and patch data descriptor index.
    d_wgt_var = d_hier_math_ops->getCellWeightVariable();
    d_wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

    // Get the volume of the physical domain.
    d_volume = d_hier_math_ops->getVolumeOfPhysicalDomain();

    // Reset the solution, rhs, and nullspace vectors.
    d_sol_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::sol_vec", d_hierarchy, 0, finest_hier_level);
    d_sol_vec->addComponent(d_U_var,d_U_scratch_idx,d_wgt_idx,d_hier_cc_data_ops);
    d_sol_vec->addComponent(d_P_var,d_P_scratch_idx,d_wgt_idx,d_hier_cc_data_ops);

    d_rhs_vec = d_sol_vec->cloneVector(d_object_name+"::rhs_vec");
    d_nul_vec = d_sol_vec->cloneVector(d_object_name+"::nul_vec");

    // Setup the nullspace object.
    PetscErrorCode ierr;
    if (d_petsc_nullsp != static_cast<MatNullSpace>(NULL))
    {
        ierr = MatNullSpaceDestroy(d_petsc_nullsp); IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(d_petsc_nullsp_vec); IBTK_CHKERRQ(ierr);
    }
    d_petsc_nullsp_vec = IBTK::PETScSAMRAIVectorReal<double>::createPETScVector(d_nul_vec, PETSC_COMM_SELF);
    Vec vecs[] = {d_petsc_nullsp_vec};
    static const PetscTruth has_cnst = PETSC_FALSE;
    ierr = MatNullSpaceCreate(PETSC_COMM_SELF, has_cnst, 1, vecs, &d_petsc_nullsp); IBTK_CHKERRQ(ierr);

    // Setup the boundary condition objects.
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* U_bc_coef = NULL;  // XXXX
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef = NULL;  // XXXX

    // Setup the patch boundary filling objects.
    typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent U_scratch_component(d_U_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, U_bc_coef);
    InterpolationTransactionComponent P_scratch_component(d_P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, P_bc_coef);

    d_U_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_U_bdry_fill_op->initializeOperatorState(U_scratch_component, d_hierarchy);

    d_P_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_P_bdry_fill_op->initializeOperatorState(P_scratch_component, d_hierarchy);

    std::vector<InterpolationTransactionComponent> U_P_transaction_comps(2);
    U_P_transaction_comps[0] = U_scratch_component;
    U_P_transaction_comps[1] = P_scratch_component;

    d_U_P_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_U_P_bdry_fill_op->initializeOperatorState(U_P_transaction_comps, d_hierarchy);

    // Setup the linear operator.
    d_linear_op = new INSNKOperator(
        d_object_name+"::linear_op",
        d_rho, d_mu, d_lambda,
        d_old_dt,
        d_integrator_time,
        d_hier_math_ops, d_U_P_bdry_fill_op);

    // Setup the linear solver.
    // XXXX extra options needed
    d_linear_solver = new IBTK::PETScKrylovLinearSolver(
        d_object_name+"::linear_solver", "ins_");
    d_linear_solver->setInitialGuessNonzero(true);
    d_linear_solver->setOperator(d_linear_op);
    d_linear_solver->initializeSolverState(*d_sol_vec,*d_rhs_vec);

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
INSNKHierarchyIntegrator::applyGradientDetector(
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

SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >
INSNKHierarchyIntegrator::getVelocityVar()
{
    return d_U_var;
}// getVelocityVar

SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >
INSNKHierarchyIntegrator::getPressureVar()
{
    return d_P_var;
}// getPressureVar

SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >
INSNKHierarchyIntegrator::getForceVar()
{
    TBOX_ASSERT(false);  // XXXX
    return NULL;
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
INSNKHierarchyIntegrator::getCurrentContext() const
{
    return d_current_context;
}// getCurrentContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
INSNKHierarchyIntegrator::getNewContext() const
{
    return d_new_context;
}// getNewContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
INSNKHierarchyIntegrator::getScratchContext() const
{
    return d_scratch_context;
}// getScratchContext

///
///  The following routines:
///
///      putToDatabase()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::tbox::Serializable abstract base class.
///

void
INSNKHierarchyIntegrator::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("INS_NK_HIERARCHY_INTEGRATOR_VERSION",
                   INS_NK_HIERARCHY_INTEGRATOR_VERSION);
    db->putDouble("d_start_time", d_start_time);
    db->putDouble("d_end_time", d_end_time);
    db->putDouble("d_grow_dt", d_grow_dt);
    db->putInteger("d_max_integrator_steps", d_max_integrator_steps);
    db->putInteger("d_num_cycles", d_num_cycles);
    db->putInteger("d_num_init_cycles", d_num_init_cycles);
    db->putInteger("d_regrid_interval", d_regrid_interval);
    db->putBool("d_using_default_tag_buffer", d_using_default_tag_buffer);
    db->putIntegerArray("d_tag_buffer", d_tag_buffer);
    db->putBool("d_using_vorticity_tagging", d_using_vorticity_tagging);
    db->putDoubleArray("d_Omega_rel_thresh", d_Omega_rel_thresh);
    db->putDoubleArray("d_Omega_abs_thresh", d_Omega_abs_thresh);
    db->putDouble("d_Omega_max", d_Omega_max);
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
INSNKHierarchyIntegrator::printClassData(
    std::ostream& os) const
{
    os << "\nINSNKHierarchyIntegrator::printClassData..." << std::endl;
    os << "this = " << const_cast<INSNKHierarchyIntegrator*>(this) << std::endl;
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
    os << "d_using_vorticity_tagging = " << d_using_vorticity_tagging << "\n"
       << "d_Omega_rel_thresh = [ ";
    std::copy(d_Omega_rel_thresh.getPointer(), d_Omega_rel_thresh.getPointer()+d_Omega_rel_thresh.size(), std::ostream_iterator<double>(os, " , "));
    os << " ]\n"
       << "d_Omega_abs_thresh = [ ";
    std::copy(d_Omega_abs_thresh.getPointer(), d_Omega_abs_thresh.getPointer()+d_Omega_abs_thresh.size(), std::ostream_iterator<double>(os, " , "));
    os << " ]\n"
       << "d_Omega_max = " << d_Omega_max << std::endl;
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
       << "d_hier_fc_data_ops = " << d_hier_fc_data_ops.getPointer() << "=n"
       << "d_hier_math_ops = " << d_hier_math_ops.getPointer() << "\n"
       << "d_is_managing_hier_math_ops = " << d_is_managing_hier_math_ops << std::endl;
    os << "d_wgt_var = " << d_wgt_var.getPointer() << "\n"
       << "d_volume = " << d_volume << std::endl;
    os << "Skipping variables, patch data descriptors, communications algorithms, etc." << std::endl;
    return;
}// printClassData

/////////////////////////////// PROTECTED ////////////////////////////////////

void
INSNKHierarchyIntegrator::registerVariable(
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
INSNKHierarchyIntegrator::registerVariable(
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
INSNKHierarchyIntegrator::integrateHierarchy_1st_order(
    const double current_time,
    const double new_time,
    const int U_current_idx,
    const int P_current_idx,
    const int U_new_idx,
    const int P_new_idx)
{
    t_integrate_hierarchy_1st_order->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(U_current_idx != d_U_scratch_idx);
    TBOX_ASSERT(P_current_idx != d_P_scratch_idx);
    TBOX_ASSERT(U_new_idx != d_U_scratch_idx);
    TBOX_ASSERT(P_new_idx != d_P_scratch_idx);
#endif

    const double dt = new_time - current_time;

    // Setup solver vectors.
    d_hier_cc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(0), U_current_idx);
    d_hier_cc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(1), P_current_idx);

    d_rhs_vec->allocateVectorData(current_time);
    computeConvectiveDerivative(
        current_time,
        d_rhs_vec->getComponentDescriptorIndex(0), d_rhs_vec->getComponentVariable(0),
        U_current_idx);
    d_hier_cc_data_ops->linearSum(d_rhs_vec->getComponentDescriptorIndex(0), d_rho/dt, U_current_idx, -d_rho, d_rhs_vec->getComponentDescriptorIndex(0));
    d_hier_cc_data_ops->setToScalar(d_rhs_vec->getComponentDescriptorIndex(1), 0.0);

    d_nul_vec->allocateVectorData(current_time);
    d_hier_cc_data_ops->setToScalar(d_nul_vec->getComponentDescriptorIndex(0), 0.0);
    d_hier_cc_data_ops->setToScalar(d_nul_vec->getComponentDescriptorIndex(1), 1.0);

    // Setup the linear operator.
    d_linear_op = new INSNKOperator(
        d_object_name+"::linear_op",
        d_rho, d_mu, d_lambda,
        dt,
        new_time,
        d_hier_math_ops, d_U_P_bdry_fill_op);
    d_linear_solver->setOperator(d_linear_op);

    // Setup the preconditioner.
    SAMRAI::solv::PoissonSpecifications poisson_spec(d_object_name+"::poisson_spec");
    poisson_spec.setCZero();
    poisson_spec.setDConstant(-1.0);

    SAMRAI::tbox::Pointer<IBTK::CCPoissonFACOperator> poisson_fac_op = new IBTK::CCPoissonFACOperator(
        d_object_name+"::FAC Op", d_poisson_fac_op_db);
    poisson_fac_op->setPoissonSpecifications(poisson_spec);

    SAMRAI::tbox::Pointer<SAMRAI::solv::FACPreconditioner<NDIM> > poisson_fac_pc = new SAMRAI::solv::FACPreconditioner<NDIM>(
        d_object_name+"::FAC Preconditioner", *poisson_fac_op, d_poisson_fac_pc_db);
    poisson_fac_op->setPreconditioner(poisson_fac_pc);

    static const bool homogeneous_bc = true;
    SAMRAI::tbox::Pointer<IBTK::CCLaplaceOperator> poisson_op = new IBTK::CCLaplaceOperator(
        d_object_name+"::Laplace Operator",
        poisson_spec, NULL, homogeneous_bc);

    SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver> poisson_solver = new IBTK::PETScKrylovLinearSolver(
        d_object_name+"::PETSc Krylov Poisson solver", "poisson_");
    poisson_solver->setInitialGuessNonzero(true);
    poisson_solver->setOperator(poisson_op);
    poisson_solver->setPreconditioner(new IBTK::FACPreconditionerLSWrapper(poisson_fac_pc, d_poisson_fac_pc_db));

    SAMRAI::solv::PoissonSpecifications helmholtz_spec(d_object_name+"::helmholtz_spec");
    helmholtz_spec.setCConstant((d_rho/dt)+d_lambda);
    helmholtz_spec.setDConstant(-d_mu);

    SAMRAI::tbox::Pointer<IBTK::CCLaplaceOperator> helmholtz_op = new IBTK::CCLaplaceOperator(
        d_object_name+"::Laplace Operator",
        helmholtz_spec, std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(NDIM,static_cast<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(NULL)), homogeneous_bc);

    SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver> helmholtz_solver = new IBTK::PETScKrylovLinearSolver(
        d_object_name+"::PETSc Krylov Helmholtz solver", "helmholtz_");
    helmholtz_solver->setInitialGuessNonzero(true);
    helmholtz_solver->setOperator(helmholtz_op);

    SAMRAI::tbox::Pointer<IBTK::LinearSolver> pc_op = new INSNKPreconditioner(
        d_object_name+"::preconditioner_op",
        d_rho, d_mu, d_lambda,
        dt,
        new_time,
        poisson_solver, helmholtz_solver,
        d_hier_math_ops, d_U_bdry_fill_op, d_P_bdry_fill_op);
    d_linear_solver->setPreconditioner(pc_op);

    // Solve system.
    d_linear_solver->solveSystem(*d_sol_vec,*d_rhs_vec);
    d_linear_solver->setPreconditioner(NULL);  // XXXX

    // Pull out solution components.
    d_hier_cc_data_ops->copyData(U_new_idx, d_sol_vec->getComponentDescriptorIndex(0));
    d_hier_cc_data_ops->copyData(P_new_idx, d_sol_vec->getComponentDescriptorIndex(1));

    // Deallocate scratch data.
    d_rhs_vec->deallocateVectorData();
    d_nul_vec->deallocateVectorData();

    t_integrate_hierarchy_1st_order->stop();
    return;
}// integrateHierarchy_1st_order

void
INSNKHierarchyIntegrator::computeConvectiveDerivative(
    const double current_time,
    const int N_current_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > N_var,
    const int U_current_idx)
{
#if (NDIM == 2)
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_V_idx, current_time);
        for (int d = 0; d < NDIM; ++d)
        {
            level->allocatePatchData(d_convective_flux_idx[d], current_time);
        }
    }

    d_hier_cc_data_ops->copyData(d_V_idx, U_current_idx);
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["SYNCH_V_DATA"][ln]->coarsenData();
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_rscheds["SYNCH_V_DATA"][ln]->fillData(current_time);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
            const SAMRAI::hier::Index<NDIM>& patch_upper = patch_box.upper();

            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U =
                patch->getPatchData(d_V_idx);

            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Flux =
                new SAMRAI::pdat::CellData<NDIM,double>(patch_box, 1, SAMRAI::hier::IntVector<NDIM>(WENO_CELLG));
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > flux_fwrd =
                new SAMRAI::pdat::FaceData<NDIM,double>(patch_box, 1, SAMRAI::hier::IntVector<NDIM>(WENO_FACEG));
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > flux_revr =
                new SAMRAI::pdat::FaceData<NDIM,double>(patch_box, 1, SAMRAI::hier::IntVector<NDIM>(WENO_FACEG));
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > flux_plus =
                new SAMRAI::pdat::FaceData<NDIM,double>(patch_box, 1, SAMRAI::hier::IntVector<NDIM>(WENO_FACEG));
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > flux_minus =
                new SAMRAI::pdat::FaceData<NDIM,double>(patch_box, 1, SAMRAI::hier::IntVector<NDIM>(WENO_FACEG));

            for (int d = 0; d < NDIM; ++d)
            {
                SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > flux =
                    patch->getPatchData(d_convective_flux_idx[d]);
                WENO_CONVECTIVE_FLUXES_F77(
                    flux      ->getPointer(0), flux      ->getPointer(1), WENO_FACEG,
                    flux_fwrd ->getPointer(0), flux_fwrd ->getPointer(1), WENO_FACEG,
                    flux_revr ->getPointer(0), flux_revr ->getPointer(1), WENO_FACEG,
                    flux_plus ->getPointer(0), flux_plus ->getPointer(1), WENO_FACEG,
                    flux_minus->getPointer(0), flux_minus->getPointer(1), WENO_FACEG,
                    Flux->getPointer(0), WENO_CELLG,
                    U   ->getPointer(d), WENO_CELLG,
                    U   ->getPointer(0), WENO_CELLG,
                    patch_lower(0), patch_upper(0),
                    patch_lower(1), patch_upper(1));
            }
        }
    }

    static const bool cf_bdry_synch = true;

    for (int d = 0; d < NDIM; ++d)
    {
        d_hier_math_ops->div(
            N_current_idx, N_var,
            1.0, d_convective_flux_idx[d], d_convective_flux_var, d_no_fill_op, current_time, cf_bdry_synch,
            0.0, -1, SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >(NULL),
            d, d);
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_V_idx);
        for (int d = 0; d < NDIM; ++d)
        {
            level->deallocatePatchData(d_convective_flux_idx[d]);
        }
    }
#endif
    return;
}// computeConvectiveDerivative

double
INSNKHierarchyIntegrator::getLevelDt(
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
INSNKHierarchyIntegrator::getPatchDt(
    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx) const
{
    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
        patch->getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const SAMRAI::hier::Index<NDIM>& ilower = patch->getBox().lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch->getBox().upper();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_data =
        patch->getPatchData(d_U_var, ctx);
    const SAMRAI::hier::IntVector<NDIM>& U_ghost_cells = U_data->getGhostCellWidth();

    double stable_dt = std::numeric_limits<double>::max();
#if (NDIM == 2)
    NAVIER_STOKES_STABLEDT_F77(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),
        U_ghost_cells(0),U_ghost_cells(1),
        U_data->getPointer(0),U_data->getPointer(1),
        stable_dt);
#endif
#if (NDIM == 3)
    NAVIER_STOKES_STABLEDT_F77(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
        U_ghost_cells(0),U_ghost_cells(1),U_ghost_cells(2),
        U_data->getPointer(0),U_data->getPointer(1),U_data->getPointer(2),
        stable_dt);
#endif
    return stable_dt;
}// getPatchDt

void
INSNKHierarchyIntegrator::getFromInput(
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

    d_output_U     = db->getBoolWithDefault("output_U"    , d_output_U     );
    d_output_P     = db->getBoolWithDefault("output_P"    , d_output_P     );
    d_output_F     = db->getBoolWithDefault("output_F"    , d_output_F     );
    d_output_Omega = db->getBoolWithDefault("output_Omega", d_output_Omega );
    d_output_Div_U = db->getBoolWithDefault("output_Div_U", d_output_Div_U );

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

    if (db->keyExists("PoissonFACOp"))
    {
        d_poisson_fac_op_db = db->getDatabase("PoissonFACOp");
    }
    if (db->keyExists("PoissonFACPreconditioner"))
    {
        d_poisson_fac_pc_db = db->getDatabase("PoissonFACPreconditioner");
    }

    if (db->keyExists("HelmholtzFACOp"))
    {
        d_helmholtz_fac_op_db = db->getDatabase("HelmholtzFACOp");
    }
    if (db->keyExists("HelmholtzFACPreconditioner"))
    {
        d_helmholtz_fac_pc_db = db->getDatabase("HelmholtzFACPreconditioner");
    }
    return;
}// getFromInput

void
INSNKHierarchyIntegrator::getFromRestart()
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

    int ver = db->getInteger("INS_NK_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != INS_NK_HIERARCHY_INTEGRATOR_VERSION)
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
    d_using_vorticity_tagging = db->getBool("d_using_vorticity_tagging");
    d_Omega_rel_thresh = db->getDoubleArray("d_Omega_rel_thresh");
    d_Omega_abs_thresh = db->getDoubleArray("d_Omega_abs_thresh");
    d_Omega_max = db->getDouble("d_Omega_max");
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
template class SAMRAI::tbox::Pointer<IBAMR::INSNKHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
