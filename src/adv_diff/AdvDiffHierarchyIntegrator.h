// Filename: AdvDiffHierarchyIntegrator.h
// Created on 16 Mar 2004 by Boyce Griffith
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

#ifndef included_AdvDiffHierarchyIntegrator
#define included_AdvDiffHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSC INCLUDES
#include <petsc.h>

// IBAMR INCLUDES
#include <ibamr/AdvDiffHypPatchOps.h>
#include <ibamr/HierarchyIntegrator.h>

// IBTK INCLUDES
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CCPoissonFACOperator.h>
#include <ibtk/KrylovLinearSolver.h>

// SAMRAI INCLUDES
#include <HyperbolicLevelIntegrator.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class AdvDiffHierarchyIntegrator manages the spatial discretization
 * and time integration of scalar- and vector-valued quantities whose dynamics
 * are governed by the advection-diffusion equation.
 *
 * Each quantity \f$ Q \f$ managed by the integrator may have a unique diffusion
 * coefficient \f$ \kappa \f$ and damping coefficient \f$ \lambda \f$, and may
 * optionally have a forcing term \f$ F \f$.  Additionally, a different
 * advection velocity may be used with each quantity registered with the
 * integrator.
 *
 * This hierarchy integrator advances all levels of the patch hierarchy
 * synchronously in time.  In particular, subcycling in time is \em not
 * performed.
 *
 * Either Crank-Nicolson (i.e., the trapezoidal rule) or backward Euler is used
 * for the linearly implicit treatment of the diffusive terms.  The advective
 * terms are discretized by the GodunovAdvector object supplied to the class
 * constructor.
 *
 * \see AdvDiffHypPatchOps
 * \see HierarchyIntegrator
 * \see GodunovAdvector
 * \see SAMRAI::algs::HyperbolicLevelIntegrator
 * \see SAMRAI::mesh::StandardTagAndInitStrategy
 * \see SAMRAI::algs::TimeRefinementIntegrator
 * \see SAMRAI::algs::TimeRefinementLevelStrategy
 */
class AdvDiffHierarchyIntegrator
    : public HierarchyIntegrator
{
public:
    /*!
     * The constructor for class AdvDiffHierarchyIntegrator sets some default
     * values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    AdvDiffHierarchyIntegrator(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geometry,
        SAMRAI::tbox::Pointer<GodunovAdvector> explicit_predictor,
        bool register_for_restart=true);

    /*!
     * The destructor for class AdvDiffHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~AdvDiffHierarchyIntegrator();

    /*!
     * Return the type of viscous time integration scheme being employed by the
     * advection-diffusion solver.
     *
     * At the present time, valid time integration schemes include:
     *
     *    - CRANK_NICOLSON
     *    - BACKWARD_EULER
     *
     * The choice of time integration scheme is set via the input database
     * provided to the class constructor.
     */
    const ViscousTimesteppingType&
    getViscousTimesteppingType() const;

    /*!
     * Register a face-centered advection velocity to be used to advect
     * cell-centered quantities by the hierarchy integrator.
     *
     * \note By default, data management for the registered advection velocity
     * will be handled by the hierarchy integrator.
     *
     * \note By default, each registered advection velocity is assumed to be
     * divergence free.
     */
    void
    registerAdvectionVelocity(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > u_var,
        const bool manage_data=true);

    /*!
     * Indicate whether a particular advection velocity is discretely divergence
     * free.
     */
    void
    setAdvectionVelocityIsDivergenceFree(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > u_var,
        const bool is_div_free);

    /*!
     * Supply an IBTK::CartGridFunction object to specify the value of a
     * particular advection velocity.
     */
    void
    setAdvectionVelocityFunction(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > u_var,
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> u_fcn);

    /*!
     * Register a cell-centered source term.
     *
     * \note By default, data management for the registered source term will be
     * handled by the hierarchy integrator.
     */
    void
    registerSourceTerm(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var,
        const bool manage_data=true);

    /*!
     * Supply an IBTK::CartGridFunction object to specify the value of a
     * particular source term.
     */
    void
    setSourceTermFunction(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var,
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> F_fcn);

    /*!
     * Register a face-centered "incompressibility fix" term.
     *
     * This term will be \em subtracted from the predicted face-centered values
     * generated by the explicit predictor.  Such terms are useful in
     * cell-centered incompressible flow solvers to account approximately for
     * the incompressibility of the velocity field.
     *
     * \note By default, data management for the registered advection velocity
     * will be handled by the hierarchy integrator.
     */
    void
    registerIncompressibilityFixTerm(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_Phi_var,
        const bool manage_data=true);

    /*!
     * Supply an IBTK::CartGridFunction object to specify the value of a
     * particular source term.
     */
    void
    setIncompressibilityFixTermFunction(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_Phi_var,
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> grad_Phi_fcn);

    /*!
     * Register a cell-centered quantity to be advected and diffused by the
     * hierarchy integrator.
     *
     * \note By default, data management for the registered quantity will be
     * handled by the hierarchy integrator.
     */
    void
    registerTransportedQuantity(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        const bool manage_data=true);

    /*!
     * Set the face-centered advection velocity to be used with a particular
     * cell-centered quantity.
     *
     * \note The specified advection velocity must have been already registered
     * with the hierarchy integrator.
     */
    void
    setAdvectionVelocity(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > u_var);

    /*!
     * Set the cell-centered source term to be used with a particular
     * cell-centered quantity.
     *
     * \note The specified source term must have been already registered with
     * the hierarchy integrator.
     */
    void
    setSourceTerm(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var);

    /*!
     * Set the face-centered "incompressibility fix" term to be used with a
     * particular cell-centered quantity.
     *
     * \note The specified "incompressibility fix" term must have been already
     * registered with the hierarchy integrator.
     */
    void
    setIncompressibilityFixTerm(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_Phi_var);

    /*!
     * Set the convective differencing form for a quantity that has been
     * registered with the hierarchy integrator.
     */
    void
    setConvectiveDifferencingType(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        const ConvectiveDifferencingType& difference_form);

    /*!
     * Set the scalar diffusion coefficient corresponding to a quantity that has
     * been registered with the hierarchy integrator.
     */
    void
    setDiffusionCoefficient(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        const double& kappa);

    /*!
     * Set the scalar linear damping coefficient corresponding to a quantity
     * that has been registered with the hierarchy integrator.
     */
    void
    setDampingCoefficient(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        const double& lambda);

    /*!
     * Set a grid function to provide initial conditions for a quantity that has
     * been registered with the hierarchy integrator.
     */
    void
    setInitialConditions(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> Q_init);

    /*!
     * Set an object to provide boundary conditions for a scalar-valued quantity
     * that has been registered with the hierarchy integrator.
     */
    void
    setPhysicalBcCoefs(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* Q_bc_coef);

    /*!
     * Set objects to provide boundary conditions for a vector-valued quantity
     * that has been registered with the hierarchy integrator.
     */
    void
    setPhysicalBcCoefs(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> Q_bc_coef);

    /*!
     * Return a pointer to the level integrator object used to integrate the
     * advective terms.
     */
    SAMRAI::tbox::Pointer<SAMRAI::algs::HyperbolicLevelIntegrator<NDIM> >
    getHyperbolicLevelIntegrator() const;

    /*!
     * Return a pointer to the patch strategy object used to specify the
     * numerical routines used to integrate the advective terms.
     */
    SAMRAI::tbox::Pointer<AdvDiffHypPatchOps>
    getHyperbolicPatchStrategy() const;

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void
    initializeHierarchyIntegrator(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void
    integrateHierarchy(
        const double current_time,
        const double new_time,
        const int cycle_num=0);

    /*!
     * Return the maximum stable time step size.
     *
     * A default implementation is provided that returns
     * min(dt_max,dt_growth_factor*dt_current).  The growth condition prevents
     * excessive changes in the time step size as the computation progresses.
     */
    double
    getStableTimestep();

    /*!
     * Reset the current data to equal the new data, update the time level of
     * the current data, and deallocate the scratch and new data.
     */
    void
    resetTimeDependentHierarchyData(
        const double new_time);

    /*!
     * Reset the hierarchy integrator to the state at the beginning of the
     * current time step.
     */
    void
    resetIntegratorToPreadvanceState();

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::initializeLevelData
     */
    void
    initializeLevelData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level=SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> >(NULL),
        const bool allocate_data=true);

    /*!
     * Reset cached hierarchy dependent data.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::resetHierarchyConfiguration
     */
    void
    resetHierarchyConfiguration(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int coarsest_level,
        const int finest_level);

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to gradient criteria specified by the
     * GodunovAdvector object.
     */
    void
    applyGradientDetector(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double error_data_time,
        const int tag_index,
        const bool initial_time,
        const bool uses_richardson_extrapolation_too);

    /*!
     * Write out specialized object state to the given database.
     */
    void
    putToDatabaseSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    /*!
     * Enum indicating the time integration employed for the implicit
     * discretization of the viscous terms.
     */
    ViscousTimesteppingType d_viscous_timestepping_type;

    /*!
     * Advection velocity data.
     */
    std::set<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > > d_u_var;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >,bool> d_manage_u_data, d_u_is_div_free;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >,SAMRAI::tbox::Pointer<IBTK::CartGridFunction> > d_u_fcn;

    /*!
     * Source term data.
     */
    std::set<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > d_F_var;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,bool> d_manage_F_data;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,SAMRAI::tbox::Pointer<IBTK::CartGridFunction> > d_F_fcn;

    /*!
     * Incompressibility fix data.
     */
    std::set<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > > d_grad_Phi_var;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >,bool> d_manage_grad_Phi_data;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >,SAMRAI::tbox::Pointer<IBTK::CartGridFunction> > d_grad_Phi_fcn;

    /*!
     * Transported quantities.
     */
    std::set<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > d_Q_var, d_Psi_var;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,bool> d_manage_Q_data;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > > d_Q_u_map, d_Q_grad_Phi_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > d_Q_F_map, d_Q_Psi_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,ConvectiveDifferencingType> d_Q_difference_form;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,double> d_Q_diffusion_coef, d_Q_damping_coef;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,SAMRAI::tbox::Pointer<IBTK::CartGridFunction> > d_Q_init;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> > d_Q_bc_coef;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffHierarchyIntegrator(
        const AdvDiffHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffHierarchyIntegrator&
    operator=(
        const AdvDiffHierarchyIntegrator& that);

    /*!
     * Read input values from a given database.
     */
    void
    getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
        bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.  The database from which the restart data are read is determined
     * by the object_name specified in the class constructor.
     */
    void
    getFromRestart();

    /*
     * The SAMRAI::algs::HyperbolicLevelIntegrator supplies generic operations
     * use to handle the explicit integration of advection terms.
     */
    SAMRAI::tbox::Pointer<SAMRAI::algs::HyperbolicLevelIntegrator<NDIM> > d_hyp_level_integrator;

    /*
     * The advection patch strategy supplies the advection-specific operations
     * needed to treat data on patches in the AMR grid hierarchy.
     */
    SAMRAI::tbox::Pointer<AdvDiffHypPatchOps> d_hyp_patch_ops;

    /*
     * Boolean value that indicates whether the integrator has been initialized.
     */
    bool d_integrator_is_initialized;

    /*
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;
    std::vector<SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> > d_hier_bdry_fill_ops;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_no_fill_op;

    /*
     * Variable context used for temporary storage.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_temp_context;

    /*
     * Linear solvers and associated data.
     */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > > d_sol_vecs, d_rhs_vecs;

    int d_max_iterations;
    double d_abs_residual_tol, d_rel_residual_tol;
    bool d_using_FAC;

    std::vector<SAMRAI::tbox::Pointer<IBTK::CCLaplaceOperator> >    d_helmholtz_ops;
    std::vector<SAMRAI::solv::PoissonSpecifications>                d_helmholtz_specs;
    std::vector<SAMRAI::tbox::Pointer<IBTK::KrylovLinearSolver> >   d_helmholtz_solvers;
    std::vector<SAMRAI::tbox::Pointer<IBTK::CCPoissonFACOperator> > d_helmholtz_fac_ops;
    std::vector<SAMRAI::tbox::Pointer<IBTK::FACPreconditioner> >    d_helmholtz_fac_pcs;
    std::vector<bool> d_helmholtz_solvers_need_init;
    int d_coarsest_reset_ln, d_finest_reset_ln;

    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_fac_op_db, d_fac_pc_db;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/AdvDiffHierarchyIntegrator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_AdvDiffHierarchyIntegrator
