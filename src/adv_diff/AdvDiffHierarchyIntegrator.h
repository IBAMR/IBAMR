// Filename: AdvDiffHierarchyIntegrator.h
// Created on 21 May 2012 by Boyce Griffith
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
#include <petscsys.h>

// IBAMR INCLUDES
#include <ibamr/ibamr_enums.h>
#include <ibamr/ibamr_utilities.h>

// IBTK INCLUDES
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CCPoissonPointRelaxationFACOperator.h>
#include <ibtk/HierarchyIntegrator.h>
#include <ibtk/KrylovLinearSolver.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class AdvDiffHierarchyIntegrator provides an abstract interface for a
 * time integrator for advection-diffusion or advection-reaction-diffusion
 * equations on an AMR grid hierarchy, along with basic data management for
 * variables defined on that hierarchy.
 */
class AdvDiffHierarchyIntegrator
    : public IBTK::HierarchyIntegrator
{
public:
    /*!
     * The destructor for class AdvDiffHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~AdvDiffHierarchyIntegrator();

    /*!
     * Return the type of diffusion time integration scheme being employed by
     * the advection-diffusion solver.
     *
     * At the present time, supported time integration schemes include:
     *
     *    - BACKWARD_EULER
     *    - TRAPEZOIDAL_RULE (same as Crank-Nicolson)
     *
     * The choice of time integration scheme is set via the input database
     * provided to the class constructor.
     */
    TimeSteppingType
    getDiffusionTimeSteppingType() const;

    /*!
     * Register a face-centered advection velocity to be used to advect
     * cell-centered quantities by the hierarchy integrator.
     *
     * Data management for the registered advection velocity will be handled by
     * the hierarchy integrator.
     *
     * \note By default, each registered advection velocity is assumed to be
     * divergence free.
     */
    void
    registerAdvectionVelocity(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > u_var);

    /*!
     * Indicate whether a particular advection velocity is discretely divergence
     * free.
     */
    void
    setAdvectionVelocityIsDivergenceFree(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > u_var,
        bool is_div_free);

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
     * Data management for the registered source term will be handled by the
     * hierarchy integrator.
     */
    void
    registerSourceTerm(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var);

    /*!
     * Supply an IBTK::CartGridFunction object to specify the value of a
     * particular source term.
     */
    void
    setSourceTermFunction(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var,
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> F_fcn);

    /*!
     * Register a cell-centered quantity to be advected and diffused by the
     * hierarchy integrator.
     *
     * Data management for the registered quantity will be handled by the
     * hierarchy integrator.
     */
    void
    registerTransportedQuantity(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var);

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
     * Set the convective differencing form for a quantity that has been
     * registered with the hierarchy integrator.
     */
    void
    setConvectiveDifferencingType(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        ConvectiveDifferencingType difference_form);

    /*!
     * Set the scalar diffusion coefficient corresponding to a quantity that has
     * been registered with the hierarchy integrator.
     */
    void
    setDiffusionCoefficient(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        double kappa);

    /*!
     * Set the scalar linear damping coefficient corresponding to a quantity
     * that has been registered with the hierarchy integrator.
     */
    void
    setDampingCoefficient(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        double lambda);

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

protected:
    /*!
     * The constructor for class AdvDiffHierarchyIntegrator sets some default
     * values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    AdvDiffHierarchyIntegrator(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        bool register_for_restart);

    /*!
     * Reset cached hierarchy dependent data.
     */
    void
    resetHierarchyConfigurationSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int coarsest_level,
        int finest_level);

    /*!
     * Write out specialized object state to the given database.
     */
    void
    putToDatabaseSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * Standard variable registration.
     */
    void
    registerVariables();

    /*
     * Boolean value that indicates whether the integrator has been initialized.
     */
    bool d_integrator_is_initialized;

    /*!
     * Enum indicating the time integration employed for the implicit
     * discretization of the diffusion terms.
     */
    TimeSteppingType d_diffusion_time_stepping_type;

    /*!
     * Advection velocity data.
     */
    std::set<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > > d_u_var;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >,bool> d_u_is_div_free;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >,SAMRAI::tbox::Pointer<IBTK::CartGridFunction> > d_u_fcn;

    /*!
     * Source term data.
     */
    std::set<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > d_F_var;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,SAMRAI::tbox::Pointer<IBTK::CartGridFunction> > d_F_fcn;

    /*!
     * Transported quantities.
     */
    std::set<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > d_Q_var, d_Psi_var;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > > d_Q_u_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > d_Q_F_map, d_Q_Psi_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,ConvectiveDifferencingType> d_Q_difference_form;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,double> d_Q_diffusion_coef, d_Q_damping_coef;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,SAMRAI::tbox::Pointer<IBTK::CartGridFunction> > d_Q_init;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> > d_Q_bc_coef;

    /*
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > d_hier_cc_data_ops;
    std::vector<SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> > d_hier_bdry_fill_ops;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_no_fill_op;

    /*
     * Linear solvers and associated data.
     */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > > d_sol_vecs, d_rhs_vecs;

    int d_max_iterations;
    double d_abs_residual_tol, d_rel_residual_tol;
    bool d_using_FAC;

    std::vector<SAMRAI::tbox::Pointer<IBTK::CCLaplaceOperator> >                   d_helmholtz_ops;
    std::vector<SAMRAI::solv::PoissonSpecifications>                               d_helmholtz_specs;
    std::vector<SAMRAI::tbox::Pointer<IBTK::KrylovLinearSolver> >                  d_helmholtz_solvers;
    std::vector<SAMRAI::tbox::Pointer<IBTK::CCPoissonPointRelaxationFACOperator> > d_helmholtz_fac_ops;
    std::vector<SAMRAI::tbox::Pointer<IBTK::FACPreconditioner> >                   d_helmholtz_fac_pcs;
    std::vector<bool> d_helmholtz_solvers_need_init;
    int d_coarsest_reset_ln, d_finest_reset_ln;

    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_fac_op_db, d_fac_pc_db;

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
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/AdvDiffHierarchyIntegrator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_AdvDiffHierarchyIntegrator
