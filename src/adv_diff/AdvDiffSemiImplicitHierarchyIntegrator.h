// Filename: AdvDiffSemiImplicitHierarchyIntegrator.h
// Created on 22 May 2012 by Boyce Griffith
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

#ifndef included_AdvDiffSemiImplicitHierarchyIntegrator
#define included_AdvDiffSemiImplicitHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSC INCLUDES
#include <petscsys.h>

// IBAMR INCLUDES
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/ConvectiveOperator.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class AdvDiffSemiImplicitHierarchyIntegrator manages the spatial
 * discretization and time integration of scalar- and vector-valued quantities
 * whose dynamics are governed by the advection-diffusion equation.
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
 * Various options are available for the spatial and temporal discretizations.
 *
 * \see HierarchyIntegrator
 * \see SAMRAI::mesh::StandardTagAndInitStrategy
 * \see SAMRAI::algs::TimeRefinementIntegrator
 * \see SAMRAI::algs::TimeRefinementLevelStrategy
 */
class AdvDiffSemiImplicitHierarchyIntegrator
    : public AdvDiffHierarchyIntegrator
{
public:
    /*!
     * The constructor for class AdvDiffSemiImplicitHierarchyIntegrator sets
     * some default values, reads in configuration information from input and
     * restart databases, and registers the integrator object with the restart
     * manager when requested.
     */
    AdvDiffSemiImplicitHierarchyIntegrator(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        bool register_for_restart=true);

    /*!
     * The destructor for class AdvDiffSemiImplicitHierarchyIntegrator
     * unregisters the integrator object with the restart manager when the
     * object is so registered.
     */
    ~AdvDiffSemiImplicitHierarchyIntegrator();

    /*!
     * Return the type of convective time integration scheme being employed by the
     * advection-diffusion solver.
     *
     * The choice of time integration scheme is set via the input database
     * provided to the class constructor.
     */
    TimeSteppingType
    getConvectiveTimeSteppingType() const;

    /*!
     * Return the type of convective time integration scheme being employed by the
     * advection-diffusion solver during the initial time step.
     *
     * The choice of time integration scheme is set via the input database
     * provided to the class constructor.
     *
     * \note This is used \em only when the basic convective time stepping
     * scheme uses a multi-step method such as Adams-Bashforth.
     */
    TimeSteppingType
    getInitialConvectiveTimeSteppingType() const;

    /*!
     * \brief Set the default convective operator type to be used by the solver.
     */
    void
    setDefaultConvectiveOperatorType(
        ConvectiveOperatorType op_type);

    /*!
     * \brief Get the default convective operator type used by the solver.
     */
    ConvectiveOperatorType
    getDefaultConvectiveOperatorType() const;

    /*!
     * \brief Set the default convective differencing form to be used by the
     * solver.
     */
    void
    setDefaultConvectiveDifferencingType(
        ConvectiveDifferencingType difference_form);

    /*!
     * \brief Get the default convective differencing form used by the solver.
     */
    ConvectiveDifferencingType
    getDefaultConvectiveDifferencingType() const;

    /*!
     * Register an operator to compute the convective derivative term u*grad Q
     * for a particular transported quantity Q.
     *
     * The boolean flag needs_reinit_when_dt_changes indicates whether the
     * operator needs to be explicitly reinitialized when the time step size
     * changes.
     */
    void
    setConvectiveOperator(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        SAMRAI::tbox::Pointer<ConvectiveOperator> convective_op,
        bool needs_reinit_when_dt_changes);

    /*!
     * Get the convective operator being used by this solver class for a
     * particular transported quantity Q.
     *
     * If the convective operator has not already been constructed, then this
     * function will initialize a default convective operator.
     */
    virtual SAMRAI::tbox::Pointer<ConvectiveOperator>
    getConvectiveOperator(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var);

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
     * Returns the number of cycles to perform for the present time step.
     */
    int
    getNumberOfCycles() const;

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void
    preprocessIntegrateHierarchy(
        double current_time,
        double new_time,
        int num_cycles=1);

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void
    integrateHierarchy(
        double current_time,
        double new_time,
        int cycle_num=0);

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void
    postprocessIntegrateHierarchy(
        double current_time,
        double new_time,
        bool skip_synchronize_new_state_data,
        int num_cycles=1);

protected:
    /*!
     * Return the maximum stable time step size.
     */
    double
    getTimeStepSizeSpecialized();

    /*!
     * Reset cached hierarchy dependent data.
     */
    void
    resetHierarchyConfigurationSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int coarsest_level,
        int finest_level);

    /*!
     * Enum indicating the time integration employed for the explicit
     * discretization of the convective terms.
     */
    TimeSteppingType d_convective_time_stepping_type;

    /*!
     * Enum indicating the time integration employed for the explicit
     * discretization of the convective terms during the \em initial time step.
     */
    TimeSteppingType d_init_convective_time_stepping_type;

    /*
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyFaceDataOpsReal<NDIM,double> > d_hier_fc_data_ops;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffSemiImplicitHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffSemiImplicitHierarchyIntegrator(
        const AdvDiffSemiImplicitHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffSemiImplicitHierarchyIntegrator&
    operator=(
        const AdvDiffSemiImplicitHierarchyIntegrator& that);

    /*!
     * Read input values from a given database.
     */
    void
    getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
        bool is_from_restart);

    /*!
     * Value indicating the number of solver cycles to be used for the present
     * time step.
     */
    int d_num_cycles_step;

    /*!
     * Advective CFL condition.
     */
    double d_cfl_max;

    /*!
     * Transported quantities.
     */
    std::set<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > d_N_var, d_N_old_var;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > d_Q_N_map, d_Q_N_old_map;

    ConvectiveOperatorType d_default_convective_op_type;
    ConvectiveDifferencingType d_default_convective_difference_form;
    std::string d_default_convective_bdry_extrap_type;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,SAMRAI::tbox::Pointer<ConvectiveOperator> > d_convective_op;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >,bool> d_convective_op_needs_init, d_convective_op_needs_reinit_when_dt_changes;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_AdvDiffSemiImplicitHierarchyIntegrator
