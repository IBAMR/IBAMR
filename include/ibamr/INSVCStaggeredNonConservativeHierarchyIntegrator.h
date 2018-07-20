// Filename: INSVCStaggeredNonConservativeHierarchyIntegrator.h
// Created on 17 May 2018 by Nishant Nangia and Amneet Bhalla
//
// Copyright (c) 2002-2018, Nishant Nangia and Amneet Bhalla
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

#ifndef included_IBAMR_INSVCStaggeredNonConservativeHierarchyIntegrator
#define included_IBAMR_INSVCStaggeredNonConservativeHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"

namespace IBAMR
{
} // namespace IBAMR
namespace IBTK
{
} // namespace IBTK
namespace SAMRAI
{
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSVCStaggeredNonConservativeHierarchyIntegrator provides a staggered-grid solver
 * for the incompressible Navier-Stokes equations on an AMR grid hierarchy, with variable
 * coefficients.
 *
 * This class always uses a non-conservative discretization of the form of the momentum equation
 * \f$\rho(\frac{\partial u}{\partial t} + N(u)) = -\nabla p + \nabla \cdot \mu (\nabla u) + (\nabla u)^T )\f$
 * where \f$ N(u) = u \cdot \nabla u \f$ for convective_difference_form = ADVECTIVE and
 * \f$ N(u) = \nabla \cdot (u u) \f$ for convective_difference_form = CONSERVATIVE.
 *
 * In other words, this class will NEVER treat the left-hand side of the momentum equation in conservative form
 * i.e. \frac{\partial \rho u}{\partial t} + \nabla \cdot (\rho u u)
 *
 * In this class, both the density and viscosity are required to be cell centered quantities, which are then
 * interpolated onto the required degrees of freedom (side-centers for rho and node(edge)-centers for mu in 2D(3D))
 *
 * Note that this class is suitable for low density ratio flows. At high density ratios, the non-conservative form
 * can lead to instabilities, and INSVCStaggeredConservativeHierarchyIntegrator should be used instead.
 *
 */
class INSVCStaggeredNonConservativeHierarchyIntegrator : public INSVCStaggeredHierarchyIntegrator
{
public:
    /*!
     * The constructor for class INSVCStaggeredNonConservativeHierarchyIntegrator sets some
     * default values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    INSVCStaggeredNonConservativeHierarchyIntegrator(const std::string& object_name,
                                                     SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                                     bool register_for_restart = true);

    /*!
     * The destructor for class INSVCStaggeredNonConservativeHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~INSVCStaggeredNonConservativeHierarchyIntegrator();

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                       SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Initialize the AMR patch hierarchy and data defined on the hierarchy at
     * the start of a computation.  If the computation is begun from a restart
     * file, the patch hierarchy and patch data are read from the hierarchy
     * database.  Otherwise, the patch hierarchy and patch data are initialized
     * by the gridding algorithm associated with the integrator object.
     *
     * The implementation of this function assumes that the hierarchy exists
     * upon entry to the function, but that it contains no patch levels.  On
     * return from this function, the state of the integrator object will be
     * such that it is possible to step through time via the advanceHierarchy()
     * function.
     */
    void initializePatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1);

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchy(double current_time, double new_time, int cycle_num = 0);

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1);

    /*!
     * Regrid the patch hierarchy.
     */
    void regridHierarchy();

    /*!
     * Explicitly remove nullspace components from a solution vector.
     */
    void removeNullSpace(const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& sol_vec);

    /*
     * \brief Supply boundary conditions for the density field, if maintained by the fluid integrator.
     *
     * \note The boundary conditions set here will be overwritten if density if being advected.
     */
    void registerMassDensityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* rho_bc_coef);

    /*
     * \brief Set the transported density variable if it is being maintained by the advection-diffusion integrator.
     *
     * \note The variable set here MUST be registered and maintained by the advection-diffusion integrator.
     */
    void setTransportedMassDensityVariable(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > rho_adv_diff_var);

protected:
    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     */
    void initializeLevelDataSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                        int level_number,
                                        double init_data_time,
                                        bool can_be_refined,
                                        bool initial_time,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
                                        bool allocate_data);

    /*!
     * Reset cached hierarchy dependent data.
     */
    void
    resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                           int coarsest_level,
                                           int finest_level);

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to the magnitude of the fluid vorticity.
     */
    void applyGradientDetectorSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                          int level_number,
                                          double error_data_time,
                                          int tag_index,
                                          bool initial_time,
                                          bool uses_richardson_extrapolation_too);

    /*!
     * Prepare variables for plotting.
     */
    void setupPlotDataSpecialized();

    /*!
     * Project the velocity field following a regridding operation.
     */
    void regridProjection();

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSVCStaggeredNonConservativeHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSVCStaggeredNonConservativeHierarchyIntegrator(const INSVCStaggeredNonConservativeHierarchyIntegrator& from) =
        delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSVCStaggeredNonConservativeHierarchyIntegrator&
    operator=(const INSVCStaggeredNonConservativeHierarchyIntegrator& that) = delete;

    /*!
     * Determine the convective time stepping type for the current time step and
     * cycle number.
     */
    TimeSteppingType getConvectiveTimeSteppingType(int cycle_num);

    /*!
     * Update the operators and solvers to account for changes due to time-dependent coefficients
     */
    void updateOperatorsAndSolvers(double current_time, double new_time);

    /*!
     * Setup solution and RHS vectors using state data maintained by the
     * integrator.
     */
    void setupSolverVectors(const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& sol_vec,
                            const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& rhs_vec,
                            double current_time,
                            double new_time,
                            int cycle_num);

    /*!
     * Copy the solution data into the state data maintained by
     * the integrator.
     */
    void resetSolverVectors(const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& sol_vec,
                            const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& rhs_vec,
                            double current_time,
                            double new_time,
                            int cycle_num);

    /*!
     * Interpolated density variable required for non-conservative discretization
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_rho_interp_var;

    /*
     * Patch data descriptor indices for all "state" variables managed by the
     * integrator.
     *
     * State variables have three contexts: current, scratch, and new.
     */
    int d_rho_current_idx, d_rho_new_idx, d_rho_scratch_idx;

    /*
     * Patch data descriptor indices for all "scratch" variables managed by the
     * integrator.
     *
     * Scratch variables have only one context: scratch.
     */
    int d_rho_interp_idx;

    /*
     * Boundary condition objects for density, which is provided by an appropriate advection-diffusion
     * integrator, or set by the fluid integrator.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_rho_bc_coef;

    /*
     * Variable to keep track of a transported density variable maintained by an advection-diffusion integrator
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_rho_adv_diff_var;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_INSVCStaggeredNonConservativeHierarchyIntegrator
