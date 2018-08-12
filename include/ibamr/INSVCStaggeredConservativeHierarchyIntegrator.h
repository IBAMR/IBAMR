// Filename: INSVCStaggeredConservativeHierarchyIntegrator.h
// Created on 15 May 2018 by Nishant Nangia and Amneet Bhalla
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

#ifndef included_IBAMR_INSVCStaggeredConservativeHierarchyIntegrator
#define included_IBAMR_INSVCStaggeredConservativeHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"

namespace IBAMR
{
class INSVCStaggeredConservativeMassMomentumIntegrator;
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
 * \brief Class INSVCStaggeredConservativeHierarchyIntegrator provides a staggered-grid solver
 * for the incompressible Navier-Stokes equations on an AMR grid hierarchy, with variable
 * coefficients.
 *
 * This class integrates the conservative form of the momentum equation
 * \f$(\frac{\partial \rho u}{\partial t} + N(\rho u)) = -\nabla p + \nabla \cdot \mu (\nabla u) + (\nabla u)^T )\f$
 * where \f$ N(u) = \nabla \cdot (\rho u u) \f$.
 *
 * In other words, the class treats the left-hand side of the momentum equation in conservative form.
 * This class is specialized to use INSVCStaggeredConservativeMassMomentumIntegrator,
 * which produces an update for the newest density by solving the mass transport equation
 * \f$ \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho u) = 0 \f$.
 * Therefore, the density variable must be registered and maintained by this class,
 * and not by any other integrator (such AdvDiffHierarchyIntegrator).
 * It is also assumed that this density variable is side-centered. In other words, given a density at the beginning
 * of the time step \f$rho^n\f$, an interpolated face density \f$rho^{n+\frac{1}{2}}\f$ is produced and used in the
 * momentum convection to obtain \f$N(\rho u)\f$ and mass advection to obtain \f$rho^{n+1}\f$. Hence, a consistent
 * momentum and mass transport is carried out, which leads to stable solutions for high-density ratio flows.
 *
 */
class INSVCStaggeredConservativeHierarchyIntegrator : public INSVCStaggeredHierarchyIntegrator
{
public:
    /*!
     * The constructor for class INSVCStaggeredConservativeHierarchyIntegrator sets some
     * default values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    INSVCStaggeredConservativeHierarchyIntegrator(const std::string& object_name,
                                                  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                                  bool register_for_restart = true);

    /*!
     * The destructor for class INSVCStaggeredConservativeHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~INSVCStaggeredConservativeHierarchyIntegrator();

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
     * \brief Supply boundary conditions for the side-centered density field, which is maintained by this integrator
     *
     */
    void registerMassDensityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* rho_bc_coef);

    /*
     * \brief Supply boundary conditions for the side-centered density field, which is maintained by this integrator
     */
    void
    registerMassDensityBoundaryConditions(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& rho_sc_bc_coefs);

    /*!
     * \brief Supply a source term for the mass update equation.
     *
     * \note Current implementation is used only to check order of accuracy via a manufactured solution.
     */
    void registerMassDensitySourceTerm(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> S_fcn);

    /*!
     * Returns the number of cycles to perform for the present time step.
     */
    int getNumberOfCycles() const;

    /*!
     * Get the convective operator being used by the integrator class.
     *
     * \note The class employs INSVCStaggeredConservativeMassMomentumIntegrator
     * to compute the conservative convective derivative. Therefore,
     * ConvectiveOperator is a NULL object.
     */
    SAMRAI::tbox::Pointer<ConvectiveOperator> getConvectiveOperator();

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
    INSVCStaggeredConservativeHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSVCStaggeredConservativeHierarchyIntegrator(const INSVCStaggeredConservativeHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSVCStaggeredConservativeHierarchyIntegrator& operator=(const INSVCStaggeredConservativeHierarchyIntegrator& that);

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
     * Side-centered density variable required for conservative discretization
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_rho_sc_var;

    /*!
     * Side-centered velocity variable maintained from previous time step
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_old_var;

    /*
     * Patch data descriptor indices for all "state" variables managed by the
     * integrator.
     *
     * State variables have three contexts: current, scratch, and new.
     */
    int d_rho_sc_current_idx, d_rho_sc_scratch_idx, d_rho_sc_new_idx;
    int d_U_old_current_idx, d_U_old_new_idx, d_U_old_scratch_idx;

    /*
     * Boundary condition object for the side-centered density variable maintained
     * by this integrator.
     */
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_rho_sc_bc_coefs;

    /*
     * Variables for plotting cell-centered density
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_rho_interp_cc_var;
    int d_rho_interp_cc_idx;

    /*
     * Source term function for the mass density update
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_S_fcn;

    /*
     * Conservative density and momentum integrator.
     */
    SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredConservativeMassMomentumIntegrator> d_rho_p_integrator;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_INSVCStaggeredConservativeHierarchyIntegrator
