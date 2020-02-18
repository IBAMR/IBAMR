// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBAMR_INSVCStaggeredConservativeHierarchyIntegrator
#define included_IBAMR_INSVCStaggeredConservativeHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/INSVCStaggeredConservativeMassMomentumIntegrator.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"

#include <string>
#include <vector>

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
 * of the time step \f$\rho^n\f$, an interpolated face density \f$\rho^{n+\frac{1}{2}}\f$ is produced and used in the
 * momentum convection to obtain \f$N(\rho u)\f$ and mass advection to obtain \f$\rho^{n+1}\f$. Hence, a consistent
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
    INSVCStaggeredConservativeHierarchyIntegrator(std::string object_name,
                                                  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                                  bool register_for_restart = true);

    /*!
     * The destructor for class INSVCStaggeredConservativeHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~INSVCStaggeredConservativeHierarchyIntegrator() = default;

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
    initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

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
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1) override;

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchy(double current_time, double new_time, int cycle_num = 0) override;

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override;

    /*!
     * Explicitly remove nullspace components from a solution vector.
     */
    void removeNullSpace(const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& sol_vec);

    /*
     * \brief Supply boundary conditions for the side-centered density field, which is maintained by this integrator
     *
     */
    void registerMassDensityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* rho_bc_coef) override;

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
    int getNumberOfCycles() const override;

    /*!
     * Get the convective operator being used by the integrator class.
     *
     * \note The class employs INSVCStaggeredConservativeMassMomentumIntegrator
     * to compute the conservative convective derivative. Therefore,
     * ConvectiveOperator is a NULL object.
     */
    SAMRAI::tbox::Pointer<ConvectiveOperator> getConvectiveOperator() override;

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
                                        bool allocate_data) override;

    /*!
     * Reset cached hierarchy dependent data.
     */
    void
    resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                           int coarsest_level,
                                           int finest_level) override;

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to the magnitude of the fluid vorticity.
     */
    void applyGradientDetectorSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                          int level_number,
                                          double error_data_time,
                                          int tag_index,
                                          bool initial_time,
                                          bool uses_richardson_extrapolation_too) override;

    /*!
     * Prepare variables for plotting.
     */
    void setupPlotDataSpecialized() override;

    /*!
     * Project the velocity field following a regridding operation.
     */
    void regridProjection() override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSVCStaggeredConservativeHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSVCStaggeredConservativeHierarchyIntegrator(const INSVCStaggeredConservativeHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSVCStaggeredConservativeHierarchyIntegrator&
    operator=(const INSVCStaggeredConservativeHierarchyIntegrator& that) = delete;

    /*!
     * Update the operators and solvers to account for changes due to time-dependent coefficients
     */
    void updateOperatorsAndSolvers(double current_time, double new_time, int cycle_num);

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

    /*
     * Patch data descriptor indices for all "state" variables managed by the
     * integrator.
     *
     * State variables have three contexts: current, scratch, and new.
     */
    int d_rho_sc_current_idx, d_rho_sc_scratch_idx, d_rho_sc_new_idx;

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
