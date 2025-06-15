// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_AcousticStreamingHierarchyIntegrator
#define included_IBAMR_AcousticStreamingHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/SOAcousticStreamingBcCoefs.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibamr/VCStaggeredStokesSpec.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/ProblemSpecification.h"
#include "ibtk/SCPoissonSolverManager.h"
#include "ibtk/SideDataSynchronization.h"
#include "ibtk/ibtk_enums.h"

#include "CellVariable.h"
#include "EdgeVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyEdgeDataOpsReal.h"
#include "HierarchyFaceDataOpsReal.h"
#include "HierarchyNodeDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "IntVector.h"
#include "NodeVariable.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

namespace IBAMR
{
class BrinkmanPenalizationStrategy;
class FOAcousticStreamingPETScLevelSolver;
} // namespace IBAMR
namespace IBTK
{
class PoissonSolver;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchLevel;
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class BasePatchHierarchy;
} // namespace hier
namespace mesh
{
template <int DIM>
class GriddingAlgorithm;
} // namespace mesh
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class AcousticStreamingHierarchyIntegrator implements time integrator for a
 * staggered-grid, surface acoustic wave (SAW)-driven (compressible) Navier-Stokes
 * solver on an AMR grid hierarchy with variable coefficients.
 *
 * The class maintains and solves both first (harmonic component) and second order (mean component)
 * equations for velocity and pressure. Note that the first order velocity, pressure, and density fields
 * are complex (i.e., have both real and imaginary components). The second-order variables are real.
 *
 * Density is related to pressure through the equation of state (EOS): \f$ p = c_0^2 \rho \f$.
 */

class AcousticStreamingHierarchyIntegrator : public IBTK::HierarchyIntegrator
{
public:
    /*!
     * The constructor for class AcousticStreamingHierarchyIntegrator sets some
     * default values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    AcousticStreamingHierarchyIntegrator(std::string object_name,
                                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                         bool register_for_restart = true);

    /*!
     * The destructor for class AcousticStreamingHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~AcousticStreamingHierarchyIntegrator();

    /*!
     * Supply a physical boundary conditions specification for the first order velocity
     * field.
     *
     * \note: For now we just impose velocity boundary conditions for the first order system
     */
    void registerFirstOrderPhysicalBoundaryConditions(
        const std::array<std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>, 2>& U_bc_coefs);

    /*!
     * Supply a physical boundary conditions specificaion for the second order fields
     */
    void registerSecondOrderPhysicalBoundaryConditions(
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs);

    /*!
     * Supply initial conditions for the first order velocity field.
     */
    void registerFirstOrderVelocityInitialConditions(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> U1_init);

    /*!
     * Supply initial conditions for the pressure.
     *
     * \note Initial conditions are not required for the pressure, but when
     * available, they can speed the convergence of the solver during the
     * initial time step.
     */
    void registerFirstOrderPressureInitialConditions(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> P1_init);

    /*!
     * Supply initial conditions for the second order velocity field.
     */
    void registerSecondOrderVelocityInitialConditions(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> U2_init);

    /*!
     * Supply initial conditions for the pressure.
     *
     * \note Initial conditions are not required for the pressure, but when
     * available, they can speed the convergence of the solver during the
     * initial time step.
     */
    void registerSecondOrderPressureInitialConditions(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> P2_init);

    /*!
     * Supply body force function for the first order momentum equation with the hierarchy integrator.
     */
    void registerFirstOrderBodyForceFunction(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> F1_fcn);

    /*!
     * Supply body force function for the second order momentum equation with the hierarchy integrator.
     */
    void registerSecondOrderBodyForceFunction(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> F2_fcn);

    /*!
     *  Supply a CartGridFunction that specifies \f$ <\nabla \rho \cdot \mathbf{u}_i + \frac{\omega}{c^2}p_r, \nabla
     * \rho \cdot \mathbf{u}_r - \frac{\omega}{c^2}p_i>\f$.
     *
     * @note If \param Q1_fcn is not specified, then the source terms are set to zero.
     */
    void registerFirstOrderVelocityDivergenceFunction(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> Q1_fcn);

    /*!
     * Supply a CartGridFunction that specifies \f$ -\nabla \rho \cdot \mathbf{u} \f$.
     *
     * @note If \param Q2_fcn is not specified, then \f$ \nabla \cdot \mathbf{u}  = 0 \f$ is imposed.
     */
    void registerSecondOrderVelocityDivergenceFunction(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> Q2_fcn);

    /*!
     * Get the subdomain solver for the second-order velocity subsystem.  Such solvers can be
     * useful in constructing block preconditioners.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> getSecondOrderVelocitySubdomainSolver();

    /*!
     * Get the subdomain solver for the second-order pressure subsystem.  Such solvers can be
     * useful in constructing block preconditioners.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> getSecondOrderPressureSubdomainSolver();

    /*!
     * Virtual method to initialize the variables, basic communications
     * algorithms, solvers, and other data structures used by a concrete time
     * integrator object.
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
     * Virtual method to initialize the AMR patch hierarchy and data defined on the hierarchy at
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
     * Virtual method to prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1) override;

    /*!
     * Virtual method to clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override;

    /*!
     * Explicitly remove nullspace components from a solution vector of the second-order system.
     */
    void
    removeSecondOrderNullSpace(const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& sol_vec);

    /*!
     * Register (zeroth order) density variable with the hierarchy integrator.
     */
    void registerMassDensityVariable(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var);

    /*!
     * Get the (zeroth order) density variable registered with the hierarchy integrator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getMassDensityVariable() const;

    /*!
     * Register shear viscosity variable with the hierarchy integrator.
     */
    void registerShearViscosityVariable(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var);

    /*!
     * Register bulk viscosity variable with the hierarchy integrator.
     */
    void registerBulkViscosityVariable(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > lambda_var);

    /*!
     * Get the shear viscosity variable registered with the hierarchy integrator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getShearViscosityVariable() const;

    /*!
     * Get the bulk viscosity variable registered with the hierarchy integrator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getBulkViscosityVariable() const;

    /*!
     * Set the interpolation type used for material properties rho
     */
    void setDensityVCInterpolationType(const IBTK::VCInterpType vc_interp_type);

    /*!
     * Set the interpolation type used for material properties mu (shear viscosity).
     */
    void setShearViscosityVCInterpolationType(const IBTK::VCInterpType vc_interp_type);

    /*!
     * \brief Function to reset fluid density or viscosity if they are
     * maintained by this integrator.
     */
    using ResetFluidPropertiesFcnPtr = void (*)(int property_idx,
                                                SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > property_var,
                                                SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                int cycle_num,
                                                double time,
                                                double current_time,
                                                double new_time,
                                                void* ctx);

    /*!
     * \brief Register function to reset fluid density.
     */
    void registerResetFluidDensityFcn(ResetFluidPropertiesFcnPtr callback, void* ctx);

    /*!
     * \brief Register function to reset fluid shear viscosity.
     */
    void registerResetFluidShearViscosityFcn(ResetFluidPropertiesFcnPtr callback, void* ctx);

    /*!
     * \brief Register function to reset fluid bulk viscosity.
     */
    void registerResetFluidBulkViscosityFcn(ResetFluidPropertiesFcnPtr callback, void* ctx);

    /*!
     * \brief Supply initial conditions for the (zeroth order) density field.
     */
    void registerMassDensityInitialConditions(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> rho_init_fcn);

    /*!
     * \brief Supply initial conditions for the shear viscosity field.
     */
    void registerShearViscosityInitialConditions(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> mu_init_fcn);

    /*!
     * \brief Supply initial conditions for the bulk viscosity field.
     */
    void registerBulkViscosityInitialConditions(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> lambda_init_fcn);

    /*
     * \brief Supply boundary conditions for the density field, if maintained by the fluid
     * integrator.
     */
    void
    registerMassDensityBoundaryConditions(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& rho_bc_coefs);

    /*
     * \brief Supply boundary conditions for the shear viscosity field, if maintained by the fluid integrator.
     */
    void registerShearViscosityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* mu_bc_coef);

    /*
     * \brief Supply boundary conditions for the bulk viscosity field, if maintained by the fluid integrator.
     */
    void registerBulkViscosityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* lambda_bc_coef);

    /*
     * \brief Get the first-order velocity variable
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getFirstOrderVelocityVariable() const
    {
        return d_U1_var;
    } // getFirstOrderVelocityVariable

    /*
     * \brief Get the first-order pressure variable
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getFirstOrderPressureVariable() const
    {
        return d_P1_var;
    } // getFirstOrderPressureVariable

    /*
     * \brief Get the second-order velocity variable
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getSecondOrderVelocityVariable() const
    {
        return d_U2_var;
    } // getSecondOrderVelocityVariable

    /*
     * \brief Get the second-order pressure variable
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getSecondOrderPressureVariable() const
    {
        return d_P2_var;
    } // getSecondOrderPressureVariable

    /*!
     * \brief Get the side-centered density patch data index, which will always be the newest one used in the linear
     * operator i.e. rho_sc in rho_sc*u^{n+1} term.
     *
     * \note These patch data will not be deallocated at the end of the time step, so they can be used for various
     * applications
     */
    inline int getLinearOperatorRhoPatchDataIndex() const
    {
        return d_rho_linear_op_idx;
    } // getLinearOperatorRhoPatchDataIndex

    /*!
     * \brief Get the cell-centered viscosity patch data index, which will always be the newest one used in the linear
     * operator.
     *
     * \note These patch data will not be deallocated at the end of the time step, so they can be used for various
     * applications
     */
    inline int getLinearOperatorMuPatchDataIndex() const
    {
        return d_mu_linear_op_idx;
    } // getLinearOperatorMuPatchDataIndex

    /*!
     * \brief Get the interpolated viscosity patch data index, which will always be the newest one used in the linear
     * operator.
     *
     * \note These patch data will not be deallocated at the end of the time step, so they can be used for various
     * applications
     */
    inline int getInterpolatedLinearOperatorMuPatchDataIndex() const
    {
        return d_mu_interp_linear_op_idx;
    } // getInterpolatedLinearOperatorMuPatchDataIndex

    /*!
     * \brief Get the shear viscosity boundary condition
     */
    inline SAMRAI::solv::RobinBcCoefStrategy<NDIM>* getShearViscosityBoundaryConditions() const
    {
        return d_mu_bc_coef;
    } // getShearViscosityBoundaryConditions

    /*!
     * \brief Get the bulk viscosity boundary condition
     */
    inline SAMRAI::solv::RobinBcCoefStrategy<NDIM>* getBulkViscosityBoundaryConditions() const
    {
        return d_lambda_bc_coef;
    } // getBulkViscosityBoundaryConditions

    /*!
     * \brief Get density boundary conditions
     */
    inline const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& getMassDensityBoundaryBoundaryConditions() const
    {
        return d_rho_bc_coefs;
    } // getMassDensityBoundaryBoundaryConditions

    /*!
     * \brief Register BrinkmanPenalizationStrategy objects to add the Brinkman penalization term
     * in the momentum equation for the first- and second-order systems.
     */
    void
    registerBrinkmanPenalizationStrategy(SAMRAI::tbox::Pointer<IBAMR::BrinkmanPenalizationStrategy> fo_brinkman_force,
                                         SAMRAI::tbox::Pointer<IBAMR::BrinkmanPenalizationStrategy> so_brinkman_force,
                                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > brinkman_var,
                                         SAMRAI::solv::RobinBcCoefStrategy<NDIM>* brinkman_bc);

    /*!
     * \brief Register level set function to perform contour integral to evaluate acoustic radiation force.
     */
    void registerContourVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > var,
                                 SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef,
                                 double val = 0.0);

    /*!
     * \brief Get the contour variables registered with this class.
     */
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > >& getContourVariables() const
    {
        return d_contour_vars;
    } // getContourVariables

    /*!
     * \brief Get the Brinkman penalization objects for the first-order velocity registered with this class.
     */
    const std::vector<SAMRAI::tbox::Pointer<IBAMR::BrinkmanPenalizationStrategy> >&
    getFOBrinkmanPenalizationStrategy() const
    {
        return d_fo_brinkman_force;
    } // getFOBrinkmanPenalizationStrategy

    /*!
     * \brief Get the Brinkman penalization objects for the second-order velocity registered with this class.
     */
    const std::vector<SAMRAI::tbox::Pointer<IBAMR::BrinkmanPenalizationStrategy> >&
    getSOBrinkmanPenalizationStrategy() const
    {
        return d_so_brinkman_force;
    } // getSOBrinkmanPenalizationStrategy

    /*!
     * \brief Get the variables associated with the Brinkman penalization objects registered with this class.
     */
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > >&
    getBrinkmanPenalizationVars() const
    {
        return d_brinkman_vars;
    } // getBrinkmanPenalizationVars

protected:
    /*!
     * Whether we need to perform a regrid projection when (re-)initializing composite hierarchy data.
     */
    bool d_do_regrid_projection = false;

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void
    integrateHierarchySpecialized(const double current_time, const double new_time, const int cycle_num = 0) override;

    /*!
     * Prepare the current hierarchy for regridding. Here we calculate the divergence.
     */
    void regridHierarchyBeginSpecialized() override;

    /*!
     * Update the current hierarchy data after regridding. Here we recalculate
     * the divergence and, if it has grown by a factor more than
     * d_regrid_max_div_growth_factor, we then project the velocity field onto
     * a divergence-free set of grid functions.
     */
    void regridHierarchyEndSpecialized() override;

    /*!
     * Perform data initialization after the entire hierarchy has been constructed.
     */
    void initializeCompositeHierarchyDataSpecialized(double init_data_time, bool initial_time) override;

    /*!
     * Virtual method to initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     */
    virtual void
    initializeLevelDataSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                   int level_number,
                                   double init_data_time,
                                   bool can_be_refined,
                                   bool initial_time,
                                   SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
                                   bool allocate_data) override;

    /*!
     * Virtual method to reset cached hierarchy dependent data.
     */
    virtual void
    resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                           int coarsest_level,
                                           int finest_level) override;

    /*!
     * Virtual method to set integer tags to "one" in cells where refinement of the given level
     * should occur according to the magnitude of the fluid vorticity.
     */
    virtual void
    applyGradientDetectorSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                     int level_number,
                                     double error_data_time,
                                     int tag_index,
                                     bool initial_time,
                                     bool uses_richardson_extrapolation_too) override;

    /*!
     * Virtual method to prepare variables for plotting.
     */
    virtual void setupPlotDataSpecialized() override;

    /*!
     * Write out specialized object state to the given database.
     */
    virtual void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * \name Parameters specific to the first order acoustic streaming system.
     */
    double d_acoustic_freq = std::numeric_limits<double>::signaling_NaN(),
           d_sound_speed = std::numeric_limits<double>::signaling_NaN();

    /*!
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM, double> > d_hier_sc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyNodeDataOpsReal<NDIM, double> > d_hier_nc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyEdgeDataOpsReal<NDIM, double> > d_hier_ec_data_ops;

    /*
     * Initial and boundary conditions and data synchronization operators.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_U1_init, d_P1_init, d_U2_init, d_P2_init;
    std::array<std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>, 2> d_U1_bc_coefs;
    std::array<IBAMR::SOAcousticStreamingBcCoefs, NDIM> d_default_so_bc_coefs;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_so_bc_coefs, d_U2_bc_coefs, d_U2_star_bc_coefs;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_P2_bc_coef;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_Phi_bc_coef;
    SAMRAI::tbox::Pointer<IBTK::SideDataSynchronization> d_side_synch1_op, d_side_synch2_op;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_rho_bdry_bc_fill_op, d_mu_bdry_bc_fill_op,
        d_lambda_bdry_bc_fill_op;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_U2_bdry_bc_fill_op, d_P2_bdry_bc_fill_op;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_no_fill_op;

    /*
     * Body force and mass source/sink terms for the first and second-order systems
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_F1_fcn, d_F2_fcn;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_Q1_fcn, d_Q2_fcn;

    /*
     * Solver for the first-order system
     */
    SAMRAI::tbox::Pointer<IBAMR::FOAcousticStreamingPETScLevelSolver> d_first_order_solver;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_first_order_solver_db;
    bool d_first_order_solver_needs_init = true;

    /*
     * Solvers for the second order system
     */
    int d_coarsest_reset_ln, d_finest_reset_ln;

    SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> d_bc_helper;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_sol1_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_rhs1_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_U2_scratch_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_U2_rhs_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_P2_scratch_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_P2_rhs_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_sol2_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_rhs2_vec;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > > d_null2_vecs;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > > d_U2_null_vecs;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > > d_U2_near_null_vecs;
    bool d_vectors_need_init, d_explicitly_remove_so_nullspace = false;

    std::string d_stokes_solver_type = IBAMR::StaggeredStokesSolverManager::UNDEFINED,
                d_stokes_precond_type = IBAMR::StaggeredStokesSolverManager::UNDEFINED;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_stokes_solver_db, d_stokes_precond_db;
    SAMRAI::tbox::Pointer<StaggeredStokesSolver> d_stokes_solver;
    bool d_stokes_solver_needs_init = true;

    std::string d_velocity_solver_type = IBTK::SCPoissonSolverManager::UNDEFINED,
                d_velocity_precond_type = IBTK::SCPoissonSolverManager::UNDEFINED;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_velocity_solver_db, d_velocity_precond_db;
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_velocity_solver;
    bool d_velocity_solver_needs_init;

    std::string d_pressure_solver_type = IBTK::CCPoissonSolverManager::UNDEFINED,
                d_pressure_precond_type = IBTK::CCPoissonSolverManager::UNDEFINED;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_pressure_solver_db, d_pressure_precond_db;
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_pressure_solver;
    bool d_pressure_solver_needs_init;

    /*!
     *
     * Boolean values indicates whether to output various quantities for
     * visualization.
     */
    bool d_output_U1 = true, d_output_U2 = true, d_output_P1 = true, d_output_P2 = true, d_output_Omega2 = true,
         d_output_rho = false, d_output_mu = false, d_output_lambda = false;

    /*!
     * Fluid solver variables.
     *
     * The 1st order variables are complex. We store them as variables with depth = 2.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U1_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U2_var;
    std::string d_U_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_U_refine_type = "CONSERVATIVE_LINEAR_REFINE";

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_P1_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_P2_var;
    std::string d_P_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_P_refine_type = "LINEAR_REFINE";

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_F1_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_F2_var;
    std::string d_F_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_F_refine_type = "CONSERVATIVE_LINEAR_REFINE";

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Q1_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Q2_var;
    std::string d_Q_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_Q_refine_type = "CONSTANT_REFINE";

    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_rho_var, d_mu_var, d_lambda_var;

    std::string d_mu_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_mu_refine_type = "CONSERVATIVE_LINEAR_REFINE";
    std::string d_mu_bdry_extrap_type = "CONSTANT";

    std::string d_rho_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_rho_refine_type = "CONSERVATIVE_LINEAR_REFINE";
    std::string d_rho_bdry_extrap_type = "CONSTANT";

    std::string d_lambda_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_lambda_refine_type = "CONSERVATIVE_LINEAR_REFINE";
    std::string d_lambda_bdry_extrap_type = "CONSTANT";

    /*!
     * Functions resetting rho, mu and lambda maintained by this integrator.
     */
    std::vector<ResetFluidPropertiesFcnPtr> d_reset_rho_fcns, d_reset_mu_fcns, d_reset_lambda_fcns;
    std::vector<void*> d_reset_rho_fcns_ctx, d_reset_mu_fcns_ctx, d_reset_lambda_fcns_ctx;

    /*!
     * Whether to couple first and second order systems.
     * By default they are coupled.
     */
    bool d_coupled_system = true;

    /*!
     * Whether to use Stokes drift formulation for the second order system.
     */
    bool d_use_stokes_drift_bc = true, d_use_stokes_drift_mass_src = true;

    /*!
     * The maximum CFL number.
     */
    double d_cfl_max = 1.0;

    /*!
     * This boolean value determines whether the pressure is normalized to have
     * zero mean (i.e., discrete integral) at the end of each timestep.
     */
    bool d_normalize_pressure = false;

    /*!
     * This boolean value determines whether the velocity is normalized to have
     * zero mean (i.e., discrete integral) at the end of each timestep.
     */
    bool d_normalize_velocity = false;

    /*!
     * This boolean value determines whether rigid body velocity is in the nullspace
     * of the velocity operator.
     */
    bool d_has_rigid_body_nullspace = false;

    /*
     * Boolean value that indicates whether the integrator has been initialized.
     */
    bool d_integrator_is_initialized = false;

    /*
     * Patch data descriptor indices for all "state" variables managed by the
     * integrator.
     *
     * State variables have three contexts: current, scratch, and new.
     */
    int d_U1_current_idx, d_U1_new_idx, d_U1_scratch_idx;
    int d_P1_current_idx, d_P1_new_idx, d_P1_scratch_idx;

    int d_F1_current_idx, d_F1_new_idx, d_F1_scratch_idx;
    int d_Q1_current_idx, d_Q1_new_idx, d_Q1_scratch_idx;

    int d_U2_current_idx, d_U2_new_idx, d_U2_scratch_idx;
    int d_P2_current_idx, d_P2_new_idx, d_P2_scratch_idx;

    int d_F2_current_idx, d_F2_new_idx, d_F2_scratch_idx;
    int d_Q2_current_idx, d_Q2_new_idx, d_Q2_scratch_idx;

    int d_mu_current_idx, d_mu_new_idx, d_mu_scratch_idx;
    int d_lambda_current_idx, d_lambda_new_idx, d_lambda_scratch_idx;
    int d_rho_current_idx, d_rho_new_idx, d_rho_scratch_idx;

    /*
     * Components of U1 and P1 stored as scratch indices.
     * These are used to compute coupling terms for the first- and
     * second-order systems, as well to fill boundary conditions for the
     * second-order velocity.
     */
    int d_U1_real_idx, d_U1_imag_idx;
    int d_p1_real_idx, d_p1_imag_idx;

    /*!
     * Cell tagging criteria based on the relative and absolute magnitudes of
     * the local vorticity.
     */
    bool d_using_vorticity_tagging = false;
    SAMRAI::tbox::Array<double> d_Omega2_rel_thresh, d_Omega2_abs_thresh;
    double d_Omega2_max = 0.0;

    /*
     * Patch data descriptor indices for all "plot" variables managed by the
     * integrator.
     *
     * Plot variables have one context: current.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_U1_plot_var, d_U2_plot_var, d_Omega2_var,
        d_rho_plot_var;
    int d_U1_plot_idx, d_U2_plot_idx, d_Omega2_idx, d_rho_plot_idx;

    /*
     * Patch data descriptor indices for all "scratch" variables managed by the
     * integrator.
     *
     * Scratch variables have only one context: scratch.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U2_regrid_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U2_src_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_indicator2_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Omega2_Norm_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_pressure_D_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_projection_D_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_velocity_C_var, d_velocity_L_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_velocity_D_cc_var;
#if (NDIM == 2)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > d_mu_interp_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > d_velocity_D_var;
#elif (NDIM == 3)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double> > d_mu_interp_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double> > d_velocity_D_var;
#endif

    int d_U2_regrid_idx, d_U2_src_idx, d_indicator2_idx;
    int d_Omega2_Norm_idx;
    int d_velocity_C_idx, d_velocity_L_idx, d_velocity_D_idx, d_velocity_D_cc_idx, d_pressure_D_idx, d_projection_D_idx;
    int d_mu_interp_idx;

    /*
     * Persistent patch data indices for the density and viscosity used in the linear operators
     */
    int d_mu_linear_op_idx;        // cell centered quantity
    int d_mu_interp_linear_op_idx; // node or edge centered
    int d_rho_linear_op_idx;       // side centered

    /*
     * Patch data index for Heaviside variable used in the projection preconditioner
     */
    int d_heaviside_cc_idx;

    /*
     * Variable to indicate the type of interpolation to be done for rho and mu.
     */
    IBTK::VCInterpType d_rho_vc_interp_type, d_mu_vc_interp_type;

    /*!
     * Brinkman force strategy objects registered with this integrator.
     */
    std::vector<SAMRAI::tbox::Pointer<IBAMR::BrinkmanPenalizationStrategy> > d_fo_brinkman_force, d_so_brinkman_force;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > > d_brinkman_vars;
    std::vector<int> d_brinkman_current_idx, d_brinkman_new_idx, d_brinkman_scratch_idx;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_brinkman_bcs;

    /*!
     * Contour variables registered with this integrator.
     */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > > d_contour_vars;
    std::vector<int> d_contour_idx;
    std::vector<double> d_contour_val;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_contour_bcs;

    /*!
     * Boolean to determine whether to write contour integrals to an output file.
     */
    bool d_write_contour_integrals = false;
    std::vector<std::unique_ptr<std::ofstream> > d_contour_integral_stream;

    /*
     * Variable to set how often the preconditioner is reinitialized.
     */
    int d_precond_reinit_interval = 1;

    /*
     * Objects to set initial condition for density and viscosity when they are maintained by the fluid integrator.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_rho_init_fcn, d_mu_init_fcn, d_lambda_init_fcn;

    /*
     * Boundary condition objects for viscosity and density.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>*d_mu_bc_coef = nullptr, *d_lambda_bc_coef = nullptr;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_rho_bc_coefs;

    /*
     * Problem specification object for VCStaggeredStokesOperator for the second order system.
     */
    VCStaggeredStokesOpSpec d_vc_stokes_op_spec;

    /*
     * Problem specification object for VCSCViscousDilatationalOp in the velocity solver.
     */
    IBTK::VCViscousDilatationalOpSpec d_vc_velocity_op_spec;
    /*
     * Problem specification object for VCStaggeredStokesProjectionPreconditioner.
     */
    VCStaggeredStokesProjectionPCSpec d_vc_projection_pc_spec;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AcousticStreamingHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AcousticStreamingHierarchyIntegrator(const AcousticStreamingHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AcousticStreamingHierarchyIntegrator& operator=(const AcousticStreamingHierarchyIntegrator& that) = delete;

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * Preprocess the operators and solvers used by the hierarchy integrator.
     */
    void preprocessOperatorsAndSolvers(double current_time, double new_time);

    /*!
     * Update the operators and solvers to account for changes due to time-dependent coefficients
     */
    void updateOperatorsAndSolvers(double current_time, double new_time, int cycle_num);

    /*!
     * Setup solution and RHS vectors using state data maintained by the
     * integrator for the first order system.
     */
    void setupSolverVectorsFOSystem(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& sol1_vec,
                                    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& rhs1_vec,
                                    double current_time,
                                    double new_time,
                                    int cycle_num);

    /*!
     * Setup solution and RHS vectors using state data maintained by the
     * integrator for the second order system.
     */
    void setupSolverVectorsSOSystem(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& sol2_vec,
                                    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& rhs2_vec,
                                    double current_time,
                                    double new_time,
                                    int cycle_num);

    /*!
     * Compute source terms for the second order system arising from the first-order solution.
     */
    void computeCouplingSourceTerms(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& sol1_vec,
                                    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& rhs2_vec,
                                    double current_time,
                                    double new_time,
                                    int cycle_num);
    /*!
     * Pull the solution out of SAMRAI vectors and reset the solver right hand side vectors.
     */
    void resetSolverVectors(const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& sol1_vec,
                            const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& rhs1_vec,
                            const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& sol2_vec,
                            const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& rhs2_vec,
                            double current_time,
                            double new_time,
                            int cycle_num);

    /*!
     * Compute contour integrals to evaluate acoustic radiation force.
     */
    void computeAcousticRadiationForce(double time);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_AcousticStreamingHierarchyIntegrator
