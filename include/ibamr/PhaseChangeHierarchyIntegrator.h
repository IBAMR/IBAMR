// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
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

#ifndef included_IBAMR_PhaseChangeHierarchyIntegrator
#define included_IBAMR_PhaseChangeHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/STSMassFluxIntegrator.h"

#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/LaplaceOperator.h"

namespace IBTK
{
class CartGridFunction;
class LaplaceOperator;
class PoissonSolver;
} // namespace IBTK

namespace IBAMR
{
class ConvectiveOperator;
} // namespace IBAMR

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchHierarchy;
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace mesh
{
template <int DIM>
class GriddingAlgorithm;
} // namespace mesh
namespace pdat
{
template <int DIM, class TYPE>
class CellVariable;
} // namespace pdat
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class PhaseChangeHierarchyIntegrator provides an abstract interface
 * for the time integration of energy and phase change (melting ad solidification)
 * equations (e.g., Allen-Cahn, enthalpy) for a phase changing material (PCM).
 * The solid and liquid phases of the PCM are allowed to have different thermophysical properties.
 *
 * Liquid fraction $\f \varphi \f$, temperature $\f T \f$, and density $\f \rho \f$ variables
 * are maintained by this integrator.
 *
 * To enable consistentcy in the mass and energy transport, an additional mass balance PDE
 * is integrated with the help of AdvDiffConservativeMassScalarTransportRKIntegrator class that
 * provides an approximation to the latest density \f$\rho^{n+1}\f$ and convective flux \f$\rho u\f$.
 * The discrete convective flux is used in the convective operator \f$\nabla \cdot (\rho u C_p T)\f$ (or \f$\nabla \cdot
 * (\rho u h)\f$). This stabilizes the numerical solution for high-density ratio flows.
 *
 * Currenlty, two concrete implementations are provided to model the phase change phenomena:
 * (1) the enthalpy approach, and (2) the Allen-Cahn phase-field method.
 *
 * Reference
 * Thirumalaisamy and Bhalla, <A HREF="https://www.sciencedirect.com/science/article/pii/S0301932223002252?via%3Dihub">
 * A low Mach enthalpy method to model non-isothermal gas–liquid–solid flows with melting and solidification</A>
 */
class PhaseChangeHierarchyIntegrator : public AdvDiffSemiImplicitHierarchyIntegrator
{
public:
    /*!
     * The constructor of the PhaseChangeHierarchyIntegrator class sets
     * some default values, reads in configuration information from input and
     * restart databases, and registers the integrator object with the restart
     * manager when requested.
     */
    PhaseChangeHierarchyIntegrator(const std::string& object_name,
                                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                   bool register_for_restart = true);

    /*!
     * The destructor of the PhaseChangeHierarchyIntegrator class
     * unregisters the integrator object with the restart manager when the
     * object is destroyed.
     */
    ~PhaseChangeHierarchyIntegrator() = default;

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
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1) override;

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    // void integrateHierarchy(double current_time, double new_time, int cycle_num = 0) override;

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override;

    /*!
     * Register the specific heat \f$ C_p \f$ variable.
     */
    void
    registerSpecificHeatVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > specific_heat_var,
                                 const bool output_Cp = false);

    /*!
     * Register the density \f$ \rho \f$ variable.
     */
    void registerDensityVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > rho_var,
                                 const bool output_rho = false);

    /*!
     * \brief Function to reset the phase density, specific heat and conductivity if they are
     * maintained by this integrator.
     */
    using ResetPhasePropertiesFcnPtr = void (*)(int property_idx,
                                                SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > property_var,
                                                SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                int cycle_num,
                                                double time,
                                                double current_time,
                                                double new_time,
                                                void* ctx);

    /*!
     * \brief Register function to reset phase density.
     */
    void registerResetDensityFcn(ResetPhasePropertiesFcnPtr callback, void* ctx);

    /*!
     * \brief Register function to reset phase specific heat.
     */
    void registerResetSpecificHeatFcn(ResetPhasePropertiesFcnPtr callback, void* ctx);

    /*!
     * \brief Register function to reset phase thermal conductivity.
     */
    void registerResetDiffusionCoefficientFcn(ResetPhasePropertiesFcnPtr callback, void* ctx);

    /*!
     * \brief Get the liquid fraction variable that is being manintained by this integrator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > getLiquidFractionVariable() const;

    /*!
     * \brief Register liquid fraction variable \f$ \varphi \f$.
     */
    virtual void registerLiquidFractionVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                                const bool output_lf_var = true);

    /*!
     * \brief Register gradient of liquid fraction variable.
     */
    virtual void registerLiquidFractionGradientVariable(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_gradient_var,
        const bool output_lf_gradient_var = true);

    /*!
     * \brief Register Heaviside variable \f$ H \f$ maintained by AdvDiffHierarchyIntegrator.
     */
    void registerHeavisideVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > H_var);

    /*!
     * \brief Register temperature variable \f$ T \f$.
     */
    void registerTemperatureVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var,
                                     const bool output_T_var = true);

    /*!
     * \brief Set a grid function to provide initial conditions for \f$ \varphi \f$ variable.
     */
    void setLiquidFractionInitialCondition(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                           SAMRAI::tbox::Pointer<IBTK::CartGridFunction> lf_init);
    /*!
     * \brief Set a grid function to provide initial conditions for  \f$ T \f$ variable.
     */
    void setTemperatureInitialCondition(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var,
                                        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> T_init);

    /*!
     * \brief Set a grid function to provide initial conditions for \f$ \rho \f$ variable.
     */
    void setDensityInitialCondition(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > rho_var,
                                    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> rho_init);

    /*!
     * \brief Set boundary conditions for \f$ T \f$ variable.
     */
    void setTemperaturePhysicalBcCoef(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var,
                                      SAMRAI::solv::RobinBcCoefStrategy<NDIM>* T_bc_coef);
    /*!
     * \brief Return \f$ T \f$ boundary condition object.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* getTemperaturePhysicalBcCoef();

    /*!
     *  \brief Add the temporal and linear terms to the RHS of the energy equation.
     */
    virtual void addTemporalAndLinearTermstoRHSOfEnergyEquation(int F_scratch_idx, const double dt) = 0;

    /*!
     * \brief Get the solver of the energy equation.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver>
    getEnergyEquationHelmholtzSolver(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > var);

    /*!
     * \brief Get the operator that is used to evaluate the right-hand side of energy equation.
     */
    SAMRAI::tbox::Pointer<IBTK::LaplaceOperator>
    getEnergyEquationHelmholtzRHSOperator(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > var);

    /*!
     * \brief Register an IBTK::CartGridFunction object to specify the value of the
     * source term for energy equation.
     */
    void setEnergyEquationSourceTermFunction(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> F_fcn);

    /*!
     * \brief Get the updated density patch data index.
     */
    inline int getUpdatedDensityIndex() const
    {
        return d_rho_new_idx;
    }

    /*!
     * \brief Get the updated specific heat patch data index.
     */
    inline int getUpdatedSpecificHeatIndex() const
    {
        return d_specific_heat_new_idx;
    }

    /*!
     * \brief Get the source term patch data index for the div U equation.
     */
    int getVelocityDivergencePatchDataIndex();

    /*!
     * \brief Compute the source term for the div U equation.
     */
    virtual void computeDivergenceVelocitySourceTerm(int Div_U_F_idx, const double new_time) = 0;

    /*!
     * \brief Register boundary conditions for the density field.
     */
    void registerMassDensityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>*& rho_bc_coefs);

    /*!
     * \brief Register the source term for the mass update equation.
     *
     * \note The current implementation is used only to check the order of accuracy via a manufactured solution.
     */
    void registerMassDensitySourceTerm(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> S_fcn);

    /*!
     * \brief Register boundary conditions for the cell-centered specific heat variable.
     */
    void registerSpecificHeatBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>*& specific_heat_bc_coefs);

    /*!
     * \brief Register boundary conditions for the cell-centered thermal conductivity variable.
     */
    void registerThermalConductivityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>*& k_bc_coefs);

    /*!
     * \brief Set the face-centered advection velocity.
     */
    void setAdvectionVelocity(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > u_var);

    /*!
     * \brief Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

protected:
    /*!
     * \brief Interpolate the cell-centered scalar variable to side-centered using simple averaging.
     */
    void interpolateCCToSCSimpleAveraging(int sc_idx, const int cc_idx);

    /*!
     * \brief Interpolate the cell-centered scalar variable to side-centered using harmonic averaging.
     */
    void interpolateCCToSCHarmonicAveraging(int sc_idx, const int cc_idx);

    /*!
     * \brief Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db, bool is_from_restart);

    /*!
     * \brief Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * Reset cached hierarchy dependent data for solvers and operators before the regridding operation.
     */
    virtual void regridHierarchyBeginSpecialized() override;

    /*!
     * Reset cached hierarchy dependent data for solvers and operators before the regridding operation.
     */
    virtual void regridHierarchyEndSpecialized() override;

    /*!
     * Bound the liquid fraction, if necessary.
     */
    void boundLiquidFraction(int lf_new_idx);

    /*!
     * Solver variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_var, d_lf_gradient_var, d_lf_pre_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_F_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_T_diffusion_coef_var, d_T_diffusion_coef_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_diffusion_coef_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_N_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_C_var, d_T_temp_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_specific_heat_var, d_rho_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_C_var, d_rho_vec_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_H_var, d_H_pre_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_D_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_updated_rho_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Div_U_F_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_u_adv_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_U_old_var;

    /*!
     * Objects to set initial condition for \f$ \varphi \f$, \f$ T \f$ and \f$ \rho \f$.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_lf_init, d_T_init, d_rho_init;

    /*!
     * Cartgrid functions to be used to set the energy equation source term.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_T_F_fcn;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_T_F_init;

    /*!
     * Boolean to output the variables in visit.
     */
    bool d_output_lf = false, d_output_lf_gradient = false, d_output_T = false, d_output_rho = false,
         d_output_Cp = false, d_output_Div_U_F = false, d_output_temp_k = false;

    /*!
     * Data shynchronization operator.
     */
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_H_bdry_bc_fill_op, d_k_bdry_bc_fill_op;

    /*!
     * Convective operator and difference type for the energy equation.
     */
    ConvectiveDifferencingType d_T_convective_difference_form = d_default_convective_difference_form;
    std::string d_T_convective_op_type = d_default_convective_op_type;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_T_convective_op_input_db = d_default_convective_op_input_db;
    SAMRAI::tbox::Pointer<ConvectiveOperator> d_T_convective_op = nullptr;
    bool d_T_convective_op_needs_init = false;

    /*!
     * Hierarchy operators, solvers and related data for the energy equation.
     */
    std::string d_T_solver_type = IBTK::CCPoissonSolverManager::UNDEFINED,
                d_T_precond_type = IBTK::CCPoissonSolverManager::UNDEFINED;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_T_solver_db, d_T_precond_db;
    SAMRAI::tbox::Pointer<IBTK::LaplaceOperator> d_T_rhs_op;
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_T_solver;
    bool d_T_solver_needs_init, d_T_rhs_op_needs_init;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_T_sol, d_T_rhs;

    TimeSteppingType d_T_diffusion_time_stepping_type = d_default_diffusion_time_stepping_type;
    TimeSteppingType d_T_convective_time_stepping_type = d_default_convective_time_stepping_type;
    TimeSteppingType d_T_init_convective_time_stepping_type = d_default_convective_time_stepping_type;

    /*!
     * Functions to reset \f$ \rho \f$, \f$ C_p \f$ and \f$ k \f$.
     */
    std::vector<ResetPhasePropertiesFcnPtr> d_reset_rho_fcns, d_reset_specific_heat_fcns, d_reset_kappa_fcns;
    std::vector<void*> d_reset_rho_fcns_ctx, d_reset_specific_heat_fcns_ctx, d_reset_kappa_fcns_ctx;

    /*!
     * Boundary conditions object for \f$ \rho \f$, \f$ C_p \f$, \f$ k \f$ and \f$ T \f$.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_rho_bc_coef = nullptr;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_k_bc_coef = nullptr;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_specific_heat_bc_coef = nullptr;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_T_bc_coef = nullptr;

    /*!
     * Boolean to identify whether the mass conservation is to solved.
     */
    bool d_solve_mass_conservation = true;

    /*
     * Pointer to mass integrator class.
     */
    SAMRAI::tbox::Pointer<IBAMR::STSMassFluxIntegrator> d_rho_p_integrator;

    /*
     * Source term function for the mass density update.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_S_fcn;

    /*!
     * Patch data descriptor indices for all "state" variables managed by the
     * integrator.
     *
     * State variables have three contexts: current, scratch, and new.
     */
    int d_lf_scratch_idx = IBTK::invalid_index, d_lf_current_idx = IBTK::invalid_index,
        d_lf_new_idx = IBTK::invalid_index;
    int d_lf_gradient_scratch_idx = IBTK::invalid_index, d_lf_gradient_current_idx = IBTK::invalid_index,
        d_lf_gradient_new_idx = IBTK::invalid_index;
    int d_T_scratch_idx = IBTK::invalid_index, d_T_current_idx = IBTK::invalid_index, d_T_new_idx = IBTK::invalid_index;
    int d_rho_scratch_idx = IBTK::invalid_index, d_rho_current_idx = IBTK::invalid_index,
        d_rho_new_idx = IBTK::invalid_index;
    int d_specific_heat_scratch_idx = IBTK::invalid_index, d_specific_heat_current_idx = IBTK::invalid_index,
        d_specific_heat_new_idx = IBTK::invalid_index;
    int d_C_scratch_idx = IBTK::invalid_index, d_C_current_idx = IBTK::invalid_index, d_C_new_idx = IBTK::invalid_index;
    int d_H_scratch_idx = IBTK::invalid_index, d_H_current_idx = IBTK::invalid_index, d_H_new_idx = IBTK::invalid_index;
    int d_U_old_current_idx = IBTK::invalid_index, d_U_old_new_idx = IBTK::invalid_index,
        d_U_old_scratch_idx = IBTK::invalid_index;
    int d_D_cc_scratch_idx = IBTK::invalid_index, d_D_cc_current_idx = IBTK::invalid_index,
        d_D_cc_new_idx = IBTK::invalid_index;
    int d_T_diffusion_coef_current_idx = IBTK::invalid_index, d_T_diffusion_coef_new_idx = IBTK::invalid_index,
        d_T_diffusion_coef_scratch_idx = IBTK::invalid_index;
    int d_T_diffusion_coef_cc_current_idx = IBTK::invalid_index, d_T_diffusion_coef_cc_new_idx = IBTK::invalid_index,
        d_T_diffusion_coef_cc_scratch_idx = IBTK::invalid_index;
    int d_u_adv_current_idx = IBTK::invalid_index, d_u_adv_new_idx = IBTK::invalid_index,
        d_u_adv_scratch_idx = IBTK::invalid_index;
    int d_T_F_current_idx = IBTK::invalid_index, d_T_F_new_idx = IBTK::invalid_index,
        d_T_F_scratch_idx = IBTK::invalid_index;

    /*!
     * Patch data descriptor indices for all "scratch" variables managed by the
     * integrator.
     *
     * Scratch variables have only one context: scratch.
     */
    int d_T_temp_rhs_idx = IBTK::invalid_index, d_T_C_idx = IBTK::invalid_index;
    int d_lf_pre_idx = IBTK::invalid_index, d_H_pre_idx = IBTK::invalid_index;
    int d_C_rhs_scratch_idx = IBTK::invalid_index;
    int d_updated_rho_idx = IBTK::invalid_index;
    int d_Div_U_F_idx = IBTK::invalid_index;
    int d_T_diffusion_coef_rhs_scratch_idx = IBTK::invalid_index;
    int d_T_rhs_scratch_idx = IBTK::invalid_index;
    int d_T_N_scratch_idx = IBTK::invalid_index;

    /*!
     * Phase change parameters.
     */
    double d_rho_liquid, d_rho_solid, d_T_melt, d_latent_heat;

    /*!
     * Variable to indicate the type of interpolation to be done for conductivity.
     */
    IBTK::VCInterpType d_k_vc_interp_type = IBTK::VC_AVERAGE_INTERP;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PhaseChangeHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PhaseChangeHierarchyIntegrator(const PhaseChangeHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PhaseChangeHierarchyIntegrator& operator=(const PhaseChangeHierarchyIntegrator& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_PhaseChangeHierarchyIntegrator
