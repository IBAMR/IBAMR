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

#ifndef included_IBAMR_IEPSemiImplicitHierarchyIntegrator
#define included_IBAMR_IEPSemiImplicitHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/AdvDiffConservativeMassTransportQuantityIntegrator.h"
#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/LaplaceOperator.h"

#include "HierarchyFaceDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <map>
#include <set>
#include <string>

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
 * \brief Class IEPSemiImplicitHierarchyIntegrator manages the spatial
 * discretization and time integration of scalar- and vector-valued quantities
 * whose dynamics are governed by the energy equation with phase-change.
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
class IEPSemiImplicitHierarchyIntegrator : public AdvDiffSemiImplicitHierarchyIntegrator
{
public:
    /*!
     * The constructor for class AdvDiffSemiImplicitHierarchyIntegrator sets
     * some default values, reads in configuration information from input and
     * restart databases, and registers the integrator object with the restart
     * manager when requested.
     */
    IEPSemiImplicitHierarchyIntegrator(const std::string& object_name,
                                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                       bool register_for_restart = true);

    /*!
     * The destructor for class AdvDiffSemiImplicitHierarchyIntegrator
     * unregisters the integrator object with the restart manager when the
     * object is so registered.
     */
    ~IEPSemiImplicitHierarchyIntegrator() = default;

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
    void integrateHierarchy(double current_time, double new_time, int cycle_num = 0) override;

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override;

    /*!
     * Register the specific heat variable.
     */
    void registerSpecificHeatVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Cp_var,
                                      const bool output_Cp = false);

    /*!
     * Register the density variable.
     */
    void registerDensityVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > rho_var,
                                 const bool output_rho = false);

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
     * \brief Register function to reset fluid viscosity.
     */
    void registerResetSpecificHeatFcn(ResetFluidPropertiesFcnPtr callback, void* ctx);

    /*!
     * \brief Register function to reset fluid viscosity.
     */
    void registerResetDiffusionCoefficientFcn(ResetFluidPropertiesFcnPtr callback, void* ctx);

    /*!
     * Register INSVCStaggeredHierarchyIntegrator class which will be
     * used to get the variables maintained by this hierarchy integrator.
     */
    //    void registerINSVCStaggeredHierarchyIntegrator(
    //        SAMRAI::tbox::Pointer<INSVCStaggeredHierarchyIntegrator> ins_hier_integrator);

    /*!
     * \brief Get the liquid fraction variable that is being manintained by an advection-diffusion integrator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > getLiquidFractionVariable() const;

    /*!
     * \brief Get the liquid fraction variable that is being manintained by an advection-diffusion integrator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > getHeavisideVariable() const;

    /*
     * \brief Register boundary condition for the density field, if maintained by an advection-diffusion
     * integrator.
     */

    void registerLiquidFractionBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* fl_bc_coef);

    /*!
     * \brief Get the liquid fraction boundary conditions.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* getLiquidFractionBoundaryConditions() const;

    /*!
     * \brief Register liquid fraction variable.
     */
    void registerLiquidFractionVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                        const bool output_lf_var = true);

    /*!
     * \brief Register liquid fraction variable.
     */
    void registerHeavisideVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > H_var,
                                   const bool output_H_var = true);

    /*!
     * \brief Register Temperature variable.
     */
    void registerTemperatureVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var,
                                     const bool output_T_var = true);

    /*!
     * Reset cached hierarchy dependent data.
     */
    void
    resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                           int coarsest_level,
                                           int finest_level) override;

    /*!
     * Set a grid function to provide initial conditions for \f$ k \f$ variable, that has
     * been registered with the hierarchy integrator.
     */
    void
    setInitialConditionsLiquidFractionEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                               SAMRAI::tbox::Pointer<IBTK::CartGridFunction> lf_init);
    /*!
     * Set a grid function to provide initial conditions for  \f$ \omega \f$ variable,
     * that has been registered with the hierarchy integrator.
     */
    void setInitialConditionsTemperatureEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var,
                                                 SAMRAI::tbox::Pointer<IBTK::CartGridFunction> T_init);

    /*!
     * Set a grid function to provide initial conditions for  \f$ \omega \f$ variable,
     * that has been registered with the hierarchy integrator.
     */
    void setDensityInitialCondition(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > rho_var,
                                    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> rho_init);

    /*!
     * Set an object to provide boundary conditions for  \f$ k \f$ variable,
     * that has been registered with the hierarchy integrator.
     */
    void
    setPhysicalBcCoefLiquidFractionEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* lf_bc_coef);

    /*!
     * Return a \f$ k \f$ boundary condition object.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* getPhysicalBcCoefLiquidFractionEquation();

    /*!
     * Set an object to provide boundary conditions for \f$ \omega \f$ variable,
     * that has been registered with the hierarchy integrator.
     */
    void setPhysicalBcCoefTemperatureEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var,
                                              SAMRAI::solv::RobinBcCoefStrategy<NDIM>* T_bc_coef);
    /*!
     * Return a \f$ \omega \f$ boundary condition object.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* getPhysicalBcCoefTemperatureEquation();

    /*!
     * Get the convective operator being used by this solver class for \f$ lf \f$ variable.
     *
     * If the convective operator has not already been constructed, then this
     * function will initialize a default convective operator.
     */
    SAMRAI::tbox::Pointer<ConvectiveOperator> getConvectiveOperatorLiquidFractionEquation(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var);

    /*!
     * Get the convective operator being used by this solver class for \f$ T \f$ variable.
     *
     * If the convective operator has not already been constructed, then this
     * function will initialize a default convective operator.
     */
    SAMRAI::tbox::Pointer<ConvectiveOperator>
    getConvectiveOperatorTemperatureEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var);

    /*!
     * Get the chemical potential patch data index.
     */
    int getChemicalPotentialIndex();

    /*
     * \brief Supply boundary conditions for the cell-centered density field, which is maintained by this integrator
     */
    void registerMassDensityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>*& rho_cc_bc_coefs);

    /*!
     * \brief Supply a source term for the mass update equation.
     *
     * \note Current implementation is used only to check order of accuracy via a manufactured solution.
     */
    void registerMassDensitySourceTerm(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> S_fcn);

    /*!
     * Set the face-centered advection velocity to be used with a cell-centered liquid fraction variable.
     *
     * \note The specified advection velocity must have been already registered
     * with the hierarchy integrator.
     */
    void
    setAdvectionVelocityLiquidFractionEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                               SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > u_var);
    /*!
     * Set the face-centered advection velocity to be used with a cell-centered \f$ \omega \f$.
     * variable.
     *
     * \note The specified advection velocity must have been already registered
     * with the hierarchy integrator.
     */
    void
    setAdvectionVelocityTemperatureEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var,
                                            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > u_var);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IEPSemiImplicitHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IEPSemiImplicitHierarchyIntegrator(const AdvDiffSemiImplicitHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IEPSemiImplicitHierarchyIntegrator& operator=(const AdvDiffSemiImplicitHierarchyIntegrator& that) = delete;

    /*!
     * \brief Compute and store the Heaviside function value on the cell-centers.
     */
    void computeHeavisideFunction(int H_idx, const int phi_idx);

    /*!
     * \brief Compute and store the double-well potential on the cell-centers based on the
     * liquid fraction value.
     */
    void computeDoubleWellPotential(int g_firstder_idx, int g_secondder_idx, const int liquid_fraction_idx);

    /*!
     * \brief Compute and store an interpolation function on the cell-centers based on the
     * liquid fraction value.
     */
    void computeInterpolationFunction(int p_firstder_idx, const int liquid_fraction_idx, const int T_idx);

    /*
     * \brief Compute source term of liquid fraction equation.
     */
    void computeLiquidFractionSourceTerm(int F_scratch_idx);

    /*
     * \brief Compute source term of liquid fraction equation.
     */
    void computeTemperatureSourceTerm(int F_scratch_idx, const double dt);

    /*
     * \brief Interpolate the cell-centered heaviside function to side-centered.
     */
    void interpolateCCHeaviside(int lf_diff_coef_idx, const int H_idx);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db, bool is_from_restart);

    /*!
     * Get the solver for the \f$ k \f$ equation.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver>
    getHelmholtzSolverLiquidFractionEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var);

    /*!
     * Get the solver for the \f$ \omega \f$ equation.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver>
    getHelmholtzSolverTemperatureEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var);

    /*!
     * Get the operator to use to evaluate the right-hand side of the \f$ k \f$ equation.
     */
    SAMRAI::tbox::Pointer<IBTK::LaplaceOperator> getHelmholtzRHSOperatorLiquidFractionEquation(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var);

    /*!
     * Get the operator to use to evaluate the right-hand side of the \f$ \omega \f$ equation.
     */
    SAMRAI::tbox::Pointer<IBTK::LaplaceOperator>
    getHelmholtzRHSOperatorTemperatureEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var);

    /*!
     * Compute the chemical potential of the Allen-Cahn equation using updated liquid fraction value at
     * the cell-centers.
     */
    void computeChemicalPotential(int chemical_potential_idx, const int H_new_idx, const double new_time);

    /*!
     * Compute LHS of AC equation.
     */
    void computeLHSOfLiquidFractionEquation(int lf_lhs_idx,
                                            const int lf_N_scratch_idx,
                                            const double dt,
                                            const int cycle_num,
                                            const double new_time,
                                            const double current_time,
                                            const double half_time);

    /*!
     * Additional variables required.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_var;

    int d_ls_scratch_idx = IBTK::invalid_index, d_ls_current_idx = IBTK::invalid_index,
        d_ls_new_idx = IBTK::invalid_index;
    int d_lf_scratch_idx = IBTK::invalid_index, d_lf_current_idx = IBTK::invalid_index,
        d_lf_new_idx = IBTK::invalid_index;
    int d_T_scratch_idx = IBTK::invalid_index, d_T_current_idx = IBTK::invalid_index, d_T_new_idx = IBTK::invalid_index;

    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_ls_init, d_lf_init, d_T_init, d_rho_init;

    /*!
     * Source function variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_F_var, d_T_F_var;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_lf_F_fcn, d_T_F_fcn;

    /*!
     * Diffusion coefficient data
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_lf_diffusion_coef_var,
        d_lf_diffusion_coef_rhs_var, d_T_diffusion_coef_var, d_T_diffusion_coef_rhs_var;

    bool d_output_ls = false, d_output_H = false, d_output_lf = false, d_output_T = false;

    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_lf_F_init, d_T_F_init;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_lf_bc_coef = nullptr;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_T_bc_coef = nullptr;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_H_bdry_bc_fill_op;

    /*!
     * Advection velocity variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_div_u_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_H_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_lf_u_var, d_T_u_var;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_lf_u_fcn, d_T_u_fcn;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_N_var, d_T_N_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_N_old_var, d_T_N_old_var;

    ConvectiveDifferencingType d_lf_convective_difference_form, d_T_convective_difference_form;
    std::string d_lf_convective_op_type, d_T_convective_op_type;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_lf_convective_op_input_db, d_T_convective_op_input_db;
    SAMRAI::tbox::Pointer<ConvectiveOperator> d_lf_convective_op, d_T_convective_op;
    bool d_lf_convective_op_needs_init, d_T_convective_op_needs_init;

    /*!
     * Solvers and related data.
     */
    std::string d_lf_solver_type, d_lf_precond_type, d_T_solver_type, d_T_precond_type;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_lf_solver_db, d_lf_precond_db, d_T_solver_db, d_T_precond_db;
    SAMRAI::tbox::Pointer<IBTK::LaplaceOperator> d_lf_rhs_op, d_T_rhs_op;
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_lf_solver, d_T_solver;
    bool d_lf_solver_needs_init, d_T_solver_needs_init, d_lf_rhs_op_needs_init, d_T_rhs_op_needs_init;

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_lf_sol, d_lf_rhs, d_T_sol, d_T_rhs;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_rhs_var, d_T_rhs_var;

    TimeSteppingType d_lf_diffusion_time_stepping_type, d_lf_convective_time_stepping_type,
        d_lf_init_convective_time_stepping_type, d_T_diffusion_time_stepping_type, d_T_convective_time_stepping_type,
        d_T_init_convective_time_stepping_type;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_C_var, d_T_C_var, d_lf_temp_rhs_var,
        d_T_temp_rhs_var;
    int d_lf_temp_rhs_idx = IBTK::invalid_index, d_T_temp_rhs_idx = IBTK::invalid_index,
        d_lf_C_idx = IBTK::invalid_index, d_T_C_idx = IBTK::invalid_index;
    double d_M_lf, d_lambda_lf, d_eta_lf;

    bool d_rho_output, d_Cp_output;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Cp_var, d_rho_var;

    std::vector<ResetFluidPropertiesFcnPtr> d_reset_rho_fcns, d_reset_Cp_fcns, d_reset_kappa_fcns;
    std::vector<void*> d_reset_rho_fcns_ctx, d_reset_Cp_fcns_ctx, d_reset_kappa_fcns_ctx;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_C_var, d_rho_vec_cc_var;

    int d_C_scratch_idx = IBTK::invalid_index, d_C_current_idx = IBTK::invalid_index, d_C_new_idx = IBTK::invalid_index,
        d_C_rhs_scratch_idx = IBTK::invalid_index, d_rho_vec_cc_current_idx = IBTK::invalid_index,
        d_rho_vec_cc_scratch_idx = IBTK::invalid_index, d_rho_vec_cc_new_idx = IBTK::invalid_index;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_H_var;
    int d_H_scratch_idx = IBTK::invalid_index, d_H_current_idx = IBTK::invalid_index, d_H_new_idx = IBTK::invalid_index;

    double d_rho_liquid, d_T_ref, d_latent_heat;

    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_rho_bc_coef;

    /*!
     * Additional variables for phase change.
     */
    bool d_phase_change = false;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_dh_var;
    int d_dh_scratch_idx = IBTK::invalid_index;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_D_cc_var;
    int d_D_cc_scratch_idx = IBTK::invalid_index, d_D_cc_current_idx = IBTK::invalid_index,
        d_D_cc_new_idx = IBTK::invalid_index;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_g_firstder_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_g_secondder_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_p_firstder_var;

    int d_g_firstder_idx = IBTK::invalid_index, d_g_secondder_idx = IBTK::invalid_index,
        d_p_firstder_idx = IBTK::invalid_index;

    // Number of interface cells to compute the Heaviside function
    double d_num_interface_cells = 2.0;

    bool d_solve_energy = false, d_solve_mass_conservation = true;

    // Variables for chemical potential.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_chemical_potential_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_grad_lf_var;
    int d_chemical_potential_idx, d_grad_lf_idx, d_H_sc_idx;

    /*
     * Conservative density and transport quantity integrator.
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffConservativeMassTransportQuantityIntegrator> d_rho_p_integrator;

    /*
     * Source term function for the mass density update
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_S_fcn;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_U_old_var;
    int d_U_old_current_idx, d_U_old_new_idx, d_U_old_scratch_idx;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_cp_old_var;
    int d_cp_old_current_idx, d_cp_old_new_idx, d_cp_old_scratch_idx;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_old_var;
    int d_T_old_current_idx, d_T_old_new_idx, d_T_old_scratch_idx;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_lhs_var, d_lf_lhs_N_var;
    int d_lf_lhs_idx, d_lf_lhs_N_scratch_idx;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_lf_N_var;
    int d_T_lf_N_scratch_idx;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IEPSemiImplicitHierarchyIntegrator
