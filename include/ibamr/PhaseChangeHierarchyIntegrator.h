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

#include <ibamr/config.h>

#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/MassIntegrator.h"

#include "ibtk/CCPoissonSolverManager.h"
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
 * \brief Class PhaseChangeHierarchyIntegrator provides an abstract interface
 * for the time integration of a cell-centered energy and Allen-Cahn equations with
 * variable material properties.
 *
 * Liquid fraction $\f \varphi \f$, Temperature $\f \varphi \f$ and the density $\f \rho \f$ are
 * maintained by this integrator.
 *
 * To have a consistency in transport of temperature (or enthalpy), conservation of mass equation is solved
 * discretly with the use of AdvDiffConservativeMassScalarTransportRKIntegrator class to get the newest
 * density \f$\rho^{n+1}\f$ and the convective flux \f$\rho u\f$. This convective flux will be used in the
 * convective operator \f$\nabla \cdot (\rho u C_p T)\f$ (or \f$\nabla \cdot (\rho u h)\f$).
 * This leads to stable solutions for high-density flows.
 *
 * Currenlty, the implementations are provided for the phase-change of a pure material using the enthalpy approach
 * and the phase-field method.
 *
 * References
 * To be added.
 */
class PhaseChangeHierarchyIntegrator : public AdvDiffSemiImplicitHierarchyIntegrator
{
public:
    /*!
     * The constructor for class PhaseChangeHierarchyIntegrator sets
     * some default values, reads in configuration information from input and
     * restart databases, and registers the integrator object with the restart
     * manager when requested.
     */
    PhaseChangeHierarchyIntegrator(const std::string& object_name,
                                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                   bool register_for_restart = true);

    /*!
     * The destructor for class PhaseChangeHierarchyIntegrator
     * unregisters the integrator object with the restart manager when the
     * object is so registered.
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
     * Register the specific heat \f$ C_p \f$.
     */
    void registerSpecificHeatVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Cp_var,
                                      const bool output_Cp = false);

    /*!
     * Register the cell-centered density \f$ \rho \f$.
     */
    void registerDensityVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > rho_var,
                                 const bool output_rho = false);

    /*!
     * \brief Function to reset fluid density, specific heat and conductivity if they are
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
     * \brief Register function to reset specific heat.
     */
    void registerResetSpecificHeatFcn(ResetFluidPropertiesFcnPtr callback, void* ctx);

    /*!
     * \brief Register function to reset thermal conductivity.
     */
    void registerResetDiffusionCoefficientFcn(ResetFluidPropertiesFcnPtr callback, void* ctx);

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
     * \brief set Heaviside variable \f$ H \f$ maintained by AdvDiffHierarchyIntegrator.
     */
    void setHeavisideVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > H_var);

    /*!
     * \brief Register Temperature variable \f$ T \f$.
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
     * Set a grid function to provide initial conditions for \f$ \varphi \f$ variable, that has
     * been registered with the hierarchy integrator.
     */
    void setLiquidFractionInitialCondition(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                           SAMRAI::tbox::Pointer<IBTK::CartGridFunction> lf_init);
    /*!
     * Set a grid function to provide initial conditions for  \f$ T \f$ variable,
     * that has been registered with the hierarchy integrator.
     */
    void setTemperatureInitialCondition(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var,
                                        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> T_init);

    /*!
     * Set a grid function to provide initial conditions for  \f$ \rho \f$ variable,
     * that has been registered with the hierarchy integrator.
     */
    void setDensityInitialCondition(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > rho_var,
                                    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> rho_init);

    /*!
     * Set an object to provide boundary conditions for \f$ T \f$ variable,
     * that has been registered with the hierarchy integrator.
     */
    void setTemperaturePhysicalBcCoef(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var,
                                      SAMRAI::solv::RobinBcCoefStrategy<NDIM>* T_bc_coef);
    /*!
     * Return a \f$ T \f$ boundary condition object.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* getTemperaturePhysicalBcCoef();

    /*!
     *  Compute the source terms of the energy equation.
     */
    virtual void computeEnergyEquationSourceTerm(int F_scratch_idx, const double dt) = 0;

    /*!
     * Get the solver for the energy equation.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver>
    getHelmholtzSolverForEnergyEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > var);

    /*!
     * Get the operator to use to evaluate the right-hand side of the energy equation.
     */
    SAMRAI::tbox::Pointer<IBTK::LaplaceOperator>
    getHelmholtzRHSOperatorForEnergyEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > var);

    /*!
     * Supply an IBTK::CartGridFunction object to specify the value of a
     * source term for the energy equation.
     */
    void setEnergyEquationSourceTermFunction(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> F_fcn);

    /*!
     * Get the updated density patch data index.
     */
    int getUpdatedDensityIndex();

    /*!
     * Get the continuity equation source term patch data index.
     */
    int getContinuityEquationSourceTermIndex();

    /*
     * \brief Supply boundary conditions for the cell-centered density field.
     */
    void registerMassDensityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>*& rho_bc_coefs);

    /*!
     * \brief Supply a source term for the mass update equation.
     *
     * \note Current implementation is used only to check order of accuracy via a manufactured solution.
     */
    void registerMassDensitySourceTerm(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> S_fcn);

    /*
     * \brief Supply boundary conditions for the cell-centered specific heat.
     */
    void registerSpecificHeatBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>*& Cp_bc_coefs);

    /*
     * \brief Supply boundary conditions for the cell-centered thermal conductivity.
     */
    void registerThermalConductivityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>*& k_bc_coefs);

    /*!
     * Set the face-centered advection velocity.
     */
    void setAdvectionVelocity(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > u_var);

    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

protected:
    /*
     * \brief Interpolate the cell-centered scalar variable to side-centered using simple averaging.
     */
    void interpolateCCToSC(int sc_idx, const int cc_idx);

    /*
     * \brief Interpolate the cell-centered scalar variable to side-centered using harmonic averaging.
     */
    void interpolateCCToSCHarmonic(int sc_idx, const int cc_idx);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

protected:
    /*!
     * Bound the liquid fraction, if necessary.
     */
    void boundLiquidFraction(int lf_new_idx);

    /*!
     * Solver variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_F_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_T_diffusion_coef_var, d_T_diffusion_coef_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_diffusion_coef_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_N_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_C_var, d_T_temp_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Cp_var, d_rho_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_C_var, d_rho_vec_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_H_var, d_H_pre_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_D_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_updated_rho_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_div_U_F_var;
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
    bool d_output_lf = false, d_output_T = false, d_output_rho = false, d_output_Cp = false;
    ;

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
     * Functions to reset the \f$ \rho \f$, \f$ C_p \f$ and \f$ k \f$.
     */
    std::vector<ResetFluidPropertiesFcnPtr> d_reset_rho_fcns, d_reset_Cp_fcns, d_reset_kappa_fcns;
    std::vector<void*> d_reset_rho_fcns_ctx, d_reset_Cp_fcns_ctx, d_reset_kappa_fcns_ctx;

    /*!
     * Boundary conditions object for \f$ \rho \f$, \f$ C_p \f$, \f$ k \f$ and \f$ T \f$.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_rho_bc_coef = nullptr;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_k_bc_coef = nullptr;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_Cp_bc_coef = nullptr;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_T_bc_coef = nullptr;

    /*!
     * Boolean to identify whether the mass conservation is to solved.
     */
    bool d_solve_mass_conservation = true;

    /*
     * Pointer to mass integrator class.
     */
    SAMRAI::tbox::Pointer<IBAMR::MassIntegrator> d_rho_p_integrator;

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
    int d_T_scratch_idx = IBTK::invalid_index, d_T_current_idx = IBTK::invalid_index, d_T_new_idx = IBTK::invalid_index;
    int d_rho_scratch_idx = IBTK::invalid_index, d_rho_current_idx = IBTK::invalid_index,
        d_rho_new_idx = IBTK::invalid_index;
    int d_Cp_scratch_idx = IBTK::invalid_index, d_Cp_current_idx = IBTK::invalid_index,
        d_Cp_new_idx = IBTK::invalid_index;
    int d_C_scratch_idx = IBTK::invalid_index, d_C_current_idx = IBTK::invalid_index, d_C_new_idx = IBTK::invalid_index;
    int d_H_scratch_idx = IBTK::invalid_index, d_H_current_idx = IBTK::invalid_index, d_H_new_idx = IBTK::invalid_index;
    int d_U_old_current_idx = IBTK::invalid_index, d_U_old_new_idx = IBTK::invalid_index,
        d_U_old_scratch_idx = IBTK::invalid_index;
    int d_D_cc_scratch_idx = IBTK::invalid_index, d_D_cc_current_idx = IBTK::invalid_index,
        d_D_cc_new_idx = IBTK::invalid_index;

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
    int d_div_U_F_idx = IBTK::invalid_index;

    /*!
     * Phase change parameters.
     */
    double d_rho_liquid, d_T_melt, d_latent_heat;

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

#endif //#ifndef included_IBAMR_PhaseChangeHierarchyIntegrator
