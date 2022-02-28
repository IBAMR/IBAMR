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
#include "ibamr/MassIntegrator.h"
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
 * discretization and time integration of the Allen-Cahn and energy equation with phase-change.
 * write more description about mass equation and consistency later.
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

    //    /*!
    //     * \brief Get the heaviside variable that is being manintained by an advection-diffusion integrator.
    //     */
    //    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > getHeavisideVariable() const;

    /*
     * \brief Register boundary condition for the liquid fraction.
     */

    void registerLiquidFractionBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* fl_bc_coef);

    /*!
     * \brief Get the liquid fraction boundary conditions.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* getLiquidFractionBoundaryConditions() const;

    /*!
     * \brief Register liquid fraction variable \f$ \varphi \f$.
     */
    void registerLiquidFractionVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                        const bool output_lf_var = true);

    /*!
     * \brief set Heaviside variable.
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
    void
    setInitialConditionsLiquidFractionEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                               SAMRAI::tbox::Pointer<IBTK::CartGridFunction> lf_init);
    /*!
     * Set a grid function to provide initial conditions for  \f$ T \f$ variable,
     * that has been registered with the hierarchy integrator.
     */
    void setInitialConditionsTemperatureEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var,
                                                 SAMRAI::tbox::Pointer<IBTK::CartGridFunction> T_init);

    /*!
     * Set a grid function to provide initial conditions for  \f$ \rho \f$ variable,
     * that has been registered with the hierarchy integrator.
     */
    void setDensityInitialCondition(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > rho_var,
                                    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> rho_init);

    /*!
     * Set an object to provide boundary conditions for  \f$ \varphi \f$ variable,
     * that has been registered with the hierarchy integrator.
     */
    void
    setPhysicalBcCoefLiquidFractionEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* lf_bc_coef);

    //    /*!
    //     * Set an object to provide boundary conditions for  \f$ H \f$ variable,
    //     * that has been registered with the hierarchy integrator.
    //     */
    //    void setPhysicalBcCoefHeavisideEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >
    //    H_var,
    //                                            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* H_bc_coef);

    /*!
     * Return a \f$ \varphi \f$ boundary condition object.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* getPhysicalBcCoefLiquidFractionEquation();

    /*!
     * Set an object to provide boundary conditions for \f$ T \f$ variable,
     * that has been registered with the hierarchy integrator.
     */
    void setPhysicalBcCoefTemperatureEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var,
                                              SAMRAI::solv::RobinBcCoefStrategy<NDIM>* T_bc_coef);
    /*!
     * Return a \f$ T \f$ boundary condition object.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* getPhysicalBcCoefTemperatureEquation();

    /*!
     * Get the convective operator being used by this solver class for \f$ \varphi \f$ variable.
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
     * Supply an IBTK::CartGridFunction object to specify the value of a
     * particular source term for temperature equation.
     */
    void setTemperatureSourceTermFunction(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> F_fcn);

    /*!
     * Get the chemical potential patch data index.
     */
    int getChemicalPotentialIndex();

    /*!
     * Get the liquid fraction material derivative patch data index.
     */
    int getLiquidFractionMaterialDerivativeIndex();

    /*!
     * Get the updated density patch data index.
     */
    int getUpdatedDensityIndex();

    /*
     * \brief Supply boundary conditions for the cell-centered density field, which is maintained by
     * AdvDiffConservativeMassTransportQuantityIntegrator.
     */
    void registerMassDensityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>*& rho_cc_bc_coefs);

    /*!
     * \brief Supply a source term for the mass update equation.
     *
     * \note Current implementation is used only to check order of accuracy via a manufactured solution.
     */
    void registerMassDensitySourceTerm(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> S_fcn);

    /*!
     * Set the face-centered advection velocity to be used with a cell-centered liquid fraction variable
     *\f$ \varphi \f$.
     *
     * \note The specified advection velocity must have been already registered
     * with the hierarchy integrator.
     */
    void
    setAdvectionVelocityLiquidFractionEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                               SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > u_var);
    /*!
     * Set the face-centered advection velocity to be used with a cell-centered \f$ T \f$.
     * variable.
     *
     * \note The specified advection velocity must have been already registered
     * with the hierarchy integrator.
     */
    void
    setAdvectionVelocityTemperatureEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var,
                                            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > u_var);

    /*!
     * \brief Function to reset the liquid fraction (fl) if they are
     * maintained by this integrator.
     */
    using ResetLiquidFractionFcnPtr = void (*)(int lf_idx,
                                               int lf_inverse_idx,
                                               int dlf_dT_idx,
                                               int T_idx,
                                               int H_idx,
                                               SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                               int cycle_num,
                                               double time,
                                               double current_time,
                                               double new_time,
                                               void* ctx);

    /*!
     * \brief Register function to reset fluid viscosity.
     */
    void registerResetLiquidFractionFcn(ResetLiquidFractionFcnPtr callback, void* ctx);

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
    void computeLiquidFractionSourceTerm(int F_scratch_idx, const double dt, const double new_time);

    /*
     * \brief Compute source term of liquid fraction equation.
     */
    void computeTemperatureSourceTerm(int F_scratch_idx, const double dt);

    /*
     * \brief Interpolate the cell-centered heaviside function to side-centered using simple averaging.
     */
    void interpolateCCToSC(int sc_idx, const int cc_idx);

    /*
     * \brief Interpolate the cell-centered heaviside function to side-centered using harmonic averaging.
     */
    void interpolateCCToSCHarmonic(int sc_idx, const int cc_idx);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db, bool is_from_restart);

    /*!
     * Get the solver for the \f$ \varphi \f$ equation.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver>
    getHelmholtzSolverLiquidFractionEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var);

    /*!
     * Get the solver for the \f$ T \f$ equation.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver>
    getHelmholtzSolverTemperatureEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var);

    /*!
     * Get the operator to use to evaluate the right-hand side of the \f$ \varphi \f$ equation.
     */
    SAMRAI::tbox::Pointer<IBTK::LaplaceOperator> getHelmholtzRHSOperatorLiquidFractionEquation(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var);

    /*!
     * Get the operator to use to evaluate the right-hand side of the \f$ T \f$ equation.
     */
    SAMRAI::tbox::Pointer<IBTK::LaplaceOperator>
    getHelmholtzRHSOperatorTemperatureEquation(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > T_var);

    /*!
     * Compute the chemical potential of the Allen-Cahn equation using updated liquid fraction value at
     * the cell-centers.
     */
    void computeChemicalPotential(int chemical_potential_idx, const int H_new_idx, const double new_time);

    /*!
     * Compute the material derivative of liquid fraction
     */
    void
    computeMaterialDerivativeOfLiquidFraction(int lf_material_derivative_idx, const double dt, const double new_time);

    //    /*!
    //     * Compute LHS of AC equation.
    //     */
    //    void computeLHSOfLiquidFractionEquation(int lf_lhs_idx,
    //                                            const int lf_N_scratch_idx,
    //                                            const double dt,
    //                                            const int cycle_num,
    //                                            const double new_time,
    //                                            const double current_time,
    //                                            const double half_time);

    /*!
     * Bound liquid fraction.
     */
    void boundLiquidFraction(int lf_new_idx);

    /*!
     * Solver variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_pre_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_F_var, d_T_F_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_lf_diffusion_coef_var,
        d_lf_diffusion_coef_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_T_diffusion_coef_var, d_T_diffusion_coef_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_diffusion_coef_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_div_u_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_H_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_lf_u_var, d_T_u_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_N_var, d_T_N_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_N_old_var, d_T_N_old_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_rhs_var, d_T_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_C_var, d_lf_temp_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_C_var, d_T_temp_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_temp_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Cp_var, d_rho_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_C_var, d_rho_vec_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_H_var, d_H_pre_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_D_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_g_firstder_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_g_secondder_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_p_firstder_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_chemical_potential_var, d_M_var,
        d_updated_rho_var, d_lf_diffusion_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_material_derivative_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_grad_lf_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_grad_H_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_U_old_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_cp_old_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_old_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_lhs_var, d_lf_lhs_N_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_lf_N_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_inverse_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_dlf_dT_var;

    /*!
     * Objects to set initial condition for \f$ \varphi \f$, \f$ T \f$ and \f$ \rho \f$.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_lf_init, d_T_init, d_rho_init;

    /*!
     * Source functions.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_lf_F_fcn, d_T_F_fcn;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_lf_F_init, d_T_F_init;

    bool d_output_ls = false, d_output_lf = false, d_output_T = false;

    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_lf_bc_coef = nullptr;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_T_bc_coef = nullptr;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_H_bdry_bc_fill_op, d_k_bdry_bc_fill_op;

    /*!
     * Advection velocity variables.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_lf_u_fcn, d_T_u_fcn;

    ConvectiveDifferencingType d_lf_convective_difference_form = d_default_convective_difference_form;
    ConvectiveDifferencingType d_T_convective_difference_form = d_default_convective_difference_form;

    std::string d_lf_convective_op_type = d_default_convective_op_type;
    std::string d_T_convective_op_type = d_default_convective_op_type;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_lf_convective_op_input_db = d_default_convective_op_input_db;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_T_convective_op_input_db = d_default_convective_op_input_db;
    SAMRAI::tbox::Pointer<ConvectiveOperator> d_lf_convective_op = nullptr, d_T_convective_op = nullptr;
    bool d_lf_convective_op_needs_init = false, d_T_convective_op_needs_init = false;

    /*!
     * Solvers and related data.
     */
    std::string d_lf_solver_type, d_lf_precond_type, d_T_solver_type, d_T_precond_type;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_lf_solver_db, d_lf_precond_db, d_T_solver_db, d_T_precond_db;
    SAMRAI::tbox::Pointer<IBTK::LaplaceOperator> d_lf_rhs_op = nullptr, d_T_rhs_op = nullptr;
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_lf_solver = nullptr, d_T_solver = nullptr;
    bool d_lf_solver_needs_init, d_T_solver_needs_init, d_lf_rhs_op_needs_init, d_T_rhs_op_needs_init;

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_lf_sol, d_lf_rhs, d_T_sol, d_T_rhs;

    TimeSteppingType d_lf_diffusion_time_stepping_type = d_default_diffusion_time_stepping_type;
    TimeSteppingType d_lf_convective_time_stepping_type = d_default_convective_time_stepping_type;
    TimeSteppingType d_lf_init_convective_time_stepping_type = d_default_convective_time_stepping_type;
    TimeSteppingType d_T_diffusion_time_stepping_type = d_default_diffusion_time_stepping_type;
    TimeSteppingType d_T_convective_time_stepping_type = d_default_convective_time_stepping_type;
    TimeSteppingType d_T_init_convective_time_stepping_type = d_default_convective_time_stepping_type;

    bool d_apply_brinkman = false, d_add_diffusion = false;

    bool d_rho_output, d_Cp_output;

    std::vector<ResetFluidPropertiesFcnPtr> d_reset_rho_fcns, d_reset_Cp_fcns, d_reset_kappa_fcns;
    std::vector<void*> d_reset_rho_fcns_ctx, d_reset_Cp_fcns_ctx, d_reset_kappa_fcns_ctx;

    std::vector<ResetLiquidFractionFcnPtr> d_reset_liquid_fraction_fcns;
    std::vector<void*> d_reset_liquid_fraction_fcns_ctx;

    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_rho_bc_coef;

    /*!
     * Boolean to identify whether phase change is involved.
     */
    bool d_phase_change = false;

    //    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_dh_var;
    //    int d_dh_scratch_idx = IBTK::invalid_index;

    bool d_solve_energy = false, d_solve_mass_conservation = true;

    /*
     * Conservative mass and transport quantity integrator.
     */
    SAMRAI::tbox::Pointer<IBAMR::MassIntegrator> d_rho_p_integrator;

    /*
     * Source term function for the mass density update
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_S_fcn;

    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_lf_brinkman_db;
    double d_lf_b = 0.0;

    /*!
     * Patch data descriptor indices for all "state" variables managed by the
     * integrator.
     *
     * State variables have three contexts: current, scratch, and new.
     */
    int d_ls_scratch_idx = IBTK::invalid_index, d_ls_current_idx = IBTK::invalid_index,
        d_ls_new_idx = IBTK::invalid_index;
    int d_lf_scratch_idx = IBTK::invalid_index, d_lf_current_idx = IBTK::invalid_index,
        d_lf_new_idx = IBTK::invalid_index, d_lf_pre_idx = IBTK::invalid_index;
    int d_T_scratch_idx = IBTK::invalid_index, d_T_current_idx = IBTK::invalid_index, d_T_new_idx = IBTK::invalid_index;
    int d_T_diff_coef_cc_current_idx = IBTK::invalid_index, d_T_diff_coef_cc_new_idx = IBTK::invalid_index,
        d_T_diff_coef_cc_scratch_idx = IBTK::invalid_index;
    int d_C_scratch_idx = IBTK::invalid_index, d_C_current_idx = IBTK::invalid_index, d_C_new_idx = IBTK::invalid_index;
    int d_H_scratch_idx = IBTK::invalid_index, d_H_current_idx = IBTK::invalid_index, d_H_new_idx = IBTK::invalid_index;
    int d_U_old_current_idx = IBTK::invalid_index, d_U_old_new_idx = IBTK::invalid_index,
        d_U_old_scratch_idx = IBTK::invalid_index;
    int d_cp_old_current_idx = IBTK::invalid_index, d_cp_old_new_idx = IBTK::invalid_index,
        d_cp_old_scratch_idx = IBTK::invalid_index;
    int d_T_old_current_idx = IBTK::invalid_index, d_T_old_new_idx = IBTK::invalid_index,
        d_T_old_scratch_idx = IBTK::invalid_index;
    int d_D_cc_scratch_idx = IBTK::invalid_index, d_D_cc_current_idx = IBTK::invalid_index,
        d_D_cc_new_idx = IBTK::invalid_index;
    int d_dlf_dT_scratch_idx = IBTK::invalid_index, d_dlf_dT_current_idx = IBTK::invalid_index,
        d_dlf_dT_new_idx = IBTK::invalid_index;

    /*!
     * Patch data descriptor indices for all "scratch" variables managed by the
     * integrator.
     *
     * Scratch variables have only one context: scratch.
     */
    int d_lf_temp_rhs_idx = IBTK::invalid_index, d_T_temp_rhs_idx = IBTK::invalid_index,
        d_T_temp_idx = IBTK::invalid_index, d_lf_C_idx = IBTK::invalid_index, d_T_C_idx = IBTK::invalid_index;
    int d_C_rhs_scratch_idx = IBTK::invalid_index;
    int d_H_pre_idx = IBTK::invalid_index;
    int d_g_firstder_idx = IBTK::invalid_index, d_g_secondder_idx = IBTK::invalid_index,
        d_p_firstder_idx = IBTK::invalid_index;
    int d_chemical_potential_idx = IBTK::invalid_index, d_grad_lf_idx = IBTK::invalid_index,
        d_lf_sc_idx = IBTK::invalid_index, d_grad_H_idx = IBTK::invalid_index, d_H_sc_idx = IBTK::invalid_index,
        d_M_idx = IBTK::invalid_index, d_lf_diffusion_idx = IBTK::invalid_index;
    int d_updated_rho_cc_idx = IBTK::invalid_index;
    int d_T_lf_N_scratch_idx = IBTK::invalid_index;
    int d_lf_inverse_scratch_idx = IBTK::invalid_index;
    int d_lf_material_derivative_idx = IBTK::invalid_index;

    /*!
     * Allen-Cahn equation parameters.
     */
    double d_M_lf, d_lambda_lf, d_eta_lf;

    /*!
     * Parameter related to Brinkman term in Allen-Cahn equation.
     */
    double d_beta = 1.0;

    /*!
     * Parameters related to additional diffusion term in Allen-Cahn equation.
     */
    double d_free_parameter = 1.0, d_eps = 1.0e-8, d_H_diffusion_coefficient = 1e-5;

    /*!
     * Energy equation parameters.
     */
    double d_rho_liquid, d_T_melt, d_latent_heat, d_latent_heat_temp;

    /*!
     * Variable to indicate the type of interpolation to be done for conductivity.
     */
    IBTK::VCInterpType d_k_vc_interp_type = IBTK::VC_AVERAGE_INTERP;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IEPSemiImplicitHierarchyIntegrator
