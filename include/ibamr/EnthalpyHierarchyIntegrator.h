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

#ifndef included_IBAMR_EnthalpyHierarchyIntegrator
#define included_IBAMR_EnthalpyHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/PhaseChangeHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

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
 * \brief Class EnthalpyHierarchyIntegrator is a concrete class that manages the
 * time integration of the energy equation using the enthalpy approach.
 *
 * In this class, implementations are provided for phase change of a pure material
 * using the phase-field method.
 *
 *
 *
 */
class EnthalpyHierarchyIntegrator : public PhaseChangeHierarchyIntegrator
{
public:
    /*!
     * The constructor for class EnthalpyHierarchyIntegrator sets
     * some default values, reads in configuration information from input and
     * restart databases, and registers the integrator object with the restart
     * manager when requested.
     */
    EnthalpyHierarchyIntegrator(const std::string& object_name,
                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                bool register_for_restart = true);

    /*!
     * The destructor for class EnthalpyHierarchyIntegrator
     * unregisters the integrator object with the restart manager when the
     * object is so registered.
     */
    ~EnthalpyHierarchyIntegrator() = default;

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
     * Reset cached hierarchy dependent data. This is not required. So does the PhaseChangeHI resetHierarchy function is
     * called?
     */
    // void
    // resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    //                                        int coarsest_level,
    //                                        int finest_level) override;

    /*!
     * Supply an IBTK::CartGridFunction object to specify the value of a
     * particular source term for energy equation.
     */
    void computeEnergyEquationSourceTerm(int F_scratch_idx, const double dt) override;

    /*!
     * Set an object to provide boundary conditions for \f$ h \f$ variable,
     * that has been registered with the hierarchy integrator.
     */
    void setEnthalpyBcCoef(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* h_bc_coef);

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
     * \brief Function to reset the dF/dT if they are
     * maintained by this integrator.
     */
    using ResetLiquidFractionDerivativeFcnPtr = void (*)(double& dlf_dT_data, const double lf_data, void* ctx);

    /*!
     * \brief Function to reset the dF/dT if they are
     * maintained by this integrator.
     */
    using ResetLiquidFractionInverseFcnPtr = void (*)(double& lf_inverse_data, const double lf_data, void* ctx);

    /*!
     * \brief Register function to reset fluid viscosity.
     */
    void registerResetLiquidFractionFcn(ResetLiquidFractionFcnPtr callback, void* ctx);

    /*!
     * \brief Register function to reset fluid viscosity.
     */
    void registerResetLiquidFractionDerivativeFcn(ResetLiquidFractionDerivativeFcnPtr callback, void* ctx);

    /*!
     * \brief Register function to reset fluid viscosity.
     */
    void registerResetLiquidFractionInverseFcn(ResetLiquidFractionInverseFcnPtr callback, void* ctx);

    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    EnthalpyHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    EnthalpyHierarchyIntegrator(const EnthalpyHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    EnthalpyHierarchyIntegrator& operator=(const EnthalpyHierarchyIntegrator& that) = delete;

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * Compute the RHS of continuity equation.
     */
    void computeContinuitySourceTerm(int div_U_F_idx,
                                     const int T_new_idx,
                                     const int h_new_idx,
                                     const int rho_new_idx,
                                     const int T_diff_coef_idx,
                                     const int H_new_idx,
                                     const double new_time);

    /*!
     * compute liquid fraction.
     */
    void computeLiquidFractionRelativeError(const int lf_new_idx, int lf_pre_idx);

    /*!
     * compute enthalpy based on h-T relation.
     */
    void computeEnthalpyBasedOnNonLinearTemperature(int h_idx,
                                                    const int T_idx,
                                                    const int rho_idx,
                                                    const int lf_idx,
                                                    const int H_idx);

    /*!
     * compute temperature based on h-T relation.
     */
    void computeTemperatureBasedOnNonLinearEnthalpy(int T_idx, const int h_idx, const int H_idx);

    /*!
     * Update enthalpy.
     */
    void updateEnthalpy(int h_idx, const int T_new_idx, const int T_pre_idx);

    /*!
     * \brief compute dh/dT based on temperature.
     */
    void computeEnthalpyDerivative(int dh_dT_data, const int T_idx, const int H_idx);

    /*!
     * Compute liquid fraction.
     */
    void computeLiquidFraction(int lf_idx, const int h_idx, const int H_idx);

    /*!
     * Solver variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_h_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_material_derivative_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_grad_T_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_pre_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_dh_dT_var;

    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_h_bc_coef = nullptr;

    /*!
     * Solvers and related data.
     */

    std::vector<ResetLiquidFractionFcnPtr> d_reset_liquid_fraction_fcns;
    std::vector<void*> d_reset_liquid_fraction_fcns_ctx;

    std::vector<ResetLiquidFractionDerivativeFcnPtr> d_reset_dlf_dT_fcns;
    std::vector<void*> d_reset_dlf_dT_fcns_ctx;

    std::vector<ResetLiquidFractionInverseFcnPtr> d_reset_lf_inverse_fcns;
    std::vector<void*> d_reset_lf_inverse_fcns_ctx;

    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_rho_bc_coef;

    /*!
     * Patch data descriptor indices for all "state" variables managed by the
     * integrator.
     *
     * State variables have three contexts: current, scratch, and new.
     */
    int d_h_scratch_idx = IBTK::invalid_index, d_h_current_idx = IBTK::invalid_index, d_h_new_idx = IBTK::invalid_index;

    /*!
     * Patch data descriptor indices for all "scratch" variables managed by the
     * integrator.
     *
     * Scratch variables have only one context: scratch.
     */
    int d_grad_T_idx = IBTK::invalid_index;
    int d_T_pre_idx = IBTK::invalid_index;
    int d_dh_dT_scratch_idx = IBTK::invalid_index;

    /*!
     * Energy equation parameters.
     */
    double d_liquidus_temperature, d_solidus_temperature, d_rho_solid, d_cp_liquid, d_cp_solid, d_cp_gas;

    /*!
     * Inner iteration parameters.
     */
    int d_max_inner_iterations = 2, d_num_cycles;
    double d_lf_iteration_error_tolerance = 1e-6;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_EnthalpyHierarchyIntegrator
