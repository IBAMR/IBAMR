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
#include "ibamr/CellConvectiveOperator.h"
#include "ibamr/PhaseChangeHierarchyIntegrator.h"

namespace IBTK
{
class CartGridFunction;
class LaplaceOperator;
class PoissonSolver;
} // namespace IBTK

namespace IBAMR
{
class CellConvectiveOperator;
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
 * time integration of energy and phase change equations using the enthalpy approach.
 */
class EnthalpyHierarchyIntegrator : public PhaseChangeHierarchyIntegrator
{
public:
    /*!
     * The constructor of EnthalpyHierarchyIntegrator class sets
     * some default values, reads in configuration information from input and
     * restart databases, and registers the integrator object with the restart
     * manager when requested.
     */
    EnthalpyHierarchyIntegrator(const std::string& object_name,
                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                bool register_for_restart = true);

    /*!
     * The destructor of EnthalpyHierarchyIntegrator class
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
     * Add the temporal and linear terms to the RHS of the energy equation.
     */
    void addTemporalAndLinearTermstoRHSOfEnergyEquation(int F_scratch_idx, const double dt) override;

    /*!
     * Compute the source term for the Div U equation.
     */
    void computeDivergenceVelocitySourceTerm(int Div_U_F_idx, const double new_time) override;

    /*!
     * \brief Register specific enthalpy variable \f$ h \f$.
     */
    virtual void
    registerSpecificEnthalpyVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > h_var,
                                     const bool output_h_var = true);

    /*!
     * Set boundary conditions for \f$ h \f$ variable.
     */
    void setEnthalpyBcCoef(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* h_bc_coef);

    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * \brief Register LevelSet variable \f$ \phi \f$ maintained by AdvDiffHierarchyIntegrator.
     */
    void registerLevelSetVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > phi_var);

    /*!
     * \brief Register liquid fraction variable for extrapolation.
     */
    void registerLiquidFractionVariableForExtrapolation(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var);

    // /*!
    //  * \brief Function to compute the advection velocity which is used for liquid fraction extrapolation to gas
    //  region.
    //  */
    // using LiquidFractionExtrapAdvectionVelocityPtr = void (*)(int F_idx,
    //                                             SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
    //                                             int cycle_num,
    //                                             double time,
    //                                             double current_time,
    //                                             double new_time,
    //                                             void* ctx);

    // /*!
    //  * \brief Register function to compute the advection velocity for liquid fraction extrapolation.
    //  */
    // void registerLiquidFractionExtrapAdvectionVelocity(LiquidFractionExtrapAdvectionVelocityPtr callback, void* ctx)

protected:
    /*!
     * Reset cached hierarchy dependent data for solvers and operators before the regridding operation.
     */
    virtual void regridHierarchyEndSpecialized() override;

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
     * \brief Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * \brief Compute enthalpy based on (nonlinear) h-T relation.
     */
    void
    computeEnthalpyBasedOnTemperature(int h_idx, const int T_idx, const int rho_idx, const int lf_idx, const int H_idx);

    /*!
     * \brief Compute temperature based on (nonlinear) h-T relation.
     */
    void computeTemperatureBasedOnEnthalpy(int T_idx, const int h_idx, const int H_idx);

    /*!
     * \brief Update enthalpy using Taylor series \f$ h_{\rm new} = h_{\rm old} + dh/dt (\Delta T) \f$.
     */
    void updateEnthalpy(int h_idx, const int T_new_idx, const int T_pre_idx);

    /*!
     * \brief compute dh/dT based on temperature.
     */
    void computeEnthalpyDerivative(int dh_dT_data, const int T_idx, const int H_idx);

    /*!
     * \brief Compute liquid fraction.
     */
    void computeLiquidFraction(int lf_idx, const int h_idx, const int H_idx);

    /*!
     * \brief Extrapolate liquid fraction from PCM to gas region using constant extrapolation.
     */
    void extrapolateLiquidFractionToGasRegion(int lf_new_idx);

    /*!
     * \brief Compute advection velocity \f$ u_{\rm adv} = H(-\phi) n \f$.
     */
    void computeAdvectionVelocityForExtrapolation(int u_adv_fc_lf_extrap_current_idx,
                                                  const double current_time,
                                                  const double new_time);

    /*!
     * Get the convective operator for solving lf extrapolation equation.
     *
     * If the convective operator has not already been constructed, then this
     * function will initialize a default convective operator.
     */
    SAMRAI::tbox::Pointer<CellConvectiveOperator> getLiquidFractionExtrapConvectiveOperator(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_extrap_var);

    /*!
     * Solver variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_h_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_grad_T_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_pre_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_dh_dT_var;

    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_h_bc_coef = nullptr;

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
     * Boolean to output the enthalpy in visit.
     */
    bool d_output_h = false;

    /*!
     * Energy equation parameters.
     */
    double d_liquidus_temperature, d_solidus_temperature, d_reference_temperature, d_Cp_liquid, d_Cp_solid, d_Cp_gas,
        d_Cp_mushy;

    /*!
     * Liquid fraction value in the gas. Default is set to be zero.
     */
    double d_gas_liquid_fraction = 0.0;

    /*!
     * Inner iteration parameters.
     */
    int d_max_inner_iterations = 5;
    double d_lf_iteration_error_tolerance = 1e-8;

    /*!
     * Machineries related to constant extrapolation.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_extrap_var, d_lf_extrap_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_phi_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_lf_extrap_interp_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_u_adv_fc_lf_extrap_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_u_adv_sc_lf_extrap_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_u_adv_cc_lf_extrap_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_normal_lf_extrap_var;

    int d_lf_extrap_current_idx = IBTK::invalid_index, d_lf_extrap_new_idx = IBTK::invalid_index,
        d_lf_extrap_scratch_idx = IBTK::invalid_index;
    int d_lf_extrap_interp_idx = IBTK::invalid_index;
    int d_u_adv_fc_lf_extrap_current_idx = IBTK::invalid_index;
    int d_u_adv_sc_lf_extrap_current_idx = IBTK::invalid_index;
    int d_u_adv_cc_lf_extrap_current_idx = IBTK::invalid_index;
    int d_normal_lf_extrap_current_idx = IBTK::invalid_index;
    int d_lf_extrap_rhs_scratch_idx = IBTK::invalid_index;

    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_lf_extrap_convective_op_input_db = d_default_convective_op_input_db;
    bool d_lf_extrap_convective_op_needs_init;
    SAMRAI::tbox::Pointer<CellConvectiveOperator> d_lf_extrap_convective_op = nullptr;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_lf_extrap_sol, d_lf_extrap_rhs;

    int d_lf_extrap_max_num_time_steps = 15;
    double d_lf_extrap_cell_size;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_EnthalpyHierarchyIntegrator
