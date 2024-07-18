// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2024 by the IBAMR developers
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

#ifndef included_IBAMR_AllenCahnHierarchyIntegrator
#define included_IBAMR_AllenCahnHierarchyIntegrator

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
 * \brief Class AllenCahnHierarchyIntegrator is a concrete class that manages the
 * time integration of the Allen-Cahn and energy equation.
 *
 * In this class, implementations are provided for phase change of a pure material
 * using the phase-field method.
 *
 */
class AllenCahnHierarchyIntegrator : public PhaseChangeHierarchyIntegrator
{
public:
    /*!
     * The constructor for class AllenCahnHierarchyIntegrator sets
     * some default values, reads in configuration information from input and
     * restart databases, and registers the integrator object with the restart
     * manager when requested.
     */
    AllenCahnHierarchyIntegrator(const std::string& object_name,
                                 IBTK::SAMRAIPointer<SAMRAI::tbox::Database> input_db,
                                 bool register_for_restart = true);

    /*!
     * The destructor for class AllenCahnHierarchyIntegrator
     * unregisters the integrator object with the restart manager when the
     * object is so registered.
     */
    ~AllenCahnHierarchyIntegrator() = default;

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void initializeHierarchyIntegrator(IBTK::SAMRAIPointer<SAMRAI::hier::PatchHierarchyNd> hierarchy,
                                       IBTK::SAMRAIPointer<SAMRAI::mesh::GriddingAlgorithmNd> gridding_alg) override;

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1) override;

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override;

    /*!
     * \brief Register liquid fraction variable \f$ \varphi \f$.
     */
    void registerLiquidFractionVariable(IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > lf_var,
                                        const bool output_lf_var = true) override;

    /*!
     *  Add the temporal terms to the RHS of the energy equation.
     */
    void addTemporalAndLinearTermstoRHSOfEnergyEquation(int F_scratch_idx, const double dt) override;

    /*!
     * Compute the source term for the Div U equation.
     */
    void computeDivergenceVelocitySourceTerm(int Div_U_F_idx, const double new_time) override;

    /*!
     * Set an object to provide boundary conditions for  \f$ \varphi \f$ variable,
     * that has been registered with the hierarchy integrator.
     */
    void setLiquidFractionPhysicalBcCoef(IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > lf_var,
                                         SAMRAI::solv::RobinBcCoefStrategyNd* lf_bc_coef);

    /*!
     * Return a \f$ \varphi \f$ boundary condition object.
     */
    SAMRAI::solv::RobinBcCoefStrategyNd* getLiquidFractionPhysicalBcCoef();

    /*!
     * Get the convective operator being used by this solver class for Allen-Cahn equation.
     *
     * If the convective operator has not already been constructed, then this
     * function will initialize a default convective operator.
     */
    IBTK::SAMRAIPointer<CellConvectiveOperator>
    getAllenCahnEquationConvectiveOperator(IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > lf_var,
                                           IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > H_var);

    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(IBTK::SAMRAIPointer<SAMRAI::tbox::Database> db) override;

protected:
    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchySpecialized(double current_time, double new_time, int cycle_num = 0) override;

    /*!
     * Reset cached hierarchy dependent data for solvers and operators before the regridding operation.
     */
    virtual void regridHierarchyBeginSpecialized() override;

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
    AllenCahnHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AllenCahnHierarchyIntegrator(const AllenCahnHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AllenCahnHierarchyIntegrator& operator=(const AllenCahnHierarchyIntegrator& that) = delete;

    /*!
     * \brief Compute the double-well potential on the cell-centers based on the
     * liquid fraction value.
     */
    void computeDoubleWellPotential(int g_firstder_idx, int g_secondder_idx, const int liquid_fraction_idx);

    /*!
     * \brief Compute an interpolation function on the cell-centers based on the
     * liquid fraction value.
     */
    void computeInterpolationFunction(int p_firstder_idx, const int liquid_fraction_idx, const int T_idx);

    /*
     * \brief Compute the source term of liquid fraction equation.
     */
    void computeLiquidFractionSourceTerm(int F_scratch_idx);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(IBTK::SAMRAIPointer<SAMRAI::tbox::Database> input_db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * Get the solver for the Allen-Cahn equation.
     */
    IBTK::SAMRAIPointer<IBTK::PoissonSolver>
    getAllenCahnEquationHelmholtzSolver(IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > lf_var);

    /*!
     * Get the operator to use to evaluate the right-hand side of the Allen-Cahn equation.
     */
    IBTK::SAMRAIPointer<IBTK::LaplaceOperator>
    getAllenCahnEquationHelmholtzRHSOperator(IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > lf_var);

    /*!
     * Solver variables.
     */
    IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > d_lf_F_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::SideVariableNd<double> > d_lf_diffusion_coef_var, d_lf_diffusion_coef_rhs_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > d_lf_H_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > d_lf_N_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > d_lf_N_old_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > d_lf_rhs_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > d_lf_C_var, d_lf_temp_rhs_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > d_g_firstder_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > d_g_secondder_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > d_q_firstder_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > d_chemical_potential_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::SideVariableNd<double> > d_grad_lf_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > d_Div_u_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > d_T_lf_N_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::FaceVariableNd<double> > d_lf_interp_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::FaceVariableNd<double> > d_H_interp_var;
    IBTK::SAMRAIPointer<SAMRAI::pdat::FaceVariableNd<double> > d_lf_flux_var;

    /*!
     * Cartgrid functions to be used to set the Allen-Cahn equation source term.
     */
    IBTK::SAMRAIPointer<IBTK::CartGridFunction> d_lf_F_fcn;
    IBTK::SAMRAIPointer<IBTK::CartGridFunction> d_lf_F_init;

    /*!
     * Boundary condition object for \f$ \varphi \f$.
     */
    SAMRAI::solv::RobinBcCoefStrategyNd* d_lf_bc_coef = nullptr;

    /*!
     * Convective operator and difference type for the Allen-Cahn equation.
     */
    ConvectiveDifferencingType d_lf_convective_difference_form = d_default_convective_difference_form;
    std::string d_lf_convective_op_type = d_default_convective_op_type;
    IBTK::SAMRAIPointer<SAMRAI::tbox::Database> d_lf_convective_op_input_db = d_default_convective_op_input_db;
    IBTK::SAMRAIPointer<CellConvectiveOperator> d_lf_convective_op = nullptr;
    bool d_lf_convective_op_needs_init = false;

    /*!
     * Hierarchy operators, solvers and related data for the Allen-Cahn equation.
     */
    std::string d_lf_solver_type = IBTK::CCPoissonSolverManager::UNDEFINED,
                d_lf_precond_type = IBTK::CCPoissonSolverManager::UNDEFINED;
    IBTK::SAMRAIPointer<SAMRAI::tbox::Database> d_lf_solver_db, d_lf_precond_db;
    IBTK::SAMRAIPointer<IBTK::LaplaceOperator> d_lf_rhs_op = nullptr;
    IBTK::SAMRAIPointer<IBTK::PoissonSolver> d_lf_solver = nullptr;
    bool d_lf_solver_needs_init, d_lf_rhs_op_needs_init;
    IBTK::SAMRAIPointer<SAMRAI::solv::SAMRAIVectorRealNd<double> > d_lf_sol, d_lf_rhs;

    TimeSteppingType d_lf_diffusion_time_stepping_type = d_default_diffusion_time_stepping_type;
    TimeSteppingType d_lf_convective_time_stepping_type = d_default_convective_time_stepping_type;
    TimeSteppingType d_lf_init_convective_time_stepping_type = d_default_convective_time_stepping_type;

    /*!
     * Boolean to identify whether the energy equation is to be solved.
     */
    bool d_solve_energy = false;

    /*!
     * Patch data descriptor indices for all "state" variables managed by the
     * integrator.
     *
     * State variables have three contexts: current, scratch, and new.
     */
    int d_lf_F_current_idx = IBTK::invalid_index, d_lf_F_scratch_idx = IBTK::invalid_index,
        d_lf_F_new_idx = IBTK::invalid_index;
    int d_lf_diff_coef_current_idx = IBTK::invalid_index, d_lf_diff_coef_scratch_idx = IBTK::invalid_index,
        d_lf_diff_coef_new_idx = IBTK::invalid_index;
    int d_lf_N_old_current_idx = IBTK::invalid_index, d_lf_N_old_new_idx = IBTK::invalid_index,
        d_lf_N_old_scratch_idx = IBTK::invalid_index;

    /*!
     * Patch data descriptor indices for all "scratch" variables managed by the
     * integrator.
     *
     * Scratch variables have only one context: scratch.
     */
    int d_lf_temp_rhs_idx = IBTK::invalid_index, d_lf_C_idx = IBTK::invalid_index;
    int d_g_firstder_idx = IBTK::invalid_index, d_g_secondder_idx = IBTK::invalid_index,
        d_q_firstder_idx = IBTK::invalid_index;
    int d_chemical_potential_idx = IBTK::invalid_index, d_grad_lf_idx = IBTK::invalid_index;
    int d_H_sc_idx = IBTK::invalid_index;
    int d_T_lf_N_scratch_idx = IBTK::invalid_index;
    int d_lf_interp_idx = IBTK::invalid_index, d_H_interp_idx = IBTK::invalid_index,
        d_lf_flux_idx = IBTK::invalid_index;
    int d_lf_diffusion_coef_rhs_scratch_idx = IBTK::invalid_index;
    int d_lf_rhs_scratch_idx = IBTK::invalid_index;
    int d_lf_H_scratch_idx = IBTK::invalid_index;
    int d_Div_u_scratch_idx = IBTK::invalid_index;
    int d_lf_N_scratch_idx = IBTK::invalid_index;

    /*!
     * Allen-Cahn equation parameters.
     */
    double d_mobility_lf, d_mixing_energy_density_lf, d_interface_thickness_lf;

    /*!
     * To smoothly extend the liquid fraction to H=0 region.
     */
    double d_numerical_diffusion = 1.0e-8;

    /*!
     * String to identify the profile that we want to use for interpolation function q'.
     */
    std::string d_interpolation_function_profile = "LINEAR_1";
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_AllenCahnHierarchyIntegrator
