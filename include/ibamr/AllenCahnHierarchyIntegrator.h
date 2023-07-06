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
                                 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
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
     * \brief Register liquid fraction variable \f$ \varphi \f$.
     */
    void registerLiquidFractionVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
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
    void setLiquidFractionPhysicalBcCoef(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                         SAMRAI::solv::RobinBcCoefStrategy<NDIM>* lf_bc_coef);

    /*!
     * Return a \f$ \varphi \f$ boundary condition object.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* getLiquidFractionPhysicalBcCoef();

    /*!
     * Get the convective operator being used by this solver class for Allen-Cahn equation.
     *
     * If the convective operator has not already been constructed, then this
     * function will initialize a default convective operator.
     */
    SAMRAI::tbox::Pointer<CellConvectiveOperator>
    getAllenCahnEquationConvectiveOperator(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > H_var);

    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

protected:
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
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * Get the solver for the Allen-Cahn equation.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver>
    getAllenCahnEquationHelmholtzSolver(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var);

    /*!
     * Get the operator to use to evaluate the right-hand side of the Allen-Cahn equation.
     */
    SAMRAI::tbox::Pointer<IBTK::LaplaceOperator>
    getAllenCahnEquationHelmholtzRHSOperator(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var);

    /*!
     * Solver variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_F_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_lf_diffusion_coef_var,
        d_lf_diffusion_coef_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_H_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_N_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_N_old_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_C_var, d_lf_temp_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_g_firstder_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_g_secondder_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_q_firstder_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_chemical_potential_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_grad_lf_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Div_u_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_lf_N_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_lf_interp_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_H_interp_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_lf_flux_var;

    /*!
     * Cartgrid functions to be used to set the Allen-Cahn equation source term.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_lf_F_fcn;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_lf_F_init;

    /*!
     * Boundary condition object for \f$ \varphi \f$.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_lf_bc_coef = nullptr;

    /*!
     * Convective operator and difference type for the Allen-Cahn equation.
     */
    ConvectiveDifferencingType d_lf_convective_difference_form = d_default_convective_difference_form;
    std::string d_lf_convective_op_type = d_default_convective_op_type;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_lf_convective_op_input_db = d_default_convective_op_input_db;
    SAMRAI::tbox::Pointer<CellConvectiveOperator> d_lf_convective_op = nullptr;
    bool d_lf_convective_op_needs_init = false;

    /*!
     * Hierarchy operators, solvers and related data for the Allen-Cahn equation.
     */
    std::string d_lf_solver_type = IBTK::CCPoissonSolverManager::UNDEFINED,
                d_lf_precond_type = IBTK::CCPoissonSolverManager::UNDEFINED;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_lf_solver_db, d_lf_precond_db;
    SAMRAI::tbox::Pointer<IBTK::LaplaceOperator> d_lf_rhs_op = nullptr;
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_lf_solver = nullptr;
    bool d_lf_solver_needs_init, d_lf_rhs_op_needs_init;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_lf_sol, d_lf_rhs;

    TimeSteppingType d_lf_diffusion_time_stepping_type = d_default_diffusion_time_stepping_type;
    TimeSteppingType d_lf_convective_time_stepping_type = d_default_convective_time_stepping_type;
    TimeSteppingType d_lf_init_convective_time_stepping_type = d_default_convective_time_stepping_type;

    /*!
     * Boolean to identify whether the energy equation is to be solved.
     */
    bool d_solve_energy = false;

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

    /*!
     * Allen-Cahn equation parameters.
     * M_lf  - Mobility
     * lambda_lf - Mixing energy density
     * eps_lf - Thickness of the liquid-solid interface
     */
    double d_M_lf, d_lambda_lf, d_eps_lf;

    /*!
     * To smoothly extend the liquid fraction to H=0 region.
     */
    double d_num_diffusion = 1.0e-8;

    /*!
     * String to identify the profile that we want to use for interpolation function q'.
     */
    std::string d_interpolation_function_profile = "LINEAR_1";
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_AllenCahnHierarchyIntegrator
