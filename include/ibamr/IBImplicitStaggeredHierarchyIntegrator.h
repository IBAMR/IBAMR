// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#ifndef included_IBAMR_IBImplicitStaggeredHierarchyIntegrator
#define included_IBAMR_IBImplicitStaggeredHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/IBImplicitStrategy.h>
#include <ibamr/StaggeredStokesIBJacobianFACPreconditioner.h>
#include <ibamr/StaggeredStokesIBJacobianOperator.h>
#include <ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesIBOperator.h>
#include <ibamr/StaggeredStokesOperator.h>

#include <ibtk/IBKernelTensorProduct.h>
#include <ibtk/IBOperatorBuilder.h>
#include <ibtk/PETScNewtonKrylovSolver.h>

#include <tbox/Pointer.h>

#include <IntVector.h>
#include <SAMRAIVectorReal.h>

#include <map>
#include <optional>
#include <string>

namespace IBAMR
{
class INSStaggeredHierarchyIntegrator;
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
namespace solv
{
class PoissonSpecifications;
} // namespace solv
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Implicit immersed boundary time integration with shared Stokes-IB operators.
 *
 * Uses the velocity-pressure formulation of StaggeredStokesIBOperator with
 * backward Euler, trapezoidal, or midpoint IB stepping. The INS integrator
 * supplies fluid time-stepping terms and boundary conditions. The matrix-free
 * hierarchy Jacobian applies the strategy's force linearization and follows the
 * Newton iterate; FAC uses separately assembled coupling, which is built once
 * per cycle at the structure configuration of the time-stepping rule's
 * evaluation time and is not refreshed between Newton iterations.
 *
 * Configure the nonlinear solver in this object's input database (PETSc prefix
 * \c ib_) and FAC in \c stokes_ib_precond_db (prefix \c stokes_ib_pc_).
 * \c jacobian_delta_fcn (default IB_4) names the delta function of the
 * interpolation matrix in the assembled coupling that FAC uses. It does not
 * affect the matrix-free Jacobian, whose action comes from the strategy, and it
 * need not equal the kernel that the strategy uses for interpolation and
 * spreading. A name selects a built-in kernel only if
 * IBTK::IBOperatorBuilder::is_built_in() accepts it (see its documentation for
 * the built-in set). registerJacobianOperatorBuilder() or setJacobianOperatorBuilder()
 * selects any other kernel.
 * The velocity DOF index
 * data has ghost width equal to the larger of the strategy's minimum ghost
 * width and the kernel's; the pressure DOF index data has no ghost cells.
 *
 * Fixed coupling is enabled on the supplied strategy at construction. Subclasses
 * overriding time-step or hierarchy hooks must call the corresponding base
 * implementation to prepare and release solver state. See
 * StaggeredStokesIBOperator::Context for the shared operator requirements.
 */
class IBImplicitStaggeredHierarchyIntegrator : public IBHierarchyIntegrator
{
public:
    /*! \brief Construct an implicit velocity-pressure integrator and enable fixed coupling. */
    IBImplicitStaggeredHierarchyIntegrator(const std::string& object_name,
                                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                           SAMRAI::tbox::Pointer<IBImplicitStrategy> ib_method_ops,
                                           SAMRAI::tbox::Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator,
                                           bool register_for_restart = true);

    /*! \brief Release the implicit solver state. */
    ~IBImplicitStaggeredHierarchyIntegrator() override;

    /*!
     * Use builder to construct the interpolation matrix in the assembled
     * coupling that FAC uses. This replaces the builder selected by the
     * jacobian_delta_fcn input, and must be called before
     * initializeHierarchyIntegrator().
     *
     * This kernel is independent of the kernels that the IB strategy uses for
     * interpolation and spreading. The builder is not written to restart files:
     * after a restart, call this function again, since the input key
     * jacobian_delta_fcn cannot name it.
     */
    void setJacobianOperatorBuilder(IBTK::IBOperatorBuilder builder);

    /*!
     * Associate builder with kernel, so that jacobian_delta_fcn = kernel's name
     * selects it. Must be called before initializeHierarchyIntegrator(). It is
     * an error to register a kernel that IBTK::IBOperatorBuilder::is_built_in()
     * already accepts, or to register the same kernel twice.
     *
     * A registration only supplies the builder for jacobian_delta_fcn if that
     * key is otherwise unresolved when initializeHierarchyIntegrator() runs;
     * an explicit call to setJacobianOperatorBuilder() always takes precedence,
     * whether it happens before or after this call. Registering a kernel that
     * jacobian_delta_fcn does not name has no other effect.
     *
     * Like setJacobianOperatorBuilder(), registrations are not written to
     * restart files: after a restart, the application must register again
     * before initializeHierarchyIntegrator().
     */
    void registerJacobianOperatorBuilder(const IBTK::IBKernelTensorProduct& kernel, IBTK::IBOperatorBuilder builder);

    /*!
     * Return the builder of the interpolation matrix in the assembled FAC
     * coupling. It is an error if jacobian_delta_fcn is not a built-in kernel
     * and no builder has been set or registered for it.
     */
    const IBTK::IBOperatorBuilder& getJacobianOperatorBuilder() const;

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

    /*! \brief Initialize coupled solver data, which requires an interpolation builder.
     *
     * The builder is selected from jacobian_delta_fcn at construction or by
     * setJacobianOperatorBuilder(); this method does not select one.
     * \see IBHierarchyIntegrator::initializeHierarchyIntegrator()
     */
    void
    initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM>> gridding_alg) override;

    /*!
     * Returns the number of cycles to perform for the present time step.
     */
    int getNumberOfCycles() const override;

protected:
    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchySpecialized(double current_time, double new_time, int cycle_num = 0) override;

    /*!
     * Reset cached hierarchy dependent data.
     */
    void resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
                                                int coarsest_level,
                                                int finest_level) override;

    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    //! Typed view of d_ib_method_ops for implicit strategy operations.
    SAMRAI::tbox::Pointer<IBImplicitStrategy> d_ib_implicit_ops;

private:
    /*! \brief Copy construction is disabled. */
    IBImplicitStaggeredHierarchyIntegrator(const IBImplicitStaggeredHierarchyIntegrator& from) = delete;
    /*! \brief Copy assignment is disabled. */
    IBImplicitStaggeredHierarchyIntegrator& operator=(const IBImplicitStaggeredHierarchyIntegrator& that) = delete;

    /*!
     * Check that the restart database exists and has this class's version.
     */
    void getFromRestart();

    /*!
     * Setup and allocate Eulerian solver vectors used in implicit solves.
     */
    void setupSolverVectors(double current_time, int coarsest_ln, int finest_ln);

    /*!
     * Setup hierarchy dependent operators and solvers used in implicit solves.
     */
    void reinitializeOperatorsAndSolvers(double current_time, double new_time);

    /*!
     * Deallocate hierarchy dependent operators and solvers used in implicit solves.
     */
    void deallocateOperatorsAndSolvers();

    // Eulerian data for storing u and p DOFs indexing.
    std::vector<int> d_num_dofs_per_proc;
    int d_u_dof_index_idx = IBTK::invalid_index, d_p_dof_index_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, int>> d_u_dof_index_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, int>> d_p_dof_index_var;

    // Solvers and associated vectors.
    std::string d_jac_delta_fcn = "IB_4";
    //! Builds the interpolation matrix used in the assembled FAC coupling.
    std::optional<IBTK::IBOperatorBuilder> d_jacobian_operator_builder;
    std::map<IBTK::IBKernelTensorProduct, IBTK::IBOperatorBuilder> d_registered_jacobian_operator_builders;
    bool d_vectors_need_init = true;
    bool d_has_velocity_nullspace = false;
    bool d_has_pressure_nullspace = true;
    //! Shared by the IB operators, which apply it; this object initializes and deallocates its state.
    SAMRAI::tbox::Pointer<StaggeredStokesOperator> d_stokes_op;
    SAMRAI::tbox::Pointer<StaggeredStokesIBOperator> d_ib_op;
    SAMRAI::tbox::Pointer<StaggeredStokesIBJacobianOperator> d_ib_jac_op;
    SAMRAI::tbox::Pointer<StaggeredStokesIBJacobianFACPreconditioner> d_ib_jac_pc;
    SAMRAI::tbox::Pointer<IBTK::PETScNewtonKrylovSolver> d_ib_solver;
    SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> d_stokes_bc_helper;
    //! References INS-owned velocity/pressure scratch components.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_sol_vec;
    //! Owns cloned components whose descriptors persist between advances.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_rhs_vec;
    //! Owned matrices, rebuilt on each cycle; FAC retains its own references.
    Mat d_ib_force_jac = nullptr;
    Mat d_ib_interp_op = nullptr;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_IBImplicitStaggeredHierarchyIntegrator
