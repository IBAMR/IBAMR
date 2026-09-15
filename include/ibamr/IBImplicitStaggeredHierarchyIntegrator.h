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

#include <ibtk/IBKernelConcepts.h>
#include <ibtk/IBKernelTensorProduct.h>
#include <ibtk/PETScNewtonKrylovSolver.h>

#include <tbox/Pointer.h>

#include <IntVector.h>
#include <SAMRAIVectorReal.h>

#include <concepts>
#include <functional>
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
 * hierarchy Jacobian applies the strategy's force linearization; FAC uses
 * separately assembled coupling.
 *
 * Configure the nonlinear solver in this object's input database (PETSc prefix
 * \c ib_) and FAC in \c stokes_ib_precond_db (prefix \c stokes_ib_pc_).
 * \c jacobian_delta_fcn selects an IBKernelTensorProduct using IB_3 through IB_6
 * or B-splines within the compiled input-selection bound documented in [CMake
 * configuration](../../doc/cmake.md#implicit-ib-interpolation-kernels). The strategy's minimum ghost width must cover
 * that kernel; with IBMethod, set \c min_ghost_cell_width when needed. Names are parsed at construction; built-in
 * availability is checked at initialization unless an explicit evaluator was supplied with
 * setJacobianInterpolationKernel().
 *
 * Fixed coupling is enabled on the supplied strategy at construction. Subclasses
 * overriding time-step or hierarchy hooks must call the corresponding base
 * implementation to prepare and release solver state. See
 * StaggeredStokesIBOperator::Context for the shared operator requirements.
 */
class IBImplicitStaggeredHierarchyIntegrator : public IBHierarchyIntegrator
{
public:
    /*! \brief Own a concrete interpolation evaluator configured before initialization.
     *
     * The evaluator is moved into owned storage and evaluated as const according
     * to IBTK::IBKernelEvaluatorCartesian and the stencil/weight conventions of
     * IBTK::PETScMatUtilities::constructPatchLevelSCInterpOp(). It must own any
     * state needed during later advances; move-only evaluators are supported.
     * This overrides any valid input kernel name, independently of the compiled
     * bound. Invalid names still fail at construction; calls after initialization
     * fail. This selects assembled FAC coupling, not live interpolation/spreading.
     * See [the configuration example](../../doc/cmake.md#implicit-ib-interpolation-kernels).
     */
    template <IBTK::IBKernelEvaluatorCartesian<double, PetscScalar> Evaluator>
    void setJacobianInterpolationKernel(Evaluator evaluator) requires(std::constructible_from<Evaluator, Evaluator&&>);

    /*! \brief Construct an implicit velocity-pressure integrator and enable fixed coupling. */
    IBImplicitStaggeredHierarchyIntegrator(const std::string& object_name,
                                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                           SAMRAI::tbox::Pointer<IBImplicitStrategy> ib_method_ops,
                                           SAMRAI::tbox::Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator,
                                           bool register_for_restart = true);

    /*! \brief Release the implicit solver state. */
    ~IBImplicitStaggeredHierarchyIntegrator() override;

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

    /*! \brief Resolve the interpolation evaluator and initialize coupled solver data.
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
    using InterpolationMatrixBuilder = std::function<
        void(Mat&, Vec, const std::vector<int>&, int, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>)>;

    /*! \brief Return a whole-matrix builder owning a concrete const evaluator. */
    template <IBTK::IBKernelEvaluatorCartesian<double, PetscScalar> Evaluator>
    static InterpolationMatrixBuilder make_matrix_builder(Evaluator&& evaluator)
        requires(std::constructible_from<Evaluator, Evaluator&&>);

    /*! \brief Select a compiled evaluator for the configured normal/tangential factors. */
    static InterpolationMatrixBuilder select_matrix_builder(const IBTK::IBKernelTensorProduct& kernel);

    /*! \brief Copy construction is disabled. */
    IBImplicitStaggeredHierarchyIntegrator(const IBImplicitStaggeredHierarchyIntegrator& from) = delete;
    /*! \brief Copy assignment is disabled. */
    IBImplicitStaggeredHierarchyIntegrator& operator=(const IBImplicitStaggeredHierarchyIntegrator& that) = delete;

    /*!
     * Read object state from the restart file and initialize class data
     * members.
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
    //! Input descriptor resolved at initialization unless an explicit evaluator is supplied.
    IBTK::IBKernelTensorProduct d_jac_kernel = IBTK::IBKernel::IB_4;
    //! Owns the operation assembling the FAC interpolation matrix.
    InterpolationMatrixBuilder d_interp_matrix_builder;
    bool d_vectors_need_init = true;
    bool d_has_velocity_nullspace = false;
    bool d_has_pressure_nullspace = true;
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
    //! Owned matrices borrowed by FAC until solver state is deallocated.
    Mat d_ib_force_jac = nullptr;
    Mat d_ib_interp_op = nullptr;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#include <ibamr/private/IBImplicitStaggeredHierarchyIntegrator-inl.h>

#endif // #ifndef included_IBAMR_IBImplicitStaggeredHierarchyIntegrator
