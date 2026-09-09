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
#include <ibtk/PETScNewtonKrylovSolver.h>

#include <tbox/Pointer.h>

#include <IntVector.h>
#include <SAMRAIVectorReal.h>

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
 * supplies fluid time-stepping terms and boundary conditions. The hierarchy
 * Jacobian applies the analytic IB linearization; FAC uses assembled coupling.
 *
 * Configure the nonlinear solver in this object's input database (PETSc prefix
 * \c ib_) and FAC in \c stokes_ib_precond_db (prefix \c stokes_ib_pc_).
 * \c jacobian_delta_fcn selects an IBTK::IBKernelTensorProduct registered for
 * interpolation-matrix construction. The strategy's minimum ghost width must
 * cover that kernel; with IBMethod, set \c min_ghost_cell_width when needed.
 * Application evaluators must be registered before hierarchy advancement.
 *
 * Fixed coupling is enabled on the supplied strategy at construction. Subclasses
 * overriding time-step or hierarchy hooks must call the corresponding base
 * implementation to prepare and release solver state. See
 * StaggeredStokesIBOperator::Context for the shared operator requirements.
 */
class IBImplicitStaggeredHierarchyIntegrator : public IBHierarchyIntegrator
{
public:
    /*!
     * The constructor for class IBImplicitStaggeredHierarchyIntegrator sets
     * some default values, reads in configuration information from input and
     * restart databases, and registers the integrator object with the restart
     * manager when requested.
     */
    IBImplicitStaggeredHierarchyIntegrator(const std::string& object_name,
                                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                           SAMRAI::tbox::Pointer<IBImplicitStrategy> ib_method_ops,
                                           SAMRAI::tbox::Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator,
                                           bool register_for_restart = true);

    /*!
     * The destructor for class IBImplicitStaggeredHierarchyIntegrator
     * unregisters the integrator object with the restart manager when the
     * object is so registered.
     */
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

    //! Strategy advanced by the coupled solve.
    SAMRAI::tbox::Pointer<IBImplicitStrategy> d_ib_implicit_ops;

private:
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
    //! Kernel used to assemble the FAC interpolation matrix.
    IBTK::IBKernelTensorProduct d_jac_kernel = IBTK::IBKernel::IB_4;
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
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_eul_sol_vec;
    //! Owns cloned components whose descriptors persist between advances.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_eul_rhs_vec;
    //! Owned matrices borrowed by FAC until solver state is deallocated.
    Mat d_ib_force_jac = nullptr;
    Mat d_ib_interp_op = nullptr;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_IBImplicitStaggeredHierarchyIntegrator
