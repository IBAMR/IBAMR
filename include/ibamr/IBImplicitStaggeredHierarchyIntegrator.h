// Filename: IBImplicitStaggeredHierarchyIntegrator.h
// Created on 07 Apr 2012 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_IBImplicitStaggeredHierarchyIntegrator
#define included_IBImplicitStaggeredHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "IntVector.h"
#include "SAMRAIVectorReal.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBImplicitStrategy.h"
#include "ibamr/StaggeredStokesFACPreconditioner.h"
#include "ibamr/StaggeredStokesIBBoxRelaxationFACOperator.h"
#include "ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h"
#include "ibamr/StaggeredStokesOperator.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscpc.h"
#include "petscsnes.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Pointer.h"

namespace IBAMR
{
class INSStaggeredHierarchyIntegrator;
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{
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
 * \brief Class IBImplicitStaggeredHierarchyIntegrator is an implementation of a
 * formally second-order accurate, nonlinearly-implicit version of the immersed
 * boundary method.
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
    ~IBImplicitStaggeredHierarchyIntegrator();

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1);

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchy(double current_time, double new_time, int cycle_num = 0);

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1);

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                       SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Returns the number of cycles to perform for the present time step.
     */
    int getNumberOfCycles() const;

protected:
    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    SAMRAI::tbox::Pointer<IBImplicitStrategy> d_ib_implicit_ops;

private:
    /*!
     * \brief Partial solver for the Stokes-IB equations.
     *
     * \note This class is designed to be used internally by class IBImplicitStaggeredHierarchyIntegrator.  It is not
     * meant to be a stand-alone solver.
     */
    class IBImplicitStaggeredStokesSolver : public IBAMR::StaggeredStokesSolver
    {
    public:
        /*!
         * \brief Class constructor.
         */
        IBImplicitStaggeredStokesSolver(const std::string& object_name,
                                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
            : StaggeredStokesSolver()
        {
            d_stokes_op = new StaggeredStokesOperator(object_name + "::stokes_op", false);
            SAMRAI::tbox::Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_op =
                new StaggeredStokesIBLevelRelaxationFACOperator(object_name + "::fac_op", input_db, "stokes_ib_pc_");
            d_stokes_fac_pc =
                new StaggeredStokesFACPreconditioner(object_name + "::fac_pc", fac_op, input_db, "stokes_ib_pc_");
            return;
        }

        /*!
         * \brief Class desctructor.
         */
        ~IBImplicitStaggeredStokesSolver()
        {
            // intentionally left blank
            return;
        }

        // \{ Implementation of IBAMR::StaggeredStokesSolver class.

        void setVelocityPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& U_problem_coefs)
        {
            StaggeredStokesSolver::setVelocityPoissonSpecifications(U_problem_coefs);
            d_stokes_op->setVelocityPoissonSpecifications(U_problem_coefs);
            d_stokes_fac_pc->setVelocityPoissonSpecifications(U_problem_coefs);
            return;
        }

        void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
                                SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef)
        {
            StaggeredStokesSolver::setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
            d_stokes_op->setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
            d_stokes_fac_pc->setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
            return;
        }

        void setPhysicalBoundaryHelper(SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
        {
            StaggeredStokesSolver::setPhysicalBoundaryHelper(bc_helper);
            d_stokes_op->setPhysicalBoundaryHelper(bc_helper);
            d_stokes_fac_pc->setPhysicalBoundaryHelper(bc_helper);
            return;
        }

        void setComponentsHaveNullspace(const bool has_velocity_nullspace, const bool has_pressure_nullspace)
        {
            StaggeredStokesSolver::setComponentsHaveNullspace(has_velocity_nullspace, has_pressure_nullspace);
            d_stokes_fac_pc->setComponentsHaveNullspace(d_has_velocity_nullspace, d_has_pressure_nullspace);
            return;
        }

        // \}

        // \{ Implementation of IBTK::GeneralSolver class.

        bool solveSystem(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& /*x*/,
                         SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& /*b*/)
        {
            TBOX_ERROR("StaggeredStokesIBSolver::solveSystem(): unimplemented.\n");
            return false;
        }

        // \}

        SAMRAI::tbox::Pointer<StaggeredStokesOperator> getStaggeredStokesOperator()
        {
            return d_stokes_op;
        }

        SAMRAI::tbox::Pointer<StaggeredStokesFACPreconditioner> getStaggeredStokesFACPreconditioner()
        {
            return d_stokes_fac_pc;
        }

    private:
        IBImplicitStaggeredStokesSolver();
        IBImplicitStaggeredStokesSolver(const IBImplicitStaggeredStokesSolver& from);
        IBImplicitStaggeredStokesSolver& operator=(const IBImplicitStaggeredStokesSolver& that);

        // Operators and solvers maintained by this class.
        SAMRAI::tbox::Pointer<StaggeredStokesOperator> d_stokes_op;
        SAMRAI::tbox::Pointer<StaggeredStokesFACPreconditioner> d_stokes_fac_pc;
    };

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBImplicitStaggeredHierarchyIntegrator(const IBImplicitStaggeredHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBImplicitStaggeredHierarchyIntegrator& operator=(const IBImplicitStaggeredHierarchyIntegrator& that);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * \brief Solve for position along with fluid variables.
     */
    void integrateHierarchy_position(double current_time, double new_time, int cycle_num);

    /*!
     * \brief Solve for fluid variables only.
     */
    void integrateHierarchy_velocity(double current_time, double new_time, int cycle_num);

    /*!
     * Static function for implicit formulation.
     */
    static PetscErrorCode IBFunction_SAMRAI(SNES snes, Vec x, Vec f, void* ctx);

    /*!
     * Function for implicit formulation that solves for u,p and X.
     */
    PetscErrorCode IBFunction_position(SNES snes, Vec x, Vec f);

    /*!
     * Function for implicit formulation that solves for u and p.
     */
    PetscErrorCode IBFunction_velocity(SNES snes, Vec x, Vec f);

    /*!
     * Static function for setting up implicit formulation Jacobian.
     */
    static PetscErrorCode IBJacobianSetup_SAMRAI(SNES snes, Vec x, Mat A, Mat B, void* p_ctx);

    /*!
     * Static function for setting up implicit formulation Jacobian that solves for u, p, and X.
     */
    PetscErrorCode IBJacobianSetup_position(SNES snes, Vec x, Mat A, Mat B);

    /*!
     * Static function for setting up implicit formulation Jacobian that solves for u and p.
     */
    PetscErrorCode IBJacobianSetup_velocity(SNES snes, Vec x, Mat A, Mat B);

    /*!
     * Static function for implicit formulation Jacobian.
     */
    static PetscErrorCode IBJacobianApply_SAMRAI(Mat A, Vec x, Vec y);

    /*!
     * Function for implicit formulation Jacobian that solves for u, p, and X.
     */
    PetscErrorCode IBJacobianApply_position(Vec x, Vec y);

    /*
     * Function for implicit formulation Jacobian that solves for u and p.
     */
    PetscErrorCode IBJacobianApply_velocity(Vec x, Vec y);

    /*!
     * Static function for implicit formulation preconditioner.
     */
    static PetscErrorCode IBPCApply_SAMRAI(PC pc, Vec x, Vec y);

    /*!
     * Function for implicit formulation preconditioner that solves for u, p, and X.
     */
    PetscErrorCode IBPCApply_position(Vec x, Vec y);

    /*!
     * Function for implicit formulation preconditioner that solves for u and p.
     */
    PetscErrorCode IBPCApply_velocity(Vec x, Vec y);

    /*!
     * Static function for implicit formulation Lagrangian Schur complement.
     */
    static PetscErrorCode lagrangianSchurApply_SAMRAI(Mat A, Vec x, Vec y);

    /*!
     * Function for implicit formulation Lagrangian Schur complement.
     */
    PetscErrorCode lagrangianSchurApply(Vec x, Vec y);

    // Eulerian data for storing u and p DOFs indexing.
    std::vector<std::vector<int> > d_num_dofs_per_proc;
    int d_u_dof_index_idx, d_p_dof_index_idx;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, int> > d_u_dof_index_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, int> > d_p_dof_index_var;

    // Solvers and associated vectors.
    bool d_solve_for_position;
    std::string d_jac_delta_fcn;
    SAMRAI::tbox::Pointer<StaggeredStokesSolver> d_stokes_solver;
    SAMRAI::tbox::Pointer<StaggeredStokesOperator> d_stokes_op;
    KSP d_schur_solver;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_u_scratch_vec, d_f_scratch_vec;
    Vec d_X_current;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBImplicitStaggeredHierarchyIntegrator
