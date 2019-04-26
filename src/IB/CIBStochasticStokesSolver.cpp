// Filename: CIBStaggeredStokesSolver.cpp
// Created on 26 Sept 2017 by Brennan Sprinkle
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of its
//      contributors may be used to endorse or promote products derived from
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/CIBStochasticStokesSolver.h"
#include "ibamr/CIBSaddlePointSolver.h"
#include "ibamr/CIBStochasticMethod.h"
#include "ibamr/CIBStrategy.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/RNG.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/LinearSolver.h"
#include "ibtk/NewtonKrylovSolver.h"
#include "ibtk/PETScKrylovLinearSolver.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "tbox/Database.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
namespace
{
// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;
} // namespace

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CIBStochasticStokesSolver::CIBStochasticStokesSolver(const std::string& object_name,
                                                     Pointer<Database> input_db,
                                                     Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
                                                     Pointer<CIBStrategy> cib_strategy,
                                                     const std::string& default_options_prefix)
    : StaggeredStokesSolver(),
      d_cib_strategy(cib_strategy, false),
      d_num_rigid_parts(d_cib_strategy->getNumberOfRigidStructures())
{
    GeneralSolver::init(object_name, /*homogeneous bcs*/ false);

    d_sp_solver = NULL;
    d_wide_u_var = NULL;
    d_wide_f_var = NULL;
    d_wide_ctx = NULL;
    d_wide_u_idx = -1;
    d_wide_f_idx = -1;
    d_x_wide = NULL;
    d_b_wide = NULL;
    d_is_initialized = false;
    d_reinitializing_solver = false;

    // Create the saddle-point solver for solving constraint problem.
    d_sp_solver = new CIBSaddlePointSolver(
        object_name, input_db, navier_stokes_integrator, d_cib_strategy, default_options_prefix);

    // Create widened variables for IB operations.
    Pointer<IBStrategy> ib_method_ops = d_cib_strategy;
    const IntVector<NDIM> ghost_width = ib_method_ops->getMinimumGhostCellWidth();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_wide_u_var = new SideVariable<NDIM, double>(d_object_name + "::wide_u_var", 1);
    d_wide_f_var = new SideVariable<NDIM, double>(d_object_name + "::wide_f_var", 1);
    d_wide_ctx = var_db->getContext(object_name + "::wide_ctx");
    d_wide_u_idx = var_db->registerVariableAndContext(d_wide_u_var, d_wide_ctx, ghost_width);
    d_wide_f_idx = var_db->registerVariableAndContext(d_wide_f_var, d_wide_ctx, ghost_width);

    return;
} // CIBStochasicStokesSolver

CIBStochasticStokesSolver::~CIBStochasticStokesSolver()
{
    d_reinitializing_solver = false;
    deallocateSolverState();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    var_db->removePatchDataIndex(d_wide_u_idx);
    var_db->removePatchDataIndex(d_wide_f_idx);

    return;
} // ~CIBStochasicStokesSolver()

void
CIBStochasticStokesSolver::setSolutionTime(double solution_time)
{
    GeneralSolver::setSolutionTime(solution_time);
    d_sp_solver->setSolutionTime(solution_time);

    return;
} // setSolutionTime

void
CIBStochasticStokesSolver::setTimeInterval(double current_time, double new_time)
{
    GeneralSolver::setTimeInterval(current_time, new_time);
    d_sp_solver->setTimeInterval(current_time, new_time);

    return;
} // setTimeInterval

Pointer<CIBSaddlePointSolver>
CIBStochasticStokesSolver::getSaddlePointSolver() const
{
    return d_sp_solver;
} // getSaddlePointSolver

void
CIBStochasticStokesSolver::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                                                 const SAMRAIVectorReal<NDIM, double>& b)
{
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized)
    {
        d_reinitializing_solver = true;
        deallocateSolverState();
    }

    // Wrap Eulerian data into PETSc Vecs.
    Pointer<PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
    const int coarsest_ln = x.getCoarsestLevelNumber();
    const int finest_ln = x.getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_wide_u_idx)) level->allocatePatchData(d_wide_u_idx);
        if (!level->checkAllocated(d_wide_f_idx)) level->allocatePatchData(d_wide_f_idx);
    }

    Pointer<CellVariable<NDIM, double> > x_p_cc_var = x.getComponentVariable(1);
    Pointer<CellVariable<NDIM, double> > b_p_cc_var = b.getComponentVariable(1);
    const int x_p_idx = x.getComponentDescriptorIndex(1);
    const int b_p_idx = b.getComponentDescriptorIndex(1);

    d_x_wide = new SAMRAIVectorReal<NDIM, double>(x.getName() + "_wide_x", hierarchy, coarsest_ln, finest_ln);
    d_b_wide = new SAMRAIVectorReal<NDIM, double>(b.getName() + "_wide_b", hierarchy, coarsest_ln, finest_ln);

    d_x_wide->addComponent(d_wide_u_var, d_wide_u_idx, x.getControlVolumeIndex(0));
    d_x_wide->addComponent(x_p_cc_var, x_p_idx, x.getControlVolumeIndex(1));
    d_b_wide->addComponent(d_wide_f_var, d_wide_f_idx, b.getControlVolumeIndex(0));
    d_b_wide->addComponent(b_p_cc_var, b_p_idx, b.getControlVolumeIndex(1));

    // Wrap SAMRAI vector into PETSc Vec
    Vec u_p = PETScSAMRAIVectorReal::createPETScVector(d_x_wide);
    Vec g_h = PETScSAMRAIVectorReal::createPETScVector(d_b_wide);

    // Get the Lagrange multiplier that maintains the rigidity constraint.
    // NOTE: The current time corresponds to the time at which solver is initialized
    // which maybe different from the current time of the timestep being integrated upon.
    Vec L;
    d_cib_strategy->getConstraintForce(&L, d_current_time);

    // Create a vector of the type imposed velocity at the material/nodal points for RHS.
    // NOTE: In the initialization stage we do not need an actual velocity vector, a reference
    // to L should suffice to know the required structure.
    Vec V = L;

    // Get the rigid body velocities that need to be solved for.
    // NOTE: The current time corresponds to the time at which solver is initialized
    // which maybe different from the current time of the timestep being integrated upon.
    Vec U;
    d_cib_strategy->getFreeRigidVelocities(&U, d_current_time);

    // Create a vector that contains net external force and torque on the body.
    // NOTE: In the initialization stage we do not need an actual force vector, a reference
    // to U should suffice to know the required structure.
    Vec F = U;

    // Create the composite vectors.
    Vec mv_x, mv_b;
    std::array<Vec,3> vx, vb;
    vx[0] = u_p;
    vx[1] = L;
    vx[2] = U;
    vb[0] = g_h;
    vb[1] = V;
    vb[2] = F;
    VecCreateNest(PETSC_COMM_WORLD, 3, NULL, &vx[0], &mv_x);
    VecCreateNest(PETSC_COMM_WORLD, 3, NULL, &vb[0], &mv_b);

    // Initialize the saddle-point solver.
    d_sp_solver->initializeSolverState(mv_x, mv_b);

    // Destroy the temporay vectors.
    PETScSAMRAIVectorReal::destroyPETScVector(u_p);
    PETScSAMRAIVectorReal::destroyPETScVector(g_h);
    VecDestroy(&mv_x);
    VecDestroy(&mv_b);

    d_is_initialized = true;
    d_reinitializing_solver = false;

    return;
} // initializeSolverState

void
CIBStochasticStokesSolver::setVelocityPoissonSpecifications(const PoissonSpecifications& u_problem_coefs)
{
    d_sp_solver->setVelocityPoissonSpecifications(u_problem_coefs);
    return;
} // setVelocityPoissonSpecifications

void
CIBStochasticStokesSolver::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                              RobinBcCoefStrategy<NDIM>* p_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(u_bc_coefs.size() == NDIM);
#endif
    d_u_bc_coefs = u_bc_coefs;
    d_sp_solver->setPhysicalBcCoefs(u_bc_coefs, p_bc_coef);

    return;
} // setPhysicalBcCoefs (Brennan Sprinkle)

void
CIBStochasticStokesSolver::setPhysicalBoundaryHelper(Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_helper);
#endif
    d_sp_solver->setPhysicalBoundaryHelper(bc_helper);
    return;
} // setPhysicalBoundaryHelper

bool
CIBStochasticStokesSolver::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
{
    // Create packaged vectors for the Saddle point solver.
    d_x_wide->copyVector(Pointer<SAMRAIVectorReal<NDIM, double> >(&x, false));
    d_b_wide->copyVector(Pointer<SAMRAIVectorReal<NDIM, double> >(&b, false));

    // Wrap SAMRAI vector into PETSc Vec
    Vec u_p = PETScSAMRAIVectorReal::createPETScVector(d_x_wide);
    Vec g_h = PETScSAMRAIVectorReal::createPETScVector(d_b_wide);

    // Get the Lagrange multiplier that maintains the rigidity constraint.
    // NOTE: We need L at new time to solve for it.
    Vec L;
    d_cib_strategy->getConstraintForce(&L, d_new_time);

    // Get the free body velocities to solve for.
    // NOTE: We need U at new time in the solver.
    Vec U;
    d_cib_strategy->getFreeRigidVelocities(&U, d_new_time);

    // Set the imposed velocity for all bodies in the RHS.
    Vec V;
    VecDuplicate(L, &V);
    for (unsigned part = 0; part < d_num_rigid_parts; ++part)
    {
        RigidDOFVector U_part;
        d_cib_strategy->getNewRigidBodyVelocity(part, U_part);

        // Zero-out free velocities.
        int num_free_dofs = 0;
        const FreeRigidDOFVector& solve_dofs = d_cib_strategy->getSolveRigidBodyVelocity(part, num_free_dofs);
        for (int k = 0; k < s_max_free_dofs; ++k)
        {
            if (solve_dofs[k]) U_part[k] = 0.0;
        }

        const double interp_scale = d_sp_solver->getInterpScale();
        U_part *= -interp_scale;
        d_cib_strategy->setRigidBodyVelocity(part, U_part, V);
    }

    // Get the net external force and torque on the bodies.
    Vec F;
    d_cib_strategy->getNetExternalForceTorque(&F, d_new_time);

    // RFD set-up
    Vec F_rfd, U_rfd;
    VecDuplicate(F, &F_rfd);
    VecDuplicate(U, &U_rfd);
    computeFUforRFD(F_rfd, U_rfd);

    // Make zero sub-vectors for RFD
    Pointer<SAMRAIVectorReal<NDIM, double> > g_rfd_wide = d_b_wide->cloneVector("");
    g_rfd_wide->allocateVectorData();
    g_rfd_wide->copyVector(d_b_wide);
    g_rfd_wide->setToScalar(0.0, false);
    Vec g_rfd = PETScSAMRAIVectorReal::createPETScVector(g_rfd_wide);

    Vec V_rfd;
    VecDuplicate(L, &V_rfd);
    VecZeroEntries(V_rfd);

    // Create multivector to pass it to the saddle point solver.
    std::vector<Vec> vx(3), vb(3), v_RFD(3);
    vx[0] = u_p;
    vx[1] = L;
    vx[2] = U;
    vb[0] = g_h;
    vb[1] = V;
    vb[2] = F;
    v_RFD[0] = g_rfd;
    v_RFD[1] = V_rfd;
    v_RFD[2] = F_rfd;

    Vec mv_x, mv_b, mv_rfd;
    VecCreateNest(PETSC_COMM_WORLD, 3, NULL, &vx[0], &mv_x);
    VecCreateNest(PETSC_COMM_WORLD, 3, NULL, &vb[0], &mv_b);
    VecCreateNest(PETSC_COMM_WORLD, 3, NULL, &v_RFD[0], &mv_rfd);

    // RFD solve with zeros RHS except random forces and torques
    bool rfd_converged = d_sp_solver->solveSystem(mv_x, mv_rfd);

    Pointer<CIBStochasticMethod> sib_method_ops = d_cib_strategy;

    // Set random rigid velocity, delta correlated with forces and torques from RFD solve
    sib_method_ops->setHalfTimeVelocity(U_rfd);

    // Modify the RHS to include RFD drift terms
    ComputeRFD_RHS(mv_b, mv_x);

    // Reset rigid velocity
    sib_method_ops->resetRFDVelocity();

    // Solve for velocity, pressure and Lagrange multipliers.
    // Notice that initial guess for U is provided by the implementation of the
    // IBAMR::CIBStrategy class. This is passed as an initial guess for the
    // next solve. We do not need to do anything special for its initial guess
    // to the Krylov solver.
    bool converged = d_sp_solver->solveSystem(mv_x, mv_b);

    converged = (converged && rfd_converged);

    // Extract solution.
    x.copyVector(d_x_wide);

    // Update free velocity DOFs.
    d_cib_strategy->updateNewRigidBodyVelocity(U,
                                               /*only_free_dofs*/ true,
                                               /*only_imposed_dofs*/ false,
                                               /*all_dofs*/ false);

    double half_time = 0.5 * (d_new_time + d_current_time);
    pout << "\n"
         << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
         << "Interpolating velocity on structure at time  " << half_time << "....\n"
         << std::endl;

    RefineAlgorithm<NDIM> ghost_fill_alg;
    ghost_fill_alg.registerRefine(d_wide_u_idx, d_wide_u_idx, d_wide_u_idx, NULL);
    Pointer<PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
    Pointer<RefineSchedule<NDIM> > ghost_fill_schd = ghost_fill_alg.createSchedule(hierarchy->getPatchLevel(0));
    ghost_fill_schd->fillData(half_time);

    d_cib_strategy->setInterpolatedVelocityVector(V, half_time);

    // Delete PETSc vectors.
    PETScSAMRAIVectorReal::destroyPETScVector(u_p);
    PETScSAMRAIVectorReal::destroyPETScVector(g_h);
    PETScSAMRAIVectorReal::destroyPETScVector(g_rfd);
    g_rfd_wide->resetLevels(
        g_rfd_wide->getCoarsestLevelNumber(),
        std::min(g_rfd_wide->getFinestLevelNumber(), g_rfd_wide->getPatchHierarchy()->getFinestLevelNumber()));
    g_rfd_wide->freeVectorComponents();
    g_rfd_wide.setNull();
    VecDestroy(&V);
    VecDestroy(&V_rfd);
    VecDestroy(&F_rfd);
    VecDestroy(&U_rfd);
    VecDestroy(&mv_x);
    VecDestroy(&mv_b);
    VecDestroy(&mv_rfd);

    return converged;

} // solveSystem

void
CIBStochasticStokesSolver::deallocateSolverState()
{
    // Deallocate the saddle-point solver if not re-initializing
    if (!d_reinitializing_solver)
    {
        d_sp_solver->deallocateSolverState();
    }

    // Deallocate widened patch data.
    Pointer<PatchHierarchy<NDIM> > hierarchy = d_x_wide->getPatchHierarchy();
    const int coarsest_ln = d_x_wide->getCoarsestLevelNumber();
    const int finest_ln = d_x_wide->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_wide_u_idx)) level->deallocatePatchData(d_wide_u_idx);
        if (level->checkAllocated(d_wide_f_idx)) level->deallocatePatchData(d_wide_f_idx);
    }

    // Free the vectors.
    d_x_wide.setNull();
    d_b_wide.setNull();

    return;
} // deallocateSolverState

void
CIBStochasticStokesSolver::computeFUforRFD(Vec F_rfd, Vec U_rfd)
{
    Pointer<CIBStochasticMethod> sib_method_ops = d_cib_strategy;
    double KbT = sib_method_ops->getkT();
    double L_scale = sib_method_ops->getLScale();

    int free_dofs_counter = 0;
    std::vector<PetscInt> indices;
    std::vector<PetscScalar> U_vec;
    std::vector<PetscScalar> F_vec;
    indices.reserve(d_num_rigid_parts * s_max_free_dofs);
    U_vec.reserve(d_num_rigid_parts * s_max_free_dofs);
    F_vec.reserve(d_num_rigid_parts * s_max_free_dofs);
    for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
    {
        int num_free_dofs;
        const FRDV& solve_dofs = d_cib_strategy->getSolveRigidBodyVelocity(part, num_free_dofs);
        if (SAMRAI_MPI::getRank() == 0)
        {
            if (num_free_dofs)
            {
                // Initialize external force and torque.
                RDV Fr;
                Eigen::Vector3d FU_rand, TW_rand, F_rand, T_rand, U_rand, W_rand;
                for (int d = 0; d < NDIM; ++d)
                {
                    RNG::genrandn(&FU_rand(d));
                    RNG::genrandn(&TW_rand(d));
                }
                F_rand = (KbT / L_scale) * FU_rand;
                T_rand = (KbT)*TW_rand;
                d_cib_strategy->eigenToRDV(F_rand, T_rand, Fr);

                U_rand = (L_scale)*FU_rand;
                W_rand = TW_rand;
                RDV Ur;
                d_cib_strategy->eigenToRDV(U_rand, W_rand, Ur);

                for (int k = 0; k < s_max_free_dofs; ++k)
                {
                    if (solve_dofs[k])
                    {
                        U_vec.push_back(Ur[k]);
                        F_vec.push_back(Fr[k]);
                        indices.push_back(free_dofs_counter);
                        ++free_dofs_counter;
                    }
                }
            }
        }
    }

    // Create PETSc Vecs for free DOFs.
    const int n = free_dofs_counter;
    if (n)
    {
        VecSetValues(U_rfd, n, &indices[0], &U_vec[0], INSERT_VALUES);
        VecSetValues(F_rfd, n, &indices[0], &F_vec[0], INSERT_VALUES);
    }

    VecAssemblyBegin(U_rfd);
    VecAssemblyEnd(U_rfd);
    VecAssemblyBegin(F_rfd);
    VecAssemblyEnd(F_rfd);

    return;
} // computeRFDforcesAndDisplacements (Brennan Sprinkle)

void
CIBStochasticStokesSolver::ComputeRFD_RHS(Vec b, Vec y)
{
    Pointer<CIBStochasticMethod> sib_method_ops = d_cib_strategy;

    // Get some constants
    const double gamma = d_sp_solver->getSpreadScale();
    const double beta = d_sp_solver->getInterpScale();
    const double half_time = 0.5 * (d_new_time + d_current_time);

    int total_comps, free_comps = 0;
    Vec *vb, *vy;
    VecNestGetSubVecs(b, NULL, &vb);
    VecNestGetSubVecs(y, &total_comps, &vy);
    VecGetSize(vb[2], &free_comps);

    Pointer<SAMRAIVectorReal<NDIM, double> > vb0, vy0;
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVector(vb[0], &vb0);
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVectorRead(vy[0], &vy0);

    // Get the individual components.
    Pointer<SAMRAIVectorReal<NDIM, double> > g_h = vb0->cloneVector("");
    g_h->allocateVectorData();
    g_h->copyVector(vb0);
    g_h->setToScalar(0.0, false); // added so we can use g_h for speading. Hopefully does not effect valure of vb

    const int g_data_idx = g_h->getComponentDescriptorIndex(0);

    Pointer<SAMRAIVectorReal<NDIM, double> > u_p = vy0;

    // Fill ghost cells of u
    int u_data_idx = u_p->getComponentDescriptorIndex(0);
    typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps;
    InterpolationTransactionComponent u_component(u_data_idx,
                                                  DATA_REFINE_TYPE,
                                                  USE_CF_INTERPOLATION,
                                                  DATA_COARSEN_TYPE,
                                                  BDRY_EXTRAP_TYPE,
                                                  CONSISTENT_TYPE_2_BDRY,
                                                  d_u_bc_coefs,
                                                  NULL);
    transaction_comps.push_back(u_component);
    // solver->d_hier_bdry_fill->initializeOperatorState(transaction_comps, u_p->getPatchHierarchy());
    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(transaction_comps, u_p->getPatchHierarchy());
    hier_bdry_fill->setHomogeneousBc(false);
    hier_bdry_fill->fillData(half_time);

    double delta = sib_method_ops->getRFdelta();
    double factor = -1.0 / delta;
    double spread_scale = gamma / delta;

    Vec F_tilde, Ds1, Ds2, Ds3;
    VecDuplicate(vb[2], &F_tilde);
    VecDuplicate(vb[1], &Ds1);
    VecDuplicate(vb[1], &Ds2);
    VecDuplicate(vb[1], &Ds3);

    for (int pm = 0; pm < 2; pm++)
    {
        // Set IB positions to RFD positions
        sib_method_ops->moveLagrangianData(delta);

        // Compute F +-= (1/delta) * K^T(Q^+-) * Lambda^RFD
        d_cib_strategy->computeNetRigidGeneralizedForce(vy[1], // maybe replace with vx[2]
                                                        F_tilde,
                                                        /*only_free_dofs*/ true,
                                                        /*only_imposed_dofs*/ false);

        VecAXPY(vb[2], factor, F_tilde);

        // Compute D^S1 +-= (1/delta) * K(Q^+-) * U^RFD
        VecZeroEntries(Ds1); // NEED THIS
        d_cib_strategy->setRigidBodyVelocity(vy[2], Ds1, /*only_free_dofs*/ true, /*only_imposed_dofs*/ false);
        VecAXPY(vb[1], factor, Ds1);

        // D^S2 -+= (1/delta) *  J(Q^+-) * u.
        VecZeroEntries(Ds2);
        d_cib_strategy->setInterpolatedVelocityVector(Ds2, half_time);
        sib_method_ops->interpolateVelocity(u_data_idx,
                                            std::vector<Pointer<CoarsenSchedule<NDIM> > >(),
                                            std::vector<Pointer<RefineSchedule<NDIM> > >(),
                                            half_time);

        d_cib_strategy->getInterpolatedVelocity(Ds2, half_time, beta);
        VecAXPY(vb[1], -1.0 * factor, Ds2);

        // g^S3 -+= (1/delta) *  S(Q^+-) * Lambda.
        d_cib_strategy->setConstraintForce(vy[1], half_time, spread_scale);
        sib_method_ops->spreadForce(g_data_idx, NULL, std::vector<Pointer<RefineSchedule<NDIM> > >(), half_time);

        // Set Factor for negative computations
        factor *= -1.0;
        delta *= -1.0;
        spread_scale *= -1.0;
    }
    // Reset IB positions
    sib_method_ops->moveLagrangianData(0);

    Pointer<StaggeredStokesSolver> LInv = d_sp_solver->getStokesSolver();
    IBTK::PETScKrylovLinearSolver* petsc_stokes_krylov_solver =
        dynamic_cast<IBTK::PETScKrylovLinearSolver*>(LInv.getPointer());

    int max_its = petsc_stokes_krylov_solver->getMaxIterations();

    petsc_stokes_krylov_solver->setKSPType("gmres");
    petsc_stokes_krylov_solver->setMaxIterations(2000);
    petsc_stokes_krylov_solver->setHomogeneousBc(false);
    petsc_stokes_krylov_solver->setLoggingEnabled(true);

    petsc_stokes_krylov_solver->solveSystem(*u_p, *g_h);

    petsc_stokes_krylov_solver->setKSPType("richardson");
    petsc_stokes_krylov_solver->setMaxIterations(max_its);
    petsc_stokes_krylov_solver->setHomogeneousBc(true);
    petsc_stokes_krylov_solver->setLoggingEnabled(false);

    u_data_idx = u_p->getComponentDescriptorIndex(0);
    transaction_comps[0] = InterpolationTransactionComponent(u_data_idx,
                                                             DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_u_bc_coefs,
                                                             NULL);

    hier_bdry_fill->initializeOperatorState(transaction_comps, u_p->getPatchHierarchy());
    hier_bdry_fill->setHomogeneousBc(false);
    hier_bdry_fill->fillData(half_time);

    // D^S3 = J(Q) * u^S3.
    VecZeroEntries(Ds3);
    d_cib_strategy->setInterpolatedVelocityVector(Ds3, half_time);
    sib_method_ops->interpolateVelocity(u_data_idx,
                                        std::vector<Pointer<CoarsenSchedule<NDIM> > >(),
                                        std::vector<Pointer<RefineSchedule<NDIM> > >(),
                                        half_time);

    d_cib_strategy->getInterpolatedVelocity(Ds3, half_time, beta);
    VecAXPY(vb[1], 1.0, Ds3);

    // Destroy temporary vectors
    g_h->resetLevels(g_h->getCoarsestLevelNumber(),
                     std::min(g_h->getFinestLevelNumber(), g_h->getPatchHierarchy()->getFinestLevelNumber()));
    g_h->freeVectorComponents();
    g_h.setNull();

    // VecDestroy(&U);
    VecDestroy(&F_tilde);
    VecDestroy(&Ds1);
    VecDestroy(&Ds2);
    VecDestroy(&Ds3);

    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVector(vb[0], &vb0);
    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(vy[0], &vy0);
}

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
