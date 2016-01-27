// Filename: INSCollocatedHierarchyIntegrator.cpp
// Created on 24 Aug 2011 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenSchedule.h"
#include "ComponentSelector.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "GriddingAlgorithm.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "HierarchyFaceDataOpsReal.h"
#include "IBAMR_config.h"
#include "Index.h"
#include "IntVector.h"
#include "LocationIndexRobinBcCoefs.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchCellDataOpsReal.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "RefinePatchStrategy.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAIVectorReal.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "VisItDataWriter.h"
#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSCollocatedConvectiveOperatorManager.h"
#include "ibamr/INSCollocatedHierarchyIntegrator.h"
#include "ibamr/INSCollocatedVelocityBcCoef.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSIntermediateVelocityBcCoef.h"
#include "ibamr/INSProjectionBcCoef.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/ibtk_enums.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define ADVECT_STABLEDT_FC IBAMR_FC_FUNC_(advect_stabledt2d, ADVECT_STABLEDT2D)
#define NAVIER_STOKES_ADV_SOURCE_FC IBAMR_FC_FUNC_(navier_stokes_adv_source2d, NAVIER_STOKES_ADV_SOURCE2D)
#define NAVIER_STOKES_CONS_SOURCE_FC IBAMR_FC_FUNC_(navier_stokes_cons_source2d, NAVIER_STOKES_CONS_SOURCE2D)
#define NAVIER_STOKES_SKEW_SYM_SOURCE_FC                                                                               \
    IBAMR_FC_FUNC_(navier_stokes_skew_sym_source2d, NAVIER_STOKES_SKEW_SYM_SOURCE2D)
#endif

#if (NDIM == 3)
#define ADVECT_STABLEDT_FC IBAMR_FC_FUNC_(advect_stabledt3d, ADVECT_STABLEDT3D)
#define NAVIER_STOKES_ADV_SOURCE_FC IBAMR_FC_FUNC_(navier_stokes_adv_source3d, NAVIER_STOKES_ADV_SOURCE3D)
#define NAVIER_STOKES_CONS_SOURCE_FC IBAMR_FC_FUNC_(navier_stokes_cons_source3d, NAVIER_STOKES_CONS_SOURCE3D)
#define NAVIER_STOKES_SKEW_SYM_SOURCE_FC                                                                               \
    IBAMR_FC_FUNC_(navier_stokes_skew_sym_source3d, NAVIER_STOKES_SKEW_SYM_SOURCE3D)
#endif

extern "C" {
void ADVECT_STABLEDT_FC(const double*,
#if (NDIM == 2)
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const double*,
                        const double*,
#endif
#if (NDIM == 3)
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const double*,
                        const double*,
                        const double*,
#endif
                        double&);

void NAVIER_STOKES_ADV_SOURCE_FC(
#if (NDIM == 2)
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
#endif
#if (NDIM == 3)
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
    const double*,
#endif
    const double*,
    double*);

void NAVIER_STOKES_CONS_SOURCE_FC(
#if (NDIM == 2)
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
#endif
#if (NDIM == 3)
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
    const double*,
#endif
    const double*,
    double*);

void NAVIER_STOKES_SKEW_SYM_SOURCE_FC(
#if (NDIM == 2)
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
#endif
#if (NDIM == 3)
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
    const double*,
#endif
    const double*,
    double*);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int FACEG = 1;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSCollocatedHierarchyIntegrator::INSCollocatedHierarchyIntegrator(const std::string& object_name,
                                                                   Pointer<Database> input_db,
                                                                   bool register_for_restart)
    : INSHierarchyIntegrator(object_name,
                             input_db,
                             new CellVariable<NDIM, double>(object_name + "::U", NDIM),
                             new CellVariable<NDIM, double>(object_name + "::P"),
                             new CellVariable<NDIM, double>(object_name + "::F", NDIM),
                             new CellVariable<NDIM, double>(object_name + "::Q"),
                             register_for_restart)
{
    // Check to make sure the time stepping type is supported.
    switch (d_viscous_time_stepping_type)
    {
    case BACKWARD_EULER:
    case FORWARD_EULER:
    case TRAPEZOIDAL_RULE:
        break;
    default:
        TBOX_ERROR(d_object_name << "::INSCollocatedHierarchyIntegrator():\n"
                                 << "  unsupported viscous time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_viscous_time_stepping_type)
                                 << " \n"
                                 << "  valid choices are: BACKWARD_EULER, FORWARD_EULER, TRAPEZOIDAL_RULE\n");
    }
    switch (d_convective_time_stepping_type)
    {
    case ADAMS_BASHFORTH:
    case FORWARD_EULER:
    case MIDPOINT_RULE:
    case TRAPEZOIDAL_RULE:
        break;
    default:
        TBOX_ERROR(d_object_name << "::INSStaggeredHierarchyIntegrator():\n"
                                 << "  unsupported convective time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_convective_time_stepping_type)
                                 << " \n"
                                 << "  valid choices are: ADAMS_BASHFORTH, FORWARD_EULER, "
                                    "MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }
    if (is_multistep_time_stepping_type(d_convective_time_stepping_type))
    {
        switch (d_init_convective_time_stepping_type)
        {
        case FORWARD_EULER:
        case MIDPOINT_RULE:
        case TRAPEZOIDAL_RULE:
            break;
        default:
            TBOX_ERROR(d_object_name << "::INSStaggeredHierarchyIntegrator():\n"
                                     << "  unsupported initial convective time stepping type: "
                                     << enum_to_string<TimeSteppingType>(d_init_convective_time_stepping_type)
                                     << " \n"
                                     << "  valid choices are: FORWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
        }
    }

    // Check to see whether the solver types have been set.
    d_velocity_solver_type = CCPoissonSolverManager::UNDEFINED;
    d_velocity_precond_type = CCPoissonSolverManager::UNDEFINED;
    if (input_db->keyExists("velocity_solver_type"))
        d_velocity_solver_type = input_db->getString("velocity_solver_type");
    if (input_db->keyExists("velocity_precond_type"))
        d_velocity_precond_type = input_db->getString("velocity_precond_type");

    d_pressure_solver_type = CCPoissonSolverManager::UNDEFINED;
    d_pressure_precond_type = CCPoissonSolverManager::UNDEFINED;
    if (input_db->keyExists("pressure_solver_type"))
        d_pressure_solver_type = input_db->getString("pressure_solver_type");
    if (input_db->keyExists("pressure_precond_type"))
        d_pressure_precond_type = input_db->getString("pressure_precond_type");

    d_regrid_projection_solver_type = CCPoissonSolverManager::UNDEFINED;
    d_regrid_projection_precond_type = CCPoissonSolverManager::UNDEFINED;
    if (input_db->keyExists("regrid_projection_solver_type"))
        d_regrid_projection_solver_type = input_db->getString("regrid_projection_solver_type");
    if (input_db->keyExists("regrid_projection_precond_type"))
        d_regrid_projection_precond_type = input_db->getString("regrid_projection_precond_type");

    // Check to see whether the convective operator type has been set.
    d_convective_op_type = INSCollocatedConvectiveOperatorManager::DEFAULT;
    if (input_db->keyExists("convective_op_type"))
        d_convective_op_type = input_db->getString("convective_op_type");
    else if (input_db->keyExists("convective_operator_type"))
        d_convective_op_type = input_db->getString("convective_operator_type");
    else if (input_db->keyExists("default_convective_op_type"))
        d_convective_op_type = input_db->getString("default_convective_op_type");
    else if (input_db->keyExists("default_convective_operator_type"))
        d_convective_op_type = input_db->getString("default_convective_operator_type");

    // Check to see what kind of projection method to use.
    d_projection_method_type = PRESSURE_INCREMENT;
    if (input_db->keyExists("proj_method_type"))
        d_projection_method_type = string_to_enum<ProjectionMethodType>(input_db->getString("proj_method_type"));
    else if (input_db->keyExists("projection_method_type"))
        d_projection_method_type = string_to_enum<ProjectionMethodType>(input_db->getString("projection_method_type"));

    d_using_2nd_order_pressure_update = true;
    if (input_db->keyExists("use_2nd_order_pressure_update"))
        d_using_2nd_order_pressure_update = input_db->getBool("use_2nd_order_pressure_update");
    else if (input_db->keyExists("using_2nd_order_pressure_update"))
        d_using_2nd_order_pressure_update = input_db->getBool("using_2nd_order_pressure_update");
    else if (input_db->keyExists("use_second_order_pressure_update"))
        d_using_2nd_order_pressure_update = input_db->getBool("use_second_order_pressure_update");
    else if (input_db->keyExists("using_second_order_pressure_update"))
        d_using_2nd_order_pressure_update = input_db->getBool("using_second_order_pressure_update");

    // Setup physical boundary conditions objects.
    d_U_bc_coefs.resize(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_U_bc_coefs[d] = new INSCollocatedVelocityBcCoef(d, this, d_bc_coefs, d_traction_bc_type);
    }

    // Initialize all variables.  The velocity, pressure, body force, and fluid
    // source variables were created above in the constructor for the
    // INSHierarchyIntegrator base class.
    d_U_var = INSHierarchyIntegrator::d_U_var;
    d_u_ADV_var = new FaceVariable<NDIM, double>(d_object_name + "::u_ADV");
    d_P_var = INSHierarchyIntegrator::d_P_var;
    d_F_var = INSHierarchyIntegrator::d_F_var;
    d_Q_var = INSHierarchyIntegrator::d_Q_var;
    d_N_old_var = new CellVariable<NDIM, double>(d_object_name + "::N_old", NDIM);

#if (NDIM == 2)
    d_Omega_var = new CellVariable<NDIM, double>(d_object_name + "::Omega");
#endif
#if (NDIM == 3)
    d_Omega_var = new CellVariable<NDIM, double>(d_object_name + "::Omega", NDIM);
#endif
    d_Div_U_var = new CellVariable<NDIM, double>(d_object_name + "::Div_U");
    d_Div_u_ADV_var = new CellVariable<NDIM, double>(d_object_name + "::Div_u_ADV");
#if (NDIM == 3)
    d_Omega_Norm_var = new CellVariable<NDIM, double>(d_object_name + "::|Omega|_2");
#endif
    d_Grad_P_var = new CellVariable<NDIM, double>(d_object_name + "::Grad_P", NDIM);
    d_Phi_var = new CellVariable<NDIM, double>(d_object_name + "::Phi");
    d_Grad_Phi_cc_var = new CellVariable<NDIM, double>(d_object_name + "::Grad_Phi_cc", NDIM);
    d_Grad_Phi_fc_var = new FaceVariable<NDIM, double>(d_object_name + "::Grad_Phi_fc");
    d_F_div_var = new CellVariable<NDIM, double>(d_object_name + "::F_div", NDIM);
    return;
} // INSCollocatedHierarchyIntegrator

INSCollocatedHierarchyIntegrator::~INSCollocatedHierarchyIntegrator()
{
    delete d_fill_after_regrid_phys_bdry_bc_op;
    d_fill_after_regrid_phys_bdry_bc_op = NULL;
    d_velocity_solver.setNull();
    d_pressure_solver.setNull();
    if (d_U_rhs_vec) d_U_rhs_vec->freeVectorComponents();
    if (d_U_adv_vec) d_U_adv_vec->freeVectorComponents();
    if (d_N_vec) d_N_vec->freeVectorComponents();
    if (d_Phi_rhs_vec) d_Phi_rhs_vec->freeVectorComponents();
    for (unsigned int k = 0; k < d_U_nul_vecs.size(); ++k)
    {
        if (d_U_nul_vecs[k]) d_U_nul_vecs[k]->freeVectorComponents();
    }
    return;
} // ~INSCollocatedHierarchyIntegrator

Pointer<ConvectiveOperator>
INSCollocatedHierarchyIntegrator::getConvectiveOperator()
{
    if (d_creeping_flow)
    {
        d_convective_op.setNull();
    }
    else if (!d_convective_op)
    {
        INSCollocatedConvectiveOperatorManager* convective_op_manager =
            INSCollocatedConvectiveOperatorManager::getManager();
        d_convective_op = convective_op_manager->allocateOperator(d_convective_op_type,
                                                                  d_object_name + "::ConvectiveOperator",
                                                                  d_convective_op_input_db,
                                                                  d_convective_difference_form,
                                                                  d_U_star_bc_coefs);
        d_convective_op_needs_init = true;
    }
    return d_convective_op;
} // getConvectiveOperator

Pointer<PoissonSolver>
INSCollocatedHierarchyIntegrator::getVelocitySubdomainSolver()
{
    if (!d_velocity_solver)
    {
        d_velocity_solver = CCPoissonSolverManager::getManager()->allocateSolver(d_velocity_solver_type,
                                                                                 d_object_name + "::velocity_solver",
                                                                                 d_velocity_solver_db,
                                                                                 "velocity_",
                                                                                 d_velocity_precond_type,
                                                                                 d_object_name + "::velocity_precond",
                                                                                 d_velocity_precond_db,
                                                                                 "velocity_pc_");
        d_velocity_solver_needs_init = true;
    }
    return d_velocity_solver;
} // getVelocitySubdomainSolver

Pointer<PoissonSolver>
INSCollocatedHierarchyIntegrator::getPressureSubdomainSolver()
{
    if (!d_pressure_solver)
    {
        d_pressure_solver = CCPoissonSolverManager::getManager()->allocateSolver(d_pressure_solver_type,
                                                                                 d_object_name + "::pressure_solver",
                                                                                 d_pressure_solver_db,
                                                                                 "pressure_",
                                                                                 d_pressure_precond_type,
                                                                                 d_object_name + "::pressure_precond",
                                                                                 d_pressure_precond_db,
                                                                                 "pressure_pc_");
        d_pressure_solver_needs_init = true;
    }
    return d_pressure_solver;
} // getPressureSubdomainSolver

void
INSCollocatedHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                                Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;

    if (d_rho_var)
    {
        TBOX_ERROR("INSCollocatedHierarchyIntegrator::initializeHierarchyIntegrator():\n"
                   << "  variable density solver not presently supported.\n");
    }

    if (d_gridding_alg->getMaxLevels() > 1 && d_projection_method_type == PRESSURE_INCREMENT)
    {
        pout << "\n"
             << "   "
                "******************************************************************************"
                "**\n"
             << "   * INSCollocatedHierarchyIntegrator:                                        "
                "    "
                "*\n"
             << "   *                                                                          "
                "    "
                "*\n"
             << "   * WARNING: PRESSURE_INCREMENT projection method exhibits very slow "
                "convergence "
                "*\n"
             << "   *          when used with locally refined grids.                           "
                "    "
                "*\n"
             << "   *                                                                          "
                "    "
                "*\n"
             << "   *          Sugest changing projection method type to PRESSURE_UPDATE.      "
                "    "
                "*\n"
             << "   "
                "******************************************************************************"
                "\n\n";
        // NOTE: It is not clear what is causing this lack of convergence.
        //
        // The basic pressure-increment projection method with is:
        //
        // (1)      rho{[u(*) - u(n)]/dt + N(n+1/2)} = - Grad_h P(n-1/2) + mu Lap_h[u(*) +
        // u(n)]/2
        //
        // (2)      rho[u(n+1) - u(*)]/dt = - Grad_h Phi
        //          Div_h u(n+1)          = 0
        //
        // (3)      P(n+1/2) = P(n-1/2) + Phi
        //
        // This formulation uses a first-order accurate pressure update formula.
        // Te second-order accurate (Brown-Cortez-Minion) pressure update is:
        //
        // (3-alt)  P(n+1/2) = P(n-1/2) + (I - mu*dt/2*rho Lap_h) Phi
        //
        // Notice that this algorithm is constructed by assuming that Grad_h and
        // Lap_h commute.  This is true only on uniform periodic grids.  On
        // locally-refined grids or in the presence of physical boundaries, this
        // is no longer the case.
        //
        // One guess is that the poor performance of the algorithm is related to
        // this lack of commutativity, and could be corrected by lagging Grad P
        // instead of P.  E.g., by keeping track of a vector-valued quantity GP,
        // as follows:
        //
        // (1)      rho{[u(*) - u(n)]/dt + N(n+1/2)} = - GP(n-1/2) + mu Lap_h[u(*) + u(n)]/2
        //
        // (2)      rho[u(n+1) - u(*)]/dt = - Grad_h Phi
        //          Div_h u(n+1)          = 0
        //
        // (3)      GP(n+1/2) = GP(n-1/2) + Grad_h Phi
        //
        // This possibility has not been investigated at this point.
    }

    // Setup solvers.
    if (d_velocity_solver_type == CCPoissonSolverManager::UNDEFINED)
    {
        d_velocity_solver_type = CCPoissonSolverManager::DEFAULT_KRYLOV_SOLVER;
    }

    if (d_velocity_precond_type == CCPoissonSolverManager::UNDEFINED)
    {
        const int max_levels = gridding_alg->getMaxLevels();
        if (max_levels == 1)
        {
            d_velocity_precond_type = CCPoissonSolverManager::DEFAULT_LEVEL_SOLVER;
        }
        else
        {
            d_velocity_precond_type = CCPoissonSolverManager::DEFAULT_FAC_PRECONDITIONER;
        }
        d_velocity_precond_db->putInteger("max_iterations", 1);
    }

    if (d_pressure_solver_type == CCPoissonSolverManager::UNDEFINED)
    {
        d_pressure_solver_type = CCPoissonSolverManager::DEFAULT_KRYLOV_SOLVER;
    }

    if (d_pressure_precond_type == CCPoissonSolverManager::UNDEFINED)
    {
        const int max_levels = gridding_alg->getMaxLevels();
        if (max_levels == 1)
        {
            d_pressure_precond_type = CCPoissonSolverManager::DEFAULT_LEVEL_SOLVER;
        }
        else
        {
            d_pressure_precond_type = CCPoissonSolverManager::DEFAULT_FAC_PRECONDITIONER;
        }
        d_pressure_precond_db->putInteger("max_iterations", 1);
    }

    if (d_regrid_projection_solver_type == CCPoissonSolverManager::UNDEFINED)
    {
        d_regrid_projection_solver_type = CCPoissonSolverManager::DEFAULT_KRYLOV_SOLVER;
    }

    if (d_regrid_projection_precond_type == CCPoissonSolverManager::UNDEFINED)
    {
        const int max_levels = gridding_alg->getMaxLevels();
        if (max_levels == 1)
        {
            d_regrid_projection_precond_type = CCPoissonSolverManager::DEFAULT_LEVEL_SOLVER;
        }
        else
        {
            d_regrid_projection_precond_type = CCPoissonSolverManager::DEFAULT_FAC_PRECONDITIONER;
        }
        d_regrid_projection_precond_db->putInteger("max_iterations", 1);
    }

    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    d_hier_cc_data_ops =
        hier_ops_manager->getOperationsDouble(new CellVariable<NDIM, double>("cc_var"), hierarchy, true);
    d_hier_fc_data_ops =
        hier_ops_manager->getOperationsDouble(new FaceVariable<NDIM, double>("fc_var"), hierarchy, true);
    d_hier_math_ops = buildHierarchyMathOps(d_hierarchy);

    // Register state variables that are maintained by the
    // INSCollocatedHierarchyIntegrator.
    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> face_ghosts = FACEG;
    const IntVector<NDIM> no_ghosts = 0;

    registerVariable(d_U_current_idx,
                     d_U_new_idx,
                     d_U_scratch_idx,
                     d_U_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE",
                     d_U_init);

    registerVariable(d_u_ADV_current_idx,
                     d_u_ADV_new_idx,
                     d_u_ADV_scratch_idx,
                     d_u_ADV_var,
                     face_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_P_current_idx,
                     d_P_new_idx,
                     d_P_scratch_idx,
                     d_P_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE",
                     d_P_init);

    if (d_F_fcn)
    {
        registerVariable(d_F_current_idx,
                         d_F_new_idx,
                         d_F_scratch_idx,
                         d_F_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE",
                         d_F_fcn);
    }
    else
    {
        d_F_current_idx = -1;
        d_F_new_idx = -1;
        d_F_scratch_idx = -1;
    }

    if (d_Q_fcn)
    {
        registerVariable(d_Q_current_idx,
                         d_Q_new_idx,
                         d_Q_scratch_idx,
                         d_Q_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSTANT_REFINE",
                         d_Q_fcn);
    }
    else
    {
        d_Q_current_idx = -1;
        d_Q_new_idx = -1;
        d_Q_scratch_idx = -1;
    }

    registerVariable(d_N_old_current_idx,
                     d_N_old_new_idx,
                     d_N_old_scratch_idx,
                     d_N_old_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    // Register plot variables that are maintained by the
    // INSCollocatedHierarchyIntegrator.
    registerVariable(d_Omega_idx, d_Omega_var, no_ghosts, getCurrentContext());
    registerVariable(d_Div_U_idx, d_Div_U_var, cell_ghosts, getCurrentContext());
    registerVariable(d_Div_u_ADV_idx, d_Div_u_ADV_var, no_ghosts, getCurrentContext());

// Register scratch variables that are maintained by the
// INSCollocatedHierarchyIntegrator.
#if (NDIM == 3)
    registerVariable(d_Omega_Norm_idx, d_Omega_Norm_var, no_ghosts);
#endif
    registerVariable(d_Grad_P_idx, d_Grad_P_var, no_ghosts);
    registerVariable(d_Phi_idx, d_Phi_var, cell_ghosts);
    registerVariable(d_Grad_Phi_cc_idx, d_Grad_Phi_cc_var, no_ghosts);
    registerVariable(d_Grad_Phi_fc_idx, d_Grad_Phi_fc_var, no_ghosts);
    if (d_Q_fcn)
    {
        registerVariable(d_F_div_idx, d_F_div_var, no_ghosts);
    }
    else
    {
        d_F_div_idx = -1;
    }

    // Register variables for plotting.
    if (d_visit_writer)
    {
        if (d_output_U)
        {
            d_visit_writer->registerPlotQuantity("U", "VECTOR", d_U_current_idx, 0, d_U_scale);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == 0) d_visit_writer->registerPlotQuantity("U_x", "SCALAR", d_U_current_idx, d, d_U_scale);
                if (d == 1) d_visit_writer->registerPlotQuantity("U_y", "SCALAR", d_U_current_idx, d, d_U_scale);
                if (d == 2) d_visit_writer->registerPlotQuantity("U_z", "SCALAR", d_U_current_idx, d, d_U_scale);
            }
        }

        if (d_output_P)
        {
            d_visit_writer->registerPlotQuantity("P", "SCALAR", d_P_current_idx, 0, d_P_scale);
        }

        if (d_F_fcn && d_output_F)
        {
            d_visit_writer->registerPlotQuantity("F", "VECTOR", d_F_current_idx, 0, d_F_scale);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == 0) d_visit_writer->registerPlotQuantity("F_x", "SCALAR", d_F_current_idx, d, d_F_scale);
                if (d == 1) d_visit_writer->registerPlotQuantity("F_y", "SCALAR", d_F_current_idx, d, d_F_scale);
                if (d == 2) d_visit_writer->registerPlotQuantity("F_z", "SCALAR", d_F_current_idx, d, d_F_scale);
            }
        }

        if (d_Q_fcn && d_output_Q)
        {
            d_visit_writer->registerPlotQuantity("Q", "SCALAR", d_Q_current_idx, 0, d_Q_scale);
        }

        if (d_output_Omega)
        {
#if (NDIM == 2)
            d_visit_writer->registerPlotQuantity("Omega", "SCALAR", d_Omega_idx, 0, d_Omega_scale);
#endif
#if (NDIM == 3)
            d_visit_writer->registerPlotQuantity("Omega", "VECTOR", d_Omega_idx, 0, d_Omega_scale);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == 0) d_visit_writer->registerPlotQuantity("Omega_x", "SCALAR", d_Omega_idx, d, d_Omega_scale);
                if (d == 1) d_visit_writer->registerPlotQuantity("Omega_y", "SCALAR", d_Omega_idx, d, d_Omega_scale);
                if (d == 2) d_visit_writer->registerPlotQuantity("Omega_z", "SCALAR", d_Omega_idx, d, d_Omega_scale);
            }
#endif
        }

        if (d_output_Div_U)
        {
            d_visit_writer->registerPlotQuantity("Div U", "SCALAR", d_Div_U_idx, 0, d_Div_U_scale);
            d_visit_writer->registerPlotQuantity("Div u_ADV", "SCALAR", d_Div_u_ADV_idx, 0, d_Div_U_scale);
        }
    }

    // Setup a specialized coarsen algorithm.
    Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<CoarsenOperator<NDIM> > coarsen_op;
    coarsen_op = grid_geom->lookupCoarsenOperator(d_U_var, "CONSERVATIVE_COARSEN");
    coarsen_alg->registerCoarsen(d_U_scratch_idx, d_U_scratch_idx, coarsen_op);
    coarsen_op = grid_geom->lookupCoarsenOperator(d_u_ADV_var, "CONSERVATIVE_COARSEN");
    coarsen_alg->registerCoarsen(d_u_ADV_scratch_idx, d_u_ADV_scratch_idx, coarsen_op);
    registerCoarsenAlgorithm(d_object_name + "::CONVECTIVE_OP", coarsen_alg);

    // Setup the component solvers.
    d_velocity_solver = getVelocitySubdomainSolver();
    d_pressure_solver = getPressureSubdomainSolver();

    // Setup the convective operator.
    d_convective_op = getConvectiveOperator();

    // Setup a boundary op to set velocity boundary conditions on regrid.
    d_fill_after_regrid_phys_bdry_bc_op = new CartCellRobinPhysBdryOp(d_U_scratch_idx, d_U_star_bc_coefs, false);

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
} // initializeHierarchyIntegrator

void
INSCollocatedHierarchyIntegrator::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                           Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    HierarchyIntegrator::initializePatchHierarchy(hierarchy, gridding_alg);

    // Project the velocity field if this is the initial time step.  Note that
    // regridProjection() also has the effect of initializing u_ADV.
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    if (initial_time)
    {
        regridProjection();
        synchronizeHierarchyData(CURRENT_DATA);
    }

    // When necessary, initialize the value of the advection velocity registered
    // with a coupled advection-diffusion solver.
    if (d_adv_diff_hier_integrator)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int U_adv_diff_current_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getCurrentContext());
        if (isAllocatedPatchData(U_adv_diff_current_idx))
        {
            d_hier_fc_data_ops->copyData(U_adv_diff_current_idx, d_u_ADV_current_idx);
        }
    }
    return;
} // initializePatchHierarhcy

void
INSCollocatedHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                               const double new_time,
                                                               const int num_cycles)
{
    INSHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Keep track of the number of cycles to be used for the present integration
    // step.
    if ((d_current_num_cycles == 1) &&
        (d_convective_time_stepping_type == MIDPOINT_RULE || d_convective_time_stepping_type == TRAPEZOIDAL_RULE))
    {
        TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                 << "  time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_convective_time_stepping_type)
                                 << " requires num_cycles > 1.\n"
                                 << "  at current time step, num_cycles = "
                                 << d_current_num_cycles
                                 << "\n");
    }

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data, new_time);
    }

    // Setup the operators and solvers.
    reinitializeOperatorsAndSolvers(current_time, new_time);

    // Allocate solver vectors.
    d_U_rhs_vec->allocateVectorData(current_time);
    d_U_rhs_vec->setToScalar(0.0);
    d_Phi_rhs_vec->allocateVectorData(current_time);
    d_Phi_rhs_vec->setToScalar(0.0);
    if (!d_creeping_flow)
    {
        d_U_adv_vec->allocateVectorData(current_time);
        d_U_adv_vec->setToScalar(0.0);
        d_N_vec->allocateVectorData(current_time);
        d_N_vec->setToScalar(0.0);
    }

    // Initialize the right-hand side terms.
    const double rho = d_problem_coefs.getRho();
    const double mu = d_problem_coefs.getMu();
    const double lambda = d_problem_coefs.getLambda();
    double K_rhs = 0.0;
    switch (d_viscous_time_stepping_type)
    {
    case BACKWARD_EULER:
        K_rhs = 0.0;
        break;
    case FORWARD_EULER:
        K_rhs = 1.0;
        break;
    case TRAPEZOIDAL_RULE:
        K_rhs = 0.5;
        break;
    default:
        TBOX_ERROR("this statment should not be reached");
    }
    PoissonSpecifications U_rhs_problem_coefs(d_object_name + "::U_rhs_problem_coefs");
    U_rhs_problem_coefs.setCConstant((rho / dt) - K_rhs * lambda);
    U_rhs_problem_coefs.setDConstant(+K_rhs * mu);
    const int U_rhs_idx = d_U_rhs_vec->getComponentDescriptorIndex(0);
    const Pointer<CellVariable<NDIM, double> > U_rhs_var = d_U_rhs_vec->getComponentVariable(0);
    d_hier_cc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        d_hier_math_ops->laplace(U_rhs_idx,
                                 U_rhs_var,
                                 U_rhs_problem_coefs,
                                 d_U_scratch_idx,
                                 d_U_var,
                                 d_U_bdry_bc_fill_op,
                                 current_time,
                                 0.0,
                                 -1,
                                 Pointer<CellVariable<NDIM, double> >(NULL),
                                 axis,
                                 axis);
    }

    // Set the initial guess.
    d_hier_cc_data_ops->copyData(d_U_new_idx, d_U_current_idx);
    d_hier_fc_data_ops->copyData(d_u_ADV_new_idx, d_u_ADV_current_idx);
    d_hier_cc_data_ops->copyData(d_P_new_idx, d_P_current_idx);

    // Initialize any registered advection-diffusion solver.
    if (d_adv_diff_hier_integrator)
    {
        const int adv_diff_num_cycles = d_adv_diff_hier_integrator->getNumberOfCycles();
        if (adv_diff_num_cycles != d_current_num_cycles && d_current_num_cycles != 1)
        {
            TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                     << "  attempting to perform "
                                     << d_current_num_cycles
                                     << " cycles of fixed point iteration.\n"
                                     << "  number of cycles required by coupled advection-diffusion solver = "
                                     << adv_diff_num_cycles
                                     << ".\n"
                                     << "  current implementation requires either that both solvers use the same "
                                        "number of cycles,\n"
                                     << "  or that the Navier-Stokes solver use only a single cycle.\n");
        }
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int U_adv_diff_current_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getCurrentContext());
        if (isAllocatedPatchData(U_adv_diff_current_idx))
        {
            d_hier_fc_data_ops->copyData(U_adv_diff_current_idx, d_u_ADV_current_idx);
        }
        d_adv_diff_hier_integrator->preprocessIntegrateHierarchy(current_time, new_time, adv_diff_num_cycles);
        const int U_adv_diff_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getScratchContext());
        if (isAllocatedPatchData(U_adv_diff_scratch_idx))
        {
            d_hier_fc_data_ops->copyData(U_adv_diff_scratch_idx, d_u_ADV_current_idx);
        }
        const int U_adv_diff_new_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getNewContext());
        if (isAllocatedPatchData(U_adv_diff_new_idx))
        {
            d_hier_fc_data_ops->copyData(U_adv_diff_new_idx, d_u_ADV_current_idx);
        }
    }

    // Account for the convective acceleration term.
    TimeSteppingType convective_time_stepping_type = d_convective_time_stepping_type;
    if (getIntegratorStep() == 0 && is_multistep_time_stepping_type(d_convective_time_stepping_type))
    {
        convective_time_stepping_type = d_init_convective_time_stepping_type;
    }
    if (!d_creeping_flow)
    {
        const int U_adv_idx = d_U_adv_vec->getComponentDescriptorIndex(0);
        d_hier_cc_data_ops->copyData(U_adv_idx, d_U_current_idx);
        d_hier_fc_data_ops->copyData(d_u_ADV_scratch_idx, d_u_ADV_current_idx);
        for (int ln = finest_ln; ln > coarsest_ln; --ln)
        {
            Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
            Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
            Pointer<CoarsenOperator<NDIM> > coarsen_op;
            coarsen_op = grid_geom->lookupCoarsenOperator(d_U_var, "CONSERVATIVE_COARSEN");
            coarsen_alg->registerCoarsen(U_adv_idx, U_adv_idx, coarsen_op);
            coarsen_op = grid_geom->lookupCoarsenOperator(d_u_ADV_var, "CONSERVATIVE_COARSEN");
            coarsen_alg->registerCoarsen(d_u_ADV_scratch_idx, d_u_ADV_scratch_idx, coarsen_op);
            coarsen_alg->resetSchedule(getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]);
            getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]->coarsenData();
            getCoarsenAlgorithm(d_object_name + "::CONVECTIVE_OP")
                ->resetSchedule(getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]);
        }
        d_convective_op->setAdvectionVelocity(d_u_ADV_scratch_idx);
        d_convective_op->setSolutionTime(current_time);
        d_convective_op->apply(*d_U_adv_vec, *d_N_vec);
        const int N_idx = d_N_vec->getComponentDescriptorIndex(0);
        d_hier_cc_data_ops->copyData(d_N_old_new_idx, N_idx);
        if (convective_time_stepping_type == FORWARD_EULER)
        {
            d_hier_cc_data_ops->axpy(
                d_U_rhs_vec->getComponentDescriptorIndex(0), -rho, N_idx, d_U_rhs_vec->getComponentDescriptorIndex(0));
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_cc_data_ops->axpy(d_U_rhs_vec->getComponentDescriptorIndex(0),
                                     -0.5 * rho,
                                     N_idx,
                                     d_U_rhs_vec->getComponentDescriptorIndex(0));
        }
    }

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
INSCollocatedHierarchyIntegrator::integrateHierarchy(const double current_time,
                                                     const double new_time,
                                                     const int cycle_num)
{
    INSHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    const double rho = d_problem_coefs.getRho();
    const double mu = d_problem_coefs.getMu();
    const double lambda = d_problem_coefs.getLambda();
    const double volume = d_hier_math_ops->getVolumeOfPhysicalDomain();
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

    // Check to make sure that the number of cycles is what we expect it to be.
    const int expected_num_cycles = getNumberOfCycles();
    if (d_current_num_cycles != expected_num_cycles)
    {
        IBAMR_DO_ONCE(
            {
                pout << "INSCollocatedHierarchyIntegrator::integrateHierarchy():\n"
                     << "  WARNING: num_cycles = " << d_current_num_cycles
                     << " but expected num_cycles = " << expected_num_cycles << ".\n";
            });
    }

    // Perform a single step of fixed point iteration.

    if (d_adv_diff_hier_integrator)
    {
        // Update the state variables maintained by the advection-diffusion
        // solver.
        d_adv_diff_hier_integrator->integrateHierarchy(current_time, new_time, cycle_num);
    }

    // Compute the current approximation to the pressure gradient.
    if (cycle_num == 0)
    {
        switch (d_projection_method_type)
        {
        case PRESSURE_INCREMENT:
            d_hier_cc_data_ops->copyData(d_P_scratch_idx, d_P_current_idx);
            break;
        case PRESSURE_UPDATE:
            d_hier_cc_data_ops->setToScalar(d_P_scratch_idx, 0.0);
            break;
        default:
            TBOX_ERROR("INSCollocatedHierarchyIntegrator::integrateHierarchy():\n"
                       << "  unsupported projection method type: "
                       << enum_to_string<ProjectionMethodType>(d_projection_method_type)
                       << " \n"
                       << "  valid choices are: PRESSURE_INCREMENT, PRESSURE_UPDATE\n");
        }
    }
    else
    {
        d_hier_cc_data_ops->copyData(d_P_scratch_idx, d_P_new_idx);
    }
    d_hier_math_ops->grad(d_Grad_P_idx, d_Grad_P_var, 1.0, d_P_scratch_idx, d_P_var, d_P_bdry_bc_fill_op, half_time);
    d_hier_cc_data_ops->subtract(
        d_U_rhs_vec->getComponentDescriptorIndex(0), d_U_rhs_vec->getComponentDescriptorIndex(0), d_Grad_P_idx);

    // Account for the convective acceleration term.
    TimeSteppingType convective_time_stepping_type = d_convective_time_stepping_type;
    if (is_multistep_time_stepping_type(convective_time_stepping_type))
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(convective_time_stepping_type == ADAMS_BASHFORTH);
#endif
        if (getIntegratorStep() == 0)
        {
            convective_time_stepping_type = d_init_convective_time_stepping_type;
        }
        else if (cycle_num > 0)
        {
            convective_time_stepping_type = MIDPOINT_RULE;
            IBAMR_DO_ONCE(
                {
                    pout << "INSStaggeredHierarchyIntegrator::integrateHierarchy():\n"
                         << "  WARNING: convective_time_stepping_type = "
                         << enum_to_string<TimeSteppingType>(d_convective_time_stepping_type)
                         << " but num_cycles = " << d_current_num_cycles << " > 1.\n"
                         << "           using " << enum_to_string<TimeSteppingType>(d_convective_time_stepping_type)
                         << " only for the first cycle in each time step;\n"
                         << "           using " << enum_to_string<TimeSteppingType>(convective_time_stepping_type)
                         << " for subsequent cycles.\n";
                });
        }
    }
    if (!d_creeping_flow && convective_time_stepping_type != FORWARD_EULER)
    {
        if (cycle_num > 0)
        {
            const int U_adv_idx = d_U_adv_vec->getComponentDescriptorIndex(0);
            double apply_time = std::numeric_limits<double>::quiet_NaN();
            if (convective_time_stepping_type == MIDPOINT_RULE)
            {
                d_hier_cc_data_ops->linearSum(U_adv_idx, 0.5, d_U_current_idx, 0.5, d_U_new_idx);
                d_hier_fc_data_ops->linearSum(d_u_ADV_scratch_idx, 0.5, d_u_ADV_current_idx, 0.5, d_u_ADV_new_idx);
                apply_time = half_time;
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                d_hier_cc_data_ops->copyData(U_adv_idx, d_U_new_idx);
                d_hier_fc_data_ops->copyData(d_u_ADV_scratch_idx, d_u_ADV_new_idx);
                apply_time = new_time;
            }
            for (int ln = finest_ln; ln > coarsest_ln; --ln)
            {
                Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
                Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
                Pointer<CoarsenOperator<NDIM> > coarsen_op;
                coarsen_op = grid_geom->lookupCoarsenOperator(d_U_var, "CONSERVATIVE_COARSEN");
                coarsen_alg->registerCoarsen(U_adv_idx, U_adv_idx, coarsen_op);
                coarsen_op = grid_geom->lookupCoarsenOperator(d_u_ADV_var, "CONSERVATIVE_COARSEN");
                coarsen_alg->registerCoarsen(d_u_ADV_scratch_idx, d_u_ADV_scratch_idx, coarsen_op);
                coarsen_alg->resetSchedule(getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]);
                getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]->coarsenData();
                getCoarsenAlgorithm(d_object_name + "::CONVECTIVE_OP")
                    ->resetSchedule(getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]);
            }
            d_convective_op->setAdvectionVelocity(d_u_ADV_scratch_idx);
            d_convective_op->setSolutionTime(apply_time);
            d_convective_op->apply(*d_U_adv_vec, *d_N_vec);
        }
        const int N_idx = d_N_vec->getComponentDescriptorIndex(0);
        if (convective_time_stepping_type == ADAMS_BASHFORTH)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(cycle_num == 0);
#endif
            const double omega = dt / d_dt_previous[0];
            d_hier_cc_data_ops->linearSum(N_idx, 1.0 + 0.5 * omega, N_idx, -0.5 * omega, d_N_old_current_idx);
        }
        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
        {
            d_hier_cc_data_ops->axpy(
                d_U_rhs_vec->getComponentDescriptorIndex(0), -rho, N_idx, d_U_rhs_vec->getComponentDescriptorIndex(0));
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_cc_data_ops->axpy(d_U_rhs_vec->getComponentDescriptorIndex(0),
                                     -0.5 * rho,
                                     N_idx,
                                     d_U_rhs_vec->getComponentDescriptorIndex(0));
        }
    }

    // Account for body forcing terms.
    if (d_F_fcn)
    {
        d_F_fcn->setDataOnPatchHierarchy(d_F_scratch_idx, d_F_var, d_hierarchy, half_time);
        d_hier_cc_data_ops->add(
            d_U_rhs_vec->getComponentDescriptorIndex(0), d_U_rhs_vec->getComponentDescriptorIndex(0), d_F_scratch_idx);
    }

    // Account for internal source/sink distributions.
    if (d_Q_fcn)
    {
        d_Q_fcn->setDataOnPatchHierarchy(d_Q_current_idx, d_Q_var, d_hierarchy, current_time);
        d_Q_fcn->setDataOnPatchHierarchy(d_Q_new_idx, d_Q_var, d_hierarchy, new_time);
        d_hier_cc_data_ops->linearSum(d_Q_scratch_idx, 0.5, d_Q_current_idx, 0.5, d_Q_new_idx);
        d_Q_bdry_bc_fill_op->fillData(half_time);
        if (!d_creeping_flow)
        {
            d_hier_cc_data_ops->linearSum(d_U_scratch_idx, 0.5, d_U_current_idx, 0.5, d_U_new_idx);
            computeDivSourceTerm(d_F_div_idx, d_Q_scratch_idx, d_U_scratch_idx);
        }
        d_hier_cc_data_ops->axpy(
            d_U_rhs_vec->getComponentDescriptorIndex(0), rho, d_F_div_idx, d_U_rhs_vec->getComponentDescriptorIndex(0));
        d_hier_cc_data_ops->subtract(
            d_Phi_rhs_vec->getComponentDescriptorIndex(0), d_Phi_rhs_vec->getComponentDescriptorIndex(0), d_Q_new_idx);
    }

    // Solve for U(*) and compute u_ADV(*).
    d_hier_cc_data_ops->copyData(d_U_scratch_idx, d_U_new_idx);
    d_velocity_solver->solveSystem(*d_U_scratch_vec, *d_U_rhs_vec);
    if (d_enable_logging)
        plog << d_object_name << "::integrateHierarchy(): velocity solve number of iterations = "
             << d_velocity_solver->getNumIterations() << "\n";
    if (d_enable_logging)
        plog << d_object_name
             << "::integrateHierarchy(): velocity solve residual norm        = " << d_velocity_solver->getResidualNorm()
             << "\n";
    d_hier_math_ops->interp(d_u_ADV_scratch_idx,
                            d_u_ADV_var,
                            /*synch_cf_bdry*/ true,
                            d_U_scratch_idx,
                            d_U_var,
                            d_U_bdry_bc_fill_op,
                            new_time);

    // Project U(*) to compute U(n+1) and u_ADV(n+1).
    const double div_fac =
        (MathUtilities<double>::equalEps(rho, 0.0) || MathUtilities<double>::equalEps(dt, 0.0) ? 1.0 : rho / dt);
    d_hier_math_ops->div(d_Phi_rhs_vec->getComponentDescriptorIndex(0),
                         d_Phi_rhs_vec->getComponentVariable(0),
                         -div_fac,
                         d_u_ADV_scratch_idx,
                         d_u_ADV_var,
                         d_no_fill_op,
                         new_time,
                         /*synch_cf_bdry*/ false,
                         +div_fac,
                         d_Q_new_idx,
                         d_Q_var);
    if (d_normalize_pressure)
    {
        const double Div_U_mean =
            (1.0 / volume) * d_hier_cc_data_ops->integral(d_Phi_rhs_vec->getComponentDescriptorIndex(0), wgt_cc_idx);
        d_hier_cc_data_ops->addScalar(
            d_Phi_rhs_vec->getComponentDescriptorIndex(0), d_Phi_rhs_vec->getComponentDescriptorIndex(0), -Div_U_mean);
    }
    if (cycle_num == 0)
    {
        d_hier_cc_data_ops->subtract(d_Phi_idx, d_P_current_idx, d_P_scratch_idx);
    }
    else
    {
        d_hier_cc_data_ops->setToScalar(d_Phi_idx, 0.0);
    }
    d_pressure_solver->solveSystem(*d_Phi_vec, *d_Phi_rhs_vec);
    if (d_enable_logging)
        plog << d_object_name << "::integrateHierarchy(): pressure solve number of iterations = "
             << d_pressure_solver->getNumIterations() << "\n";
    if (d_enable_logging)
        plog << d_object_name
             << "::integrateHierarchy(): pressure solve residual norm        = " << d_pressure_solver->getResidualNorm()
             << "\n";
    d_hier_math_ops->grad(d_Grad_Phi_fc_idx,
                          d_Grad_Phi_fc_var,
                          /*synch_cf_bdry*/ true,
                          1.0,
                          d_Phi_idx,
                          d_Phi_var,
                          d_Phi_bdry_bc_fill_op,
                          half_time);
    d_hier_fc_data_ops->axpy(d_u_ADV_new_idx, -1.0 / div_fac, d_Grad_Phi_fc_idx, d_u_ADV_scratch_idx);
    d_hier_math_ops->interp(d_Grad_Phi_cc_idx,
                            d_Grad_Phi_cc_var,
                            d_Grad_Phi_fc_idx,
                            d_Grad_Phi_fc_var,
                            d_no_fill_op,
                            half_time,
                            /*synch_cf_bdry*/ false);
    d_hier_cc_data_ops->axpy(d_U_new_idx, -1.0 / div_fac, d_Grad_Phi_cc_idx, d_U_scratch_idx);

    // Determine P(n+1/2).
    double K = 0.0;
    switch (d_viscous_time_stepping_type)
    {
    case BACKWARD_EULER:
        K = 1.0;
        break;
    case FORWARD_EULER:
        K = 0.0;
        break;
    case TRAPEZOIDAL_RULE:
        K = 0.5;
        break;
    default:
        TBOX_ERROR("this statment should not be reached");
    }
    PoissonSpecifications helmholtz_spec(d_object_name + "::helmholtz_spec");
    if (MathUtilities<double>::equalEps(rho, 0.0))
    {
        helmholtz_spec.setCConstant(0.0);
        helmholtz_spec.setDConstant(-K * mu);
    }
    else if (d_using_2nd_order_pressure_update)
    {
        helmholtz_spec.setCConstant(1.0 + K * dt * lambda / rho);
        helmholtz_spec.setDConstant(-K * dt * mu / rho);
    }
    else
    {
        helmholtz_spec.setCConstant(1.0);
        helmholtz_spec.setDConstant(0.0);
    }
    d_hier_math_ops->laplace(d_P_new_idx,
                             d_P_var,
                             helmholtz_spec,
                             d_Phi_idx,
                             d_Phi_var,
                             d_no_fill_op,
                             half_time,
                             1.0,
                             d_P_scratch_idx,
                             d_P_var);
    if (d_normalize_pressure)
    {
        const double P_mean = (1.0 / volume) * d_hier_cc_data_ops->integral(d_P_new_idx, wgt_cc_idx);
        d_hier_cc_data_ops->addScalar(d_P_new_idx, d_P_new_idx, -P_mean);
    }

    // Reset the right-hand side vector.
    d_hier_cc_data_ops->add(
        d_U_rhs_vec->getComponentDescriptorIndex(0), d_U_rhs_vec->getComponentDescriptorIndex(0), d_Grad_P_idx);
    if (!d_creeping_flow)
    {
        const int N_idx = d_N_vec->getComponentDescriptorIndex(0);
        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
        {
            d_hier_cc_data_ops->axpy(
                d_U_rhs_vec->getComponentDescriptorIndex(0), +rho, N_idx, d_U_rhs_vec->getComponentDescriptorIndex(0));
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_cc_data_ops->axpy(d_U_rhs_vec->getComponentDescriptorIndex(0),
                                     +0.5 * rho,
                                     N_idx,
                                     d_U_rhs_vec->getComponentDescriptorIndex(0));
        }
    }
    if (d_F_fcn)
    {
        d_hier_cc_data_ops->subtract(
            d_U_rhs_vec->getComponentDescriptorIndex(0), d_U_rhs_vec->getComponentDescriptorIndex(0), d_F_scratch_idx);
        d_hier_cc_data_ops->copyData(d_F_new_idx, d_F_scratch_idx);
    }
    if (d_Q_fcn)
    {
        d_hier_cc_data_ops->axpy(d_U_rhs_vec->getComponentDescriptorIndex(0),
                                 -rho,
                                 d_F_div_idx,
                                 d_U_rhs_vec->getComponentDescriptorIndex(0));
        d_hier_cc_data_ops->add(
            d_Phi_rhs_vec->getComponentDescriptorIndex(0), d_Phi_rhs_vec->getComponentDescriptorIndex(0), d_Q_new_idx);
    }

    if (d_adv_diff_hier_integrator)
    {
        // Update the advection velocities used by the advection-diffusion
        // solver.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int U_adv_diff_new_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getNewContext());
        if (isAllocatedPatchData(U_adv_diff_new_idx))
        {
            d_hier_fc_data_ops->copyData(U_adv_diff_new_idx, d_u_ADV_new_idx);
        }
        const int U_adv_diff_current_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getCurrentContext());
        const int U_adv_diff_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getScratchContext());
        if (isAllocatedPatchData(U_adv_diff_scratch_idx))
        {
            d_hier_fc_data_ops->linearSum(U_adv_diff_scratch_idx, 0.5, U_adv_diff_current_idx, 0.5, U_adv_diff_new_idx);
        }

        // Update the state variables maintained by the advection-diffusion
        // solver.
        //
        // NOTE: We already performed cycle 0 above.
        const int adv_diff_num_cycles = d_adv_diff_hier_integrator->getNumberOfCycles();
        if (d_current_num_cycles != adv_diff_num_cycles)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(d_current_num_cycles == 1);
#endif
            for (int adv_diff_cycle_num = 1; adv_diff_cycle_num < adv_diff_num_cycles; ++adv_diff_cycle_num)
            {
                d_adv_diff_hier_integrator->integrateHierarchy(current_time, new_time, adv_diff_cycle_num);
            }
        }
    }

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
    return;
} // integrateHierarchy

void
INSCollocatedHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                const double new_time,
                                                                const bool skip_synchronize_new_state_data,
                                                                const int num_cycles)
{
    INSHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Synchronize new state data.
    if (!skip_synchronize_new_state_data)
    {
        if (d_enable_logging)
            plog << d_object_name << "::postprocessIntegrateHierarchy(): synchronizing updated data\n";
        synchronizeHierarchyData(NEW_DATA);
    }

    // Determine the CFL number.
    if (!d_parent_integrator)
    {
        double cfl_max = 0.0;
        PatchCellDataOpsReal<NDIM, double> patch_cc_ops;
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double dx_min = *(std::min_element(dx, dx + NDIM));
                Pointer<CellData<NDIM, double> > u_cc_new_data = patch->getPatchData(d_U_new_idx);
                double u_max = 0.0;
                u_max = patch_cc_ops.maxNorm(u_cc_new_data, patch_box);
                cfl_max = std::max(cfl_max, u_max * dt / dx_min);
            }
        }
        cfl_max = SAMRAI_MPI::maxReduction(cfl_max);
        if (d_enable_logging)
            plog << d_object_name << "::postprocessIntegrateHierarchy(): CFL number = " << cfl_max << "\n";
    }

    // Compute max |Omega|_2.
    if (d_using_vorticity_tagging)
    {
        d_hier_cc_data_ops->copyData(d_U_scratch_idx, d_U_new_idx);
        d_hier_math_ops->curl(d_Omega_idx, d_Omega_var, d_U_scratch_idx, d_U_var, d_U_bdry_bc_fill_op, new_time);
#if (NDIM == 3)
        d_hier_math_ops->pointwiseL2Norm(d_Omega_Norm_idx, d_Omega_Norm_var, d_Omega_idx, d_Omega_var);
#endif
        const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
#if (NDIM == 2)
        d_Omega_max = d_hier_cc_data_ops->maxNorm(d_Omega_idx, wgt_cc_idx);
#endif
#if (NDIM == 3)
        d_Omega_max = d_hier_cc_data_ops->max(d_Omega_Norm_idx, wgt_cc_idx);
#endif
    }

    // Deallocate scratch data.
    d_U_rhs_vec->deallocateVectorData();
    d_Phi_rhs_vec->deallocateVectorData();
    if (!d_creeping_flow)
    {
        d_U_adv_vec->deallocateVectorData();
        d_N_vec->deallocateVectorData();
    }

    // Deallocate any registered advection-diffusion solver.
    if (d_adv_diff_hier_integrator)
    {
        const int adv_diff_num_cycles = d_adv_diff_hier_integrator->getNumberOfCycles();
        d_adv_diff_hier_integrator->postprocessIntegrateHierarchy(
            current_time, new_time, skip_synchronize_new_state_data, adv_diff_num_cycles);
    }

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

void
INSCollocatedHierarchyIntegrator::regridHierarchy()
{
    const int coarsest_ln = 0;

    // Regrid the hierarchy.
    switch (d_regrid_mode)
    {
    case STANDARD:
        d_gridding_alg->regridAllFinerLevels(d_hierarchy, coarsest_ln, d_integrator_time, d_tag_buffer);
        break;
    case AGGRESSIVE:
        for (int k = 0; k < d_gridding_alg->getMaxLevels(); ++k)
        {
            d_gridding_alg->regridAllFinerLevels(d_hierarchy, coarsest_ln, d_integrator_time, d_tag_buffer);
        }
        break;
    default:
        TBOX_ERROR(d_object_name << "::regridHierarchy():\n"
                                 << "  unrecognized regrid mode: "
                                 << IBTK::enum_to_string<RegridMode>(d_regrid_mode)
                                 << "."
                                 << std::endl);
    }

    // Project the interpolated velocity.
    regridProjection();

    // Synchronize the state data on the patch hierarchy.
    synchronizeHierarchyData(CURRENT_DATA);
    return;
} // regridHierarchy

/////////////////////////////// PROTECTED ////////////////////////////////////

void
INSCollocatedHierarchyIntegrator::initializeLevelDataSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int level_number,
    const double init_data_time,
    const bool /*can_be_refined*/,
    const bool initial_time,
    const Pointer<BasePatchLevel<NDIM> > base_old_level,
    const bool /*allocate_data*/)
{
    const Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    const Pointer<PatchLevel<NDIM> > old_level = base_old_level;
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    if (old_level)
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
#endif
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Initialize level data at the initial time.
    if (initial_time)
    {
        // Allocate scratch data.
        for (int ln = 0; ln <= level_number; ++ln)
        {
            hierarchy->getPatchLevel(ln)->allocatePatchData(d_U_scratch_idx, init_data_time);
#if (NDIM == 3)
            hierarchy->getPatchLevel(ln)->allocatePatchData(d_Omega_Norm_idx, init_data_time);
#endif
        }

        // Fill ghost cells.
        HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
        Pointer<HierarchyCellDataOpsReal<NDIM, double> > hier_cc_data_ops =
            hier_ops_manager->getOperationsDouble(d_U_var, d_hierarchy, true);
        hier_cc_data_ops->resetLevels(0, level_number);
        hier_cc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        InterpolationTransactionComponent U_bc_component(d_U_scratch_idx,
                                                         DATA_REFINE_TYPE,
                                                         USE_CF_INTERPOLATION,
                                                         DATA_COARSEN_TYPE,
                                                         d_bdry_extrap_type,
                                                         CONSISTENT_TYPE_2_BDRY,
                                                         d_U_star_bc_coefs);
        HierarchyGhostCellInterpolation U_bdry_bc_fill_op;
        U_bdry_bc_fill_op.initializeOperatorState(U_bc_component, d_hierarchy, 0, level_number);
        U_bdry_bc_fill_op.fillData(init_data_time);
        HierarchyMathOps hier_math_ops(d_object_name + "::HierarchyLevelMathOps", d_hierarchy, 0, level_number);

        // Interpolate U to u_ADV.
        hier_math_ops.interp(d_u_ADV_current_idx,
                             d_u_ADV_var,
                             /*synch_cf_bdry*/ false,
                             d_U_scratch_idx,
                             d_U_var,
                             d_no_fill_op,
                             init_data_time);

        // Initialize the maximum value of |Omega|_2 on the grid.
        if (d_using_vorticity_tagging)
        {
            if (level_number == 0) d_Omega_max = 0.0;

            // Compute max |Omega|_2.
            hier_math_ops.curl(d_Omega_idx, d_Omega_var, d_U_scratch_idx, d_U_var, d_U_bdry_bc_fill_op, init_data_time);
#if (NDIM == 3)
            hier_math_ops.pointwiseL2Norm(d_Omega_Norm_idx, d_Omega_Norm_var, d_Omega_idx, d_Omega_var);
#endif
            const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
#if (NDIM == 2)
            d_Omega_max = hier_cc_data_ops->maxNorm(d_Omega_idx, wgt_cc_idx);
#endif
#if (NDIM == 3)
            d_Omega_max = hier_cc_data_ops->max(d_Omega_Norm_idx, wgt_cc_idx);
#endif
        }

        // Deallocate scratch data.
        for (int ln = 0; ln <= level_number; ++ln)
        {
            hierarchy->getPatchLevel(ln)->deallocatePatchData(d_U_scratch_idx);
#if (NDIM == 3)
            hierarchy->getPatchLevel(ln)->deallocatePatchData(d_Omega_Norm_idx);
#endif
        }
    }
    return;
} // initializeLevelDataSpecialized

void
INSCollocatedHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    const Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((coarsest_level >= 0) && (coarsest_level <= finest_level) &&
                (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(hierarchy->getPatchLevel(ln));
    }
#else
    NULL_USE(coarsest_level);
    NULL_USE(finest_level);
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // Reset the hierarchy operations objects for the new hierarchy configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);

    d_hier_fc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_fc_data_ops->resetLevels(0, finest_hier_level);

    // Setup the patch boundary filling objects.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent U_bc_component(d_U_scratch_idx,
                                                     DATA_REFINE_TYPE,
                                                     USE_CF_INTERPOLATION,
                                                     DATA_COARSEN_TYPE,
                                                     d_bdry_extrap_type,
                                                     CONSISTENT_TYPE_2_BDRY,
                                                     d_U_star_bc_coefs);
    d_U_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_U_bdry_bc_fill_op->initializeOperatorState(U_bc_component, d_hierarchy);

    InterpolationTransactionComponent P_bc_component(d_P_scratch_idx,
                                                     DATA_REFINE_TYPE,
                                                     USE_CF_INTERPOLATION,
                                                     DATA_COARSEN_TYPE,
                                                     d_bdry_extrap_type,
                                                     CONSISTENT_TYPE_2_BDRY);
    d_P_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_P_bdry_bc_fill_op->initializeOperatorState(P_bc_component, d_hierarchy);

    InterpolationTransactionComponent Phi_bc_component(d_Phi_idx,
                                                       DATA_REFINE_TYPE,
                                                       USE_CF_INTERPOLATION,
                                                       DATA_COARSEN_TYPE,
                                                       d_bdry_extrap_type,
                                                       CONSISTENT_TYPE_2_BDRY,
                                                       d_Phi_bc_coef);
    d_Phi_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_Phi_bdry_bc_fill_op->initializeOperatorState(Phi_bc_component, d_hierarchy);

    if (d_Q_fcn)
    {
        InterpolationTransactionComponent Q_bc_component(d_Q_scratch_idx,
                                                         DATA_REFINE_TYPE,
                                                         USE_CF_INTERPOLATION,
                                                         DATA_COARSEN_TYPE,
                                                         d_bdry_extrap_type,
                                                         CONSISTENT_TYPE_2_BDRY);
        d_Q_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
        d_Q_bdry_bc_fill_op->initializeOperatorState(Q_bc_component, d_hierarchy);
    }

    // Indicate that vectors and solvers need to be re-initialized.
    d_coarsest_reset_ln = coarsest_level;
    d_finest_reset_ln = finest_level;
    d_vectors_need_init = true;
    d_convective_op_needs_init = true;
    d_velocity_solver_needs_init = true;
    d_pressure_solver_needs_init = true;
    return;
} // resetHierarchyConfigurationSpecialized

void
INSCollocatedHierarchyIntegrator::applyGradientDetectorSpecialized(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                                                   const int level_number,
                                                                   const double /*error_data_time*/,
                                                                   const int tag_index,
                                                                   const bool /*initial_time*/,
                                                                   const bool /*uses_richardson_extrapolation_too*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
#endif
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Tag cells based on the magnitude of the vorticity.
    //
    // Note that if either the relative or absolute threshold is zero for a
    // particular level, no tagging is performed on that level.
    if (d_using_vorticity_tagging)
    {
        double Omega_rel_thresh = 0.0;
        if (d_Omega_rel_thresh.size() > 0)
        {
            Omega_rel_thresh = d_Omega_rel_thresh[std::max(std::min(level_number, d_Omega_rel_thresh.size() - 1), 0)];
        }
        double Omega_abs_thresh = 0.0;
        if (d_Omega_abs_thresh.size() > 0)
        {
            Omega_abs_thresh = d_Omega_abs_thresh[std::max(std::min(level_number, d_Omega_abs_thresh.size() - 1), 0)];
        }
        if (Omega_rel_thresh > 0.0 || Omega_abs_thresh > 0.0)
        {
            double thresh = std::numeric_limits<double>::max();
            if (Omega_rel_thresh > 0.0) thresh = std::min(thresh, Omega_rel_thresh * d_Omega_max);
            if (Omega_abs_thresh > 0.0) thresh = std::min(thresh, Omega_abs_thresh);
            thresh += sqrt(std::numeric_limits<double>::epsilon());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, int> > tags_data = patch->getPatchData(tag_index);
                Pointer<CellData<NDIM, double> > Omega_data = patch->getPatchData(d_Omega_idx);
                for (CellIterator<NDIM> ic(patch_box); ic; ic++)
                {
                    const Index<NDIM>& i = ic();
#if (NDIM == 2)
                    if (std::abs((*Omega_data)(i)) > thresh)
                    {
                        (*tags_data)(i) = 1;
                    }
#endif
#if (NDIM == 3)
                    double norm_Omega_sq = 0.0;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        norm_Omega_sq += (*Omega_data)(i, d) * (*Omega_data)(i, d);
                    }
                    const double norm_Omega = sqrt(norm_Omega_sq);
                    if (norm_Omega > thresh)
                    {
                        (*tags_data)(i) = 1;
                    }
#endif
                }
            }
        }
    }
    return;
} // applyGradientDetectorSpecialized

void
INSCollocatedHierarchyIntegrator::setupPlotDataSpecialized()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Allocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_U_scratch_idx, d_integrator_time);
    }

    // Fill U ghost cells.
    d_hier_cc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
    d_U_bdry_bc_fill_op->fillData(d_integrator_time);

    // Compute Omega = curl U.
    if (d_output_Omega)
        d_hier_math_ops->curl(d_Omega_idx, d_Omega_var, d_U_scratch_idx, d_U_var, d_no_fill_op, d_integrator_time);

    // Compute Div U.
    if (d_output_Div_U)
    {
        d_hier_math_ops->div(d_Div_U_idx, d_Div_U_var, 1.0, d_U_scratch_idx, d_U_var, d_no_fill_op, d_integrator_time);
        d_hier_math_ops->div(d_Div_u_ADV_idx,
                             d_Div_u_ADV_var,
                             +1.0,
                             d_u_ADV_current_idx,
                             d_u_ADV_var,
                             d_no_fill_op,
                             d_integrator_time,
                             /*synch_cf_bdry*/ false);
    }

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_U_scratch_idx);
    }
    synchronizeHierarchyData(CURRENT_DATA);
    return;
} // setupPlotDataSpecialized

void
INSCollocatedHierarchyIntegrator::regridProjection()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const double volume = d_hier_math_ops->getVolumeOfPhysicalDomain();

    // Setup the solver vectors.
    SAMRAIVectorReal<NDIM, double> sol_vec(d_object_name + "::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec.addComponent(d_Phi_var, d_Phi_idx, wgt_cc_idx, d_hier_cc_data_ops);

    SAMRAIVectorReal<NDIM, double> rhs_vec(d_object_name + "::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    rhs_vec.addComponent(d_Div_U_var, d_Div_U_idx, wgt_cc_idx, d_hier_cc_data_ops);

    // Setup the regrid Poisson solver.
    Pointer<PoissonSolver> regrid_projection_solver =
        CCPoissonSolverManager::getManager()->allocateSolver(d_regrid_projection_solver_type,
                                                             d_object_name + "::regrid_projection_solver",
                                                             d_regrid_projection_solver_db,
                                                             "regrid_projection_",
                                                             d_regrid_projection_precond_type,
                                                             d_object_name + "::regrid_projection_precond",
                                                             d_regrid_projection_precond_db,
                                                             "regrid_projection_pc_");
    PoissonSpecifications regrid_projection_spec(d_object_name + "::regrid_projection_spec");
    regrid_projection_spec.setCZero();
    regrid_projection_spec.setDConstant(-1.0);
    LocationIndexRobinBcCoefs<NDIM> Phi_bc_coef;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        Phi_bc_coef.setBoundarySlope(2 * d, 0.0);
        Phi_bc_coef.setBoundarySlope(2 * d + 1, 0.0);
    }
    regrid_projection_solver->setPoissonSpecifications(regrid_projection_spec);
    regrid_projection_solver->setPhysicalBcCoef(&Phi_bc_coef);
    regrid_projection_solver->setHomogeneousBc(true);
    regrid_projection_solver->setSolutionTime(d_integrator_time);
    regrid_projection_solver->setTimeInterval(d_integrator_time, d_integrator_time);
    LinearSolver* p_regrid_projection_solver = dynamic_cast<LinearSolver*>(regrid_projection_solver.getPointer());
    if (p_regrid_projection_solver)
    {
        p_regrid_projection_solver->setInitialGuessNonzero(false);
        p_regrid_projection_solver->setNullspace(true);
    }

    // Allocate temporary data.
    ComponentSelector scratch_idxs;
    scratch_idxs.setFlag(d_U_scratch_idx);
    scratch_idxs.setFlag(d_Phi_idx);
    scratch_idxs.setFlag(d_Grad_Phi_cc_idx);
    scratch_idxs.setFlag(d_Grad_Phi_fc_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(scratch_idxs, d_integrator_time);
    }

    // Interpolate U to u_ADV.
    d_hier_cc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
    d_hier_math_ops->interp(d_u_ADV_current_idx,
                            d_u_ADV_var,
                            /*synch_cf_bdry*/ true,
                            d_U_scratch_idx,
                            d_U_var,
                            d_U_bdry_bc_fill_op,
                            d_integrator_time);

    // Setup the right-hand side vector for the projection-Poisson solve.
    d_hier_math_ops->div(d_Div_U_idx,
                         d_Div_U_var,
                         -1.0,
                         d_u_ADV_current_idx,
                         d_u_ADV_var,
                         d_no_fill_op,
                         d_integrator_time,
                         /*synch_cf_bdry*/ false,
                         +1.0,
                         d_Q_current_idx,
                         d_Q_var);
    const double Div_U_mean = (1.0 / volume) * d_hier_cc_data_ops->integral(d_Div_U_idx, wgt_cc_idx);
    d_hier_cc_data_ops->addScalar(d_Div_U_idx, d_Div_U_idx, -Div_U_mean);

    // Solve the projection pressure-Poisson problem.
    regrid_projection_solver->solveSystem(sol_vec, rhs_vec);
    if (d_enable_logging)
        plog << d_object_name << "::regridProjection(): projection solve number of iterations = "
             << regrid_projection_solver->getNumIterations() << "\n";
    if (d_enable_logging)
        plog << d_object_name << "::regridProjection(): projection solve residual norm        = "
             << regrid_projection_solver->getResidualNorm() << "\n";

    // Fill ghost cells for Phi, compute Grad Phi, and set U := U - Grad Phi.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent Phi_bc_component(d_Phi_idx,
                                                       DATA_REFINE_TYPE,
                                                       USE_CF_INTERPOLATION,
                                                       DATA_COARSEN_TYPE,
                                                       d_bdry_extrap_type,
                                                       CONSISTENT_TYPE_2_BDRY,
                                                       &Phi_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> Phi_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    Phi_bdry_bc_fill_op->initializeOperatorState(Phi_bc_component, d_hierarchy);
    Phi_bdry_bc_fill_op->setHomogeneousBc(true);
    Phi_bdry_bc_fill_op->fillData(d_integrator_time);
    d_hier_math_ops->grad(d_Grad_Phi_fc_idx,
                          d_Grad_Phi_fc_var,
                          /*synch_cf_bdry*/ true,
                          1.0,
                          d_Phi_idx,
                          d_Phi_var,
                          d_no_fill_op,
                          d_integrator_time);
    d_hier_fc_data_ops->subtract(d_u_ADV_current_idx, d_u_ADV_current_idx, d_Grad_Phi_fc_idx);
    d_hier_math_ops->interp(d_Grad_Phi_cc_idx,
                            d_Grad_Phi_cc_var,
                            d_Grad_Phi_fc_idx,
                            d_Grad_Phi_fc_var,
                            d_no_fill_op,
                            d_integrator_time,
                            /*synch_cf_bdry*/ false);
    d_hier_cc_data_ops->subtract(d_U_current_idx, d_U_current_idx, d_Grad_Phi_cc_idx);

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(scratch_idxs);
    }
    return;
} // regridProjection

double
INSCollocatedHierarchyIntegrator::getStableTimestep(Pointer<Patch<NDIM> > patch) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const Index<NDIM>& ilower = patch->getBox().lower();
    const Index<NDIM>& iupper = patch->getBox().upper();

    Pointer<FaceData<NDIM, double> > u_ADV_data = patch->getPatchData(d_u_ADV_var, getCurrentContext());
    const IntVector<NDIM>& u_ADV_ghost_cells = u_ADV_data->getGhostCellWidth();

    double stable_dt = std::numeric_limits<double>::max();
    ADVECT_STABLEDT_FC(dx,
#if (NDIM == 2)
                       ilower(0),
                       iupper(0),
                       ilower(1),
                       iupper(1),
                       u_ADV_ghost_cells(0),
                       u_ADV_ghost_cells(1),
                       u_ADV_data->getPointer(0),
                       u_ADV_data->getPointer(1),
#endif
#if (NDIM == 3)
                       ilower(0),
                       iupper(0),
                       ilower(1),
                       iupper(1),
                       ilower(2),
                       iupper(2),
                       u_ADV_ghost_cells(0),
                       u_ADV_ghost_cells(1),
                       u_ADV_ghost_cells(2),
                       u_ADV_data->getPointer(0),
                       u_ADV_data->getPointer(1),
                       u_ADV_data->getPointer(2),
#endif
                       stable_dt);
    return stable_dt;
} // getStableTimestep

/////////////////////////////// PRIVATE //////////////////////////////////////

void
INSCollocatedHierarchyIntegrator::reinitializeOperatorsAndSolvers(const double current_time, const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const double rho = d_problem_coefs.getRho();
    const double mu = d_problem_coefs.getMu();
    const double lambda = d_problem_coefs.getLambda();
    double K = 0.0;
    switch (d_viscous_time_stepping_type)
    {
    case BACKWARD_EULER:
        K = 1.0;
        break;
    case FORWARD_EULER:
        K = 0.0;
        break;
    case TRAPEZOIDAL_RULE:
        K = 0.5;
        break;
    default:
        TBOX_ERROR("this statment should not be reached");
    }
    PoissonSpecifications U_problem_coefs(d_object_name + "::U_problem_coefs");
    U_problem_coefs.setCConstant((rho / dt) + K * lambda);
    U_problem_coefs.setDConstant(-K * mu);
    PoissonSpecifications P_problem_coefs(d_object_name + "::P_problem_coefs");
    P_problem_coefs.setCZero();
    P_problem_coefs.setDConstant(-1.0);

    // Ensure that solver components are appropriately reinitialized when the
    // time step size changes.
    const bool dt_change = initial_time || !MathUtilities<double>::equalEps(dt, d_dt_previous[0]);
    if (dt_change)
    {
        d_velocity_solver_needs_init = true;
    }

    // Setup solver vectors.
    const bool has_velocity_nullspace = d_normalize_velocity && MathUtilities<double>::equalEps(rho, 0.0);
    const bool has_pressure_nullspace = d_normalize_pressure;
    if (d_vectors_need_init)
    {
        d_U_scratch_vec =
            new SAMRAIVectorReal<NDIM, double>(d_object_name + "::U_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_U_scratch_vec->addComponent(d_U_var, d_U_scratch_idx, wgt_cc_idx, d_hier_cc_data_ops);

        d_Phi_vec =
            new SAMRAIVectorReal<NDIM, double>(d_object_name + "::Phi_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_Phi_vec->addComponent(d_Phi_var, d_Phi_idx, wgt_cc_idx, d_hier_cc_data_ops);

        if (d_U_rhs_vec) d_U_rhs_vec->freeVectorComponents();
        if (d_U_adv_vec) d_U_adv_vec->freeVectorComponents();
        if (d_N_vec) d_N_vec->freeVectorComponents();
        if (d_Phi_rhs_vec) d_Phi_rhs_vec->freeVectorComponents();

        d_U_rhs_vec = d_U_scratch_vec->cloneVector(d_object_name + "::U_rhs_vec");
        d_U_adv_vec = d_U_scratch_vec->cloneVector(d_object_name + "::U_adv_vec");
        d_N_vec = d_U_scratch_vec->cloneVector(d_object_name + "::N_vec");
        d_Phi_rhs_vec = d_Phi_vec->cloneVector(d_object_name + "::Phi_rhs_vec");

        for (unsigned int k = 0; k < d_U_nul_vecs.size(); ++k)
        {
            if (d_U_nul_vecs[k]) d_U_nul_vecs[k]->freeVectorComponents();
        }
        const int n_U_nul_vecs = (has_velocity_nullspace ? NDIM : 0);
        d_U_nul_vecs.resize(n_U_nul_vecs);

        if (has_velocity_nullspace)
        {
            for (unsigned int k = 0; k < NDIM; ++k)
            {
                std::ostringstream stream;
                stream << k;
                d_U_nul_vecs[k] = d_U_scratch_vec->cloneVector(d_object_name + "::U_nul_vec_U_" + stream.str());
                d_U_nul_vecs[k]->allocateVectorData(current_time);
                for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
                {
                    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        Pointer<CellData<NDIM, double> > U_nul_data =
                            patch->getPatchData(d_U_nul_vecs[k]->getComponentDescriptorIndex(0));
                        U_nul_data->fillAll(0.0);
                        U_nul_data->fill(1.0, k);
                    }
                }
            }
        }

        d_vectors_need_init = false;
    }

    // Setup boundary conditions objects.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        INSIntermediateVelocityBcCoef* U_star_bc_coef =
            dynamic_cast<INSIntermediateVelocityBcCoef*>(d_U_star_bc_coefs[d]);
        U_star_bc_coef->setPhysicalBcCoefs(d_bc_coefs);
        U_star_bc_coef->setSolutionTime(new_time);
        U_star_bc_coef->setTimeInterval(current_time, new_time);
    }
    INSProjectionBcCoef* Phi_bc_coef = dynamic_cast<INSProjectionBcCoef*>(d_Phi_bc_coef);
    Phi_bc_coef->setPhysicalBcCoefs(d_bc_coefs);
    Phi_bc_coef->setSolutionTime(0.5 * (current_time + new_time));
    Phi_bc_coef->setTimeInterval(current_time, new_time);

    // Setup convective operator.
    if (d_convective_op && d_convective_op_needs_init)
    {
        if (d_enable_logging)
            plog << d_object_name << "::preprocessIntegrateHierarchy(): initializing convective operator" << std::endl;
        d_convective_op->setAdvectionVelocity(d_U_scratch_idx);
        d_convective_op->initializeOperatorState(*d_U_scratch_vec, *d_U_rhs_vec);
        d_convective_op_needs_init = false;
    }

    // Setup subdomain solvers.
    if (d_velocity_solver)
    {
        d_velocity_solver->setPoissonSpecifications(U_problem_coefs);
        d_velocity_solver->setPhysicalBcCoefs(d_U_star_bc_coefs);
        d_velocity_solver->setSolutionTime(new_time);
        d_velocity_solver->setTimeInterval(current_time, new_time);
        LinearSolver* p_velocity_solver = dynamic_cast<LinearSolver*>(d_velocity_solver.getPointer());
        if (d_velocity_solver_needs_init)
        {
            if (d_enable_logging)
                plog << d_object_name << "::preprocessIntegrateHierarchy(): initializing "
                                         "velocity subdomain solver"
                     << std::endl;
            if (p_velocity_solver)
            {
                p_velocity_solver->setInitialGuessNonzero(true);
                if (has_velocity_nullspace) p_velocity_solver->setNullspace(false, d_U_nul_vecs);
            }
            d_velocity_solver->initializeSolverState(*d_U_scratch_vec, *d_U_rhs_vec);
            d_velocity_solver_needs_init = false;
        }
    }

    if (d_pressure_solver)
    {
        d_pressure_solver->setPoissonSpecifications(P_problem_coefs);
        d_pressure_solver->setPhysicalBcCoef(d_Phi_bc_coef);
        d_pressure_solver->setSolutionTime(half_time);
        d_pressure_solver->setTimeInterval(current_time, new_time);
        LinearSolver* p_pressure_solver = dynamic_cast<LinearSolver*>(d_pressure_solver.getPointer());
        if (d_pressure_solver_needs_init)
        {
            if (d_enable_logging)
                plog << d_object_name << "::preprocessIntegrateHierarchy(): initializing "
                                         "pressure subdomain solver"
                     << std::endl;
            if (p_pressure_solver)
            {
                p_pressure_solver->setInitialGuessNonzero(true);
                if (has_pressure_nullspace) p_pressure_solver->setNullspace(true);
            }
            d_pressure_solver->initializeSolverState(*d_Phi_vec, *d_Phi_rhs_vec);
            d_pressure_solver_needs_init = false;
        }
    }
    return;
} // reinitializeOperatorsAndSolvers

void
INSCollocatedHierarchyIntegrator::computeDivSourceTerm(const int F_idx, const int Q_idx, const int u_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            const Index<NDIM>& ilower = patch->getBox().lower();
            const Index<NDIM>& iupper = patch->getBox().upper();

            Pointer<FaceData<NDIM, double> > u_data = patch->getPatchData(u_idx);
            Pointer<CellData<NDIM, double> > Q_data = patch->getPatchData(Q_idx);
            Pointer<CellData<NDIM, double> > F_data = patch->getPatchData(F_idx);

            const IntVector<NDIM>& u_data_gc = u_data->getGhostCellWidth();
            const IntVector<NDIM>& Q_data_gc = Q_data->getGhostCellWidth();
            const IntVector<NDIM>& F_data_gc = F_data->getGhostCellWidth();

            switch (d_convective_op->getConvectiveDifferencingType())
            {
            case CONSERVATIVE:
                NAVIER_STOKES_CONS_SOURCE_FC(
#if (NDIM == 2)
                    ilower(0),
                    iupper(0),
                    ilower(1),
                    iupper(1),
                    u_data_gc(0),
                    u_data_gc(1),
                    Q_data_gc(0),
                    Q_data_gc(1),
                    F_data_gc(0),
                    F_data_gc(1),
                    u_data->getPointer(0),
                    u_data->getPointer(1),
#endif
#if (NDIM == 3)
                    ilower(0),
                    iupper(0),
                    ilower(1),
                    iupper(1),
                    ilower(2),
                    iupper(2),
                    u_data_gc(0),
                    u_data_gc(1),
                    u_data_gc(2),
                    Q_data_gc(0),
                    Q_data_gc(1),
                    Q_data_gc(2),
                    F_data_gc(0),
                    F_data_gc(1),
                    F_data_gc(2),
                    u_data->getPointer(0),
                    u_data->getPointer(1),
                    u_data->getPointer(2),
#endif
                    Q_data->getPointer(),
                    F_data->getPointer());
                break;
            case ADVECTIVE:
                NAVIER_STOKES_ADV_SOURCE_FC(
#if (NDIM == 2)
                    ilower(0),
                    iupper(0),
                    ilower(1),
                    iupper(1),
                    u_data_gc(0),
                    u_data_gc(1),
                    Q_data_gc(0),
                    Q_data_gc(1),
                    F_data_gc(0),
                    F_data_gc(1),
                    u_data->getPointer(0),
                    u_data->getPointer(1),
#endif
#if (NDIM == 3)
                    ilower(0),
                    iupper(0),
                    ilower(1),
                    iupper(1),
                    ilower(2),
                    iupper(2),
                    u_data_gc(0),
                    u_data_gc(1),
                    u_data_gc(2),
                    Q_data_gc(0),
                    Q_data_gc(1),
                    Q_data_gc(2),
                    F_data_gc(0),
                    F_data_gc(1),
                    F_data_gc(2),
                    u_data->getPointer(0),
                    u_data->getPointer(1),
                    u_data->getPointer(2),
#endif
                    Q_data->getPointer(),
                    F_data->getPointer());
                break;
            case SKEW_SYMMETRIC:
                NAVIER_STOKES_SKEW_SYM_SOURCE_FC(
#if (NDIM == 2)
                    ilower(0),
                    iupper(0),
                    ilower(1),
                    iupper(1),
                    u_data_gc(0),
                    u_data_gc(1),
                    Q_data_gc(0),
                    Q_data_gc(1),
                    F_data_gc(0),
                    F_data_gc(1),
                    u_data->getPointer(0),
                    u_data->getPointer(1),
#endif
#if (NDIM == 3)
                    ilower(0),
                    iupper(0),
                    ilower(1),
                    iupper(1),
                    ilower(2),
                    iupper(2),
                    u_data_gc(0),
                    u_data_gc(1),
                    u_data_gc(2),
                    Q_data_gc(0),
                    Q_data_gc(1),
                    Q_data_gc(2),
                    F_data_gc(0),
                    F_data_gc(1),
                    F_data_gc(2),
                    u_data->getPointer(0),
                    u_data->getPointer(1),
                    u_data->getPointer(2),
#endif
                    Q_data->getPointer(),
                    F_data->getPointer());
                break;
            default:
                TBOX_ERROR(
                    "INSCollocatedHierarchyIntegrator::computeDivSourceTerm():\n"
                    << "  unsupported differencing form: "
                    << enum_to_string<ConvectiveDifferencingType>(d_convective_op->getConvectiveDifferencingType())
                    << " \n"
                    << "  valid choices are: ADVECTIVE, CONSERVATIVE, SKEW_SYMMETRIC\n");
            }
        }
    }
    return;
} // computeDivSourceTerm

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
