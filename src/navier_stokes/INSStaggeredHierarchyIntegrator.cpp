// Filename: INSStaggeredHierarchyIntegrator.cpp
// Created on 20 Mar 2008 by Boyce Griffith
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

#include "ArrayData.h"
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
#include "HierarchySideDataOpsReal.h"
#include "IBAMR_config.h"
#include "Index.h"
#include "IntVector.h"
#include "LocationIndexRobinBcCoefs.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchSideDataOpsReal.h"
#include "PoissonSpecifications.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefinePatchStrategy.h"
#include "RefineSchedule.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "VisItDataWriter.h"
#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSIntermediateVelocityBcCoef.h"
#include "ibamr/INSProjectionBcCoef.h"
#include "ibamr/INSStaggeredConvectiveOperatorManager.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/INSStaggeredPressureBcCoef.h"
#include "ibamr/INSStaggeredVelocityBcCoef.h"
#include "ibamr/StaggeredStokesBlockPreconditioner.h"
#include "ibamr/StaggeredStokesFACPreconditioner.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/CartSideDoubleDivPreservingRefine.h"
#include "ibtk/CartSideDoubleSpecializedConstantRefine.h"
#include "ibtk/CartSideDoubleSpecializedLinearRefine.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/NewtonKrylovSolver.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/SCPoissonSolverManager.h"
#include "ibtk/SideDataSynchronization.h"
#include "ibtk/ibtk_enums.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define NAVIER_STOKES_SC_STABLEDT_FC IBAMR_FC_FUNC_(navier_stokes_sc_stabledt2d, NAVIER_STOKES_SC_STABLEDT2D)
#define NAVIER_STOKES_SIDE_TO_FACE_FC IBAMR_FC_FUNC_(navier_stokes_side_to_face2d, NAVIER_STOKES_SIDE_TO_FACE2D)
#define NAVIER_STOKES_STAGGERED_ADV_SOURCE_FC                                                                          \
    IBAMR_FC_FUNC_(navier_stokes_staggered_adv_source2d, NAVIER_STOKES_STAGGERED_ADV_SOURCE2D)
#define NAVIER_STOKES_STAGGERED_CONS_SOURCE_FC                                                                         \
    IBAMR_FC_FUNC_(navier_stokes_staggered_cons_source2d, NAVIER_STOKES_STAGGERED_CONS_SOURCE2D)
#define NAVIER_STOKES_STAGGERED_SKEW_SYM_SOURCE_FC                                                                     \
    IBAMR_FC_FUNC_(navier_stokes_staggered_skew_sym_source2d, NAVIER_STOKES_STAGGERED_SKEW_SYM_SOURCE2D)
#endif

#if (NDIM == 3)
#define NAVIER_STOKES_SC_STABLEDT_FC IBAMR_FC_FUNC_(navier_stokes_sc_stabledt3d, NAVIER_STOKES_SC_STABLEDT3D)
#define NAVIER_STOKES_SIDE_TO_FACE_FC IBAMR_FC_FUNC_(navier_stokes_side_to_face3d, NAVIER_STOKES_SIDE_TO_FACE3D)
#define NAVIER_STOKES_STAGGERED_ADV_SOURCE_FC                                                                          \
    IBAMR_FC_FUNC_(navier_stokes_staggered_adv_source3d, NAVIER_STOKES_STAGGERED_ADV_SOURCE3D)
#define NAVIER_STOKES_STAGGERED_CONS_SOURCE_FC                                                                         \
    IBAMR_FC_FUNC_(navier_stokes_staggered_cons_source3d, NAVIER_STOKES_STAGGERED_CONS_SOURCE3D)
#define NAVIER_STOKES_STAGGERED_SKEW_SYM_SOURCE_FC                                                                     \
    IBAMR_FC_FUNC_(navier_stokes_staggered_skew_sym_source3d, NAVIER_STOKES_STAGGERED_SKEW_SYM_SOURCE3D)
#endif

extern "C" {
void NAVIER_STOKES_SC_STABLEDT_FC(const double*,
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

void NAVIER_STOKES_SIDE_TO_FACE_FC(
#if (NDIM == 2)
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
    const int&,
    double*,
    double*,
    const int&
#endif
#if (NDIM == 3)
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
    const double*,
    const int&,
    double*,
    double*,
    double*,
    const int&
#endif
    );

void NAVIER_STOKES_STAGGERED_ADV_SOURCE_FC(
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
    const double*,
    double*,
    double*
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
    const double*,
    double*,
    double*,
    double*
#endif
    );

void NAVIER_STOKES_STAGGERED_CONS_SOURCE_FC(
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
    const double*,
    double*,
    double*
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
    const double*,
    double*,
    double*,
    double*
#endif
    );

void NAVIER_STOKES_STAGGERED_SKEW_SYM_SOURCE_FC(
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
    const double*,
    double*,
    double*
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
    const double*,
    double*,
    double*,
    double*
#endif
    );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Copy data from a side-centered variable to a face-centered variable.
void
copy_side_to_face(const int U_fc_idx, const int U_sc_idx, Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Index<NDIM>& ilower = patch->getBox().lower();
            const Index<NDIM>& iupper = patch->getBox().upper();
            Pointer<SideData<NDIM, double> > U_sc_data = patch->getPatchData(U_sc_idx);
            Pointer<FaceData<NDIM, double> > U_fc_data = patch->getPatchData(U_fc_idx);
#if !defined(NDEBUG)
            TBOX_ASSERT(U_sc_data->getGhostCellWidth().min() == U_sc_data->getGhostCellWidth().max());
            TBOX_ASSERT(U_fc_data->getGhostCellWidth().min() == U_fc_data->getGhostCellWidth().max());
#endif
            const int U_sc_gcw = U_sc_data->getGhostCellWidth().max();
            const int U_fc_gcw = U_fc_data->getGhostCellWidth().max();
            NAVIER_STOKES_SIDE_TO_FACE_FC(ilower(0),
                                          iupper(0),
                                          ilower(1),
                                          iupper(1),
#if (NDIM == 3)
                                          ilower(2),
                                          iupper(2),
#endif
                                          U_sc_data->getPointer(0),
                                          U_sc_data->getPointer(1),
#if (NDIM == 3)
                                          U_sc_data->getPointer(2),
#endif
                                          U_sc_gcw,
                                          U_fc_data->getPointer(0),
                                          U_fc_data->getPointer(1),
#if (NDIM == 3)
                                          U_fc_data->getPointer(2),
#endif
                                          U_fc_gcw);
        }
    }
    return;
} // copy_side_to_face
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredHierarchyIntegrator::INSStaggeredHierarchyIntegrator(const std::string& object_name,
                                                                 Pointer<Database> input_db,
                                                                 bool register_for_restart)
    : INSHierarchyIntegrator(object_name,
                             input_db,
                             new SideVariable<NDIM, double>(object_name + "::U"),
                             new CellVariable<NDIM, double>(object_name + "::P"),
                             new SideVariable<NDIM, double>(object_name + "::F"),
                             new CellVariable<NDIM, double>(object_name + "::Q"),
                             register_for_restart)
{
    // Check to see whether the solver types have been set.
    d_stokes_solver_type = StaggeredStokesSolverManager::UNDEFINED;
    d_stokes_precond_type = StaggeredStokesSolverManager::UNDEFINED;
    if (input_db->keyExists("stokes_solver_type")) d_stokes_solver_type = input_db->getString("stokes_solver_type");
    if (input_db->keyExists("stokes_precond_type")) d_stokes_precond_type = input_db->getString("stokes_precond_type");

    d_velocity_solver_type = SCPoissonSolverManager::UNDEFINED;
    d_velocity_precond_type = SCPoissonSolverManager::UNDEFINED;
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

    // Check to make sure the time stepping types are supported.
    switch (d_viscous_time_stepping_type)
    {
    case BACKWARD_EULER:
    case FORWARD_EULER:
    case TRAPEZOIDAL_RULE:
        break;
    default:
        TBOX_ERROR(d_object_name << "::INSStaggeredHierarchyIntegrator():\n"
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

    // Check to see whether the convective operator type has been set.
    d_convective_op_type = INSStaggeredConvectiveOperatorManager::DEFAULT;
    if (input_db->keyExists("convective_op_type"))
        d_convective_op_type = input_db->getString("convective_op_type");
    else if (input_db->keyExists("convective_operator_type"))
        d_convective_op_type = input_db->getString("convective_operator_type");
    else if (input_db->keyExists("default_convective_op_type"))
        d_convective_op_type = input_db->getString("default_convective_op_type");
    else if (input_db->keyExists("default_convective_operator_type"))
        d_convective_op_type = input_db->getString("default_convective_operator_type");

    // Setup Stokes solver options.
    if (input_db->keyExists("stokes_solver_type"))
    {
        d_stokes_solver_type = input_db->getString("stokes_solver_type");
        if (input_db->keyExists("stokes_solver_db")) d_stokes_solver_db = input_db->getDatabase("stokes_solver_db");
    }
    if (!d_stokes_solver_db) d_stokes_solver_db = new MemoryDatabase("stokes_solver_db");

    if (input_db->keyExists("stokes_precond_type"))
    {
        d_stokes_precond_type = input_db->getString("stokes_precond_type");
        if (input_db->keyExists("stokes_precond_db")) d_stokes_precond_db = input_db->getDatabase("stokes_precond_db");
    }
    if (!d_stokes_precond_db) d_stokes_precond_db = new MemoryDatabase("stokes_precond_db");

    // Setup physical boundary conditions objects.
    d_bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
    d_U_bc_coefs.resize(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_U_bc_coefs[d] = new INSStaggeredVelocityBcCoef(d, this, d_bc_coefs, d_traction_bc_type);
    }
    d_P_bc_coef = new INSStaggeredPressureBcCoef(this, d_bc_coefs, d_traction_bc_type);

    // Initialize all variables.  The velocity, pressure, body force, and fluid
    // source variables were created above in the constructor for the
    // INSHierarchyIntegrator base class.
    d_U_var = INSHierarchyIntegrator::d_U_var;
    d_P_var = INSHierarchyIntegrator::d_P_var;
    d_F_var = INSHierarchyIntegrator::d_F_var;
    d_Q_var = INSHierarchyIntegrator::d_Q_var;
    d_N_old_var = new SideVariable<NDIM, double>(d_object_name + "::N_old");

    d_U_cc_var = new CellVariable<NDIM, double>(d_object_name + "::U_cc", NDIM);
    d_F_cc_var = new CellVariable<NDIM, double>(d_object_name + "::F_cc", NDIM);
#if (NDIM == 2)
    d_Omega_var = new CellVariable<NDIM, double>(d_object_name + "::Omega");
#endif
#if (NDIM == 3)
    d_Omega_var = new CellVariable<NDIM, double>(d_object_name + "::Omega", NDIM);
#endif
    d_Div_U_var = new CellVariable<NDIM, double>(d_object_name + "::Div_U");

#if (NDIM == 3)
    d_Omega_Norm_var = new CellVariable<NDIM, double>(d_object_name + "::|Omega|_2");
#endif
    d_U_regrid_var = new SideVariable<NDIM, double>(d_object_name + "::U_regrid");
    d_U_src_var = new SideVariable<NDIM, double>(d_object_name + "::U_src");
    d_indicator_var = new SideVariable<NDIM, double>(d_object_name + "::indicator");
    d_F_div_var = new SideVariable<NDIM, double>(d_object_name + "::F_div");
    return;
} // INSStaggeredHierarchyIntegrator

INSStaggeredHierarchyIntegrator::~INSStaggeredHierarchyIntegrator()
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        delete d_U_bc_coefs[d];
        d_U_bc_coefs[d] = NULL;
    }
    delete d_P_bc_coef;
    d_P_bc_coef = NULL;
    delete d_fill_after_regrid_phys_bdry_bc_op;
    d_fill_after_regrid_phys_bdry_bc_op = NULL;
    d_velocity_solver.setNull();
    d_pressure_solver.setNull();
    if (d_U_rhs_vec) d_U_rhs_vec->freeVectorComponents();
    if (d_U_adv_vec) d_U_adv_vec->freeVectorComponents();
    if (d_N_vec) d_N_vec->freeVectorComponents();
    if (d_P_rhs_vec) d_P_rhs_vec->freeVectorComponents();
    for (unsigned int k = 0; k < d_nul_vecs.size(); ++k)
    {
        if (d_nul_vecs[k]) d_nul_vecs[k]->freeVectorComponents();
    }
    for (unsigned int k = 0; k < d_U_nul_vecs.size(); ++k)
    {
        if (d_U_nul_vecs[k]) d_U_nul_vecs[k]->freeVectorComponents();
    }
    return;
} // ~INSStaggeredHierarchyIntegrator

Pointer<ConvectiveOperator>
INSStaggeredHierarchyIntegrator::getConvectiveOperator()
{
    if (d_creeping_flow)
    {
        d_convective_op.setNull();
    }
    else if (!d_convective_op)
    {
        INSStaggeredConvectiveOperatorManager* convective_op_manager =
            INSStaggeredConvectiveOperatorManager::getManager();
        d_convective_op = convective_op_manager->allocateOperator(d_convective_op_type,
                                                                  d_object_name + "::ConvectiveOperator",
                                                                  d_convective_op_input_db,
                                                                  d_convective_difference_form,
                                                                  d_U_bc_coefs);
        d_convective_op_needs_init = true;
    }
    return d_convective_op;
} // getConvectiveOperator

Pointer<PoissonSolver>
INSStaggeredHierarchyIntegrator::getVelocitySubdomainSolver()
{
    if (!d_velocity_solver)
    {
        d_velocity_solver = SCPoissonSolverManager::getManager()->allocateSolver(d_velocity_solver_type,
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
INSStaggeredHierarchyIntegrator::getPressureSubdomainSolver()
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
INSStaggeredHierarchyIntegrator::setStokesSolver(Pointer<StaggeredStokesSolver> stokes_solver)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_stokes_solver);
#endif
    d_stokes_solver = stokes_solver;
    d_stokes_solver_needs_init = true;
    return;
} // setStokesSolver

Pointer<StaggeredStokesSolver>
INSStaggeredHierarchyIntegrator::getStokesSolver()
{
    if (!d_stokes_solver)
    {
        d_stokes_solver = StaggeredStokesSolverManager::getManager()->allocateSolver(d_stokes_solver_type,
                                                                                     d_object_name + "::stokes_solver",
                                                                                     d_stokes_solver_db,
                                                                                     "stokes_",
                                                                                     d_stokes_precond_type,
                                                                                     d_object_name + "::stokes_precond",
                                                                                     d_stokes_precond_db,
                                                                                     "stokes_pc_");
        d_stokes_solver_needs_init = true;
    }
    return d_stokes_solver;
} // getStokesSolver

void
INSStaggeredHierarchyIntegrator::setStokesSolverNeedsInit()
{
    d_stokes_solver_needs_init = true;
    return;
}

void
INSStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                               Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;

    // Setup solvers.
    if (d_stokes_solver_type == StaggeredStokesSolverManager::UNDEFINED)
    {
        d_stokes_solver_type = StaggeredStokesSolverManager::PETSC_KRYLOV_SOLVER;
        d_stokes_solver_db->putString("ksp_type", "fgmres");
    }

    if (d_stokes_precond_type == StaggeredStokesSolverManager::UNDEFINED)
    {
        d_stokes_precond_type = StaggeredStokesSolverManager::DEFAULT_BLOCK_PRECONDITIONER;
        d_stokes_precond_db->putInteger("max_iterations", 1);
    }

    if (d_velocity_solver_type == SCPoissonSolverManager::UNDEFINED)
    {
        d_velocity_solver_type = SCPoissonSolverManager::PETSC_KRYLOV_SOLVER;
        d_velocity_solver_db->putString("ksp_type", "richardson");
        d_velocity_solver_db->putInteger("max_iterations", 10);
        d_velocity_solver_db->putDouble("rel_residual_tol", 1.0e-1);
    }

    if (d_velocity_precond_type == SCPoissonSolverManager::UNDEFINED)
    {
        const int max_levels = gridding_alg->getMaxLevels();
        if (max_levels == 1)
        {
            d_velocity_precond_type = SCPoissonSolverManager::DEFAULT_LEVEL_SOLVER;
        }
        else
        {
            d_velocity_precond_type = SCPoissonSolverManager::DEFAULT_FAC_PRECONDITIONER;
        }
        d_velocity_precond_db->putInteger("max_iterations", 1);
    }

    if (d_pressure_solver_type == CCPoissonSolverManager::UNDEFINED)
    {
        d_pressure_solver_type = CCPoissonSolverManager::PETSC_KRYLOV_SOLVER;
        d_pressure_solver_db->putString("ksp_type", "richardson");
        d_pressure_solver_db->putInteger("max_iterations", 10);
        d_pressure_solver_db->putDouble("rel_residual_tol", 1.0e-1);
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
    d_hier_sc_data_ops =
        hier_ops_manager->getOperationsDouble(new SideVariable<NDIM, double>("sc_var"), hierarchy, true);
    d_hier_math_ops = buildHierarchyMathOps(d_hierarchy);

    // Register state variables that are maintained by the
    // INSStaggeredHierarchyIntegrator.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    grid_geom->addSpatialRefineOperator(new CartSideDoubleSpecializedConstantRefine());
    grid_geom->addSpatialRefineOperator(new CartSideDoubleSpecializedLinearRefine());

    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> side_ghosts = SIDEG;
    const IntVector<NDIM> no_ghosts = 0;

    registerVariable(d_U_current_idx,
                     d_U_new_idx,
                     d_U_scratch_idx,
                     d_U_var,
                     side_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE",
                     d_U_init);

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
                         side_ghosts,
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
                     side_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    d_rho_var = INSHierarchyIntegrator::d_rho_var;
    if (INSHierarchyIntegrator::d_rho_var && !d_rho_var)
    {
        TBOX_ERROR("INSStaggeredHierarchyIntegrator::initializeHierarchyIntegrator():\n"
                   << "  mass density variable must be side-centered.\n");
    }
    if (d_rho_var)
    {
        registerVariable(d_rho_current_idx,
                         d_rho_new_idx,
                         d_rho_scratch_idx,
                         d_rho_var,
                         side_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE",
                         d_rho_fcn);
    }
    else
    {
        d_rho_current_idx = -1;
        d_rho_new_idx = -1;
        d_rho_scratch_idx = -1;
    }

    // Register plot variables that are maintained by the
    // INSCollocatedHierarchyIntegrator.
    registerVariable(d_U_cc_idx, d_U_cc_var, no_ghosts, getCurrentContext());
    if (d_F_fcn)
    {
        registerVariable(d_F_cc_idx, d_F_cc_var, no_ghosts, getCurrentContext());
    }
    else
    {
        d_F_cc_idx = -1;
    }
    registerVariable(d_Omega_idx, d_Omega_var, no_ghosts, getCurrentContext());
    registerVariable(d_Div_U_idx, d_Div_U_var, cell_ghosts, getCurrentContext());
    if (d_rho_var)
    {
        d_rho_cc_var = new CellVariable<NDIM, double>(d_object_name + "::rho_cc");
        registerVariable(d_rho_cc_idx, d_rho_cc_var, no_ghosts, getCurrentContext());
    }

// Register scratch variables that are maintained by the
// INSStaggeredHierarchyIntegrator.
#if (NDIM == 3)
    registerVariable(d_Omega_Norm_idx, d_Omega_Norm_var, no_ghosts);
#endif
    registerVariable(d_U_regrid_idx, d_U_regrid_var, CartSideDoubleDivPreservingRefine::REFINE_OP_STENCIL_WIDTH);
    registerVariable(d_U_src_idx, d_U_src_var, CartSideDoubleDivPreservingRefine::REFINE_OP_STENCIL_WIDTH);
    registerVariable(d_indicator_idx, d_indicator_var, CartSideDoubleDivPreservingRefine::REFINE_OP_STENCIL_WIDTH);
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
            d_visit_writer->registerPlotQuantity("U", "VECTOR", d_U_cc_idx, 0, d_U_scale);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == 0) d_visit_writer->registerPlotQuantity("U_x", "SCALAR", d_U_cc_idx, d, d_U_scale);
                if (d == 1) d_visit_writer->registerPlotQuantity("U_y", "SCALAR", d_U_cc_idx, d, d_U_scale);
                if (d == 2) d_visit_writer->registerPlotQuantity("U_z", "SCALAR", d_U_cc_idx, d, d_U_scale);
            }
        }

        if (d_output_P)
        {
            d_visit_writer->registerPlotQuantity("P", "SCALAR", d_P_current_idx, 0, d_P_scale);
        }

        if (d_F_fcn && d_output_F)
        {
            d_visit_writer->registerPlotQuantity("F", "VECTOR", d_F_cc_idx, 0, d_F_scale);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == 0) d_visit_writer->registerPlotQuantity("F_x", "SCALAR", d_F_cc_idx, d, d_F_scale);
                if (d == 1) d_visit_writer->registerPlotQuantity("F_y", "SCALAR", d_F_cc_idx, d, d_F_scale);
                if (d == 2) d_visit_writer->registerPlotQuantity("F_z", "SCALAR", d_F_cc_idx, d, d_F_scale);
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
        }
    }

    // Setup a specialized coarsen algorithm.
    Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
    Pointer<CoarsenOperator<NDIM> > coarsen_op = grid_geom->lookupCoarsenOperator(d_U_var, "CONSERVATIVE_COARSEN");
    coarsen_alg->registerCoarsen(d_U_scratch_idx, d_U_scratch_idx, coarsen_op);
    registerCoarsenAlgorithm(d_object_name + "::CONVECTIVE_OP", coarsen_alg);

    // Setup the Stokes solver.
    d_stokes_solver = getStokesSolver();
    Pointer<LinearSolver> p_stokes_linear_solver = d_stokes_solver;
    if (!p_stokes_linear_solver)
    {
        Pointer<NewtonKrylovSolver> p_stokes_newton_solver = d_stokes_solver;
        if (p_stokes_newton_solver) p_stokes_linear_solver = p_stokes_newton_solver->getLinearSolver();
    }
    if (p_stokes_linear_solver)
    {
        Pointer<StaggeredStokesBlockPreconditioner> p_stokes_block_pc = p_stokes_linear_solver;
        if (!p_stokes_block_pc)
        {
            Pointer<KrylovLinearSolver> p_stokes_krylov_solver = p_stokes_linear_solver;
            if (p_stokes_krylov_solver) p_stokes_block_pc = p_stokes_krylov_solver->getPreconditioner();
        }
        if (p_stokes_block_pc)
        {
            if (p_stokes_block_pc->needsVelocitySubdomainSolver())
            {
                p_stokes_block_pc->setVelocitySubdomainSolver(getVelocitySubdomainSolver());
            }
            if (p_stokes_block_pc->needsPressureSubdomainSolver())
            {
                p_stokes_block_pc->setPressureSubdomainSolver(getPressureSubdomainSolver());
            }
        }
    }

    // Setup the convective operator.
    d_convective_op = getConvectiveOperator();

    // Setup a boundary op to set velocity boundary conditions on regrid.
    d_fill_after_regrid_phys_bdry_bc_op = new CartSideRobinPhysBdryOp(d_U_scratch_idx, d_U_bc_coefs, false);

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
} // initializeHierarchyIntegrator

void
INSStaggeredHierarchyIntegrator::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                          Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    HierarchyIntegrator::initializePatchHierarchy(hierarchy, gridding_alg);

    // Project the velocity field if this is the initial time step.
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
            copy_side_to_face(U_adv_diff_current_idx, d_U_current_idx, d_hierarchy);
        }
    }
    return;
} // initializePatchHierarhcy

void
INSStaggeredHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                              const double new_time,
                                                              const int num_cycles)
{
    INSHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Keep track of the number of cycles to be used for the present integration
    // step.
    if (!d_creeping_flow && (d_current_num_cycles == 1) &&
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
    d_P_rhs_vec->allocateVectorData(current_time);
    d_P_rhs_vec->setToScalar(0.0);
    if (!d_creeping_flow)
    {
        d_U_adv_vec->allocateVectorData(current_time);
        d_U_adv_vec->setToScalar(0.0);
        d_N_vec->allocateVectorData(current_time);
        d_N_vec->setToScalar(0.0);
    }

    // Cache BC data.
    d_bc_helper->cacheBcCoefData(d_bc_coefs, new_time, d_hierarchy);

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
    const Pointer<SideVariable<NDIM, double> > U_rhs_var = d_U_rhs_vec->getComponentVariable(0);
    d_hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(d_U_bc_coefs,
                                                              /*P_bc_coef*/ NULL,
                                                              d_U_scratch_idx,
                                                              /*P_data_idx*/ -1,
                                                              /*homogeneous_bc*/ false);
    d_U_bdry_bc_fill_op->fillData(current_time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs,
                                                              /*P_bc_coef*/ NULL);
    d_hier_math_ops->laplace(
        U_rhs_idx, U_rhs_var, U_rhs_problem_coefs, d_U_scratch_idx, d_U_var, d_no_fill_op, current_time);
    d_hier_sc_data_ops->copyData(d_U_src_idx, d_U_scratch_idx, /*interior_only*/ false);

    // Set the initial guess.
    d_hier_sc_data_ops->copyData(d_U_new_idx, d_U_current_idx);
    d_hier_cc_data_ops->copyData(d_P_new_idx, d_P_current_idx);

    // Set up inhomogeneous BCs.
    d_stokes_solver->setHomogeneousBc(false);

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
            copy_side_to_face(U_adv_diff_current_idx, d_U_current_idx, d_hierarchy);
        }
        d_adv_diff_hier_integrator->preprocessIntegrateHierarchy(current_time, new_time, adv_diff_num_cycles);
        const int U_adv_diff_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getScratchContext());
        if (isAllocatedPatchData(U_adv_diff_scratch_idx))
        {
            d_hier_fc_data_ops->copyData(U_adv_diff_scratch_idx, U_adv_diff_current_idx);
        }
        const int U_adv_diff_new_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getNewContext());
        if (isAllocatedPatchData(U_adv_diff_new_idx))
        {
            d_hier_fc_data_ops->copyData(U_adv_diff_new_idx, U_adv_diff_current_idx);
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
        d_hier_sc_data_ops->copyData(U_adv_idx, d_U_current_idx);
        for (int ln = finest_ln; ln > coarsest_ln; --ln)
        {
            Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
            Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
            Pointer<CoarsenOperator<NDIM> > coarsen_op =
                grid_geom->lookupCoarsenOperator(d_U_var, "CONSERVATIVE_COARSEN");
            coarsen_alg->registerCoarsen(U_adv_idx, U_adv_idx, coarsen_op);
            coarsen_alg->resetSchedule(getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]);
            getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]->coarsenData();
            getCoarsenAlgorithm(d_object_name + "::CONVECTIVE_OP")
                ->resetSchedule(getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]);
        }
        d_convective_op->setAdvectionVelocity(d_U_adv_vec->getComponentDescriptorIndex(0));
        d_convective_op->setSolutionTime(current_time);
        d_convective_op->apply(*d_U_adv_vec, *d_N_vec);
        const int N_idx = d_N_vec->getComponentDescriptorIndex(0);
        d_hier_sc_data_ops->copyData(d_N_old_new_idx, N_idx);
        if (convective_time_stepping_type == FORWARD_EULER)
        {
            d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0),
                                     -1.0 * rho,
                                     N_idx,
                                     d_rhs_vec->getComponentDescriptorIndex(0));
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0),
                                     -0.5 * rho,
                                     N_idx,
                                     d_rhs_vec->getComponentDescriptorIndex(0));
        }
    }

    // Compute variable mass density.
    if (d_rho_var)
    {
        if (d_rho_fcn)
        {
            d_rho_fcn->setDataOnPatchHierarchy(d_rho_current_idx, d_rho_var, d_hierarchy, current_time);
            d_rho_fcn->setDataOnPatchHierarchy(d_rho_new_idx, d_rho_var, d_hierarchy, new_time);
        }
        else
        {
            d_hier_sc_data_ops->copyData(d_rho_new_idx, d_rho_current_idx);
        }
        d_hier_sc_data_ops->linearSum(d_rho_scratch_idx, 0.5, d_rho_current_idx, 0.5, d_rho_new_idx);
    }

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
INSStaggeredHierarchyIntegrator::integrateHierarchy(const double current_time,
                                                    const double new_time,
                                                    const int cycle_num)
{
    INSHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);

    // Check to make sure that the number of cycles is what we expect it to be.
    const int expected_num_cycles = getNumberOfCycles();
    if (d_current_num_cycles != expected_num_cycles)
    {
        IBAMR_DO_ONCE(
            {
                pout << "INSStaggeredHierarchyIntegrator::integrateHierarchy():\n"
                     << "  WARNING: num_cycles = " << d_current_num_cycles
                     << " but expected num_cycles = " << expected_num_cycles << ".\n";
            });
    }

    // Update the state variables of any linked advection-diffusion solver.
    if (d_adv_diff_hier_integrator)
    {
        d_adv_diff_hier_integrator->integrateHierarchy(current_time, new_time, cycle_num);
    }

    // Setup the solution and right-hand-side vectors.
    setupSolverVectors(d_sol_vec, d_rhs_vec, current_time, new_time, cycle_num);

    // Solve for u(n+1), p(n+1/2).
    d_stokes_solver->solveSystem(*d_sol_vec, *d_rhs_vec);
    if (d_enable_logging)
        plog << d_object_name
             << "::integrateHierarchy(): stokes solve number of iterations = " << d_stokes_solver->getNumIterations()
             << "\n";
    if (d_enable_logging)
        plog << d_object_name
             << "::integrateHierarchy(): stokes solve residual norm        = " << d_stokes_solver->getResidualNorm()
             << "\n";

    // Reset the solution and right-hand-side vectors.
    resetSolverVectors(d_sol_vec, d_rhs_vec, current_time, new_time, cycle_num);

    // Update the state variables of any linked advection-diffusion solver.
    if (d_adv_diff_hier_integrator)
    {
        // Update the advection velocities used by the advection-diffusion
        // solver.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int U_adv_diff_new_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getNewContext());
        if (isAllocatedPatchData(U_adv_diff_new_idx))
        {
            copy_side_to_face(U_adv_diff_new_idx, d_U_new_idx, d_hierarchy);
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
INSStaggeredHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
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
        PatchSideDataOpsReal<NDIM, double> patch_sc_ops;
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
                Pointer<SideData<NDIM, double> > u_sc_new_data = patch->getPatchData(d_U_new_idx);
                double u_max = 0.0;
                u_max = patch_sc_ops.maxNorm(u_sc_new_data, patch_box);
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
        d_hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_new_idx);
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
    d_P_rhs_vec->deallocateVectorData();
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
INSStaggeredHierarchyIntegrator::regridHierarchy()
{
    const int coarsest_ln = 0;

    // Determine the divergence of the velocity field before regridding.
    d_hier_math_ops->div(d_Div_U_idx,
                         d_Div_U_var,
                         1.0,
                         d_U_current_idx,
                         d_U_var,
                         d_no_fill_op,
                         d_integrator_time,
                         /*synch_cf_bdry*/ false,
                         -1.0,
                         d_Q_current_idx,
                         d_Q_var);
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const double Div_U_norm_1_pre = d_hier_cc_data_ops->L1Norm(d_Div_U_idx, wgt_cc_idx);
    const double Div_U_norm_2_pre = d_hier_cc_data_ops->L2Norm(d_Div_U_idx, wgt_cc_idx);
    const double Div_U_norm_oo_pre = d_hier_cc_data_ops->maxNorm(d_Div_U_idx, wgt_cc_idx);

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

    // Determine the divergence of the velocity field after regridding.
    d_hier_math_ops->div(d_Div_U_idx,
                         d_Div_U_var,
                         1.0,
                         d_U_current_idx,
                         d_U_var,
                         d_no_fill_op,
                         d_integrator_time,
                         /*synch_cf_bdry*/ true,
                         -1.0,
                         d_Q_current_idx,
                         d_Q_var);
    const double Div_U_norm_1_post = d_hier_cc_data_ops->L1Norm(d_Div_U_idx, wgt_cc_idx);
    const double Div_U_norm_2_post = d_hier_cc_data_ops->L2Norm(d_Div_U_idx, wgt_cc_idx);
    const double Div_U_norm_oo_post = d_hier_cc_data_ops->maxNorm(d_Div_U_idx, wgt_cc_idx);

    // Project the interpolated velocity if needed.
    if (Div_U_norm_1_post > d_regrid_max_div_growth_factor * Div_U_norm_1_pre ||
        Div_U_norm_2_post > d_regrid_max_div_growth_factor * Div_U_norm_2_pre ||
        Div_U_norm_oo_post > d_regrid_max_div_growth_factor * Div_U_norm_oo_pre)
    {
        pout << d_object_name << "::regridHierarchy():\n"
             << "  WARNING: projecting the interpolated velocity field\n";
        regridProjection();
    }

    // Synchronize the state data on the patch hierarchy.
    synchronizeHierarchyData(CURRENT_DATA);
    return;
} // regridHierarchy

void
INSStaggeredHierarchyIntegrator::setupSolverVectors(const Pointer<SAMRAIVectorReal<NDIM, double> >& sol_vec,
                                                    const Pointer<SAMRAIVectorReal<NDIM, double> >& rhs_vec,
                                                    const double current_time,
                                                    const double new_time,
                                                    const int cycle_num)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    const double rho = d_problem_coefs.getRho();

    if (rhs_vec->getComponentDescriptorIndex(0) != d_U_rhs_vec->getComponentDescriptorIndex(0))
    {
        d_hier_sc_data_ops->copyData(rhs_vec->getComponentDescriptorIndex(0),
                                     d_U_rhs_vec->getComponentDescriptorIndex(0));
    }
    if (rhs_vec->getComponentDescriptorIndex(1) != d_P_rhs_vec->getComponentDescriptorIndex(0))
    {
        d_hier_cc_data_ops->copyData(rhs_vec->getComponentDescriptorIndex(1),
                                     d_P_rhs_vec->getComponentDescriptorIndex(0));
    }

    // Account for the convective acceleration term.
    if (!d_creeping_flow)
    {
        const TimeSteppingType convective_time_stepping_type = getConvectiveTimeSteppingType(cycle_num);
        if (cycle_num > 0)
        {
            const int U_adv_idx = d_U_adv_vec->getComponentDescriptorIndex(0);
            double apply_time = std::numeric_limits<double>::quiet_NaN();
            if (convective_time_stepping_type == MIDPOINT_RULE)
            {
                d_hier_sc_data_ops->linearSum(U_adv_idx, 0.5, d_U_current_idx, 0.5, d_U_new_idx);
                apply_time = half_time;
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                d_hier_sc_data_ops->copyData(U_adv_idx, d_U_new_idx);
                apply_time = new_time;
            }
            for (int ln = finest_ln; ln > coarsest_ln; --ln)
            {
                Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
                Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
                Pointer<CoarsenOperator<NDIM> > coarsen_op =
                    grid_geom->lookupCoarsenOperator(d_U_var, "CONSERVATIVE_COARSEN");
                coarsen_alg->registerCoarsen(U_adv_idx, U_adv_idx, coarsen_op);
                coarsen_alg->resetSchedule(getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]);
                getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]->coarsenData();
                getCoarsenAlgorithm(d_object_name + "::CONVECTIVE_OP")
                    ->resetSchedule(getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]);
            }
            d_convective_op->setAdvectionVelocity(d_U_adv_vec->getComponentDescriptorIndex(0));
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
            d_hier_sc_data_ops->linearSum(N_idx, 1.0 + 0.5 * omega, N_idx, -0.5 * omega, d_N_old_current_idx);
        }
        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
        {
            d_hier_sc_data_ops->axpy(
                rhs_vec->getComponentDescriptorIndex(0), -1.0 * rho, N_idx, rhs_vec->getComponentDescriptorIndex(0));
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_sc_data_ops->axpy(
                rhs_vec->getComponentDescriptorIndex(0), -0.5 * rho, N_idx, rhs_vec->getComponentDescriptorIndex(0));
        }
    }

    // Account for body forcing terms.
    if (d_F_fcn)
    {
        d_F_fcn->setDataOnPatchHierarchy(d_F_scratch_idx, d_F_var, d_hierarchy, half_time);
        d_hier_sc_data_ops->add(
            rhs_vec->getComponentDescriptorIndex(0), rhs_vec->getComponentDescriptorIndex(0), d_F_scratch_idx);
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
            d_hier_sc_data_ops->linearSum(d_U_scratch_idx, 0.5, d_U_current_idx, 0.5, d_U_new_idx);
            computeDivSourceTerm(d_F_div_idx, d_Q_scratch_idx, d_U_scratch_idx);
        }
        d_hier_sc_data_ops->axpy(
            rhs_vec->getComponentDescriptorIndex(0), rho, d_F_div_idx, rhs_vec->getComponentDescriptorIndex(0));
        d_hier_cc_data_ops->subtract(
            rhs_vec->getComponentDescriptorIndex(1), rhs_vec->getComponentDescriptorIndex(1), d_Q_new_idx);
    }

    // Set solution components to equal most recent approximations to u(n+1) and
    // p(n+1/2).
    d_hier_sc_data_ops->copyData(sol_vec->getComponentDescriptorIndex(0), d_U_new_idx);
    d_hier_cc_data_ops->copyData(sol_vec->getComponentDescriptorIndex(1), d_P_new_idx);

    // Synchronize solution and right-hand-side data before solve.
    typedef SideDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent;
    SynchronizationTransactionComponent sol_synch_transaction =
        SynchronizationTransactionComponent(sol_vec->getComponentDescriptorIndex(0), "CONSERVATIVE_COARSEN");
    d_side_synch_op->resetTransactionComponent(sol_synch_transaction);
    d_side_synch_op->synchronizeData(current_time);
    SynchronizationTransactionComponent rhs_synch_transaction =
        SynchronizationTransactionComponent(rhs_vec->getComponentDescriptorIndex(0), "CONSERVATIVE_COARSEN");
    d_side_synch_op->resetTransactionComponent(rhs_synch_transaction);
    d_side_synch_op->synchronizeData(current_time);
    SynchronizationTransactionComponent default_synch_transaction =
        SynchronizationTransactionComponent(d_U_scratch_idx, "CONSERVATIVE_COARSEN");
    d_side_synch_op->resetTransactionComponent(default_synch_transaction);
    return;
} // setupSolverVectors

void
INSStaggeredHierarchyIntegrator::resetSolverVectors(const Pointer<SAMRAIVectorReal<NDIM, double> >& sol_vec,
                                                    const Pointer<SAMRAIVectorReal<NDIM, double> >& rhs_vec,
                                                    const double current_time,
                                                    const double /*new_time*/,
                                                    const int cycle_num)
{
    // Synchronize solution data after solve.
    typedef SideDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent;
    SynchronizationTransactionComponent sol_synch_transaction =
        SynchronizationTransactionComponent(sol_vec->getComponentDescriptorIndex(0), "CONSERVATIVE_COARSEN");
    d_side_synch_op->synchronizeData(current_time);
    SynchronizationTransactionComponent default_synch_transaction =
        SynchronizationTransactionComponent(d_U_scratch_idx, "CONSERVATIVE_COARSEN");
    d_side_synch_op->resetTransactionComponent(default_synch_transaction);

    // Pull out solution components.
    d_hier_sc_data_ops->copyData(d_U_new_idx, sol_vec->getComponentDescriptorIndex(0));
    d_hier_cc_data_ops->copyData(d_P_new_idx, sol_vec->getComponentDescriptorIndex(1));

    // Reset the right-hand side vector.
    const double rho = d_problem_coefs.getRho();
    if (!d_creeping_flow)
    {
        const TimeSteppingType convective_time_stepping_type = getConvectiveTimeSteppingType(cycle_num);
        const int N_idx = d_N_vec->getComponentDescriptorIndex(0);
        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
        {
            d_hier_sc_data_ops->axpy(
                rhs_vec->getComponentDescriptorIndex(0), +1.0 * rho, N_idx, rhs_vec->getComponentDescriptorIndex(0));
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_sc_data_ops->axpy(
                rhs_vec->getComponentDescriptorIndex(0), +0.5 * rho, N_idx, rhs_vec->getComponentDescriptorIndex(0));
        }
    }
    if (d_F_fcn)
    {
        d_hier_sc_data_ops->subtract(
            rhs_vec->getComponentDescriptorIndex(0), rhs_vec->getComponentDescriptorIndex(0), d_F_scratch_idx);
        d_hier_sc_data_ops->copyData(d_F_new_idx, d_F_scratch_idx);
    }
    if (d_Q_fcn)
    {
        d_hier_sc_data_ops->axpy(
            rhs_vec->getComponentDescriptorIndex(0), -rho, d_F_div_idx, rhs_vec->getComponentDescriptorIndex(0));
        d_hier_cc_data_ops->add(
            rhs_vec->getComponentDescriptorIndex(1), rhs_vec->getComponentDescriptorIndex(1), d_Q_new_idx);
    }
    return;
} // resetSolverVectors

/////////////////////////////// PROTECTED ////////////////////////////////////

void
INSStaggeredHierarchyIntegrator::initializeLevelDataSpecialized(const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
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

    // Correct the divergence of the interpolated velocity data.
    if (!initial_time && level_number > 0)
    {
        // Allocate scratch data.
        ComponentSelector scratch_data;
        scratch_data.setFlag(d_U_regrid_idx);
        scratch_data.setFlag(d_U_src_idx);
        scratch_data.setFlag(d_indicator_idx);
        level->allocatePatchData(scratch_data, init_data_time);
        if (old_level) old_level->allocatePatchData(scratch_data, init_data_time);

        // Set the indicator data to equal "0" in each patch of the new patch
        // level, and initialize values of U to cause floating point errors if
        // we fail to re-initialize it properly.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > indicator_data = patch->getPatchData(d_indicator_idx);
            indicator_data->fillAll(0.0);

            Pointer<SideData<NDIM, double> > U_current_data = patch->getPatchData(d_U_current_idx);
            Pointer<SideData<NDIM, double> > U_regrid_data = patch->getPatchData(d_U_regrid_idx);
            Pointer<SideData<NDIM, double> > U_src_data = patch->getPatchData(d_U_src_idx);
            U_current_data->fillAll(std::numeric_limits<double>::quiet_NaN());
            U_regrid_data->fillAll(std::numeric_limits<double>::quiet_NaN());
            U_src_data->fillAll(std::numeric_limits<double>::quiet_NaN());
        }

        if (old_level)
        {
            // Set the indicator data to equal "1" on each patch of the old
            // patch level and reset U.
            for (PatchLevel<NDIM>::Iterator p(old_level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = old_level->getPatch(p());

                Pointer<SideData<NDIM, double> > indicator_data = patch->getPatchData(d_indicator_idx);
                indicator_data->fillAll(1.0);

                Pointer<SideData<NDIM, double> > U_current_data = patch->getPatchData(d_U_current_idx);
                Pointer<SideData<NDIM, double> > U_regrid_data = patch->getPatchData(d_U_regrid_idx);
                Pointer<SideData<NDIM, double> > U_src_data = patch->getPatchData(d_U_src_idx);
                U_regrid_data->copy(*U_current_data);
                U_src_data->copy(*U_current_data);
            }

            // Create a communications schedule to copy data from the old patch
            // level to the new patch level.
            //
            // Note that this will set the indicator data to equal "1" at each
            // location in the new patch level that is a copy of a location from
            // the old patch level.
            RefineAlgorithm<NDIM> copy_data;
            copy_data.registerRefine(d_U_regrid_idx, d_U_regrid_idx, d_U_regrid_idx, NULL);
            copy_data.registerRefine(d_U_src_idx, d_U_src_idx, d_U_src_idx, NULL);
            copy_data.registerRefine(d_indicator_idx, d_indicator_idx, d_indicator_idx, NULL);
            ComponentSelector bc_fill_data;
            bc_fill_data.setFlag(d_U_regrid_idx);
            bc_fill_data.setFlag(d_U_src_idx);
            CartSideRobinPhysBdryOp phys_bdry_bc_op(bc_fill_data, d_U_bc_coefs, false);
            copy_data.createSchedule(level, old_level, &phys_bdry_bc_op)->fillData(init_data_time);
        }

        // Setup the divergence- and curl-preserving prolongation refine
        // algorithm and refine the velocity data.
        RefineAlgorithm<NDIM> fill_div_free_prolongation;
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        fill_div_free_prolongation.registerRefine(d_U_current_idx, d_U_current_idx, d_U_regrid_idx, NULL);
        Pointer<RefineOperator<NDIM> > refine_op =
            grid_geom->lookupRefineOperator(d_U_var, "CONSERVATIVE_LINEAR_REFINE");
        Pointer<CoarsenOperator<NDIM> > coarsen_op = grid_geom->lookupCoarsenOperator(d_U_var, "CONSERVATIVE_COARSEN");
        CartSideRobinPhysBdryOp phys_bdry_bc_op(d_U_regrid_idx, d_U_bc_coefs, false);
        CartSideDoubleDivPreservingRefine div_preserving_op(
            d_U_regrid_idx, d_U_src_idx, d_indicator_idx, refine_op, coarsen_op, init_data_time, &phys_bdry_bc_op);
        fill_div_free_prolongation.createSchedule(level, old_level, level_number - 1, hierarchy, &div_preserving_op)
            ->fillData(init_data_time);

        // Free scratch data.
        level->deallocatePatchData(scratch_data);
        if (old_level) old_level->deallocatePatchData(scratch_data);
    }

    // Initialize level data at the initial time.
    if (initial_time)
    {
        // Initialize the maximum value of |Omega|_2 on the grid.
        if (d_using_vorticity_tagging)
        {
            if (level_number == 0) d_Omega_max = 0.0;

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
                hier_ops_manager->getOperationsDouble(d_U_cc_var, d_hierarchy, true);
            Pointer<HierarchySideDataOpsReal<NDIM, double> > hier_sc_data_ops =
                hier_ops_manager->getOperationsDouble(d_U_var, d_hierarchy, true);
            hier_sc_data_ops->resetLevels(0, level_number);
            hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
            typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent
                InterpolationTransactionComponent;
            InterpolationTransactionComponent U_bc_component(d_U_scratch_idx,
                                                             DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             d_bdry_extrap_type,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_U_bc_coefs);
            HierarchyGhostCellInterpolation U_bdry_bc_fill_op;
            U_bdry_bc_fill_op.initializeOperatorState(U_bc_component, d_hierarchy, 0, level_number);
            U_bdry_bc_fill_op.fillData(init_data_time);

            // Compute max |Omega|_2.
            HierarchyMathOps hier_math_ops(d_object_name + "::HierarchyLevelMathOps", d_hierarchy, 0, level_number);
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

            // Deallocate scratch data.
            for (int ln = 0; ln <= level_number; ++ln)
            {
                hierarchy->getPatchLevel(ln)->deallocatePatchData(d_U_scratch_idx);
#if (NDIM == 3)
                hierarchy->getPatchLevel(ln)->deallocatePatchData(d_Omega_Norm_idx);
#endif
            }
        }
    }
    return;
} // initializeLevelDataSpecialized

void
INSStaggeredHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
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

    d_hier_sc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_sc_data_ops->resetLevels(0, finest_hier_level);

    // Setup the patch boundary filling objects.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent U_bc_component(d_U_scratch_idx,
                                                     DATA_REFINE_TYPE,
                                                     USE_CF_INTERPOLATION,
                                                     DATA_COARSEN_TYPE,
                                                     d_bdry_extrap_type,
                                                     CONSISTENT_TYPE_2_BDRY,
                                                     d_U_bc_coefs);
    d_U_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_U_bdry_bc_fill_op->initializeOperatorState(U_bc_component, d_hierarchy);

    InterpolationTransactionComponent P_bc_component(d_P_scratch_idx,
                                                     DATA_REFINE_TYPE,
                                                     USE_CF_INTERPOLATION,
                                                     DATA_COARSEN_TYPE,
                                                     d_bdry_extrap_type,
                                                     CONSISTENT_TYPE_2_BDRY,
                                                     d_P_bc_coef);
    d_P_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_P_bdry_bc_fill_op->initializeOperatorState(P_bc_component, d_hierarchy);

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

    // Setup the patch boundary synchronization objects.
    typedef SideDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent;
    SynchronizationTransactionComponent synch_transaction =
        SynchronizationTransactionComponent(d_U_scratch_idx, "CONSERVATIVE_COARSEN");
    d_side_synch_op = new SideDataSynchronization();
    d_side_synch_op->initializeOperatorState(synch_transaction, d_hierarchy);

    // Indicate that vectors and solvers need to be re-initialized.
    d_coarsest_reset_ln = coarsest_level;
    d_finest_reset_ln = finest_level;
    d_vectors_need_init = true;
    d_convective_op_needs_init = true;
    d_velocity_solver_needs_init = true;
    d_pressure_solver_needs_init = true;
    d_stokes_solver_needs_init = true;
    return;
} // resetHierarchyConfigurationSpecialized

void
INSStaggeredHierarchyIntegrator::applyGradientDetectorSpecialized(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
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
INSStaggeredHierarchyIntegrator::setupPlotDataSpecialized()
{
    Pointer<VariableContext> ctx = getCurrentContext();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    static const bool synch_cf_interface = true;

    // Interpolate u to cell centers.
    if (d_output_U)
    {
        const int U_sc_idx = var_db->mapVariableAndContextToIndex(d_U_var, ctx);
        const int U_cc_idx = var_db->mapVariableAndContextToIndex(d_U_cc_var, ctx);
        d_hier_math_ops->interp(
            U_cc_idx, d_U_cc_var, U_sc_idx, d_U_var, d_no_fill_op, d_integrator_time, synch_cf_interface);
    }

    // Interpolate f to cell centers.
    if (d_F_fcn && d_output_F)
    {
        const int F_sc_idx = var_db->mapVariableAndContextToIndex(d_F_var, ctx);
        const int F_cc_idx = var_db->mapVariableAndContextToIndex(d_F_cc_var, ctx);
        d_hier_math_ops->interp(
            F_cc_idx, d_F_cc_var, F_sc_idx, d_F_var, d_no_fill_op, d_integrator_time, synch_cf_interface);
    }

    // Compute Omega = curl U.
    if (d_output_Omega)
    {
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(d_U_scratch_idx, d_integrator_time);
        }
        d_hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
        d_U_bdry_bc_fill_op->fillData(d_integrator_time);
        d_hier_math_ops->curl(d_Omega_idx, d_Omega_var, d_U_scratch_idx, d_U_var, d_no_fill_op, d_integrator_time);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(d_U_scratch_idx);
        }
    }

    // Compute Div U.
    if (d_output_Div_U)
    {
        d_hier_math_ops->div(
            d_Div_U_idx, d_Div_U_var, 1.0, d_U_current_idx, d_U_var, d_no_fill_op, d_integrator_time, false);
    }
    return;
} // setupPlotDataSpecialized

void
INSStaggeredHierarchyIntegrator::regridProjection()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const double volume = d_hier_math_ops->getVolumeOfPhysicalDomain();

    // Setup the solver vectors.
    SAMRAIVectorReal<NDIM, double> sol_vec(d_object_name + "::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec.addComponent(d_P_var, d_P_scratch_idx, wgt_cc_idx, d_hier_cc_data_ops);

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
    scratch_idxs.setFlag(d_P_scratch_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(scratch_idxs, d_integrator_time);
    }

    // Setup the right-hand-side vector for the projection-Poisson solve.
    d_hier_math_ops->div(d_Div_U_idx,
                         d_Div_U_var,
                         -1.0,
                         d_U_current_idx,
                         d_U_var,
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
        plog << d_object_name << "::regridProjection(): regrid projection solve number of iterations = "
             << regrid_projection_solver->getNumIterations() << "\n";
    if (d_enable_logging)
        plog << d_object_name << "::regridProjection(): regrid projection solve residual norm        = "
             << regrid_projection_solver->getResidualNorm() << "\n";

    // Fill ghost cells for Phi, compute Grad Phi, and set U := U - Grad Phi
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent Phi_bc_component(d_P_scratch_idx,
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
    d_hier_math_ops->grad(d_U_current_idx,
                          d_U_var,
                          /*synch_cf_bdry*/ true,
                          -1.0,
                          d_P_scratch_idx,
                          d_P_var,
                          d_no_fill_op,
                          d_integrator_time,
                          +1.0,
                          d_U_current_idx,
                          d_U_var);

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(scratch_idxs);
    }
    return;
} // regridProjection

double
INSStaggeredHierarchyIntegrator::getStableTimestep(Pointer<Patch<NDIM> > patch) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const Index<NDIM>& ilower = patch->getBox().lower();
    const Index<NDIM>& iupper = patch->getBox().upper();

    Pointer<SideData<NDIM, double> > U_data = patch->getPatchData(d_U_var, getCurrentContext());
    const IntVector<NDIM>& U_ghost_cells = U_data->getGhostCellWidth();

    double stable_dt = std::numeric_limits<double>::max();
    NAVIER_STOKES_SC_STABLEDT_FC(dx,
#if (NDIM == 2)
                                 ilower(0),
                                 iupper(0),
                                 ilower(1),
                                 iupper(1),
                                 U_ghost_cells(0),
                                 U_ghost_cells(1),
                                 U_data->getPointer(0),
                                 U_data->getPointer(1),
#endif
#if (NDIM == 3)
                                 ilower(0),
                                 iupper(0),
                                 ilower(1),
                                 iupper(1),
                                 ilower(2),
                                 iupper(2),
                                 U_ghost_cells(0),
                                 U_ghost_cells(1),
                                 U_ghost_cells(2),
                                 U_data->getPointer(0),
                                 U_data->getPointer(1),
                                 U_data->getPointer(2),
#endif
                                 stable_dt);
    return stable_dt;
} // getStableTimestep

/////////////////////////////// PRIVATE //////////////////////////////////////

void
INSStaggeredHierarchyIntegrator::reinitializeOperatorsAndSolvers(const double current_time, const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const int wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
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
    P_problem_coefs.setDConstant(rho == 0.0 ? -1.0 : -1.0 / rho);

    // Ensure that solver components are appropriately reinitialized when the
    // time step size changes.
    const bool dt_change = initial_time || !MathUtilities<double>::equalEps(dt, d_dt_previous[0]);
    if (dt_change)
    {
        d_velocity_solver_needs_init = true;
        d_stokes_solver_needs_init = true;
    }

    // Setup solver vectors.
    const bool has_velocity_nullspace = d_normalize_velocity && MathUtilities<double>::equalEps(rho, 0.0);
    const bool has_pressure_nullspace = d_normalize_pressure;
    if (d_vectors_need_init)
    {
        d_U_scratch_vec =
            new SAMRAIVectorReal<NDIM, double>(d_object_name + "::U_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_U_scratch_vec->addComponent(d_U_var, d_U_scratch_idx, wgt_sc_idx, d_hier_sc_data_ops);

        d_P_scratch_vec =
            new SAMRAIVectorReal<NDIM, double>(d_object_name + "::P_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_P_scratch_vec->addComponent(d_P_var, d_P_scratch_idx, wgt_cc_idx, d_hier_cc_data_ops);

        if (d_U_rhs_vec) d_U_rhs_vec->freeVectorComponents();
        if (d_U_adv_vec) d_U_adv_vec->freeVectorComponents();
        if (d_N_vec) d_N_vec->freeVectorComponents();
        if (d_P_rhs_vec) d_P_rhs_vec->freeVectorComponents();

        d_U_rhs_vec = d_U_scratch_vec->cloneVector(d_object_name + "::U_rhs_vec");
        d_U_adv_vec = d_U_scratch_vec->cloneVector(d_object_name + "::U_adv_vec");
        d_N_vec = d_U_scratch_vec->cloneVector(d_object_name + "::N_vec");
        d_P_rhs_vec = d_P_scratch_vec->cloneVector(d_object_name + "::P_rhs_vec");

        d_sol_vec =
            new SAMRAIVectorReal<NDIM, double>(d_object_name + "::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_sol_vec->addComponent(d_U_var, d_U_scratch_idx, wgt_sc_idx, d_hier_sc_data_ops);
        d_sol_vec->addComponent(d_P_var, d_P_scratch_idx, wgt_cc_idx, d_hier_cc_data_ops);

        d_rhs_vec =
            new SAMRAIVectorReal<NDIM, double>(d_object_name + "::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
        const int U_rhs_idx = d_U_rhs_vec->getComponentDescriptorIndex(0);
        d_rhs_vec->addComponent(d_U_var, U_rhs_idx, wgt_sc_idx, d_hier_sc_data_ops);
        const int P_rhs_idx = d_P_rhs_vec->getComponentDescriptorIndex(0);
        d_rhs_vec->addComponent(d_P_var, P_rhs_idx, wgt_cc_idx, d_hier_cc_data_ops);

        for (unsigned int k = 0; k < d_nul_vecs.size(); ++k)
        {
            if (d_nul_vecs[k]) d_nul_vecs[k]->freeVectorComponents();
        }
        const int n_nul_vecs = (has_pressure_nullspace ? 1 : 0) + (has_velocity_nullspace ? NDIM : 0);
        d_nul_vecs.resize(n_nul_vecs);

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
                d_nul_vecs[k] = d_sol_vec->cloneVector(d_object_name + "::nul_vec_U_" + stream.str());
                d_nul_vecs[k]->allocateVectorData(current_time);
                d_nul_vecs[k]->setToScalar(0.0);
                d_U_nul_vecs[k] = d_U_scratch_vec->cloneVector(d_object_name + "::U_nul_vec_U_" + stream.str());
                d_U_nul_vecs[k]->allocateVectorData(current_time);
                d_U_nul_vecs[k]->setToScalar(0.0);
                for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
                {
                    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        Pointer<SideData<NDIM, double> > nul_data =
                            patch->getPatchData(d_nul_vecs[k]->getComponentDescriptorIndex(0));
                        nul_data->getArrayData(k).fillAll(1.0);
                        Pointer<SideData<NDIM, double> > U_nul_data =
                            patch->getPatchData(d_U_nul_vecs[k]->getComponentDescriptorIndex(0));
                        U_nul_data->getArrayData(k).fillAll(1.0);
                    }
                }
            }
        }

        if (has_pressure_nullspace)
        {
            d_nul_vecs.back() = d_sol_vec->cloneVector(d_object_name + "::nul_vec_p");
            d_nul_vecs.back()->allocateVectorData(current_time);
            d_hier_sc_data_ops->setToScalar(d_nul_vecs.back()->getComponentDescriptorIndex(0), 0.0);
            d_hier_cc_data_ops->setToScalar(d_nul_vecs.back()->getComponentDescriptorIndex(1), 1.0);
        }

        d_vectors_need_init = false;
    }

    // Setup boundary conditions objects.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        INSStaggeredVelocityBcCoef* U_bc_coef = dynamic_cast<INSStaggeredVelocityBcCoef*>(d_U_bc_coefs[d]);
        U_bc_coef->setStokesSpecifications(&d_problem_coefs);
        U_bc_coef->setPhysicalBcCoefs(d_bc_coefs);
        U_bc_coef->setSolutionTime(new_time);
        U_bc_coef->setTimeInterval(current_time, new_time);
    }
    INSStaggeredPressureBcCoef* P_bc_coef = dynamic_cast<INSStaggeredPressureBcCoef*>(d_P_bc_coef);
    P_bc_coef->setStokesSpecifications(&d_problem_coefs);
    P_bc_coef->setPhysicalBcCoefs(d_bc_coefs);
    P_bc_coef->setSolutionTime(new_time);
    P_bc_coef->setTimeInterval(current_time, new_time);
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
        d_convective_op->setSolutionTime(d_integrator_time);
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
        if (d_velocity_solver_needs_init)
        {
            if (d_enable_logging)
                plog << d_object_name << "::preprocessIntegrateHierarchy(): initializing "
                                         "velocity subdomain solver"
                     << std::endl;
            LinearSolver* p_velocity_solver = dynamic_cast<LinearSolver*>(d_velocity_solver.getPointer());
            if (p_velocity_solver)
            {
                p_velocity_solver->setInitialGuessNonzero(false);
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
        if (d_pressure_solver_needs_init)
        {
            if (d_enable_logging)
                plog << d_object_name << "::preprocessIntegrateHierarchy(): initializing "
                                         "pressure subdomain solver"
                     << std::endl;
            LinearSolver* p_pressure_solver = dynamic_cast<LinearSolver*>(d_pressure_solver.getPointer());
            if (p_pressure_solver)
            {
                p_pressure_solver->setInitialGuessNonzero(false);
                if (has_pressure_nullspace) p_pressure_solver->setNullspace(true);
            }
            d_pressure_solver->initializeSolverState(*d_P_scratch_vec, *d_P_rhs_vec);
            d_pressure_solver_needs_init = false;
        }
    }

    // Setup Stokes solver.
    d_stokes_solver->setVelocityPoissonSpecifications(U_problem_coefs);
    d_stokes_solver->setPhysicalBcCoefs(d_U_bc_coefs, d_P_bc_coef);
    d_stokes_solver->setPhysicalBoundaryHelper(d_bc_helper);
    d_stokes_solver->setSolutionTime(new_time);
    d_stokes_solver->setTimeInterval(current_time, new_time);
    d_stokes_solver->setComponentsHaveNullspace(has_velocity_nullspace, has_pressure_nullspace);
    LinearSolver* p_stokes_linear_solver = dynamic_cast<LinearSolver*>(d_stokes_solver.getPointer());
    if (!p_stokes_linear_solver)
    {
        NewtonKrylovSolver* p_stokes_newton_solver = dynamic_cast<NewtonKrylovSolver*>(d_stokes_solver.getPointer());
        if (p_stokes_newton_solver) p_stokes_linear_solver = p_stokes_newton_solver->getLinearSolver().getPointer();
    }
    if (p_stokes_linear_solver)
    {
        StaggeredStokesBlockPreconditioner* p_stokes_block_pc =
            dynamic_cast<StaggeredStokesBlockPreconditioner*>(p_stokes_linear_solver);
        StaggeredStokesFACPreconditioner* p_stokes_fac_pc =
            dynamic_cast<StaggeredStokesFACPreconditioner*>(p_stokes_linear_solver);
        if (!(p_stokes_block_pc || p_stokes_fac_pc))
        {
            KrylovLinearSolver* p_stokes_krylov_solver = dynamic_cast<KrylovLinearSolver*>(p_stokes_linear_solver);
            if (p_stokes_krylov_solver)
            {
                p_stokes_block_pc = dynamic_cast<StaggeredStokesBlockPreconditioner*>(
                    p_stokes_krylov_solver->getPreconditioner().getPointer());

                p_stokes_fac_pc = dynamic_cast<StaggeredStokesFACPreconditioner*>(
                    p_stokes_krylov_solver->getPreconditioner().getPointer());
            }
        }
        if (p_stokes_block_pc)
        {
            p_stokes_block_pc->setPressurePoissonSpecifications(P_problem_coefs);
            p_stokes_block_pc->setPhysicalBcCoefs(d_U_star_bc_coefs, d_Phi_bc_coef);
            p_stokes_block_pc->setComponentsHaveNullspace(has_velocity_nullspace, has_pressure_nullspace);
        }
        else if (p_stokes_fac_pc)
        {
            p_stokes_fac_pc->setPhysicalBcCoefs(d_U_star_bc_coefs, d_Phi_bc_coef);
            p_stokes_fac_pc->setComponentsHaveNullspace(has_velocity_nullspace, has_pressure_nullspace);
        }
        else
        {
            TBOX_WARNING("No special BCs set for the preconditioner \n");
        }
    }
    if (d_stokes_solver_needs_init)
    {
        if (d_enable_logging)
            plog << d_object_name << "::preprocessIntegrateHierarchy(): initializing "
                                     "incompressible Stokes solver"
                 << std::endl;
        if (p_stokes_linear_solver)
        {
            p_stokes_linear_solver->setInitialGuessNonzero(true);
            if (has_velocity_nullspace || has_pressure_nullspace)
                p_stokes_linear_solver->setNullspace(false, d_nul_vecs);
        }
        d_stokes_solver->initializeSolverState(*d_sol_vec, *d_rhs_vec);
        d_stokes_solver_needs_init = false;
    }
    return;
} // reinitializeOperatorsAndSolvers

void
INSStaggeredHierarchyIntegrator::computeDivSourceTerm(const int F_idx, const int Q_idx, const int U_idx)
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

            Pointer<SideData<NDIM, double> > U_data = patch->getPatchData(U_idx);
            Pointer<CellData<NDIM, double> > Q_data = patch->getPatchData(Q_idx);
            Pointer<SideData<NDIM, double> > F_data = patch->getPatchData(F_idx);

            const IntVector<NDIM>& U_data_gc = U_data->getGhostCellWidth();
            const IntVector<NDIM>& Q_data_gc = Q_data->getGhostCellWidth();
            const IntVector<NDIM>& F_data_gc = F_data->getGhostCellWidth();

            switch (d_convective_op->getConvectiveDifferencingType())
            {
            case CONSERVATIVE:
                NAVIER_STOKES_STAGGERED_CONS_SOURCE_FC(
#if (NDIM == 2)
                    ilower(0),
                    iupper(0),
                    ilower(1),
                    iupper(1),
                    U_data_gc(0),
                    U_data_gc(1),
                    Q_data_gc(0),
                    Q_data_gc(1),
                    F_data_gc(0),
                    F_data_gc(1),
                    U_data->getPointer(0),
                    U_data->getPointer(1),
                    Q_data->getPointer(),
                    F_data->getPointer(0),
                    F_data->getPointer(1)
#endif
#if (NDIM == 3)
                        ilower(0),
                    iupper(0),
                    ilower(1),
                    iupper(1),
                    ilower(2),
                    iupper(2),
                    U_data_gc(0),
                    U_data_gc(1),
                    U_data_gc(2),
                    Q_data_gc(0),
                    Q_data_gc(1),
                    Q_data_gc(2),
                    F_data_gc(0),
                    F_data_gc(1),
                    F_data_gc(2),
                    U_data->getPointer(0),
                    U_data->getPointer(1),
                    U_data->getPointer(2),
                    Q_data->getPointer(),
                    F_data->getPointer(0),
                    F_data->getPointer(1),
                    F_data->getPointer(2)
#endif
                        );
                break;
            case ADVECTIVE:
                NAVIER_STOKES_STAGGERED_ADV_SOURCE_FC(
#if (NDIM == 2)
                    ilower(0),
                    iupper(0),
                    ilower(1),
                    iupper(1),
                    U_data_gc(0),
                    U_data_gc(1),
                    Q_data_gc(0),
                    Q_data_gc(1),
                    F_data_gc(0),
                    F_data_gc(1),
                    U_data->getPointer(0),
                    U_data->getPointer(1),
                    Q_data->getPointer(),
                    F_data->getPointer(0),
                    F_data->getPointer(1)
#endif
#if (NDIM == 3)
                        ilower(0),
                    iupper(0),
                    ilower(1),
                    iupper(1),
                    ilower(2),
                    iupper(2),
                    U_data_gc(0),
                    U_data_gc(1),
                    U_data_gc(2),
                    Q_data_gc(0),
                    Q_data_gc(1),
                    Q_data_gc(2),
                    F_data_gc(0),
                    F_data_gc(1),
                    F_data_gc(2),
                    U_data->getPointer(0),
                    U_data->getPointer(1),
                    U_data->getPointer(2),
                    Q_data->getPointer(),
                    F_data->getPointer(0),
                    F_data->getPointer(1),
                    F_data->getPointer(2)
#endif
                        );
                break;
            case SKEW_SYMMETRIC:
                NAVIER_STOKES_STAGGERED_SKEW_SYM_SOURCE_FC(
#if (NDIM == 2)
                    ilower(0),
                    iupper(0),
                    ilower(1),
                    iupper(1),
                    U_data_gc(0),
                    U_data_gc(1),
                    Q_data_gc(0),
                    Q_data_gc(1),
                    F_data_gc(0),
                    F_data_gc(1),
                    U_data->getPointer(0),
                    U_data->getPointer(1),
                    Q_data->getPointer(),
                    F_data->getPointer(0),
                    F_data->getPointer(1)
#endif
#if (NDIM == 3)
                        ilower(0),
                    iupper(0),
                    ilower(1),
                    iupper(1),
                    ilower(2),
                    iupper(2),
                    U_data_gc(0),
                    U_data_gc(1),
                    U_data_gc(2),
                    Q_data_gc(0),
                    Q_data_gc(1),
                    Q_data_gc(2),
                    F_data_gc(0),
                    F_data_gc(1),
                    F_data_gc(2),
                    U_data->getPointer(0),
                    U_data->getPointer(1),
                    U_data->getPointer(2),
                    Q_data->getPointer(),
                    F_data->getPointer(0),
                    F_data->getPointer(1),
                    F_data->getPointer(2)
#endif
                        );
                break;
            default:
                TBOX_ERROR(
                    "INSStaggeredHierarchyIntegrator::computeDivSourceTerm():\n"
                    << "  unsupported differencing form: "
                    << enum_to_string<ConvectiveDifferencingType>(d_convective_op->getConvectiveDifferencingType())
                    << " \n"
                    << "  valid choices are: ADVECTIVE, CONSERVATIVE, SKEW_SYMMETRIC\n");
            }
        }
    }
    return;
} // computeDivSourceTerm

TimeSteppingType
INSStaggeredHierarchyIntegrator::getConvectiveTimeSteppingType(const int cycle_num)
{
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
    return convective_time_stepping_type;
} // getConvectiveTimeSteppingType

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
