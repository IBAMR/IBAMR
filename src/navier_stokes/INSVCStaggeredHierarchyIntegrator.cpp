// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/BrinkmanPenalizationStrategy.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSIntermediateVelocityBcCoef.h"
#include "ibamr/INSProjectionBcCoef.h"
#include "ibamr/INSStaggeredConvectiveOperatorManager.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibamr/INSVCStaggeredPressureBcCoef.h"
#include "ibamr/INSVCStaggeredVelocityBcCoef.h"
#include "ibamr/PETScKrylovStaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesBlockPreconditioner.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/VCStaggeredStokesOperator.h"
#include "ibamr/VCStaggeredStokesProjectionPreconditioner.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/CartSideDoubleDivPreservingRefine.h"
#include "ibtk/CartSideDoubleRT0Refine.h"
#include "ibtk/CartSideDoubleSpecializedLinearRefine.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/NewtonKrylovSolver.h"
#include "ibtk/PETScKrylovLinearSolver.h"
#include "ibtk/PETScKrylovPoissonSolver.h"
#include "ibtk/PoissonFACPreconditioner.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/SCPoissonSolverManager.h"
#include "ibtk/SideDataSynchronization.h"
#include "ibtk/VCSCViscousOpPointRelaxationFACOperator.h"
#include "ibtk/VCSCViscousOperator.h"
#include "ibtk/VCSCViscousPETScLevelSolver.h"
#include "ibtk/ibtk_enums.h"

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
#include "ComponentSelector.h"
#include "EdgeVariable.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "GriddingAlgorithm.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "HierarchyEdgeDataOpsReal.h"
#include "HierarchyFaceDataOpsReal.h"
#include "HierarchyNodeDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "NodeVariable.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchSideDataOpsReal.h"
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
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

// FORTRAN ROUTINES
#if (NDIM == 2)
#define NAVIER_STOKES_SC_STABLEDT_FC IBAMR_FC_FUNC_(navier_stokes_sc_stabledt2d, NAVIER_STOKES_SC_STABLEDT2D)
#define NAVIER_STOKES_SIDE_TO_FACE_FC IBAMR_FC_FUNC_(navier_stokes_side_to_face2d, NAVIER_STOKES_SIDE_TO_FACE2D)
#endif

#if (NDIM == 3)
#define NAVIER_STOKES_SC_STABLEDT_FC IBAMR_FC_FUNC_(navier_stokes_sc_stabledt3d, NAVIER_STOKES_SC_STABLEDT3D)
#define NAVIER_STOKES_SIDE_TO_FACE_FC IBAMR_FC_FUNC_(navier_stokes_side_to_face3d, NAVIER_STOKES_SIDE_TO_FACE3D)
#endif

extern "C"
{
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
static const int NODEG = 1;
static const int EDGEG = 1;
static const int MUCELLG = 2;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Default solvers and preconditioners to register with the solver managers
static const std::string DEFAULT_VC_STAGGERED_STOKES_SOLVER = "VC_STAGGERED_STOKES_PETSC_KRYLOV_SOLVER";
static const std::string DEFAULT_VC_STAGGERED_STOKES_PRECOND = "VC_STAGGERED_STOKES_PROJECTION_PRECONDITIONER";
static const std::string DEFAULT_VC_VELOCITY_SOLVER = "VC_VELOCITY_PETSC_KRYLOV_SOLVER";
static const std::string DEFAULT_VC_VELOCITY_PRECOND = "VC_VELOCITY_POINT_RELAXATION_FAC_PRECONDITIONER";
static const std::string DEFAULT_VC_VELOCITY_LEVEL_SOLVER = "VC_VELOCITY_PETSC_LEVEL_SOLVER";

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
            const hier::Index<NDIM>& ilower = patch->getBox().lower();
            const hier::Index<NDIM>& iupper = patch->getBox().upper();
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

Pointer<StaggeredStokesSolver>
allocate_vc_stokes_krylov_solver(const std::string& solver_object_name,
                                 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> solver_input_db,
                                 const std::string& solver_default_options_prefix)
{
    Pointer<PETScKrylovStaggeredStokesSolver> krylov_solver =
        new PETScKrylovStaggeredStokesSolver(solver_object_name, solver_input_db, solver_default_options_prefix);
    krylov_solver->setOperator(new VCStaggeredStokesOperator(solver_object_name + "::vc_staggered_stokes_operator"));
    return krylov_solver;
}

Pointer<PoissonSolver>
allocate_vc_velocity_krylov_solver(const std::string& solver_object_name,
                                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> solver_input_db,
                                   const std::string& solver_default_options_prefix)
{
    Pointer<PETScKrylovPoissonSolver> krylov_solver =
        new PETScKrylovPoissonSolver(solver_object_name, solver_input_db, solver_default_options_prefix);
    krylov_solver->setOperator(new VCSCViscousOperator(solver_object_name + "::vc_viscous_operator"));
    return krylov_solver;
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSVCStaggeredHierarchyIntegrator::INSVCStaggeredHierarchyIntegrator(std::string object_name,
                                                                     Pointer<Database> input_db,
                                                                     bool register_for_restart)
    : INSHierarchyIntegrator(std::move(object_name),
                             input_db,
                             new SideVariable<NDIM, double>(object_name + "::U"),
                             "CONSERVATIVE_COARSEN",
                             "CONSERVATIVE_LINEAR_REFINE",
                             new CellVariable<NDIM, double>(object_name + "::P"),
                             "CONSERVATIVE_COARSEN",
                             "LINEAR_REFINE",
                             new SideVariable<NDIM, double>(object_name + "::F"),
                             "CONSERVATIVE_COARSEN",
                             "CONSERVATIVE_LINEAR_REFINE",
                             new CellVariable<NDIM, double>(object_name + "::Q"),
                             "CONSERVATIVE_COARSEN",
                             "CONSTANT_REFINE",
                             register_for_restart),
      d_rho_vc_interp_type(VC_HARMONIC_INTERP),
      d_mu_vc_interp_type(VC_HARMONIC_INTERP)
{
    // Get plotting options from database
    if (input_db->keyExists("rho_scale")) d_rho_scale = input_db->getDouble("rho_scale");
    if (input_db->keyExists("mu_scale")) d_mu_scale = input_db->getDouble("mu_scale");
    if (input_db->keyExists("output_rho")) d_output_rho = input_db->getBool("output_rho");
    if (input_db->keyExists("output_mu")) d_output_mu = input_db->getBool("output_mu");

    // Register solver factory functions for variable coefficient Stokes and
    // viscous solvers
    StaggeredStokesSolverManager::getManager()->registerSolverFactoryFunction(DEFAULT_VC_STAGGERED_STOKES_SOLVER,
                                                                              allocate_vc_stokes_krylov_solver);
    StaggeredStokesSolverManager::getManager()->registerSolverFactoryFunction(
        DEFAULT_VC_STAGGERED_STOKES_PRECOND, VCStaggeredStokesProjectionPreconditioner::allocate_solver);
    SCPoissonSolverManager::getManager()->registerSolverFactoryFunction(DEFAULT_VC_VELOCITY_SOLVER,
                                                                        allocate_vc_velocity_krylov_solver);
    SCPoissonSolverManager::getManager()->registerSolverFactoryFunction(
        DEFAULT_VC_VELOCITY_PRECOND, VCSCViscousOpPointRelaxationFACOperator::allocate_solver);
    SCPoissonSolverManager::getManager()->registerSolverFactoryFunction(DEFAULT_VC_VELOCITY_LEVEL_SOLVER,
                                                                        VCSCViscousPETScLevelSolver::allocate_solver);

    // Check to see whether the solver types have been set.
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

    // Check if density and/or viscosity is constant
    if (input_db->keyExists("rho_is_const")) d_rho_is_const = input_db->getBool("rho_is_const");
    if (input_db->keyExists("mu_is_const")) d_mu_is_const = input_db->getBool("mu_is_const");
    if (d_rho_is_const && d_mu_is_const)
    {
        TBOX_ERROR(d_object_name << "::INSVCStaggeredHierarchyIntegrator():\n"
                                 << " for constant coefficient problems,\n"
                                 << " variable coefficient integrator is less\n"
                                 << " efficient than constant coefficient version.\n"
                                 << " Use INSStaggeredHierarchyIntegrator instead.");
    }

    // Get the interpolation type for the material properties
    if (input_db->keyExists("vc_interpolation_type"))
    {
        d_rho_vc_interp_type = IBTK::string_to_enum<VCInterpType>(input_db->getString("vc_interpolation_type"));
        d_mu_vc_interp_type = IBTK::string_to_enum<VCInterpType>(input_db->getString("vc_interpolation_type"));
    }
    if (input_db->keyExists("rho_vc_interpolation_type"))
    {
        d_rho_vc_interp_type = IBTK::string_to_enum<VCInterpType>(input_db->getString("rho_vc_interpolation_type"));
    }
    if (input_db->keyExists("mu_vc_interpolation_type"))
    {
        d_mu_vc_interp_type = IBTK::string_to_enum<VCInterpType>(input_db->getString("mu_vc_interpolation_type"));
    }
    switch (d_rho_vc_interp_type)
    {
    case VC_HARMONIC_INTERP:
    case VC_AVERAGE_INTERP:
        break;
    default:
        TBOX_ERROR(d_object_name << "::INSVCStaggeredHierarchyIntegrator():\n"
                                 << "  unsupported density interpolation type: "
                                 << IBTK::enum_to_string<VCInterpType>(d_rho_vc_interp_type) << " \n"
                                 << "  valid choices are: VC_HARMONIC_INTERP, VC_AVERAGE_INTERP\n");
    }
    switch (d_mu_vc_interp_type)
    {
    case VC_HARMONIC_INTERP:
    case VC_AVERAGE_INTERP:
        break;
    default:
        TBOX_ERROR(d_object_name << "::INSVCStaggeredHierarchyIntegrator():\n"
                                 << "  unsupported viscosity interpolation type: "
                                 << IBTK::enum_to_string<VCInterpType>(d_mu_vc_interp_type) << " \n"
                                 << "  valid choices are: VC_HARMONIC_INTERP, VC_AVERAGE_INTERP\n");
    }

    // Get the scaling coefficients
    if (input_db->keyExists("operator_scale_factors"))
    {
        d_A_scale = input_db->getDoubleArray("operator_scale_factors");
        for (int k = 0; k < d_A_scale.size(); ++k)
        {
            if (d_A_scale[k] <= 0.0)
            {
                TBOX_ERROR(d_object_name << "::INSVCStaggeredHierarchyIntegrator():\n"
                                         << " scaling for improving the condition number of\n"
                                         << " the Stokes system must be positive.\n"
                                         << " operator_scale_factors[" << k << "] = " << d_A_scale[k] << "\n");
            }
        }
    }
    else
    {
        // The default is no scaling, which will be applied to all levels of the
        // patch hierarchy
        d_A_scale.resizeArray(1);
        d_A_scale[0] = 1.0;
    }

    // By default, reinitialize the preconditioner every time step
    if (input_db->keyExists("precond_reinit_interval"))
        d_precond_reinit_interval = input_db->getInteger("precond_reinit_interval");
    if (d_precond_reinit_interval <= 0)
    {
        TBOX_ERROR(d_object_name << "::INSVCStaggeredHierarchyIntegrator():\n"
                                 << " preconditioner reinitialization interval\n"
                                 << " must be a positive integer");
    }

    // Check to make sure the time stepping types are supported.
    switch (d_viscous_time_stepping_type)
    {
    case BACKWARD_EULER:
    case FORWARD_EULER:
    case TRAPEZOIDAL_RULE:
        break;
    default:
        TBOX_ERROR(d_object_name << "::INSVCStaggeredHierarchyIntegrator():\n"
                                 << "  unsupported viscous time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_viscous_time_stepping_type) << " \n"
                                 << "  valid choices are: BACKWARD_EULER, FORWARD_EULER, "
                                    "TRAPEZOIDAL_RULE\n");
    }

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

    // Flag to determine whether we explicitly remove any null space components.
    if (input_db->keyExists("explicitly_remove_nullspace"))
        d_explicitly_remove_nullspace = input_db->getBool("explicitly_remove_nullspace");

    // Setup physical boundary conditions objects.
    d_bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
    d_U_bc_coefs.resize(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_U_bc_coefs[d] = new INSVCStaggeredVelocityBcCoef(d, this, d_bc_coefs, d_traction_bc_type);
    }
    d_P_bc_coef = new INSVCStaggeredPressureBcCoef(this, d_bc_coefs, d_traction_bc_type);

    // Get coarsen and refine operator types.
    if (input_db->keyExists("N_coarsen_type")) d_N_coarsen_type = input_db->getString("N_coarsen_type");
    if (input_db->keyExists("N_refine_type")) d_N_refine_type = input_db->getString("N_refine_type");

    if (input_db->keyExists("mu_coarsen_type")) d_mu_coarsen_type = input_db->getString("mu_coarsen_type");
    if (input_db->keyExists("mu_refine_type")) d_mu_refine_type = input_db->getString("mu_refine_type");
    if (input_db->keyExists("mu_bdry_extrap_type")) d_mu_bdry_extrap_type = input_db->getString("mu_bdry_extrap_type");

    if (input_db->keyExists("rho_coarsen_type")) d_rho_coarsen_type = input_db->getString("rho_coarsen_type");
    if (input_db->keyExists("rho_refine_type")) d_rho_refine_type = input_db->getString("rho_refine_type");
    if (input_db->keyExists("rho_bdry_extrap_type"))
        d_rho_bdry_extrap_type = input_db->getString("rho_bdry_extrap_type");

    // Initialize all variables.  The velocity, pressure, body force, and fluid
    // source variables were created above in the constructor for the
    // INSHierarchyIntegrator base class.
    d_U_var = INSHierarchyIntegrator::d_U_var;
    d_P_var = INSHierarchyIntegrator::d_P_var;
    d_F_var = INSHierarchyIntegrator::d_F_var;
    d_Q_var = INSHierarchyIntegrator::d_Q_var;
    d_N_old_var = new SideVariable<NDIM, double>(d_object_name + "::N_old");
    d_U_old_var = new SideVariable<NDIM, double>(d_object_name + "::U_old");

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
    d_EE_var = new CellVariable<NDIM, double>(d_object_name + "::EE", NDIM * NDIM);

    return;
} // INSVCStaggeredHierarchyIntegrator

INSVCStaggeredHierarchyIntegrator::~INSVCStaggeredHierarchyIntegrator()
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        delete d_U_bc_coefs[d];
        d_U_bc_coefs[d] = nullptr;
    }
    delete d_P_bc_coef;
    d_P_bc_coef = nullptr;
    d_velocity_solver.setNull();
    d_pressure_solver.setNull();
    if (d_U_rhs_vec) d_U_rhs_vec->freeVectorComponents();
    if (d_U_adv_vec) d_U_adv_vec->freeVectorComponents();
    if (d_N_vec) d_N_vec->freeVectorComponents();
    if (d_P_rhs_vec) d_P_rhs_vec->freeVectorComponents();
    for (const auto& nul_vec : d_nul_vecs)
    {
        if (nul_vec) nul_vec->freeVectorComponents();
    }
    for (const auto& U_nul_vec : d_U_nul_vecs)
    {
        if (U_nul_vec) U_nul_vec->freeVectorComponents();
    }

    return;
} // ~INSVCStaggeredHierarchyIntegrator

Pointer<ConvectiveOperator>
INSVCStaggeredHierarchyIntegrator::getConvectiveOperator()
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
INSVCStaggeredHierarchyIntegrator::getVelocitySubdomainSolver()
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
INSVCStaggeredHierarchyIntegrator::getPressureSubdomainSolver()
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
INSVCStaggeredHierarchyIntegrator::setStokesSolver(Pointer<StaggeredStokesSolver> stokes_solver)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_stokes_solver);
#endif
    d_stokes_solver = stokes_solver;
    d_stokes_solver_needs_init = true;
    return;
} // setStokesSolver

Pointer<StaggeredStokesSolver>
INSVCStaggeredHierarchyIntegrator::getStokesSolver()
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
INSVCStaggeredHierarchyIntegrator::setStokesSolverNeedsInit()
{
    d_stokes_solver_needs_init = true;
    return;
}

void
INSVCStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                                 Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;
    const int max_levels = gridding_alg->getMaxLevels();

    // Setup the condition number scaling array.
    // Extra entries will be stripped, while unset values will be set to that of
    // the next coarsest level.
    const int orig_size = d_A_scale.size();
    const int size_difference = max_levels - orig_size;
    d_A_scale.resizeArray(max_levels);
    if (size_difference > 0)
    {
        const double last_scale = d_A_scale[orig_size - 1];
        for (int k = orig_size; k < max_levels; ++k)
        {
            d_A_scale[k] = last_scale;
        }
    }

    // Setup solvers.
    if (d_stokes_solver_type == StaggeredStokesSolverManager::UNDEFINED)
    {
        d_stokes_solver_type = DEFAULT_VC_STAGGERED_STOKES_SOLVER;
        d_stokes_solver_db->putString("ksp_type", "fgmres");
    }

    if (d_stokes_precond_type == StaggeredStokesSolverManager::UNDEFINED)
    {
        d_stokes_precond_type = DEFAULT_VC_STAGGERED_STOKES_PRECOND;
        d_stokes_precond_db->putInteger("max_iterations", 1);
    }

    if (d_velocity_solver_type == SCPoissonSolverManager::UNDEFINED)
    {
        d_velocity_solver_type = DEFAULT_VC_VELOCITY_SOLVER;
        d_velocity_solver_db->putString("ksp_type", "richardson");
        d_velocity_solver_db->putInteger("max_iterations", 10);
        d_velocity_solver_db->putDouble("rel_residual_tol", 1.0e-1);
    }

    if (d_velocity_precond_type == SCPoissonSolverManager::UNDEFINED)
    {
        if (max_levels == 1)
        {
            d_velocity_precond_type = DEFAULT_VC_VELOCITY_LEVEL_SOLVER;
        }
        else
        {
            d_velocity_precond_type = DEFAULT_VC_VELOCITY_PRECOND;
            d_velocity_precond_db->putString("coarse_solver_type", DEFAULT_VC_VELOCITY_LEVEL_SOLVER);
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
    d_hier_nc_data_ops =
        hier_ops_manager->getOperationsDouble(new NodeVariable<NDIM, double>("nc_var"), hierarchy, true);
    d_hier_ec_data_ops =
        hier_ops_manager->getOperationsDouble(new EdgeVariable<NDIM, double>("ec_var"), hierarchy, true);
    d_hier_math_ops = buildHierarchyMathOps(d_hierarchy);

    // Register state variables that are maintained by the
    // INSVCStaggeredHierarchyIntegrator.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    grid_geom->addSpatialRefineOperator(new CartSideDoubleRT0Refine());
    grid_geom->addSpatialRefineOperator(new CartSideDoubleSpecializedLinearRefine());

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> side_ghosts = SIDEG;
    const IntVector<NDIM> node_ghosts = NODEG;
    const IntVector<NDIM> edge_ghosts = EDGEG;
    const IntVector<NDIM> no_ghosts = 0;
    const IntVector<NDIM> mu_cell_ghosts = MUCELLG;

    registerVariable(d_U_old_current_idx,
                     d_U_old_new_idx,
                     d_U_old_scratch_idx,
                     d_U_old_var,
                     side_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE",
                     d_U_init);

    registerVariable(d_U_current_idx,
                     d_U_new_idx,
                     d_U_scratch_idx,
                     d_U_var,
                     side_ghosts,
                     d_U_coarsen_type,
                     d_U_refine_type,
                     d_U_init);

    registerVariable(d_P_current_idx,
                     d_P_new_idx,
                     d_P_scratch_idx,
                     d_P_var,
                     cell_ghosts,
                     d_P_coarsen_type,
                     d_P_refine_type,
                     d_P_init);

    if (d_F_fcn)
    {
        registerVariable(d_F_current_idx,
                         d_F_new_idx,
                         d_F_scratch_idx,
                         d_F_var,
                         side_ghosts,
                         d_F_coarsen_type,
                         d_F_refine_type,
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
                         d_Q_coarsen_type,
                         d_Q_refine_type,
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
                     d_N_coarsen_type,
                     d_N_refine_type);

    // Get the viscosity variable, which can either be an advected field
    // maintained by an appropriate advection-diffusion integrator, or a set
    // field with some functional form maintained by the INS integrator
    if (!d_mu_is_const)
    {
        if (d_adv_diff_hier_integrator && d_mu_adv_diff_var)
        {
#if !defined(NDEBUG)
            // AdvDiffHierarchyIntegrator should initialize and maintain the viscosity
            // variable.
            TBOX_ASSERT(!d_mu_var);
            TBOX_ASSERT(!d_mu_init_fcn);
#endif
            d_mu_var = Pointer<CellVariable<NDIM, double> >(nullptr);
            // Ensure that boundary conditions are provided by the advection-diffusion
            // integrator
            d_mu_bc_coef = (d_adv_diff_hier_integrator->getPhysicalBcCoefs(d_mu_adv_diff_var)).front();
        }
        else if (d_mu_var)
        {
            Pointer<CellVariable<NDIM, double> > cc_var = d_mu_var;
            if (!cc_var)
            {
                TBOX_ERROR(
                    "INSVCStaggeredHierarchyIntegrator::"
                    "initializeHierarchyIntegrator():\n"
                    << " registered viscosity variable must be cell centered");
            }
        }
        else
        {
            TBOX_ERROR(
                "INSVCStaggeredHierarchyIntegrator::"
                "initializeHierarchyIntegrator():\n"
                << "  mu_is_const == false but no viscosity variable has been "
                   "registered.\n");
        }
    }

    if (d_mu_var)
    {
#if !defined(NDEBUG)
        // INSVCStaggeredHierarchyIntegrator should initialize the viscosity
        // variable.
        TBOX_ASSERT(d_mu_init_fcn || d_reset_mu_fcns.size() > 0);
#endif
        registerVariable(d_mu_current_idx,
                         d_mu_new_idx,
                         d_mu_scratch_idx,
                         d_mu_var,
                         mu_cell_ghosts,
                         d_mu_coarsen_type,
                         d_mu_refine_type,
                         d_mu_init_fcn);
    }
    else
    {
        d_mu_current_idx = -1;
        d_mu_new_idx = -1;
        d_mu_init_fcn = nullptr;

        Pointer<CellVariable<NDIM, double> > mu_cc_scratch_var =
            new CellVariable<NDIM, double>(d_object_name + "_mu_cc_scratch_var",
                                           /*depth*/ 1);
        d_mu_scratch_idx = var_db->registerVariableAndContext(mu_cc_scratch_var, getScratchContext(), mu_cell_ghosts);
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

// Register scratch variables that are maintained by the
// INSVCStaggeredHierarchyIntegrator.
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

        if (d_output_mu && !d_mu_is_const && d_mu_var)
        {
            d_visit_writer->registerPlotQuantity("mu_ins", "SCALAR", d_mu_current_idx, 0, d_mu_scale);
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

        if (d_output_EE)
        {
            registerVariable(d_EE_idx, d_EE_var, no_ghosts, getCurrentContext());
            d_visit_writer->registerPlotQuantity("EE", "TENSOR", d_EE_idx);
        }
    }

    // Initialize and register variables used to compute INS coefficients
    d_pressure_D_var = new SideVariable<NDIM, double>(d_object_name + "::pressure_D");
    d_pressure_D_idx = var_db->registerVariableAndContext(d_pressure_D_var, getCurrentContext(), side_ghosts);

    d_pressure_rhs_D_var = new SideVariable<NDIM, double>(d_object_name + "::pressure_rhs_D");
    d_pressure_rhs_D_idx = var_db->registerVariableAndContext(d_pressure_rhs_D_var, getCurrentContext(), side_ghosts);

    d_velocity_C_var = new SideVariable<NDIM, double>(d_object_name + "::velocity_C");
    d_velocity_C_idx = var_db->registerVariableAndContext(d_velocity_C_var, getCurrentContext(), no_ghosts);

    d_velocity_L_var = new SideVariable<NDIM, double>(d_object_name + "::velocity_L");
    d_velocity_L_idx = var_db->registerVariableAndContext(d_velocity_L_var, getCurrentContext(), side_ghosts);

    d_velocity_rhs_C_var = new SideVariable<NDIM, double>(d_object_name + "::velocity_rhs_C");
    d_velocity_rhs_C_idx = var_db->registerVariableAndContext(d_velocity_rhs_C_var, getCurrentContext(), no_ghosts);
#if (NDIM == 2)
    d_velocity_D_var = new NodeVariable<NDIM, double>(d_object_name + "::velocity_D");
    d_velocity_rhs_D_var = new NodeVariable<NDIM, double>(d_object_name + "::velocity_rhs_D");
#elif (NDIM == 3)
    d_velocity_D_var = new EdgeVariable<NDIM, double>(d_object_name + "::velocity_D");
    d_velocity_rhs_D_var = new EdgeVariable<NDIM, double>(d_object_name + "::velocity_rhs_D");
#endif
    d_velocity_D_idx = var_db->registerVariableAndContext(
        d_velocity_D_var, getCurrentContext(), (NDIM == 2 ? node_ghosts : edge_ghosts));
    d_velocity_rhs_D_idx = var_db->registerVariableAndContext(
        d_velocity_rhs_D_var, getCurrentContext(), (NDIM == 2 ? node_ghosts : edge_ghosts));

    d_velocity_D_cc_var = new CellVariable<NDIM, double>(d_object_name + "::velocity_D_cc");
    d_velocity_D_cc_idx = var_db->registerVariableAndContext(d_velocity_D_cc_var, getCurrentContext(), no_ghosts);

    d_temp_sc_var = new SideVariable<NDIM, double>(d_object_name + "::temp_sc");
    d_temp_sc_idx = var_db->registerVariableAndContext(d_temp_sc_var, getCurrentContext(), no_ghosts);
    d_temp_cc_var = new CellVariable<NDIM, double>(d_object_name + ":temp_cc",
                                                   /*depth*/ NDIM);
    d_temp_cc_idx = var_db->registerVariableAndContext(d_temp_cc_var, getCurrentContext(), cell_ghosts);

#if (NDIM == 2)
    d_mu_interp_var = new NodeVariable<NDIM, double>(d_object_name + "::mu_interp");
#elif (NDIM == 3)
    d_mu_interp_var = new EdgeVariable<NDIM, double>(d_object_name + "::mu_interp");
#endif
    d_mu_interp_idx =
        var_db->registerVariableAndContext(d_mu_interp_var, getCurrentContext(), NDIM == 2 ? node_ghosts : edge_ghosts);

    d_N_full_var = new SideVariable<NDIM, double>(d_object_name + "N_full");
    d_N_full_idx = var_db->registerVariableAndContext(d_N_full_var, getCurrentContext(), no_ghosts);

    // Register persistent variables to be used for boundary conditions and other
    // applications.
    // Note: these will not be deallocated.
    Pointer<CellVariable<NDIM, double> > mu_cc_linear_op_var =
        new CellVariable<NDIM, double>(d_object_name + "_mu_cc_linear_op_var",
                                       /*depth*/ 1);
    d_mu_linear_op_idx = var_db->registerVariableAndContext(
        mu_cc_linear_op_var, var_db->getContext(d_object_name + "::mu_linear_op"), mu_cell_ghosts);
    d_mu_interp_linear_op_idx =
        var_db->registerVariableAndContext(d_mu_interp_var,
                                           var_db->getContext(d_object_name + "::mu_interp_linear_op"),
                                           NDIM == 2 ? node_ghosts : edge_ghosts);
    Pointer<SideVariable<NDIM, double> > rho_sc_linear_op_var =
        new SideVariable<NDIM, double>(d_object_name + "_rho_sc_linear_op_var",
                                       /*depth*/ 1);
    d_rho_linear_op_idx = var_db->registerVariableAndContext(
        rho_sc_linear_op_var, var_db->getContext(d_object_name + "::rho_linear_op_var"), no_ghosts);

    // Setup a specialized coarsen algorithm.
    Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
    Pointer<CoarsenOperator<NDIM> > coarsen_op = grid_geom->lookupCoarsenOperator(d_U_var, d_U_coarsen_type);
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
            if (p_stokes_krylov_solver)
            {
                p_stokes_block_pc = p_stokes_krylov_solver->getPreconditioner();
                Pointer<VCStaggeredStokesOperator> p_vc_stokes_op = p_stokes_krylov_solver->getOperator();
                if (p_vc_stokes_op) p_vc_stokes_op->setDPatchDataInterpolationType(d_mu_vc_interp_type);
            }
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

    // Set the velocity subdomain solver interpolation type if necessary
    auto p_velocity_solver = dynamic_cast<IBTK::PETScKrylovLinearSolver*>(d_velocity_solver.getPointer());
    if (p_velocity_solver)
    {
        Pointer<PoissonFACPreconditioner> p_poisson_fac_pc = p_velocity_solver->getPreconditioner();
        Pointer<VCSCViscousOpPointRelaxationFACOperator> p_vc_point_fac_op =
            p_poisson_fac_pc->getFACPreconditionerStrategy();
        if (p_vc_point_fac_op) p_vc_point_fac_op->setDPatchDataInterpolationType(d_mu_vc_interp_type);
        if (p_vc_point_fac_op) p_vc_point_fac_op->setOperatorScaling(d_A_scale);

        Pointer<VCSCViscousOperator> p_velocity_op = p_velocity_solver->getOperator();
        if (p_velocity_op) p_velocity_op->setDPatchDataInterpolationType(d_mu_vc_interp_type);

        if (p_vc_point_fac_op)
        {
            Pointer<VCSCViscousPETScLevelSolver> p_vc_level_solver = p_vc_point_fac_op->getCoarseSolver();
            if (p_poisson_fac_pc) p_vc_level_solver->setViscosityInterpolationType(d_mu_vc_interp_type);
        }
    }

    // Setup the convective operator.
    d_convective_op = getConvectiveOperator();

    // Setup a boundary op to set velocity boundary conditions on regrid.
    d_fill_after_regrid_phys_bdry_bc_op.reset(new CartSideRobinPhysBdryOp(d_U_scratch_idx, d_U_bc_coefs, false));

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
} // initializeHierarchyIntegrator

void
INSVCStaggeredHierarchyIntegrator::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                            Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    HierarchyIntegrator::initializePatchHierarchy(hierarchy, gridding_alg);

    // When necessary, initialize the value of the advection velocity registered
    // with a coupled advection-diffusion solver.
    if (d_adv_diff_hier_integrator)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int U_adv_diff_current_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getCurrentContext());
        if (isAllocatedPatchData(U_adv_diff_current_idx))
        {
            copySideToFace(U_adv_diff_current_idx, d_U_current_idx, d_hierarchy);
        }
    }
    return;
} // initializePatchHierarhcy

void
INSVCStaggeredHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                                const double new_time,
                                                                const int num_cycles)
{
    INSHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data, new_time);
        level->allocatePatchData(d_velocity_C_idx, current_time);
        level->allocatePatchData(d_velocity_L_idx, current_time);
        level->allocatePatchData(d_velocity_rhs_C_idx, current_time);
        level->allocatePatchData(d_velocity_D_idx, current_time);
        level->allocatePatchData(d_velocity_D_cc_idx, current_time);
        level->allocatePatchData(d_velocity_rhs_D_idx, current_time);
        level->allocatePatchData(d_pressure_D_idx, current_time);
        level->allocatePatchData(d_pressure_rhs_D_idx, current_time);
        level->allocatePatchData(d_temp_sc_idx, current_time);
        level->allocatePatchData(d_mu_interp_idx, current_time);
        level->allocatePatchData(d_N_full_idx, current_time);
        if (d_mu_var.isNull()) level->allocatePatchData(d_mu_scratch_idx, current_time);
        if (!level->checkAllocated(d_mu_linear_op_idx)) level->allocatePatchData(d_mu_linear_op_idx, current_time);
        if (!level->checkAllocated(d_mu_interp_linear_op_idx))
            level->allocatePatchData(d_mu_interp_linear_op_idx, current_time);
        if (!level->checkAllocated(d_rho_linear_op_idx)) level->allocatePatchData(d_rho_linear_op_idx, current_time);
    }

    // Preprocess the operators and solvers
    preprocessOperatorsAndSolvers(current_time, new_time);

    // Preprocess Brinkman penalization objects.
    for (auto& brinkman_force : d_brinkman_force)
    {
        brinkman_force->setTimeInterval(current_time, new_time);
        brinkman_force->preprocessComputeBrinkmanPenalization(current_time, new_time, num_cycles);
    }

    // Keep track of the time-lagged velocity
    d_hier_sc_data_ops->copyData(d_U_old_new_idx, d_U_current_idx);

    return;
} // preprocessIntegrateHierarchy

void
INSVCStaggeredHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
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
        cfl_max = IBTK_MPI::maxReduction(cfl_max);
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

    // Deallocate any temporary data used to compute coefficients
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_velocity_C_idx);
        level->deallocatePatchData(d_velocity_L_idx);
        level->deallocatePatchData(d_velocity_D_idx);
        level->deallocatePatchData(d_velocity_D_cc_idx);
        level->deallocatePatchData(d_pressure_D_idx);
        level->deallocatePatchData(d_velocity_rhs_C_idx);
        level->deallocatePatchData(d_velocity_rhs_D_idx);
        level->deallocatePatchData(d_pressure_rhs_D_idx);
        level->deallocatePatchData(d_temp_sc_idx);
        level->deallocatePatchData(d_mu_interp_idx);
        level->deallocatePatchData(d_N_full_idx);
        if (d_mu_var.isNull()) level->deallocatePatchData(d_mu_scratch_idx);
    }

    // Postprocess Brinkman penalization objects.
    for (auto& brinkman_force : d_brinkman_force)
    {
        brinkman_force->postprocessComputeBrinkmanPenalization(current_time, new_time, num_cycles);
    }

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

void
INSVCStaggeredHierarchyIntegrator::removeNullSpace(const Pointer<SAMRAIVectorReal<NDIM, double> >& sol_vec)
{
    if (d_nul_vecs.empty()) return;
    for (const auto& nul_vec : d_nul_vecs)
    {
        const double sol_dot_nul = sol_vec->dot(nul_vec);
        const double nul_L2_norm = std::sqrt(nul_vec->dot(nul_vec));
        sol_vec->axpy(-sol_dot_nul / nul_L2_norm, nul_vec, sol_vec);
    }
    return;
} // removeNullSpace

void
INSVCStaggeredHierarchyIntegrator::registerMassDensityVariable(Pointer<Variable<NDIM> > rho_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_rho_var);
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_rho_var = rho_var;
    return;
} // registerMassDensityVariable

Pointer<Variable<NDIM> >
INSVCStaggeredHierarchyIntegrator::getMassDensityVariable() const
{
    return d_rho_var;
} // getMassDensityVariable

void
INSVCStaggeredHierarchyIntegrator::registerViscosityVariable(Pointer<Variable<NDIM> > mu_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_mu_var);
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_mu_var = mu_var;
    return;
} // registerViscosityVariable

Pointer<Variable<NDIM> >
INSVCStaggeredHierarchyIntegrator::getViscosityVariable() const
{
    return d_mu_var;
} // getViscosityVariable

void
INSVCStaggeredHierarchyIntegrator::setDensityVCInterpolationType(const IBTK::VCInterpType vc_interp_type)
{
    d_rho_vc_interp_type = vc_interp_type;
    return;
} // setDensityVCInterpolationType

void
INSVCStaggeredHierarchyIntegrator::setViscosityVCInterpolationType(const IBTK::VCInterpType vc_interp_type)
{
    d_mu_vc_interp_type = vc_interp_type;
    return;
} // setViscosityVCInterpolationType

void
INSVCStaggeredHierarchyIntegrator::registerResetFluidDensityFcn(ResetFluidPropertiesFcnPtr callback, void* ctx)
{
    d_reset_rho_fcns.push_back(callback);
    d_reset_rho_fcns_ctx.push_back(ctx);
    return;
} // registerResetFluidDensityFcn

void
INSVCStaggeredHierarchyIntegrator::registerResetFluidViscosityFcn(ResetFluidPropertiesFcnPtr callback, void* ctx)
{
    d_reset_mu_fcns.push_back(callback);
    d_reset_mu_fcns_ctx.push_back(ctx);
    return;
} // registerResetFluidViscosityFcn

void
INSVCStaggeredHierarchyIntegrator::registerBrinkmanPenalizationStrategy(
    Pointer<BrinkmanPenalizationStrategy> brinkman_force)
{
    d_brinkman_force.push_back(brinkman_force);
    return;
} // registerBrinkmanPenalizationStrategy

void
INSVCStaggeredHierarchyIntegrator::registerMassDensityInitialConditions(const Pointer<CartGridFunction> rho_init_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_rho_init_fcn = rho_init_fcn;
    return;
} // registerMassDensityInitialConditions

void
INSVCStaggeredHierarchyIntegrator::registerViscosityInitialConditions(const Pointer<CartGridFunction> mu_init_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_mu_init_fcn = mu_init_fcn;
    return;
} // registerViscosityInitialConditions

void
INSVCStaggeredHierarchyIntegrator::registerViscosityBoundaryConditions(RobinBcCoefStrategy<NDIM>* mu_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_mu_bc_coef = mu_bc_coef;
    return;
} // registerViscosityBoundaryConditions

void
INSVCStaggeredHierarchyIntegrator::setTransportedViscosityVariable(Pointer<CellVariable<NDIM, double> > mu_adv_diff_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_adv_diff_hier_integrator);
#endif
    d_mu_adv_diff_var = mu_adv_diff_var;
    return;
} // setTransportedViscosityVariable

Pointer<CellVariable<NDIM, double> >
INSVCStaggeredHierarchyIntegrator::getTransportedViscosityVariable() const
{
    return d_mu_adv_diff_var;
} // getTransportedViscosityVariable

/////////////////////////////// PROTECTED ////////////////////////////////////

double
INSVCStaggeredHierarchyIntegrator::getStableTimestep(Pointer<Patch<NDIM> > patch) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const hier::Index<NDIM>& ilower = patch->getBox().lower();
    const hier::Index<NDIM>& iupper = patch->getBox().upper();

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

void
INSVCStaggeredHierarchyIntegrator::regridHierarchyBeginSpecialized()
{
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
    d_div_U_norm_1_pre = d_hier_cc_data_ops->L1Norm(d_Div_U_idx, wgt_cc_idx);
    d_div_U_norm_2_pre = d_hier_cc_data_ops->L2Norm(d_Div_U_idx, wgt_cc_idx);
    d_div_U_norm_oo_pre = d_hier_cc_data_ops->maxNorm(d_Div_U_idx, wgt_cc_idx);

    return;
} // regridHierarchyBeginSpecialized

void
INSVCStaggeredHierarchyIntegrator::regridHierarchyEndSpecialized()
{
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
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
    d_div_U_norm_1_post = d_hier_cc_data_ops->L1Norm(d_Div_U_idx, wgt_cc_idx);
    d_div_U_norm_2_post = d_hier_cc_data_ops->L2Norm(d_Div_U_idx, wgt_cc_idx);
    d_div_U_norm_oo_post = d_hier_cc_data_ops->maxNorm(d_Div_U_idx, wgt_cc_idx);
    d_do_regrid_projection = d_div_U_norm_1_post > d_regrid_max_div_growth_factor * d_div_U_norm_1_pre ||
                             d_div_U_norm_2_post > d_regrid_max_div_growth_factor * d_div_U_norm_2_pre ||
                             d_div_U_norm_oo_post > d_regrid_max_div_growth_factor * d_div_U_norm_oo_pre;
    return;
} // regridHierarchyEndSpecialized

void
INSVCStaggeredHierarchyIntegrator::initializeCompositeHierarchyDataSpecialized(const double /*init_data_time*/,
                                                                               const bool initial_time)
{
    // Project the interpolated velocity if needed.
    if (initial_time || d_do_regrid_projection)
    {
        plog << d_object_name << "::initializeCompositeHierarchyData():\n"
             << "  projecting the interpolated velocity field\n";
        regridProjection();
        d_do_regrid_projection = false;
    }
    return;
} // initializeCompositeHierarchyDataSpecialized

void
INSVCStaggeredHierarchyIntegrator::initializeLevelDataSpecialized(
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
            copy_data.registerRefine(d_U_regrid_idx, d_U_regrid_idx, d_U_regrid_idx, nullptr);
            copy_data.registerRefine(d_U_src_idx, d_U_src_idx, d_U_src_idx, nullptr);
            copy_data.registerRefine(d_indicator_idx, d_indicator_idx, d_indicator_idx, nullptr);
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
        fill_div_free_prolongation.registerRefine(d_U_current_idx, d_U_current_idx, d_U_regrid_idx, nullptr);
        Pointer<RefineOperator<NDIM> > refine_op = grid_geom->lookupRefineOperator(d_U_var, d_U_refine_type);
        Pointer<CoarsenOperator<NDIM> > coarsen_op = grid_geom->lookupCoarsenOperator(d_U_var, d_U_coarsen_type);
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
            using InterpolationTransactionComponent =
                HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
            InterpolationTransactionComponent U_bc_component(d_U_scratch_idx,
                                                             DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             d_bdry_extrap_type, // TODO: update variable name
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
INSVCStaggeredHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
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

    d_hier_nc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_nc_data_ops->resetLevels(0, finest_hier_level);

    d_hier_ec_data_ops->setPatchHierarchy(hierarchy);
    d_hier_ec_data_ops->resetLevels(0, finest_hier_level);

    // Setup the patch boundary filling objects.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent U_bc_component(d_U_scratch_idx,
                                                     DATA_REFINE_TYPE,
                                                     USE_CF_INTERPOLATION,
                                                     DATA_COARSEN_TYPE,
                                                     d_bdry_extrap_type, // TODO: update variable name
                                                     CONSISTENT_TYPE_2_BDRY,
                                                     d_U_bc_coefs);
    d_U_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_U_bdry_bc_fill_op->initializeOperatorState(U_bc_component, d_hierarchy);

    InterpolationTransactionComponent P_bc_component(d_P_scratch_idx,
                                                     DATA_REFINE_TYPE,
                                                     USE_CF_INTERPOLATION,
                                                     DATA_COARSEN_TYPE,
                                                     d_bdry_extrap_type, // TODO: update variable name
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
                                                         d_bdry_extrap_type, // TODO: update variable name
                                                         CONSISTENT_TYPE_2_BDRY);
        d_Q_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
        d_Q_bdry_bc_fill_op->initializeOperatorState(Q_bc_component, d_hierarchy);
    }

    if (!d_mu_is_const)
    {
        // These options are chosen to ensure that information is propagated
        // conservatively from the coarse cells only
        InterpolationTransactionComponent mu_bc_component(
            d_mu_scratch_idx, d_mu_refine_type, false, d_mu_coarsen_type, d_mu_bdry_extrap_type, false, d_mu_bc_coef);
        d_mu_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
        d_mu_bdry_bc_fill_op->initializeOperatorState(mu_bc_component, d_hierarchy);
    }

    // Setup the patch boundary synchronization objects.
    using SynchronizationTransactionComponent = SideDataSynchronization::SynchronizationTransactionComponent;
    SynchronizationTransactionComponent synch_transaction =
        SynchronizationTransactionComponent(d_U_scratch_idx, d_U_coarsen_type);
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
INSVCStaggeredHierarchyIntegrator::applyGradientDetectorSpecialized(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
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
            thresh += std::sqrt(std::numeric_limits<double>::epsilon());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, int> > tags_data = patch->getPatchData(tag_index);
                Pointer<CellData<NDIM, double> > Omega_data = patch->getPatchData(d_Omega_idx);
                for (CellIterator<NDIM> ic(patch_box); ic; ic++)
                {
                    const hier::Index<NDIM>& i = ic();
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
                    const double norm_Omega = std::sqrt(norm_Omega_sq);
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
INSVCStaggeredHierarchyIntegrator::setupPlotDataSpecialized()
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

    // Compute EE = 0.5*(grad u + grad u^T).
    if (d_output_EE)
    {
        const int EE_idx = var_db->mapVariableAndContextToIndex(d_EE_var, ctx);
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(d_U_scratch_idx, d_integrator_time);
        }
        d_hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
        d_U_bdry_bc_fill_op->fillData(d_integrator_time);
        d_hier_math_ops->strain_rate(EE_idx, d_EE_var, d_U_scratch_idx, d_U_var, d_no_fill_op, d_integrator_time);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(d_U_scratch_idx);
        }
    }
    return;
} // setupPlotDataSpecialized

void
INSVCStaggeredHierarchyIntegrator::copySideToFace(const int U_fc_idx,
                                                  const int U_sc_idx,
                                                  Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    copy_side_to_face(U_fc_idx, U_sc_idx, hierarchy);
    return;
} // copySideToFace

/////////////////////////////// PRIVATE //////////////////////////////////////

void
INSVCStaggeredHierarchyIntegrator::preprocessOperatorsAndSolvers(const double current_time, const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double rho = d_rho_is_const ? d_problem_coefs.getRho() : -1.0;
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const int wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();

    // Setup solver vectors.
    const bool has_velocity_nullspace =
        d_normalize_velocity && (d_rho_is_const && MathUtilities<double>::equalEps(rho, 0.0));
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

        for (const auto& nul_vec : d_nul_vecs)
        {
            if (nul_vec) nul_vec->freeVectorComponents();
        }
        const int n_nul_vecs = (has_pressure_nullspace ? 1 : 0) + (has_velocity_nullspace ? NDIM : 0);
        d_nul_vecs.resize(n_nul_vecs);

        for (const auto& U_nul_vec : d_U_nul_vecs)
        {
            if (U_nul_vec) U_nul_vec->freeVectorComponents();
        }
        const int n_U_nul_vecs = (has_velocity_nullspace ? NDIM : 0);
        d_U_nul_vecs.resize(n_U_nul_vecs);

        if (has_velocity_nullspace)
        {
            for (unsigned int k = 0; k < NDIM; ++k)
            {
                d_nul_vecs[k] = d_sol_vec->cloneVector(d_object_name + "::nul_vec_U_" + std::to_string(k));
                d_nul_vecs[k]->allocateVectorData(current_time);
                d_nul_vecs[k]->setToScalar(0.0);
                d_U_nul_vecs[k] = d_U_scratch_vec->cloneVector(d_object_name + "::U_nul_vec_U_" + std::to_string(k));
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
        auto U_bc_coef = dynamic_cast<INSVCStaggeredVelocityBcCoef*>(d_U_bc_coefs[d]);
        U_bc_coef->setStokesSpecifications(&d_problem_coefs);
        U_bc_coef->setPhysicalBcCoefs(d_bc_coefs);
        U_bc_coef->setSolutionTime(new_time);
        U_bc_coef->setTimeInterval(current_time, new_time);
    }
    auto P_bc_coef = dynamic_cast<INSVCStaggeredPressureBcCoef*>(d_P_bc_coef);
    P_bc_coef->setStokesSpecifications(&d_problem_coefs);
    P_bc_coef->setPhysicalBcCoefs(d_bc_coefs);
    P_bc_coef->setSolutionTime(new_time);
    P_bc_coef->setTimeInterval(current_time, new_time);
    P_bc_coef->setViscosityInterpolationType(d_mu_vc_interp_type);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto U_star_bc_coef = dynamic_cast<INSIntermediateVelocityBcCoef*>(d_U_star_bc_coefs[d]);
        U_star_bc_coef->setPhysicalBcCoefs(d_bc_coefs);
        U_star_bc_coef->setSolutionTime(new_time);
        U_star_bc_coef->setTimeInterval(current_time, new_time);
    }
    auto Phi_bc_coef = dynamic_cast<INSProjectionBcCoef*>(d_Phi_bc_coef.get());
    Phi_bc_coef->setPhysicalBcCoefs(d_bc_coefs);
    Phi_bc_coef->setSolutionTime(0.5 * (current_time + new_time));
    Phi_bc_coef->setTimeInterval(current_time, new_time);

    // Setup convective operator.
    if (d_convective_op && d_convective_op_needs_init)
    {
        if (d_enable_logging)
            plog << d_object_name
                 << "::preprocessIntegrateHierarchy(): initializing "
                    "convective operator"
                 << std::endl;
        d_convective_op->setAdvectionVelocity(d_U_scratch_idx);
        d_convective_op->setSolutionTime(d_integrator_time);
        d_convective_op->setTimeInterval(current_time, new_time);
        d_convective_op->initializeOperatorState(*d_U_scratch_vec, *d_U_rhs_vec);
        d_convective_op_needs_init = false;
    }
    return;
} // preprocessOperatorsAndSolvers
//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
