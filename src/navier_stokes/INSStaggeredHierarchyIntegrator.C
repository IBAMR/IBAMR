// Filename: INSStaggeredHierarchyIntegrator.C
// Created on 20 Mar 2008 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
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

#include "INSStaggeredHierarchyIntegrator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/INSIntermediateVelocityBcCoef.h>
#include <ibamr/INSProjectionBcCoef.h>
#include <ibamr/INSStaggeredConvectiveOperatorManager.h>
#include <ibamr/INSStaggeredPressureBcCoef.h>
#include <ibamr/INSStaggeredVelocityBcCoef.h>
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/CCPoissonSolverManager.h>
#include <ibtk/CartSideDoubleDivPreservingRefine.h>
#include <ibtk/CartSideDoubleSpecializedConstantRefine.h>
#include <ibtk/CartSideDoubleSpecializedLinearRefine.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/NewtonKrylovSolver.h>
#include <ibtk/RefinePatchStrategySet.h>
#include <ibtk/SCPoissonSolverManager.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define NAVIER_STOKES_SC_STABLEDT_FC FC_FUNC_(navier_stokes_sc_stabledt2d, NAVIER_STOKES_SC_STABLEDT2D)
#define NAVIER_STOKES_SIDE_TO_FACE_FC FC_FUNC_(navier_stokes_side_to_face2d, NAVIER_STOKES_SIDE_TO_FACE2D)
#define NAVIER_STOKES_STAGGERED_ADV_SOURCE_FC FC_FUNC_(navier_stokes_staggered_adv_source2d, NAVIER_STOKES_STAGGERED_ADV_SOURCE2D)
#define NAVIER_STOKES_STAGGERED_CONS_SOURCE_FC FC_FUNC_(navier_stokes_staggered_cons_source2d, NAVIER_STOKES_STAGGERED_CONS_SOURCE2D)
#define NAVIER_STOKES_STAGGERED_SKEW_SYM_SOURCE_FC FC_FUNC_(navier_stokes_staggered_skew_sym_source2d, NAVIER_STOKES_STAGGERED_SKEW_SYM_SOURCE2D)
#endif

#if (NDIM == 3)
#define NAVIER_STOKES_SC_STABLEDT_FC FC_FUNC_(navier_stokes_sc_stabledt3d, NAVIER_STOKES_SC_STABLEDT3D)
#define NAVIER_STOKES_SIDE_TO_FACE_FC FC_FUNC_(navier_stokes_side_to_face3d, NAVIER_STOKES_SIDE_TO_FACE3D)
#define NAVIER_STOKES_STAGGERED_ADV_SOURCE_FC FC_FUNC_(navier_stokes_staggered_adv_source3d, NAVIER_STOKES_STAGGERED_ADV_SOURCE3D)
#define NAVIER_STOKES_STAGGERED_CONS_SOURCE_FC FC_FUNC_(navier_stokes_staggered_cons_source3d, NAVIER_STOKES_STAGGERED_CONS_SOURCE3D)
#define NAVIER_STOKES_STAGGERED_SKEW_SYM_SOURCE_FC FC_FUNC_(navier_stokes_staggered_skew_sym_source3d, NAVIER_STOKES_STAGGERED_SKEW_SYM_SOURCE3D)
#endif

extern "C"
{
    void
    NAVIER_STOKES_SC_STABLEDT_FC(
        const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
#endif
        double&
                                 );

    void
    NAVIER_STOKES_SIDE_TO_FACE_FC(
#if (NDIM == 2)
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* , const int& ,
        double* , double* , const int&
#endif
#if (NDIM == 3)
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* , const int& ,
        double* , double* , double* , const int&
#endif
                                  );

    void
    NAVIER_STOKES_STAGGERED_ADV_SOURCE_FC(
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const double* ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* ,
        double* , double* , double*
#endif
                                          );

    void
    NAVIER_STOKES_STAGGERED_CONS_SOURCE_FC(
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const double* ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* ,
        double* , double* , double*
#endif
                                           );

    void
    NAVIER_STOKES_STAGGERED_SKEW_SYM_SOURCE_FC(
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const double* ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* ,
        double* , double* , double*
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
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
static const int SIDEG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string CELL_DATA_COARSEN_TYPE = "CUBIC_COARSEN";
static const std::string SIDE_DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Copy data from a side-centered variable to a face-centered variable.
void
copy_side_to_face(
    const int U_fc_idx,
    const int U_sc_idx,
    Pointer<PatchHierarchy<NDIM> > hierarchy)
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
            Pointer<SideData<NDIM,double> > U_sc_data = patch->getPatchData(U_sc_idx);
            Pointer<FaceData<NDIM,double> > U_fc_data = patch->getPatchData(U_fc_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(U_sc_data->getGhostCellWidth().min() == U_sc_data->getGhostCellWidth().max());
                TBOX_ASSERT(U_fc_data->getGhostCellWidth().min() == U_fc_data->getGhostCellWidth().max());
#endif
                const int U_sc_gcw = U_sc_data->getGhostCellWidth().max();
                const int U_fc_gcw = U_fc_data->getGhostCellWidth().max();
                NAVIER_STOKES_SIDE_TO_FACE_FC(
                    ilower(0), iupper(0),
                    ilower(1), iupper(1),
#if (NDIM == 3)
                    ilower(2), iupper(2),
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
}// copy_side_to_face
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredHierarchyIntegrator::INSStaggeredHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    bool register_for_restart)
    : INSHierarchyIntegrator(
        object_name, input_db,
        new SideVariable<NDIM,double>(object_name+"::U"),
        new CellVariable<NDIM,double>(object_name+"::P"),
        new SideVariable<NDIM,double>(object_name+"::F"),
        new CellVariable<NDIM,double>(object_name+"::Q"),
        register_for_restart)
{
    // Check to make sure the time stepping types are supported.
    switch (d_viscous_time_stepping_type)
    {
        case BACKWARD_EULER:
        case FORWARD_EULER:
        case TRAPEZOIDAL_RULE:
            break;
        default:
            TBOX_ERROR(d_object_name << "::INSStaggeredHierarchyIntegrator():\n"
                       << "  unsupported viscous time stepping type: " << enum_to_string<TimeSteppingType>(d_viscous_time_stepping_type) << " \n"
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
                       << "  unsupported convective time stepping type: " << enum_to_string<TimeSteppingType>(d_convective_time_stepping_type) << " \n"
                       << "  valid choices are: ADAMS_BASHFORTH, FORWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
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
                           << "  unsupported initial convective time stepping type: " << enum_to_string<TimeSteppingType>(d_init_convective_time_stepping_type) << " \n"
                           << "  valid choices are: FORWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
        }
    }

    // Check to see whether the convective operator type has been set.
    d_default_convective_op_type = INSStaggeredConvectiveOperatorManager::DEFAULT;
    if      (input_db->keyExists("convective_op_type"))               d_default_convective_op_type = input_db->getString("convective_op_type");
    else if (input_db->keyExists("convective_operator_type"))         d_default_convective_op_type = input_db->getString("convective_operator_type");
    else if (input_db->keyExists("default_convective_op_type"))       d_default_convective_op_type = input_db->getString("default_convective_op_type");
    else if (input_db->keyExists("default_convective_operator_type")) d_default_convective_op_type = input_db->getString("default_convective_operator_type");

    // Setup physical boundary conditions objects.
    d_bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
    d_P_bc_coef = new INSStaggeredPressureBcCoef(&d_problem_coefs,d_bc_coefs);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_U_bc_coefs[d] = new INSStaggeredVelocityBcCoef(d,&d_problem_coefs,d_bc_coefs,dynamic_cast<INSStaggeredPressureBcCoef*>(d_P_bc_coef));
    }

    // Initialize all variables.  The velocity, pressure, body force, and fluid
    // source variables were created above in the constructor for the
    // INSHierarchyIntegrator base class.
    d_U_var          = INSHierarchyIntegrator::d_U_var;
    d_P_var          = INSHierarchyIntegrator::d_P_var;
    d_F_var          = INSHierarchyIntegrator::d_F_var;
    d_Q_var          = INSHierarchyIntegrator::d_Q_var;
    d_N_old_var      = new SideVariable<NDIM,double>(d_object_name+"::N_old"      );

    d_U_cc_var       = new CellVariable<NDIM,double>(d_object_name+"::U_cc",  NDIM);
    d_F_cc_var       = new CellVariable<NDIM,double>(d_object_name+"::F_cc",  NDIM);
#if (NDIM == 2)
    d_Omega_var      = new CellVariable<NDIM,double>(d_object_name+"::Omega"      );
#endif
#if (NDIM == 3)
    d_Omega_var      = new CellVariable<NDIM,double>(d_object_name+"::Omega", NDIM);
#endif
    d_Div_U_var      = new CellVariable<NDIM,double>(d_object_name+"::Div_U"      );

#if (NDIM == 3)
    d_Omega_Norm_var = new CellVariable<NDIM,double>(d_object_name+"::|Omega|_2"  );
#endif
    d_U_regrid_var   = new SideVariable<NDIM,double>(d_object_name+"::U_regrid"   );
    d_U_src_var      = new SideVariable<NDIM,double>(d_object_name+"::U_src"      );
    d_indicator_var  = new SideVariable<NDIM,double>(d_object_name+"::indicator"  );
    d_F_div_var      = new SideVariable<NDIM,double>(d_object_name+"::F_div"      );
    return;
}// INSStaggeredHierarchyIntegrator

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
    if (!d_U_rhs_vec.isNull()) d_U_rhs_vec->freeVectorComponents();
    if (!d_U_adv_vec.isNull()) d_U_adv_vec->freeVectorComponents();
    if (!d_N_vec    .isNull()) d_N_vec    ->freeVectorComponents();
    if (!d_P_rhs_vec.isNull()) d_P_rhs_vec->freeVectorComponents();
    for (unsigned int k = 0; k < d_nul_vecs.size(); ++k)
    {
        if (!d_nul_vecs[k].isNull()) d_nul_vecs[k]->freeVectorComponents();
    }
    for (unsigned int k = 0; k < d_U_nul_vecs.size(); ++k)
    {
        if (!d_U_nul_vecs[k].isNull()) d_U_nul_vecs[k]->freeVectorComponents();
    }
    return;
}// ~INSStaggeredHierarchyIntegrator

const std::vector<RobinBcCoefStrategy<NDIM>*>&
INSStaggeredHierarchyIntegrator::getVelocityBoundaryConditions() const
{
    return d_U_bc_coefs;
}// getVelocityBoundaryConditions

RobinBcCoefStrategy<NDIM>*
INSStaggeredHierarchyIntegrator::getPressureBoundaryConditions() const
{
    return d_P_bc_coef;
}// getPressureBoundaryConditions

Pointer<ConvectiveOperator>
INSStaggeredHierarchyIntegrator::getConvectiveOperator()
{
    if (d_creeping_flow)
    {
        d_convective_op.setNull();
    }
    else if (d_convective_op.isNull())
    {
        INSStaggeredConvectiveOperatorManager* convective_op_manager = INSStaggeredConvectiveOperatorManager::getManager();
        d_convective_op = convective_op_manager->allocateOperator(
            d_default_convective_op_type, d_object_name+"::ConvectiveOperator", d_default_convective_difference_form, d_default_convective_bdry_extrap_type);
        d_convective_op_needs_init = true;
    }
    return d_convective_op;
}// getConvectiveOperator

Pointer<PoissonSolver>
INSStaggeredHierarchyIntegrator::getVelocitySubdomainSolver()
{
    if (d_velocity_solver.isNull())
    {
        d_velocity_solver = SCPoissonSolverManager::getManager()->allocateSolver(
            d_velocity_solver_type , d_object_name+"::velocity_solver" , d_velocity_solver_db ,
            d_velocity_precond_type, d_object_name+"::velocity_precond", d_velocity_precond_db);
    }
    return d_velocity_solver;
}// getVelocitySubdomainSolver

Pointer<PoissonSolver>
INSStaggeredHierarchyIntegrator::getPressureSubdomainSolver()
{
    if (d_pressure_solver.isNull())
    {
        d_pressure_solver = CCPoissonSolverManager::getManager()->allocateSolver(
            d_pressure_solver_type , d_object_name+"::pressure_solver" , d_pressure_solver_db ,
            d_pressure_precond_type, d_object_name+"::pressure_precond", d_pressure_precond_db);
    }
    return d_pressure_solver;
}// getPressureSubdomainSolver

void
INSStaggeredHierarchyIntegrator::setStokesSolver(
    Pointer<StaggeredStokesSolver> stokes_solver)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_stokes_solver.isNull());
#endif
    d_stokes_solver = stokes_solver;
    d_stokes_solver_needs_init = true;
    return;
}// setStokesSolver

Pointer<StaggeredStokesSolver>
INSStaggeredHierarchyIntegrator::getStokesSolver()
{
    if (d_stokes_solver.isNull())
    {
        TBOX_ERROR("ack!\n");
    }
    return d_stokes_solver;
}// getStokesSolver

void
INSStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;

    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(new CellVariable<NDIM,double>("cc_var"), hierarchy, true);
    d_hier_fc_data_ops = hier_ops_manager->getOperationsDouble(new FaceVariable<NDIM,double>("fc_var"), hierarchy, true);
    d_hier_sc_data_ops = hier_ops_manager->getOperationsDouble(new SideVariable<NDIM,double>("sc_var"), hierarchy, true);
    d_hier_math_ops = buildHierarchyMathOps(d_hierarchy);

    // Register state variables that are maintained by the
    // INSStaggeredHierarchyIntegrator.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    grid_geom->addSpatialRefineOperator(new CartSideDoubleSpecializedConstantRefine());
    grid_geom->addSpatialRefineOperator(new CartSideDoubleSpecializedLinearRefine());

    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> side_ghosts = SIDEG;
    const IntVector<NDIM> no_ghosts = 0;

    registerVariable(d_U_current_idx, d_U_new_idx, d_U_scratch_idx, d_U_var, side_ghosts, "CONSERVATIVE_COARSEN", "SPECIALIZED_LINEAR_REFINE", d_U_init);

    registerVariable(d_P_current_idx, d_P_new_idx, d_P_scratch_idx, d_P_var, cell_ghosts, "CONSERVATIVE_COARSEN", "LINEAR_REFINE", d_P_init);

    if (!d_F_fcn.isNull())
    {
        registerVariable(d_F_current_idx, d_F_new_idx, d_F_scratch_idx, d_F_var, side_ghosts, "CONSERVATIVE_COARSEN", "SPECIALIZED_LINEAR_REFINE", d_F_fcn);
    }
    else
    {
        d_F_current_idx = -1;
        d_F_new_idx     = -1;
        d_F_scratch_idx = -1;
    }

    if (!d_Q_fcn.isNull())
    {
        registerVariable(d_Q_current_idx, d_Q_new_idx, d_Q_scratch_idx, d_Q_var, cell_ghosts, "CONSERVATIVE_COARSEN", "CONSTANT_REFINE", d_Q_fcn);
    }
    else
    {
        d_Q_current_idx = -1;
        d_Q_new_idx     = -1;
        d_Q_scratch_idx = -1;
    }

    registerVariable(d_N_old_current_idx, d_N_old_new_idx, d_N_old_scratch_idx, d_N_old_var, side_ghosts, "CONSERVATIVE_COARSEN", "SPECIALIZED_LINEAR_REFINE");

    // Register plot variables that are maintained by the
    // INSCollocatedHierarchyIntegrator.
    registerVariable(d_U_cc_idx, d_U_cc_var, no_ghosts, getCurrentContext());
    if (!d_F_fcn.isNull())
    {
        registerVariable(d_F_cc_idx, d_F_cc_var, no_ghosts, getCurrentContext());
    }
    else
    {
        d_F_cc_idx = -1;
    }
    registerVariable(d_Omega_idx, d_Omega_var,   no_ghosts, getCurrentContext());
    registerVariable(d_Div_U_idx, d_Div_U_var, cell_ghosts, getCurrentContext());

    // Register scratch variables that are maintained by the
    // INSStaggeredHierarchyIntegrator.
#if (NDIM == 3)
    registerVariable(d_Omega_Norm_idx, d_Omega_Norm_var, no_ghosts);
#endif
    registerVariable( d_U_regrid_idx,  d_U_regrid_var, CartSideDoubleDivPreservingRefine::REFINE_OP_STENCIL_WIDTH);
    registerVariable(    d_U_src_idx,     d_U_src_var, CartSideDoubleDivPreservingRefine::REFINE_OP_STENCIL_WIDTH);
    registerVariable(d_indicator_idx, d_indicator_var, CartSideDoubleDivPreservingRefine::REFINE_OP_STENCIL_WIDTH);
    if (!d_Q_fcn.isNull())
    {
        registerVariable(d_F_div_idx, d_F_div_var, no_ghosts);
    }
    else
    {
        d_F_div_idx = -1;
    }

    // Register variables for plotting.
    if (!d_visit_writer.isNull())
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

        if (!d_F_fcn.isNull() && d_output_F)
        {
            d_visit_writer->registerPlotQuantity("F", "VECTOR", d_F_cc_idx, 0, d_F_scale);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == 0) d_visit_writer->registerPlotQuantity("F_x", "SCALAR", d_F_cc_idx, d, d_F_scale);
                if (d == 1) d_visit_writer->registerPlotQuantity("F_y", "SCALAR", d_F_cc_idx, d, d_F_scale);
                if (d == 2) d_visit_writer->registerPlotQuantity("F_z", "SCALAR", d_F_cc_idx, d, d_F_scale);
            }
        }

        if (!d_Q_fcn.isNull() && d_output_Q)
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
    registerCoarsenAlgorithm(d_object_name+"::CONVECTIVE_OP", coarsen_alg);

    // Setup the Stokes solver.
    d_stokes_solver = getStokesSolver();

    // Setup the convective operator.
    d_convective_op = getConvectiveOperator();

    // Setup a boundary op to set velocity boundary conditions on regrid.
    d_fill_after_regrid_phys_bdry_bc_op = new CartSideRobinPhysBdryOp(d_U_scratch_idx, d_U_bc_coefs, false);

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
}// initializeHierarchyIntegrator

void
INSStaggeredHierarchyIntegrator::initializePatchHierarchy(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
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
    if (!d_adv_diff_hier_integrator.isNull())
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int U_adv_diff_current_idx = var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getCurrentContext());
        copy_side_to_face(U_adv_diff_current_idx, d_U_current_idx, d_hierarchy);
    }
    return;
}// initializePatchHierarhcy

void
INSStaggeredHierarchyIntegrator::preprocessIntegrateHierarchy(
    const double current_time,
    const double new_time,
    const int num_cycles)
{
    INSHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;

    // Keep track of the number of cycles to be used for the present integration
    // step.
    d_current_num_cycles = num_cycles;
    if (!d_creeping_flow && (d_current_num_cycles == 1) && (d_convective_time_stepping_type == MIDPOINT_RULE || d_convective_time_stepping_type == TRAPEZOIDAL_RULE))
    {
        TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                   << "  time stepping type: " << enum_to_string<TimeSteppingType>(d_convective_time_stepping_type) << " requires num_cycles > 1.\n"
                   << "  at current time step, num_cycles = " << d_current_num_cycles << "\n");
    }

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data    ,     new_time);
    }

    // Setup the operators and solvers.
    reinitializeOperatorsAndSolvers(current_time, new_time);

    // Allocate solver vectors.
    d_U_rhs_vec->allocateVectorData(current_time);  d_U_rhs_vec->setToScalar(0.0);
    d_U_adv_vec->allocateVectorData(current_time);  d_U_adv_vec->setToScalar(0.0);
    d_N_vec    ->allocateVectorData(current_time);  d_N_vec    ->setToScalar(0.0);
    d_P_rhs_vec->allocateVectorData(current_time);  d_P_rhs_vec->setToScalar(0.0);

    // Initialize the right-hand side terms.
    const double rho    = d_problem_coefs.getRho();
    const double mu     = d_problem_coefs.getMu();
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
    PoissonSpecifications rhs_spec(d_object_name+"::rhs_spec");
    rhs_spec.setCConstant((rho/dt)-K_rhs*lambda);
    rhs_spec.setDConstant(        +K_rhs*mu    );
    const int U_rhs_idx = d_U_rhs_vec->getComponentDescriptorIndex(0);
    const Pointer<SideVariable<NDIM,double> > U_rhs_var = d_U_rhs_vec->getComponentVariable(0);
    d_hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
    d_hier_math_ops->laplace(U_rhs_idx, U_rhs_var, rhs_spec, d_U_scratch_idx, d_U_var, d_U_bdry_bc_fill_op, current_time);

    // Set the initial guess.
    d_hier_sc_data_ops->copyData(d_U_new_idx, d_U_current_idx);
    d_hier_cc_data_ops->copyData(d_P_new_idx, d_P_current_idx);

    // Setup inhomogeneous boundary conditions.
    d_bc_helper->clearBcCoefData();
    d_bc_helper->cacheBcCoefData(d_U_scratch_idx, d_U_var, d_U_bc_coefs, new_time, IntVector<NDIM>(SIDEG), d_hierarchy);
    d_stokes_solver->setHomogeneousBc(false);
    Pointer<KrylovLinearSolver> p_stokes_solver = d_stokes_solver;
    if (!p_stokes_solver.isNull())
    {
        p_stokes_solver->getOperator()->modifyRhsForInhomogeneousBc(*d_U_rhs_vec);
        p_stokes_solver->setHomogeneousBc(true);
    }

    // Initialize any registered advection-diffusion solver.
    if (!d_adv_diff_hier_integrator.isNull())
    {
        const int adv_diff_num_cycles = d_adv_diff_hier_integrator->getNumberOfCycles();
        if (adv_diff_num_cycles != d_current_num_cycles && d_current_num_cycles != 1)
        {
            TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                       << "  attempting to perform " << d_current_num_cycles << " cycles of fixed point iteration.\n"
                       << "  number of cycles required by coupled advection-diffusion solver = " << adv_diff_num_cycles << ".\n"
                       << "  current implementation requires either that both solvers use the same number of cycles,\n"
                       << "  or that the Navier-Stokes solver use only a single cycle.\n");
        }
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int U_adv_diff_current_idx = var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getCurrentContext());
        copy_side_to_face(U_adv_diff_current_idx, d_U_current_idx, d_hierarchy);
        d_adv_diff_hier_integrator->preprocessIntegrateHierarchy(current_time, new_time, adv_diff_num_cycles);
        const int U_adv_diff_scratch_idx = var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getScratchContext());
        const int U_adv_diff_new_idx     = var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getNewContext()    );
        d_hier_fc_data_ops->copyData(U_adv_diff_scratch_idx, U_adv_diff_current_idx);
        d_hier_fc_data_ops->copyData(U_adv_diff_new_idx    , U_adv_diff_current_idx);
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
            Pointer<CoarsenOperator<NDIM> > coarsen_op = grid_geom->lookupCoarsenOperator(d_U_var, "CONSERVATIVE_COARSEN");
            coarsen_alg->registerCoarsen(U_adv_idx, U_adv_idx, coarsen_op);
            coarsen_alg->resetSchedule(getCoarsenSchedules(d_object_name+"::CONVECTIVE_OP")[ln]);
            getCoarsenSchedules(d_object_name+"::CONVECTIVE_OP")[ln]->coarsenData();
            getCoarsenAlgorithm(d_object_name+"::CONVECTIVE_OP")->resetSchedule(getCoarsenSchedules(d_object_name+"::CONVECTIVE_OP")[ln]);
        }
        d_convective_op->setAdvectionVelocity(d_U_adv_vec->getComponentDescriptorIndex(0));
        d_convective_op->apply(*d_U_adv_vec, *d_N_vec);
        const int N_idx = d_N_vec->getComponentDescriptorIndex(0);
        d_hier_sc_data_ops->copyData(d_N_old_new_idx, N_idx);
        if (convective_time_stepping_type == FORWARD_EULER)
        {
            d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0), -1.0*rho, N_idx, d_rhs_vec->getComponentDescriptorIndex(0));
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0), -0.5*rho, N_idx, d_rhs_vec->getComponentDescriptorIndex(0));
        }
    }
    return;
}// preprocessIntegrateHierarchy

void
INSStaggeredHierarchyIntegrator::integrateHierarchy(
    const double current_time,
    const double new_time,
    const int cycle_num)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_current_dt = new_time-current_time);
    TBOX_ASSERT(cycle_num < d_current_num_cycles);
#endif
    d_current_cycle_num = cycle_num;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt  = new_time-current_time;
    const double rho = d_problem_coefs.getRho();

    // Check to make sure that the number of cycles is what we expect it to be.
    const int expected_num_cycles = getNumberOfCycles();
    if (d_current_num_cycles != expected_num_cycles)
    {
        IBAMR_DO_ONCE(
            {
                pout << "INSStaggeredHierarchyIntegrator::integrateHierarchy():\n"
                     << "  WARNING: num_cycles = " << d_current_num_cycles << " but expected num_cycles = " << expected_num_cycles << ".\n";
            }
                      );
    }

    // Perform a single step of fixed point iteration.

    if (!d_adv_diff_hier_integrator.isNull())
    {
        // Update the state variables maintained by the advection-diffusion
        // solver.
        d_adv_diff_hier_integrator->integrateHierarchy(current_time, new_time, cycle_num);
    }

    // Account for the convective acceleration term.
    TimeSteppingType convective_time_stepping_type = d_convective_time_stepping_type;
    if (is_multistep_time_stepping_type(convective_time_stepping_type))
    {
#ifdef DEBUG_CHECK_ASSERTIONS
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
                         << "  WARNING: convective_time_stepping_type = " << enum_to_string<TimeSteppingType>(d_convective_time_stepping_type) << " but num_cycles = " << d_current_num_cycles << " > 1.\n"
                         << "           using " << enum_to_string<TimeSteppingType>(d_convective_time_stepping_type) << " only for the first cycle in each time step;\n"
                         << "           using " << enum_to_string<TimeSteppingType>(  convective_time_stepping_type) << " for subsequent cycles.\n";
                }
                          );
        }
    }
    if (!d_creeping_flow && convective_time_stepping_type != FORWARD_EULER)
    {
        if (cycle_num > 0)
        {
            const int U_adv_idx = d_U_adv_vec->getComponentDescriptorIndex(0);
            if (convective_time_stepping_type == MIDPOINT_RULE)
            {
                d_hier_sc_data_ops->linearSum(U_adv_idx, 0.5, d_U_current_idx, 0.5, d_U_new_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                d_hier_sc_data_ops->copyData(U_adv_idx, d_U_new_idx);
            }
            for (int ln = finest_ln; ln > coarsest_ln; --ln)
            {
                Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
                Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
                Pointer<CoarsenOperator<NDIM> > coarsen_op = grid_geom->lookupCoarsenOperator(d_U_var, "CONSERVATIVE_COARSEN");
                coarsen_alg->registerCoarsen(U_adv_idx, U_adv_idx, coarsen_op);
                coarsen_alg->resetSchedule(getCoarsenSchedules(d_object_name+"::CONVECTIVE_OP")[ln]);
                getCoarsenSchedules(d_object_name+"::CONVECTIVE_OP")[ln]->coarsenData();
                getCoarsenAlgorithm(d_object_name+"::CONVECTIVE_OP")->resetSchedule(getCoarsenSchedules(d_object_name+"::CONVECTIVE_OP")[ln]);
            }
            d_convective_op->setAdvectionVelocity(d_U_adv_vec->getComponentDescriptorIndex(0));
            d_convective_op->apply(*d_U_adv_vec, *d_N_vec);
        }
        const int N_idx = d_N_vec->getComponentDescriptorIndex(0);
        if (convective_time_stepping_type == ADAMS_BASHFORTH)
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(cycle_num == 0);
#endif
            const double omega = dt / d_dt_previous[0];
            d_hier_sc_data_ops->linearSum(N_idx, 1.0 + 0.5*omega, N_idx, -0.5*omega, d_N_old_current_idx);
        }
        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
        {
            d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0), -1.0*rho, N_idx, d_rhs_vec->getComponentDescriptorIndex(0));
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0), -0.5*rho, N_idx, d_rhs_vec->getComponentDescriptorIndex(0));
        }
    }

    // Account for body forcing terms.
    if (!d_F_fcn.isNull())
    {
        d_F_fcn->setDataOnPatchHierarchy(d_F_scratch_idx, d_F_var, d_hierarchy, current_time+0.5*dt);
        d_hier_sc_data_ops->add(d_rhs_vec->getComponentDescriptorIndex(0), d_rhs_vec->getComponentDescriptorIndex(0), d_F_scratch_idx);
    }

    // Account for internal source/sink distributions.
    if (!d_Q_fcn.isNull())
    {
        d_Q_fcn->setDataOnPatchHierarchy(d_Q_current_idx, d_Q_var, d_hierarchy, current_time);
        d_Q_fcn->setDataOnPatchHierarchy(d_Q_new_idx    , d_Q_var, d_hierarchy, new_time    );
        d_hier_cc_data_ops->linearSum(d_Q_scratch_idx, 0.5, d_Q_current_idx, 0.5, d_Q_new_idx);
        d_Q_bdry_bc_fill_op->fillData(current_time+0.5*dt);
        if (!d_creeping_flow)
        {
            d_hier_sc_data_ops->linearSum(d_U_scratch_idx, 0.5, d_U_current_idx, 0.5, d_U_new_idx);
            computeDivSourceTerm(d_F_div_idx, d_Q_scratch_idx, d_U_scratch_idx);
        }
        d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0), rho, d_F_div_idx, d_rhs_vec->getComponentDescriptorIndex(0));
        d_hier_cc_data_ops->subtract(d_rhs_vec->getComponentDescriptorIndex(1), d_rhs_vec->getComponentDescriptorIndex(1), d_Q_new_idx);
    }

    // Synchronize solution and right-hand-side data before solve.
    typedef SideDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent;
    SynchronizationTransactionComponent sol_synch_transaction = SynchronizationTransactionComponent(d_sol_vec->getComponentDescriptorIndex(0), "CONSERVATIVE_COARSEN");
    d_side_synch_op->resetTransactionComponent(sol_synch_transaction);
    d_side_synch_op->synchronizeData(current_time);
    SynchronizationTransactionComponent rhs_synch_transaction = SynchronizationTransactionComponent(d_rhs_vec->getComponentDescriptorIndex(0), "CONSERVATIVE_COARSEN");
    d_side_synch_op->resetTransactionComponent(rhs_synch_transaction);
    d_side_synch_op->synchronizeData(current_time);

    // Set solution components to equal most recent approximations to u(n+1) and
    // p(n+1/2).
    d_hier_sc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(0), d_U_new_idx);
    d_hier_cc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(1), d_P_new_idx);

    // Ensure there is no forcing at Dirichlet boundaries (the Dirichlet
    // boundary condition takes precedence).
    d_bc_helper->enforceDirichletBcs(d_rhs_vec->getComponentDescriptorIndex(0), /*homogeneous_bcs*/ true);

    // Solve for u(n+1), p(n+1/2).
    d_stokes_solver->solveSystem(*d_sol_vec,*d_rhs_vec);

    // Synchronize solution data after solve.
    d_side_synch_op->resetTransactionComponent(sol_synch_transaction);
    d_side_synch_op->synchronizeData(current_time);

    // Enforce Dirichlet boundary conditions.
    d_bc_helper->enforceDirichletBcs(d_sol_vec->getComponentDescriptorIndex(0), /*homogeneous_bcs*/ false);

    // Pull out solution components.
    d_hier_sc_data_ops->copyData(d_U_new_idx, d_sol_vec->getComponentDescriptorIndex(0));
    d_hier_cc_data_ops->copyData(d_P_new_idx, d_sol_vec->getComponentDescriptorIndex(1));

    // Reset the right-hand side vector.
    if (!d_creeping_flow)
    {
        const int N_idx = d_N_vec->getComponentDescriptorIndex(0);
        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
        {
            d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0), +1.0*rho, N_idx, d_rhs_vec->getComponentDescriptorIndex(0));
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0), +0.5*rho, N_idx, d_rhs_vec->getComponentDescriptorIndex(0));
        }
    }
    if (!d_F_fcn.isNull())
    {
        d_hier_sc_data_ops->subtract(d_rhs_vec->getComponentDescriptorIndex(0), d_rhs_vec->getComponentDescriptorIndex(0), d_F_scratch_idx);
        d_hier_sc_data_ops->copyData(d_F_new_idx, d_F_scratch_idx);
    }
    if (!d_Q_fcn.isNull())
    {
        d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0), -rho, d_F_div_idx, d_rhs_vec->getComponentDescriptorIndex(0));
        d_hier_cc_data_ops->add(d_rhs_vec->getComponentDescriptorIndex(1), d_rhs_vec->getComponentDescriptorIndex(1), d_Q_new_idx);
    }

    if (!d_adv_diff_hier_integrator.isNull())
    {
        // Update the advection velocities used by the advection-diffusion
        // solver.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int U_adv_diff_new_idx     = var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getNewContext()    );
        copy_side_to_face(U_adv_diff_new_idx, d_U_new_idx, d_hierarchy);
        const int U_adv_diff_current_idx = var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getCurrentContext());
        const int U_adv_diff_scratch_idx = var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getScratchContext());
        d_hier_fc_data_ops->linearSum(U_adv_diff_scratch_idx, 0.5, U_adv_diff_current_idx, 0.5, U_adv_diff_new_idx);

        // Update the state variables maintained by the advection-diffusion
        // solver.
        //
        // NOTE: We already performed cycle 0 above.
        const int adv_diff_num_cycles = d_adv_diff_hier_integrator->getNumberOfCycles();
        if (d_current_num_cycles != adv_diff_num_cycles)
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_current_num_cycles == 1);
#endif
            for (int cycle_num = 1; cycle_num < adv_diff_num_cycles; ++cycle_num)
            {
                d_adv_diff_hier_integrator->integrateHierarchy(current_time, new_time, cycle_num);
            }
        }
    }
    return;
}// integrateHierarchy

void
INSStaggeredHierarchyIntegrator::postprocessIntegrateHierarchy(
    const double current_time,
    const double new_time,
    const bool skip_synchronize_new_state_data,
    const int num_cycles)
{
    INSHierarchyIntegrator::postprocessIntegrateHierarchy(current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Synchronize new state data.
    if (!skip_synchronize_new_state_data)
    {
        if (d_enable_logging) plog << d_object_name << "::postprocessIntegrateHierarchy(): synchronizing updated data\n";
        synchronizeHierarchyData(NEW_DATA);
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
    d_U_adv_vec->deallocateVectorData();
    d_N_vec    ->deallocateVectorData();
    d_P_rhs_vec->deallocateVectorData();

    // Deallocate any registered advection-diffusion solver.
    if (!d_adv_diff_hier_integrator.isNull())
    {
        const int adv_diff_num_cycles = d_adv_diff_hier_integrator->getNumberOfCycles();
        d_adv_diff_hier_integrator->postprocessIntegrateHierarchy(current_time, new_time, skip_synchronize_new_state_data, adv_diff_num_cycles);
    }
    return;
}// postprocessIntegrateHierarchy

void
INSStaggeredHierarchyIntegrator::regridHierarchy()
{
    const int coarsest_ln = 0;

    // Determine the divergence of the velocity field before regridding.
    d_hier_math_ops->div(d_Div_U_idx, d_Div_U_var, 1.0, d_U_current_idx, d_U_var, d_no_fill_op, d_integrator_time, /*synch_cf_bdry*/ false, -1.0, d_Q_current_idx, d_Q_var);
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const double Div_U_norm_1_pre  = d_hier_cc_data_ops->L1Norm( d_Div_U_idx, wgt_cc_idx);
    const double Div_U_norm_2_pre  = d_hier_cc_data_ops->L2Norm( d_Div_U_idx, wgt_cc_idx);
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
                       << "  unrecognized regrid mode: " << IBTK::enum_to_string<RegridMode>(d_regrid_mode) << "." << std::endl);
    }

    // Determine the divergence of the velocity field after regridding.
    d_hier_math_ops->div(d_Div_U_idx, d_Div_U_var, 1.0, d_U_current_idx, d_U_var, d_no_fill_op, d_integrator_time, /*synch_cf_bdry*/ true, -1.0, d_Q_current_idx, d_Q_var);
    const double Div_U_norm_1_post  = d_hier_cc_data_ops->L1Norm( d_Div_U_idx, wgt_cc_idx);
    const double Div_U_norm_2_post  = d_hier_cc_data_ops->L2Norm( d_Div_U_idx, wgt_cc_idx);
    const double Div_U_norm_oo_post = d_hier_cc_data_ops->maxNorm(d_Div_U_idx, wgt_cc_idx);

    // Project the interpolated velocity if needed.
    if (Div_U_norm_1_post  > d_regrid_max_div_growth_factor*Div_U_norm_1_pre ||
        Div_U_norm_2_post  > d_regrid_max_div_growth_factor*Div_U_norm_2_pre ||
        Div_U_norm_oo_post > d_regrid_max_div_growth_factor*Div_U_norm_oo_pre)
    {
        pout << d_object_name << "::regridHierarchy():\n"
             << "  WARNING: projecting the interpolated velocity field\n";
        regridProjection();
    }

    // Synchronize the state data on the patch hierarchy.
    synchronizeHierarchyData(CURRENT_DATA);
    return;
}// regridHierarchy

/////////////////////////////// PROTECTED ////////////////////////////////////

void
INSStaggeredHierarchyIntegrator::initializeLevelDataSpecialized(
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
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Correct the divergence of the interpolated velocity data.
    if (!initial_time && level_number > 0)
    {
        // Allocate scratch data.
        ComponentSelector scratch_data;
        scratch_data.setFlag( d_U_regrid_idx);
        scratch_data.setFlag(    d_U_src_idx);
        scratch_data.setFlag(d_indicator_idx);
        level->allocatePatchData(scratch_data, init_data_time);
        if (!old_level.isNull()) old_level->allocatePatchData(scratch_data, init_data_time);

        // Set the indicator data to equal "0" in each patch of the new patch
        // level, and initialize values of U to cause floating point errors if
        // we fail to re-initialize it properly.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM,double> > indicator_data = patch->getPatchData(d_indicator_idx);
            indicator_data->fillAll(0.0);

            Pointer<SideData<NDIM,double> > U_current_data = patch->getPatchData(d_U_current_idx);
            Pointer<SideData<NDIM,double> >  U_regrid_data = patch->getPatchData( d_U_regrid_idx);
            Pointer<SideData<NDIM,double> >     U_src_data = patch->getPatchData(    d_U_src_idx);
            U_current_data->fillAll(std::numeric_limits<double>::quiet_NaN());
            U_regrid_data ->fillAll(std::numeric_limits<double>::quiet_NaN());
            U_src_data    ->fillAll(std::numeric_limits<double>::quiet_NaN());
        }

        if (!old_level.isNull())
        {
            // Set the indicator data to equal "1" on each patch of the old
            // patch level and reset U.
            for (PatchLevel<NDIM>::Iterator p(old_level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = old_level->getPatch(p());

                Pointer<SideData<NDIM,double> > indicator_data = patch->getPatchData(d_indicator_idx);
                indicator_data->fillAll(1.0);

                Pointer<SideData<NDIM,double> > U_current_data = patch->getPatchData(d_U_current_idx);
                Pointer<SideData<NDIM,double> >  U_regrid_data = patch->getPatchData( d_U_regrid_idx);
                Pointer<SideData<NDIM,double> >     U_src_data = patch->getPatchData(    d_U_src_idx);
                U_regrid_data->copy(*U_current_data);
                U_src_data   ->copy(*U_current_data);
            }

            // Create a communications schedule to copy data from the old patch
            // level to the new patch level.
            //
            // Note that this will set the indicator data to equal "1" at each
            // location in the new patch level that is a copy of a location from
            // the old patch level.
            RefineAlgorithm<NDIM> copy_data;
            copy_data.registerRefine( d_U_regrid_idx,  d_U_regrid_idx,  d_U_regrid_idx, NULL);
            copy_data.registerRefine(    d_U_src_idx,     d_U_src_idx,     d_U_src_idx, NULL);
            copy_data.registerRefine(d_indicator_idx, d_indicator_idx, d_indicator_idx, NULL);
            ComponentSelector bc_fill_data;
            bc_fill_data.setFlag(d_U_regrid_idx);
            bc_fill_data.setFlag(   d_U_src_idx);
            CartSideRobinPhysBdryOp phys_bdry_bc_op(bc_fill_data, d_U_bc_coefs, false);
            copy_data.createSchedule(level, old_level, &phys_bdry_bc_op)->fillData(init_data_time);
        }

        // Setup the divergence- and curl-preserving prolongation refine
        // algorithm and refine the velocity data.
        RefineAlgorithm<NDIM> fill_div_free_prolongation;
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        fill_div_free_prolongation.registerRefine(d_U_current_idx, d_U_current_idx, d_U_regrid_idx, NULL);
        Pointer<RefineOperator<NDIM> > refine_op = grid_geom->lookupRefineOperator(d_U_var, "SPECIALIZED_LINEAR_REFINE");
        Pointer<CoarsenOperator<NDIM> > coarsen_op = grid_geom->lookupCoarsenOperator(d_U_var, "CONSERVATIVE_COARSEN");
        CartSideRobinPhysBdryOp phys_bdry_bc_op(d_U_regrid_idx, d_U_bc_coefs, false);
        CartSideDoubleDivPreservingRefine div_preserving_op(d_U_regrid_idx, d_U_src_idx, d_indicator_idx, refine_op, coarsen_op, init_data_time, &phys_bdry_bc_op);
        fill_div_free_prolongation.createSchedule(level, old_level, level_number-1, hierarchy, &div_preserving_op)->fillData(init_data_time);

        // Free scratch data.
        level->deallocatePatchData(scratch_data);
        if (!old_level.isNull()) old_level->deallocatePatchData(scratch_data);
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
            Pointer<HierarchyCellDataOpsReal<NDIM,double> > hier_cc_data_ops = hier_ops_manager->getOperationsDouble(d_U_cc_var, d_hierarchy, true);
            Pointer<HierarchySideDataOpsReal<NDIM,double> > hier_sc_data_ops = hier_ops_manager->getOperationsDouble(d_U_var, d_hierarchy, true);
            hier_sc_data_ops->resetLevels(0, level_number);
            hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
            typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
            InterpolationTransactionComponent U_bc_component(d_U_scratch_idx, SIDE_DATA_COARSEN_TYPE, d_bdry_extrap_type, CONSISTENT_TYPE_2_BDRY, d_U_bc_coefs);
            HierarchyGhostCellInterpolation U_bdry_bc_fill_op;
            U_bdry_bc_fill_op.initializeOperatorState(U_bc_component, d_hierarchy, 0, level_number);
            U_bdry_bc_fill_op.fillData(init_data_time);

            // Compute max |Omega|_2.
            HierarchyMathOps hier_math_ops(d_object_name+"::HierarchyLevelMathOps", d_hierarchy, 0, level_number);
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
}// initializeLevelDataSpecialized

void
INSStaggeredHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    const Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0) && (coarsest_level <= finest_level) && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
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
    InterpolationTransactionComponent U_bc_component(d_U_scratch_idx, SIDE_DATA_COARSEN_TYPE, d_bdry_extrap_type, CONSISTENT_TYPE_2_BDRY, d_U_bc_coefs);
    d_U_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_U_bdry_bc_fill_op->initializeOperatorState(U_bc_component, d_hierarchy);

    InterpolationTransactionComponent P_bc_component(d_P_scratch_idx, CELL_DATA_COARSEN_TYPE, d_bdry_extrap_type, CONSISTENT_TYPE_2_BDRY);
    d_P_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_P_bdry_bc_fill_op->initializeOperatorState(P_bc_component, d_hierarchy);

    if (!d_Q_fcn.isNull())
    {
        InterpolationTransactionComponent Q_bc_component(d_Q_scratch_idx, CELL_DATA_COARSEN_TYPE, d_bdry_extrap_type, CONSISTENT_TYPE_2_BDRY);
        d_Q_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
        d_Q_bdry_bc_fill_op->initializeOperatorState(Q_bc_component, d_hierarchy);
    }

    // Setup the patch boundary synchronization objects.
    typedef SideDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent;
    SynchronizationTransactionComponent synch_transaction = SynchronizationTransactionComponent(d_U_scratch_idx, "CONSERVATIVE_COARSEN");
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
}// resetHierarchyConfigurationSpecialized

void
INSStaggeredHierarchyIntegrator::applyGradientDetectorSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double /*error_data_time*/,
    const int tag_index,
    const bool /*initial_time*/,
    const bool /*uses_richardson_extrapolation_too*/)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
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
            Omega_rel_thresh = d_Omega_rel_thresh[std::max(std::min(level_number,d_Omega_rel_thresh.size()-1),0)];
        }
        double Omega_abs_thresh = 0.0;
        if (d_Omega_abs_thresh.size() > 0)
        {
            Omega_abs_thresh = d_Omega_abs_thresh[std::max(std::min(level_number,d_Omega_abs_thresh.size()-1),0)];
        }
        if (Omega_rel_thresh > 0.0 || Omega_abs_thresh > 0.0)
        {
            double thresh = std::numeric_limits<double>::max();
            if (Omega_rel_thresh > 0.0) thresh = std::min(thresh, Omega_rel_thresh*d_Omega_max);
            if (Omega_abs_thresh > 0.0) thresh = std::min(thresh, Omega_abs_thresh            );
            thresh += sqrt(std::numeric_limits<double>::epsilon());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
                Pointer<CellData<NDIM,double> > Omega_data = patch->getPatchData(d_Omega_idx);
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
                        norm_Omega_sq += (*Omega_data)(i,d)*(*Omega_data)(i,d);
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
}// applyGradientDetectorSpecialized

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
        d_hier_math_ops->interp(U_cc_idx, d_U_cc_var, U_sc_idx, d_U_var, d_no_fill_op, d_integrator_time, synch_cf_interface);
    }

    // Interpolate f to cell centers.
    if (d_output_F && !d_F_fcn.isNull())
    {
        const int F_sc_idx = var_db->mapVariableAndContextToIndex(d_F_var, ctx);
        const int F_cc_idx = var_db->mapVariableAndContextToIndex(d_F_cc_var, ctx);
        d_hier_math_ops->interp(F_cc_idx, d_F_cc_var, F_sc_idx, d_F_var, d_no_fill_op, d_integrator_time, synch_cf_interface);
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
        d_hier_math_ops->curl(d_Omega_idx, d_Omega_var, d_U_scratch_idx, d_U_var, d_U_bdry_bc_fill_op, d_integrator_time);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(d_U_scratch_idx);
        }
    }

    // Compute Div U.
    if (d_output_Div_U)
    {
        d_hier_math_ops->div(d_Div_U_idx, d_Div_U_var, 1.0, d_U_current_idx, d_U_var, d_no_fill_op, d_integrator_time, false);
    }
    return;
}// setupPlotDataSpecialized

void
INSStaggeredHierarchyIntegrator::regridProjection()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const double volume = d_hier_math_ops->getVolumeOfPhysicalDomain();

    // Setup the solver vectors.
    SAMRAIVectorReal<NDIM,double> sol_vec(d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec.addComponent(d_P_var, d_P_scratch_idx, wgt_cc_idx, d_hier_cc_data_ops);

    SAMRAIVectorReal<NDIM,double> rhs_vec(d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    rhs_vec.addComponent(d_Div_U_var, d_Div_U_idx, wgt_cc_idx, d_hier_cc_data_ops);

    // Setup the regrid Poisson solver.
    Pointer<PoissonSolver> regrid_projection_solver = CCPoissonSolverManager::getManager()->allocateSolver(
            d_regrid_projection_solver_type , d_object_name+"::regrid_projection_solver" , d_regrid_projection_solver_db ,
            d_regrid_projection_precond_type, d_object_name+"::regrid_projection_precond", d_regrid_projection_precond_db);
    PoissonSpecifications regrid_projection_spec(d_object_name+"::regrid_projection_spec");
    regrid_projection_spec.setCZero();
    regrid_projection_spec.setDConstant(-1.0);
    LocationIndexRobinBcCoefs<NDIM> Phi_bc_coef;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        Phi_bc_coef.setBoundarySlope(2*d  ,0.0);
        Phi_bc_coef.setBoundarySlope(2*d+1,0.0);
    }
    regrid_projection_solver->setPoissonSpecifications(regrid_projection_spec);
    regrid_projection_solver->setPhysicalBcCoef(&Phi_bc_coef);
    regrid_projection_solver->setHomogeneousBc(true);
    regrid_projection_solver->setSolutionTime(d_integrator_time);
    regrid_projection_solver->setTimeInterval(d_integrator_time, d_integrator_time);
    regrid_projection_solver->setNullspace(true);
    regrid_projection_solver->setInitialGuessNonzero(false);

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
    d_hier_math_ops->div(d_Div_U_idx, d_Div_U_var, -1.0, d_U_current_idx, d_U_var, d_no_fill_op, d_integrator_time, /*synch_cf_bdry*/ false, +1.0, d_Q_current_idx, d_Q_var);
    const double Div_U_mean = (1.0/volume)*d_hier_cc_data_ops->integral(d_Div_U_idx, wgt_cc_idx);
    d_hier_cc_data_ops->addScalar(d_Div_U_idx, d_Div_U_idx, -Div_U_mean);

    // Solve the projection pressure-Poisson problem.
    regrid_projection_solver->solveSystem(sol_vec,rhs_vec);

    // Fill ghost cells for Phi, compute Grad Phi, and set U := U - Grad Phi
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent Phi_bc_component(d_P_scratch_idx, CELL_DATA_COARSEN_TYPE, d_bdry_extrap_type, CONSISTENT_TYPE_2_BDRY, &Phi_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> Phi_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    Phi_bdry_bc_fill_op->initializeOperatorState(Phi_bc_component, d_hierarchy);
    Phi_bdry_bc_fill_op->setHomogeneousBc(true);
    Phi_bdry_bc_fill_op->fillData(d_integrator_time);
    d_hier_math_ops->grad(d_U_current_idx, d_U_var, /*synch_cf_bdry*/ true, -1.0, d_P_scratch_idx, d_P_var, d_no_fill_op, d_integrator_time, +1.0, d_U_current_idx, d_U_var);

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(scratch_idxs);
    }
    return;
}// regridProjection

double
INSStaggeredHierarchyIntegrator::getStableTimestep(
    Pointer<Patch<NDIM> > patch) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const Index<NDIM>& ilower = patch->getBox().lower();
    const Index<NDIM>& iupper = patch->getBox().upper();

    Pointer<SideData<NDIM,double> > U_data = patch->getPatchData(d_U_var, getCurrentContext());
    const IntVector<NDIM>& U_ghost_cells = U_data->getGhostCellWidth();

    double stable_dt = std::numeric_limits<double>::max();
    NAVIER_STOKES_SC_STABLEDT_FC(
        dx,
#if (NDIM == 2)
        ilower(0),iupper(0),ilower(1),iupper(1),
        U_ghost_cells(0),U_ghost_cells(1),
        U_data->getPointer(0),U_data->getPointer(1),
#endif
#if (NDIM == 3)
        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
        U_ghost_cells(0),U_ghost_cells(1),U_ghost_cells(2),
        U_data->getPointer(0),U_data->getPointer(1),U_data->getPointer(2),
#endif
        stable_dt);
    return stable_dt;
}// getStableTimestep

/////////////////////////////// PRIVATE //////////////////////////////////////

void
INSStaggeredHierarchyIntegrator::reinitializeOperatorsAndSolvers(
    const double current_time,
    const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    const double dt = new_time-current_time;
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const int wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
    const double rho    = d_problem_coefs.getRho();
    const double mu     = d_problem_coefs.getMu();
    const double lambda = d_problem_coefs.getLambda();

    // Ensure that solver components are appropriately reinitialized when the
    // time step size changes.
    const bool dt_change = initial_time || !MathUtilities<double>::equalEps(dt,d_dt_previous[0]);
    if (dt_change)
    {
        d_velocity_solver_needs_init = true;
        d_stokes_solver_needs_init = true;
    }

    // Setup solver vectors.
    if (d_vectors_need_init)
    {
        d_U_scratch_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::U_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_U_scratch_vec->addComponent(d_U_var, d_U_scratch_idx, wgt_sc_idx, d_hier_sc_data_ops);

        d_P_scratch_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::P_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_P_scratch_vec->addComponent(d_P_var, d_P_scratch_idx, wgt_cc_idx, d_hier_cc_data_ops);

        if (!d_U_rhs_vec.isNull()) d_U_rhs_vec->freeVectorComponents();
        if (!d_U_adv_vec.isNull()) d_U_adv_vec->freeVectorComponents();
        if (!d_N_vec    .isNull()) d_N_vec    ->freeVectorComponents();
        if (!d_P_rhs_vec.isNull()) d_P_rhs_vec->freeVectorComponents();

        d_U_rhs_vec = d_U_scratch_vec->cloneVector(d_object_name+"::U_rhs_vec");
        d_U_adv_vec = d_U_scratch_vec->cloneVector(d_object_name+"::U_adv_vec");
        d_N_vec     = d_U_scratch_vec->cloneVector(d_object_name+"::N_vec"    );
        d_P_rhs_vec = d_P_scratch_vec->cloneVector(d_object_name+"::P_rhs_vec");

        d_sol_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_sol_vec->addComponent(d_U_var,d_U_scratch_idx,wgt_sc_idx,d_hier_sc_data_ops);
        d_sol_vec->addComponent(d_P_var,d_P_scratch_idx,wgt_cc_idx,d_hier_cc_data_ops);

        d_rhs_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
        const int U_rhs_idx = d_U_rhs_vec->getComponentDescriptorIndex(0);
        d_rhs_vec->addComponent(d_U_var,U_rhs_idx,wgt_sc_idx,d_hier_sc_data_ops);
        const int P_rhs_idx = d_P_rhs_vec->getComponentDescriptorIndex(0);
        d_rhs_vec->addComponent(d_P_var,P_rhs_idx,wgt_cc_idx,d_hier_cc_data_ops);

        d_nul_vecs.resize((d_normalize_pressure ? 1 : 0) + (MathUtilities<double>::equalEps(rho, 0.0) ? NDIM : 0));
        d_U_nul_vecs.resize(MathUtilities<double>::equalEps(rho, 0.0) ? NDIM : 0);
        if (d_normalize_pressure)
        {
            if (!d_nul_vecs[0].isNull()) d_nul_vecs[0]->freeVectorComponents();
            d_nul_vecs[0] = d_sol_vec->cloneVector(d_object_name+"::nul_vec_p");
            d_nul_vecs[0]->allocateVectorData(current_time);
            d_hier_sc_data_ops->setToScalar(d_nul_vecs[0]->getComponentDescriptorIndex(0), 0.0);
            d_hier_cc_data_ops->setToScalar(d_nul_vecs[0]->getComponentDescriptorIndex(1), 1.0);
        }
        if (MathUtilities<double>::equalEps(rho, 0.0))
        {
            const int offset = (d_normalize_pressure ? 1 : 0);
            for (unsigned int k = 0; k < NDIM; ++k)
            {
                if (!d_nul_vecs[k+offset].isNull()) d_nul_vecs[k+offset]->freeVectorComponents();
                std::ostringstream stream;
                stream << k;
                d_nul_vecs[k+offset] = d_sol_vec->cloneVector(d_object_name+"::nul_vec_U_"+stream.str());
                d_nul_vecs[k+offset]->allocateVectorData(current_time);
                for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
                {
                    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        Pointer<SideData<NDIM,double> > nul_data = patch->getPatchData(d_nul_vecs[k+offset]->getComponentDescriptorIndex(0));
                        nul_data->fillAll(0.0);
                        nul_data->getArrayData(k).fillAll(1.0);
                    }
                }
            }
            for (unsigned int k = 0; k < NDIM; ++k)
            {
                if (!d_U_nul_vecs[k].isNull()) d_U_nul_vecs[k]->freeVectorComponents();
                std::ostringstream stream;
                stream << k;
                d_U_nul_vecs[k] = d_U_scratch_vec->cloneVector(d_object_name+"::U_nul_vec_U"+stream.str());
                d_U_nul_vecs[k]->allocateVectorData(current_time);
                for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
                {
                    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        Pointer<SideData<NDIM,double> > U_nul_data = patch->getPatchData(d_U_nul_vecs[k]->getComponentDescriptorIndex(0));
                        U_nul_data->fillAll(0.0);
                        U_nul_data->getArrayData(k).fillAll(1.0);
                    }
                }
            }
        }

        d_vectors_need_init = false;
    }

    // Setup boundary conditions objects.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        INSStaggeredVelocityBcCoef* U_bc_coef = dynamic_cast<INSStaggeredVelocityBcCoef*>(d_U_bc_coefs[d]);
        U_bc_coef->setPhysicalBoundaryConditions(d_bc_coefs);
        U_bc_coef->setTimeInterval(current_time,new_time);
    }
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        INSIntermediateVelocityBcCoef* U_star_bc_coef = dynamic_cast<INSIntermediateVelocityBcCoef*>(d_U_star_bc_coefs[d]);
        U_star_bc_coef->setPhysicalBoundaryConditions(d_bc_coefs);
        U_star_bc_coef->setTimeInterval(current_time,new_time);
    }
    INSStaggeredPressureBcCoef* P_bc_coef = dynamic_cast<INSStaggeredPressureBcCoef*>(d_P_bc_coef);
    P_bc_coef->setPhysicalBoundaryConditions(d_bc_coefs);
    P_bc_coef->setVelocityCurrentPatchDataIndex(d_U_current_idx);
    P_bc_coef->setVelocityNewPatchDataIndex(d_U_new_idx);
    P_bc_coef->setTimeInterval(current_time,new_time);
    INSProjectionBcCoef* Phi_bc_coef = dynamic_cast<INSProjectionBcCoef*>(d_Phi_bc_coef);
    Phi_bc_coef->setPhysicalBoundaryConditions(d_bc_coefs);
    Phi_bc_coef->setTimeInterval(current_time,new_time);

    // Setup convective operator.
    if (!d_convective_op.isNull() && d_convective_op_needs_init)
    {
        if (d_enable_logging) plog << d_object_name << "::preprocessIntegrateHierarchy(): initializing convective operator" << std::endl;
        d_convective_op->setAdvectionVelocity(d_U_scratch_idx);
        d_convective_op->initializeOperatorState(*d_U_scratch_vec,*d_U_rhs_vec);
        d_convective_op_needs_init = false;
    }

    // Setup subdomain solvers.
    if (!d_velocity_solver.isNull())
    {
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
        PoissonSpecifications velocity_solver_spec(d_object_name+"::velocity_solver_spec");
        velocity_solver_spec.setCConstant((rho/dt)+K*lambda);
        velocity_solver_spec.setDConstant(        -K*mu    );
        d_velocity_solver->setPoissonSpecifications(velocity_solver_spec);
        d_velocity_solver->setPhysicalBcCoefs(d_U_star_bc_coefs);
        d_velocity_solver->setSolutionTime(new_time);
        d_velocity_solver->setTimeInterval(current_time, new_time);
        d_velocity_solver->setInitialGuessNonzero(true);
        if (!d_U_nul_vecs.empty())
        {
            d_velocity_solver->setNullspace(false, d_U_nul_vecs);
        }
        if (d_velocity_solver_needs_init)
        {
            if (d_enable_logging) plog << d_object_name << "::preprocessIntegrateHierarchy(): initializing velocity subdomain solver" << std::endl;
            d_velocity_solver->initializeSolverState(*d_U_scratch_vec,*d_U_rhs_vec);
            d_velocity_solver_needs_init = false;
        }
    }

    if (!d_pressure_solver.isNull())
    {
        PoissonSpecifications pressure_solver_spec(d_object_name+"::pressure_solver_spec");
        pressure_solver_spec.setCZero();
        pressure_solver_spec.setDConstant(-1.0);
        d_pressure_solver->setPoissonSpecifications(pressure_solver_spec);
        d_pressure_solver->setPhysicalBcCoef(d_Phi_bc_coef);
        d_pressure_solver->setSolutionTime(current_time+0.5*dt);
        d_pressure_solver->setTimeInterval(current_time, new_time);
        d_pressure_solver->setInitialGuessNonzero(true);
        if (d_normalize_pressure)
        {
            d_pressure_solver->setNullspace(true);
        }
        if (d_pressure_solver_needs_init)
        {
            if (d_enable_logging) plog << d_object_name << "::preprocessIntegrateHierarchy(): Initializing pressure subdomain solver" << std::endl;
            d_pressure_solver->initializeSolverState(*d_P_scratch_vec,*d_P_rhs_vec);
            d_pressure_solver_needs_init = false;
        }
    }

    // Setup Stokes solver.
    d_stokes_solver->setSolutionTime(new_time);
    d_stokes_solver->setTimeInterval(current_time,new_time);
    if (!d_nul_vecs.empty())
    {
        Pointer<LinearSolver> p_stokes_linear_solver = d_stokes_solver;
        if (!p_stokes_linear_solver.isNull())
        {
            Pointer<NewtonKrylovSolver> p_stokes_newton_solver = d_stokes_solver;
            if (!p_stokes_newton_solver.isNull()) p_stokes_linear_solver = p_stokes_newton_solver->getLinearSolver();
        }
        if (!p_stokes_linear_solver.isNull())
        {
            p_stokes_linear_solver->setNullspace(false, d_nul_vecs);
        }
    }
    if (d_stokes_solver_needs_init)
    {
        if (d_enable_logging) plog << d_object_name << "::preprocessIntegrateHierarchy(): initializing incompressible Stokes solver" << std::endl;
        d_stokes_solver->initializeSolverState(*d_sol_vec,*d_rhs_vec);
        d_stokes_solver_needs_init = false;
    }
    return;
}// reinitializeOperatorsAndSolvers

void
INSStaggeredHierarchyIntegrator::computeDivSourceTerm(
    const int F_idx,
    const int Q_idx,
    const int U_idx)
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

            Pointer<SideData<NDIM,double> > U_data = patch->getPatchData(U_idx);
            Pointer<CellData<NDIM,double> > Q_data = patch->getPatchData(Q_idx);
            Pointer<SideData<NDIM,double> > F_data = patch->getPatchData(F_idx);

            const IntVector<NDIM>& U_data_gc = U_data->getGhostCellWidth();
            const IntVector<NDIM>& Q_data_gc = Q_data->getGhostCellWidth();
            const IntVector<NDIM>& F_data_gc = F_data->getGhostCellWidth();

            switch (d_convective_op->getConvectiveDifferencingType())
            {
                case CONSERVATIVE:
                    NAVIER_STOKES_STAGGERED_CONS_SOURCE_FC(
#if (NDIM == 2)
                        ilower(0),iupper(0),ilower(1),iupper(1),
                        U_data_gc(0),U_data_gc(1),
                        Q_data_gc(0),Q_data_gc(1),
                        F_data_gc(0),F_data_gc(1),
                        U_data->getPointer(0),U_data->getPointer(1),
                        Q_data->getPointer(),
                        F_data->getPointer(0),F_data->getPointer(1)
#endif
#if (NDIM == 3)
                        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                        U_data_gc(0),U_data_gc(1),U_data_gc(2),
                        Q_data_gc(0),Q_data_gc(1),Q_data_gc(2),
                        F_data_gc(0),F_data_gc(1),F_data_gc(2),
                        U_data->getPointer(0),U_data->getPointer(1),U_data->getPointer(2),
                        Q_data->getPointer(),
                        F_data->getPointer(0),F_data->getPointer(1),F_data->getPointer(2)
#endif
                                                           );
                    break;
                case ADVECTIVE:
                    NAVIER_STOKES_STAGGERED_ADV_SOURCE_FC(
#if (NDIM == 2)
                        ilower(0),iupper(0),ilower(1),iupper(1),
                        U_data_gc(0),U_data_gc(1),
                        Q_data_gc(0),Q_data_gc(1),
                        F_data_gc(0),F_data_gc(1),
                        U_data->getPointer(0),U_data->getPointer(1),
                        Q_data->getPointer(),
                        F_data->getPointer(0),F_data->getPointer(1)
#endif
#if (NDIM == 3)
                        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                        U_data_gc(0),U_data_gc(1),U_data_gc(2),
                        Q_data_gc(0),Q_data_gc(1),Q_data_gc(2),
                        F_data_gc(0),F_data_gc(1),F_data_gc(2),
                        U_data->getPointer(0),U_data->getPointer(1),U_data->getPointer(2),
                        Q_data->getPointer(),
                        F_data->getPointer(0),F_data->getPointer(1),F_data->getPointer(2)
#endif
                                                          );
                    break;
                case SKEW_SYMMETRIC:
                    NAVIER_STOKES_STAGGERED_SKEW_SYM_SOURCE_FC(
#if (NDIM == 2)
                        ilower(0),iupper(0),ilower(1),iupper(1),
                        U_data_gc(0),U_data_gc(1),
                        Q_data_gc(0),Q_data_gc(1),
                        F_data_gc(0),F_data_gc(1),
                        U_data->getPointer(0),U_data->getPointer(1),
                        Q_data->getPointer(),
                        F_data->getPointer(0),F_data->getPointer(1)
#endif
#if (NDIM == 3)
                        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                        U_data_gc(0),U_data_gc(1),U_data_gc(2),
                        Q_data_gc(0),Q_data_gc(1),Q_data_gc(2),
                        F_data_gc(0),F_data_gc(1),F_data_gc(2),
                        U_data->getPointer(0),U_data->getPointer(1),U_data->getPointer(2),
                        Q_data->getPointer(),
                        F_data->getPointer(0),F_data->getPointer(1),F_data->getPointer(2)
#endif
                                                               );
                    break;
                default:
                    TBOX_ERROR("INSStaggeredHierarchyIntegrator::computeDivSourceTerm():\n"
                               << "  unsupported differencing form: " << enum_to_string<ConvectiveDifferencingType>(d_convective_op->getConvectiveDifferencingType()) << " \n"
                               << "  valid choices are: ADVECTIVE, CONSERVATIVE, SKEW_SYMMETRIC\n");
            }
        }
    }
    return;
}// computeDivSourceTerm

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
