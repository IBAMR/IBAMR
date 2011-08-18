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
#include <ibamr/INSStaggeredIntermediateVelocityBcCoef.h>
#include <ibamr/INSStaggeredBlockFactorizationPreconditioner.h>
#include <ibamr/INSStaggeredPPMConvectiveOperator.h>
#include <ibamr/INSStaggeredProjectionBcCoef.h>
#include <ibamr/INSStaggeredProjectionPreconditioner.h>
#include <ibamr/INSStaggeredPressureBcCoef.h>
#include <ibamr/INSStaggeredVelocityBcCoef.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/CartSideDoubleDivPreservingRefine.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/RefinePatchStrategySet.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

// C++ STDLIB INCLUDES
#include <limits>

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

// Version of INSStaggeredHierarchyIntegrator restart file data.
static const int INS_STAGGERED_HIERARCHY_INTEGRATOR_VERSION = 2;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredHierarchyIntegrator::INSStaggeredHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    bool register_for_restart)
    : INSHierarchyIntegrator(object_name, input_db, register_for_restart)
{
    // Set some default values.
    d_integrator_is_initialized = false;
    d_num_cycles = 3;
    d_cfl_max = 1.0;
    d_using_vorticity_tagging = false;
    d_Omega_max = 0.0;
    d_normalize_pressure = false;
    d_convective_difference_form = ADVECTIVE;
    d_creeping_flow = false;
    d_regrid_max_div_growth_factor = 1.1;
    d_U_scale = 1.0;
    d_P_scale = 1.0;
    d_F_scale = 1.0;
    d_Q_scale = 1.0;
    d_Omega_scale = 1.0;
    d_Div_U_scale = 1.0;
    d_output_U = true;
    d_output_P = true;
    d_output_F = false;
    d_output_Q = false;
    d_output_Omega = true;
    d_output_Div_U = true;

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_U_bc_coefs[d] = NULL;
        d_U_star_bc_coefs[d] = NULL;
    }
    d_P_bc_coef = NULL;
    d_Phi_bc_coef = NULL;

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);
    if (from_restart)
    {
        getFromRestart();
    }
    else
    {
        d_integrator_time = d_start_time;
        d_integrator_step = 0;
    }

    // Initialize all variables.
    d_U_sc_var       = new SideVariable<NDIM,double>(d_object_name+"::U"          );
    d_U_cc_var       = new CellVariable<NDIM,double>(d_object_name+"::U_cc"  ,NDIM);
    d_P_cc_var       = new CellVariable<NDIM,double>(d_object_name+"::P"          );
    d_F_sc_var       = new SideVariable<NDIM,double>(d_object_name+"::F"          );
    d_F_cc_var       = new CellVariable<NDIM,double>(d_object_name+"::F_cc"  ,NDIM);
    d_Q_cc_var       = new CellVariable<NDIM,double>(d_object_name+"::Q"          );
#if (NDIM == 2)
    d_Omega_var      = new CellVariable<NDIM,double>(d_object_name+"::Omega"      );
#endif
#if (NDIM == 3)
    d_Omega_var      = new CellVariable<NDIM,double>(d_object_name+"::Omega" ,NDIM);
    d_Omega_Norm_var = new CellVariable<NDIM,double>(d_object_name+"::||Omega||_2");
#endif
    d_Div_U_var      = new CellVariable<NDIM,double>(d_object_name+"::Div_U"      );
    d_U_regrid_var   = new SideVariable<NDIM,double>(d_object_name+"::U_regrid"   );
    d_U_src_var      = new SideVariable<NDIM,double>(d_object_name+"::U_src"      );
    d_indicator_var  = new SideVariable<NDIM,double>(d_object_name+"::indicator"  );
    d_F_div_var      = new SideVariable<NDIM,double>(d_object_name+"::F_div"      );

    // Initialize variables maintained by the INSHierarchyIntegrator base class.
    d_U_var = d_U_sc_var;
    d_P_var = d_P_cc_var;
    d_F_var = d_F_sc_var;
    d_Q_var = d_Q_cc_var;
    return;
}// INSStaggeredHierarchyIntegrator

INSStaggeredHierarchyIntegrator::~INSStaggeredHierarchyIntegrator()
{
    if (d_helmholtz_spec != NULL) delete d_helmholtz_spec;
    if (d_poisson_spec   != NULL) delete d_poisson_spec;
    if (!d_U_rhs_vec .isNull()) d_U_rhs_vec ->freeVectorComponents();
    if (!d_U_half_vec.isNull()) d_U_half_vec->freeVectorComponents();
    if (!d_N_vec     .isNull()) d_N_vec     ->freeVectorComponents();
    if (!d_P_rhs_vec .isNull()) d_P_rhs_vec ->freeVectorComponents();
    if (!d_nul_vec   .isNull()) d_nul_vec   ->freeVectorComponents();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (d_U_bc_coefs[d] != NULL) delete d_U_bc_coefs[d];
        if (d_U_star_bc_coefs[d] != NULL) delete d_U_star_bc_coefs[d];
    }
    if (d_P_bc_coef != NULL) delete d_P_bc_coef;
    if (d_Phi_bc_coef != NULL) delete d_Phi_bc_coef;
    return;
}// ~INSStaggeredHierarchyIntegrator

void
INSStaggeredHierarchyIntegrator::setConvectiveOperator(
    Pointer<GeneralOperator> convective_op)
{
    d_convective_op = convective_op;
    if (d_convective_op.isNull()) d_creeping_flow = true;
    return;
}// setConvectiveOperator

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

    Pointer<CellVariable<NDIM,double> > cc_var = new CellVariable<NDIM,double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, hierarchy, true);

    Pointer<SideVariable<NDIM,double> > sc_var = new SideVariable<NDIM,double>("sc_var");
    d_hier_sc_data_ops = hier_ops_manager->getOperationsDouble(sc_var, hierarchy, true);

    d_hier_math_ops = buildHierarchyMathOps(d_hierarchy);

    // Register state variables that are maintained by the
    // INSStaggeredHierarchyIntegrator.
    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> side_ghosts = SIDEG;
    const IntVector<NDIM> no_ghosts = 0;

    registerVariable(d_U_current_idx, d_U_new_idx, d_U_scratch_idx,
                     d_U_sc_var, side_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE",
                     d_U_init);

    registerVariable(d_U_cc_current_idx, d_U_cc_new_idx, d_U_cc_scratch_idx,
                     d_U_cc_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_P_current_idx, d_P_new_idx, d_P_scratch_idx,
                     d_P_cc_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE",
                     d_P_init);

    if (!d_F_fcn.isNull())
    {
        registerVariable(d_F_current_idx, d_F_new_idx, d_F_scratch_idx,
                         d_F_sc_var, side_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE",
                         d_F_fcn);

        registerVariable(d_F_cc_current_idx, d_F_cc_new_idx, d_F_cc_scratch_idx,
                         d_F_cc_var, cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");
    }
    else
    {
        d_F_current_idx = -1;
        d_F_new_idx     = -1;
        d_F_scratch_idx = -1;
        d_F_cc_current_idx = -1;
        d_F_cc_new_idx     = -1;
        d_F_cc_scratch_idx = -1;
    }

    if (!d_Q_fcn.isNull())
    {
        registerVariable(d_Q_current_idx, d_Q_new_idx, d_Q_scratch_idx,
                         d_Q_cc_var, cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSTANT_REFINE",
                         d_Q_fcn);
    }
    else
    {
        d_Q_current_idx = -1;
        d_Q_new_idx     = -1;
        d_Q_scratch_idx = -1;
    }

    registerVariable(d_Omega_current_idx, d_Omega_new_idx, d_Omega_scratch_idx,
                     d_Omega_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE");
#if (NDIM == 3)
    registerVariable(d_Omega_Norm_current_idx, d_Omega_Norm_new_idx, d_Omega_Norm_scratch_idx,
                     d_Omega_Norm_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "LINEAR_REFINE");
#endif
    registerVariable(d_Div_U_current_idx, d_Div_U_new_idx, d_Div_U_scratch_idx,
                     d_Div_U_var, cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSTANT_REFINE");

    // Register scratch variables that are maintained by the
    // INSStaggeredHierarchyIntegrator.
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
            d_visit_writer->registerPlotQuantity(d_U_sc_var->getName(), "VECTOR", d_U_cc_current_idx, 0, d_U_scale);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(d_U_sc_var->getName()+stream.str(), "SCALAR", d_U_cc_current_idx, d, d_U_scale);
            }
        }

        if (d_output_P)
        {
            d_visit_writer->registerPlotQuantity(d_P_cc_var->getName(), "SCALAR", d_P_current_idx, 0, d_P_scale);
        }

        if (!d_F_fcn.isNull() && d_output_F)
        {
            d_visit_writer->registerPlotQuantity(d_F_sc_var->getName(), "VECTOR", d_F_cc_current_idx, 0, d_F_scale);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(d_F_sc_var->getName()+stream.str(), "SCALAR", d_F_cc_current_idx, d, d_F_scale);
            }
        }

        if (!d_Q_fcn.isNull() && d_output_Q)
        {
            d_visit_writer->registerPlotQuantity(d_Q_cc_var->getName(), "SCALAR", d_Q_current_idx, 0, d_Q_scale);
        }

        if (d_output_Omega)
        {
            d_visit_writer->registerPlotQuantity(d_Omega_var->getName(), (NDIM == 2) ? "SCALAR" : "VECTOR", d_Omega_current_idx);
#if (NDIM == 3)
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                d_visit_writer->registerPlotQuantity(d_Omega_var->getName()+stream.str(), "SCALAR", d_Omega_current_idx, d, d_Omega_scale);
            }
#endif
        }

        if (d_output_Div_U)
        {
            d_visit_writer->registerPlotQuantity(d_Div_U_var->getName(), "SCALAR", d_Div_U_current_idx, 0, d_Div_U_scale);
        }
    }

    // Set the current integration time.
    if (!RestartManager::getManager()->isFromRestart())
    {
        d_integrator_time = d_start_time;
        d_integrator_step = 0;
    }

    // Initialize operator and solver objects.
    initializeOperatorsAndSolvers();

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
}// initializeHierarchyIntegrator

void
INSStaggeredHierarchyIntegrator::preprocessIntegrateHierarchy(
    const double current_time,
    const double new_time,
    const int /*num_cycles*/)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;
    if (d_do_log) plog << d_object_name << "::integrateHierarchy_initialize(): current_time = " << current_time << ", new_time = " << new_time << ", dt = " << dt << "\n";

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
    d_U_rhs_vec ->allocateVectorData(current_time);  d_U_rhs_vec ->setToScalar(0.0);
    d_U_half_vec->allocateVectorData(current_time);  d_U_half_vec->setToScalar(0.0);
    d_N_vec     ->allocateVectorData(current_time);  d_N_vec     ->setToScalar(0.0);
    d_P_rhs_vec ->allocateVectorData(current_time);  d_P_rhs_vec ->setToScalar(0.0);

    // Initialize the right-hand side terms.
    const double rho    = d_problem_coefs.getRho();
    const double mu     = d_problem_coefs.getMu();
    const double lambda = d_problem_coefs.getLambda();
    PoissonSpecifications rhs_spec(d_object_name+"::rhs_spec");
    rhs_spec.setCConstant((rho/dt)-0.5*lambda);
    rhs_spec.setDConstant(        +0.5*mu    );
    const int U_rhs_idx = d_U_rhs_vec->getComponentDescriptorIndex(0);
    const Pointer<SideVariable<NDIM,double> > U_rhs_var = d_U_rhs_vec->getComponentVariable(0);
    d_hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
    d_hier_math_ops->laplace(U_rhs_idx, U_rhs_var, rhs_spec, d_U_scratch_idx, d_U_sc_var, d_U_bdry_bc_fill_op, current_time);

    // Set the initial guess.
    d_hier_sc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(0), d_U_current_idx);
    d_hier_cc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(1), d_P_current_idx);
    d_hier_sc_data_ops->copyData(d_U_new_idx, d_sol_vec->getComponentDescriptorIndex(0));
    d_hier_cc_data_ops->copyData(d_P_new_idx, d_sol_vec->getComponentDescriptorIndex(1));

    // Setup inhomogeneous boundary conditions.
    d_U_bc_helper->clearBcCoefData();
    d_U_bc_helper->cacheBcCoefData(d_U_scratch_idx, d_U_sc_var, d_U_bc_coefs, new_time, IntVector<NDIM>(SIDEG), d_hierarchy);
    d_stokes_op->setHomogeneousBc(false);
    d_stokes_op->modifyRhsForInhomogeneousBc(*d_rhs_vec);
    d_stokes_op->setHomogeneousBc(true);
    return;
}// preprocessIntegrateHierarchy

void
INSStaggeredHierarchyIntegrator::integrateHierarchy(
    const double current_time,
    const double new_time,
    const int /*cycle_num*/)
{
    const double dt  = new_time-current_time;
    const double rho = d_problem_coefs.getRho();

    // Perform a single step of fixed point iteration.

    // Compute U_half := 0.5*(u(n)+u(n+1)).
    const int U_half_idx = d_U_half_vec->getComponentDescriptorIndex(0);
    d_hier_sc_data_ops->linearSum(U_half_idx, 0.5, d_U_current_idx, 0.5, d_U_new_idx);

    // Setup the right-hand side vector.
    if (!d_creeping_flow && !d_convective_op.isNull())
    {
        d_convective_op->apply(*d_U_half_vec, *d_N_vec);
        const int N_idx = d_N_vec->getComponentDescriptorIndex(0);
        d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0), -rho, N_idx, d_rhs_vec->getComponentDescriptorIndex(0));
    }
    if (!d_F_fcn.isNull())
    {
        d_F_fcn->setDataOnPatchHierarchy(d_F_scratch_idx, d_F_sc_var, d_hierarchy, current_time+0.5*dt);
        d_hier_sc_data_ops->add(d_rhs_vec->getComponentDescriptorIndex(0), d_rhs_vec->getComponentDescriptorIndex(0), d_F_scratch_idx);
    }
    if (!d_Q_fcn.isNull())
    {
        d_Q_fcn->setDataOnPatchHierarchy(d_Q_current_idx, d_Q_cc_var, d_hierarchy, current_time);
        d_Q_fcn->setDataOnPatchHierarchy(d_Q_new_idx    , d_Q_cc_var, d_hierarchy, new_time    );
        d_hier_cc_data_ops->linearSum(d_Q_scratch_idx, 0.5, d_Q_current_idx, 0.5, d_Q_new_idx);
        d_Q_bdry_bc_fill_op->fillData(current_time+0.5*dt);
        if (!d_creeping_flow) computeDivSourceTerm(d_F_div_idx, d_Q_scratch_idx, U_half_idx);
        d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0), rho, d_F_div_idx, d_rhs_vec->getComponentDescriptorIndex(0));
        d_hier_cc_data_ops->subtract(d_rhs_vec->getComponentDescriptorIndex(1), d_rhs_vec->getComponentDescriptorIndex(1), d_Q_new_idx);
    }

    // Ensure there is no forcing at Dirichlet boundaries (the Dirichlet
    // boundary condition takes precedence).
    d_U_bc_helper->zeroValuesAtDirichletBoundaries(d_rhs_vec->getComponentDescriptorIndex(0));

    // Synchronize solution and right-hand-side data before solve.
    typedef SideDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent;
    SynchronizationTransactionComponent sol_synch_transaction = SynchronizationTransactionComponent(d_sol_vec->getComponentDescriptorIndex(0), "CONSERVATIVE_COARSEN");
    d_side_synch_op->resetTransactionComponent(sol_synch_transaction);
    d_side_synch_op->synchronizeData(current_time);
    SynchronizationTransactionComponent rhs_synch_transaction = SynchronizationTransactionComponent(d_rhs_vec->getComponentDescriptorIndex(0), "CONSERVATIVE_COARSEN");
    d_side_synch_op->resetTransactionComponent(rhs_synch_transaction);
    d_side_synch_op->synchronizeData(current_time);

    // Solve for u(n+1), p(n+1/2).
    d_stokes_solver->solveSystem(*d_sol_vec,*d_rhs_vec);

    // Synchronize solution data after solve.
    d_side_synch_op->resetTransactionComponent(sol_synch_transaction);
    d_side_synch_op->synchronizeData(current_time);

    // Enforce Dirichlet boundary conditions.
    d_U_bc_helper->resetValuesAtDirichletBoundaries(d_sol_vec->getComponentDescriptorIndex(0));

    // Pull out solution components.
    d_hier_sc_data_ops->copyData(d_U_new_idx, d_sol_vec->getComponentDescriptorIndex(0));
    d_hier_cc_data_ops->copyData(d_P_new_idx, d_sol_vec->getComponentDescriptorIndex(1));

    // Reset the right-hand side vector.
    if (!d_creeping_flow && !d_convective_op.isNull())
    {
        const int N_idx = d_N_vec->getComponentDescriptorIndex(0);
        d_hier_sc_data_ops->axpy(d_rhs_vec->getComponentDescriptorIndex(0), +rho, N_idx, d_rhs_vec->getComponentDescriptorIndex(0));
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
    return;
}// integrateHierarchy

void
INSStaggeredHierarchyIntegrator::postprocessIntegrateHierarchy(
    const double /*current_time*/,
    const double new_time,
    const bool skip_synchronize_new_state_data,
    const int /*num_cycles*/)
{
    // Synchronize new state data.
    if (!skip_synchronize_new_state_data)
    {
        if (d_do_log) plog << d_object_name << "::advanceHierarchy(): synchronizing updated data\n";
        synchronizeHierarchyData(d_new_context);
    }

    if (d_using_vorticity_tagging)
    {
        // Compute max ||Omega||_2.
        d_hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_new_idx);
        d_hier_math_ops->curl(d_Omega_scratch_idx, d_Omega_var, d_U_scratch_idx, d_U_sc_var, d_U_bdry_bc_fill_op, new_time);
#if (NDIM == 3)
        d_hier_math_ops->pointwiseL2Norm(d_Omega_Norm_scratch_idx, d_Omega_Norm_var, d_Omega_scratch_idx, d_Omega_var);
#endif

        const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
#if (NDIM == 2)
        d_Omega_max = d_hier_cc_data_ops->maxNorm(d_Omega_scratch_idx, wgt_cc_idx);
#endif
#if (NDIM == 3)
        d_Omega_max = d_hier_cc_data_ops->max(d_Omega_Norm_scratch_idx, wgt_cc_idx);
#endif
    }

    // Deallocate scratch data.
    d_U_rhs_vec ->deallocateVectorData();
    d_U_half_vec->deallocateVectorData();
    d_N_vec     ->deallocateVectorData();
    d_P_rhs_vec ->deallocateVectorData();
    return;
}// postprocessIntegrateHierarchy

void
INSStaggeredHierarchyIntegrator::regridHierarchy()
{
    const int coarsest_ln = 0;

    // Determine the divergence of the velocity field before regridding.
    d_hier_math_ops->div(d_Div_U_current_idx, d_Div_U_var, 1.0, d_U_current_idx, d_U_sc_var, d_no_fill_op, d_integrator_time, false, -1.0, d_Q_current_idx, d_Q_cc_var);
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const double Div_U_norm_1_pre  = d_hier_cc_data_ops->L1Norm( d_Div_U_current_idx, wgt_cc_idx);
    const double Div_U_norm_2_pre  = d_hier_cc_data_ops->L2Norm( d_Div_U_current_idx, wgt_cc_idx);
    const double Div_U_norm_oo_pre = d_hier_cc_data_ops->maxNorm(d_Div_U_current_idx, wgt_cc_idx);

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
                       << "  unrecognized regrid mode: " << enum_to_string<RegridMode>(d_regrid_mode) << "." << std::endl);
    }

    // Determine the divergence of the velocity field after regridding.
    d_hier_math_ops->div(d_Div_U_current_idx, d_Div_U_var, 1.0, d_U_current_idx, d_U_sc_var, d_no_fill_op, d_integrator_time, true, -1.0, d_Q_current_idx, d_Q_cc_var);
    const double Div_U_norm_1_post  = d_hier_cc_data_ops->L1Norm( d_Div_U_current_idx, wgt_cc_idx);
    const double Div_U_norm_2_post  = d_hier_cc_data_ops->L2Norm( d_Div_U_current_idx, wgt_cc_idx);
    const double Div_U_norm_oo_post = d_hier_cc_data_ops->maxNorm(d_Div_U_current_idx, wgt_cc_idx);

    // Project the interpolated velocity if needed.
    if (Div_U_norm_1_post  > d_regrid_max_div_growth_factor*Div_U_norm_1_pre ||
        Div_U_norm_2_post  > d_regrid_max_div_growth_factor*Div_U_norm_2_pre ||
        Div_U_norm_oo_post > d_regrid_max_div_growth_factor*Div_U_norm_oo_pre)
    {
        regridProjection();
    }

    // Synchronize the state data on the patch hierarchy.
    synchronizeHierarchyData(getCurrentContext());
    return;
}// regridHierarchy

/////////////////////////////// PROTECTED ////////////////////////////////////

double
INSStaggeredHierarchyIntegrator::getTimeStepSizeSpecialized()
{
    double dt = d_dt_max;
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        dt = std::min(dt, d_cfl_max*getStableTimestep(level));
    }
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    if (!initial_time && d_dt_growth_factor >= 1.0)
    {
        dt = std::min(dt,d_dt_growth_factor*d_dt_previous);
    }
    return dt;
}// getTimeStepSizeSpecialized

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

    // Use divergence- and curl-preserving prolongation to re-set the velocity
    // data.
    if (!initial_time && (level_number > 0 || !old_level.isNull()))
    {
        // Allocate scratch data.
        ComponentSelector scratch_data;
        scratch_data.setFlag( d_U_regrid_idx);
        scratch_data.setFlag(    d_U_src_idx);
        scratch_data.setFlag(d_indicator_idx);
        level->allocatePatchData(d_scratch_data, init_data_time);
        if (!old_level.isNull()) old_level->allocatePatchData(scratch_data, init_data_time);

        // Set the indicator data to equal "0" in each patch of the new
        // patch level and initialize values of U to cause a floating point
        // error if used incorrectly.
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
            // location in the new patch level which is a copy of a location
            // from the old patch level.
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
        fill_div_free_prolongation.registerRefine(d_U_current_idx, d_U_current_idx, d_U_regrid_idx, NULL);
        CartSideRobinPhysBdryOp phys_bdry_bc_op(d_U_regrid_idx, d_U_bc_coefs, false);
        CartSideDoubleDivPreservingRefine div_preserving_op(d_U_regrid_idx, d_U_src_idx, d_indicator_idx, init_data_time, Pointer<RefinePatchStrategy<NDIM> >(&phys_bdry_bc_op,false));
        fill_div_free_prolongation.createSchedule(level, old_level, level_number-1, hierarchy, &div_preserving_op)->fillData(init_data_time);

        // Free scratch data.
        level->deallocatePatchData(scratch_data);
        if (!old_level.isNull()) old_level->deallocatePatchData(scratch_data);
    }

    // Initialize level data at the initial time.
    if (initial_time)
    {
        // Initialize the maximum value of ||Omega||_2 on the grid.
        if (d_using_vorticity_tagging)
        {
            if (level_number == 0) d_Omega_max = 0.0;

            // Allocate scratch data.
            level->allocatePatchData(d_U_scratch_idx, init_data_time);

            // Fill ghost cells.
            HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
            Pointer<HierarchySideDataOpsReal<NDIM,double> > hier_sc_data_ops = hier_ops_manager->getOperationsDouble(d_U_sc_var, d_hierarchy, true);
            hier_sc_data_ops->resetLevels(0, level_number);
            hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
            typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
            InterpolationTransactionComponent U_bc_component(d_U_scratch_idx, SIDE_DATA_COARSEN_TYPE, d_bdry_extrap_type, CONSISTENT_TYPE_2_BDRY, d_U_bc_coefs);
            HierarchyGhostCellInterpolation U_bdry_bc_fill_op;
            U_bdry_bc_fill_op.initializeOperatorState(U_bc_component, d_hierarchy, 0, level_number);
            U_bdry_bc_fill_op.fillData(init_data_time);

            // Compute Omega = curl U on the grid and determine max ||Omega||_2.
            PatchCellDataOpsReal<NDIM,double> patch_cc_data_ops;
            PatchMathOps patch_math_ops;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<SideData<NDIM,double> > U_scratch_data = patch->getPatchData(d_U_scratch_idx);
                Pointer<CellData<NDIM,double> > Omega_current_data = patch->getPatchData(d_Omega_current_idx);
                patch_math_ops.curl(Omega_current_data, U_scratch_data, patch);
#if (NDIM == 2)
                d_Omega_max = std::max(d_Omega_max, +patch_cc_data_ops.max(Omega_current_data, patch_box));
                d_Omega_max = std::max(d_Omega_max, -patch_cc_data_ops.min(Omega_current_data, patch_box));
#endif
#if (NDIM == 3)
                Pointer<CellData<NDIM,double> > Omega_Norm_current_data = patch->getPatchData(d_Omega_Norm_current_idx);
                patch_math_ops.pointwiseL2Norm(Omega_Norm_current_data, Omega_current_data, patch);
                d_Omega_max = std::max(d_Omega_max, patch_cc_data_ops.max(Omega_Norm_current_data, patch_box));
#endif
            }

            // Deallocate scratch data.
            level->deallocatePatchData(d_U_scratch_idx);
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
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // Reset the hierarchy operations objects for the new hierarchy configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);

    d_hier_sc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_sc_data_ops->resetLevels(0, finest_hier_level);

    // Setup the patch boundary filling objects.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent U_bc_component(d_U_scratch_idx, SIDE_DATA_COARSEN_TYPE, d_bdry_extrap_type, CONSISTENT_TYPE_2_BDRY, d_U_bc_coefs);
    d_U_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_U_bdry_bc_fill_op->initializeOperatorState(U_bc_component, d_hierarchy);

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
    d_vectors_need_init = true;
    d_stokes_op_needs_init = true;
    d_convective_op_needs_init = true;
    d_helmholtz_solver_needs_init = true;
    d_poisson_solver_needs_init = true;
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

    // Untag all cells.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
        tags_data->fillAll(0);
    }

    // Tag cells based on the magnitude of the vorticity.
    //
    // Note that if either the relative or absolute threshold is zero for a
    // particular level, no tagging is performed on that level.
    if (d_using_vorticity_tagging)
    {
        const double Omega_rel_thresh =
            (level_number >= 0 && level_number < d_Omega_rel_thresh.getSize()
             ? d_Omega_rel_thresh[level_number]
             : (level_number < 0
                ? d_Omega_rel_thresh[0]
                : d_Omega_rel_thresh[d_Omega_rel_thresh.size()-1]));
        const double Omega_abs_thresh =
            (level_number >= 0 && level_number < d_Omega_abs_thresh.getSize()
             ? d_Omega_abs_thresh[level_number]
             : (level_number < 0
                ? d_Omega_abs_thresh[0]
                : d_Omega_abs_thresh[d_Omega_abs_thresh.size()-1]));
        if (Omega_rel_thresh > 0.0 && Omega_abs_thresh > 0.0)
        {
            const double thresh = sqrt(std::numeric_limits<double>::epsilon()) +
                std::min(Omega_rel_thresh*d_Omega_max, Omega_abs_thresh);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
                Pointer<CellData<NDIM,double> > Omega_current_data = patch->getPatchData(d_Omega_current_idx);
                for (CellIterator<NDIM> ic(patch_box); ic; ic++)
                {
                    const Index<NDIM>& i = ic();
#if (NDIM == 2)
                    if (std::abs((*Omega_current_data)(i)) > thresh)
                    {
                        (*tags_data)(i) = 1;
                    }
#endif
#if (NDIM == 3)
                    double norm_Omega_sq = 0.0;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        norm_Omega_sq += (*Omega_current_data)(i,d)*(*Omega_current_data)(i,d);
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
        const int U_sc_idx = var_db->mapVariableAndContextToIndex(d_U_sc_var, ctx);
        const int U_cc_idx = var_db->mapVariableAndContextToIndex(d_U_cc_var, ctx);
        d_hier_math_ops->interp(U_cc_idx, d_U_cc_var, U_sc_idx, d_U_sc_var, d_no_fill_op, d_integrator_time, synch_cf_interface);
    }

    // Interpolate f to cell centers.
    if (d_output_F && !d_F_fcn.isNull())
    {
        const int F_sc_idx = var_db->mapVariableAndContextToIndex(d_F_sc_var, ctx);
        const int F_cc_idx = var_db->mapVariableAndContextToIndex(d_F_cc_var, ctx);
        d_hier_math_ops->interp(F_cc_idx, d_F_cc_var, F_sc_idx, d_F_sc_var, d_no_fill_op, d_integrator_time, synch_cf_interface);
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
        d_hier_math_ops->curl(d_Omega_current_idx, d_Omega_var, d_U_scratch_idx, d_U_sc_var, d_U_bdry_bc_fill_op, d_integrator_time);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(d_U_scratch_idx);
        }
    }

    // Compute Div U.
    if (d_output_Div_U)
    {
        d_hier_math_ops->div(d_Div_U_current_idx, d_Div_U_var, 1.0, d_U_current_idx, d_U_sc_var, d_no_fill_op, d_integrator_time, false);
    }
    return;
}// setupPlotDataSpecialized

void
INSStaggeredHierarchyIntegrator::putToDatabaseSpecialized(
    Pointer<Database> db)
{
    db->putInteger("INS_STAGGERED_HIERARCHY_INTEGRATOR_VERSION",INS_STAGGERED_HIERARCHY_INTEGRATOR_VERSION);
    db->putDouble("d_rho",d_problem_coefs.getRho());
    db->putDouble("d_mu",d_problem_coefs.getMu());
    db->putDouble("d_lambda",d_problem_coefs.getLambda());
    db->putInteger("d_num_cycles",d_num_cycles);
    db->putDouble("d_cfl_max",d_cfl_max);
    db->putBool("d_using_vorticity_tagging",d_using_vorticity_tagging);
    db->putDoubleArray("d_Omega_rel_thresh",d_Omega_rel_thresh);
    db->putDoubleArray("d_Omega_abs_thresh",d_Omega_abs_thresh);
    db->putDouble("d_Omega_max",d_Omega_max);
    db->putBool("d_normalize_pressure",d_normalize_pressure);
    db->putString("d_convective_difference_form",enum_to_string<ConvectiveDifferencingType>(d_convective_difference_form));
    db->putBool("d_creeping_flow",d_creeping_flow);
    db->putDouble("d_regrid_max_div_growth_factor",d_regrid_max_div_growth_factor);
    db->putDouble("d_U_scale",d_U_scale);
    db->putDouble("d_P_scale",d_P_scale);
    db->putDouble("d_F_scale",d_F_scale);
    db->putDouble("d_Q_scale",d_Q_scale);
    db->putDouble("d_Omega_scale",d_Omega_scale);
    db->putDouble("d_Div_U_scale",d_Div_U_scale);
    db->putBool("d_output_U",d_output_U);
    db->putBool("d_output_P",d_output_P);
    db->putBool("d_output_F",d_output_F);
    db->putBool("d_output_Q",d_output_Q);
    db->putBool("d_output_Omega",d_output_Omega);
    db->putBool("d_output_Div_U",d_output_Div_U);
    return;
}// putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

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

            switch (d_convective_difference_form)
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
                               << "  unsupported differencing form: " << enum_to_string<ConvectiveDifferencingType>(d_convective_difference_form) << " \n"
                               << "  valid choices are: ADVECTIVE, CONSERVATIVE, SKEW_SYMMETRIC\n");
            }
        }
    }
    return;
}// computeDivSourceTerm

void
INSStaggeredHierarchyIntegrator::regridProjection()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const double volume = d_hier_math_ops->getVolumeOfPhysicalDomain();

    // Allocate temporary data.
    ComponentSelector scratch_idxs;
    scratch_idxs.setFlag(d_U_scratch_idx);
    scratch_idxs.setFlag(d_P_scratch_idx);
    scratch_idxs.setFlag(d_Div_U_scratch_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(scratch_idxs, d_integrator_time);
    }

    // Compute div U before applying the projection operator.
    const bool U_current_cf_bdry_synch = true;
    d_hier_math_ops->div(
        d_Div_U_scratch_idx, d_Div_U_var, // dst
        +1.0,                             // alpha
        d_U_current_idx, d_U_sc_var,      // src1
        d_no_fill_op,                     // src1_bdry_fill
        d_integrator_time,                // src1_bdry_fill_time
        U_current_cf_bdry_synch,          // src1_cf_bdry_synch
        -1.0,                             // beta
        d_Q_current_idx, d_Q_cc_var);     // src2

    // Write out the divergence of the velocity field to the log file.
    if (d_do_log)
    {
        const double Div_U_norm_1  = d_hier_cc_data_ops->L1Norm( d_Div_U_scratch_idx, wgt_cc_idx);
        const double Div_U_norm_2  = d_hier_cc_data_ops->L2Norm( d_Div_U_scratch_idx, wgt_cc_idx);
        const double Div_U_norm_oo = d_hier_cc_data_ops->maxNorm(d_Div_U_scratch_idx, wgt_cc_idx);
        if (d_Q_fcn.isNull())
        {
            plog << d_object_name << "::regridProjection():\n"
                 << "  performing regrid projection\n"
                 << "  before projection:\n"
                 << "    ||Div U||_1  = " << Div_U_norm_1  << "\n"
                 << "    ||Div U||_2  = " << Div_U_norm_2  << "\n"
                 << "    ||Div U||_oo = " << Div_U_norm_oo << "\n";
        }
        else
        {
            plog << d_object_name << "::regridProjection():\n"
                 << "  performing regrid projection\n"
                 << "  before projection:\n"
                 << "    ||Div U - Q||_1  = " << Div_U_norm_1  << "\n"
                 << "    ||Div U - Q||_2  = " << Div_U_norm_2  << "\n"
                 << "    ||Div U - Q||_oo = " << Div_U_norm_oo << "\n";
        }
    }

    // Setup the solver vectors.
    d_hier_cc_data_ops->setToScalar(d_P_scratch_idx, 0.0, false);
    d_hier_cc_data_ops->scale(d_Div_U_scratch_idx, -1.0, d_Div_U_scratch_idx);
    const double Div_U_mean = (1.0/volume)*d_hier_cc_data_ops->integral(d_Div_U_scratch_idx, wgt_cc_idx);
    d_hier_cc_data_ops->addScalar(d_Div_U_scratch_idx, d_Div_U_scratch_idx, -Div_U_mean);

    SAMRAIVectorReal<NDIM,double> sol_vec(d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec.addComponent(d_P_var, d_P_scratch_idx, wgt_cc_idx, d_hier_cc_data_ops);

    SAMRAIVectorReal<NDIM,double> rhs_vec(d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    rhs_vec.addComponent(d_Div_U_var, d_Div_U_scratch_idx, wgt_cc_idx, d_hier_cc_data_ops);

    // Setup the Poisson solver.
    const std::string regrid_projection_prefix = "regrid_projection_";
    LocationIndexRobinBcCoefs<NDIM> Phi_bc_coef;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        Phi_bc_coef.setBoundarySlope(2*d  ,0.0);
        Phi_bc_coef.setBoundarySlope(2*d+1,0.0);
    }

    PoissonSpecifications regrid_projection_spec(d_object_name+"::regrid_projection_spec");
    regrid_projection_spec.setCZero();
    regrid_projection_spec.setDConstant(-1.0);

    Pointer<CCLaplaceOperator> regrid_projection_op = new CCLaplaceOperator(d_object_name+"::Regrid Projection Poisson Operator", regrid_projection_spec, &Phi_bc_coef, true);
    regrid_projection_op->setHierarchyMathOps(d_hier_math_ops);
    regrid_projection_op->setPoissonSpecifications(regrid_projection_spec);
    regrid_projection_op->setPhysicalBcCoef(&Phi_bc_coef);
    regrid_projection_op->setHomogeneousBc(true);
    regrid_projection_op->setTime(d_integrator_time);
    regrid_projection_op->setHierarchyMathOps(d_hier_math_ops);

    if (d_regrid_projection_fac_pc_db.isNull())
    {
        TBOX_WARNING(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                     "  regrid projection poisson fac pc solver database is null." << std::endl);
    }
    Pointer<CCPoissonFACOperator> regrid_projection_fac_op = new CCPoissonFACOperator(d_object_name+"::Regrid Projection Poisson FAC Operator", d_regrid_projection_fac_pc_db);
    regrid_projection_fac_op->setPoissonSpecifications(regrid_projection_spec);
    regrid_projection_fac_op->setPoissonSpecifications(regrid_projection_spec);
    regrid_projection_fac_op->setPhysicalBcCoef(&Phi_bc_coef);
    regrid_projection_fac_op->setTime(d_integrator_time);
    Pointer<FACPreconditioner> regrid_projection_fac_pc = new FACPreconditioner(d_object_name+"::Regrid Projection Poisson Preconditioner", *regrid_projection_fac_op, d_regrid_projection_fac_pc_db);

    PETScKrylovLinearSolver regrid_projection_solver(d_object_name+"::Regrid Projection Poisson Krylov Solver", regrid_projection_prefix);
    regrid_projection_solver.setOperator(regrid_projection_op);
    regrid_projection_solver.setPreconditioner(regrid_projection_fac_pc);

    // Set some default options.
    regrid_projection_solver.setKSPType("gmres");
    regrid_projection_solver.setAbsoluteTolerance(1.0e-12);
    regrid_projection_solver.setRelativeTolerance(1.0e-08);
    regrid_projection_solver.setMaxIterations(25);
    regrid_projection_solver.setNullspace(true, NULL);
    regrid_projection_solver.setInitialGuessNonzero(false);

    // Solve the projection Poisson problem.
    regrid_projection_solver.solveSystem(sol_vec,rhs_vec);

    // Setup the interpolation transaction information.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent Phi_bc_component(d_P_scratch_idx, CELL_DATA_COARSEN_TYPE, d_bdry_extrap_type, CONSISTENT_TYPE_2_BDRY, &Phi_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> Phi_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    Phi_bdry_bc_fill_op->initializeOperatorState(Phi_bc_component, d_hierarchy);

    // Fill the physical boundary conditions for Phi.
    Phi_bdry_bc_fill_op->setHomogeneousBc(true);
    Phi_bdry_bc_fill_op->fillData(d_integrator_time);

    // Set U := U - grad Phi.
    const bool U_scratch_cf_bdry_synch = true;
    d_hier_math_ops->grad(
        d_U_scratch_idx, d_U_sc_var,  // dst
        U_scratch_cf_bdry_synch,      // dst_cf_bdry_synch
        1.0,                          // alpha
        d_P_scratch_idx, d_P_var,     // src
        d_no_fill_op,                 // src_bdry_fill
        d_integrator_time);           // src_bdry_fill_time
    d_hier_sc_data_ops->axpy(d_U_current_idx, -1.0, d_U_scratch_idx, d_U_current_idx);

    // Write out the divergence of the velocity field to the log file.
    if (d_do_log)
    {
        const bool U_current_cf_bdry_synch = true;
        d_hier_math_ops->div(
            d_Div_U_scratch_idx, d_Div_U_var, // dst
            +1.0,                             // alpha
            d_U_current_idx, d_U_sc_var,      // src1
            d_no_fill_op,                     // src1_bdry_fill
            d_integrator_time,                // src1_bdry_fill_time
            U_current_cf_bdry_synch,          // src1_cf_bdry_synch
            -1.0,                             // beta
            d_Q_current_idx, d_Q_cc_var);     // src2
        const double Div_U_norm_1  = d_hier_cc_data_ops->L1Norm( d_Div_U_scratch_idx, wgt_cc_idx);
        const double Div_U_norm_2  = d_hier_cc_data_ops->L2Norm( d_Div_U_scratch_idx, wgt_cc_idx);
        const double Div_U_norm_oo = d_hier_cc_data_ops->maxNorm(d_Div_U_scratch_idx, wgt_cc_idx);
        if (d_Q_fcn.isNull())
        {
            plog << "  after projection:\n"
                 << "    ||Div U||_1  = " << Div_U_norm_1  << "\n"
                 << "    ||Div U||_2  = " << Div_U_norm_2  << "\n"
                 << "    ||Div U||_oo = " << Div_U_norm_oo << "\n";
        }
        else
        {
            plog << "  after projection:\n"
                 << "    ||Div U - Q||_1  = " << Div_U_norm_1  << "\n"
                 << "    ||Div U - Q||_2  = " << Div_U_norm_2  << "\n"
                 << "    ||Div U - Q||_oo = " << Div_U_norm_oo << "\n";
        }
    }

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(scratch_idxs);
    }
    return;
}// regridProjection

void
INSStaggeredHierarchyIntegrator::initializeOperatorsAndSolvers()
{
    // Setup physical boundary conditions objects.
    d_U_bc_helper = new INSStaggeredPhysicalBoundaryHelper();
    blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM> bc_coefs = d_U_bc_coefs;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_U_bc_coefs[d] = new INSStaggeredVelocityBcCoef(d,d_problem_coefs,bc_coefs);
        d_U_star_bc_coefs[d] = new INSStaggeredIntermediateVelocityBcCoef(d,bc_coefs);
    }
    d_P_bc_coef = new INSStaggeredPressureBcCoef(d_problem_coefs,bc_coefs);
    d_Phi_bc_coef = new INSStaggeredProjectionBcCoef(bc_coefs);

    // Setup the Stokes operator.
    d_stokes_op_needs_init = true;
    d_stokes_op = new INSStaggeredStokesOperator(d_problem_coefs, d_U_bc_coefs, d_U_bc_helper, d_P_bc_coef, d_hier_math_ops);

    // Setup the convective operator.
    d_convective_op_needs_init = true;
    if (d_convective_op.isNull() && !d_creeping_flow)
    {
        d_convective_op = new INSStaggeredPPMConvectiveOperator(d_convective_difference_form);
    }

    // Setup the linear solver.
    const std::string stokes_prefix = "stokes_";
    d_stokes_solver_needs_init = true;
    d_stokes_solver = new PETScKrylovLinearSolver(d_object_name+"::stokes_solver", stokes_prefix);
    d_stokes_solver->setInitialGuessNonzero(true);
    d_stokes_solver->setOperator(d_stokes_op);
    d_stokes_solver->setKSPType("fgmres");

    // Setup the preconditioner and preconditioner sub-solvers.
    std::vector<std::string> pc_shell_types(4);
    pc_shell_types[0] = "projection";
    pc_shell_types[1] = "vanka";
    pc_shell_types[2] = "block_factorization";
    pc_shell_types[3] = "none";
    d_stokes_solver->setValidPCShellTypes(pc_shell_types);

    int ierr;
    PetscTruth flg;
    static const size_t len = 255;
    char stokes_pc_type_str[len];
    ierr = PetscOptionsGetString("stokes_", "-pc_type", stokes_pc_type_str, len, &flg);  IBTK_CHKERRQ(ierr);
    std::string stokes_pc_type = "shell";
    if (flg)
    {
        stokes_pc_type = std::string(stokes_pc_type_str);
    }

    if (!(stokes_pc_type == "none" || stokes_pc_type == "shell"))
    {
        TBOX_ERROR(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                   "  invalid stokes preconditioner type: " << stokes_pc_type << "\n"
                   "  valid stokes preconditioner types: shell, none" << std::endl);
    }

    char helmholtz_pc_type_str[len];
    ierr = PetscOptionsGetString("helmholtz_", "-pc_type", helmholtz_pc_type_str, len, &flg);  IBTK_CHKERRQ(ierr);
    std::string helmholtz_pc_type = "shell";
    if (flg)
    {
        helmholtz_pc_type = std::string(helmholtz_pc_type_str);
    }

    char poisson_pc_type_str[len];
    ierr = PetscOptionsGetString("poisson_", "-pc_type", poisson_pc_type_str, len, &flg);  IBTK_CHKERRQ(ierr);
    std::string poisson_pc_type = "shell";
    if (flg)
    {
        poisson_pc_type = std::string(poisson_pc_type_str);
    }

    std::string stokes_pc_shell_type = "none";
    if (stokes_pc_type == "shell")
    {
        char stokes_pc_shell_type_str[len];
        ierr = PetscOptionsGetString("stokes_", "-pc_shell_type", stokes_pc_shell_type_str, len, &flg);  IBTK_CHKERRQ(ierr);
        stokes_pc_shell_type = "projection";
        if (flg)
        {
            stokes_pc_shell_type = std::string(stokes_pc_shell_type_str);
        }

        if (!(stokes_pc_shell_type == "none" || stokes_pc_shell_type == "projection" || stokes_pc_shell_type == "vanka" || stokes_pc_shell_type == "block_factorization"))
        {
            TBOX_ERROR(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                       "  invalid stokes shell preconditioner type: " << stokes_pc_shell_type << "\n"
                       "  valid stokes shell preconditioner types: projection, vanka, block_factorization, none" << std::endl);
        }
    }

    // Setup the velocity subdomain solver.
    const bool needs_helmholtz_solver = stokes_pc_type == "shell" && (stokes_pc_shell_type == "projection" || stokes_pc_shell_type == "block_factorization");
    if (needs_helmholtz_solver)
    {
        const std::string helmholtz_prefix = "helmholtz_";

        // Setup the various solver components.
        d_helmholtz_spec = new PoissonSpecifications(d_object_name+"::helmholtz_spec");
        d_helmholtz_op = new SCLaplaceOperator(d_object_name+"::Helmholtz Operator", *d_helmholtz_spec, d_U_star_bc_coefs, true);
        d_helmholtz_op->setHierarchyMathOps(d_hier_math_ops);

        d_helmholtz_solver_needs_init = true;
        d_helmholtz_solver = new PETScKrylovLinearSolver(d_object_name+"::Helmholtz Krylov Solver", helmholtz_prefix);
        d_helmholtz_solver->setInitialGuessNonzero(false);
        d_helmholtz_solver->setOperator(d_helmholtz_op);

        if (helmholtz_pc_type != "none")
        {
            if (d_gridding_alg->getMaxLevels() == 1)
            {
                if (d_helmholtz_hypre_pc_db.isNull())
                {
                    TBOX_WARNING(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                                 "  helmholtz hypre pc solver database is null." << std::endl);
                }
                d_helmholtz_hypre_pc = new SCPoissonHypreLevelSolver(d_object_name+"::Helmholtz Preconditioner", d_helmholtz_hypre_pc_db);
                d_helmholtz_hypre_pc->setPoissonSpecifications(*d_helmholtz_spec);
                d_helmholtz_solver->setPreconditioner(d_helmholtz_hypre_pc);
            }
            else
            {
                if (d_helmholtz_fac_pc_db.isNull())
                {
                    TBOX_WARNING(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                                 "  helmholtz fac pc solver database is null." << std::endl);
                }
                d_helmholtz_fac_op = new SCPoissonFACOperator(d_object_name+"::Helmholtz FAC Operator", d_helmholtz_fac_pc_db);
                d_helmholtz_fac_op->setPoissonSpecifications(*d_helmholtz_spec);
                d_helmholtz_fac_pc = new IBTK::FACPreconditioner(d_object_name+"::Helmholtz Preconditioner", *d_helmholtz_fac_op, d_helmholtz_fac_pc_db);
                d_helmholtz_solver->setPreconditioner(d_helmholtz_fac_pc);
            }
        }

        // Set some default options.
        d_helmholtz_solver->setKSPType(d_gridding_alg->getMaxLevels() == 1 ? "preonly" : "gmres");
        d_helmholtz_solver->setAbsoluteTolerance(1.0e-30);
        d_helmholtz_solver->setRelativeTolerance(1.0e-02);
        d_helmholtz_solver->setMaxIterations(25);
    }
    else
    {
        d_helmholtz_spec = NULL;
        d_helmholtz_op = NULL;
        d_helmholtz_hypre_pc = NULL;
        d_helmholtz_fac_op = NULL;
        d_helmholtz_fac_pc = NULL;
        d_helmholtz_solver = NULL;
    }

    // Setup the pressure subdomain solver.
    const bool needs_poisson_solver = stokes_pc_type == "shell" && (stokes_pc_shell_type == "projection" || stokes_pc_shell_type == "block_factorization");
    if (needs_poisson_solver)
    {
        const std::string poisson_prefix = "poisson_";

        // Setup the various solver components.
        d_poisson_spec = new PoissonSpecifications(d_object_name+"::poisson_spec");
        d_poisson_op = new CCLaplaceOperator(d_object_name+"::Poisson Operator", *d_poisson_spec, d_Phi_bc_coef, true);
        d_poisson_op->setHierarchyMathOps(d_hier_math_ops);

        d_poisson_solver_needs_init = true;
        d_poisson_solver = new PETScKrylovLinearSolver(d_object_name+"::Poisson Krylov Solver", poisson_prefix);
        d_poisson_solver->setInitialGuessNonzero(false);
        d_poisson_solver->setOperator(d_poisson_op);

        if (poisson_pc_type != "none")
        {
            if (d_gridding_alg->getMaxLevels() == 1)
            {
                if (d_poisson_hypre_pc_db.isNull())
                {
                    TBOX_WARNING(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                                 "  Poisson hypre PC solver database is null." << std::endl);
                }
                d_poisson_hypre_pc = new CCPoissonHypreLevelSolver(d_object_name+"::Poisson Preconditioner", d_poisson_hypre_pc_db);
                d_poisson_hypre_pc->setPoissonSpecifications(*d_poisson_spec);
                d_poisson_solver->setPreconditioner(d_poisson_hypre_pc);
            }
            else
            {
                if (d_poisson_fac_pc_db.isNull())
                {
                    TBOX_WARNING(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                                 "  Poisson FAC PC solver database is null." << std::endl);
                }
                d_poisson_fac_op = new CCPoissonFACOperator(d_object_name+"::Poisson FAC Operator", d_poisson_fac_pc_db);
                d_poisson_fac_op->setPoissonSpecifications(*d_poisson_spec);
                d_poisson_fac_pc = new IBTK::FACPreconditioner(d_object_name+"::Poisson Preconditioner", *d_poisson_fac_op, d_poisson_fac_pc_db);
                d_poisson_solver->setPreconditioner(d_poisson_fac_pc);
            }
        }

        // Set some default options.
        d_poisson_solver->setKSPType(d_gridding_alg->getMaxLevels() == 1 ? "preonly" : "gmres");
        d_poisson_solver->setAbsoluteTolerance(1.0e-30);
        d_poisson_solver->setRelativeTolerance(1.0e-02);
        d_poisson_solver->setMaxIterations(25);
        if (d_normalize_pressure)
        {
            d_poisson_solver->setNullspace(d_normalize_pressure, NULL);
        }
    }
    else
    {
        d_poisson_spec = NULL;
        d_poisson_op = NULL;
        d_poisson_hypre_pc = NULL;
        d_poisson_fac_op = NULL;
        d_poisson_fac_pc = NULL;
        d_poisson_solver = NULL;
    }

    // Setup the Stokes preconditioner.
    if (stokes_pc_type == "shell")
    {
        if (stokes_pc_shell_type == "projection")
        {
            d_stokes_pc = new INSStaggeredProjectionPreconditioner(d_problem_coefs, d_Phi_bc_coef, d_normalize_pressure, d_helmholtz_solver, d_poisson_solver, d_hier_cc_data_ops, d_hier_sc_data_ops, d_hier_math_ops);
        }
        else if (stokes_pc_shell_type == "vanka")
        {
            if (d_vanka_fac_pc_db.isNull())
            {
                TBOX_WARNING(d_object_name << "::initializeHierarchyIntegrator():\n" <<
                             "  Vanka FAC PC solver database is null." << std::endl);
            }
            d_vanka_fac_op = new INSStaggeredBoxRelaxationFACOperator(d_object_name+"::Vanka FAC Operator", d_vanka_fac_pc_db);
            d_stokes_pc = new IBTK::FACPreconditioner(d_object_name+"::Vanka Preconditioner", *d_vanka_fac_op, d_vanka_fac_pc_db);
        }
        else if (stokes_pc_shell_type == "block_factorization")
        {
            d_stokes_pc = new INSStaggeredBlockFactorizationPreconditioner(d_problem_coefs, d_Phi_bc_coef, d_normalize_pressure, d_helmholtz_solver, d_poisson_solver, d_hier_cc_data_ops, d_hier_sc_data_ops, d_hier_math_ops);
        }
        d_stokes_solver->setPreconditioner(d_stokes_pc);
    }
    else
    {
        d_stokes_pc = NULL;
    }
    return;
}// initializeOperatorsAndSolvers

void
INSStaggeredHierarchyIntegrator::reinitializeOperatorsAndSolvers(
    const double current_time,
    const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const int wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
    const double rho    = d_problem_coefs.getRho();
    const double mu     = d_problem_coefs.getMu();
    const double lambda = d_problem_coefs.getLambda();

    if (d_vectors_need_init)
    {
        d_U_scratch_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::U_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_U_scratch_vec->addComponent(d_U_sc_var, d_U_scratch_idx, wgt_sc_idx, d_hier_sc_data_ops);

        d_P_scratch_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::P_scratch_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_P_scratch_vec->addComponent(d_P_cc_var, d_P_scratch_idx, wgt_cc_idx, d_hier_cc_data_ops);

        if (!d_U_rhs_vec .isNull()) d_U_rhs_vec ->freeVectorComponents();
        if (!d_U_half_vec.isNull()) d_U_half_vec->freeVectorComponents();
        if (!d_N_vec     .isNull()) d_N_vec     ->freeVectorComponents();
        if (!d_P_rhs_vec .isNull()) d_P_rhs_vec ->freeVectorComponents();

        d_U_rhs_vec  = d_U_scratch_vec->cloneVector(d_object_name+"::U_rhs_vec" );
        d_U_half_vec = d_U_scratch_vec->cloneVector(d_object_name+"::U_half_vec");
        d_N_vec      = d_U_scratch_vec->cloneVector(d_object_name+"::N_vec"     );
        d_P_rhs_vec  = d_P_scratch_vec->cloneVector(d_object_name+"::P_rhs_vec" );

        d_sol_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_sol_vec->addComponent(d_U_sc_var,d_U_scratch_idx,wgt_sc_idx,d_hier_sc_data_ops);
        d_sol_vec->addComponent(d_P_cc_var,d_P_scratch_idx,wgt_cc_idx,d_hier_cc_data_ops);

        d_rhs_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
        const int U_rhs_idx = d_U_rhs_vec->getComponentDescriptorIndex(0);
        d_rhs_vec->addComponent(d_U_sc_var,U_rhs_idx,wgt_sc_idx,d_hier_sc_data_ops);
        const int P_rhs_idx = d_P_rhs_vec->getComponentDescriptorIndex(0);
        d_rhs_vec->addComponent(d_P_cc_var,P_rhs_idx,wgt_cc_idx,d_hier_cc_data_ops);

        if (d_normalize_pressure)
        {
            if (!d_nul_vec.isNull()) d_nul_vec->freeVectorComponents();
            d_nul_vec = d_sol_vec->cloneVector(d_object_name+"::nul_vec");
        }

        d_vectors_need_init = false;
    }

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        INSStaggeredVelocityBcCoef* U_bc_coef = dynamic_cast<INSStaggeredVelocityBcCoef*>(d_U_bc_coefs[d]);
        U_bc_coef->setTimeInterval(current_time,new_time);
    }
    INSStaggeredPressureBcCoef* P_bc_coef = dynamic_cast<INSStaggeredPressureBcCoef*>(d_P_bc_coef);
    P_bc_coef->setTimeInterval(current_time,new_time);
    P_bc_coef->setVelocityCurrentPatchDataIndex(d_U_current_idx);
    P_bc_coef->setVelocityNewPatchDataIndex(d_U_new_idx);

    if (!d_helmholtz_solver.isNull())
    {
        d_helmholtz_spec->setCConstant((rho/dt)+0.5*lambda);
        d_helmholtz_spec->setDConstant(        -0.5*mu    );

        d_helmholtz_op->setPoissonSpecifications(*d_helmholtz_spec);
        d_helmholtz_op->setPhysicalBcCoefs(d_U_star_bc_coefs);
        d_helmholtz_op->setHomogeneousBc(true);
        d_helmholtz_op->setTime(new_time);
        d_helmholtz_op->setHierarchyMathOps(d_hier_math_ops);

        if (!d_helmholtz_hypre_pc.isNull())
        {
            d_helmholtz_hypre_pc->setPoissonSpecifications(*d_helmholtz_spec);
            d_helmholtz_hypre_pc->setPhysicalBcCoefs(d_U_star_bc_coefs);
            d_helmholtz_hypre_pc->setHomogeneousBc(true);
            d_helmholtz_hypre_pc->setTime(new_time);
        }
        else if (!d_helmholtz_fac_op.isNull())
        {
            d_helmholtz_fac_op->setPoissonSpecifications(*d_helmholtz_spec);
            d_helmholtz_fac_op->setPhysicalBcCoefs(d_U_star_bc_coefs);
            d_helmholtz_fac_op->setTime(new_time);
        }

        d_helmholtz_solver->setInitialGuessNonzero(false);
        d_helmholtz_solver->setOperator(d_helmholtz_op);
        if (d_helmholtz_solver_needs_init || !MathUtilities<double>::equalEps(dt,d_dt_previous))
        {
            if (d_do_log) plog << d_object_name << "::integrateHierarchy(): Initializing Helmholtz solver" << std::endl;
            d_helmholtz_solver->initializeSolverState(*d_U_scratch_vec,*d_U_rhs_vec);
        }
        d_helmholtz_solver_needs_init = false;
    }

    if (!d_poisson_solver.isNull())
    {
        d_poisson_spec->setCZero();
        d_poisson_spec->setDConstant(-1.0);

        d_poisson_op->setPoissonSpecifications(*d_poisson_spec);
        d_poisson_op->setPhysicalBcCoef(d_Phi_bc_coef);
        d_poisson_op->setHomogeneousBc(true);
        d_poisson_op->setTime(current_time+0.5*dt);
        d_poisson_op->setHierarchyMathOps(d_hier_math_ops);

        if (!d_poisson_hypre_pc.isNull())
        {
            d_poisson_hypre_pc->setPoissonSpecifications(*d_poisson_spec);
            d_poisson_hypre_pc->setPhysicalBcCoef(d_Phi_bc_coef);
            d_poisson_hypre_pc->setHomogeneousBc(true);
            d_poisson_hypre_pc->setTime(current_time+0.5*dt);
        }
        else if (!d_poisson_fac_op.isNull())
        {
            d_poisson_fac_op->setPoissonSpecifications(*d_poisson_spec);
            d_poisson_fac_op->setPhysicalBcCoef(d_Phi_bc_coef);
            d_poisson_fac_op->setTime(current_time+0.5*dt);
        }

        d_poisson_solver->setInitialGuessNonzero(false);
        d_poisson_solver->setOperator(d_poisson_op);
        if (d_poisson_solver_needs_init)
        {
            if (d_do_log) plog << d_object_name << "::integrateHierarchy(): Initializing Poisson solver" << std::endl;
            d_poisson_solver->initializeSolverState(*d_P_scratch_vec,*d_P_rhs_vec);
        }
        d_poisson_solver_needs_init = false;
    }

    if (!d_stokes_op.isNull())
    {
        if (d_stokes_op_needs_init && !d_stokes_solver_needs_init)
        {
            if (d_do_log) plog << d_object_name << "::integrateHierarchy(): Initializing incompressible Stokes operator" << std::endl;
            d_stokes_op->initializeOperatorState(*d_U_scratch_vec,*d_U_rhs_vec);
        }
        d_stokes_op_needs_init = false;
    }

    if (!d_stokes_solver.isNull())
    {
        d_stokes_op->setTimeInterval(current_time,new_time);
        d_stokes_solver->setOperator(d_stokes_op);
        if (d_stokes_solver_needs_init)
        {
            if (d_do_log) plog << d_object_name << "::integrateHierarchy(): Initializing incompressible Stokes solver" << std::endl;
            if (!d_vanka_fac_op.isNull())
            {
                d_vanka_fac_op->setProblemCoefficients(d_problem_coefs,dt);
                d_vanka_fac_op->setTimeInterval(current_time,new_time);
                d_vanka_fac_op->setPhysicalBcCoefs(d_U_star_bc_coefs,d_Phi_bc_coef);
            }
            if (!d_stokes_pc.isNull())
            {
                d_stokes_pc->setTimeInterval(current_time,new_time);
            }
            d_stokes_solver->initializeSolverState(*d_sol_vec,*d_rhs_vec);
            if (d_normalize_pressure)
            {
                d_nul_vec->allocateVectorData(current_time);
                d_hier_sc_data_ops->setToScalar(d_nul_vec->getComponentDescriptorIndex(0), 0.0);
                d_hier_cc_data_ops->setToScalar(d_nul_vec->getComponentDescriptorIndex(1), 1.0);
                d_stokes_solver->setNullspace(false, d_nul_vec);
            }
        }
        d_stokes_solver_needs_init = false;
    }

    if (!d_convective_op.isNull())
    {
        if (d_convective_op_needs_init)
        {
            if (d_do_log) plog << d_object_name << "::integrateHierarchy(): Initializing convective operator" << std::endl;
            d_convective_op->initializeOperatorState(*d_U_scratch_vec,*d_U_rhs_vec);
        }
        d_convective_op_needs_init = false;
    }
    return;
}// reinitializeOperatorsAndSolvers

double
INSStaggeredHierarchyIntegrator::getStableTimestep(
    Pointer<PatchLevel<NDIM> > level) const
{
    double stable_dt = std::numeric_limits<double>::max();
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        stable_dt = std::min(stable_dt,getStableTimestep(patch));
    }
    stable_dt = SAMRAI_MPI::minReduction(stable_dt);
    return stable_dt;
}// getStableTimestep

double
INSStaggeredHierarchyIntegrator::getStableTimestep(
    Pointer<Patch<NDIM> > patch) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const Index<NDIM>& ilower = patch->getBox().lower();
    const Index<NDIM>& iupper = patch->getBox().upper();

    Pointer<SideData<NDIM,double> > U_data = patch->getPatchData(d_U_sc_var, getCurrentContext());
    const IntVector<NDIM>& U_ghost_cells = U_data->getGhostCellWidth();

    double stable_dt = std::numeric_limits<double>::max();
    NAVIER_STOKES_SC_STABLEDT_FC(
#if (NDIM == 2)
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),
        U_ghost_cells(0),U_ghost_cells(1),
        U_data->getPointer(0),U_data->getPointer(1),
        stable_dt
#endif
#if (NDIM == 3)
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
        U_ghost_cells(0),U_ghost_cells(1),U_ghost_cells(2),
        U_data->getPointer(0),U_data->getPointer(1),U_data->getPointer(2),
        stable_dt
#endif
                                 );
    return stable_dt;
}// getStableTimestep

void
INSStaggeredHierarchyIntegrator::getFromInput(
    Pointer<Database> db,
    const bool is_from_restart)
{
    if (!is_from_restart)
    {
        if (db->keyExists("rho"))
        {
            d_problem_coefs.setRho(db->getDouble("rho"));
        }
        else
        {
            TBOX_ERROR(d_object_name << ":  "
                       << "Key data `rho' not found in input.");
        }

        if (db->keyExists("mu"))
        {
            d_problem_coefs.setMu(db->getDouble("mu"));
        }
        else
        {
            TBOX_ERROR(d_object_name << ":  "
                       << "Key data `mu' not found in input.");
        }

        if (db->keyExists("lambda"))
        {
            d_problem_coefs.setLambda(db->getDouble("lambda"));
        }
        else
        {
            d_problem_coefs.setLambda(0.0);
        }
    }
    if (db->keyExists("num_cycles")) d_num_cycles = db->getInteger("num_cycles");
    if (db->keyExists("cfl")) d_cfl_max = db->getDouble("cfl");
    else if (db->keyExists("cfl_max")) d_cfl_max = db->getDouble("cfl_max");
    else if (db->keyExists("CFL")) d_cfl_max = db->getDouble("CFL");
    else if (db->keyExists("CFL_max")) d_cfl_max = db->getDouble("CFL_max");
    if (db->keyExists("using_vorticity_tagging")) d_using_vorticity_tagging = db->getBool("using_vorticity_tagging");
    if (db->keyExists("Omega_rel_thresh")) db->getDoubleArray("Omega_rel_thresh");
    if (db->keyExists("Omega_abs_thresh")) db->getDoubleArray("Omega_abs_thresh");
    if (db->keyExists("normalize_pressure")) d_normalize_pressure = db->getBool("normalize_pressure");
    if (db->keyExists("convective_difference_form")) d_convective_difference_form = string_to_enum<ConvectiveDifferencingType>(db->getString("convective_difference_form"));
    if (db->keyExists("creeping_flow")) d_creeping_flow = db->getBool("creeping_flow");
    if (db->keyExists("regrid_max_div_growth_factor")) d_regrid_max_div_growth_factor = db->getDouble("regrid_max_div_growth_factor");
    if (db->keyExists("U_scale")) d_U_scale = db->getDouble("U_scale");
    if (db->keyExists("P_scale")) d_P_scale = db->getDouble("P_scale");
    if (db->keyExists("F_scale")) d_F_scale = db->getDouble("F_scale");
    if (db->keyExists("Q_scale")) d_Q_scale = db->getDouble("Q_scale");
    if (db->keyExists("Omega_scale")) d_Omega_scale = db->getDouble("Omega_scale");
    if (db->keyExists("Div_U_scale")) d_Div_U_scale = db->getDouble("Div_U_scale");
    if (db->keyExists("output_U")) d_output_U = db->getBool("output_U");
    if (db->keyExists("output_P")) d_output_P = db->getBool("output_P");
    if (db->keyExists("output_F")) d_output_F = db->getBool("output_F");
    if (db->keyExists("output_Q")) d_output_Q = db->getBool("output_Q");
    if (db->keyExists("output_Omega")) d_output_Omega = db->getBool("output_Omega");
    if (db->keyExists("output_Div_U")) d_output_Div_U = db->getBool("output_Div_U");
    if (db->keyExists("HelmholtzHypreSolver")) d_helmholtz_hypre_pc_db = db->getDatabase("HelmholtzHypreSolver");
    else if (db->keyExists("HelmholtzHyprePreconditioner")) d_helmholtz_hypre_pc_db = db->getDatabase("HelmholtzHyprePreconditioner");
    if (db->keyExists("HelmholtzFACSolver")) d_helmholtz_fac_pc_db = db->getDatabase("HelmholtzFACSolver");
    else if (db->keyExists("HelmholtzFACPreconditioner")) d_helmholtz_fac_pc_db = db->getDatabase("HelmholtzFACPreconditioner");
    if (db->keyExists("PoissonHypreSolver")) d_poisson_hypre_pc_db = db->getDatabase("PoissonHypreSolver");
    else if (db->keyExists("PoissonHyprePreconditioner")) d_poisson_hypre_pc_db = db->getDatabase("PoissonHyprePreconditioner");
    if (db->keyExists("PoissonFACSolver")) d_poisson_fac_pc_db = db->getDatabase("PoissonFACSolver");
    else if (db->keyExists("PoissonFACPreconditioner")) d_poisson_fac_pc_db = db->getDatabase("PoissonFACPreconditioner");
    if (db->keyExists("VankaFACSolver")) d_vanka_fac_pc_db = db->getDatabase("VankaFACSolver");
    else if (db->keyExists("VankaFACPreconditioner")) d_vanka_fac_pc_db = db->getDatabase("VankaFACPreconditioner");
    if (db->keyExists("RegridProjectionFACSolver")) d_regrid_projection_fac_pc_db = db->getDatabase("RegridProjectionFACSolver");
    else if (db->keyExists("RegridProjectionFACPreconditioner")) d_regrid_projection_fac_pc_db = db->getDatabase("RegridProjectionFACPreconditioner");
    return;
}// getFromInput

void
INSStaggeredHierarchyIntegrator::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to "
                   << d_object_name << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("INS_STAGGERED_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != INS_STAGGERED_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    d_problem_coefs.setRho(db->getDouble("d_rho"));
    d_problem_coefs.setMu(db->getDouble("d_mu"));
    d_problem_coefs.setLambda(db->getDouble("d_lambda"));
    d_num_cycles = db->getInteger("d_num_cycles");
    d_cfl_max = db->getDouble("d_cfl_max");
    d_using_vorticity_tagging = db->getBool("d_using_vorticity_tagging");
    d_Omega_rel_thresh = db->getDoubleArray("d_Omega_rel_thresh");
    d_Omega_abs_thresh = db->getDoubleArray("d_Omega_abs_thresh");
    d_Omega_max = db->getDouble("d_Omega_max");
    d_normalize_pressure = db->getBool("d_normalize_pressure");
    d_convective_difference_form = string_to_enum<ConvectiveDifferencingType>(db->getString("d_convective_difference_form"));
    d_creeping_flow = db->getBool("d_creeping_flow");
    d_regrid_max_div_growth_factor = db->getDouble("d_regrid_max_div_growth_factor");
    d_U_scale = db->getDouble("d_U_scale");
    d_P_scale = db->getDouble("d_P_scale");
    d_F_scale = db->getDouble("d_F_scale");
    d_Q_scale = db->getDouble("d_Q_scale");
    d_Omega_scale = db->getDouble("d_Omega_scale");
    d_Div_U_scale = db->getDouble("d_Div_U_scale");
    d_output_U = db->getBool("d_output_U");
    d_output_P = db->getBool("d_output_P");
    d_output_F = db->getBool("d_output_F");
    d_output_Q = db->getBool("d_output_Q");
    d_output_Omega = db->getBool("d_output_Omega");
    d_output_Div_U = db->getBool("d_output_Div_U");
    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::INSStaggeredHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
