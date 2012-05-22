// Filename: AdvDiffHierarchyIntegrator.C
// Created on 21 May 2012 by Boyce Griffith
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

#include "AdvDiffHierarchyIntegrator.h"

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
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>
#include <tbox/NullDatabase.h>

// C++ STDLIB INCLUDES
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Version of AdvDiffHierarchyIntegrator restart file data.
static const int ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION = 2;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffHierarchyIntegrator::~AdvDiffHierarchyIntegrator()
{
    // Deallocate all solver components.
    //
    // NOTE: The following code ensures that the solver components are
    // deallocated in the correct order.
    d_helmholtz_solvers.clear();
    d_helmholtz_fac_pcs.clear();
    d_helmholtz_fac_ops.clear();
    d_helmholtz_ops.clear();
    d_helmholtz_specs.clear();
    return;
}// ~AdvDiffHierarchyIntegrator

ViscousTimesteppingType
AdvDiffHierarchyIntegrator::getViscousTimesteppingType() const
{
    return d_viscous_timestepping_type;
}// getViscousTimesteppingType

void
AdvDiffHierarchyIntegrator::registerAdvectionVelocity(
    Pointer<FaceVariable<NDIM,double> > u_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!u_var.isNull());
#endif
    d_u_var.insert(u_var);
    d_u_is_div_free[u_var] = true;
    return;
}// registerAdvectionVelocity

void
AdvDiffHierarchyIntegrator::setAdvectionVelocityIsDivergenceFree(
    Pointer<FaceVariable<NDIM,double> > u_var,
    const bool is_div_free)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_u_var.find(u_var) != d_u_var.end());
#endif
    d_u_is_div_free[u_var] = is_div_free;
    return;
}// setAdvectionVelocityIsDivergenceFree

void
AdvDiffHierarchyIntegrator::setAdvectionVelocityFunction(
    Pointer<FaceVariable<NDIM,double> > u_var,
    Pointer<IBTK::CartGridFunction> u_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_u_var.find(u_var) != d_u_var.end());
#endif
    d_u_fcn[u_var] = u_fcn;
    return;
}// setAdvectionVelocityFunction

void
AdvDiffHierarchyIntegrator::registerSourceTerm(
    Pointer<CellVariable<NDIM,double> > F_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!F_var.isNull());
#endif
    d_F_var.insert(F_var);
    return;
}// registerSourceTerm

void
AdvDiffHierarchyIntegrator::setSourceTermFunction(
    Pointer<CellVariable<NDIM,double> > F_var,
    Pointer<IBTK::CartGridFunction> F_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_F_var.find(F_var) != d_F_var.end());
#endif
    d_F_fcn[F_var] = F_fcn;
    return;
}// setSourceTermFunction

void
AdvDiffHierarchyIntegrator::registerTransportedQuantity(
    Pointer<CellVariable<NDIM,double> > Q_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!Q_var.isNull());
#endif
    d_Q_var.insert(Q_var);
    d_Q_difference_form[Q_var] = CONSERVATIVE;
    d_Q_diffusion_coef[Q_var] = 0.0;
    d_Q_damping_coef[Q_var] = 0.0;

    Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
    const int Q_depth = Q_factory->getDefaultDepth();
    d_Q_bc_coef[Q_var] = std::vector<RobinBcCoefStrategy<NDIM>*>(Q_depth,static_cast<RobinBcCoefStrategy<NDIM>*>(NULL));

    Pointer<CellVariable<NDIM,double> > Psi_var = new CellVariable<NDIM,double>(Q_var->getName()+"::Psi",Q_depth);
    d_Psi_var.insert(Psi_var);
    d_Q_Psi_map[Q_var] = Psi_var;
    return;
}// registerTransportedQuantity

void
AdvDiffHierarchyIntegrator::setAdvectionVelocity(
    Pointer<CellVariable<NDIM,double> > Q_var,
    Pointer<FaceVariable<NDIM,double> > u_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    TBOX_ASSERT(d_u_var.find(u_var) != d_u_var.end());
#endif
    d_Q_u_map[Q_var] = u_var;
    return;
}// setAdvectionVelocity

void
AdvDiffHierarchyIntegrator::setSourceTerm(
    Pointer<CellVariable<NDIM,double> > Q_var,
    Pointer<CellVariable<NDIM,double> > F_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    TBOX_ASSERT(d_F_var.find(F_var) != d_F_var.end());
#endif
    d_Q_F_map[Q_var] = F_var;
    return;
}// setSourceTerm

void
AdvDiffHierarchyIntegrator::setConvectiveDifferencingType(
    Pointer<CellVariable<NDIM,double> > Q_var,
    const ConvectiveDifferencingType difference_form)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_difference_form[Q_var] = difference_form;
    return;
}// setConvectiveDifferencingType

void
AdvDiffHierarchyIntegrator::setDiffusionCoefficient(
    Pointer<CellVariable<NDIM,double> > Q_var,
    const double kappa)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_diffusion_coef[Q_var] = kappa;
    return;
}// setDiffusionCoefficient

void
AdvDiffHierarchyIntegrator::setDampingCoefficient(
    Pointer<CellVariable<NDIM,double> > Q_var,
    const double lambda)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_damping_coef[Q_var] = lambda;
    return;
}// setDampingCoefficient

void
AdvDiffHierarchyIntegrator::setInitialConditions(
    Pointer<CellVariable<NDIM,double> > Q_var,
    Pointer<IBTK::CartGridFunction> Q_init)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_init[Q_var] = Q_init;
    return;
}// setInitialConditions

void
AdvDiffHierarchyIntegrator::setPhysicalBcCoefs(
    Pointer<CellVariable<NDIM,double> > Q_var,
    RobinBcCoefStrategy<NDIM>* Q_bc_coef)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
    const unsigned int Q_depth = Q_factory->getDefaultDepth();
    TBOX_ASSERT(Q_depth == 1);
#endif
    d_Q_bc_coef[Q_var] = std::vector<RobinBcCoefStrategy<NDIM>*>(1,Q_bc_coef);
    return;
}// setPhysicalBcCoefs

void
AdvDiffHierarchyIntegrator::setPhysicalBcCoefs(
    Pointer<CellVariable<NDIM,double> > Q_var,
    std::vector<RobinBcCoefStrategy<NDIM>*> Q_bc_coef)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
    const unsigned int Q_depth = Q_factory->getDefaultDepth();
    TBOX_ASSERT(Q_depth == Q_bc_coef.size());
#endif
    d_Q_bc_coef[Q_var] = Q_bc_coef;
    return;
}// setPhysicalBcCoefs

/////////////////////////////// PROTECTED ////////////////////////////////////

AdvDiffHierarchyIntegrator::AdvDiffHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    bool register_for_restart)
    : HierarchyIntegrator(object_name, input_db, register_for_restart),
      d_integrator_is_initialized(false),
      d_viscous_timestepping_type(CRANK_NICOLSON),
      d_u_var(),
      d_u_is_div_free(),
      d_u_fcn(),
      d_F_var(),
      d_F_fcn(),
      d_Q_var(),
      d_Psi_var(),
      d_Q_u_map(),
      d_Q_F_map(),
      d_Q_Psi_map(),
      d_Q_difference_form(),
      d_Q_diffusion_coef(),
      d_Q_damping_coef(),
      d_Q_init(),
      d_Q_bc_coef(),
      d_hier_cc_data_ops(NULL),
      d_sol_vecs(),
      d_rhs_vecs(),
      d_max_iterations(25),
      d_abs_residual_tol(1.0e-30),
      d_rel_residual_tol(1.0e-8),
      d_using_FAC(true),
      d_helmholtz_ops(),
      d_helmholtz_specs(),
      d_helmholtz_solvers(),
      d_helmholtz_fac_ops(),
      d_helmholtz_fac_pcs(),
      d_helmholtz_solvers_need_init(),
      d_coarsest_reset_ln(-1),
      d_finest_reset_ln(-1),
      d_fac_op_db(NULL),
      d_fac_pc_db(NULL)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
#endif
    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);

    // Get initialization data for the FAC ops and FAC preconditioners.
    if (d_using_FAC)
    {
        if (input_db->keyExists("FACOp"))
        {
            d_fac_op_db = input_db->getDatabase("FACOp");
        }
        else if (input_db->keyExists("FACOps"))
        {
            d_fac_op_db = input_db->getDatabase("FACOps");
        }
        else
        {
            d_fac_op_db = new NullDatabase();
        }

        if (input_db->keyExists("FACPreconditioner"))
        {
            d_fac_pc_db = input_db->getDatabase("FACPreconditioner");
        }
        else if (input_db->keyExists("FACPreconditioners"))
        {
            d_fac_pc_db = input_db->getDatabase("FACPreconditioners");
        }
        else
        {
            d_fac_pc_db = new NullDatabase();
        }
    }
    return;
}// AdvDiffHierarchyIntegrator

void
AdvDiffHierarchyIntegrator::putToDatabaseSpecialized(
    Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION", ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION);
    db->putString("d_viscous_timestepping_type", enum_to_string<ViscousTimesteppingType>(d_viscous_timestepping_type));
    return;
}// putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
AdvDiffHierarchyIntegrator::getFromInput(
    Pointer<Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    // Read in data members from input database.
    if (!is_from_restart)
    {
        if (db->keyExists("viscous_timestepping_type")) d_viscous_timestepping_type = string_to_enum<ViscousTimesteppingType>(db->getString("viscous_timestepping_type"));
    }
    if (db->keyExists("max_iterations")) d_max_iterations = db->getInteger("max_iterations");
    if (db->keyExists("abs_residual_tol")) d_abs_residual_tol = db->getDouble("abs_residual_tol");
    if (db->keyExists("rel_residual_tol")) d_rel_residual_tol = db->getDouble("rel_residual_tol");
    if (db->keyExists("using_FAC")) d_using_FAC = db->getBool("using_FAC");
    return;
}// getFromInput

void
AdvDiffHierarchyIntegrator::getFromRestart()
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
    int ver = db->getInteger("ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    d_viscous_timestepping_type = string_to_enum<ViscousTimesteppingType>(db->getString("d_viscous_timestepping_type"));
    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
