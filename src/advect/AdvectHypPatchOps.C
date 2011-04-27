// Filename: AdvectHypPatchOps.C
// Created on 12 Mar 2004 by Boyce Griffith
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

#include "AdvectHypPatchOps.h"

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
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/PhysicalBoundaryUtilities.h>

// SAMRAI INCLUDES
#include <ArrayDataBasicOps.h>
#include <BoundaryBox.h>
#include <BoxArray.h>
#include <CartesianPatchGeometry.h>
#include <CellDataFactory.h>
#include <CellIndex.h>
#include <CellIterator.h>
#include <CellVariable.h>
#include <FaceIndex.h>
#include <FaceVariable.h>
#include <Index.h>
#include <LoadBalancer.h>
#include <PatchCellDataOpsReal.h>
#include <PatchFaceDataOpsReal.h>
#include <PatchData.h>
#include <VariableDatabase.h>
#include <tbox/Array.h>
#include <tbox/RestartManager.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <limits>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define ADVECT_CONSDIFF_FC FC_FUNC_(advect_consdiff2d, ADVECT_CONSDIFF2D)
#define ADVECT_CONSDIFFWITHDIVSOURCE_FC FC_FUNC_(advect_consdiffwithdivsource2d, ADVECT_CONSDIFFWITHDIVSOURCE2D)
#define ADVECT_DETECTGRAD_FC FC_FUNC_(advect_detectgrad2d, ADVECT_DETECTGRAD2D)
#endif

#if (NDIM == 3)
#define ADVECT_CONSDIFF_FC FC_FUNC_(advect_consdiff3d, ADVECT_CONSDIFF3D)
#define ADVECT_CONSDIFFWITHDIVSOURCE_FC FC_FUNC_(advect_consdiffwithdivsource3d, ADVECT_CONSDIFFWITHDIVSOURCE3D)
#define ADVECT_DETECTGRAD_FC FC_FUNC_(advect_detectgrad3d, ADVECT_DETECTGRAD3D)
#endif

extern "C"
{
    void
    ADVECT_CONSDIFF_FC(
        const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
#endif
        double*
                       );

    void
    ADVECT_CONSDIFFWITHDIVSOURCE_FC(
        const double* , const double& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const double* , const double* ,
        const double* , const double* ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
#endif
        double*);

    void
    ADVECT_DETECTGRAD_FC(
#if (NDIM == 2)
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& , const int&,
        const int& , const int& , const int&,
#endif
#if (NDIM == 3)
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& , const int&,
        const int& , const int& , const int&,
        const int& , const int& , const int&,
#endif
        const double* ,
        const double& ,
        const int&,
        const double*,
        int* , int*);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 4;
static const int FLUXG = 1;

// Values for cell tagging routines.
static const int TRUE_VAL  = 1;
static const int FALSE_VAL = 0;

// Version of AdvectHypPatchOps restart file data.
static const int ADVECT_HYP_PATCH_OPS_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvectHypPatchOps::AdvectHypPatchOps(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<GodunovAdvector> godunov_advector,
    Pointer<CartesianGridGeometry<NDIM> > grid_geom,
    bool register_for_restart)
    : d_integrator(NULL),
      d_godunov_advector(godunov_advector),
      d_u_var(),
      d_manage_u_data(),
      d_u_is_div_free(),
      d_u_fcn(),
      d_compute_init_velocity(true),
      d_compute_half_velocity(true),
      d_compute_final_velocity(true),
      d_F_var(),
      d_manage_F_data(),
      d_F_fcn(),
      d_grad_Phi_var(),
      d_manage_grad_Phi_data(),
      d_grad_Phi_fcn(),
      d_Q_var(),
      d_manage_Q_data(),
      d_Q_u_map(),
      d_Q_grad_Phi_map(),
      d_Q_F_map(),
      d_Q_difference_form(),
      d_Q_init(),
      d_Q_bc_coef(),
      d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_grid_geometry(grid_geom),
      d_visit_writer(NULL),
      d_extrap_bc_helper(),
      d_ghosts(CELLG),
      d_flux_ghosts(FLUXG),
      d_extrap_type("CONSTANT"),
      d_refinement_criteria(),
      d_dev_tol(),
      d_dev(),
      d_dev_time_max(),
      d_dev_time_min(),
      d_grad_tol(),
      d_grad_time_max(),
      d_grad_time_min()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
    TBOX_ASSERT(!godunov_advector.isNull());
    TBOX_ASSERT(!grid_geom.isNull());
#endif

    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from given input/restart databases.
    bool is_from_restart = RestartManager::getManager()->isFromRestart();
    if (is_from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, is_from_restart);
    return;
}// AdvectHypPatchOps

AdvectHypPatchOps::~AdvectHypPatchOps()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
    return;
}// ~AdvectHypPatchOps

const std::string&
AdvectHypPatchOps::getName() const
{
    return d_object_name;
}// getName

void
AdvectHypPatchOps::registerVisItDataWriter(
    Pointer<VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!visit_writer.isNull());
#endif
    d_visit_writer = visit_writer;
    return;
}// registerVisItDataWriter

void
AdvectHypPatchOps::registerAdvectionVelocity(
    Pointer<FaceVariable<NDIM,double> > u_var,
    const bool manage_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!u_var.isNull());
#endif
    d_u_var.insert(u_var);
    d_manage_u_data[u_var] = manage_data;
    d_u_is_div_free[u_var] = true;
    return;
}// registerAdvectionVelocity

void
AdvectHypPatchOps::setAdvectionVelocityIsDivergenceFree(
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
AdvectHypPatchOps::setAdvectionVelocityFunction(
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
AdvectHypPatchOps::registerSourceTerm(
    Pointer<CellVariable<NDIM,double> > F_var,
    const bool manage_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!F_var.isNull());
#endif
    d_F_var.insert(F_var);
    d_manage_F_data[F_var] = manage_data;
    return;
}// registerSourceTerm

void
AdvectHypPatchOps::setSourceTermFunction(
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
AdvectHypPatchOps::registerIncompressibilityFixTerm(
    Pointer<FaceVariable<NDIM,double> > grad_Phi_var,
    const bool manage_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!grad_Phi_var.isNull());
#endif
    d_grad_Phi_var.insert(grad_Phi_var);
    d_manage_grad_Phi_data[grad_Phi_var] = manage_data;
    return;
}// registerIncompressibilityFixTerm

void
AdvectHypPatchOps::setIncompressibilityFixTermFunction(
    Pointer<FaceVariable<NDIM,double> > grad_Phi_var,
    Pointer<IBTK::CartGridFunction> grad_Phi_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_grad_Phi_var.find(grad_Phi_var) != d_grad_Phi_var.end());
#endif
    d_grad_Phi_fcn[grad_Phi_var] = grad_Phi_fcn;
    return;
}// setIncompressibilityFixTermFunction

void
AdvectHypPatchOps::registerTransportedQuantity(
    Pointer<CellVariable<NDIM,double> > Q_var,
    const bool manage_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!Q_var.isNull());
#endif
    d_Q_var.insert(Q_var);
    d_manage_Q_data[Q_var] = manage_data;
    d_Q_difference_form[Q_var] = CONSERVATIVE;
    return;
}// registerTransportedQuantity

void
AdvectHypPatchOps::setAdvectionVelocity(
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
AdvectHypPatchOps::setSourceTerm(
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
AdvectHypPatchOps::setIncompressibilityFixTerm(
    Pointer<CellVariable<NDIM,double> > Q_var,
    Pointer<FaceVariable<NDIM,double> > grad_Phi_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    TBOX_ASSERT(d_grad_Phi_var.find(grad_Phi_var) != d_grad_Phi_var.end());
#endif
    d_Q_grad_Phi_map[Q_var] = grad_Phi_var;
    return;
}// setIncompressibilityFixTerm

void
AdvectHypPatchOps::setConvectiveDifferencingType(
    Pointer<CellVariable<NDIM,double> > Q_var,
    const ConvectiveDifferencingType& difference_form)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_difference_form[Q_var] = difference_form;
    return;
}// setConvectiveDifferencingType

void
AdvectHypPatchOps::setInitialConditions(
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
AdvectHypPatchOps::setBoundaryConditions(
    Pointer<CellVariable<NDIM,double> > Q_var,
    RobinBcCoefStrategy<NDIM>* Q_bc_coef)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
    const unsigned Q_depth = Q_factory->getDefaultDepth();
    TBOX_ASSERT(Q_depth == 1);
#endif
    d_Q_bc_coef[Q_var] = std::vector<RobinBcCoefStrategy<NDIM>*>(1,Q_bc_coef);
    return;
}// setBoundaryConditions

void
AdvectHypPatchOps::setBoundaryConditions(
    Pointer<CellVariable<NDIM,double> > Q_var,
    std::vector<RobinBcCoefStrategy<NDIM>*> Q_bc_coef)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
    const unsigned Q_depth = Q_factory->getDefaultDepth();
    TBOX_ASSERT(Q_depth == Q_bc_coef.size());
#endif
    d_Q_bc_coef[Q_var] = Q_bc_coef;
    return;
}// setBoundaryConditions

///
///  The following routines:
///
///      registerModelVariables(),
///      initializeDataOnPatch(),
///      computeStableDtOnPatch(),
///      computeFluxesOnPatch(),
///      conservativeDifferenceOnPatch(),
///      preprocessAdvanceLevelState(),
///      postprocessAdvanceLevelState(),
///      tagGradientDetectorCells()
///
///  are concrete implementations of functions declared in the
///  HyperbolicPatchStrategy abstract base class.
///

void
AdvectHypPatchOps::registerModelVariables(
    HyperbolicLevelIntegrator<NDIM>* integrator)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(integrator != static_cast<HyperbolicLevelIntegrator<NDIM>*>(NULL));
#endif
    d_integrator = integrator;

    for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin();
         cit != d_u_var.end(); ++cit)
    {
        Pointer<FaceVariable<NDIM,double> > u_var = *cit;
        if (d_manage_u_data[u_var])
        {
            d_integrator->registerVariable(
                u_var, d_ghosts,
                HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
        }
    }

    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_F_var.begin();
         cit != d_F_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > F_var = *cit;
        if (d_manage_F_data[F_var])
        {
            d_integrator->registerVariable(
                F_var, d_ghosts,
                HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
        }
    }

    for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_grad_Phi_var.begin();
         cit != d_grad_Phi_var.end(); ++cit)
    {
        Pointer<FaceVariable<NDIM,double> > grad_Phi_var = *cit;
        if (d_manage_grad_Phi_data[grad_Phi_var])
        {
            d_integrator->registerVariable(
                grad_Phi_var, d_ghosts,
                HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
        }
    }

    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
        if (d_manage_Q_data[Q_var])
        {
            d_integrator->registerVariable(
                Q_var, d_ghosts,
                HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
        }

        if (!d_visit_writer.isNull())
        {
            const int Q_idx = VariableDatabase<NDIM>::getDatabase()->mapVariableAndContextToIndex(Q_var, d_integrator->getPlotContext());
            const int depth = Q_factory->getDefaultDepth();
            if (depth == 1)
            {
                d_visit_writer->registerPlotQuantity(Q_var->getName(), "SCALAR", Q_idx);
            }
            else
            {
                if (depth == NDIM)
                {
                    d_visit_writer->registerPlotQuantity(Q_var->getName(), "VECTOR", Q_idx);
                }
                for (int d = 0; d < depth; ++d)
                {
                    std::ostringstream stream;
                    stream << d;
                    d_visit_writer->registerPlotQuantity(Q_var->getName()+"_"+stream.str(), "SCALAR", Q_idx, d);
                }
            }
        }

        const bool conservation_form = d_Q_difference_form[Q_var] == CONSERVATIVE;
        if (conservation_form)
        {
            d_flux_integral_var[Q_var] = new FaceVariable<NDIM,double>(
                d_object_name+"::"+Q_var->getName()+" advective flux time integral",
                Q_factory->getDefaultDepth());
            d_integrator->registerVariable(
                d_flux_integral_var[Q_var], d_flux_ghosts,
                HyperbolicLevelIntegrator<NDIM>::FLUX,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "NO_REFINE");
        }

        Pointer<FaceVariable<NDIM,double> > u_var = d_Q_u_map[Q_var];
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!u_var.isNull());
#endif
        const bool u_is_div_free = d_u_is_div_free[u_var];
        if (!conservation_form || !u_is_div_free)
        {
            d_q_integral_var[Q_var] = new FaceVariable<NDIM,double>(
                d_object_name+"::"+Q_var->getName()+" time integral",
                Q_factory->getDefaultDepth());
            d_integrator->registerVariable(
                d_q_integral_var[Q_var], d_flux_ghosts,
                HyperbolicLevelIntegrator<NDIM>::FLUX,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "NO_REFINE");
            if (d_u_integral_var[u_var].isNull())
            {
                d_u_integral_var[u_var] = new FaceVariable<NDIM,double>(
                    d_object_name+"::"+u_var->getName()+" time integral");
                d_integrator->registerVariable(
                    d_u_integral_var[u_var], d_flux_ghosts,
                    HyperbolicLevelIntegrator<NDIM>::FLUX,
                    d_grid_geometry,
                    "CONSERVATIVE_COARSEN",
                    "NO_REFINE");
            }
        }
    }
    return;
}// registerModelVariables

void
AdvectHypPatchOps::initializeDataOnPatch(
    Patch<NDIM>& patch,
    const double data_time,
    const bool initial_time)
{
    if (initial_time)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin();
             cit != d_u_var.end(); ++cit)
        {
            Pointer<FaceVariable<NDIM,double> > u_var = *cit;
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, getDataContext());
            if (d_u_fcn[u_var].isNull())
            {
                Pointer<FaceData<NDIM,double> > u_data = patch.getPatchData(u_idx);
                u_data->fillAll(0.0);
            }
            else
            {
                d_u_fcn[u_var]->setDataOnPatch(u_idx, u_var, Pointer<Patch<NDIM> >(&patch,false), data_time, initial_time);
            }
        }

        for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_F_var.begin();
             cit != d_F_var.end(); ++cit)
        {
            Pointer<CellVariable<NDIM,double> > F_var = *cit;
            const int F_idx = var_db->mapVariableAndContextToIndex(F_var, getDataContext());
            if (d_F_fcn[F_var].isNull())
            {
                Pointer<CellData<NDIM,double> > F_data = patch.getPatchData(F_idx);
                F_data->fillAll(0.0);
            }
            else
            {
                d_F_fcn[F_var]->setDataOnPatch(F_idx, F_var, Pointer<Patch<NDIM> >(&patch,false), data_time, initial_time);
            }
        }

        for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_grad_Phi_var.begin();
             cit != d_grad_Phi_var.end(); ++cit)
        {
            Pointer<FaceVariable<NDIM,double> > grad_Phi_var = *cit;
            const int grad_Phi_idx = var_db->mapVariableAndContextToIndex(grad_Phi_var, getDataContext());
            if (d_grad_Phi_fcn[grad_Phi_var].isNull())
            {
                Pointer<FaceData<NDIM,double> > grad_Phi_data = patch.getPatchData(grad_Phi_idx);
                grad_Phi_data->fillAll(0.0);
            }
            else
            {
                d_grad_Phi_fcn[grad_Phi_var]->setDataOnPatch(grad_Phi_idx, grad_Phi_var, Pointer<Patch<NDIM> >(&patch,false), data_time, initial_time);
            }
        }

        for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
             cit != d_Q_var.end(); ++cit)
        {
            Pointer<CellVariable<NDIM,double> > Q_var = *cit;
            const int Q_idx = var_db->mapVariableAndContextToIndex(Q_var, getDataContext());
            if (d_Q_init[Q_var].isNull())
            {
                Pointer<CellData<NDIM,double> > Q_data = patch.getPatchData(Q_var, getDataContext());
                Q_data->fillAll(0.0);
            }
            else
            {
                d_Q_init[Q_var]->setDataOnPatch(Q_idx, Q_var, Pointer<Patch<NDIM> >(&patch,false), data_time, initial_time);
            }
        }
    }
    return;
}// initializeDataOnPatch

double
AdvectHypPatchOps::computeStableDtOnPatch(
    Patch<NDIM>& patch,
    const bool initial_time,
    const double dt_time)
{
    double stable_dt = std::numeric_limits<double>::max();
    for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin();
         cit != d_u_var.end(); ++cit)
    {
        Pointer<FaceVariable<NDIM,double> > u_var = *cit;
        Pointer<FaceData<NDIM,double> > u_data = patch.getPatchData(u_var, getDataContext());
        stable_dt = std::min(stable_dt,d_godunov_advector->computeStableDtOnPatch(*u_data,patch));
    }
    return stable_dt;
}// computeStableDtOnPatch

void
AdvectHypPatchOps::computeFluxesOnPatch(
    Patch<NDIM>& patch,
    const double time,
    const double dt)
{
    Pointer<PatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const Box<NDIM>& patch_box = patch.getBox();

    PatchFaceDataOpsReal<NDIM,double> patch_fc_data_ops;

    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        Pointer<FaceVariable<NDIM,double> > u_var = d_Q_u_map[Q_var];

        Pointer<CellData<NDIM,double> > Q_data = patch.getPatchData(Q_var, getDataContext());
        Pointer<FaceData<NDIM,double> > u_data = patch.getPatchData(u_var, getDataContext());
        Pointer<FaceData<NDIM,double> > q_integral_data = getQIntegralData(Q_var, patch, getDataContext());

        // Predict time- and face-centered values.
        if (!d_Q_F_map[Q_var].isNull())
        {
            Pointer<CellData<NDIM,double> > F_data = patch.getPatchData(d_Q_F_map[Q_var], getDataContext());
            d_godunov_advector->predictValueWithSourceTerm(*q_integral_data, *u_data, *Q_data, *F_data, patch, dt);
        }
        else
        {
            d_godunov_advector->predictValue(*q_integral_data, *u_data, *Q_data, patch, dt);
        }

        // For incompressible flow problems, we allow for the specification of
        // an auxiliary gradient that is used to enforce the incompressibility
        // constraint in an extremely approximate manner.
        Pointer<FaceVariable<NDIM,double> > grad_Phi_var = d_Q_grad_Phi_map[Q_var];
        if (!grad_Phi_var.isNull())
        {
            Pointer<FaceData<NDIM,double> > grad_Phi_data = patch.getPatchData(grad_Phi_var, getDataContext());
            if (!d_grad_Phi_fcn[grad_Phi_var].isNull() && d_grad_Phi_fcn[grad_Phi_var]->isTimeDependent())
            {
                VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
                const int grad_Phi_idx = var_db->mapVariableAndContextToIndex(grad_Phi_var, getDataContext());
                d_grad_Phi_fcn[grad_Phi_var]->setDataOnPatch(grad_Phi_idx, grad_Phi_var, Pointer<Patch<NDIM> >(&patch,false), time+0.5*dt);
            }
            d_godunov_advector->enforceIncompressibility(*q_integral_data, *u_data, *grad_Phi_data, patch);
        }
    }

    // Set physical boundary conditions for the face-centered predicted values
    // at inflow boundaries only.
    if (pgeom->getTouchesRegularBoundary())
    {
        setInflowBoundaryConditions(patch, time+0.5*dt);
    }

    // Update the advection velocity.
    if (d_compute_half_velocity)
    {
        for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin();
             cit != d_u_var.end(); ++cit)
        {
            Pointer<FaceVariable<NDIM,double> > u_var = *cit;
            if (!d_u_fcn[u_var].isNull() && d_u_fcn[u_var]->isTimeDependent())
            {
                VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
                const int u_idx = var_db->mapVariableAndContextToIndex(u_var, getDataContext());
                d_u_fcn[u_var]->setDataOnPatch(u_idx, u_var, Pointer<Patch<NDIM> >(&patch,false), time+0.5*dt);
            }
        }
    }

    // Compute fluxes and other face-centered quantities.
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        Pointer<FaceVariable<NDIM,double> > u_var = d_Q_u_map[Q_var];

        Pointer<FaceData<NDIM,double> > u_data = patch.getPatchData(u_var, getDataContext());

        const bool conservation_form = d_Q_difference_form[Q_var] == CONSERVATIVE;
        const bool u_is_div_free = d_u_is_div_free[u_var];

        if (conservation_form)
        {
            Pointer<FaceData<NDIM,double> > flux_integral_data = getFluxIntegralData(Q_var, patch, getDataContext());
            Pointer<FaceData<NDIM,double> > q_integral_data = getQIntegralData(Q_var, patch, getDataContext());
            d_godunov_advector->computeFlux(*flux_integral_data, *u_data, *q_integral_data, patch, dt);
        }
        else if (!conservation_form || !u_is_div_free)
        {
            Pointer<FaceData<NDIM,double> > q_integral_data = getQIntegralData(Q_var, patch, getDataContext());
            patch_fc_data_ops.scale(q_integral_data, dt, q_integral_data, patch_box);

            Pointer<FaceData<NDIM,double> > u_integral_data = getUIntegralData(Q_var, patch, getDataContext());
            patch_fc_data_ops.scale(u_integral_data, dt, u_data, patch_box);
        }
    }
    return;
}// computeFluxesOnPatch

void
AdvectHypPatchOps::conservativeDifferenceOnPatch(
    Patch<NDIM>& patch,
    const double time,
    const double dt,
    bool at_synchronization)
{
    const Box<NDIM>& patch_box = patch.getBox();
    const Index<NDIM>& ilower = patch_box.lower();
    const Index<NDIM>& iupper = patch_box.upper();

    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        Pointer<FaceVariable<NDIM,double> > u_var = d_Q_u_map[Q_var];

        Pointer<CellData<NDIM,double> > Q_data = patch.getPatchData(Q_var, getDataContext());
        Pointer<FaceData<NDIM,double> > flux_integral_data = getFluxIntegralData(Q_var, patch, getDataContext());
        Pointer<FaceData<NDIM,double> > q_integral_data = getQIntegralData(Q_var, patch, getDataContext());
        Pointer<FaceData<NDIM,double> > u_integral_data = getUIntegralData(Q_var, patch, getDataContext());

        const IntVector<NDIM>& Q_data_ghost_cells = Q_data->getGhostCellWidth();
        const IntVector<NDIM>& flux_integral_data_ghost_cells =
            (!flux_integral_data.isNull()
             ? flux_integral_data->getGhostCellWidth()
             : 0);
        const IntVector<NDIM>& q_integral_data_ghost_cells =
            (!q_integral_data.isNull()
             ? q_integral_data->getGhostCellWidth()
             : 0);
        const IntVector<NDIM>& u_integral_data_ghost_cells =
            (!u_integral_data.isNull()
             ? u_integral_data->getGhostCellWidth()
             : 0);

        const bool u_is_div_free = d_u_is_div_free[u_var];

        switch (d_Q_difference_form[Q_var])
        {
            case CONSERVATIVE:
            {
                for (int depth = 0; depth < Q_data->getDepth(); ++depth)
                {
                    if (u_is_div_free)
                    {
#if (NDIM == 2)
                        ADVECT_CONSDIFF_FC(
                            dx,
                            ilower(0),iupper(0),ilower(1),iupper(1),
                            flux_integral_data_ghost_cells(0),flux_integral_data_ghost_cells(1),
                            Q_data_ghost_cells(0),Q_data_ghost_cells(1),
                            flux_integral_data->getPointer(0,depth),
                            flux_integral_data->getPointer(1,depth),
                            Q_data->getPointer(depth));
#endif
#if (NDIM == 3)
                        ADVECT_CONSDIFF_FC(
                            dx,
                            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                            flux_integral_data_ghost_cells(0),flux_integral_data_ghost_cells(1),flux_integral_data_ghost_cells(2),
                            Q_data_ghost_cells(0),Q_data_ghost_cells(1),Q_data_ghost_cells(2),
                            flux_integral_data->getPointer(0,depth),
                            flux_integral_data->getPointer(1,depth),
                            flux_integral_data->getPointer(2,depth),
                            Q_data->getPointer(depth));
#endif
                    }
                    else
                    {
#if (NDIM == 2)
                        ADVECT_CONSDIFFWITHDIVSOURCE_FC(
                            dx,dt,
                            ilower(0),iupper(0),ilower(1),iupper(1),
                            flux_integral_data_ghost_cells(0),flux_integral_data_ghost_cells(1),
                            q_integral_data_ghost_cells(0),q_integral_data_ghost_cells(1),
                            u_integral_data_ghost_cells(0),u_integral_data_ghost_cells(1),
                            Q_data_ghost_cells(0),Q_data_ghost_cells(1),
                            flux_integral_data->getPointer(0,depth),
                            flux_integral_data->getPointer(1,depth),
                            q_integral_data->getPointer(0,depth),
                            q_integral_data->getPointer(1,depth),
                            u_integral_data->getPointer(0),
                            u_integral_data->getPointer(1),
                            Q_data->getPointer(depth));
#endif
#if (NDIM == 3)
                        ADVECT_CONSDIFFWITHDIVSOURCE_FC(
                            dx,dt,
                            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                            flux_integral_data_ghost_cells(0),flux_integral_data_ghost_cells(1),flux_integral_data_ghost_cells(2),
                            q_integral_data_ghost_cells(0),q_integral_data_ghost_cells(1),q_integral_data_ghost_cells(2),
                            u_integral_data_ghost_cells(0),u_integral_data_ghost_cells(1),u_integral_data_ghost_cells(2),
                            Q_data_ghost_cells(0),Q_data_ghost_cells(1),Q_data_ghost_cells(2),
                            flux_integral_data->getPointer(0,depth),
                            flux_integral_data->getPointer(1,depth),
                            flux_integral_data->getPointer(2,depth),
                            q_integral_data->getPointer(0,depth),
                            q_integral_data->getPointer(1,depth),
                            q_integral_data->getPointer(2,depth),
                            u_integral_data->getPointer(0),
                            u_integral_data->getPointer(1),
                            u_integral_data->getPointer(2),
                            Q_data->getPointer(depth));
#endif
                    }
                }
                break;
            }
            case ADVECTIVE:
            {
                CellData<NDIM,double> N_data(patch_box,Q_data->getDepth(),0);
                d_godunov_advector->computeAdvectiveDerivative(N_data, *u_integral_data, *q_integral_data, patch);
                PatchCellDataOpsReal<NDIM,double> patch_cc_data_ops;
                patch_cc_data_ops.axpy(Q_data, -1.0/dt, Pointer<CellData<NDIM,double> >(&N_data,false), Q_data, patch_box);
                break;
            }
            default:
            {
                TBOX_ERROR("AdvectHypPatchOps::conservativeDifferenceOnPatch():\n"
                           << "  unsupported differencing form: " << enum_to_string<ConvectiveDifferencingType>(d_Q_difference_form[Q_var]) << " \n"
                           << "  valid choices are: ADVECTIVE, CONSERVATIVE\n");
            }
        }
    }
    return;
}// conservativeDifferenceOnPatch

void
AdvectHypPatchOps::preprocessAdvanceLevelState(
    const Pointer<PatchLevel<NDIM> >& level,
    double current_time,
    double dt,
    bool first_step,
    bool last_step,
    bool regrid_advance)
{
    if (!d_compute_init_velocity) return;

    // Update the advection velocity.
    for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin();
         cit != d_u_var.end(); ++cit)
    {
        Pointer<FaceVariable<NDIM,double> > u_var = *cit;
        if (!d_u_fcn[u_var].isNull() && d_u_fcn[u_var]->isTimeDependent())
        {
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, d_integrator->getScratchContext());
            d_u_fcn[u_var]->setDataOnPatchLevel(u_idx, u_var, level, current_time);
        }
    }
    return;
}// preprocessAdvanceLevelState

void
AdvectHypPatchOps::postprocessAdvanceLevelState(
    const Pointer<PatchLevel<NDIM> >& level,
    double current_time,
    double dt,
    bool first_step,
    bool last_step,
    bool regrid_advance)
{
    PatchCellDataOpsReal<NDIM,double> patch_cc_data_ops;

    Pointer<VariableContext> new_context     = d_integrator->getNewContext();
    Pointer<VariableContext> scratch_context = d_integrator->getScratchContext();

    // Update the values of any time-dependent source terms and add the values
    // of all source terms to the advected quantities.
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        Pointer<CellVariable<NDIM,double> > F_var = d_Q_F_map[Q_var];
        if (F_var.isNull()) continue;
        Pointer<CartGridFunction> F_fcn = d_F_fcn[F_var];
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM,double> > Q_data = patch->getPatchData(Q_var, new_context);
            Pointer<CellData<NDIM,double> > F_data = patch->getPatchData(F_var, scratch_context);

            if (!F_fcn.isNull())
            {
                patch_cc_data_ops.axpy(Q_data, 0.5*dt, F_data, Q_data, patch_box);
            }
            else
            {
                patch_cc_data_ops.axpy(Q_data, dt, F_data, Q_data, patch_box);
            }
        }
    }

    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_F_var.begin();
         cit != d_F_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > F_var = *cit;
        Pointer<CartGridFunction> F_fcn = d_F_fcn[F_var];
        if (!F_fcn.isNull() && F_fcn->isTimeDependent())
        {
            F_fcn->setDataOnPatchLevel(F_var, new_context, level, current_time+dt);
        }
    }

    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        Pointer<CellVariable<NDIM,double> > F_var = d_Q_F_map[Q_var];
        if (F_var.isNull()) continue;
        Pointer<CartGridFunction> F_fcn = d_F_fcn[F_var];
        if (!F_fcn.isNull())
        {
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();

                Pointer<CellData<NDIM,double> > Q_data = patch->getPatchData(Q_var, new_context);
                Pointer<CellData<NDIM,double> > F_data = patch->getPatchData(F_var, scratch_context);

                patch_cc_data_ops.axpy(Q_data, 0.5*dt, F_data, Q_data, patch_box);
            }
        }
    }

    if (!d_compute_final_velocity) return;

    // Update the advection velocity (or velocities).
    for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin();
         cit != d_u_var.end(); ++cit)
    {
        Pointer<FaceVariable<NDIM,double> > u_var = *cit;
        if (!d_u_fcn[u_var].isNull() && d_u_fcn[u_var]->isTimeDependent())
        {
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, d_integrator->getScratchContext());
            d_u_fcn[u_var]->setDataOnPatchLevel(u_idx, u_var, level, current_time+dt);
        }
    }
    return;
}// postprocessAdvanceLevelState

void
AdvectHypPatchOps::tagGradientDetectorCells(
    Patch<NDIM>& patch,
    const double regrid_time,
    const bool initial_error,
    const int tag_indx,
    const bool uses_richardson_extrapolation_too)
{
    const int error_level_number = patch.getPatchLevelNumber();

    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const Box<NDIM>& patch_box = patch.getBox();
    const Index<NDIM>& ilower = patch.getBox().lower();
    const Index<NDIM>& iupper = patch.getBox().upper();

    Pointer<CellData<NDIM,int> > tags = patch.getPatchData(tag_indx);

    const int not_refine_tag_val = FALSE_VAL;
    const int refine_tag_val = TRUE_VAL;

    // Create a set of temporary tags and set to untagged value.
    CellData<NDIM,int> temp_tags(patch_box, 1, d_ghosts);
    temp_tags.fillAll(not_refine_tag_val);

    // Possible tagging criteria include: QVAL_DEVIATION, QVAL_GRADIENT.  The
    // criteria are specified over a time interval.
    //
    // Loop over criteria provided and check to make sure we are in the
    // specified time interval.  If so, apply appropriate tagging for the level.
    for (int ncrit = 0; ncrit < d_refinement_criteria.getSize(); ++ncrit)
    {
        std::string ref = d_refinement_criteria[ncrit];
        const IntVector<NDIM>& tagghost = tags->getGhostCellWidth();

        int size = 0;
        double tol = 0.0;
        bool time_allowed = false;

        if (ref == "QVAL_DEVIATION")
        {
            size = d_dev_tol.getSize();
            tol = (error_level_number < size
                   ? d_dev_tol[error_level_number]
                   : d_dev_tol[size-1]);
            size = d_dev.getSize();
            double dev = (error_level_number < size
                          ? d_dev[error_level_number]
                          : d_dev[size-1]);
            size = d_dev_time_min.getSize();
            double time_min = (error_level_number < size
                               ? d_dev_time_min[error_level_number]
                               : d_dev_time_min[size-1]);
            size = d_dev_time_max.getSize();
            double time_max = (error_level_number < size
                               ? d_dev_time_max[error_level_number]
                               : d_dev_time_max[size-1]);
            time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);

            if (time_allowed)
            {
                // Check for tags that have already been set in a previous step.
                for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
                     cit != d_Q_var.end(); ++cit)
                {
                    Pointer<CellVariable<NDIM,double> > Q_var = *cit;
                    Pointer<CellData<NDIM,double> > Q_data = patch.getPatchData(Q_var, getDataContext());
                    for (int depth = 0; depth < Q_data->getDepth(); ++depth)
                    {
                        for (CellIterator<NDIM> ic(patch_box); ic; ic++)
                        {
                            double locden = tol;
                            int tag_val = (*tags)(ic(),0);
                            if (tag_val)
                            {
                                locden *= 0.75;
                            }
                            if (std::abs((*Q_data)(ic(),depth)-dev) > locden)
                            {
                                temp_tags(ic(),0) = refine_tag_val;
                            }
                        }
                    }
                }
            }
        }

        if (ref == "QVAL_GRADIENT")
        {
            size = d_grad_tol.getSize();
            tol = (error_level_number < size
                   ? d_grad_tol[error_level_number]
                   : d_grad_tol[size-1]);
            size = d_grad_time_min.getSize();
            double time_min = (error_level_number < size
                               ? d_grad_time_min[error_level_number]
                               : d_grad_time_min[size-1]);
            size = d_grad_time_max.getSize();
            double time_max = (error_level_number < size
                               ? d_grad_time_max[error_level_number]
                               : d_grad_time_max[size-1]);
            time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);

            if (time_allowed)
            {
                for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
                     cit != d_Q_var.end(); ++cit)
                {
                    Pointer<CellVariable<NDIM,double> > Q_var = *cit;
                    Pointer<CellData<NDIM,double> > Q_data = patch.getPatchData(Q_var, getDataContext());
                    const IntVector<NDIM>& Q_ghost = Q_data->getGhostCellWidth();
                    for (int depth = 0; depth < Q_data->getDepth(); ++depth)
                    {
                        ADVECT_DETECTGRAD_FC(
#if (NDIM == 2)
                            ilower(0),iupper(0), ilower(1),iupper(1),
                            Q_ghost(0),tagghost(0),d_ghosts(0),
                            Q_ghost(1),tagghost(1),d_ghosts(1),
#endif
#if (NDIM == 3)
                            ilower(0),iupper(0), ilower(1),iupper(1),ilower(2),iupper(2),
                            Q_ghost(0),tagghost(0),d_ghosts(0),
                            Q_ghost(1),tagghost(1),d_ghosts(1),
                            Q_ghost(2),tagghost(2),d_ghosts(2),
#endif
                            dx,
                            tol,
                            refine_tag_val,
                            Q_data->getPointer(depth),
                            tags->getPointer(),temp_tags.getPointer());
                    }
                }
            }
        }
    } // loop over criteria

    // Update tags.
    for (CellIterator<NDIM> ic(patch_box); ic; ic++)
    {
        (*tags)(ic(),0) = temp_tags(ic(),0);
    }
    return;
}// tagGradientDetectorCells

///
///  The following routine:
///
///      setPhysicalBoundaryConditions()
///
///  is a concrete implementations of functions declared in the
///  RefinePatchStrategy abstract base class.
///

void
AdvectHypPatchOps::setPhysicalBoundaryConditions(
    Patch<NDIM>& patch,
    const double fill_time,
    const IntVector<NDIM>& ghost_width_to_fill)
{
    // Extrapolate the interior data to set the ghost cell values for the state
    // variables and for any forcing terms.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    ComponentSelector u_patch_data_indices;
    for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin();
         cit != d_u_var.end(); ++cit)
    {
        Pointer<FaceVariable<NDIM,double> > u_var = *cit;
        const int u_data_idx = var_db->mapVariableAndContextToIndex(u_var, d_integrator->getScratchContext());
        u_patch_data_indices.setFlag(u_data_idx);
    }
    d_extrap_bc_helper.setExtrapolationType(d_extrap_type);
    d_extrap_bc_helper.setPatchDataIndices(u_patch_data_indices);
    d_extrap_bc_helper.setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);

    ComponentSelector patch_data_indices;
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_F_var.begin();
         cit != d_F_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > F_var = *cit;
        const int F_data_idx = var_db->mapVariableAndContextToIndex(F_var, d_integrator->getScratchContext());
        patch_data_indices.setFlag(F_data_idx);
    }
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        const int Q_data_idx = var_db->mapVariableAndContextToIndex(Q_var, d_integrator->getScratchContext());
        patch_data_indices.setFlag(Q_data_idx);
    }
    d_extrap_bc_helper.setExtrapolationType(d_extrap_type);
    d_extrap_bc_helper.setPatchDataIndices(patch_data_indices);
    d_extrap_bc_helper.setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);
    return;
}// setPhysicalBoundaryConditions

///
///  The following routines:
///
///      putToDatabase()
///
///  are concrete implementations of functions declared in the
///  Serializable abstract base class.
///

void
AdvectHypPatchOps::putToDatabase(
    Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif

    db->putInteger("ADVECT_HYP_PATCH_OPS_VERSION", ADVECT_HYP_PATCH_OPS_VERSION);

    db->putIntegerArray("d_ghosts", static_cast<int*>(d_ghosts), NDIM);
    db->putIntegerArray("d_flux_ghosts", static_cast<int*>(d_flux_ghosts), NDIM);

    if (d_refinement_criteria.getSize() > 0)
    {
        db->putStringArray("d_refinement_criteria", d_refinement_criteria);
    }
    for (int i = 0; i < d_refinement_criteria.getSize(); ++i)
    {
        if (d_refinement_criteria[i] == "QVAL_DEVIATION")
        {
            db->putDoubleArray("d_dev_tol", d_dev_tol);
            db->putDoubleArray("d_dev", d_dev);
            db->putDoubleArray("d_dev_time_max", d_dev_time_max);
            db->putDoubleArray("d_dev_time_min", d_dev_time_min);
        }
        else if (d_refinement_criteria[i] == "QVAL_GRADIENT")
        {
            db->putDoubleArray("d_grad_tol", d_grad_tol);
            db->putDoubleArray("d_grad_time_max", d_grad_time_max);
            db->putDoubleArray("d_grad_time_min", d_grad_time_min);
        }
    }
    return;
}// putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

Pointer<FaceData<NDIM,double> >
AdvectHypPatchOps::getFluxIntegralData(
    Pointer<CellVariable<NDIM,double> > Q_var,
    Patch<NDIM>& patch,
    Pointer<VariableContext> context)
{
    if (d_Q_difference_form[Q_var] == CONSERVATIVE)
    {
        return patch.getPatchData(d_flux_integral_var[Q_var], context);
    }
    else
    {
        return Pointer<FaceData<NDIM,double> >(NULL);
    }
}// getFluxIntegralData

Pointer<FaceData<NDIM,double> >
AdvectHypPatchOps::getQIntegralData(
    Pointer<CellVariable<NDIM,double> > Q_var,
    Patch<NDIM>& patch,
    Pointer<VariableContext> context)
{
    if (!d_u_is_div_free[d_Q_u_map[Q_var]] || d_Q_difference_form[Q_var] != CONSERVATIVE)
    {
        return patch.getPatchData(d_q_integral_var[Q_var], context);
    }
    else
    {
        return getFluxIntegralData(Q_var, patch, context);
    }
}// getQIntegralData

Pointer<FaceData<NDIM,double> >
AdvectHypPatchOps::getUIntegralData(
    Pointer<CellVariable<NDIM,double> > Q_var,
    Patch<NDIM>& patch,
    Pointer<VariableContext> context)
{
    Pointer<FaceVariable<NDIM,double> > u_var = d_Q_u_map[Q_var];
    if (!d_u_is_div_free[u_var] || d_Q_difference_form[Q_var] != CONSERVATIVE)
    {
        return patch.getPatchData(d_u_integral_var[u_var], context);
    }
    else
    {
        return Pointer<FaceData<NDIM,double> >(NULL);
    }
}// getUIntegralData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
AdvectHypPatchOps::setInflowBoundaryConditions(
    Patch<NDIM>& patch,
    const double fill_time)
{
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

    // There is nothing to do if the patch does not touch a regular (physical)
    // boundary.
    if (!pgeom->getTouchesRegularBoundary()) return;

    // Compute the co-dimension one boundary boxes.
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);

    // There is nothing to do if the patch does not have any co-dimension one
    // boundary boxes.
    if (physical_codim1_boxes.size() == 0) return;

    // Loop over the boundary fill boxes and set boundary conditions at inflow
    // boundaries only.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const Box<NDIM>& patch_box = patch.getBox();
    const double* const dx = pgeom->getDx();
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        Pointer<FaceVariable<NDIM,double> > u_var = d_Q_u_map[Q_var];

        const int Q_data_idx = var_db->mapVariableAndContextToIndex(Q_var, d_integrator->getScratchContext());

        Pointer<FaceData<NDIM,double> > u_data = patch.getPatchData(u_var, getDataContext());
        Pointer<FaceData<NDIM,double> > q_integral_data = getQIntegralData(Q_var, patch, getDataContext());

        // Setup any extended Robin BC coef objects.
        for (int depth = 0; depth < q_integral_data->getDepth(); ++depth)
        {
            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_Q_bc_coef[Q_var][depth]);
            if (extended_bc_coef != NULL)
            {
                extended_bc_coef->setHomogeneousBc(false);
                extended_bc_coef->setTargetPatchDataIndex(Q_data_idx);
            }
        }

        // Set the boundary conditions.
        for (int n = 0; n < physical_codim1_boxes.size(); ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
            const int location_index   = bdry_box.getLocationIndex();
            const int bdry_normal_axis = location_index/2;
            const bool is_lower        = location_index%2 == 0;

            static const IntVector<NDIM> gcw_to_fill = 1;
            const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
            const BoundaryBox<NDIM> trimmed_bdry_box(bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), bdry_box.getLocationIndex());
            const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            // Loop over the boundary box indices and compute the nearest
            // interior index.
            for (int depth = 0; depth < q_integral_data->getDepth(); ++depth)
            {
                ArrayData<NDIM,double> acoef_data(bc_coef_box, 1);
                ArrayData<NDIM,double> bcoef_data(bc_coef_box, 1);
                ArrayData<NDIM,double> gcoef_data(bc_coef_box, 1);

                Pointer<ArrayData<NDIM,double> > acoef_data_ptr(&acoef_data,false);
                Pointer<ArrayData<NDIM,double> > bcoef_data_ptr(&bcoef_data,false);
                Pointer<ArrayData<NDIM,double> > gcoef_data_ptr(&gcoef_data,false);

                d_Q_bc_coef[Q_var][depth]->setBcCoefs(
                    acoef_data_ptr, bcoef_data_ptr, gcoef_data_ptr,
                    Q_var, patch, trimmed_bdry_box, fill_time);

                for (Box<NDIM>::Iterator b(bc_coef_box); b; b++)
                {
                    const Index<NDIM>& i = b();
                    const FaceIndex<NDIM> i_f(i, bdry_normal_axis, FaceIndex<NDIM>::Lower);

                    bool is_inflow_bdry = (is_lower && (*u_data)(i_f) > 0.0) || (!is_lower && (*u_data)(i_f) < 0.0);
                    if (is_inflow_bdry)
                    {
                        const double& a = acoef_data(i,0);
                        const double& b = bcoef_data(i,0);
                        const double& g = gcoef_data(i,0);
                        const double& h = dx[bdry_normal_axis];

                        Index<NDIM> i_intr(i);
                        if (is_lower)
                        {
                            // intentionally blank
                        }
                        else
                        {
                            i_intr(bdry_normal_axis) -= 1;
                        }

                        const FaceIndex<NDIM> i_f_bdry(i_intr, bdry_normal_axis, (is_lower ? FaceIndex<NDIM>::Lower : FaceIndex<NDIM>::Upper));
                        const FaceIndex<NDIM> i_f_intr(i_intr, bdry_normal_axis, (is_lower ? FaceIndex<NDIM>::Upper : FaceIndex<NDIM>::Lower));
                        const double& q_i = (*q_integral_data)(i_f_intr,depth);
                        const double q_b = (b*q_i + g*h)/(a*h + b);
                        (*q_integral_data)(i_f_bdry,depth) = q_b;
                    }
                }
            }
        }
    }
    return;
}// setInflowBoundaryConditions

void
AdvectHypPatchOps::getFromInput(
    Pointer<Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    d_compute_init_velocity  = db->getBoolWithDefault("compute_init_velocity" , d_compute_init_velocity);
    d_compute_half_velocity  = db->getBoolWithDefault("compute_half_velocity" , d_compute_half_velocity);
    d_compute_final_velocity = db->getBoolWithDefault("compute_final_velocity", d_compute_final_velocity);

    d_extrap_type = db->getStringWithDefault("extrap_type", d_extrap_type);
    if (!(d_extrap_type == "CONSTANT" || d_extrap_type == "LINEAR" || d_extrap_type == "QUADRATIC"))
    {
        TBOX_ERROR("AdvectHypPatchOps::getFromInput():\n"
                   << "  unknown extrapolation type: " << d_extrap_type << "\n"
                   << "  valid selections are: CONSTANT, LINEAR, or QUADRATIC" << std::endl);
    }

    if (db->keyExists("Refinement_data"))
    {
        Pointer<Database> refine_db = db->getDatabase("Refinement_data");
        Array<std::string> refinement_keys = refine_db->getAllKeys();
        int num_keys = refinement_keys.getSize();

        if (refine_db->keyExists("refine_criteria"))
        {
            d_refinement_criteria = refine_db->getStringArray("refine_criteria");
        }
        else
        {
            TBOX_WARNING(d_object_name << ":\n"
                         << "  No key `refine_criteria' found in data for"
                         << " Refinement_data. No refinement will occur.\n");
        }

        Array<std::string> ref_keys_defined(num_keys);
        int def_key_cnt = 0;
        Pointer<Database> error_db;

        for (int i = 0; i < refinement_keys.getSize(); ++i)
        {
            std::string error_key = refinement_keys[i];
            error_db.setNull();

            if (error_key != "refine_criteria")
            {
                if (!(error_key == "QVAL_DEVIATION" ||
                      error_key == "QVAL_GRADIENT"  ))
                {
                    TBOX_ERROR(d_object_name << ":\n"
                               << "  Unknown refinement criteria: " << error_key
                               << "\nin input.\n");
                }
                else
                {
                    error_db = refine_db->getDatabase(error_key);
                    ref_keys_defined[def_key_cnt] = error_key;
                    def_key_cnt++;
                }

                if (!error_db.isNull() && error_key == "QVAL_DEVIATION")
                {
                    if (error_db->keyExists("dev_tol"))
                    {
                        d_dev_tol = error_db->getDoubleArray("dev_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name << ":\n"
                                   << "  No key `dev_tol' found in data for "
                                   << error_key << "\n");
                    }

                    if (error_db->keyExists("qval_dev"))
                    {
                        d_dev = error_db->getDoubleArray("qval_dev");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name << ":\n"
                                   << "  No key `qval_dev' found in data for "
                                   << error_key << "\n");
                    }

                    if (error_db->keyExists("time_max"))
                    {
                        d_dev_time_max = error_db->getDoubleArray("time_max");
                    }
                    else
                    {
                        d_dev_time_max.resizeArray(1);
                        d_dev_time_max[0] = std::numeric_limits<double>::max();
                    }

                    if (error_db->keyExists("time_min"))
                    {
                        d_dev_time_min = error_db->getDoubleArray("time_min");
                    }
                    else
                    {
                        d_dev_time_min.resizeArray(1);
                        d_dev_time_min[0] = 0.0;
                    }
                }

                if (!error_db.isNull() && error_key == "QVAL_GRADIENT")
                {
                    if (error_db->keyExists("grad_tol"))
                    {
                        d_grad_tol = error_db->getDoubleArray("grad_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name << ":\n"
                                   << "  No key `grad_tol' found in data for "
                                   << error_key << "\n");
                    }

                    if (error_db->keyExists("time_max"))
                    {
                        d_grad_time_max = error_db->getDoubleArray("time_max");
                    }
                    else
                    {
                        d_grad_time_max.resizeArray(1);
                        d_grad_time_max[0] = std::numeric_limits<double>::max();
                    }

                    if (error_db->keyExists("time_min"))
                    {
                        d_grad_time_min = error_db->getDoubleArray("time_min");
                    }
                    else
                    {
                        d_grad_time_min.resizeArray(1);
                        d_grad_time_min[0] = 0.0;
                    }
                }

            }
        } // loop over refine criteria

        // Check that input is found for each string identifier in key list.
        for (int k0 = 0; k0 < d_refinement_criteria.getSize(); ++k0)
        {
            std::string use_key = d_refinement_criteria[k0];
            bool key_found = false;
            for (int k1 = 0; k1 < def_key_cnt; ++k1)
            {
                std::string def_key = ref_keys_defined[k1];
                if (def_key == use_key)
                {
                    key_found = true;
                }
            }

            if (!key_found)
            {
                TBOX_ERROR(d_object_name << ":\n"
                           << "  No input found for specified refine criteria: "
                           << d_refinement_criteria[k0] << "\n");
            }
        }
    } // refine db entry exists
    return;
}// getFromInput

void
AdvectHypPatchOps::getFromRestart()
{
    Pointer<Database> root_db =
        RestartManager::getManager()->getRootDatabase();

    Pointer<Database> db;

    if (root_db->isDatabase(d_object_name))
    {
        db = root_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart database corresponding to "
                   << d_object_name << " not found in restart file.");
    }

    int ver = db->getInteger("ADVECT_HYP_PATCH_OPS_VERSION");
    if (ver != ADVECT_HYP_PATCH_OPS_VERSION)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  Restart file version different than class version.");
    }

    db->getIntegerArray("d_ghosts", static_cast<int*>(d_ghosts), NDIM);
    if (d_ghosts != IntVector<NDIM>(CELLG))
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  Key data `d_ghosts' in restart file != CELLG.\n");
    }

    db->getIntegerArray("d_flux_ghosts", static_cast<int*>(d_flux_ghosts), NDIM);
    if (d_flux_ghosts != IntVector<NDIM>(FLUXG))
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  Key data `d_flux_ghosts' in restart file != FLUXG.\n");
    }

    if (db->keyExists("d_refinement_criteria"))
    {
        d_refinement_criteria = db->getStringArray("d_refinement_criteria");
    }
    for (int i = 0; i < d_refinement_criteria.getSize(); ++i)
    {
        if (d_refinement_criteria[i] == "QVAL_DEVIATION")
        {
            d_dev_tol = db->getDoubleArray("d_dev_tol");
            d_dev = db->getDoubleArray("d_dev");
            d_dev_time_max = db->getDoubleArray("d_dev_time_max");
            d_dev_time_min = db->getDoubleArray("d_dev_time_min");
        }
        else if (d_refinement_criteria[i] == "QVAL_GRADIENT")
        {
            d_grad_tol = db->getDoubleArray("d_grad_tol");
            d_grad_time_max = db->getDoubleArray("d_grad_time_max");
            d_grad_time_min = db->getDoubleArray("d_grad_time_min");
        }
    }
    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::AdvectHypPatchOps>;

//////////////////////////////////////////////////////////////////////////////
