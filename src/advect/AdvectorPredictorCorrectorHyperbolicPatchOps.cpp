// Filename: AdvectorPredictorCorrectorHyperbolicPatchOps.cpp
// Created on 12 Mar 2004 by Boyce Griffith
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
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/xfer/Geometry.h"
#include "SAMRAI/algs/HyperbolicLevelIntegrator.h"
#include "IBAMR_config.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/math/PatchFaceDataOpsReal.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/solv/RobinBcCoefStrategy.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "ibamr/AdvectorExplicitPredictorPatchOps.h"
#include "ibamr/AdvectorPredictorCorrectorHyperbolicPatchOps.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartExtrapPhysBdryOp.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/CartGridFunctionSet.h"
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/PIO.h"

#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define ADVECT_CONSDIFF_FC IBAMR_FC_FUNC_(advect_consdiff2d, ADVECT_CONSDIFF2D)
#define ADVECT_CONSDIFFWITHDIVSOURCE_FC IBAMR_FC_FUNC_(advect_consdiffwithdivsource2d, ADVECT_CONSDIFFWITHDIVSOURCE2D)
#define ADVECT_DETECTGRAD_FC IBAMR_FC_FUNC_(advect_detectgrad2d, ADVECT_DETECTGRAD2D)
#endif

#if (NDIM == 3)
#define ADVECT_CONSDIFF_FC IBAMR_FC_FUNC_(advect_consdiff3d, ADVECT_CONSDIFF3D)
#define ADVECT_CONSDIFFWITHDIVSOURCE_FC IBAMR_FC_FUNC_(advect_consdiffwithdivsource3d, ADVECT_CONSDIFFWITHDIVSOURCE3D)
#define ADVECT_DETECTGRAD_FC IBAMR_FC_FUNC_(advect_detectgrad3d, ADVECT_DETECTGRAD3D)
#endif

extern "C" {
void ADVECT_CONSDIFF_FC(const double*,
#if (NDIM == 2)
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
                        const double*,
                        const double*,
                        const double*,
#endif
                        double*);

void ADVECT_CONSDIFFWITHDIVSOURCE_FC(const double*,
                                     const double&,
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
                                     const int&,
                                     const int&,
                                     const double*,
                                     const double*,
                                     const double*,
                                     const double*,
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
                                     const int&,
                                     const int&,
                                     const int&,
                                     const double*,
                                     const double*,
                                     const double*,
                                     const double*,
                                     const double*,
                                     const double*,
                                     const double*,
                                     const double*,
                                     const double*,
#endif
                                     double*);

void ADVECT_DETECTGRAD_FC(
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
#endif
    const double*,
    const double&,
    const int&,
    const double*,
    int*,
    int*);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Values for cell tagging routines.
static const int TRUE_VAL = 1;
static const int FALSE_VAL = 0;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvectorPredictorCorrectorHyperbolicPatchOps::AdvectorPredictorCorrectorHyperbolicPatchOps(
    const std::string& object_name,
    boost::shared_ptr<Database> input_db,
    boost::shared_ptr<AdvectorExplicitPredictorPatchOps> explicit_predictor,
    boost::shared_ptr<CartesianGridGeometry> grid_geom,
    bool register_for_restart)
    : HyperbolicPatchStrategy(DIM), d_integrator(NULL), d_explicit_predictor(explicit_predictor), d_u_var(),
      d_u_is_div_free(), d_u_fcn(), d_compute_init_velocity(true), d_compute_half_velocity(true),
      d_compute_final_velocity(true), d_F_var(), d_F_fcn(), d_Q_var(), d_Q_u_map(), d_Q_F_map(), d_Q_difference_form(),
      d_Q_init(), d_Q_bc_coef(), d_overwrite_tags(true), d_object_name(object_name),
      d_registered_for_restart(register_for_restart), d_grid_geometry(grid_geom), d_visit_writer(NULL),
      d_extrap_bc_helper(), d_ghosts(DIM), d_flux_ghosts(DIM), d_extrap_type("CONSTANT"), d_refinement_criteria(),
      d_dev_tol(), d_dev(), d_dev_time_max(), d_dev_time_min(), d_grad_tol(), d_grad_time_max(), d_grad_time_min()
{
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(d_explicit_predictor);
    TBOX_ASSERT(grid_geom);

    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from given input/restart databases.
    bool is_from_restart = RestartManager::getManager()->isFromRestart();
    if (is_from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, is_from_restart);
    // Get number of ghost cells from the explicit predictor.
    d_ghosts = IntVector(DIM, explicit_predictor->getNumberCellGhosts());
    d_flux_ghosts = IntVector(DIM, explicit_predictor->getNumberFluxGhosts());
    return;
} // AdvectorPredictorCorrectorHyperbolicPatchOps

AdvectorPredictorCorrectorHyperbolicPatchOps::~AdvectorPredictorCorrectorHyperbolicPatchOps()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
    return;
} // ~AdvectorPredictorCorrectorHyperbolicPatchOps

const std::string& AdvectorPredictorCorrectorHyperbolicPatchOps::getName() const
{
    return d_object_name;
} // getName

void
AdvectorPredictorCorrectorHyperbolicPatchOps::registerVisItDataWriter(boost::shared_ptr<VisItDataWriter> visit_writer)
{
    TBOX_ASSERT(visit_writer);
    d_visit_writer = visit_writer;
    return;
} // registerVisItDataWriter

void
AdvectorPredictorCorrectorHyperbolicPatchOps::registerAdvectionVelocity(boost::shared_ptr<FaceVariable<double> > u_var)
{
    TBOX_ASSERT(u_var);
    d_u_var.insert(u_var);
    d_u_is_div_free[u_var] = true;
    return;
} // registerAdvectionVelocity

void AdvectorPredictorCorrectorHyperbolicPatchOps::setAdvectionVelocityIsDivergenceFree(
    boost::shared_ptr<FaceVariable<double> > u_var,
    const bool is_div_free)
{
    TBOX_ASSERT(d_u_var.find(u_var) != d_u_var.end());
    d_u_is_div_free[u_var] = is_div_free;
    return;
} // setAdvectionVelocityIsDivergenceFree

void AdvectorPredictorCorrectorHyperbolicPatchOps::setAdvectionVelocityFunction(
    boost::shared_ptr<FaceVariable<double> > u_var,
    boost::shared_ptr<CartGridFunction> u_fcn)
{
    TBOX_ASSERT(d_u_var.find(u_var) != d_u_var.end());
    d_u_fcn[u_var] = u_fcn;
    return;
} // setAdvectionVelocityFunction

void AdvectorPredictorCorrectorHyperbolicPatchOps::registerSourceTerm(boost::shared_ptr<CellVariable<double> > F_var)
{
    TBOX_ASSERT(F_var);
    d_F_var.insert(F_var);
    return;
} // registerSourceTerm

void AdvectorPredictorCorrectorHyperbolicPatchOps::setSourceTermFunction(boost::shared_ptr<CellVariable<double> > F_var,
                                                                         boost::shared_ptr<CartGridFunction> F_fcn)
{
    TBOX_ASSERT(d_F_var.find(F_var) != d_F_var.end());
    if (d_F_fcn[F_var])
    {
        const std::string& F_var_name = F_var->getName();
        boost::shared_ptr<CartGridFunctionSet> p_F_fcn = d_F_fcn[F_var];
        if (!p_F_fcn)
        {
            pout << d_object_name << "::setSourceTermFunction(): WARNING:\n"
                 << "  source term function for source term variable " << F_var_name << " has already been set.\n"
                 << "  functions will be evaluated in the order in which they were registered "
                    "with "
                    "the solver\n"
                 << "  when evaluating the source term value.\n";
            p_F_fcn =
                boost::make_shared<CartGridFunctionSet>(d_object_name + "::" + F_var_name + "::source_function_set");
            p_F_fcn->addFunction(d_F_fcn[F_var]);
        }
        p_F_fcn->addFunction(F_fcn);
    }
    else
    {
        d_F_fcn[F_var] = F_fcn;
    }
    return;
} // setSourceTermFunction

void AdvectorPredictorCorrectorHyperbolicPatchOps::registerTransportedQuantity(
    boost::shared_ptr<CellVariable<double> > Q_var)
{
    TBOX_ASSERT(Q_var);
    d_Q_var.insert(Q_var);
    d_Q_difference_form[Q_var] = CONSERVATIVE;
    return;
} // registerTransportedQuantity

void AdvectorPredictorCorrectorHyperbolicPatchOps::setAdvectionVelocity(boost::shared_ptr<CellVariable<double> > Q_var,
                                                                        boost::shared_ptr<FaceVariable<double> > u_var)
{
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    TBOX_ASSERT(d_u_var.find(u_var) != d_u_var.end());
    d_Q_u_map[Q_var] = u_var;
    return;
} // setAdvectionVelocity

void AdvectorPredictorCorrectorHyperbolicPatchOps::setSourceTerm(boost::shared_ptr<CellVariable<double> > Q_var,
                                                                 boost::shared_ptr<CellVariable<double> > F_var)
{
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    TBOX_ASSERT(d_F_var.find(F_var) != d_F_var.end());
    d_Q_F_map[Q_var] = F_var;
    return;
} // setSourceTerm

void AdvectorPredictorCorrectorHyperbolicPatchOps::setConvectiveDifferencingType(
    boost::shared_ptr<CellVariable<double> > Q_var,
    const ConvectiveDifferencingType difference_form)
{
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    d_Q_difference_form[Q_var] = difference_form;
    return;
} // setConvectiveDifferencingType

void AdvectorPredictorCorrectorHyperbolicPatchOps::setInitialConditions(boost::shared_ptr<CellVariable<double> > Q_var,
                                                                        boost::shared_ptr<CartGridFunction> Q_init)
{
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    d_Q_init[Q_var] = Q_init;
    return;
} // setInitialConditions

void AdvectorPredictorCorrectorHyperbolicPatchOps::setPhysicalBcCoefs(boost::shared_ptr<CellVariable<double> > Q_var,
                                                                      RobinBcCoefStrategy* Q_bc_coef)
{
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    const unsigned int Q_depth = Q_var->getDepth();
    TBOX_ASSERT(Q_depth == 1);
    d_Q_bc_coef[Q_var] = std::vector<RobinBcCoefStrategy*>(1, Q_bc_coef);
    return;
} // setPhysicalBcCoefs

void AdvectorPredictorCorrectorHyperbolicPatchOps::setPhysicalBcCoefs(boost::shared_ptr<CellVariable<double> > Q_var,
                                                                      std::vector<RobinBcCoefStrategy*> Q_bc_coef)
{
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    const unsigned int Q_depth = Q_var->getDepth();
    TBOX_ASSERT(Q_depth == Q_bc_coef.size());
    d_Q_bc_coef[Q_var] = Q_bc_coef;
    return;
} // setPhysicalBcCoefs

void AdvectorPredictorCorrectorHyperbolicPatchOps::registerModelVariables(HyperbolicLevelIntegrator* integrator)
{
    TBOX_ASSERT(integrator);
    d_integrator = integrator;

    for (auto cit = d_u_var.begin();
         cit != d_u_var.end();
         ++cit)
    {
        auto u_var = *cit;
        d_integrator->registerVariable(u_var,
                                       d_ghosts,
                                       HyperbolicLevelIntegrator::TIME_DEP,
                                       d_grid_geometry,
                                       "CONSERVATIVE_COARSEN",
                                       "CONSERVATIVE_LINEAR_REFINE");
    }

    for (auto cit = d_F_var.begin();
         cit != d_F_var.end();
         ++cit)
    {
        auto F_var = *cit;
        d_integrator->registerVariable(F_var,
                                       d_ghosts,
                                       HyperbolicLevelIntegrator::TIME_DEP,
                                       d_grid_geometry,
                                       "CONSERVATIVE_COARSEN",
                                       "CONSERVATIVE_LINEAR_REFINE");
    }

    for (auto cit = d_Q_var.begin();
         cit != d_Q_var.end();
         ++cit)
    {
        auto Q_var = *cit;
        boost::shared_ptr<CellDataFactory<double> > Q_factory = Q_var->getPatchDataFactory();
        d_integrator->registerVariable(Q_var,
                                       d_ghosts,
                                       HyperbolicLevelIntegrator::TIME_DEP,
                                       d_grid_geometry,
                                       "CONSERVATIVE_COARSEN",
                                       "CONSERVATIVE_LINEAR_REFINE");

        if (d_visit_writer)
        {
            const int Q_idx =
                VariableDatabase::getDatabase()->mapVariableAndContextToIndex(Q_var, d_integrator->getPlotContext());
            const int depth = Q_var->getDepth();
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
                    d_visit_writer->registerPlotQuantity(Q_var->getName() + "_" + stream.str(), "SCALAR", Q_idx, d);
                }
            }
        }

        const bool conservation_form = d_Q_difference_form[Q_var] == CONSERVATIVE;
        if (conservation_form)
        {
            d_flux_integral_var[Q_var] = boost::make_shared<FaceVariable<double> >(
                DIM, d_object_name + "::" + Q_var->getName() + " advective flux time integral", Q_var->getDepth());
            d_integrator->registerVariable(d_flux_integral_var[Q_var],
                                           d_flux_ghosts,
                                           HyperbolicLevelIntegrator::FLUX,
                                           d_grid_geometry,
                                           "CONSERVATIVE_COARSEN",
                                           "NO_REFINE");
        }

        auto u_var = d_Q_u_map[Q_var];
        const bool u_is_div_free = d_u_is_div_free[u_var];
        if (!conservation_form || !u_is_div_free)
        {
            d_q_integral_var[Q_var] = boost::make_shared<FaceVariable<double> >(
                DIM, d_object_name + "::" + Q_var->getName() + " time integral", Q_var->getDepth());
            d_integrator->registerVariable(d_q_integral_var[Q_var],
                                           d_flux_ghosts,
                                           HyperbolicLevelIntegrator::FLUX,
                                           d_grid_geometry,
                                           "CONSERVATIVE_COARSEN",
                                           "NO_REFINE");
            if (u_var && !d_u_integral_var[u_var])
            {
                d_u_integral_var[u_var] = boost::make_shared<FaceVariable<double> >(
                    DIM, d_object_name + "::" + u_var->getName() + " time integral");
                d_integrator->registerVariable(d_u_integral_var[u_var],
                                               d_flux_ghosts,
                                               HyperbolicLevelIntegrator::FLUX,
                                               d_grid_geometry,
                                               "CONSERVATIVE_COARSEN",
                                               "NO_REFINE");
            }
        }
    }
    return;
} // registerModelVariables

void AdvectorPredictorCorrectorHyperbolicPatchOps::initializeDataOnPatch(Patch& patch,
                                                                         const double data_time,
                                                                         const bool initial_time)
{
    if (initial_time)
    {
        VariableDatabase* var_db = VariableDatabase::getDatabase();
        for (auto cit = d_u_var.begin();
             cit != d_u_var.end();
             ++cit)
        {
            auto u_var = *cit;
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, getDataContext());
            if (d_u_fcn[u_var])
            {
                d_u_fcn[u_var]->setDataOnPatch(
                    u_idx, u_var, boost::shared_ptr<Patch>(&patch, NullDeleter()), data_time, initial_time);
            }
            else
            {
                boost::shared_ptr<FaceData<double> > u_data = patch.getPatchData(u_idx);
                u_data->fillAll(0.0);
            }
        }

        for (auto cit = d_F_var.begin();
             cit != d_F_var.end();
             ++cit)
        {
            auto F_var = *cit;
            const int F_idx = var_db->mapVariableAndContextToIndex(F_var, getDataContext());
            if (d_F_fcn[F_var])
            {
                d_F_fcn[F_var]->setDataOnPatch(
                    F_idx, F_var, boost::shared_ptr<Patch>(&patch, NullDeleter()), data_time, initial_time);
            }
            else
            {
                boost::shared_ptr<CellData<double> > F_data = patch.getPatchData(F_idx);
                F_data->fillAll(0.0);
            }
        }

        for (auto cit = d_Q_var.begin();
             cit != d_Q_var.end();
             ++cit)
        {
            auto Q_var = *cit;
            const int Q_idx = var_db->mapVariableAndContextToIndex(Q_var, getDataContext());
            if (d_Q_init[Q_var])
            {
                d_Q_init[Q_var]->setDataOnPatch(
                    Q_idx, Q_var, boost::shared_ptr<Patch>(&patch, NullDeleter()), data_time, initial_time);
            }
            else
            {
                boost::shared_ptr<CellData<double> > Q_data = patch.getPatchData(Q_var, getDataContext());
                Q_data->fillAll(0.0);
            }
        }
    }
    return;
} // initializeDataOnPatch

double AdvectorPredictorCorrectorHyperbolicPatchOps::computeStableDtOnPatch(Patch& patch,
                                                                            const bool /*initial_time*/,
                                                                            const double /*dt_time*/)
{
    double stable_dt = std::numeric_limits<double>::max();
    for (auto cit = d_u_var.begin();
         cit != d_u_var.end();
         ++cit)
    {
        auto u_var = *cit;
        boost::shared_ptr<FaceData<double> > u_data = patch.getPatchData(u_var, getDataContext());
        stable_dt = std::min(stable_dt, d_explicit_predictor->computeStableDtOnPatch(*u_data, patch));
    }
    return stable_dt;
} // computeStableDtOnPatch

void
AdvectorPredictorCorrectorHyperbolicPatchOps::computeFluxesOnPatch(Patch& patch, const double time, const double dt)
{
    auto pgeom = patch.getPatchGeometry();
    const Box& patch_box = patch.getBox();

    PatchFaceDataOpsReal<double> patch_fc_data_ops;

    for (auto cit = d_Q_var.begin();
         cit != d_Q_var.end();
         ++cit)
    {
        auto Q_var = *cit;
        auto u_var = d_Q_u_map[Q_var];
        boost::shared_ptr<FaceData<double> > q_integral_data = getQIntegralData(Q_var, patch, getDataContext());
        if (!u_var)
        {
            q_integral_data->fillAll(0.0);
            continue;
        }

        // Predict time- and face-centered values.
        boost::shared_ptr<CellData<double> > Q_data = patch.getPatchData(Q_var, getDataContext());
        boost::shared_ptr<FaceData<double> > u_data = patch.getPatchData(u_var, getDataContext());
        auto F_var = d_Q_F_map[Q_var];
        if (F_var)
        {
            boost::shared_ptr<CellData<double> > F_data = patch.getPatchData(F_var, getDataContext());
            d_explicit_predictor->predictValueWithSourceTerm(*q_integral_data, *u_data, *Q_data, *F_data, patch, dt);
        }
        else
        {
            d_explicit_predictor->predictValue(*q_integral_data, *u_data, *Q_data, patch, dt);
        }
    }

    // Set physical boundary conditions for the face-centered predicted values
    // at inflow boundaries only.
    if (pgeom->getTouchesRegularBoundary())
    {
        setInflowBoundaryConditions(patch, time + 0.5 * dt);
    }

    // Update the advection velocity.
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    if (d_compute_half_velocity)
    {
        for (auto cit = d_u_var.begin();
             cit != d_u_var.end();
             ++cit)
        {
            auto u_var = *cit;
            if (d_u_fcn[u_var] && d_u_fcn[u_var]->isTimeDependent())
            {
                const int u_idx = var_db->mapVariableAndContextToIndex(u_var, getDataContext());
                d_u_fcn[u_var]->setDataOnPatch(
                    u_idx, u_var, boost::shared_ptr<Patch>(&patch, NullDeleter()), time + 0.5 * dt);
            }
        }
    }

    // Compute fluxes and other face-centered quantities.
    for (auto cit = d_Q_var.begin();
         cit != d_Q_var.end();
         ++cit)
    {
        auto Q_var = *cit;
        auto u_var = d_Q_u_map[Q_var];

        if (!u_var) continue;

        boost::shared_ptr<FaceData<double> > u_data = patch.getPatchData(u_var, getDataContext());
        const bool conservation_form = d_Q_difference_form[Q_var] == CONSERVATIVE;
        const bool u_is_div_free = d_u_is_div_free[u_var];

        if (conservation_form)
        {
            boost::shared_ptr<FaceData<double> > flux_integral_data =
                getFluxIntegralData(Q_var, patch, getDataContext());
            boost::shared_ptr<FaceData<double> > q_integral_data = getQIntegralData(Q_var, patch, getDataContext());
            d_explicit_predictor->computeFlux(*flux_integral_data, *u_data, *q_integral_data, patch, dt);
        }

        if (!conservation_form || !u_is_div_free)
        {
            boost::shared_ptr<FaceData<double> > q_integral_data = getQIntegralData(Q_var, patch, getDataContext());
            patch_fc_data_ops.scale(q_integral_data, dt, q_integral_data, patch_box);

            boost::shared_ptr<FaceData<double> > u_integral_data = getUIntegralData(Q_var, patch, getDataContext());
            patch_fc_data_ops.scale(u_integral_data, dt, u_data, patch_box);
        }
    }
    return;
} // computeFluxesOnPatch

void AdvectorPredictorCorrectorHyperbolicPatchOps::conservativeDifferenceOnPatch(Patch& patch,
                                                                                 const double /*time*/,
                                                                                 const double dt,
                                                                                 bool /*at_synchronization*/)
{
    const Box& patch_box = patch.getBox();
    const Index& ilower = patch_box.lower();
    const Index& iupper = patch_box.upper();

    auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch.getPatchGeometry());
    const double* const dx = pgeom->getDx();

    for (auto cit = d_Q_var.begin();
         cit != d_Q_var.end();
         ++cit)
    {
        auto Q_var = *cit;
        auto u_var = d_Q_u_map[Q_var];

        if (!u_var) continue;

        boost::shared_ptr<CellData<double> > Q_data = patch.getPatchData(Q_var, getDataContext());
        boost::shared_ptr<FaceData<double> > flux_integral_data = getFluxIntegralData(Q_var, patch, getDataContext());
        boost::shared_ptr<FaceData<double> > q_integral_data = getQIntegralData(Q_var, patch, getDataContext());
        boost::shared_ptr<FaceData<double> > u_integral_data = getUIntegralData(Q_var, patch, getDataContext());

        const IntVector no_ghosts = IntVector::getZero(DIM);
        const IntVector& Q_data_ghost_cells = Q_data->getGhostCellWidth();
        const IntVector& flux_integral_data_ghost_cells =
            (flux_integral_data ? flux_integral_data->getGhostCellWidth() : no_ghosts);
        const IntVector& q_integral_data_ghost_cells =
            (q_integral_data ? q_integral_data->getGhostCellWidth() : no_ghosts);
        const IntVector& u_integral_data_ghost_cells =
            (u_integral_data ? u_integral_data->getGhostCellWidth() : no_ghosts);

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
                    ADVECT_CONSDIFF_FC(dx,
                                       ilower(0),
                                       iupper(0),
                                       ilower(1),
                                       iupper(1),
                                       flux_integral_data_ghost_cells(0),
                                       flux_integral_data_ghost_cells(1),
                                       Q_data_ghost_cells(0),
                                       Q_data_ghost_cells(1),
                                       flux_integral_data->getPointer(0, depth),
                                       flux_integral_data->getPointer(1, depth),
                                       Q_data->getPointer(depth));
#endif
#if (NDIM == 3)
                    ADVECT_CONSDIFF_FC(dx,
                                       ilower(0),
                                       iupper(0),
                                       ilower(1),
                                       iupper(1),
                                       ilower(2),
                                       iupper(2),
                                       flux_integral_data_ghost_cells(0),
                                       flux_integral_data_ghost_cells(1),
                                       flux_integral_data_ghost_cells(2),
                                       Q_data_ghost_cells(0),
                                       Q_data_ghost_cells(1),
                                       Q_data_ghost_cells(2),
                                       flux_integral_data->getPointer(0, depth),
                                       flux_integral_data->getPointer(1, depth),
                                       flux_integral_data->getPointer(2, depth),
                                       Q_data->getPointer(depth));
#endif
                }
                else
                {
#if (NDIM == 2)
                    ADVECT_CONSDIFFWITHDIVSOURCE_FC(dx,
                                                    dt,
                                                    ilower(0),
                                                    iupper(0),
                                                    ilower(1),
                                                    iupper(1),
                                                    flux_integral_data_ghost_cells(0),
                                                    flux_integral_data_ghost_cells(1),
                                                    q_integral_data_ghost_cells(0),
                                                    q_integral_data_ghost_cells(1),
                                                    u_integral_data_ghost_cells(0),
                                                    u_integral_data_ghost_cells(1),
                                                    Q_data_ghost_cells(0),
                                                    Q_data_ghost_cells(1),
                                                    flux_integral_data->getPointer(0, depth),
                                                    flux_integral_data->getPointer(1, depth),
                                                    q_integral_data->getPointer(0, depth),
                                                    q_integral_data->getPointer(1, depth),
                                                    u_integral_data->getPointer(0),
                                                    u_integral_data->getPointer(1),
                                                    Q_data->getPointer(depth));
#endif
#if (NDIM == 3)
                    ADVECT_CONSDIFFWITHDIVSOURCE_FC(dx,
                                                    dt,
                                                    ilower(0),
                                                    iupper(0),
                                                    ilower(1),
                                                    iupper(1),
                                                    ilower(2),
                                                    iupper(2),
                                                    flux_integral_data_ghost_cells(0),
                                                    flux_integral_data_ghost_cells(1),
                                                    flux_integral_data_ghost_cells(2),
                                                    q_integral_data_ghost_cells(0),
                                                    q_integral_data_ghost_cells(1),
                                                    q_integral_data_ghost_cells(2),
                                                    u_integral_data_ghost_cells(0),
                                                    u_integral_data_ghost_cells(1),
                                                    u_integral_data_ghost_cells(2),
                                                    Q_data_ghost_cells(0),
                                                    Q_data_ghost_cells(1),
                                                    Q_data_ghost_cells(2),
                                                    flux_integral_data->getPointer(0, depth),
                                                    flux_integral_data->getPointer(1, depth),
                                                    flux_integral_data->getPointer(2, depth),
                                                    q_integral_data->getPointer(0, depth),
                                                    q_integral_data->getPointer(1, depth),
                                                    q_integral_data->getPointer(2, depth),
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
            CellData<double> N_data(patch_box, Q_data->getDepth(), IntVector::getZero(DIM));
            d_explicit_predictor->computeAdvectiveDerivative(N_data, *u_integral_data, *q_integral_data, patch);
            PatchCellDataOpsReal<double> patch_cc_data_ops;
            patch_cc_data_ops.axpy(
                Q_data, -1.0 / dt, boost::shared_ptr<CellData<double> >(&N_data, NullDeleter()), Q_data, patch_box);
            break;
        }
        default:
        {
            TBOX_ERROR(
                "AdvectorPredictorCorrectorHyperbolicPatchOps::conservativeDifferenceOnPatch()"
                ":\n"
                << "  unsupported differencing form: "
                << enum_to_string<ConvectiveDifferencingType>(d_Q_difference_form[Q_var]) << " \n"
                << "  valid choices are: ADVECTIVE, CONSERVATIVE\n");
        }
        }
    }
    return;
} // conservativeDifferenceOnPatch

void
AdvectorPredictorCorrectorHyperbolicPatchOps::preprocessAdvanceLevelState(const boost::shared_ptr<PatchLevel>& level,
                                                                          double current_time,
                                                                          double /*dt*/,
                                                                          bool /*first_step*/,
                                                                          bool /*last_step*/,
                                                                          bool /*regrid_advance*/)
{
    VariableDatabase* var_db = VariableDatabase::getDatabase();

    // Update the source term.
    for (auto cit = d_F_var.begin();
         cit != d_F_var.end();
         ++cit)
    {
        auto F_var = *cit;
        if (d_F_fcn[F_var] && d_F_fcn[F_var]->isTimeDependent())
        {
            const int F_idx = var_db->mapVariableAndContextToIndex(F_var, d_integrator->getScratchContext());
            d_F_fcn[F_var]->setDataOnPatchLevel(F_idx, F_var, level, current_time);
        }
    }

    if (!d_compute_init_velocity) return;

    // Update the advection velocity.
    for (auto cit = d_u_var.begin();
         cit != d_u_var.end();
         ++cit)
    {
        auto u_var = *cit;
        if (d_u_fcn[u_var] && d_u_fcn[u_var]->isTimeDependent())
        {
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, d_integrator->getScratchContext());
            d_u_fcn[u_var]->setDataOnPatchLevel(u_idx, u_var, level, current_time);
        }
    }
    return;
} // preprocessAdvanceLevelState

void
AdvectorPredictorCorrectorHyperbolicPatchOps::postprocessAdvanceLevelState(const boost::shared_ptr<PatchLevel>& level,
                                                                           double current_time,
                                                                           double dt,
                                                                           bool /*first_step*/,
                                                                           bool /*last_step*/,
                                                                           bool /*regrid_advance*/)
{
    VariableDatabase* var_db = VariableDatabase::getDatabase();

    PatchCellDataOpsReal<double> patch_cc_data_ops;

    auto new_context = d_integrator->getNewContext();
    auto scratch_context = d_integrator->getScratchContext();

    // Update the values of any time-dependent source terms and add the values
    // of all source terms to the advected quantities.
    for (auto cit = d_Q_var.begin();
         cit != d_Q_var.end();
         ++cit)
    {
        auto Q_var = *cit;
        auto F_var = d_Q_F_map[Q_var];
        if (!F_var) continue;
        auto F_fcn = d_F_fcn[F_var];
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch =*p;
            const Box& patch_box = patch->getBox();

            boost::shared_ptr<CellData<double> > Q_data = patch->getPatchData(Q_var, new_context);
            if (F_fcn)
            {
                const int F_scratch_idx = var_db->mapVariableAndContextToIndex(F_var, scratch_context);
                const int F_new_idx = var_db->mapVariableAndContextToIndex(F_var, new_context);
                boost::shared_ptr<CellData<double> > F_scratch_data = patch->getPatchData(F_scratch_idx);
                boost::shared_ptr<CellData<double> > F_new_data = patch->getPatchData(F_new_idx);
                F_fcn->setDataOnPatchLevel(F_new_idx, F_var, level, current_time + dt);
                patch_cc_data_ops.axpy(Q_data, 0.5 * dt, F_scratch_data, Q_data, patch_box);
                patch_cc_data_ops.axpy(Q_data, 0.5 * dt, F_new_data, Q_data, patch_box);
            }
            else
            {
                const int F_scratch_idx = var_db->mapVariableAndContextToIndex(F_var, scratch_context);
                boost::shared_ptr<CellData<double> > F_scratch_data = patch->getPatchData(F_scratch_idx);
                patch_cc_data_ops.axpy(Q_data, dt, F_scratch_data, Q_data, patch_box);
            }
        }
    }

    if (!d_compute_final_velocity) return;

    // Update the advection velocity.
    for (auto cit = d_u_var.begin();
         cit != d_u_var.end();
         ++cit)
    {
        auto u_var = *cit;
        if (d_u_fcn[u_var] && d_u_fcn[u_var]->isTimeDependent())
        {
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, d_integrator->getNewContext());
            d_u_fcn[u_var]->setDataOnPatchLevel(u_idx, u_var, level, current_time + dt);
        }
    }
    return;
} // postprocessAdvanceLevelState

void
AdvectorPredictorCorrectorHyperbolicPatchOps::tagGradientDetectorCells(Patch& patch,
                                                                       const double regrid_time,
                                                                       const bool /*initial_error*/,
                                                                       const int tag_indx,
                                                                       const bool /*uses_richardson_extrapolation_too*/)
{
    const int error_level_number = patch.getPatchLevelNumber();

    auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch.getPatchGeometry());
    const double* const dx = pgeom->getDx();

    const Box& patch_box = patch.getBox();
    const Index& ilower = patch.getBox().lower();
    const Index& iupper = patch.getBox().upper();

    boost::shared_ptr<CellData<int> > tags = patch.getPatchData(tag_indx);

    const int not_refine_tag_val = FALSE_VAL;
    const int refine_tag_val = TRUE_VAL;

    // Create a set of temporary tags and set to untagged value.
    CellData<int> temp_tags(patch_box, 1, d_ghosts);
    temp_tags.fillAll(not_refine_tag_val);

    // Possible tagging criteria include: QVAL_DEVIATION, QVAL_GRADIENT.  The
    // criteria are specified over a time interval.
    //
    // Loop over criteria provided and check to make sure we are in the
    // specified time interval.  If so, apply appropriate tagging for the level.
    for (int ncrit = 0; ncrit < d_refinement_criteria.getSize(); ++ncrit)
    {
        std::string ref = d_refinement_criteria[ncrit];
        const IntVector& tagghost = tags->getGhostCellWidth();

        int size = 0;
        double tol = 0.0;
        bool time_allowed = false;

        if (ref == "QVAL_DEVIATION")
        {
            size = d_dev_tol.getSize();
            tol = (error_level_number < size ? d_dev_tol[error_level_number] : d_dev_tol[size - 1]);
            size = d_dev.getSize();
            double dev = (error_level_number < size ? d_dev[error_level_number] : d_dev[size - 1]);
            size = d_dev_time_min.getSize();
            double time_min =
                (error_level_number < size ? d_dev_time_min[error_level_number] : d_dev_time_min[size - 1]);
            size = d_dev_time_max.getSize();
            double time_max =
                (error_level_number < size ? d_dev_time_max[error_level_number] : d_dev_time_max[size - 1]);
            time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);

            if (time_allowed)
            {
                // Check for tags that have already been set in a previous step.
                for (auto cit = d_Q_var.begin();
                     cit != d_Q_var.end();
                     ++cit)
                {
                    auto Q_var = *cit;
                    boost::shared_ptr<CellData<double> > Q_data = patch.getPatchData(Q_var, getDataContext());
                    for (int depth = 0; depth < Q_data->getDepth(); ++depth)
                    {
                        for (CellIterator ic(patch_box); ic; ic++)
                        {
                            double locden = tol;
                            int tag_val = (*tags)(ic(), 0);
                            if (tag_val)
                            {
                                locden *= 0.75;
                            }
                            if (std::abs((*Q_data)(ic(), depth) - dev) > locden)
                            {
                                temp_tags(ic(), 0) = refine_tag_val;
                            }
                        }
                    }
                }
            }
        }

        if (ref == "QVAL_GRADIENT")
        {
            size = d_grad_tol.getSize();
            tol = (error_level_number < size ? d_grad_tol[error_level_number] : d_grad_tol[size - 1]);
            size = d_grad_time_min.getSize();
            double time_min =
                (error_level_number < size ? d_grad_time_min[error_level_number] : d_grad_time_min[size - 1]);
            size = d_grad_time_max.getSize();
            double time_max =
                (error_level_number < size ? d_grad_time_max[error_level_number] : d_grad_time_max[size - 1]);
            time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);

            if (time_allowed)
            {
                for (auto cit = d_Q_var.begin();
                     cit != d_Q_var.end();
                     ++cit)
                {
                    auto Q_var = *cit;
                    boost::shared_ptr<CellData<double> > Q_data = patch.getPatchData(Q_var, getDataContext());
                    const IntVector& Q_ghost = Q_data->getGhostCellWidth();
                    for (int depth = 0; depth < Q_data->getDepth(); ++depth)
                    {
                        ADVECT_DETECTGRAD_FC(
#if (NDIM == 2)
                            ilower(0),
                            iupper(0),
                            ilower(1),
                            iupper(1),
                            Q_ghost(0),
                            tagghost(0),
                            d_ghosts(0),
                            Q_ghost(1),
                            tagghost(1),
                            d_ghosts(1),
#endif
#if (NDIM == 3)
                            ilower(0),
                            iupper(0),
                            ilower(1),
                            iupper(1),
                            ilower(2),
                            iupper(2),
                            Q_ghost(0),
                            tagghost(0),
                            d_ghosts(0),
                            Q_ghost(1),
                            tagghost(1),
                            d_ghosts(1),
                            Q_ghost(2),
                            tagghost(2),
                            d_ghosts(2),
#endif
                            dx,
                            tol,
                            refine_tag_val,
                            Q_data->getPointer(depth),
                            tags->getPointer(),
                            temp_tags.getPointer());
                    }
                }
            }
        }
    } // loop over criteria

    // Update tags.
    if (d_overwrite_tags)
    {
        for (CellIterator ic(patch_box); ic; ic++)
        {
            (*tags)(ic(), 0) = temp_tags(ic(), 0);
        }
    }
    else
    {
        for (CellIterator ic(patch_box); ic; ic++)
        {
            (*tags)(ic(), 0) = (*tags)(ic(), 0) || temp_tags(ic(), 0);
        }
    }
    return;
} // tagGradientDetectorCells

void AdvectorPredictorCorrectorHyperbolicPatchOps::setPhysicalBoundaryConditions(Patch& patch,
                                                                                 const double fill_time,
                                                                                 const IntVector& ghost_width_to_fill)
{
    VariableDatabase* var_db = VariableDatabase::getDatabase();

    // Extrapolate the interior data to set the ghost cell values for the state
    // variables and for any forcing terms.
    ComponentSelector u_patch_data_indices;
    for (auto cit = d_u_var.begin();
         cit != d_u_var.end();
         ++cit)
    {
        auto u_var = *cit;
        const int u_data_idx = var_db->mapVariableAndContextToIndex(u_var, d_integrator->getScratchContext());
        u_patch_data_indices.setFlag(u_data_idx);
    }
    d_extrap_bc_helper.setExtrapolationType(d_extrap_type);
    d_extrap_bc_helper.setPatchDataIndices(u_patch_data_indices);
    d_extrap_bc_helper.setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);

    ComponentSelector patch_data_indices;
    for (auto cit = d_F_var.begin();
         cit != d_F_var.end();
         ++cit)
    {
        auto F_var = *cit;
        const int F_data_idx = var_db->mapVariableAndContextToIndex(F_var, d_integrator->getScratchContext());
        patch_data_indices.setFlag(F_data_idx);
    }
    for (auto cit = d_Q_var.begin();
         cit != d_Q_var.end();
         ++cit)
    {
        auto Q_var = *cit;
        const int Q_data_idx = var_db->mapVariableAndContextToIndex(Q_var, d_integrator->getScratchContext());
        patch_data_indices.setFlag(Q_data_idx);
    }
    d_extrap_bc_helper.setExtrapolationType(d_extrap_type);
    d_extrap_bc_helper.setPatchDataIndices(patch_data_indices);
    d_extrap_bc_helper.setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);
    return;
} // setPhysicalBoundaryConditions

void AdvectorPredictorCorrectorHyperbolicPatchOps::putToRestart(const boost::shared_ptr<Database>& db) const
{
    TBOX_ASSERT(db);
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
} // putToRestart

/////////////////////////////// PROTECTED ////////////////////////////////////

boost::shared_ptr<FaceData<double> >
AdvectorPredictorCorrectorHyperbolicPatchOps::getFluxIntegralData(boost::shared_ptr<CellVariable<double> > Q_var,
                                                                  Patch& patch,
                                                                  boost::shared_ptr<VariableContext> context)
{
    TBOX_ASSERT(Q_var);
    if (d_Q_difference_form[Q_var] == CONSERVATIVE)
    {
        return patch.getPatchData(d_flux_integral_var[Q_var], context);
    }
    else
    {
        return NULL;
    }
} // getFluxIntegralData

boost::shared_ptr<FaceData<double> >
AdvectorPredictorCorrectorHyperbolicPatchOps::getQIntegralData(boost::shared_ptr<CellVariable<double> > Q_var,
                                                               Patch& patch,
                                                               boost::shared_ptr<VariableContext> context)
{
    TBOX_ASSERT(Q_var);
    if (!d_u_is_div_free[d_Q_u_map[Q_var]] || d_Q_difference_form[Q_var] != CONSERVATIVE)
    {
        return patch.getPatchData(d_q_integral_var[Q_var], context);
    }
    else
    {
        return getFluxIntegralData(Q_var, patch, context);
    }
} // getQIntegralData

boost::shared_ptr<FaceData<double> >
AdvectorPredictorCorrectorHyperbolicPatchOps::getUIntegralData(boost::shared_ptr<CellVariable<double> > Q_var,
                                                               Patch& patch,
                                                               boost::shared_ptr<VariableContext> context)
{
    TBOX_ASSERT(Q_var);
    auto u_var = d_Q_u_map[Q_var];
    if (u_var && (!d_u_is_div_free[u_var] || d_Q_difference_form[Q_var] != CONSERVATIVE))
    {
        return patch.getPatchData(d_u_integral_var[u_var], context);
    }
    else
    {
        return NULL;
    }
} // getUIntegralData

/////////////////////////////// PRIVATE //////////////////////////////////////

void AdvectorPredictorCorrectorHyperbolicPatchOps::setInflowBoundaryConditions(Patch& patch, const double fill_time)
{
    VariableDatabase* var_db = VariableDatabase::getDatabase();

    auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch.getPatchGeometry());

    // There is nothing to do if the patch does not touch a regular (physical)
    // boundary.
    if (!pgeom->getTouchesRegularBoundary()) return;

    // Compute the co-dimension one boundary boxes.
    const std::vector<BoundaryBox> physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);

    // There is nothing to do if the patch does not have any co-dimension one
    // boundary boxes.
    if (physical_codim1_boxes.size() == 0) return;

    // Loop over the boundary fill boxes and set boundary conditions at inflow
    // boundaries only.
    const Box& patch_box = patch.getBox();
    const double* const dx = pgeom->getDx();
    for (auto cit = d_Q_var.begin();
         cit != d_Q_var.end();
         ++cit)
    {
        auto Q_var = *cit;
        auto u_var = d_Q_u_map[Q_var];

        if (!u_var) continue;

        const int Q_data_idx = var_db->mapVariableAndContextToIndex(Q_var, d_integrator->getScratchContext());

        boost::shared_ptr<FaceData<double> > u_data = patch.getPatchData(u_var, getDataContext());
        boost::shared_ptr<FaceData<double> > q_integral_data = getQIntegralData(Q_var, patch, getDataContext());

        // Setup any extended Robin BC coef objects.
        for (int depth = 0; depth < q_integral_data->getDepth(); ++depth)
        {
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_Q_bc_coef[Q_var][depth]);
            if (extended_bc_coef)
            {
                extended_bc_coef->setTargetPatchDataIndex(Q_data_idx);
                extended_bc_coef->setHomogeneousBc(false);
            }
        }

        // Set the boundary conditions.
        for (int n = 0; n < physical_codim1_boxes.size(); ++n)
        {
            const BoundaryBox& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = location_index % 2 == 0;

            static const IntVector gcw_to_fill = IntVector::getOne(DIM);
            const Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
            const BoundaryBox trimmed_bdry_box(
                bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), bdry_box.getLocationIndex());
            const Box bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            // Loop over the boundary box indices and compute the nearest
            // interior index.
            for (int depth = 0; depth < q_integral_data->getDepth(); ++depth)
            {
                auto acoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
                auto bcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
                auto gcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
                d_Q_bc_coef[Q_var][depth]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, Q_var, patch, trimmed_bdry_box, fill_time);
                for (CellIterator b = CellGeometry::begin(bc_coef_box); b != CellGeometry::end(bc_coef_box); ++b)
                {
                    const CellIndex& i = *b;
                    const FaceIndex i_f(i, bdry_normal_axis, FaceIndex::Lower);

                    bool is_inflow_bdry = (is_lower && (*u_data)(i_f) > 0.0) || (!is_lower && (*u_data)(i_f) < 0.0);
                    if (is_inflow_bdry)
                    {
                        const double& a = (*acoef_data)(i, 0);
                        const double& b = (*bcoef_data)(i, 0);
                        const double& g = (*gcoef_data)(i, 0);
                        const double& h = dx[bdry_normal_axis];

                        Index i_intr(i);
                        if (is_lower)
                        {
                            // intentionally blank
                        }
                        else
                        {
                            i_intr(bdry_normal_axis) -= 1;
                        }

                        const FaceIndex i_f_bdry(
                            i_intr, bdry_normal_axis, (is_lower ? FaceIndex::Lower : FaceIndex::Upper));
                        const FaceIndex i_f_intr(
                            i_intr, bdry_normal_axis, (is_lower ? FaceIndex::Upper : FaceIndex::Lower));
                        const double& q_i = (*q_integral_data)(i_f_intr, depth);
                        const double q_b = (b * q_i + g * h) / (a * h + b);
                        (*q_integral_data)(i_f_bdry, depth) = q_b;
                    }
                }
            }
        }
    }
    return;
} // setInflowBoundaryConditions

void AdvectorPredictorCorrectorHyperbolicPatchOps::getFromInput(boost::shared_ptr<Database> db,
                                                                bool /*is_from_restart*/)
{
    TBOX_ASSERT(db);

    if (db->keyExists("compute_init_velocity")) d_compute_init_velocity = db->getBool("compute_init_velocity");
    if (db->keyExists("compute_half_velocity")) d_compute_half_velocity = db->getBool("compute_half_velocity");
    if (db->keyExists("compute_final_velocity")) d_compute_final_velocity = db->getBool("compute_final_velocity");

    if (db->keyExists("extrap_type")) d_extrap_type = db->getString("extrap_type");
    if (!(d_extrap_type == "CONSTANT" || d_extrap_type == "LINEAR" || d_extrap_type == "QUADRATIC"))
    {
        TBOX_ERROR("AdvectorPredictorCorrectorHyperbolicPatchOps::getFromInput():\n"
                   << "  unknown extrapolation type: " << d_extrap_type << "\n"
                   << "  valid selections are: CONSTANT, LINEAR, or QUADRATIC" << std::endl);
    }

    if (db->keyExists("Refinement_data"))
    {
        auto refine_db = db->getDatabase("Refinement_data");
        std::vector<std::string> refinement_keys = refine_db->getAllKeys();
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

        std::vector<std::string> ref_keys_defined(num_keys);
        int def_key_cnt = 0;
        boost::shared_ptr<Database> error_db;

        for (int i = 0; i < refinement_keys.getSize(); ++i)
        {
            std::string error_key = refinement_keys[i];
            error_db.reset();

            if (error_key != "refine_criteria")
            {
                if (!(error_key == "QVAL_DEVIATION" || error_key == "QVAL_GRADIENT"))
                {
                    TBOX_ERROR(d_object_name << ":\n"
                                             << "  Unknown refinement criteria: " << error_key << "\nin input.\n");
                }
                else
                {
                    error_db = refine_db->getDatabase(error_key);
                    ref_keys_defined[def_key_cnt] = error_key;
                    def_key_cnt++;
                }

                if (error_db && error_key == "QVAL_DEVIATION")
                {
                    if (error_db->keyExists("dev_tol"))
                    {
                        d_dev_tol = error_db->getDoubleArray("dev_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name << ":\n"
                                                 << "  No key `dev_tol' found in data for " << error_key << "\n");
                    }

                    if (error_db->keyExists("qval_dev"))
                    {
                        d_dev = error_db->getDoubleArray("qval_dev");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name << ":\n"
                                                 << "  No key `qval_dev' found in data for " << error_key << "\n");
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

                if (error_db && error_key == "QVAL_GRADIENT")
                {
                    if (error_db->keyExists("grad_tol"))
                    {
                        d_grad_tol = error_db->getDoubleArray("grad_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name << ":\n"
                                                 << "  No key `grad_tol' found in data for " << error_key << "\n");
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
} // getFromInput

void AdvectorPredictorCorrectorHyperbolicPatchOps::getFromRestart()
{
    auto root_db = RestartManager::getManager()->getRootDatabase();

    boost::shared_ptr<Database> db;

    if (root_db->isDatabase(d_object_name))
    {
        db = root_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  "
                                 << "Restart database corresponding to " << d_object_name
                                 << " not found in restart file.");
    }

    // db->getIntegerArray("d_flux_ghosts", &d_flux_ghosts[0], NDIM);
    // if (d_flux_ghosts != IntVector(FLUXG))
    // {
    // TBOX_ERROR(d_object_name << ":\n"
    // << "  Key data `d_flux_ghosts' in restart file != FLUXG.\n");
    // }

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
} // getFromRestart

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
