// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#include "ibamr/AdvectorExplicitPredictorPatchOps.h"
#include "ibamr/AdvectorPredictorCorrectorHyperbolicPatchOps.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CartExtrapPhysBdryOp.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/CartGridFunctionSet.h"
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/PhysicalBoundaryUtilities.h"

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellDataFactory.h"
#include "CellIterator.h"
#include "CellVariable.h"
#include "ComponentSelector.h"
#include "FaceData.h"
#include "FaceIndex.h"
#include "FaceVariable.h"
#include "Geometry.h"
#include "HyperbolicLevelIntegrator.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchCellDataOpsReal.h"
#include "PatchFaceDataOpsReal.h"
#include "PatchGeometry.h"
#include "PatchLevel.h"
#include "RobinBcCoefStrategy.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "VisItDataWriter.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

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

extern "C"
{
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
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvectorPredictorCorrectorHyperbolicPatchOps::AdvectorPredictorCorrectorHyperbolicPatchOps(
    const std::string& object_name,
    SAMRAIPointer<Database> input_db,
    SAMRAIPointer<AdvectorExplicitPredictorPatchOps> explicit_predictor,
    SAMRAIPointer<CartesianGridGeometryNd> grid_geom,
    bool register_for_restart)
    : d_explicit_predictor(explicit_predictor),
      d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_grid_geometry(grid_geom)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(d_explicit_predictor);
    TBOX_ASSERT(grid_geom);
#endif

    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from given input/restart databases.
    bool is_from_restart = RestartManager::getManager()->isFromRestart();
    if (is_from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, is_from_restart);
    // Get number of ghost cells from the explicit predictor.
    d_ghosts = IntVectorNd(explicit_predictor->getNumberCellGhosts());
    d_flux_ghosts = IntVectorNd(explicit_predictor->getNumberFluxGhosts());
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

const std::string&
AdvectorPredictorCorrectorHyperbolicPatchOps::getName() const
{
    return d_object_name;
} // getName

void
AdvectorPredictorCorrectorHyperbolicPatchOps::registerVisItDataWriter(SAMRAIPointer<VisItDataWriterNd> visit_writer)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(visit_writer);
#endif
    d_visit_writer = visit_writer;
    return;
} // registerVisItDataWriter

void
AdvectorPredictorCorrectorHyperbolicPatchOps::registerAdvectionVelocity(SAMRAIPointer<FaceVariableNd<double> > u_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(u_var);
#endif
    d_u_var.insert(u_var);
    d_u_is_div_free[u_var] = true;
    return;
} // registerAdvectionVelocity

void
AdvectorPredictorCorrectorHyperbolicPatchOps::setAdvectionVelocityIsDivergenceFree(
    SAMRAIPointer<FaceVariableNd<double> > u_var,
    const bool is_div_free)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_u_var.find(u_var) != d_u_var.end());
#endif
    d_u_is_div_free[u_var] = is_div_free;
    return;
} // setAdvectionVelocityIsDivergenceFree

void
AdvectorPredictorCorrectorHyperbolicPatchOps::setAdvectionVelocityFunction(SAMRAIPointer<FaceVariableNd<double> > u_var,
                                                                           SAMRAIPointer<CartGridFunction> u_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_u_var.find(u_var) != d_u_var.end());
#endif
    d_u_fcn[u_var] = u_fcn;
    return;
} // setAdvectionVelocityFunction

void
AdvectorPredictorCorrectorHyperbolicPatchOps::registerSourceTerm(SAMRAIPointer<CellVariableNd<double> > F_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(F_var);
#endif
    d_F_var.insert(F_var);
    return;
} // registerSourceTerm

void
AdvectorPredictorCorrectorHyperbolicPatchOps::setSourceTermFunction(SAMRAIPointer<CellVariableNd<double> > F_var,
                                                                    SAMRAIPointer<CartGridFunction> F_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_F_var.find(F_var) != d_F_var.end());
#endif
    if (d_F_fcn[F_var])
    {
        const std::string& F_var_name = F_var->getName();
        SAMRAIPointer<CartGridFunctionSet> p_F_fcn = d_F_fcn[F_var];
        if (!p_F_fcn)
        {
            pout << d_object_name << "::setSourceTermFunction(): WARNING:\n"
                 << "  source term function for source term variable " << F_var_name << " has already been set.\n"
                 << "  functions will be evaluated in the order in which they were "
                    "registered with the solver\n"
                 << "  when evaluating the source term value.\n";
            p_F_fcn = new CartGridFunctionSet(d_object_name + "::" + F_var_name + "::source_function_set");
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

void
AdvectorPredictorCorrectorHyperbolicPatchOps::registerTransportedQuantity(SAMRAIPointer<CellVariableNd<double> > Q_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_var);
#endif
    d_Q_var.insert(Q_var);
    d_Q_difference_form[Q_var] = CONSERVATIVE;
    return;
} // registerTransportedQuantity

void
AdvectorPredictorCorrectorHyperbolicPatchOps::setAdvectionVelocity(SAMRAIPointer<CellVariableNd<double> > Q_var,
                                                                   SAMRAIPointer<FaceVariableNd<double> > u_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    TBOX_ASSERT(d_u_var.find(u_var) != d_u_var.end());
#endif
    d_Q_u_map[Q_var] = u_var;
    return;
} // setAdvectionVelocity

void
AdvectorPredictorCorrectorHyperbolicPatchOps::setSourceTerm(SAMRAIPointer<CellVariableNd<double> > Q_var,
                                                            SAMRAIPointer<CellVariableNd<double> > F_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    TBOX_ASSERT(d_F_var.find(F_var) != d_F_var.end());
#endif
    d_Q_F_map[Q_var] = F_var;
    return;
} // setSourceTerm

void
AdvectorPredictorCorrectorHyperbolicPatchOps::setConvectiveDifferencingType(
    SAMRAIPointer<CellVariableNd<double> > Q_var,
    const ConvectiveDifferencingType difference_form)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_difference_form[Q_var] = difference_form;
    return;
} // setConvectiveDifferencingType

void
AdvectorPredictorCorrectorHyperbolicPatchOps::setInitialConditions(SAMRAIPointer<CellVariableNd<double> > Q_var,
                                                                   SAMRAIPointer<CartGridFunction> Q_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_init[Q_var] = Q_init;
    return;
} // setInitialConditions

void
AdvectorPredictorCorrectorHyperbolicPatchOps::setPhysicalBcCoefs(SAMRAIPointer<CellVariableNd<double> > Q_var,
                                                                 RobinBcCoefStrategyNd* Q_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    SAMRAIPointer<CellDataFactoryNd<double> > Q_factory = Q_var->getPatchDataFactory();
    const unsigned int Q_depth = Q_factory->getDefaultDepth();
    TBOX_ASSERT(Q_depth == 1);
#endif
    d_Q_bc_coef[Q_var] = std::vector<RobinBcCoefStrategyNd*>(1, Q_bc_coef);
    return;
} // setPhysicalBcCoefs

void
AdvectorPredictorCorrectorHyperbolicPatchOps::setPhysicalBcCoefs(SAMRAIPointer<CellVariableNd<double> > Q_var,
                                                                 std::vector<RobinBcCoefStrategyNd*> Q_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    SAMRAIPointer<CellDataFactoryNd<double> > Q_factory = Q_var->getPatchDataFactory();
    const unsigned int Q_depth = Q_factory->getDefaultDepth();
    TBOX_ASSERT(Q_depth == Q_bc_coef.size());
#endif
    d_Q_bc_coef[Q_var] = Q_bc_coef;
    return;
} // setPhysicalBcCoefs

void
AdvectorPredictorCorrectorHyperbolicPatchOps::registerModelVariables(HyperbolicLevelIntegratorNd* integrator)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(integrator);
#endif
    d_integrator = integrator;

    for (const auto& u_var : d_u_var)
    {
        d_integrator->registerVariable(u_var,
                                       d_ghosts,
                                       HyperbolicLevelIntegratorNd::TIME_DEP,
                                       d_grid_geometry,
                                       "CONSERVATIVE_COARSEN",
                                       "CONSERVATIVE_LINEAR_REFINE");
    }

    for (const auto& F_var : d_F_var)
    {
        d_integrator->registerVariable(F_var,
                                       d_ghosts,
                                       HyperbolicLevelIntegratorNd::TIME_DEP,
                                       d_grid_geometry,
                                       "CONSERVATIVE_COARSEN",
                                       "CONSERVATIVE_LINEAR_REFINE");
    }

    for (const auto& Q_var : d_Q_var)
    {
        SAMRAIPointer<CellDataFactoryNd<double> > Q_factory = Q_var->getPatchDataFactory();
        d_integrator->registerVariable(Q_var,
                                       d_ghosts,
                                       HyperbolicLevelIntegratorNd::TIME_DEP,
                                       d_grid_geometry,
                                       "CONSERVATIVE_COARSEN",
                                       "CONSERVATIVE_LINEAR_REFINE");

        if (d_visit_writer)
        {
            const int Q_idx =
                VariableDatabaseNd::getDatabase()->mapVariableAndContextToIndex(Q_var, d_integrator->getPlotContext());
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
                    d_visit_writer->registerPlotQuantity(
                        Q_var->getName() + "_" + std::to_string(d), "SCALAR", Q_idx, d);
                }
            }
        }

        const bool conservation_form = d_Q_difference_form[Q_var] == CONSERVATIVE;
        if (conservation_form)
        {
            d_flux_integral_var[Q_var] =
                new FaceVariableNd<double>(d_object_name + "::" + Q_var->getName() + " advective flux time integral",
                                           Q_factory->getDefaultDepth());
            d_integrator->registerVariable(d_flux_integral_var[Q_var],
                                           d_flux_ghosts,
                                           HyperbolicLevelIntegratorNd::FLUX,
                                           d_grid_geometry,
                                           "CONSERVATIVE_COARSEN",
                                           "NO_REFINE");
        }

        SAMRAIPointer<FaceVariableNd<double> > u_var = d_Q_u_map[Q_var];
        const bool u_is_div_free = d_u_is_div_free[u_var];
        if (!conservation_form || !u_is_div_free)
        {
            d_q_integral_var[Q_var] = new FaceVariableNd<double>(
                d_object_name + "::" + Q_var->getName() + " time integral", Q_factory->getDefaultDepth());
            d_integrator->registerVariable(d_q_integral_var[Q_var],
                                           d_flux_ghosts,
                                           HyperbolicLevelIntegratorNd::FLUX,
                                           d_grid_geometry,
                                           "CONSERVATIVE_COARSEN",
                                           "NO_REFINE");
            if (u_var && !d_u_integral_var[u_var])
            {
                d_u_integral_var[u_var] =
                    new FaceVariableNd<double>(d_object_name + "::" + u_var->getName() + " time integral");
                d_integrator->registerVariable(d_u_integral_var[u_var],
                                               d_flux_ghosts,
                                               HyperbolicLevelIntegratorNd::FLUX,
                                               d_grid_geometry,
                                               "CONSERVATIVE_COARSEN",
                                               "NO_REFINE");
            }
        }
    }
    return;
} // registerModelVariables

void
AdvectorPredictorCorrectorHyperbolicPatchOps::initializeDataOnPatch(PatchNd& patch,
                                                                    const double data_time,
                                                                    const bool initial_time)
{
    if (initial_time)
    {
        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
        for (const auto& u_var : d_u_var)
        {
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, getDataContext());
            if (d_u_fcn[u_var])
            {
                d_u_fcn[u_var]->setDataOnPatch(
                    u_idx, u_var, SAMRAIPointer<PatchNd>(&patch, false), data_time, initial_time);
            }
            else
            {
                SAMRAIPointer<FaceDataNd<double> > u_data = patch.getPatchData(u_idx);
                u_data->fillAll(0.0);
            }
        }

        for (const auto& F_var : d_F_var)
        {
            const int F_idx = var_db->mapVariableAndContextToIndex(F_var, getDataContext());
            if (d_F_fcn[F_var])
            {
                d_F_fcn[F_var]->setDataOnPatch(
                    F_idx, F_var, SAMRAIPointer<PatchNd>(&patch, false), data_time, initial_time);
            }
            else
            {
                SAMRAIPointer<CellDataNd<double> > F_data = patch.getPatchData(F_idx);
                F_data->fillAll(0.0);
            }
        }

        for (const auto& Q_var : d_Q_var)
        {
            const int Q_idx = var_db->mapVariableAndContextToIndex(Q_var, getDataContext());
            if (d_Q_init[Q_var])
            {
                d_Q_init[Q_var]->setDataOnPatch(
                    Q_idx, Q_var, SAMRAIPointer<PatchNd>(&patch, false), data_time, initial_time);
            }
            else
            {
                SAMRAIPointer<CellDataNd<double> > Q_data = patch.getPatchData(Q_var, getDataContext());
                Q_data->fillAll(0.0);
            }
        }
    }
    return;
} // initializeDataOnPatch

double
AdvectorPredictorCorrectorHyperbolicPatchOps::computeStableDtOnPatch(PatchNd& patch,
                                                                     const bool /*initial_time*/,
                                                                     const double /*dt_time*/)
{
    double stable_dt = std::numeric_limits<double>::max();
    for (const auto& u_var : d_u_var)
    {
        SAMRAIPointer<FaceDataNd<double> > u_data = patch.getPatchData(u_var, getDataContext());
        stable_dt = std::min(stable_dt, d_explicit_predictor->computeStableDtOnPatch(*u_data, patch));
    }
    return stable_dt;
} // computeStableDtOnPatch

void
AdvectorPredictorCorrectorHyperbolicPatchOps::computeFluxesOnPatch(PatchNd& patch, const double time, const double dt)
{
    SAMRAIPointer<PatchGeometryNd> pgeom = patch.getPatchGeometry();
    const BoxNd& patch_box = patch.getBox();

    PatchFaceDataOpsRealNd<double> patch_fc_data_ops;

    for (const auto& Q_var : d_Q_var)
    {
        SAMRAIPointer<FaceVariableNd<double> > u_var = d_Q_u_map[Q_var];
        SAMRAIPointer<FaceDataNd<double> > q_integral_data = getQIntegralData(Q_var, patch, getDataContext());
        if (!u_var)
        {
            q_integral_data->fillAll(0.0);
            continue;
        }

        // Predict time- and face-centered values.
        SAMRAIPointer<CellDataNd<double> > Q_data = patch.getPatchData(Q_var, getDataContext());
        SAMRAIPointer<FaceDataNd<double> > u_data = patch.getPatchData(u_var, getDataContext());
        SAMRAIPointer<CellVariableNd<double> > F_var = d_Q_F_map[Q_var];
        if (F_var)
        {
            SAMRAIPointer<CellDataNd<double> > F_data = patch.getPatchData(F_var, getDataContext());
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
    VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
    if (d_compute_half_velocity)
    {
        for (const auto& u_var : d_u_var)
        {
            if (d_u_fcn[u_var] && d_u_fcn[u_var]->isTimeDependent())
            {
                const int u_idx = var_db->mapVariableAndContextToIndex(u_var, getDataContext());
                d_u_fcn[u_var]->setDataOnPatch(u_idx, u_var, SAMRAIPointer<PatchNd>(&patch, false), time + 0.5 * dt);
            }
        }
    }

    // Compute fluxes and other face-centered quantities.
    for (const auto& Q_var : d_Q_var)
    {
        SAMRAIPointer<FaceVariableNd<double> > u_var = d_Q_u_map[Q_var];

        if (!u_var) continue;

        SAMRAIPointer<FaceDataNd<double> > u_data = patch.getPatchData(u_var, getDataContext());
        const bool conservation_form = d_Q_difference_form[Q_var] == CONSERVATIVE;
        const bool u_is_div_free = d_u_is_div_free[u_var];

        if (conservation_form)
        {
            SAMRAIPointer<FaceDataNd<double> > flux_integral_data = getFluxIntegralData(Q_var, patch, getDataContext());
            SAMRAIPointer<FaceDataNd<double> > q_integral_data = getQIntegralData(Q_var, patch, getDataContext());
            d_explicit_predictor->computeFlux(*flux_integral_data, *u_data, *q_integral_data, patch, dt);
        }

        if (!conservation_form || !u_is_div_free)
        {
            SAMRAIPointer<FaceDataNd<double> > q_integral_data = getQIntegralData(Q_var, patch, getDataContext());
            patch_fc_data_ops.scale(q_integral_data, dt, q_integral_data, patch_box);

            SAMRAIPointer<FaceDataNd<double> > u_integral_data = getUIntegralData(Q_var, patch, getDataContext());
            patch_fc_data_ops.scale(u_integral_data, dt, u_data, patch_box);
        }
    }
    return;
} // computeFluxesOnPatch

void
AdvectorPredictorCorrectorHyperbolicPatchOps::conservativeDifferenceOnPatch(PatchNd& patch,
                                                                            const double /*time*/,
                                                                            const double dt,
                                                                            bool /*at_synchronization*/)
{
    const BoxNd& patch_box = patch.getBox();
    const hier::IndexNd& ilower = patch_box.lower();
    const hier::IndexNd& iupper = patch_box.upper();

    const SAMRAIPointer<CartesianPatchGeometryNd> patch_geom = patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    for (const auto& Q_var : d_Q_var)
    {
        SAMRAIPointer<FaceVariableNd<double> > u_var = d_Q_u_map[Q_var];

        if (!u_var) continue;

        SAMRAIPointer<CellDataNd<double> > Q_data = patch.getPatchData(Q_var, getDataContext());
        SAMRAIPointer<FaceDataNd<double> > flux_integral_data = getFluxIntegralData(Q_var, patch, getDataContext());
        SAMRAIPointer<FaceDataNd<double> > q_integral_data = getQIntegralData(Q_var, patch, getDataContext());
        SAMRAIPointer<FaceDataNd<double> > u_integral_data = getUIntegralData(Q_var, patch, getDataContext());

        const IntVectorNd& Q_data_ghost_cells = Q_data->getGhostCellWidth();
        const IntVectorNd& flux_integral_data_ghost_cells =
            (flux_integral_data ? flux_integral_data->getGhostCellWidth() : 0);
        const IntVectorNd& q_integral_data_ghost_cells = (q_integral_data ? q_integral_data->getGhostCellWidth() : 0);
        const IntVectorNd& u_integral_data_ghost_cells = (u_integral_data ? u_integral_data->getGhostCellWidth() : 0);

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
            CellDataNd<double> N_data(patch_box, Q_data->getDepth(), 0);
            d_explicit_predictor->computeAdvectiveDerivative(N_data, *u_integral_data, *q_integral_data, patch);
            PatchCellDataOpsRealNd<double> patch_cc_data_ops;
            patch_cc_data_ops.axpy(
                Q_data, -1.0 / dt, SAMRAIPointer<CellDataNd<double> >(&N_data, false), Q_data, patch_box);
            break;
        }
        default:
        {
            TBOX_ERROR(
                "AdvectorPredictorCorrectorHyperbolicPatchOps::"
                "conservativeDifferenceOnPatch():\n"
                << "  unsupported differencing form: "
                << enum_to_string<ConvectiveDifferencingType>(d_Q_difference_form[Q_var]) << " \n"
                << "  valid choices are: ADVECTIVE, CONSERVATIVE\n");
        }
        }
    }
    return;
} // conservativeDifferenceOnPatch

void
AdvectorPredictorCorrectorHyperbolicPatchOps::preprocessAdvanceLevelState(const SAMRAIPointer<PatchLevelNd>& level,
                                                                          double current_time,
                                                                          double /*dt*/,
                                                                          bool /*first_step*/,
                                                                          bool /*last_step*/,
                                                                          bool /*regrid_advance*/)
{
    VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();

    // Update the source term.
    for (const auto& F_var : d_F_var)
    {
        if (d_F_fcn[F_var] && d_F_fcn[F_var]->isTimeDependent())
        {
            const int F_idx = var_db->mapVariableAndContextToIndex(F_var, d_integrator->getScratchContext());
            d_F_fcn[F_var]->setDataOnPatchLevel(F_idx, F_var, level, current_time);
        }
    }

    if (!d_compute_init_velocity) return;

    // Update the advection velocity.
    for (const auto& u_var : d_u_var)
    {
        if (d_u_fcn[u_var] && d_u_fcn[u_var]->isTimeDependent())
        {
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, d_integrator->getScratchContext());
            d_u_fcn[u_var]->setDataOnPatchLevel(u_idx, u_var, level, current_time);
        }
    }
    return;
} // preprocessAdvanceLevelState

void
AdvectorPredictorCorrectorHyperbolicPatchOps::postprocessAdvanceLevelState(const SAMRAIPointer<PatchLevelNd>& level,
                                                                           double current_time,
                                                                           double dt,
                                                                           bool /*first_step*/,
                                                                           bool /*last_step*/,
                                                                           bool /*regrid_advance*/)
{
    VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();

    PatchCellDataOpsRealNd<double> patch_cc_data_ops;

    SAMRAIPointer<VariableContext> new_context = d_integrator->getNewContext();
    SAMRAIPointer<VariableContext> scratch_context = d_integrator->getScratchContext();

    // Update the values of any time-dependent source terms and add the values
    // of all source terms to the advected quantities.
    for (const auto& Q_var : d_Q_var)
    {
        SAMRAIPointer<CellVariableNd<double> > F_var = d_Q_F_map[Q_var];
        if (!F_var) continue;
        SAMRAIPointer<CartGridFunction> F_fcn = d_F_fcn[F_var];
        for (PatchLevelNd::Iterator p(level); p; p++)
        {
            SAMRAIPointer<PatchNd> patch = level->getPatch(p());
            const BoxNd& patch_box = patch->getBox();

            SAMRAIPointer<CellDataNd<double> > Q_data = patch->getPatchData(Q_var, new_context);
            if (F_fcn)
            {
                const int F_scratch_idx = var_db->mapVariableAndContextToIndex(F_var, scratch_context);
                const int F_new_idx = var_db->mapVariableAndContextToIndex(F_var, new_context);
                SAMRAIPointer<CellDataNd<double> > F_scratch_data = patch->getPatchData(F_scratch_idx);
                SAMRAIPointer<CellDataNd<double> > F_new_data = patch->getPatchData(F_new_idx);
                F_fcn->setDataOnPatchLevel(F_new_idx, F_var, level, current_time + dt);
                patch_cc_data_ops.axpy(Q_data, 0.5 * dt, F_scratch_data, Q_data, patch_box);
                patch_cc_data_ops.axpy(Q_data, 0.5 * dt, F_new_data, Q_data, patch_box);
            }
            else
            {
                const int F_scratch_idx = var_db->mapVariableAndContextToIndex(F_var, scratch_context);
                SAMRAIPointer<CellDataNd<double> > F_scratch_data = patch->getPatchData(F_scratch_idx);
                patch_cc_data_ops.axpy(Q_data, dt, F_scratch_data, Q_data, patch_box);
            }
        }
    }

    if (!d_compute_final_velocity) return;

    // Update the advection velocity.
    for (const auto& u_var : d_u_var)
    {
        if (d_u_fcn[u_var] && d_u_fcn[u_var]->isTimeDependent())
        {
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, d_integrator->getNewContext());
            d_u_fcn[u_var]->setDataOnPatchLevel(u_idx, u_var, level, current_time + dt);
        }
    }
    return;
} // postprocessAdvanceLevelState

void
AdvectorPredictorCorrectorHyperbolicPatchOps::tagGradientDetectorCells(PatchNd& patch,
                                                                       const double regrid_time,
                                                                       const bool /*initial_error*/,
                                                                       const int tag_indx,
                                                                       const bool /*uses_richardson_extrapolation_too*/)
{
    const int error_level_number = patch.getPatchLevelNumber();

    const SAMRAIPointer<CartesianPatchGeometryNd> patch_geom = patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const BoxNd& patch_box = patch.getBox();
    const hier::IndexNd& ilower = patch.getBox().lower();
    const hier::IndexNd& iupper = patch.getBox().upper();

    SAMRAIPointer<CellDataNd<int> > tags = patch.getPatchData(tag_indx);

    const int not_refine_tag_val = FALSE_VAL;
    const int refine_tag_val = TRUE_VAL;

    // Create a set of temporary tags and set to untagged value.
    CellDataNd<int> temp_tags(patch_box, 1, d_ghosts);
    temp_tags.fillAll(not_refine_tag_val);

    // Possible tagging criteria include: QVAL_DEVIATION, QVAL_GRADIENT.  The
    // criteria are specified over a time interval.
    //
    // Loop over criteria provided and check to make sure we are in the
    // specified time interval.  If so, apply appropriate tagging for the level.
    for (int ncrit = 0; ncrit < d_refinement_criteria.getSize(); ++ncrit)
    {
        std::string ref = d_refinement_criteria[ncrit];
        const IntVectorNd& tagghost = tags->getGhostCellWidth();

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
                for (const auto& Q_var : d_Q_var)
                {
                    SAMRAIPointer<CellDataNd<double> > Q_data = patch.getPatchData(Q_var, getDataContext());
                    for (int depth = 0; depth < Q_data->getDepth(); ++depth)
                    {
                        for (CellIteratorNd ic(patch_box); ic; ic++)
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
                for (const auto& Q_var : d_Q_var)
                {
                    SAMRAIPointer<CellDataNd<double> > Q_data = patch.getPatchData(Q_var, getDataContext());
                    const IntVectorNd& Q_ghost = Q_data->getGhostCellWidth();
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
        for (CellIteratorNd ic(patch_box); ic; ic++)
        {
            (*tags)(ic(), 0) = temp_tags(ic(), 0);
        }
    }
    else
    {
        for (CellIteratorNd ic(patch_box); ic; ic++)
        {
            (*tags)(ic(), 0) = (*tags)(ic(), 0) || temp_tags(ic(), 0);
        }
    }
    return;
} // tagGradientDetectorCells

void
AdvectorPredictorCorrectorHyperbolicPatchOps::setPhysicalBoundaryConditions(PatchNd& patch,
                                                                            const double fill_time,
                                                                            const IntVectorNd& ghost_width_to_fill)
{
    VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();

    // Extrapolate the interior data to set the ghost cell values for the state
    // variables and for any forcing terms.
    ComponentSelector u_patch_data_indices;
    for (const auto& u_var : d_u_var)
    {
        const int u_data_idx = var_db->mapVariableAndContextToIndex(u_var, d_integrator->getScratchContext());
        u_patch_data_indices.setFlag(u_data_idx);
    }
    d_extrap_bc_helper.setExtrapolationType(d_extrap_type);
    d_extrap_bc_helper.setPatchDataIndices(u_patch_data_indices);
    d_extrap_bc_helper.setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);

    ComponentSelector patch_data_indices;
    for (const auto& F_var : d_F_var)
    {
        const int F_data_idx = var_db->mapVariableAndContextToIndex(F_var, d_integrator->getScratchContext());
        patch_data_indices.setFlag(F_data_idx);
    }
    for (const auto& Q_var : d_Q_var)
    {
        const int Q_data_idx = var_db->mapVariableAndContextToIndex(Q_var, d_integrator->getScratchContext());
        patch_data_indices.setFlag(Q_data_idx);
    }
    d_extrap_bc_helper.setExtrapolationType(d_extrap_type);
    d_extrap_bc_helper.setPatchDataIndices(patch_data_indices);
    d_extrap_bc_helper.setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);
    return;
} // setPhysicalBoundaryConditions

void
AdvectorPredictorCorrectorHyperbolicPatchOps::putToDatabase(SAMRAIPointer<Database> db)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif

    // db->putIntegerArray("d_flux_ghosts", &d_flux_ghosts[0], NDIM);

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
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

SAMRAIPointer<FaceDataNd<double> >
AdvectorPredictorCorrectorHyperbolicPatchOps::getFluxIntegralData(SAMRAIPointer<CellVariableNd<double> > Q_var,
                                                                  PatchNd& patch,
                                                                  SAMRAIPointer<VariableContext> context)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_var);
#endif
    if (d_Q_difference_form[Q_var] == CONSERVATIVE)
    {
        return patch.getPatchData(d_flux_integral_var[Q_var], context);
    }
    else
    {
        return SAMRAIPointer<FaceDataNd<double> >(nullptr);
    }
} // getFluxIntegralData

SAMRAIPointer<FaceDataNd<double> >
AdvectorPredictorCorrectorHyperbolicPatchOps::getQIntegralData(SAMRAIPointer<CellVariableNd<double> > Q_var,
                                                               PatchNd& patch,
                                                               SAMRAIPointer<VariableContext> context)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_var);
#endif
    if (!d_u_is_div_free[d_Q_u_map[Q_var]] || d_Q_difference_form[Q_var] != CONSERVATIVE)
    {
        return patch.getPatchData(d_q_integral_var[Q_var], context);
    }
    else
    {
        return getFluxIntegralData(Q_var, patch, context);
    }
} // getQIntegralData

SAMRAIPointer<FaceDataNd<double> >
AdvectorPredictorCorrectorHyperbolicPatchOps::getUIntegralData(SAMRAIPointer<CellVariableNd<double> > Q_var,
                                                               PatchNd& patch,
                                                               SAMRAIPointer<VariableContext> context)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_var);
#endif
    SAMRAIPointer<FaceVariableNd<double> > u_var = d_Q_u_map[Q_var];
    if (u_var && (!d_u_is_div_free[u_var] || d_Q_difference_form[Q_var] != CONSERVATIVE))
    {
        return patch.getPatchData(d_u_integral_var[u_var], context);
    }
    else
    {
        return SAMRAIPointer<FaceDataNd<double> >(nullptr);
    }
} // getUIntegralData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
AdvectorPredictorCorrectorHyperbolicPatchOps::setInflowBoundaryConditions(PatchNd& patch, const double fill_time)
{
    VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();

    SAMRAIPointer<CartesianPatchGeometryNd> pgeom = patch.getPatchGeometry();

    // There is nothing to do if the patch does not touch a regular (physical)
    // boundary.
    if (!pgeom->getTouchesRegularBoundary()) return;

    // Compute the co-dimension one boundary boxes.
    const Array<BoundaryBoxNd> physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);

    // There is nothing to do if the patch does not have any co-dimension one
    // boundary boxes.
    if (physical_codim1_boxes.size() == 0) return;

    // Loop over the boundary fill boxes and set boundary conditions at inflow
    // boundaries only.
    const BoxNd& patch_box = patch.getBox();
    const double* const dx = pgeom->getDx();
    for (const auto& Q_var : d_Q_var)
    {
        SAMRAIPointer<FaceVariableNd<double> > u_var = d_Q_u_map[Q_var];

        if (!u_var) continue;

        const int Q_data_idx = var_db->mapVariableAndContextToIndex(Q_var, d_integrator->getScratchContext());

        SAMRAIPointer<FaceDataNd<double> > u_data = patch.getPatchData(u_var, getDataContext());
        SAMRAIPointer<FaceDataNd<double> > q_integral_data = getQIntegralData(Q_var, patch, getDataContext());

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
            const BoundaryBoxNd& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = location_index % 2 == 0;

            static const IntVectorNd gcw_to_fill = 1;
            const BoxNd bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
            const BoundaryBoxNd trimmed_bdry_box(
                bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), bdry_box.getLocationIndex());
            const BoxNd bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            // Loop over the boundary box indices and compute the nearest
            // interior index.
            for (int depth = 0; depth < q_integral_data->getDepth(); ++depth)
            {
                SAMRAIPointer<ArrayDataNd<double> > acoef_data = new ArrayDataNd<double>(bc_coef_box, 1);
                SAMRAIPointer<ArrayDataNd<double> > bcoef_data = new ArrayDataNd<double>(bc_coef_box, 1);
                SAMRAIPointer<ArrayDataNd<double> > gcoef_data = new ArrayDataNd<double>(bc_coef_box, 1);
                d_Q_bc_coef[Q_var][depth]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, Q_var, patch, trimmed_bdry_box, fill_time);
                for (BoxNd::Iterator b(bc_coef_box); b; b++)
                {
                    const hier::IndexNd& i = b();
                    const FaceIndexNd i_f(i, bdry_normal_axis, FaceIndexNd::Lower);

                    bool is_inflow_bdry = (is_lower && (*u_data)(i_f) > 0.0) || (!is_lower && (*u_data)(i_f) < 0.0);
                    if (is_inflow_bdry)
                    {
                        const double& a = (*acoef_data)(i, 0);
                        const double& b = (*bcoef_data)(i, 0);
                        const double& g = (*gcoef_data)(i, 0);
                        const double& h = dx[bdry_normal_axis];

                        hier::IndexNd i_intr(i);
                        if (is_lower)
                        {
                            // intentionally blank
                        }
                        else
                        {
                            i_intr(bdry_normal_axis) -= 1;
                        }

                        const FaceIndexNd i_f_bdry(
                            i_intr, bdry_normal_axis, (is_lower ? FaceIndexNd::Lower : FaceIndexNd::Upper));
                        const FaceIndexNd i_f_intr(
                            i_intr, bdry_normal_axis, (is_lower ? FaceIndexNd::Upper : FaceIndexNd::Lower));
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

void
AdvectorPredictorCorrectorHyperbolicPatchOps::getFromInput(SAMRAIPointer<Database> db, bool /*is_from_restart*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif
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
        SAMRAIPointer<Database> refine_db = db->getDatabase("Refinement_data");
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
        SAMRAIPointer<Database> error_db;

        for (int i = 0; i < refinement_keys.getSize(); ++i)
        {
            std::string error_key = refinement_keys[i];
            error_db.setNull();

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

void
AdvectorPredictorCorrectorHyperbolicPatchOps::getFromRestart()
{
    SAMRAIPointer<Database> root_db = RestartManager::getManager()->getRootDatabase();

    SAMRAIPointer<Database> db;

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
    // if (d_flux_ghosts != IntVectorNd(FLUXG))
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
