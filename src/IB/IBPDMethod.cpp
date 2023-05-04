// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2020 by the IBAMR developers
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

#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBPDForceGen.h"
#include "ibamr/IBPDMethod.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LInitStrategy.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/ibtk_utilities.h"

#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include "petscvec.h"

#include <math.h>
#include <stddef.h>

#include <limits>
#include <ostream>
#include <string>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of IBPDMethod restart file data.
static const int IB_PD_METHOD_VERSION = 1;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBPDMethod::IBPDMethod(std::string object_name, Pointer<Database> input_db, bool register_for_restart)
    : IBMethod(std::move(object_name), input_db, register_for_restart)
{
    // NOTE: Parent class constructor registers class with the restart manager, sets object
    // name.

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    return;
} // IBPDMethod

void
IBPDMethod::registerInitialCoordinateMappingFunction(CoordinateMappingFcnPtr fcn, void* ctx)
{
    registerInitialCoordinateMappingFunction(CoordinateMappingFcnData(fcn, ctx));
    return;
} // registerInitialCoordinateMappingFunction

void
IBPDMethod::registerInitialCoordinateMappingFunction(const CoordinateMappingFcnData& data)
{
    d_coordinate_mapping_fcn_data = data;
    return;
} // registerInitialCoordinateMappingFunction

void
IBPDMethod::registerIBPDForceGen(Pointer<IBPDForceGen> ib_pd_force_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(ib_pd_force_fcn);
#endif
    d_ib_pd_force_fcn = ib_pd_force_fcn;
    return;
} // registerIBPDForceGen

void
IBPDMethod::setIBPDForceGenNeedsInit()
{
    d_ib_pd_force_fcn_needs_init = true;
    return;
} // setIBPDForceGenNeedsInit()

void
IBPDMethod::registerEulerianVariables()
{
    IBMethod::registerEulerianVariables();
    return;
} // registerEulerianVariables

void
IBPDMethod::registerEulerianCommunicationAlgorithms()
{
    IBMethod::registerEulerianCommunicationAlgorithms();
    return;
} // registerEulerianCommunicationAlgorithms

void
IBPDMethod::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    IBMethod::preprocessIntegrateData(current_time, new_time, num_cycles);
    const double start_time = d_ib_solver->getStartTime();
    if (d_ib_pd_force_fcn)
    {
        if (d_ib_pd_force_fcn_needs_init)
        {
            const bool initial_time = MathUtilities<double>::equalEps(current_time, start_time);
            resetLagrangianPDForceFunction(current_time, initial_time);
            d_ib_pd_force_fcn_needs_init = false;
        }
        d_ib_pd_force_fcn->setTimeInterval(current_time, new_time);
    }

    return;
} // preprocessIntegrateData

void
IBPDMethod::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    IBMethod::postprocessIntegrateData(current_time, new_time, num_cycles);
    return;
} // postprocessIntegrateData

void
IBPDMethod::interpolateVelocity(const int u_data_idx,
                                const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                const double data_time)
{
    // Interpolate fluid velocity.
    IBMethod::interpolateVelocity(u_data_idx, u_synch_scheds, u_ghost_fill_scheds, data_time);

    return;
} // interpolateVelocity

void
IBPDMethod::forwardEulerStep(const double current_time, const double new_time)
{
    IBMethod::forwardEulerStep(current_time, new_time);

    return;
} // forwardEulerStep

void
IBPDMethod::midpointStep(const double current_time, const double new_time)
{
    IBMethod::midpointStep(current_time, new_time);
    return;

} // midpointStep

void
IBPDMethod::trapezoidalStep(const double /*current_time*/, const double /*new_time*/)
{
    TBOX_ERROR(d_object_name << "::trapezoidalStep():\n"
                             << "  time-stepping type TRAPEZOIDAL_RULE not supported by class IBPDMethod;\n"
                             << "  use MIDPOINT_RULE instead.\n");
    return;
} // trapezoidalStep

void
IBPDMethod::computeLagrangianForce(const double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    std::vector<Pointer<LData> >* F_data = nullptr;
    std::vector<Pointer<LData> >* X_data = nullptr;
    std::vector<Pointer<LData> >* U_data = nullptr;
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        d_F_current_needs_ghost_fill = true;
        F_data = &d_F_current_data;
        X_data = &d_X_current_data;
        U_data = &d_U_current_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        d_F_half_needs_ghost_fill = true;
        F_data = &d_F_half_data;
        X_data = &d_X_half_data;
        U_data = &d_U_half_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        d_F_new_needs_ghost_fill = true;
        F_data = &d_F_new_data;
        X_data = &d_X_new_data;
        U_data = &d_U_new_data;
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        if (d_ib_pd_force_fcn)
        {
            d_ib_pd_force_fcn->computeLagrangianForceAndDamage((*F_data)[ln],
                                                               d_l_data_manager->getLData("damage", ln),
                                                               d_l_data_manager->getLData("jacobian", ln),
                                                               (*X_data)[ln],
                                                               (*U_data)[ln],
                                                               d_hierarchy,
                                                               ln,
                                                               data_time,
                                                               d_l_data_manager);
        }
    }

    return;
} // computeLagrangianForce

void
IBPDMethod::spreadForce(const int f_data_idx,
                        RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                        const std::vector<Pointer<RefineSchedule<NDIM> > >& f_prolongation_scheds,
                        const double data_time)
{
    IBMethod::spreadForce(f_data_idx, f_phys_bdry_op, f_prolongation_scheds, data_time);

    return;
} // spreadForce

void
IBPDMethod::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                     Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
                                     int u_data_idx,
                                     const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                     const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                     int integrator_step,
                                     double init_data_time,
                                     bool initial_time)
{
    // Initialize various Lagrangian data objects required by the conventional
    // IB method.
    IBMethod::initializePatchHierarchy(hierarchy,
                                       gridding_alg,
                                       u_data_idx,
                                       u_synch_scheds,
                                       u_ghost_fill_scheds,
                                       integrator_step,
                                       init_data_time,
                                       initial_time);

    // Indicate that the force strategy needs to be re-initialized.
    d_ib_pd_force_fcn_needs_init = true;

    const int struct_ln = d_hierarchy->getFinestLevelNumber();
    if (initial_time)
    {
        VecSet(d_l_data_manager->getLData("damage", struct_ln)->getVec(), 0.0);
        VecSet(d_l_data_manager->getLData("jacobian", struct_ln)->getVec(), 1.0);
    }

    // Register plot quantities.
    if (d_silo_writer)
    {
        Pointer<LData> dmg_data = d_l_data_manager->getLData("damage", struct_ln);
        d_silo_writer->registerVariableData("damage", dmg_data, struct_ln);

        Pointer<LData> jacobian_data = d_l_data_manager->getLData("jacobian", struct_ln);
        d_silo_writer->registerVariableData("jacobian", jacobian_data, struct_ln);
    }

    // Initialize unshifted X0 data.
    {
        Pointer<LData> X0_unshifted_data = d_l_data_manager->getLData("X0_unshifted", struct_ln);
        Pointer<LData> X0_data = d_l_data_manager->getLData("X0", struct_ln);

        boost::multi_array_ref<double, 2>& X0_unshifted_data_array = *X0_unshifted_data->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& X0_data_array = *X0_data->getLocalFormVecArray();

        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(struct_ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        for (const auto& node : local_nodes)
        {
            const int local_idx = node->getLocalPETScIndex();
            const Vector& displacement_0 = node->getInitialPeriodicDisplacement();
            double* const X0_unshifted = &X0_unshifted_data_array[local_idx][0];
            const double* const X0 = &X0_data_array[local_idx][0];

            for (int d = 0; d < NDIM; ++d)
            {
                X0_unshifted[d] = X0[d] + displacement_0[d];
            }
        }

        X0_unshifted_data->restoreArrays();
        X0_data->restoreArrays();
    }

    return;
} // initializePatchHierarchy

void
IBPDMethod::initializePDData()
{
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) return;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        const bool identity_mapping = !d_coordinate_mapping_fcn_data.fcn;
        if (identity_mapping) break;

        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        Pointer<LData> X0_data = d_l_data_manager->getLData("X0_unshifted", ln);
        Pointer<LData> X_data = d_l_data_manager->getLData("X", ln);
        boost::multi_array_ref<double, 2>& X0_data_array = *X0_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& X_data_array = *X_data->getLocalFormVecArray();

        for (const auto& node : local_nodes)
        {
            const int lag_idx = node->getLagrangianIndex();
            const int local_idx = node->getLocalPETScIndex();
            const double* X0 = &X0_data_array[local_idx][0];
            double* X = &X_data_array[local_idx][0];
            Eigen::Map<const IBTK::Point> eig_X0(X0);
            Eigen::Map<IBTK::Point> eig_X(X);
            d_coordinate_mapping_fcn_data.fcn(eig_X, eig_X0, lag_idx, ln, d_coordinate_mapping_fcn_data.ctx);
        }
        X0_data->restoreArrays();
        X_data->restoreArrays();

        // Update ghost DOFs
        X_data->beginGhostUpdate();
        X_data->endGhostUpdate();

        // Redistribute nodes.
        d_l_data_manager->beginDataRedistribution(ln, ln);
        d_l_data_manager->endDataRedistribution(ln, ln);
    }

    d_ib_pd_force_fcn_needs_init = true;

    return;

} // initializePDData

void
IBPDMethod::registerLoadBalancer(Pointer<LoadBalancer<NDIM> > load_balancer, int workload_data_idx)
{
    IBMethod::registerLoadBalancer(load_balancer,workload_data_idx);
    return;
} // registerLoadBalancer

void
IBPDMethod::addWorkloadEstimate(Pointer<PatchHierarchy<NDIM> > hierarchy, const int workload_data_idx)
{
    IBMethod::addWorkloadEstimate(hierarchy,workload_data_idx);
    return;
} // addWorkloadEstimate

void IBPDMethod::beginDataRedistribution(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                       Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    IBMethod::beginDataRedistribution(hierarchy,gridding_alg);
    return;
} // beginDataRedistribution

void IBPDMethod::endDataRedistribution(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                       Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    IBMethod::endDataRedistribution(hierarchy,gridding_alg);

    d_ib_pd_force_fcn_needs_init = true;

    return;
} // endDataRedistribution

void
IBPDMethod::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                int level_number,
                                double init_data_time,
                                bool can_be_refined,
                                bool initial_time,
                                Pointer<BasePatchLevel<NDIM> > old_level,
                                bool allocate_data)
{
    IBMethod::initializeLevelData(
        hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);

    // Allocate LData corresponding to the Lagrange multiplier.
    if (initial_time && d_l_data_manager->levelContainsLagrangianData(level_number))
    {   
        // Create LData for damage variable.
        d_l_data_manager->createLData("damage", level_number, 1, /*manage_data*/ true);
        d_l_data_manager->createLData("jacobian", level_number, 1, /*manage_data*/ true);

        // Create unshifted initial position of the structure.
        Pointer<IBTK::LData> X0_unshifted_data = d_l_data_manager->createLData("X0_unshifted",
                                                                               level_number,
                                                                               NDIM,
                                                                               /*manage_data*/ true);
    }

    return;
} // initializeLevelData

void
IBPDMethod::resetHierarchyConfiguration(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                      int coarsest_level,
                                      int finest_level)
{
    IBMethod::resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    return;
} // resetHierarchyConfiguration

void
IBPDMethod::putToDatabase(Pointer<Database> db)
{
    IBMethod::putToDatabase(db);
    db->putInteger("IB_PD_METHOD_VERSION", IB_PD_METHOD_VERSION);
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBPDMethod::resetLagrangianPDForceFunction(const double init_data_time, const bool initial_time)
{
    if (!d_ib_pd_force_fcn) return;
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        d_ib_pd_force_fcn->initializeLevelData(d_hierarchy, ln, init_data_time, initial_time, d_l_data_manager);
    }
    return;
} // resetLagrangianPDForceFunction

void
IBPDMethod::getFromInput(Pointer<Database> /*db*/, bool /*is_from_restart*/)
{
    // intentionally blank
    return;
} // getFromInput

void
IBPDMethod::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("IB_PD_METHOD_VERSION");
    if (ver != IB_PD_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
