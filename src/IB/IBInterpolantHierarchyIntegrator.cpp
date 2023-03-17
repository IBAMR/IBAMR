// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
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

#include "ibamr/BrinkmanPenalizationRigidBodyDynamics.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBInterpolantHierarchyIntegrator.h"
#include "ibamr/IBInterpolantMethod.h"
#include "ibamr/IBLevelSetMethod.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/ibtk_enums.h"

#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchCellDataOpsReal.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchSideDataOpsReal.h"
#include "SideData.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include "Eigen/Core"

#include <algorithm>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Box;
} // namespace hier
} // namespace SAMRAI

namespace IBAMR
{
class BrinkmanPenalizationStrategy;
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of IBInterpolantHierarchyIntegrator restart file data.
static const int IB_INTERPOLANT_HIERARCHY_INTEGRATOR_VERSION = 1;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBInterpolantHierarchyIntegrator::IBInterpolantHierarchyIntegrator(std::string object_name,
                                                                   Pointer<Database> input_db,
                                                                   Pointer<IBLevelSetMethod> ib_ls_method_ops,
                                                                   Pointer<INSHierarchyIntegrator> ins_hier_integrator,
                                                                   bool register_for_restart)
    : IBHierarchyIntegrator(std::move(object_name),
                            input_db,
                            ib_ls_method_ops,
                            ins_hier_integrator,
                            register_for_restart),
      d_ib_ls_method_ops(ib_ls_method_ops)
{
    d_ib_interpolant_method_ops = d_ib_ls_method_ops->getIBMethodOps();

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    return;
} // IBInterpolantHierarchyIntegrator

void
IBInterpolantHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                               const double new_time,
                                                               const int num_cycles)
{
    // preprocess our dependencies...
    IBHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    // ... and preprocess objects owned by this class.
    bool initial_time = IBTK::abs_equal_eps(current_time, 0.0);
    if (initial_time) d_ib_interpolant_method_ops->interpolateQ();

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);

    return;
} // preprocessIntegrateHierarchy

void
IBInterpolantHierarchyIntegrator::integrateHierarchy(const double current_time,
                                                     const double new_time,
                                                     const int cycle_num)
{
    IBHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);

    // Here we implement the following time integration scheme:
    //
    // (1) Update the interpolation mesh position to X^{n+1}.
    //
    // (2) Spread the scalar to new position.
    //
    // (3) Solve the fluid equations using the updated value of the scalar.

    // Move the mesh to new location.
    Pointer<INSVCStaggeredHierarchyIntegrator> vc_ins_integrator = d_ins_hier_integrator;
    const std::vector<Pointer<BrinkmanPenalizationStrategy> >& brinkman_force =
        vc_ins_integrator->getBrinkmanPenalizationStrategy();
    const std::size_t num_objects = brinkman_force.size();
    EigenAlignedVector<Eigen::Vector3d> U(num_objects), W(num_objects);
    for (std::size_t k = 0; k < num_objects; ++k)
    {
        Pointer<BrinkmanPenalizationRigidBodyDynamics> rbd = brinkman_force[k];
        U[k] = rbd->getNewCOMTransVelocity();
        W[k] = rbd->getNewCOMRotVelocity();
    }
    d_ib_interpolant_method_ops->updateMeshPosition(current_time, new_time, U, W);

    // Spread the scalar to the grid.
    d_ib_interpolant_method_ops->spreadQ(new_time);

    // Solve the INS equations.
    d_ins_hier_integrator->integrateHierarchy(current_time, new_time, cycle_num);

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);

    return;
} // integrateHierarchy

void
IBInterpolantHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                const double new_time,
                                                                const bool skip_synchronize_new_state_data,
                                                                const int num_cycles)
{
    // We don't have any data ourselves that needs to be postprocessed so defer
    // immediately to the base class:
    IBHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    return;
} // postprocessIntegrateHierarchy

void
IBInterpolantHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                                Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    // Pass on body force function to INS integrator.
    if (d_body_force_fcn)
    {
        d_ins_hier_integrator->registerBodyForceFunction(d_body_force_fcn);
    }

    // Finish initializing the hierarchy integrator.
    IBHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);
    return;
} // initializeHierarchyIntegrator

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBInterpolantHierarchyIntegrator::putToDatabaseSpecialized(Pointer<Database> db)
{
    IBHierarchyIntegrator::putToDatabaseSpecialized(db);
    db->putInteger("IB_INTERPOLANT_HIERARCHY_INTEGRATOR_VERSION", IB_INTERPOLANT_HIERARCHY_INTEGRATOR_VERSION);
    return;
} // putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBInterpolantHierarchyIntegrator::getFromRestart()
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
    int ver = db->getInteger("IB_INTERPOLANT_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_INTERPOLANT_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
