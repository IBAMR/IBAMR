// Filename: IBInterpolantHierarchyIntegrator.cpp
// Created on 20 Nov 2018 by Amneet Bhalla
//
// Copyright (c) 2002-2018, Amneet Bhalla
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

#include "ibamr/IBInterpolantHierarchyIntegrator.h"
#include "ibamr/BrinkmanPenalizationRigidBodyDynamics.h"
#include "ibamr/IBInterpolantMethod.h"
#include "ibamr/IBLevelSetMethod.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of IBInterpolantHierarchyIntegrator restart file data.
static const int IB_INTERPOLANT_HIERARCHY_INTEGRATOR_VERSION = 1;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBInterpolantHierarchyIntegrator::IBInterpolantHierarchyIntegrator(const std::string& object_name,
                                                                   Pointer<Database> input_db,
                                                                   Pointer<IBLevelSetMethod> ib_ls_method_ops,
                                                                   Pointer<INSHierarchyIntegrator> ins_hier_integrator,
                                                                   bool register_for_restart)
    : IBHierarchyIntegrator(object_name, input_db, ib_ls_method_ops, ins_hier_integrator, register_for_restart),
      d_ib_ls_method_ops(ib_ls_method_ops)
{
    d_ib_interpolant_method_ops = d_ib_ls_method_ops->getIBMethodOps();

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    return;

    // intentionally blank
    return;
} // IBInterpolantHierarchyIntegrator

IBInterpolantHierarchyIntegrator::~IBInterpolantHierarchyIntegrator()
{
    // intentionally blank
    return;
} // ~IBInterpolantHierarchyIntegrator

void
IBInterpolantHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                               const double new_time,
                                                               const int num_cycles)
{
    IBHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_level_num = 0;
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();

    // Allocate Eulerian scratch and new data.
    for (int level_num = coarsest_level_num; level_num <= finest_level_num; ++level_num)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data, new_time);
    }

    // Initialize IB data.
    d_ib_method_ops->preprocessIntegrateData(current_time, new_time, num_cycles);

    // Initialize the fluid solver.
    d_ins_hier_integrator->preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    // At initial time interpolate Q.
    bool initial_time = MathUtilities<double>::equalEps(current_time, 0.0);
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
    std::vector<Eigen::Vector3d> U(num_objects), W(num_objects);
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
    IBHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Synchronize new state data.
    if (!skip_synchronize_new_state_data)
    {
        if (d_enable_logging)
            plog << d_object_name << "::postprocessIntegrateHierarchy(): synchronizing updated data\n";
        synchronizeHierarchyData(NEW_DATA);
    }

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                               d_ins_hier_integrator->getNewContext());

    // Determine the CFL number.
    double cfl_max = 0.0;
    PatchCellDataOpsReal<NDIM, double> patch_cc_ops;
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
            Pointer<CellData<NDIM, double> > u_cc_new_data = patch->getPatchData(u_new_idx);
            Pointer<SideData<NDIM, double> > u_sc_new_data = patch->getPatchData(u_new_idx);
            double u_max = 0.0;
            if (u_cc_new_data) u_max = patch_cc_ops.maxNorm(u_cc_new_data, patch_box);
            if (u_sc_new_data) u_max = patch_sc_ops.maxNorm(u_sc_new_data, patch_box);
            cfl_max = std::max(cfl_max, u_max * dt / dx_min);
        }
    }
    cfl_max = SAMRAI_MPI::maxReduction(cfl_max);
    d_regrid_cfl_estimate += cfl_max;
    if (d_enable_logging)
        plog << d_object_name << "::postprocessIntegrateHierarchy(): CFL number = " << cfl_max << "\n";
    if (d_enable_logging)
        plog << d_object_name
             << "::postprocessIntegrateHierarchy(): estimated upper bound on IB "
                "point displacement since last regrid = "
             << d_regrid_cfl_estimate << "\n";

    // Deallocate the fluid solver.
    d_ins_hier_integrator->postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Deallocate IB data.
    d_ib_method_ops->postprocessIntegrateData(current_time, new_time, num_cycles);

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

    // Finish initializing the hierarchy integrator.  This function call should
    // come at the end of this function.
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
