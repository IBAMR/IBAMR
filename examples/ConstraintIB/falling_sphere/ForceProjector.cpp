// Filename: ForceProjector.cpp
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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
// POSSIBILITY OF SUCH DAMAGE

/////////////////////////////////////// INCLUDES ///////////////////////////////////////////

#include "ForceProjector.h"

#include "Box.h"
#include "Patch.h"
#include "PatchLevel.h"
#include "VariableDatabase.h"
#include "ibamr/namespaces.h"
#include "tbox/Array.h"
#include "tbox/PIO.h"

namespace IBTK
{
void
callForceProjectorCallBackFunction(const double current_time, const double new_time, const int cycle_num, void* ctx)
{
    static ForceProjector* ptr_forceprojector = static_cast<ForceProjector*>(ctx);
    if (cycle_num == 0)
    {
        ptr_forceprojector->calculateLagrangianBodyForce(new_time, current_time);
        ptr_forceprojector->calculateEulerianBodyForce(new_time, current_time);
    }

    return;

} // callForceProjectorCallBackFunction

ForceProjector::ForceProjector(const std::string& object_name,
                               LDataManager* lag_data_manager,
                               Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                               Pointer<Database> input_db,
                               const std::string solver_type)
    : d_object_name(object_name),
      d_lag_data_manager(lag_data_manager),
      d_patch_hierarchy(patch_hierarchy),
      d_solver_type(solver_type),
      d_grav_const(NDIM)
{
    // put some default values.
    d_rho_fluid = 1.0;
    d_rho_body = 1.0;

    // Initialize  variables & variable contexts associated with Eulerian forces.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_body_force_context = var_db->getContext(d_object_name + "::BODYFORCE");
    if (d_solver_type == "STAGGERED")
        d_body_force_var = new SideVariable<NDIM, double>(d_object_name + "::BodyForce_sc_var");
    if (d_solver_type == "COLLOCATED")
        d_body_force_var = new CellVariable<NDIM, double>(d_object_name + "::BodyForce_cc_var", NDIM);
    d_body_force_idx = var_db->registerVariableAndContext(d_body_force_var, d_body_force_context, 0);

    getFromInput(input_db);

    return;

} // ForceProjector

ForceProjector::~ForceProjector()
{
    // intentionally left blank
    return;

} //~ForceProjector

void
ForceProjector::getFromInput(Pointer<Database> input_db)
{
    d_rho_fluid = input_db->getDoubleWithDefault("rho_fluid", d_rho_fluid);
    d_rho_body = input_db->getDoubleWithDefault("rho_body", d_rho_body);
    d_grav_const = input_db->getDoubleArray("gravitational_constant");

    return;

} // getFromInput

void
ForceProjector::registerLagrangianQuantityName(const std::string& lag_quantity_name)
{
    registerLagrangianQuantitiesName(std::vector<std::string>(1, lag_quantity_name));

    return;

} // registerLagrangianQuantityName

void
ForceProjector::registerLagrangianQuantitiesName(const std::vector<std::string>& lag_quantities_name)
{
    const unsigned size = lag_quantities_name.size();
    for (unsigned i = 0; i < size; ++i)
    {
        d_lag_quantities_name.push_back(lag_quantities_name[i]);
    }

    return;

} // registerLagrangianQuantitiesName

void
ForceProjector::associateVolumeElement(const double vol_lag_pt)
{
    d_vol_lag_pt = vol_lag_pt;

    return;

} // associateVolumeElement

void
ForceProjector::calculateLagrangianBodyForce(const double /*new_time*/, const double /*current_time*/)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    d_lag_force.clear();
    d_lag_force.resize(finest_ln + 1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_lag_data_manager->levelContainsLagrangianData(ln)) continue;
        d_lag_force[ln] = d_lag_data_manager->createLData(d_object_name + "::lag_force_data", ln, NDIM, false);
    }

    // Set Lagrangian gravitational force.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_lag_data_manager->levelContainsLagrangianData(ln)) continue;

        // Get ponter to LData corresponding to lagrangian force.
        boost::multi_array_ref<double, 2>& F_data = *d_lag_force[ln]->getLocalFormVecArray();
        const Pointer<LMesh> mesh = d_lag_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int local_idx = node_idx->getLocalPETScIndex();
            double* const F = &F_data[local_idx][0];

            for (int d = 0; d < 3; ++d) F[d] = d_rho_fluid * d_grav_const[d] * d_vol_lag_pt;
        }
        d_lag_force[ln]->restoreArrays();

    } // all levels

    return;

} // calculateLagrangianBodyForce

void
ForceProjector::calculateEulerianBodyForce(const double /*new_time*/, const double current_time)
{
    // allocate patch data for Eulerian forcing.
    const int coarsest_ln = 0;
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_body_force_idx)) level->deallocatePatchData(d_body_force_idx);
        level->allocatePatchData(d_body_force_idx, current_time);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            if (d_solver_type == "STAGGERED")
            {
                Pointer<SideData<NDIM, double> > body_force_data = patch->getPatchData(d_body_force_idx);
                body_force_data->fill(0.0);
            }
            else if (d_solver_type == "COLLOCATED")
            {
                Pointer<CellData<NDIM, double> > body_force_data = patch->getPatchData(d_body_force_idx);
                body_force_data->fill(0.0);
            }
            else
            {
                TBOX_ERROR("ForceProjector::calculateEulerianBodyForce() "
                           << "UNKNOWN SOLVER ENCOUNTERED" << std::endl);
            }

        } // iterate over patches

    } // all levels.

    // spread the lagrangian force from finest level to the finest level.
    std::vector<Pointer<LData> > F_data(finest_ln + 1, Pointer<LData>(NULL));
    std::vector<Pointer<LData> > X_data(finest_ln + 1, Pointer<LData>(NULL));

    // Fill in the above vectors at the finest level.
    F_data[finest_ln] = d_lag_force[finest_ln];
    X_data[finest_ln] = d_lag_data_manager->getLData("X", finest_ln);

    // Spread the deformation velocities.
    d_lag_data_manager->spread(d_body_force_idx, F_data, X_data, (RobinPhysBdryPatchStrategy*)NULL);

    return;

} // calculateEulerianBodyForce

} // namespace IBTK
