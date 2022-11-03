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

#include "ibamr/IBFEDirectForcingKinematics.h"
#include "ibamr/IBFEPostProcessor.h"

#include "ibtk/FEDataManager.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/libmesh_utilities.h"

#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/equation_systems.h"
#include "libmesh/system.h"

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

namespace libMesh
{
template <typename T>
class NumericVector;
} // namespace libMesh

using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEPostProcessor::IBFEPostProcessor(std::string name, FEDataManager* fe_data_manager)
    : d_name(std::move(name)),
      d_mesh(&fe_data_manager->getEquationSystems()->get_mesh()),
      d_fe_data_manager(fe_data_manager)
{
    // intentionally blank
    return;
} // IBFEPostProcessor

void
IBFEPostProcessor::registerScalarVariable(const std::string& name,
                                          libMesh::FEFamily fe_family,
                                          libMesh::Order fe_order,
                                          ScalarMeshFcnPtr fcn,
                                          const std::vector<SystemData>& system_data,
                                          void* fcn_ctx)
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    auto& system = equation_systems->add_system<System>(name + " reconstruction system");
    RestartManager* restart_manager = RestartManager::getManager();
    const bool is_from_restart = restart_manager->isFromRestart();
    if (!is_from_restart) system.add_variable(name, fe_order, fe_family);
    d_scalar_var_systems.push_back(&system);
    d_scalar_var_fcns.push_back(fcn);
    d_scalar_var_system_data.push_back(system_data);
    d_scalar_var_fcn_ctxs.push_back(fcn_ctx);
    d_var_systems.push_back(&system);
    return;
} // registerScalarVariable

void
IBFEPostProcessor::registerVectorVariable(const std::string& name,
                                          libMesh::FEFamily fe_family,
                                          libMesh::Order fe_order,
                                          VectorMeshFcnPtr fcn,
                                          const std::vector<SystemData>& system_data,
                                          void* fcn_ctx,
                                          unsigned int dim)
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    auto& system = equation_systems->add_system<System>(name + " reconstruction system");
    RestartManager* restart_manager = RestartManager::getManager();
    const bool is_from_restart = restart_manager->isFromRestart();
    for (unsigned int i = 0; i < dim; ++i)
    {
        if (!is_from_restart) system.add_variable(name + "_" + std::to_string(i), fe_order, fe_family);
    }
    d_vector_var_systems.push_back(&system);
    d_vector_var_fcns.push_back(fcn);
    d_vector_var_system_data.push_back(system_data);
    d_vector_var_fcn_ctxs.push_back(fcn_ctx);
    d_vector_var_dims.push_back(dim);
    d_var_systems.push_back(&system);
    return;
} // registerVectorVariable

void
IBFEPostProcessor::registerTensorVariable(const std::string& var_name,
                                          libMesh::FEFamily var_fe_family,
                                          libMesh::Order var_fe_order,
                                          TensorMeshFcnPtr var_fcn,
                                          const std::vector<SystemData>& system_data,
                                          void* var_fcn_ctx,
                                          unsigned int var_dim)
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    auto& system = equation_systems->add_system<System>(var_name + " reconstruction system");
    RestartManager* restart_manager = RestartManager::getManager();
    const bool is_from_restart = restart_manager->isFromRestart();
    for (unsigned int i = 0; i < var_dim; ++i)
    {
        for (unsigned int j = 0; j < var_dim; ++j)
        {
            if (!is_from_restart)
                system.add_variable(
                    var_name + "_" + std::to_string(i) + std::to_string(j), var_fe_order, var_fe_family);
        }
    }
    d_tensor_var_systems.push_back(&system);
    d_tensor_var_fcns.push_back(var_fcn);
    d_tensor_var_system_data.push_back(system_data);
    d_tensor_var_fcn_ctxs.push_back(var_fcn_ctx);
    d_tensor_var_dims.push_back(var_dim);
    d_var_systems.push_back(&system);
    return;
} // registerTensorVariable

void
IBFEPostProcessor::registerInterpolatedScalarEulerianVariable(
    const std::string& var_name,
    libMesh::FEFamily var_fe_family,
    libMesh::Order var_fe_order,
    Pointer<hier::Variable<NDIM> > var,
    Pointer<VariableContext> ctx,
    const HierarchyGhostCellInterpolation::InterpolationTransactionComponent& ghost_fill_transaction)
{
    registerInterpolatedScalarEulerianVariable(var_name,
                                               var_fe_family,
                                               var_fe_order,
                                               var,
                                               ctx,
                                               ghost_fill_transaction,
                                               d_fe_data_manager->getDefaultInterpSpec());
    return;
} // registerInterpolatedScalarEulerianVariable

void
IBFEPostProcessor::registerInterpolatedScalarEulerianVariable(
    const std::string& var_name,
    libMesh::FEFamily var_fe_family,
    libMesh::Order var_fe_order,
    Pointer<hier::Variable<NDIM> > var,
    Pointer<VariableContext> ctx,
    const HierarchyGhostCellInterpolation::InterpolationTransactionComponent& ghost_fill_transaction,
    const FEDataManager::InterpSpec& interp_spec)
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    auto& system = equation_systems->add_system<System>(var_name + " interpolation system");
    RestartManager* restart_manager = RestartManager::getManager();
    const bool is_from_restart = restart_manager->isFromRestart();
    if (!is_from_restart) system.add_variable(var_name, var_fe_order, var_fe_family);
    d_scalar_interp_var_systems.push_back(&system);
    d_scalar_interp_vars.push_back(var);
    d_scalar_interp_ctxs.push_back(ctx);
    d_scalar_interp_data_idxs.push_back(-1); // These must be set up just before they are used.
    d_scalar_interp_scratch_idxs.push_back(-1);
    d_scalar_interp_fill_transactions.push_back(ghost_fill_transaction);
    d_scalar_interp_specs.push_back(interp_spec);
    d_var_systems.push_back(&system);
} // registerInterpolatedEulerianScalarVariable

void
IBFEPostProcessor::initializeFEData()
{
    if (d_fe_data_initialized) return;
    for (const auto& var_system : d_var_systems)
    {
        System& system = *var_system;
        system.assemble_before_solve = false;
        system.assemble();
    }
    d_fe_data_initialized = true;
    return;
} // initializeFEData

void
IBFEPostProcessor::postProcessData(const double data_time)
{
    // First interpolate variables from the Eulerian grid, then reconstruct
    // variables on the Lagrangian mesh.
    interpolateVariables(data_time);
    reconstructVariables(data_time);
    return;
} // postProcessData

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBFEPostProcessor::interpolateVariables(const double data_time)
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = d_fe_data_manager->getPatchHierarchy();
    const int coarsest_ln = d_fe_data_manager->getCoarsestPatchLevelNumber();
    const int finest_ln = d_fe_data_manager->getFinestPatchLevelNumber();

    const size_t num_eulerian_vars = d_scalar_interp_var_systems.size();

    // Set up Eulerian scratch space and fill ghost cell values.
    std::set<int> scratch_idxs;
    for (unsigned int k = 0; k < num_eulerian_vars; ++k)
    {
        int& data_idx = d_scalar_interp_data_idxs[k];
        int& scratch_idx = d_scalar_interp_scratch_idxs[k];
        if (data_idx < 0 || scratch_idx < 0)
        {
            TBOX_ASSERT(data_idx < 0 || scratch_idx < 0);
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            Pointer<hier::Variable<NDIM> > data_var = d_scalar_interp_vars[k];
            Pointer<VariableContext> data_ctx = d_scalar_interp_ctxs[k];
            data_idx = var_db->mapVariableAndContextToIndex(data_var, data_ctx);
            TBOX_ASSERT(data_idx >= 0);
            Pointer<VariableContext> scratch_ctx = var_db->getContext(d_name + "::SCRATCH");
            const FEDataManager::InterpSpec& interp_spec = d_scalar_interp_specs[k];
            const int ghost_width = LEInteractor::getMinimumGhostWidth(interp_spec.kernel_fcn) + 1;
            scratch_idx = var_db->registerVariableAndContext(data_var, scratch_ctx, ghost_width);
            scratch_idxs.insert(scratch_idx);
            d_scalar_interp_fill_transactions[k].d_src_data_idx = data_idx;
            d_scalar_interp_fill_transactions[k].d_dst_data_idx = scratch_idx;
        }
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (unsigned int k = 0; k < num_eulerian_vars; ++k)
        {
            const int scratch_idx = d_scalar_interp_scratch_idxs[k];
            if (!level->checkAllocated(scratch_idx)) level->allocatePatchData(scratch_idx, data_time);
        }
    }

    HierarchyGhostCellInterpolation ghost_fill_op;
    ghost_fill_op.initializeOperatorState(d_scalar_interp_fill_transactions, hierarchy);
    ghost_fill_op.fillData(data_time);

    // Interpolate variables.
    std::unique_ptr<libMesh::PetscVector<double> > X_ghost_vec =
        d_fe_data_manager->buildIBGhostedVector(d_fe_data_manager->getCurrentCoordinatesSystemName());
    copy_and_synch(*d_fe_data_manager->getSolutionVector(d_fe_data_manager->getCurrentCoordinatesSystemName()),
                   *X_ghost_vec);
    for (unsigned int k = 0; k < num_eulerian_vars; ++k)
    {
        System* system = d_scalar_interp_var_systems[k];
        const std::string& system_name = system->name();
        const int scratch_idx = d_scalar_interp_scratch_idxs[k];
        d_fe_data_manager->interp(scratch_idx, *system->solution, *X_ghost_vec, system_name, d_scalar_interp_specs[k]);
    }

    // Deallocate Eulerian scratch space.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (unsigned int k = 0; k < num_eulerian_vars; ++k)
        {
            const int scratch_idx = d_scalar_interp_scratch_idxs[k];
            level->deallocatePatchData(scratch_idx);
        }
    }
    return;
} // interpolateVariables

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
