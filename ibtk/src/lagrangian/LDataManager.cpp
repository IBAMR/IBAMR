// Filename: LDataManager.cpp
// Created on 01 Mar 2004 by Boyce Griffith
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

#include <math.h>
#include <stddef.h>
#include <algorithm>
#include <limits>
#include <map>
#include <numeric>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxArray.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxTree.h"
#include "SAMRAI/geom/CartesianCellDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/HierarchyDataOpsManager.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"

#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/ProcessorMapping.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "boost/array.hpp"
#include "boost/multi_array.hpp"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/LIndexSetData.h"
#include "ibtk/LInitStrategy.h"
#include "ibtk/LMesh.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeIndex.h"
#include "ibtk/LNodeSet.h"
#include "ibtk/LNodeSetData.h"
#include "ibtk/LNodeSetVariable.h"
#include "ibtk/LNodeTransaction.h"
#include "ibtk/LSet.h"
#include "ibtk/LSetData.h"
#include "ibtk/LSetDataIterator.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/LTransaction.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/compiler_hints.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "petscao.h"
#include "petscis.h"
#include "petscsys.h"
#include "petscvec.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Schedule.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Transaction.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{

class Index;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static boost::shared_ptr<Timer> t_spread;
static boost::shared_ptr<Timer> t_interp;
static boost::shared_ptr<Timer> t_map_lagrangian_to_petsc;
static boost::shared_ptr<Timer> t_map_petsc_to_lagrangian;
static boost::shared_ptr<Timer> t_begin_data_redistribution;
static boost::shared_ptr<Timer> t_end_data_redistribution;
static boost::shared_ptr<Timer> t_update_workload_estimates;
static boost::shared_ptr<Timer> t_update_node_count_data;
static boost::shared_ptr<Timer> t_initialize_level_data;
static boost::shared_ptr<Timer> t_reset_hierarchy_configuration;
static boost::shared_ptr<Timer> t_apply_gradient_detector;
static boost::shared_ptr<Timer> t_put_to_database;
static boost::shared_ptr<Timer> t_begin_nonlocal_data_fill;
static boost::shared_ptr<Timer> t_end_nonlocal_data_fill;
static boost::shared_ptr<Timer> t_compute_node_distribution;
static boost::shared_ptr<Timer> t_compute_node_offsets;

// Assume max(U)dt/dx <= 2.
static const int CFL_WIDTH = 2;

// Version of LDataManager restart file data.
static const int LDATA_MANAGER_VERSION = 1;

inline int round(double x)
{
    return floor(x + 0.5);
} // round
}

const std::string LDataManager::POSN_DATA_NAME = "X";
const std::string LDataManager::INIT_POSN_DATA_NAME = "X0";
const std::string LDataManager::VEL_DATA_NAME = "U";
std::map<std::string, LDataManager*> LDataManager::s_data_manager_instances;
bool LDataManager::s_registered_callback = false;
unsigned char LDataManager::s_shutdown_priority = 200;
std::vector<int> LDataManager::s_ao_dummy(1, -1);

LDataManager* LDataManager::getManager(const std::string& name,
                                       const std::string& default_interp_kernel_fcn,
                                       const std::string& default_spread_kernel_fcn,
                                       const IntVector& min_ghost_width,
                                       bool register_for_restart)
{
    if (s_data_manager_instances.find(name) == s_data_manager_instances.end())
    {
        const IntVector ghost_width =
            IntVector::max(min_ghost_width,
                           IntVector(DIM,
                                     std::max(LEInteractor::getMinimumGhostWidth(default_interp_kernel_fcn),
                                              LEInteractor::getMinimumGhostWidth(default_spread_kernel_fcn))));
        s_data_manager_instances[name] = boost::make_shared<LDataManager>(
            name, default_interp_kernel_fcn, default_spread_kernel_fcn, ghost_width, register_for_restart);
    }
    if (!s_registered_callback)
    {
        static StartupShutdownManager::Handler handler(NULL, NULL, freeAllManagers, NULL, s_shutdown_priority);
        StartupShutdownManager::registerHandler(&handler);
        s_registered_callback = true;
    }
    return s_data_manager_instances[name];
} // getManager

void LDataManager::freeAllManagers()
{
    for (auto it = s_data_manager_instances.begin();
         it != s_data_manager_instances.end();
         ++it)
    {
        if (it->second)
        {
            delete it->second;
        }
        it->second = NULL;
    }
    return;
} // freeManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

void LDataManager::setPatchHierarchy(boost::shared_ptr<PatchHierarchy> hierarchy)
{
    TBOX_ASSERT(hierarchy);

    // Reset the hierarchy.
    d_hierarchy = hierarchy;
    d_grid_geom = hierarchy->getGridGeometry();
    return;
} // setPatchHierarchy

boost::shared_ptr<PatchHierarchy> LDataManager::getPatchHierarchy() const
{
    return d_hierarchy;
} // getPatchHierarchy

void LDataManager::setPatchLevels(const int coarsest_ln, const int finest_ln)
{
    TBOX_ASSERT(d_hierarchy);
    TBOX_ASSERT((coarsest_ln >= 0) && (finest_ln >= coarsest_ln) && (finest_ln <= d_hierarchy->getFinestLevelNumber()));

    // Destroy any unneeded AO objects.
    int ierr;
    for (int level_number = std::max(d_coarsest_ln, 0); (level_number <= d_finest_ln) && (level_number < coarsest_ln);
         ++level_number)
    {
        if (d_ao[level_number])
        {
            ierr = AODestroy(&d_ao[level_number]);
            IBTK_CHKERRQ(ierr);
        }
    }
    for (int level_number = finest_ln + 1; level_number <= d_finest_ln; ++level_number)
    {
        if (d_ao[level_number])
        {
            ierr = AODestroy(&d_ao[level_number]);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Reset the level numbers.
    d_coarsest_ln = coarsest_ln;
    d_finest_ln = finest_ln;

    // Resize some arrays.
    d_level_contains_lag_data.resize(d_finest_ln + 1);
    d_strct_name_to_strct_id_map.resize(d_finest_ln + 1);
    d_strct_id_to_strct_name_map.resize(d_finest_ln + 1);
    d_strct_id_to_lag_idx_range_map.resize(d_finest_ln + 1);
    d_last_lag_idx_to_strct_id_map.resize(d_finest_ln + 1);
    d_inactive_strcts.resize(d_finest_ln + 1);
    d_lag_mesh.resize(d_finest_ln + 1);
    d_lag_mesh_data.resize(d_finest_ln + 1);
    d_needs_synch.resize(d_finest_ln + 1, false);
    d_ao.resize(d_finest_ln + 1);
    d_num_nodes.resize(d_finest_ln + 1);
    d_node_offset.resize(d_finest_ln + 1);
    d_local_lag_indices.resize(d_finest_ln + 1);
    d_nonlocal_lag_indices.resize(d_finest_ln + 1);
    d_local_petsc_indices.resize(d_finest_ln + 1);
    d_nonlocal_petsc_indices.resize(d_finest_ln + 1);
    return;
} // setPatchLevels

std::pair<int, int> LDataManager::getPatchLevels() const
{
    return std::make_pair(d_coarsest_ln, d_finest_ln + 1);
} // getPatchLevels

void LDataManager::spread(const int f_data_idx,
                          boost::shared_ptr<LData> F_data,
                          boost::shared_ptr<LData> X_data,
                          boost::shared_ptr<LData> ds_data,
                          RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                          const int level_num,
                          const std::vector<boost::shared_ptr<RefineSchedule> >& f_prolongation_scheds,
                          const double fill_data_time,
                          const bool F_data_ghost_node_update,
                          const bool X_data_ghost_node_update,
                          const bool ds_data_ghost_node_update)
{
    spread(f_data_idx,
           F_data,
           X_data,
           ds_data,
           d_default_spread_kernel_fcn,
           f_phys_bdry_op,
           level_num,
           f_prolongation_scheds,
           fill_data_time,
           F_data_ghost_node_update,
           X_data_ghost_node_update,
           ds_data_ghost_node_update);
    return;
} // spread

void LDataManager::spread(const int f_data_idx,
                          boost::shared_ptr<LData> F_data,
                          boost::shared_ptr<LData> X_data,
                          boost::shared_ptr<LData> ds_data,
                          const std::string& spread_kernel_fcn,
                          RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                          const int level_num,
                          const std::vector<boost::shared_ptr<RefineSchedule> >& f_prolongation_scheds,
                          const double fill_data_time,
                          const bool F_data_ghost_node_update,
                          const bool X_data_ghost_node_update,
                          const bool ds_data_ghost_node_update)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    TBOX_ASSERT(coarsest_ln <= level_num && level_num <= finest_ln);
    std::vector<boost::shared_ptr<LData> > F_data_vec(finest_ln + 1);
    std::vector<boost::shared_ptr<LData> > X_data_vec(finest_ln + 1);
    std::vector<boost::shared_ptr<LData> > ds_data_vec(finest_ln + 1);
    F_data_vec[level_num] = F_data;
    X_data_vec[level_num] = X_data;
    ds_data_vec[level_num] = ds_data;
    spread(f_data_idx,
           F_data_vec,
           X_data_vec,
           ds_data_vec,
           spread_kernel_fcn,
           f_phys_bdry_op,
           f_prolongation_scheds,
           fill_data_time,
           F_data_ghost_node_update,
           X_data_ghost_node_update,
           ds_data_ghost_node_update,
           coarsest_ln,
           finest_ln);
    return;
} // spread

void LDataManager::spread(const int f_data_idx,
                          std::vector<boost::shared_ptr<LData> >& F_data,
                          std::vector<boost::shared_ptr<LData> >& X_data,
                          std::vector<boost::shared_ptr<LData> >& ds_data,
                          RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                          const std::vector<boost::shared_ptr<RefineSchedule> >& f_prolongation_scheds,
                          const double fill_data_time,
                          const bool F_data_ghost_node_update,
                          const bool X_data_ghost_node_update,
                          const bool ds_data_ghost_node_update,
                          const int coarsest_ln,
                          const int finest_ln)
{
    spread(f_data_idx,
           F_data,
           X_data,
           ds_data,
           d_default_spread_kernel_fcn,
           f_phys_bdry_op,
           f_prolongation_scheds,
           fill_data_time,
           F_data_ghost_node_update,
           X_data_ghost_node_update,
           ds_data_ghost_node_update,
           coarsest_ln,
           finest_ln);
    return;
} // spread

void LDataManager::spread(const int f_data_idx,
                          std::vector<boost::shared_ptr<LData> >& F_data,
                          std::vector<boost::shared_ptr<LData> >& X_data,
                          std::vector<boost::shared_ptr<LData> >& ds_data,
                          const std::string& spread_kernel_fcn,
                          RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                          const std::vector<boost::shared_ptr<RefineSchedule> >& f_prolongation_scheds,
                          const double fill_data_time,
                          const bool F_data_ghost_node_update,
                          const bool X_data_ghost_node_update,
                          const bool ds_data_ghost_node_update,
                          const int coarsest_ln_in,
                          const int finest_ln_in)
{
    IBTK_TIMER_START(t_spread);

    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? d_hierarchy->getFinestLevelNumber() : finest_ln_in);

    // Zero inactivated components.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        zeroInactivatedComponents(F_data[ln], ln);
    }

    // Compute F*ds.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!levelContainsLagrangianData(ln)) continue;

        if (F_data_ghost_node_update) F_data[ln]->beginGhostUpdate();
        if (ds_data_ghost_node_update) ds_data[ln]->beginGhostUpdate();
    }
    std::vector<boost::shared_ptr<LData> > F_ds_data(F_data.size());
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!levelContainsLagrangianData(ln)) continue;

        if (F_data_ghost_node_update) F_data[ln]->endGhostUpdate();
        if (ds_data_ghost_node_update) ds_data[ln]->endGhostUpdate();

        const int depth = F_data[ln]->getDepth();
        F_ds_data[ln] = boost::make_shared<LData("", getNumberOfLocalNodes>(ln), depth, d_nonlocal_petsc_indices[ln]);
        boost::multi_array_ref<double, 2>& F_ds_arr = *F_ds_data[ln]->getGhostedLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& F_arr = *F_data[ln]->getGhostedLocalFormVecArray();
        const boost::multi_array_ref<double, 1>& ds_arr = *ds_data[ln]->getGhostedLocalFormArray();
        for (int k = 0; k < static_cast<int>(F_data[ln]->getLocalNodeCount() + F_data[ln]->getGhostNodeCount()); ++k)
        {
            for (int d = 0; d < depth; ++d)
            {
                F_ds_arr[k][d] = F_arr[k][d] * ds_arr[k];
            }
        }
        F_ds_data[ln]->restoreArrays();
        F_data[ln]->restoreArrays();
        ds_data[ln]->restoreArrays();
    }

    IBTK_TIMER_STOP(t_spread);

    // Spread data from the Lagrangian mesh to the Eulerian grid.
    spread(f_data_idx,
           F_ds_data,
           X_data,
           spread_kernel_fcn,
           f_phys_bdry_op,
           f_prolongation_scheds,
           fill_data_time,
           /*F_data_ghost_node_update*/ false,
           X_data_ghost_node_update,
           coarsest_ln,
           finest_ln);
    return;
} // spread

void LDataManager::spread(const int f_data_idx,
                          boost::shared_ptr<LData> F_data,
                          boost::shared_ptr<LData> X_data,
                          RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                          const int level_num,
                          const std::vector<boost::shared_ptr<RefineSchedule> >& f_prolongation_scheds,
                          const double fill_data_time,
                          const bool F_data_ghost_node_update,
                          const bool X_data_ghost_node_update)
{
    spread(f_data_idx,
           F_data,
           X_data,
           d_default_spread_kernel_fcn,
           f_phys_bdry_op,
           level_num,
           f_prolongation_scheds,
           fill_data_time,
           F_data_ghost_node_update,
           X_data_ghost_node_update);
    return;
} // spread

void LDataManager::spread(const int f_data_idx,
                          boost::shared_ptr<LData> F_data,
                          boost::shared_ptr<LData> X_data,
                          const std::string& spread_kernel_fcn,
                          RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                          const int level_num,
                          const std::vector<boost::shared_ptr<RefineSchedule> >& f_prolongation_scheds,
                          const double fill_data_time,
                          const bool F_data_ghost_node_update,
                          const bool X_data_ghost_node_update)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    TBOX_ASSERT(coarsest_ln <= level_num && level_num <= finest_ln);
    std::vector<boost::shared_ptr<LData> > F_data_vec(finest_ln + 1);
    std::vector<boost::shared_ptr<LData> > X_data_vec(finest_ln + 1);
    F_data_vec[level_num] = F_data;
    X_data_vec[level_num] = X_data;
    spread(f_data_idx,
           F_data_vec,
           X_data_vec,
           spread_kernel_fcn,
           f_phys_bdry_op,
           f_prolongation_scheds,
           fill_data_time,
           F_data_ghost_node_update,
           X_data_ghost_node_update,
           coarsest_ln,
           finest_ln);
    return;
} // spread

void LDataManager::spread(const int f_data_idx,
                          std::vector<boost::shared_ptr<LData> >& F_data,
                          std::vector<boost::shared_ptr<LData> >& X_data,
                          RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                          const std::vector<boost::shared_ptr<RefineSchedule> >& f_prolongation_scheds,
                          const double fill_data_time,
                          const bool F_data_ghost_node_update,
                          const bool X_data_ghost_node_update,
                          const int coarsest_ln,
                          const int finest_ln)
{
    spread(f_data_idx,
           F_data,
           X_data,
           d_default_spread_kernel_fcn,
           f_phys_bdry_op,
           f_prolongation_scheds,
           fill_data_time,
           F_data_ghost_node_update,
           X_data_ghost_node_update,
           coarsest_ln,
           finest_ln);
    return;
} // spread

void LDataManager::spread(const int f_data_idx,
                          std::vector<boost::shared_ptr<LData> >& F_data,
                          std::vector<boost::shared_ptr<LData> >& X_data,
                          const std::string& spread_kernel_fcn,
                          RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                          const std::vector<boost::shared_ptr<RefineSchedule> >& f_prolongation_scheds,
                          const double fill_data_time,
                          const bool F_data_ghost_node_update,
                          const bool X_data_ghost_node_update,
                          const int coarsest_ln_in,
                          const int finest_ln_in)
{
    IBTK_TIMER_START(t_spread);

    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? d_hierarchy->getFinestLevelNumber() : finest_ln_in);
    VariableDatabase* var_db = VariableDatabase::getDatabase();

    // Determine the type of data centering.
    boost::shared_ptr<Variable> f_var;
    var_db->mapIndexToVariable(f_data_idx, f_var);
    auto f_cc_var = boost::dynamic_pointer_cast<CellVariable<double> >(f_var);
    auto f_ec_var = boost::dynamic_pointer_cast<EdgeVariable<double> >(f_var);
    auto f_nc_var = boost::dynamic_pointer_cast<NodeVariable<double> >(f_var);
    auto f_sc_var = boost::dynamic_pointer_cast<SideVariable<double> >(f_var);
    const bool cc_data = f_cc_var;
    const bool ec_data = f_ec_var;
    const bool nc_data = f_nc_var;
    const bool sc_data = f_sc_var;
    TBOX_ASSERT(cc_data || ec_data || nc_data || sc_data);

    // Make a copy of the Eulerian data.
    const int f_copy_data_idx = var_db->registerClonedPatchDataIndex(f_var, f_data_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        auto level =d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(f_copy_data_idx);
    }
    boost::shared_ptr<HierarchyDataOpsReal<double> > f_data_ops =
        HierarchyDataOpsManager::getManager()->getOperationsDouble(f_var, d_hierarchy, true);
    f_data_ops->swapData(f_copy_data_idx, f_data_idx);
    f_data_ops->setToScalar(f_data_idx, 0.0, /*interior_only*/ false);

    // Start filling Lagrangian ghost node values.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!levelContainsLagrangianData(ln)) continue;

        if (F_data_ghost_node_update) F_data[ln]->beginGhostUpdate();
        if (X_data_ghost_node_update) X_data[ln]->beginGhostUpdate();
    }

    // Spread data from the Lagrangian mesh to the Eulerian grid.
    auto grid_geom = BOOST_CAST<CartesianGridGeometry>(d_hierarchy->getGridGeometry());
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        // If there are coarser levels in the patch hierarchy, prolong data from
        // the coarser levels before spreading data on this level.
        if (ln > coarsest_ln && ln < static_cast<int>(f_prolongation_scheds.size()) && f_prolongation_scheds[ln])
        {
            f_prolongation_scheds[ln]->fillData(fill_data_time);
        }

        if (!levelContainsLagrangianData(ln)) continue;

        // Spread data onto the grid.
        if (F_data_ghost_node_update) F_data[ln]->endGhostUpdate();
        if (X_data_ghost_node_update) X_data[ln]->endGhostUpdate();
        auto level =d_hierarchy->getPatchLevel(ln);
        const IntVector& periodic_shift = grid_geom->getPeriodicShift(level->getRatioToLevelZero());
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch =*p;
            boost::shared_ptr<PatchData> f_data = patch->getPatchData(f_data_idx);
            boost::shared_ptr<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
            const Box& box = idx_data->getGhostBox();
            if (cc_data)
            {
                boost::shared_ptr<CellData<double> > f_cc_data = f_data;
                LEInteractor::spread(
                    f_cc_data, F_data[ln], X_data[ln], idx_data, patch, box, periodic_shift, spread_kernel_fcn);
            }
            if (ec_data)
            {
                boost::shared_ptr<EdgeData<double> > f_ec_data = f_data;
                LEInteractor::spread(
                    f_ec_data, F_data[ln], X_data[ln], idx_data, patch, box, periodic_shift, spread_kernel_fcn);
            }
            if (nc_data)
            {
                boost::shared_ptr<NodeData<double> > f_nc_data = f_data;
                LEInteractor::spread(
                    f_nc_data, F_data[ln], X_data[ln], idx_data, patch, box, periodic_shift, spread_kernel_fcn);
            }
            if (sc_data)
            {
                boost::shared_ptr<SideData<double> > f_sc_data = f_data;
                LEInteractor::spread(
                    f_sc_data, F_data[ln], X_data[ln], idx_data, patch, box, periodic_shift, spread_kernel_fcn);
            }
            if (f_phys_bdry_op)
            {
                f_phys_bdry_op->setPatchDataIndex(f_data_idx);
                f_phys_bdry_op->accumulateFromPhysicalBoundaryData(*patch, fill_data_time, f_data->getGhostCellWidth());
            }
        }
    }

    // Accumulate data.
    f_data_ops->swapData(f_copy_data_idx, f_data_idx);
    f_data_ops->add(f_data_idx, f_data_idx, f_copy_data_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        auto level =d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(f_copy_data_idx);
    }
    var_db->removePatchDataIndex(f_copy_data_idx);

    IBTK_TIMER_STOP(t_spread);
    return;
} // spread

void LDataManager::interp(const int f_data_idx,
                          boost::shared_ptr<LData> F_data,
                          boost::shared_ptr<LData> X_data,
                          const int level_num,
                          const std::vector<boost::shared_ptr<CoarsenSchedule> >& f_synch_scheds,
                          const std::vector<boost::shared_ptr<RefineSchedule> >& f_ghost_fill_scheds,
                          const double fill_data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    TBOX_ASSERT(coarsest_ln <= level_num && level_num <= finest_ln);
    std::vector<boost::shared_ptr<LData> > F_data_vec(finest_ln + 1);
    std::vector<boost::shared_ptr<LData> > X_data_vec(finest_ln + 1);
    F_data_vec[level_num] = F_data;
    X_data_vec[level_num] = X_data;
    interp(f_data_idx,
           F_data_vec,
           X_data_vec,
           f_synch_scheds,
           f_ghost_fill_scheds,
           fill_data_time,
           coarsest_ln,
           finest_ln);
    return;
} // interp

void LDataManager::interp(const int f_data_idx,
                          std::vector<boost::shared_ptr<LData> >& F_data,
                          std::vector<boost::shared_ptr<LData> >& X_data,
                          const std::vector<boost::shared_ptr<CoarsenSchedule> >& f_synch_scheds,
                          const std::vector<boost::shared_ptr<RefineSchedule> >& f_ghost_fill_scheds,
                          const double fill_data_time,
                          const int coarsest_ln_in,
                          const int finest_ln_in)
{
    IBTK_TIMER_START(t_interp);

    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? d_hierarchy->getFinestLevelNumber() : finest_ln_in);
    VariableDatabase* var_db = VariableDatabase::getDatabase();

    // Determine the type of data centering.
    boost::shared_ptr<Variable> f_var;
    var_db->mapIndexToVariable(f_data_idx, f_var);
    auto f_cc_var = boost::dynamic_pointer_cast<CellVariable<double> >(f_var);
    auto f_ec_var = boost::dynamic_pointer_cast<EdgeVariable<double> >(f_var);
    auto f_nc_var = boost::dynamic_pointer_cast<NodeVariable<double> >(f_var);
    auto f_sc_var = boost::dynamic_pointer_cast<SideVariable<double> >(f_var);
    const bool cc_data = f_cc_var;
    const bool ec_data = f_ec_var;
    const bool nc_data = f_nc_var;
    const bool sc_data = f_sc_var;
    TBOX_ASSERT(cc_data || ec_data || nc_data || sc_data);

    // Synchronize Eulerian values.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        if (ln < static_cast<int>(f_synch_scheds.size()) && f_synch_scheds[ln])
        {
            f_synch_scheds[ln]->coarsenData();
        }
    }

    // Interpolate data from the Eulerian grid to the Lagrangian mesh.
    auto grid_geom = BOOST_CAST<CartesianGridGeometry>(d_hierarchy->getGridGeometry());
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!levelContainsLagrangianData(ln)) continue;

        if (ln < static_cast<int>(f_ghost_fill_scheds.size()) && f_ghost_fill_scheds[ln])
        {
            f_ghost_fill_scheds[ln]->fillData(fill_data_time);
        }
        auto level =d_hierarchy->getPatchLevel(ln);
        const IntVector& periodic_shift = grid_geom->getPeriodicShift(level->getRatioToLevelZero());
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch =*p;
            boost::shared_ptr<PatchData> f_data = patch->getPatchData(f_data_idx);
            boost::shared_ptr<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
            const Box& box = idx_data->getBox();
            if (cc_data)
            {
                boost::shared_ptr<CellData<double> > f_cc_data = f_data;
                LEInteractor::interpolate(F_data[ln],
                                          X_data[ln],
                                          idx_data,
                                          f_cc_data,
                                          patch,
                                          box,
                                          periodic_shift,
                                          d_default_interp_kernel_fcn);
            }
            if (ec_data)
            {
                boost::shared_ptr<EdgeData<double> > f_ec_data = f_data;
                LEInteractor::interpolate(F_data[ln],
                                          X_data[ln],
                                          idx_data,
                                          f_ec_data,
                                          patch,
                                          box,
                                          periodic_shift,
                                          d_default_interp_kernel_fcn);
            }
            if (nc_data)
            {
                boost::shared_ptr<NodeData<double> > f_nc_data = f_data;
                LEInteractor::interpolate(F_data[ln],
                                          X_data[ln],
                                          idx_data,
                                          f_nc_data,
                                          patch,
                                          box,
                                          periodic_shift,
                                          d_default_interp_kernel_fcn);
            }
            if (sc_data)
            {
                boost::shared_ptr<SideData<double> > f_sc_data = f_data;
                LEInteractor::interpolate(F_data[ln],
                                          X_data[ln],
                                          idx_data,
                                          f_sc_data,
                                          patch,
                                          box,
                                          periodic_shift,
                                          d_default_interp_kernel_fcn);
            }
        }
    }

    // Zero inactivated components.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        zeroInactivatedComponents(F_data[ln], ln);
    }

    IBTK_TIMER_STOP(t_interp);
    return;
} // interp

void LDataManager::registerLInitStrategy(boost::shared_ptr<LInitStrategy> lag_init)
{
    TBOX_ASSERT(lag_init);
    d_lag_init = lag_init;
    return;
} // registerLInitStrategy

void LDataManager::freeLInitStrategy()
{
    d_lag_init.reset();
    return;
} // freeLInitStrategy

void LDataManager::registerVisItDataWriter(boost::shared_ptr<VisItDataWriter> visit_writer)
{
    TBOX_ASSERT(visit_writer);
    d_visit_writer = visit_writer;
    if (d_output_workload)
    {
        d_visit_writer->registerPlotQuantity("workload estimate", "SCALAR", d_workload_idx);
    }
    if (d_output_node_count)
    {
        d_visit_writer->registerPlotQuantity("node count", "SCALAR", d_node_count_idx);
    }
    return;
} // registerVisItDataWriter

void LDataManager::registerLSiloDataWriter(boost::shared_ptr<LSiloDataWriter> silo_writer)
{
    TBOX_ASSERT(silo_writer);
    d_silo_writer = silo_writer;
    return;
} // registerLSiloDataWriter

void LDataManager::registerLoadBalancer(boost::shared_ptr<ChopAndPackLoadBalancer> load_balancer, int workload_idx)
{
    TBOX_ASSERT(load_balancer);
    d_load_balancer = load_balancer;
    d_workload_idx = workload_idx;
    boost::shared_ptr<Variable> workload_var;
    VariableDatabase::getDatabase()->mapIndexToVariable(d_workload_idx, workload_var);
    d_workload_var = workload_var;
    TBOX_ASSERT(d_workload_var);
    return;
} // return

boost::shared_ptr<LData> LDataManager::createLData(const std::string& quantity_name,
                                         const int level_number,
                                         const unsigned int depth,
                                         const bool maintain_data)
{
    TBOX_ASSERT(!maintain_data ||
                (d_lag_mesh_data[level_number].find(quantity_name) == d_lag_mesh_data[level_number].end()));
    TBOX_ASSERT(level_number >= 0);
    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);
    TBOX_ASSERT(depth > 0);
    auto ret_val = boost::make_shared<LData>(quantity_name, getNumberOfLocalNodes(level_number), depth, d_nonlocal_petsc_indices[level_number]);
    if (maintain_data)
    {
        d_lag_mesh_data[level_number][quantity_name] = ret_val;
    }
    return ret_val;
} // createLData

Point LDataManager::computeLagrangianStructureCenterOfMass(const int structure_id, const int level_number)
{
    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);
    int node_counter = 0;
    Point X_com(Point::Zero());
    std::pair<int, int> lag_idx_range = getLagrangianStructureIndexRange(structure_id, level_number);

    const boost::multi_array_ref<double, 2>& X_data =
        *d_lag_mesh_data[level_number][POSN_DATA_NAME]->getLocalFormVecArray();
    const auto mesh = getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    for (auto cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const int lag_idx = node_idx->getLagrangianIndex();
        if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
        {
            ++node_counter;
            const int local_idx = node_idx->getLocalPETScIndex();
            const double* const X = &X_data[local_idx][0];
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_com[d] += X[d];
            }
        }
    }
    d_lag_mesh_data[level_number][POSN_DATA_NAME]->restoreArrays();

    tbox::SAMRAI_MPI comm(MPI_COMM_WORLD);
    comm.AllReduce(&X_com[0], NDIM, MPI_SUM);
    comm.AllReduce(&node_counter, 1, MPI_SUM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        X_com[d] /= static_cast<double>(node_counter);
    }
    return X_com;
} // computeLagrangianStructureCenterOfMass

std::pair<Point, Point> LDataManager::computeLagrangianStructureBoundingBox(const int structure_id,
                                                                            const int level_number)
{
    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);
    Point X_lower(Point::Constant((std::numeric_limits<double>::max() - sqrt(std::numeric_limits<double>::epsilon()))));
    Point X_upper(
        Point::Constant(-(std::numeric_limits<double>::max() - sqrt(std::numeric_limits<double>::epsilon()))));
    std::pair<int, int> lag_idx_range = getLagrangianStructureIndexRange(structure_id, level_number);

    const boost::multi_array_ref<double, 2>& X_data =
        *d_lag_mesh_data[level_number][POSN_DATA_NAME]->getLocalFormVecArray();
    const auto mesh = getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    for (auto cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const int lag_idx = node_idx->getLagrangianIndex();
        if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
        {
            const int local_idx = node_idx->getLocalPETScIndex();
            const double* const X = &X_data[local_idx][0];
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_lower[d] = std::min(X_lower[d], X[d]);
                X_upper[d] = std::max(X_upper[d], X[d]);
            }
        }
    }
    d_lag_mesh_data[level_number][POSN_DATA_NAME]->restoreArrays();

    tbox::SAMRAI_MPI comm(MPI_COMM_WORLD);
    comm.AllReduce(&X_lower[0], NDIM, MPI_MIN);
    comm.AllReduce(&X_upper[0], NDIM, MPI_MAX);
    return std::make_pair(X_lower, X_upper);
} // computeLagrangianStructureBoundingBox

void LDataManager::activateLagrangianStructures(const std::vector<int>& structure_ids, const int level_number)
{
    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);
#ifndef NDEBUG
    tbox::SAMRAI_MPI comm(MPI_COMM_WORLD);
    comm.Barrier();
#endif
    for (auto cit = structure_ids.begin(); cit != structure_ids.end(); ++cit)
    {
        d_inactive_strcts[level_number].erase(*cit);
    }
    return;
} // activateLagrangianStructures

void LDataManager::inactivateLagrangianStructures(const std::vector<int>& structure_ids, const int level_number)
{
    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);
#ifndef NDEBUG
    tbox::SAMRAI_MPI comm(MPI_COMM_WORLD);
    comm.Barrier();
#endif
    for (auto cit = structure_ids.begin(); cit != structure_ids.end(); ++cit)
    {
        d_inactive_strcts[level_number].insert(*cit);
    }
    return;
} // inactivateLagrangianStructures

void LDataManager::zeroInactivatedComponents(boost::shared_ptr<LData> lag_data, const int level_number) const
{
    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);

    if (LIKELY(d_inactive_strcts[level_number].empty())) return;

    // Construct a list of all of the inactivated Lagrangian indices.
    //
    // NOTE: The vector idxs is NOT a list of the LOCAL inactivated indices.
    // Instead, it is a list of ALL of the inactivated indices.  Thus, idxs will
    // have the same contents for all MPI processes.
    std::vector<int> idxs;
    for (auto cit = d_inactive_strcts[level_number].begin();
         cit != d_inactive_strcts[level_number].end();
         ++cit)
    {
        const int strct_id = *cit;
        const std::pair<int, int>& lag_index_range =
            d_strct_id_to_lag_idx_range_map[level_number].find(strct_id)->second;
        for (int l = lag_index_range.first; l < lag_index_range.second; ++l)
        {
            idxs.push_back(l);
        }
    }

    // There is nothing left to do if there are no inactivated indices.
    if (LIKELY(idxs.empty())) return;

    // Map the Lagrangian indices to PETSc indices.
    mapLagrangianToPETSc(idxs, level_number);

    // Determine the local extents of the global PETSc Vec.
    int ierr;
    Vec lag_data_vec = lag_data->getVec();
    int lo, hi, bs;
    ierr = VecGetOwnershipRange(lag_data_vec, &lo, &hi);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetBlockSize(lag_data_vec, &bs);
    IBTK_CHKERRQ(ierr);

    // Zero-out all local inactivated components.
    std::vector<int> ix(bs);
    std::vector<double> y(bs, 0.0);
    std::sort(idxs.begin(), idxs.end());
    for (auto cit = idxs.begin(); cit != idxs.end(); ++cit)
    {
        const int l = *cit;
        if (l * bs >= lo && l * bs < hi)
        {
            for (int k = 0; k < bs; ++k)
            {
                ix[k] = bs * l + k;
            }
            ierr = VecSetValues(lag_data_vec, bs, &ix[0], &y[0], INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }
    ierr = VecAssemblyBegin(lag_data_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(lag_data_vec);
    IBTK_CHKERRQ(ierr);
    lag_data->beginGhostUpdate();
    lag_data->endGhostUpdate();
    return;
} // zeroInactivatedComponents

void LDataManager::mapLagrangianToPETSc(std::vector<int>& inds, const int level_number) const
{
    IBTK_TIMER_START(t_map_lagrangian_to_petsc);

    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);
    const int ierr =
        AOApplicationToPetsc(d_ao[level_number],
                             (!inds.empty() ? static_cast<int>(inds.size()) : static_cast<int>(s_ao_dummy.size())),
                             (!inds.empty() ? &inds[0] : &s_ao_dummy[0]));
    IBTK_CHKERRQ(ierr);

    IBTK_TIMER_STOP(t_map_lagrangian_to_petsc);
    return;
} // mapLagrangianToPETSc

void LDataManager::mapPETScToLagrangian(std::vector<int>& inds, const int level_number) const
{
    IBTK_TIMER_START(t_map_petsc_to_lagrangian);

    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);
    const int ierr =
        AOPetscToApplication(d_ao[level_number],
                             (!inds.empty() ? static_cast<int>(inds.size()) : static_cast<int>(s_ao_dummy.size())),
                             (!inds.empty() ? &inds[0] : &s_ao_dummy[0]));
    IBTK_CHKERRQ(ierr);

    IBTK_TIMER_STOP(t_map_petsc_to_lagrangian);
    return;
} // mapPETScToLagrangian

void LDataManager::scatterLagrangianToPETSc(Vec& lagrangian_vec, Vec& petsc_vec, const int level_number) const
{
    scatterData(petsc_vec, lagrangian_vec, level_number, SCATTER_REVERSE);
    return;
} // scatterLagrangianToPETSc

void LDataManager::scatterPETScToLagrangian(Vec& petsc_vec, Vec& lagrangian_vec, const int level_number) const
{
    scatterData(lagrangian_vec, petsc_vec, level_number, SCATTER_FORWARD);
    return;
} // scatterPETScToLagrangian

void LDataManager::scatterToAll(Vec& parallel_vec, Vec& sequential_vec) const
{
    int ierr;
    const bool create_vout = !sequential_vec;
    VecScatter ctx;
    ierr = VecScatterCreateToAll(parallel_vec, &ctx, (create_vout ? &sequential_vec : NULL));
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, parallel_vec, sequential_vec, INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, parallel_vec, sequential_vec, INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx);
    IBTK_CHKERRQ(ierr);
    return;
} // scatterToAll

void LDataManager::scatterToZero(Vec& parallel_vec, Vec& sequential_vec) const
{
    int ierr;
    const bool create_vout = !sequential_vec;
    VecScatter ctx;
    ierr = VecScatterCreateToZero(parallel_vec, &ctx, (create_vout ? &sequential_vec : NULL));
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, parallel_vec, sequential_vec, INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, parallel_vec, sequential_vec, INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx);
    IBTK_CHKERRQ(ierr);
    return;
} // scatterToZero

void LDataManager::beginDataRedistribution(const int coarsest_ln_in, const int finest_ln_in)
{
    IBTK_TIMER_START(t_begin_data_redistribution);

    const int coarsest_ln = (coarsest_ln_in == -1) ? d_coarsest_ln : coarsest_ln_in;
    const int finest_ln = (finest_ln_in == -1) ? d_finest_ln : finest_ln_in;

    TBOX_ASSERT(coarsest_ln >= d_coarsest_ln && coarsest_ln <= d_finest_ln);
    TBOX_ASSERT(finest_ln >= d_coarsest_ln && finest_ln <= d_finest_ln);

    // Emit warnings if things seem to be out of synch.
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (!d_level_contains_lag_data[level_number]) continue;
        if (d_needs_synch[level_number])
        {
            TBOX_WARNING("LDataManager::beginDataRedistribution():\n"
                         << "\tLData is not synchronized with LNodeSetData.\n"
                         << "\tLagrangian node position data is probably invalid!\n");
        }
    }

    // Ensure that no IB points manage to escape the computational domain.
    const double* const domain_x_lower = d_grid_geom->getXLower();
    const double* const domain_x_upper = d_grid_geom->getXUpper();
    const double* const domain_dx = d_grid_geom->getDx();
    std::vector<std::map<int, IntVector> > periodic_offset_data(finest_ln + 1);
    std::vector<std::map<int, Vector> > periodic_displacement_data(finest_ln + 1);
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (!d_level_contains_lag_data[level_number]) continue;
        d_lag_mesh_data[level_number][POSN_DATA_NAME]->beginGhostUpdate();
    }
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (!d_level_contains_lag_data[level_number]) continue;
        d_lag_mesh_data[level_number][POSN_DATA_NAME]->endGhostUpdate();
        // NOTE: We cannot use LMesh data structures here because they have not
        // been (re-)initialized yet.
        boost::multi_array_ref<double, 2>& X_data =
            *d_lag_mesh_data[level_number][POSN_DATA_NAME]->getGhostedLocalFormVecArray();
        const size_t num_nodes = X_data.shape()[0];
        const IntVector& ratio = d_hierarchy->getPatchLevel(level_number)->getRatioToLevelZero();
        const IntVector& periodic_shift = d_grid_geom->getPeriodicShift(ratio);
        double level_dx[NDIM];
        for (unsigned int d = 0; d < NDIM; ++d) level_dx[d] = domain_dx[d] / ratio[d];
        for (unsigned int local_idx = 0; local_idx < num_nodes; ++local_idx)
        {
            Eigen::Map<Vector> X(&X_data[local_idx][0], NDIM);
            Vector X_real = X;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (periodic_shift[d])
                {
                    double domain_length = domain_x_upper[d] - domain_x_lower[d];
                    while (X[d] < domain_x_lower[d]) X[d] += domain_length;
                    while (X[d] >= domain_x_upper[d]) X[d] -= domain_length;
                    TBOX_ASSERT(X[d] >= domain_x_lower[d] && X[d] < domain_x_upper[d]);
                }
                X[d] = std::max(X[d], domain_x_lower[d]);
                X[d] = std::min(X[d], domain_x_upper[d] - std::numeric_limits<double>::epsilon());
            }
            Vector periodic_displacement = X_real - X;
            IntVector periodic_offset(DIM);
            for (int d = 0; d < NDIM; ++d)
            {
                periodic_offset[d] = round(periodic_displacement[d] / level_dx[d]);
            }
            if (periodic_offset != IntVector::getZero(DIM))
            {
                periodic_offset_data[level_number].insert(std::make_pair(local_idx, periodic_offset));
                periodic_displacement_data[level_number][local_idx] = periodic_displacement;
            }
        }
        d_lag_mesh_data[level_number][POSN_DATA_NAME]->restoreArrays();
    }

    // Update the index patch data.
    //
    // We only keep nodes whose new locations are in the patch interior.  That
    // is to say, we only keep the nodes that the patch will own after
    // redistribution.
    //
    // Notice that it is possible for periodic copies of nodes to appear in
    // multiple grid cells within the ghost cell region.  We must therefore
    // ensure that nodes passing through periodic boundaries are added to the
    // patch only once.
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (!d_level_contains_lag_data[level_number]) continue;
        boost::multi_array_ref<double, 2>& X_data =
            *d_lag_mesh_data[level_number][POSN_DATA_NAME]->getGhostedLocalFormVecArray();
        auto level =d_hierarchy->getPatchLevel(level_number);
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch =*p;
            boost::shared_ptr<LNodeSetData> current_idx_data = patch->getPatchData(d_lag_node_index_current_idx);
            auto new_idx_data = boost::make_shared<LNodeSetData>(current_idx_data->getBox(), current_idx_data->getGhostCellWidth());
            const Box& patch_box = patch->getBox();
            auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch->getPatchGeometry());
            const Index& patch_lower = patch_box.lower();
            const Index& patch_upper = patch_box.upper();
            const double* const patch_x_lower = pgeom->getXLower();
            const double* const patch_x_upper = pgeom->getXUpper();
            const double* const patch_dx = pgeom->getDx();
            std::set<int> registered_periodic_idx;
            for (LNodeSetData::CellIterator it(Box::grow(patch_box, IntVector(DIM, CFL_WIDTH))); it; it++)
            {
                const Index& old_cell_idx = *it;
                LNodeSet* const old_node_set = current_idx_data->getItem(old_cell_idx);
                if (old_node_set)
                {
                    for (auto n = old_node_set->begin(); n != old_node_set->end(); ++n)
                    {
                        LNodeSet::value_type& node_idx = *n;
                        const int local_idx = node_idx->getLocalPETScIndex();
                        double* const X = &X_data[local_idx][0];
                        const CellIndex new_cell_idx(IndexUtilities::getCellIndex(
                            X, patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper));
                        if (patch_box.contains(new_cell_idx))
                        {
                            std::map<int, IntVector>::const_iterator it_offset =
                                periodic_offset_data[level_number].find(local_idx);
                            const bool periodic_node = it_offset != periodic_offset_data[level_number].end();
                            const bool unregistered_periodic_node =
                                periodic_node &&
                                registered_periodic_idx.find(local_idx) == registered_periodic_idx.end();
                            if (!periodic_node || unregistered_periodic_node)
                            {
                                if (unregistered_periodic_node)
                                {
                                    const IntVector& periodic_offset = it_offset->second;
                                    std::map<int, Vector>::const_iterator it_displacement =
                                        periodic_displacement_data[level_number].find(local_idx);
                                    const Vector& periodic_displacement = it_displacement->second;
                                    node_idx->registerPeriodicShift(periodic_offset, periodic_displacement);
                                    registered_periodic_idx.insert(local_idx);
                                }
                                if (!new_idx_data->isElement(new_cell_idx))
                                    new_idx_data->appendItemPointer(new_cell_idx, new LNodeSet());
                                LNodeSet* const new_node_set = new_idx_data->getItem(new_cell_idx);
                                new_node_set->push_back(node_idx);
                            }
                        }
                    }
                }
            }
            for (LNodeSetData::SetIterator it(*new_idx_data); it; it++)
            {
                LNodeSet::DataSet& node_set = (*it).getDataSet();
                std::sort(node_set.begin(), node_set.end(), LNodeIndexLagrangianIndexComp());
                node_set.erase(std::unique(node_set.begin(), node_set.end(), LNodeIndexLagrangianIndexEqual()),
                               node_set.end());
            }
            patch->setPatchData(d_lag_node_index_current_idx, new_idx_data);
        }
        d_lag_mesh_data[level_number][POSN_DATA_NAME]->restoreArrays();
        d_needs_synch[level_number] = true;
    }

    IBTK_TIMER_STOP(t_begin_data_redistribution);
    return;
} // beginDataRedistribution

void LDataManager::endDataRedistribution(const int coarsest_ln_in, const int finest_ln_in)
{
    IBTK_TIMER_START(t_end_data_redistribution);

    const int coarsest_ln = (coarsest_ln_in == -1) ? d_coarsest_ln : coarsest_ln_in;
    const int finest_ln = (finest_ln_in == -1) ? d_finest_ln : finest_ln_in;

    TBOX_ASSERT(coarsest_ln >= d_coarsest_ln && coarsest_ln <= d_finest_ln);
    TBOX_ASSERT(finest_ln >= d_coarsest_ln && finest_ln <= d_finest_ln);

    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (d_level_contains_lag_data[level_number] && (!d_needs_synch[level_number]))
        {
            TBOX_WARNING("LDataManager::endDataRedistribution():\n"
                         << "\tLData is already synchronized with LNodeSetData.\n"
                         << "\tlevel = " << level_number << "\n");
        }
    }

    // Fill the ghost cells of each level.
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (!d_level_contains_lag_data[level_number]) continue;

        auto level =d_hierarchy->getPatchLevel(level_number);
        level->allocatePatchData(d_scratch_data);

        const double current_time = 0.0;
        level->setTime(current_time, d_current_data);
        level->setTime(current_time, d_scratch_data);

        d_lag_node_index_bdry_fill_scheds[level_number]->fillData(current_time);

        level->deallocatePatchData(d_scratch_data);
    }

    // Define the PETSc data needed to communicate the LData from its
    // old configuration to its new configuration.
    int ierr;

    std::vector<AO> new_ao(finest_ln + 1);

    std::vector<std::vector<Vec> > src_vec(finest_ln + 1);
    std::vector<std::vector<Vec> > dst_vec(finest_ln + 1);
    std::vector<std::vector<VecScatter> > scatter(finest_ln + 1);
    std::vector<std::map<int, IS> > src_IS(finest_ln + 1);
    std::vector<std::map<int, IS> > dst_IS(finest_ln + 1);
    std::vector<std::map<int, VecScatter> > scatter_template(finest_ln + 1);

    // The number of all local (e.g., on processor) and ghost (e.g., off
    // processor) nodes.
    //
    // NOTE:  num_local_nodes   [ln] == d_local_lag_indices   [ln].size()
    //        num_nonlocal_nodes[ln] == d_nonlocal_lag_indices[ln].size()
    std::vector<int> num_local_nodes(finest_ln + 1);
    std::vector<int> num_nonlocal_nodes(finest_ln + 1);

    // Setup maps from patch numbers to the nodes indexed in the patch interior
    // and the patch ghost cell region.
    //
    // NOTE 1: The ghost cell region used for each patch is defined by the ghost
    // cell width of the indexing variable.
    //
    // NOTE 2: These indices are in the local PETSc ordering.  Hence they can be
    // used to access elements in the local form of ghosted parallel PETSc Vec
    // objects.
    //
    // NOTE 3: The PETSc ordering is maintained so that the data corresponding
    // to patch interiors is contiguous (as long as there are no overlapping
    // patches).  Nodes in the ghost region of a patch will not in general be
    // stored as contiguous data, and no attempt is made to do so.

    // In the following loop over patch levels, we first compute the new
    // distribution data (e.g., all of these indices).
    //
    // Next, we use the old and new PETSc AO (application ordering) objects to
    // define a mapping from the old distribution to the new distribution.
    //
    // Finally, we create the new PETSc Vec (vector) objects that are used to
    // store the Lagrangian data in the new distribution.
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (!d_level_contains_lag_data[level_number]) continue;

        std::map<std::string, boost::shared_ptr<LData> >& level_data = d_lag_mesh_data[level_number];
        const std::vector<int>::size_type num_data = level_data.size();
        src_vec[level_number].resize(num_data);
        dst_vec[level_number].resize(num_data);
        scatter[level_number].resize(num_data);

        // Get the new distribution of nodes for the level.
        //
        // NOTE: This process updates the local PETSc indices of the LNodeSet
        // objects contained in the current patch.
        computeNodeDistribution(new_ao[level_number],
                                d_local_lag_indices[level_number],
                                d_nonlocal_lag_indices[level_number],
                                d_local_petsc_indices[level_number],
                                d_nonlocal_petsc_indices[level_number],
                                d_num_nodes[level_number],
                                d_node_offset[level_number],
                                level_number);
        num_local_nodes[level_number] = static_cast<int>(d_local_lag_indices[level_number].size());
        num_nonlocal_nodes[level_number] = static_cast<int>(d_nonlocal_lag_indices[level_number].size());

        // Setup src indices.
        std::vector<int> src_inds(num_local_nodes[level_number]);
        for (int k = 0; k < num_local_nodes[level_number]; ++k)
        {
            src_inds[k] = d_node_offset[level_number] + k;
        }

        // Convert dst indices from the old ordering to the new ordering.
        std::vector<int> dst_inds = d_local_petsc_indices[level_number];
        ierr = AOPetscToApplication(
            d_ao[level_number],
            static_cast<int>(num_local_nodes[level_number] > 0 ? num_local_nodes[level_number] : s_ao_dummy.size()),
            (num_local_nodes[level_number] > 0 ? &dst_inds[0] : &s_ao_dummy[0]));
        IBTK_CHKERRQ(ierr);
        ierr = AOApplicationToPetsc(
            new_ao[level_number],
            static_cast<int>(num_local_nodes[level_number] > 0 ? num_local_nodes[level_number] : s_ao_dummy.size()),
            (num_local_nodes[level_number] > 0 ? &dst_inds[0] : &s_ao_dummy[0]));
        IBTK_CHKERRQ(ierr);

        // Setup VecScatter objects for each LData object and start scattering
        // data.
        std::map<std::string, boost::shared_ptr<LData> >::iterator it;
        int i;
        for (it = level_data.begin(), i = 0; it != level_data.end(); ++it, ++i)
        {
            boost::shared_ptr<LData> data = it->second;
            TBOX_ASSERT(data);
            const int depth = data->getDepth();

            // Determine the PETSc indices of the source nodes for use when
            // scattering values from the old configuration to the new
            // configuration.  Notice that a different IS object must be used
            // for each unique data depth.
            if (src_IS[level_number].find(depth) == src_IS[level_number].end())
            {
                ierr = ISCreateBlock(PETSC_COMM_WORLD,
                                     depth,
                                     num_local_nodes[level_number],
                                     num_local_nodes[level_number] > 0 ? &src_inds[0] : NULL,
                                     PETSC_COPY_VALUES,
                                     &src_IS[level_number][depth]);
                IBTK_CHKERRQ(ierr);
            }

            // Determine the PETSc indices of the destination nodes for use when
            // scattering values from the old configuration to the new
            // configuration.  Notice that a different IS object must be used
            // for each unique data depth.
            if (dst_IS[level_number].find(depth) == dst_IS[level_number].end())
            {
                ierr = ISCreateBlock(PETSC_COMM_WORLD,
                                     depth,
                                     num_local_nodes[level_number],
                                     num_local_nodes[level_number] > 0 ? &dst_inds[0] : NULL,
                                     PETSC_COPY_VALUES,
                                     &dst_IS[level_number][depth]);
                IBTK_CHKERRQ(ierr);
            }

            // Create the destination Vec.
            src_vec[level_number][i] = data->getVec();
            ierr = VecCreateGhostBlock(
                PETSC_COMM_WORLD,
                depth,
                depth * num_local_nodes[level_number],
                PETSC_DECIDE,
                num_nonlocal_nodes[level_number],
                num_nonlocal_nodes[level_number] > 0 ? &d_nonlocal_petsc_indices[level_number][0] : NULL,
                &dst_vec[level_number][i]);
            IBTK_CHKERRQ(ierr);

            // Create the VecScatter.
            ierr = VecScatterCreate(src_vec[level_number][i],
                                    src_IS[level_number][depth],
                                    dst_vec[level_number][i],
                                    dst_IS[level_number][depth],
                                    &scatter[level_number][i]);
            IBTK_CHKERRQ(ierr);

            // Begin scattering data.
            ierr = VecScatterBegin(scatter[level_number][i],
                                   src_vec[level_number][i],
                                   dst_vec[level_number][i],
                                   INSERT_VALUES,
                                   SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Update cached indexing information on each grid patch and setup new LMesh
    // data structures.
    auto grid_geom = BOOST_CAST<CartesianGridGeometry>(d_hierarchy->getGridGeometry());
    d_local_and_ghost_nodes.resize(finest_ln + 1);
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (!d_level_contains_lag_data[level_number]) continue;

        const int num_local_nodes = getNumberOfLocalNodes(level_number);
        const int num_ghost_nodes = getNumberOfGhostNodes(level_number);
        const int num_local_and_ghost_nodes = num_local_nodes + num_ghost_nodes;
        auto new_local_and_ghost_nodes = boost::make_shared<std::vector<LNode>>(num_local_and_ghost_nodes);
        std::vector<boost::shared_ptr<LNode> > local_and_ghost_node_ptrs(num_local_and_ghost_nodes);
        for (int k = 0; k < num_local_and_ghost_nodes; ++k)
        {
            local_and_ghost_node_ptrs[k] = boost::shared_ptr<LNode>(&(*new_local_and_ghost_nodes)[k], NullDeleter());
        }
        auto level =d_hierarchy->getPatchLevel(level_number);
        const IntVector& periodic_shift = grid_geom->getPeriodicShift(level->getRatioToLevelZero());
        std::set<int> local_petsc_idxs;
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch =*p;
            boost::shared_ptr<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
            idx_data->cacheLocalIndices(patch, periodic_shift);
            const Box& ghost_box = idx_data->getGhostBox();
            for (LNodeSetData::DataIterator it = idx_data->data_begin(ghost_box); it != idx_data->data_end(); ++it)
            {
                LNode* const tmp_node_idx = *it;
                const int local_petsc_idx = tmp_node_idx->getLocalPETScIndex();
                TBOX_ASSERT(0 <= local_petsc_idx && local_petsc_idx < num_local_and_ghost_nodes);
                auto node_idx = local_and_ghost_node_ptrs[local_petsc_idx];
                *node_idx = *tmp_node_idx;
                local_petsc_idxs.insert(local_petsc_idx);
                *it = node_idx;
                TBOX_ASSERT((*it).getPointer() == node_idx.getPointer());
            }
        }
        // Swap out the LNode vectors now --- doing so earlier would potentially
        // reset existing nodes on the patch.
        d_local_and_ghost_nodes[level_number] = new_local_and_ghost_nodes;
        std::vector<LNode*> local_nodes(num_local_nodes);
        for (int k = 0; k < num_local_nodes; ++k)
        {
            local_nodes[k] = &(*d_local_and_ghost_nodes[level_number])[k];
        }
        std::vector<LNode*> ghost_nodes(num_ghost_nodes);
        for (int k = 0; k < num_ghost_nodes; ++k)
        {
            ghost_nodes[k] = &(*d_local_and_ghost_nodes[level_number])[num_local_nodes + k];
        }
        std::ostringstream name_stream;
        name_stream << d_object_name << "::mesh::level_" << level_number;
        d_lag_mesh[level_number] = boost::make_shared<LMesh(name_stream.str>(), local_nodes, ghost_nodes);
    }

    // End scattering data, reset LData objects, and destroy the VecScatter
    // contexts.
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (!d_level_contains_lag_data[level_number]) continue;

        std::map<std::string, boost::shared_ptr<LData> >& level_data = d_lag_mesh_data[level_number];
        std::map<std::string, boost::shared_ptr<LData> >::iterator it;
        int i;
        for (it = level_data.begin(), i = 0; it != level_data.end(); ++it, ++i)
        {
            ierr = VecScatterEnd(scatter[level_number][i],
                                 src_vec[level_number][i],
                                 dst_vec[level_number][i],
                                 INSERT_VALUES,
                                 SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterDestroy(&scatter[level_number][i]);
            IBTK_CHKERRQ(ierr);
            boost::shared_ptr<LData> data = it->second;
            data->resetData(dst_vec[level_number][i], d_nonlocal_petsc_indices[level_number]);
        }
    }

    // Distribute nonlocal data to the new configuration.
    beginNonlocalDataFill(coarsest_ln, finest_ln);
    endNonlocalDataFill(coarsest_ln, finest_ln);

    // Indicate that the levels have been synchronized and destroy unneeded
    // ordering and indexing objects.
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        d_needs_synch[level_number] = false;

        if (d_ao[level_number])
        {
            ierr = AODestroy(&d_ao[level_number]);
            IBTK_CHKERRQ(ierr);
        }
        d_ao[level_number] = new_ao[level_number];

        for (auto it = src_IS[level_number].begin(); it != src_IS[level_number].end(); ++it)
        {
            ierr = ISDestroy(&it->second);
            IBTK_CHKERRQ(ierr);
        }

        for (auto it = dst_IS[level_number].begin(); it != dst_IS[level_number].end(); ++it)
        {
            ierr = ISDestroy(&it->second);
            IBTK_CHKERRQ(ierr);
        }
    }

    // If a Silo data writer is registered with the manager, give it access to
    // the new application orderings.
    if (d_silo_writer)
    {
        d_silo_writer->registerLagrangianAO(d_ao, coarsest_ln, finest_ln);
    }

    IBTK_TIMER_STOP(t_end_data_redistribution);
    return;
} // endDataRedistribution

void LDataManager::updateWorkloadEstimates(const int coarsest_ln_in, const int finest_ln_in)
{
    if (!d_load_balancer) return;

    IBTK_TIMER_START(t_update_workload_estimates);

    const int coarsest_ln = (coarsest_ln_in == -1) ? d_coarsest_ln : coarsest_ln_in;
    const int finest_ln = (finest_ln_in == -1) ? d_finest_ln : finest_ln_in;

    TBOX_ASSERT(coarsest_ln >= d_coarsest_ln && coarsest_ln <= d_finest_ln);
    TBOX_ASSERT(finest_ln >= d_coarsest_ln && finest_ln <= d_finest_ln);

    updateNodeCountData(coarsest_ln, finest_ln);
    HierarchyCellDataOpsReal<double> hier_cc_data_ops(d_hierarchy, coarsest_ln, finest_ln);
    hier_cc_data_ops.axpy(d_workload_idx, d_beta_work, d_node_count_idx, d_workload_idx);

    IBTK_TIMER_STOP(t_update_workload_estimates);
    return;
} // updateWorkloadEstimates

void LDataManager::updateNodeCountData(const int coarsest_ln_in, const int finest_ln_in)
{
    IBTK_TIMER_START(t_update_node_count_data);

    const int coarsest_ln = (coarsest_ln_in == -1) ? d_coarsest_ln : coarsest_ln_in;
    const int finest_ln = (finest_ln_in == -1) ? d_finest_ln : finest_ln_in;

    TBOX_ASSERT(coarsest_ln >= d_coarsest_ln && coarsest_ln <= d_finest_ln);
    TBOX_ASSERT(finest_ln >= d_coarsest_ln && finest_ln <= d_finest_ln);

    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        auto level =d_hierarchy->getPatchLevel(level_number);
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch =*p;
            const Box& patch_box = patch->getBox();
            const boost::shared_ptr<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
            boost::shared_ptr<CellData<double> > node_count_data = patch->getPatchData(d_node_count_idx);
            node_count_data->fillAll(0.0);
            for (LNodeSetData::SetIterator it(*idx_data); it; it++)
            {
                const Index& i = it.getIndex();
                if (patch_box.contains(i))
                {
                    const LNodeSet& node_set = *it;
                    (*node_count_data)(CellIndex(i)) = node_set.size();
                }
            }
        }
    }

    IBTK_TIMER_STOP(t_update_node_count_data);
    return;
} // updateNodeCountData

void LDataManager::initializeLevelData(const boost::shared_ptr<PatchHierarchy> hierarchy,
                                       const int level_number,
                                       const double init_data_time,
                                       const bool can_be_refined,
                                       const bool initial_time,
                                       const boost::shared_ptr<PatchLevel> old_level,
                                       const bool allocate_data)
{
    IBTK_TIMER_START(t_initialize_level_data);

    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    if (old_level)
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));

    auto level =hierarchy->getPatchLevel(level_number);

    // Allocate storage needed to initialize the level and fill data from
    // coarser levels in AMR hierarchy, if any.
    //
    // Since time gets set when we allocate data, re-stamp it to current time if
    // we don't need to allocate.
    if (allocate_data && (level_number >= d_coarsest_ln) && (level_number <= d_finest_ln))
    {
        level->allocatePatchData(d_current_data, init_data_time);
    }
    else if (level_number >= d_coarsest_ln && level_number <= d_finest_ln)
    {
        level->setTime(init_data_time, d_current_data);
    }

    // Fill data from the old level when available.
    if (old_level && d_level_contains_lag_data[level_number])
    {
        level->allocatePatchData(d_scratch_data, init_data_time);
        d_lag_node_index_bdry_fill_alg->createSchedule(level, old_level)->fillData(init_data_time);
        level->deallocatePatchData(d_scratch_data);
    }

    // Initialize the data on the level and, when appropriate, move data from
    // coarser levels to finer levels.
    if (initial_time)
    {
        // Resize some arrays.
        d_level_contains_lag_data.resize(level_number + 1);
        d_strct_name_to_strct_id_map.resize(level_number + 1);
        d_strct_id_to_strct_name_map.resize(level_number + 1);
        d_strct_id_to_lag_idx_range_map.resize(level_number + 1);
        d_last_lag_idx_to_strct_id_map.resize(level_number + 1);
        d_inactive_strcts.resize(level_number + 1);
        d_lag_mesh.resize(level_number + 1);
        d_lag_mesh_data.resize(level_number + 1);
        d_needs_synch.resize(level_number + 1, false);
        d_ao.resize(level_number + 1);
        d_num_nodes.resize(level_number + 1);
        d_node_offset.resize(level_number + 1);
        d_local_lag_indices.resize(level_number + 1);
        d_nonlocal_lag_indices.resize(level_number + 1);
        d_local_petsc_indices.resize(level_number + 1);
        d_nonlocal_petsc_indices.resize(level_number + 1);
        TBOX_ASSERT(d_lag_init);
        d_level_contains_lag_data[level_number] = d_lag_init->getLevelHasLagrangianData(level_number, can_be_refined);
    }
    if (initial_time && d_level_contains_lag_data[level_number])
    {
        int ierr;

        // 1. Determine the number of local (on processor) nodes to be allocated
        //    on the patch level and allocate space for the local and non-local
        //    index data.
        const unsigned int num_global_nodes = d_lag_init->computeGlobalNodeCountOnPatchLevel(
            hierarchy, level_number, init_data_time, can_be_refined, initial_time);
        const unsigned int num_local_nodes = d_lag_init->computeLocalNodeCountOnPatchLevel(
            hierarchy, level_number, init_data_time, can_be_refined, initial_time);
        tbox::SAMRAI_MPI comm(MPI_COMM_WORLD);
        int sum_num_local_nodes = num_local_nodes;
        comm.AllReduce(&sum_num_local_nodes, 1, MPI_SUM);
        if (num_global_nodes != static_cast<unsigned int>(sum_num_local_nodes))
        {
            TBOX_ERROR("LDataManager::initializeLevelData()"
                       << "\n"
                       << "  num_global_nodes    = " << num_global_nodes << "\n"
                       << "  sum num_local_nodes = " << sum_num_local_nodes << "\n");
        }

        d_local_lag_indices[level_number].resize(num_local_nodes, -1);
        d_local_petsc_indices[level_number].resize(num_local_nodes, -1);

        d_nonlocal_lag_indices[level_number].clear();
        d_nonlocal_petsc_indices[level_number].clear();

        computeNodeOffsets(d_num_nodes[level_number], d_node_offset[level_number], num_local_nodes);
        TBOX_ASSERT(d_num_nodes[level_number] == num_global_nodes);

        // 2. Allocate LData corresponding to the curvilinear mesh node
        //    positions and velocities.
        static const bool maintain_data = true;
        createLData(POSN_DATA_NAME, level_number, NDIM, maintain_data);
        createLData(INIT_POSN_DATA_NAME, level_number, NDIM, maintain_data);
        createLData(VEL_DATA_NAME, level_number, NDIM, maintain_data);

        // 3. Initialize the Lagrangian data.
        d_lag_init->initializeStructureIndexingOnPatchLevel(d_strct_id_to_strct_name_map[level_number],
                                                            d_strct_id_to_lag_idx_range_map[level_number],
                                                            level_number,
                                                            init_data_time,
                                                            can_be_refined,
                                                            initial_time,
                                                            this);

        for (auto cit(d_strct_id_to_strct_name_map[level_number].begin());
             cit != d_strct_id_to_strct_name_map[level_number].end();
             ++cit)
        {
            d_strct_name_to_strct_id_map[level_number][cit->second] = cit->first;
        }

        for (auto cit(
                 d_strct_id_to_lag_idx_range_map[level_number].begin());
             cit != d_strct_id_to_lag_idx_range_map[level_number].end();
             ++cit)
        {
            d_last_lag_idx_to_strct_id_map[level_number][cit->second.second - 1] = cit->first;
        }

        // WARNING: If either of the following offsets is ever nonzero, note
        // that it may be necessary to modify IBHierarchyIntegrator, in
        // particular the code where the data related to the implementation of
        // the penalty IB method are initialized.
        static const unsigned int global_index_offset = 0;
        static const unsigned int local_index_offset = 0;
        const unsigned int num_initialized_local_nodes =
            d_lag_init->initializeDataOnPatchLevel(d_lag_node_index_current_idx,
                                                   global_index_offset,
                                                   local_index_offset,
                                                   d_lag_mesh_data[level_number][POSN_DATA_NAME],
                                                   d_lag_mesh_data[level_number][VEL_DATA_NAME],
                                                   hierarchy,
                                                   level_number,
                                                   init_data_time,
                                                   can_be_refined,
                                                   initial_time,
                                                   this);

        ierr = VecCopy(d_lag_mesh_data[level_number][POSN_DATA_NAME]->getVec(),
                       d_lag_mesh_data[level_number][INIT_POSN_DATA_NAME]->getVec());
        IBTK_CHKERRQ(ierr);

        if (num_local_nodes != num_initialized_local_nodes)
        {
            TBOX_ERROR("LDataManager::initializeLevelData()"
                       << "\n"
                       << "  num_local_nodes             = " << num_local_nodes << "\n"
                       << "  num_initialized_local_nodes = " << num_initialized_local_nodes << "\n");
        }

        // 4. Compute the initial distribution (indexing) data.
        auto grid_geom = BOOST_CAST<CartesianGridGeometry>(d_hierarchy->getGridGeometry());
        const IntVector& periodic_shift = grid_geom->getPeriodicShift(level->getRatioToLevelZero());
        std::set<LNode*, LNodeIndexLocalPETScIndexComp> local_nodes, ghost_nodes;
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch =*p;
            const Box& patch_box = patch->getBox();

            boost::shared_ptr<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
            boost::shared_ptr<CellData<double> > node_count_data = patch->getPatchData(d_node_count_idx);

            node_count_data->fillAll(0.0);

            idx_data->cacheLocalIndices(patch, periodic_shift);
            for (LNodeSetData::SetIterator it(*idx_data); it; it++)
            {
                const Index& i = it.getIndex();
                LNodeSet& node_set = *it;
                const bool patch_owns_idx_set = patch_box.contains(i);
                if (patch_owns_idx_set)
                {
                    (*node_count_data)(CellIndex(i)) = node_set.size();
                }

                for (auto n = node_set.begin(); n != node_set.end(); ++n)
                {
                    LNode* const node_idx = *n;
                    const int lag_idx = node_idx->getLagrangianIndex();
                    const int local_idx = node_idx->getLocalPETScIndex();
                    if (!(0 <= local_idx && local_idx < static_cast<int>(num_local_nodes)))
                    {
                        TBOX_ERROR("LDataManager::initializeLevelData()"
                                   << "\n"
                                   << "  local_idx       = " << local_idx << "\n"
                                   << "  num_local_nodes = " << num_local_nodes << "\n");
                    }
                    d_local_lag_indices[level_number][local_idx] = lag_idx;
                    d_local_petsc_indices[level_number][local_idx] = local_idx + d_node_offset[level_number];
                    local_nodes.insert(node_idx);
                }
            }
        }
        int num_initialized_global_nodes = static_cast<int>(local_nodes.size());
        comm.AllReduce(&num_initialized_global_nodes, 1, MPI_SUM);
        if (d_num_nodes[level_number] != static_cast<unsigned int>(num_initialized_global_nodes))
        {
            TBOX_ERROR("LDataManager::initializeLevelData()"
                       << "\n"
                       << "  num_nodes[level_number] = " << d_num_nodes[level_number] << "\n"
                       << "  num_initialized_global_nodes = " << num_initialized_global_nodes << "\n");
        }

        std::ostringstream name_stream;
        name_stream << d_object_name << "::mesh::level_" << level_number;
        d_lag_mesh[level_number] = boost::make_shared<LMesh(name_stream.str>(),
                                             std::vector<LNode*>(local_nodes.begin(), local_nodes.end()),
                                             std::vector<LNode*>(ghost_nodes.begin(), ghost_nodes.end()));

        // 5. The AO (application order) is determined by the initial values of
        //    the local Lagrangian indices.
        if (d_ao[level_number])
        {
            ierr = AODestroy(&d_ao[level_number]);
            IBTK_CHKERRQ(ierr);
        }

        ierr = AOCreateMapping(PETSC_COMM_WORLD,
                               num_local_nodes,
                               num_local_nodes > 0 ? &d_local_lag_indices[level_number][0] : NULL,
                               num_local_nodes > 0 ? &d_local_petsc_indices[level_number][0] : NULL,
                               &d_ao[level_number]);
        IBTK_CHKERRQ(ierr);
    }

    // If a Silo data writer is registered with the manager, give it access to
    // the new application ordering.
    if (d_silo_writer && d_level_contains_lag_data[level_number])
    {
        d_silo_writer->registerCoordsData(d_lag_mesh_data[level_number][POSN_DATA_NAME], level_number);
        d_silo_writer->registerLagrangianAO(d_ao[level_number], level_number);
    }

    IBTK_TIMER_STOP(t_initialize_level_data);
    return;
} // initializeLevelData

void LDataManager::resetHierarchyConfiguration(const boost::shared_ptr<PatchHierarchy> hierarchy,
                                               const int coarsest_ln,
                                               const int finest_ln)
{
    IBTK_TIMER_START(t_reset_hierarchy_configuration);

    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((coarsest_ln >= 0) && (coarsest_ln <= finest_ln) && (finest_ln <= hierarchy->getFinestLevelNumber()));
    for (int level_number = 0; level_number <= finest_ln; ++level_number)
    {
        TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
    }
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // Reset the patch hierarchy and levels.
    setPatchHierarchy(hierarchy);
    setPatchLevels(0, finest_hier_level);

    // Reset the Silo data writer.
    if (d_silo_writer)
    {
        d_silo_writer->setPatchHierarchy(hierarchy);
        d_silo_writer->resetLevels(d_coarsest_ln, d_finest_ln);
        for (int level_number = d_coarsest_ln; level_number <= d_finest_ln; ++level_number)
        {
            if (!d_level_contains_lag_data[level_number]) continue;
            d_silo_writer->registerCoordsData(d_lag_mesh_data[level_number][POSN_DATA_NAME], level_number);
        }
    }

    // If we have added or removed a level, resize the schedule vectors.
    d_lag_node_index_bdry_fill_scheds.resize(finest_hier_level + 1);
    d_node_count_coarsen_scheds.resize(finest_hier_level + 1);

    // (Re)build refine communication schedules.  These are created for only the
    // specified levels in the hierarchy.
    //
    // NOTE: These schedules do not fill from coarser levels in the hierarchy.
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        auto level =hierarchy->getPatchLevel(level_number);
        d_lag_node_index_bdry_fill_scheds[level_number] = d_lag_node_index_bdry_fill_alg->createSchedule(level);
    }

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (int level_number = std::max(coarsest_ln, 1); level_number <= finest_hier_level; ++level_number)
    {
        auto level =hierarchy->getPatchLevel(level_number);
        auto coarser_level = hierarchy->getPatchLevel(level_number - 1);
        d_node_count_coarsen_scheds[level_number] = d_node_count_coarsen_alg->createSchedule(coarser_level, level);
    }

    IBTK_TIMER_STOP(t_reset_hierarchy_configuration);
    return;
} // resetHierarchyConfiguration

void LDataManager::applyGradientDetector(const boost::shared_ptr<PatchHierarchy> hierarchy,
                                         const int level_number,
                                         const double error_data_time,
                                         const int tag_index,
                                         const bool initial_time,
                                         const bool /*uses_richardson_extrapolation_too*/)
{
    IBTK_TIMER_START(t_apply_gradient_detector);

    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));

    if (initial_time)
    {
        // Tag cells for refinement based on the initial configuration of the
        // Lagrangian structure.
        d_lag_init->tagCellsForInitialRefinement(hierarchy, level_number, error_data_time, tag_index);
    }
    else if (hierarchy->finerLevelExists(level_number))
    {
        auto level =hierarchy->getPatchLevel(level_number);
        auto finer_level = hierarchy->getPatchLevel(level_number + 1);

        // Zero out the node count data on the current level.
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch =*p;
            boost::shared_ptr<CellData<double> > node_count_data = patch->getPatchData(d_node_count_idx);
            node_count_data->fillAll(0.0);
        }

        // Compute the node count data on the next finer level of the patch
        // hierarchy.
        updateNodeCountData(level_number + 1, level_number + 1);

        // Coarsen the node count data from the next finer level of the patch
        // hierarchy.
        d_node_count_coarsen_scheds[level_number + 1]->coarsenData();

        // Tag cells for refinement wherever there exist nodes on the next finer
        // level of the Cartesian grid.
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            const auto patch =*p;
            const Box& patch_box = patch->getBox();

            boost::shared_ptr<CellData<int> > tag_data = patch->getPatchData(tag_index);
            const boost::shared_ptr<CellData<double> > node_count_data = patch->getPatchData(d_node_count_idx);

            for (CellIterator ic(patch_box); ic; ic++)
            {
                const CellIndex& i = ic();
                if (!MathUtilities<double>::equalEps((*node_count_data)(i), 0.0))
                {
                    (*tag_data)(i) = 1;
                }
            }
        }

        // Re-compute the node count data on the present level of the patch
        // hierarchy (since it was invalidated above).
        updateNodeCountData(level_number, level_number);
    }

    IBTK_TIMER_STOP(t_apply_gradient_detector);
    return;
} // applyGradientDetector

void LDataManager::putToRestart(const boost::shared_ptr<Database>& db) const
{
    IBTK_TIMER_START(t_put_to_database);

    TBOX_ASSERT(db);

    db->putInteger("LDATA_MANAGER_VERSION", LDATA_MANAGER_VERSION);
    db->putInteger("d_coarsest_ln", d_coarsest_ln);
    db->putInteger("d_finest_ln", d_finest_ln);
    db->putDouble("d_beta_work", d_beta_work);

    // Write out data that is stored on a level-by-level basis.
    for (int level_number = d_coarsest_ln; level_number <= d_finest_ln; ++level_number)
    {
        std::ostringstream stream;
        stream << "level_" << level_number;
        const std::string level_db_name = stream.str();
        auto level_db = db->putDatabase(level_db_name);

        level_db->putBool("d_level_contains_lag_data", d_level_contains_lag_data[level_number]);

        if (!d_level_contains_lag_data[level_number]) continue;

        std::vector<int> lstruct_ids, lstruct_lag_idx_range_first, lstruct_lag_idx_range_second, lstruct_activation;
        std::vector<std::string> lstruct_names;
        for (auto it = d_strct_id_to_strct_name_map[level_number].begin();
             it != d_strct_id_to_strct_name_map[level_number].end();
             ++it)
        {
            const int id = it->first;
            lstruct_ids.push_back(id);
            lstruct_names.push_back(it->second);
            lstruct_lag_idx_range_first.push_back(d_strct_id_to_lag_idx_range_map[level_number].find(id)->second.first);
            lstruct_lag_idx_range_second.push_back(
                d_strct_id_to_lag_idx_range_map[level_number].find(id)->second.second);
            lstruct_activation.push_back(
                static_cast<int>(d_inactive_strcts[level_number].find(id) != d_inactive_strcts[level_number].end()));
        }
        level_db->putInteger("n_lstructs", static_cast<int>(lstruct_ids.size()));
        if (!lstruct_ids.empty())
        {
            level_db->putIntegerArray("lstruct_ids", &lstruct_ids[0], static_cast<int>(lstruct_ids.size()));
            level_db->putIntegerArray("lstruct_lag_idx_range_first",
                                      &lstruct_lag_idx_range_first[0],
                                      static_cast<int>(lstruct_lag_idx_range_first.size()));
            level_db->putIntegerArray("lstruct_lag_idx_range_second",
                                      &lstruct_lag_idx_range_second[0],
                                      static_cast<int>(lstruct_lag_idx_range_second.size()));
            level_db->putIntegerArray(
                "lstruct_activation", &lstruct_activation[0], static_cast<int>(lstruct_activation.size()));
            level_db->putStringArray("lstruct_names", &lstruct_names[0], static_cast<int>(lstruct_names.size()));
        }

        std::vector<std::string> ldata_names;
        for (auto it = d_lag_mesh_data[level_number].begin();
             it != d_lag_mesh_data[level_number].end();
             ++it)
        {
            ldata_names.push_back(it->first);
            it->second->putToRestart(level_db->putDatabase(ldata_names.back()));
        }
        level_db->putInteger("n_ldata_names", static_cast<int>(ldata_names.size()));
        if (!ldata_names.empty())
        {
            level_db->putStringArray("ldata_names", &ldata_names[0], static_cast<int>(ldata_names.size()));
        }

        level_db->putInteger("d_num_nodes", d_num_nodes[level_number]);
        level_db->putInteger("d_node_offset", d_node_offset[level_number]);

        level_db->putInteger("n_local_lag_indices", static_cast<int>(d_local_lag_indices[level_number].size()));
        if (!d_local_lag_indices[level_number].empty())
        {
            level_db->putIntegerArray("d_local_lag_indices",
                                      &d_local_lag_indices[level_number][0],
                                      static_cast<int>(d_local_lag_indices[level_number].size()));
        }
        level_db->putInteger("n_nonlocal_lag_indices", static_cast<int>(d_nonlocal_lag_indices[level_number].size()));
        if (!d_nonlocal_lag_indices[level_number].empty())
        {
            level_db->putIntegerArray("d_nonlocal_lag_indices",
                                      &d_nonlocal_lag_indices[level_number][0],
                                      static_cast<int>(d_nonlocal_lag_indices[level_number].size()));
        }
        level_db->putInteger("n_local_petsc_indices", static_cast<int>(d_local_petsc_indices[level_number].size()));
        if (!d_local_petsc_indices[level_number].empty())
        {
            level_db->putIntegerArray("d_local_petsc_indices",
                                      &d_local_petsc_indices[level_number][0],
                                      static_cast<int>(d_local_petsc_indices[level_number].size()));
        }
        // NOTE: d_nonlocal_petsc_indices[level_number] is a map from the data
        // depth to the nonlocal petsc indices for that particular depth.  We
        // only serialize the indices corresponding to a data depth of 1.
        level_db->putInteger("n_nonlocal_petsc_indices",
                             static_cast<int>(d_nonlocal_petsc_indices[level_number].size()));
        if (!d_nonlocal_petsc_indices[level_number].empty())
        {
            level_db->putIntegerArray("d_nonlocal_petsc_indices",
                                      &d_nonlocal_petsc_indices[level_number][0],
                                      static_cast<int>(d_nonlocal_petsc_indices[level_number].size()));
        }
    }

    IBTK_TIMER_STOP(t_put_to_database);
    return;
} // putToRestart

/////////////////////////////// PROTECTED ////////////////////////////////////

LDataManager::LDataManager(const std::string& object_name,
                           const std::string& default_interp_kernel_fcn,
                           const std::string& default_spread_kernel_fcn,
                           const IntVector& ghost_width,
                           bool register_for_restart)
    : d_object_name(object_name), d_registered_for_restart(register_for_restart), d_hierarchy(NULL), d_grid_geom(NULL),
      d_coarsest_ln(-1), d_finest_ln(-1), d_visit_writer(NULL), d_silo_writer(NULL), d_load_balancer(NULL),
      d_lag_init(NULL), d_level_contains_lag_data(), d_lag_node_index_var(NULL), d_lag_node_index_current_idx(-1),
      d_lag_node_index_scratch_idx(-1), d_beta_work(1.0), d_workload_var(NULL), d_workload_idx(-1),
      d_output_workload(false), d_node_count_var(NULL), d_node_count_idx(-1), d_output_node_count(false),
      d_default_interp_kernel_fcn(default_interp_kernel_fcn), d_default_spread_kernel_fcn(default_spread_kernel_fcn),
      d_ghost_width(ghost_width), d_lag_node_index_bdry_fill_alg(NULL), d_lag_node_index_bdry_fill_scheds(),
      d_node_count_coarsen_alg(NULL), d_node_count_coarsen_scheds(), d_current_context(NULL), d_scratch_context(NULL),
      d_current_data(), d_scratch_data(), d_lag_mesh(), d_lag_mesh_data(), d_needs_synch(true), d_ao(), d_num_nodes(),
      d_node_offset(), d_local_lag_indices(), d_nonlocal_lag_indices(), d_local_petsc_indices(),
      d_nonlocal_petsc_indices()
{
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(ghost_width.min() >= 0);

    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }

    // Create/look up the variable contexts.
    VariableDatabase* var_db = VariableDatabase::getDatabase();

    d_current_context = var_db->getContext(d_object_name + "::CURRENT");
    d_scratch_context = var_db->getContext(d_object_name + "::SCRATCH");

    // Register the SAMRAI variables with the VariableDatabase.
    d_lag_node_index_var = boost::make_shared<LNodeSetVariable>(d_object_name + "::lag_node_index");

    // Setup the current context.
    d_lag_node_index_current_idx =
        var_db->registerVariableAndContext(d_lag_node_index_var, d_current_context, d_ghost_width);
    d_current_data.setFlag(d_lag_node_index_current_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_lag_node_index_current_idx);
    }

    // Setup the scratch context.
    d_lag_node_index_scratch_idx =
        var_db->registerVariableAndContext(d_lag_node_index_var, d_scratch_context, d_ghost_width);
    d_scratch_data.setFlag(d_lag_node_index_scratch_idx);

    // Setup a refine algorithm, used to fill LNode boundary data.
    boost::shared_ptr<RefineOperator> lag_node_index_bdry_fill_op;
    d_lag_node_index_bdry_fill_alg = boost::make_shared<RefineAlgorithm>();
    d_lag_node_index_bdry_fill_alg->registerRefine(d_lag_node_index_current_idx,
                                                   d_lag_node_index_current_idx,
                                                   d_lag_node_index_scratch_idx,
                                                   lag_node_index_bdry_fill_op);

    // Register the node count variable with the VariableDatabase.
    d_node_count_var = boost::make_shared<CellVariable<double> >(DIM, d_object_name + "::node_count");
    d_node_count_idx = var_db->registerVariableAndContext(d_node_count_var, d_current_context, IntVector::getZero(DIM));
    d_current_data.setFlag(d_node_count_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_node_count_idx);
    }

    auto node_count_coarsen_op = boost::make_shared<CartesianCellDoubleWeightedAverage>(DIM);
    d_node_count_coarsen_alg = boost::make_shared<CoarsenAlgorithm>(DIM);
    d_node_count_coarsen_alg->registerCoarsen(d_node_count_idx, d_node_count_idx, node_count_coarsen_op);

    // Setup Timers.
    IBTK_DO_ONCE(
        t_spread = TimerManager::getManager()->getTimer("IBTK::LDataManager::spread()");
        t_interp = TimerManager::getManager()->getTimer("IBTK::LDataManager::interp()");
        t_map_lagrangian_to_petsc = TimerManager::getManager()->getTimer("IBTK::LDataManager::mapLagrangianToPETSc()");
        t_map_petsc_to_lagrangian = TimerManager::getManager()->getTimer("IBTK::LDataManager::mapPETScToLagrangian()");
        t_begin_data_redistribution =
            TimerManager::getManager()->getTimer("IBTK::LDataManager::beginDataRedistribution()");
        t_end_data_redistribution = TimerManager::getManager()->getTimer("IBTK::LDataManager::endDataRedistribution()");
        t_update_workload_estimates =
            TimerManager::getManager()->getTimer("IBTK::LDataManager::updateWorkloadEstimates()");
        t_update_node_count_data = TimerManager::getManager()->getTimer("IBTK::LDataManager::updateNodeCountData()");
        t_initialize_level_data = TimerManager::getManager()->getTimer("IBTK::LDataManager::initializeLevelData()");
        t_reset_hierarchy_configuration =
            TimerManager::getManager()->getTimer("IBTK::LDataManager::resetHierarchyConfiguration()");
        t_apply_gradient_detector = TimerManager::getManager()->getTimer("IBTK::LDataManager::applyGradientDetector()");
        t_put_to_database = TimerManager::getManager()->getTimer("IBTK::LDataManager::putToRestart()");
        t_begin_nonlocal_data_fill =
            TimerManager::getManager()->getTimer("IBTK::LDataManager::beginNonlocalDataFill()");
        t_end_nonlocal_data_fill = TimerManager::getManager()->getTimer("IBTK::LDataManager::endNonlocalDataFill()");
        t_compute_node_distribution =
            TimerManager::getManager()->getTimer("IBTK::LDataManager::computeNodeDistribution()");
        t_compute_node_offsets = TimerManager::getManager()->getTimer("IBTK::LDataManager::computeNodeOffsets()"););
    return;
} // LDataManager

LDataManager::~LDataManager()
{
    // Destroy any remaining AO objects.
    int ierr;
    for (int level_number = d_coarsest_ln; level_number <= d_finest_ln; ++level_number)
    {
        if (d_ao[level_number])
        {
            ierr = AODestroy(&d_ao[level_number]);
            IBTK_CHKERRQ(ierr);
        }
    }
    return;
} // ~LDataManager

/////////////////////////////// PRIVATE //////////////////////////////////////

void LDataManager::scatterData(Vec& lagrangian_vec, Vec& petsc_vec, const int level_number, ScatterMode mode) const
{
    int ierr;

    // Get the vector sizes.
    int petsc_size, lagrangian_size;
    ierr = VecGetSize(petsc_vec, &petsc_size);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetSize(lagrangian_vec, &lagrangian_size);
    IBTK_CHKERRQ(ierr);
    TBOX_ASSERT(petsc_size == lagrangian_size);
    int petsc_bs, lagrangian_bs;
    ierr = VecGetBlockSize(petsc_vec, &petsc_bs);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetBlockSize(lagrangian_vec, &lagrangian_bs);
    IBTK_CHKERRQ(ierr);
    TBOX_ASSERT(petsc_bs == lagrangian_bs);
    const int depth = petsc_bs;

    // Determine the application indices corresponding to the local PETSc
    // indices.
    int local_sz;
    ierr = VecGetLocalSize(lagrangian_vec, &local_sz);
    IBTK_CHKERRQ(ierr);
    local_sz /= depth;
    std::vector<int> local_lag_idxs(local_sz, -1);

    int ilo, ihi;
    ierr = VecGetOwnershipRange(lagrangian_vec, &ilo, &ihi);
    IBTK_CHKERRQ(ierr);
    ilo /= depth;
    ihi /= depth;
    for (int k = 0; k < local_sz; ++k)
    {
        local_lag_idxs[k] = ilo + k;
    }
    mapLagrangianToPETSc(local_lag_idxs, level_number);

    IS lag_is;
    ierr = ISCreateBlock(PETSC_COMM_WORLD,
                         depth,
                         static_cast<int>(local_lag_idxs.size()),
                         local_lag_idxs.empty() ? NULL : &local_lag_idxs[0],
                         PETSC_COPY_VALUES,
                         &lag_is);
    IBTK_CHKERRQ(ierr);

    // Create a VecScatter to scatter data from the distributed PETSc
    // representation to the distributed Lagrangian representation.
    VecScatter vec_scatter;
    ierr = VecScatterCreate(petsc_vec, lag_is, lagrangian_vec, NULL, &vec_scatter);
    IBTK_CHKERRQ(ierr);

    // Scatter the values.
    ierr = VecScatterBegin(vec_scatter, petsc_vec, lagrangian_vec, INSERT_VALUES, mode);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterEnd(vec_scatter, petsc_vec, lagrangian_vec, INSERT_VALUES, mode);
    IBTK_CHKERRQ(ierr);

    // Cleanup allocated data.
    ierr = ISDestroy(&lag_is);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterDestroy(&vec_scatter);
    IBTK_CHKERRQ(ierr);
    return;
} // scatterData

void LDataManager::beginNonlocalDataFill(const int coarsest_ln_in, const int finest_ln_in)
{
    IBTK_TIMER_START(t_begin_nonlocal_data_fill);

    const int coarsest_ln = (coarsest_ln_in == -1) ? d_coarsest_ln : coarsest_ln_in;
    const int finest_ln = (finest_ln_in == -1) ? d_finest_ln : finest_ln_in;

    TBOX_ASSERT(coarsest_ln >= d_coarsest_ln && coarsest_ln <= d_finest_ln);
    TBOX_ASSERT(finest_ln >= d_coarsest_ln && finest_ln <= d_finest_ln);

    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        std::map<std::string, boost::shared_ptr<LData> >& level_data = d_lag_mesh_data[level_number];
        for (auto it = level_data.begin(); it != level_data.end(); ++it)
        {
            it->second->beginGhostUpdate();
        }
    }

    IBTK_TIMER_STOP(t_begin_nonlocal_data_fill);
    return;
} // beginNonlocalDataFill

void LDataManager::endNonlocalDataFill(const int coarsest_ln_in, const int finest_ln_in)
{
    IBTK_TIMER_START(t_end_nonlocal_data_fill);

    const int coarsest_ln = (coarsest_ln_in == -1) ? d_coarsest_ln : coarsest_ln_in;
    const int finest_ln = (finest_ln_in == -1) ? d_finest_ln : finest_ln_in;

    TBOX_ASSERT(coarsest_ln >= d_coarsest_ln && coarsest_ln <= d_finest_ln);
    TBOX_ASSERT(finest_ln >= d_coarsest_ln && finest_ln <= d_finest_ln);

    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        std::map<std::string, boost::shared_ptr<LData> >& level_data = d_lag_mesh_data[level_number];
        for (auto it = level_data.begin(); it != level_data.end(); ++it)
        {
            it->second->endGhostUpdate();
        }
    }

    IBTK_TIMER_STOP(t_end_nonlocal_data_fill);
    return;
} // endNonlocalDataFill

void LDataManager::computeNodeDistribution(AO& ao,
                                           std::vector<int>& local_lag_indices,
                                           std::vector<int>& nonlocal_lag_indices,
                                           std::vector<int>& local_petsc_indices,
                                           std::vector<int>& nonlocal_petsc_indices,
                                           unsigned int& num_nodes,
                                           unsigned int& node_offset,
                                           const int level_number)
{
    IBTK_TIMER_START(t_compute_node_distribution);

    TBOX_ASSERT(level_number >= d_coarsest_ln && level_number <= d_finest_ln);

    local_lag_indices.clear();
    nonlocal_lag_indices.clear();
    local_petsc_indices.clear();
    nonlocal_petsc_indices.clear();

    // Determine the Lagrangian IDs of all of the Lagrangian nodes on the
    // specified level of the patch hierarchy.
    //
    // We differentiate between nodes that are local to the processor
    // (i.e. nodes that live in the interior of some patch owned by the
    // processor) and nodes that are non-local (i.e. nodes that live in the
    // interior of a patch owned by a different processor).
    //
    // It is important to emphasize that while a local node by definition lives
    // on the interior of some patch on this processor, it may also live in the
    // ghost cell regions of other patches owned by this processor.
    //
    // Non-local nodes ONLY appear in ghost cells for on processor patches.
    auto level =d_hierarchy->getPatchLevel(level_number);

    // Collect the local nodes and assign local indices to the local nodes.
    unsigned int local_offset = 0;
    std::map<int, int> lag_idx_to_petsc_idx;
    for (auto p = level->begin(); p != level->end(); ++p)
    {
        const auto patch =*p;
        const Box& patch_box = patch->getBox();
        const boost::shared_ptr<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
        for (LNodeSetData::DataIterator it = idx_data->data_begin(patch_box); it != idx_data->data_end(); ++it)
        {
            LNode* const node_idx = *it;
            const int lag_idx = node_idx->getLagrangianIndex();
            local_lag_indices.push_back(lag_idx);
            const int petsc_idx = local_offset++;
            node_idx->setLocalPETScIndex(petsc_idx);
            lag_idx_to_petsc_idx[lag_idx] = petsc_idx;
        }
    }

    // Determine the Lagrangian indices of the nonlocal nodes.
    for (auto p = level->begin(); p != level->end(); ++p)
    {
        const auto patch =*p;
        const Box& patch_box = patch->getBox();
        const boost::shared_ptr<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
        BoxContainer ghost_boxes(idx_data->getGhostBox());
        ghost_boxes.removeIntersections(patch_box);
        for (auto bl(ghost_boxes); bl; bl++)
        {
            for (LNodeSetData::DataIterator it = idx_data->data_begin(bl()); it != idx_data->data_end(); ++it)
            {
                LNode* const node_idx = *it;
                const int lag_idx = node_idx->getLagrangianIndex();
                std::map<int, int>::const_iterator idx_it = lag_idx_to_petsc_idx.find(lag_idx);
                if (idx_it == lag_idx_to_petsc_idx.end())
                {
                    // This is the first time we have encountered this index; it
                    // must be a nonlocal index.
                    nonlocal_lag_indices.push_back(lag_idx);
                    const int petsc_idx = local_offset++;
                    node_idx->setLocalPETScIndex(petsc_idx);
                    lag_idx_to_petsc_idx[lag_idx] = petsc_idx;
                }
                else
                {
                    node_idx->setLocalPETScIndex(idx_it->second);
                }
            }
        }
    }

    // Compute the new PETSc global ordering and initialize the AO object.
    int ierr;

    // Determine how many nodes are on each processor to calculate the PETSc
    // indexing scheme.
    const unsigned int num_local_nodes = static_cast<unsigned int>(local_lag_indices.size());
    const unsigned int num_nonlocal_nodes = static_cast<unsigned int>(nonlocal_lag_indices.size());

    if (local_offset != (num_local_nodes + num_nonlocal_nodes))
    {
        TBOX_ERROR("LDataManager::computeNodeDistribution()"
                   << "\n"
                   << "  local_offset       = " << local_offset << "\n"
                   << "  num_local_nodes    = " << num_local_nodes << "\n"
                   << "  num_nonlocal_nodes = " << num_nonlocal_nodes << "\n");
    }

    computeNodeOffsets(num_nodes, node_offset, num_local_nodes);

    // Determine the PETSc ordering and setup the new AO object.
    const int num_proc_nodes = num_local_nodes + num_nonlocal_nodes;

    std::vector<int> node_indices;
    node_indices.reserve(num_proc_nodes);
    node_indices.insert(node_indices.end(), local_lag_indices.begin(), local_lag_indices.end());

    local_petsc_indices.resize(num_local_nodes);
    for (unsigned int k = 0; k < num_local_nodes; ++k)
    {
        local_petsc_indices[k] = node_offset + k;
    }

    if (ao)
    {
        ierr = AODestroy(&ao);
        IBTK_CHKERRQ(ierr);
    }

    ierr = AOCreateMapping(PETSC_COMM_WORLD,
                           num_local_nodes,
                           num_local_nodes > 0 ? &node_indices[0] : NULL,
                           num_local_nodes > 0 ? &local_petsc_indices[0] : NULL,
                           &ao);
    IBTK_CHKERRQ(ierr);

    // Determine the PETSc local to global mapping (including PETSc Vec ghost
    // indices).
    //
    // NOTE: After this operation, data stored in node_indices are in the global
    // PETSc ordering.
    node_indices.reserve(node_indices.size() + nonlocal_lag_indices.size());
    node_indices.insert(node_indices.end(), nonlocal_lag_indices.begin(), nonlocal_lag_indices.end());
    ierr = AOApplicationToPetsc(ao,
                                (num_proc_nodes > 0 ? num_proc_nodes : static_cast<int>(s_ao_dummy.size())),
                                (num_proc_nodes > 0 ? &node_indices[0] : &s_ao_dummy[0]));
    IBTK_CHKERRQ(ierr);

    // Keep track of the global PETSc indices of the ghost nodes.
    nonlocal_petsc_indices.clear();
    nonlocal_petsc_indices.reserve(num_nonlocal_nodes);
    nonlocal_petsc_indices.insert(
        nonlocal_petsc_indices.end(), node_indices.begin() + num_local_nodes, node_indices.end());

    // Store the global PETSc index in the local LNode objects.
    for (auto p = level->begin(); p != level->end(); ++p)
    {
        const auto patch =*p;
        const boost::shared_ptr<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
        const Box& ghost_box = idx_data->getGhostBox();
        for (LNodeSetData::DataIterator it = idx_data->data_begin(ghost_box); it != idx_data->data_end(); ++it)
        {
            LNode* const node_idx = *it;
            node_idx->setGlobalPETScIndex(node_indices[node_idx->getLocalPETScIndex()]);
        }
    }

    IBTK_TIMER_STOP(t_compute_node_distribution);
    return;
} // computeNodeDistribution

void
LDataManager::computeNodeOffsets(unsigned int& num_nodes, unsigned int& node_offset, const unsigned int num_local_nodes)
{
    IBTK_TIMER_START(t_compute_node_offsets);

    tbox::SAMRAI_MPI comm(MPI_COMM_WORLD);
    const int mpi_size = comm.getSize();
    const int mpi_rank = comm.getRank();

    std::vector<int> num_nodes_proc(mpi_size, 0);

    int num_local_nodes_int = num_local_nodes;
    comm.Allgather(&num_local_nodes_int, 1, MPI_INT, &num_nodes_proc[0], mpi_size, MPI_INT);

    node_offset = std::accumulate(num_nodes_proc.begin(), num_nodes_proc.begin() + mpi_rank, 0);

    num_nodes = std::accumulate(num_nodes_proc.begin() + mpi_rank, num_nodes_proc.end(), node_offset);

    IBTK_TIMER_STOP(t_compute_node_offsets);
    return;
} // computeNodeOffsets

void LDataManager::getFromRestart()
{
    auto restart_db = RestartManager::getManager()->getRootDatabase();

    boost::shared_ptr<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to " << d_object_name << " not found in restart file.");
    }

    int ver = db->getInteger("LDATA_MANAGER_VERSION");
    if (ver != LDATA_MANAGER_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                                 << "Restart file version different than class version.");
    }

    d_coarsest_ln = db->getInteger("d_coarsest_ln");
    d_finest_ln = db->getInteger("d_finest_ln");
    d_beta_work = db->getDouble("d_beta_work");

    // Resize some arrays.
    d_level_contains_lag_data.resize(d_finest_ln + 1, false);
    d_strct_name_to_strct_id_map.resize(d_finest_ln + 1);
    d_strct_id_to_strct_name_map.resize(d_finest_ln + 1);
    d_strct_id_to_lag_idx_range_map.resize(d_finest_ln + 1);
    d_last_lag_idx_to_strct_id_map.resize(d_finest_ln + 1);
    d_inactive_strcts.resize(d_finest_ln + 1);
    d_lag_mesh.resize(d_finest_ln + 1);
    d_lag_mesh_data.resize(d_finest_ln + 1);
    d_needs_synch.resize(d_finest_ln + 1, false);
    d_ao.resize(d_finest_ln + 1);
    d_num_nodes.resize(d_finest_ln + 1);
    d_node_offset.resize(d_finest_ln + 1);
    d_local_lag_indices.resize(d_finest_ln + 1);
    d_nonlocal_lag_indices.resize(d_finest_ln + 1);
    d_local_petsc_indices.resize(d_finest_ln + 1);
    d_nonlocal_petsc_indices.resize(d_finest_ln + 1);

    // Read in data that is stored on a level-by-level basis.
    for (int level_number = d_coarsest_ln; level_number <= d_finest_ln; ++level_number)
    {
        std::ostringstream stream;
        stream << "level_" << level_number;
        const std::string level_db_name = stream.str();
        auto level_db = db->getDatabase(level_db_name);

        d_level_contains_lag_data[level_number] = level_db->getBool("d_level_contains_lag_data");

        if (!d_level_contains_lag_data[level_number]) continue;

        const int n_lstructs = level_db->getInteger("n_lstructs");
        std::vector<int> lstruct_ids(n_lstructs), lstruct_lag_idx_range_first(n_lstructs),
            lstruct_lag_idx_range_second(n_lstructs), lstruct_activation(n_lstructs);
        std::vector<std::string> lstruct_names(n_lstructs);
        if (n_lstructs > 0)
        {
            level_db->getIntegerArray("lstruct_ids", &lstruct_ids[0], static_cast<int>(lstruct_ids.size()));
            level_db->getIntegerArray("lstruct_lag_idx_range_first",
                                      &lstruct_lag_idx_range_first[0],
                                      static_cast<int>(lstruct_lag_idx_range_first.size()));
            level_db->getIntegerArray("lstruct_lag_idx_range_second",
                                      &lstruct_lag_idx_range_second[0],
                                      static_cast<int>(lstruct_lag_idx_range_second.size()));
            level_db->getIntegerArray(
                "lstruct_activation", &lstruct_activation[0], static_cast<int>(lstruct_activation.size()));
            level_db->getStringArray("lstruct_names", &lstruct_names[0], static_cast<int>(lstruct_names.size()));
        }
        for (int k = 0; k < n_lstructs; ++k)
        {
            d_strct_id_to_strct_name_map[level_number][lstruct_ids[k]] = lstruct_names[k];
            d_strct_id_to_lag_idx_range_map[level_number][lstruct_ids[k]] =
                std::make_pair(lstruct_lag_idx_range_first[k], lstruct_lag_idx_range_second[k]);
            if (lstruct_activation[k] == 1)
            {
                d_inactive_strcts[level_number].insert(k);
            }
        }

        for (auto cit(d_strct_id_to_strct_name_map[level_number].begin());
             cit != d_strct_id_to_strct_name_map[level_number].end();
             ++cit)
        {
            d_strct_name_to_strct_id_map[level_number][cit->second] = cit->first;
        }

        for (auto cit(
                 d_strct_id_to_lag_idx_range_map[level_number].begin());
             cit != d_strct_id_to_lag_idx_range_map[level_number].end();
             ++cit)
        {
            d_last_lag_idx_to_strct_id_map[level_number][cit->second.second - 1] = cit->first;
        }

        const int n_ldata_names = level_db->getInteger("n_ldata_names");
        std::vector<std::string> ldata_names(n_ldata_names);
        if (!ldata_names.empty())
        {
            level_db->getStringArray("ldata_names", &ldata_names[0], n_ldata_names);
        }

        std::set<int> data_depths;
        for (auto it = ldata_names.begin(); it != ldata_names.end(); ++it)
        {
            const std::string& ldata_name = *it;
            d_lag_mesh_data[level_number][ldata_name] = boost::make_shared<LData(level_db->getDatabase>(ldata_name));
            data_depths.insert(d_lag_mesh_data[level_number][ldata_name]->getDepth());
        }

        d_num_nodes[level_number] = level_db->getInteger("d_num_nodes");
        d_node_offset[level_number] = level_db->getInteger("d_node_offset");

        const int n_local_lag_indices = level_db->getInteger("n_local_lag_indices");
        if (n_local_lag_indices > 0)
        {
            d_local_lag_indices[level_number].resize(n_local_lag_indices);
            level_db->getIntegerArray(
                "d_local_lag_indices", &d_local_lag_indices[level_number][0], n_local_lag_indices);
        }
        const int n_nonlocal_lag_indices = level_db->getInteger("n_nonlocal_lag_indices");
        if (n_nonlocal_lag_indices > 0)
        {
            d_nonlocal_lag_indices[level_number].resize(n_nonlocal_lag_indices);
            level_db->getIntegerArray(
                "d_nonlocal_lag_indices", &d_nonlocal_lag_indices[level_number][0], n_nonlocal_lag_indices);
        }
        const int n_local_petsc_indices = level_db->getInteger("n_local_petsc_indices");
        if (n_local_petsc_indices > 0)
        {
            d_local_petsc_indices[level_number].resize(n_local_petsc_indices);
            level_db->getIntegerArray(
                "d_local_petsc_indices", &d_local_petsc_indices[level_number][0], n_local_petsc_indices);
        }
        const int n_nonlocal_petsc_indices = level_db->getInteger("n_nonlocal_petsc_indices");
        if (n_nonlocal_petsc_indices > 0)
        {
            d_nonlocal_petsc_indices[level_number].resize(n_nonlocal_petsc_indices);
            level_db->getIntegerArray(
                "d_nonlocal_petsc_indices", &d_nonlocal_petsc_indices[level_number][0], n_nonlocal_petsc_indices);
        }

        // Rebuild the application ordering.
        int ierr;
        ierr = AOCreateMapping(PETSC_COMM_WORLD,
                               n_local_lag_indices,
                               n_local_lag_indices > 0 ? &d_local_lag_indices[level_number][0] : NULL,
                               n_local_lag_indices > 0 ? &d_local_petsc_indices[level_number][0] : NULL,
                               &d_ao[level_number]);
        IBTK_CHKERRQ(ierr);
    }
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
