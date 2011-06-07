// Filename: LDataManager.C
// Created on 01 Mar 2004 by Boyce Griffith
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

#include "LDataManager.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LData.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/LNodeSet.h>
#include <ibtk/LNodeSetData.h>
#include <ibtk/LNodeTransaction.h>
#include <ibtk/LSiloDataWriter.h>
#include <ibtk/namespaces.h>

#if (NDIM == 3)
#include <ibtk/LM3DDataWriter.h>
#endif

// SAMRAI INCLUDES
#include <Box.h>
#include <BoxList.h>
#include <CartesianCellDoubleWeightedAverage.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <CoarsenOperator.h>
#include <HierarchyDataOpsManager.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <RefineOperator.h>
#include <VariableDatabase.h>
#include <tbox/MathUtilities.h>
#include <tbox/RestartManager.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/ShutdownRegistry.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>
#include <set>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_spread;
static Timer* t_interp;
static Timer* t_map_lagrangian_to_petsc;
static Timer* t_map_petsc_to_lagrangian;
static Timer* t_begin_data_redistribution;
static Timer* t_end_data_redistribution;
static Timer* t_update_workload_data;
static Timer* t_initialize_level_data;
static Timer* t_reset_hierarchy_configuration;
static Timer* t_apply_gradient_detector;
static Timer* t_put_to_database;
static Timer* t_begin_nonlocal_data_fill;
static Timer* t_end_nonlocal_data_fill;
static Timer* t_compute_node_distribution;
static Timer* t_compute_node_offsets;

// Assume max(U)dt/dx <= 2.
static const int CFL_WIDTH = 2;

// Version of LDataManager restart file data.
static const int LDATA_MANAGER_VERSION = 1;

inline CellIndex<NDIM>
get_canonical_cell_index(
    const CellIndex<NDIM>& cell_idx,
    const Box<NDIM>& domain_box,
    const IntVector<NDIM>& periodic_shift)
{
    if (periodic_shift == IntVector<NDIM>(0)) return cell_idx;

    if (domain_box.contains(cell_idx)) return cell_idx;

    CellIndex<NDIM> shifted_idx = cell_idx;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if      (shifted_idx(d) < domain_box.lower()(d)) shifted_idx(d) += periodic_shift(d);
        else if (shifted_idx(d) > domain_box.upper()(d)) shifted_idx(d) -= periodic_shift(d);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(domain_box.contains(shifted_idx));
#endif
    return shifted_idx;
}// get_canonical_cell_index
}

const std::string LDataManager::     POSN_DATA_NAME = "X";
const std::string LDataManager::INIT_POSN_DATA_NAME = "X0";
const std::string LDataManager::      VEL_DATA_NAME = "U";
std::map<std::string,LDataManager*> LDataManager::s_data_manager_instances;
bool LDataManager::s_registered_callback;
unsigned char LDataManager::s_shutdown_priority = 200;
std::vector<int> LDataManager::s_ao_dummy(1,-1);

LDataManager*
LDataManager::getManager(
    const std::string& name,
    const std::string& interp_weighting_fcn,
    const std::string& spread_weighting_fcn,
    const IntVector<NDIM>& ghost_cell_width,
    bool register_for_restart)
{
    if (s_data_manager_instances.find(name) == s_data_manager_instances.end())
    {
        const int stencil_size = std::max(LEInteractor::getStencilSize(interp_weighting_fcn),LEInteractor::getStencilSize(interp_weighting_fcn));
        const IntVector<NDIM> gcw = IntVector<NDIM>::max(IntVector<NDIM>(static_cast<int>(floor(0.5*static_cast<double>(stencil_size)))+1),ghost_cell_width);
        s_data_manager_instances[name] = new LDataManager(name, interp_weighting_fcn, spread_weighting_fcn, gcw, register_for_restart);
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeAllManagers, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_data_manager_instances[name];
}// getManager

void
LDataManager::freeAllManagers()
{
    for (std::map<std::string,LDataManager*>::iterator it = s_data_manager_instances.begin();
         it != s_data_manager_instances.end(); ++it)
    {
        if (it->second)
        {
            delete it->second;
        }
        it->second = NULL;
    }
    return;
}// freeManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
LDataManager::setPatchHierarchy(
    Pointer<PatchHierarchy<NDIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
#endif

    // Reset the hierarchy.
    d_hierarchy = hierarchy;
    d_grid_geom = hierarchy->getGridGeometry();
    return;
}// setPatchHierarchy

void
LDataManager::resetLevels(
    const int coarsest_ln,
    const int finest_ln)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((coarsest_ln >= 0) &&
                (finest_ln >= coarsest_ln) &&
                (finest_ln <= d_hierarchy->getFinestLevelNumber()));
#endif
    // Destroy any unneeded AO objects.
    int ierr;
    for (int level_number = std::max(d_coarsest_ln,0); (level_number <= d_finest_ln) && (level_number < coarsest_ln); ++level_number)
    {
        if (d_ao[level_number])
        {
            ierr = AODestroy(d_ao[level_number]);
            IBTK_CHKERRQ(ierr);
        }
    }
    for (int level_number = finest_ln+1; level_number <= d_finest_ln; ++level_number)
    {
        if (d_ao[level_number])
        {
            ierr = AODestroy(d_ao[level_number]);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Reset the level numbers.
    d_coarsest_ln = coarsest_ln;
    d_finest_ln   = finest_ln;

    // Resize some arrays.
    d_level_contains_lag_data       .resize(d_finest_ln+1);
    d_strct_name_to_strct_id_map    .resize(d_finest_ln+1);
    d_strct_id_to_strct_name_map    .resize(d_finest_ln+1);
    d_strct_id_to_lag_idx_range_map .resize(d_finest_ln+1);
    d_last_lag_idx_to_strct_id_map  .resize(d_finest_ln+1);
    d_inactive_strcts               .resize(d_finest_ln+1);
    d_displaced_strct_ids           .resize(d_finest_ln+1);
    d_displaced_strct_bounding_boxes.resize(d_finest_ln+1);
    d_displaced_strct_lnode_idxs    .resize(d_finest_ln+1);
    d_displaced_strct_lnode_posns   .resize(d_finest_ln+1);
    d_lag_mesh                      .resize(d_finest_ln+1);
    d_lag_mesh_data                 .resize(d_finest_ln+1);
    d_needs_synch                   .resize(d_finest_ln+1,false);
    d_ao                            .resize(d_finest_ln+1);
    d_num_nodes                     .resize(d_finest_ln+1);
    d_node_offset                   .resize(d_finest_ln+1);
    d_local_lag_indices             .resize(d_finest_ln+1);
    d_nonlocal_lag_indices          .resize(d_finest_ln+1);
    d_local_petsc_indices           .resize(d_finest_ln+1);
    d_nonlocal_petsc_indices        .resize(d_finest_ln+1);
    return;
}// resetLevels

void
LDataManager::spread(
    const int f_data_idx,
    std::vector<Pointer<LData> >& F_data,
    std::vector<Pointer<LData> >& X_data,
    std::vector<Pointer<LData> >& ds_data,
    std::vector<Pointer<RefineSchedule<NDIM> > > f_prolongation_scheds,
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
    std::vector<Pointer<LData> > F_ds_data(F_data.size());
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (levelContainsLagrangianData(ln))
        {
            const int depth = F_data[ln]->getDepth();
            F_ds_data[ln] = new LData("", getNumberOfLocalNodes(ln), depth, d_nonlocal_petsc_indices[ln][depth]);
            blitz::Array<double,2>&       F_ds_arr = *F_ds_data[ln]->getGhostedLocalFormVecArray();
            const blitz::Array<double,2>&    F_arr = *   F_data[ln]->getGhostedLocalFormVecArray();
            const blitz::Array<double,1>&   ds_arr = *  ds_data[ln]->getGhostedLocalFormArray();
            for (int k = 0; k < static_cast<int>(F_data[ln]->getLocalNodeCount() + F_data[ln]->getGhostNodeCount()); ++k)
            {
                for (int d = 0; d < depth; ++d)
                {
                    F_ds_arr(k,d) = F_arr(k,d)*ds_arr(k);
                }
            }
            F_ds_data[ln]->restoreArrays();
            F_data   [ln]->restoreArrays();
            ds_data  [ln]->restoreArrays();
        }
    }

    IBTK_TIMER_STOP(t_spread);

    // Spread data from the Lagrangian mesh to the Eulerian grid.
    spread(f_data_idx, F_ds_data, X_data, f_prolongation_scheds, F_data_ghost_node_update || ds_data_ghost_node_update, X_data_ghost_node_update, coarsest_ln, finest_ln);
    return;
}// spread

void
LDataManager::spread(
    const int f_data_idx,
    std::vector<Pointer<LData> >& F_data,
    std::vector<Pointer<LData> >& X_data,
    std::vector<Pointer<RefineSchedule<NDIM> > > f_prolongation_scheds,
    const bool F_data_ghost_node_update,
    const bool X_data_ghost_node_update,
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    IBTK_TIMER_START(t_spread);

    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? d_hierarchy->getFinestLevelNumber() : finest_ln_in);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Used specialized spreading routine if possible.
    bool use_spread_specialized = true;
    for (int ln = coarsest_ln; ln < finest_ln && use_spread_specialized; ++ln)
    {
        if (levelContainsLagrangianData(ln)) use_spread_specialized = false;
    }
    if (use_spread_specialized)
    {
        spread_specialized(f_data_idx, F_data, X_data, F_data_ghost_node_update, X_data_ghost_node_update, coarsest_ln, finest_ln);
    }
    else
    {
        // Make a copy of the Eulerian data.
        Pointer<Variable<NDIM> > f_var;
        var_db->mapIndexToVariable(f_data_idx, f_var);
        const int f_copy_data_idx = var_db->registerClonedPatchDataIndex(f_var, f_data_idx);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(f_copy_data_idx);
        }
        Pointer<HierarchyDataOpsReal<NDIM,double> > f_data_ops = HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(f_var, d_hierarchy, true);
        f_data_ops->resetLevels(coarsest_ln, finest_ln);
        f_data_ops->copyData(f_copy_data_idx, f_data_idx);

        // Start filling Lagrangian ghost node values.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (levelContainsLagrangianData(ln))
            {
                if (F_data_ghost_node_update) F_data[ln]->beginGhostUpdate();
                if (X_data_ghost_node_update) X_data[ln]->beginGhostUpdate();
            }
        }

        // Spread data from the Lagrangian mesh to the Eulerian grid.
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

            // On the coarsest level of the patch hierarchy, initialize the
            // Eulerian data to equal zero before spreading.
            //
            // For each of the finer levels in the patch hierarchy, prolong the
            // Eulerian data from coarser levels in the patch hierarchy before
            // spreading.
            if (ln > coarsest_ln && ln < static_cast<int>(f_prolongation_scheds.size()) && !f_prolongation_scheds[ln].isNull())
            {
                f_prolongation_scheds[ln]->fillData(0.0);
            }
            else
            {
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_data_idx);
                    Pointer<CellData<NDIM,double> > f_cc_data = f_data;
                    Pointer<SideData<NDIM,double> > f_sc_data = f_data;
                    const bool is_cc_data = !f_cc_data.isNull();
                    const bool is_sc_data = !f_sc_data.isNull();
                    if (is_cc_data) f_cc_data->fillAll(0.0);
                    if (is_sc_data) f_sc_data->fillAll(0.0);
                }
            }

            // Spread data onto the grid.
            if (levelContainsLagrangianData(ln))
            {
                if (F_data_ghost_node_update) F_data[ln]->endGhostUpdate();
                if (X_data_ghost_node_update) X_data[ln]->endGhostUpdate();
                const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_data_idx);
                    Pointer<CellData<NDIM,double> > f_cc_data = f_data;
                    Pointer<SideData<NDIM,double> > f_sc_data = f_data;
                    const bool is_cc_data = !f_cc_data.isNull();
                    const bool is_sc_data = !f_sc_data.isNull();
                    const Pointer<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
                    const Box<NDIM>& box = idx_data->getGhostBox();
                    if (is_cc_data) LEInteractor::spread(f_cc_data, F_data[ln], X_data[ln], idx_data, patch, box, periodic_shift, d_spread_weighting_fcn);
                    if (is_sc_data) LEInteractor::spread(f_sc_data, F_data[ln], X_data[ln], idx_data, patch, box, periodic_shift, d_spread_weighting_fcn);
                }
            }
        }

        // Accumulate data.
        f_data_ops->add(f_data_idx, f_data_idx, f_copy_data_idx);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(f_copy_data_idx);
        }
        var_db->removePatchDataIndex(f_copy_data_idx);
    }

    IBTK_TIMER_STOP(t_spread);
    return;
}// spread

void
LDataManager::interp(
    const int f_data_idx,
    std::vector<Pointer<LData> >& F_data,
    std::vector<Pointer<LData> >& X_data,
    std::vector<Pointer<CoarsenSchedule<NDIM> > > f_synch_scheds,
    std::vector<Pointer<RefineSchedule<NDIM> > > f_ghost_fill_scheds,
    const double fill_data_time,
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    IBTK_TIMER_START(t_interp);

    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? d_hierarchy->getFinestLevelNumber() : finest_ln_in);

    // Used specialized spreading routine if possible.
    bool use_interp_specialized = true;
    for (int ln = coarsest_ln; ln < finest_ln && use_interp_specialized; ++ln)
    {
        if (levelContainsLagrangianData(ln))
        {
            use_interp_specialized = false;
        }
    }
    if (use_interp_specialized)
    {
        interp_specialized(f_data_idx, F_data, X_data, f_ghost_fill_scheds, fill_data_time, coarsest_ln, finest_ln);
    }
    else
    {
        // Synchronize Eulerian values.
        for (int ln = finest_ln; ln > coarsest_ln; --ln)
        {
            if (ln < static_cast<int>(f_synch_scheds.size()) && !f_synch_scheds[ln].isNull())
            {
                f_synch_scheds[ln]->coarsenData();
            }
        }

        // Interpolate data from the Eulerian grid to the Lagrangian mesh.
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
            if (levelContainsLagrangianData(ln))
            {
                if (ln < static_cast<int>(f_ghost_fill_scheds.size()) && !f_ghost_fill_scheds[ln].isNull())
                {
                    f_ghost_fill_scheds[ln]->fillData(fill_data_time);
                }
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_data_idx);
                    Pointer<CellData<NDIM,double> > f_cc_data = f_data;
                    Pointer<SideData<NDIM,double> > f_sc_data = f_data;
                    const bool is_cc_data = !f_cc_data.isNull();
                    const bool is_sc_data = !f_sc_data.isNull();
                    const Pointer<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
                    const Box<NDIM>& box = idx_data->getBox();
                    if (is_cc_data) LEInteractor::interpolate(F_data[ln], X_data[ln], idx_data, f_cc_data, patch, box, periodic_shift, d_interp_weighting_fcn);
                    if (is_sc_data) LEInteractor::interpolate(F_data[ln], X_data[ln], idx_data, f_sc_data, patch, box, periodic_shift, d_interp_weighting_fcn);
                }
            }
        }

        // Zero inactivated components.
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            zeroInactivatedComponents(F_data[ln], ln);
        }
    }

    IBTK_TIMER_STOP(t_interp);
    return;
}// interp

void
LDataManager::registerLInitStrategy(
    Pointer<LInitStrategy> lag_init)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!lag_init.isNull());
#endif
    d_lag_init = lag_init;
    return;
}// registerLInitStrategy

void
LDataManager::freeLInitStrategy()
{
    d_lag_init.setNull();
    return;
}// freeLInitStrategy

void
LDataManager::registerVisItDataWriter(
    Pointer<VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!visit_writer.isNull());
#endif
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
}// registerVisItDataWriter

void
LDataManager::registerLSiloDataWriter(
    Pointer<LSiloDataWriter> silo_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!silo_writer.isNull());
#endif
    d_silo_writer = silo_writer;
    return;
}// registerLSiloDataWriter

#if (NDIM == 3)
void
LDataManager::registerLM3DDataWriter(
    Pointer<LM3DDataWriter> m3D_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!m3D_writer.isNull());
#endif
    d_m3D_writer = m3D_writer;
    return;
}// registerLM3DDataWriter
#endif

void
LDataManager::registerLoadBalancer(
    Pointer<LoadBalancer<NDIM> > load_balancer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!load_balancer.isNull());
#endif
    d_load_balancer = load_balancer;
    return;
}// return

Pointer<LData>
LDataManager::createLData(
    const std::string& quantity_name,
    const int level_number,
    const unsigned int depth,
    const bool maintain_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!maintain_data || (d_lag_mesh_data[level_number].find(quantity_name) == d_lag_mesh_data[level_number].end()));
    TBOX_ASSERT(level_number >= 0);
    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);
    TBOX_ASSERT(depth > 0);
#endif

    if (d_nonlocal_petsc_indices[level_number].find(depth) == d_nonlocal_petsc_indices[level_number].end())
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(depth != 1);
#endif
        d_nonlocal_petsc_indices[level_number][depth].resize(d_nonlocal_petsc_indices[level_number][1].size());
        std::transform(d_nonlocal_petsc_indices[level_number][    1].begin(),
                       d_nonlocal_petsc_indices[level_number][    1].end(),
                       d_nonlocal_petsc_indices[level_number][depth].begin(),
                       std::bind2nd(std::multiplies<int>(),depth));
    }

    Pointer<LData> ret_val = new LData(quantity_name, getNumberOfLocalNodes(level_number), depth, d_nonlocal_petsc_indices[level_number][depth]);
    if (maintain_data)
    {
        d_lag_mesh_data[level_number][quantity_name] = ret_val;
    }
    return ret_val;
}// createLData

blitz::TinyVector<double,NDIM>
LDataManager::computeLagrangianStructureCenterOfMass(
    const int structure_id,
    const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif
    int node_counter = 0;
    blitz::TinyVector<double,NDIM> X_com(0.0);
    std::pair<int,int> lag_idx_range = getLagrangianStructureIndexRange(structure_id, level_number);

    const blitz::Array<double,2>& X_data = *d_lag_mesh_data[level_number][POSN_DATA_NAME]->getLocalFormVecArray();
    const Pointer<LMesh> mesh = getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getNodes();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const int lag_idx = node_idx->getLagrangianIndex();
        if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
        {
            ++node_counter;
            const int local_idx = node_idx->getLocalPETScIndex();
            const double* const X = &X_data(local_idx,0);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_com[d] += X[d];
            }
        }
    }
    d_lag_mesh_data[level_number][POSN_DATA_NAME]->restoreArrays();

    SAMRAI_MPI::sumReduction(&X_com[0],NDIM);
    node_counter = SAMRAI_MPI::sumReduction(node_counter);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        X_com[d] /= static_cast<double>(node_counter);
    }
    return X_com;
}// computeLagrangianStructureCenterOfMass

std::pair<blitz::TinyVector<double,NDIM>,blitz::TinyVector<double,NDIM> >
LDataManager::computeLagrangianStructureBoundingBox(
    const int structure_id,
    const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif
    blitz::TinyVector<double,NDIM> X_lower( (std::numeric_limits<double>::max()-sqrt(std::numeric_limits<double>::epsilon())));
    blitz::TinyVector<double,NDIM> X_upper(-(std::numeric_limits<double>::max()-sqrt(std::numeric_limits<double>::epsilon())));
    std::pair<int,int> lag_idx_range = getLagrangianStructureIndexRange(structure_id, level_number);

    const blitz::Array<double,2>& X_data = *d_lag_mesh_data[level_number][POSN_DATA_NAME]->getLocalFormVecArray();
    const Pointer<LMesh> mesh = getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getNodes();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const int lag_idx = node_idx->getLagrangianIndex();
        if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
        {
            const int local_idx = node_idx->getLocalPETScIndex();
            const double* const X = &X_data(local_idx,0);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_lower[d] = std::min(X_lower[d],X[d]);
                X_upper[d] = std::max(X_upper[d],X[d]);
            }
        }
    }
    d_lag_mesh_data[level_number][POSN_DATA_NAME]->restoreArrays();

    SAMRAI_MPI::minReduction(&X_lower[0],NDIM);
    SAMRAI_MPI::maxReduction(&X_upper[0],NDIM);
    return std::make_pair(X_lower,X_upper);
}// computeLagrangianStructureBoundingBox

void
LDataManager::reinitLagrangianStructure(
    const blitz::TinyVector<double,NDIM>& X_center,
    const int structure_id,
    const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif
    d_displaced_strct_ids[level_number].push_back(structure_id);

    // Compute the bounding box of the structure in its reference configuration.
    const blitz::Array<double,2>& X0_data = *d_lag_mesh_data[level_number][INIT_POSN_DATA_NAME]->getLocalFormVecArray();
    blitz::TinyVector<double,NDIM> X_lower( (std::numeric_limits<double>::max()-sqrt(std::numeric_limits<double>::epsilon())));
    blitz::TinyVector<double,NDIM> X_upper(-(std::numeric_limits<double>::max()-sqrt(std::numeric_limits<double>::epsilon())));
    std::pair<int,int> lag_idx_range = getLagrangianStructureIndexRange(structure_id, level_number);

    const Pointer<LMesh> mesh = getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getNodes();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const int lag_idx = node_idx->getLagrangianIndex();
        if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
        {
            const int local_idx = node_idx->getLocalPETScIndex();
            const double* const X0 = &X0_data(local_idx,0);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_lower[d] = std::min(X_lower[d],X0[d]);
                X_upper[d] = std::max(X_upper[d],X0[d]);
            }
        }
    }
    SAMRAI_MPI::minReduction(&X_lower[0],NDIM);
    SAMRAI_MPI::maxReduction(&X_upper[0],NDIM);
    std::pair<blitz::TinyVector<double,NDIM>,blitz::TinyVector<double,NDIM> > bounding_box = std::make_pair(X_lower,X_upper);

    // Compute the displacement.
    blitz::TinyVector<double,NDIM> dX(0.0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        dX[d] = X_center[d] - 0.5*(X_upper[d]+X_lower[d]);
    }

    // Compute the shifted bounding box.
    std::pair<blitz::TinyVector<double,NDIM>,blitz::TinyVector<double,NDIM> > shifted_bounding_box = bounding_box;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        shifted_bounding_box.first [d] += dX[d];
        shifted_bounding_box.second[d] += dX[d];
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    const double* const gridXLower = d_grid_geom->getXLower();
    const double* const gridXUpper = d_grid_geom->getXUpper();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(gridXLower[d] <= shifted_bounding_box.first [d]);
        TBOX_ASSERT(gridXUpper[d] >= shifted_bounding_box.second[d]);
    }
#endif
    d_displaced_strct_bounding_boxes[level_number].push_back(shifted_bounding_box);

    // For each node within the shifted structure: update the position of the
    // node; excise that node from the LNodeSetData patch data structure; and
    // collect that node into a set of shifted Lagrangian nodes.
    blitz::Array<double,2>& X_data = *d_lag_mesh_data[level_number][POSN_DATA_NAME]->getLocalFormVecArray();
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
        for (LNodeSetData::CellIterator it(patch_box); it; it++)
        {
            const Index<NDIM>& i = *it;
            LNodeSet* const node_set = idx_data->getItem(i);
            if (node_set != NULL)
            {
                LNodeSet::DataSet unmoved_idxs;
                unmoved_idxs.reserve(node_set->size());
                for (LNodeSet::iterator it = node_set->begin(); it != node_set->end(); ++it)
                {
                    LNodeSet::value_type& node_idx = *it;
                    const int lag_idx = node_idx->getLagrangianIndex();
                    if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                    {
                        const int local_idx = node_idx->getLocalPETScIndex();
                        double* const X = &X_data(local_idx,0);
                        const double* const X0 = &X0_data(local_idx,0);
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            X[d] = X0[d] + dX[d];
                        }
                        d_displaced_strct_lnode_idxs [level_number].push_back(node_idx);
                        blitz::TinyVector<double,NDIM> X_displaced;
                        for (unsigned int d = 0; d < NDIM; ++d) X_displaced[d] = X[d];
                        d_displaced_strct_lnode_posns[level_number].push_back(X_displaced);
                    }
                    else
                    {
                        unmoved_idxs.push_back(node_idx);
                    }
                }
                node_set->setDataSet(unmoved_idxs);
                if (node_set->empty()) idx_data->removeItem(i);
            }
        }
    }
    d_lag_mesh_data[level_number][     POSN_DATA_NAME]->restoreArrays();
    d_lag_mesh_data[level_number][INIT_POSN_DATA_NAME]->restoreArrays();
    return;
}// reinitLagrangianStructure

void
LDataManager::displaceLagrangianStructure(
    const blitz::TinyVector<double,NDIM>& dX,
    const int structure_id,
    const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif
    d_displaced_strct_ids[level_number].push_back(structure_id);

    // Compute the shifted bounding box.
    std::pair<blitz::TinyVector<double,NDIM>,blitz::TinyVector<double,NDIM> > bounding_box = computeLagrangianStructureBoundingBox(structure_id, level_number);
    std::pair<blitz::TinyVector<double,NDIM>,blitz::TinyVector<double,NDIM> > shifted_bounding_box = bounding_box;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        shifted_bounding_box.first [d] += dX[d];
        shifted_bounding_box.second[d] += dX[d];
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    const double* const gridXLower = d_grid_geom->getXLower();
    const double* const gridXUpper = d_grid_geom->getXUpper();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(gridXLower[d] <= shifted_bounding_box.first [d]);
        TBOX_ASSERT(gridXUpper[d] >= shifted_bounding_box.second[d]);
    }
#endif
    d_displaced_strct_bounding_boxes[level_number].push_back(shifted_bounding_box);

    // For each node within the shifted structure: update the position of the
    // node; excise that node from the LNodeSetData patch data structure; and
    // collect that node into a set of shifted Lagrangian nodes.
    blitz::Array<double,2>& X_data = *d_lag_mesh_data[level_number][POSN_DATA_NAME]->getLocalFormVecArray();
    std::pair<int,int> lag_idx_range = getLagrangianStructureIndexRange(structure_id, level_number);
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
        for (LNodeSetData::CellIterator it(patch_box); it; it++)
        {
            const Index<NDIM>& i = *it;
            LNodeSet* const node_set = idx_data->getItem(i);
            if (node_set != NULL)
            {
                LNodeSet::DataSet unmoved_idxs;
                unmoved_idxs.reserve(node_set->size());
                for (LNodeSet::iterator it = node_set->begin(); it != node_set->end(); ++it)
                {
                    LNodeSet::value_type& node_idx = *it;
                    const int lag_idx = node_idx->getLagrangianIndex();
                    if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                    {
                        const int local_idx = node_idx->getLocalPETScIndex();
                        double* const X = &X_data(local_idx,0);
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            X[d] += dX[d];
                        }
                        d_displaced_strct_lnode_idxs [level_number].push_back(node_idx);
                        blitz::TinyVector<double,NDIM> X_displaced;
                        for (unsigned int d = 0; d < NDIM; ++d) X_displaced[d] = X[d];
                        d_displaced_strct_lnode_posns[level_number].push_back(X_displaced);
                    }
                    else
                    {
                        unmoved_idxs.push_back(node_idx);
                    }
                }
                node_set->setDataSet(unmoved_idxs);
                if (node_set->empty()) idx_data->removeItem(i);
            }
        }
    }
    d_lag_mesh_data[level_number][POSN_DATA_NAME]->restoreArrays();
    return;
}// displaceLagrangianStructure

void
LDataManager::activateLagrangianStructures(
    const std::vector<int>& structure_ids,
    const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif
    for (std::vector<int>::const_iterator cit = structure_ids.begin();
         cit != structure_ids.end(); ++cit)
    {
        d_inactive_strcts[level_number].removeItem(*cit);
    }
    d_inactive_strcts[level_number].communicateData();
    return;
}// activateLagrangianStructures

void
LDataManager::inactivateLagrangianStructures(
    const std::vector<int>& structure_ids,
    const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif
    for (std::vector<int>::const_iterator cit = structure_ids.begin();
         cit != structure_ids.end(); ++cit)
    {
        d_inactive_strcts[level_number].addItem(*cit);
    }
    d_inactive_strcts[level_number].communicateData();
    return;
}// inactivateLagrangianStructures

void
LDataManager::zeroInactivatedComponents(
    Pointer<LData> lag_data,
    const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif

    if (LIKELY(d_inactive_strcts[level_number].getSet().empty())) return;

    // Construct a list of all of the inactivated Lagrangian indices.
    //
    // NOTE: The vector idxs is NOT a list of the LOCAL inactivated indices.
    // Instead, it is a list of ALL of the inactivated indices.  Thus, idxs will
    // have the same contents for all MPI processes.
    std::vector<int> idxs;
    for (std::set<int>::const_iterator cit = d_inactive_strcts[level_number].getSet().begin();
         cit != d_inactive_strcts[level_number].getSet().end(); ++cit)
    {
        const int strct_id = *cit;
        const std::pair<int,int>& lag_index_range = d_strct_id_to_lag_idx_range_map[level_number].find(strct_id)->second;
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
    ierr = VecGetOwnershipRange(lag_data_vec, &lo, &hi);  IBTK_CHKERRQ(ierr);
    ierr = VecGetBlockSize(lag_data_vec, &bs);  IBTK_CHKERRQ(ierr);

    // Zero-out all local inactivated components.
    std::vector<int> ix(bs);
    std::vector<double> y(bs,0.0);
    std::sort(idxs.begin(), idxs.end());
    for (std::vector<int>::const_iterator cit = idxs.begin(); cit != idxs.end(); ++cit)
    {
        const int l = *cit;
        if (l*bs >= lo && l*bs < hi)
        {
            for (int k = 0; k < bs; ++k)
            {
                ix[k] = bs*l + k;
            }
            ierr = VecSetValues(lag_data_vec, bs, &ix[0], &y[0], INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        }
    }
    ierr = VecAssemblyBegin(lag_data_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(  lag_data_vec);  IBTK_CHKERRQ(ierr);
    lag_data->beginGhostUpdate();
    lag_data->  endGhostUpdate();
    return;
}// zeroInactivatedComponents

void
LDataManager::mapLagrangianToPETSc(
    std::vector<int>& inds,
    const int level_number) const
{
    IBTK_TIMER_START(t_map_lagrangian_to_petsc);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif

    const int ierr = AOApplicationToPetsc(
        d_ao[level_number],
        (!inds.empty() ? static_cast<int>(inds.size()) : static_cast<int>(s_ao_dummy.size())),
        (!inds.empty() ? &inds[0]                      : &s_ao_dummy[0]));
    IBTK_CHKERRQ(ierr);

    IBTK_TIMER_STOP(t_map_lagrangian_to_petsc);
    return;
}// mapLagrangianToPETSc

void
LDataManager::mapLagrangianToPETSc(
    blitz::Array<int,1>& inds,
    const int level_number) const
{
    IBTK_TIMER_START(t_map_lagrangian_to_petsc);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif

    const int ierr = AOApplicationToPetsc(
        d_ao[level_number],
        (inds.size() > 0 ? static_cast<int>(inds.size()) : static_cast<int>(s_ao_dummy.size())),
        (inds.size() > 0 ? inds.data()                   : &s_ao_dummy[0]));
    IBTK_CHKERRQ(ierr);

    IBTK_TIMER_STOP(t_map_lagrangian_to_petsc);
    return;
}// mapLagrangianToPETSc

void
LDataManager::mapPETScToLagrangian(
    std::vector<int>& inds,
    const int level_number) const
{
    IBTK_TIMER_START(t_map_petsc_to_lagrangian);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif

    const int ierr = AOPetscToApplication(
        d_ao[level_number],
        (!inds.empty() ? static_cast<int>(inds.size()) : static_cast<int>(s_ao_dummy.size())),
        (!inds.empty() ? &inds[0]                      : &s_ao_dummy[0]));
    IBTK_CHKERRQ(ierr);

    IBTK_TIMER_STOP(t_map_petsc_to_lagrangian);
    return;
}// mapPETScToLagrangian

void
LDataManager::mapPETScToLagrangian(
    blitz::Array<int,1>& inds,
    const int level_number) const
{
    IBTK_TIMER_START(t_map_petsc_to_lagrangian);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif

    const int ierr = AOPetscToApplication(
        d_ao[level_number],
        (inds.size() > 0 ? static_cast<int>(inds.size()) : static_cast<int>(s_ao_dummy.size())),
        (inds.size() > 0 ? inds.data()                   : &s_ao_dummy[0]));
    IBTK_CHKERRQ(ierr);

    IBTK_TIMER_STOP(t_map_petsc_to_lagrangian);
    return;
}// mapPETScToLagrangian

void
LDataManager::scatterLagrangianToPETSc(
    Vec& lagrangian_vec,
    Vec& petsc_vec,
    const int level_number) const
{
    scatterData(lagrangian_vec,petsc_vec,level_number,SCATTER_REVERSE);
    return;
}// scatterLagrangianToPETSc

void
LDataManager::scatterPETScToLagrangian(
    Vec& petsc_vec,
    Vec& lagrangian_vec,
    const int level_number) const
{
    scatterData(lagrangian_vec,petsc_vec,level_number,SCATTER_FORWARD);
    return;
}// scatterPETScToLagrangian

void
LDataManager::scatterToAll(
    Vec& parallel_vec,
    Vec& sequential_vec) const
{
    int ierr;
    const bool create_vout = sequential_vec != PETSC_NULL;
    VecScatter ctx;
    ierr = VecScatterCreateToAll(parallel_vec, &ctx, (create_vout ? &sequential_vec : PETSC_NULL));  IBTK_CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, parallel_vec, sequential_vec, INSERT_VALUES, SCATTER_FORWARD);  IBTK_CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, parallel_vec, sequential_vec, INSERT_VALUES, SCATTER_FORWARD);  IBTK_CHKERRQ(ierr);
    ierr = VecScatterDestroy(ctx);  IBTK_CHKERRQ(ierr);
    return;
}// scatterToAll

void
LDataManager::scatterToZero(
    Vec& parallel_vec,
    Vec& sequential_vec) const
{
    int ierr;
    const bool create_vout = sequential_vec != PETSC_NULL;
    VecScatter ctx;
    ierr = VecScatterCreateToZero(parallel_vec, &ctx, (create_vout ? &sequential_vec : PETSC_NULL));  IBTK_CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, parallel_vec, sequential_vec, INSERT_VALUES, SCATTER_FORWARD);  IBTK_CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, parallel_vec, sequential_vec, INSERT_VALUES, SCATTER_FORWARD);  IBTK_CHKERRQ(ierr);
    ierr = VecScatterDestroy(ctx);  IBTK_CHKERRQ(ierr);
    return;
}// scatterToZero

void
LDataManager::beginDataRedistribution(
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    IBTK_TIMER_START(t_begin_data_redistribution);

    const int coarsest_ln = (coarsest_ln_in == -1) ? d_coarsest_ln : coarsest_ln_in;
    const int finest_ln = (finest_ln_in == -1) ? d_finest_ln : finest_ln_in;

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(coarsest_ln >= d_coarsest_ln &&
                coarsest_ln <= d_finest_ln);
    TBOX_ASSERT(finest_ln   >= d_coarsest_ln &&
                finest_ln   <= d_finest_ln);
#endif

    const double* const gridXLower = d_grid_geom->getXLower();
    const double* const gridXUpper = d_grid_geom->getXUpper();
    blitz::TinyVector<double,NDIM> gridXLength;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        gridXLength[d] = gridXUpper[d] - gridXLower[d];
    }
    const IntVector<NDIM>& periodic_shift = d_grid_geom->getPeriodicShift();

    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (!d_level_contains_lag_data[level_number]) continue;

        if (d_needs_synch[level_number])
        {
            TBOX_WARNING("LDataManager::beginDataRedistribution():\n" <<
                         "\tLData is not synchronized with LNodeSetData.\n" <<
                         "\tLagrangian node position data is probably invalid!\n");
        }

        // Ensure that no IB points manage to escape the computational domain.
        //
        // NOTE: We cannot use the LMesh data structure here because it has not
        // been initialized.
        blitz::Array<double,2>& X_data = *d_lag_mesh_data[level_number][POSN_DATA_NAME]->getLocalFormVecArray();
        static const double edge_tol = sqrt(std::numeric_limits<double>::epsilon());
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM> patch_box = patch->getBox();
            Pointer<LNodeSetData> lag_node_idx_data = patch->getPatchData(d_lag_node_index_current_idx);
            for (LNodeSetData::DataIterator it = lag_node_idx_data->data_begin(patch_box);
                 it != lag_node_idx_data->data_end(); ++it)
            {
                LNodeSet::value_type& node_idx = *it;
                const int local_idx = node_idx->getLocalPETScIndex();
                double* const X = &X_data(local_idx,0);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    if (periodic_shift[d] == 0)
                    {
                        X[d] = std::max(X[d],gridXLower[d]+edge_tol);
                        X[d] = std::min(X[d],gridXUpper[d]-edge_tol);
                    }
                }
            }
        }
        d_lag_mesh_data[level_number][POSN_DATA_NAME]->restoreArrays();

        // Begin updating the ghost values of the Lagrangian nodal positions.
        d_lag_mesh_data[level_number][POSN_DATA_NAME]->beginGhostUpdate();
    }

    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (!d_level_contains_lag_data[level_number]) continue;

        // Finish updating the ghost values of the Lagrangian nodal positions.
        d_lag_mesh_data[level_number][POSN_DATA_NAME]->endGhostUpdate();

        // Update the index patch data on the level.
        blitz::Array<double,2>& X_data = *d_lag_mesh_data[level_number][POSN_DATA_NAME]->getGhostedLocalFormVecArray();
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const CellIndex<NDIM>& patch_lower = patch_box.lower();
            const CellIndex<NDIM>& patch_upper = patch_box.upper();
            const double* const patchXLower = patch_geom->getXLower();
            const double* const patchXUpper = patch_geom->getXUpper();
            const double* const patchDx = patch_geom->getDx();
            const bool touches_periodic_bdry = patch_geom->getTouchesPeriodicBoundary();

            Pointer<LNodeSetData> current_idx_data = patch->getPatchData(d_lag_node_index_current_idx);
            Pointer<LNodeSetData> new_idx_data = new LNodeSetData(current_idx_data->getBox(), current_idx_data->getGhostCellWidth());

            // Initialize LNodeSet objects for each cell index in the patch
            // interior that will contain the LNode objects AFTER
            // redistribution.
            //
            // We only keep nodes whose new locations are in the patch interior.
            // That is to say, we only keep the nodes that the patch will own
            // after redistribution.
            const Box<NDIM> cfl_box = Box<NDIM>::grow(patch_box, IntVector<NDIM>(CFL_WIDTH));
            for (LNodeSetData::CellIterator it(cfl_box); it; it++)
            {
                const Index<NDIM>& old_cell_idx = *it;
                LNodeSet* const old_node_set = current_idx_data->getItem(old_cell_idx);
                if (old_node_set != NULL)
                {
                    const bool patch_owns_node_at_old_loc = patch_box.contains(old_cell_idx);
                    const IntVector<NDIM>& periodic_offset = old_node_set->getPeriodicOffset();
                    blitz::TinyVector<double,NDIM> X_shifted;
                    for (LNodeSet::iterator n = old_node_set->begin(); n != old_node_set->end(); ++n)
                    {
                        LNodeSet::value_type& node_idx = *n;
                        const int local_idx = node_idx->getLocalPETScIndex();
                        double* const X = &X_data(local_idx,0);
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            X_shifted[d] = X[d] + static_cast<double>(periodic_offset(d))*patchDx[d];
                        }

                        const bool patch_owns_node_at_new_loc =
                            ((  patchXLower[0] <= X_shifted[0])&&(X_shifted[0] < patchXUpper[0]))
#if (NDIM > 1)
                            &&((patchXLower[1] <= X_shifted[1])&&(X_shifted[1] < patchXUpper[1]))
#if (NDIM > 2)
                            &&((patchXLower[2] <= X_shifted[2])&&(X_shifted[2] < patchXUpper[2]))
#endif
#endif
                            ;

                        if (patch_owns_node_at_new_loc)
                        {
                            if (periodic_offset != IntVector<NDIM>(0))
                            {
                                blitz::TinyVector<double,NDIM> displacement(0.0);
                                for (unsigned int d = 0; d < NDIM; ++d)
                                {
                                    displacement[d] = static_cast<double>(periodic_offset(d))*patchDx[d];
                                }
                                node_idx->registerPeriodicShift(periodic_offset,displacement);
                            }

                            const CellIndex<NDIM> new_cell_idx = IndexUtilities::getCellIndex(X_shifted, patchXLower, patchXUpper, patchDx, patch_lower, patch_upper);
                            if (!new_idx_data->isElement(new_cell_idx))
                            {
                                new_idx_data->appendItemPointer(new_cell_idx, new LNodeSet());
                            }
                            LNodeSet* const new_node_set = new_idx_data->getItem(new_cell_idx);
                            new_node_set->push_back(node_idx);
                        }

                        // If a node leaves the patch via a periodic boundary,
                        // we have to make sure that its location is properly
                        // updated.
                        //
                        // IMPORTANT NOTE: The location of a node is itself a
                        // *quantity* living on the Lagrangian mesh.  Until the
                        // Lagrangian data are redistributed, the patch still
                        // owns the LData associated with the old LNode
                        // distribution.  Thus it is the responsibility of this
                        // patch to update the location of the node if the node
                        // leaves via a periodic boundary.
                        if (touches_periodic_bdry && (patch_owns_node_at_old_loc || patch_owns_node_at_new_loc))
                        {
                            static const int lower = 0;
                            static const int upper = 1;
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                if      (patch_geom->getTouchesPeriodicBoundary(d,lower) && X[d] < gridXLower[d])
                                {
                                    X[d] += gridXLength[d];
                                }
                                else if (patch_geom->getTouchesPeriodicBoundary(d,upper) && X[d] >= gridXUpper[d])
                                {
                                    X[d] -= gridXLength[d];
                                }
                            }
                        }
                    }
                }
            }

            // Sort the LNode objects according to their Lagrangian indices.
            for (LNodeSetData::SetIterator it(*new_idx_data); it; it++)
            {
                LNodeSet& node_set = *it;
                std::sort(node_set.begin(), node_set.end(), LNodeIndexLagrangianIndexComp());
            }

            // Swap the old and new patch data pointers.
            patch->setPatchData(d_lag_node_index_current_idx, new_idx_data);
        }
        d_lag_mesh_data[level_number][POSN_DATA_NAME]->restoreArrays();
        d_needs_synch[level_number] = true;
    }

    IBTK_TIMER_STOP(t_begin_data_redistribution);
    return;
}// beginDataRedistribution

void
LDataManager::endDataRedistribution(
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    IBTK_TIMER_START(t_end_data_redistribution);

    const int coarsest_ln = (coarsest_ln_in == -1) ? d_coarsest_ln : coarsest_ln_in;
    const int finest_ln = (finest_ln_in == -1) ? d_finest_ln : finest_ln_in;

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(coarsest_ln >= d_coarsest_ln &&
                coarsest_ln <= d_finest_ln);
    TBOX_ASSERT(finest_ln   >= d_coarsest_ln &&
                finest_ln   <= d_finest_ln);
#endif

    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (d_level_contains_lag_data[level_number] && (!d_needs_synch[level_number]))
        {
            TBOX_WARNING("LDataManager::endDataRedistribution():\n" <<
                         "\tLData is already synchronized with LNodeSetData.\n" <<
                         "\tlevel = " << level_number << "\n");
        }
    }

    // Update parallel data structures to account for any displaced nodes.
    const double* const dx0 = d_grid_geom->getDx();
    const double* const gridXLower = d_grid_geom->getXLower();
    const double* const gridXUpper = d_grid_geom->getXUpper();
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (!d_level_contains_lag_data[level_number] || d_displaced_strct_ids[level_number].empty()) continue;

        const int num_procs = SAMRAI_MPI::getNodes();
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
        Pointer<BoxTree<NDIM> > box_tree = level->getBoxTree();
        const ProcessorMapping& processor_mapping = level->getProcessorMapping();
        const Box<NDIM>& domain_box = level->getPhysicalDomain()[0];
        const CellIndex<NDIM>& domain_lower = domain_box.lower();
        const CellIndex<NDIM>& domain_upper = domain_box.upper();
        const IntVector<NDIM>& ratio = level->getRatio();
        blitz::TinyVector<double,NDIM> dx;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dx[d] = dx0[d]/static_cast<double>(ratio(d));
        }

        // Determine which processor owns each of the local displaced nodes.
        size_t num_nodes = d_displaced_strct_lnode_idxs[level_number].size();
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(d_displaced_strct_lnode_posns[level_number].size() == num_nodes);
#endif
        typedef LNodeTransaction::LTransactionComponent LNodeTransactionComponent;
        std::vector<std::vector<LNodeTransactionComponent> > src_index_set(num_procs);
        for (size_t k = 0; k < num_nodes; ++k)
        {
            LNodeSet::value_type& lag_idx = d_displaced_strct_lnode_idxs[level_number][k];
            const blitz::TinyVector<double,NDIM>& posn = d_displaced_strct_lnode_posns[level_number][k];
            const CellIndex<NDIM> cell_idx = IndexUtilities::getCellIndex(
                posn, gridXLower, gridXUpper, dx.data(), domain_lower, domain_upper);

            Array<int> indices;
            box_tree->findOverlapIndices(indices, Box<NDIM>(cell_idx,cell_idx));
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(indices.getSize() == 1);
#endif
            const int patch_num = indices[0];
#if DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(patch_num >= 0 && patch_num < level->getNumberOfPatches());
#endif
            const int dst_proc = processor_mapping.getProcessorAssignment(patch_num);
            LNodeTransactionComponent component(lag_idx, posn);
            src_index_set[dst_proc].push_back(component);
        }

        // Setup communication transactions between each pair of processors.
        Schedule lnode_idx_data_mover;
        std::vector<std::vector<Pointer<Transaction> > > transactions(num_procs, std::vector<Pointer<Transaction> >(num_procs));
        for (int src_proc = 0; src_proc < num_procs; ++src_proc)
        {
            for (int dst_proc = 0; dst_proc < num_procs; ++dst_proc)
            {
                if (src_proc == SAMRAI_MPI::getRank())
                {
                    transactions[src_proc][dst_proc] = new LNodeTransaction(src_proc, dst_proc, src_index_set[dst_proc]);
                }
                else
                {
                    transactions[src_proc][dst_proc] = new LNodeTransaction(src_proc, dst_proc);
                }
                lnode_idx_data_mover.appendTransaction(transactions[src_proc][dst_proc]);
            }
        }

        // Communicate the data.
        lnode_idx_data_mover.communicate();

        // Clear the cached displaced nodes.
        d_displaced_strct_lnode_idxs [level_number].clear();
        d_displaced_strct_lnode_posns[level_number].clear();

        // Retrieve the communicated values.
        for (int src_proc = 0; src_proc < num_procs; ++src_proc)
        {
            for (int dst_proc = 0; dst_proc < num_procs; ++dst_proc)
            {
                if (dst_proc == SAMRAI_MPI::getRank())
                {
                    Pointer<LNodeTransaction> transaction = transactions[src_proc][dst_proc];
                    const std::vector<LNodeTransactionComponent>& dst_index_set = transaction->getDestinationData();
                    for (std::vector<LNodeTransactionComponent>::const_iterator
                             cit = dst_index_set.begin(); cit != dst_index_set.end(); ++cit)
                    {
                        d_displaced_strct_lnode_idxs [level_number].push_back(cit->item);
                        d_displaced_strct_lnode_posns[level_number].push_back(cit->posn);
                    }
                }
            }
        }

        // Determine which patch owns each of the local displaced nodes.
        num_nodes = d_displaced_strct_lnode_idxs[level_number].size();
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(d_displaced_strct_lnode_posns[level_number].size() == num_nodes);
#endif
        for (size_t k = 0; k < num_nodes; ++k)
        {
            const LNodeSet::value_type& lag_idx = d_displaced_strct_lnode_idxs[level_number][k];
            const blitz::TinyVector<double,NDIM>& posn = d_displaced_strct_lnode_posns[level_number][k];
            const CellIndex<NDIM> cell_idx = IndexUtilities::getCellIndex(
                posn, gridXLower, gridXUpper, dx.data(), domain_lower, domain_upper);

            Array<int> indices;
            box_tree->findOverlapIndices(indices, Box<NDIM>(cell_idx,cell_idx));
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(indices.getSize() == 1);
#endif
            const int patch_num = indices[0];
#if DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(patch_num >= 0 && patch_num < level->getNumberOfPatches());
#endif
            Pointer<Patch<NDIM> > patch = level->getPatch(patch_num);
            Pointer<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
            if (!idx_data->isElement(cell_idx))
            {
                idx_data->appendItemPointer(cell_idx, new LNodeSet());
            }
            LNodeSet* const node_set = idx_data->getItem(cell_idx);
            node_set->push_back(lag_idx);
        }

        // Clear all cached data associated with this patch level.
        d_displaced_strct_ids           [level_number].clear();
        d_displaced_strct_bounding_boxes[level_number].clear();
        d_displaced_strct_lnode_idxs    [level_number].clear();
        d_displaced_strct_lnode_posns   [level_number].clear();
    }

    // Fill the ghost cells of each level.
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (!d_level_contains_lag_data[level_number]) continue;

        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
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

    std::vector<AO> new_ao(finest_ln+1);

    std::vector<std::vector<Vec> > src_vec(finest_ln+1);
    std::vector<std::vector<Vec> > dst_vec(finest_ln+1);
    std::vector<std::vector<VecScatter> > scatter(finest_ln+1);
    std::vector<std::map<int,IS> > src_IS(finest_ln+1);
    std::vector<std::map<int,IS> > dst_IS(finest_ln+1);
    std::vector<std::map<int,VecScatter> > scatter_template(finest_ln+1);

    // The number of all local (e.g., on processor) and ghost (e.g., off
    // processor) nodes.
    //
    // NOTE:  num_local_nodes   [ln] == d_local_lag_indices   [ln].size()
    //        num_nonlocal_nodes[ln] == d_nonlocal_lag_indices[ln].size()
    std::vector<int> num_local_nodes   (finest_ln+1);
    std::vector<int> num_nonlocal_nodes(finest_ln+1);

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

        // Reset the nonlocal PETSc indices.
        d_nonlocal_petsc_indices[level_number].clear();

        // The destination indices.
        std::vector<int> dst_inds;
        bool dst_inds_set = false;

        // Compute the new data distribution and start scattering.
        std::map<std::string,Pointer<LData> >& level_data = d_lag_mesh_data[level_number];
        const std::vector<int>::size_type num_data = level_data.size();

        src_vec[level_number].resize(num_data);
        dst_vec[level_number].resize(num_data);
        scatter[level_number].resize(num_data);

        // Get the new distribution of nodes for the level.
        //
        // NOTE: This process updates the local PETSc indices of the LNodeSet
        // objects contained in the current patch.
        ierr = computeNodeDistribution(d_local_lag_indices   [level_number],
                                       d_nonlocal_lag_indices[level_number],
                                       new_ao[level_number],
                                       d_local_petsc_indices   [level_number],
                                       d_nonlocal_petsc_indices[level_number][1],
                                       d_num_nodes[level_number],
                                       d_node_offset[level_number],
                                       level_number);
        IBTK_CHKERRQ(ierr);

        num_local_nodes   [level_number] = d_local_lag_indices   [level_number].size();
        num_nonlocal_nodes[level_number] = d_nonlocal_lag_indices[level_number].size();

        std::map<std::string, Pointer<LData> >::iterator it;
        int i;
        for (it = level_data.begin(), i = 0; it != level_data.end(); ++it, ++i)
        {
            Pointer<LData> data = it->second;
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!data.isNull());
#endif
            const int depth = data->getDepth();

            // Determine the PETSc indices of the ghost nodes required in the
            // destination Vec.
            //
            // (This is only computed once for each data depth encountered,
            // including depth==1.)
            if (d_nonlocal_petsc_indices[level_number].find(depth) == d_nonlocal_petsc_indices[level_number].end())
            {
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(depth != 1);
#endif
                d_nonlocal_petsc_indices[level_number][depth].resize(num_nonlocal_nodes[level_number]);

                std::transform(d_nonlocal_petsc_indices[level_number][    1].begin(),
                               d_nonlocal_petsc_indices[level_number][    1].end(),
                               d_nonlocal_petsc_indices[level_number][depth].begin(),
                               std::bind2nd(std::multiplies<int>(),depth));
            }

            // Determine the PETSc indices of the source nodes for use when
            // scattering values from the old configuration to the new
            // configuration.
            //
            // (This is only computed once for each data depth encountered,
            // including depth==1.)
            if (src_IS[level_number].find(depth) == src_IS[level_number].end())
            {
                ierr = ISCreateStride(PETSC_COMM_WORLD,
                                      depth*num_local_nodes[level_number],
                                      depth*d_node_offset[level_number],
                                      1, &src_IS[level_number][depth]);
                IBTK_CHKERRQ(ierr);
            }

            // Determine the PETSc indices of the destination nodes for use when
            // scattering values from the old configuration to the new
            // configuration.
            //
            // (This is only computed once for each data depth encountered,
            // including depth==1.)
            if (!dst_inds_set)
            {
                dst_inds = d_local_petsc_indices[level_number];

                ierr = AOPetscToApplication(
                    d_ao[level_number], // the old AO
                    (num_local_nodes[level_number] > 0 ? num_local_nodes[level_number] : static_cast<int>(s_ao_dummy.size())),
                    (num_local_nodes[level_number] > 0 ? &dst_inds[0]                  : &s_ao_dummy[0]));
                IBTK_CHKERRQ(ierr);
                ierr = AOApplicationToPetsc(
                    new_ao[level_number],
                    (num_local_nodes[level_number] > 0 ? num_local_nodes[level_number] : static_cast<int>(s_ao_dummy.size())),
                    (num_local_nodes[level_number] > 0 ? &dst_inds[0]                  : &s_ao_dummy[0]));
                IBTK_CHKERRQ(ierr);

                dst_inds_set = true;
            }

            if (dst_IS[level_number].find(depth) == dst_IS[level_number].end())
            {
                if (depth == 1)
                {
                    ierr = ISCreateGeneral(PETSC_COMM_WORLD,
                                           num_local_nodes[level_number],
                                           num_local_nodes[level_number] > 0 ? &dst_inds[0] : NULL,
                                           &dst_IS[level_number][depth]);
                    IBTK_CHKERRQ(ierr);
                }
                else
                {
                    std::vector<int> scaled_dst_inds(dst_inds.size());
                    std::transform(dst_inds.begin(), dst_inds.end(),
                                   scaled_dst_inds.begin(),
                                   std::bind2nd(std::multiplies<int>(),depth));

                    ierr = ISCreateBlock(PETSC_COMM_WORLD,
                                         depth,
                                         num_local_nodes[level_number],
                                         num_local_nodes[level_number] > 0 ? &scaled_dst_inds[0] : NULL,
                                         &dst_IS[level_number][depth]);
                    IBTK_CHKERRQ(ierr);
                }
            }

            // Create the destination Vec and VecScatter contexts.
            src_vec[level_number][i] = data->getVec();

            if (depth == 1)
            {
                ierr = VecCreateGhost(PETSC_COMM_WORLD,
                                      num_local_nodes[level_number], PETSC_DECIDE,
                                      num_nonlocal_nodes[level_number],
                                      num_nonlocal_nodes[level_number] > 0 ? &d_nonlocal_petsc_indices[level_number][depth][0] : NULL,
                                      &dst_vec[level_number][i]);
                IBTK_CHKERRQ(ierr);
            }
            else
            {
                ierr = VecCreateGhostBlock(PETSC_COMM_WORLD, depth,
                                           depth*num_local_nodes[level_number], PETSC_DECIDE,
                                           num_nonlocal_nodes[level_number],
                                           num_nonlocal_nodes[level_number] > 0 ? &d_nonlocal_petsc_indices[level_number][depth][0] : NULL,
                                           &dst_vec[level_number][i]);
                IBTK_CHKERRQ(ierr);
            }
            ierr = VecSetBlockSize(dst_vec[level_number][i], depth);
            IBTK_CHKERRQ(ierr);

            ierr = VecScatterCreate(src_vec[level_number][i], src_IS[level_number][depth],
                                    dst_vec[level_number][i], dst_IS[level_number][depth],
                                    &scatter[level_number][i]);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Begin scattering data.
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (!d_level_contains_lag_data[level_number]) continue;

        std::map<std::string,Pointer<LData> >& level_data = d_lag_mesh_data[level_number];
        std::map<std::string,Pointer<LData> >::iterator it;
        int i;
        for (it = level_data.begin(), i = 0; it != level_data.end(); ++it, ++i)
        {
            ierr = VecScatterBegin(scatter[level_number][i], src_vec[level_number][i], dst_vec[level_number][i], INSERT_VALUES, SCATTER_FORWARD); IBTK_CHKERRQ(ierr);
        }
    }

    // Update cached indexing information on each grid patch and setup LMesh
    // data structures.
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (!d_level_contains_lag_data[level_number]) continue;

        std::set<SAMRAI::tbox::Pointer<LNode>,LNodeIndexLocalPETScIndexComp> local_nodes, ghost_nodes;
        const int num_local_nodes = getNumberOfLocalNodes(level_number);
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<LNodeSetData> lag_node_idx_data = patch->getPatchData(d_lag_node_index_current_idx);
            lag_node_idx_data->cacheLocalIndices();
            const Box<NDIM>& ghost_box = lag_node_idx_data->getGhostBox();
            for (LNodeSetData::DataIterator it = lag_node_idx_data->data_begin(ghost_box);
                 it != lag_node_idx_data->data_end(); ++it)
            {
                LNodeSet::value_type& node_idx = *it;
                if (node_idx->getLocalPETScIndex() < num_local_nodes)
                {
                    local_nodes.insert(node_idx);
                }
                else
                {
                    ghost_nodes.insert(node_idx);
                }
            }
        }
        std::ostringstream name_stream;
        name_stream << d_object_name << "::mesh::level_" << level_number;
        d_lag_mesh[level_number] = new LMesh(name_stream.str());
        static const int sorted = true;
        d_lag_mesh[level_number]->setNodes(std::vector<LNode*>(local_nodes.begin(), local_nodes.end()), sorted);
        d_lag_mesh[level_number]->setGhostNodes(std::vector<LNode*>(ghost_nodes.begin(), ghost_nodes.end()), sorted);
    }

    // End scattering data, reset data, and destroy the VecScatter contexts.
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        if (!d_level_contains_lag_data[level_number]) continue;

        std::map<std::string,Pointer<LData> >& level_data = d_lag_mesh_data[level_number];
        std::map<std::string,Pointer<LData> >::iterator it;
        int i;
        for (it = level_data.begin(), i = 0; it != level_data.end(); ++it, ++i)
        {
            ierr = VecScatterEnd(scatter[level_number][i], src_vec[level_number][i], dst_vec[level_number][i], INSERT_VALUES, SCATTER_FORWARD); IBTK_CHKERRQ(ierr);
            ierr = VecScatterDestroy(scatter[level_number][i]); IBTK_CHKERRQ(ierr);
            Pointer<LData> data = it->second;
            const int depth = data->getDepth();
            data->resetData(dst_vec[level_number][i], d_nonlocal_petsc_indices[level_number][depth]);
        }
    }

    // Distribute nonlocal data to the new configuration.
    beginNonlocalDataFill(coarsest_ln,finest_ln);
    endNonlocalDataFill(  coarsest_ln,finest_ln);

    // Indicate that the levels have been synchronized and destroy unneeded
    // ordering and indexing objects.
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        d_needs_synch[level_number] = false;

        if (d_ao[level_number])
        {
            ierr = AODestroy(d_ao[level_number]);
            IBTK_CHKERRQ(ierr);
        }
        d_ao[level_number] = new_ao[level_number];

        for (std::map<int,IS>::iterator it = src_IS[level_number].begin();
             it != src_IS[level_number].end(); ++it)
        {
            ierr = ISDestroy(it->second);
            IBTK_CHKERRQ(ierr);
        }

        for (std::map<int,IS>::iterator it = dst_IS[level_number].begin();
             it != dst_IS[level_number].end(); ++it)
        {
            ierr = ISDestroy(it->second);
            IBTK_CHKERRQ(ierr);
        }
    }

    // If a Silo data writer is registered with the manager, give it access to
    // the new application orderings.
    if (!d_silo_writer.isNull())
    {
        d_silo_writer->registerLagrangianAO(d_ao, coarsest_ln, finest_ln);
    }
#if (NDIM == 3)
    // If a myocardial3D data writer is registered with the manager, give it
    // access to the new application orderings.
    if (!d_m3D_writer.isNull())
    {
        d_m3D_writer->registerLagrangianAO(d_ao, coarsest_ln, finest_ln);
    }
#endif
    IBTK_TIMER_STOP(t_end_data_redistribution);
    return;
}// endDataRedistribution

void
LDataManager::updateWorkloadData(
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    IBTK_TIMER_START(t_update_workload_data);

    const int coarsest_ln = (coarsest_ln_in == -1) ? d_coarsest_ln : coarsest_ln_in;
    const int finest_ln = (finest_ln_in == -1) ? d_finest_ln : finest_ln_in;

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(coarsest_ln >= d_coarsest_ln && coarsest_ln <= d_finest_ln);
    TBOX_ASSERT(finest_ln   >= d_coarsest_ln && finest_ln   <= d_finest_ln);
#endif

    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
            Pointer<CellData<NDIM,double> > workload_data = patch->getPatchData(d_workload_idx);
            Pointer<CellData<NDIM,double> > node_count_data = patch->getPatchData(d_node_count_idx);
            workload_data->fillAll(d_alpha_work);
            node_count_data->fillAll(0.0);
            for (LNodeSetData::SetIterator it(*idx_data); it; it++)
            {
                const Index<NDIM>& i = it.getIndex();
                if (patch_box.contains(i))
                {
                    const LNodeSet& node_set = *it;
                    (*node_count_data)(i) = node_set.size();
                    (*workload_data)(i) += d_beta_work*(*node_count_data)(i);
                }
            }
        }
    }

    IBTK_TIMER_STOP(t_update_workload_data);
    return;
}// updateWorkloadData

///
///  The following routines:
///
///      initializeLevelData(),
///      resetHierarchyConfiguration(),
///      applyGradientDetector()
///
///  are concrete implementations of functions declared in the
///  StandardTagAndInitStrategy<NDIM> abstract base class.
///

void
LDataManager::initializeLevelData(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const Pointer<BasePatchLevel<NDIM> > old_level,
    const bool allocate_data)
{
    IBTK_TIMER_START(t_initialize_level_data);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

#ifdef DEBUG_CHECK_ASSERTIONS
    // Check for overlapping boxes on this level.
    //
    // (This is potentially fairly expensive and hence is only done when
    // assertion checking is active.)
    BoxList<NDIM> boxes(level->getBoxes());

    std::vector<bool> patch_overlaps(boxes.getNumberOfItems());
    std::vector<bool>::size_type j, k;

    for (k = 0; k < patch_overlaps.size(); ++k)
    {
        patch_overlaps[k] = false;
    }
    k = 0;
    while (!boxes.isEmpty())
    {
        j = k+1;
        Box<NDIM> tryme = boxes.getFirstItem();
        boxes.removeFirstItem();

        for (BoxList<NDIM>::Iterator ib(boxes); ib; ib++)
        {
            if (tryme.intersects(ib()))
            {
                patch_overlaps[k] = true;
                patch_overlaps[j] = true;
            }
            ++j;
        }
        ++k;
    }

    for (k = 0; k < patch_overlaps.size(); ++k)
    {
        if (patch_overlaps[k])
        {
            TBOX_ERROR(d_object_name << "::initializeLevelData()\n"
                       << "  patch " << k << " overlaps another patch!\n");
        }
    }
#endif

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
    if (!old_level.isNull() && d_level_contains_lag_data[level_number])
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
        d_level_contains_lag_data       .resize(level_number+1);
        d_strct_name_to_strct_id_map    .resize(level_number+1);
        d_strct_id_to_strct_name_map    .resize(level_number+1);
        d_strct_id_to_lag_idx_range_map .resize(level_number+1);
        d_last_lag_idx_to_strct_id_map  .resize(level_number+1);
        d_inactive_strcts               .resize(level_number+1);
        d_displaced_strct_ids           .resize(d_finest_ln+1);
        d_displaced_strct_bounding_boxes.resize(d_finest_ln+1);
        d_displaced_strct_lnode_idxs    .resize(d_finest_ln+1);
        d_displaced_strct_lnode_posns   .resize(d_finest_ln+1);
        d_lag_mesh                      .resize(level_number+1);
        d_lag_mesh_data                 .resize(level_number+1);
        d_needs_synch                   .resize(level_number+1,false);
        d_ao                            .resize(level_number+1);
        d_num_nodes                     .resize(level_number+1);
        d_node_offset                   .resize(level_number+1);
        d_local_lag_indices             .resize(level_number+1);
        d_nonlocal_lag_indices          .resize(level_number+1);
        d_local_petsc_indices           .resize(level_number+1);
        d_nonlocal_petsc_indices        .resize(level_number+1);

#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!d_lag_init.isNull());
#endif
        d_level_contains_lag_data[level_number] = d_lag_init->getLevelHasLagrangianData(level_number, can_be_refined);
    }
    if (initial_time && d_level_contains_lag_data[level_number])
    {
        int ierr;

        // 1. Determine the number of local (on processor) nodes to be allocated
        //    on the patch level and allocate space for the local and non-local
        //    index data.
        const unsigned int num_local_nodes = d_lag_init->computeLocalNodeCountOnPatchLevel(
            hierarchy, level_number, init_data_time, can_be_refined, initial_time);

        d_local_lag_indices  [level_number].resize(num_local_nodes,-1);
        d_local_petsc_indices[level_number].resize(num_local_nodes,-1);

        d_nonlocal_lag_indices  [level_number]   .clear();
        d_nonlocal_petsc_indices[level_number]   .clear();
        d_nonlocal_petsc_indices[level_number][1].clear();

        computeNodeOffsets(d_num_nodes[level_number], d_node_offset[level_number], num_local_nodes);

        // 2. Allocate LData corresponding to the curvilinear mesh node
        //    positions and velocities.
        static const bool maintain_data = true;
        createLData(     POSN_DATA_NAME, level_number, NDIM, maintain_data);
        createLData(INIT_POSN_DATA_NAME, level_number, NDIM, maintain_data);
        createLData(      VEL_DATA_NAME, level_number, NDIM, maintain_data);

        // 3. Initialize the Lagrangian data.
        d_lag_init->initializeStructureIndexingOnPatchLevel(
            d_strct_id_to_strct_name_map   [level_number],
            d_strct_id_to_lag_idx_range_map[level_number],
            level_number,
            init_data_time, can_be_refined, initial_time, this);

        for (std::map<int,std::string>::const_iterator cit(d_strct_id_to_strct_name_map[level_number].begin());
             cit != d_strct_id_to_strct_name_map[level_number].end(); ++cit)
        {
            d_strct_name_to_strct_id_map[level_number][cit->second] = cit->first;
        }

        for (std::map<int,std::pair<int,int> >::const_iterator cit(d_strct_id_to_lag_idx_range_map[level_number].begin());
             cit != d_strct_id_to_lag_idx_range_map[level_number].end(); ++cit)
        {
            d_last_lag_idx_to_strct_id_map[level_number][cit->second.second-1] = cit->first;
        }

        // WARNING: If either of the following offsets is ever nonzero, note
        // that it may be necessary to modify IBHierarchyIntegrator, in
        // particular the code where the data related to the implementation of
        // the penalty IB method are initialized.
        static const unsigned int global_index_offset = 0;
        static const unsigned int local_index_offset = 0;
        const unsigned int num_initialized_local_nodes = d_lag_init->initializeDataOnPatchLevel(
            d_lag_node_index_current_idx,
            global_index_offset, local_index_offset,
            d_lag_mesh_data[level_number][POSN_DATA_NAME],
            d_lag_mesh_data[level_number][ VEL_DATA_NAME],
            hierarchy, level_number,
            init_data_time, can_be_refined, initial_time, this);

        ierr = VecCopy(d_lag_mesh_data[level_number][     POSN_DATA_NAME]->getVec(),
                       d_lag_mesh_data[level_number][INIT_POSN_DATA_NAME]->getVec());
        IBTK_CHKERRQ(ierr);

        if (num_local_nodes != num_initialized_local_nodes)
        {
            TBOX_ERROR("LDataManager::initializeLevelData()"                             << "\n" <<
                       "  num_local_nodes             = " << num_local_nodes             << "\n" <<
                       "  num_initialized_local_nodes = " << num_initialized_local_nodes << "\n");
        }

        // 4. Compute the initial distribution (indexing) data.
        std::set<SAMRAI::tbox::Pointer<LNode>,LNodeIndexLocalPETScIndexComp> local_nodes, ghost_nodes;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
            Pointer<CellData<NDIM,double> > node_count_data = patch->getPatchData(d_node_count_idx);
            Pointer<CellData<NDIM,double> >   workload_data = patch->getPatchData(  d_workload_idx);

            node_count_data->fillAll(0.0);
            workload_data  ->fillAll(0.0);

            idx_data->cacheLocalIndices();
            for (LNodeSetData::SetIterator it(*idx_data); it; it++)
            {
                const CellIndex<NDIM>& i = it.getIndex();
                LNodeSet& node_set = *it;
                const bool patch_owns_idx_set = patch_box.contains(i);
                if (patch_owns_idx_set)
                {
                    (*node_count_data)(i) = node_set.size();
                    (*  workload_data)(i) = d_alpha_work + d_beta_work*(*node_count_data)(i);
                }

                for (LNodeSet::iterator n = node_set.begin(); n != node_set.end(); ++n)
                {
                    LNodeSet::value_type& node_idx = *n;
                    const int lag_idx   = node_idx->getLagrangianIndex();
                    const int local_idx = node_idx->getLocalPETScIndex();
                    if (!(0 <= local_idx && local_idx < static_cast<int>(num_local_nodes)))
                    {
                        TBOX_ERROR("LDataManager::initializeLevelData()"     << "\n" <<
                                   "  local_idx       = " << local_idx       << "\n" <<
                                   "  num_local_nodes = " << num_local_nodes << "\n");
                    }
                    d_local_lag_indices  [level_number][local_idx] = lag_idx;
                    d_local_petsc_indices[level_number][local_idx] = local_idx + d_node_offset[level_number];
                    local_nodes.insert(node_idx);
                }
            }
        }
        std::ostringstream name_stream;
        name_stream << d_object_name << "::mesh::level_" << level_number;
        d_lag_mesh[level_number] = new LMesh(name_stream.str());
        static const int sorted = true;
        d_lag_mesh[level_number]->setNodes(std::vector<LNode*>(local_nodes.begin(), local_nodes.end()),sorted);
        d_lag_mesh[level_number]->setGhostNodes(std::vector<LNode*>(ghost_nodes.begin(), ghost_nodes.end()),sorted);

        // 5. The AO (application order) is determined by the initial values of
        //    the local Lagrangian indices.
        if (d_ao[level_number])
        {
            ierr = AODestroy(d_ao[level_number]);
            IBTK_CHKERRQ(ierr);
        }

        ierr = AOCreateBasic(PETSC_COMM_WORLD,
                             num_local_nodes,
                             num_local_nodes > 0 ? &d_local_lag_indices  [level_number][0] : NULL,
                             num_local_nodes > 0 ? &d_local_petsc_indices[level_number][0] : NULL,
                             &d_ao[level_number]);
        IBTK_CHKERRQ(ierr);
    }

    // Initialize workload data and setup the load balancer.
    if (!d_load_balancer.isNull())
    {
        HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(hierarchy,level_number,level_number);
        hier_cc_data_ops.setToScalar(d_workload_idx, d_alpha_work);
        if (!old_level.isNull() && d_level_contains_lag_data[level_number])
        {
            Pointer<RefineOperator<NDIM> > regrid_fill_op = Pointer<RefineOperator<NDIM> >(NULL);
            Pointer<RefineAlgorithm<NDIM> > regrid_fill_alg = new RefineAlgorithm<NDIM>();
            regrid_fill_alg->registerRefine(d_workload_idx, d_workload_idx, d_workload_idx, regrid_fill_op);
            regrid_fill_alg->createSchedule(level, old_level)->fillData(init_data_time);
        }
        if (d_level_contains_lag_data[level_number])
        {
            d_load_balancer->setWorkloadPatchDataIndex(d_workload_idx, level_number);
        }
        else
        {
            d_load_balancer->setUniformWorkload(level_number);
        }

        // Ensure that workload and related data are allocated.
        for (int ln = 0; ln <= level_number; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            if (!level->checkAllocated(d_workload_idx)) level->allocatePatchData(d_workload_idx, init_data_time);
            if (!level->checkAllocated(d_node_count_idx)) level->allocatePatchData(d_node_count_idx, init_data_time);
        }
    }

    // If a Silo data writer is registered with the manager, give it access to
    // the new application ordering.
    if (!d_silo_writer.isNull() && d_level_contains_lag_data[level_number])
    {
        d_silo_writer->registerCoordsData(d_lag_mesh_data[level_number][POSN_DATA_NAME], level_number);
        d_silo_writer->registerLagrangianAO(d_ao[level_number], level_number);
    }
#if (NDIM == 3)
    // If a myocardial3D data writer is registered with the manager, give it
    // access to the new application ordering.
    if (!d_m3D_writer.isNull() && d_level_contains_lag_data[level_number])
    {
        d_m3D_writer->registerCoordsData(d_lag_mesh_data[level_number][POSN_DATA_NAME], level_number);
        d_m3D_writer->registerLagrangianAO(d_ao[level_number], level_number);
    }
#endif
    IBTK_TIMER_STOP(t_initialize_level_data);
    return;
}// initializeLevelData

void
LDataManager::resetHierarchyConfiguration(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_ln,
    const int finest_ln)
{
    IBTK_TIMER_START(t_reset_hierarchy_configuration);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_ln >= 0)
                && (coarsest_ln <= finest_ln)
                && (finest_ln <= hierarchy->getFinestLevelNumber()));
    for (int level_number = 0; level_number <= finest_ln; ++level_number)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
    }
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // Reset the patch hierarchy and levels.
    setPatchHierarchy(hierarchy);
    resetLevels(0,finest_hier_level);

    // Reset the Silo data writer.
    if (!d_silo_writer.isNull())
    {
        d_silo_writer->setPatchHierarchy(hierarchy);
        d_silo_writer->resetLevels(d_coarsest_ln, d_finest_ln);
        for (int level_number = d_coarsest_ln; level_number <= d_finest_ln; ++level_number)
        {
            if (!d_level_contains_lag_data[level_number]) continue;
            d_silo_writer->registerCoordsData(d_lag_mesh_data[level_number][POSN_DATA_NAME], level_number);
        }
    }
#if (NDIM == 3)
    // Reset the myocardial3D data writer.
    if (!d_m3D_writer.isNull())
    {
        d_m3D_writer->setPatchHierarchy(hierarchy);
        d_m3D_writer->resetLevels(d_coarsest_ln, d_finest_ln);
        for (int level_number = d_coarsest_ln; level_number <= d_finest_ln; ++level_number)
        {
            if (!d_level_contains_lag_data[level_number]) continue;
            d_m3D_writer->registerCoordsData(d_lag_mesh_data[level_number][POSN_DATA_NAME], level_number);
        }
    }
#endif
    // If we have added or removed a level, resize the schedule vectors.
    d_lag_node_index_bdry_fill_scheds.resize(finest_hier_level+1);
    d_node_count_coarsen_scheds      .resize(finest_hier_level+1);

    // (Re)build refine communication schedules.  These are created for only the
    // specified levels in the hierarchy.
    //
    // NOTE: These schedules do not fill from coarser levels in the hierarchy.
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
        d_lag_node_index_bdry_fill_scheds[level_number] = d_lag_node_index_bdry_fill_alg->createSchedule(level);
    }

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (int level_number = std::max(coarsest_ln,1); level_number <= finest_hier_level; ++level_number)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
        Pointer<PatchLevel<NDIM> > coarser_level = hierarchy->getPatchLevel(level_number-1);
        d_node_count_coarsen_scheds[level_number] = d_node_count_coarsen_alg->createSchedule(coarser_level, level);
    }

    IBTK_TIMER_STOP(t_reset_hierarchy_configuration);
    return;
}// resetHierarchyConfiguration

void
LDataManager::applyGradientDetector(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    IBTK_TIMER_START(t_apply_gradient_detector);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    if (initial_time)
    {
        // Tag cells for refinement based on the initial configuration of the
        // Lagrangian structure.
        d_lag_init->tagCellsForInitialRefinement(hierarchy, level_number, error_data_time, tag_index);
    }
    else if (hierarchy->finerLevelExists(level_number))
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
        Pointer<PatchLevel<NDIM> > finer_level = hierarchy->getPatchLevel(level_number+1);

        // Zero out the node count data on the current level.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM,double> > node_count_data = patch->getPatchData(d_node_count_idx);
            node_count_data->fillAll(0.0);
        }

        // Compute the workload and node count data on the next finer level of
        // the patch hierarchy.
        updateWorkloadData(level_number+1);

        // Coarsen the node count data from the next finer level of the patch
        // hierarchy.
        d_node_count_coarsen_scheds[level_number+1]->coarsenData();

        // Tag cells for refinement wherever there exist nodes on the next finer
        // level of the Cartesian grid.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM,int> > tag_data = patch->getPatchData(tag_index);
            const Pointer<CellData<NDIM,double> > node_count_data = patch->getPatchData(d_node_count_idx);

            for (CellIterator<NDIM> ic(patch_box); ic; ic++)
            {
                const CellIndex<NDIM>& i = ic();
                if (!MathUtilities<double>::equalEps((*node_count_data)(i),0.0))
                {
                    (*tag_data)(i) = 1;
                }
            }
        }

        // Tag cells for refinement within the bounding boxes of any displaced
        // structures on finer levels of the patch hierarchy.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const CellIndex<NDIM>& patch_lower = patch_box.lower();
            const CellIndex<NDIM>& patch_upper = patch_box.upper();
            const double* const patchXLower = patch_geom->getXLower();
            const double* const patchXUpper = patch_geom->getXUpper();
            const double* const patchDx = patch_geom->getDx();

            Pointer<CellData<NDIM,int> > tag_data = patch->getPatchData(tag_index);

            for (int ln = level_number+1; ln <= d_finest_ln; ++ln)
            {
                for (std::vector<std::pair<blitz::TinyVector<double,NDIM>,blitz::TinyVector<double,NDIM> > >::const_iterator
                         cit = d_displaced_strct_bounding_boxes[ln].begin();
                     cit != d_displaced_strct_bounding_boxes[ln].end(); ++cit)
                {
                    const std::pair<blitz::TinyVector<double,NDIM>,blitz::TinyVector<double,NDIM> >& bounding_box = *cit;

                    // Determine the region of index space covered by the
                    // displaced structure bounding box.
                    const CellIndex<NDIM> bbox_lower = IndexUtilities::getCellIndex(bounding_box.first , patchXLower, patchXUpper, patchDx, patch_lower, patch_upper);
                    const CellIndex<NDIM> bbox_upper = IndexUtilities::getCellIndex(bounding_box.second, patchXLower, patchXUpper, patchDx, patch_lower, patch_upper);
                    const Box<NDIM> tag_box(bbox_lower,bbox_upper);
                    tag_data->fillAll(1,tag_box);
                }
            }
        }

        // Compute the workload and node count data on the present level of the
        // patch hierarchy (since it was invalidated above).
        updateWorkloadData(level_number);
    }

    IBTK_TIMER_STOP(t_apply_gradient_detector);
    return;
}// applyGradientDetector

void
LDataManager::putToDatabase(
    Pointer<Database> db)
{
    IBTK_TIMER_START(t_put_to_database);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("LDATA_MANAGER_VERSION", LDATA_MANAGER_VERSION);

    db->putInteger("d_coarsest_ln", d_coarsest_ln);
    db->putInteger("d_finest_ln"  , d_finest_ln  );
    db->putDouble ("d_alpha_work" , d_alpha_work );
    db->putDouble ("d_beta_work"  , d_beta_work  );

    // Write out data that is stored on a level-by-level basis.
    for (int level_number = d_coarsest_ln; level_number <= d_finest_ln; ++level_number)
    {
        std::ostringstream stream;
        stream << "level_" << level_number;
        const std::string level_db_name = stream.str();
        Pointer<Database> level_db = db->putDatabase(level_db_name);

        level_db->putBool("d_level_contains_lag_data", d_level_contains_lag_data[level_number]);

        if (!d_level_contains_lag_data[level_number]) continue;

        std::vector<int> lstruct_ids, lstruct_lag_idx_range_first, lstruct_lag_idx_range_second, lstruct_activation;
        std::vector<std::string> lstruct_names;
        for (std::map<int,std::string>::iterator it = d_strct_id_to_strct_name_map[level_number].begin();
             it != d_strct_id_to_strct_name_map[level_number].end(); ++it)
        {
            const int id = it->first;
            lstruct_ids.push_back(id);
            lstruct_names.push_back(it->second);
            lstruct_lag_idx_range_first.push_back(d_strct_id_to_lag_idx_range_map[level_number].find(id)->second.first);
            lstruct_lag_idx_range_second.push_back(d_strct_id_to_lag_idx_range_map[level_number].find(id)->second.second);
            lstruct_activation.push_back(static_cast<int>(d_inactive_strcts[level_number].getSet().find(id) != d_inactive_strcts[level_number].getSet().end()));
        }
        level_db->putInteger("n_lstructs", lstruct_ids.size());
        if (!lstruct_ids.empty())
        {
            level_db->putIntegerArray("lstruct_ids", &lstruct_ids[0], lstruct_ids.size());
            level_db->putIntegerArray("lstruct_lag_idx_range_first", &lstruct_lag_idx_range_first[0], lstruct_lag_idx_range_first.size());
            level_db->putIntegerArray("lstruct_lag_idx_range_second", &lstruct_lag_idx_range_second[0], lstruct_lag_idx_range_second.size());
            level_db->putIntegerArray("lstruct_activation", &lstruct_activation[0], lstruct_activation.size());
            level_db->putStringArray("lstruct_names", &lstruct_names[0], lstruct_names.size());
        }

        std::vector<std::string> ldata_names;
        for (std::map<std::string,Pointer<LData> >::iterator it = d_lag_mesh_data[level_number].begin();
             it != d_lag_mesh_data[level_number].end(); ++it)
        {
            ldata_names.push_back(it->first);
            it->second->putToDatabase(level_db->putDatabase(ldata_names.back()));
        }
        level_db->putInteger("n_ldata_names", ldata_names.size());
        if (!ldata_names.empty())
        {
            level_db->putStringArray("ldata_names",
                                     &ldata_names[0],
                                     ldata_names.size());
        }

        level_db->putInteger("d_num_nodes"  , d_num_nodes  [level_number]);
        level_db->putInteger("d_node_offset", d_node_offset[level_number]);

        level_db->putInteger("n_local_lag_indices",
                             d_local_lag_indices[level_number].size());
        if (!d_local_lag_indices[level_number].empty())
        {
            level_db->putIntegerArray("d_local_lag_indices",
                                      &d_local_lag_indices[level_number][0],
                                      d_local_lag_indices [level_number].size());
        }
        level_db->putInteger("n_nonlocal_lag_indices",
                             d_nonlocal_lag_indices[level_number].size());
        if (!d_nonlocal_lag_indices[level_number].empty())
        {
            level_db->putIntegerArray("d_nonlocal_lag_indices",
                                      &d_nonlocal_lag_indices[level_number][0],
                                      d_nonlocal_lag_indices [level_number].size());
        }
        level_db->putInteger("n_local_petsc_indices",
                             d_local_petsc_indices[level_number].size());
        if (!d_local_petsc_indices[level_number].empty())
        {
            level_db->putIntegerArray("d_local_petsc_indices",
                                      &d_local_petsc_indices[level_number][0],
                                      d_local_petsc_indices [level_number].size());
        }
        // NOTE: d_nonlocal_petsc_indices[level_number] is a map from the data
        // depth to the nonlocal petsc indices for that particular depth.  We
        // only serialize the indices corresponding to a data depth of 1.
        level_db->putInteger("n_nonlocal_petsc_indices",
                             d_nonlocal_petsc_indices[level_number][1].size());
        if (!d_nonlocal_petsc_indices[level_number][1].empty())
        {
            level_db->putIntegerArray("d_nonlocal_petsc_indices",
                                      &d_nonlocal_petsc_indices[level_number][1][0],
                                      d_nonlocal_petsc_indices [level_number][1].size());
        }
    }

    IBTK_TIMER_STOP(t_put_to_database);
    return;
}// putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

LDataManager::LDataManager(
    const std::string& object_name,
    const std::string& interp_weighting_fcn,
    const std::string& spread_weighting_fcn,
    const IntVector<NDIM>& ghost_width,
    bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_hierarchy(NULL),
      d_grid_geom(NULL),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_visit_writer(NULL),
      d_silo_writer(NULL),
#if (NDIM == 3)
      d_m3D_writer(NULL),
#endif
      d_load_balancer(NULL),
      d_lag_init(NULL),
      d_level_contains_lag_data(),
      d_lag_node_index_var(NULL),
      d_lag_node_index_current_idx(-1),
      d_lag_node_index_scratch_idx(-1),
      d_alpha_work(1.0),
      d_beta_work(1.0),
      d_workload_var(NULL),
      d_workload_idx(-1),
      d_output_workload(false),
      d_node_count_var(NULL),
      d_node_count_idx(-1),
      d_output_node_count(false),
      d_interp_weighting_fcn(interp_weighting_fcn),
      d_spread_weighting_fcn(spread_weighting_fcn),
      d_ghost_width(ghost_width),
      d_lag_node_index_bdry_fill_alg(NULL),
      d_lag_node_index_bdry_fill_scheds(),
      d_node_count_coarsen_alg(NULL),
      d_node_count_coarsen_scheds(),
      d_current_context(NULL),
      d_scratch_context(NULL),
      d_current_data(),
      d_scratch_data(),
      d_lag_mesh(),
      d_lag_mesh_data(),
      d_needs_synch(true),
      d_ao(),
      d_num_nodes(),
      d_node_offset(),
      d_local_lag_indices(),
      d_nonlocal_lag_indices(),
      d_local_petsc_indices(),
      d_nonlocal_petsc_indices()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(ghost_width.min() >= 0);
#endif

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
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    d_current_context = var_db->getContext(d_object_name+"::CURRENT");
    d_scratch_context = var_db->getContext(d_object_name+"::SCRATCH");

    // Register the SAMRAI variables with the VariableDatabase.
    d_lag_node_index_var = new LNodeSetVariable(d_object_name+"::lag_node_index");

    // Setup the current context.
    d_lag_node_index_current_idx = var_db->registerVariableAndContext(d_lag_node_index_var, d_current_context, d_ghost_width);
    d_current_data.setFlag(d_lag_node_index_current_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_lag_node_index_current_idx);
    }

    // Setup the scratch context.
    d_lag_node_index_scratch_idx = var_db->registerVariableAndContext(d_lag_node_index_var, d_scratch_context, d_ghost_width);
    d_scratch_data.setFlag(d_lag_node_index_scratch_idx);

    // Setup a refine algorithm, used to fill LNode boundary data.
    Pointer<RefineOperator<NDIM> > lag_node_index_bdry_fill_op = Pointer<RefineOperator<NDIM> >(NULL);
    d_lag_node_index_bdry_fill_alg = new RefineAlgorithm<NDIM>();
    d_lag_node_index_bdry_fill_alg->registerRefine(d_lag_node_index_current_idx, d_lag_node_index_current_idx, d_lag_node_index_scratch_idx, lag_node_index_bdry_fill_op);

    // Register the workload variable with the VariableDatabase.
    d_workload_var = new CellVariable<NDIM,double>(d_object_name+"::workload");
    d_workload_idx = var_db->registerVariableAndContext(d_workload_var, d_current_context, 0);
    d_current_data.setFlag(d_workload_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_workload_idx);
    }

    // Register the node count variable with the VariableDatabase.
    d_node_count_var = new CellVariable<NDIM,double>(d_object_name+"::node_count");
    d_node_count_idx = var_db->registerVariableAndContext(d_node_count_var, d_current_context, 0);
    d_current_data.setFlag(d_node_count_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_node_count_idx);
    }

    Pointer<CoarsenOperator<NDIM> > node_count_coarsen_op = new CartesianCellDoubleWeightedAverage<NDIM>();
    d_node_count_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    d_node_count_coarsen_alg->registerCoarsen(d_node_count_idx, d_node_count_idx, node_count_coarsen_op);

    // Setup Timers.
    IBTK_DO_ONCE(
        t_spread = TimerManager::getManager()->getTimer("IBTK::LDataManager::spread()");
        t_interp = TimerManager::getManager()->getTimer("IBTK::LDataManager::interp()");
        t_map_lagrangian_to_petsc = TimerManager::getManager()->getTimer("IBTK::LDataManager::mapLagrangianToPETSc()");
        t_map_petsc_to_lagrangian = TimerManager::getManager()->getTimer("IBTK::LDataManager::mapPETScToLagrangian()");
        t_begin_data_redistribution = TimerManager::getManager()->getTimer("IBTK::LDataManager::beginDataRedistribution()");
        t_end_data_redistribution = TimerManager::getManager()->getTimer("IBTK::LDataManager::endDataRedistribution()");
        t_update_workload_data = TimerManager::getManager()->getTimer("IBTK::LDataManager::updateWorkloadData()");
        t_initialize_level_data = TimerManager::getManager()->getTimer("IBTK::LDataManager::initializeLevelData()");
        t_reset_hierarchy_configuration = TimerManager::getManager()->getTimer("IBTK::LDataManager::resetHierarchyConfiguration()");
        t_apply_gradient_detector = TimerManager::getManager()->getTimer("IBTK::LDataManager::applyGradientDetector()");
        t_put_to_database = TimerManager::getManager()->getTimer("IBTK::LDataManager::putToDatabase()");
        t_begin_nonlocal_data_fill = TimerManager::getManager()->getTimer("IBTK::LDataManager::beginNonlocalDataFill()");
        t_end_nonlocal_data_fill = TimerManager::getManager()->getTimer("IBTK::LDataManager::endNonlocalDataFill()");
        t_compute_node_distribution = TimerManager::getManager()->getTimer("IBTK::LDataManager::computeNodeDistribution()");
        t_compute_node_offsets = TimerManager::getManager()->getTimer("IBTK::LDataManager::computeNodeOffsets()");
                 );
    return;
}// LDataManager

LDataManager::~LDataManager()
{
    // Destroy any remaining AO objects.
    int ierr;
    for (int level_number = d_coarsest_ln; level_number <= d_finest_ln; ++level_number)
    {
        if (d_ao[level_number])
        {
            ierr = AODestroy(d_ao[level_number]);
            IBTK_CHKERRQ(ierr);
        }
    }
    return;
}// ~LDataManager

/////////////////////////////// PRIVATE //////////////////////////////////////

void
LDataManager::spread_specialized(
    const int f_data_idx,
    std::vector<Pointer<LData> >& F_data,
    std::vector<Pointer<LData> >& X_data,
    const bool F_data_ghost_node_update,
    const bool X_data_ghost_node_update,
    const int coarsest_ln,
    const int finest_ln)
{
    // Zero inactivated components.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        zeroInactivatedComponents(F_data[ln], ln);
    }

    // Spread data from the Lagrangian mesh to the Eulerian grid.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (levelContainsLagrangianData(ln))
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(ln == finest_ln);
#endif
            if (F_data_ghost_node_update) F_data[ln]->beginGhostUpdate();
            if (X_data_ghost_node_update) X_data[ln]->beginGhostUpdate();
            if (F_data_ghost_node_update) F_data[ln]->endGhostUpdate();
            if (X_data_ghost_node_update) X_data[ln]->endGhostUpdate();
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_data_idx);
                Pointer<CellData<NDIM,double> > f_cc_data = f_data;
                Pointer<SideData<NDIM,double> > f_sc_data = f_data;
                const bool is_cc_data = !f_cc_data.isNull();
                const bool is_sc_data = !f_sc_data.isNull();
                const Pointer<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
                const Box<NDIM>& box = idx_data->getGhostBox();
                if (is_cc_data) LEInteractor::spread(f_cc_data, F_data[ln], X_data[ln], idx_data, patch, box, periodic_shift, d_spread_weighting_fcn);
                if (is_sc_data) LEInteractor::spread(f_sc_data, F_data[ln], X_data[ln], idx_data, patch, box, periodic_shift, d_spread_weighting_fcn);
            }
        }
    }
    return;
}// spread_specialized

void
LDataManager::interp_specialized(
    const int f_data_idx,
    std::vector<Pointer<LData> >& F_data,
    std::vector<Pointer<LData> >& X_data,
    std::vector<Pointer<RefineSchedule<NDIM> > > f_ghost_fill_scheds,
    const double fill_data_time,
    const int coarsest_ln,
    const int finest_ln)
{
    // Interpolate data from the Eulerian grid to the Lagrangian mesh.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
        if (levelContainsLagrangianData(ln))
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(ln == finest_ln);
#endif
            if (ln < static_cast<int>(f_ghost_fill_scheds.size()) && !f_ghost_fill_scheds[ln].isNull()) f_ghost_fill_scheds[ln]->fillData(fill_data_time);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_data_idx);
                Pointer<CellData<NDIM,double> > f_cc_data = f_data;
                Pointer<SideData<NDIM,double> > f_sc_data = f_data;
                const bool is_cc_data = !f_cc_data.isNull();
                const bool is_sc_data = !f_sc_data.isNull();
                const Pointer<LNodeSetData> idx_data = patch->getPatchData(d_lag_node_index_current_idx);
                const Box<NDIM>& box = idx_data->getBox();
                if (is_cc_data) LEInteractor::interpolate(F_data[ln], X_data[ln], idx_data, f_cc_data, patch, box, periodic_shift, d_interp_weighting_fcn);
                if (is_sc_data) LEInteractor::interpolate(F_data[ln], X_data[ln], idx_data, f_sc_data, patch, box, periodic_shift, d_interp_weighting_fcn);
            }
        }
    }

    // Zero inactivated components.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        zeroInactivatedComponents(F_data[ln], ln);
    }
    return;
}// interp_specialized

void
LDataManager::scatterData(
    Vec& lagrangian_vec,
    Vec& petsc_vec,
    const int level_number,
    ScatterMode mode) const
{
    int ierr;

    // Get the vector sizes.
    int petsc_size, lagrangian_size;
    ierr = VecGetSize(petsc_vec, &petsc_size);  IBTK_CHKERRQ(ierr);
    ierr = VecGetSize(lagrangian_vec, &lagrangian_size);  IBTK_CHKERRQ(ierr);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(petsc_size == lagrangian_size);
#endif
    int petsc_bs, lagrangian_bs;
    ierr = VecGetBlockSize(petsc_vec, &petsc_bs);  IBTK_CHKERRQ(ierr);
    ierr = VecGetBlockSize(lagrangian_vec, &lagrangian_bs);  IBTK_CHKERRQ(ierr);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(petsc_bs == lagrangian_bs);
#endif
    const int depth = petsc_bs;

    // Determine the application indices corresponding to the local PETSc
    // indices.
    int local_sz;
    ierr = VecGetLocalSize(lagrangian_vec, &local_sz);  IBTK_CHKERRQ(ierr);
    local_sz /= depth;
    std::vector<int> local_lag_idxs(local_sz,-1);

    int ilo, ihi;
    ierr = VecGetOwnershipRange(lagrangian_vec, &ilo, &ihi);  IBTK_CHKERRQ(ierr);
    ilo /= depth;
    ihi /= depth;
    for (int k = 0; k < local_sz; ++k)
    {
        local_lag_idxs[k] = ilo+k;
    }
    mapLagrangianToPETSc(local_lag_idxs, level_number);
    for (int k = 0; k < local_sz; ++k)
    {
        local_lag_idxs[k] *= depth;
    }

    IS lag_is;
    ierr = ISCreateBlock(PETSC_COMM_WORLD, depth,
                         local_lag_idxs.size(),
                         local_lag_idxs.empty() ? NULL : &local_lag_idxs[0],
                         &lag_is);  IBTK_CHKERRQ(ierr);

    // Create a VecScatter to scatter data from the distributed PETSc
    // representation to the distributed Lagrangian representation.
    VecScatter vec_scatter;
    ierr = VecScatterCreate(petsc_vec, lag_is, lagrangian_vec,
                            PETSC_NULL, &vec_scatter);  IBTK_CHKERRQ(ierr);

    // Scatter the values.
    ierr = VecScatterBegin(vec_scatter, petsc_vec, lagrangian_vec,
                           INSERT_VALUES, mode);  IBTK_CHKERRQ(ierr);
    ierr = VecScatterEnd(vec_scatter, petsc_vec, lagrangian_vec,
                         INSERT_VALUES, mode);  IBTK_CHKERRQ(ierr);

    // Cleanup allocated data.
    ierr = ISDestroy(lag_is);  IBTK_CHKERRQ(ierr);
    ierr = VecScatterDestroy(vec_scatter);  IBTK_CHKERRQ(ierr);
    return;
}// scatterData

void
LDataManager::beginNonlocalDataFill(
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    IBTK_TIMER_START(t_begin_nonlocal_data_fill);

    const int coarsest_ln = (coarsest_ln_in == -1) ? d_coarsest_ln : coarsest_ln_in;
    const int finest_ln = (finest_ln_in == -1) ? d_finest_ln : finest_ln_in;

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(coarsest_ln >= d_coarsest_ln &&
                coarsest_ln <= d_finest_ln);
    TBOX_ASSERT(finest_ln   >= d_coarsest_ln &&
                finest_ln   <= d_finest_ln);
#endif

    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        std::map<std::string,Pointer<LData> >& level_data = d_lag_mesh_data[level_number];
        for (std::map<std::string,Pointer<LData> >::iterator it = level_data.begin();
             it != level_data.end(); ++it)
        {
            it->second->beginGhostUpdate();
        }
    }

    IBTK_TIMER_STOP(t_begin_nonlocal_data_fill);
    return;
}// beginNonlocalDataFill

void
LDataManager::endNonlocalDataFill(
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    IBTK_TIMER_START(t_end_nonlocal_data_fill);

    const int coarsest_ln = (coarsest_ln_in == -1) ? d_coarsest_ln : coarsest_ln_in;
    const int finest_ln = (finest_ln_in == -1) ? d_finest_ln : finest_ln_in;

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(coarsest_ln >= d_coarsest_ln && coarsest_ln <= d_finest_ln);
    TBOX_ASSERT(finest_ln   >= d_coarsest_ln && finest_ln   <= d_finest_ln);
#endif

    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        std::map<std::string,Pointer<LData> >& level_data = d_lag_mesh_data[level_number];
        for (std::map<std::string,Pointer<LData> >::iterator it = level_data.begin();
             it != level_data.end(); ++it)
        {
            it->second->endGhostUpdate();
        }
    }

    IBTK_TIMER_STOP(t_end_nonlocal_data_fill);
    return;
}// endNonlocalDataFill

int
LDataManager::computeNodeDistribution(
    std::vector<int>& local_lag_indices,
    std::vector<int>& nonlocal_lag_indices,
    AO& ao,
    std::vector<int>& local_petsc_indices,
    std::vector<int>& nonlocal_petsc_indices,
    unsigned int& num_nodes,
    unsigned int& node_offset,
    const int level_number)
{
    IBTK_TIMER_START(t_compute_node_distribution);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(level_number >= d_coarsest_ln &&
                level_number <= d_finest_ln);
#endif

    local_lag_indices.clear();
    nonlocal_lag_indices.clear();
    local_petsc_indices.clear();
    nonlocal_petsc_indices.clear();

    // Determine the Lagrangian IDs of all of the Lagrangian nodes on the
    // specified level of the patch hierarchy.
    //
    // We differentiate between nodes that are local to the processor
    // (i.e. nodes that live in the interior of a patch owned by the processor)
    // and nodes that are non-local (i.e. nodes that live in the interior of a
    // patch owned by a different processor).
    //
    // It is important to emphasize that while a local node by definition lives
    // on the interior of some patch on this processor, it may also live in the
    // ghost cell regions of other patches owned by this processor.
    //
    // Non-local nodes ONLY appear in ghost cells for on processor patches.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);

    // Collect the local nodes.
    std::set<SAMRAI::tbox::Pointer<LNode>,LNodeIndexLagrangianIndexComp> local_node_idxs;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<LNodeSetData> lag_node_index_data = patch->getPatchData(d_lag_node_index_current_idx);
        for (LNodeSetData::DataIterator it = lag_node_index_data->data_begin(patch_box);
             it != lag_node_index_data->data_end(); ++it)
        {
            LNodeSet::value_type& node_idx = *it;
            local_node_idxs.insert(node_idx);
        }
    }

    // Assign local indices to the local nodes.
    unsigned int local_offset = 0;
    std::map<int,int> lag_idx_to_petsc_idx;
    for (std::set<SAMRAI::tbox::Pointer<LNode>,LNodeIndexLagrangianIndexComp>::const_iterator cit = local_node_idxs.begin();
         cit != local_node_idxs.end(); ++cit)
    {
        LNode* const node_idx = *cit;
        const int lag_idx = node_idx->getLagrangianIndex();
        const int petsc_idx = local_offset++;
        node_idx->setLocalPETScIndex(petsc_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(lag_idx_to_petsc_idx.find(lag_idx) == lag_idx_to_petsc_idx.end());
#endif
        lag_idx_to_petsc_idx[lag_idx] = petsc_idx;
        local_lag_indices.push_back(lag_idx);
    }

    // Determine the Lagrangian indices of the nonlocal nodes.
    std::set<int> nonlocal_lag_idx_set;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<LNodeSetData> lag_node_index_data = patch->getPatchData(d_lag_node_index_current_idx);
        const Box<NDIM>& ghost_box = lag_node_index_data->getGhostBox();
        for (LNodeSetData::DataIterator it = lag_node_index_data->data_begin(ghost_box);
             it != lag_node_index_data->data_end(); ++it)
        {
            if (patch_box.contains(it.getCellIndex())) continue;

            LNodeSet::value_type& node_idx = *it;
            const int lag_idx = node_idx->getLagrangianIndex();
            if (lag_idx_to_petsc_idx.find(lag_idx) == lag_idx_to_petsc_idx.end())
            {
                // If we have not already assigned a local index to this
                // Lagrangian index, then it must be a nonlocal index.
                nonlocal_lag_idx_set.insert(lag_idx);
            }
        }
    }

    // Assign local indices to the nonlocal nodes.
    for (std::set<int>::const_iterator cit = nonlocal_lag_idx_set.begin();
         cit != nonlocal_lag_idx_set.end(); ++cit)
    {
        const int lag_idx = *cit;
        lag_idx_to_petsc_idx[lag_idx] = local_offset++;
        nonlocal_lag_indices.push_back(lag_idx);
    }
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<LNodeSetData> lag_node_index_data = patch->getPatchData(d_lag_node_index_current_idx);
        const Box<NDIM>& ghost_box = lag_node_index_data->getGhostBox();
        for (LNodeSetData::DataIterator it = lag_node_index_data->data_begin(ghost_box);
             it != lag_node_index_data->data_end(); ++it)
        {
            if (patch_box.contains(it.getCellIndex())) continue;

            LNodeSet::value_type& node_idx = *it;
            const int lag_idx = node_idx->getLagrangianIndex();
            node_idx->setLocalPETScIndex(lag_idx_to_petsc_idx[lag_idx]);
        }
    }

    // Compute the new PETSc global ordering and initialize the AO object.
    int ierr;

    // Determine how many nodes are on each processor to calculate the PETSc
    // indexing scheme.
    const unsigned int    num_local_nodes =    local_lag_indices.size();
    const unsigned int num_nonlocal_nodes = nonlocal_lag_indices.size();

    if (local_offset != (num_local_nodes+num_nonlocal_nodes))
    {
        TBOX_ERROR("LDataManager::computeNodeDistribution()"       << "\n" <<
                   "  local_offset       = " << local_offset       << "\n" <<
                   "  num_local_nodes    = " << num_local_nodes    << "\n" <<
                   "  num_nonlocal_nodes = " << num_nonlocal_nodes << "\n");
    }

    computeNodeOffsets(num_nodes, node_offset, num_local_nodes);

    // Determine the PETSc ordering and setup the new AO object.
    const int num_proc_nodes = num_local_nodes + num_nonlocal_nodes;

    std::vector<int> node_indices;
    node_indices.reserve(num_proc_nodes);
    node_indices.insert(node_indices.end(),
                        local_lag_indices.begin(),
                        local_lag_indices.end());

    local_petsc_indices.resize(num_local_nodes);
    for (unsigned int k = 0; k < num_local_nodes; ++k)
    {
        local_petsc_indices[k] = node_offset+k;
    }

    if (ao)
    {
        ierr = AODestroy(ao);
        IBTK_CHKERRQ(ierr);
    }

    ierr = AOCreateBasic(PETSC_COMM_WORLD,
                         num_local_nodes,
                         num_local_nodes > 0 ? &       node_indices[0] : NULL,
                         num_local_nodes > 0 ? &local_petsc_indices[0] : NULL, &ao);
    IBTK_CHKERRQ(ierr);

    // Determine the PETSc local to global mapping (including PETSc Vec ghost
    // indices).
    //
    // NOTE: After this operation, node_indices are in the global PETSc
    // ordering.
    node_indices.insert(node_indices.end(),
                        nonlocal_lag_indices.begin(),
                        nonlocal_lag_indices.end());

    ierr = AOApplicationToPetsc(
        ao,
        (num_proc_nodes > 0 ? num_proc_nodes   : static_cast<int>(s_ao_dummy.size())),
        (num_proc_nodes > 0 ? &node_indices[0] : &s_ao_dummy[0]));
    IBTK_CHKERRQ(ierr);

    // Keep track of the global PETSc indices of the ghost nodes.
    nonlocal_petsc_indices.clear();
    nonlocal_petsc_indices.reserve(num_nonlocal_nodes);

    nonlocal_petsc_indices.insert(nonlocal_petsc_indices.end(),
                                  node_indices.begin()+num_local_nodes,
                                  node_indices.end());

    // Store the global PETSc index in the local LNode objects.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<LNodeSetData> lag_node_index_data = patch->getPatchData(d_lag_node_index_current_idx);
        const Box<NDIM>& ghost_box = lag_node_index_data->getGhostBox();
        for (LNodeSetData::DataIterator it = lag_node_index_data->data_begin(ghost_box);
             it != lag_node_index_data->data_end(); ++it)
        {
            LNodeSet::value_type& node_idx = *it;
            node_idx->setGlobalPETScIndex(node_indices[node_idx->getLocalPETScIndex()]);
        }
    }

    IBTK_TIMER_STOP(t_compute_node_distribution);
    return 0;
}// computeNodeDistribution

void
LDataManager::computeNodeOffsets(
    unsigned int& num_nodes,
    unsigned int& node_offset,
    const unsigned int& num_local_nodes)
{
    IBTK_TIMER_START(t_compute_node_offsets);

    const int mpi_size = SAMRAI_MPI::getNodes();
    const int mpi_rank = SAMRAI_MPI::getRank();

    std::vector<int> num_nodes_proc(mpi_size,0);

    SAMRAI_MPI::allGather(num_local_nodes, &num_nodes_proc[0]);

    node_offset = std::accumulate(num_nodes_proc.begin(),
                                  num_nodes_proc.begin()+mpi_rank, 0);

    num_nodes = std::accumulate(num_nodes_proc.begin()+mpi_rank,
                                num_nodes_proc.end(), node_offset);

    IBTK_TIMER_STOP(t_compute_node_offsets);
    return;
}// computeNodeOffsets

void
LDataManager::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();

    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to "
                   << d_object_name << " not found in restart file.");
    }

    int ver = db->getInteger("LDATA_MANAGER_VERSION");
    if (ver != LDATA_MANAGER_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different than class version.");
    }

    d_coarsest_ln = db->getInteger("d_coarsest_ln");
    d_finest_ln   = db->getInteger("d_finest_ln"  );
    d_alpha_work  = db->getDouble ("d_alpha_work" );
    d_beta_work   = db->getDouble ("d_beta_work"  );

    // Resize some arrays.
    d_level_contains_lag_data       .resize(d_finest_ln+1,false);
    d_strct_name_to_strct_id_map    .resize(d_finest_ln+1);
    d_strct_id_to_strct_name_map    .resize(d_finest_ln+1);
    d_strct_id_to_lag_idx_range_map .resize(d_finest_ln+1);
    d_last_lag_idx_to_strct_id_map  .resize(d_finest_ln+1);
    d_inactive_strcts               .resize(d_finest_ln+1);
    d_displaced_strct_ids           .resize(d_finest_ln+1);
    d_displaced_strct_bounding_boxes.resize(d_finest_ln+1);
    d_displaced_strct_lnode_idxs    .resize(d_finest_ln+1);
    d_displaced_strct_lnode_posns   .resize(d_finest_ln+1);
    d_lag_mesh                      .resize(d_finest_ln+1);
    d_lag_mesh_data                 .resize(d_finest_ln+1);
    d_needs_synch                   .resize(d_finest_ln+1,false);
    d_ao                            .resize(d_finest_ln+1);
    d_num_nodes                     .resize(d_finest_ln+1);
    d_node_offset                   .resize(d_finest_ln+1);
    d_local_lag_indices             .resize(d_finest_ln+1);
    d_nonlocal_lag_indices          .resize(d_finest_ln+1);
    d_local_petsc_indices           .resize(d_finest_ln+1);
    d_nonlocal_petsc_indices        .resize(d_finest_ln+1);

    // Read in data that is stored on a level-by-level basis.
    for (int level_number = d_coarsest_ln; level_number <= d_finest_ln; ++level_number)
    {
        std::ostringstream stream;
        stream << "level_" << level_number;
        const std::string level_db_name = stream.str();
        Pointer<Database> level_db = db->getDatabase(level_db_name);

        d_level_contains_lag_data[level_number] = level_db->getBool("d_level_contains_lag_data");

        if (!d_level_contains_lag_data[level_number]) continue;

        const int n_lstructs = level_db->getIntegerWithDefault("n_lstructs",0);
        std::vector<int> lstruct_ids(n_lstructs), lstruct_lag_idx_range_first(n_lstructs), lstruct_lag_idx_range_second(n_lstructs), lstruct_activation(n_lstructs);
        std::vector<std::string> lstruct_names(n_lstructs);
        if (n_lstructs > 0)
        {
            level_db->getIntegerArray("lstruct_ids", &lstruct_ids[0], lstruct_ids.size());
            level_db->getIntegerArray("lstruct_lag_idx_range_first", &lstruct_lag_idx_range_first[0], lstruct_lag_idx_range_first.size());
            level_db->getIntegerArray("lstruct_lag_idx_range_second", &lstruct_lag_idx_range_second[0], lstruct_lag_idx_range_second.size());
            level_db->getIntegerArray("lstruct_activation", &lstruct_activation[0], lstruct_activation.size());
            level_db->getStringArray("lstruct_names", &lstruct_names[0], lstruct_names.size());
        }
        for (int k = 0; k < n_lstructs; ++k)
        {
            d_strct_id_to_strct_name_map   [level_number][lstruct_ids[k]] = lstruct_names[k];
            d_strct_id_to_lag_idx_range_map[level_number][lstruct_ids[k]] = std::make_pair(lstruct_lag_idx_range_first[k], lstruct_lag_idx_range_second[k]);
            if (lstruct_activation[k] == 1)
            {
                d_inactive_strcts[level_number].addItem(k);
            }
        }
        d_inactive_strcts[level_number].communicateData();

        for (std::map<int,std::string>::const_iterator cit(d_strct_id_to_strct_name_map[level_number].begin());
             cit != d_strct_id_to_strct_name_map[level_number].end(); ++cit)
        {
            d_strct_name_to_strct_id_map[level_number][cit->second] = cit->first;
        }

        for (std::map<int,std::pair<int,int> >::const_iterator cit(d_strct_id_to_lag_idx_range_map[level_number].begin());
             cit != d_strct_id_to_lag_idx_range_map[level_number].end(); ++cit)
        {
            d_last_lag_idx_to_strct_id_map[level_number][cit->second.second-1] = cit->first;
        }

        const int n_ldata_names = level_db->getIntegerWithDefault("n_ldata_names",0);
        std::vector<std::string> ldata_names(n_ldata_names);
        if (!ldata_names.empty())
        {
            level_db->getStringArray("ldata_names",
                                     &ldata_names[0],
                                     n_ldata_names);
        }

        std::set<int> data_depths;
        for (std::vector<std::string>::iterator it = ldata_names.begin();
             it != ldata_names.end(); ++it)
        {
            const std::string& ldata_name = *it;
            d_lag_mesh_data[level_number][ldata_name] = new LData(level_db->getDatabase(ldata_name));
            data_depths.insert(d_lag_mesh_data[level_number][ldata_name]->getDepth());
        }

        d_num_nodes  [level_number] = level_db->getInteger("d_num_nodes"  );
        d_node_offset[level_number] = level_db->getInteger("d_node_offset");

        const int n_local_lag_indices = level_db->getInteger("n_local_lag_indices");
        if (n_local_lag_indices > 0)
        {
            d_local_lag_indices[level_number].resize(n_local_lag_indices);
            level_db->getIntegerArray("d_local_lag_indices",
                                      &d_local_lag_indices[level_number][0],
                                      n_local_lag_indices);
        }
        const int n_nonlocal_lag_indices = level_db->getInteger("n_nonlocal_lag_indices");
        if (n_nonlocal_lag_indices > 0)
        {
            d_nonlocal_lag_indices[level_number].resize(n_nonlocal_lag_indices);
            level_db->getIntegerArray("d_nonlocal_lag_indices",
                                      &d_nonlocal_lag_indices[level_number][0],
                                      n_nonlocal_lag_indices);
        }
        const int n_local_petsc_indices = level_db->getInteger("n_local_petsc_indices");
        if (n_local_petsc_indices > 0)
        {
            d_local_petsc_indices[level_number].resize(n_local_petsc_indices);
            level_db->getIntegerArray("d_local_petsc_indices",
                                      &d_local_petsc_indices[level_number][0],
                                      n_local_petsc_indices);
        }
        const int n_nonlocal_petsc_indices = level_db->getInteger("n_nonlocal_petsc_indices");
        if (n_nonlocal_petsc_indices > 0)
        {
            d_nonlocal_petsc_indices[level_number][1].resize(n_nonlocal_petsc_indices);
            level_db->getIntegerArray("d_nonlocal_petsc_indices",
                                      &d_nonlocal_petsc_indices[level_number][1][0],
                                      n_nonlocal_petsc_indices);

            // Rebuild the nonlocal PETSc indices for the other data depths.
            for (std::set<int>::const_iterator it = data_depths.begin();
                 it != data_depths.end(); ++it)
            {
                const int depth = *it;

                d_nonlocal_petsc_indices[level_number][depth].
                    resize(d_nonlocal_petsc_indices[level_number][1].size());

                std::transform(d_nonlocal_petsc_indices[level_number][    1].begin(),
                               d_nonlocal_petsc_indices[level_number][    1].end  (),
                               d_nonlocal_petsc_indices[level_number][depth].begin(),
                               std::bind2nd(std::multiplies<int>(),depth));
            }
        }

        // Rebuild the application ordering.
        int ierr;
        ierr = AOCreateBasic(PETSC_COMM_WORLD,
                             n_local_lag_indices,
                             n_local_lag_indices > 0 ? &d_local_lag_indices  [level_number][0] : NULL,
                             n_local_lag_indices > 0 ? &d_local_petsc_indices[level_number][0] : NULL,
                             &d_ao[level_number]);
        IBTK_CHKERRQ(ierr);
    }
    return;
}// getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
