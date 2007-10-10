// Filename: LDataManager.C
// Last modified: <10.Oct.2007 18:13:19 griffith@box221.cims.nyu.edu>
// Created on 01 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

#include "LDataManager.h"

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
#include <ibamr/LEInteractor.h>
#include <ibamr/LNodeIndexSet.h>
#include <ibamr/LNodeIndexData2.h>
#include <ibamr/LNodeLevelData.h>
#include <ibamr/LagSiloDataWriter.h>

#if (NDIM == 3)
#include <ibamr/LagM3DDataWriter.h>
#endif

// STOOLS INCLUDES
#include <stools/PETSC_SAMRAI_ERROR.h>
#include <stools/STOOLS_Utilities.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <BoxList.h>
#include <CartesianCellDoubleWeightedAverage.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <CoarsenOperator.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <RefineOperator.h>
#include <VariableDatabase.h>
#include <tbox/MPI.h>
#include <tbox/RestartManager.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/ShutdownRegistry.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <functional>
#include <numeric>
#include <set>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_map_lagrangian_to_petsc;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_map_petsc_to_lagrangian;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_begin_data_redistribution;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_end_data_redistribution;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_update_workload_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_update_irregular_cell_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_restore_location_pointers;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_invalidate_location_pointers;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_hierarchy_configuration;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_apply_gradient_detector;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_put_to_database;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_begin_nonlocal_data_fill;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_end_nonlocal_data_fill;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_node_distribution;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_node_offsets;

// Assume max(U)dt/dx <= 1.
static const int CFL_WIDTH = 1;

// Version of LDataManager restart file data.
static const int LDATA_MANAGER_VERSION = 1;

inline SAMRAI::pdat::CellIndex<NDIM>
get_canonical_cell_index(
    const SAMRAI::pdat::CellIndex<NDIM>& cell_idx,
    const SAMRAI::hier::Box<NDIM>& domain_box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift)
{
    if (periodic_shift == SAMRAI::hier::IntVector<NDIM>(0)) return cell_idx;

    if (domain_box.contains(cell_idx)) return cell_idx;

    SAMRAI::pdat::CellIndex<NDIM> shifted_idx = cell_idx;
    for (int d = 0; d < NDIM; ++d)
    {
        if      (shifted_idx(d) < domain_box.lower()(d)) shifted_idx(d) += periodic_shift(d);
        else if (shifted_idx(d) > domain_box.upper()(d)) shifted_idx(d) -= periodic_shift(d);
    }

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(domain_box.contains(shifted_idx));
#endif

    return shifted_idx;
}// get_canonical_cell_index

struct LNodeIndexLessThan
    : std::binary_function<SAMRAI::tbox::Pointer<LNodeIndex>,SAMRAI::tbox::Pointer<LNodeIndex>,bool>
{
    inline bool
    operator()(
        const SAMRAI::tbox::Pointer<LNodeIndex>& lhs,
        const SAMRAI::tbox::Pointer<LNodeIndex>& rhs) const
        {
            return *lhs < *rhs;
        }
};

struct CellIndexFortranOrder
    : std::binary_function<SAMRAI::pdat::CellIndex<NDIM>,SAMRAI::pdat::CellIndex<NDIM>,bool>
{
    inline bool
    operator()(
        const SAMRAI::pdat::CellIndex<NDIM>& lhs,
        const SAMRAI::pdat::CellIndex<NDIM>& rhs) const
        {
            return (lhs(0) < rhs(0)
#if (NDIM>1)
                    || (lhs(0) == rhs(0) && lhs(1) < rhs(1))
#if (NDIM>2)
                    || (lhs(0) == rhs(0) && lhs(1) == rhs(1) && lhs(2) < rhs(2))
#endif
#endif
                    );
        }
};

struct BeginLNodeLevelDataNonlocalFill
    : std::unary_function<pair<std::string,SAMRAI::tbox::Pointer<LNodeLevelData> >,void>
{
    inline void
    operator()(
        const pair<std::string,SAMRAI::tbox::Pointer<LNodeLevelData> >& data) const
        {
            data.second->beginGhostUpdate();
            return;
        }
};

class RestoreLNodeIndexLocationPointers
    : std::unary_function<SAMRAI::tbox::Pointer<LNodeIndex>,void>
{
public:
    inline
    RestoreLNodeIndexLocationPointers(
        double* const X_arr)
        : d_X_arr(X_arr)
        {
            return;
        }

    inline void
    operator()(
        const SAMRAI::tbox::Pointer<LNodeIndex>& data) const
        {
            data->setNodeLocation(
                &(d_X_arr[NDIM*data->getLocalPETScIndex()]));
            return;
        }

private:
    double* const d_X_arr;
};

struct InvalidateLNodeIndexLocationPointers
    : std::unary_function<SAMRAI::tbox::Pointer<LNodeIndex>,void>
{
    inline void
    operator()(
        const SAMRAI::tbox::Pointer<LNodeIndex>& data) const
        {
            data->setNodeLocation(NULL);
            return;
        }
};

struct EndLNodeLevelDataNonlocalFill
    : std::unary_function<pair<std::string,SAMRAI::tbox::Pointer<LNodeLevelData> >,void>
{
    inline void
    operator()(
        const pair<std::string,SAMRAI::tbox::Pointer<LNodeLevelData> >& data) const
        {
            data.second->endGhostUpdate();
            return;
        }
};

struct GetLagrangianIndex
    : std::unary_function<SAMRAI::tbox::Pointer<LNodeIndex>,int>
{
    inline int
    operator()(
        const SAMRAI::tbox::Pointer<LNodeIndex>& index) const
        {
            return index->getLagrangianIndex();
        }
};

class SetLocalPETScIndex
    : public std::unary_function<SAMRAI::tbox::Pointer<LNodeIndex>,void>,
      public std::unary_function<void,int>
{
public:
    inline
    SetLocalPETScIndex(
        const int offset)
        : d_current_index(offset)
        {
            return;
        }

    inline void
    operator()(
        const SAMRAI::tbox::Pointer<LNodeIndex>& index)
        {
            index->setLocalPETScIndex(d_current_index++);
            return;
        }

    inline int
    operator()()
        {
            return d_current_index++;
        }

private:
    int d_current_index;
};

class GetLocalPETScIndexFromIDSet
    : public std::unary_function<SAMRAI::tbox::Pointer<LNodeIndex>,void>,
      public std::unary_function<void,int>
{
public:
    inline
    GetLocalPETScIndexFromIDSet(
        const LNodeIndexSet::const_iterator& begin)
        : d_index(begin)
        {
            return;
        }

    inline void
    operator()(
        const SAMRAI::tbox::Pointer<LNodeIndex>& index)
        {
            index->setLocalPETScIndex((*(d_index++))->getLocalPETScIndex());
            return;
        }

    inline int
    operator()()
        {
            return (*(d_index++))->getLocalPETScIndex();
        }

private:
    LNodeIndexSet::const_iterator d_index;
};
}

const std::string LDataManager::POSN_DATA_NAME = "X";
const std::string LDataManager::VEL_DATA_NAME = "U";
std::map<std::string,LDataManager*> LDataManager::s_data_manager_instances;
bool LDataManager::s_registered_callback;
unsigned char LDataManager::s_shutdown_priority = 200;
std::vector<int> LDataManager::s_ao_dummy(1,-1);

LDataManager*
LDataManager::getManager(
    const std::string& name,
    const SAMRAI::hier::IntVector<NDIM>& ghosts,
    bool register_for_restart)
{
    if (s_data_manager_instances.find(name) ==
        s_data_manager_instances.end())
    {
        s_data_manager_instances[name] = new LDataManager(
            name, ghosts, register_for_restart);
    }
    if (!s_registered_callback)
    {
        SAMRAI::tbox::ShutdownRegistry::registerShutdownRoutine(
            freeAllManagers,s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_data_manager_instances[name];
}// getManager

void
LDataManager::freeAllManagers()
{
    for (std::map<std::string,LDataManager*>::iterator it =
             s_data_manager_instances.begin();
         it != s_data_manager_instances.end();
         ++it)
    {
        if ((*it).second)
        {
            delete (*it).second;
        }
        (*it).second = NULL;
    }
    return;
}// freeManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
LDataManager::setPatchHierarchy(
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
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
    assert(!d_hierarchy.isNull());
    assert((coarsest_ln >= 0) &&
           (finest_ln >= coarsest_ln) &&
           (finest_ln <= d_hierarchy->getFinestLevelNumber()));
#endif
    // Destroy any un-needed AO objects.
    int ierr;
    for (int ln = SAMRAI::tbox::Utilities::imax(d_coarsest_ln,0);
         (ln <= d_finest_ln) && (ln < coarsest_ln); ++ln)
    {
        if (d_ao[ln])
        {
            ierr = AODestroy(d_ao[ln]);
            PETSC_SAMRAI_ERROR(ierr);
        }
    }
    for (int ln = finest_ln+1; ln <= d_finest_ln; ++ln)
    {
        if (d_ao[ln])
        {
            ierr = AODestroy(d_ao[ln]);
            PETSC_SAMRAI_ERROR(ierr);
        }
    }

    // Reset the level numbers.
    d_coarsest_ln = coarsest_ln;
    d_finest_ln   = finest_ln;

    // Resize some arrays.
    d_level_contains_lag_data.resize(d_finest_ln+1);
    d_lag_quantity_data      .resize(d_finest_ln+1);
    d_needs_synch            .resize(d_finest_ln+1,false);
    d_ao                     .resize(d_finest_ln+1);
    d_num_nodes              .resize(d_finest_ln+1);
    d_node_offset            .resize(d_finest_ln+1);
    d_local_lag_indices      .resize(d_finest_ln+1);
    d_nonlocal_lag_indices   .resize(d_finest_ln+1);
    d_local_petsc_indices    .resize(d_finest_ln+1);
    d_nonlocal_petsc_indices .resize(d_finest_ln+1);

    return;
}// resetLevels

void
LDataManager::registerLNodeInitStrategy(
    SAMRAI::tbox::Pointer<LNodeInitStrategy> lag_init)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!lag_init.isNull());
#endif
    d_lag_init = lag_init;
    return;
}// registerLNodeInitStrategy

void
LDataManager::freeLNodeInitStrategy()
{
    d_lag_init.setNull();
    return;
}// freeLNodeInitStrategy

void
LDataManager::registerVisItDataWriter(
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!visit_writer.isNull());
#endif
    d_visit_writer = visit_writer;
    if (d_output_workload)
    {
        d_visit_writer->registerPlotQuantity(
            "workload estimate", "SCALAR", d_workload_idx);
    }
    if (d_output_node_count)
    {
        d_visit_writer->registerPlotQuantity(
            "node count", "SCALAR", d_node_count_idx);
    }
    if (d_output_mpi_proc)
    {
        d_visit_writer->registerPlotQuantity(
            "MPI process", "SCALAR", d_mpi_proc_idx);
    }
    return;
}// registerVisItDataWriter

void
LDataManager::registerLagSiloDataWriter(
    SAMRAI::tbox::Pointer<LagSiloDataWriter> silo_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!silo_writer.isNull());
#endif
    d_silo_writer = silo_writer;
    return;
}// registerLagSiloDataWriter

#if (NDIM == 3)
void
LDataManager::registerLagM3DDataWriter(
    SAMRAI::tbox::Pointer<LagM3DDataWriter> m3D_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!m3D_writer.isNull());
#endif
    d_m3D_writer = m3D_writer;
    return;
}// registerLagM3DDataWriter
#endif

void
LDataManager::registerLoadBalancer(
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!load_balancer.isNull());
#endif
    d_load_balancer = load_balancer;
    return;
}// return

bool
LDataManager::levelContainsLagrangianData(
    const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(level_number >= 0);
#endif
    if (!(d_coarsest_ln <= level_number &&
          d_finest_ln   >= level_number))
    {
        return false;
    }
    else
    {
        return d_level_contains_lag_data[level_number];
    }
}// levelContainsLagrangianData

int
LDataManager::getNumberOfNodes(
    const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(level_number >= 0);
    assert(d_coarsest_ln <= level_number &&
           d_finest_ln   >= level_number);
#endif
    return d_num_nodes[level_number];
}// getNumberOfNodes

int
LDataManager::getNumberOfLocalNodes(
    const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(level_number >= 0);
    assert(d_coarsest_ln <= level_number &&
           d_finest_ln   >= level_number);
#endif
    return d_local_lag_indices[level_number].size();
}// getNumberOfLocalNodes

int
LDataManager::getGlobalNodeOffset(
    const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(level_number >= 0);
    assert(d_coarsest_ln <= level_number &&
           d_finest_ln   >= level_number);
#endif
    return d_node_offset[level_number];
}// getGlobalNodeOffset

SAMRAI::tbox::Pointer<LNodeLevelData>
LDataManager::getLNodeLevelData(
    const std::string& quantity_name,
    const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_lag_quantity_data[level_number].find(quantity_name) !=
           d_lag_quantity_data[level_number].end());
    assert(level_number >= 0);
    assert(d_coarsest_ln <= level_number &&
           d_finest_ln   >= level_number);
#endif
    return d_lag_quantity_data[level_number][quantity_name];
}// getLNodeLevelData

SAMRAI::tbox::Pointer<LNodeLevelData>
LDataManager::createLNodeLevelData(
    const std::string& quantity_name,
    const int level_number,
    const int depth,
    const bool maintain_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!maintain_data ||
           (d_lag_quantity_data[level_number].find(quantity_name) ==
            d_lag_quantity_data[level_number].end()));
    assert(level_number >= 0);
    assert(d_coarsest_ln <= level_number &&
           d_finest_ln   >= level_number);
    assert(depth > 0);
#endif

    bool created_nonlocal_petsc_indices = false;

    if (d_nonlocal_petsc_indices[level_number].find(depth) ==
        d_nonlocal_petsc_indices[level_number].end())
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(depth != 1);
#endif
        d_nonlocal_petsc_indices[level_number][depth].resize(
            d_nonlocal_petsc_indices[level_number][1].size());

        std::transform(d_nonlocal_petsc_indices[level_number][    1].begin(),
                       d_nonlocal_petsc_indices[level_number][    1].end(),
                       d_nonlocal_petsc_indices[level_number][depth].begin(),
                       std::bind2nd(std::multiplies<int>(),depth));

        created_nonlocal_petsc_indices = true;
    }

    SAMRAI::tbox::Pointer<LNodeLevelData> ret_val = new LNodeLevelData(
        quantity_name,getNumberOfLocalNodes(level_number),depth,
        d_nonlocal_petsc_indices[level_number][depth]);

    if (maintain_data)
    {
        d_lag_quantity_data[level_number][quantity_name] = ret_val;
    }

    if (!maintain_data && created_nonlocal_petsc_indices)
    {
        d_nonlocal_petsc_indices[level_number].erase(depth);
    }

    return ret_val;
}// createLNodeLevelData

int
LDataManager::getLNodeIndexPatchDescriptorIndex() const
{
    return d_lag_node_index_current_idx;
}// getLNodeIndexPatchDescriptorIndex

int
LDataManager::getWorkloadPatchDescriptorIndex() const
{
    return d_workload_idx;
}// getWorkloadPatchDescriptorIndex

int
LDataManager::getNodeCountPatchDescriptorIndex() const
{
    return d_node_count_idx;
}// getNodeCountPatchDescriptorIndex

int
LDataManager::getIrregularCellPatchDescriptorIndex() const
{
    return d_irregular_cell_idx;
}// getIrregularCellPatchDescriptorIndex

int
LDataManager::getProcMappingPatchDescriptorIndex() const
{
    return d_mpi_proc_idx;
}// getProcMappingPatchDescriptorIndex

void
LDataManager::mapLagrangianToPETSc(
    std::vector<int>& inds,
    const int level_number) const
{
    t_map_lagrangian_to_petsc->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_coarsest_ln <= level_number &&
           d_finest_ln   >= level_number);
#endif

    const int ierr = AOApplicationToPetsc(
        d_ao[level_number],
        (!inds.empty() ? static_cast<int>(inds.size()) : static_cast<int>(s_ao_dummy.size())),
        (!inds.empty() ? &inds[0]                      : &s_ao_dummy[0]));
    PETSC_SAMRAI_ERROR(ierr);

    t_map_lagrangian_to_petsc->stop();
    return;
}// mapLagrangianToPETSc

void
LDataManager::mapPETScToLagrangian(
    std::vector<int>& inds,
    const int level_number) const
{
    t_map_petsc_to_lagrangian->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_coarsest_ln <= level_number &&
           d_finest_ln   >= level_number);
#endif

    const int ierr = AOPetscToApplication(
        d_ao[level_number],
        (!inds.empty() ? static_cast<int>(inds.size()) : static_cast<int>(s_ao_dummy.size())),
        (!inds.empty() ? &inds[0]                      : &s_ao_dummy[0]));
    PETSC_SAMRAI_ERROR(ierr);

    t_map_petsc_to_lagrangian->stop();
    return;
}// mapPETScToLagrangian

void
LDataManager::scatterLagrangianToPETSc(
    Vec& lagrangian_vec,
    Vec& petsc_vec,
    const int level_number) const
{
    scatterData(lagrangian_vec,petsc_vec,level_number,SCATTER_FORWARD);
    return;
}// scatterLagrangianToPETSc

void
LDataManager::scatterPETScToLagrangian(
    Vec& petsc_vec,
    Vec& lagrangian_vec,
    const int level_number) const
{
    scatterData(lagrangian_vec,petsc_vec,level_number,SCATTER_REVERSE);
    return;
}// scatterPETScToLagrangian

void
LDataManager::beginDataRedistribution(
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    t_begin_data_redistribution->start();

    const int coarsest_ln =
        (coarsest_ln_in == -1)
        ? d_coarsest_ln
        : coarsest_ln_in;
    const int finest_ln =
        (finest_ln_in == -1)
        ? d_finest_ln
        : finest_ln_in;

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(coarsest_ln >= d_coarsest_ln &&
           coarsest_ln <= d_finest_ln);
    assert(finest_ln   >= d_coarsest_ln &&
           finest_ln   <= d_finest_ln);
#endif

    const double* const gridXLower = d_grid_geom->getXLower();
    const double* const gridXUpper = d_grid_geom->getXUpper();
    double gridXLength[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
        gridXLength[d] = gridXUpper[d] - gridXLower[d];
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_level_contains_lag_data[ln])
        {
            if (d_needs_synch[ln])
            {
                TBOX_WARNING("LDataManager::beginDataRedistribution():\n" <<
                             "\tLNodeLevelData is not synchronized with LNodeIndexData2.\n" <<
                             "\tLagrangian node position data is probably invalid!\n");
            }

            // Update the ghost values of the Lagrangian nodal positions.
            d_lag_quantity_data[ln][POSN_DATA_NAME]->beginGhostUpdate();
            d_lag_quantity_data[ln][POSN_DATA_NAME]->endGhostUpdate();

            // Make sure that the location pointers are properly set for each
            // LNodeIndex.  They are directly used below to locate the
            // Lagrangian nodes.  They are also used to re-sort the node index
            // sets in an attempt to maximize data locality.
            restoreLocationPointers(ln,ln);
        }
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_level_contains_lag_data[ln])
        {
            // Update the index patch data on the level.
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            const SAMRAI::hier::IntVector<NDIM>& periodic_shift = d_grid_geom->getPeriodicShift(
                level->getRatio());
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
                    patch->getPatchGeometry();
                const SAMRAI::pdat::CellIndex<NDIM>& patch_lower = patch_box.lower();
                const SAMRAI::pdat::CellIndex<NDIM>& patch_upper = patch_box.upper();
                const double* const patchXLower = patch_geom->getXLower();
                const double* const patchXUpper = patch_geom->getXUpper();
                const double* const patchDx = patch_geom->getDx();
                const bool touches_periodic_bdry =
                    patch_geom->getTouchesPeriodicBoundary();

                SAMRAI::tbox::Pointer<LNodeIndexData2> current_idx_data =
                    patch->getPatchData(d_lag_node_index_current_idx);
                SAMRAI::tbox::Pointer<LNodeIndexData2> new_idx_data = new LNodeIndexData2(
                    current_idx_data->getBox(), current_idx_data->getGhostCellWidth());

                // Initialize LNodeIndexSet objects for each cell index in the
                // patch interior which will contain LNodeIndex objects AFTER
                // redistribution.
                //
                // We only keep nodes whose new locations are in the patch
                // interior.  That is to say, we only keep the nodes which the
                // patch will own after redistribution.
                const SAMRAI::hier::Box<NDIM> cfl_box = SAMRAI::hier::Box<NDIM>::grow(
                    patch_box, SAMRAI::hier::IntVector<NDIM>(CFL_WIDTH));

                for (LNodeIndexData2::Iterator it(cfl_box); it; it++)
                {
                    const SAMRAI::pdat::CellIndex<NDIM>& old_cell_idx = *it;
                    const LNodeIndexSet& old_node_set = (*current_idx_data)(old_cell_idx);
                    const bool patch_owns_node_at_old_loc = patch_box.contains(old_cell_idx);
                    const SAMRAI::hier::IntVector<NDIM>& offset = old_node_set.getPeriodicOffset();
                    double X_shifted[NDIM];

                    for (LNodeIndexSet::const_iterator n = old_node_set.begin(); n != old_node_set.end(); ++n)
                    {
                        const LNodeIndexSet::value_type& node_idx = *n;
                        double* const X = node_idx->getNodeLocation();
                        for (int d = 0; d < NDIM; ++d)
                        {
                            X_shifted[d] = X[d] + double(offset(d))*patchDx[d];
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
                            const SAMRAI::pdat::CellIndex<NDIM> new_cell_idx = STOOLS::STOOLS_Utilities::getCellIndex(
                                X_shifted, patchXLower, patchXUpper, patchDx, patch_lower, patch_upper);
                            LNodeIndexSet& new_node_set = (*new_idx_data)(new_cell_idx);
                            const LNodeIndexSet::iterator pos = lower_bound(
                                new_node_set.begin(), new_node_set.end(), node_idx, LNodeIndexLessThan());
                            new_node_set.d_set.insert(pos, node_idx);
                        }
                        else
                        {
                            // If a node leaves the patch via a periodic
                            // boundary, we have to make sure that its location
                            // is properly updated.
                            //
                            // NOTE: The location of a node is itself a
                            // *quantity* living on the Lagrangian mesh.  Until
                            // the Lagrangian data is redistributed, the patch
                            // still owns the LNodeLevelData associated with the
                            // old LNodeIndex distribution.  Thus it is the
                            // responsibility of this patch to update the
                            // location of the node if the node leaves via a
                            // periodic boundary.
                            if (touches_periodic_bdry && patch_owns_node_at_old_loc)
                            {
                                bool shifted = false;
                                SAMRAI::hier::IntVector<NDIM> offset = 0;
                                std::vector<double> D(NDIM,0.0);
                                static const int lower = 0;
                                static const int upper = 1;
                                for (int d = 0; d < NDIM; ++d)
                                {
                                    if      (patch_geom->getTouchesPeriodicBoundary(d,lower) &&
                                             X[d] < gridXLower[d])
                                    {
                                        shifted = true;
                                        offset(d) = +periodic_shift(d);
                                        D[d] = +gridXLength[d];
                                    }
                                    else if (patch_geom->getTouchesPeriodicBoundary(d,upper) &&
                                             X[d] >= gridXUpper[d])
                                    {
                                        shifted = true;
                                        offset(d) = -periodic_shift(d);
                                        D[d] = -gridXLength[d];
                                    }
                                }
                                if (shifted)
                                {
                                    for (int d = 0; d < NDIM; ++d)
                                    {
                                        X[d] += D[d];
                                    }
                                    node_idx->registerPeriodicShift(offset, D);
                                }
                            }
                        }
                    }
                }

                // Swap the old and new patch data pointers.
                patch->setPatchData(d_lag_node_index_current_idx, new_idx_data);
            }
        }
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_level_contains_lag_data[ln])
        {
            d_needs_synch[ln] = true;
        }
    }

    t_begin_data_redistribution->stop();
    return;
}// beginDataRedistribution

void
LDataManager::endDataRedistribution(
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    t_end_data_redistribution->start();

    const int coarsest_ln =
        (coarsest_ln_in == -1)
        ? d_coarsest_ln
        : coarsest_ln_in;
    const int finest_ln =
        (finest_ln_in == -1)
        ? d_finest_ln
        : finest_ln_in;

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(coarsest_ln >= d_coarsest_ln &&
           coarsest_ln <= d_finest_ln);
    assert(finest_ln   >= d_coarsest_ln &&
           finest_ln   <= d_finest_ln);
#endif

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_level_contains_lag_data[ln] && (!d_needs_synch[ln]))
        {
            TBOX_WARNING("LDataManager::endDataRedistribution():\n" <<
                         "\tLNodeLevelData is already synchronized with LNodeIndexData2.\n" <<
                         "\tlevel = " << ln << "\n");
        }
    }

    // Fill the ghost cells of each level.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_level_contains_lag_data[ln])
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

            level->allocatePatchData(d_scratch_data);

            const double current_time = 0.0;
            level->setTime(current_time, d_current_data);
            level->setTime(current_time, d_scratch_data);

            d_lag_node_index_bdry_fill_scheds[ln]->fillData(current_time);

            level->deallocatePatchData(d_scratch_data);
        }
    }

    // Define the PETSc data needed to communicate the LNodeLevelData from its
    // old configuration to its new configuration.
    int ierr;

    std::vector<AO> new_ao(finest_ln+1);

    std::vector<std::vector<Vec> > src(finest_ln+1);
    std::vector<std::vector<Vec> > dst(finest_ln+1);
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
    // patches---but overlapping patches are the work of the devil).  Nodes in
    // the ghost region of a patch will not in general be stored as contiguous
    // data, and no attempt is made to do so.
    std::vector<std::map<int,std::vector<int>*> > patch_interior_local_indices(finest_ln+1);
    std::vector<std::map<int,std::vector<int>*> > patch_ghost_local_indices   (finest_ln+1);

    // In the following loop over patch levels, we first compute the new
    // distribution data (e.g., all of these indices).
    //
    // Next, we use the old and new PETSc AO (application ordering) objects to
    // define a mapping from the old distribution to the new distribution.
    //
    // Finally, we create the new PETSc Vec (vector) objects which are used to
    // store the Lagrangian data in the new distribution and begin the process
    // of scattering data from the old configuration into the new one.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        // Reset the nonlocal PETSc indices.
        d_nonlocal_petsc_indices[ln].clear();

        // The destination indices.
        std::vector<int> dst_inds;
        bool dst_inds_set = false;

        // Compute the new data distribution and start scattering.
        if (d_level_contains_lag_data[ln])
        {
            std::map<std::string,SAMRAI::tbox::Pointer<LNodeLevelData> >& level_data =
                d_lag_quantity_data[ln];
            const std::vector<int>::size_type num_data = level_data.size();

            src[    ln].resize(num_data);
            dst[    ln].resize(num_data);
            scatter[ln].resize(num_data);

            // Obtain pointers to the patch local and ghost index sets.
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const int patch_num = patch->getPatchNumber();
                SAMRAI::tbox::Pointer<LNodeIndexData2> lag_node_idx_data =
                    patch->getPatchData(d_lag_node_index_current_idx);

                patch_interior_local_indices[ln][patch_num] =
                    &(lag_node_idx_data->d_interior_local_indices);
                patch_ghost_local_indices[ln][patch_num] =
                    &(lag_node_idx_data->d_ghost_local_indices);
            }

            // Get the new distribution of nodes for the level.
            //
            // NOTE: This process updates the local PETSc indices of the
            // LNodeIndexSet objects contained in the current patch.
            ierr = computeNodeDistribution(d_local_lag_indices   [ln],
                                           d_nonlocal_lag_indices[ln],
                                           new_ao[ln],
                                           d_local_petsc_indices   [ln],
                                           d_nonlocal_petsc_indices[ln][1],
                                           d_num_nodes[ln],
                                           d_node_offset[ln],
                                           patch_interior_local_indices[ln],
                                           patch_ghost_local_indices   [ln],
                                           ln);
            PETSC_SAMRAI_ERROR(ierr);

            num_local_nodes   [ln] = static_cast<int>(d_local_lag_indices   [ln].size());
            num_nonlocal_nodes[ln] = static_cast<int>(d_nonlocal_lag_indices[ln].size());

            std::map<std::string, SAMRAI::tbox::Pointer<LNodeLevelData> >::iterator it;
            int i;
            for (it = level_data.begin(), i = 0; it != level_data.end();
                 ++it, ++i)
            {
                SAMRAI::tbox::Pointer<LNodeLevelData> data = (*it).second;
#ifdef DEBUG_CHECK_ASSERTIONS
                assert(!data.isNull());
#endif
                const int depth = data->getDepth();

                // Determine the PETSc indices of the ghost nodes required in
                // the destination Vec.
                //
                // (This is only computed once for each data depth encountered,
                // including depth==1.)
                if (d_nonlocal_petsc_indices[ln].find(depth) ==
                    d_nonlocal_petsc_indices[ln].end())
                {
#ifdef DEBUG_CHECK_ASSERTIONS
                    assert(depth != 1);
#endif
                    d_nonlocal_petsc_indices[ln][depth].
                        resize(num_nonlocal_nodes[ln]);

                    std::transform(d_nonlocal_petsc_indices[ln][    1].begin(),
                                   d_nonlocal_petsc_indices[ln][    1].end(),
                                   d_nonlocal_petsc_indices[ln][depth].begin(),
                                   std::bind2nd(std::multiplies<int>(),depth));
                }

                // Determine the PETSc indices of the source nodes for use when
                // scattering values from the old configuration to the new
                // configuration.
                //
                // (This is only computed once for each data depth encountered,
                // including depth==1.)
                if (src_IS[ln].find(depth) == src_IS[ln].end())
                {
                    ierr = ISCreateStride(PETSC_COMM_WORLD,
                                          depth*num_local_nodes[ln],
                                          depth*d_node_offset[ln],
                                          1, &src_IS[ln][depth]);
                    PETSC_SAMRAI_ERROR(ierr);
                }

                // Determine the PETSc indices of the destination nodes for use
                // when scattering values from the old configuration to the new
                // configuration.
                //
                // (This is only computed once for each data depth encountered,
                // including depth==1.)
                if (!dst_inds_set)
                {
                    dst_inds = d_local_petsc_indices[ln];

                    ierr = AOPetscToApplication(
                        d_ao[ln], // the old AO
                        (num_local_nodes[ln] > 0 ? num_local_nodes[ln] : static_cast<int>(s_ao_dummy.size())),
                        (num_local_nodes[ln] > 0 ? &dst_inds[0]        : &s_ao_dummy[0]));
                    PETSC_SAMRAI_ERROR(ierr);
                    ierr = AOApplicationToPetsc(
                        new_ao[ln],
                        (num_local_nodes[ln] > 0 ? num_local_nodes[ln] : static_cast<int>(s_ao_dummy.size())),
                        (num_local_nodes[ln] > 0 ? &dst_inds[0]        : &s_ao_dummy[0]));
                    PETSC_SAMRAI_ERROR(ierr);

                    dst_inds_set = true;
                }

                if (dst_IS[ln].find(depth) == dst_IS[ln].end())
                {
                    if (depth == 1)
                    {
                        ierr = ISCreateGeneral(PETSC_COMM_WORLD,
                                               num_local_nodes[ln],
                                               &dst_inds[0],
                                               &dst_IS[ln][depth]);
                        PETSC_SAMRAI_ERROR(ierr);
                    }
                    else
                    {
                        std::vector<int> scaled_dst_inds(dst_inds.size());
                        std::transform(dst_inds.begin(), dst_inds.end(),
                                       scaled_dst_inds.begin(),
                                       std::bind2nd(std::multiplies<int>(),depth));

                        ierr = ISCreateBlock(PETSC_COMM_WORLD,
                                             depth, num_local_nodes[ln],
                                             &scaled_dst_inds[0],
                                             &dst_IS[ln][depth]);
                        PETSC_SAMRAI_ERROR(ierr);
                    }
                }

                // Create the destination Vec and start scattering data from the
                // old configuration to the new configuration.
                src[ln][i] = data->getGlobalVec();

                if (depth == 1)
                {
                    ierr = VecCreateGhost(PETSC_COMM_WORLD,
                                          num_local_nodes[ln], PETSC_DECIDE,
                                          num_nonlocal_nodes[ln],
                                          &d_nonlocal_petsc_indices[ln][depth][0],
                                          &dst[ln][i]);
                    PETSC_SAMRAI_ERROR(ierr);
                }
                else
                {
                    ierr = VecCreateGhostBlock(PETSC_COMM_WORLD, depth,
                                               depth*num_local_nodes[ln], PETSC_DECIDE,
                                               num_nonlocal_nodes[ln],
                                               &d_nonlocal_petsc_indices[ln][depth][0],
                                               &dst[ln][i]);
                    PETSC_SAMRAI_ERROR(ierr);
                }
                ierr = VecSetBlockSize(dst[ln][i], depth);
                PETSC_SAMRAI_ERROR(ierr);

                if (scatter_template[ln].find(depth) ==
                    scatter_template[ln].end())
                {
                    ierr = VecScatterCreate(src[ln][i], src_IS[ln][depth],
                                            dst[ln][i], dst_IS[ln][depth],
                                            &scatter_template[ln][depth]);
                    PETSC_SAMRAI_ERROR(ierr);
                    scatter[ln][i] = scatter_template[ln][depth];
                }
                else if (scatter_template[ln][depth]->copy != NULL)
                {
                    ierr = VecScatterCopy(scatter_template[ln][depth],
                                          &scatter[ln][i]);
                    PETSC_SAMRAI_ERROR(ierr);
                }
                else
                {
                    ierr = VecScatterCreate(src[ln][i], src_IS[ln][depth],
                                            dst[ln][i], dst_IS[ln][depth],
                                            &scatter[ln][i]);
                    PETSC_SAMRAI_ERROR(ierr);
                }

                ierr = VecScatterBegin(scatter[ln][i], src[ln][i], dst[ln][i],
                                       INSERT_VALUES, SCATTER_FORWARD);
                PETSC_SAMRAI_ERROR(ierr);
            }
        }// if (d_level_contains_lag_data[ln])
    }// for (int ln = coarsest_ln; ln <= finest_ln; ++ln)

    // Complete the data scattering process, destroy the source Vec objects, and
    // distribute nonlocal data to the new configuration.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_level_contains_lag_data[ln])
        {
            std::map<std::string,SAMRAI::tbox::Pointer<LNodeLevelData> >& level_data =
                d_lag_quantity_data[ln];

            std::map<std::string,SAMRAI::tbox::Pointer<LNodeLevelData> >::iterator it;
            int i;
            for (it = level_data.begin(), i = 0; it != level_data.end();
                 ++it, ++i)
            {
                SAMRAI::tbox::Pointer<LNodeLevelData> data = (*it).second;

                ierr = VecScatterEnd(scatter[ln][i], src[ln][i], dst[ln][i],
                                     INSERT_VALUES, SCATTER_FORWARD);
                PETSC_SAMRAI_ERROR(ierr);
                ierr = VecScatterDestroy(scatter[ln][i]);
                PETSC_SAMRAI_ERROR(ierr);

                ierr = VecDestroy(src[ln][i]);
                PETSC_SAMRAI_ERROR(ierr);

                const int depth = data->getDepth();
                data->resetData(
                    dst[ln][i], d_nonlocal_petsc_indices[ln][depth]);
            }
            beginNonlocalDataFill(ln,ln);
        }// if (d_level_contains_lag_data[ln])
    }// for (int ln = coarsest_ln; ln <= finest_ln; ++ln)

    // Finish distributing nonlocal data to the new configuration.
    endNonlocalDataFill(coarsest_ln,finest_ln);

    // Indicate that the levels have been synchronized and destroy unneeded
    // ordering and indexing objects.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_needs_synch[ln] = false;

        if (d_ao[ln])
        {
            ierr = AODestroy(d_ao[ln]);
            PETSC_SAMRAI_ERROR(ierr);
        }
        d_ao[ln] = new_ao[ln];

        for (std::map<int,IS>::iterator it = src_IS[ln].begin();
             it != src_IS[ln].end(); ++it)
        {
            ierr = ISDestroy((*it).second);
            PETSC_SAMRAI_ERROR(ierr);
        }

        for (std::map<int,IS>::iterator it = dst_IS[ln].begin();
             it != dst_IS[ln].end(); ++it)
        {
            ierr = ISDestroy((*it).second);
            PETSC_SAMRAI_ERROR(ierr);
        }
    }

    // Restore the position data pointers for the LNodeIndex objects.
    restoreLocationPointers(coarsest_ln, finest_ln);

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
    t_end_data_redistribution->stop();
    return;
}// endDataRedistribution

void
LDataManager::updateWorkloadData(
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    t_update_workload_data->start();

    const int coarsest_ln =
        (coarsest_ln_in == -1)
        ? d_coarsest_ln
        : coarsest_ln_in;
    const int finest_ln =
        (finest_ln_in == -1)
        ? d_finest_ln
        : finest_ln_in;

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(coarsest_ln >= d_coarsest_ln &&
           coarsest_ln <= d_finest_ln);
    assert(finest_ln   >= d_coarsest_ln &&
           finest_ln   <= d_finest_ln);
#endif

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

            const SAMRAI::tbox::Pointer<LNodeIndexData2> idx_data =
                patch->getPatchData(d_lag_node_index_current_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > workload_data =
                patch->getPatchData(d_workload_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > node_count_data =
                patch->getPatchData(d_node_count_idx);

            workload_data->fillAll(d_alpha_work);
            node_count_data->fillAll(0.0);

            for (LNodeIndexData2::Iterator it(patch_box); it; it++)
            {
                const SAMRAI::pdat::CellIndex<NDIM>& i = *it;
                const LNodeIndexSet& node_set = (*idx_data)(i);
                (*node_count_data)(i) = node_set.size();
                (*workload_data)(i) += d_beta_work*(*node_count_data)(i);
            }
        }
    }

    t_update_workload_data->stop();
    return;
}// updateWorkloadData

void
LDataManager::updateIrregularCellData(
    const int stencil_size,
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    t_update_irregular_cell_data->start();

    const int coarsest_ln =
        (coarsest_ln_in == -1)
        ? d_coarsest_ln
        : coarsest_ln_in;
    const int finest_ln =
        (finest_ln_in == -1)
        ? d_finest_ln
        : finest_ln_in;

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(coarsest_ln >= d_coarsest_ln &&
           coarsest_ln <= d_finest_ln);
    assert(finest_ln   >= d_coarsest_ln &&
           finest_ln   <= d_finest_ln);
#endif

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::tbox::Pointer<LNodeIndexData2> idx_data =
                patch->getPatchData(d_lag_node_index_current_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > irregular_cell_data =
                patch->getPatchData(d_irregular_cell_idx);
            irregular_cell_data->fillAll(0.0);
            for (LNodeIndexData2::Iterator it(patch_box); it; it++)
            {
                const SAMRAI::pdat::CellIndex<NDIM>& i = *it;
                const LNodeIndexSet& node_set = (*idx_data)(i);
                if (!node_set.empty())
                {
                    const SAMRAI::hier::Box<NDIM> stencil_box =
                        SAMRAI::hier::Box<NDIM>::grow(SAMRAI::hier::Box<NDIM>(i,i),
                                                      SAMRAI::hier::IntVector<NDIM>(stencil_size/2));
                    irregular_cell_data->fillAll(1.0,stencil_box);
                }
            }
        }
    }

    t_update_irregular_cell_data->stop();
    return;
}// updateIrregularCellData

void
LDataManager::restoreLocationPointers(
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    t_restore_location_pointers->start();

    const int coarsest_ln =
        (coarsest_ln_in == -1)
        ? d_coarsest_ln
        : coarsest_ln_in;
    const int finest_ln =
        (finest_ln_in == -1)
        ? d_finest_ln
        : finest_ln_in;

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(coarsest_ln >= d_coarsest_ln &&
           coarsest_ln <= d_finest_ln);
    assert(finest_ln   >= d_coarsest_ln &&
           finest_ln   <= d_finest_ln);
#endif

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<LNodeIndexData2> idx_data =
                patch->getPatchData(d_lag_node_index_current_idx);
            const SAMRAI::hier::Box<NDIM>& ghost_box = idx_data->getGhostBox();

            for (LNodeIndexData2::Iterator it(ghost_box); it; it++)
            {
                const SAMRAI::pdat::CellIndex<NDIM>& i = *it;
                const LNodeIndexSet& node_set = (*idx_data)(i);
                if (!node_set.empty())
                {
                    for_each(node_set.begin(), node_set.end(),
                             RestoreLNodeIndexLocationPointers(
                                 d_lag_quantity_data[ln][POSN_DATA_NAME]->
                                 getLocalFormArray()));
                }
            }
        }
    }

    t_restore_location_pointers->stop();
    return;
}// restoreLocationPointers

void
LDataManager::invalidateLocationPointers(
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    t_invalidate_location_pointers->start();

    const int coarsest_ln =
        (coarsest_ln_in == -1)
        ? d_coarsest_ln
        : coarsest_ln_in;
    const int finest_ln =
        (finest_ln_in == -1)
        ? d_finest_ln
        : finest_ln_in;

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(coarsest_ln >= d_coarsest_ln &&
           coarsest_ln <= d_finest_ln);
    assert(finest_ln   >= d_coarsest_ln &&
           finest_ln   <= d_finest_ln);
#endif

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<LNodeIndexData2> idx_data =
                patch->getPatchData(d_lag_node_index_current_idx);
            const SAMRAI::hier::Box<NDIM>& ghost_box = idx_data->getGhostBox();

            for (LNodeIndexData2::Iterator it(ghost_box); it; it++)
            {
                const SAMRAI::pdat::CellIndex<NDIM>& i = *it;
                const LNodeIndexSet& node_set = (*idx_data)(i);
                if (!node_set.empty())
                {
                    for_each(node_set.begin(), node_set.end(),
                             InvalidateLNodeIndexLocationPointers());
                }
            }
        }
    }

    t_invalidate_location_pointers->stop();
    return;
}// invalidateLocationPointers

///
///  The following routines:
///
///      initializeLevelData(),
///      resetHierarchyConfiguration(),
///      applyGradientDetector()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::mesh::StandardTagAndInitStrategy<NDIM> abstract base class.
///

void
LDataManager::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
    const bool allocate_data)
{
    t_initialize_level_data->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((level_number >= 0)
           && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        assert(level_number == old_level->getLevelNumber());
    }
    assert(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

#ifdef DEBUG_CHECK_ASSERTIONS
    // Check for overlapping boxes on this level.
    //
    // (This is potentially fairly expensive and hence is only done when
    // assertion checking is active.)
    SAMRAI::hier::BoxList<NDIM> boxes(level->getBoxes());

    std::vector<bool> patch_overlaps(boxes.getNumberItems());
    std::vector<bool>::size_type j, k;

    for (k = 0; k < patch_overlaps.size(); ++k)
    {
        patch_overlaps[k] = false;
    }
    k = 0;
    while (!boxes.isEmpty())
    {
        j = k+1;
        SAMRAI::hier::Box<NDIM> tryme = boxes.getFirstItem();
        boxes.removeFirstItem();

        for (SAMRAI::hier::BoxList<NDIM>::Iterator ib(boxes); ib; ib++)
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
    if (allocate_data &&
        (level_number >= d_coarsest_ln) && (level_number <= d_finest_ln))
    {
        level->allocatePatchData(d_current_data, init_data_time);
    }
    else if (level_number >= d_coarsest_ln &&
             level_number <= d_finest_ln)
    {
        level->setTime(init_data_time, d_current_data);
    }

    // Fill data from the old level when available.
    if (!old_level.isNull())
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(old_level->getLevelNumber() == level_number);
#endif
        if (d_level_contains_lag_data[level_number])
        {
            level->allocatePatchData(d_scratch_data, init_data_time);

            d_lag_node_index_bdry_fill_alg->
                createSchedule(level, old_level)->fillData(init_data_time);

            level->deallocatePatchData(d_scratch_data);
        }
    }

    // Initialize the data on the level and, when appropriate, move data from
    // coarser levels to finer levels.
    if (initial_time)
    {
        d_level_contains_lag_data.resize(level_number+1);
        d_lag_quantity_data      .resize(level_number+1);
        d_needs_synch            .resize(level_number+1,false);
        d_ao                     .resize(level_number+1);
        d_num_nodes              .resize(level_number+1);
        d_node_offset            .resize(level_number+1);
        d_local_lag_indices      .resize(level_number+1);
        d_nonlocal_lag_indices   .resize(level_number+1);

#ifdef DEBUG_CHECK_ASSERTIONS
        assert(!d_lag_init.isNull());
#endif
        d_level_contains_lag_data[level_number] = d_lag_init->
            getLevelHasLagrangianData(level_number, can_be_refined);

        if (d_level_contains_lag_data[level_number])
        {
            // 1. Determine the number of local (on processor) nodes to be
            //    allocated on the patch level and allocate space for the local
            //    and non-local index data.
            const int num_local_nodes = d_lag_init->getLocalNodeCountOnPatchLevel(
                hierarchy, level_number,
                init_data_time, can_be_refined, initial_time);

            d_local_lag_indices  [level_number].resize(num_local_nodes,-1);
            d_local_petsc_indices[level_number].resize(num_local_nodes,-1);

            d_nonlocal_lag_indices  [level_number]   .clear();
            d_nonlocal_petsc_indices[level_number]   .clear();
            d_nonlocal_petsc_indices[level_number][1].clear();

            computeNodeOffsets(d_num_nodes[level_number],
                               d_node_offset[level_number], num_local_nodes);

            // 2. Allocate LNodeLevelData corresponding to the curvilinear mesh
            //    node positions and velocities.
            static const bool maintain_data = true;
            createLNodeLevelData(POSN_DATA_NAME, level_number, NDIM, maintain_data);
            createLNodeLevelData(VEL_DATA_NAME, level_number, NDIM, maintain_data);

            // 3. Initialize the Lagrangian data.

            // WARNING: If either of the following offsets is ever nonzero, note
            // that it may be necessary to modify IBHierarchyIntegrator, in
            // particular the code where the data related to the implementation
            // of the penalty IB method are initialized.
            static const int global_index_offset = 0;
            static const int local_index_offset = 0;

            const int num_initialized_local_nodes = d_lag_init->initializeDataOnPatchLevel(
                d_lag_node_index_current_idx,
                global_index_offset, local_index_offset,
                d_lag_quantity_data[level_number][POSN_DATA_NAME],
                d_lag_quantity_data[level_number][VEL_DATA_NAME],
                hierarchy, level_number,
                init_data_time, can_be_refined, initial_time, this);

            if (num_local_nodes != num_initialized_local_nodes)
            {
                TBOX_ERROR("LDataManager::initializeLevelData()"                             << "\n" <<
                           "  num_local_nodes             = " << num_local_nodes             << "\n" <<
                           "  num_initialized_local_nodes = " << num_initialized_local_nodes << "\n");
            }

            // 4. Compute the initial distribution (indexing) data.
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

                SAMRAI::tbox::Pointer<LNodeIndexData2> idx_data =
                    patch->getPatchData(d_lag_node_index_current_idx);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > workload_data =
                    patch->getPatchData(d_workload_idx);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > node_count_data =
                    patch->getPatchData(d_node_count_idx);
                const SAMRAI::hier::Box<NDIM>& ghost_box = idx_data->getGhostBox();

                idx_data->d_interior_local_indices.clear();
                idx_data->d_ghost_local_indices.clear();

                node_count_data->fillAll(0.0);

                for (LNodeIndexData2::Iterator it(ghost_box); it; it++)
                {
                    const SAMRAI::pdat::CellIndex<NDIM>& i = *it;
                    const LNodeIndexSet& node_set = (*idx_data)(i);
                    const bool patch_owns_idx = patch_box.contains(i);

                    if (patch_owns_idx)
                    {
                        (*node_count_data)(i) = node_set.size();
                        (*workload_data)(i) = d_alpha_work + d_beta_work*(*node_count_data)(i);
                    }

                    for (LNodeIndexSet::const_iterator n = node_set.begin();
                         n != node_set.end(); ++n)
                    {
                        const LNodeIndexSet::value_type& node_idx = *n;
                        const int lag_idx   = node_idx->getLagrangianIndex();
                        const int local_idx = node_idx->getLocalPETScIndex();
                        if (!(0 <= local_idx && local_idx < num_local_nodes))
                        {
                            TBOX_ERROR("LDataManager::initializeLevelData()"     << "\n" <<
                                       "  local_idx       = " << local_idx       << "\n" <<
                                       "  num_local_nodes = " << num_local_nodes << "\n");
                        }
                        d_local_lag_indices  [level_number][local_idx] = lag_idx;
                        d_local_petsc_indices[level_number][local_idx] = local_idx + d_node_offset[level_number];

                        if (patch_owns_idx)
                        {
                            idx_data->d_interior_local_indices.push_back(local_idx);
                        }
                        else
                        {
                            idx_data->d_ghost_local_indices.push_back(local_idx);
                        }
                    }
                }
            }

            // 5. The AO (application order) is determined by the initial values
            //    of the local Lagrangian indices.
            int ierr;

            if (d_ao[level_number])
            {
                ierr = AODestroy(d_ao[level_number]);
                PETSC_SAMRAI_ERROR(ierr);
            }

            ierr = AOCreateBasic(PETSC_COMM_WORLD, num_local_nodes,
                                 &d_local_lag_indices  [level_number][0],
                                 &d_local_petsc_indices[level_number][0],
                                 &d_ao[level_number]);
            PETSC_SAMRAI_ERROR(ierr);

            // 6. Restore the position data pointers for the LNodeIndex objects.
            restoreLocationPointers(level_number, level_number);

            // 7. Indicate to use nonuniform workload estimate on the specified
            //    level.
            if (!d_load_balancer.isNull())
            {
                d_load_balancer->setWorkloadPatchDataIndex(
                    d_workload_idx, level_number);
            }
        }
        else
        {
            // Indicate to use uniform workload estimate on the specified level.
            if (!d_load_balancer.isNull())
            {
                d_load_balancer->setUniformWorkload(level_number);
            }
        }
    }

    // Update the MPI process mapping data.
    const int mpi_process = SAMRAI::tbox::MPI::getRank();
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > mpi_proc_data =
            patch->getPatchData(d_mpi_proc_idx);
        mpi_proc_data->fillAll(mpi_process);
    }

    // If a Silo data writer is registered with the manager, give it access to
    // the new application ordering.
    if (!d_silo_writer.isNull() && d_level_contains_lag_data[level_number])
    {
        d_silo_writer->registerCoordsData(
            d_lag_quantity_data[level_number][POSN_DATA_NAME], level_number);
        d_silo_writer->registerLagrangianAO(
            d_ao[level_number], level_number);
    }
#if (NDIM == 3)
    // If a myocardial3D data writer is registered with the manager, give it
    // access to the new application ordering.
    if (!d_m3D_writer.isNull() && d_level_contains_lag_data[level_number])
    {
        d_m3D_writer->registerCoordsData(
            d_lag_quantity_data[level_number][POSN_DATA_NAME], level_number);
        d_m3D_writer->registerLagrangianAO(
            d_ao[level_number], level_number);
    }
#endif
    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
LDataManager::resetHierarchyConfiguration(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    t_reset_hierarchy_configuration->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((coarsest_level >= 0)
           && (coarsest_level <= finest_level)
           && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        assert(!(hierarchy->getPatchLevel(ln)).isNull());
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
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            if (d_level_contains_lag_data[ln])
            {
                d_silo_writer->registerCoordsData(
                    d_lag_quantity_data[ln][POSN_DATA_NAME], ln);
            }
        }
    }
#if (NDIM == 3)
    // Reset the myocardial3D data writer.
    if (!d_m3D_writer.isNull())
    {
        d_m3D_writer->setPatchHierarchy(hierarchy);
        d_m3D_writer->resetLevels(d_coarsest_ln, d_finest_ln);
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            if (d_level_contains_lag_data[ln])
            {
                d_m3D_writer->registerCoordsData(
                    d_lag_quantity_data[ln][POSN_DATA_NAME], ln);
            }
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
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        d_lag_node_index_bdry_fill_scheds[ln] =
            d_lag_node_index_bdry_fill_alg->createSchedule(level);
    }

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (int ln = SAMRAI::tbox::Utilities::imax(coarsest_level,1);
         ln <= finest_hier_level; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarser_level =
            hierarchy->getPatchLevel(ln-1);

        d_node_count_coarsen_scheds[ln] = d_node_count_coarsen_alg->
            createSchedule(coarser_level, level);
    }

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
LDataManager::applyGradientDetector(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    t_apply_gradient_detector->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((level_number >= 0)
           && (level_number <= hierarchy->getFinestLevelNumber()));
    assert(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    if (initial_time)
    {
        // Tag cells for refinement based on the initial configuration of the
        // Lagrangian structure.
        d_lag_init->tagCellsForInitialRefinement(
            hierarchy, level_number, error_data_time, tag_index);
    }
    else if (hierarchy->finerLevelExists(level_number))
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
            hierarchy->getPatchLevel(level_number);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > finer_level =
            hierarchy->getPatchLevel(level_number+1);

        // Zero out the node count data on the current level.
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > node_count_data =
                patch->getPatchData(d_node_count_idx);
            node_count_data->fillAll(0.0);
        }

        // Compute the workload and node count data on the next finer level of
        // the patch hierarchy.
        updateWorkloadData(level_number+1);

        // Coarsen the node count data from the next finer level of the patch
        // hierarchy.
        d_node_count_coarsen_scheds[level_number+1]->coarsenData();

        // Tag cells for refinement wherever there exist nodes on the the next
        // finer level of the Cartesian grid.
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > tag_data =
                patch->getPatchData(tag_index);
            const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > node_count_data =
                patch->getPatchData(d_node_count_idx);

            for (SAMRAI::pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
            {
                const SAMRAI::pdat::CellIndex<NDIM>& i = ic();
                if (!SAMRAI::tbox::Utilities::deq((*node_count_data)(i),0.0))
                {
                    (*tag_data)(i) = 1;
                }
            }
        }

        // Compute the workload and node count data on the present level of the
        // patch hierarchy (since it was invalidated above).
        updateWorkloadData(level_number);
    }

    t_apply_gradient_detector->stop();
    return;
}// applyGradientDetector

void
LDataManager::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif
    db->putInteger("LDATA_MANAGER_VERSION", LDATA_MANAGER_VERSION);

    db->putInteger("d_coarsest_ln", d_coarsest_ln);
    db->putInteger("d_finest_ln"  , d_finest_ln  );
    db->putDouble ("d_alpha_work" , d_alpha_work );
    db->putDouble ("d_beta_work"  , d_beta_work  );

    // Write out data that is stored on a level-by-level basis.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        std::ostringstream stream;
        stream << "level_" << ln;
        const std::string level_db_name = stream.str();
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> level_db = db->putDatabase(level_db_name);

        level_db->putBool("d_level_contains_lag_data", d_level_contains_lag_data[ln]);

        if (d_level_contains_lag_data[ln])
        {
            std::vector<std::string> ldata_names;
            for (std::map<std::string,SAMRAI::tbox::Pointer<LNodeLevelData> >::iterator it =
                     d_lag_quantity_data[ln].begin();
                 it != d_lag_quantity_data[ln].end(); ++it)
            {
                ldata_names.push_back((*it).first);
                (*it).second->putToDatabase(level_db->putDatabase(ldata_names.back()));
            }
            level_db->putInteger("n_ldata_names", ldata_names.size());
            level_db->putStringArray("ldata_names",
                                     &ldata_names[0], ldata_names.size());

            level_db->putInteger("d_num_nodes"  , d_num_nodes  [ln]);
            level_db->putInteger("d_node_offset", d_node_offset[ln]);

            level_db->putInteger("n_local_lag_indices",
                                 d_local_lag_indices[ln].size());
            if (!d_local_lag_indices[ln].empty())
            {
                level_db->putIntegerArray("d_local_lag_indices",
                                          &d_local_lag_indices[ln][0],
                                          d_local_lag_indices [ln].size());
            }
            level_db->putInteger("n_nonlocal_lag_indices",
                                 d_nonlocal_lag_indices[ln].size());
            if (!d_nonlocal_lag_indices[ln].empty())
            {
                level_db->putIntegerArray("d_nonlocal_lag_indices",
                                          &d_nonlocal_lag_indices[ln][0],
                                          d_nonlocal_lag_indices [ln].size());
            }
            level_db->putInteger("n_local_petsc_indices",
                                 d_local_petsc_indices[ln].size());
            if (!d_local_petsc_indices[ln].empty())
            {
                level_db->putIntegerArray("d_local_petsc_indices",
                                          &d_local_petsc_indices[ln][0],
                                          d_local_petsc_indices [ln].size());
            }
            // NOTE: d_nonlocal_petsc_indices[ln] is a map from the data depth
            // to the nonlocal petsc indices for that particular depth.  We only
            // serialize the indices corresponding to a data depth of 1.
            level_db->putInteger("n_nonlocal_petsc_indices",
                                 d_nonlocal_petsc_indices[ln][1].size());
            if (!d_nonlocal_petsc_indices[ln][1].empty())
            {
                level_db->putIntegerArray("d_nonlocal_petsc_indices",
                                          &d_nonlocal_petsc_indices[ln][1][0],
                                          d_nonlocal_petsc_indices [ln][1].size());
            }
        }
    }

    t_put_to_database->stop();
    return;
}// putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

LDataManager::LDataManager(
    const std::string& object_name,
    const SAMRAI::hier::IntVector<NDIM>& ghosts,
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
      d_irregular_cell_var(NULL),
      d_irregular_cell_idx(-1),
      d_mpi_proc_var(NULL),
      d_mpi_proc_idx(-1),
      d_output_mpi_proc(false),
      d_ghosts(ghosts),
      d_lag_node_index_bdry_fill_alg(NULL),
      d_lag_node_index_bdry_fill_scheds(),
      d_node_count_coarsen_alg(NULL),
      d_node_count_coarsen_scheds(),
      d_current_context(NULL),
      d_scratch_context(NULL),
      d_current_data(),
      d_scratch_data(),
      d_lag_quantity_data(),
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
    assert(!object_name.empty());
    assert(ghosts.min() >= 0);
#endif

    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->
            registerRestartItem(d_object_name, this);
    }

    bool from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }

    // Create/lookup the variable contexts.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    d_current_context = var_db->getContext(d_object_name+"::CURRENT");
    d_scratch_context = var_db->getContext(d_object_name+"::SCRATCH");

    // Register the SAMRAI variables with the VariableDatabase.
    d_lag_node_index_var = new LNodeIndexVariable2(
        d_object_name+"::lag_node_index");

    // Setup the current context.
    d_lag_node_index_current_idx = var_db->registerVariableAndContext(
        d_lag_node_index_var, d_current_context, d_ghosts);
    d_current_data.setFlag(d_lag_node_index_current_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_lag_node_index_current_idx);
    }

    // Setup the scratch context.
    d_lag_node_index_scratch_idx = var_db->registerVariableAndContext(
        d_lag_node_index_var, d_scratch_context, d_ghosts);
    d_scratch_data.setFlag(d_lag_node_index_scratch_idx);

    // Setup a refine algorithm, used to fill LNodeIndex boundary data.
    d_lag_node_index_bdry_fill_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    d_lag_node_index_bdry_fill_alg->registerRefine(
        d_lag_node_index_current_idx, // destination
        d_lag_node_index_current_idx, // source
        d_lag_node_index_scratch_idx, // temporary work space
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> >(NULL));

    // Register the workload variable with the VariableDatabase.
    d_workload_var = new SAMRAI::pdat::CellVariable<NDIM,double>(
        d_object_name+"::workload");

    d_workload_idx = var_db->registerVariableAndContext(
        d_workload_var, d_current_context, 0);
    d_current_data.setFlag(d_workload_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_workload_idx);
    }

    // Register the node count variable with the VariableDatabase.
    d_node_count_var = new SAMRAI::pdat::CellVariable<NDIM,double>(
        d_object_name+"::node_count");

    d_node_count_idx = var_db->registerVariableAndContext(
        d_node_count_var, d_current_context, 0);
    d_current_data.setFlag(d_node_count_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_node_count_idx);
    }

    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > coarsen_op =
        new SAMRAI::geom::CartesianCellDoubleWeightedAverage<NDIM>();

    d_node_count_coarsen_alg = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    d_node_count_coarsen_alg->
        registerCoarsen(d_node_count_idx, // destination
                        d_node_count_idx, // source
                        coarsen_op);

    // Register the irregular cell variable with the VariableDatabase.
    d_irregular_cell_var = new SAMRAI::pdat::CellVariable<NDIM,double>(
        d_object_name+"::irregular_cell");

    d_irregular_cell_idx = var_db->registerVariableAndContext(
        d_irregular_cell_var, d_current_context, 0);
    d_current_data.setFlag(d_irregular_cell_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_irregular_cell_idx);
    }

    // Register the MPI process mapping variable with the VariableDatabase.
    d_mpi_proc_var = new SAMRAI::pdat::CellVariable<NDIM,int>(
        d_object_name+"::mpi_proc");

    d_mpi_proc_idx = var_db->registerVariableAndContext(
        d_mpi_proc_var, d_current_context, 0);
    d_current_data.setFlag(d_mpi_proc_idx);

    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_mpi_proc_idx);
    }

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_map_lagrangian_to_petsc = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::mapLagrangianToPETSc()");
        t_map_petsc_to_lagrangian = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::mapPETScToLagrangian()");
        t_begin_data_redistribution = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::beginDataRedistribution()");
        t_end_data_redistribution = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::endDataRedistribution()");
        t_update_workload_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::updateWorkloadData()");
        t_update_irregular_cell_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::updateIrregularCellData()");
        t_restore_location_pointers = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::restoreLocationPointers()");
        t_invalidate_location_pointers = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::invalidateLocationPointers()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::initializeLevelData()");
        t_reset_hierarchy_configuration = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::resetHierarchyConfiguration()");
        t_apply_gradient_detector = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::applyGradientDetector()");
        t_put_to_database = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::putToDatabase()");
        t_begin_nonlocal_data_fill = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::beginNonlocalDataFill()");
        t_end_nonlocal_data_fill = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::endNonlocalDataFill()");
        t_compute_node_distribution = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::computeNodeDistribution()");
        t_compute_node_offsets = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LDataManager::computeNodeOffsets()");
        LEInteractor::initializeTimers();
        timers_need_init = false;
    }
    return;
}// LDataManager

LDataManager::~LDataManager()
{
    // Destroy any remaining AO objects.
    int ierr;
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        if (d_ao[ln])
        {
            ierr = AODestroy(d_ao[ln]);
            PETSC_SAMRAI_ERROR(ierr);
        }
    }
    return;
}// ~LDataManager

/////////////////////////////// PRIVATE //////////////////////////////////////

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
    ierr = VecGetSize(petsc_vec, &petsc_size);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecGetSize(lagrangian_vec, &lagrangian_size);  PETSC_SAMRAI_ERROR(ierr);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(petsc_size == lagrangian_size);
#endif
    int petsc_bs, lagrangian_bs;
    ierr = VecGetBlockSize(petsc_vec, &petsc_bs);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecGetBlockSize(lagrangian_vec, &lagrangian_bs);  PETSC_SAMRAI_ERROR(ierr);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(petsc_bs == lagrangian_bs);
#endif
    const int depth = petsc_bs;

    int lagrangian_lo, lagrangian_hi;
    ierr = VecGetOwnershipRange(lagrangian_vec, &lagrangian_lo, &lagrangian_hi);  PETSC_SAMRAI_ERROR(ierr);

    // Setup the IS data used to generate the VecScatters that redistribute the
    // data into the Lagrangian indexing scheme.
    std::vector<int> ref_is_idxs(lagrangian_hi-lagrangian_lo);
    for (size_t k = 0; k < ref_is_idxs.size(); ++k)
    {
        ref_is_idxs[k] = lagrangian_lo+k;
    }

    // Map Lagrangian indices to PETSc indices.
    std::vector<int> ao_dummy(1,-1);
    ierr = AOApplicationToPetsc(
        d_ao[level_number],
        (!ref_is_idxs.empty() ? static_cast<int>(ref_is_idxs.size()) : static_cast<int>(ao_dummy.size())),
        (!ref_is_idxs.empty() ? &ref_is_idxs[0]                      : &ao_dummy[0]));  PETSC_SAMRAI_ERROR(ierr);

    // Create the VecScatter object.
    std::vector<int> petsc_is_idxs(ref_is_idxs.size());
    std::transform(ref_is_idxs.begin(), ref_is_idxs.end(),
                   petsc_is_idxs.begin(),
                   std::bind2nd(std::multiplies<int>(),depth));

    IS petsc_is;
    ierr = ISCreateBlock(PETSC_COMM_WORLD, depth, petsc_is_idxs.size(),
                         &petsc_is_idxs[0], &petsc_is);  PETSC_SAMRAI_ERROR(ierr);

    VecScatter vec_scatter;
    ierr = VecScatterCreate(petsc_vec, petsc_is, lagrangian_vec,
                            PETSC_NULL, &vec_scatter);  PETSC_SAMRAI_ERROR(ierr);

    // Scatter the values.
    ierr = VecScatterBegin(vec_scatter, petsc_vec, lagrangian_vec,
                           INSERT_VALUES, mode);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecScatterEnd(vec_scatter, petsc_vec, lagrangian_vec,
                         INSERT_VALUES, mode);  PETSC_SAMRAI_ERROR(ierr);

    // Cleanup allocated data.
    ierr = ISDestroy(petsc_is);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecScatterDestroy(vec_scatter);  PETSC_SAMRAI_ERROR(ierr);

    return;
}// scatterData

void
LDataManager::beginNonlocalDataFill(
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    t_begin_nonlocal_data_fill->start();

    const int coarsest_ln =
        (coarsest_ln_in == -1)
        ? d_coarsest_ln
        : coarsest_ln_in;
    const int finest_ln =
        (finest_ln_in == -1)
        ? d_finest_ln
        : finest_ln_in;

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(coarsest_ln >= d_coarsest_ln &&
           coarsest_ln <= d_finest_ln);
    assert(finest_ln   >= d_coarsest_ln &&
           finest_ln   <= d_finest_ln);
#endif

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        std::map<std::string,SAMRAI::tbox::Pointer<LNodeLevelData> >& level_data =
            d_lag_quantity_data[ln];
        for_each(level_data.begin(), level_data.end(),
                 BeginLNodeLevelDataNonlocalFill());
    }

    t_begin_nonlocal_data_fill->stop();
    return;
}// beginNonlocalDataFill

void
LDataManager::endNonlocalDataFill(
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    t_end_nonlocal_data_fill->start();

    const int coarsest_ln =
        (coarsest_ln_in == -1)
        ? d_coarsest_ln
        : coarsest_ln_in;
    const int finest_ln =
        (finest_ln_in == -1)
        ? d_finest_ln
        : finest_ln_in;

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(coarsest_ln >= d_coarsest_ln &&
           coarsest_ln <= d_finest_ln);
    assert(finest_ln   >= d_coarsest_ln &&
           finest_ln   <= d_finest_ln);
#endif

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        std::map<std::string,SAMRAI::tbox::Pointer<LNodeLevelData> >& level_data =
            d_lag_quantity_data[ln];
        for_each(level_data.begin(), level_data.end(),
                 EndLNodeLevelDataNonlocalFill());
    }

    t_end_nonlocal_data_fill->stop();
    return;
}// endNonlocalDataFill

int
LDataManager::computeNodeDistribution(
    std::vector<int>& local_lag_indices,
    std::vector<int>& nonlocal_lag_indices,
    AO& ao,
    std::vector<int>& local_petsc_indices,
    std::vector<int>& nonlocal_petsc_indices,
    int& num_nodes,
    int& node_offset,
    std::map<int,std::vector<int>*>& patch_interior_local_indices,
    std::map<int,std::vector<int>*>& patch_ghost_local_indices,
    const int ln)
{
    t_compute_node_distribution->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(ln >= d_coarsest_ln &&
           ln <= d_finest_ln);
#endif

    typedef std::map<SAMRAI::pdat::CellIndex<NDIM>,const LNodeIndexSet*,CellIndexFortranOrder> IndexSetMap;

    // Collect the Lagrangian IDs of all of the Lagrangian nodes on the
    // specified level of the patch hierarchy.
    //
    // We differentiate between nodes which are local to the processor
    // (i.e. nodes which live in the interior of a patch owned by the processor)
    // and nodes which are non-local (i.e. nodes which live in the interior of a
    // patch owned by a different processor).
    //
    // It is important to emphasize that while a local node by definition lives
    // on the interior of some patch on this processor, it may also live in the
    // ghost cell regions of other patches owned by this processor.
    //
    // Non-local nodes ONLY appear in ghost cells for on processor patches.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

    const SAMRAI::hier::ProcessorMapping& proc_mapping = level->getProcessorMapping();
    const SAMRAI::hier::BoxArray<NDIM>& global_boxes = level->getBoxes();

    const int num_local_boxes = proc_mapping.getNumberOfLocalIndices();
    const SAMRAI::tbox::Array<int>& local_box_indices = proc_mapping.getLocalIndices();

    SAMRAI::hier::BoxArray<NDIM> local_boxes(num_local_boxes);

    for (int i = 0; i < num_local_boxes; ++i)
    {
        local_boxes(i) = global_boxes(local_box_indices[i]);
    }

    // Wipe out the old level and patch index data.
    local_lag_indices.clear();
    nonlocal_lag_indices.clear();

    // Keep track of the number of local indices we've encountered so far.
    int local_offset = 0;

    // Set all patch interior local indices on the level.
    IndexSetMap ghost_cell_local_map;
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::tbox::Pointer<LNodeIndexData2> lag_node_index_data =
            patch->getPatchData(d_lag_node_index_current_idx);
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const int& patch_num = patch->getPatchNumber();

        std::vector<int>& patch_interior_indices = *patch_interior_local_indices[patch_num];
        patch_interior_indices.clear();

        for (LNodeIndexData2::Iterator it(patch_box); it; it++)
        {
            const SAMRAI::pdat::CellIndex<NDIM>& cell_idx = *it;
            const LNodeIndexSet* const id_set = &((*lag_node_index_data)(cell_idx));
            const LNodeIndexSet::size_type& num_ids = id_set->size();
            if (num_ids > 0)
            {
                // All nodes located in this cell are owned by the patch and
                // consequently by the processor.
                //
                // We can immediately determine the local PETSc index for such
                // nodes.
                std::vector<int> cell_lag_ids(num_ids);
                std::transform(id_set->begin(), id_set->end(), cell_lag_ids.begin(),
                               GetLagrangianIndex());

                local_lag_indices.reserve(
                    local_lag_indices.size()+num_ids);
                local_lag_indices.insert(
                    local_lag_indices.end(),
                    cell_lag_ids.begin(), cell_lag_ids.end());

                for_each(id_set->begin(), id_set->end(),
                         SetLocalPETScIndex(local_offset));

                patch_interior_indices.resize(
                    patch_interior_indices.size()+num_ids);
                generate(patch_interior_indices.end()-num_ids,
                         patch_interior_indices.end(),
                         SetLocalPETScIndex(local_offset));

                local_offset += static_cast<int>(num_ids);

                ghost_cell_local_map[cell_idx] = id_set;
            }
        }
    }

    // Get the periodic shift vector and the refined domain box.
    const SAMRAI::hier::IntVector<NDIM>& ratio = level->getRatio();
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift = d_grid_geom->getPeriodicShift(ratio);
    if (periodic_shift != SAMRAI::hier::IntVector<NDIM>(0)) assert(d_grid_geom->getDomainIsSingleBox());
    const SAMRAI::hier::Box<NDIM> domain_box =
        SAMRAI::hier::Box<NDIM>::refine(d_grid_geom->getPhysicalDomain()(0),ratio);

    // Set all remaining local and nonlocal indices on the level.
    IndexSetMap ghost_cell_nonlocal_map;
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::tbox::Pointer<LNodeIndexData2> lag_node_index_data =
            patch->getPatchData(d_lag_node_index_current_idx);
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const SAMRAI::hier::Box<NDIM>& ghost_box = lag_node_index_data->getGhostBox();
        const int& patch_num = patch->getPatchNumber();

        std::vector<int>& patch_ghost_indices = *patch_ghost_local_indices[patch_num];
        patch_ghost_indices.clear();

        SAMRAI::hier::BoxList<NDIM> ghost_boxes;
        ghost_boxes.removeIntersections(ghost_box,patch_box);

        for (SAMRAI::hier::BoxList<NDIM>::Iterator b(ghost_boxes); b; b++)
        {
            for (LNodeIndexData2::Iterator it(*b); it; it++)
            {
                const SAMRAI::pdat::CellIndex<NDIM>& cell_idx = *it;
                const LNodeIndexSet* const id_set = &((*lag_node_index_data)(cell_idx));
                const LNodeIndexSet::size_type& num_ids = id_set->size();
                if (num_ids > 0)
                {
                    const SAMRAI::pdat::CellIndex<NDIM> canonical_cell_idx =
                        get_canonical_cell_index(cell_idx, domain_box, periodic_shift);
                    if (local_boxes.contains(canonical_cell_idx))
                    {
#ifdef DEBUG_CHECK_ASSERTIONS
                        assert(ghost_cell_local_map.count(canonical_cell_idx) == 1);
#endif
                        // The nodes are local nodes, so we just look-up their
                        // local IDs.
                        for_each(id_set->begin(), id_set->end(),
                                 GetLocalPETScIndexFromIDSet(
                                     (ghost_cell_local_map[canonical_cell_idx])->begin()));

                        patch_ghost_indices.resize(
                            patch_ghost_indices.size()+num_ids);
                        generate(patch_ghost_indices.end()-num_ids,
                                 patch_ghost_indices.end(),
                                 GetLocalPETScIndexFromIDSet(
                                     (ghost_cell_local_map[canonical_cell_idx])->begin()));
                    }
                    else
                    {
                        // The nodes are not local to the processor, so we must
                        // either assign or lookup their local IDs.
                        if (ghost_cell_nonlocal_map.count(canonical_cell_idx) == 0)
                        {
                            // We have not set the local IDs for this cell
                            // index, so we must assign them.
                            std::vector<int> cell_lag_ids(num_ids);
                            std::transform(id_set->begin(), id_set->end(), cell_lag_ids.begin(),
                                           GetLagrangianIndex());

                            nonlocal_lag_indices.reserve(
                                nonlocal_lag_indices.size()+num_ids);
                            nonlocal_lag_indices.insert(
                                nonlocal_lag_indices.end(),
                                cell_lag_ids.begin(), cell_lag_ids.end());

                            for_each(id_set->begin(), id_set->end(),
                                     SetLocalPETScIndex(local_offset));

                            local_offset += static_cast<int>(num_ids);

                            ghost_cell_nonlocal_map[canonical_cell_idx] = id_set;

                            patch_ghost_indices.resize(
                                patch_ghost_indices.size()+num_ids);
                            generate(patch_ghost_indices.end()-num_ids,
                                     patch_ghost_indices.end(),
                                     GetLocalPETScIndexFromIDSet(id_set->begin()));
                        }
                        else
                        {
                            // We have already set the local IDs for this cell
                            // index, so we just look-up their local IDs.
                            for_each(id_set->begin(), id_set->end(),
                                     GetLocalPETScIndexFromIDSet(
                                         (ghost_cell_nonlocal_map[canonical_cell_idx])->begin()));

                            patch_ghost_indices.resize(
                                patch_ghost_indices.size()+num_ids);
                            generate(patch_ghost_indices.end()-num_ids,
                                     patch_ghost_indices.end(),
                                     GetLocalPETScIndexFromIDSet(
                                         (ghost_cell_nonlocal_map[canonical_cell_idx])->begin()));
                        }
                    }
                }
            }
        }
    }

    // Trim-to-fit vectors to get rid of any excess capacity.
    if (local_lag_indices.size() != local_lag_indices.capacity())
    {
        std::vector<int>(local_lag_indices).swap(local_lag_indices);
    }
    if (nonlocal_lag_indices.size() != nonlocal_lag_indices.capacity())
    {
        std::vector<int>(nonlocal_lag_indices).swap(nonlocal_lag_indices);
    }

    // We now compute the new PETSc global ordering and initialize the AO
    // object.
    int ierr;

    // Determine how many nodes are on each processor to calculate the PETSc
    // indexing scheme.
    const int    num_local_nodes =    local_lag_indices.size();
    const int num_nonlocal_nodes = nonlocal_lag_indices.size();

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
    generate(local_petsc_indices.begin(),local_petsc_indices.end(),
             SetLocalPETScIndex(node_offset));

    if (ao)
    {
        ierr = AODestroy(ao);
        PETSC_SAMRAI_ERROR(ierr);
    }

    ierr = AOCreateBasic(PETSC_COMM_WORLD, num_local_nodes,
                         &node_indices[0], &local_petsc_indices[0], &ao);
    PETSC_SAMRAI_ERROR(ierr);

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
    PETSC_SAMRAI_ERROR(ierr);

//  If desired, we can build an ISLocalToGlobalMapping at this point.  It is not
//  necessary to do so, since all of the information this mapping provides is
//  encoded in d_node_offset and the two vectors d_local_lag_indices and
//  d_nonlocal_lag_indices.
//
//  However, here is the code that you would use to build such a mapping:
//
//  // Construct a local to global mapping for the new ordering.
//  ISLocalToGlobalMapping local_to_global_map;
//
//  if (local_to_global_map)
//  {
//      ierr = ISLocalToGlobalMappingDestroy(local_to_global_map);
//      PETSC_SAMRAI_ERROR(ierr);
//  }
//  ierr = ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, num_proc_nodes,
//                                      &node_indices[0],
//                                      &local_to_global_map);
//  PETSC_SAMRAI_ERROR(ierr);

    // Keep track of the global PETSc indices of the ghost nodes.
    nonlocal_petsc_indices.clear();
    nonlocal_petsc_indices.reserve(num_nonlocal_nodes);

    nonlocal_petsc_indices.insert(nonlocal_petsc_indices.end(),
                                  node_indices.begin()+num_local_nodes,
                                  node_indices.end());

    t_compute_node_distribution->stop();
    return 0;
}// computeNodeDistribution

void
LDataManager::computeNodeOffsets(
    int& num_nodes,
    int& node_offset,
    const int& num_local_nodes)
{
    t_compute_node_offsets->start();

    const int mpi_size = SAMRAI::tbox::MPI::getNodes();
    const int mpi_rank = SAMRAI::tbox::MPI::getRank();

    std::vector<int> num_nodes_proc(mpi_size,0);

    SAMRAI::tbox::MPI::allGather(num_local_nodes, &num_nodes_proc[0]);

    node_offset = std::accumulate(num_nodes_proc.begin(),
                                  num_nodes_proc.begin()+mpi_rank, 0);

    num_nodes = std::accumulate(num_nodes_proc.begin()+mpi_rank,
                                num_nodes_proc.end(), node_offset);

    t_compute_node_offsets->stop();
    return;
}// computeNodeOffsets

// bool
// LDataManager::checkAllPeriodicShifts(
//     SAMRAI::pdat::CellIndex<NDIM>& shifted_idx,
//     const SAMRAI::pdat::CellIndex<NDIM>& idx,
//     const SAMRAI::hier::BoxArray<NDIM>& boxes,
//     const SAMRAI::hier::IntVector<NDIM>& periodic_shift)
// {
// #if (NDIM > 3)
//     TBOX_ERROR("NDIM > 3 NOT SUPPORTED!\n");
// #endif
//     // Check all periodic shifts.
//     //
//     // NOTE: This could be made more efficient, but hopefully this
//     // check will occur infrequently.
//     bool boxes_contains_shifted_idx = false;
//
//     for (int d = 0; d < NDIM && !boxes_contains_shifted_idx; ++d)
//     {
//         if (!boxes_contains_shifted_idx)
//         {
//             shifted_idx = idx;
//             shifted_idx(d) += periodic_shift(d);
//             if (boxes.contains(shifted_idx))
//             {
//                 boxes_contains_shifted_idx = true;
//             }
//         }
//
//         if (!boxes_contains_shifted_idx)
//         {
//             shifted_idx = idx;
//             shifted_idx(d) -= periodic_shift(d);
//             if (boxes.contains(shifted_idx))
//             {
//                 boxes_contains_shifted_idx = true;
//             }
//         }
//     }
//
//     for (int d0 = 0; d0 < NDIM && !boxes_contains_shifted_idx; ++d0)
//     {
//         for (int d1 = 0; d1 < NDIM && !boxes_contains_shifted_idx; ++d1)
//         {
//             if (d1 == d0) break;
//
//             if (!boxes_contains_shifted_idx)
//             {
//                 shifted_idx = idx;
//                 shifted_idx(d0) += periodic_shift(d0);
//                 shifted_idx(d1) += periodic_shift(d1);
//                 if (boxes.contains(shifted_idx))
//                 {
//                     boxes_contains_shifted_idx = true;
//                 }
//             }
//
//             if (!boxes_contains_shifted_idx)
//             {
//                 shifted_idx = idx;
//                 shifted_idx(d0) += periodic_shift(d0);
//                 shifted_idx(d1) -= periodic_shift(d1);
//                 if (boxes.contains(shifted_idx))
//                 {
//                     boxes_contains_shifted_idx = true;
//                 }
//             }
//
//             if (!boxes_contains_shifted_idx)
//             {
//                 shifted_idx = idx;
//                 shifted_idx(d0) -= periodic_shift(d0);
//                 shifted_idx(d1) += periodic_shift(d1);
//                 if (boxes.contains(shifted_idx))
//                 {
//                     boxes_contains_shifted_idx = true;
//                 }
//             }
//
//             if (!boxes_contains_shifted_idx)
//             {
//                 shifted_idx = idx;
//                 shifted_idx(d0) -= periodic_shift(d0);
//                 shifted_idx(d1) -= periodic_shift(d1);
//                 if (boxes.contains(shifted_idx))
//                 {
//                     boxes_contains_shifted_idx = true;
//                 }
//             }
//         }
//     }
//
//     for (int d0 = 0; d0 < NDIM && !boxes_contains_shifted_idx; ++d0)
//     {
//         for (int d1 = 0; d1 < NDIM && !boxes_contains_shifted_idx; ++d1)
//         {
//             if (d1 == d0) break;
//             for (int d2 = 0; d2 < NDIM && !boxes_contains_shifted_idx; ++d2)
//             {
//                 if (d2 == d0) break;
//                 if (d2 == d1) break;
//
//                 if (!boxes_contains_shifted_idx)
//                 {
//                     shifted_idx = idx;
//                     shifted_idx(d0) += periodic_shift(d0);
//                     shifted_idx(d1) += periodic_shift(d1);
//                     shifted_idx(d2) += periodic_shift(d2);
//                     if (boxes.contains(shifted_idx))
//                     {
//                         boxes_contains_shifted_idx = true;
//                     }
//                 }
//
//                 if (!boxes_contains_shifted_idx)
//                 {
//                     shifted_idx = idx;
//                     shifted_idx(d0) += periodic_shift(d0);
//                     shifted_idx(d1) += periodic_shift(d1);
//                     shifted_idx(d2) -= periodic_shift(d2);
//                     if (boxes.contains(shifted_idx))
//                     {
//                         boxes_contains_shifted_idx = true;
//                     }
//                 }
//
//                 if (!boxes_contains_shifted_idx)
//                 {
//                     shifted_idx = idx;
//                     shifted_idx(d0) += periodic_shift(d0);
//                     shifted_idx(d1) -= periodic_shift(d1);
//                     shifted_idx(d2) += periodic_shift(d2);
//                     if (boxes.contains(shifted_idx))
//                     {
//                         boxes_contains_shifted_idx = true;
//                     }
//                 }
//
//                 if (!boxes_contains_shifted_idx)
//                 {
//                     shifted_idx = idx;
//                     shifted_idx(d0) += periodic_shift(d0);
//                     shifted_idx(d1) -= periodic_shift(d1);
//                     shifted_idx(d2) -= periodic_shift(d2);
//                     if (boxes.contains(shifted_idx))
//                     {
//                         boxes_contains_shifted_idx = true;
//                     }
//                 }
//
//                 if (!boxes_contains_shifted_idx)
//                 {
//                     shifted_idx = idx;
//                     shifted_idx(d0) -= periodic_shift(d0);
//                     shifted_idx(d1) += periodic_shift(d1);
//                     shifted_idx(d2) += periodic_shift(d2);
//                     if (boxes.contains(shifted_idx))
//                     {
//                         boxes_contains_shifted_idx = true;
//                     }
//                 }
//
//                 if (!boxes_contains_shifted_idx)
//                 {
//                     shifted_idx = idx;
//                     shifted_idx(d0) -= periodic_shift(d0);
//                     shifted_idx(d1) += periodic_shift(d1);
//                     shifted_idx(d2) -= periodic_shift(d2);
//                     if (boxes.contains(shifted_idx))
//                     {
//                         boxes_contains_shifted_idx = true;
//                     }
//                 }
//
//                 if (!boxes_contains_shifted_idx)
//                 {
//                     shifted_idx = idx;
//                     shifted_idx(d0) -= periodic_shift(d0);
//                     shifted_idx(d1) -= periodic_shift(d1);
//                     shifted_idx(d2) += periodic_shift(d2);
//                     if (boxes.contains(shifted_idx))
//                     {
//                         boxes_contains_shifted_idx = true;
//                     }
//                 }
//
//                 if (!boxes_contains_shifted_idx)
//                 {
//                     shifted_idx = idx;
//                     shifted_idx(d0) -= periodic_shift(d0);
//                     shifted_idx(d1) -= periodic_shift(d1);
//                     shifted_idx(d2) -= periodic_shift(d2);
//                     if (boxes.contains(shifted_idx))
//                     {
//                         boxes_contains_shifted_idx = true;
//                     }
//                 }
//             }
//         }
//     }
//     return boxes_contains_shifted_idx;
// }// checkAllPeriodicShifts

void
LDataManager::getFromRestart()
{
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> restart_db =
        SAMRAI::tbox::RestartManager::getManager()->getRootDatabase();

    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db;
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
    d_level_contains_lag_data.resize(d_finest_ln+1);
    d_lag_quantity_data      .resize(d_finest_ln+1);
    d_needs_synch            .resize(d_finest_ln+1,false);
    d_ao                     .resize(d_finest_ln+1);
    d_num_nodes              .resize(d_finest_ln+1);
    d_node_offset            .resize(d_finest_ln+1);
    d_local_lag_indices      .resize(d_finest_ln+1);
    d_nonlocal_lag_indices   .resize(d_finest_ln+1);
    d_local_petsc_indices    .resize(d_finest_ln+1);
    d_nonlocal_petsc_indices .resize(d_finest_ln+1);

    // Read in data that is stored on a level-by-level basis.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        std::ostringstream stream;
        stream << "level_" << ln;
        const std::string level_db_name = stream.str();
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> level_db = db->getDatabase(level_db_name);

        d_level_contains_lag_data[ln] = level_db->getBool("d_level_contains_lag_data");

        if (d_level_contains_lag_data[ln])
        {
            const int n_ldata_names = level_db->getInteger("n_ldata_names");
            std::vector<std::string> ldata_names(n_ldata_names);
            level_db->getStringArray("ldata_names", &ldata_names[0], n_ldata_names);

            std::set<int> data_depths;
            for (std::vector<std::string>::iterator it = ldata_names.begin();
                 it != ldata_names.end(); ++it)
            {
                const std::string& ldata_name = *it;
                d_lag_quantity_data[ln][ldata_name] =
                    new LNodeLevelData(level_db->getDatabase(ldata_name));
                data_depths.insert(
                    d_lag_quantity_data[ln][ldata_name]->getDepth());
            }

            d_num_nodes  [ln] = level_db->getInteger("d_num_nodes"  );
            d_node_offset[ln] = level_db->getInteger("d_node_offset");

            const int n_local_lag_indices = level_db->getInteger("n_local_lag_indices");
            if (n_local_lag_indices > 0)
            {
                d_local_lag_indices[ln].resize(n_local_lag_indices);
                level_db->getIntegerArray("d_local_lag_indices",
                                          &d_local_lag_indices[ln][0],
                                          n_local_lag_indices);
            }
            const int n_nonlocal_lag_indices = level_db->getInteger("n_nonlocal_lag_indices");
            if (n_nonlocal_lag_indices > 0)
            {
                d_nonlocal_lag_indices[ln].resize(n_nonlocal_lag_indices);
                level_db->getIntegerArray("d_nonlocal_lag_indices",
                                          &d_nonlocal_lag_indices[ln][0],
                                          n_nonlocal_lag_indices);
            }
            const int n_local_petsc_indices = level_db->getInteger("n_local_petsc_indices");
            if (n_local_petsc_indices > 0)
            {
                d_local_petsc_indices[ln].resize(n_local_petsc_indices);
                level_db->getIntegerArray("d_local_petsc_indices",
                                          &d_local_petsc_indices[ln][0],
                                          n_local_petsc_indices);
            }
            const int n_nonlocal_petsc_indices = level_db->getInteger("n_nonlocal_petsc_indices");
            if (n_nonlocal_petsc_indices > 0)
            {
                d_nonlocal_petsc_indices[ln][1].resize(n_nonlocal_petsc_indices);
                level_db->getIntegerArray("d_nonlocal_petsc_indices",
                                          &d_nonlocal_petsc_indices[ln][1][0],
                                          n_nonlocal_petsc_indices);

                // Rebuild the nonlocal PETSc indices for the other data depths.
                for (std::set<int>::const_iterator it = data_depths.begin();
                     it != data_depths.end(); ++it)
                {
                    const int depth = *it;

                    d_nonlocal_petsc_indices[ln][depth].
                        resize(d_nonlocal_petsc_indices[ln][1].size());

                    std::transform(d_nonlocal_petsc_indices[ln][    1].begin(),
                                   d_nonlocal_petsc_indices[ln][    1].end  (),
                                   d_nonlocal_petsc_indices[ln][depth].begin(),
                                   std::bind2nd(std::multiplies<int>(),depth));
                }
            }

            // Rebuild the application ordering.
            int ierr;
            ierr = AOCreateBasic(PETSC_COMM_WORLD, n_local_lag_indices,
                                 &d_local_lag_indices  [ln][0],
                                 &d_local_petsc_indices[ln][0], &d_ao[ln]);
            PETSC_SAMRAI_ERROR(ierr);
        }
    }

    return;
}// getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
