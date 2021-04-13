// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2021 by the IBAMR developers
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

#include "ibtk/FECache.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/FEMappingCache.h"
#include "ibtk/FEProjector.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/QuadratureCache.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/SAMRAIDataCache.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"

#include "BasePatchHierarchy.h"
#include "Box.h"
#include "CartesianCellDoubleWeightedAverage.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchData.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineSchedule.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "SideIterator.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

#include "libmesh/boundary_info.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_parallel_type.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_type.h"
#include "libmesh/id_types.h"
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_version.h"
#include "libmesh/linear_solver.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"
#include "libmesh/quadrature_grid.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_vector.h"
#include "libmesh/variant_filter_iterator.h"

#include "petscvec.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include "boost/multi_array.hpp"
IBTK_ENABLE_EXTRA_WARNINGS

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

namespace libMesh
{
namespace Parallel
{
class Communicator;
} // namespace Parallel
} // namespace libMesh

namespace libMesh
{
template <typename T>
class VectorValue;
} // namespace libMesh

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_reinit_element_mappings;
static Timer* t_build_ghosted_solution_vector;
static Timer* t_build_ghosted_vector;
static Timer* t_spread;
static Timer* t_prolong_data;
static Timer* t_interp;
static Timer* t_interp_weighted;
static Timer* t_restrict_data;
static Timer* t_update_workload_estimates;
static Timer* t_initialize_level_data;
static Timer* t_reset_hierarchy_configuration;
static Timer* t_apply_gradient_detector;
static Timer* t_put_to_database;

// Version of FEDataManager restart file data.
static const int FE_DATA_MANAGER_VERSION = 4;

// Local helper functions.
struct ElemComp
{
    inline bool operator()(const Elem* const x, const Elem* const y) const
    {
        return x->id() < y->id();
    } // operator()
};

template <class ContainerOfContainers>
inline void
collect_unique_elems(std::vector<Elem*>& elems, const ContainerOfContainers& elem_patch_map)
{
    std::set<Elem*, ElemComp> elem_set;
    for (auto it = elem_patch_map.begin(); it != elem_patch_map.end(); ++it)
    {
        elem_set.insert(it->begin(), it->end());
    }
    elems.assign(elem_set.begin(), elem_set.end());
    return;
} // collect_unique_elems

std::set<libMesh::subdomain_id_type>
collect_subdomain_ids(const libMesh::MeshBase& mesh)
{
    std::set<libMesh::subdomain_id_type> subdomain_ids;
    // Get all subdomain ids present, not just local elements.
    const auto el_begin = mesh.elements_begin();
    const auto el_end = mesh.elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        subdomain_ids.insert((*el_it)->subdomain_id());
    }
    return subdomain_ids;
}

/**
 * A difficulty with FEDataManager is that it needs to work with both volumetric
 * and surface meshes, even though for performance reasons we template FEMapping
 * on the current topological dimension of the mesh. Fortunately we do not need
 * the contravariants or other topological dimension-dependent values: we just
 * need the JxW values. Hence this little helper function lets one get the
 * correct JxW values from the correct mapping object based on the runtime
 * dimensionality of the mesh.
 */
const std::vector<double>&
get_JxW(quadrature_key_type key,
        const Elem* elem,
        const bool is_volume_mesh,
        FEMappingCache<NDIM, NDIM>& volume_mapping_cache,
        FEMappingCache<NDIM - 1, NDIM>& surface_mapping_cache)
{
    if (is_volume_mesh)
    {
        auto& mapping = volume_mapping_cache[key];
        mapping.reinit(elem);
        return mapping.getJxW();
    }
    else
    {
        auto& mapping = surface_mapping_cache[key];
        mapping.reinit(elem);
        return mapping.getJxW();
    }
}

#if LIBMESH_VERSION_LESS_THAN(1, 6, 0)
// libMesh's box intersection code is slow and not in a header (i.e., cannot be
// inlined). This is problematic for us since we presently call this function
// for every element on every patch: i.e., for N patches on the current
// processor and K *total* (i.e., on all processors) elements we call this
// function K*N times, which takes up a lot of time in regrids. For example:
// switching to this function lowers the time required to get to the end of the
// first time step in the TAVR model by 20%.
inline bool
bbox_intersects(const libMeshWrappers::BoundingBox& a, const libMeshWrappers::BoundingBox& b)
{
    const libMesh::Point& a_lower = a.first;
    const libMesh::Point& a_upper = a.second;

    const libMesh::Point& b_lower = b.first;
    const libMesh::Point& b_upper = b.second;

    // Since boxes are tensor products of line intervals it suffices to check
    // that the line segments for each coordinate axis overlap.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        // Line segments can intersect in two ways:
        // 1. They can overlap.
        // 2. One can be inside the other.
        //
        // In the first case we want to see if either end point of the second
        // line segment lies within the first. In the second case we can simply
        // check that one end point of the first line segment lies in the second
        // line segment. Note that we don't need, in the second case, to do two
        // checks since that case is already covered by the first.
        if (!((a_lower(d) <= b_lower(d) && b_lower(d) <= a_upper(d)) ||
              (a_lower(d) <= b_upper(d) && b_upper(d) <= a_upper(d))) &&
            !((b_lower(d) <= a_lower(d) && a_lower(d) <= b_upper(d))))
        {
            return false;
        }
    }

    return true;
}
#endif
} // namespace

FEData::FEData(std::string object_name, EquationSystems& equation_systems, const bool register_for_restart)
    : d_object_name(std::move(object_name)),
      d_registered_for_restart(register_for_restart),
      d_es(&equation_systems),
      d_quadrature_cache(d_es->get_mesh().mesh_dimension())
{
    TBOX_ASSERT(!d_object_name.empty());

    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        FEData::getFromRestart();
    }
} // FEData

FEData::~FEData()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
} // ~FEData

void
FEData::getFromRestart()
{
} // getFromRestart

void FEData::putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> /*db*/)
{
} // putToDatabase

void
FEData::setEquationSystems(EquationSystems* const equation_systems, const int /*level_number*/)
{
    d_es = equation_systems;

    // Now that we have the EquationSystems object we know the dimensionality
    // of the mesh.
    d_quadrature_cache = QuadratureCache(d_es->get_mesh().mesh_dimension());
    return;
} // setEquationSystems

EquationSystems*
FEData::getEquationSystems() const
{
    return d_es;
} // getEquationSystems

FEData::SystemDofMapCache*
FEData::getDofMapCache(const std::string& system_name)
{
    TBOX_ASSERT(d_es);
    return getDofMapCache(d_es->get_system(system_name).number());
}

FEData::SystemDofMapCache*
FEData::getDofMapCache(unsigned int system_num)
{
    System& system = d_es->get_system(system_num);
    const unsigned int n_vars = system.n_vars();
    TBOX_ASSERT(n_vars > 0);
    const FEType& fe_type = system.variable_type(0);

    // we assume that all variables have the same FEType inside IBAMR
    for (unsigned int var_n = 1; var_n < n_vars; ++var_n)
    {
        TBOX_ASSERT(fe_type == system.variable_type(var_n));
    }

    // We assume that the dofs have the same layout as long as they have the
    // same fe type and number of variables
    std::pair<unsigned int, FEType> key(n_vars, fe_type);

    std::unique_ptr<FEData::SystemDofMapCache>& dof_map_cache = d_system_dof_map_cache[key];
    if (dof_map_cache == nullptr)
    {
        dof_map_cache.reset(new FEData::SystemDofMapCache(system));
    }
    return dof_map_cache.get();
} // getDofMapCache

void
FEData::clearPatchHierarchyDependentData()
{
    d_system_dof_map_cache.clear();
    d_quadrature_cache.clear();
}

SubdomainToPatchLevelTranslation::SubdomainToPatchLevelTranslation(
    const int max_level_number,
    const std::set<libMesh::subdomain_id_type>& subdomain_ids,
    const Pointer<Database>& input_db)
    : d_max_level_number(max_level_number)
{
    std::fill(d_fixed.begin(), d_fixed.end(), IBTK::invalid_level_number);

    // Overwrite all mesh subdomain ids to be on the finest level. Any that are
    // specified in the input database will be corrected when we read it.
    for (const libMesh::subdomain_id_type id : subdomain_ids)
    {
        get(id) = max_level_number;
    }

    // check that no id is set more than once
    std::set<libMesh::subdomain_id_type> known_subdomain_ids;
    if (input_db)
    {
        for (int ln = 0; ln <= d_max_level_number; ++ln)
        {
            const std::string key = "level_" + std::to_string(ln);
            if (input_db->keyExists(key))
            {
                const Array<int> subdomain_ids = input_db->getIntegerArray(key);
                for (int i = 0; i < subdomain_ids.getSize(); ++i)
                {
                    const auto id = subdomain_ids[i];
                    get(id) = ln;

                    if (known_subdomain_ids.count(id) == 0)
                        known_subdomain_ids.insert(id);
                    else
                    {
                        TBOX_ERROR("The input database lists the subdomain id "
                                   << id
                                   << " more than once. Since each subdomain can only "
                                      "be on exactly one level this is not permitted.");
                    }
                }
            }
        }
    }
}

const boundary_id_type FEDataManager::ZERO_DISPLACEMENT_X_BDRY_ID = 0x100;
const boundary_id_type FEDataManager::ZERO_DISPLACEMENT_Y_BDRY_ID = 0x200;
const boundary_id_type FEDataManager::ZERO_DISPLACEMENT_Z_BDRY_ID = 0x400;
const boundary_id_type FEDataManager::ZERO_DISPLACEMENT_XY_BDRY_ID = 0x100 | 0x200;
const boundary_id_type FEDataManager::ZERO_DISPLACEMENT_XZ_BDRY_ID = 0x100 | 0x400;
const boundary_id_type FEDataManager::ZERO_DISPLACEMENT_YZ_BDRY_ID = 0x200 | 0x400;
const boundary_id_type FEDataManager::ZERO_DISPLACEMENT_XYZ_BDRY_ID = 0x100 | 0x200 | 0x400;
std::map<std::string, FEDataManager*> FEDataManager::s_data_manager_instances;
bool FEDataManager::s_registered_callback = false;
unsigned char FEDataManager::s_shutdown_priority = 200;

FEDataManager*
FEDataManager::getManager(std::shared_ptr<FEData> fe_data,
                          const std::string& name,
                          const Pointer<Database>& input_db,
                          const int max_levels,
                          const FEDataManager::InterpSpec& default_interp_spec,
                          const FEDataManager::SpreadSpec& default_spread_spec,
                          const FEDataManager::WorkloadSpec& default_workload_spec,
                          const IntVector<NDIM>& min_ghost_width,
                          std::shared_ptr<SAMRAIDataCache> eulerian_data_cache,
                          bool register_for_restart)
{
    if (s_data_manager_instances.find(name) == s_data_manager_instances.end())
    {
        const IntVector<NDIM> ghost_width = IntVector<NDIM>::max(
            min_ghost_width,
            IntVector<NDIM>(std::max(LEInteractor::getMinimumGhostWidth(default_interp_spec.kernel_fcn),
                                     LEInteractor::getMinimumGhostWidth(default_spread_spec.kernel_fcn))));
        s_data_manager_instances[name] = new FEDataManager(name,
                                                           input_db,
                                                           max_levels,
                                                           default_interp_spec,
                                                           default_spread_spec,
                                                           default_workload_spec,
                                                           ghost_width,
                                                           eulerian_data_cache,
                                                           fe_data,
                                                           register_for_restart);
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeAllManagers, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_data_manager_instances[name];
} // getManager

void
FEDataManager::freeAllManagers()
{
    for (auto& data_manager_instance : s_data_manager_instances)
    {
        if (data_manager_instance.second)
        {
            delete data_manager_instance.second;
        }
        data_manager_instance.second = nullptr;
    }
    return;
} // freeAllManagers

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
FEDataManager::setEquationSystems(EquationSystems* const equation_systems, const int level_number)
{
    d_fe_data->setEquationSystems(equation_systems, level_number);
    return;
} // setEquationSystems

EquationSystems*
FEDataManager::getEquationSystems() const
{
    return d_fe_data->getEquationSystems();
} // getEquationSystems

FEDataManager::SystemDofMapCache*
FEDataManager::getDofMapCache(const std::string& system_name)
{
    return d_fe_data->getDofMapCache(system_name);
} // getDofMapCache

FEDataManager::SystemDofMapCache*
FEDataManager::getDofMapCache(unsigned int system_num)
{
    return d_fe_data->getDofMapCache(system_num);
} // getDofMapCache

void
FEDataManager::setPatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    // Reset the hierarchy.
    TBOX_ASSERT(hierarchy);
    d_hierarchy = hierarchy;
    if (!d_eulerian_data_cache)
    {
        d_eulerian_data_cache = std::make_shared<SAMRAIDataCache>();
        d_eulerian_data_cache->setPatchHierarchy(hierarchy);
        d_eulerian_data_cache->resetLevels(0, hierarchy->getFinestLevelNumber());
    }
    else
    {
        d_eulerian_data_cache->setPatchHierarchy(hierarchy);
    }
    return;
} // setPatchHierarchy

Pointer<PatchHierarchy<NDIM> >
FEDataManager::getPatchHierarchy() const
{
    return d_hierarchy;
} // getPatchHierarchy

int
FEDataManager::getCoarsestPatchLevelNumber() const
{
    for (int ln = 0; ln <= d_max_level_number; ++ln)
    {
        if (d_level_lookup.levelHasElements(ln)) return ln;
    }
    TBOX_ERROR("There must be a coarsest patch level number.");
    return IBTK::invalid_level_number;
}

int
FEDataManager::getFinestPatchLevelNumber() const
{
    for (int ln = d_max_level_number; 0 <= ln; --ln)
    {
        if (d_level_lookup.levelHasElements(ln)) return ln;
    }
    TBOX_ERROR("There must be a finest patch level number.");
    return IBTK::invalid_level_number;
}

const IntVector<NDIM>&
FEDataManager::getGhostCellWidth() const
{
    return d_ghost_width;
} // getGhostCellWidth

const FEDataManager::InterpSpec&
FEDataManager::getDefaultInterpSpec() const
{
    return d_default_interp_spec;
} // getDefaultInterpSpec

const FEDataManager::SpreadSpec&
FEDataManager::getDefaultSpreadSpec() const
{
    return d_default_spread_spec;
} // getDefaultSpreadSpec

const std::vector<std::vector<Elem*> >&
FEDataManager::getActivePatchElementMap() const
{
    return d_active_patch_elem_map.back();
} // getActivePatchElementMap

const std::vector<std::vector<Node*> >&
FEDataManager::getActivePatchNodeMap() const
{
    return d_active_patch_node_map.back();
} // getActivePatchNodeMap

void
FEDataManager::reinitElementMappings()
{
    IBTK_TIMER_START(t_reinit_element_mappings);

    // We reinitialize mappings after repartitioning, so clear the cache since
    // its content is no longer relevant:
    d_fe_data->clearPatchHierarchyDependentData();

    // The ghosted diagonal mass matrix needs to be rebuilt after repartitioning.
    d_L2_proj_matrix_diag_ghost.clear();

    // Delete cached hierarchy-dependent data.
    d_active_patch_elem_map.clear();
    d_active_patch_elem_map.resize(d_max_level_number + 1);
    d_active_patch_node_map.clear();
    d_active_patch_node_map.resize(d_max_level_number + 1);
    d_active_patch_ghost_dofs.clear();
    d_active_elems.clear();
    d_system_ghost_vec.clear();
    d_system_ib_ghost_vec.clear();

    // Reset the mappings between grid patches and active mesh
    // elements.
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        collectActivePatchElements(d_active_patch_elem_map[ln], ln, ln, ln);
        collectActivePatchNodes(d_active_patch_node_map[ln], d_active_patch_elem_map[ln]);
    }

    std::set<Elem*> elem_set;
    for (const std::vector<std::vector<Elem*> >& level_elems : d_active_patch_elem_map)
    {
        for (const std::vector<Elem*>& patch_elems : level_elems)
        {
            elem_set.insert(patch_elems.begin(), patch_elems.end());
        }
    }
    d_active_elems.assign(elem_set.begin(), elem_set.end());

    // If we are not regridding in the usual way (i.e., if
    // IBHierarchyIntegrator::d_regrid_cfl_interval > 1) then it is possible
    // that an element has traveled outside of it's assigned patch level. If
    // this happens then IBFE won't work - we cannot correctly interpolate
    // velocity at that point. Hence try to detect it by checking that all nodes
    // are on the interior of some patch (or outside the domain) at the moment.
    {
        const Pointer<CartesianGridGeometry<NDIM> > hier_geom = d_hierarchy->getGridGeometry();
        // TODO - we only support single box geometries right now
        TBOX_ASSERT(hier_geom);
        const double* const hier_x_lower = hier_geom->getXLower();
        const double* const hier_x_upper = hier_geom->getXUpper();
        const int rank = IBTK_MPI::getRank();
        const int n_procs = IBTK_MPI::getNodes();
        const MeshBase& mesh = getEquationSystems()->get_mesh();
        std::unique_ptr<PetscVector<double> > X_petsc_vec = buildIBGhostedVector(COORDINATES_SYSTEM_NAME);
        *X_petsc_vec = *getCoordsVector();
        X_petsc_vec->close();
        const double* const X_local_soln = X_petsc_vec->get_array_read();
        const DofMap& X_dof_map = d_fe_data->d_es->get_system(COORDINATES_SYSTEM_NAME).get_dof_map();

        std::vector<int> node_ranks(mesh.parallel_n_nodes());
        std::vector<dof_id_type> X_idxs;
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            int local_patch_num = 0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
            {
                const Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* const patch_x_lower = patch_geom->getXLower();
                const double* const patch_x_upper = patch_geom->getXUpper();

                for (const Node* n : d_active_patch_node_map[ln][local_patch_num])
                {
                    IBTK::Point X;
                    bool inside_patch = true;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        IBTK::get_nodal_dof_indices(X_dof_map, n, d, X_idxs);
                        X[d] = X_local_soln[X_petsc_vec->map_global_to_local_index(X_idxs[0])];
                        // Due to how SAMRAI computes patch boundaries, even if the
                        // patch's domain is [0, 1]^2 the patches on the boundary
                        // may not actually end at 1.0. Hence allow a small
                        // tolerance here to account for the case where two patches
                        // are adjacent to each-other but their patch boundaries
                        // don't quite line up:
                        const double x_lower = patch_x_lower[d] - std::max(1.0, std::abs(patch_x_lower[d])) * 1e-14;
                        const double x_upper = patch_x_upper[d] + std::max(1.0, std::abs(patch_x_upper[d])) * 1e-14;
                        inside_patch = inside_patch && (x_lower <= X[d] && X[d] <= x_upper);
                    }
                    if (inside_patch)
                    {
                        node_ranks[n->id()] = rank + 1;
                    }
                    else
                    {
                        // Points are allowed to be outside the domain - they simply
                        // are no longer used for IB calculations.
                        bool inside_hier = true;
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            // Like above - permit points very close to the boundary
                            // to pass the check and treat them as being outside
                            const double x_lower = hier_x_lower[d] + std::max(1.0, std::abs(hier_x_lower[d])) * 1e-14;
                            const double x_upper = hier_x_upper[d] - std::max(1.0, std::abs(hier_x_upper[d])) * 1e-14;
                            inside_hier = inside_hier && (x_lower <= X[d] && X[d] <= x_upper);
                        }
                        if (!inside_hier)
                        {
                            node_ranks[n->id()] = n_procs + 1;
                        }
                    }
                }
            }
        }
        X_petsc_vec->restore_array();

        // send everything to rank 0 instead of doing an all-to-all:
        IBTK_MPI::allToOneSumReduction(node_ranks.data(), node_ranks.size());
        if (rank == 0)
        {
            const std::string message =
                "At least one node in the current mesh is inside the fluid domain and not associated with any "
                "patch. This class currently assumes that all elements are on the finest level and will not "
                "work correctly if this assumption does not hold. This usually happens when you use multiple "
                "patch levels and set the regrid CFL interval to a value larger than one. To change this check "
                "set node_outside_patch_check to a different value in the input database: see the documentation "
                "of FEDataManager for more information.";
            for (const int node_rank : node_ranks)
            {
                if (node_rank == 0)
                {
                    switch (d_node_patch_check)
                    {
                    case NODE_OUTSIDE_PERMIT:
                        break;
                    case NODE_OUTSIDE_WARN:
                        TBOX_WARNING(message);
                        break;
                    case NODE_OUTSIDE_ERROR:
                        TBOX_ERROR(message);
                        break;
                    default:
                        // we shouldn't get here
                        TBOX_ERROR("unrecognized value for d_node_patch_check");
                        break;
                    }
                    // no need to check more nodes if we already found one outside
                    break;
                }
            }
        }
    }

    IBTK_TIMER_STOP(t_reinit_element_mappings);
    return;
} // reinitElementMappings

NumericVector<double>*
FEDataManager::getSolutionVector(const std::string& system_name) const
{
    return d_fe_data->d_es->get_system(system_name).solution.get();
} // getSolutionVector

NumericVector<double>*
FEDataManager::buildGhostedSolutionVector(const std::string& system_name, const bool localize_data)
{
    IBTK_TIMER_START(t_build_ghosted_solution_vector);

    reinitializeIBGhostedDOFs(system_name);
    NumericVector<double>* sol_vec = getSolutionVector(system_name);
    TBOX_ASSERT(sol_vec);
    if (!d_system_ghost_vec.count(system_name))
    {
        if (d_enable_logging)
        {
            plog << "FEDataManager::buildGhostedSolutionVector(): building ghosted solution vector for system: "
                 << system_name << "\n";
        }
        TBOX_ASSERT(d_active_patch_ghost_dofs.count(system_name));
        std::unique_ptr<NumericVector<double> > sol_ghost_vec = NumericVector<double>::build(sol_vec->comm());
        sol_ghost_vec->init(
            sol_vec->size(), sol_vec->local_size(), d_active_patch_ghost_dofs[system_name], true, GHOSTED);
        d_system_ghost_vec[system_name] = std::move(sol_ghost_vec);
    }
    NumericVector<double>* sol_ghost_vec = d_system_ghost_vec[system_name].get();
    if (localize_data) copy_and_synch(*sol_vec, *sol_ghost_vec, /*close_v_in*/ false);

    IBTK_TIMER_STOP(t_build_ghosted_solution_vector);
    return sol_ghost_vec;
} // buildGhostedSolutionVector

std::unique_ptr<PetscVector<double> >
FEDataManager::buildIBGhostedVector(const std::string& system_name)
{
    IBTK_TIMER_START(t_build_ghosted_vector);

    reinitializeIBGhostedDOFs(system_name);
    TBOX_ASSERT(d_system_ib_ghost_vec.find(system_name) != d_system_ib_ghost_vec.end());
    const std::unique_ptr<PetscVector<double> >& exemplar_ib_vector = d_system_ib_ghost_vec.at(system_name);
    TBOX_ASSERT(exemplar_ib_vector);
    // Sanity check:
    const System& system = d_fe_data->d_es->get_system(system_name);
    if (system.solution->local_size() != exemplar_ib_vector->local_size())
    {
        TBOX_ERROR(
            "The locally stored example IB vector does not have the same local size as the relevant system's solution "
            "vector. This usually occurs when the finite element data has been modified (e.g., the mesh has been "
            "refined) outside of one of IBAMR's classes (like FEDataManager or IBFEMethod).");
    }

    std::unique_ptr<NumericVector<double> > clone = exemplar_ib_vector->zero_clone();
    auto ptr = dynamic_cast<PetscVector<double>*>(clone.release());
    TBOX_ASSERT(ptr);

    IBTK_TIMER_STOP(t_build_ghosted_vector);
    return std::unique_ptr<PetscVector<double> >(ptr);
}

NumericVector<double>*
FEDataManager::getCoordsVector() const
{
    return getSolutionVector(COORDINATES_SYSTEM_NAME);
} // getCoordsVector

NumericVector<double>*
FEDataManager::buildGhostedCoordsVector(const bool localize_data)
{
    return buildGhostedSolutionVector(COORDINATES_SYSTEM_NAME, localize_data);
} // buildGhostedCoordsVector

std::shared_ptr<FEData>&
FEDataManager::getFEData()
{
    return d_fe_data;
} // getFEData

const std::shared_ptr<FEData>&
FEDataManager::getFEData() const
{
    return d_fe_data;
} // getFEData

/**
 * @brief Compute the product of the finite element representation of the
 * force and (optionally) a set of weights (e.g., JxW values) at all
 * quadrature points on an element. See integrate_elem_rhs for details on the
 * first two template arguments.
 *
 * @tparam weights_are_unity Assume that all values in @p weights are 1. This
 * value should be true when we are only interpolating values at quadrature
 * points and otherwise false.
 *
 * @param[in] qp_offset Offset used to calculate an index into @p F_w_qp. The
 * relevant part of the array is assumed to start at <code>n_vars *
 * qp_offset</code>.
 *
 * @param[in] phi_F Values of test functions evaluated at quadrature points,
 * indexed by test function number and then quadrature point number.
 *
 * @param[in] weights Vector containing weights at quadrature points: in the
 * case of force spreading this is the standard JxW array. If
 * <code>weights_are_unity</code> is true then its values are assumed to be 1.
 *
 * @param[in] F_node Values of the finite element solution (i.e., the
 * multipliers on each trial function) on the current element.
 *
 * @param[out] F_w_qp Globally indexed array containing the product of the
 * finite element representation of the force and weights at the quadrature
 * points on the current element. The array is indexed by quadrature point and
 * then by variable: i.e., if @p n_vars is greater than one then variables of
 * the vector-valued force being evaluated at a certain quadrature point are
 * contiguous.
 */
template <int n_vars, int n_basis, bool weights_are_unity = false>
void
sum_weighted_elem_solution_n_vars_n_basis(const int qp_offset,
                                          const std::vector<std::vector<double> >& phi_F,
                                          const std::vector<double>& weights,
                                          const boost::multi_array<double, 2>& F_node,
                                          std::vector<double>& F_w_qp)
{
    const int n_qp = phi_F[0].size();
    if (n_vars == -1 || n_basis == -1)
    {
        const int n_basis_ = phi_F.size();
        const int n_vars_ = F_node.shape()[1];
        for (int qp = 0; qp < n_qp; ++qp)
        {
            const int idx = n_vars_ * (qp_offset + qp);
            for (int k = 0; k < n_basis_; ++k)
                for (int i = 0; i < n_vars_; ++i) F_w_qp[idx + i] += F_node[k][i] * phi_F[k][qp];
            if (!weights_are_unity)
                for (int i = 0; i < n_vars_; ++i) F_w_qp[idx + i] *= weights[qp];
        }
    }
    else
    {
        for (int qp = 0; qp < n_qp; ++qp)
        {
            const int idx = n_vars * (qp_offset + qp);
            for (int k = 0; k < n_basis; ++k)
                for (int i = 0; i < n_vars; ++i) F_w_qp[idx + i] += F_node[k][i] * phi_F[k][qp];
            if (!weights_are_unity)
                for (int i = 0; i < n_vars; ++i) F_w_qp[idx + i] *= weights[qp];
        }
    }
}

// Dispatch to the above function based on the number of variables and number
// of basis functions.
template <int n_vars, bool weights_are_unity = false>
void
sum_weighted_elem_solution_n_vars(const int n_basis,
                                  const int qp_offset,
                                  const std::vector<std::vector<double> >& phi_F,
                                  const std::vector<double>& weights,
                                  const boost::multi_array<double, 2>& F_node,
                                  std::vector<double>& F_w_qp)
{
    switch (n_basis)
    {
#define SWITCH_CASE(k)                                                                                                 \
    case k:                                                                                                            \
    {                                                                                                                  \
        sum_weighted_elem_solution_n_vars_n_basis<n_vars, k, weights_are_unity>(                                       \
            qp_offset, phi_F, weights, F_node, F_w_qp);                                                                \
        return;                                                                                                        \
    }
        SWITCH_CASE(1)
        SWITCH_CASE(2)
        SWITCH_CASE(3)
        SWITCH_CASE(4)
        SWITCH_CASE(5)
        SWITCH_CASE(6)
        SWITCH_CASE(8)
        SWITCH_CASE(9)
        SWITCH_CASE(10)
        SWITCH_CASE(27)
#undef SWITCH_CASE
    default:
        sum_weighted_elem_solution_n_vars_n_basis<n_vars, -1, weights_are_unity>(
            qp_offset, phi_F, weights, F_node, F_w_qp);
        return;
    }
    return;
}

// Dispatch to the above function based on the number of variables.
template <bool weights_are_unity = false>
void
sum_weighted_elem_solution(const int n_vars,
                           const int n_basis,
                           const int qp_offset,
                           const std::vector<std::vector<double> >& phi_F,
                           const std::vector<double>& weights,
                           const boost::multi_array<double, 2>& F_node,
                           std::vector<double>& F_w_qp)
{
    switch (n_vars)
    {
#define SWITCH_CASE(k)                                                                                                 \
    case k:                                                                                                            \
    {                                                                                                                  \
        sum_weighted_elem_solution_n_vars<k, weights_are_unity>(n_basis, qp_offset, phi_F, weights, F_node, F_w_qp);   \
        return;                                                                                                        \
    }
        SWITCH_CASE(1)
        SWITCH_CASE(2)
        SWITCH_CASE(3)
        SWITCH_CASE(4)
        SWITCH_CASE(5)
#undef SWITCH_CASE
    default:
        sum_weighted_elem_solution_n_vars<-1, weights_are_unity>(n_basis, qp_offset, phi_F, weights, F_node, F_w_qp);
    }
    return;
}

void
FEDataManager::spread(const int f_data_idx,
                      NumericVector<double>& F_vec,
                      NumericVector<double>& X_vec,
                      const std::string& system_name)
{
    spread(f_data_idx, F_vec, X_vec, system_name, d_default_spread_spec);

    return;
} // spread

void
FEDataManager::spread(const int f_data_idx,
                      NumericVector<double>& F_vec,
                      NumericVector<double>& X_vec,
                      const std::string& system_name,
                      const FEDataManager::SpreadSpec& spread_spec)
{
    IBTK_TIMER_START(t_spread);

    // Determine the type of data centering.
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<hier::Variable<NDIM> > f_var;
    var_db->mapIndexToVariable(f_data_idx, f_var);
    Pointer<CellVariable<NDIM, double> > f_cc_var = f_var;
    Pointer<SideVariable<NDIM, double> > f_sc_var = f_var;
    const bool cc_data = f_cc_var;
    const bool sc_data = f_sc_var;
    TBOX_ASSERT(cc_data || sc_data);

    // Extract the mesh.
    const MeshBase& mesh = d_fe_data->d_es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // Extract the FE systems and DOF maps, and setup the FECache objects.
    System& F_system = d_fe_data->d_es->get_system(system_name);
    const unsigned int n_vars = F_system.n_vars();
    const DofMap& F_dof_map = F_system.get_dof_map();
    FEData::SystemDofMapCache& F_dof_map_cache = *getDofMapCache(system_name);
    System& X_system = d_fe_data->d_es->get_system(COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    FEData::SystemDofMapCache& X_dof_map_cache = *getDofMapCache(COORDINATES_SYSTEM_NAME);
    FEType F_fe_type = F_dof_map.variable_type(0);
    Order F_order = F_dof_map.variable_order(0);
    for (unsigned i = 0; i < n_vars; ++i)
    {
        TBOX_ASSERT(F_dof_map.variable_type(i) == F_fe_type);
        TBOX_ASSERT(F_dof_map.variable_order(i) == F_order);
    }
    FEType X_fe_type = X_dof_map.variable_type(0);
    Order X_order = X_dof_map.variable_order(0);
    for (unsigned d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
        TBOX_ASSERT(X_dof_map.variable_order(d) == X_order);
    }

    // convenience alias for the quadrature key type used by FECache and FEMappingCache
    using quad_key_type = quadrature_key_type;
    FECache F_fe_cache(dim, F_fe_type, FEUpdateFlags::update_phi);
    FECache X_fe_cache(dim, X_fe_type, FEUpdateFlags::update_phi);

    // We have to support both volumetric and surface meshes based on runtime data
    FEMappingCache<NDIM, NDIM> volume_mapping_cache(FEUpdateFlags::update_JxW);
    FEMappingCache<NDIM - 1, NDIM> surface_mapping_cache(FEUpdateFlags::update_JxW);
    const bool is_volume_mesh = dim == NDIM;

    // Check to see if we are using nodal quadrature.
    const bool use_nodal_quadrature = spread_spec.use_nodal_quadrature;
    if (use_nodal_quadrature) TBOX_ASSERT(F_fe_type == X_fe_type && F_order == X_order);

    // This will break, at some point in the future, if we ever use
    // nonnodal-interpolating finite elements. TODO: more finite elements
    // will probably work.
    std::vector<FEFamily> fe_family_whitelist{ LAGRANGE, L2_LAGRANGE };
    TBOX_ASSERT(std::find(fe_family_whitelist.begin(), fe_family_whitelist.end(), F_fe_type.family) !=
                fe_family_whitelist.end());
    TBOX_ASSERT(std::find(fe_family_whitelist.begin(), fe_family_whitelist.end(), X_fe_type.family) !=
                fe_family_whitelist.end());

    if (use_nodal_quadrature)
    {
        // Multiply by the nodal volume fractions (to convert densities into
        // values).
        PetscVector<double>* dX_vec = buildIBGhostedDiagonalL2MassMatrix(system_name);
        std::unique_ptr<NumericVector<double> > F_x_dX_vec = F_vec.clone();
        F_x_dX_vec->pointwise_mult(F_vec, *dX_vec);

        // Extract local form vectors.
        auto F_x_dX_petsc_vec = static_cast<PetscVector<double>*>(F_x_dX_vec.get());
        const double* const F_x_dX_local_soln = F_x_dX_petsc_vec->get_array_read();
        auto X_petsc_vec = static_cast<PetscVector<double>*>(&X_vec);
        const double* const X_local_soln = X_petsc_vec->get_array_read();

        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            int local_patch_num = 0;
            // Spread from the nodes.
            for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
            {
                // The relevant collection of nodes.
                const std::vector<Node*>& patch_nodes = d_active_patch_node_map[ln][local_patch_num];
                const size_t num_active_patch_nodes = patch_nodes.size();
                if (!num_active_patch_nodes) continue;

                const Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* const patch_x_lower = patch_geom->getXLower();
                const double* const patch_x_upper = patch_geom->getXUpper();
                std::array<bool, NDIM> touches_upper_regular_bdry;
                for (unsigned int d = 0; d < NDIM; ++d)
                    touches_upper_regular_bdry[d] = patch_geom->getTouchesRegularBoundary(d, 1);

                // Store the values of F_JxW and X at the nodes inside the patch.
                std::vector<double> F_x_dX_node, X_node;
                F_x_dX_node.reserve(n_vars * num_active_patch_nodes);
                X_node.reserve(NDIM * num_active_patch_nodes);
                std::vector<dof_id_type> F_idxs, X_idxs;
                IBTK::Point X;
                for (unsigned int k = 0; k < num_active_patch_nodes; ++k)
                {
                    const Node* const n = patch_nodes[k];
                    bool inside_patch = true;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        IBTK::get_nodal_dof_indices(X_dof_map, n, d, X_idxs);
                        X[d] = X_local_soln[X_petsc_vec->map_global_to_local_index(X_idxs[0])];
                        inside_patch =
                            inside_patch && (X[d] >= patch_x_lower[d]) &&
                            ((X[d] < patch_x_upper[d]) || (touches_upper_regular_bdry[d] && X[d] <= patch_x_upper[d]));
                    }
                    if (inside_patch)
                    {
                        for (unsigned int i = 0; i < n_vars; ++i)
                        {
                            IBTK::get_nodal_dof_indices(F_dof_map, n, i, F_idxs);
                            for (const auto& F_idx : F_idxs)
                            {
                                F_x_dX_node.push_back(
                                    F_x_dX_local_soln[F_x_dX_petsc_vec->map_global_to_local_index(F_idx)]);
                            }
                        }
                        X_node.insert(X_node.end(), &X[0], &X[0] + NDIM);
                    }
                }
                TBOX_ASSERT(F_x_dX_node.size() <= n_vars * num_active_patch_nodes);
                TBOX_ASSERT(X_node.size() <= NDIM * num_active_patch_nodes);

                // Spread values from the nodes to the Cartesian grid patch.
                //
                // \todo Fix this for FE structures with periodic boundaries.  Nodes
                // on periodic boundaries will be "double counted".
                //
                // \todo Add warnings for FE structures with periodic boundaries.
                const Box<NDIM> spread_box = patch->getBox();
                Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_data_idx);
                if (cc_data)
                {
                    Pointer<CellData<NDIM, double> > f_cc_data = f_data;
                    LEInteractor::spread(
                        f_cc_data, F_x_dX_node, n_vars, X_node, NDIM, patch, spread_box, spread_spec.kernel_fcn);
                }
                if (sc_data)
                {
                    Pointer<SideData<NDIM, double> > f_sc_data = f_data;
                    LEInteractor::spread(
                        f_sc_data, F_x_dX_node, n_vars, X_node, NDIM, patch, spread_box, spread_spec.kernel_fcn);
                }
            }
        }

        // Restore local form vectors.
        F_x_dX_petsc_vec->restore_array();
        X_petsc_vec->restore_array();
    }
    else
    {
        // Extract local form vectors.
        auto F_petsc_vec = static_cast<PetscVector<double>*>(&F_vec);
        const double* const F_local_soln = F_petsc_vec->get_array_read();

        auto X_petsc_vec = static_cast<PetscVector<double>*>(&X_vec);
        const double* const X_local_soln = X_petsc_vec->get_array_read();

        // Loop over the patches to interpolate nodal values on the FE mesh to
        // the element quadrature points, then spread those values onto the
        // Eulerian grid.
        boost::multi_array<double, 2> F_node;
        std::vector<double> F_JxW_qp, X_qp;
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            int local_patch_num = 0;
            // Spread from the nodes.
            for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
            {
                // The relevant collection of elements.
                const std::vector<Elem*>& patch_elems = d_active_patch_elem_map[ln][local_patch_num];
                const size_t num_active_patch_elems = patch_elems.size();
                if (!num_active_patch_elems) continue;

                const Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* const patch_dx = patch_geom->getDx();
                const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);

                // Determining which quadrature rule should be used on which
                // processor is surprisingly expensive, so cache the keys:
                std::vector<quad_key_type> quad_keys(num_active_patch_elems);

                // Cache interpolated positions too:
                std::vector<boost::multi_array<double, 2> > X_nodes(num_active_patch_elems);

                // Setup vectors to store the values of F_JxW and X at the
                // quadrature points.
                unsigned int n_qp_patch = 0;
                for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
                {
                    Elem* const elem = patch_elems[e_idx];
#ifndef NDEBUG
                    TBOX_ASSERT(getPatchLevel(elem) == ln);
#endif // ifndef NDEBUG
                    const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                    get_values_for_interpolation(X_nodes[e_idx], *X_petsc_vec, X_local_soln, X_dof_indices);
                    const quad_key_type key = getQuadratureKey(spread_spec.quad_type,
                                                               spread_spec.quad_order,
                                                               spread_spec.use_adaptive_quadrature,
                                                               spread_spec.point_density,
                                                               spread_spec.allow_rules_with_negative_weights,
                                                               elem,
                                                               X_nodes[e_idx],
                                                               patch_dx_min);
                    quad_keys[e_idx] = key;
                    QBase& qrule = d_fe_data->d_quadrature_cache[key];
                    n_qp_patch += qrule.n_points();
                }
                if (!n_qp_patch) continue;
                F_JxW_qp.resize(n_vars * n_qp_patch);
                X_qp.resize(NDIM * n_qp_patch);

                // Loop over the elements and compute the values to be spread and
                // the positions of the quadrature points.
                int qp_offset = 0;
                for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
                {
                    Elem* const elem = patch_elems[e_idx];
                    const auto& F_dof_indices = F_dof_map_cache.dof_indices(elem);
                    get_values_for_interpolation(F_node, *F_petsc_vec, F_local_soln, F_dof_indices);
                    const quad_key_type& key = quad_keys[e_idx];
                    const FEBase& X_fe = X_fe_cache(key, elem);
                    const FEBase& F_fe = F_fe_cache(key, elem);
                    const QBase& qrule = d_fe_data->d_quadrature_cache[key];

                    // JxW depends on the element
                    const std::vector<double>& JxW_F =
                        get_JxW(key, elem, is_volume_mesh, volume_mapping_cache, surface_mapping_cache);
                    const std::vector<std::vector<double> >& phi_F = F_fe.get_phi();
                    const std::vector<std::vector<double> >& phi_X = X_fe.get_phi();

                    const unsigned int n_qp = qrule.n_points();
                    TBOX_ASSERT(n_qp == phi_F[0].size());
                    TBOX_ASSERT(n_qp == phi_X[0].size());
                    TBOX_ASSERT(n_qp == JxW_F.size());
                    double* F_begin = &F_JxW_qp[n_vars * qp_offset];
                    double* X_begin = &X_qp[NDIM * qp_offset];
                    std::fill(F_begin, F_begin + n_vars * n_qp, 0.0);
                    std::fill(X_begin, X_begin + NDIM * n_qp, 0.0);

                    sum_weighted_elem_solution</*weights_are_unity*/ false>(
                        n_vars, F_dof_indices[0].size(), qp_offset, phi_F, JxW_F, F_node, F_JxW_qp);
                    sum_weighted_elem_solution</*weights_are_unity*/ true>(
                        NDIM, phi_X.size(), qp_offset, phi_X, {}, X_nodes[e_idx], X_qp);
                    qp_offset += n_qp;
                }

                zeroExteriorValues(*patch_geom, X_qp, F_JxW_qp, n_vars);

                // Spread values from the quadrature points to the Cartesian grid
                // patch.
                const Box<NDIM> spread_box = patch->getBox();
                Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_data_idx);
                if (cc_data)
                {
                    Pointer<CellData<NDIM, double> > f_cc_data = f_data;
                    LEInteractor::spread(
                        f_cc_data, F_JxW_qp, n_vars, X_qp, NDIM, patch, spread_box, spread_spec.kernel_fcn);
                }
                if (sc_data)
                {
                    Pointer<SideData<NDIM, double> > f_sc_data = f_data;
                    LEInteractor::spread(
                        f_sc_data, F_JxW_qp, n_vars, X_qp, NDIM, patch, spread_box, spread_spec.kernel_fcn);
                }
            }
        }

        // Restore local form vectors.
        F_petsc_vec->restore_array();
        X_petsc_vec->restore_array();
    }

    IBTK_TIMER_STOP(t_spread);
    return;
} // spread

void
FEDataManager::spread(const int f_data_idx,
                      NumericVector<double>& F_vec,
                      NumericVector<double>& X_vec,
                      const std::string& system_name,
                      RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                      const double fill_data_time,
                      const bool close_F,
                      const bool close_X)
{
    spread(
        f_data_idx, F_vec, X_vec, system_name, d_default_spread_spec, f_phys_bdry_op, fill_data_time, close_X, close_F);
    return;
} // spread

void
FEDataManager::spread(const int f_data_idx,
                      NumericVector<double>& F_vec,
                      NumericVector<double>& X_vec,
                      const std::string& system_name,
                      const FEDataManager::SpreadSpec& spread_spec,
                      RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                      const double fill_data_time,
                      const bool close_F,
                      const bool close_X)
{
    // Determine the type of data centering.
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<hier::Variable<NDIM> > f_var;
    var_db->mapIndexToVariable(f_data_idx, f_var);
    Pointer<CellVariable<NDIM, double> > f_cc_var = f_var;
    Pointer<SideVariable<NDIM, double> > f_sc_var = f_var;
    const bool cc_data = f_cc_var;
    const bool sc_data = f_sc_var;
    TBOX_ASSERT(cc_data || sc_data);

    // Make a copy of the Eulerian data.
    const auto f_copy_data_idx = d_eulerian_data_cache->getCachedPatchDataIndex(f_data_idx);
    Pointer<HierarchyDataOpsReal<NDIM, double> > f_data_ops =
        HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(f_var, d_hierarchy, true);
    f_data_ops->resetLevels(0, d_hierarchy->getFinestLevelNumber());
    f_data_ops->swapData(f_copy_data_idx, f_data_idx);
    f_data_ops->setToScalar(f_data_idx, 0.0, /*interior_only*/ false);

    // Communicate any unsynchronized ghost data.
    if (close_F) F_vec.close();
    if (close_X) X_vec.close();

    // Spread values.
    spread(f_data_idx, F_vec, X_vec, system_name, spread_spec);

    // sum values outside the physical domain into the physical domain.
    if (f_phys_bdry_op)
    {
        f_phys_bdry_op->setPatchDataIndex(f_data_idx);
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            int local_patch_num = 0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
            {
                // If there are no elements on this patch then we can skip it
                if (d_active_patch_elem_map[ln][local_patch_num].size() == 0) continue;
                const Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_data_idx);
                f_phys_bdry_op->accumulateFromPhysicalBoundaryData(*patch, fill_data_time, f_data->getGhostCellWidth());
            }
        }
    }

    // Accumulate data.
    f_data_ops->swapData(f_copy_data_idx, f_data_idx);
    f_data_ops->add(f_data_idx, f_data_idx, f_copy_data_idx, /*interior_only*/ false);

    return;
} // spread

void
FEDataManager::prolongData(const int f_data_idx,
                           NumericVector<double>& F_vec,
                           NumericVector<double>& X_vec,
                           const std::string& system_name,
                           const bool is_density,
                           const bool accumulate_on_grid,
                           const bool close_F,
                           const bool close_X)
{
    IBTK_TIMER_START(t_prolong_data);

    // NOTE #1: This routine is specialized for a staggered-grid Eulerian
    // discretization.  It should be straightforward to generalize it to work
    // with other data centerings.
    //
    // NOTE #2: This code is specialized for isoparametric elements.  It is less
    // clear how to relax this assumption.
    //
    // NOTE #3: This implementation optionally uses the pointwise value of J =
    // det(dX/ds) to convert a Lagrangian density into an Eulerian density.  We
    // should investigate whether there is any advantage to using a projection
    // of J onto a (discontinuous) FE basis instead of evaluating J directly
    // from the discrete deformation.

    // Extract the mesh.
    const MeshBase& mesh = d_fe_data->d_es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    TBOX_ASSERT(dim == NDIM);

    // Extract the FE systems and DOF maps, and setup the FE object.
    System& F_system = d_fe_data->d_es->get_system(system_name);
    const unsigned int n_vars = F_system.n_vars();
    TBOX_ASSERT(n_vars == NDIM); // specialized to side-centered data
    const DofMap& F_dof_map = F_system.get_dof_map();
    FEData::SystemDofMapCache& F_dof_map_cache = *getDofMapCache(system_name);
    System& X_system = d_fe_data->d_es->get_system(system_name);
    const DofMap& X_dof_map = X_system.get_dof_map();
    FEData::SystemDofMapCache& X_dof_map_cache = *getDofMapCache(COORDINATES_SYSTEM_NAME);
    FEType F_fe_type = F_dof_map.variable_type(0);
    for (unsigned i = 0; i < n_vars; ++i)
    {
        TBOX_ASSERT(F_dof_map.variable_type(i) == F_fe_type);
    }
    FEType X_fe_type = X_dof_map.variable_type(0);
    for (unsigned d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
    }
    std::unique_ptr<FEBase> F_fe_autoptr(FEBase::build(dim, F_fe_type)), X_fe_autoptr;
    if (F_fe_type != X_fe_type)
    {
        X_fe_autoptr = std::unique_ptr<FEBase>(FEBase::build(dim, X_fe_type));
    }
    FEBase* F_fe = F_fe_autoptr.get();
    FEBase* X_fe = X_fe_autoptr.get() ? X_fe_autoptr.get() : F_fe_autoptr.get();
    const std::vector<std::vector<double> >& phi_F = F_fe->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& dphi_X = X_fe->get_dphi();

    // Communicate any unsynchronized ghost data and extract the underlying
    // solution data.
    if (close_F) F_vec.close();
    auto F_petsc_vec = static_cast<PetscVector<double>*>(&F_vec);
    const double* const F_local_soln = F_petsc_vec->get_array_read();

    if (close_X) X_vec.close();
    auto X_petsc_vec = static_cast<PetscVector<double>*>(&X_vec);
    const double* const X_local_soln = X_petsc_vec->get_array_read();

    // Loop over the patches to interpolate nodal values on the FE mesh to the
    // points of the Eulerian grid.
    TensorValue<double> dX_ds;
    boost::multi_array<double, 2> F_node, X_node;
    std::vector<libMesh::Point> s_node_cache, X_node_cache;
    Point X_min, X_max;
    std::vector<libMesh::Point> intersection_ref_coords;
    std::vector<SideIndex<NDIM> > intersection_indices;
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const IntVector<NDIM>& ratio = level->getRatio();
        const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            // The relevant collection of elements.
            const std::vector<Elem*>& patch_elems = d_active_patch_elem_map[ln][local_patch_num];
            const size_t num_active_patch_elems = patch_elems.size();
            if (!num_active_patch_elems) continue;

            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_data_idx);
            if (!accumulate_on_grid) f_data->fillAll(0.0);
            const Box<NDIM>& patch_box = patch->getBox();
            const CellIndex<NDIM>& patch_lower = patch_box.lower();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_dx = patch_geom->getDx();

            std::array<Box<NDIM>, NDIM> side_boxes;
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                side_boxes[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
            }

            SideData<NDIM, int> num_intersections(patch_box, 1, IntVector<NDIM>(0));
            num_intersections.fillAll(0);

            // Loop over the elements and compute the values to be prolonged.
            for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems[e_idx];
                const unsigned int n_node = elem->n_nodes();
                const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);

                // Cache the nodal and physical coordinates of the element,
                // determine the bounding box of the current configuration of the
                // element, and set the nodal coordinates of the element to
                // correspond to the physical coordinates.
                s_node_cache.resize(n_node);
                X_node_cache.resize(n_node);
                X_min = Point::Constant(std::numeric_limits<double>::max());
                X_max = Point::Constant(-std::numeric_limits<double>::max());
                for (unsigned int k = 0; k < n_node; ++k)
                {
                    s_node_cache[k] = elem->point(k);
                    libMesh::Point& X = X_node_cache[k];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        X(d) = X_vec(X_dof_indices[d][k]);
                        X_min[d] = std::min(X_min[d], X(d));
                        X_max[d] = std::max(X_max[d], X(d));
                    }
                    elem->point(k) = X;
                }
                Box<NDIM> box(IndexUtilities::getCellIndex(&X_min[0], grid_geom, ratio),
                              IndexUtilities::getCellIndex(&X_max[0], grid_geom, ratio));
                box.grow(IntVector<NDIM>(1));
                box = box * patch_box;

                // Loop over coordinate directions and look for Eulerian grid points
                // that are covered by the element.
                intersection_ref_coords.clear();
                intersection_indices.clear();
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    // Loop over the relevant range of indices.
                    for (SideIterator<NDIM> b(box, axis); b; b++)
                    {
                        const SideIndex<NDIM>& i_s = b();
                        if (num_intersections(i_s) == 0 && side_boxes[axis].contains(i_s))
                        {
                            libMesh::Point p;
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                p(d) = patch_x_lower[d] + patch_dx[d] * (static_cast<double>(i_s(d) - patch_lower[d]) +
                                                                         (d == axis ? 0.0 : 0.5));
                            }
                            static const double TOL = std::sqrt(std::numeric_limits<double>::epsilon());
                            const libMesh::Point ref_coords =
                                FEInterface::inverse_map(dim, X_fe_type, elem, p, TOL, false);
                            if (FEInterface::on_reference_element(ref_coords, elem->type(), TOL))
                            {
                                intersection_ref_coords.push_back(ref_coords);
                                intersection_indices.push_back(i_s);
                                num_intersections(i_s) += 1;
                            }
                        }
                    }
                }

                // Restore the nodal coordinates.
                for (unsigned int k = 0; k < n_node; ++k)
                {
                    elem->point(k) = s_node_cache[k];
                }

                // If there are no intersection points, then continue on to the next
                // element.
                if (intersection_ref_coords.empty()) continue;

                // Evaluate the Lagrangian quantity at the Eulerian grid point and
                // update the data on the grid.
                const auto& F_dof_indices = F_dof_map_cache.dof_indices(elem);
                get_values_for_interpolation(F_node, *F_petsc_vec, F_local_soln, F_dof_indices);
                if (is_density) get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);
                F_fe->reinit(elem, &intersection_ref_coords);
                if (X_fe != F_fe) X_fe->reinit(elem, &intersection_ref_coords);
                for (unsigned int qp = 0; qp < intersection_ref_coords.size(); ++qp)
                {
                    const SideIndex<NDIM>& i_s = intersection_indices[qp];
                    const int axis = i_s.getAxis();
                    using range = boost::multi_array_types::index_range;
                    double F_qp = interpolate(qp, F_node[boost::indices[range(0, n_node)][axis]], phi_F);
                    if (is_density)
                    {
                        jacobian(dX_ds, qp, X_node, dphi_X);
                        F_qp /= std::abs(dX_ds.det());
                    }
                    (*f_data)(i_s) += F_qp / static_cast<double>(num_intersections(i_s));
                }
            }
        }
    }

    F_petsc_vec->restore_array();
    X_petsc_vec->restore_array();

    IBTK_TIMER_STOP(t_prolong_data);
    return;
} // prolongData

/**
 * @brief Assemble the element contribution to a load vector.
 *
 * The three functions below allow compile-time selection of the number of
 * basis functions and number of variables. Fixing the loop bounds via
 * templates and switch statements makes most IBFE codes a few percent faster.
 *
 * This function (the last of the three) takes the number of variables and
 * number of basis functions as template argument which allows the compiler to
 * unroll the two inner loops when appropriate. Setting either value to -1
 * results in run-time selection of loop bounds (which is considerably
 * slower).
 *
 * @param[in] qp_offset Offset used to calculate an index into @p F_qp. The
 * relevant part of the array is assumed to start at <code>n_vars *
 * qp_offset</code>.
 *
 * @param[in] phi_F Values of test functions evaluated at quadrature points,
 * indexed by test function number and then quadrature point number.
 *
 * @param[in] JxW_F Products of Jacobian and quadrature weight at each
 * quadrature point.
 *
 * @param[in] F_qp Globally indexed array containing function values at
 * quadrature points. The array is indexed by quadrature point and then by
 * variable: i.e., if @p n_vars is greater than one then components of the
 * vector-valued function being projected at a certain quadrature point are
 * contiguous.
 *
 * @param[out] F_rhs_concatenated Vector containing element integrals of
 * products of test functions and values of the interpolated function. Unlike
 * @p F_qp, the values in this vector are first indexed by variable and then
 * by quadrature point.
 */
template <int n_vars, int n_basis>
void
integrate_elem_rhs_n_vars_n_basis(const int qp_offset,
                                  const std::vector<std::vector<double> >& phi_F,
                                  const std::vector<double>& JxW_F,
                                  const std::vector<double>& F_qp,
                                  std::vector<double>& F_rhs_concatenated)
{
    const int n_qp = phi_F[0].size();
    if (n_vars == -1 || n_basis == -1)
    {
        const int n_basis_ = phi_F.size();
        const int n_vars_ = F_rhs_concatenated.size() / n_basis_;
        for (int qp = 0; qp < n_qp; ++qp)
        {
            const int idx = n_vars_ * (qp_offset + qp);
            for (int k = 0; k < n_basis_; ++k)
            {
                const double p_JxW_F = phi_F[k][qp] * JxW_F[qp];
                for (int i = 0; i < n_vars_; ++i) F_rhs_concatenated[n_basis_ * i + k] += F_qp[idx + i] * p_JxW_F;
            }
        }
    }
    else
    {
        for (int qp = 0; qp < n_qp; ++qp)
        {
            const int idx = n_vars * (qp_offset + qp);
            for (int k = 0; k < n_basis; ++k)
            {
                const double p_JxW_F = phi_F[k][qp] * JxW_F[qp];
                for (int i = 0; i < n_vars; ++i) F_rhs_concatenated[n_basis * i + k] += F_qp[idx + i] * p_JxW_F;
            }
        }
    }
}

// Dispatch to the above function based on the number of variables and number
// of basis functions.
template <int n_vars>
void
integrate_elem_rhs_n_vars(const int n_basis,
                          const int qp_offset,
                          const std::vector<std::vector<double> >& phi_F,
                          const std::vector<double>& JxW_F,
                          const std::vector<double>& F_qp,
                          std::vector<double>& F_rhs_concatenated)
{
    switch (n_basis)
    {
#define SWITCH_CASE(k)                                                                                                 \
    case k:                                                                                                            \
    {                                                                                                                  \
        integrate_elem_rhs_n_vars_n_basis<n_vars, k>(qp_offset, phi_F, JxW_F, F_qp, F_rhs_concatenated);               \
        return;                                                                                                        \
    }
        SWITCH_CASE(1)
        SWITCH_CASE(2)
        SWITCH_CASE(3)
        SWITCH_CASE(4)
        SWITCH_CASE(5)
        SWITCH_CASE(6)
        SWITCH_CASE(8)
        SWITCH_CASE(9)
        SWITCH_CASE(10)
        SWITCH_CASE(27)
#undef SWITCH_CASE
    default:
        integrate_elem_rhs_n_vars_n_basis<n_vars, -1>(qp_offset, phi_F, JxW_F, F_qp, F_rhs_concatenated);
        return;
    }
    return;
}

// Dispatch to the above function based on the number of variables.
void
integrate_elem_rhs(const int n_vars,
                   const int n_basis,
                   const int qp_offset,
                   const std::vector<std::vector<double> >& phi_F,
                   const std::vector<double>& JxW_F,
                   const std::vector<double>& F_qp,
                   std::vector<double>& F_rhs_concatenated)
{
    switch (n_vars)
    {
#define SWITCH_CASE(k)                                                                                                 \
    case k:                                                                                                            \
    {                                                                                                                  \
        integrate_elem_rhs_n_vars<k>(n_basis, qp_offset, phi_F, JxW_F, F_qp, F_rhs_concatenated);                      \
        return;                                                                                                        \
    }
        SWITCH_CASE(1)
        SWITCH_CASE(2)
        SWITCH_CASE(3)
        SWITCH_CASE(4)
        SWITCH_CASE(5)
#undef SWITCH_CASE
    default:
        integrate_elem_rhs_n_vars<-1>(n_basis, qp_offset, phi_F, JxW_F, F_qp, F_rhs_concatenated);
    }
    return;
}

void
FEDataManager::interpWeighted(const int f_data_idx,
                              NumericVector<double>& F_vec,
                              NumericVector<double>& X_vec,
                              const std::string& system_name,
                              const std::vector<Pointer<RefineSchedule<NDIM> > >& f_refine_scheds,
                              const double fill_data_time,
                              const bool close_F,
                              const bool close_X)
{
    interpWeighted(f_data_idx,
                   F_vec,
                   X_vec,
                   system_name,
                   d_default_interp_spec,
                   f_refine_scheds,
                   fill_data_time,
                   close_F,
                   close_X);
    return;
} // interpWeighted

void
FEDataManager::interpWeighted(const int f_data_idx,
                              NumericVector<double>& F_vec,
                              NumericVector<double>& X_vec,
                              const std::string& system_name,
                              const FEDataManager::InterpSpec& interp_spec,
                              const std::vector<Pointer<RefineSchedule<NDIM> > >& f_refine_scheds,
                              const double fill_data_time,
                              const bool close_F,
                              const bool close_X)
{
    IBTK_TIMER_START(t_interp_weighted);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Determine the type of data centering.
    Pointer<hier::Variable<NDIM> > f_var;
    var_db->mapIndexToVariable(f_data_idx, f_var);
    Pointer<CellVariable<NDIM, double> > f_cc_var = f_var;
    Pointer<SideVariable<NDIM, double> > f_sc_var = f_var;
    const bool cc_data = f_cc_var;
    const bool sc_data = f_sc_var;
    TBOX_ASSERT(cc_data || sc_data);

    // Extract the mesh.
    const MeshBase& mesh = d_fe_data->d_es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // Extract the FE systems and DOF maps, and setup the FECache objects.
    System& F_system = d_fe_data->d_es->get_system(system_name);
    const unsigned int n_vars = F_system.n_vars();
    const DofMap& F_dof_map = F_system.get_dof_map();
    FEData::SystemDofMapCache& F_dof_map_cache = *getDofMapCache(system_name);
    System& X_system = d_fe_data->d_es->get_system(COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    FEData::SystemDofMapCache& X_dof_map_cache = *getDofMapCache(COORDINATES_SYSTEM_NAME);
    FEType F_fe_type = F_dof_map.variable_type(0);
    Order F_order = F_dof_map.variable_order(0);
    for (unsigned i = 0; i < n_vars; ++i)
    {
        TBOX_ASSERT(F_dof_map.variable_type(i) == F_fe_type);
        TBOX_ASSERT(F_dof_map.variable_order(i) == F_order);
    }
    FEType X_fe_type = X_dof_map.variable_type(0);
    Order X_order = X_dof_map.variable_order(0);
    for (unsigned d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
        TBOX_ASSERT(X_dof_map.variable_order(d) == X_order);
    }

    // convenience alias for the quadrature key type used by FECache and FEMappingCache
    using quad_key_type = quadrature_key_type;
    FECache F_fe_cache(dim, F_fe_type, FEUpdateFlags::update_phi);
    FECache X_fe_cache(dim, X_fe_type, FEUpdateFlags::update_phi);
    FEMappingCache<NDIM, NDIM> volume_mapping_cache(FEUpdateFlags::update_JxW);
    FEMappingCache<NDIM - 1, NDIM> surface_mapping_cache(FEUpdateFlags::update_JxW);
    const bool is_volume_mesh = dim == NDIM;

    // Communicate any unsynchronized ghost data.
    for (const auto& f_refine_sched : f_refine_scheds)
    {
        if (f_refine_sched) f_refine_sched->fillData(fill_data_time);
    }

    if (close_X) X_vec.close();

    // Check to see if we are using nodal quadrature.
    const bool use_nodal_quadrature = interp_spec.use_nodal_quadrature;
    if (use_nodal_quadrature) TBOX_ASSERT(F_fe_type == X_fe_type && F_order == X_order);

    // This will break, at some point in the future, if we ever use
    // nonnodal-interpolating finite elements. TODO: more finite elements
    // will probably work.
    std::vector<FEFamily> fe_family_whitelist{ LAGRANGE, L2_LAGRANGE };
    TBOX_ASSERT(std::find(fe_family_whitelist.begin(), fe_family_whitelist.end(), F_fe_type.family) !=
                fe_family_whitelist.end());
    TBOX_ASSERT(std::find(fe_family_whitelist.begin(), fe_family_whitelist.end(), X_fe_type.family) !=
                fe_family_whitelist.end());

    // Extract local form vectors.
    auto X_petsc_vec = dynamic_cast<PetscVector<double>*>(&X_vec);
    TBOX_ASSERT(X_petsc_vec != nullptr);
    const double* const X_local_soln = X_petsc_vec->get_array_read();
    // Since we do a lot of assembly in this routine into off-processor
    // entries we will directly insert into the ghost values (and then
    // scatter in the calling function with the usual batch function).
    F_vec.zero();
    auto F_petsc_vec = dynamic_cast<PetscVector<double>*>(&F_vec);
    Vec F_local_form = nullptr;
    double* F_local_soln = nullptr;
    PetscInt F_local_size = -1;
    const bool is_ghosted = F_vec.type() == GHOSTED;
    if (is_ghosted)
    {
        TBOX_ASSERT(F_petsc_vec != nullptr);
        int ierr = VecGhostGetLocalForm(F_petsc_vec->vec(), &F_local_form);
        IBTK_CHKERRQ(ierr);
        ierr = VecGetArray(F_local_form, &F_local_soln);
        IBTK_CHKERRQ(ierr);
        ierr = VecGetSize(F_local_form, &F_local_size);
        IBTK_CHKERRQ(ierr);
    }

    if (use_nodal_quadrature)
    {
        // Extract local form vectors.
        PetscVector<double>* dX_vec = buildIBGhostedDiagonalL2MassMatrix(system_name);
        const double* const dX_local_soln = dX_vec->get_array_read();

        // Loop over the patches to interpolate values to the nodes from the grid, then use these values to
        // compute the projection of the interpolated velocity field onto the FE basis functions.
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            int local_patch_num = 0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
            {
                // The relevant collection of nodes.
                const std::vector<Node*>& patch_nodes = d_active_patch_node_map[ln][local_patch_num];
                const size_t num_active_patch_nodes = patch_nodes.size();
                if (!num_active_patch_nodes) continue;

                const Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* const patch_x_lower = patch_geom->getXLower();
                const double* const patch_x_upper = patch_geom->getXUpper();
                std::array<bool, NDIM> touches_upper_regular_bdry;
                for (unsigned int d = 0; d < NDIM; ++d)
                    touches_upper_regular_bdry[d] = patch_geom->getTouchesRegularBoundary(d, 1);

                // Store the value of X at the nodes that are inside the current
                // patch.
                std::vector<dof_id_type> F_node_idxs;
                std::vector<double> F_node, X_node;
                F_node_idxs.reserve(n_vars * num_active_patch_nodes);
                F_node.reserve(n_vars * num_active_patch_nodes);
                X_node.reserve(NDIM * num_active_patch_nodes);
                std::vector<dof_id_type> F_idxs, X_idxs;
                IBTK::Point X;
                for (unsigned int k = 0; k < num_active_patch_nodes; ++k)
                {
                    const Node* const n = patch_nodes[k];
                    bool inside_patch = true;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        IBTK::get_nodal_dof_indices(X_dof_map, n, d, X_idxs);
                        X[d] = X_local_soln[X_petsc_vec->map_global_to_local_index(X_idxs[0])];
                        inside_patch =
                            inside_patch && (X[d] >= patch_x_lower[d]) &&
                            ((X[d] < patch_x_upper[d]) || (touches_upper_regular_bdry[d] && X[d] <= patch_x_upper[d]));
                    }
                    if (inside_patch)
                    {
                        F_node.resize(F_node.size() + n_vars, 0.0);
                        X_node.insert(X_node.end(), &X[0], &X[0] + NDIM);
                        for (unsigned int i = 0; i < n_vars; ++i)
                        {
                            IBTK::get_nodal_dof_indices(F_dof_map, n, i, F_idxs);
                            F_node_idxs.insert(F_node_idxs.end(), F_idxs.begin(), F_idxs.end());
                        }
                    }
                }
                TBOX_ASSERT(F_node.size() <= n_vars * num_active_patch_nodes);
                TBOX_ASSERT(X_node.size() <= NDIM * num_active_patch_nodes);
                TBOX_ASSERT(F_node_idxs.size() <= n_vars * num_active_patch_nodes);

                if (F_node.empty()) continue;

                // Interpolate values from the Cartesian grid patch to the nodes.
                //
                // NOTE: The only nodes that are treated here are the nodes that are
                // inside of the patch.  The interpolation box is grown by a ghost
                // cell width of 1 to ensure that roundoff errors do not
                // inadvertently exclude the selected points.
                const Box<NDIM>& interp_box = Box<NDIM>::grow(patch->getBox(), IntVector<NDIM>(1));
                Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_data_idx);
                if (cc_data)
                {
                    Pointer<CellData<NDIM, double> > f_cc_data = f_data;
                    LEInteractor::interpolate(
                        F_node, n_vars, X_node, NDIM, f_cc_data, patch, interp_box, interp_spec.kernel_fcn);
                }
                if (sc_data)
                {
                    Pointer<SideData<NDIM, double> > f_sc_data = f_data;
                    LEInteractor::interpolate(
                        F_node, n_vars, X_node, NDIM, f_sc_data, patch, interp_box, interp_spec.kernel_fcn);
                }

                // Scale by the diagonal mass matrix.
                std::vector<dof_id_type> F_local_idxs(F_node_idxs.size());
                for (unsigned int i = 0; i < F_node_idxs.size(); ++i)
                {
                    F_local_idxs[i] = F_petsc_vec->map_global_to_local_index(F_node_idxs[i]);
                    TBOX_ASSERT(F_local_idxs[i] == dX_vec->map_global_to_local_index(F_node_idxs[i]));
                    F_node[i] *= dX_local_soln[F_local_idxs[i]];
                }

                // Insert the values into the global array.
                if (is_ghosted)
                {
                    for (unsigned int i = 0; i < F_node_idxs.size(); ++i)
                    {
                        F_local_soln[F_local_idxs[i]] += F_node[i];
                    }
                }
                else
                {
                    F_vec.add_vector(F_node, F_node_idxs);
                }
            }
        }

        // Restore local form vectors.
        dX_vec->restore_array();
    }
    else
    {
        // Loop over the patches to interpolate values to the element quadrature
        // points from the grid, then use these values to compute the projection
        // of the interpolated velocity field onto the FE basis functions.
        DenseVector<double> F_rhs;
        // Assemble F_rhs_e's vectors in an interleaved format (see the implementation):
        std::vector<double> F_rhs_concatenated;
        std::vector<double> F_qp, X_qp;
        std::vector<libMesh::dof_id_type> dof_id_scratch;
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            int local_patch_num = 0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
            {
                // The relevant collection of elements.
                const std::vector<Elem*>& patch_elems = d_active_patch_elem_map[ln][local_patch_num];
                const size_t num_active_patch_elems = patch_elems.size();
                if (!num_active_patch_elems) continue;

                const Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* const patch_dx = patch_geom->getDx();
                const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);

                // Determining which quadrature rule should be used on which
                // processor is surprisingly expensive, so cache the keys:
                std::vector<quad_key_type> quad_keys(num_active_patch_elems);

                // Cache interpolated positions too:
                std::vector<boost::multi_array<double, 2> > X_nodes(num_active_patch_elems);

                // Setup vectors to store the values of F and X at the quadrature
                // points.
                unsigned int n_qp_patch = 0;
                for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
                {
                    Elem* const elem = patch_elems[e_idx];
                    const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                    get_values_for_interpolation(X_nodes[e_idx], *X_petsc_vec, X_local_soln, X_dof_indices);
                    const quad_key_type key = getQuadratureKey(interp_spec.quad_type,
                                                               interp_spec.quad_order,
                                                               interp_spec.use_adaptive_quadrature,
                                                               interp_spec.point_density,
                                                               interp_spec.allow_rules_with_negative_weights,
                                                               elem,
                                                               X_nodes[e_idx],
                                                               patch_dx_min);
                    QBase& qrule = d_fe_data->d_quadrature_cache[key];
                    n_qp_patch += qrule.n_points();
                    quad_keys[e_idx] = key;
                }
                if (!n_qp_patch) continue;
                F_qp.resize(n_vars * n_qp_patch);
                X_qp.resize(NDIM * n_qp_patch);
                std::fill(F_qp.begin(), F_qp.end(), 0.0);

                // Loop over the elements and compute the positions of the
                // quadrature points.
                int qp_offset = 0;
                for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
                {
                    Elem* const elem = patch_elems[e_idx];
                    TBOX_ASSERT(elem->active());
                    const quad_key_type& key = quad_keys[e_idx];
                    const QBase& qrule = d_fe_data->d_quadrature_cache[key];
                    const FEBase& X_fe = X_fe_cache(key, elem);
                    const std::vector<std::vector<double> >& phi_X = X_fe.get_phi();

                    const unsigned int n_node = elem->n_nodes();
                    const unsigned int n_qp = qrule.n_points();
                    TBOX_ASSERT(n_qp == phi_X[0].size());
                    double* X_begin = &X_qp[NDIM * qp_offset];
                    std::fill(X_begin, X_begin + NDIM * n_qp, 0.0);
                    sum_weighted_elem_solution<true>(NDIM, n_node, qp_offset, phi_X, {}, X_nodes[e_idx], X_qp);
                    qp_offset += n_qp;
                }

                // Interpolate values from the Cartesian grid patch to the
                // quadrature points.
                //
                // NOTE: Values are interpolated only to those quadrature points
                // that are within the patch interior.
                const Box<NDIM>& interp_box = patch->getBox();
                Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_data_idx);
                if (cc_data)
                {
                    Pointer<CellData<NDIM, double> > f_cc_data = f_data;
                    LEInteractor::interpolate(
                        F_qp, n_vars, X_qp, NDIM, f_cc_data, patch, interp_box, interp_spec.kernel_fcn);
                }
                if (sc_data)
                {
                    Pointer<SideData<NDIM, double> > f_sc_data = f_data;
                    LEInteractor::interpolate(
                        F_qp, n_vars, X_qp, NDIM, f_sc_data, patch, interp_box, interp_spec.kernel_fcn);
                }

                // Loop over the elements and accumulate the right-hand-side values.
                qp_offset = 0;
                for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
                {
                    Elem* const elem = patch_elems[e_idx];
                    const auto& F_dof_indices = F_dof_map_cache.dof_indices(elem);
                    // check the concatenation assumption
#ifndef NDEBUG
                    for (unsigned int i = 0; i < n_vars; ++i)
                    {
                        TBOX_ASSERT(F_dof_indices[i].size() == F_dof_indices[0].size());
                    }
#endif
                    F_rhs_concatenated.resize(n_vars * F_dof_indices[0].size());
                    std::fill(F_rhs_concatenated.begin(), F_rhs_concatenated.end(), 0.0);
                    const quad_key_type& key = quad_keys[e_idx];
                    const FEBase& F_fe = F_fe_cache(key, elem);
                    const QBase& qrule = d_fe_data->d_quadrature_cache[key];

                    // JxW depends on the element
                    const std::vector<double>& JxW_F =
                        get_JxW(key, elem, is_volume_mesh, volume_mapping_cache, surface_mapping_cache);
                    const std::vector<std::vector<double> >& phi_F = F_fe.get_phi();

                    const unsigned int n_qp = qrule.n_points();
                    TBOX_ASSERT(n_qp == phi_F[0].size());
                    TBOX_ASSERT(n_qp == JxW_F.size());
                    const size_t n_basis = F_dof_indices[0].size();
                    integrate_elem_rhs(n_vars, n_basis, qp_offset, phi_F, JxW_F, F_qp, F_rhs_concatenated);

                    for (unsigned int var_n = 0; var_n < n_vars; ++var_n)
                    {
                        F_rhs.resize(F_dof_indices[var_n].size());
                        std::copy(F_rhs_concatenated.begin() + var_n * n_basis,
                                  F_rhs_concatenated.begin() + (var_n + 1) * n_basis,
                                  F_rhs.get_values().begin());

                        // We do *not* apply constraints here. See the note in the
                        // documentation of this function for an explanation.
                        if (is_ghosted)
                        {
                            for (unsigned int i = 0; i < F_dof_indices[var_n].size(); ++i)
                            {
                                const PetscInt index = F_petsc_vec->map_global_to_local_index(F_dof_indices[var_n][i]);
#ifndef NDEBUG
                                TBOX_ASSERT(0 <= index);
                                TBOX_ASSERT(index < F_local_size);
#endif
                                F_local_soln[index] += F_rhs(i);
                            }
                        }
                        else
                        {
                            copy_dof_ids_to_vector(var_n, F_dof_indices, dof_id_scratch);
                            F_vec.add_vector(F_rhs, dof_id_scratch);
                        }
                    }
                    qp_offset += n_qp;
                }
            }
        }
    }

    // Restore local form vectors.
    X_petsc_vec->restore_array();
    if (is_ghosted)
    {
        int ierr = VecRestoreArray(F_local_form, &F_local_soln);
        IBTK_CHKERRQ(ierr);
        ierr = VecGhostRestoreLocalForm(F_petsc_vec->vec(), &F_local_form);
        IBTK_CHKERRQ(ierr);

        if (close_F)
        {
            ierr = VecGhostUpdateBegin(F_petsc_vec->vec(), ADD_VALUES, SCATTER_REVERSE);
            IBTK_CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(F_petsc_vec->vec(), ADD_VALUES, SCATTER_REVERSE);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Accumulate data.
    if (close_F) F_vec.close();

    IBTK_TIMER_STOP(t_interp_weighted);
    return;
} // interpWeighted

void
FEDataManager::interp(const int f_data_idx,
                      NumericVector<double>& F_vec,
                      NumericVector<double>& X_vec,
                      const std::string& system_name,
                      const std::vector<Pointer<RefineSchedule<NDIM> > >& f_refine_scheds,
                      const double fill_data_time,
                      const bool close_X)
{
    interp(f_data_idx, F_vec, X_vec, system_name, d_default_interp_spec, f_refine_scheds, fill_data_time, close_X);
    return;
} // interp

void
FEDataManager::interp(const int f_data_idx,
                      NumericVector<double>& F_vec,
                      NumericVector<double>& X_vec,
                      const std::string& system_name,
                      const FEDataManager::InterpSpec& interp_spec,
                      const std::vector<Pointer<RefineSchedule<NDIM> > >& f_refine_scheds,
                      const double fill_data_time,
                      const bool close_X)
{
    IBTK_TIMER_START(t_interp);

    // Interpolate quantity at quadrature points and filter it to nodal points.
    std::unique_ptr<NumericVector<double> > F_rhs_vec = F_vec.zero_clone();
    interpWeighted(f_data_idx,
                   *F_rhs_vec,
                   X_vec,
                   system_name,
                   interp_spec,
                   f_refine_scheds,
                   fill_data_time,
                   /*close_F*/ true,
                   close_X);

    // Solve for the nodal values.
    computeL2Projection(
        F_vec, *F_rhs_vec, system_name, interp_spec.use_consistent_mass_matrix, true /*close_U*/, false /*close_F*/);

    IBTK_TIMER_STOP(t_interp);
    return;
} // interp

void
FEDataManager::restrictData(const int f_data_idx,
                            NumericVector<double>& F_vec,
                            NumericVector<double>& X_vec,
                            const std::string& system_name,
                            const bool use_consistent_mass_matrix,
                            const bool close_X)
{
    IBTK_TIMER_START(t_restrict_data);

    // NOTE #1: This routine is sepcialized for a staggered-grid Eulerian
    // discretization.  It should be straightforward to generalize it to work
    // with other data centerings.
    //
    // NOTE #2: This code is specialized for isoparametric elements.  It is less
    // clear how to relax this assumption.
    //
    // NOTE #3: This implementation uses the pointwise value of J = det(dX/ds)
    // to convert a Lagrangian density into an Eulerian density.  We should
    // investigate whether there is any advantage to using a projection of J
    // onto a (discontinuous) FE basis instead of evaluating J directly from the
    // discrete deformation.

    // Extract the mesh.
    const MeshBase& mesh = d_fe_data->d_es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    TBOX_ASSERT(dim == NDIM);

    // Extract the FE systems and DOF maps, and setup the FE object.
    System& F_system = d_fe_data->d_es->get_system(system_name);
    const unsigned int n_vars = F_system.n_vars();
    const DofMap& F_dof_map = F_system.get_dof_map();
    FEData::SystemDofMapCache& F_dof_map_cache = *getDofMapCache(system_name);
    System& X_system = d_fe_data->d_es->get_system(COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    FEData::SystemDofMapCache& X_dof_map_cache = *getDofMapCache(COORDINATES_SYSTEM_NAME);
    FEType F_fe_type = F_dof_map.variable_type(0);
    for (unsigned i = 0; i < n_vars; ++i)
    {
        TBOX_ASSERT(F_dof_map.variable_type(i) == F_fe_type);
    }
    FEType X_fe_type = X_dof_map.variable_type(0);
    for (unsigned d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
    }
    std::unique_ptr<FEBase> F_fe_autoptr(FEBase::build(dim, F_fe_type)), X_fe_autoptr;
    if (F_fe_type != X_fe_type)
    {
        X_fe_autoptr = std::unique_ptr<FEBase>(FEBase::build(dim, X_fe_type));
    }
    FEBase* F_fe = F_fe_autoptr.get();
    FEBase* X_fe = X_fe_autoptr.get() ? X_fe_autoptr.get() : F_fe_autoptr.get();
    const std::vector<std::vector<double> >& phi_F = F_fe->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& dphi_X = X_fe->get_dphi();

    // Communicate any unsynchronized ghost data and extract the underlying
    // solution data.
    if (close_X) X_vec.close();
    auto X_petsc_vec = static_cast<PetscVector<double>*>(&X_vec);
    const double* const X_local_soln = X_petsc_vec->get_array_read();

    // Loop over the patches to assemble the right-hand-side vector used to
    // solve for F.
    std::unique_ptr<NumericVector<double> > F_rhs_vec = F_vec.zero_clone();
    std::vector<DenseVector<double> > F_rhs_e(n_vars);
    TensorValue<double> dX_ds;
    boost::multi_array<double, 2> X_node;
    std::vector<libMesh::Point> s_node_cache, X_node_cache;
    Point X_min, X_max;
    std::vector<libMesh::Point> intersection_ref_coords;
    std::vector<SideIndex<NDIM> > intersection_indices;
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const IntVector<NDIM>& ratio = level->getRatio();
        const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            // The relevant collection of elements.
            const std::vector<Elem*>& patch_elems = d_active_patch_elem_map[ln][local_patch_num];
            const size_t num_active_patch_elems = patch_elems.size();
            if (!num_active_patch_elems) continue;

            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_data_idx);
            const Box<NDIM>& patch_box = patch->getBox();
            const CellIndex<NDIM>& patch_lower = patch_box.lower();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_dx = patch_geom->getDx();
            double dV = 1.0;
            for (unsigned int d = 0; d < NDIM; ++d) dV *= patch_dx[d];

            std::array<Box<NDIM>, NDIM> side_boxes;
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                side_boxes[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                if (!patch_geom->getTouchesRegularBoundary(axis, 1)) side_boxes[axis].growUpper(axis, -1);
            }

            SideData<NDIM, bool> interpolated_value_at_loc(patch_box, 1, IntVector<NDIM>(0));
            interpolated_value_at_loc.fillAll(false);

            // Loop over the elements.
            std::vector<libMesh::dof_id_type> dof_id_scratch;
            for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems[e_idx];
                const unsigned int n_node = elem->n_nodes();
                const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);

                // Cache the nodal and physical coordinates of the element,
                // determine the bounding box of the current configuration of the
                // element, and set the nodal coordinates of the element to
                // correspond to the physical coordinates.
                s_node_cache.resize(n_node);
                X_node_cache.resize(n_node);
                X_min = Point::Constant(std::numeric_limits<double>::max());
                X_max = Point::Constant(-std::numeric_limits<double>::max());
                for (unsigned int k = 0; k < n_node; ++k)
                {
                    s_node_cache[k] = elem->point(k);
                    for (int d = 0; d < NDIM; ++d)
                    {
                        X_node_cache[k](d) = X_vec(X_dof_indices[d][k]);
                        X_min[d] = std::min(X_min[d], X_node_cache[k](d));
                        X_max[d] = std::max(X_max[d], X_node_cache[k](d));
                    }
                    elem->point(k) = X_node_cache[k];
                }

                Box<NDIM> box(IndexUtilities::getCellIndex(&X_min[0], grid_geom, ratio),
                              IndexUtilities::getCellIndex(&X_max[0], grid_geom, ratio));
                box.grow(IntVector<NDIM>(1));
                box = box * patch_box;

                // Loop over coordinate directions and look for Eulerian grid points
                // that are covered by the element.
                intersection_ref_coords.clear();
                intersection_indices.clear();
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    // Loop over the relevant range of indices.
                    for (SideIterator<NDIM> b(box, axis); b; b++)
                    {
                        const SideIndex<NDIM>& i_s = b();
                        if (side_boxes[axis].contains(i_s) && !interpolated_value_at_loc(i_s))
                        {
                            libMesh::Point p;
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                p(d) = patch_x_lower[d] + patch_dx[d] * (static_cast<double>(i_s(d) - patch_lower[d]) +
                                                                         (d == axis ? 0.0 : 0.5));
                            }
                            static const double TOL = std::sqrt(std::numeric_limits<double>::epsilon());
                            const libMesh::Point ref_coords =
                                FEInterface::inverse_map(dim, X_fe_type, elem, p, TOL, false);
                            if (FEInterface::on_reference_element(ref_coords, elem->type(), TOL))
                            {
                                intersection_ref_coords.push_back(ref_coords);
                                intersection_indices.push_back(i_s);
                                interpolated_value_at_loc(i_s) = true;
                            }
                        }
                    }
                }

                // Restore the nodal coordinates.
                for (unsigned int k = 0; k < n_node; ++k)
                {
                    elem->point(k) = s_node_cache[k];
                }

                // If there are no intersection points, then continue on to the next
                // element.
                if (intersection_ref_coords.empty()) continue;

                // Evaluate the Eulerian value and rescale it by 1.0/det(dX/ds).
                F_fe->reinit(elem, &intersection_ref_coords);
                if (X_fe != F_fe) X_fe->reinit(elem, &intersection_ref_coords);
                const auto& F_dof_indices = F_dof_map_cache.dof_indices(elem);
                for (unsigned int i = 0; i < n_vars; ++i)
                {
                    F_rhs_e[i].resize(static_cast<int>(F_dof_indices[i].size()));
                }
                get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);
                const size_t n_basis = F_dof_indices[0].size();
                for (unsigned int qp = 0; qp < intersection_ref_coords.size(); ++qp)
                {
                    const SideIndex<NDIM>& i_s = intersection_indices[qp];
                    const int axis = i_s.getAxis();
                    jacobian(dX_ds, qp, X_node, dphi_X);
                    const double J = std::abs(dX_ds.det());
                    const double F_qp = (*f_data)(i_s)*dV / J;
                    for (unsigned int k = 0; k < n_basis; ++k)
                    {
                        F_rhs_e[axis](k) += F_qp * phi_F[k][qp];
                    }
                }
                for (unsigned int var_n = 0; var_n < n_vars; ++var_n)
                {
                    copy_dof_ids_to_vector(var_n, F_dof_indices, dof_id_scratch);
                    F_dof_map.constrain_element_vector(F_rhs_e[var_n], dof_id_scratch);
                    F_rhs_vec->add_vector(F_rhs_e[var_n], dof_id_scratch);
                }
            }
        }
    }

    X_petsc_vec->restore_array();

    // Solve for the nodal values.
    F_rhs_vec->close();
    computeL2Projection(F_vec, *F_rhs_vec, system_name, use_consistent_mass_matrix, true /*close_U*/, true /*close_F*/);

    IBTK_TIMER_STOP(t_restrict_data);
    return;
} // restrictData

NumericVector<double>*
FEDataManager::buildDiagonalL2MassMatrix(const std::string& system_name)
{
    return d_fe_projector->buildDiagonalL2MassMatrix(system_name);
} // buildDiagonalL2MassMatrix

PetscVector<double>*
FEDataManager::buildIBGhostedDiagonalL2MassMatrix(const std::string& system_name)
{
    if (!d_L2_proj_matrix_diag_ghost.count(system_name))
    {
        std::unique_ptr<PetscVector<double> > M_vec = buildIBGhostedVector(system_name);
        *M_vec = *d_fe_projector->buildDiagonalL2MassMatrix(system_name);
        M_vec->close();
        d_L2_proj_matrix_diag_ghost[system_name] = std::move(M_vec);
    }
    return d_L2_proj_matrix_diag_ghost[system_name].get();
}

bool
FEDataManager::computeL2Projection(NumericVector<double>& U_vec,
                                   NumericVector<double>& F_vec,
                                   const std::string& system_name,
                                   const bool consistent_mass_matrix,
                                   const bool close_U,
                                   const bool close_F,
                                   const double tol,
                                   const unsigned int max_its)
{
    return d_fe_projector->computeL2Projection(*static_cast<PetscVector<double>*>(&U_vec),
                                               *static_cast<PetscVector<double>*>(&F_vec),
                                               system_name,
                                               consistent_mass_matrix,
                                               close_U,
                                               close_F,
                                               tol,
                                               max_its);
} // computeL2Projection

bool
FEDataManager::updateQuadratureRule(std::unique_ptr<QBase>& qrule,
                                    QuadratureType type,
                                    Order order,
                                    bool use_adaptive_quadrature,
                                    double point_density,
                                    bool allow_rules_with_negative_weights,
                                    const Elem* const elem,
                                    const boost::multi_array<double, 2>& X_node,
                                    const double dx_min)
{
    unsigned int elem_dim = elem->dim();
    const unsigned int elem_p_level = elem->p_level();
    bool qrule_updated = false;

    ElemType elem_type;
    std::tie(elem_type, type, order, allow_rules_with_negative_weights) = getQuadratureKey(
        type, order, use_adaptive_quadrature, point_density, allow_rules_with_negative_weights, elem, X_node, dx_min);

    if (!qrule || qrule->type() != type || qrule->get_dim() != elem_dim || qrule->get_order() != order ||
        qrule->allow_rules_with_negative_weights != allow_rules_with_negative_weights ||
        qrule->get_elem_type() != elem_type || qrule->get_p_level() != elem_p_level)
    {
        qrule =
            (type == QGRID ? std::unique_ptr<QBase>(new QGrid(elem_dim, order)) : QBase::build(type, elem_dim, order));
        qrule->allow_rules_with_negative_weights = allow_rules_with_negative_weights;
        qrule->init(elem_type, elem_p_level);
        qrule_updated = true;
    }
    return qrule_updated;
}

bool
FEDataManager::updateInterpQuadratureRule(std::unique_ptr<QBase>& qrule,
                                          const FEDataManager::InterpSpec& spec,
                                          const Elem* const elem,
                                          const boost::multi_array<double, 2>& X_node,
                                          const double dx_min)
{
    return updateQuadratureRule(qrule,
                                spec.quad_type,
                                spec.quad_order,
                                spec.use_adaptive_quadrature,
                                spec.point_density,
                                spec.allow_rules_with_negative_weights,
                                elem,
                                X_node,
                                dx_min);
}

bool
FEDataManager::updateSpreadQuadratureRule(std::unique_ptr<QBase>& qrule,
                                          const FEDataManager::SpreadSpec& spec,
                                          const Elem* const elem,
                                          const boost::multi_array<double, 2>& X_node,
                                          const double dx_min)
{
    return updateQuadratureRule(qrule,
                                spec.quad_type,
                                spec.quad_order,
                                spec.allow_rules_with_negative_weights,
                                spec.point_density,
                                spec.use_adaptive_quadrature,
                                elem,
                                X_node,
                                dx_min);
}

void
FEDataManager::addWorkloadEstimate(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                   const int workload_data_idx,
                                   const int coarsest_ln_in,
                                   const int finest_ln_in)
{
    IBTK_TIMER_START(t_update_workload_estimates);

    const int coarsest_ln = (coarsest_ln_in == -1) ? getCoarsestPatchLevelNumber() : coarsest_ln_in;
    const int finest_ln = (finest_ln_in == -1) ? getFinestPatchLevelNumber() : finest_ln_in;
    TBOX_ASSERT(coarsest_ln >= getCoarsestPatchLevelNumber() && coarsest_ln <= getFinestPatchLevelNumber());
    TBOX_ASSERT(finest_ln >= getCoarsestPatchLevelNumber() && finest_ln <= getFinestPatchLevelNumber());

    {
        updateQuadPointCountData(0, hierarchy->getFinestLevelNumber());
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);
        hier_cc_data_ops.axpy(
            workload_data_idx, d_default_workload_spec.q_point_weight, d_qp_count_idx, workload_data_idx);
    }

    // Add work estimates from duplicated nodes.
    if (d_default_workload_spec.duplicated_node_weight != 0.0)
    {
        FEData::SystemDofMapCache& X_dof_map_cache = *getDofMapCache(COORDINATES_SYSTEM_NAME);

        // Extract the underlying solution data.
        NumericVector<double>* X_ghost_vec = buildGhostedCoordsVector();
        auto X_petsc_vec = static_cast<PetscVector<double>*>(X_ghost_vec);
        const double* const X_local_soln = X_petsc_vec->get_array_read();

        // Determine the position of the nodes of each element.
        boost::multi_array<double, 2> X_node;
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            const IntVector<NDIM>& ratio = level->getRatio();
            const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
            int local_patch_num = 0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
            {
                const std::vector<Elem*>& patch_elems = d_active_patch_elem_map[ln][local_patch_num];
                const size_t num_active_patch_elems = patch_elems.size();
                if (!num_active_patch_elems) continue;

                const Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, double> > workload_data = patch->getPatchData(workload_data_idx);

                for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
                {
                    Elem* const elem = patch_elems[e_idx];
                    if (!(getPatchLevel(elem) == ln)) continue;
                    const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                    get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);

                    const auto n_nodes = elem->n_nodes();
                    for (unsigned int node_n = 0; node_n < n_nodes; ++node_n)
                    {
                        Point node;
                        for (unsigned int d = 0; d < NDIM; ++d) node(d) = X_node[node_n][d];
                        const hier::Index<NDIM> i = IndexUtilities::getCellIndex(node, grid_geom, ratio);
                        if (patch_box.contains(i))
                        {
                            (*workload_data)(i) += d_default_workload_spec.duplicated_node_weight;
                        }
                    }
                }
            }
        }

        X_petsc_vec->restore_array();
    }

    IBTK_TIMER_STOP(t_update_workload_estimates);
    return;
} // addWorkloadEstimate

void
FEDataManager::applyGradientDetector(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                     const int level_number,
                                     const double /*error_data_time*/,
                                     const int tag_index,
                                     const bool initial_time,
                                     const bool /*uses_richardson_extrapolation_too*/)
{
    if (level_number >= d_max_level_number) return;

    IBTK_TIMER_START(t_apply_gradient_detector);

    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
    TBOX_ASSERT(hierarchy == d_hierarchy);

    if (initial_time)
    {
        // Determine the active elements associated with the prescribed patch
        // level.
        std::vector<std::vector<Elem*> > active_level_elem_map;
        const IntVector<NDIM> ghost_width = 1;
        collectActivePatchElements(active_level_elem_map, level_number, level_number + 1, d_max_level_number);
        std::vector<unsigned int> X_ghost_dofs;
        std::vector<Elem*> active_level_elems;
        collect_unique_elems(active_level_elems, active_level_elem_map);
        collectGhostDOFIndices(X_ghost_dofs, active_level_elems, COORDINATES_SYSTEM_NAME);

        // Extract the mesh.
        const MeshBase& mesh = d_fe_data->d_es->get_mesh();
        const Parallel::Communicator& comm = mesh.comm();
        const unsigned int dim = mesh.mesh_dimension();

        // Extract the FE system and DOF map, and setup the FE object.
        System& X_system = d_fe_data->d_es->get_system(COORDINATES_SYSTEM_NAME);
        const DofMap& X_dof_map = X_system.get_dof_map();
        FEData::SystemDofMapCache& X_dof_map_cache = *getDofMapCache(COORDINATES_SYSTEM_NAME);
        FEType fe_type = X_dof_map.variable_type(0);
        for (unsigned d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(X_dof_map.variable_type(d) == fe_type);
        }
        using quad_key_type = quadrature_key_type;
        FECache X_fe_cache(dim, fe_type, FEUpdateFlags::update_phi);

        // Setup and extract the underlying solution data.
        NumericVector<double>* X_vec = getCoordsVector();
        std::unique_ptr<NumericVector<double> > X_ghost_vec = NumericVector<double>::build(comm);
        X_ghost_vec->init(X_vec->size(), X_vec->local_size(), X_ghost_dofs, true, GHOSTED);
        copy_and_synch(*X_vec, *X_ghost_vec, /*close_v_in*/ false);
        auto X_petsc_vec = static_cast<PetscVector<double>*>(X_ghost_vec.get());
        const double* const X_local_soln = X_petsc_vec->get_array_read();

        // Tag cells for refinement whenever they contain active element
        // quadrature points.
        boost::multi_array<double, 2> X_node;
        Point X_qp;
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
        const IntVector<NDIM>& ratio = level->getRatio();
        const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            // The relevant collection of elements.
            const std::vector<Elem*>& patch_elems = active_level_elem_map[local_patch_num];
            const size_t num_active_patch_elems = patch_elems.size();
            if (!num_active_patch_elems) continue;

            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);

            Pointer<CellData<NDIM, int> > tag_data = patch->getPatchData(tag_index);

            for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems[e_idx];
                const int elem_ln = getPatchLevel(elem);
                TBOX_ASSERT(level_number + 1 <= elem_ln && elem_ln <= d_max_level_number);
                const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);

                const quad_key_type key = getQuadratureKey(d_default_interp_spec.quad_type,
                                                           d_default_interp_spec.quad_order,
                                                           d_default_interp_spec.use_adaptive_quadrature,
                                                           d_default_interp_spec.point_density,
                                                           d_default_interp_spec.allow_rules_with_negative_weights,
                                                           elem,
                                                           X_node,
                                                           patch_dx_min);
                const FEBase& X_fe = X_fe_cache(key, elem);
                const QBase& qrule = d_fe_data->d_quadrature_cache[key];
                const std::vector<std::vector<double> >& X_phi = X_fe.get_phi();
                TBOX_ASSERT(qrule.n_points() == X_phi[0].size());

                for (unsigned int qp = 0; qp < qrule.n_points(); ++qp)
                {
                    interpolate(&X_qp[0], qp, X_node, X_phi);

                    const hier::Index<NDIM> i = IndexUtilities::getCellIndex(X_qp, grid_geom, ratio);
                    tag_data->fill(1, Box<NDIM>(i - hier::Index<NDIM>(1), i + hier::Index<NDIM>(1)));
                }
            }
        }

        X_petsc_vec->restore_array();
    }
    else
    {
        // Tag all cells that correspond to elements that should be on finer
        // element levels.
        updateQuadPointCountData(level_number + 1, d_hierarchy->getFinestLevelNumber());
        const auto qp_scratch_idx = d_eulerian_data_cache->getCachedPatchDataIndex(d_qp_count_idx);
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(
            d_hierarchy, level_number, d_hierarchy->getFinestLevelNumber());
        hier_cc_data_ops.setToScalar(qp_scratch_idx, 0.0);

        // The semantics of d_qp_count_idx are slightly abused here: since we
        // need to tag cells for refinement whenever they contain a quadrature
        // point on a finer level, the values (after this loop) in that index
        // correspond to all quadrature points on finer levels.
        for (int finer_ln = d_hierarchy->getFinestLevelNumber(); finer_ln > level_number; --finer_ln)
        {
            Pointer<PatchLevel<NDIM> > finer_level = d_hierarchy->getPatchLevel(finer_ln);
            const int coarser_ln = finer_ln - 1;
            TBOX_ASSERT(coarser_ln >= level_number);
            Pointer<PatchLevel<NDIM> > coarser_level = d_hierarchy->getPatchLevel(coarser_ln);
            Pointer<CoarsenOperator<NDIM> > coarsen_op = new CartesianCellDoubleWeightedAverage<NDIM>();
            Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
            // Coarsen into the scratch index and then add that to the quadrature count.
            coarsen_alg->registerCoarsen(qp_scratch_idx, d_qp_count_idx, coarsen_op);
            coarsen_alg->createSchedule(coarser_level, finer_level)->coarsenData();
            hier_cc_data_ops.resetLevels(coarser_ln, coarser_ln);
            hier_cc_data_ops.add(d_qp_count_idx, d_qp_count_idx, qp_scratch_idx);
        }

        // Tag cells for refinement whenever they contain element quadrature
        // points.
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, int> > tag_data = patch->getPatchData(tag_index);
            Pointer<CellData<NDIM, double> > qp_count_data = patch->getPatchData(d_qp_count_idx);
            for (CellIterator<NDIM> b(patch_box); b; b++)
            {
                const CellIndex<NDIM>& i_c = b();
                if ((*qp_count_data)(i_c) > 0.0)
                {
                    (*tag_data)(i_c) = 1;
                }
            }
        }
    }

    IBTK_TIMER_STOP(t_apply_gradient_detector);
    return;
} // applyGradientDetector

void
FEDataManager::putToDatabase(Pointer<Database> db)
{
    IBTK_TIMER_START(t_put_to_database);
    d_fe_data->putToDatabase(db);

    TBOX_ASSERT(db);
    db->putInteger("FE_DATA_MANAGER_VERSION", FE_DATA_MANAGER_VERSION);

    IBTK_TIMER_STOP(t_put_to_database);
    return;
} // putToDatabase

void
FEDataManager::zeroExteriorValues(const CartesianPatchGeometry<NDIM>& patch_geom,
                                  const std::vector<double>& X_qp,
                                  std::vector<double>& F_qp,
                                  const int n_vars)
{
    // Check that we actually have input in the packed format
    std::size_t n_qp = X_qp.size() / NDIM;
    TBOX_ASSERT(n_qp * NDIM == X_qp.size());
    TBOX_ASSERT(F_qp.size() == n_vars * n_qp);
    libMeshWrappers::BoundingBox patch_box;

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        double min = patch_geom.getXLower()[d];
        double max = patch_geom.getXUpper()[d];
        patch_box.min()(d) = min;
        patch_box.max()(d) = max;
    }
    // fill extra dimension with 0.0, which is libMesh's convention
    for (unsigned int d = NDIM; d < LIBMESH_DIM; ++d)
    {
        patch_box.min()(d) = 0.0;
        patch_box.max()(d) = 0.0;
    }

    for (unsigned int i = 0; i < n_qp; ++i)
    {
        libMesh::Point qp;
        for (unsigned int d = 0; d < NDIM; ++d) qp(d) = X_qp[i * NDIM + d];

        bool contains_point = true;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            // note that this is a half-open interval
            if (qp(d) < patch_box.min()(d) || qp(d) >= patch_box.max()(d))
            {
                contains_point = false;
                break;
            }
        }

        if (!contains_point)
        {
            std::fill(F_qp.begin() + i * n_vars, F_qp.begin() + (i + 1) * n_vars, 0.0);
        }
    }
}

/////////////////////////////// PROTECTED ////////////////////////////////////

namespace
{
Pointer<Database>
setup_fe_projector_db(const Pointer<Database>& input_db)
{
    Pointer<Database> db;
    if (input_db->keyExists("FEProjector"))
    {
        db = input_db->getDatabase("FEProjector");
    }
    else
    {
        db = new InputDatabase("FEProjector");
    }
    return db;
}
} // namespace

FEDataManager::FEDataManager(std::string object_name,
                             const Pointer<Database>& input_db,
                             const int max_levels,
                             FEDataManager::InterpSpec default_interp_spec,
                             FEDataManager::SpreadSpec default_spread_spec,
                             FEDataManager::WorkloadSpec default_workload_spec,
                             IntVector<NDIM> ghost_width,
                             std::shared_ptr<SAMRAIDataCache> eulerian_data_cache,
                             std::shared_ptr<FEData> fe_data,
                             bool register_for_restart)
    : d_fe_data(fe_data),
      d_fe_projector(new FEProjector(d_fe_data, setup_fe_projector_db(input_db))),
      COORDINATES_SYSTEM_NAME(d_fe_data->d_coordinates_system_name),
      d_level_lookup(max_levels - 1,
                     collect_subdomain_ids(d_fe_data->getEquationSystems()->get_mesh()),
                     input_db ? (input_db->keyExists("subdomain_ids_on_levels") ?
                                     input_db->getDatabase("subdomain_ids_on_levels") :
                                     nullptr) :
                                nullptr),
      d_object_name(std::move(object_name)),
      d_registered_for_restart(register_for_restart),
      d_max_level_number(max_levels - 1),
      d_eulerian_data_cache(eulerian_data_cache),
      d_default_workload_spec(default_workload_spec),
      d_default_interp_spec(default_interp_spec),
      d_default_spread_spec(default_spread_spec),
      d_ghost_width(std::move(ghost_width))
{
    TBOX_ASSERT(!d_object_name.empty());

    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        FEDataManager::getFromRestart();
    }

    // Create/look up the variable context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name + "::CONTEXT");

    // Register the node count variable with the VariableDatabase.
    d_qp_count_var = new CellVariable<NDIM, double>(d_object_name + "::qp_count");
    d_qp_count_idx = var_db->registerVariableAndContext(d_qp_count_var, d_context, 0);

    // We cannot set up or check the Eulerian data cache yet since we don't have
    // the hierarchy - that will come later.

    // Setup the FE projector logging state.
    d_fe_projector->setLoggingEnabled(getLoggingEnabled());

    d_node_patch_check = string_to_enum<NodeOutsidePatchCheckType>(
        input_db->getStringWithDefault("node_outside_patch_check", "NODE_OUTSIDE_ERROR"));
    if (d_node_patch_check == UNKNOWN_NODE_OUTSIDE_PATCH_CHECK_TYPE)
    {
        TBOX_ERROR("unrecognized value " << input_db->getString("node_outside_patch_check")
                                         << "for input entry 'node_outside_patch_check'.");
    }

    // Setup Timers.
    IBTK_DO_ONCE(
        t_reinit_element_mappings =
            TimerManager::getManager()->getTimer("IBTK::FEDataManager::reinitElementMappings()");
        t_build_ghosted_solution_vector =
            TimerManager::getManager()->getTimer("IBTK::FEDataManager::buildGhostedSolutionVector()");
        t_build_ghosted_vector = TimerManager::getManager()->getTimer("IBTK::FEDataManager::buildGhostedVector()");
        t_spread = TimerManager::getManager()->getTimer("IBTK::FEDataManager::spread()");
        t_prolong_data = TimerManager::getManager()->getTimer("IBTK::FEDataManager::prolongData()");
        t_interp_weighted = TimerManager::getManager()->getTimer("IBTK::FEDataManager::interpWeighted()");
        t_interp = TimerManager::getManager()->getTimer("IBTK::FEDataManager::interp()");
        t_restrict_data = TimerManager::getManager()->getTimer("IBTK::FEDataManager::restrictData()");
        t_update_workload_estimates =
            TimerManager::getManager()->getTimer("IBTK::FEDataManager::updateWorkloadEstimates()");
        t_initialize_level_data = TimerManager::getManager()->getTimer("IBTK::FEDataManager::initializeLevelData()");
        t_reset_hierarchy_configuration =
            TimerManager::getManager()->getTimer("IBTK::FEDataManager::resetHierarchyConfiguration()");
        t_apply_gradient_detector =
            TimerManager::getManager()->getTimer("IBTK::FEDataManager::applyGradientDetector()");
        t_put_to_database = TimerManager::getManager()->getTimer("IBTK::FEDataManager::putToDatabase()");)
    return;
} // FEDataManager

FEDataManager::~FEDataManager()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
} // ~FEDataManager

void
FEDataManager::setLoggingEnabled(bool enable_logging)
{
    d_enable_logging = enable_logging;
    d_fe_projector->setLoggingEnabled(enable_logging);
    return;
} // setLoggingEnabled

bool
FEDataManager::getLoggingEnabled() const
{
    return d_enable_logging;
} // getLoggingEnabled

/////////////////////////////// PRIVATE //////////////////////////////////////

void
FEDataManager::updateQuadPointCountData(const int coarsest_ln, const int finest_ln)
{
    // The patch hierarchy should be completely set up at this point so we can
    // 1. assert that the Eulerian data cache is ready
    TBOX_ASSERT(d_eulerian_data_cache);
    TBOX_ASSERT(d_eulerian_data_cache->getCoarsestLevelNumber() == 0);
    TBOX_ASSERT(d_eulerian_data_cache->getFinestLevelNumber() >= getFinestPatchLevelNumber());

    // 2. actually allocate the Eulerian data we keep
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_qp_count_idx)) level->allocatePatchData(d_qp_count_idx);
    }

    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
    hier_cc_data_ops.setToScalar(d_qp_count_idx, 0.0);
    unsigned long n_local_q_points = 0;

    // Extract the mesh.
    const MeshBase& mesh = d_fe_data->d_es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // Extract the FE system and DOF map, and setup the FE object.
    System& X_system = d_fe_data->d_es->get_system(COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    FEData::SystemDofMapCache& X_dof_map_cache = *getDofMapCache(COORDINATES_SYSTEM_NAME);
    FEType fe_type = X_dof_map.variable_type(0);
    for (unsigned d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == fe_type);
    }

    // convenience alias for the quadrature key type used by FECache and MappingCache
    using quad_key_type = quadrature_key_type;
    FECache X_fe_cache(dim, fe_type, FEUpdateFlags::update_phi);

    // Extract the underlying solution data.
    NumericVector<double>* X_ghost_vec = buildGhostedCoordsVector();
    auto X_petsc_vec = static_cast<PetscVector<double>*>(X_ghost_vec);
    const double* const X_local_soln = X_petsc_vec->get_array_read();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const IntVector<NDIM>& ratio = level->getRatio();
        const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();

        // Determine the number of element quadrature points associated with
        // each Cartesian grid cell.
        boost::multi_array<double, 2> X_node;
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);

            Pointer<CellData<NDIM, double> > qp_count_data = patch->getPatchData(d_qp_count_idx);

            // nodal case:
            if (d_default_interp_spec.use_nodal_quadrature || d_default_spread_spec.use_nodal_quadrature)
            {
                // The relevant collection of nodes.
                const std::vector<Node*>& patch_nodes = d_active_patch_node_map[ln][local_patch_num];
                const size_t num_active_patch_nodes = patch_nodes.size();

                std::vector<dof_id_type> X_idxs;
                IBTK::Point X_qp;
                for (unsigned int k = 0; k < num_active_patch_nodes; ++k)
                {
                    const Node* const n = patch_nodes[k];
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        IBTK::get_nodal_dof_indices(X_dof_map, n, d, X_idxs);
                        X_qp[d] = X_local_soln[X_petsc_vec->map_global_to_local_index(X_idxs[0])];
                    }
                    const hier::Index<NDIM> i = IndexUtilities::getCellIndex(X_qp, grid_geom, ratio);
                    if (patch_box.contains(i))
                    {
                        (*qp_count_data)(i) += 1.0;
                        ++n_local_q_points;
                    }
                }
            }
            // non-nodal case:
            else if (!d_default_interp_spec.use_nodal_quadrature || !d_default_spread_spec.use_nodal_quadrature)
            {
                const std::vector<Elem*>& patch_elems = d_active_patch_elem_map[ln][local_patch_num];
                const size_t num_active_patch_elems = patch_elems.size();
                for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
                {
                    Elem* const elem = patch_elems[e_idx];
                    TBOX_ASSERT(getPatchLevel(elem) == ln);
                    const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                    get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);

                    const quad_key_type key = getQuadratureKey(d_default_interp_spec.quad_type,
                                                               d_default_interp_spec.quad_order,
                                                               d_default_interp_spec.use_adaptive_quadrature,
                                                               d_default_interp_spec.point_density,
                                                               d_default_interp_spec.allow_rules_with_negative_weights,
                                                               elem,
                                                               X_node,
                                                               patch_dx_min);
                    const FEBase& X_fe = X_fe_cache(key, elem);
                    const QBase& qrule = d_fe_data->d_quadrature_cache[key];
                    const std::vector<std::vector<double> >& X_phi = X_fe.get_phi();
                    TBOX_ASSERT(qrule.n_points() == X_phi[0].size());

                    Point X_qp;
                    for (unsigned int qp = 0; qp < qrule.n_points(); ++qp)
                    {
                        interpolate(&X_qp[0], qp, X_node, X_phi);
                        const hier::Index<NDIM> i = IndexUtilities::getCellIndex(X_qp, grid_geom, ratio);
                        if (patch_box.contains(i))
                        {
                            (*qp_count_data)(i) += 1.0;
                            ++n_local_q_points;
                        }
                    }
                }
            }
        }
    }

    if (d_enable_logging)
    {
        const int n_processes = IBTK_MPI::getNodes();
        const int current_rank = IBTK_MPI::getRank();
        const auto right_padding = std::size_t(std::log10(n_processes)) + 1;

        std::vector<unsigned long> n_q_points_on_processors(n_processes);
        n_q_points_on_processors[current_rank] = n_local_q_points;

        const int ierr = MPI_Allreduce(MPI_IN_PLACE,
                                       n_q_points_on_processors.data(),
                                       n_q_points_on_processors.size(),
                                       MPI_UNSIGNED_LONG,
                                       MPI_SUM,
                                       IBTK_MPI::getCommunicator());
        TBOX_ASSERT(ierr == 0);
        if (current_rank == 0)
        {
            for (int rank = 0; rank < n_processes; ++rank)
            {
                SAMRAI::tbox::plog << "quadrature points on processor " << std::setw(right_padding) << std::left << rank
                                   << " = " << n_q_points_on_processors[rank] << '\n';
            }
        }
    }

    X_petsc_vec->restore_array();
    return;
} // updateQuadPointCountData

void
FEDataManager::collectActivePatchElements(std::vector<std::vector<Elem*> >& active_patch_elems,
                                          const int level_number,
                                          const int coarsest_elem_ln,
                                          const int finest_elem_ln)
{
    // Get the necessary FE data.
    const MeshBase& mesh = d_fe_data->d_es->get_mesh();
    System& X_system = d_fe_data->d_es->get_system(COORDINATES_SYSTEM_NAME);

    // Setup data structures used to assign elements to patches.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    const int num_local_patches = level->getProcessorMapping().getNumberOfLocalIndices();
    std::vector<std::set<Elem*> > local_patch_elems(num_local_patches);
    active_patch_elems.resize(num_local_patches);

    // Try to exit quickly if no patches will actually have elements (i.e., if
    // there are no elements assigned to the range of levels we are working
    // with)
    bool levels_have_elements = false;
    for (int ln = coarsest_elem_ln; ln <= finest_elem_ln; ++ln)
    {
        levels_have_elements = levels_have_elements || d_level_lookup.levelHasElements(ln);
    }
    if (!levels_have_elements) return;

    // We associate an element with a Cartesian grid patch if the element's
    // bounding box (which is computed based on the bounds of quadrature
    // points) intersects the patch interior grown by
    // d_associated_elem_ghost_width (which is presently assumed to be 1).
    double dx_0 = std::numeric_limits<double>::max();
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        dx_0 = std::min(dx_0, *std::min_element(pgeom->getDx(), pgeom->getDx() + NDIM));
    }
    dx_0 = IBTK_MPI::minReduction(dx_0);
    TBOX_ASSERT(dx_0 != std::numeric_limits<double>::max());

    // be a bit paranoid by computing bounding boxes for elements as the union
    // of the bounding box of the nodes and the bounding box of the quadrature
    // points:
    const std::vector<libMeshWrappers::BoundingBox> local_nodal_bboxes =
        get_local_element_bounding_boxes(mesh, X_system);
    const std::vector<libMeshWrappers::BoundingBox> local_qp_bboxes =
        get_local_element_bounding_boxes(mesh,
                                         X_system,
                                         d_default_interp_spec.quad_type,
                                         d_default_interp_spec.quad_order,
                                         d_default_interp_spec.use_adaptive_quadrature,
                                         d_default_interp_spec.point_density,
                                         d_default_interp_spec.allow_rules_with_negative_weights,
                                         dx_0);
    TBOX_ASSERT(local_nodal_bboxes.size() == local_qp_bboxes.size());
    std::vector<libMeshWrappers::BoundingBox> local_bboxes;
    for (std::size_t box_n = 0; box_n < local_nodal_bboxes.size(); ++box_n)
    {
        local_bboxes.emplace_back(local_nodal_bboxes[box_n]);
#if LIBMESH_VERSION_LESS_THAN(1, 2, 0)
        // no BoundingBox::union_with in libMesh 1.1
        auto& box = local_bboxes.back();
        auto& other_box = local_qp_bboxes[box_n];
        for (unsigned int d = 0; d < LIBMESH_DIM; ++d)
        {
            box.first(d) = std::min(box.first(d), other_box.first(d));
            box.second(d) = std::max(box.second(d), other_box.second(d));
        }
#else
        local_bboxes.back().union_with(local_qp_bboxes[box_n]);
#endif
    }
    const std::vector<libMeshWrappers::BoundingBox> global_bboxes =
        get_global_element_bounding_boxes(mesh, local_bboxes);

    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        std::set<Elem*>& elems = local_patch_elems[local_patch_num];
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        // TODO: reimplement this with an rtree description of SAMRAI's patches
        libMeshWrappers::BoundingBox patch_bbox;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            patch_bbox.first(d) = pgeom->getXLower()[d] - dx[d] * d_associated_elem_ghost_width(d);
            patch_bbox.second(d) = pgeom->getXUpper()[d] + dx[d] * d_associated_elem_ghost_width(d);
        }
        for (unsigned int d = NDIM; d < LIBMESH_DIM; ++d)
        {
            patch_bbox.first(d) = 0.0;
            patch_bbox.second(d) = 0.0;
        }

        auto el_it = mesh.elements_begin();
        for (const libMeshWrappers::BoundingBox& bbox : global_bboxes)
        {
            if ((*el_it)->active())
            {
                const int elem_ln = getPatchLevel(*el_it);
                if (coarsest_elem_ln <= elem_ln && elem_ln <= finest_elem_ln)
                {
#if LIBMESH_VERSION_LESS_THAN(1, 6, 0)
                    if (bbox_intersects(bbox, patch_bbox)) elems.insert(*el_it);
#else
                    // New versions of libMesh have this function's performance
                    // problems fixed
                    if (bbox.intersects(patch_bbox)) elems.insert(*el_it);
#endif
                }
            }
            ++el_it;
        }
    }

    // Set the active patch element data.
    local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        const std::set<Elem*>& local_elems = local_patch_elems[local_patch_num];
        std::vector<Elem*>& active_elems = active_patch_elems[local_patch_num];
        active_elems.resize(local_elems.size());
        std::copy(local_elems.begin(), local_elems.end(), active_elems.begin());
    }
    return;
} // collectActivePatchElements

void
FEDataManager::collectActivePatchNodes(std::vector<std::vector<Node*> >& active_patch_nodes,
                                       const std::vector<std::vector<Elem*> >& active_patch_elems)
{
    const MeshBase& mesh = d_fe_data->d_es->get_mesh();
    const unsigned int num_local_patches = active_patch_elems.size();
    active_patch_nodes.resize(num_local_patches);
    for (unsigned int k = 0; k < num_local_patches; ++k)
    {
        std::set<dof_id_type> active_node_ids;
        for (const auto& elem : active_patch_elems[k])
        {
            for (unsigned int n = 0; n < elem->n_nodes(); ++n)
            {
                active_node_ids.insert(elem->node_id(n));
            }
        }
        const unsigned int num_active_nodes = active_node_ids.size();
        active_patch_nodes[k].reserve(num_active_nodes);
        for (const auto& active_node_id : active_node_ids)
        {
            active_patch_nodes[k].push_back(const_cast<Node*>(mesh.node_ptr(active_node_id)));
        }
    }
    return;
}

int
FEDataManager::getPatchLevel(const Elem* elem) const
{
    return d_level_lookup[elem->subdomain_id()];
}

void
FEDataManager::collectGhostDOFIndices(std::vector<unsigned int>& ghost_dofs,
                                      const std::vector<Elem*>& active_elems,
                                      const std::string& system_name)
{
    System& system = d_fe_data->d_es->get_system(system_name);
    const unsigned int sys_num = system.number();
    const DofMap& dof_map = system.get_dof_map();
    const unsigned int first_local_dof = dof_map.first_dof();
    const unsigned int end_local_dof = dof_map.end_dof();

    // Include non-local DOF constraint dependencies for local DOFs in the list
    // of ghost DOFs.
    std::vector<unsigned int> constraint_dependency_dof_list;
    for (auto i = dof_map.constraint_rows_begin(); i != dof_map.constraint_rows_end(); ++i)
    {
        const unsigned int constrained_dof = i->first;
        if (constrained_dof >= first_local_dof && constrained_dof < end_local_dof)
        {
            const DofConstraintRow& constraint_row = i->second;
            for (const auto& elem : constraint_row)
            {
                const unsigned int constraint_dependency = elem.first;
                if (constraint_dependency < first_local_dof || constraint_dependency >= end_local_dof)
                {
                    constraint_dependency_dof_list.push_back(constraint_dependency);
                }
            }
        }
    }

    // Record the local DOFs associated with the active local elements.
    std::set<unsigned int> ghost_dof_set(constraint_dependency_dof_list.begin(), constraint_dependency_dof_list.end());
    std::vector<unsigned int> dof_indices;
    for (unsigned int e = 0; e < active_elems.size(); ++e)
    {
        const Elem* const elem = active_elems[e];
        for (unsigned int var_num = 0; var_num < elem->n_vars(sys_num); ++var_num)
        {
            dof_map.dof_indices(elem, dof_indices, var_num);
            for (std::vector<unsigned int>::iterator it = dof_indices.begin(); it != dof_indices.end(); ++it)
            {
                const unsigned int dof_index = *it;
                if (dof_index < first_local_dof || dof_index >= end_local_dof)
                {
                    ghost_dof_set.insert(dof_index);
                }
            }
        }
    }
    ghost_dofs.clear();
    ghost_dofs.insert(ghost_dofs.end(), ghost_dof_set.begin(), ghost_dof_set.end());
    return;
} // collectGhostDOFIndices

void
FEDataManager::reinitializeIBGhostedDOFs(const std::string& system_name)
{
    // Reset the sets of dofs corresponding to IB ghost data. This is usually
    // a superset of the standard (i.e., all dofs on unowned cells adjacent to
    // locally owned cells) ghost data.
    if (!d_active_patch_ghost_dofs.count(system_name))
    {
        const System& system = d_fe_data->d_es->get_system(system_name);
        std::vector<libMesh::dof_id_type> ib_ghost_dofs;
        collectGhostDOFIndices(ib_ghost_dofs, d_active_elems, system_name);

        // Match the expected vector sizes by using the solution for non-ghost
        // sizes:
        const NumericVector<double>& solution = *system.solution;
        std::unique_ptr<PetscVector<double> > exemplar_ib_vector(new PetscVector<double>(
            system.comm(), solution.size(), solution.local_size(), ib_ghost_dofs, libMesh::GHOSTED));
        d_active_patch_ghost_dofs[system_name] = std::move(ib_ghost_dofs);
        d_system_ib_ghost_vec[system_name] = std::move(exemplar_ib_vector);
    }
}

void
FEDataManager::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();

    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to " << d_object_name << " not found in restart file.");
    }

    int ver = db->getInteger("FE_DATA_MANAGER_VERSION");
    if (ver != FE_DATA_MANAGER_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                                 << "Restart file version different than class version.");
    }

    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
