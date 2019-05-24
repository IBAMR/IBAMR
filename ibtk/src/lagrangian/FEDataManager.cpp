// Filename: FEDataManager.cpp
// Created on 19 Apr 2010 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianCellDoubleWeightedAverage.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenSchedule.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchData.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "ProcessorMapping.h"
#include "RefineSchedule.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "SideIterator.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "boost/multi_array.hpp"
#include "ibtk/FECache.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/FEMapCache.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_parallel_type.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_type.h"
#include "libmesh/fem_context.h"
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
#include "libmesh/string_to_enum.h"
#include "libmesh/system.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_vector.h"
#include "libmesh/variant_filter_iterator.h"
#include "petscksp.h"
#include "petscoptions.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

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
static Timer* t_spread;
static Timer* t_prolong_data;
static Timer* t_interp;
static Timer* t_interp_weighted;
static Timer* t_restrict_data;
static Timer* t_build_l2_projection_solver;
static Timer* t_build_diagonal_l2_mass_matrix;
static Timer* t_compute_l2_projection;
static Timer* t_update_workload_estimates;
static Timer* t_initialize_level_data;
static Timer* t_reset_hierarchy_configuration;
static Timer* t_apply_gradient_detector;
static Timer* t_put_to_database;

// Version of FEDataManager restart file data.
static const int FE_DATA_MANAGER_VERSION = 1;

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

inline boundary_id_type
get_dirichlet_bdry_ids(const std::vector<boundary_id_type>& bdry_ids)
{
    boundary_id_type dirichlet_bdry_ids = 0;
    for (const auto& bdry_id : bdry_ids)
    {
        if (bdry_id == FEDataManager::ZERO_DISPLACEMENT_X_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_Y_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_Z_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_XY_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_XZ_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_YZ_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_XYZ_BDRY_ID)
        {
            dirichlet_bdry_ids |= bdry_id;
        }
    }
    return dirichlet_bdry_ids;
} // get_dirichlet_bdry_ids
} // namespace

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
FEDataManager::getManager(const std::string& name,
                          const FEDataManager::InterpSpec& default_interp_spec,
                          const FEDataManager::SpreadSpec& default_spread_spec,
                          const FEDataManager::WorkloadSpec& default_workload_spec,
                          const IntVector<NDIM>& min_ghost_width,
                          bool register_for_restart)
{
    if (s_data_manager_instances.find(name) == s_data_manager_instances.end())
    {
        const IntVector<NDIM> ghost_width = IntVector<NDIM>::max(
            min_ghost_width,
            IntVector<NDIM>(std::max(LEInteractor::getMinimumGhostWidth(default_interp_spec.kernel_fcn),
                                     LEInteractor::getMinimumGhostWidth(default_spread_spec.kernel_fcn))));
        s_data_manager_instances[name] = new FEDataManager(
            name, default_interp_spec, default_spread_spec, default_workload_spec, ghost_width, register_for_restart);
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeAllManagers, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_data_manager_instances[name];
} // getManager

FEDataManager*
FEDataManager::getManager(const std::string& name,
                          const FEDataManager::InterpSpec& default_interp_spec,
                          const FEDataManager::SpreadSpec& default_spread_spec,
                          const IntVector<NDIM>& min_ghost_width,
                          bool register_for_restart)
{
    return getManager(name, default_interp_spec, default_spread_spec, {}, min_ghost_width, register_for_restart);
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
FEDataManager::registerLoadBalancer(Pointer<LoadBalancer<NDIM> > load_balancer, int workload_data_idx)
{
    IBTK_DEPRECATED_MEMBER_FUNCTION1("FEDataManager", "registerLoadBalancer");
    TBOX_ASSERT(load_balancer);
    d_load_balancer = load_balancer;
    d_workload_idx = workload_data_idx;
    return;
} // registerLoadBalancer

void
FEDataManager::setPatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    // Reset the hierarchy.
    TBOX_ASSERT(hierarchy);
    d_hierarchy = hierarchy;
    d_cached_eulerian_data.setPatchHierarchy(hierarchy);
    return;
} // setPatchHierarchy

Pointer<PatchHierarchy<NDIM> >
FEDataManager::getPatchHierarchy() const
{
    return d_hierarchy;
} // getPatchHierarchy

void
FEDataManager::setPatchLevels(const int coarsest_ln, const int finest_ln)
{
    // Reset the level numbers.
    TBOX_ASSERT(d_hierarchy);
    TBOX_ASSERT((coarsest_ln >= 0) && (finest_ln >= coarsest_ln) && (finest_ln <= d_hierarchy->getFinestLevelNumber()));
    d_coarsest_ln = coarsest_ln;
    d_finest_ln = finest_ln;
    d_cached_eulerian_data.resetLevels(coarsest_ln, finest_ln);
    return;
} // setPatchLevels

std::pair<int, int>
FEDataManager::getPatchLevels() const
{
    return std::make_pair(d_coarsest_ln, d_finest_ln + 1);
} // getPatchLevels

void
FEDataManager::setEquationSystems(EquationSystems* const equation_systems, const int level_number)
{
    d_es = equation_systems;
    d_level_number = level_number;

    // Now that we have the EquationSystems object we know the dimensionality
    // of the mesh.
    d_quadrature_cache = QuadratureCache(d_es->get_mesh().mesh_dimension());
    return;
} // setEquationSystems

EquationSystems*
FEDataManager::getEquationSystems() const
{
    return d_es;
} // getEquationSystems

FEDataManager::SystemDofMapCache*
FEDataManager::getDofMapCache(const std::string& system_name)
{
    TBOX_ASSERT(d_es);
    return getDofMapCache(d_es->get_system(system_name).number());
}

FEDataManager::SystemDofMapCache*
FEDataManager::getDofMapCache(unsigned int system_num)
{
    std::unique_ptr<SystemDofMapCache>& dof_map_cache = d_system_dof_map_cache[system_num];
    if (dof_map_cache == nullptr)
    {
        dof_map_cache.reset(new SystemDofMapCache(d_es->get_system(system_num)));
    }
    return dof_map_cache.get();
} // getDofMapCache

int
FEDataManager::getLevelNumber() const
{
    return d_level_number;
} // getLevelNumber

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
    return d_active_patch_elem_map;
} // getActivePatchElementMap

const std::vector<std::vector<Node*> >&
FEDataManager::getActivePatchNodeMap() const
{
    return d_active_patch_node_map;
} // getActivePatchNodeMap

void
FEDataManager::reinitElementMappings()
{
    IBTK_TIMER_START(t_reinit_element_mappings);

    // We reinitialize mappings after repartitioning, so clear the cache since
    // its content is no longer relevant:
    d_system_dof_map_cache.clear();

    // Delete cached hierarchy-dependent data.
    d_active_patch_elem_map.clear();
    d_active_patch_node_map.clear();
    d_active_patch_ghost_dofs.clear();
    d_active_elem_bboxes.clear();
    d_system_ghost_vec.clear();
    d_system_ib_ghost_vec.clear();

    // also clear any stored matrices:
    d_L2_proj_solver.clear();
    d_L2_proj_matrix.clear();
    d_L2_proj_matrix_diag.clear();

    // Reset the mappings between grid patches and active mesh
    // elements. collectActivePatchElements will populate d_active_elem_bboxes
    // and use it.
    collectActivePatchElements(d_active_patch_elem_map, d_level_number, d_ghost_width);
    collectActivePatchNodes(d_active_patch_node_map, d_active_patch_elem_map);

    // Reset the sets of dofs corresponding to IB ghost data. This is usually
    // a superset of the standard (i.e., all dofs on unowned cells adjacent to
    // locally owned cells) ghost data.
    std::vector<Elem*> active_elems;
    collect_unique_elems(active_elems, d_active_patch_elem_map);
    for (unsigned int system_n = 0; system_n < d_es->n_systems(); ++system_n)
    {
        const System& system = d_es->get_system(system_n);
        std::vector<libMesh::dof_id_type> ib_ghost_dofs;
        collectGhostDOFIndices(ib_ghost_dofs, active_elems, system.name());

        // Match the expected vector sizes by using the solution for non-ghost
        // sizes:
        const NumericVector<double>& solution = *system.solution;
        std::unique_ptr<PetscVector<double> > exemplar_ib_vector(new PetscVector<double>(
            system.comm(), solution.size(), solution.local_size(), ib_ghost_dofs, libMesh::GHOSTED));
        d_active_patch_ghost_dofs[system.name()] = std::move(ib_ghost_dofs);
        d_system_ib_ghost_vec[system.name()] = std::move(exemplar_ib_vector);
    }

    IBTK_TIMER_STOP(t_reinit_element_mappings);
    return;
} // reinitElementMappings

NumericVector<double>*
FEDataManager::getSolutionVector(const std::string& system_name) const
{
    return d_es->get_system(system_name).solution.get();
} // getSolutionVector

NumericVector<double>*
FEDataManager::buildGhostedSolutionVector(const std::string& system_name, const bool localize_data)
{
    IBTK_TIMER_START(t_build_ghosted_solution_vector);

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
FEDataManager::buildIBGhostedVector(const std::string& system_name) const
{
    TBOX_ASSERT(d_system_ib_ghost_vec.find(system_name) != d_system_ib_ghost_vec.end());
    const std::unique_ptr<PetscVector<double> >& exemplar_ib_vector = d_system_ib_ghost_vec.at(system_name);
    TBOX_ASSERT(exemplar_ib_vector);
    std::unique_ptr<NumericVector<double> > clone = exemplar_ib_vector->zero_clone();
    PetscVector<double>* ptr = dynamic_cast<PetscVector<double>*>(clone.release());
    TBOX_ASSERT(ptr);
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

    // Make a copy of the Eulerian data.
    const auto f_copy_data_idx = d_cached_eulerian_data.getCachedPatchDataIndex(f_data_idx);
    Pointer<HierarchyDataOpsReal<NDIM, double> > f_data_ops =
        HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(f_var, d_hierarchy, true);
    f_data_ops->swapData(f_copy_data_idx, f_data_idx);
    f_data_ops->setToScalar(f_data_idx, 0.0, /*interior_only*/ false);

    // We spread directly to the finest level of the patch hierarchy.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_number);

    // Extract the mesh.
    const MeshBase& mesh = d_es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // Extract the FE systems and DOF maps, and setup the FECache objects.
    System& F_system = d_es->get_system(system_name);
    const unsigned int n_vars = F_system.n_vars();
    const DofMap& F_dof_map = F_system.get_dof_map();
    SystemDofMapCache& F_dof_map_cache = *getDofMapCache(system_name);
    System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    SystemDofMapCache& X_dof_map_cache = *getDofMapCache(COORDINATES_SYSTEM_NAME);
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

    // convenience alias for the quadrature key type used by FECache and FEMapCache
    using quad_key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;
    FECache F_fe_cache(dim, F_fe_type);
    FECache X_fe_cache(dim, X_fe_type);
    FEMapCache fe_map_cache(dim);

    // Check to see if we are using nodal quadrature.
    const bool use_nodal_quadrature =
        spread_spec.use_nodal_quadrature && (F_fe_type == X_fe_type && F_order == X_order);

    // Communicate any unsynchronized ghost data.
    if (close_F) F_vec.close();
    if (close_X) X_vec.close();

    // We only use the FECache objects if we do *not* use nodal quadrature so
    // only perform sanity checks in that case:
    if (!use_nodal_quadrature)
    {
        // This will break, at some point in the future, if we ever use
        // nonnodal-interpolating finite elements. TODO: more finite elements
        // will probably work.
        std::vector<FEFamily> fe_family_whitelist{ LAGRANGE, L2_LAGRANGE };
        TBOX_ASSERT(std::find(fe_family_whitelist.begin(), fe_family_whitelist.end(), F_fe_type.family) !=
                    fe_family_whitelist.end());
        TBOX_ASSERT(std::find(fe_family_whitelist.begin(), fe_family_whitelist.end(), X_fe_type.family) !=
                    fe_family_whitelist.end());
    }

    if (use_nodal_quadrature)
    {
        // Store a copy of the F vector.
        std::unique_ptr<NumericVector<double> > F_vec_bak = F_vec.clone();
        *F_vec_bak = F_vec;

        // Multiply by the nodal volume fractions (to convert densities into
        // values).
        NumericVector<double>& F_dX_vec = F_vec;
        const NumericVector<double>& dX_vec = *buildDiagonalL2MassMatrix(system_name);
        F_dX_vec.pointwise_mult(F_vec, dX_vec);
        F_dX_vec.close();

        // Extract local form vectors.
        auto F_dX_petsc_vec = static_cast<PetscVector<double>*>(&F_dX_vec);
        const double* const F_dX_local_soln = F_dX_petsc_vec->get_array_read();

        auto X_petsc_vec = static_cast<PetscVector<double>*>(&X_vec);
        const double* const X_local_soln = X_petsc_vec->get_array_read();

        // Spread from the nodes.
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            // The relevant collection of nodes.
            const std::vector<Node*>& patch_nodes = d_active_patch_node_map[local_patch_num];
            const size_t num_active_patch_nodes = patch_nodes.size();
            if (!num_active_patch_nodes) continue;

            // Store the values of F_JxW and X at the nodes.
            std::vector<double> F_dX_node, X_node;
            F_dX_node.reserve(n_vars * num_active_patch_nodes);
            X_node.reserve(NDIM * num_active_patch_nodes);
            std::vector<dof_id_type> F_idxs, X_idxs;
            for (unsigned int k = 0; k < num_active_patch_nodes; ++k)
            {
                const Node* const n = patch_nodes[k];
                for (unsigned int i = 0; i < n_vars; ++i)
                {
                    IBTK::get_nodal_dof_indices(F_dof_map, n, i, F_idxs);
                    for (const auto& F_idx : F_idxs)
                    {
                        F_dX_node.push_back(F_dX_local_soln[F_dX_petsc_vec->map_global_to_local_index(F_idx)]);
                    }
                }
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    IBTK::get_nodal_dof_indices(X_dof_map, n, d, X_idxs);
                    for (const auto& X_idx : X_idxs)
                    {
                        X_node.push_back(X_local_soln[X_petsc_vec->map_global_to_local_index(X_idx)]);
                    }
                }
            }
            TBOX_ASSERT(F_dX_node.size() == n_vars * num_active_patch_nodes);
            TBOX_ASSERT(X_node.size() == NDIM * num_active_patch_nodes);

            // Spread values from the nodes to the Cartesian grid patch.
            //
            // NOTE: Values are spread only from those nodes that are within the
            // ghost cell width of the patch interior.
            //
            // \todo Fix this for FE structures with periodic boundaries.  Nodes
            // on periodic boundaries will be "double counted".
            //
            // \todo Add warnings for FE structures with periodic boundaries.
            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM> spread_box = Box<NDIM>::grow(patch->getBox(), d_ghost_width);
            Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_data_idx);
            if (cc_data)
            {
                Pointer<CellData<NDIM, double> > f_cc_data = f_data;
                LEInteractor::spread(
                    f_cc_data, F_dX_node, n_vars, X_node, NDIM, patch, spread_box, spread_spec.kernel_fcn);
            }
            if (sc_data)
            {
                Pointer<SideData<NDIM, double> > f_sc_data = f_data;
                LEInteractor::spread(
                    f_sc_data, F_dX_node, n_vars, X_node, NDIM, patch, spread_box, spread_spec.kernel_fcn);
            }
            if (f_phys_bdry_op)
            {
                f_phys_bdry_op->setPatchDataIndex(f_data_idx);
                f_phys_bdry_op->accumulateFromPhysicalBoundaryData(*patch, fill_data_time, f_data->getGhostCellWidth());
            }
        }

        // Restore local form vectors.
        F_dX_petsc_vec->restore_array();
        X_petsc_vec->restore_array();

        // Restore the value of the F vector.
        F_vec = *F_vec_bak;
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
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            // The relevant collection of elements.
            const std::vector<Elem*>& patch_elems = d_active_patch_elem_map[local_patch_num];
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
                const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                get_values_for_interpolation(X_nodes[e_idx], *X_petsc_vec, X_local_soln, X_dof_indices);
                const quad_key_type key = getQuadratureKey(spread_spec.quad_type,
                                                           spread_spec.quad_order,
                                                           spread_spec.use_adaptive_quadrature,
                                                           spread_spec.point_density,
                                                           elem,
                                                           X_nodes[e_idx],
                                                           patch_dx_min);
                quad_keys[e_idx] = key;
                QBase& qrule = d_quadrature_cache[key];
                n_qp_patch += qrule.n_points();
            }
            if (!n_qp_patch) continue;
            F_JxW_qp.resize(n_vars * n_qp_patch);
            X_qp.resize(NDIM * n_qp_patch);

            // Loop over the elements and compute the values to be spread and
            // the positions of the quadrature points.
            int qp_offset = 0;
            std::set<quad_key_type> used_quadratures;
            for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems[e_idx];
                const auto& F_dof_indices = F_dof_map_cache.dof_indices(elem);
                get_values_for_interpolation(F_node, *F_petsc_vec, F_local_soln, F_dof_indices);
                const quad_key_type& key = quad_keys[e_idx];
                FEBase& X_fe = X_fe_cache[key];
                FEBase& F_fe = F_fe_cache[key];
                FEMap& fe_map = fe_map_cache[key];
                QBase& qrule = d_quadrature_cache[key];

                // See the note in interpWeighted to explain why we override
                // libMesh's reinit logic here
                if (used_quadratures.find(key) == used_quadratures.end())
                {
                    F_fe.get_phi();
                    F_fe.reinit(elem);
                    X_fe.get_phi();
                    X_fe.reinit(elem);
                    used_quadratures.insert(key);
                }

                // JxW depends on the element
                fe_map.compute_map(dim, qrule.get_weights(), elem, /*second derivatives*/ false);
                const std::vector<double>& JxW_F = fe_map.get_JxW();
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

            // Spread values from the quadrature points to the Cartesian grid
            // patch.
            //
            // NOTE: Values are spread only from those quadrature points that
            // are within the ghost cell width of the patch interior.
            const Box<NDIM> spread_box = Box<NDIM>::grow(patch->getBox(), d_ghost_width);
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
            if (f_phys_bdry_op)
            {
                f_phys_bdry_op->setPatchDataIndex(f_data_idx);
                f_phys_bdry_op->accumulateFromPhysicalBoundaryData(*patch, fill_data_time, f_data->getGhostCellWidth());
            }
        }

        // Restore local form vectors.
        F_petsc_vec->restore_array();
        X_petsc_vec->restore_array();
    }

    // Accumulate data.
    f_data_ops->swapData(f_copy_data_idx, f_data_idx);
    f_data_ops->add(f_data_idx, f_data_idx, f_copy_data_idx);

    IBTK_TIMER_STOP(t_spread);
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

    // NOTE #1: This routine is sepcialized for a staggered-grid Eulerian
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
    const MeshBase& mesh = d_es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    TBOX_ASSERT(dim == NDIM);

    // Extract the FE systems and DOF maps, and setup the FE object.
    System& F_system = d_es->get_system(system_name);
    const unsigned int n_vars = F_system.n_vars();
    TBOX_ASSERT(n_vars == NDIM); // specialized to side-centered data
    const DofMap& F_dof_map = F_system.get_dof_map();
    SystemDofMapCache& F_dof_map_cache = *getDofMapCache(system_name);
    System& X_system = d_es->get_system(system_name);
    const DofMap& X_dof_map = X_system.get_dof_map();
    SystemDofMapCache& X_dof_map_cache = *getDofMapCache(COORDINATES_SYSTEM_NAME);
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
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_number);
    const IntVector<NDIM>& ratio = level->getRatio();
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const std::vector<Elem*>& patch_elems = d_active_patch_elem_map[local_patch_num];
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
                        const libMesh::Point ref_coords = FEInterface::inverse_map(dim, X_fe_type, elem, p, TOL, false);
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

    // We interpolate directly from the finest level of the patch hierarchy.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_number);

    // Extract the mesh.
    const MeshBase& mesh = d_es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // Extract the FE systems and DOF maps, and setup the FECache objects.
    System& F_system = d_es->get_system(system_name);
    const unsigned int n_vars = F_system.n_vars();
    const DofMap& F_dof_map = F_system.get_dof_map();
    SystemDofMapCache& F_dof_map_cache = *getDofMapCache(system_name);
    System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    SystemDofMapCache& X_dof_map_cache = *getDofMapCache(COORDINATES_SYSTEM_NAME);
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

    // convenience alias for the quadrature key type used by FECache and FEMapCache
    using quad_key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;
    FECache F_fe_cache(dim, F_fe_type);
    FECache X_fe_cache(dim, X_fe_type);
    FEMapCache fe_map_cache(dim);

    // Communicate any unsynchronized ghost data.
    for (const auto& f_refine_sched : f_refine_scheds)
    {
        if (f_refine_sched) f_refine_sched->fillData(fill_data_time);
    }

    if (close_X) X_vec.close();

    // Check to see if we are using nodal quadrature.
    const bool use_nodal_quadrature =
        interp_spec.use_nodal_quadrature && (F_fe_type == X_fe_type && F_order == X_order);

    // We only use the FECache objects if we do *not* use nodal quadrature so
    // only perform sanity checks in that case:
    if (!use_nodal_quadrature)
    {
        // This will break, at some point in the future, if we ever use
        // nonnodal-interpolating finite elements. TODO: more finite elements
        // will probably work.
        std::vector<FEFamily> fe_family_whitelist{ LAGRANGE, L2_LAGRANGE };
        TBOX_ASSERT(std::find(fe_family_whitelist.begin(), fe_family_whitelist.end(), F_fe_type.family) !=
                    fe_family_whitelist.end());
        TBOX_ASSERT(std::find(fe_family_whitelist.begin(), fe_family_whitelist.end(), X_fe_type.family) !=
                    fe_family_whitelist.end());
    }

    if (use_nodal_quadrature)
    {
        // Extract the local form vectors.
        auto X_petsc_vec = static_cast<PetscVector<double>*>(&X_vec);
        const double* const X_local_soln = X_petsc_vec->get_array_read();

        // Interpolate to the nodes.
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            // The relevant collection of nodes.
            const std::vector<Node*>& patch_nodes = d_active_patch_node_map[local_patch_num];
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
            for (unsigned int k = 0; k < num_active_patch_nodes; ++k)
            {
                const Node* const n = patch_nodes[k];
                IBTK::Point X;
                bool inside_patch = true;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    IBTK::get_nodal_dof_indices(X_dof_map, n, d, X_idxs);
                    X[d] = X_local_soln[X_petsc_vec->map_global_to_local_index(X_idxs[0])];
                    inside_patch = inside_patch && (X[d] >= patch_x_lower[d]) &&
                                   ((touches_upper_regular_bdry[d] && X[d] <= patch_x_upper[d]) ||
                                    (!touches_upper_regular_bdry[d] && X[d] < patch_x_upper[d]));
                }
                if (inside_patch)
                {
                    F_node.resize(F_node.size() + n_vars);
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

            // Insesrt the values of F at the nodes.
            F_vec.insert(F_node, F_node_idxs);
        }

        // Restore local form vectors.
        X_petsc_vec->restore_array();

        // Scale by the diagonal mass matrix.
        F_vec.close();
        const NumericVector<double>& dX_vec = *buildDiagonalL2MassMatrix(system_name);
        F_vec.pointwise_mult(F_vec, dX_vec);
        if (close_F) F_vec.close();
    }
    else
    {
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
        const bool is_ghosted = F_vec.type() == GHOSTED;
        if (is_ghosted)
        {
            TBOX_ASSERT(F_petsc_vec != nullptr);
            int ierr = VecGhostGetLocalForm(F_petsc_vec->vec(), &F_local_form);
            IBTK_CHKERRQ(ierr);
            ierr = VecGetArray(F_local_form, &F_local_soln);
            IBTK_CHKERRQ(ierr);
        }

        // Loop over the patches to interpolate values to the element quadrature
        // points from the grid, then use these values to compute the projection
        // of the interpolated velocity field onto the FE basis functions.
        DenseVector<double> F_rhs;
        // Assemble F_rhs_e's vectors in an interleaved format (see the implementation):
        std::vector<double> F_rhs_concatenated;
        std::vector<double> F_qp, X_qp;
        int local_patch_num = 0;
        std::vector<libMesh::dof_id_type> dof_id_scratch;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            // The relevant collection of elements.
            const std::vector<Elem*>& patch_elems = d_active_patch_elem_map[local_patch_num];
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
                                                           elem,
                                                           X_nodes[e_idx],
                                                           patch_dx_min);
                QBase& qrule = d_quadrature_cache[key];
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
            std::set<quad_key_type> used_X_quadratures;
            for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems[e_idx];
                const quad_key_type& key = quad_keys[e_idx];
                QBase& qrule = d_quadrature_cache[key];
                FEBase& X_fe = X_fe_cache[key];

                // libMesh::FE defaults to recalculating *everything* when we
                // call reinit unless we first request things (see
                // FEGenericBase<>::determine_calculations() for the
                // reinitialization logic). This only works for Lagrange
                // elements (and similar nodal value elements).
                //
                // Furthermore: we can skip calling reinit() completely as
                // long as we have called it on this FE object at least once.
                if (used_X_quadratures.find(key) == used_X_quadratures.end())
                {
                    // get_phi changes the state of the object if called
                    // before reinit
                    X_fe.get_phi();
                    X_fe.reinit(elem);
                    used_X_quadratures.insert(key);
                }

                const std::vector<std::vector<double> >& phi_X = X_fe.get_phi();

                const unsigned int n_node = elem->n_nodes();
                const unsigned int n_qp = qrule.n_points();
                TBOX_ASSERT(n_qp == phi_X[0].size());
                double* X_begin = &X_qp[NDIM * qp_offset];
                std::fill(X_begin, X_begin + NDIM * n_qp, 0.0);
                for (unsigned int k = 0; k < n_node; ++k)
                {
                    for (unsigned int qp = 0; qp < n_qp; ++qp)
                    {
                        const double& p_X = phi_X[k][qp];
                        for (unsigned int i = 0; i < NDIM; ++i)
                        {
                            X_qp[NDIM * (qp_offset + qp) + i] += X_nodes[e_idx][k][i] * p_X;
                        }
                    }
                }
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
            std::set<quad_key_type> used_F_quadratures;
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
                FEBase& F_fe = F_fe_cache[key];
                FEMap& fe_map = fe_map_cache[key];
                QBase& qrule = d_quadrature_cache[key];

                // Like above: conditionally initialize the FE object if it is
                // new
                if (used_F_quadratures.find(key) == used_F_quadratures.end())
                {
                    // Same as above
                    F_fe.get_phi();
                    F_fe.reinit(elem);
                    used_F_quadratures.insert(key);
                }

                // JxW depends on the element
                fe_map.compute_map(dim, qrule.get_weights(), elem, /*second derivatives*/ false);
                const std::vector<double>& JxW_F = fe_map.get_JxW();
                const std::vector<std::vector<double> >& phi_F = F_fe.get_phi();

                const unsigned int n_qp = qrule.n_points();
                TBOX_ASSERT(n_qp == phi_F[0].size());
                TBOX_ASSERT(n_qp == JxW_F.size());
                const size_t n_basis = F_dof_indices[0].size();
                integrate_elem_rhs(n_vars, n_basis, qp_offset, phi_F, JxW_F, F_qp, F_rhs_concatenated);

                for (unsigned int i = 0; i < n_vars; ++i)
                {
                    // libMesh sometimes resizes F_rhs inside
                    // constrain_element_vector, so ensure it has the right size:
                    F_rhs.resize(F_dof_indices[i].size());
                    std::copy(F_rhs_concatenated.begin() + i * n_basis,
                              F_rhs_concatenated.begin() + (i + 1) * n_basis,
                              F_rhs.get_values().begin());

                    dof_id_scratch = F_dof_indices[i];
                    F_dof_map.constrain_element_vector(F_rhs, dof_id_scratch);
                    if (is_ghosted)
                    {
                        for (unsigned int i = 0; i < dof_id_scratch.size(); ++i)
                        {
                            F_local_soln[F_petsc_vec->map_global_to_local_index(dof_id_scratch[i])] += F_rhs(i);
                        }
                    }
                    else
                    {
                        F_vec.add_vector(F_rhs, dof_id_scratch);
                    }
                }
                qp_offset += n_qp;
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
    const MeshBase& mesh = d_es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    TBOX_ASSERT(dim == NDIM);

    // Extract the FE systems and DOF maps, and setup the FE object.
    System& F_system = d_es->get_system(system_name);
    const unsigned int n_vars = F_system.n_vars();
    const DofMap& F_dof_map = F_system.get_dof_map();
    SystemDofMapCache& F_dof_map_cache = *getDofMapCache(system_name);
    System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    SystemDofMapCache& X_dof_map_cache = *getDofMapCache(COORDINATES_SYSTEM_NAME);
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
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_number);
    const IntVector<NDIM>& ratio = level->getRatio();
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const std::vector<Elem*>& patch_elems = d_active_patch_elem_map[local_patch_num];
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
                        const libMesh::Point ref_coords = FEInterface::inverse_map(dim, X_fe_type, elem, p, TOL, false);
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
            for (unsigned int i = 0; i < n_vars; ++i)
            {
                dof_id_scratch = F_dof_indices[i];
                F_dof_map.constrain_element_vector(F_rhs_e[i], dof_id_scratch);
                F_rhs_vec->add_vector(F_rhs_e[i], dof_id_scratch);
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

std::pair<LinearSolver<double>*, SparseMatrix<double>*>
FEDataManager::buildL2ProjectionSolver(const std::string& system_name)
{
    IBTK_TIMER_START(t_build_l2_projection_solver);

    if (!d_L2_proj_solver.count(system_name) || !d_L2_proj_matrix.count(system_name))
    {
        if (d_enable_logging)
        {
            plog << "FEDataManager::buildL2ProjectionSolver(): building L2 projection solver for system: "
                 << system_name << "\n";
        }

        // Extract the mesh.
        const MeshBase& mesh = d_es->get_mesh();
        const Parallel::Communicator& comm = mesh.comm();
        const unsigned int dim = mesh.mesh_dimension();

        // Extract the FE system and DOF map, and setup the FE object.
        System& system = d_es->get_system(system_name);
        const int sys_num = system.number();
        DofMap& dof_map = system.get_dof_map();
        SystemDofMapCache& dof_map_cache = *getDofMapCache(system_name);
        dof_map.compute_sparsity(mesh);
        FEType fe_type = dof_map.variable_type(0);
        std::unique_ptr<QBase> qrule = fe_type.default_quadrature_rule(dim);
        std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<double>& JxW = fe->get_JxW();
        const std::vector<std::vector<double> >& phi = fe->get_phi();

        // Build solver components.
        std::unique_ptr<LinearSolver<double> > solver = LinearSolver<double>::build(comm);
        solver->init();

        std::unique_ptr<SparseMatrix<double> > M_mat = SparseMatrix<double>::build(comm);
        M_mat->attach_dof_map(dof_map);
        M_mat->init();

        // Loop over the mesh to construct the system matrix.
        DenseMatrix<double> M_e;
        std::vector<libMesh::dof_id_type> dof_id_scratch;
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            fe->reinit(elem);

            const auto& dof_indices = dof_map_cache.dof_indices(elem);
            for (unsigned int var_num = 0; var_num < dof_map.n_variables(); ++var_num)
            {
                const auto& dof_indices_var = dof_indices[var_num];
                const auto dof_indices_sz = static_cast<unsigned int>(dof_indices_var.size());
                M_e.resize(dof_indices_sz, dof_indices_sz);
                const size_t n_basis = dof_indices_var.size();
                const unsigned int n_qp = qrule->n_points();
                for (unsigned int i = 0; i < n_basis; ++i)
                {
                    for (unsigned int j = 0; j < n_basis; ++j)
                    {
                        for (unsigned int qp = 0; qp < n_qp; ++qp)
                        {
                            M_e(i, j) += (phi[i][qp] * phi[j][qp]) * JxW[qp];
                        }
                    }
                }
                dof_id_scratch = dof_indices_var;
                dof_map.constrain_element_matrix(M_e, dof_id_scratch);
                M_mat->add_matrix(M_e, dof_id_scratch);
            }
        }

        // Flush assemble the matrix.
        M_mat->close();

        // Reset values at Dirichlet boundaries.
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                if (elem->neighbor_ptr(side)) continue;
                static const boundary_id_type dirichlet_bdry_id_set[3] = { ZERO_DISPLACEMENT_X_BDRY_ID,
                                                                           ZERO_DISPLACEMENT_Y_BDRY_ID,
                                                                           ZERO_DISPLACEMENT_Z_BDRY_ID };
                std::vector<boundary_id_type> bdry_ids;
                mesh.boundary_info->boundary_ids(elem, side, bdry_ids);
                const boundary_id_type dirichlet_bdry_ids = get_dirichlet_bdry_ids(bdry_ids);
                if (!dirichlet_bdry_ids) continue;
                fe->reinit(elem);
                for (unsigned int n = 0; n < elem->n_nodes(); ++n)
                {
                    if (elem->is_node_on_side(n, side))
                    {
                        const Node* const node = elem->node_ptr(n);
                        const auto& dof_indices = dof_map_cache.dof_indices(elem);
                        for (unsigned int var_num = 0; var_num < dof_map.n_variables(); ++var_num)
                        {
                            const unsigned int n_comp = node->n_comp(sys_num, var_num);
                            for (unsigned int comp = 0; comp < n_comp; ++comp)
                            {
                                if (!(dirichlet_bdry_ids & dirichlet_bdry_id_set[comp])) continue;
                                const unsigned int node_dof_index = node->dof_number(sys_num, var_num, comp);
                                if (!dof_map.is_constrained_dof(node_dof_index)) continue;
                                for (const auto& idx : dof_indices[var_num])
                                {
                                    M_mat->set(node_dof_index, idx, (node_dof_index == idx ? 1.0 : 0.0));
                                    M_mat->set(idx, node_dof_index, (node_dof_index == idx ? 1.0 : 0.0));
                                }
                            }
                        }
                    }
                }
            }
        }

        // Assemble the matrix.
        M_mat->close();

        // Setup the solver.
        solver->reuse_preconditioner(true);

        // Store the solver, mass matrix, and configuration options.
        d_L2_proj_solver[system_name] = std::move(solver);
        d_L2_proj_matrix[system_name] = std::move(M_mat);
    }

    IBTK_TIMER_STOP(t_build_l2_projection_solver);
    return std::make_pair(d_L2_proj_solver[system_name].get(), d_L2_proj_matrix[system_name].get());
} // buildL2ProjectionSolver

NumericVector<double>*
FEDataManager::buildDiagonalL2MassMatrix(const std::string& system_name)
{
    IBTK_TIMER_START(t_build_diagonal_l2_mass_matrix);

    if (!d_L2_proj_matrix_diag.count(system_name))
    {
        if (d_enable_logging)
        {
            plog << "FEDataManager::buildDiagonalL2MassMatrix(): building diagonal L2 mass matrix for system: "
                 << system_name << "\n";
        }

        // Extract the mesh.
        const MeshBase& mesh = d_es->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        // Extract the FE system and DOF map, and setup the FE object.
        System& system = d_es->get_system(system_name);
        const int sys_num = system.number();
        DofMap& dof_map = system.get_dof_map();
        SystemDofMapCache& dof_map_cache = *getDofMapCache(system_name);
        dof_map.compute_sparsity(mesh);
        FEType fe_type = dof_map.variable_type(0);
        std::unique_ptr<QBase> qrule = fe_type.default_quadrature_rule(dim);
        std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<double>& JxW = fe->get_JxW();
        const std::vector<std::vector<double> >& phi = fe->get_phi();

        // Build solver components.
        std::unique_ptr<NumericVector<double> > M_vec = system.solution->zero_clone();

        // Loop over the mesh to construct the system matrix.
        DenseMatrix<double> M_e;
        DenseVector<double> M_e_vec;
        std::vector<libMesh::dof_id_type> dof_id_scratch;
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            fe->reinit(elem);
            const auto& dof_indices = dof_map_cache.dof_indices(elem);
            for (unsigned int var_num = 0; var_num < dof_map.n_variables(); ++var_num)
            {
                const auto dof_indices_sz = static_cast<unsigned int>(dof_indices[var_num].size());
                M_e.resize(dof_indices_sz, dof_indices_sz);
                M_e_vec.resize(dof_indices_sz);
                const size_t n_basis = dof_indices[var_num].size();
                const unsigned int n_qp = qrule->n_points();
                for (unsigned int i = 0; i < n_basis; ++i)
                {
                    for (unsigned int j = 0; j < n_basis; ++j)
                    {
                        for (unsigned int qp = 0; qp < n_qp; ++qp)
                        {
                            M_e(i, j) += (phi[i][qp] * phi[j][qp]) * JxW[qp];
                        }
                    }
                }

                const double vol = elem->volume();
                double tr_M = 0.0;
                for (unsigned int i = 0; i < n_basis; ++i) tr_M += M_e(i, i);
                for (unsigned int i = 0; i < n_basis; ++i)
                {
                    M_e_vec(i) = vol * M_e(i, i) / tr_M;
                }

                dof_id_scratch = dof_indices[var_num];
                dof_map.constrain_element_vector(M_e_vec, dof_id_scratch);
                M_vec->add_vector(M_e_vec, dof_id_scratch);
            }
        }

        // Flush assemble the matrix.
        M_vec->close();

        // Reset values at Dirichlet boundaries.
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                if (elem->neighbor_ptr(side)) continue;
                static const boundary_id_type dirichlet_bdry_id_set[3] = { ZERO_DISPLACEMENT_X_BDRY_ID,
                                                                           ZERO_DISPLACEMENT_Y_BDRY_ID,
                                                                           ZERO_DISPLACEMENT_Z_BDRY_ID };
                std::vector<boundary_id_type> bdry_ids;
                mesh.boundary_info->boundary_ids(elem, side, bdry_ids);
                const boundary_id_type dirichlet_bdry_ids = get_dirichlet_bdry_ids(bdry_ids);
                if (!dirichlet_bdry_ids) continue;
                fe->reinit(elem);
                for (unsigned int n = 0; n < elem->n_nodes(); ++n)
                {
                    if (elem->is_node_on_side(n, side))
                    {
                        const Node* const node = elem->node_ptr(n);
                        for (unsigned int var_num = 0; var_num < dof_map.n_variables(); ++var_num)
                        {
                            const unsigned int n_comp = node->n_comp(sys_num, var_num);
                            for (unsigned int comp = 0; comp < n_comp; ++comp)
                            {
                                if (!(dirichlet_bdry_ids & dirichlet_bdry_id_set[comp])) continue;
                                const unsigned int node_dof_index = node->dof_number(sys_num, var_num, comp);
                                if (!dof_map.is_constrained_dof(node_dof_index)) continue;
                                M_vec->set(node_dof_index, 1.0);
                            }
                        }
                    }
                }
            }
        }

        // Assemble the vector.
        M_vec->close();

        // Store the diagonal mass matrix.
        d_L2_proj_matrix_diag[system_name] = std::move(M_vec);
    }

    IBTK_TIMER_STOP(t_build_diagonal_l2_mass_matrix);
    return d_L2_proj_matrix_diag[system_name].get();
} // buildDiagonalL2MassMatrix

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
    IBTK_TIMER_START(t_compute_l2_projection);

    int ierr;
    bool converged = false;

    if (close_F) F_vec.close();
    const System& system = d_es->get_system(system_name);
    const DofMap& dof_map = system.get_dof_map();
    if (consistent_mass_matrix)
    {
        std::pair<libMesh::LinearSolver<double>*, SparseMatrix<double>*> proj_solver_components =
            buildL2ProjectionSolver(system_name);
        auto solver = static_cast<PetscLinearSolver<double>*>(proj_solver_components.first);
        auto M_mat = static_cast<PetscMatrix<double>*>(proj_solver_components.second);
        PetscBool rtol_set;
        double runtime_rtol;
        ierr = PetscOptionsGetReal(nullptr, "", "-ksp_rtol", &runtime_rtol, &rtol_set);
        IBTK_CHKERRQ(ierr);
        PetscBool max_it_set;
        int runtime_max_it;
        ierr = PetscOptionsGetInt(nullptr, "", "-ksp_max_it", &runtime_max_it, &max_it_set);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetFromOptions(solver->ksp());
        IBTK_CHKERRQ(ierr);
        solver->solve(
            *M_mat, *M_mat, U_vec, F_vec, rtol_set ? runtime_rtol : tol, max_it_set ? runtime_max_it : max_its);
        KSPConvergedReason reason;
        ierr = KSPGetConvergedReason(solver->ksp(), &reason);
        IBTK_CHKERRQ(ierr);
        converged = reason > 0;
    }
    else
    {
        auto M_diag_vec = static_cast<PetscVector<double>*>(buildDiagonalL2MassMatrix(system_name));
        Vec M_diag_petsc_vec = M_diag_vec->vec();
        Vec U_petsc_vec = static_cast<PetscVector<double>*>(&U_vec)->vec();
        Vec F_petsc_vec = static_cast<PetscVector<double>*>(&F_vec)->vec();
        ierr = VecPointwiseDivide(U_petsc_vec, F_petsc_vec, M_diag_petsc_vec);
        IBTK_CHKERRQ(ierr);
        converged = true;
    }
    if (close_U) U_vec.close();
    dof_map.enforce_constraints_exactly(system, &U_vec);

    IBTK_TIMER_STOP(t_compute_l2_projection);
    return converged;
} // computeL2Projection

bool
FEDataManager::updateQuadratureRule(std::unique_ptr<QBase>& qrule,
                                    QuadratureType type,
                                    Order order,
                                    bool use_adaptive_quadrature,
                                    double point_density,
                                    const Elem* const elem,
                                    const boost::multi_array<double, 2>& X_node,
                                    const double dx_min)
{
    unsigned int elem_dim = elem->dim();
    const unsigned int elem_p_level = elem->p_level();
    bool qrule_updated = false;

    ElemType elem_type;
    std::tie(elem_type, type, order) =
        getQuadratureKey(type, order, use_adaptive_quadrature, point_density, elem, X_node, dx_min);

    if (!qrule || qrule->type() != type || qrule->get_dim() != elem_dim || qrule->get_order() != order ||
        qrule->get_elem_type() != elem_type || qrule->get_p_level() != elem_p_level)
    {
        qrule =
            (type == QGRID ? std::unique_ptr<QBase>(new QGrid(elem_dim, order)) : QBase::build(type, elem_dim, order));
        // qrule->allow_rules_with_negative_weights = false;
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
    return updateQuadratureRule(
        qrule, spec.quad_type, spec.quad_order, spec.use_adaptive_quadrature, spec.point_density, elem, X_node, dx_min);
}

bool
FEDataManager::updateSpreadQuadratureRule(std::unique_ptr<QBase>& qrule,
                                          const FEDataManager::SpreadSpec& spec,
                                          const Elem* const elem,
                                          const boost::multi_array<double, 2>& X_node,
                                          const double dx_min)
{
    return updateQuadratureRule(
        qrule, spec.quad_type, spec.quad_order, spec.use_adaptive_quadrature, spec.point_density, elem, X_node, dx_min);
}

void
FEDataManager::addWorkloadEstimate(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                   const int workload_data_idx,
                                   const int coarsest_ln_in,
                                   const int finest_ln_in)
{
    IBTK_TIMER_START(t_update_workload_estimates);

    const int coarsest_ln = (coarsest_ln_in == -1) ? d_coarsest_ln : coarsest_ln_in;
    const int finest_ln = (finest_ln_in == -1) ? d_finest_ln : finest_ln_in;
    TBOX_ASSERT(coarsest_ln >= d_coarsest_ln && coarsest_ln <= d_finest_ln);
    TBOX_ASSERT(finest_ln >= d_coarsest_ln && finest_ln <= d_finest_ln);

    // Workload estimates are computed only on the level to which the FE mesh
    // has been assigned.
    const int ln = d_level_number;
    if (coarsest_ln <= ln && ln <= finest_ln)
    {
        updateQuadPointCountData(ln, ln);
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, ln, ln);
        hier_cc_data_ops.axpy(
            workload_data_idx, d_default_workload_spec.q_point_weight, d_qp_count_idx, workload_data_idx);
    }

    IBTK_TIMER_STOP(t_update_workload_estimates);
    return;
} // addWorkloadEstimate

void
FEDataManager::initializeLevelData(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                   const int level_number,
                                   const double /*init_data_time*/,
                                   const bool /*can_be_refined*/,
                                   const bool /*initial_time*/,
                                   const Pointer<BasePatchLevel<NDIM> > old_level,
                                   const bool /*allocate_data*/)
{
    IBTK_TIMER_START(t_initialize_level_data);

    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    if (old_level)
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));

    IBTK_TIMER_STOP(t_initialize_level_data);
    return;
} // initializeLevelData

void
FEDataManager::resetHierarchyConfiguration(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                           const int coarsest_ln,
                                           const int finest_ln)
{
    IBTK_TIMER_START(t_reset_hierarchy_configuration);

    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((coarsest_ln >= 0) && (coarsest_ln <= finest_ln) && (finest_ln <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        TBOX_ASSERT(hierarchy->getPatchLevel(ln));
    }

    // Reset the patch hierarchy and levels.
    setPatchHierarchy(hierarchy);
    setPatchLevels(0, d_hierarchy->getFinestLevelNumber());

    IBTK_TIMER_STOP(t_reset_hierarchy_configuration);
    return;
} // resetHierarchyConfiguration

void
FEDataManager::applyGradientDetector(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                     const int level_number,
                                     const double /*error_data_time*/,
                                     const int tag_index,
                                     const bool initial_time,
                                     const bool /*uses_richardson_extrapolation_too*/)
{
    if (level_number >= d_level_number) return;

    IBTK_TIMER_START(t_apply_gradient_detector);

    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));

    if (initial_time)
    {
        // Determine the active elements associated with the prescribed patch
        // level.
        std::vector<std::vector<Elem*> > active_level_elem_map;
        const IntVector<NDIM> ghost_width = 1;
        collectActivePatchElements(active_level_elem_map, level_number, ghost_width);
        std::vector<unsigned int> X_ghost_dofs;
        std::vector<Elem*> active_level_elems;
        collect_unique_elems(active_level_elems, active_level_elem_map);
        collectGhostDOFIndices(X_ghost_dofs, active_level_elems, COORDINATES_SYSTEM_NAME);

        // Extract the mesh.
        const MeshBase& mesh = d_es->get_mesh();
        const Parallel::Communicator& comm = mesh.comm();
        const unsigned int dim = mesh.mesh_dimension();

        // Extract the FE system and DOF map, and setup the FE object.
        System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
        const DofMap& X_dof_map = X_system.get_dof_map();
        SystemDofMapCache& X_dof_map_cache = *getDofMapCache(COORDINATES_SYSTEM_NAME);
        FEType fe_type = X_dof_map.variable_type(0);
        for (unsigned d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(X_dof_map.variable_type(d) == fe_type);
        }
        using quad_key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;
        FECache X_fe_cache(dim, fe_type);

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

            std::set<quad_key_type> used_quadratures;
            for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems[e_idx];
                const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);

                const quad_key_type key = getQuadratureKey(d_default_interp_spec.quad_type,
                                                           d_default_interp_spec.quad_order,
                                                           d_default_interp_spec.use_adaptive_quadrature,
                                                           d_default_interp_spec.point_density,
                                                           elem,
                                                           X_node,
                                                           patch_dx_min);
                FEBase& X_fe = X_fe_cache[key];
                QBase& qrule = d_quadrature_cache[key];

                // See the note in interpWeighted to explain why we override
                // libMesh's reinit logic here
                if (used_quadratures.find(key) == used_quadratures.end())
                {
                    X_fe.get_phi();
                    X_fe.reinit(elem);
                    used_quadratures.insert(key);
                }
                const std::vector<std::vector<double> >& X_phi = X_fe.get_phi();

                for (unsigned int qp = 0; qp < qrule.n_points(); ++qp)
                {
                    interpolate(&X_qp[0], qp, X_node, X_phi);

                    const Index<NDIM> i = IndexUtilities::getCellIndex(X_qp, grid_geom, ratio);
                    tag_data->fill(1, Box<NDIM>(i - Index<NDIM>(1), i + Index<NDIM>(1)));
                }
            }
        }

        X_petsc_vec->restore_array();
    }
    else if (level_number + 1 == d_level_number && level_number < d_hierarchy->getFinestLevelNumber())
    {
        Pointer<PatchLevel<NDIM> > finer_level = d_hierarchy->getPatchLevel(level_number + 1);
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);

        // Update the node count data and coarsen it from the finer level.
        updateQuadPointCountData(level_number, level_number + 1);
        Pointer<CoarsenOperator<NDIM> > coarsen_op = new CartesianCellDoubleWeightedAverage<NDIM>();
        Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
        coarsen_alg->registerCoarsen(d_qp_count_idx, d_qp_count_idx, coarsen_op);
        coarsen_alg->createSchedule(level, finer_level)->coarsenData();

        // Tag cells for refinement whenever they contain element quadrature
        // points.
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

    TBOX_ASSERT(db);
    db->putInteger("FE_DATA_MANAGER_VERSION", FE_DATA_MANAGER_VERSION);
    db->putInteger("d_coarsest_ln", d_coarsest_ln);
    db->putInteger("d_finest_ln", d_finest_ln);

    IBTK_TIMER_STOP(t_put_to_database);
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

FEDataManager::FEDataManager(std::string object_name,
                             FEDataManager::InterpSpec default_interp_spec,
                             FEDataManager::SpreadSpec default_spread_spec,
                             FEDataManager::WorkloadSpec default_workload_spec,
                             IntVector<NDIM> ghost_width,
                             bool register_for_restart)
    : d_object_name(std::move(object_name)),
      d_registered_for_restart(register_for_restart),
      d_default_workload_spec(default_workload_spec),
      d_default_interp_spec(default_interp_spec),
      d_default_spread_spec(default_spread_spec),
      d_ghost_width(std::move(ghost_width)),
      d_quadrature_cache(NDIM)
{
    TBOX_ASSERT(!d_object_name.empty());

    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }

    // Create/look up the variable context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name + "::CONTEXT");

    // Register the node count variable with the VariableDatabase.
    d_qp_count_var = new CellVariable<NDIM, double>(d_object_name + "::qp_count");
    d_qp_count_idx = var_db->registerVariableAndContext(d_qp_count_var, d_context, 0);

    // Setup Timers.
    IBTK_DO_ONCE(
        t_reinit_element_mappings =
            TimerManager::getManager()->getTimer("IBTK::FEDataManager::reinitElementMappings()");
        t_build_ghosted_solution_vector =
            TimerManager::getManager()->getTimer("IBTK::FEDataManager::buildGhostedSolutionVector()");
        t_spread = TimerManager::getManager()->getTimer("IBTK::FEDataManager::spread()");
        t_prolong_data = TimerManager::getManager()->getTimer("IBTK::FEDataManager::prolongData()");
        t_interp_weighted = TimerManager::getManager()->getTimer("IBTK::FEDataManager::interpWeighted()");
        t_interp = TimerManager::getManager()->getTimer("IBTK::FEDataManager::interp()");
        t_restrict_data = TimerManager::getManager()->getTimer("IBTK::FEDataManager::restrictData()");
        t_build_l2_projection_solver =
            TimerManager::getManager()->getTimer("IBTK::FEDataManager::buildL2ProjectionSolver()");
        t_build_diagonal_l2_mass_matrix =
            TimerManager::getManager()->getTimer("IBTK::FEDataManager::buildDiagonalL2MassMatrix()");
        t_compute_l2_projection = TimerManager::getManager()->getTimer("IBTK::FEDataManager::computeL2Projection()");
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

void
FEDataManager::setLoggingEnabled(bool enable_logging)
{
    d_enable_logging = enable_logging;
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
    // Set the node count data on the specified range of levels of the
    // hierarchy.
    unsigned long n_local_q_points = 0;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const IntVector<NDIM>& ratio = level->getRatio();
        const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
        if (!level->checkAllocated(d_qp_count_idx)) level->allocatePatchData(d_qp_count_idx);
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(d_hierarchy, ln, ln);
        hier_cc_data_ops.setToScalar(d_qp_count_idx, 0.0);
        if (ln != d_level_number) continue;

        // Extract the mesh.
        const MeshBase& mesh = d_es->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        // Extract the FE system and DOF map, and setup the FE object.
        System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
        const DofMap& X_dof_map = X_system.get_dof_map();
        SystemDofMapCache& X_dof_map_cache = *getDofMapCache(COORDINATES_SYSTEM_NAME);
        FEType fe_type = X_dof_map.variable_type(0);
        for (unsigned d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(X_dof_map.variable_type(d) == fe_type);
        }

        // convenience alias for the quadrature key type used by FECache and FEMapCache
        using quad_key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;
        FECache X_fe_cache(dim, fe_type);

        // Extract the underlying solution data.
        NumericVector<double>* X_ghost_vec = buildGhostedCoordsVector();
        auto X_petsc_vec = static_cast<PetscVector<double>*>(X_ghost_vec);
        const double* const X_local_soln = X_petsc_vec->get_array_read();

        // Determine the number of element quadrature points associated with
        // each Cartesian grid cell.
        boost::multi_array<double, 2> X_node;
        Point X_qp;
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            const std::vector<Elem*>& patch_elems = d_active_patch_elem_map[local_patch_num];
            const size_t num_active_patch_elems = patch_elems.size();
            if (!num_active_patch_elems) continue;

            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);

            Pointer<CellData<NDIM, double> > qp_count_data = patch->getPatchData(d_qp_count_idx);

            std::set<quad_key_type> used_quadratures;
            for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems[e_idx];
                const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);

                const quad_key_type key = getQuadratureKey(d_default_interp_spec.quad_type,
                                                           d_default_interp_spec.quad_order,
                                                           d_default_interp_spec.use_adaptive_quadrature,
                                                           d_default_interp_spec.point_density,
                                                           elem,
                                                           X_node,
                                                           patch_dx_min);
                FEBase& X_fe = X_fe_cache[key];
                QBase& qrule = d_quadrature_cache[key];

                // See the note in interpWeighted to explain why we override
                // libMesh's reinit logic here
                if (used_quadratures.find(key) == used_quadratures.end())
                {
                    X_fe.get_phi();
                    X_fe.reinit(elem);
                    used_quadratures.insert(key);
                }
                const std::vector<std::vector<double> >& X_phi = X_fe.get_phi();

                for (unsigned int qp = 0; qp < qrule.n_points(); ++qp)
                {
                    interpolate(&X_qp[0], qp, X_node, X_phi);
                    const Index<NDIM> i = IndexUtilities::getCellIndex(X_qp, grid_geom, ratio);
                    if (patch_box.contains(i))
                    {
                        (*qp_count_data)(i) += 1.0;
                        ++n_local_q_points;
                    }
                }
            }
        }

        if (d_enable_logging)
        {
            const int n_processes = SAMRAI::tbox::SAMRAI_MPI::getNodes();
            const int current_rank = SAMRAI::tbox::SAMRAI_MPI::getRank();
            const auto right_padding = std::size_t(std::log10(n_processes)) + 1;

            std::vector<unsigned long> n_q_points_on_processors(n_processes);
            n_q_points_on_processors[current_rank] = n_local_q_points;

            const int ierr = MPI_Allreduce(MPI_IN_PLACE,
                                           n_q_points_on_processors.data(),
                                           n_q_points_on_processors.size(),
                                           MPI_UNSIGNED_LONG,
                                           MPI_SUM,
                                           SAMRAI::tbox::SAMRAI_MPI::commWorld);
            TBOX_ASSERT(ierr == 0);
            if (current_rank == 0)
            {
                for (int rank = 0; rank < n_processes; ++rank)
                {
                    SAMRAI::tbox::plog << "quadrature points on processor " << std::setw(right_padding) << std::left
                                       << rank << " = " << n_q_points_on_processors[rank] << '\n';
                }
            }
        }

        X_petsc_vec->restore_array();
    }
    return;
} // updateQuadPointCountData

std::vector<std::pair<Point, Point> >*
FEDataManager::computeActiveElementBoundingBoxes()
{
    // Get the necessary FE data.
    const MeshBase& mesh = d_es->get_mesh();
    const unsigned int n_elem = mesh.max_elem_id() + 1;
    System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_vec = *X_system.solution;
    NumericVector<double>& X_ghost_vec = *X_system.current_local_solution;
    copy_and_synch(X_vec, X_ghost_vec, /*close_v_in*/ false);

    // Compute the lower and upper bounds of all active local elements in the
    // mesh.  Assumes nodal basis functions.
    d_active_elem_bboxes.resize(n_elem);
    std::fill(d_active_elem_bboxes.begin(), d_active_elem_bboxes.end(), std::make_pair(Point::Zero(), Point::Zero()));
    std::vector<unsigned int> dof_indices;
    MeshBase::const_element_iterator el_it = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const unsigned int elem_id = elem->id();
        Point& elem_lower_bound = d_active_elem_bboxes[elem_id].first;
        Point& elem_upper_bound = d_active_elem_bboxes[elem_id].second;
        elem_lower_bound = Point::Constant(std::numeric_limits<double>::max());
        elem_upper_bound = Point::Constant(-std::numeric_limits<double>::max());

        const unsigned int n_nodes = elem->n_nodes();
        dof_indices.clear();
        for (unsigned int k = 0; k < n_nodes; ++k)
        {
            const Node* const node = elem->node_ptr(k);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                TBOX_ASSERT(node->n_dofs(X_sys_num, d) == 1);
                dof_indices.push_back(node->dof_number(X_sys_num, d, 0));
            }
        }
        std::vector<double> X_node;
        X_ghost_vec.get(dof_indices, X_node);
        for (unsigned int k = 0; k < n_nodes; ++k)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const double& X = X_node[k * NDIM + d];
                elem_lower_bound[d] = std::min(elem_lower_bound[d], X);
                elem_upper_bound[d] = std::max(elem_upper_bound[d], X);
            }
        }
    }

    // Parallel sum elem_lower_bound and elem_upper_bound so that each process
    // has access to the bounding box data for each active element in the mesh.
    std::vector<double> d_active_elem_bboxes_flattened(2 * NDIM * n_elem);
    for (unsigned int e = 0; e < n_elem; ++e)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            d_active_elem_bboxes_flattened[2 * e * NDIM + d] = d_active_elem_bboxes[e].first[d];
            d_active_elem_bboxes_flattened[(2 * e + 1) * NDIM + d] = d_active_elem_bboxes[e].second[d];
        }
    }
    SAMRAI_MPI::sumReduction(&d_active_elem_bboxes_flattened[0],
                             static_cast<int>(d_active_elem_bboxes_flattened.size()));
    for (unsigned int e = 0; e < n_elem; ++e)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            d_active_elem_bboxes[e].first[d] = d_active_elem_bboxes_flattened[2 * e * NDIM + d];
            d_active_elem_bboxes[e].second[d] = d_active_elem_bboxes_flattened[(2 * e + 1) * NDIM + d];
        }
    }
    return &d_active_elem_bboxes;
} // computeActiveElementBoundingBoxes

void
FEDataManager::collectActivePatchElements(std::vector<std::vector<Elem*> >& active_patch_elems,
                                          const int level_number,
                                          const IntVector<NDIM>& ghost_width)
{
    // Get the necessary FE data.
    const MeshBase& mesh = d_es->get_mesh();
    const Parallel::Communicator& comm = mesh.comm();
    const unsigned int dim = mesh.mesh_dimension();
    System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    SystemDofMapCache& X_dof_map_cache = *getDofMapCache(COORDINATES_SYSTEM_NAME);
    FEType fe_type = X_dof_map.variable_type(0);
    for (unsigned d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == fe_type);
    }
    NumericVector<double>* X_vec = getCoordsVector();
    std::unique_ptr<NumericVector<double> > X_ghost_vec = NumericVector<double>::build(comm);

    // convenience alias for the quadrature key type used by FECache and FEMapCache
    using quad_key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;
    FECache X_fe_cache(dim, fe_type);

    // Setup data structures used to assign elements to patches.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
    const IntVector<NDIM>& ratio = level->getRatio();
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    const int num_local_patches = level->getProcessorMapping().getNumberOfLocalIndices();
    std::vector<std::set<Elem*> > local_patch_elems(num_local_patches);
    std::vector<std::set<Elem*> > nonlocal_patch_elems(num_local_patches);
    std::vector<std::set<Elem*> > frontier_patch_elems(num_local_patches);

    // We provisionally associate an element with a Cartesian grid patch if the
    // element's bounding box intersects the patch interior grown by the
    // specified ghost cell width.
    //
    // NOTE: Following the call to computeActiveElementBoundingBoxes, each
    // processor will have access to all of the element bounding boxes.  This is
    // not a scalable approach, but we won't worry about this until it becomes
    // an actual issue.
    computeActiveElementBoundingBoxes();
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        std::set<Elem*>& frontier_elems = frontier_patch_elems[local_patch_num];
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        Point x_lower;
        Point x_upper;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            x_lower[d] = pgeom->getXLower()[d] - dx[d] * ghost_width[d];
            x_upper[d] = pgeom->getXUpper()[d] + dx[d] * ghost_width[d];
        }

        MeshBase::const_element_iterator el_it = mesh.active_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_elements_end();
        for (; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            const unsigned int elem_id = elem->id();
            const Point& elem_lower_bound = d_active_elem_bboxes[elem_id].first;
            const Point& elem_upper_bound = d_active_elem_bboxes[elem_id].second;
            bool in_patch = true;
            for (unsigned int d = 0; d < NDIM && in_patch; ++d)
            {
                in_patch = in_patch && ((elem_upper_bound[d] >= x_lower[d] && elem_upper_bound[d] <= x_upper[d]) ||
                                        (elem_lower_bound[d] >= x_lower[d] && elem_lower_bound[d] <= x_upper[d]));
            }
            if (in_patch)
            {
                frontier_elems.insert(elem);
            }
        }
    }

    // Recursively add/remove elements from the active sets that were generated
    // via the bounding box method.
    bool done = false;
    while (!done)
    {
        // Setup an appropriately ghosted temporary coordinates vector.
        std::vector<unsigned int> X_ghost_dofs;
        std::vector<Elem*> frontier_elems;
        collect_unique_elems(frontier_elems, frontier_patch_elems);
        collectGhostDOFIndices(X_ghost_dofs, frontier_elems, COORDINATES_SYSTEM_NAME);
        X_ghost_vec->init(X_vec->size(), X_vec->local_size(), X_ghost_dofs, true, GHOSTED);
        copy_and_synch(*X_vec, *X_ghost_vec, /*close_v_in*/ false);
        auto X_petsc_vec = static_cast<PetscVector<double>*>(X_ghost_vec.get());
        const double* const X_local_soln = X_petsc_vec->get_array_read();

        // Keep only those elements that have a quadrature point on the local
        // patch.
        boost::multi_array<double, 2> X_node;
        Point X_qp;
        int local_patch_num = 0;
        std::set<quad_key_type> used_quadratures;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            const std::set<Elem*>& frontier_elems = frontier_patch_elems[local_patch_num];
            std::set<Elem*>& local_elems = local_patch_elems[local_patch_num];
            std::set<Elem*>& nonlocal_elems = nonlocal_patch_elems[local_patch_num];
            if (frontier_elems.empty()) continue;

            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Box<NDIM> ghost_box = Box<NDIM>::grow(patch_box, ghost_width);
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);

            auto el_it = frontier_elems.begin();
            const std::set<Elem*>::const_iterator el_end = frontier_elems.end();
            for (; el_it != el_end; ++el_it)
            {
                Elem* const elem = *el_it;
                const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);

                const quad_key_type key = getQuadratureKey(d_default_interp_spec.quad_type,
                                                           d_default_interp_spec.quad_order,
                                                           d_default_interp_spec.use_adaptive_quadrature,
                                                           d_default_interp_spec.point_density,
                                                           elem,
                                                           X_node,
                                                           patch_dx_min);
                FEBase& X_fe = X_fe_cache[key];
                QBase& qrule = d_quadrature_cache[key];

                // See the note in interpWeighted to explain why we override
                // libMesh's reinit logic here
                if (used_quadratures.find(key) == used_quadratures.end())
                {
                    X_fe.get_phi();
                    X_fe.reinit(elem);
                    used_quadratures.insert(key);
                }
                const std::vector<std::vector<double> >& X_phi = X_fe.get_phi();

                bool found_qp = false;
                for (unsigned int qp = 0; qp < qrule.n_points() && !found_qp; ++qp)
                {
                    interpolate(&X_qp[0], qp, X_node, X_phi);
                    const Index<NDIM> i = IndexUtilities::getCellIndex(X_qp, grid_geom, ratio);
                    if (ghost_box.contains(i))
                    {
                        local_elems.insert(elem);
                        found_qp = true;
                    }
                }
                if (!found_qp) nonlocal_elems.insert(elem);
            }
        }

        X_petsc_vec->restore_array();

        // Rebuild the set of frontier elements, which are any neighbors of a
        // local element that has not already been determined to be either a
        // local or a nonlocal element.
        bool new_frontier = false;
        local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            std::set<Elem*>& frontier_elems = frontier_patch_elems[local_patch_num];
            const std::set<Elem*>& local_elems = local_patch_elems[local_patch_num];
            const std::set<Elem*>& nonlocal_elems = nonlocal_patch_elems[local_patch_num];
            frontier_elems.clear();
            if (local_elems.empty()) continue;

            for (const auto& elem : local_elems)
            {
                for (unsigned int n = 0; n < elem->n_neighbors(); ++n)
                {
                    Elem* const nghbr_elem = elem->neighbor_ptr(n);
                    if (nghbr_elem)
                    {
                        const bool is_local_elem = local_elems.find(nghbr_elem) != local_elems.end();
                        const bool is_nonlocal_elem = nonlocal_elems.find(nghbr_elem) != nonlocal_elems.end();
                        if (!(is_local_elem || is_nonlocal_elem))
                        {
                            frontier_elems.insert(nghbr_elem);
                            new_frontier = true;
                        }
                    }
                }
            }
        }

        // Check to see if we are done.
        done = SAMRAI_MPI::sumReduction(new_frontier ? 1 : 0) == 0;
    }

    // Set the active patch element data.
    active_patch_elems.resize(num_local_patches);
    local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        std::vector<Elem*>& active_elems = active_patch_elems[local_patch_num];
        const std::set<Elem*>& local_elems = local_patch_elems[local_patch_num];
        active_elems.resize(local_elems.size());
        int k = 0;
        for (auto cit = local_elems.begin(); cit != local_elems.end(); ++cit, ++k)
        {
            active_elems[k] = *cit;
        }
    }
    return;
} // collectActivePatchElements

void
FEDataManager::collectActivePatchNodes(std::vector<std::vector<Node*> >& active_patch_nodes,
                                       const std::vector<std::vector<Elem*> >& active_patch_elems)
{
    const MeshBase& mesh = d_es->get_mesh();
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

void
FEDataManager::collectGhostDOFIndices(std::vector<unsigned int>& ghost_dofs,
                                      const std::vector<Elem*>& active_elems,
                                      const std::string& system_name)
{
    System& system = d_es->get_system(system_name);
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
    for (const auto& elem : active_elems)
    {
        // DOFs associated with the element.
        for (unsigned int var_num = 0; var_num < elem->n_vars(sys_num); ++var_num)
        {
            if (elem->n_dofs(sys_num, var_num) > 0)
            {
                const unsigned int dof_index = elem->dof_number(sys_num, var_num, 0);
                if (dof_index < first_local_dof || dof_index >= end_local_dof)
                {
                    ghost_dof_set.insert(dof_index);
                }
            }
        }

        // DOFs associated with the nodes of the element.
        for (unsigned int k = 0; k < elem->n_nodes(); ++k)
        {
            const Node* const node = elem->node_ptr(k);
            for (unsigned int var_num = 0; var_num < node->n_vars(sys_num); ++var_num)
            {
                if (node->n_dofs(sys_num, var_num) > 0)
                {
                    const unsigned int dof_index = node->dof_number(sys_num, var_num, 0);
                    if (dof_index < first_local_dof || dof_index >= end_local_dof)
                    {
                        ghost_dof_set.insert(dof_index);
                    }
                }
            }
        }
    }
    ghost_dofs.clear();
    ghost_dofs.insert(ghost_dofs.end(), ghost_dof_set.begin(), ghost_dof_set.end());
    return;
} // collectGhostDOFIndices

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

    d_coarsest_ln = db->getInteger("d_coarsest_ln");
    d_finest_ln = db->getInteger("d_finest_ln");
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
