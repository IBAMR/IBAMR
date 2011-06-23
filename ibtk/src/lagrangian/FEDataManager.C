// Filename: FEDataManager.C
// Created on 19 Apr 2010 by Boyce Griffith
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

#include "FEDataManager.h"

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
#include <ibtk/IndexUtilities.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/compiler_hints.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/namespaces.h>

// LIBMESH INCLUDES
#include <boundary_info.h>
#include <fe.h>
#include <fe_interface.h>
#include <dense_matrix.h>
#include <dense_vector.h>
#include <dof_map.h>
#include <explicit_system.h>
#include <mesh.h>
#include <numeric_vector.h>
#include <parallel.h>
#include <petsc_linear_solver.h>
#include <petsc_matrix.h>
#include <petsc_vector.h>
#include <quadrature_gauss.h>
using namespace libMesh;

// SAMRAI INCLUDES
#include <CartesianCellDoubleWeightedAverage.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CoarsenAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>

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
static Timer* t_prolong_value;
static Timer* t_prolong_density;
static Timer* t_interp;
static Timer* t_restrict_value;
static Timer* t_build_l2_projection_solver;
static Timer* t_build_diagonal_l2_mass_matrix;
static Timer* t_compute_l2_projection;
static Timer* t_update_workload_data;
static Timer* t_initialize_level_data;
static Timer* t_reset_hierarchy_configuration;
static Timer* t_apply_gradient_detector;
static Timer* t_put_to_database;

// Version of FEDataManager restart file data.
static const int FE_DATA_MANAGER_VERSION = 1;

// Local helper functions.
template<class T>
inline void
flatten(
    blitz::Array<Elem*,1>& elems,
    const T& elem_patch_map)
{
    std::set<Elem*> elem_set;
    for (int k = 0; k < elem_patch_map.size(); ++k)
    {
        elem_set.insert(elem_patch_map(k).begin(),elem_patch_map(k).end());
    }

    elems.resize(elem_set.size());
    int k = 0;
    for (std::set<Elem*>::const_iterator cit = elem_set.begin(); cit != elem_set.end(); ++cit, ++k)
    {
        elems(k) = *cit;
    }
    return;
}// flatten
}

const short int FEDataManager::DIRICHLET_BDRY_ID;
std::map<std::string,FEDataManager*> FEDataManager::s_data_manager_instances;
bool FEDataManager::s_registered_callback;
unsigned char FEDataManager::s_shutdown_priority;

FEDataManager*
FEDataManager::getManager(
    const std::string& name,
    const std::string& interp_weighting_fcn,
    const std::string& spread_weighting_fcn,
    QBase* const qrule,
    const bool interp_uses_consistent_mass_matrix,
    bool register_for_restart)
{
    if (s_data_manager_instances.find(name) == s_data_manager_instances.end())
    {
        const int stencil_size = std::max(LEInteractor::getStencilSize(interp_weighting_fcn),LEInteractor::getStencilSize(spread_weighting_fcn));
        const IntVector<NDIM> ghost_width = static_cast<int>(floor(0.5*static_cast<double>(stencil_size)))+1;
        s_data_manager_instances[name] = new FEDataManager(name, interp_weighting_fcn, spread_weighting_fcn, qrule, interp_uses_consistent_mass_matrix, ghost_width, register_for_restart);
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeAllManagers, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_data_manager_instances[name];
}// getManager

void
FEDataManager::freeAllManagers()
{
    for (std::map<std::string,FEDataManager*>::iterator it = s_data_manager_instances.begin();
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
FEDataManager::registerLoadBalancer(
    Pointer<LoadBalancer<NDIM> > load_balancer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!load_balancer.isNull());
#endif
    d_load_balancer = load_balancer;
    return;
}// return

void
FEDataManager::setPatchHierarchy(
    Pointer<PatchHierarchy<NDIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
#endif
    // Reset the hierarchy.
    d_hierarchy = hierarchy;
    return;
}// setPatchHierarchy

void
FEDataManager::resetLevels(
    const int coarsest_ln,
    const int finest_ln)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((coarsest_ln >= 0) &&
                (finest_ln >= coarsest_ln) &&
                (finest_ln <= d_hierarchy->getFinestLevelNumber()));
#endif
    // Reset the level numbers.
    d_coarsest_ln = coarsest_ln;
    d_finest_ln   = finest_ln;
    return;
}// resetLevels

void
FEDataManager::setEquationSystems(
    EquationSystems* const equation_systems,
    const int level_number)
{
    d_es = equation_systems;
    d_level_number = level_number;
    return;
}// setEquationSystems

EquationSystems*
FEDataManager::getEquationSystems() const
{
    return d_es;
}// getEquationSystems

int
FEDataManager::getLevelNumber() const
{
    return d_level_number;
}// getLevelNumber

const IntVector<NDIM>&
FEDataManager::getGhostCellWidth() const
{
    return d_ghost_width;
}// getGhostCellWidth

const std::string&
FEDataManager::getInterpWeightingFunction() const
{
    return d_interp_weighting_fcn;
}// getInterpWeightingFunction

const std::string&
FEDataManager::getSpreadWeightingFunction() const
{
    return d_spread_weighting_fcn;
}// getSpreadWeightingFunction

QBase*
FEDataManager::getQuadratureRule() const
{
    return d_qrule;
}// getQuadratureRule

bool
FEDataManager::getInterpUsesConsistentMassMatrix() const
{
    return d_interp_uses_consistent_mass_matrix;
}// getInterpUsesConsistentMassMatrix

const blitz::Array<blitz::Array<Elem*,1>,1>&
FEDataManager::getActivePatchElementMap() const
{
    return d_active_patch_elem_map;
}// getActivePatchElementMap

void
FEDataManager::reinitElementMappings()
{
    IBTK_TIMER_START(t_reinit_element_mappings);

    // Delete cached hierarchy-dependent data.
    d_active_patch_elem_map  .free();
    d_active_patch_ghost_dofs.clear();
    for (std::map<std::string,NumericVector<double>*>::iterator it = d_system_ghost_vec.begin();
         it != d_system_ghost_vec.end(); ++it)
    {
        delete it->second;
    }
    d_system_ghost_vec.clear();

    // Reset the mappings between grid patches and active mesh elements.
    collectActivePatchElements(d_active_patch_elem_map, d_level_number, d_ghost_width);

    IBTK_TIMER_STOP(t_reinit_element_mappings);
    return;
}// reinitElementMappings

NumericVector<double>*
FEDataManager::getSolutionVector(
    const std::string& system_name) const
{
    System& system = d_es->get_system(system_name);
    return system.solution.get();
}// getSolutionVector

NumericVector<double>*
FEDataManager::buildGhostedSolutionVector(
    const std::string& system_name)
{
    IBTK_TIMER_START(t_build_ghosted_solution_vector);

    NumericVector<double>* sol_vec = getSolutionVector(system_name);
    if (d_system_ghost_vec.count(system_name) == 0)
    {
        if (d_active_patch_ghost_dofs.count(system_name) == 0)
        {
            blitz::Array<Elem*,1> active_elems;
            flatten(active_elems, d_active_patch_elem_map);
            collectGhostDOFIndices(d_active_patch_ghost_dofs[system_name], active_elems, system_name);
        }
        AutoPtr<NumericVector<double> > sol_ghost_vec = NumericVector<double>::build();
        sol_ghost_vec->init(sol_vec->size(), sol_vec->local_size(), d_active_patch_ghost_dofs[system_name], true, GHOSTED);
        d_system_ghost_vec[system_name] = sol_ghost_vec.release();
    }
    System& system = d_es->get_system(system_name);
    NumericVector<double>* sol_ghost_vec = d_system_ghost_vec[system_name];
    sol_vec->localize(*sol_ghost_vec);
    system.get_dof_map().enforce_constraints_exactly(system, sol_ghost_vec);

    IBTK_TIMER_STOP(t_build_ghosted_solution_vector);
    return sol_ghost_vec;
}// buildGhostedSolutionVector

NumericVector<double>*
FEDataManager::getCoordsVector() const
{
    return getSolutionVector(COORDINATES_SYSTEM_NAME);
}// getCoordsVector

NumericVector<double>*
FEDataManager::buildGhostedCoordsVector()
{
    return buildGhostedSolutionVector(COORDINATES_SYSTEM_NAME);
}// buildGhostedCoordsVector

void
FEDataManager::spread(
    const int f_data_idx,
    NumericVector<double>& F_vec,
    NumericVector<double>& X_vec,
    const std::string& system_name,
    const bool close_F,
    const bool close_X)
{
    IBTK_TIMER_START(t_spread);

    // Extract the mesh.
    const MeshBase& mesh = d_es->get_mesh();
    const int dim = mesh.mesh_dimension();

    // Extract the FE systems and DOF maps, and setup the FE objects.
    System& F_system = d_es->get_system(system_name);
    const unsigned int n_vars = F_system.n_vars();
    const DofMap& F_dof_map = F_system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned i = 0; i < n_vars; ++i) TBOX_ASSERT(F_dof_map.variable_type(i) == F_dof_map.variable_type(0));
#endif
    blitz::Array<std::vector<unsigned int>,1> F_dof_indices(n_vars);
    for (unsigned int i = 0; i < n_vars; ++i) F_dof_indices(i).reserve(NDIM == 2 ? 9 : 27);
    AutoPtr<FEBase> F_fe(FEBase::build(dim, F_dof_map.variable_type(0)));
    F_fe->attach_quadrature_rule(d_qrule);
    const std::vector<double>& JxW_F = F_fe->get_JxW();
    const std::vector<std::vector<double> >& phi_F = F_fe->get_phi();

    System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_dof_map.variable_type(0));
#endif
    blitz::Array<std::vector<unsigned int>,1> X_dof_indices(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) X_dof_indices(d).reserve(NDIM == 2 ? 9 : 27);
    AutoPtr<FEBase> X_fe(FEBase::build(dim, X_dof_map.variable_type(0)));
    X_fe->attach_quadrature_rule(d_qrule);
    const std::vector<std::vector<double> >& phi_X = X_fe->get_phi();

    // Communicate any unsynchronized ghost data and enforce any constraints.
    if (close_F) F_vec.close();
    F_dof_map.enforce_constraints_exactly(F_system, &F_vec);

    if (close_X) X_vec.close();
    X_dof_map.enforce_constraints_exactly(X_system, &X_vec);

    // Loop over the patches to interpolate nodal values on the FE mesh to the
    // element quadrature points, then spread thost values onto the Eulerian
    // grid.
    blitz::Array<double,2> F_node, X_node;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_number);
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const blitz::Array<Elem*,1>& patch_elems = d_active_patch_elem_map(local_patch_num);
        const unsigned int num_active_patch_elems = patch_elems.size();
        if (num_active_patch_elems == 0) continue;

        // Setup vectors to store the values of F_JxW and X at the quadrature
        // points.  We compute a conservative upper bound on the number of
        // quadrature points to try to avoid unnecessary reallocations.
        if (dim == 2)
        {
            d_qrule->init(QUAD9);
        }
        if (dim == 3)
        {
            d_qrule->init(HEX27);
        }
        const unsigned int n_qp_estimate = d_qrule->n_points();
        static const unsigned int safety_factor = 2;
        std::vector<double> F_JxW_qp(safety_factor*n_vars*n_qp_estimate*num_active_patch_elems);
        std::vector<double>     X_qp(safety_factor*NDIM  *n_qp_estimate*num_active_patch_elems);

        // Loop over the elements and compute the values to be spread and the
        // positions of the quadrature points.
        int qp_offset = 0;
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            const Elem* const elem = patch_elems(e_idx);

            F_fe->reinit(elem);
            for (unsigned int i = 0; i < n_vars; ++i)
            {
                F_dof_map.dof_indices(elem, F_dof_indices(i), i);
            }

            X_fe->reinit(elem);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(elem, X_dof_indices(d), d);
            }

            const unsigned int n_qp = d_qrule->n_points();
            if (UNLIKELY(F_JxW_qp.size() < n_vars*(qp_offset+n_qp))) F_JxW_qp.resize(n_vars*(qp_offset+n_qp));
            if (UNLIKELY(    X_qp.size() < NDIM  *(qp_offset+n_qp)))     X_qp.resize(NDIM  *(qp_offset+n_qp));

            get_values_for_interpolation(F_node, F_vec, F_dof_indices);
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const int idx = n_vars*(qp+qp_offset);
                interpolate(&F_JxW_qp[idx],qp,F_node,phi_F);
                for (unsigned int i = 0; i < n_vars; ++i)
                {
                    F_JxW_qp[idx+i] *= JxW_F[qp];
                }
            }

            get_values_for_interpolation(X_node, X_vec, X_dof_indices);
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const int idx = NDIM*(qp+qp_offset);
                interpolate(&X_qp[idx],qp,X_node,phi_X);
            }

            qp_offset += n_qp;
        }

        if (qp_offset == 0) continue;

        F_JxW_qp.resize(n_vars*qp_offset);
        X_qp    .resize(NDIM  *qp_offset);

        // Spread values from the quadrature points to the Cartesian grid patch.
        //
        // NOTE: Values are spread only from those quadrature points that are
        // within the ghost cell width of the patch interior.
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM> spread_box = Box<NDIM>::grow(patch->getBox(), d_ghost_width);
        Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_data_idx);
        Pointer<CellData<NDIM,double> > f_cc_data = f_data;
        Pointer<SideData<NDIM,double> > f_sc_data = f_data;
        const bool is_cc_data = !f_cc_data.isNull();
        const bool is_sc_data = !f_sc_data.isNull();
        if (is_cc_data) LEInteractor::spread(f_cc_data, F_JxW_qp, n_vars, X_qp, NDIM, patch, spread_box, d_spread_weighting_fcn);
        if (is_sc_data) LEInteractor::spread(f_sc_data, F_JxW_qp, n_vars, X_qp, NDIM, patch, spread_box, d_spread_weighting_fcn);
    }

    IBTK_TIMER_STOP(t_spread);
    return;
}// spread

void
FEDataManager::prolongValue(
    const int f_data_idx,
    NumericVector<double>& F_vec,
    NumericVector<double>& X_vec,
    const std::string& system_name,
    const bool close_F,
    const bool close_X)
{
    IBTK_TIMER_START(t_prolong_value);

    // NOTE #1: This routine is sepcialized for a staggered-grid Eulerian
    // discretization.  It should be straightforward to generalize it to work
    // with other data centerings.
    //
    // NOTE #2: This code is specialized for isoparametric elements.  It is less
    // clear how to relax this assumption.

    // Extract the mesh.
    const MeshBase& mesh = d_es->get_mesh();
    const int dim = mesh.mesh_dimension();

    // Extract the FE systems and DOF maps, and setup the FE objects.
    System& F_system = d_es->get_system(system_name);
    const unsigned int n_vars = F_system.n_vars();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(n_vars == NDIM);  // specialized to side-centered data
#endif
    const DofMap& F_dof_map = F_system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned i = 0; i < n_vars; ++i) TBOX_ASSERT(F_dof_map.variable_type(i) == F_dof_map.variable_type(0));
#endif
    blitz::Array<std::vector<unsigned int>,1> F_dof_indices(n_vars);
    for (unsigned int i = 0; i < n_vars; ++i) F_dof_indices(i).reserve(NDIM == 2 ? 9 : 27);
    AutoPtr<FEBase> F_fe(FEBase::build(dim, F_dof_map.variable_type(0)));
    const std::vector<std::vector<double> >& phi_F = F_fe->get_phi();

    System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_dof_map.variable_type(0));
#endif
    FEType X_fe_type = X_dof_map.variable_type(0);
    blitz::Array<std::vector<unsigned int>,1> X_dof_indices(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) X_dof_indices(d).reserve(NDIM == 2 ? 9 : 27);

    // Communicate any unsynchronized ghost data and enforce any constraints.
    if (close_F) F_vec.close();
    F_dof_map.enforce_constraints_exactly(F_system, &F_vec);

    if (close_X) X_vec.close();
    X_dof_map.enforce_constraints_exactly(X_system, &X_vec);

    // Loop over the patches to interpolate nodal values on the FE mesh to the
    // the points of the Eulerian grid.
    blitz::Array<double,2> F_node;
    static const unsigned int MAX_NODES = (NDIM == 2 ? 9 : 27);
    Point s_node_cache[MAX_NODES], X_node_cache[MAX_NODES];
    blitz::TinyVector<double,NDIM> X_min, X_max;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_number);
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const blitz::Array<Elem*,1>& patch_elems = d_active_patch_elem_map(local_patch_num);
        const unsigned int num_active_patch_elems = patch_elems.size();
        if (num_active_patch_elems == 0) continue;

        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<SideData<NDIM,double> > f_data = patch->getPatchData(f_data_idx);
        const Box<NDIM>& patch_box = patch->getBox();
        const CellIndex<NDIM>& patch_lower = patch_box.lower();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_x_lower = patch_geom->getXLower();
        const double* const patch_dx = patch_geom->getDx();

        blitz::TinyVector<Box<NDIM>,NDIM> side_boxes;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            side_boxes[axis] = SideGeometry<NDIM>::toSideBox(patch_box,axis);
        }

        SideData<NDIM,int> spread_value_at_loc(patch_box, 1, IntVector<NDIM>(0));
        spread_value_at_loc.fillAll(0);

        // Loop over the elements and compute the values to be prolonged.
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems(e_idx);
            const unsigned int n_node = elem->n_nodes();
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(elem, X_dof_indices(d), d);
            }

            // Cache the nodal and physical coordinates of the element,
            // determine the bounding box of the current configuration of the
            // element, and set the nodal coordinates of the element to
            // correspond to the physical coordinates.
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(n_node <= MAX_NODES);
#endif
            X_min =  0.5*std::numeric_limits<double>::max();
            X_max = -0.5*std::numeric_limits<double>::max();
            for (unsigned int k = 0; k < n_node; ++k)
            {
                s_node_cache[k] = elem->point(k);
                for (int d = 0; d < NDIM; ++d)
                {
                    X_node_cache[k](d) = X_vec(X_dof_indices(d)[k]);
                    X_min[d] = std::min(X_min[d],X_node_cache[k](d));
                    X_max[d] = std::max(X_max[d],X_node_cache[k](d));
                }
                elem->point(k) = X_node_cache[k];
            }

            // Loop over coordinate directions and look for Eulerian grid points
            // that are covered by the element.
            std::vector<Point>            intersection_master_coords;
            std::vector<SideIndex<NDIM> > intersection_indices;
            static const int estimated_max_size = (NDIM == 2 ? 64 : 512);
            intersection_master_coords.reserve(estimated_max_size);
            intersection_indices      .reserve(estimated_max_size);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                // Loop over the relevant range of indices.
                blitz::TinyVector<int,NDIM> i_begin, i_end, ic;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    if (d == axis)
                    {
                        i_begin[d] = std::ceil((X_min[d]-patch_x_lower[d])/patch_dx[d]) + patch_lower[d];
                        i_end  [d] = std::ceil((X_max[d]-patch_x_lower[d])/patch_dx[d]) + patch_lower[d];
                    }
                    else
                    {
                        i_begin[d] = std::ceil((X_min[d]-patch_x_lower[d])/patch_dx[d] - 0.5) + patch_lower[d];
                        i_end  [d] = std::ceil((X_max[d]-patch_x_lower[d])/patch_dx[d] - 0.5) + patch_lower[d];
                    }
                }
#if (NDIM == 3)
                for (ic[2] = i_begin[2]; ic[2] < i_end[2]; ++ic[2])
                {
#endif
                    for (ic[1] = i_begin[1]; ic[1] < i_end[1]; ++ic[1])
                    {
                        for (ic[0] = i_begin[0]; ic[0] < i_end[0]; ++ic[0])
                        {
                            Point p;
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                p(d) = patch_x_lower[d] + patch_dx[d]*(static_cast<double>(ic[d]-patch_lower[d])+(d == axis ? 0.0 : 0.5));
                            }
                            const Point master_coords = FEInterface::inverse_map(dim, X_fe_type, elem, p, TOLERANCE, false);
                            if (FEInterface::on_reference_element(master_coords,elem->type()))
                            {
                                intersection_master_coords.push_back(master_coords);
#if (NDIM == 2)
                                SideIndex<NDIM> s_i(Index<NDIM>(ic[0],ic[1]),axis,0);
#endif
#if (NDIM == 3)
                                SideIndex<NDIM> s_i(Index<NDIM>(ic[0],ic[1],ic[2]),axis,0);
#endif
                                intersection_indices.push_back(s_i);
                            }
                        }
                    }
#if (NDIM == 3)
                }
#endif
            }

            // Restore the nodal coordinates.
            for (unsigned int k = 0; k < n_node; ++k)
            {
                elem->point(k) = s_node_cache[k];
            }

            // If there are no intersection points, then continue on to the next
            // element.
            if (intersection_master_coords.empty()) continue;

            // Evaluate the Lagrangian value at the Eulerian grid point, and
            // set the value on the Eulerian grid to be F/det(dX/ds).
            F_fe->reinit(elem, &intersection_master_coords);
            for (unsigned int i = 0; i < n_vars; ++i)
            {
                F_dof_map.dof_indices(elem, F_dof_indices(i), i);
            }

            get_values_for_interpolation(F_node, F_vec, F_dof_indices);
            for (unsigned int qp = 0; qp < intersection_master_coords.size(); ++qp)
            {
                const SideIndex<NDIM>& s_i = intersection_indices[qp];
                const int axis = s_i.getAxis();
                if (!side_boxes[axis].contains(s_i)) continue;
                if (spread_value_at_loc(s_i) != 0) continue;  // each value may be spread to only once
                spread_value_at_loc(s_i) = 1;
                const double F_qp = interpolate(qp,F_node(blitz::Range::all(),axis),phi_F);
                (*f_data)(s_i) += F_qp;
            }
        }
    }

    IBTK_TIMER_STOP(t_prolong_value);
    return;
}// prolongValue

void
FEDataManager::prolongDensity(
    const int f_data_idx,
    NumericVector<double>& F_vec,
    NumericVector<double>& X_vec,
    const std::string& system_name,
    const bool close_F,
    const bool close_X)
{
    IBTK_TIMER_START(t_prolong_density);

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
    // onto a (possibly discontinuous) FE basis instead of evaluating J directly
    // from the discrete deformation.

    // Extract the mesh.
    const MeshBase& mesh = d_es->get_mesh();
    const int dim = mesh.mesh_dimension();

    // Extract the FE systems and DOF maps, and setup the FE objects.
    System& F_system = d_es->get_system(system_name);
    const unsigned int n_vars = F_system.n_vars();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(n_vars == NDIM);  // specialized to side-centered data
#endif
    const DofMap& F_dof_map = F_system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned i = 0; i < n_vars; ++i) TBOX_ASSERT(F_dof_map.variable_type(i) == F_dof_map.variable_type(0));
#endif
    blitz::Array<std::vector<unsigned int>,1> F_dof_indices(n_vars);
    for (unsigned int i = 0; i < n_vars; ++i) F_dof_indices(i).reserve(NDIM == 2 ? 9 : 27);
    AutoPtr<FEBase> F_fe(FEBase::build(dim, F_dof_map.variable_type(0)));
    const std::vector<std::vector<double> >& phi_F = F_fe->get_phi();

    System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_dof_map.variable_type(0));
#endif
    FEType X_fe_type = X_dof_map.variable_type(0);
    blitz::Array<std::vector<unsigned int>,1> X_dof_indices(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) X_dof_indices(d).reserve(NDIM == 2 ? 9 : 27);
    AutoPtr<FEBase> X_fe(FEBase::build(dim, X_dof_map.variable_type(0)));
    const std::vector<std::vector<VectorValue<double> > >& dphi_X = X_fe->get_dphi();

    // Communicate any unsynchronized ghost data and enforce any constraints.
    if (close_F) F_vec.close();
    F_dof_map.enforce_constraints_exactly(F_system, &F_vec);

    if (close_X) X_vec.close();
    X_dof_map.enforce_constraints_exactly(X_system, &X_vec);

    // Loop over the patches to interpolate nodal values on the FE mesh to the
    // the points of the Eulerian grid.
    TensorValue<double> dX_ds;
    blitz::Array<double,2> F_node, X_node;
    static const unsigned int MAX_NODES = (NDIM == 2 ? 9 : 27);
    Point s_node_cache[MAX_NODES], X_node_cache[MAX_NODES];
    blitz::TinyVector<double,NDIM> X_min, X_max;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_number);
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const blitz::Array<Elem*,1>& patch_elems = d_active_patch_elem_map(local_patch_num);
        const unsigned int num_active_patch_elems = patch_elems.size();
        if (num_active_patch_elems == 0) continue;

        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<SideData<NDIM,double> > f_data = patch->getPatchData(f_data_idx);
        const Box<NDIM>& patch_box = patch->getBox();
        const CellIndex<NDIM>& patch_lower = patch_box.lower();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_x_lower = patch_geom->getXLower();
        const double* const patch_dx = patch_geom->getDx();

        blitz::TinyVector<Box<NDIM>,NDIM> side_boxes;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            side_boxes[axis] = SideGeometry<NDIM>::toSideBox(patch_box,axis);
        }

        SideData<NDIM,int> spread_value_at_loc(patch_box, 1, IntVector<NDIM>(0));
        spread_value_at_loc.fillAll(0);

        // Loop over the elements and compute the values to be prolonged.
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems(e_idx);
            const unsigned int n_node = elem->n_nodes();
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(elem, X_dof_indices(d), d);
            }

            // Cache the nodal and physical coordinates of the element,
            // determine the bounding box of the current configuration of the
            // element, and set the nodal coordinates of the element to
            // correspond to the physical coordinates.
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(n_node <= MAX_NODES);
#endif
            X_min =  0.5*std::numeric_limits<double>::max();
            X_max = -0.5*std::numeric_limits<double>::max();
            for (unsigned int k = 0; k < n_node; ++k)
            {
                s_node_cache[k] = elem->point(k);
                for (int d = 0; d < NDIM; ++d)
                {
                    X_node_cache[k](d) = X_vec(X_dof_indices(d)[k]);
                    X_min[d] = std::min(X_min[d],X_node_cache[k](d));
                    X_max[d] = std::max(X_max[d],X_node_cache[k](d));
                }
                elem->point(k) = X_node_cache[k];
            }

            // Loop over coordinate directions and look for Eulerian grid points
            // that are covered by the element.
            std::vector<Point>            intersection_master_coords;
            std::vector<SideIndex<NDIM> > intersection_indices;
            static const int estimated_max_size = (NDIM == 2 ? 64 : 512);
            intersection_master_coords.reserve(estimated_max_size);
            intersection_indices      .reserve(estimated_max_size);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                // Loop over the relevant range of indices.
                blitz::TinyVector<int,NDIM> i_begin, i_end, ic;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    if (d == axis)
                    {
                        i_begin[d] = std::ceil((X_min[d]-patch_x_lower[d])/patch_dx[d]) + patch_lower[d];
                        i_end  [d] = std::ceil((X_max[d]-patch_x_lower[d])/patch_dx[d]) + patch_lower[d];
                    }
                    else
                    {
                        i_begin[d] = std::ceil((X_min[d]-patch_x_lower[d])/patch_dx[d] - 0.5) + patch_lower[d];
                        i_end  [d] = std::ceil((X_max[d]-patch_x_lower[d])/patch_dx[d] - 0.5) + patch_lower[d];
                    }
                }
#if (NDIM == 3)
                for (ic[2] = i_begin[2]; ic[2] < i_end[2]; ++ic[2])
                {
#endif
                    for (ic[1] = i_begin[1]; ic[1] < i_end[1]; ++ic[1])
                    {
                        for (ic[0] = i_begin[0]; ic[0] < i_end[0]; ++ic[0])
                        {
                            Point p;
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                p(d) = patch_x_lower[d] + patch_dx[d]*(static_cast<double>(ic[d]-patch_lower[d])+(d == axis ? 0.0 : 0.5));
                            }
                            const Point master_coords = FEInterface::inverse_map(dim, X_fe_type, elem, p, TOLERANCE, false);
                            if (FEInterface::on_reference_element(master_coords,elem->type()))
                            {
                                intersection_master_coords.push_back(master_coords);
#if (NDIM == 2)
                                SideIndex<NDIM> s_i(Index<NDIM>(ic[0],ic[1]),axis,0);
#endif
#if (NDIM == 3)
                                SideIndex<NDIM> s_i(Index<NDIM>(ic[0],ic[1],ic[2]),axis,0);
#endif
                                intersection_indices.push_back(s_i);
                            }
                        }
                    }
#if (NDIM == 3)
                }
#endif
            }

            // Restore the nodal coordinates.
            for (unsigned int k = 0; k < n_node; ++k)
            {
                elem->point(k) = s_node_cache[k];
            }

            // If there are no intersection points, then continue on to the next
            // element.
            if (intersection_master_coords.empty()) continue;

            // Evaluate the Lagrangian density at the Eulerian grid point, and
            // set the value on the Eulerian grid to be F/det(dX/ds).
            F_fe->reinit(elem, &intersection_master_coords);
            X_fe->reinit(elem, &intersection_master_coords);
            for (unsigned int i = 0; i < n_vars; ++i)
            {
                F_dof_map.dof_indices(elem, F_dof_indices(i), i);
            }

            get_values_for_interpolation(F_node, F_vec, F_dof_indices);
            get_values_for_interpolation(X_node, X_vec, X_dof_indices);
            for (unsigned int qp = 0; qp < intersection_master_coords.size(); ++qp)
            {
                const SideIndex<NDIM>& s_i = intersection_indices[qp];
                const int axis = s_i.getAxis();
                if (!side_boxes[axis].contains(s_i)) continue;
                if (spread_value_at_loc(s_i) != 0) continue;  // each value may be spread to only once
                spread_value_at_loc(s_i) = 1;
                jacobian(dX_ds,qp,X_node,dphi_X);
                const double J = std::abs(dX_ds.det());
                const double F_qp = interpolate(qp,F_node(blitz::Range::all(),axis),phi_F)/J;
                (*f_data)(s_i) += F_qp;
            }
        }
    }

    IBTK_TIMER_STOP(t_prolong_density);
    return;
}// prolongDensity

void
FEDataManager::interp(
    const int f_data_idx,
    NumericVector<double>& F_vec,
    NumericVector<double>& X_vec,
    const std::string& system_name,
    std::vector<Pointer<RefineSchedule<NDIM> > > f_refine_scheds,
    const double fill_data_time,
    const bool close_X)
{
    IBTK_TIMER_START(t_interp);

    // Extract the mesh.
    const MeshBase& mesh = d_es->get_mesh();
    const int dim = mesh.mesh_dimension();

    // Extract the FE systems and DOF maps, and setup the FE objects.
    System& F_system = d_es->get_system(system_name);
    const unsigned int n_vars = F_system.n_vars();
    const DofMap& F_dof_map = F_system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned i = 0; i < n_vars; ++i) TBOX_ASSERT(F_dof_map.variable_type(i) == F_dof_map.variable_type(0));
#endif
    blitz::Array<std::vector<unsigned int>,1> F_dof_indices(n_vars);
    for (unsigned int i = 0; i < n_vars; ++i) F_dof_indices(i).reserve(NDIM == 2 ? 9 : 27);
    AutoPtr<FEBase> F_fe(FEBase::build(dim, F_dof_map.variable_type(0)));
    F_fe->attach_quadrature_rule(d_qrule);
    const std::vector<double>& JxW_F = F_fe->get_JxW();
    const std::vector<std::vector<double> >& phi_F = F_fe->get_phi();

    System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_dof_map.variable_type(0));
#endif
    blitz::Array<std::vector<unsigned int>,1> X_dof_indices(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) X_dof_indices(d).reserve(NDIM == 2 ? 9 : 27);
    AutoPtr<FEBase> X_fe(FEBase::build(dim, X_dof_map.variable_type(0)));
    X_fe->attach_quadrature_rule(d_qrule);
    const std::vector<std::vector<double> >& phi_X = X_fe->get_phi();

    // Communicate any unsynchronized ghost data and enforce any constraints.
    for (unsigned int k = 0; k < f_refine_scheds.size(); ++k)
    {
        if (!f_refine_scheds[k].isNull()) f_refine_scheds[k]->fillData(fill_data_time);
    }

    if (close_X) X_vec.close();
    X_dof_map.enforce_constraints_exactly(X_system, &X_vec);

    // Loop over the patches to interpolate values to the element quadrature
    // points from the grid, then use these values to compute the projection of
    // the interpolated velocity field onto the FE basis functions.
    AutoPtr<NumericVector<double> > F_rhs_vec = F_vec.zero_clone();
    std::vector<DenseVector<double> > F_rhs_e(n_vars);
    blitz::Array<double,2> X_node;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_number);
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const blitz::Array<Elem*,1>& patch_elems = d_active_patch_elem_map(local_patch_num);
        const unsigned int num_active_patch_elems = patch_elems.size();
        if (num_active_patch_elems == 0) continue;

        // Setup vectors to store the values of F and X at the quadrature
        // points.  We compute a conservative upper bound on the number of
        // quadrature points to try to avoid unnecessary reallocations.
        if (dim == 2)
        {
            d_qrule->init(QUAD9);
        }
        if (dim == 3)
        {
            d_qrule->init(HEX27);
        }
        const unsigned int n_qp_estimate = d_qrule->n_points();
        static const unsigned int safety_factor = 2;
        std::vector<double> F_qp(safety_factor*n_vars*n_qp_estimate*num_active_patch_elems);
        std::vector<double> X_qp(safety_factor*NDIM  *n_qp_estimate*num_active_patch_elems);

        // Loop over the elements and compute the positions of the quadrature points.
        int qp_offset = 0;
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            const Elem* const elem = patch_elems(e_idx);

            X_fe->reinit(elem);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(elem, X_dof_indices(d), d);
            }

            const unsigned int n_qp = d_qrule->n_points();
            if (UNLIKELY(F_qp.size() < n_vars*(qp_offset+n_qp))) F_qp.resize(n_vars*(qp_offset+n_qp));
            if (UNLIKELY(X_qp.size() < NDIM  *(qp_offset+n_qp))) X_qp.resize(NDIM  *(qp_offset+n_qp));

            get_values_for_interpolation(X_node, X_vec, X_dof_indices);
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const int idx = NDIM*(qp+qp_offset);
                interpolate(&X_qp[idx],qp,X_node,phi_X);
            }

            qp_offset += n_qp;
        }

        if (qp_offset == 0) continue;

        F_qp.resize(n_vars*qp_offset);
        X_qp.resize(NDIM  *qp_offset);

        // Interpolate values from the Cartesian grid patch to the quadrature
        // points.
        //
        // NOTE: Values are interpolated only to those quadrature points that
        // are within the patch interior.
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& interp_box = patch->getBox();
        Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_data_idx);
        Pointer<CellData<NDIM,double> > f_cc_data = f_data;
        Pointer<SideData<NDIM,double> > f_sc_data = f_data;
        const bool is_cc_data = !f_cc_data.isNull();
        const bool is_sc_data = !f_sc_data.isNull();
        if (is_cc_data) LEInteractor::interpolate(F_qp, n_vars, X_qp, NDIM, f_cc_data, patch, interp_box, d_interp_weighting_fcn);
        if (is_sc_data) LEInteractor::interpolate(F_qp, n_vars, X_qp, NDIM, f_sc_data, patch, interp_box, d_interp_weighting_fcn);

        // Loop over the elements and accumulate the right-hand-side values.
        qp_offset = 0;
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            const Elem* const elem = patch_elems(e_idx);

            F_fe->reinit(elem);
            for (unsigned int i = 0; i < n_vars; ++i)
            {
                F_dof_map.dof_indices(elem, F_dof_indices(i), i);
                if (F_rhs_e[i].size() != F_dof_indices(i).size())
                {
                    F_rhs_e[i].resize(F_dof_indices(i).size());  // NOTE: DenseVector::resize() automatically zeroes the vector contents.
                }
                else
                {
                    F_rhs_e[i].zero();
                }
            }

            const unsigned int n_qp = d_qrule->n_points();
            const unsigned int n_basis = F_dof_indices(0).size();

            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const int idx = n_vars*(qp+qp_offset);
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    const double phi_JxW_F = phi_F[k][qp]*JxW_F[qp];
                    for (unsigned int i = 0; i < n_vars; ++i)
                    {
                        F_rhs_e[i](k) += F_qp[idx+i]*phi_JxW_F;
                    }
                }
            }

            for (unsigned int i = 0; i < n_vars; ++i)
            {
                F_dof_map.constrain_element_vector(F_rhs_e[i], F_dof_indices(i));
                F_rhs_vec->add_vector(F_rhs_e[i], F_dof_indices(i));
            }

            qp_offset += n_qp;
        }
    }

    // Solve for the nodal values.
    computeL2Projection(F_vec, *F_rhs_vec, system_name, d_interp_uses_consistent_mass_matrix);

    IBTK_TIMER_STOP(t_interp);
    return;
}// interp

void
FEDataManager::restrictValue(
    const int f_data_idx,
    NumericVector<double>& F_vec,
    NumericVector<double>& X_vec,
    const std::string& system_name,
    const bool close_X)
{
    IBTK_TIMER_START(t_restrict_value);

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
    const int dim = mesh.mesh_dimension();

    // Extract the FE systems and DOF maps, and setup the FE objects.
    System& F_system = d_es->get_system(system_name);
    const unsigned int n_vars = F_system.n_vars();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(n_vars == NDIM);  // specialized to side-centered data
#endif
    const DofMap& F_dof_map = F_system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned i = 0; i < n_vars; ++i) TBOX_ASSERT(F_dof_map.variable_type(i) == F_dof_map.variable_type(0));
#endif
    blitz::Array<std::vector<unsigned int>,1> F_dof_indices(n_vars);
    for (unsigned int i = 0; i < n_vars; ++i) F_dof_indices(i).reserve(NDIM == 2 ? 9 : 27);
    AutoPtr<FEBase> F_fe(FEBase::build(dim, F_dof_map.variable_type(0)));
    const std::vector<std::vector<double> >& phi_F = F_fe->get_phi();

    System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_dof_map.variable_type(0));
#endif
    FEType X_fe_type = X_dof_map.variable_type(0);
    blitz::Array<std::vector<unsigned int>,1> X_dof_indices(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) X_dof_indices(d).reserve(NDIM == 2 ? 9 : 27);
    AutoPtr<FEBase> X_fe(FEBase::build(dim, X_dof_map.variable_type(0)));
    const std::vector<std::vector<VectorValue<double> > >& dphi_X = X_fe->get_dphi();

    // Communicate any unsynchronized ghost data and enforce any constraints.
    if (close_X) X_vec.close();
    X_dof_map.enforce_constraints_exactly(X_system, &X_vec);

    // Loop over the patches to assemble the right-hand-side vector used to
    // solve for F.
    AutoPtr<NumericVector<double> > F_rhs_vec = F_vec.zero_clone();
    std::vector<DenseVector<double> > F_rhs_e(n_vars);
    TensorValue<double> dX_ds;
    blitz::Array<double,2> X_node;
    static const unsigned int MAX_NODES = (NDIM == 2 ? 9 : 27);
    Point s_node_cache[MAX_NODES], X_node_cache[MAX_NODES];
    blitz::TinyVector<double,NDIM> X_min, X_max;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_number);
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const blitz::Array<Elem*,1>& patch_elems = d_active_patch_elem_map(local_patch_num);
        const unsigned int num_active_patch_elems = patch_elems.size();
        if (num_active_patch_elems == 0) continue;

        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<SideData<NDIM,double> > f_data = patch->getPatchData(f_data_idx);
        const Box<NDIM>& patch_box = patch->getBox();
        const CellIndex<NDIM>& patch_lower = patch_box.lower();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_x_lower = patch_geom->getXLower();
        const double* const patch_dx = patch_geom->getDx();
        double dx = 1.0;
        for (unsigned int d = 0; d < NDIM; ++d) dx *= patch_dx[d];

        blitz::TinyVector<Box<NDIM>,NDIM> side_boxes;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            side_boxes[axis] = SideGeometry<NDIM>::toSideBox(patch_box,axis);
            if (!patch_geom->getTouchesRegularBoundary(axis,1)) side_boxes[axis].growUpper(axis,-1);
        }

        SideData<NDIM,int> interpolated_value_at_loc(patch_box, 1, IntVector<NDIM>(0));
        interpolated_value_at_loc.fillAll(0);

        // Loop over the elements.
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems(e_idx);
            const unsigned int n_node = elem->n_nodes();
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(elem, X_dof_indices(d), d);
            }

            // Cache the nodal and physical coordinates of the element,
            // determine the bounding box of the current configuration of the
            // element, and set the nodal coordinates of the element to
            // correspond to the physical coordinates.
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(n_node <= MAX_NODES);
#endif
            X_min =  0.5*std::numeric_limits<double>::max();
            X_max = -0.5*std::numeric_limits<double>::max();
            for (unsigned int k = 0; k < n_node; ++k)
            {
                s_node_cache[k] = elem->point(k);
                for (int d = 0; d < NDIM; ++d)
                {
                    X_node_cache[k](d) = X_vec(X_dof_indices(d)[k]);
                    X_min[d] = std::min(X_min[d],X_node_cache[k](d));
                    X_max[d] = std::max(X_max[d],X_node_cache[k](d));
                }
                elem->point(k) = X_node_cache[k];
            }

            // Loop over coordinate directions and look for Eulerian grid points
            // that are covered by the element.
            std::vector<Point>            intersection_master_coords;
            std::vector<SideIndex<NDIM> > intersection_indices;
            static const int estimated_max_size = (NDIM == 2 ? 64 : 512);
            intersection_master_coords.reserve(estimated_max_size);
            intersection_indices      .reserve(estimated_max_size);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                // Loop over the relevant range of indices.
                blitz::TinyVector<int,NDIM> i_begin, i_end, ic;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    if (d == axis)
                    {
                        i_begin[d] = std::ceil((X_min[d]-patch_x_lower[d])/patch_dx[d]) + patch_lower[d];
                        i_end  [d] = std::ceil((X_max[d]-patch_x_lower[d])/patch_dx[d]) + patch_lower[d];
                    }
                    else
                    {
                        i_begin[d] = std::ceil((X_min[d]-patch_x_lower[d])/patch_dx[d] - 0.5) + patch_lower[d];
                        i_end  [d] = std::ceil((X_max[d]-patch_x_lower[d])/patch_dx[d] - 0.5) + patch_lower[d];
                    }
                }
#if (NDIM == 3)
                for (ic[2] = i_begin[2]; ic[2] < i_end[2]; ++ic[2])
                {
#endif
                    for (ic[1] = i_begin[1]; ic[1] < i_end[1]; ++ic[1])
                    {
                        for (ic[0] = i_begin[0]; ic[0] < i_end[0]; ++ic[0])
                        {
                            Point p;
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                p(d) = patch_x_lower[d] + patch_dx[d]*(static_cast<double>(ic[d]-patch_lower[d])+(d == axis ? 0.0 : 0.5));
                            }
                            const Point master_coords = FEInterface::inverse_map(dim, X_fe_type, elem, p, TOLERANCE, false);
                            if (FEInterface::on_reference_element(master_coords,elem->type()))
                            {
                                intersection_master_coords.push_back(master_coords);
#if (NDIM == 2)
                                SideIndex<NDIM> s_i(Index<NDIM>(ic[0],ic[1]),axis,0);
#endif
#if (NDIM == 3)
                                SideIndex<NDIM> s_i(Index<NDIM>(ic[0],ic[1],ic[2]),axis,0);
#endif
                                intersection_indices.push_back(s_i);
                            }
                        }
                    }
#if (NDIM == 3)
                }
#endif
            }

            // Restore the nodal coordinates.
            for (unsigned int k = 0; k < n_node; ++k)
            {
                elem->point(k) = s_node_cache[k];
            }

            // If there are no intersection points, then continue on to the next
            // element.
            if (intersection_master_coords.empty()) continue;

            // Evaluate the Eulerian value and rescale it by 1.0/det(dX/ds).
            F_fe->reinit(elem, &intersection_master_coords);
            X_fe->reinit(elem, &intersection_master_coords);
            for (unsigned int i = 0; i < n_vars; ++i)
            {
                F_dof_map.dof_indices(elem, F_dof_indices(i), i);
                if (F_rhs_e[i].size() != F_dof_indices(i).size())
                {
                    F_rhs_e[i].resize(F_dof_indices(i).size());  // NOTE: DenseVector::resize() automatically zeroes the vector contents.
                }
                else
                {
                    F_rhs_e[i].zero();
                }
            }

            const unsigned int n_basis = F_dof_indices(0).size();

            get_values_for_interpolation(X_node, X_vec, X_dof_indices);
            for (unsigned int qp = 0; qp < intersection_master_coords.size(); ++qp)
            {
                const SideIndex<NDIM>& s_i = intersection_indices[qp];
                const int axis = s_i.getAxis();
                if (!side_boxes[axis].contains(s_i)) continue;
                if (interpolated_value_at_loc(s_i) != 0) continue;  // each value may be interpolated only once
                interpolated_value_at_loc(s_i) = 1;
                jacobian(dX_ds,qp,X_node,dphi_X);
                const double J = std::abs(dX_ds.det());
                const double F_qp = (*f_data)(s_i)*dx/J;
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    F_rhs_e[axis](k) += F_qp*phi_F[k][qp];
                }
            }

            for (unsigned int i = 0; i < n_vars; ++i)
            {
                F_dof_map.constrain_element_vector(F_rhs_e[i], F_dof_indices(i));
                F_rhs_vec->add_vector(F_rhs_e[i], F_dof_indices(i));
            }
        }
    }

    // Solve for the nodal values.
    computeL2Projection(F_vec, *F_rhs_vec, system_name, d_interp_uses_consistent_mass_matrix);

    IBTK_TIMER_STOP(t_restrict_value);
    return;
}// restrictValue

std::pair<LinearSolver<double>*,SparseMatrix<double>*>
FEDataManager::buildL2ProjectionSolver(
    const std::string& system_name,
    const QuadratureType quad_type,
    const Order quad_order)
{
    IBTK_TIMER_START(t_build_l2_projection_solver);

    if ((d_L2_proj_solver.count(system_name) == 0 || d_L2_proj_matrix.count(system_name) == 0) ||
        (d_L2_proj_quad_type[system_name] != quad_type) || (d_L2_proj_quad_order[system_name] != quad_order))
    {
        const MeshBase& mesh = d_es->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();
        AutoPtr<QBase> qrule = QBase::build(quad_type, dim, quad_order);

        System& system = d_es->get_system(system_name);
        const int sys_num = system.number();
        DofMap& dof_map = system.get_dof_map();
        // Sparsity patterns are not automatically computed for all system
        // types.  If one has not been computed for this system, we compute it
        // now.
        if (dof_map.get_n_nz().size() != dof_map.n_local_dofs())
        {
            dof_map.compute_sparsity(mesh);
        }
        std::vector<unsigned int> dof_indices;
        AutoPtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<double>& JxW = fe->get_JxW();
        const std::vector<std::vector<double> >& phi = fe->get_phi();

        LinearSolver<double>* solver = LinearSolver<double>::build().release();
        solver->init();

        SparseMatrix<double>* M_mat = SparseMatrix<double>::build().release();
        M_mat->attach_dof_map(dof_map);
        M_mat->init();

        DenseMatrix<double> M_e;

        // Loop over the mesh to construct the system matrix.
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end   = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            fe->reinit(elem);
            for (unsigned int var_num = 0; var_num < dof_map.n_variables(); ++var_num)
            {
                dof_map.dof_indices(elem, dof_indices, var_num);
                M_e.resize(dof_indices.size(), dof_indices.size());
                const unsigned int n_basis = dof_indices.size();
                const unsigned int n_qp = qrule->n_points();
                for (unsigned int i = 0; i < n_basis; ++i)
                {
                    for (unsigned int j = 0; j < n_basis; ++j)
                    {
                        for (unsigned int qp = 0; qp < n_qp; ++qp)
                        {
                            M_e(i,j) += (phi[i][qp]*phi[j][qp])*JxW[qp];
                        }
                    }
                }
                dof_map.constrain_element_matrix(M_e, dof_indices);
                M_mat->add_matrix(M_e, dof_indices);
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
                if (elem->neighbor(side) != NULL) continue;
                const std::vector<short int>& bdry_ids = mesh.boundary_info->boundary_ids(elem, side);
                const bool at_dirichlet_bdry = std::find(bdry_ids.begin(), bdry_ids.end(), DIRICHLET_BDRY_ID) != bdry_ids.end();
                if (!at_dirichlet_bdry) continue;
                fe->reinit(elem);
                for (unsigned int n = 0; n < elem->n_nodes(); ++n)
                {
                    if (elem->is_node_on_side(n, side))
                    {
                        Node* node = elem->get_node(n);
                        for (unsigned int var_num = 0; var_num < dof_map.n_variables(); ++var_num)
                        {
                            dof_map.dof_indices(elem, dof_indices, var_num);
                            const unsigned int n_comp = node->n_comp(sys_num, var_num);
                            for (unsigned int comp = 0; comp < n_comp; ++comp)
                            {
                                const unsigned int node_dof_index = node->dof_number(sys_num, var_num, comp);
                                if (dof_map.is_constrained_dof(node_dof_index))
                                {
                                    for (std::vector<unsigned int>::const_iterator cit = dof_indices.begin(); cit != dof_indices.end(); ++cit)
                                    {
                                        const unsigned int k = *cit;
                                        M_mat->set(node_dof_index, k, (node_dof_index == k ? 1.0 : 0.0));
                                        M_mat->set(k, node_dof_index, (node_dof_index == k ? 1.0 : 0.0));
                                    }
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
        d_L2_proj_solver[system_name] = solver;
        d_L2_proj_matrix[system_name] = M_mat;
        d_L2_proj_quad_type[system_name] = quad_type;
        d_L2_proj_quad_order[system_name] = quad_order;
    }

    IBTK_TIMER_STOP(t_build_l2_projection_solver);
    return std::make_pair(d_L2_proj_solver[system_name], d_L2_proj_matrix[system_name]);
}// buildL2ProjectionSolver

NumericVector<double>*
FEDataManager::buildDiagonalL2MassMatrix(
    const std::string& system_name,
    const QuadratureType quad_type,
    const Order quad_order)
{
    IBTK_TIMER_START(t_build_diagonal_l2_mass_matrix);

    if (d_L2_proj_matrix_diag.count(system_name) == 0)
    {
        const MeshBase& mesh = d_es->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();
        AutoPtr<QBase> qrule = QBase::build(quad_type, dim, quad_order);

        System& system = d_es->get_system(system_name);
        const int sys_num = system.number();
        DofMap& dof_map = system.get_dof_map();
        std::vector<unsigned int> dof_indices;
        AutoPtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<double>& JxW = fe->get_JxW();
        const std::vector<std::vector<double> >& phi = fe->get_phi();

        NumericVector<double>* M_vec = system.solution->zero_clone().release();
        DenseVector<double> M_diag_e;

        // Loop over the mesh to construct the (diagonal) system matrix.
        //
        // We construct diagonal elemental mass matrices by taking the diagonal
        // part of the consistent elemental mass matrix and rescaling it so that
        // it has the correct total elemental mass.
        //
        // Ref: E. Hinton, T. Rock and O.C. Zienkiewicz.  A note on mass lumping
        // and related processes in the finite element method.  Earthquake Eng
        // Struct Dyn.  4:245--249 (1976).
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end   = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            if (elem->default_order() != FIRST)
            {
                IBTK_DO_ONCE(
                        pout << "WARNING: use of diagonal mass matrices is not recommended for higher-order elements.\n"
                             << "         suggest using consistent mass matrix instead.\n"
                             );
            }
            fe->reinit(elem);
            const unsigned int n_qp = qrule->n_points();
            double elem_mass = 0.0;
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                elem_mass += JxW[qp];
            }
            for (unsigned int var_num = 0; var_num < dof_map.n_variables(); ++var_num)
            {
                dof_map.dof_indices(elem, dof_indices, var_num);
                M_diag_e.resize(dof_indices.size());
                const unsigned int n_basis = dof_indices.size();
                double diagonal_mass = 0.0;
                for (unsigned int i = 0; i < n_basis; ++i)
                {
                    for (unsigned int qp = 0; qp < n_qp; ++qp)
                    {
                        M_diag_e(i) += (phi[i][qp]*phi[i][qp])*JxW[qp];
                    }
                    diagonal_mass += M_diag_e(i);
                }
                M_diag_e *= elem_mass/diagonal_mass;
                dof_map.constrain_element_vector(M_diag_e, dof_indices);
                M_vec->add_vector(M_diag_e, dof_indices);
            }
        }

        // Flush assemble the vector representation of the diagonal matrix.
        M_vec->close();

        // Reset values at Dirichlet boundaries.
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                if (elem->neighbor(side) != NULL) continue;
                const std::vector<short int>& bdry_ids = mesh.boundary_info->boundary_ids(elem, side);
                const bool at_dirichlet_bdry = std::find(bdry_ids.begin(), bdry_ids.end(), DIRICHLET_BDRY_ID) != bdry_ids.end();
                if (!at_dirichlet_bdry) continue;
                for (unsigned int n = 0; n < elem->n_nodes(); ++n)
                {
                    if (elem->is_node_on_side(n, side))
                    {
                        Node* node = elem->get_node(n);
                        for (unsigned int var_num = 0; var_num < dof_map.n_variables(); ++var_num)
                        {
                            const unsigned int n_comp = node->n_comp(sys_num, var_num);
                            for (unsigned int comp = 0; comp < n_comp; ++comp)
                            {
                                const unsigned int node_dof_index = node->dof_number(sys_num, var_num, comp);
                                if (dof_map.is_constrained_dof(node_dof_index)) M_vec->set(node_dof_index, 1.0);
                            }
                        }
                    }
                }
            }
        }

        // Assemble the vector representation of the diagonal matrix.
        M_vec->close();

        // Store the diagonal mass matrix.
        d_L2_proj_matrix_diag[system_name] = M_vec;
    }

    IBTK_TIMER_STOP(t_build_diagonal_l2_mass_matrix);
    return d_L2_proj_matrix_diag[system_name];
}// buildDiagonalL2MassMatrix

bool
FEDataManager::computeL2Projection(
    NumericVector<double>& U_vec,
    NumericVector<double>& F_vec,
    const std::string& system_name,
    const bool consistent_mass_matrix,
    const QuadratureType quad_type,
    const Order quad_order,
    const double tol,
    const unsigned int max_its)
{
    IBTK_TIMER_START(t_compute_l2_projection);

    int ierr;
    bool converged = false;

    F_vec.close();
    const System& system = d_es->get_system(system_name);
    const DofMap& dof_map = system.get_dof_map();
    dof_map.enforce_constraints_exactly(system, &F_vec);
    if (consistent_mass_matrix)
    {
        std::pair<libMesh::LinearSolver<double>*,SparseMatrix<double>*> proj_solver_components = buildL2ProjectionSolver(system_name, quad_type, quad_order);
        PetscLinearSolver<double>* solver = dynamic_cast<PetscLinearSolver<double>*>(proj_solver_components.first);
        PetscMatrix<double>* M_mat = dynamic_cast<PetscMatrix<double>*>(proj_solver_components.second);
        solver->solve(*M_mat, *M_mat, U_vec, F_vec, tol, max_its);
        KSPConvergedReason reason;
        ierr = KSPGetConvergedReason(solver->ksp(), &reason); IBTK_CHKERRQ(ierr);
        converged = reason > 0;
    }
    else
    {
        PetscVector<double>* M_diag_vec = dynamic_cast<PetscVector<double>*>(buildDiagonalL2MassMatrix(system_name,quad_type,quad_order));
        Vec M_diag_petsc_vec = M_diag_vec->vec();
        Vec U_petsc_vec = dynamic_cast<PetscVector<double>*>(&U_vec)->vec();
        Vec F_petsc_vec = dynamic_cast<PetscVector<double>*>(&F_vec)->vec();
        ierr = VecPointwiseDivide(U_petsc_vec, F_petsc_vec, M_diag_petsc_vec); IBTK_CHKERRQ(ierr);
        converged = true;
    }
    dof_map.enforce_constraints_exactly(system, &U_vec);

    IBTK_TIMER_STOP(t_compute_l2_projection);
    return converged;
}// computeL2Projection

void
FEDataManager::updateWorkloadData(
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

    // Workload estimates are computed only on the level to which the FE mesh
    // has been assigned.
    if (!d_load_balancer.isNull())
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (ln == d_level_number)
            {
                updateQuadPointCountData(ln,ln);
                HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(d_hierarchy,ln,ln);
                hier_cc_data_ops.setToScalar(d_workload_idx, 1.0);
                hier_cc_data_ops.add(d_workload_idx, d_qp_count_idx, d_workload_idx);
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
FEDataManager::initializeLevelData(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool /*can_be_refined*/,
    const bool /*initial_time*/,
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

    // Allocate storage needed to initialize the level.
    //
    // Since time gets set when we allocate data, re-stamp it to current time if
    // we don't need to allocate.
    if (allocate_data)
    {
        level->allocatePatchData(d_workload_idx, init_data_time);
    }
    else
    {
        level->setTime(d_workload_idx, init_data_time);
    }

    // Initialize workload data and setup the load balancer.
    if (!d_load_balancer.isNull())
    {
        HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(hierarchy,level_number,level_number);
        hier_cc_data_ops.setToScalar(d_workload_idx, 1.0);
        if (!old_level.isNull() && level_number == d_level_number)
        {
            Pointer<RefineOperator<NDIM> > regrid_fill_op = Pointer<RefineOperator<NDIM> >(NULL);
            Pointer<RefineAlgorithm<NDIM> > regrid_fill_alg = new RefineAlgorithm<NDIM>();
            regrid_fill_alg->registerRefine(d_workload_idx, d_workload_idx, d_workload_idx, regrid_fill_op);
            regrid_fill_alg->createSchedule(level, old_level)->fillData(init_data_time);
        }
        if (level_number == d_level_number)
        {
            d_load_balancer->setWorkloadPatchDataIndex(d_workload_idx, level_number);
        }
        else
        {
            d_load_balancer->setUniformWorkload(level_number);
        }
    }

    IBTK_TIMER_STOP(t_initialize_level_data);
    return;
}// initializeLevelData

void
FEDataManager::resetHierarchyConfiguration(
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
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif

    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();

    // Reset the patch hierarchy and levels.
    setPatchHierarchy(hierarchy);
    resetLevels(0,finest_hier_level);

    IBTK_TIMER_STOP(t_reset_hierarchy_configuration);
    return;
}// resetHierarchyConfiguration

void
FEDataManager::applyGradientDetector(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double /*error_data_time*/,
    const int tag_index,
    const bool initial_time,
    const bool /*uses_richardson_extrapolation_too*/)
{
    if (level_number >= d_level_number) return;

    IBTK_TIMER_START(t_apply_gradient_detector);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    if (initial_time)
    {
        // Determine the active elements associated with the prescribed patch
        // level.
        blitz::Array<blitz::Array<Elem*,1>,1> active_level_elem_map;
        const IntVector<NDIM> ghost_width = 1;
        collectActivePatchElements(active_level_elem_map, level_number, ghost_width);
        std::vector<unsigned int> X_ghost_dofs;
        blitz::Array<Elem*,1> active_level_elems;
        flatten(active_level_elems, active_level_elem_map);
        collectGhostDOFIndices(X_ghost_dofs, active_level_elems, COORDINATES_SYSTEM_NAME);

        const MeshBase& mesh = d_es->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
        const DofMap& X_dof_map = X_system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
        for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_dof_map.variable_type(0));
#endif
        blitz::Array<std::vector<unsigned int>,1> X_dof_indices(dim);
        AutoPtr<FEBase> X_fe(FEBase::build(dim, X_dof_map.variable_type(0)));
        X_fe->attach_quadrature_rule(d_qrule);
        const std::vector<std::vector<double> >& phi_X = X_fe->get_phi();
        NumericVector<double>* X_vec = getCoordsVector();
        AutoPtr<NumericVector<double> > X_ghost_vec = NumericVector<double>::build();
        X_ghost_vec->init(X_vec->size(), X_vec->local_size(), X_ghost_dofs, true, GHOSTED);
        X_vec->localize(*X_ghost_vec);
        X_dof_map.enforce_constraints_exactly(X_system, X_ghost_vec.get());

        // Tag cells for refinement whenever they contain active element
        // quadrature points.
        blitz::Array<double,2> X_node;
        blitz::TinyVector<double,NDIM> X_qp;
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            // The relevant collection of elements.
            const blitz::Array<Elem*,1>& patch_elems = active_level_elem_map(local_patch_num);
            const unsigned int num_active_patch_elems = patch_elems.size();
            if (num_active_patch_elems == 0) continue;

            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const CellIndex<NDIM>& patch_lower = patch_box.lower();
            const CellIndex<NDIM>& patch_upper = patch_box.upper();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();
            const double* const patch_dx = patch_geom->getDx();

            Pointer<CellData<NDIM,int> > tag_data = patch->getPatchData(tag_index);

            for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                const Elem* const elem = patch_elems(e_idx);
                X_fe->reinit(elem);
                for (unsigned int d = 0; d < dim; ++d)
                {
                    X_dof_map.dof_indices(elem, X_dof_indices(d), d);
                }
                get_values_for_interpolation(X_node, *X_ghost_vec.get(), X_dof_indices);
                for (unsigned int qp = 0; qp < d_qrule->n_points(); ++qp)
                {
                    interpolate(&X_qp[0], qp, X_node, phi_X);
                    const Index<NDIM> i = IndexUtilities::getCellIndex(X_qp, patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper);
                    tag_data->fill(1,Box<NDIM>(i-Index<NDIM>(1),i+Index<NDIM>(1)));
                }
            }
        }
    }
    else if (level_number+1 == d_level_number && level_number < d_hierarchy->getFinestLevelNumber())
    {
        Pointer<PatchLevel<NDIM> > finer_level = d_hierarchy->getPatchLevel(level_number+1);
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);

        // Update the node count data and coarsen it from the finer level.
        updateQuadPointCountData(level_number,level_number+1);
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
            Pointer<CellData<NDIM,int> > tag_data = patch->getPatchData(tag_index);
            Pointer<CellData<NDIM,double> > qp_count_data = patch->getPatchData(d_qp_count_idx);
            for (CellIterator<NDIM> ic(patch_box); ic; ic++)
            {
                const CellIndex<NDIM>& i = ic();
                if ((*qp_count_data)(i) > 0.0)
                {
                    (*tag_data)(i) = 1;
                }
            }
        }
    }

    IBTK_TIMER_STOP(t_apply_gradient_detector);
    return;
}// applyGradientDetector

void
FEDataManager::putToDatabase(
    Pointer<Database> db)
{
    IBTK_TIMER_START(t_put_to_database);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("FE_DATA_MANAGER_VERSION", FE_DATA_MANAGER_VERSION);

    db->putInteger("d_coarsest_ln", d_coarsest_ln);
    db->putInteger("d_finest_ln"  , d_finest_ln  );

    IBTK_TIMER_STOP(t_put_to_database);
    return;
}// putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

FEDataManager::FEDataManager(
    const std::string& object_name,
    const std::string& interp_weighting_fcn,
    const std::string& spread_weighting_fcn,
    QBase* const qrule,
    const bool interp_uses_consistent_mass_matrix,
    const IntVector<NDIM>& ghost_width,
    bool register_for_restart)
    : COORDINATES_SYSTEM_NAME("coordinates system"),
      d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_load_balancer(NULL),
      d_hierarchy(NULL),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_interp_weighting_fcn(interp_weighting_fcn),
      d_spread_weighting_fcn(spread_weighting_fcn),
      d_qrule(qrule),
      d_interp_uses_consistent_mass_matrix(interp_uses_consistent_mass_matrix),
      d_ghost_width(ghost_width),
      d_es(NULL),
      d_level_number(-1),
      d_active_patch_ghost_dofs(),
      d_L2_proj_solver(),
      d_L2_proj_matrix(),
      d_L2_proj_matrix_diag(),
      d_L2_proj_quad_type(),
      d_L2_proj_quad_order()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
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

    // Create/look up the variable context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name+"::CONTEXT");

    // Register the node count and workload variables with the VariableDatabase.
    d_qp_count_var = new CellVariable<NDIM,double>(d_object_name+"::qp_count");
    d_qp_count_idx = var_db->registerVariableAndContext(d_qp_count_var, d_context, 0);

    d_workload_var = new CellVariable<NDIM,double>(d_object_name+"::workload");
    d_workload_idx = var_db->registerVariableAndContext(d_workload_var, d_context, 0);

    // Setup Timers.
    IBTK_DO_ONCE(
        t_reinit_element_mappings = TimerManager::getManager()->getTimer("IBTK::FEDataManager::reinitElementMappings()");
        t_build_ghosted_solution_vector = TimerManager::getManager()->getTimer("IBTK::FEDataManager::buildGhostedSolutionVector()");
        t_spread = TimerManager::getManager()->getTimer("IBTK::FEDataManager::spread()");
        t_prolong_value = TimerManager::getManager()->getTimer("IBTK::FEDataManager::prolongValue()");
        t_prolong_density = TimerManager::getManager()->getTimer("IBTK::FEDataManager::prolongDensity()");
        t_interp = TimerManager::getManager()->getTimer("IBTK::FEDataManager::interp()");
        t_restrict_value = TimerManager::getManager()->getTimer("IBTK::FEDataManager::restrictValue()");
        t_build_l2_projection_solver = TimerManager::getManager()->getTimer("IBTK::FEDataManager::buildL2ProjectionSolver()");
        t_build_diagonal_l2_mass_matrix = TimerManager::getManager()->getTimer("IBTK::FEDataManager::buildDiagonalL2MassMatrix()");
        t_compute_l2_projection = TimerManager::getManager()->getTimer("IBTK::FEDataManager::computeL2Projection()");
        t_update_workload_data = TimerManager::getManager()->getTimer("IBTK::FEDataManager::updateWorkloadData()");
        t_initialize_level_data = TimerManager::getManager()->getTimer("IBTK::FEDataManager::initializeLevelData()");
        t_reset_hierarchy_configuration = TimerManager::getManager()->getTimer("IBTK::FEDataManager::resetHierarchyConfiguration()");
        t_apply_gradient_detector = TimerManager::getManager()->getTimer("IBTK::FEDataManager::applyGradientDetector()");
        t_put_to_database = TimerManager::getManager()->getTimer("IBTK::FEDataManager::putToDatabase()");
                 );
    return;
}// FEDataManager

FEDataManager::~FEDataManager()
{
    for (std::map<std::string,NumericVector<double>*>::iterator it = d_system_ghost_vec.begin();
         it != d_system_ghost_vec.end(); ++it)
    {
        delete it->second;
    }
    for (std::map<std::string,LinearSolver<double>*>::iterator it = d_L2_proj_solver.begin();
         it != d_L2_proj_solver.end(); ++it)
    {
        delete it->second;
    }
    for (std::map<std::string,SparseMatrix<double>*>::iterator it = d_L2_proj_matrix.begin();
         it != d_L2_proj_matrix.end(); ++it)
    {
        delete it->second;
    }
    for (std::map<std::string,NumericVector<double>*>::iterator it = d_L2_proj_matrix_diag.begin();
         it != d_L2_proj_matrix_diag.end(); ++it)
    {
        delete it->second;
    }
    return;
}// ~FEDataManager

/////////////////////////////// PRIVATE //////////////////////////////////////

void
FEDataManager::updateQuadPointCountData(
    const int coarsest_ln,
    const int finest_ln)
{
    // Set the node count data on the specified range of levels of the
    // hierarchy.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_qp_count_idx)) level->allocatePatchData(d_qp_count_idx);
        HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(d_hierarchy,ln,ln);
        hier_cc_data_ops.setToScalar(d_qp_count_idx, 0.0);
        if (ln != d_level_number) continue;

        const MeshBase& mesh = d_es->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
        const DofMap& X_dof_map = X_system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
        for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_dof_map.variable_type(0));
#endif
        blitz::Array<std::vector<unsigned int>,1> X_dof_indices(dim);
        AutoPtr<FEBase> X_fe(FEBase::build(dim, X_dof_map.variable_type(0)));
        X_fe->attach_quadrature_rule(d_qrule);
        const std::vector<std::vector<double> >& phi_X = X_fe->get_phi();
        NumericVector<double>* X_vec = getCoordsVector();
        NumericVector<double>* X_ghost_vec = buildGhostedCoordsVector();
        X_vec->localize(*X_ghost_vec);
        X_dof_map.enforce_constraints_exactly(X_system, X_ghost_vec);
        blitz::Array<double,2> X_node;
        blitz::TinyVector<double,NDIM> X_qp;
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            const blitz::Array<Elem*,1>& patch_elems = d_active_patch_elem_map(local_patch_num);
            const unsigned int num_active_patch_elems = patch_elems.size();
            if (num_active_patch_elems == 0) continue;

            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const CellIndex<NDIM>& patch_lower = patch_box.lower();
            const CellIndex<NDIM>& patch_upper = patch_box.upper();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();
            const double* const patch_dx = patch_geom->getDx();

            Pointer<CellData<NDIM,double> > qp_count_data = patch->getPatchData(d_qp_count_idx);
            for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                const Elem* const elem = patch_elems(e_idx);
                X_fe->reinit(elem);
                for (unsigned int d = 0; d < dim; ++d)
                {
                    X_dof_map.dof_indices(elem, X_dof_indices(d), d);
                }
                get_values_for_interpolation(X_node, *X_ghost_vec, X_dof_indices);
                for (unsigned int qp = 0; qp < d_qrule->n_points(); ++qp)
                {
                    interpolate(&X_qp[0], qp, X_node, phi_X);
                    const Index<NDIM> i = IndexUtilities::getCellIndex(X_qp, patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper);
                    if (patch_box.contains(i)) (*qp_count_data)(i) += 1.0;
                }
            }
        }
    }
    return;
}// updateQuadPointCountData

blitz::Array<std::pair<blitz::TinyVector<double,NDIM>,blitz::TinyVector<double,NDIM> >,1>*
FEDataManager::computeActiveElementBoundingBoxes()
{
    const MeshBase& mesh = d_es->get_mesh();
    const unsigned int n_elem = mesh.max_elem_id()+1;
    System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    const DofMap& X_dof_map = X_system.get_dof_map();
    NumericVector<double>& X_vec = *X_system.solution;
    NumericVector<double>& X_ghost_vec = *X_system.current_local_solution;
    X_vec.localize(X_ghost_vec);
    X_dof_map.enforce_constraints_exactly(X_system, &X_ghost_vec);

    // Compute the lower and upper bounds of all active local elements in the
    // mesh.  Assumes nodal basis functions.
    d_active_elem_bboxes.resize(n_elem);
    d_active_elem_bboxes = std::make_pair(blitz::TinyVector<double,NDIM>(0.0),blitz::TinyVector<double,NDIM>(0.0));
    MeshBase::const_element_iterator       el_it  = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for ( ; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const unsigned int elem_id = elem->id();
        blitz::TinyVector<double,NDIM>& elem_lower_bound = d_active_elem_bboxes(elem_id).first;
        blitz::TinyVector<double,NDIM>& elem_upper_bound = d_active_elem_bboxes(elem_id).second;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            elem_lower_bound[d] =  0.5*std::numeric_limits<double>::max();
            elem_upper_bound[d] = -0.5*std::numeric_limits<double>::max();
        }

        const unsigned int n_nodes = elem->n_nodes();
        std::vector<unsigned int> dof_indices;
        dof_indices.reserve(NDIM*n_nodes);
        for (unsigned int k = 0; k < n_nodes; ++k)
        {
            Node* node = elem->get_node(k);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(node->n_dofs(X_sys_num,0) > 0);
#endif
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dof_indices.push_back(node->dof_number(X_sys_num,d,0));
            }
        }
        std::vector<double> X_node;
        X_ghost_vec.get(dof_indices, X_node);
        for (unsigned int k = 0; k < n_nodes; ++k)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const double& X = X_node[k*NDIM+d];
                elem_lower_bound[d] = std::min(elem_lower_bound[d], X);
                elem_upper_bound[d] = std::max(elem_upper_bound[d], X);
            }
        }
    }

    // Parallel sum elem_lower_bound and elem_upper_bound so that each process
    // has access to the bounding box data for each active element in the mesh.
    std::vector<double> d_active_elem_bboxes_flattened(2*NDIM*n_elem);
    for (unsigned int e = 0; e < n_elem; ++e)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            d_active_elem_bboxes_flattened[ 2*e   *NDIM+d] = d_active_elem_bboxes(e).first [d];
            d_active_elem_bboxes_flattened[(2*e+1)*NDIM+d] = d_active_elem_bboxes(e).second[d];
        }
    }
    SAMRAI_MPI::sumReduction(&d_active_elem_bboxes_flattened[0], d_active_elem_bboxes_flattened.size());
    for (unsigned int e = 0; e < n_elem; ++e)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            d_active_elem_bboxes(e).first [d] = d_active_elem_bboxes_flattened[ 2*e   *NDIM+d];
            d_active_elem_bboxes(e).second[d] = d_active_elem_bboxes_flattened[(2*e+1)*NDIM+d];
        }
    }
    return &d_active_elem_bboxes;
}// computeActiveElementBoundingBoxes

void
FEDataManager::collectActivePatchElements(
    blitz::Array<blitz::Array<Elem*,1>,1>& active_patch_elems,
    const int level_number,
    const IntVector<NDIM>& ghost_width)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif

    // Get the necessary FE data.
    const MeshBase& mesh = d_es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_dof_map.variable_type(0));
#endif
    blitz::Array<std::vector<unsigned int>,1> X_dof_indices(dim);
    AutoPtr<FEBase> X_fe(FEBase::build(dim, X_dof_map.variable_type(0)));
    X_fe->attach_quadrature_rule(d_qrule);
    const std::vector<std::vector<double> >& phi_X = X_fe->get_phi();
    NumericVector<double>* X_vec = getCoordsVector();
    AutoPtr<NumericVector<double> > X_ghost_vec = NumericVector<double>::build();

    // Setup data structures used to assign elements to patches.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
    const int num_local_patches = level->getProcessorMapping().getNumberOfLocalIndices();
    blitz::Array<std::set<Elem*>,1>    local_patch_elems(num_local_patches);
    blitz::Array<std::set<Elem*>,1> nonlocal_patch_elems(num_local_patches);
    blitz::Array<std::set<Elem*>,1> frontier_patch_elems(num_local_patches);

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
        std::set<Elem*>& frontier_elems = frontier_patch_elems(local_patch_num);
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        blitz::TinyVector<double,NDIM> x_lower;  for (unsigned int d = 0; d < NDIM; ++d) x_lower[d] = pgeom->getXLower()[d];
        blitz::TinyVector<double,NDIM> x_upper;  for (unsigned int d = 0; d < NDIM; ++d) x_upper[d] = pgeom->getXUpper()[d];
        const double* const dx = pgeom->getDx();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            x_lower[d] -= dx[d]*ghost_width[d];
            x_upper[d] += dx[d]*ghost_width[d];
        }

        MeshBase::const_element_iterator       el_it  = mesh.active_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_elements_end();
        for ( ; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            const unsigned int elem_id = elem->id();
            const blitz::TinyVector<double,NDIM>& elem_lower_bound = d_active_elem_bboxes(elem_id).first;
            const blitz::TinyVector<double,NDIM>& elem_upper_bound = d_active_elem_bboxes(elem_id).second;
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
        blitz::Array<Elem*,1> frontier_elems;
        flatten(frontier_elems, frontier_patch_elems);
        collectGhostDOFIndices(X_ghost_dofs, frontier_elems, COORDINATES_SYSTEM_NAME);
        X_ghost_vec->init(X_vec->size(), X_vec->local_size(), X_ghost_dofs, true, GHOSTED);
        X_vec->localize(*X_ghost_vec);
        X_dof_map.enforce_constraints_exactly(X_system, X_ghost_vec.get());

        // Keep only those elements that have a quadrature point on the local
        // patch.
        blitz::Array<double,2> X_node;
        blitz::TinyVector<double,NDIM> X_qp;
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            const std::set<Elem*>& frontier_elems = frontier_patch_elems(local_patch_num);
            std::set<Elem*>&          local_elems =    local_patch_elems(local_patch_num);
            std::set<Elem*>&       nonlocal_elems = nonlocal_patch_elems(local_patch_num);
            if (frontier_elems.empty()) continue;

            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Box<NDIM>  ghost_box = Box<NDIM>::grow(patch_box, ghost_width);
            const CellIndex<NDIM>& patch_lower = patch_box.lower();
            const CellIndex<NDIM>& patch_upper = patch_box.upper();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();
            const double* const patch_dx = patch_geom->getDx();

            std::set<Elem*>::const_iterator       el_it  = frontier_elems.begin();
            const std::set<Elem*>::const_iterator el_end = frontier_elems.end();
            for ( ; el_it != el_end; ++el_it)
            {
                Elem* const elem = *el_it;
                X_fe->reinit(elem);
                for (unsigned int d = 0; d < dim; ++d)
                {
                    X_dof_map.dof_indices(elem, X_dof_indices(d), d);
                }

                get_values_for_interpolation(X_node, *X_ghost_vec, X_dof_indices);
                bool found_qp = false;
                for (unsigned int qp = 0; qp < d_qrule->n_points() && !found_qp; ++qp)
                {
                    interpolate(&X_qp[0], qp, X_node, phi_X);
                    const Index<NDIM> i = IndexUtilities::getCellIndex(X_qp, patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper);
                    if (ghost_box.contains(i))
                    {
                        local_elems.insert(elem);
                        found_qp = true;
                    }
                }
                if (!found_qp) nonlocal_elems.insert(elem);
            }
        }

        // Rebuild the set of frontier elements, which are any neighbors of a
        // local element that has not already been determined to be either a
        // local or a nonlocal element.
        bool new_frontier = false;
        local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            std::set<Elem*>&       frontier_elems = frontier_patch_elems(local_patch_num);
            const std::set<Elem*>&    local_elems =    local_patch_elems(local_patch_num);
            const std::set<Elem*>& nonlocal_elems = nonlocal_patch_elems(local_patch_num);
            frontier_elems.clear();
            if (local_elems.empty()) continue;

            for (std::set<Elem*>::const_iterator cit = local_elems.begin();
                 cit != local_elems.end(); ++cit)
            {
                const Elem* const elem = *cit;
                for (unsigned int n = 0; n < elem->n_neighbors(); ++n)
                {
                    Elem* const nghbr_elem = elem->neighbor(n);
                    if (nghbr_elem == NULL) continue;
                    const bool    is_local_elem =    local_elems.find(nghbr_elem) !=    local_elems.end();
                    const bool is_nonlocal_elem = nonlocal_elems.find(nghbr_elem) != nonlocal_elems.end();
                    if (!(is_local_elem || is_nonlocal_elem))
                    {
                        frontier_elems.insert(nghbr_elem);
                        new_frontier = true;
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
        blitz::Array<Elem*,1>& active_elems = active_patch_elems(local_patch_num);
        const std::set<Elem*>&  local_elems =  local_patch_elems(local_patch_num);
        active_elems.resize(local_elems.size());
        int k = 0;
        for (std::set<Elem*>::const_iterator cit = local_elems.begin();
             cit != local_elems.end(); ++cit, ++k)
        {
            active_elems(k) = *cit;
        }
    }
    return;
}// collectActivePatchElements

void
FEDataManager::collectGhostDOFIndices(
    std::vector<unsigned int>& ghost_dofs,
    const blitz::Array<Elem*,1>& active_elems,
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
    for (DofConstraints::const_iterator i = dof_map.constraint_rows_begin();
         i != dof_map.constraint_rows_end(); ++i)
    {
        const unsigned int constrained_dof = i->first;
        if (constrained_dof >= first_local_dof && constrained_dof < end_local_dof)
        {
            const DofConstraintRow& constraint_row = i->second;
            for (DofConstraintRow::const_iterator j = constraint_row.begin();
                 j != constraint_row.end(); ++j)
            {
                const unsigned int constraint_dependency = j->first;
                if (constraint_dependency < first_local_dof || constraint_dependency >= end_local_dof)
                {
                    constraint_dependency_dof_list.push_back(constraint_dependency);
                }
            }
        }
    }

    // Record the local DOFs associated with the active local elements.
    std::set<unsigned int> ghost_dof_set(constraint_dependency_dof_list.begin(), constraint_dependency_dof_list.end());
    for (int e = 0; e < active_elems.size(); ++e)
    {
        const Elem* const elem = active_elems(e);

        // DOFs associated with the element.
        for (unsigned int var_num = 0; var_num < elem->n_vars(sys_num); ++var_num)
        {
            if (elem->n_dofs(sys_num, var_num) > 0)
            {
                const unsigned int dof_index = elem->dof_number(sys_num,var_num,0);
                if (dof_index < first_local_dof || dof_index >= end_local_dof)
                {
                    ghost_dof_set.insert(dof_index);
                }
            }
        }

        // DOFs associated with the nodes of the element.
        for (unsigned int k = 0; k < elem->n_nodes(); ++k)
        {
            Node* node = elem->get_node(k);
            for (unsigned int var_num = 0; var_num < node->n_vars(sys_num); ++var_num)
            {
                if (node->n_dofs(sys_num, var_num) > 0)
                {
                    const unsigned int dof_index = node->dof_number(sys_num,var_num,0);
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
}// collectGhostDOFIndices

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
        TBOX_ERROR("Restart database corresponding to "
                   << d_object_name << " not found in restart file.");
    }

    int ver = db->getInteger("FE_DATA_MANAGER_VERSION");
    if (ver != FE_DATA_MANAGER_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different than class version.");
    }

    d_coarsest_ln = db->getInteger("d_coarsest_ln");
    d_finest_ln   = db->getInteger("d_finest_ln"  );
    return;
}// getFromRestart

void
FEDataManager::do_partition(
    MeshBase& mesh,
    const unsigned int /*n*/)
{
    // Compute the centroids of the elements.
    const unsigned int n_elem = mesh.max_elem_id()+1;
    System& X_system = d_es->get_system(COORDINATES_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    const DofMap& X_dof_map = X_system.get_dof_map();
    NumericVector<double>& X_vec = *X_system.solution;
    NumericVector<double>& X_ghost_vec = *X_system.current_local_solution;
    X_vec.localize(X_ghost_vec);
    X_dof_map.enforce_constraints_exactly(X_system, &X_ghost_vec);

    // Compute the lower and upper bounds of all local elements in the mesh.
    //
    // NOTE: Assuming nodal basis functions.
    std::vector<blitz::TinyVector<double,NDIM> > elem_centroids(n_elem, blitz::TinyVector<double,NDIM>(0.0));
    MeshBase::element_iterator       el_it  = mesh.local_elements_begin();
    const MeshBase::element_iterator el_end = mesh.local_elements_end();
    for ( ; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const unsigned int elem_id = elem->id();
        const unsigned int n_nodes = elem->n_nodes();
        std::vector<unsigned int> dof_indices;
        dof_indices.reserve(NDIM*n_nodes);
        for (unsigned int k = 0; k < n_nodes; ++k)
        {
            Node* node = elem->get_node(k);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(node->n_dofs(X_sys_num,0) > 0);
#endif
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dof_indices.push_back(node->dof_number(X_sys_num,d,0));
            }
        }
        std::vector<double> X_node;
        X_ghost_vec.get(dof_indices, X_node);
        for (unsigned int k = 0; k < n_nodes; ++k)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const double& X = X_node[k*NDIM+d];
                elem_centroids[elem_id][d] += X;
            }
        }
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            elem_centroids[elem_id][d] /= static_cast<double>(n_nodes);
        }
    }

    // Parallel sum elem_centroids so that each process has access to the
    // centroid data for each active element in the mesh.
    std::vector<double> elem_centroids_flattened(NDIM*n_elem);
    for (unsigned int e = 0; e < n_elem; ++e)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            elem_centroids_flattened[e*NDIM+d] = elem_centroids[e][d];
        }
    }
    SAMRAI_MPI::sumReduction(&elem_centroids_flattened[0], elem_centroids_flattened.size());
    for (unsigned int e = 0; e < n_elem; ++e)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            elem_centroids[e][d] = elem_centroids_flattened[e*NDIM+d];
        }
    }

    // Find which patch contains each element centroid.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_number);
    const IntVector<NDIM>& ratio = level->getRatio();
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
#endif
    const Box<NDIM>& domain_box0 = grid_geom->getPhysicalDomain()[0];
    const Box<NDIM>  domain_box  = Box<NDIM>::refine(domain_box0, ratio);
    const Index<NDIM>& lower = domain_box.lower();
    const Index<NDIM>& upper = domain_box.upper();
    const double* const x_lower = grid_geom->getXLower();
    const double* const x_upper = grid_geom->getXUpper();
    const double* const dx0 = grid_geom->getDx();
    double dx[NDIM];
    for (unsigned int d = 0; d < NDIM; ++d) dx[d] = dx0[d]/static_cast<double>(ratio(d));
    Pointer<BoxTree<NDIM> > box_tree = level->getBoxTree();
    const ProcessorMapping& proc_map = level->getProcessorMapping();
    for (el_it = mesh.elements_begin(); el_it != mesh.elements_end(); ++el_it)
    {
        Elem* const elem = *el_it;
        const unsigned int elem_id = elem->id();
        const double* const X = elem_centroids[elem_id].data();
        const Index<NDIM> i = IndexUtilities::getCellIndex(X, x_lower, x_upper, dx, lower, upper);
        Array<int> indices;
        box_tree->findOverlapIndices(indices, Box<NDIM>(i,i));
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(indices.size() == 1);
#endif
        elem->processor_id() = proc_map.getProcessorAssignment(indices[0]);
    }
    return;
}// do_partition

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
