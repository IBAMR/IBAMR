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
#include <ibtk/ibtk_utilities.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/namespaces.h>

// LIBMESH INCLUDES
#include <boundary_info.h>
#include <fe.h>
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
#include <CartesianPatchGeometry.h>
#include <CoarsenAlgorithm.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Pointer<Timer> t_initialize_level_data;
static Pointer<Timer> t_reset_hierarchy_configuration;
static Pointer<Timer> t_apply_gradient_detector;
static Pointer<Timer> t_put_to_database;

// Version of FEDataManager restart file data.
static const int FE_DATA_MANAGER_VERSION = 1;
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
        const IntVector<NDIM> ghost_width = int(floor(0.5*double(stencil_size)))+1;
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

std::pair<const std::vector<std::vector<unsigned int> >&,const std::vector<Elem*>&>
FEDataManager::getActivePatchElements()
{
    return std::make_pair(d_active_patch_elem_map,d_active_patch_elems);
}// getActivePatchElements

void
FEDataManager::reinitElementMappings()
{
    // Reset the mappings between mesh elements and grid patches.
    d_active_patch_elem_map.clear();
    d_active_patch_elems.clear();
    d_active_patch_ghost_dofs.clear();
    collectActivePatchElements(d_active_patch_elem_map, d_active_patch_elems, d_level_number, d_ghost_width);
    clearCachedLEInteractionFEData();
    return;
}// reinitElementMappings

NumericVector<double>*
FEDataManager::getSolutionVector(
    const std::string& system_name)
{
    System& system = d_es->get_system<System>(system_name);
    return system.solution.get();
}// getSolutionVector

NumericVector<double>*
FEDataManager::getGhostedSolutionVector(
    const std::string& system_name)
{
    NumericVector<double>* solution = getSolutionVector(system_name);
    if (d_system_ghost_vec.count(system_name) == 0)
    {
        if (d_active_patch_ghost_dofs.count(system_name) == 0)
        {
            collectGhostDOFIndices(d_active_patch_ghost_dofs[system_name], d_active_patch_elems, system_name);
        }
        AutoPtr<NumericVector<double> > ghosted_solution = NumericVector<double>::build();
        ghosted_solution->init(solution->size(), solution->local_size(), d_active_patch_ghost_dofs[system_name], true, GHOSTED);
        d_system_ghost_vec[system_name] = ghosted_solution.release();
    }
    System& system = d_es->get_system<System>(system_name);
    NumericVector<double>* ghosted_solution = d_system_ghost_vec[system_name];
    solution->localize(*ghosted_solution);
    system.get_dof_map().enforce_constraints_exactly(system, ghosted_solution);
    return ghosted_solution;
}// getGhostedSolutionVector

NumericVector<double>*
FEDataManager::getCoordsVector()
{
    return getSolutionVector(COORDINATES_SYSTEM_NAME);
}// getCoordsVector

NumericVector<double>*
FEDataManager::getGhostedCoordsVector()
{
    return getGhostedSolutionVector(COORDINATES_SYSTEM_NAME);
}// getGhostedCoordsVector

void
FEDataManager::spread(
    const int f_data_idx,
    NumericVector<double>& F,
    NumericVector<double>& X,
    const std::string& system_name,
    const bool close_F,
    const bool close_X)
{
    computeCachedLEInteractionFEData(system_name);
    System& system = d_es->get_system<System>(system_name);
    const unsigned int n_vars = system.n_vars();
    const DofMap& dof_map = system.get_dof_map();
    if (close_F) F.close();
    dof_map.enforce_constraints_exactly(system, &F);
    const blitz::Array<blitz::Array<std::vector<unsigned int>,1>,1>& cached_dof_indices = d_cached_fe_system_data[system_name]->dof_indices();
    const blitz::Array<blitz::Array<double,2>,1>& cached_phi_JxW = d_cached_fe_system_data[system_name]->phi_JxW();

    System& coords_system = d_es->get_system<System>(COORDINATES_SYSTEM_NAME);
    const DofMap& coords_dof_map = coords_system.get_dof_map();
    if (close_X) X.close();
    coords_dof_map.enforce_constraints_exactly(coords_system, &X);
    const blitz::Array<blitz::Array<std::vector<unsigned int>,1>,1>& cached_coords_dof_indices = d_cached_fe_system_data[COORDINATES_SYSTEM_NAME]->dof_indices();
    const blitz::Array<blitz::Array<double,2>,1>& cached_coords_phi = d_cached_fe_system_data[COORDINATES_SYSTEM_NAME]->phi();

    // Loop over the patches to interpolate nodal values on the FE mesh to the
    // element quadrature points, then spread values from the element quadrature
    // points onto the grid to obtain the values on the Eulerian grid.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_number);
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const std::vector<unsigned int>& active_patch_elem_map = d_active_patch_elem_map[local_patch_num];
        const int num_active_patch_elems = active_patch_elem_map.size();
        if (num_active_patch_elems == 0) continue;

        // Loop over the elements and compute the values to be spread and the
        // positions of the quadrature points.
        int qp_offset = 0;
        std::vector<double> F_JxW_qp, X_qp;
        for (int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            const unsigned int e = active_patch_elem_map[e_idx];
            const blitz::Array<std::vector<unsigned int>,1>& dof_indices = cached_dof_indices(e);
            const blitz::Array<double,2>& phi_JxW = cached_phi_JxW(e);
            const int n_qp = phi_JxW.extent(blitz::firstDim);
            const int n_basis = phi_JxW.extent(blitz::secondDim);
            F_JxW_qp.resize(F_JxW_qp.size()+n_vars*n_qp,0.0);
            for (int qp = 0; qp < n_qp; ++qp)
            {
                const int idx = n_vars*(qp+qp_offset);
                for (int k = 0; k < n_basis; ++k)
                {
                    const double& p_JxW = phi_JxW(qp,k);
                    for (unsigned int i = 0; i < n_vars; ++i)
                    {
                        F_JxW_qp[idx+i] += F(dof_indices(i)[k])*p_JxW;
                    }
                }
            }

            const blitz::Array<std::vector<unsigned int>,1>& coords_dof_indices = cached_coords_dof_indices(e);
            const blitz::Array<double,2>& coords_phi = cached_coords_phi(e);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(coords_phi.extent(blitz::firstDim) == n_qp);
#endif
            const int n_coords_basis = coords_phi.extent(blitz::secondDim);
            X_qp.resize(X_qp.size()+NDIM*n_qp,0.0);
            for (int qp = 0; qp < n_qp; ++qp)
            {
                const int idx = NDIM*(qp+qp_offset);
                for (int k = 0; k < n_coords_basis; ++k)
                {
                    const double& p = coords_phi(qp,k);
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        X_qp[idx+i] += X(coords_dof_indices(i)[k])*p;
                    }
                }
            }

            qp_offset += n_qp;
        }

        if (qp_offset == 0) continue;

        // Spread the values to the Cartesian grid patch.
        //
        // NOTE: Values are spread only from those quadrature points that are
        // within the ghost cell width of the patch interior.
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const int stencil_size = LEInteractor::getStencilSize(d_spread_weighting_fcn);
        const hier::IntVector<NDIM> ghost_width = int(ceil(0.5*double(stencil_size)));
        const Box<NDIM> spread_box = Box<NDIM>::grow(patch->getBox(), ghost_width);
        Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_data_idx);
        Pointer<CellData<NDIM,double> > f_cc_data = f_data;
        Pointer<SideData<NDIM,double> > f_sc_data = f_data;
        const bool is_cc_data = !f_cc_data.isNull();
        const bool is_sc_data = !f_sc_data.isNull();
        if (is_cc_data) LEInteractor::spread(f_cc_data, F_JxW_qp, n_vars, X_qp, NDIM, patch, spread_box, d_spread_weighting_fcn);
        if (is_sc_data) LEInteractor::spread(f_sc_data, F_JxW_qp, n_vars, X_qp, NDIM, patch, spread_box, d_spread_weighting_fcn);
    }
    return;
}// spread

void
FEDataManager::interp(
    const int f_data_idx,
    NumericVector<double>& F,
    NumericVector<double>& X,
    const std::string& system_name,
    std::vector<Pointer<RefineSchedule<NDIM> > > f_refine_scheds,
    const double fill_data_time,
    const bool close_X)
{
    computeCachedLEInteractionFEData(system_name);
    ExplicitSystem& system = d_es->get_system<ExplicitSystem>(system_name);
    const unsigned int n_vars = system.n_vars();
    const DofMap& dof_map = system.get_dof_map();
    AutoPtr<NumericVector<double> > rhs = F.clone();
    rhs->zero();
    const blitz::Array<blitz::Array<std::vector<unsigned int>,1>,1>& cached_dof_indices = d_cached_fe_system_data[system_name]->dof_indices();
    const blitz::Array<blitz::Array<double,2>,1>& cached_phi_JxW = d_cached_fe_system_data[system_name]->phi_JxW();

    System& coords_system = d_es->get_system<System>(COORDINATES_SYSTEM_NAME);
    const DofMap& coords_dof_map = coords_system.get_dof_map();
    if (close_X) X.close();
    coords_dof_map.enforce_constraints_exactly(coords_system, &X);
    const blitz::Array<blitz::Array<std::vector<unsigned int>,1>,1>& cached_coords_dof_indices = d_cached_fe_system_data[COORDINATES_SYSTEM_NAME]->dof_indices();
    const blitz::Array<blitz::Array<double,2>,1>& cached_coords_phi = d_cached_fe_system_data[COORDINATES_SYSTEM_NAME]->phi();

    for (unsigned int k = 0; k < f_refine_scheds.size(); ++k)
    {
        if (!f_refine_scheds[k].isNull()) f_refine_scheds[k]->fillData(fill_data_time);
    }

    // Loop over the patches to interpolate values to the element quadrature
    // points from the grid, then integrate the values at the element quadrature
    // points to obtain the nodal values on the Lagrangian mesh.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_number);
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const std::vector<unsigned int>& active_patch_elem_map = d_active_patch_elem_map[local_patch_num];
        const int num_active_patch_elems = active_patch_elem_map.size();
        if (num_active_patch_elems == 0) continue;

        // Loop over the elements and compute the positions of the quadrature points.
        int qp_offset = 0;
        std::vector<double> F_qp, X_qp;
        for (int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            const unsigned int e = active_patch_elem_map[e_idx];
            const blitz::Array<std::vector<unsigned int>,1>& coords_dof_indices = cached_coords_dof_indices(e);
            const blitz::Array<double,2>& coords_phi = cached_coords_phi(e);
            const int n_qp = coords_phi.extent(blitz::firstDim);
            const int n_basis = coords_phi.extent(blitz::secondDim);
            F_qp.resize(F_qp.size()+n_vars*n_qp,0.0);
            X_qp.resize(X_qp.size()+NDIM  *n_qp,0.0);
            for (int qp = 0; qp < n_qp; ++qp)
            {
                const int idx = NDIM*(qp+qp_offset);
                for (int k = 0; k < n_basis; ++k)
                {
                    const double& p = coords_phi(qp,k);
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        X_qp[idx+i] += X(coords_dof_indices(i)[k])*p;
                    }
                }
            }
            qp_offset += n_qp;
        }

        if (qp_offset == 0) continue;

        // Interpolate the values from the Cartesian grid patch to the
        // quadrature points.
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
        std::vector<DenseVector<double> > rhs_e(n_vars);
        for (int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            const unsigned int e = active_patch_elem_map[e_idx];
            blitz::Array<std::vector<unsigned int>,1> dof_indices = cached_dof_indices(e);
            for (unsigned int i = 0; i < n_vars; ++i)
            {
                rhs_e[i].resize(dof_indices(i).size());
            }
            const blitz::Array<double,2>& phi_JxW = cached_phi_JxW(e);
            const int n_qp = phi_JxW.extent(blitz::firstDim);
            const int n_basis = phi_JxW.extent(blitz::secondDim);
            for (int qp = 0; qp < n_qp; ++qp)
            {
                const int idx = n_vars*(qp+qp_offset);
                for (int k = 0; k < n_basis; ++k)
                {
                    const double& p_JxW = phi_JxW(qp,k);
                    for (unsigned int i = 0; i < n_vars; ++i)
                    {
                        rhs_e[i](k) += F_qp[idx+i]*p_JxW;
                    }
                }
            }
            for (unsigned int i = 0; i < n_vars; ++i)
            {
                dof_map.constrain_element_vector(rhs_e[i], dof_indices(i));
                rhs->add_vector(rhs_e[i], dof_indices(i));
            }
            qp_offset += n_qp;
        }
    }

    // Solve for the nodal values.
    computeL2Projection(F, *rhs, system_name, d_interp_uses_consistent_mass_matrix);
    return;
}// interp

std::pair<LinearSolver<double>*,SparseMatrix<double>*>
FEDataManager::getL2ProjectionSolver(
    const std::string& system_name,
    const bool consistent_mass_matrix,
    const QuadratureType quad_type,
    const Order quad_order)
{
    if ((d_L2_proj_solver.count(system_name) == 0 || d_L2_proj_matrix.count(system_name) == 0) ||
        (d_L2_proj_consistent_mass_matrix[system_name] != consistent_mass_matrix) ||
        (d_L2_proj_quad_type[system_name] != quad_type) || (d_L2_proj_quad_order[system_name] != quad_order))
    {
        const MeshBase& mesh = d_es->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();
        AutoPtr<QBase> qrule = QBase::build(quad_type, dim, quad_order);

        System& system = d_es->get_system<System>(system_name);
        DofMap& dof_map = system.get_dof_map();
        // NOTE: Sparsity patterns are not automatically computed for all system
        // types; if it has not been computed for this system, we must compute
        // it now.
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

        SparseMatrix<double>* M = SparseMatrix<double>::build().release();
        M->attach_dof_map(dof_map);
        M->init();

        DenseMatrix<double> M_e;

        // Loop over the mesh to construct the system matrix.
        MeshBase::const_element_iterator       el_it  = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for ( ; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            fe->reinit(elem);
            for (unsigned int var_num = 0; var_num < dof_map.n_variables(); ++var_num)
            {
                dof_map.dof_indices(elem, dof_indices, var_num);
                M_e.resize(dof_indices.size(), dof_indices.size());
                for (unsigned int i = 0; i < phi.size(); ++i)
                {
                    for (unsigned int j = 0; j < phi.size(); ++j)
                    {
                        for (unsigned int qp = 0; qp < qrule->n_points(); ++qp)
                        {
                            if (consistent_mass_matrix)
                            {
                                M_e(i,j) += (phi[i][qp]*phi[j][qp])*JxW[qp];
                            }
                            else
                            {
                                M_e(j,j) += (phi[i][qp]*phi[j][qp])*JxW[qp];
                            }
                        }
                    }
                }
                dof_map.constrain_element_matrix(M_e, dof_indices);
                M->add_matrix(M_e, dof_indices);
            }
        }

        // Assemble the matrix.
        M->close();

        // Setup the solver.
        solver->same_preconditioner = true;

        // Store the solver, mass matrix, and configuration options.
        d_L2_proj_solver[system_name] = solver;
        d_L2_proj_matrix[system_name] = M;
        d_L2_proj_consistent_mass_matrix[system_name] = consistent_mass_matrix;
        d_L2_proj_quad_type[system_name] = quad_type;
        d_L2_proj_quad_order[system_name] = quad_order;
    }
    return std::make_pair(d_L2_proj_solver[system_name], d_L2_proj_matrix[system_name]);
}// getL2ProjectionSolver

NumericVector<double>*
FEDataManager::getDiagonalL2MassMatrix(
    const std::string& system_name,
    const QuadratureType quad_type,
    const Order quad_order)
{
    if (d_L2_proj_matrix_diag.count(system_name) == 0)
    {
        const MeshBase& mesh = d_es->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();
        AutoPtr<QBase> qrule = QBase::build(quad_type, dim, quad_order);

        System& system = d_es->get_system<System>(system_name);
        DofMap& dof_map = system.get_dof_map();
        // NOTE: Sparsity patterns are not automatically computed for all system
        // types; if it has not been computed for this system, we must compute
        // it now.
        if (dof_map.get_n_nz().size() != dof_map.n_local_dofs())
        {
            dof_map.compute_sparsity(mesh);
        }
        std::vector<unsigned int> dof_indices;
        AutoPtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<double>& JxW = fe->get_JxW();
        const std::vector<std::vector<double> >& phi = fe->get_phi();

        NumericVector<double>* M = system.solution->clone().release();
        DenseVector<double> M_e;

        // Loop over the mesh to construct the system matrix.
        MeshBase::const_element_iterator       el_it  = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for ( ; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            fe->reinit(elem);
            for (unsigned int var_num = 0; var_num < dof_map.n_variables(); ++var_num)
            {
                dof_map.dof_indices(elem, dof_indices, var_num);
                M_e.resize(dof_indices.size());
                for (unsigned int i = 0; i < phi.size(); ++i)
                {
                    for (unsigned int j = 0; j < phi.size(); ++j)
                    {
                        for (unsigned int qp = 0; qp < qrule->n_points(); ++qp)
                        {
                            M_e(j) += (phi[i][qp]*phi[j][qp])*JxW[qp];
                        }
                    }
                }
                dof_map.constrain_element_vector(M_e, dof_indices);
                M->add_vector(M_e, dof_indices);
            }
        }

        // Assemble the vector representation of the diagonal matrix.
        M->close();

        // Store the diagonal mass matrix.
        d_L2_proj_matrix_diag[system_name] = M;
    }
    return d_L2_proj_matrix_diag[system_name];
}// getDiagonalL2MassMatrix

bool
FEDataManager::computeL2Projection(
    NumericVector<double>& U,
    NumericVector<double>& F,
    const std::string& system_name,
    const bool consistent_mass_matrix,
    const QuadratureType quad_type,
    const Order quad_order,
    const double tol,
    const unsigned int max_its)
{
    int ierr;
    bool converged = false;

    F.close();
    const System& system = d_es->get_system<System>(system_name);
    const DofMap& dof_map = system.get_dof_map();
    dof_map.enforce_constraints_exactly(system, &F);
    if (consistent_mass_matrix)
    {
        std::pair<libMesh::LinearSolver<double>*,SparseMatrix<double>*> proj_solver_components =
            getL2ProjectionSolver(system_name, consistent_mass_matrix, quad_type, quad_order);
        PetscLinearSolver<double>* solver = dynamic_cast<PetscLinearSolver<double>*>(proj_solver_components.first);
        PetscMatrix<double>* M = dynamic_cast<PetscMatrix<double>*>(proj_solver_components.second);
        solver->solve(*M, *M, U, F, tol, max_its);
        KSPConvergedReason reason;
        ierr = KSPGetConvergedReason(solver->ksp(), &reason); IBTK_CHKERRQ(ierr);
        converged = reason > 0;
    }
    else
    {
        PetscVector<double>* M_diag = dynamic_cast<PetscVector<double>*>(getDiagonalL2MassMatrix(system_name,quad_type,quad_order));
        Vec M_diag_vec = M_diag->vec();
        Vec U_vec = dynamic_cast<PetscVector<double>*>(&U)->vec();
        Vec F_vec = dynamic_cast<PetscVector<double>*>(&F)->vec();
        ierr = VecPointwiseDivide(U_vec, F_vec, M_diag_vec); IBTK_CHKERRQ(ierr);
        converged = true;
    }
    dof_map.enforce_constraints_exactly(system, &U);
    return converged;
}// computeL2Projection

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
    const bool can_be_refined,
    const bool initial_time,
    const Pointer<BasePatchLevel<NDIM> > old_level,
    const bool allocate_data)
{
    t_initialize_level_data->start();

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

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
FEDataManager::resetHierarchyConfiguration(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_ln,
    const int finest_ln)
{
    t_reset_hierarchy_configuration->start();

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

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
FEDataManager::applyGradientDetector(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    if (level_number >= d_level_number) return;

    t_apply_gradient_detector->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    if (initial_time)
    {
        // Compute element bounding boxes and determine the active elements
        // associated with the prescribed patch level.
        std::vector<std::vector<unsigned int> > active_level_elem_map;
        std::vector<Elem*> active_level_elems;
        const IntVector<NDIM> ghost_width = 1;
        collectActivePatchElements(active_level_elem_map, active_level_elems, level_number, ghost_width);
        std::vector<unsigned int> coords_ghost_dofs;
        collectGhostDOFIndices(coords_ghost_dofs, active_level_elems, COORDINATES_SYSTEM_NAME);

        const MeshBase& mesh = d_es->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        System& coords_system = d_es->get_system<System>(COORDINATES_SYSTEM_NAME);
        const DofMap& coords_dof_map = coords_system.get_dof_map();
        blitz::Array<std::vector<unsigned int>,1> coords_dof_indices(dim);
        AutoPtr<FEBase> coords_fe(FEBase::build(dim, coords_dof_map.variable_type(0)));
        coords_fe->attach_quadrature_rule(d_qrule);
        const std::vector<std::vector<double> >& coords_phi = coords_fe->get_phi();
        NumericVector<double>* X = getCoordsVector();
        AutoPtr<NumericVector<double> > X_ghost = NumericVector<double>::build();
        X_ghost->init(X->size(), X->local_size(), coords_ghost_dofs, true, GHOSTED);
        X->localize(*X_ghost);
        coords_dof_map.enforce_constraints_exactly(coords_system, X_ghost.get());

        // Tag cells for refinement whenever they contain element quadrature
        // points.
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            const Pointer<Patch<NDIM> > patch = level->getPatch(p());

            const Box<NDIM>& patch_box = patch->getBox();
            const CellIndex<NDIM>& patch_lower = patch_box.lower();
            const CellIndex<NDIM>& patch_upper = patch_box.upper();

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();
            const double* const patch_dx = patch_geom->getDx();

            Pointer<CellData<NDIM,int> > tag_data = patch->getPatchData(tag_index);

            const std::vector<unsigned int>& active_patch_elem_map = active_level_elem_map[local_patch_num];
            const int num_active_patch_elems = active_patch_elem_map.size();
            if (num_active_patch_elems == 0) continue;
            for (int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                const unsigned int e = active_patch_elem_map[e_idx];
                const Elem* const elem = active_level_elems[e];
                coords_fe->reinit(elem);
                for (unsigned int d = 0; d < dim; ++d)
                {
                    coords_dof_map.dof_indices(elem, coords_dof_indices(d), d);
                }
                for (unsigned int qp = 0; qp < d_qrule->n_points(); ++qp)
                {
                    std::vector<double> X_qp(NDIM,0.0);
                    for (unsigned int k = 0; k < coords_phi.size(); ++k)
                    {
                        for (unsigned int i = 0; i < NDIM; ++i)
                        {
                            X_qp[i] += (*X_ghost)(coords_dof_indices(i)[k])*coords_phi[k][qp];
                        }
                    }
                    const Index<NDIM> i = IndexUtilities::getCellIndex(X_qp, patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper);
                    tag_data->fill(1,Box<NDIM>::Box(i-Index<NDIM>(1),i+Index<NDIM>(1)));
                }
            }
        }
    }
    else if (level_number+1 == d_level_number && level_number < d_hierarchy->getFinestLevelNumber())
    {
        Pointer<PatchLevel<NDIM> > finer_level = d_hierarchy->getPatchLevel(level_number+1);
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);

        // Update the node count data and coarsen it from the finer level.
        if (!      level->checkAllocated(d_qp_count_idx))       level->allocatePatchData(d_qp_count_idx);
        if (!finer_level->checkAllocated(d_qp_count_idx)) finer_level->allocatePatchData(d_qp_count_idx);
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
                if (!MathUtilities<double>::equalEps((*qp_count_data)(i),0.0))
                {
                    (*tag_data)(i) = 1;
                }
            }
        }

        // Deallocate the node count data.
        level->deallocatePatchData(d_qp_count_idx);
        finer_level->deallocatePatchData(d_qp_count_idx);
    }
    t_apply_gradient_detector->stop();
    return;
}// applyGradientDetector

void
FEDataManager::putToDatabase(
    Pointer<Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("FE_DATA_MANAGER_VERSION", FE_DATA_MANAGER_VERSION);

    db->putInteger("d_coarsest_ln", d_coarsest_ln);
    db->putInteger("d_finest_ln"  , d_finest_ln  );

    t_put_to_database->stop();
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
      d_active_patch_elems(),
      d_active_patch_ghost_dofs(),
      d_L2_proj_solver(),
      d_L2_proj_matrix(),
      d_L2_proj_matrix_diag(),
      d_L2_proj_consistent_mass_matrix(),
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

    // Register the node count variable with the VariableDatabase.
    d_qp_count_var = new CellVariable<NDIM,double>(d_object_name+"::qp_count");
    d_qp_count_idx = var_db->registerVariableAndContext(d_qp_count_var, d_context, 0);

    // Setup Timers.
    IBTK_DO_ONCE(
        t_initialize_level_data = TimerManager::getManager()->getTimer("IBTK::FEDataManager::initializeLevelData()");
        t_reset_hierarchy_configuration = TimerManager::getManager()->getTimer("IBTK::FEDataManager::resetHierarchyConfiguration()");
        t_apply_gradient_detector = TimerManager::getManager()->getTimer("IBTK::FEDataManager::applyGradientDetector()");
        t_put_to_database = TimerManager::getManager()->getTimer("IBTK::FEDataManager::putToDatabase()");
        LEInteractor::initializeTimers();
                 );
    return;
}// FEDataManager

FEDataManager::~FEDataManager()
{
    for (std::map<std::string,NumericVector<double>*>::iterator it = d_system_ghost_vec.begin();
         it != d_system_ghost_vec.end(); ++it)
    {
        delete (*it).second;
    }
    for (std::map<std::string,LinearSolver<double>*>::iterator it = d_L2_proj_solver.begin();
         it != d_L2_proj_solver.end(); ++it)
    {
        delete (*it).second;
    }
    for (std::map<std::string,SparseMatrix<double>*>::iterator it = d_L2_proj_matrix.begin();
         it != d_L2_proj_matrix.end(); ++it)
    {
        delete (*it).second;
    }
    for (std::map<std::string,NumericVector<double>*>::iterator it = d_L2_proj_matrix_diag.begin();
         it != d_L2_proj_matrix_diag.end(); ++it)
    {
        delete (*it).second;
    }
    return;
}// ~FEDataManager

/////////////////////////////// PRIVATE //////////////////////////////////////

void
FEDataManager::updateQuadPointCountData(
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1) ? d_coarsest_ln : coarsest_ln_in;
    const int finest_ln = (finest_ln_in == -1) ? d_finest_ln : finest_ln_in;

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(coarsest_ln >= d_coarsest_ln && coarsest_ln <= d_finest_ln);
    TBOX_ASSERT(finest_ln   >= d_coarsest_ln && finest_ln   <= d_finest_ln);
#endif

    const MeshBase& mesh = d_es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    System& coords_system = d_es->get_system<System>(COORDINATES_SYSTEM_NAME);
    const DofMap& coords_dof_map = coords_system.get_dof_map();
    blitz::Array<std::vector<unsigned int>,1> coords_dof_indices(dim);
    AutoPtr<FEBase> coords_fe(FEBase::build(dim, coords_dof_map.variable_type(0)));
    coords_fe->attach_quadrature_rule(d_qrule);
    const std::vector<std::vector<double> >& coords_phi = coords_fe->get_phi();
    NumericVector<double>* X = getCoordsVector();
    NumericVector<double>* X_ghost = getGhostedCoordsVector();
    X->localize(*X_ghost);
    coords_dof_map.enforce_constraints_exactly(coords_system, X_ghost);

    // Set the node count data on the specified range of levels of the
    // hierarchy.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_qp_count_idx)) level->allocatePatchData(d_qp_count_idx);
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM,double> > qp_count_data = patch->getPatchData(d_qp_count_idx);
            qp_count_data->fillAll(0.0);
            if (ln == d_level_number)
            {
                const Box<NDIM>& patch_box = patch->getBox();
                const CellIndex<NDIM>& patch_lower = patch_box.lower();
                const CellIndex<NDIM>& patch_upper = patch_box.upper();

                const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* const patch_x_lower = patch_geom->getXLower();
                const double* const patch_x_upper = patch_geom->getXUpper();
                const double* const patch_dx = patch_geom->getDx();

                // Keep track of the number of quadrature points in each Cartesian grid cell.
                const std::vector<unsigned int>& active_patch_elem_map = d_active_patch_elem_map[local_patch_num];
                const int num_active_patch_elems = active_patch_elem_map.size();
                if (num_active_patch_elems == 0) continue;
                for (int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
                {
                    const unsigned int e = active_patch_elem_map[e_idx];
                    const Elem* const elem = d_active_patch_elems[e];
                    coords_fe->reinit(elem);
                    for (unsigned int d = 0; d < dim; ++d)
                    {
                        coords_dof_map.dof_indices(elem, coords_dof_indices(d), d);
                    }
                    for (unsigned int qp = 0; qp < d_qrule->n_points(); ++qp)
                    {
                        std::vector<double> X_qp(NDIM,0.0);
                        for (unsigned int k = 0; k < coords_phi.size(); ++k)
                        {
                            for (unsigned int i = 0; i < NDIM; ++i)
                            {
                                X_qp[i] += (*X_ghost)(coords_dof_indices(i)[k])*coords_phi[k][qp];
                            }
                        }
                        const Index<NDIM> i = IndexUtilities::getCellIndex(X_qp, patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper);
                        if (patch_box.contains(i)) (*qp_count_data)(i) += 1.0;
                    }
                }
            }
        }
    }
    return;
}// updateQuadPointCountData

void
FEDataManager::computeActiveElementBoundingBoxes(
    std::vector<double>& elem_bounds)
{
    elem_bounds.clear();

    const MeshBase& mesh = d_es->get_mesh();
    const unsigned int n_elem = mesh.n_elem();
    System& coords_system = d_es->get_system<System>(COORDINATES_SYSTEM_NAME);
    const DofMap& coords_dof_map = coords_system.get_dof_map();
    NumericVector<double>& X = *coords_system.solution;
    NumericVector<double>& X_ghost = *coords_system.current_local_solution;
    X.localize(X_ghost);
    coords_dof_map.enforce_constraints_exactly(coords_system, &X_ghost);
    const unsigned int coords_sys_num = coords_system.number();

    // Compute the lower and upper bounds of all active local elements in the
    // mesh.  Assumes nodal basis functions.
    elem_bounds.resize(2*n_elem*NDIM,0.0);
    MeshBase::const_element_iterator       el_it  = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for ( ; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const unsigned int elem_id = elem->id();
        double* const elem_lower_bound = &elem_bounds[2*NDIM*elem_id     ];
        double* const elem_upper_bound = &elem_bounds[2*NDIM*elem_id+NDIM];
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            elem_lower_bound[d] =  0.5*std::numeric_limits<double>::max();
            elem_upper_bound[d] = -0.5*std::numeric_limits<double>::max();
        }

        for (unsigned int k = 0; k < elem->n_nodes(); ++k)
        {
            Node* node = elem->get_node(k);
            if (node->n_dofs(coords_sys_num,0) > 0)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    const int dof_index = node->dof_number(coords_sys_num,d,0);
                    elem_lower_bound[d] = std::min(elem_lower_bound[d], X_ghost(dof_index));
                    elem_upper_bound[d] = std::max(elem_upper_bound[d], X_ghost(dof_index));
                }
            }
        }
    }

    // Parallel sum elem_lower_bound and elem_upper_bound so that each process
    // has access to the bounding box data for each active element in the mesh.
    SAMRAI_MPI::sumReduction(&elem_bounds[0], elem_bounds.size());
    return;
}// computeActiveElementBoundingBoxes

void
FEDataManager::collectActivePatchElements(
    std::vector<std::vector<unsigned int> >& active_patch_elem_map,
    std::vector<Elem*>& active_patch_elems,
    const int level_number,
    const IntVector<NDIM>& ghost_width)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif

    active_patch_elem_map.clear();
    active_patch_elems.clear();

    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
    const int num_local_patches = level->getProcessorMapping().getNumberOfLocalIndices();

    std::vector<std::vector<Elem*> > active_patch_elem_vec;
    collectActivePatchElements_helper(active_patch_elem_vec, level_number, ghost_width);
    std::set<Elem*> active_patch_elems_set;
    for (int local_patch_num = 0; local_patch_num < num_local_patches; ++local_patch_num)
    {
        active_patch_elems_set.insert(active_patch_elem_vec[local_patch_num].begin(),active_patch_elem_vec[local_patch_num].end());
    }
    active_patch_elems.insert(active_patch_elems.end(), active_patch_elems_set.begin(), active_patch_elems_set.end());

    active_patch_elem_map.resize(num_local_patches);
    const std::vector<Elem*>::iterator el_begin = active_patch_elems.begin();
    const std::vector<Elem*>::iterator el_end   = active_patch_elems.end();
    for (int local_patch_num = 0; local_patch_num < num_local_patches; ++local_patch_num)
    {
        const int num_elems = active_patch_elem_vec[local_patch_num].size();
        active_patch_elem_map[local_patch_num].resize(num_elems);
        for (int k = 0; k < num_elems; ++k)
        {
            const Elem* const elem = active_patch_elem_vec[local_patch_num][k];
            active_patch_elem_map[local_patch_num][k] = std::distance(el_begin,std::equal_range(el_begin,el_end,elem).first);
        }
    }
    return;
}// collectActivePatchElements

void
FEDataManager::collectActivePatchElements_helper(
    std::vector<std::vector<Elem*> >& active_patch_elems,
    const int level_number,
    const IntVector<NDIM>& ghost_width)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif

    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
    active_patch_elems.clear();
    const int num_local_patches = level->getProcessorMapping().getNumberOfLocalIndices();
    active_patch_elems.resize(num_local_patches);

    // We initially associate an element with a Cartesian grid patch if the
    // element's bounding box intersects the patch interior grown by the
    // specified ghost cell width.
    std::vector<double> elem_bounds;
    computeActiveElementBoundingBoxes(elem_bounds);
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

        std::vector<double> xLower(pgeom->getXLower(),pgeom->getXLower()+NDIM);
        std::vector<double> xUpper(pgeom->getXUpper(),pgeom->getXUpper()+NDIM);
        const double* const dx = pgeom->getDx();
        for (int d = 0; d < NDIM; ++d)
        {
            xLower[d] -= dx[d]*ghost_width[d];
            xUpper[d] += dx[d]*ghost_width[d];
        }

        const MeshBase& mesh = d_es->get_mesh();
        MeshBase::const_element_iterator       el_it  = mesh.active_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_elements_end();
        for ( ; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            const unsigned int elem_id = elem->id();
            const double* const elem_lower_bound = &elem_bounds[2*NDIM*elem_id     ];
            const double* const elem_upper_bound = &elem_bounds[2*NDIM*elem_id+NDIM];
            bool in_patch = true;
            for (int d = 0; d < NDIM && in_patch; ++d)
            {
                in_patch = in_patch && ((elem_upper_bound[d] >= xLower[d] && elem_upper_bound[d] <= xUpper[d]) ||
                                        (elem_lower_bound[d] >= xLower[d] && elem_lower_bound[d] <= xUpper[d]));
            }
            if (in_patch)
            {
                active_patch_elems[local_patch_num].push_back(elem);
            }
        }
    }

    // Recursively add/remove elements from the active sets that were generated
    // via the bounding box method.
    bool done = false;
    while (!done)
    {
        std::vector<std::set<Elem*> > neighbor_elems(num_local_patches);
        for (int local_patch_num = 0; local_patch_num < num_local_patches; ++local_patch_num)
        {
            std::vector<Elem*>& patch_elems = active_patch_elems[local_patch_num];
            if (patch_elems.empty()) continue;

            // Collect all active elements that neighbor the set of elements
            // currently associated with the patch.
            std::set<Elem*> patch_elems_set(patch_elems.begin(), patch_elems.end());
            for (std::vector<Elem*>::const_iterator cit = patch_elems.begin(); cit != patch_elems.end(); ++cit)
            {
                const Elem* const elem = *cit;
                for (unsigned n = 0; n < elem->n_neighbors(); ++n)
                {
                    Elem* const nghbr_elem = elem->neighbor(n);
                    if (nghbr_elem != NULL && nghbr_elem->active() &&
                        patch_elems_set.count(nghbr_elem) == 0)
                    {
                        neighbor_elems[local_patch_num].insert(nghbr_elem);
                    }
                }
            }

            // Merge the set of neighboring elements with the set of elements
            // currently associated with the patch.
            patch_elems.insert(patch_elems.end(),
                               neighbor_elems[local_patch_num].begin(),
                               neighbor_elems[local_patch_num].end());
        }

        // Setup an appropriately ghosted temporary coordinates vector.
        const MeshBase& mesh = d_es->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        System& coords_system = d_es->get_system<System>(COORDINATES_SYSTEM_NAME);
        const DofMap& coords_dof_map = coords_system.get_dof_map();
        std::vector<std::vector<unsigned int> > coords_dof_indices(dim);
        AutoPtr<FEBase> coords_fe(FEBase::build(dim, coords_dof_map.variable_type(0)));
        coords_fe->attach_quadrature_rule(d_qrule);
        const std::vector<std::vector<double> >& coords_phi = coords_fe->get_phi();
        NumericVector<double>* X = getCoordsVector();
        AutoPtr<NumericVector<double> > X_ghost = NumericVector<double>::build();
        std::vector<unsigned int> coords_ghost_dofs;
        collectGhostDOFIndices(coords_ghost_dofs, active_patch_elems, COORDINATES_SYSTEM_NAME);
        X_ghost->init(X->size(), X->local_size(), coords_ghost_dofs, true, GHOSTED);
        X->localize(*X_ghost);
        coords_dof_map.enforce_constraints_exactly(coords_system, X_ghost.get());

        // Keep only those elements that have a quadrature point on the local
        // patch.  We also keep track of whether we added any new elements to
        // the list of active elements; if not, then we have found all of the
        // elements associated with the local patches.
        std::vector<std::vector<Elem*> > new_active_patch_elems(num_local_patches);
        bool inserted_ngbhr_elem = false;
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            const std::vector<Elem*>& patch_elems = active_patch_elems[local_patch_num];
            if (patch_elems.empty()) continue;

            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Box<NDIM>  ghost_box = Box<NDIM>::grow(patch_box, ghost_width);
            const CellIndex<NDIM>& patch_lower = patch_box.lower();
            const CellIndex<NDIM>& patch_upper = patch_box.upper();

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();
            const double* const patch_dx = patch_geom->getDx();

            std::vector<Elem*>::const_iterator       el_it  = patch_elems.begin();
            const std::vector<Elem*>::const_iterator el_end = patch_elems.end();
            for ( ; el_it != el_end; ++el_it)
            {
                Elem* const elem = *el_it;
                coords_fe->reinit(elem);
                for (unsigned int d = 0; d < dim; ++d)
                {
                    coords_dof_map.dof_indices(elem, coords_dof_indices[d], d);
                }

                bool found_qp = false;
                for (unsigned int qp = 0; qp < d_qrule->n_points() && !found_qp; ++qp)
                {
                    std::vector<double> X_qp(NDIM,0.0);
                    for (unsigned int k = 0; k < coords_phi.size(); ++k)
                    {
                        for (unsigned int i = 0; i < NDIM; ++i)
                        {
                            X_qp[i] += (*X_ghost)(coords_dof_indices[i][k])*coords_phi[k][qp];
                        }
                    }
                    const Index<NDIM> i = IndexUtilities::getCellIndex(X_qp, patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper);
                    if (ghost_box.contains(i))
                    {
                        new_active_patch_elems[local_patch_num].push_back(elem);
                        if (neighbor_elems[local_patch_num].count(elem) > 0) inserted_ngbhr_elem = true;
                        found_qp = true;
                    }
                }
            }
        }
        active_patch_elems = new_active_patch_elems;

        // Check to see if any processors have added a neighbor element; if so,
        // we are not yet done.
        done = SAMRAI_MPI::sumReduction(inserted_ngbhr_elem ? 1 : 0) == 0;
    }

    // Sort the element pointers so that they are in ascending order, in an
    // attempt to improve cache performance.  This probably doesn't make any
    // real difference in practice.
    for (unsigned k = 0; k < active_patch_elems.size(); ++k)
    {
        std::sort(active_patch_elems[k].begin(), active_patch_elems[k].end());
    }
    return;
}// collectActivePatchElements_helper

void
FEDataManager::collectGhostDOFIndices(
    std::vector<unsigned int>& ghost_dofs,
    std::vector<Elem*>& active_patch_elems,
    const std::string& system_name)
{
    ghost_dofs.clear();

    System& system = d_es->get_system<System>(system_name);
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
    for (std::vector<Elem*>::const_iterator cit = active_patch_elems.begin(); cit != active_patch_elems.end(); ++cit)
    {
        const Elem* const elem = *cit;

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
    ghost_dofs.insert(ghost_dofs.end(), ghost_dof_set.begin(), ghost_dof_set.end());
    return;
}// collectGhostDOFIndices

void
FEDataManager::clearCachedLEInteractionFEData()
{
    d_cached_fe_system_data.clear();
    return;
}// clearCachedLEInteractionFEData

void
FEDataManager::computeCachedLEInteractionFEData(
    const std::string& system_name)
{
    if (d_cached_fe_system_data[system_name] != NULL) return;
    TBOX_ASSERT(false);
    return;
}// computeCachedLEInteractionFEData

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

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
