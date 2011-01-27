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

const std::map<int,std::set<Elem*> >&
FEDataManager::getActivePatchElements()
{
    return d_active_patch_elems;
}// getActivePatchElements

void
FEDataManager::reinitElementMappings()
{
    // Reset the mappings between mesh elements and grid patches.
    d_active_patch_elems.clear();
    d_active_patch_ghost_dofs.clear();
    collectActivePatchElements(d_active_patch_elems, d_level_number, d_ghost_width);

    // Delete cached vectors.
    for (std::map<std::string,NumericVector<double>*>::iterator it = d_system_ghost_vec.begin();
         it != d_system_ghost_vec.end(); ++it)
    {
        delete (*it).second;
    }
    d_system_ghost_vec.clear();

    for (std::map<std::string,LinearSolver<double>*>::iterator it = d_L2_projection_solvers.begin();
         it != d_L2_projection_solvers.end(); ++it)
    {
        delete (*it).second;
    }
    d_L2_projection_solvers.clear();

    for (std::map<std::string,SparseMatrix<double>*>::iterator it = d_L2_mass_matrices.begin();
         it != d_L2_mass_matrices.end(); ++it)
    {
        delete (*it).second;
    }
    d_L2_mass_matrices.clear();
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
        ghosted_solution->init(solution->size(), solution->local_size(), d_active_patch_ghost_dofs[system_name], GHOSTED);
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
    const MeshBase& mesh = d_es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    System& system = d_es->get_system<System>(system_name);
    const unsigned int n_vars = system.n_vars();
    const DofMap& dof_map = system.get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(n_vars);
    AutoPtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    fe->attach_quadrature_rule(d_qrule);
    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<std::vector<double> >& phi = fe->get_phi();
    if (close_F) F.close();
    dof_map.enforce_constraints_exactly(system, &F);

    System& coords_system = d_es->get_system<System>(COORDINATES_SYSTEM_NAME);
    const DofMap& coords_dof_map = coords_system.get_dof_map();
    std::vector<std::vector<unsigned int> > coords_dof_indices(dim);
    AutoPtr<FEBase> coords_fe(FEBase::build(dim, coords_dof_map.variable_type(0)));
    coords_fe->attach_quadrature_rule(d_qrule);
    const std::vector<std::vector<double> >& coords_phi = coords_fe->get_phi();
    if (close_X) X.close();
    coords_dof_map.enforce_constraints_exactly(coords_system, &X);

    // Loop over the patches to interpolate nodal values on the FE mesh to the
    // element quadrature points, then spread values from the element quadrature
    // points onto the grid to obtain the values on the Eulerian grid.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const int patch_num = patch->getPatchNumber();

        // The relevant collection of elements.
        std::map<int,std::set<Elem*> >::const_iterator pos = d_active_patch_elems.find(patch_num);
        if (pos == d_active_patch_elems.end()) continue;
        std::set<Elem*>::const_iterator       el_it  = pos->second.begin();
        const std::set<Elem*>::const_iterator el_end = pos->second.end();

        // Loop over the elements and compute the values to be spread and the
        // positions of the quadrature points.
        int qp_offset = 0;
        std::vector<double> F_JxW_qp, X_qp;
        F_JxW_qp.reserve(n_vars*d_qrule->n_points()*pos->second.size());
        X_qp.reserve(NDIM*d_qrule->n_points()*pos->second.size());
        for ( ; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;

            fe->reinit(elem);
            for (unsigned int i = 0; i < n_vars; ++i)
            {
                dof_map.dof_indices(elem, dof_indices[i], i);
            }

            F_JxW_qp.resize(F_JxW_qp.size()+n_vars*d_qrule->n_points(),0.0);
            for (unsigned int qp = 0; qp < d_qrule->n_points(); ++qp)
            {
                for (unsigned int k = 0; k < phi.size(); ++k)
                {
                    for (unsigned int i = 0; i < n_vars; ++i)
                    {
                        F_JxW_qp[n_vars*(qp+qp_offset)+i] += F(dof_indices[i][k])*phi[k][qp]*JxW[qp];
                    }
                }
            }

            coords_fe->reinit(elem);
            for (unsigned int d = 0; d < dim; ++d)
            {
                coords_dof_map.dof_indices(elem, coords_dof_indices[d], d);
            }

            X_qp.resize(X_qp.size()+NDIM*d_qrule->n_points(),0.0);
            for (unsigned int qp = 0; qp < d_qrule->n_points(); ++qp)
            {
                for (unsigned int k = 0; k < coords_phi.size(); ++k)
                {
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        X_qp[NDIM*(qp+qp_offset)+i] += X(coords_dof_indices[i][k])*coords_phi[k][qp];
                    }
                }
            }

            qp_offset += d_qrule->n_points();
        }

        if (qp_offset == 0) continue;

        // Spread the values to the Cartesian grid patch.
        //
        // NOTE: Values are spread only from those quadrature points that are
        // within the ghost cell width of the patch interior.
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
    const MeshBase& mesh = d_es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    ExplicitSystem& system = d_es->get_system<ExplicitSystem>(system_name);
    const unsigned int n_vars = system.n_vars();
    const DofMap& dof_map = system.get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(n_vars);
    AutoPtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    fe->attach_quadrature_rule(d_qrule);
    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<std::vector<double> >& phi = fe->get_phi();
    AutoPtr<NumericVector<double> > rhs = F.clone();
    rhs->zero();

    System& coords_system = d_es->get_system<System>(COORDINATES_SYSTEM_NAME);
    const DofMap& coords_dof_map = coords_system.get_dof_map();
    std::vector<std::vector<unsigned int> > coords_dof_indices(dim);
    AutoPtr<FEBase> coords_fe(FEBase::build(dim, coords_dof_map.variable_type(0)));
    coords_fe->attach_quadrature_rule(d_qrule);
    const std::vector<std::vector<double> >& coords_phi = coords_fe->get_phi();
    if (close_X) X.close();
    coords_dof_map.enforce_constraints_exactly(coords_system, &X);

    for (unsigned int k = 0; k < f_refine_scheds.size(); ++k)
    {
        if (!f_refine_scheds[k].isNull()) f_refine_scheds[k]->fillData(fill_data_time);
    }

    // Loop over the patches to interpolate values to the element quadrature
    // points from the grid, then integrate the values at the element quadrature
    // points to obtain the nodal values on the Lagrangian mesh.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const int patch_num = patch->getPatchNumber();

        // The relevant collection of elements.
        std::map<int,std::set<Elem*> >::const_iterator pos = d_active_patch_elems.find(patch_num);
        if (pos == d_active_patch_elems.end()) continue;
        const std::set<Elem*>::const_iterator el_begin = pos->second.begin();
        const std::set<Elem*>::const_iterator el_end   = pos->second.end();

        // Loop over the elements and compute the positions of the quadrature points.
        int qp_offset = 0;
        std::vector<double> F_qp, X_qp;
        F_qp.reserve(n_vars*d_qrule->n_points()*pos->second.size());
        X_qp.reserve(NDIM*d_qrule->n_points()*pos->second.size());
        for (std::set<Elem*>::const_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;

            fe->reinit(elem);
            for (unsigned int i = 0; i < n_vars; ++i)
            {
                dof_map.dof_indices(elem, dof_indices[i], i);
            }

            F_qp.resize(F_qp.size()+n_vars*d_qrule->n_points(),0.0);

            coords_fe->reinit(elem);
            for (unsigned int d = 0; d < dim; ++d)
            {
                coords_dof_map.dof_indices(elem, coords_dof_indices[d], d);
            }

            X_qp.resize(X_qp.size()+NDIM*d_qrule->n_points(),0.0);
            for (unsigned int qp = 0; qp < d_qrule->n_points(); ++qp)
            {
                for (unsigned int k = 0; k < coords_phi.size(); ++k)
                {
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        X_qp[NDIM*(qp+qp_offset)+i] += X(coords_dof_indices[i][k])*coords_phi[k][qp];
                    }
                }
            }

            qp_offset += d_qrule->n_points();
        }

        if (qp_offset == 0) continue;

        // Interpolate the values from the Cartesian grid patch to the
        // quadrature points.
        //
        // NOTE: Values are interpolated only to those quadrature points that
        // are within the patch interior.
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
        for (std::set<Elem*>::const_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;

            fe->reinit(elem);
            for (unsigned int i = 0; i < n_vars; ++i)
            {
                dof_map.dof_indices(elem, dof_indices[i], i);
                rhs_e[i].resize(dof_indices[i].size());
            }

            for (unsigned int qp = 0; qp < d_qrule->n_points(); ++qp)
            {
                for (unsigned int k = 0; k < phi.size(); ++k)
                {
                    for (unsigned int i = 0; i < n_vars; ++i)
                    {
                        rhs_e[i](k) += F_qp[n_vars*(qp+qp_offset)+i]*phi[k][qp]*JxW[qp];
                    }
                }
            }

            for (unsigned int i = 0; i < n_vars; ++i)
            {
                dof_map.constrain_element_vector(rhs_e[i], dof_indices[i]);
                rhs->add_vector(rhs_e[i], dof_indices[i]);
            }

            qp_offset += d_qrule->n_points();
        }
    }

    // Assemble the rhs vector.
    rhs->close();

    // Solve for the nodal values.
    std::pair<LinearSolver<double>*,SparseMatrix<double>*> proj_solver_components = getL2ProjectionSolver(system_name, d_interp_uses_consistent_mass_matrix);
    LinearSolver<double>* solver = proj_solver_components.first;
    SparseMatrix<double>* M = proj_solver_components.second;
    const double tol = 1.0e-10;
    const unsigned int max_its = 100;
    dof_map.enforce_constraints_exactly(system, rhs.get());
    solver->solve(*M, *M, F, *rhs, tol, max_its);
    dof_map.enforce_constraints_exactly(system, &F);
    return;
}// interp

std::pair<LinearSolver<double>*,SparseMatrix<double>*>
FEDataManager::getL2ProjectionSolver(
    const std::string& system_name,
    const bool consistent_mass_matrix,
    const QuadratureType quad_type,
    const Order quad_order)
{
    if ((d_L2_projection_solvers.count(system_name) == 0 && d_L2_mass_matrices.count(system_name) == 0) ||
        (d_L2_consistent_mass_matrix[system_name] != consistent_mass_matrix))
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

        DenseMatrix<double> Me;

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
                Me.resize(dof_indices.size(), dof_indices.size());
                for (unsigned int i = 0; i < phi.size(); ++i)
                {
                    for (unsigned int j = 0; j < phi.size(); ++j)
                    {
                        for (unsigned int qp = 0; qp < qrule->n_points(); ++qp)
                        {
                            if (consistent_mass_matrix)
                            {
                                Me(i,j) += (phi[i][qp]*phi[j][qp])*JxW[qp];
                            }
                            else
                            {
                                Me(j,j) += (phi[i][qp]*phi[j][qp])*JxW[qp];
                            }
                        }
                    }
                }
                dof_map.constrain_element_matrix(Me, dof_indices);
                M->add_matrix(Me, dof_indices);
            }
        }

        // Assemble the matrix.
        M->close();

        // Setup the solver.
        solver->same_preconditioner = true;

        d_L2_projection_solvers[system_name] = solver;
        d_L2_mass_matrices[system_name] = M;
        d_L2_consistent_mass_matrix[system_name] = consistent_mass_matrix;
    }
    return std::make_pair(d_L2_projection_solvers[system_name], d_L2_mass_matrices[system_name]);
}// getL2ProjectionSolver

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
        std::map<int,std::set<Elem*> > active_level_elems;
        const IntVector<NDIM> ghost_width = 1;
        collectActivePatchElements(active_level_elems, level_number, ghost_width);
        std::vector<unsigned int> coords_ghost_dofs;
        collectGhostDOFIndices(coords_ghost_dofs, active_level_elems, COORDINATES_SYSTEM_NAME);

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
        X_ghost->init(X->size(), X->local_size(), coords_ghost_dofs, GHOSTED);
        X->localize(*X_ghost);
        coords_dof_map.enforce_constraints_exactly(coords_system, X_ghost.get());

        // Tag cells for refinement whenever they contain element quadrature
        // points.
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const int patch_num = patch->getPatchNumber();

            const Box<NDIM>& patch_box = patch->getBox();
            const CellIndex<NDIM>& patch_lower = patch_box.lower();
            const CellIndex<NDIM>& patch_upper = patch_box.upper();

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();
            const double* const patch_dx = patch_geom->getDx();

            Pointer<CellData<NDIM,int> > tag_data = patch->getPatchData(tag_index);

            std::set<Elem*>::const_iterator       el_it  = active_level_elems[patch_num].begin();
            const std::set<Elem*>::const_iterator el_end = active_level_elems[patch_num].end();
            for ( ; el_it != el_end; ++el_it)
            {
                const Elem* const elem = *el_it;

                coords_fe->reinit(elem);
                for (unsigned int d = 0; d < dim; ++d)
                {
                    coords_dof_map.dof_indices(elem, coords_dof_indices[d], d);
                }

                for (unsigned int qp = 0; qp < d_qrule->n_points(); ++qp)
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
                    tag_data->fill(1,Box<NDIM>::Box(i-Index<NDIM>(1),i+Index<NDIM>(1)));
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
      d_L2_projection_solvers(),
      d_L2_mass_matrices()
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
    for (std::map<std::string,LinearSolver<double>*>::iterator it = d_L2_projection_solvers.begin();
         it != d_L2_projection_solvers.end(); ++it)
    {
        delete (*it).second;
    }
    for (std::map<std::string,SparseMatrix<double>*>::iterator it = d_L2_mass_matrices.begin();
         it != d_L2_mass_matrices.end(); ++it)
    {
        delete (*it).second;
    }
    for (std::map<std::string,NumericVector<double>*>::iterator it = d_system_ghost_vec.begin();
         it != d_system_ghost_vec.end(); ++it)
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
    std::vector<std::vector<unsigned int> > coords_dof_indices(dim);
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
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const int patch_num = patch->getPatchNumber();
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
                std::map<int,std::set<Elem*> >::const_iterator pos = d_active_patch_elems.find(patch_num);
                if (pos == d_active_patch_elems.end()) continue;
                std::set<Elem*>::const_iterator       el_it  = pos->second.begin();
                const std::set<Elem*>::const_iterator el_end = pos->second.end();
                for ( ; el_it != el_end; ++el_it)
                {
                    const Elem* const elem = *el_it;

                    coords_fe->reinit(elem);
                    for (unsigned int d = 0; d < dim; ++d)
                    {
                        coords_dof_map.dof_indices(elem, coords_dof_indices[d], d);
                    }

                    for (unsigned int qp = 0; qp < d_qrule->n_points(); ++qp)
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
    std::map<int,std::set<Elem*> >& active_patch_elems,
    const int level_number,
    const IntVector<NDIM>& ghost_width)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif

    active_patch_elems.clear();

    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);

    // We initially associate an element with a Cartesian grid patch if the
    // element's bounding box intersects the patch interior grown by the
    // specified ghost cell width.
    std::vector<double> elem_bounds;
    computeActiveElementBoundingBoxes(elem_bounds);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
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
                active_patch_elems[patch_num].insert(elem);
            }
        }
    }

    // Recursively add/remove elements from the active sets that were generated
    // via the bounding box method.
    bool done = false;
    while (!done)
    {
        std::map<int,std::set<Elem*> > neighbor_elems;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const int patch_num = patch->getPatchNumber();

            std::map<int,std::set<Elem*> >::iterator pos = active_patch_elems.find(patch_num);
            if (pos == active_patch_elems.end()) continue;

            // Make a set of all active elements that neighbor the set of
            // elements currently associated with the patch.
            for (std::set<Elem*>::const_iterator cit = pos->second.begin(); cit != pos->second.end(); ++cit)
            {
                const Elem* const elem = *cit;
                for (unsigned n = 0; n < elem->n_neighbors(); ++n)
                {
                    Elem* const nghbr_elem = elem->neighbor(n);
                    if (nghbr_elem != NULL && nghbr_elem->active() &&
                        pos->second.count(nghbr_elem) == 0)
                    {
                        neighbor_elems[patch_num].insert(nghbr_elem);
                    }
                }
            }

            // Merge the set of neighboring elements with the set of elements
            // currently associated with the patch.
            pos->second.insert(neighbor_elems[patch_num].begin(), neighbor_elems[patch_num].end());
        }

        // Setup an appropriately ghosted coordinates vector.
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
        X_ghost->init(X->size(), X->local_size(), coords_ghost_dofs, GHOSTED);
        X->localize(*X_ghost);
        coords_dof_map.enforce_constraints_exactly(coords_system, X_ghost.get());

        // Keep only those elements that have a quadrature point on the local
        // patch.
        std::map<int,std::set<Elem*> > new_active_patch_elems;
        bool inserted_ngbhr_elem = false;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const int patch_num = patch->getPatchNumber();

            std::map<int,std::set<Elem*> >::const_iterator pos = active_patch_elems.find(patch_num);
            if (pos == active_patch_elems.end()) continue;

            const Box<NDIM>& patch_box = patch->getBox();
            const Box<NDIM>  ghost_box = Box<NDIM>::grow(patch_box, ghost_width);
            const CellIndex<NDIM>& patch_lower = patch_box.lower();
            const CellIndex<NDIM>& patch_upper = patch_box.upper();

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();
            const double* const patch_dx = patch_geom->getDx();

            std::set<Elem*>::const_iterator       el_it  = pos->second.begin();
            const std::set<Elem*>::const_iterator el_end = pos->second.end();
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
                        new_active_patch_elems[patch_num].insert(elem);
                        if (neighbor_elems[patch_num].count(elem) > 0) inserted_ngbhr_elem = true;
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
    return;
}// collectActivePatchElements

void
FEDataManager::collectGhostDOFIndices(
    std::vector<unsigned int>& ghost_dofs,
    const std::map<int,std::set<Elem*> >& active_patch_elems,
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
    for (std::map<int,std::set<Elem*> >::const_iterator cit1 = active_patch_elems.begin(); cit1 != active_patch_elems.end(); ++cit1)
    {
        const std::set<Elem*>& elems = (*cit1).second;
        for (std::set<Elem*>::const_iterator cit2 = elems.begin(); cit2 != elems.end(); ++cit2)
        {
            const Elem* const elem = *cit2;

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
    }
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

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
