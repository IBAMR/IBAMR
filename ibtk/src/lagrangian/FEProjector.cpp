// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
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

#include <IBTK_config.h>

#include <ibtk/FEDataManager.h>
#include <ibtk/FEProjector.h>
#include <ibtk/namespaces.h> // IWYU pragma: keep

#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

#include <libmesh/boundary_info.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_build_l2_projection_solver;
static Timer* t_build_diagonal_l2_mass_matrix;
static Timer* t_compute_l2_projection;

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
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

FEProjector::FEProjector(EquationSystems* equation_systems, const bool enable_logging)
    : d_fe_data(new FEData("FEProjector", /*register_for_restart*/ false)), d_enable_logging(enable_logging)
{
    d_fe_data->setEquationSystems(equation_systems, /*level_number*/ 0);
    IBTK_DO_ONCE(t_build_l2_projection_solver =
                     TimerManager::getManager()->getTimer("IBTK::FEProjector::buildL2ProjectionSolver()");
                 t_build_diagonal_l2_mass_matrix =
                     TimerManager::getManager()->getTimer("IBTK::FEProjector::buildDiagonalL2MassMatrix()");
                 t_compute_l2_projection =
                     TimerManager::getManager()->getTimer("IBTK::FEProjector::computeL2Projection()");)
}

FEProjector::FEProjector(std::shared_ptr<FEData> fe_data, const bool enable_logging)
    : d_fe_data(fe_data), d_enable_logging(enable_logging)
{
    TBOX_ASSERT(d_fe_data);
    IBTK_DO_ONCE(t_build_l2_projection_solver =
                     TimerManager::getManager()->getTimer("IBTK::FEProjector::buildL2ProjectionSolver()");
                 t_build_diagonal_l2_mass_matrix =
                     TimerManager::getManager()->getTimer("IBTK::FEProjector::buildDiagonalL2MassMatrix()");
                 t_compute_l2_projection =
                     TimerManager::getManager()->getTimer("IBTK::FEProjector::computeL2Projection()");)
}

std::pair<PetscLinearSolver<double>*, PetscMatrix<double>*>
FEProjector::buildL2ProjectionSolver(const std::string& system_name)
{
    IBTK_TIMER_START(t_build_l2_projection_solver);

    if (!d_L2_proj_solver.count(system_name) || !d_L2_proj_matrix.count(system_name))
    {
        if (d_enable_logging)
        {
            plog << "FEProjector::buildL2ProjectionSolver(): building L2 projection solver for system: " << system_name
                 << "\n";
        }

        // Extract the mesh.
        const MeshBase& mesh = d_fe_data->getEquationSystems()->get_mesh();
        const Parallel::Communicator& comm = mesh.comm();
        const unsigned int dim = mesh.mesh_dimension();

        // Extract the FE system and DOF map, and setup the FE object.
        System& system = d_fe_data->getEquationSystems()->get_system(system_name);
        const int sys_num = system.number();
        DofMap& dof_map = system.get_dof_map();
        FEData::SystemDofMapCache& dof_map_cache = *d_fe_data->getDofMapCache(system_name);
        dof_map.compute_sparsity(mesh);
        FEType fe_type = dof_map.variable_type(0);
        std::unique_ptr<QBase> qrule = fe_type.default_quadrature_rule(dim);
        std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<double>& JxW = fe->get_JxW();
        const std::vector<std::vector<double> >& phi = fe->get_phi();

        // Build solver components.
        std::unique_ptr<PetscLinearSolver<double> > solver(new PetscLinearSolver<double>(comm));
        solver->init();

        std::unique_ptr<PetscMatrix<double> > M_mat(new PetscMatrix<double>(comm));
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
                static const std::array<boundary_id_type, 3> dirichlet_bdry_id_set = {
                    FEDataManager::ZERO_DISPLACEMENT_X_BDRY_ID,
                    FEDataManager::ZERO_DISPLACEMENT_Y_BDRY_ID,
                    FEDataManager::ZERO_DISPLACEMENT_Z_BDRY_ID
                };
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
}

PetscVector<double>*
FEProjector::buildDiagonalL2MassMatrix(const std::string& system_name)
{
    IBTK_TIMER_START(t_build_diagonal_l2_mass_matrix);

    if (!d_L2_proj_matrix_diag.count(system_name))
    {
        if (d_enable_logging)
        {
            plog << "FEProjector::buildDiagonalL2MassMatrix(): building diagonal L2 mass matrix for system: "
                 << system_name << "\n";
        }

        // Extract the mesh.
        const MeshBase& mesh = d_fe_data->getEquationSystems()->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        // Extract the FE system and DOF map, and setup the FE object.
        System& system = d_fe_data->getEquationSystems()->get_system(system_name);
        const int sys_num = system.number();
        DofMap& dof_map = system.get_dof_map();
        FEData::SystemDofMapCache& dof_map_cache = *d_fe_data->getDofMapCache(system_name);
        dof_map.compute_sparsity(mesh);
        FEType fe_type = dof_map.variable_type(0);
        std::unique_ptr<QBase> qrule = fe_type.default_quadrature_rule(dim);
        std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<double>& JxW = fe->get_JxW();
        const std::vector<std::vector<double> >& phi = fe->get_phi();

        // Build solver components.
        std::unique_ptr<PetscVector<double> > M_vec(
            static_cast<PetscVector<double>*>(system.solution->zero_clone().release()));

        // Loop over the mesh to construct the system matrix.
        DenseMatrix<double> M_e;
        DenseVector<double> M_e_vec;
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

                // We explicitly do *not* apply constraints in this function
                // for two reasons:
                //
                // 1. Since we have a diagonal matrix no DoF is coupled to any
                //    other DoF: therefore constraints aren't necessary.
                // 2. We want to explicitly ignore constraints in the case
                //    of periodic boundary conditions. This is because we want
                //    to calculate volume of the elements containing a
                //    node. For example: if we had a periodic dof X which
                //    should be equal to a periodic dof Y with the following
                //    geometry:
                //
                //    +---+                 +----------+
                //    |   |                 |          |
                //    X---+   (etc)         +----------Y
                //    |   |                 |          |
                //    +---+                 +----------+
                //
                //    then we would still want to have the entry in the mass
                //    matrix corresponding to X be the volume of the two cells
                //    on which its shape function has support, and we would
                //    want the entry in the mass matrix corresponding to Y to
                //    have the same property. Clearly, these values should not
                //    be equal - in fact, the periodic constraints are totally
                //    irrelevant to what we are trying to compute in the first
                //    place.
                M_vec->add_vector(M_e_vec, dof_indices[var_num]);
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
                static const std::array<boundary_id_type, 3> dirichlet_bdry_id_set = {
                    FEDataManager::ZERO_DISPLACEMENT_X_BDRY_ID,
                    FEDataManager::ZERO_DISPLACEMENT_Y_BDRY_ID,
                    FEDataManager::ZERO_DISPLACEMENT_Z_BDRY_ID
                };
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
}

bool
FEProjector::computeL2Projection(PetscVector<double>& U_vec,
                                 PetscVector<double>& F_vec,
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
    const System& system = d_fe_data->getEquationSystems()->get_system(system_name);
    const DofMap& dof_map = system.get_dof_map();
    if (consistent_mass_matrix)
    {
        std::pair<PetscLinearSolver<double>*, PetscMatrix<double>*> proj_solver_components =
            buildL2ProjectionSolver(system_name);
        PetscLinearSolver<double>* solver = proj_solver_components.first;
        PetscMatrix<double>* M_mat = proj_solver_components.second;
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
        PetscVector<double>* M_diag_vec = buildDiagonalL2MassMatrix(system_name);
        ierr = VecPointwiseDivide(U_vec.vec(), F_vec.vec(), M_diag_vec->vec());
        IBTK_CHKERRQ(ierr);
        converged = true;
    }
    if (close_U) U_vec.close();
    dof_map.enforce_constraints_exactly(system, &U_vec);

    IBTK_TIMER_STOP(t_compute_l2_projection);
    return converged;
}

void
FEProjector::setLoggingEnabled(const bool enable_logging)
{
    d_enable_logging = enable_logging;
}

bool
FEProjector::getLoggingEnabled() const
{
    return d_enable_logging;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
