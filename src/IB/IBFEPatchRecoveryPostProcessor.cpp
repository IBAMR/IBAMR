// Filename: IBFEPatchRecoveryPostProcessor.cpp
// Created on 2 Jul 2013 by Boyce Griffith
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

#include "IBAMR_config.h"
#include "Eigen/Dense"
#include "ibamr/IBFEPatchRecoveryPostProcessor.h"
#include "SAMRAI_config.h"
#include "boost/multi_array.hpp"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/libmesh_utilities.h"
#include "libmesh/boundary_info.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/mesh.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/quadrature.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/string_to_enum.h"

using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
unsigned int
num_polynomial_basis_fcns(const unsigned int dim, const unsigned int order)
{
    unsigned int num_basis_fcns = 0;
    unsigned int order_p_1 = order + 1;
    switch (dim)
    {
    case 1:
        num_basis_fcns = order_p_1;
        break;
    case 2:
        num_basis_fcns = order_p_1 * order_p_1 + order_p_1;
        num_basis_fcns /= 2;
        break;
    case 3:
        num_basis_fcns = order_p_1 * order_p_1 * order_p_1 + 3 * order_p_1 * order_p_1 + 2 * order_p_1;
        num_basis_fcns /= 6;
        break;
    default:
        TBOX_ERROR("only supports dim = 1, 2, or 3\n");
    }
    return num_basis_fcns;
} // num_polynomial_basis_fcns

void
evaluate_polynomial_basis_fcns(Eigen::VectorXd& P,
                               const libMesh::Point& x_center,
                               const libMesh::Point& x_eval,
                               const unsigned int dim,
                               const unsigned int order)
{
    TBOX_ASSERT(static_cast<unsigned int>(P.size()) == num_polynomial_basis_fcns(dim, order));

    // Compute powers of the components of x up to the specified order.
    libMesh::Point x = x_center - x_eval;
    boost::multi_array<double, 2> x_pow(boost::extents[dim][order + 1]);
    for (unsigned int d = 0; d < dim; ++d)
    {
        x_pow[d][0] = 1.0;
        for (unsigned int k = 1; k <= order; ++k)
        {
            x_pow[d][k] = x(d) * x_pow[d][k - 1];
        }
    }

    // Evaluate the complete polynomial basis functions of the specified order.
    static const unsigned int X_IDX = 0;
    static const unsigned int Y_IDX = 1;
    static const unsigned int Z_IDX = 2;
    switch (dim)
    {
    case 1:
        for (unsigned int total_pow = 0, k = 0; total_pow <= order; ++total_pow, ++k)
        {
            unsigned int x_exp = total_pow;
            P(k) = x_pow[X_IDX][x_exp];
        }
        break;
    case 2:
        for (unsigned int total_pow = 0, k = 0; total_pow <= order; ++total_pow)
        {
            for (unsigned int x_exp = 0; x_exp <= total_pow; ++x_exp, ++k)
            {
                unsigned int y_exp = total_pow - x_exp;
                P(k) = x_pow[X_IDX][x_exp] * x_pow[Y_IDX][y_exp];
            }
        }
        break;
    case 3:
        for (unsigned int total_pow = 0, k = 0; total_pow <= order; ++total_pow)
        {
            for (unsigned int x_exp = 0; x_exp <= total_pow; ++x_exp)
            {
                for (unsigned int y_exp = 0; y_exp <= total_pow - x_exp; ++y_exp, ++k)
                {
                    unsigned int z_exp = total_pow - (x_exp + y_exp);
                    P(k) = x_pow[X_IDX][x_exp] * x_pow[Y_IDX][y_exp] * x_pow[Z_IDX][z_exp];
                }
            }
        }
        break;
    default:
        TBOX_ERROR("only supports dim = 1, 2, or 3\n");
    }
    return;
} // evaluate_polynomial_basis_fcns
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEPatchRecoveryPostProcessor::IBFEPatchRecoveryPostProcessor(MeshBase* mesh, FEDataManager* fe_data_manager)
    : d_mesh(mesh),
      d_fe_data_manager(fe_data_manager),
      d_periodic_boundaries(NULL),
      d_interp_order(INVALID_ORDER),
      d_quad_order(INVALID_ORDER)
{
    // Active local elements.
    const MeshBase::const_element_iterator el_begin = d_mesh->active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = d_mesh->active_local_elements_end();

    // Determine the number of quadrature/interpolation points in each element.
    //
    // We use full-order Gaussian quadrature rules (i.e. third-order Gauss
    // quadrature for first-order elements and fifth-order Gauss quadrature for
    // second-order elements) in all elements to avoid special treatment at
    // boundary nodes.
    bool first_order_elems = false;
    bool second_order_elems = false;
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const unsigned int dim = elem->dim();
        TBOX_ASSERT(dim == d_mesh->mesh_dimension());
        Order elem_order = elem->default_order();
        if (elem_order == FIRST) first_order_elems = true;
        if (elem_order == SECOND) second_order_elems = true;
        TBOX_ASSERT(elem_order == FIRST || elem_order == SECOND);
    }
    if (first_order_elems && second_order_elems)
    {
        TBOX_ERROR("cannot have both first- and second-order elements in the same mesh.\n");
    }
    d_interp_order = first_order_elems ? FIRST : SECOND;
    d_quad_order = first_order_elems ? THIRD : FIFTH;
    return;
} // IBFEPatchRecoveryPostProcessor

IBFEPatchRecoveryPostProcessor::~IBFEPatchRecoveryPostProcessor()
{
    // intentionally blank
    return;
} // ~IBFEPatchRecoveryPostProcessor

void
IBFEPatchRecoveryPostProcessor::initializeFEData(const PeriodicBoundaries* const periodic_boundaries)
{
    d_periodic_boundaries = periodic_boundaries;

    const Parallel::Communicator& comm = d_mesh->comm();
    const int mpi_rank = comm.rank();
    const int mpi_size = comm.size();

    // Active local elements.
    const MeshBase::const_element_iterator el_begin = d_mesh->active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = d_mesh->active_local_elements_end();

    // Determine the element patches associated with each node N, which is
    // defined to be the collection of elements that contain node N.
    //
    // Unlike the standard Z-Z patch recovery algorithm, we use "tight" element
    // patches for non-vertex nodes.
    AutoPtr<PointLocatorBase> point_locator = PointLocatorBase::build(TREE, *d_mesh);
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        for (unsigned int n = 0; n < elem->n_nodes(); ++n)
        {
            // Only set up patches for local nodes.
            const Node* const node = elem->get_node(n);
            if (node->processor_id() != mpi_rank) continue;

            // Only set up patches once for each node.
            const dof_id_type node_id = node->id();
            if (d_local_elem_patches.find(node_id) != d_local_elem_patches.end()) continue;

            // Find the elements that touch this node.
            ElemPatch& elem_patch = d_local_elem_patches[node_id];
            std::set<const Elem*> elems;
            elem->find_point_neighbors(*node, elems);
            for (std::set<const Elem*>::const_iterator it = elems.begin(); it != elems.end(); ++it)
            {
                elem_patch.insert(boost::make_tuple(*it, CompositePeriodicMapping(), CompositePeriodicMapping()));
            }

            // Account for periodic boundaries.
            bool done = false || !d_periodic_boundaries;
            while (!done)
            {
                ElemPatch periodic_neighbors;
                for (ElemPatch::const_iterator it = elem_patch.begin(); it != elem_patch.end(); ++it)
                {
                    const Elem* const elem = it->get<0>();
                    const CompositePeriodicMapping& forward_mapping = it->get<1>();
                    const CompositePeriodicMapping& inverse_mapping = it->get<2>();
                    const libMesh::Point p = apply_composite_periodic_mapping(forward_mapping, *node);
                    TBOX_ASSERT(elem->contains_point(p));
                    for (unsigned int i = 0; i < elem->n_neighbors(); ++i)
                    {
                        if (!elem->neighbor(i))
                        {
                            const std::vector<boundary_id_type>& boundary_ids =
                                d_mesh->boundary_info->boundary_ids(elem, i);
                            for (std::vector<boundary_id_type>::const_iterator j = boundary_ids.begin();
                                 j != boundary_ids.end();
                                 ++j)
                            {
                                const boundary_id_type boundary_id = *j;
                                const PeriodicBoundaryBase* const periodic_boundary =
                                    d_periodic_boundaries->boundary(boundary_id);
                                if (periodic_boundary)
                                {
                                    const libMesh::Point periodic_image = periodic_boundary->get_corresponding_pos(p);
                                    const Elem* neighbor =
                                        d_periodic_boundaries->neighbor(boundary_id, *point_locator, elem, i);
                                    if (elem->level() < neighbor->level())
                                    {
                                        neighbor = neighbor->parent();
                                    }
                                    if (neighbor->contains_point(periodic_image))
                                    {
                                        std::set<const Elem*> elems;
                                        neighbor->find_point_neighbors(periodic_image, elems);
                                        for (std::set<const Elem*>::const_iterator k = elems.begin(); k != elems.end();
                                             ++k)
                                        {
                                            const Elem* const elem = *k;
                                            if (elem_patch.find(elem) == elem_patch.end() &&
                                                periodic_neighbors.find(elem) == periodic_neighbors.end())
                                            {
                                                CompositePeriodicMapping forward = forward_mapping;
                                                CompositePeriodicMapping inverse = inverse_mapping;
                                                forward.insert(
                                                    forward.end(),
                                                    periodic_boundary->clone(PeriodicBoundaryBase::FORWARD).release());
                                                inverse.insert(
                                                    inverse.begin(),
                                                    periodic_boundary->clone(PeriodicBoundaryBase::INVERSE).release());
                                                periodic_neighbors.insert(boost::make_tuple(elem, forward, inverse));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                elem_patch.insert(periodic_neighbors.begin(), periodic_neighbors.end());
                periodic_neighbors.clear();
                done = periodic_neighbors.empty();
            }
        }
    }

    // Setup mappings used to fill global indexing data structures.
    //
    // We use full-order Gaussian quadrature rules (i.e. third-order Gauss
    // quadrature for first-order elements and fifth-order Gauss quadrature for
    // second-order elements) in all elements to avoid special treatment at
    // boundary nodes.
    d_n_qp_global = 0;
    d_n_qp_local = 0;
    d_qp_global_offset = 0;
    const unsigned int n_elem = d_mesh->n_elem();
    d_elem_n_qp.resize(n_elem, 0);
    d_elem_qp_global_offset.resize(n_elem, 0);
    d_elem_qp_local_offset.resize(n_elem, 0);
    AutoPtr<QBase> qrule;
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const unsigned int dim = elem->dim();
        bool reinit_qrule = false;
        if (!qrule.get() || qrule->get_dim() != dim || qrule->get_order() != d_quad_order)
        {
            qrule = QBase::build(QGAUSS, dim, d_quad_order);
            reinit_qrule = true;
        }
        else if (qrule->get_elem_type() != elem->type() || qrule->get_p_level() != elem->p_level())
        {
            reinit_qrule = true;
        }
        if (reinit_qrule) qrule->init(elem->type(), elem->p_level());
        unsigned int n_qp = qrule->n_points();
        const dof_id_type elem_id = elem->id();
        d_elem_n_qp[elem_id] = n_qp;
        d_elem_qp_local_offset[elem_id] = d_n_qp_local;
        d_n_qp_local += n_qp;
        d_elem_sigma[elem->id()].resize(n_qp);
        d_elem_pressure[elem->id()].resize(n_qp);
    }
    std::vector<int> n_qp_per_proc(mpi_size);
    n_qp_per_proc[mpi_rank] = d_n_qp_local;
    comm.sum(n_qp_per_proc);
    d_qp_global_offset = std::accumulate(n_qp_per_proc.begin(), n_qp_per_proc.begin() + mpi_rank, 0);
    d_n_qp_global = std::accumulate(n_qp_per_proc.begin() + mpi_rank, n_qp_per_proc.end(), d_qp_global_offset);
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const dof_id_type elem_id = elem->id();
        d_elem_qp_global_offset[elem_id] = d_elem_qp_local_offset[elem_id] + d_qp_global_offset;
    }
    comm.sum(d_elem_qp_global_offset);
    comm.sum(d_elem_qp_local_offset);

    // Set up element patch L2 projection matrices.
    unsigned int dim = d_mesh->mesh_dimension();
    const unsigned int num_basis_fcns = num_polynomial_basis_fcns(dim, d_interp_order);
    Eigen::MatrixXd M(num_basis_fcns, num_basis_fcns);
    Eigen::VectorXd P(num_basis_fcns);
    AutoPtr<FEBase> fe(FEBase::build(dim, FEType(d_interp_order, LAGRANGE)));
    const std::vector<libMesh::Point>& q_point = fe->get_xyz();
    qrule = QBase::build(QGAUSS, dim, d_quad_order);
    fe->attach_quadrature_rule(qrule.get());
    d_local_patch_proj_solver.resize(d_local_elem_patches.size());
    unsigned int k = 0;
    for (std::map<dof_id_type, ElemPatch>::iterator it = d_local_elem_patches.begin(); it != d_local_elem_patches.end();
         ++it, ++k)
    {
        const dof_id_type node_id = it->first;
        const Node& node = d_mesh->node(node_id);
        ElemPatch& elem_patch = it->second;
        M.setZero();
        for (ElemPatch::const_iterator el_it = elem_patch.begin(); el_it != elem_patch.end(); ++el_it)
        {
            const Elem* const elem = el_it->get<0>();
            const CompositePeriodicMapping& inverse_mapping = el_it->get<2>();
            fe->reinit(elem);
            for (unsigned int qp = 0; qp < qrule->n_points(); ++qp)
            {
                evaluate_polynomial_basis_fcns(
                    P, node, apply_composite_periodic_mapping(inverse_mapping, q_point[qp]), dim, d_interp_order);
                M += P * P.transpose();
            }
        }
        d_local_patch_proj_solver[k] = M.colPivHouseholderQr();
        if (!d_local_patch_proj_solver[k].isInvertible())
        {
            TBOX_ERROR(
                "IBFEPatchRecoveryPostProcessor could not construct L2 reconstruction for "
                "element "
                "patch associated with node "
                << node_id
                << "\n");
        }
    }
    return;
} // initializeFEData

System*
IBFEPatchRecoveryPostProcessor::initializeCauchyStressSystem()
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    System* sigma_system = &equation_systems->add_system<System>("CAUCHY_STRESS_RECOVERY_SYSTEM");
    for (unsigned int i = 0; i < NDIM; ++i)
    {
        for (unsigned int j = i; j < NDIM; ++j)
        {
            std::ostringstream os;
            os << "sigma_" << (i == 0 ? "x" : i == 1 ? "y" : "z") << (j == 0 ? "x" : j == 1 ? "y" : "z");
            sigma_system->add_variable(os.str(), d_interp_order, LAGRANGE);
        }
    }
    return sigma_system;
} // initializeCauchyStressSystem

System*
IBFEPatchRecoveryPostProcessor::initializePressureSystem()
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    System* p_system = &equation_systems->add_system<System>("PRESSURE_RECOVERY_SYSTEM");
    p_system->add_variable("p", d_interp_order, LAGRANGE);
    return p_system;
} // initializePressureSystem

void
IBFEPatchRecoveryPostProcessor::registerCauchyStressValue(const Elem* const elem,
                                                          const QBase* const qrule,
                                                          const unsigned int qp,
                                                          const TensorValue<double>& sigma)
{
    const Parallel::Communicator& comm = d_mesh->comm();
    if (elem->processor_id() != comm.rank() || !elem->active())
    {
        TBOX_ERROR("must register stresses only for active local elements\n");
    }
    TBOX_ASSERT(elem->default_order() == d_interp_order);
    TBOX_ASSERT(qrule->type() == QGAUSS);
    TBOX_ASSERT(qrule->get_order() == d_quad_order);
    TBOX_ASSERT(qrule->get_elem_type() == elem->type());
    TBOX_ASSERT(qrule->get_p_level() == elem->p_level());
    TBOX_ASSERT(qp < d_elem_n_qp[elem->id()]);
    d_elem_sigma[elem->id()][qp] = sigma;
    return;
} // registerCauchyStressValue

void
IBFEPatchRecoveryPostProcessor::registerPressureValue(const Elem* const elem,
                                                      const QBase* const qrule,
                                                      const unsigned int qp,
                                                      const double p)
{
    const Parallel::Communicator& comm = d_mesh->comm();
    if (elem->processor_id() != comm.rank() || !elem->active())
    {
        TBOX_ERROR("must register pressures only for active local elements\n");
    }
    TBOX_ASSERT(elem->default_order() == d_interp_order);
    TBOX_ASSERT(qrule->type() == QGAUSS);
    TBOX_ASSERT(qrule->get_order() == d_quad_order);
    TBOX_ASSERT(qrule->get_elem_type() == elem->type());
    TBOX_ASSERT(qrule->get_p_level() == elem->p_level());
    TBOX_ASSERT(qp < d_elem_n_qp[elem->id()]);
    d_elem_pressure[elem->id()][qp] = p;
    return;
} // registerPressureValue

void
IBFEPatchRecoveryPostProcessor::reconstructCauchyStress(System& sigma_system)
{
    const unsigned int sigma_sys_num = sigma_system.number();
    NumericVector<double>& sigma_vec = *sigma_system.solution;

    // Communicate the stored values of the Cauchy stress.
    static const unsigned int NVARS = (NDIM * (NDIM + 1)) / 2;
    std::vector<double> sigma_vals(NVARS * d_n_qp_global, 0.0);
    const MeshBase::const_element_iterator el_begin = d_mesh->active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = d_mesh->active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const dof_id_type elem_id = elem->id();
        const int global_offset = d_elem_qp_global_offset[elem_id];
        for (unsigned int qp = 0; qp < d_elem_n_qp[elem_id]; ++qp)
        {
            const TensorValue<double>& stress = d_elem_sigma[elem_id][qp];
            for (unsigned int i = 0, k = 0; i < NDIM; ++i)
            {
                for (unsigned int j = i; j < NDIM; ++j, ++k)
                {
                    sigma_vals[NVARS * (global_offset + qp) + k] = stress(i, j);
                }
            }
        }
    }
    const Parallel::Communicator& comm = d_mesh->comm();
    comm.sum(sigma_vals);

    // Perform element patch L2 projections.
    const unsigned int dim = d_mesh->mesh_dimension();
    const unsigned int num_basis_fcns = num_polynomial_basis_fcns(dim, d_interp_order);
    Eigen::VectorXd P(num_basis_fcns), a(num_basis_fcns), f(num_basis_fcns);
    AutoPtr<FEBase> fe(FEBase::build(dim, FEType(d_interp_order, LAGRANGE)));
    const std::vector<libMesh::Point>& q_point = fe->get_xyz();
    AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, d_quad_order);
    fe->attach_quadrature_rule(qrule.get());
    unsigned int k = 0;
    for (std::map<dof_id_type, ElemPatch>::const_iterator it = d_local_elem_patches.begin();
         it != d_local_elem_patches.end();
         ++it, ++k)
    {
        const dof_id_type node_id = it->first;
        const Node& node = d_mesh->node(node_id);
        const ElemPatch& elem_patch = it->second;
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& patch_proj_solver = d_local_patch_proj_solver[k];
        for (unsigned int var = 0; var < NVARS; ++var)
        {
            // Solve for the coefficients of the reconstruction.
            f.setZero();
            for (ElemPatch::const_iterator el_it = elem_patch.begin(); el_it != elem_patch.end(); ++el_it)
            {
                const Elem* const elem = el_it->get<0>();
                const CompositePeriodicMapping& inverse_mapping = el_it->get<2>();
                const dof_id_type elem_id = elem->id();
                const int global_offset = d_elem_qp_global_offset[elem_id];
                fe->reinit(elem);
                for (unsigned int qp = 0; qp < qrule->n_points(); ++qp)
                {
                    evaluate_polynomial_basis_fcns(
                        P, node, apply_composite_periodic_mapping(inverse_mapping, q_point[qp]), dim, d_interp_order);
                    f += P * sigma_vals[NVARS * (global_offset + qp) + var];
                }
            }
            a = patch_proj_solver.solve(f);

            // Evaluate the reconstruction at the node.
            const int dof_index = node.dof_number(sigma_sys_num, var, 0);
            sigma_vec.set(dof_index, a(0));
        }
    }
    return;
} // reconstructCauchyStress

void
IBFEPatchRecoveryPostProcessor::reconstructPressure(System& p_system)
{
    const unsigned int p_sys_num = p_system.number();
    NumericVector<double>& p_vec = *p_system.solution;

    // Communicate the stored values of the Cauchy stress.
    std::vector<double> pressure_vals(d_n_qp_global, 0.0);
    const MeshBase::const_element_iterator el_begin = d_mesh->active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = d_mesh->active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const dof_id_type elem_id = elem->id();
        const int global_offset = d_elem_qp_global_offset[elem_id];
        for (unsigned int qp = 0; qp < d_elem_n_qp[elem_id]; ++qp)
        {
            pressure_vals[global_offset + qp] = d_elem_pressure[elem_id][qp];
        }
    }
    const Parallel::Communicator& comm = d_mesh->comm();
    comm.sum(pressure_vals);

    // Perform element patch L2 projections.
    const unsigned int dim = d_mesh->mesh_dimension();
    const unsigned int num_basis_fcns = num_polynomial_basis_fcns(dim, d_interp_order);
    Eigen::VectorXd P(num_basis_fcns), a(num_basis_fcns), f(num_basis_fcns);
    AutoPtr<FEBase> fe(FEBase::build(dim, FEType(d_interp_order, LAGRANGE)));
    const std::vector<libMesh::Point>& q_point = fe->get_xyz();
    AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, d_quad_order);
    fe->attach_quadrature_rule(qrule.get());
    unsigned int k = 0;
    for (std::map<dof_id_type, ElemPatch>::const_iterator it = d_local_elem_patches.begin();
         it != d_local_elem_patches.end();
         ++it, ++k)
    {
        const dof_id_type node_id = it->first;
        const Node& node = d_mesh->node(node_id);
        const ElemPatch& elem_patch = it->second;
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& patch_proj_solver = d_local_patch_proj_solver[k];

        // Solve for the coefficients of the reconstruction.
        f.setZero();
        for (ElemPatch::const_iterator el_it = elem_patch.begin(); el_it != elem_patch.end(); ++el_it)
        {
            const Elem* const elem = el_it->get<0>();
            const CompositePeriodicMapping& inverse_mapping = el_it->get<2>();
            const dof_id_type elem_id = elem->id();
            const int global_offset = d_elem_qp_global_offset[elem_id];
            fe->reinit(elem);
            for (unsigned int qp = 0; qp < qrule->n_points(); ++qp)
            {
                evaluate_polynomial_basis_fcns(
                    P, node, apply_composite_periodic_mapping(inverse_mapping, q_point[qp]), dim, d_interp_order);
                f += P * pressure_vals[global_offset + qp];
            }
        }
        a = patch_proj_solver.solve(f);

        // Evaluate the reconstruction at the node.
        const unsigned int var = 0;
        const int dof_index = node.dof_number(p_sys_num, var, 0);
        p_vec.set(dof_index, a(0));
    }
    return;
} // reconstructPressure

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
