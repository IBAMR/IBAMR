// Filename: IBFEPatchRecoveryPostProcessor.cpp
// Created on 2 Jul 2013 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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
#include "libmesh/string_to_enum.h"

using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
unsigned int
num_basis_functions(
    const unsigned int dim,
    const unsigned int order)
{
    unsigned int num_basis = 0;
    unsigned int order_p_1 = order+1;
    switch (dim)
    {
        case 1:
            num_basis = order_p_1;
            break;
        case 2:
            num_basis = order_p_1*order_p_1+order_p_1;
            num_basis /= 2;
            break;
        case 3:
            num_basis = order_p_1*order_p_1*order_p_1+3*order_p_1*order_p_1+2*order_p_1;
            num_basis /= 6;
            break;
        default:
            TBOX_ERROR("only supports dim = 1, 2, or 3\n");
    }
    return num_basis;
}// num_basis_functions

void
evaluate_basis_functions(
    Eigen::VectorXd& P,
    const libMesh::Point& x_center,
    const libMesh::Point& x_eval,
    const unsigned int dim,
    const unsigned int order)
{
    TBOX_ASSERT(static_cast<unsigned int>(P.size()) == num_basis_functions(dim,order));

    // Compute powers of the components of x up to the specified order.
    libMesh::Point x = x_center-x_eval;
    boost::multi_array<double,2> x_pow(boost::extents[dim][order+1]);
    for (unsigned int d = 0; d < dim; ++d)
    {
        x_pow[d][0] = 1.0;
        for (unsigned int k = 1; k <= order; ++k)
        {
            x_pow[d][k] = x(d)*x_pow[d][k-1];
        }
    }

    // Evaluate complete polynomials up to the specified order.
    switch (dim)
    {
        case 1:
            for (unsigned int total_pow = 0, k = 0; total_pow <= order; ++total_pow, ++k)
            {
                unsigned int x_exp = total_pow;
                P(k) = x_pow[0][x_exp];
            }
            break;
        case 2:
            for (unsigned int total_pow = 0, k = 0; total_pow <= order; ++total_pow)
            {
                for (unsigned int x_exp = 0; x_exp <= total_pow; ++x_exp, ++k)
                {
                    unsigned int y_exp = total_pow-x_exp;
                    P(k) = x_pow[0][x_exp]*x_pow[1][y_exp];
                }
            }
            break;
        case 3:
            for (unsigned int total_pow = 0, k = 0; total_pow <= order; ++total_pow)
            {
                for (unsigned int x_exp = 0; x_exp <= total_pow; ++x_exp)
                {
                    for (unsigned int y_exp = 0; y_exp <= total_pow-x_exp; ++y_exp, ++k)
                    {
                        unsigned int z_exp = total_pow-(x_exp+y_exp);
                        P(k) = x_pow[0][x_exp]*x_pow[1][y_exp]*x_pow[2][z_exp];
                    }
                }
            }
            break;
        default:
            TBOX_ERROR("only supports dim = 1, 2, or 3\n");
    }
    return;
}// evaluate_basis_functions
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEPatchRecoveryPostProcessor::IBFEPatchRecoveryPostProcessor(
    MeshBase* mesh,
    FEDataManager* fe_data_manager)
    : d_mesh(mesh),
      d_fe_data_manager(fe_data_manager)
{
    const int rank = libMesh::processor_id();
    const int size = libMesh::n_processors();

    // Active local elements.
    const MeshBase::const_element_iterator el_begin = d_mesh->active_local_elements_begin();
    const MeshBase::const_element_iterator el_end   = d_mesh->active_local_elements_end();

    // Determine the element patches associated with each node N, which is
    // defined to be the collection of elements that contain node N.
    //
    // Unlike the standard Z-Z patch recovery algorithm, we use "tight" element
    // patches for non-vertex nodes.
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        for (unsigned int n = 0; n < elem->n_nodes(); ++n)
        {
            // Only set up patches for local nodes.
            const Node* const node = elem->get_node(n);
            if (node->processor_id() != rank) continue;

            // Only set up patches once for each node.
            const dof_id_type node_id = node->id();
            if (d_local_elem_patches.find(node_id) != d_local_elem_patches.end()) continue;

            // Find the elements that touch this node.
            std::set<const Elem*>& elem_patch = d_local_elem_patches[node_id];
            elem->find_point_neighbors(*node, elem_patch);
            d_local_elems.insert(d_local_elems.end(), elem_patch.begin(), elem_patch.end());
        }
    }

    // Determine the number of quadrature/interpolation points in each element
    // and setup mappings used to fill global data structures.
    //
    // We use full Gaussian quadrature rules (i.e. third-order Gauss quadrature
    // for first-order elements and fifth-order Gauss quadrature for
    // second-order elements) in all elements to avoid special treatment at
    // boundary nodes.
    d_n_qp_global = 0;
    d_n_qp_local  = 0;
    d_qp_global_offset = 0;
    const unsigned int n_elem = d_mesh->n_elem();
    d_elem_n_qp.resize(n_elem,0);
    d_elem_qp_global_offset.resize(n_elem,0);
    d_elem_qp_local_offset .resize(n_elem,0);
    AutoPtr<QBase> qrule;
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
        const Order quad_order = (elem_order == FIRST ? THIRD : FIFTH);
        bool reinit_qrule = false;
        if (!qrule.get() || qrule->get_dim() != dim || qrule->get_order() != quad_order)
        {
            qrule = QBase::build(QGAUSS, dim, quad_order);
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
        d_elem_sigma   [elem_id].resize(n_qp);
        d_elem_pressure[elem_id].resize(n_qp);
    }
    if (first_order_elems && second_order_elems)
    {
        TBOX_ERROR("cannot have both first- and second-order elements in the same mesh.\n");
    }
    std::vector<int> n_qp_per_proc(size);
    n_qp_per_proc[rank] = d_n_qp_local;
    Parallel::sum(n_qp_per_proc);
    d_qp_global_offset = std::accumulate(n_qp_per_proc.begin(), n_qp_per_proc.begin()+rank, 0);
    d_n_qp_global = std::accumulate(n_qp_per_proc.begin()+rank, n_qp_per_proc.end(), d_qp_global_offset);
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const dof_id_type elem_id = elem->id();
        d_elem_qp_global_offset[elem_id] = d_elem_qp_local_offset[elem_id]+d_qp_global_offset;
    }
    Parallel::sum(d_elem_qp_global_offset);
    Parallel::sum(d_elem_qp_local_offset);

    // Set up element patch L2 projection matrices.
    unsigned int dim = d_mesh->mesh_dimension();
    d_interp_order = first_order_elems ? FIRST : SECOND;
    d_quad_order = first_order_elems ? THIRD : FIFTH;
    const unsigned int num_basis = num_basis_functions(dim, d_interp_order);
    Eigen::MatrixXd M(num_basis, num_basis);
    Eigen::VectorXd P(num_basis);
    AutoPtr<FEBase> fe(FEBase::build(dim, FEType(d_interp_order, LAGRANGE)));
    const std::vector<libMesh::Point>& q_point = fe->get_xyz();
    qrule = QBase::build(QGAUSS, dim, d_quad_order);
    fe->attach_quadrature_rule(qrule.get());
    d_local_patch_proj_solver.resize(d_local_elem_patches.size());
    unsigned int k = 0;
    for (std::map<libMesh::dof_id_type,ElemPatch>::const_iterator it = d_local_elem_patches.begin();
         it != d_local_elem_patches.end(); ++it, ++k)
    {
        const dof_id_type node_id = it->first;
        const Node& node = d_mesh->node(node_id);
        const ElemPatch& elem_patch = it->second;
        M.setZero();
        for (ElemPatch::const_iterator el_it = elem_patch.begin(); el_it != elem_patch.end(); ++el_it)
        {
            const Elem* const elem = *el_it;
            fe->reinit(elem);
            for (unsigned int qp = 0; qp < qrule->n_points(); ++qp)
            {
                evaluate_basis_functions(P, node, q_point[qp], dim, d_interp_order);
                M += P*P.transpose();
            }
        }
        d_local_patch_proj_solver[k] = M.colPivHouseholderQr();
        if (!d_local_patch_proj_solver[k].isInvertible())
        {
            // TODO: We should try enlarging the element patch before emitting
            // an error message.
            TBOX_ERROR("encountered singular patch L2 projection matrix\n");
        }
    }
    return;
}// IBFEPatchRecoveryPostProcessor

IBFEPatchRecoveryPostProcessor::~IBFEPatchRecoveryPostProcessor()
{
    // intentionally blank
    return;
}// IBFEPatchRecoveryPostProcessor

void
IBFEPatchRecoveryPostProcessor::registerCauchyStressValue(
    const Elem* const elem,
    const QBase* const qrule,
    const unsigned int qp,
    const TensorValue<double>& sigma)
{
    if (elem->processor_id() != libMesh::processor_id())
    {
        TBOX_ERROR("must register stresses only for local elements\n");
    }
    TBOX_ASSERT(elem->default_order() == d_interp_order);
    TBOX_ASSERT(qrule->type() == QGAUSS);
    TBOX_ASSERT(qrule->get_order() == d_quad_order);
    TBOX_ASSERT(0 >= qp && qp < qrule->n_points());
    d_elem_sigma[elem->id()][qp] = sigma;
    return;
}// registerCauchyStressValue

void
IBFEPatchRecoveryPostProcessor::registerPressureValue(
    const Elem* const elem,
    const QBase* const qrule,
    const unsigned int qp,
    const double p)
{
    if (elem->processor_id() != libMesh::processor_id())
    {
        TBOX_ERROR("must register pressures only for local elements\n");
    }
    TBOX_ASSERT(elem->default_order() == d_interp_order);
    TBOX_ASSERT(qrule->type() == QGAUSS);
    TBOX_ASSERT(qrule->get_order() == d_quad_order);
    TBOX_ASSERT(0 >= qp && qp < qrule->n_points());
    d_elem_pressure[elem->id()][qp] = p;
    return;
}// registerPressureValue

void
IBFEPatchRecoveryPostProcessor::reconstructCauchyStress()
{
    // Communicate the stored values of the Cauchy stress.
    static const unsigned int NVARS = 2*(NDIM+1);
    std::vector<double> sigma_vals(NVARS*d_n_qp_global,0.0);
    const MeshBase::const_element_iterator el_begin = d_mesh->active_local_elements_begin();
    const MeshBase::const_element_iterator el_end   = d_mesh->active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const dof_id_type elem_id = elem->id();
        const int global_offset = d_elem_qp_global_offset[elem_id];
        for (int qp = 0; qp < d_elem_n_qp[elem_id]; ++qp)
        {
            const TensorValue<double>& stress = d_elem_sigma[elem_id][qp];
            for (unsigned int i = 0, k = 0; i < NDIM; ++i)
            {
                for (unsigned int j = i; j < NDIM; ++j, ++k)
                {
                    sigma_vals[NVARS*(global_offset+qp)+k] = stress(i,j);
                }
            }
        }
    }
    Parallel::sum(sigma_vals);

    // Perform element patch L2 projections.
    const unsigned int dim = d_mesh->mesh_dimension();
    const unsigned int num_basis = num_basis_functions(dim, d_interp_order);
    Eigen::VectorXd P(num_basis), a(num_basis), f(num_basis);
    AutoPtr<FEBase> fe(FEBase::build(dim, FEType(d_interp_order, LAGRANGE)));
    const std::vector<libMesh::Point>& q_point = fe->get_xyz();
    AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, d_quad_order);
    fe->attach_quadrature_rule(qrule.get());
    unsigned int k = 0;
    for (std::map<libMesh::dof_id_type,ElemPatch>::const_iterator it = d_local_elem_patches.begin();
         it != d_local_elem_patches.end(); ++it, ++k)
    {
        const dof_id_type node_id = it->first;
        const Node& node = d_mesh->node(node_id);
        const ElemPatch& elem_patch = it->second;
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& patch_proj_solver = d_local_patch_proj_solver[k];
        for (unsigned int var = 0; var < NVARS; ++var)
        {
            f.setZero();
            for (ElemPatch::const_iterator el_it = elem_patch.begin(); el_it != elem_patch.end(); ++el_it)
            {
                const Elem* const elem = *el_it;
                const dof_id_type elem_id = elem->id();
                const int global_offset = d_elem_qp_global_offset[elem_id];
                fe->reinit(elem);
                for (unsigned int qp = 0; qp < qrule->n_points(); ++qp)
                {
                    evaluate_basis_functions(P, node, q_point[qp], dim, d_interp_order);
                    f += P*sigma_vals[NVARS*(global_offset+qp)+var];
                }
            }
            a = patch_proj_solver.solve(f);
        }
    }
    return;
}// reconstructCauchyStress

void
IBFEPatchRecoveryPostProcessor::reconstructPressure()
{
    // Communicate the stored values of the Cauchy stress.
    std::vector<double> pressure_vals(d_n_qp_global,0.0);
    const MeshBase::const_element_iterator el_begin = d_mesh->active_local_elements_begin();
    const MeshBase::const_element_iterator el_end   = d_mesh->active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const dof_id_type elem_id = elem->id();
        const int global_offset = d_elem_qp_global_offset[elem_id];
        for (int qp = 0; qp < d_elem_n_qp[elem_id]; ++qp)
        {
            pressure_vals[global_offset+qp] = d_elem_pressure[elem_id][qp];
        }
    }
    Parallel::sum(pressure_vals);

    // Perform element patch L2 projections.
    const unsigned int dim = d_mesh->mesh_dimension();
    const unsigned int num_basis = num_basis_functions(dim, d_interp_order);
    Eigen::VectorXd P(num_basis), a(num_basis), f(num_basis);
    AutoPtr<FEBase> fe(FEBase::build(dim, FEType(d_interp_order, LAGRANGE)));
    const std::vector<libMesh::Point>& q_point = fe->get_xyz();
    AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, d_quad_order);
    fe->attach_quadrature_rule(qrule.get());
    unsigned int k = 0;
    for (std::map<libMesh::dof_id_type,ElemPatch>::const_iterator it = d_local_elem_patches.begin();
         it != d_local_elem_patches.end(); ++it, ++k)
    {
        const dof_id_type node_id = it->first;
        const Node& node = d_mesh->node(node_id);
        const ElemPatch& elem_patch = it->second;
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& patch_proj_solver = d_local_patch_proj_solver[k];
        f.setZero();
        for (ElemPatch::const_iterator el_it = elem_patch.begin(); el_it != elem_patch.end(); ++el_it)
        {
            const Elem* const elem = *el_it;
            const dof_id_type elem_id = elem->id();
            const int global_offset = d_elem_qp_global_offset[elem_id];
            fe->reinit(elem);
            for (unsigned int qp = 0; qp < qrule->n_points(); ++qp)
            {
                evaluate_basis_functions(P, node, q_point[qp], dim, d_interp_order);
                f += P*pressure_vals[global_offset+qp];
            }
        }
        a = patch_proj_solver.solve(f);
    }
    return;
}// reconstructPressure

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
