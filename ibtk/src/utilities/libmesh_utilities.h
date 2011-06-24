// Filename: libmesh_utilities.h
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

#ifndef included_libmesh_utilities
#define included_libmesh_utilities

/////////////////////////////// INCLUDES /////////////////////////////////////

// BLITZ++ INCLUDES
#include <blitz/array.h>

// LIBMESH INCLUDES
#define LIBMESH_REQUIRE_SEPARATE_NAMESPACE
#include <dof_object.h>
#include <petsc_vector.h>
#include <point.h>
#include <type_tensor.h>
#include <type_vector.h>
#include <vector_value.h>

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////

namespace IBTK
{
inline void
get_values_for_interpolation(
    blitz::Array<double,1>& U_node,
    libMesh::NumericVector<double>& U_vec,
    const std::vector<unsigned int>& dof_indices)
{
    const int n_nodes = dof_indices.size();
    if (U_node.extent(0) != n_nodes) U_node.resize(n_nodes);
    libMesh::PetscVector<double>* U_petsc_vec = dynamic_cast<libMesh::PetscVector<double>*>(&U_vec);
    Vec U_global_vec = U_petsc_vec->vec();
    Vec U_local_vec;
    VecGhostGetLocalForm(U_global_vec,&U_local_vec);
    double* values;
    VecGetArray(U_local_vec, &values);
    for (int k = 0; k < n_nodes; ++k)
    {
        const unsigned int local_index = U_petsc_vec->map_global_to_local_index(dof_indices[k]);
        U_node(k) = values[local_index];
    }
    VecRestoreArray(U_local_vec, &values);
    VecGhostRestoreLocalForm(U_global_vec,&U_local_vec);
    return;
}// get_values_for_interpolation

inline void
get_values_for_interpolation(
    blitz::Array<double,2>& U_node,
    libMesh::NumericVector<double>& U_vec,
    const blitz::Array<std::vector<unsigned int>,1>& dof_indices)
{
    const int n_vars = dof_indices.extent(0);
    const int n_nodes = dof_indices(0).size();
    if (U_node.extent(0) != n_nodes || U_node.extent(1) != n_vars) U_node.resize(n_nodes,n_vars);
    libMesh::PetscVector<double>* U_petsc_vec = dynamic_cast<libMesh::PetscVector<double>*>(&U_vec);
    Vec U_global_vec = U_petsc_vec->vec();
    Vec U_local_vec;
    VecGhostGetLocalForm(U_global_vec,&U_local_vec);
    double* values;
    VecGetArray(U_local_vec, &values);
    for (int k = 0; k < n_nodes; ++k)
    {
        for (int i = 0; i < n_vars; ++i)
        {
            const unsigned int local_index = U_petsc_vec->map_global_to_local_index(dof_indices(i)[k]);
            U_node(k,i) = values[local_index];
        }
    }
    VecRestoreArray(U_local_vec, &values);
    VecGhostRestoreLocalForm(U_global_vec,&U_local_vec);
    return;
}// get_values_for_interpolation

inline void
interpolate(
    double& U,
    const int qp,
    const blitz::Array<double,1>& U_node,
    const blitz::Array<double,2>& phi)
{
    const int n_nodes = U_node.extent(blitz::firstDim);
    U = 0.0;
    for (int k = 0; k < n_nodes; ++k)
    {
        U += U_node(k)*phi(qp,k);
    }
    return;
}// interpolate

inline void
interpolate(
    double& U,
    const int qp,
    const blitz::Array<double,1>& U_node,
    const std::vector<std::vector<double> >& phi)
{
    const int n_nodes = U_node.extent(blitz::firstDim);
    U = 0.0;
    for (int k = 0; k < n_nodes; ++k)
    {
        U += U_node(k)*phi[k][qp];
    }
    return;
}// interpolate

inline double
interpolate(
    const int qp,
    const blitz::Array<double,1>& U_node,
    const blitz::Array<double,2>& phi)
{
    const int n_nodes = U_node.extent(blitz::firstDim);
    double U = 0.0;
    for (int k = 0; k < n_nodes; ++k)
    {
        U += U_node(k)*phi(qp,k);
    }
    return U;
}// interpolate

inline double
interpolate(
    const int qp,
    const blitz::Array<double,1>& U_node,
    const std::vector<std::vector<double> >& phi)
{
    const int n_nodes = U_node.extent(blitz::firstDim);
    double U = 0.0;
    for (int k = 0; k < n_nodes; ++k)
    {
        U += U_node(k)*phi[k][qp];
    }
    return U;
}// interpolate

inline void
interpolate(
    double* const U,
    const int qp,
    const blitz::Array<double,2>& U_node,
    const blitz::Array<double,2>& phi)
{
    const int n_nodes = U_node.extent(blitz::firstDim);
    const int n_vars  = U_node.extent(blitz::secondDim);
    std::fill(U, U+n_vars, 0.0);
    for (int k = 0; k < n_nodes; ++k)
    {
        const double& p = phi(qp,k);
        for (int i = 0; i < n_vars; ++i)
        {
            U[i] += U_node(k,i)*p;
        }
    }
    return;
}// interpolate

inline void
interpolate(
    double* const U,
    const int qp,
    const blitz::Array<double,2>& U_node,
    const std::vector<std::vector<double> >& phi)
{
    const int n_nodes = U_node.extent(blitz::firstDim);
    const int n_vars  = U_node.extent(blitz::secondDim);
    std::fill(U, U+n_vars, 0.0);
    for (int k = 0; k < n_nodes; ++k)
    {
        const double& p = phi[k][qp];
        for (int i = 0; i < n_vars; ++i)
        {
            U[i] += U_node(k,i)*p;
        }
    }
    return;
}// interpolate

inline void
interpolate(
    libMesh::TypeVector<double>& U,
    const int qp,
    const blitz::Array<double,2>& U_node,
    const blitz::Array<double,2>& phi)
{
    const int n_nodes = U_node.extent(blitz::firstDim);
    const int n_vars  = U_node.extent(blitz::secondDim);
    U.zero();
    for (int k = 0; k < n_nodes; ++k)
    {
        const double& p = phi(qp,k);
        for (int i = 0; i < n_vars; ++i)
        {
            U(i) += U_node(k,i)*p;
        }
    }
    return;
}// interpolate

inline void
interpolate(
    libMesh::TypeVector<double>& U,
    const int qp,
    const blitz::Array<double,2>& U_node,
    const std::vector<std::vector<double> >& phi)
{
    const int n_nodes = U_node.extent(blitz::firstDim);
    const int n_vars  = U_node.extent(blitz::secondDim);
    U.zero();
    for (int k = 0; k < n_nodes; ++k)
    {
        const double& p = phi[k][qp];
        for (int i = 0; i < n_vars; ++i)
        {
            U(i) += U_node(k,i)*p;
        }
    }
    return;
}// interpolate

inline void
jacobian(
    libMesh::TypeTensor<double>& dX_ds,
    const int qp,
    const blitz::Array<double,2>& X_node,
    const blitz::Array<libMesh::VectorValue<double>,2>& dphi)
{
    const int n_nodes = X_node.extent(blitz::firstDim);
    const int dim     = X_node.extent(blitz::secondDim);
    dX_ds.zero();
    for (int k = 0; k < n_nodes; ++k)
    {
        const libMesh::VectorValue<double>& dphi_ds = dphi(qp,k);
        for (int i = 0; i < dim; ++i)
        {
            const double& X = X_node(k,i);
            for (int j = 0; j < dim; ++j)
            {
                dX_ds(i,j) += X*dphi_ds(j);
            }
        }
    }
    if (dim == 2)
    {
        dX_ds(2,2) = 1.0;
    }
    return;
}// jacobian

inline void
jacobian(
    libMesh::TypeTensor<double>& dX_ds,
    const int qp,
    const blitz::Array<double,2>& X_node,
    const std::vector<std::vector<libMesh::VectorValue<double> > >& dphi)
{
    const int n_nodes = X_node.extent(blitz::firstDim);
    const int dim     = X_node.extent(blitz::secondDim);
    dX_ds.zero();
    for (int k = 0; k < n_nodes; ++k)
    {
        const libMesh::VectorValue<double>& dphi_ds = dphi[k][qp];
        for (int i = 0; i < dim; ++i)
        {
            const double& X = X_node(k,i);
            for (int j = 0; j < dim; ++j)
            {
                dX_ds(i,j) += X*dphi_ds(j);
            }
        }
    }
    if (dim == 2)
    {
        dX_ds(2,2) = 1.0;
    }
    return;
}// jacobian

inline void
jacobian(
    libMesh::TypeTensor<double>& dX_ds,
    const int qp,
    const blitz::Array<double,2>& X_node,
    const blitz::Array<libMesh::VectorValue<double>,2>& dphi,
    const blitz::Array<double,1>* const strain_J_bar_node,
    const blitz::Array<double,2>* const strain_phi)
{
    jacobian(dX_ds, qp, X_node, dphi);
    if (strain_J_bar_node == NULL || strain_phi == NULL) return;
    const int dim = X_node.extent(blitz::secondDim);
    const double J = dX_ds.det();
    const double J_bar = interpolate(qp,*strain_J_bar_node,*strain_phi);
    const double alpha = pow(J_bar/J,1.0/static_cast<double>(dim));
    dX_ds *= alpha;
    if (dim == 2) dX_ds(2,2) = 1.0;
    return;
}// jacobian

inline void
jacobian(
    libMesh::TypeTensor<double>& dX_ds,
    const int qp,
    const blitz::Array<double,2>& X_node,
    const std::vector<std::vector<libMesh::VectorValue<double> > >& dphi,
    const blitz::Array<double,1>* const strain_J_bar_node,
    const std::vector<std::vector<double> >* const strain_phi)
{
    jacobian(dX_ds, qp, X_node, dphi);
    if (strain_J_bar_node == NULL || strain_phi == NULL) return;
    const int dim = X_node.extent(blitz::secondDim);
    const double J = dX_ds.det();
    const double J_bar = interpolate(qp,*strain_J_bar_node,*strain_phi);
    const double alpha = pow(J_bar/J,1.0/static_cast<double>(dim));
    dX_ds *= alpha;
    if (dim == 2) dX_ds(2,2) = 1.0;
    return;
}// jacobian

inline void
tensor_inverse(
    libMesh::TensorValue<double>& A_inv,
    const libMesh::TensorValue<double>& A,
    const int dim)
{
    const double det_A = A.det();
    if (dim == 2)
    {
        A_inv(0,0) =  A(1,1)/det_A;
        A_inv(0,1) = -A(0,1)/det_A;
        A_inv(0,2) = 0.0;
        A_inv(1,0) = -A(1,0)/det_A;
        A_inv(1,1) =  A(0,0)/det_A;
        A_inv(1,2) = 0.0;
        A_inv(2,0) = 0.0;
        A_inv(2,1) = 0.0;
        A_inv(2,2) = 1.0;
    }
    else
    {
        A_inv(0,0) =  (A(2,2)*A(1,1)-A(2,1)*A(1,2))/det_A;
        A_inv(0,1) = -(A(2,2)*A(0,1)-A(2,1)*A(0,2))/det_A;
        A_inv(0,2) =  (A(1,2)*A(0,1)-A(1,1)*A(0,2))/det_A;
        A_inv(1,0) = -(A(2,2)*A(1,0)-A(2,0)*A(1,2))/det_A;
        A_inv(1,1) =  (A(2,2)*A(0,0)-A(2,0)*A(0,2))/det_A;
        A_inv(1,2) = -(A(1,2)*A(0,0)-A(1,0)*A(0,2))/det_A;
        A_inv(2,0) =  (A(2,1)*A(1,0)-A(2,0)*A(1,1))/det_A;
        A_inv(2,1) = -(A(2,1)*A(0,0)-A(2,0)*A(0,1))/det_A;
        A_inv(2,2) =  (A(1,1)*A(0,0)-A(1,0)*A(0,1))/det_A;
    }
    return;
}// tensor_inverse

inline libMesh::TensorValue<double>
tensor_inverse(
    const libMesh::TensorValue<double>& A,
    const int dim)
{
    libMesh::TensorValue<double> A_inv;
    const double det_A = A.det();
    if (dim == 2)
    {
        A_inv(0,0) =  A(1,1)/det_A;
        A_inv(0,1) = -A(0,1)/det_A;
        A_inv(0,2) = 0.0;
        A_inv(1,0) = -A(1,0)/det_A;
        A_inv(1,1) =  A(0,0)/det_A;
        A_inv(1,2) = 0.0;
        A_inv(2,0) = 0.0;
        A_inv(2,1) = 0.0;
        A_inv(2,2) = 1.0;
    }
    else
    {
        A_inv(0,0) =  (A(2,2)*A(1,1)-A(2,1)*A(1,2))/det_A;
        A_inv(0,1) = -(A(2,2)*A(0,1)-A(2,1)*A(0,2))/det_A;
        A_inv(0,2) =  (A(1,2)*A(0,1)-A(1,1)*A(0,2))/det_A;
        A_inv(1,0) = -(A(2,2)*A(1,0)-A(2,0)*A(1,2))/det_A;
        A_inv(1,1) =  (A(2,2)*A(0,0)-A(2,0)*A(0,2))/det_A;
        A_inv(1,2) = -(A(1,2)*A(0,0)-A(1,0)*A(0,2))/det_A;
        A_inv(2,0) =  (A(2,1)*A(1,0)-A(2,0)*A(1,1))/det_A;
        A_inv(2,1) = -(A(2,1)*A(0,0)-A(2,0)*A(0,1))/det_A;
        A_inv(2,2) =  (A(1,1)*A(0,0)-A(1,0)*A(0,1))/det_A;
    }
    return A_inv;
}// tensor_inverse

inline void
tensor_inverse_transpose(
    libMesh::TensorValue<double>& A_inv_trans,
    const libMesh::TensorValue<double>& A,
    const int dim)
{
    const double det_A = A.det();
    if (dim == 2)
    {
        A_inv_trans(0,0) =  A(1,1)/det_A;
        A_inv_trans(0,1) = -A(1,0)/det_A;
        A_inv_trans(0,2) = 0.0;
        A_inv_trans(1,0) = -A(0,1)/det_A;
        A_inv_trans(1,1) =  A(0,0)/det_A;
        A_inv_trans(1,2) = 0.0;
        A_inv_trans(2,0) = 0.0;
        A_inv_trans(2,1) = 0.0;
        A_inv_trans(2,2) = 1.0;
    }
    else
    {
        A_inv_trans(0,0) =  (A(2,2)*A(1,1)-A(2,1)*A(1,2))/det_A;
        A_inv_trans(0,1) = -(A(2,2)*A(1,0)-A(2,0)*A(1,2))/det_A;
        A_inv_trans(0,2) =  (A(2,1)*A(1,0)-A(2,0)*A(1,1))/det_A;
        A_inv_trans(1,0) = -(A(2,2)*A(0,1)-A(2,1)*A(0,2))/det_A;
        A_inv_trans(1,1) =  (A(2,2)*A(0,0)-A(2,0)*A(0,2))/det_A;
        A_inv_trans(1,2) = -(A(2,1)*A(0,0)-A(2,0)*A(0,1))/det_A;
        A_inv_trans(2,0) =  (A(1,2)*A(0,1)-A(1,1)*A(0,2))/det_A;
        A_inv_trans(2,1) = -(A(1,2)*A(0,0)-A(1,0)*A(0,2))/det_A;
        A_inv_trans(2,2) =  (A(1,1)*A(0,0)-A(1,0)*A(0,1))/det_A;
    }
    return;
}// tensor_inverse_transpose

inline libMesh::TensorValue<double>
tensor_inverse_transpose(
    const libMesh::TensorValue<double>& A,
    const int dim)
{
    libMesh::TensorValue<double> A_inv_trans;
    const double det_A = A.det();
    if (dim == 2)
    {
        A_inv_trans(0,0) =  A(1,1)/det_A;
        A_inv_trans(0,1) = -A(1,0)/det_A;
        A_inv_trans(0,2) = 0.0;
        A_inv_trans(1,0) = -A(0,1)/det_A;
        A_inv_trans(1,1) =  A(0,0)/det_A;
        A_inv_trans(1,2) = 0.0;
        A_inv_trans(2,0) = 0.0;
        A_inv_trans(2,1) = 0.0;
        A_inv_trans(2,2) = 1.0;
    }
    else
    {
        A_inv_trans(0,0) =  (A(2,2)*A(1,1)-A(2,1)*A(1,2))/det_A;
        A_inv_trans(0,1) = -(A(2,2)*A(1,0)-A(2,0)*A(1,2))/det_A;
        A_inv_trans(0,2) =  (A(2,1)*A(1,0)-A(2,0)*A(1,1))/det_A;
        A_inv_trans(1,0) = -(A(2,2)*A(0,1)-A(2,1)*A(0,2))/det_A;
        A_inv_trans(1,1) =  (A(2,2)*A(0,0)-A(2,0)*A(0,2))/det_A;
        A_inv_trans(1,2) = -(A(2,1)*A(0,0)-A(2,0)*A(0,1))/det_A;
        A_inv_trans(2,0) =  (A(1,2)*A(0,1)-A(1,1)*A(0,2))/det_A;
        A_inv_trans(2,1) = -(A(1,2)*A(0,0)-A(1,0)*A(0,2))/det_A;
        A_inv_trans(2,2) =  (A(1,1)*A(0,0)-A(1,0)*A(0,1))/det_A;
    }
    return A_inv_trans;
}// tensor_inverse_transpose

inline void
outer_product(
    libMesh::TensorValue<double>& u_prod_v,
    const libMesh::TypeVector<double>& u,
    const libMesh::TypeVector<double>& v)
{
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    {
        for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
        {
            u_prod_v(i,j) = u(i)*v(j);
        }
    }
    return;
}// outer_product

inline libMesh::TensorValue<double>
outer_product(
    const libMesh::TypeVector<double>& u,
    const libMesh::TypeVector<double>& v)
{
    libMesh::TensorValue<double> u_prod_v;
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    {
        for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
        {
            u_prod_v(i,j) = u(i)*v(j);
        }
    }
    return u_prod_v;
}// outer_product

struct DofObjectComp
    : std::binary_function<const libMesh::DofObject* const,const libMesh::DofObject* const,bool>
{
    inline bool
    operator()(
        const libMesh::DofObject* const lhs,
        const libMesh::DofObject* const rhs)
        {
            return lhs->id() < rhs->id();
        }
};
}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_libmesh_utilities
