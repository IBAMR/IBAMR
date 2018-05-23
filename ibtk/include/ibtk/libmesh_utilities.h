// Filename: libmesh_utilities.h
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

#ifndef included_IBTK_libmesh_utilities
#define included_IBTK_libmesh_utilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "libmesh/dof_map.h"
#include "libmesh/dof_object.h"
#include "libmesh/edge.h"
#include "libmesh/face.h"
#include "libmesh/fe.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/point.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/type_tensor.h"
#include "libmesh/type_vector.h"
#include "libmesh/vector_value.h"
#include "tbox/Utilities.h"

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////

namespace IBTK
{
/*!
 * Struct allowing for the specification of system variables / gradients and the NumericVector used to evaluate
 * those quantities.
 */
struct SystemData
{
    SystemData(const std::string& system_name = "",
               const std::vector<int>& vars = std::vector<int>(),
               const std::vector<int>& grad_vars = std::vector<int>(),
               libMesh::NumericVector<double>* const system_vec = NULL)
        : system_name(system_name), vars(vars), grad_vars(grad_vars), system_vec(system_vec)
    {
    }

    std::string system_name;
    std::vector<int> vars, grad_vars;
    libMesh::NumericVector<double>* system_vec;
};

typedef void (*ScalarMeshFcnPtr)(
    double& F,
    const libMesh::TensorValue<double>& FF,
    const libMesh::Point& x,
    const libMesh::Point& X,
    libMesh::Elem* elem,
    const std::vector<const std::vector<double>*>& system_var_data,
    const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
    double data_time,
    void* ctx);

typedef void (*VectorMeshFcnPtr)(
    libMesh::VectorValue<double>& F,
    const libMesh::TensorValue<double>& FF,
    const libMesh::Point& x,
    const libMesh::Point& X,
    libMesh::Elem* elem,
    const std::vector<const std::vector<double>*>& system_var_data,
    const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
    double data_time,
    void* ctx);

typedef void (*TensorMeshFcnPtr)(
    libMesh::TensorValue<double>& F,
    const libMesh::TensorValue<double>& FF,
    const libMesh::Point& x,
    const libMesh::Point& X,
    libMesh::Elem* elem,
    const std::vector<const std::vector<double>*>& system_var_data,
    const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
    double data_time,
    void* ctx);

typedef void (*ScalarSurfaceFcnPtr)(
    double& F,
    const libMesh::VectorValue<double>& n,
    const libMesh::VectorValue<double>& N,
    const libMesh::TensorValue<double>& FF,
    const libMesh::Point& x,
    const libMesh::Point& X,
    libMesh::Elem* elem,
    unsigned short int side,
    const std::vector<const std::vector<double>*>& system_var_data,
    const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
    double data_time,
    void* ctx);

typedef void (*VectorSurfaceFcnPtr)(
    libMesh::VectorValue<double>& F,
    const libMesh::VectorValue<double>& n,
    const libMesh::VectorValue<double>& N,
    const libMesh::TensorValue<double>& FF,
    const libMesh::Point& x,
    const libMesh::Point& X,
    libMesh::Elem* elem,
    unsigned short int side,
    const std::vector<const std::vector<double>*>& system_var_data,
    const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
    double data_time,
    void* ctx);

typedef void (*TensorSurfaceFcnPtr)(
    libMesh::TensorValue<double>& F,
    const libMesh::VectorValue<double>& n,
    const libMesh::VectorValue<double>& N,
    const libMesh::TensorValue<double>& FF,
    const libMesh::Point& x,
    const libMesh::Point& X,
    libMesh::Elem* elem,
    unsigned short int side,
    const std::vector<const std::vector<double>*>& system_var_data,
    const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
    double data_time,
    void* ctx);

template <class MultiArray, class Array>
inline void
get_values_for_interpolation(MultiArray& U_node,
                             const libMesh::PetscVector<double>& U_petsc_vec,
                             const Array& U_local_soln,
                             const std::vector<unsigned int>& dof_indices)
{
    const std::size_t n_nodes = dof_indices.size();
    if (U_node.shape()[0] != n_nodes)
    {
        typename MultiArray::extent_gen extents;
        U_node.resize(extents[n_nodes]);
    }
    for (std::size_t k = 0; k < n_nodes; ++k)
    {
        U_node[k] = U_local_soln[U_petsc_vec.map_global_to_local_index(dof_indices[k])];
    }
    return;
} // get_values_for_interpolation

template <class MultiArray, class Array>
inline void
get_values_for_interpolation(MultiArray& U_node,
                             const libMesh::PetscVector<double>& U_petsc_vec,
                             const Array& U_local_soln,
                             const std::vector<std::vector<unsigned int> >& dof_indices)
{
    const std::size_t n_vars = dof_indices.size();
    const std::size_t n_nodes = dof_indices[0].size();
    if (U_node.shape()[0] != n_nodes || U_node.shape()[1] != n_vars)
    {
        typename MultiArray::extent_gen extents;
        U_node.resize(extents[n_nodes][n_vars]);
    }
    for (std::size_t k = 0; k < n_nodes; ++k)
    {
        for (std::size_t i = 0; i < n_vars; ++i)
        {
            U_node[k][i] = U_local_soln[U_petsc_vec.map_global_to_local_index(dof_indices[i][k])];
        }
    }
    return;
} // get_values_for_interpolation

template <class MultiArray>
inline void
get_values_for_interpolation(MultiArray& U_node,
                             libMesh::NumericVector<double>& U_vec,
                             const std::vector<unsigned int>& dof_indices)
{
    libMesh::PetscVector<double>* U_petsc_vec = dynamic_cast<libMesh::PetscVector<double>*>(&U_vec);
    Vec U_global_vec = U_petsc_vec->vec();
    Vec U_local_vec;
    VecGhostGetLocalForm(U_global_vec, &U_local_vec);
    double* U_local_soln;
    VecGetArray(U_local_vec, &U_local_soln);
    get_values_for_interpolation(U_node, *U_petsc_vec, U_local_soln, dof_indices);
    VecRestoreArray(U_local_vec, &U_local_soln);
    VecGhostRestoreLocalForm(U_global_vec, &U_local_vec);
    return;
} // get_values_for_interpolation

template <class MultiArray>
inline void
get_values_for_interpolation(MultiArray& U_node,
                             libMesh::NumericVector<double>& U_vec,
                             const std::vector<std::vector<unsigned int> >& dof_indices)
{
    libMesh::PetscVector<double>* U_petsc_vec = dynamic_cast<libMesh::PetscVector<double>*>(&U_vec);
    Vec U_global_vec = U_petsc_vec->vec();
    Vec U_local_vec;
    VecGhostGetLocalForm(U_global_vec, &U_local_vec);
    double* U_local_soln;
    VecGetArray(U_local_vec, &U_local_soln);
    get_values_for_interpolation(U_node, *U_petsc_vec, U_local_soln, dof_indices);
    VecRestoreArray(U_local_vec, &U_local_soln);
    VecGhostRestoreLocalForm(U_global_vec, &U_local_vec);
    return;
} // get_values_for_interpolation

template <class MultiArray>
inline void
interpolate(double& U, const int qp, const MultiArray& U_node, const std::vector<std::vector<double> >& phi)
{
    const int n_nodes = static_cast<int>(U_node.shape()[0]);
    U = 0.0;
    for (int k = 0; k < n_nodes; ++k)
    {
        U += U_node[k] * phi[k][qp];
    }
    return;
} // interpolate

template <class MultiArray>
inline double
interpolate(const int qp, const MultiArray& U_node, const std::vector<std::vector<double> >& phi)
{
    const int n_nodes = static_cast<int>(U_node.shape()[0]);
    double U = 0.0;
    for (int k = 0; k < n_nodes; ++k)
    {
        U += U_node[k] * phi[k][qp];
    }
    return U;
} // interpolate

template <class MultiArray>
inline void
interpolate(double* const U, const int qp, const MultiArray& U_node, const std::vector<std::vector<double> >& phi)
{
    const int n_nodes = static_cast<int>(U_node.shape()[0]);
    const int n_vars = static_cast<int>(U_node.shape()[1]);
    std::fill(U, U + n_vars, 0.0);
    for (int k = 0; k < n_nodes; ++k)
    {
        const double& p = phi[k][qp];
        for (int i = 0; i < n_vars; ++i)
        {
            U[i] += U_node[k][i] * p;
        }
    }
    return;
} // interpolate

template <class MultiArray>
inline void
interpolate(libMesh::TypeVector<double>& U,
            const int qp,
            const MultiArray& U_node,
            const std::vector<std::vector<double> >& phi)
{
    const int n_nodes = static_cast<int>(U_node.shape()[0]);
    const int n_vars = static_cast<int>(U_node.shape()[1]);
    U.zero();
    for (int k = 0; k < n_nodes; ++k)
    {
        const double& p = phi[k][qp];
        for (int i = 0; i < n_vars; ++i)
        {
            U(i) += U_node[k][i] * p;
        }
    }
    return;
} // interpolate

template <class MultiArray>
inline void
jacobian(libMesh::TypeTensor<double>& dX_ds,
         const int qp,
         const MultiArray& X_node,
         const std::vector<std::vector<libMesh::VectorValue<double> > >& dphi)
{
    const int n_nodes = static_cast<int>(X_node.shape()[0]);
    const int dim = static_cast<int>(X_node.shape()[1]);
    dX_ds.zero();
    for (int k = 0; k < n_nodes; ++k)
    {
        const libMesh::VectorValue<double>& dphi_ds = dphi[k][qp];
        for (int i = 0; i < dim; ++i)
        {
            const double& X = X_node[k][i];
            for (int j = 0; j < dim; ++j)
            {
                dX_ds(i, j) += X * dphi_ds(j);
            }
        }
    }
    if (dim == 2)
    {
        dX_ds(2, 2) = 1.0;
    }
    return;
} // jacobian

inline void
tensor_inverse(libMesh::TensorValue<double>& A_inv, const libMesh::TensorValue<double>& A, const int dim = NDIM)
{
    const double det_A = A.det();
    if (dim == 2)
    {
        A_inv(0, 0) = +A(1, 1) / det_A;
        A_inv(0, 1) = -A(0, 1) / det_A;
        A_inv(0, 2) = 0.0;
        A_inv(1, 0) = -A(1, 0) / det_A;
        A_inv(1, 1) = +A(0, 0) / det_A;
        A_inv(1, 2) = 0.0;
        A_inv(2, 0) = 0.0;
        A_inv(2, 1) = 0.0;
        A_inv(2, 2) = 1.0;
    }
    else
    {
        A_inv(0, 0) = +(A(2, 2) * A(1, 1) - A(2, 1) * A(1, 2)) / det_A;
        A_inv(0, 1) = -(A(2, 2) * A(0, 1) - A(2, 1) * A(0, 2)) / det_A;
        A_inv(0, 2) = +(A(1, 2) * A(0, 1) - A(1, 1) * A(0, 2)) / det_A;
        A_inv(1, 0) = -(A(2, 2) * A(1, 0) - A(2, 0) * A(1, 2)) / det_A;
        A_inv(1, 1) = +(A(2, 2) * A(0, 0) - A(2, 0) * A(0, 2)) / det_A;
        A_inv(1, 2) = -(A(1, 2) * A(0, 0) - A(1, 0) * A(0, 2)) / det_A;
        A_inv(2, 0) = +(A(2, 1) * A(1, 0) - A(2, 0) * A(1, 1)) / det_A;
        A_inv(2, 1) = -(A(2, 1) * A(0, 0) - A(2, 0) * A(0, 1)) / det_A;
        A_inv(2, 2) = +(A(1, 1) * A(0, 0) - A(1, 0) * A(0, 1)) / det_A;
    }
    return;
} // tensor_inverse

inline libMesh::TensorValue<double>
tensor_inverse(const libMesh::TensorValue<double>& A, const int dim = NDIM)
{
    libMesh::TensorValue<double> A_inv;
    const double det_A = A.det();
    if (dim == 2)
    {
        A_inv(0, 0) = +A(1, 1) / det_A;
        A_inv(0, 1) = -A(0, 1) / det_A;
        A_inv(0, 2) = 0.0;
        A_inv(1, 0) = -A(1, 0) / det_A;
        A_inv(1, 1) = +A(0, 0) / det_A;
        A_inv(1, 2) = 0.0;
        A_inv(2, 0) = 0.0;
        A_inv(2, 1) = 0.0;
        A_inv(2, 2) = 1.0;
    }
    else
    {
        A_inv(0, 0) = +(A(2, 2) * A(1, 1) - A(2, 1) * A(1, 2)) / det_A;
        A_inv(0, 1) = -(A(2, 2) * A(0, 1) - A(2, 1) * A(0, 2)) / det_A;
        A_inv(0, 2) = +(A(1, 2) * A(0, 1) - A(1, 1) * A(0, 2)) / det_A;
        A_inv(1, 0) = -(A(2, 2) * A(1, 0) - A(2, 0) * A(1, 2)) / det_A;
        A_inv(1, 1) = +(A(2, 2) * A(0, 0) - A(2, 0) * A(0, 2)) / det_A;
        A_inv(1, 2) = -(A(1, 2) * A(0, 0) - A(1, 0) * A(0, 2)) / det_A;
        A_inv(2, 0) = +(A(2, 1) * A(1, 0) - A(2, 0) * A(1, 1)) / det_A;
        A_inv(2, 1) = -(A(2, 1) * A(0, 0) - A(2, 0) * A(0, 1)) / det_A;
        A_inv(2, 2) = +(A(1, 1) * A(0, 0) - A(1, 0) * A(0, 1)) / det_A;
    }
    return A_inv;
} // tensor_inverse

inline void
tensor_inverse_transpose(libMesh::TensorValue<double>& A_inv_trans,
                         const libMesh::TensorValue<double>& A,
                         const int dim = NDIM)
{
    const double det_A = A.det();
    if (dim == 2)
    {
        A_inv_trans(0, 0) = +A(1, 1) / det_A;
        A_inv_trans(0, 1) = -A(1, 0) / det_A;
        A_inv_trans(0, 2) = 0.0;
        A_inv_trans(1, 0) = -A(0, 1) / det_A;
        A_inv_trans(1, 1) = +A(0, 0) / det_A;
        A_inv_trans(1, 2) = 0.0;
        A_inv_trans(2, 0) = 0.0;
        A_inv_trans(2, 1) = 0.0;
        A_inv_trans(2, 2) = 1.0;
    }
    else
    {
        A_inv_trans(0, 0) = +(A(2, 2) * A(1, 1) - A(2, 1) * A(1, 2)) / det_A;
        A_inv_trans(0, 1) = -(A(2, 2) * A(1, 0) - A(2, 0) * A(1, 2)) / det_A;
        A_inv_trans(0, 2) = +(A(2, 1) * A(1, 0) - A(2, 0) * A(1, 1)) / det_A;
        A_inv_trans(1, 0) = -(A(2, 2) * A(0, 1) - A(2, 1) * A(0, 2)) / det_A;
        A_inv_trans(1, 1) = +(A(2, 2) * A(0, 0) - A(2, 0) * A(0, 2)) / det_A;
        A_inv_trans(1, 2) = -(A(2, 1) * A(0, 0) - A(2, 0) * A(0, 1)) / det_A;
        A_inv_trans(2, 0) = +(A(1, 2) * A(0, 1) - A(1, 1) * A(0, 2)) / det_A;
        A_inv_trans(2, 1) = -(A(1, 2) * A(0, 0) - A(1, 0) * A(0, 2)) / det_A;
        A_inv_trans(2, 2) = +(A(1, 1) * A(0, 0) - A(1, 0) * A(0, 1)) / det_A;
    }
    return;
} // tensor_inverse_transpose

inline libMesh::TensorValue<double>
tensor_inverse_transpose(const libMesh::TensorValue<double>& A, const int dim = NDIM)
{
    libMesh::TensorValue<double> A_inv_trans;
    const double det_A = A.det();
    if (dim == 2)
    {
        A_inv_trans(0, 0) = +A(1, 1) / det_A;
        A_inv_trans(0, 1) = -A(1, 0) / det_A;
        A_inv_trans(0, 2) = 0.0;
        A_inv_trans(1, 0) = -A(0, 1) / det_A;
        A_inv_trans(1, 1) = +A(0, 0) / det_A;
        A_inv_trans(1, 2) = 0.0;
        A_inv_trans(2, 0) = 0.0;
        A_inv_trans(2, 1) = 0.0;
        A_inv_trans(2, 2) = 1.0;
    }
    else
    {
        A_inv_trans(0, 0) = +(A(2, 2) * A(1, 1) - A(2, 1) * A(1, 2)) / det_A;
        A_inv_trans(0, 1) = -(A(2, 2) * A(1, 0) - A(2, 0) * A(1, 2)) / det_A;
        A_inv_trans(0, 2) = +(A(2, 1) * A(1, 0) - A(2, 0) * A(1, 1)) / det_A;
        A_inv_trans(1, 0) = -(A(2, 2) * A(0, 1) - A(2, 1) * A(0, 2)) / det_A;
        A_inv_trans(1, 1) = +(A(2, 2) * A(0, 0) - A(2, 0) * A(0, 2)) / det_A;
        A_inv_trans(1, 2) = -(A(2, 1) * A(0, 0) - A(2, 0) * A(0, 1)) / det_A;
        A_inv_trans(2, 0) = +(A(1, 2) * A(0, 1) - A(1, 1) * A(0, 2)) / det_A;
        A_inv_trans(2, 1) = -(A(1, 2) * A(0, 0) - A(1, 0) * A(0, 2)) / det_A;
        A_inv_trans(2, 2) = +(A(1, 1) * A(0, 0) - A(1, 0) * A(0, 1)) / det_A;
    }
    return A_inv_trans;
} // tensor_inverse_transpose

inline void
outer_product(libMesh::TensorValue<double>& u_prod_v,
              const libMesh::TypeVector<double>& u,
              const libMesh::TypeVector<double>& v,
              const int dim = NDIM)
{
    u_prod_v.zero();
    for (int i = 0; i < dim; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            u_prod_v(i, j) = u(i) * v(j);
        }
    }
    return;
} // outer_product

inline libMesh::TensorValue<double>
outer_product(const libMesh::TypeVector<double>& u, const libMesh::TypeVector<double>& v, const int dim = NDIM)
{
    libMesh::TensorValue<double> u_prod_v;
    for (int i = 0; i < dim; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            u_prod_v(i, j) = u(i) * v(j);
        }
    }
    return u_prod_v;
} // outer_product

// WARNING: This code is specialized to the case in which q is a unit vector
// aligned with the coordinate axes.
inline void
intersect_line_with_edge(std::vector<std::pair<double, libMesh::Point> >& t_vals,
                         libMesh::Edge* elem,
                         libMesh::Point r,
                         libMesh::VectorValue<double> q)
{
    t_vals.resize(0);
    switch (elem->type())
    {
    case libMesh::EDGE2:
    {
        // Linear interpolation:
        //
        //    0.5*(1-u)*p0 + 0.5*(1+u)*p1 = r + t*q
        //
        // Factor the interpolation formula:
        //
        //    0.5*(p1-p0)*u+0.5*(p1+p0) = r + t*q
        //
        // Solve for u:
        //
        //    a*u + b = 0
        //
        // with:
        //
        //    a = 0.5*(-p0+p1)
        //    b = 0.5*(p0+p1) - r
        const libMesh::Point& p0 = *elem->get_node(0);
        const libMesh::Point& p1 = *elem->get_node(1);
        double a, b;
        if (q(0) != 0.0)
        {
            a = 0.5 * (p1(1) - p0(1));
            b = 0.5 * (p1(1) + p0(1)) - r(1);
        }
        else
        {
            a = 0.5 * (p1(0) - p0(0));
            b = 0.5 * (p1(0) + p0(0)) - r(0);
        }
        const double u = -b / a;

        // Look for intersections within the element interior.
        if (u >= -1.0 && u <= 1.0)
        {
            double t;
            if (std::abs(q(0)) >= std::abs(q(1)))
            {
                const double p = p0(0) * 0.5 * (1.0 - u) + p1(0) * 0.5 * (1.0 + u);
                t = (p - r(0)) / q(0);
            }
            else
            {
                const double p = p0(1) * 0.5 * (1.0 - u) + p1(1) * 0.5 * (1.0 + u);
                t = (p - r(1)) / q(1);
            }
            t_vals.push_back(std::make_pair(t, libMesh::Point(u, 0.0, 0.0)));
        }
        break;
    }
    case libMesh::EDGE3:
    {
        // Quadratic interpolation:
        //
        //    0.5*u*(u-1)*p0 + 0.5*u*(u+1)*p1 + (1-u*u)*p2 = r + t*q
        //
        // Factor the interpolation formula:
        //
        //    (0.5*p0+0.5*p1-p2)*u^2 + 0.5*(p1-p0)*u + p2 = r + t*q
        //
        // Solve for u:
        //
        //    a*u^2 + b*u + c = 0
        //
        // with:
        //
        //    a = (0.5*p0+0.5*p1-p2)
        //    b = 0.5*(p1-p0)
        //    c = p2-r
        const libMesh::Point& p0 = *elem->get_node(0);
        const libMesh::Point& p1 = *elem->get_node(1);
        const libMesh::Point& p2 = *elem->get_node(2);
        double a, b, c;
        if (q(0) != 0.0)
        {
            a = (0.5 * p0(1) + 0.5 * p1(1) - p2(1));
            b = 0.5 * (p1(1) - p0(1));
            c = p2(1) - r(1);
        }
        else
        {
            a = (0.5 * p0(0) + 0.5 * p1(0) - p2(0));
            b = 0.5 * (p1(0) - p0(0));
            c = p2(0) - r(0);
        }
        const double disc = b * b - 4.0 * a * c;
        std::vector<double> u_vals;
        if (disc > 0.0)
        {
            const double q = -0.5 * (b + (b > 0.0 ? 1.0 : -1.0) * sqrt(disc));
            const double u0 = q / a;
            u_vals.push_back(u0);
            const double u1 = c / q;
            if (std::abs(u0 - u1) > std::numeric_limits<double>::epsilon())
            {
                u_vals.push_back(u1);
            }
        }

        // Look for intersections within the element interior.
        for (unsigned int k = 0; k < u_vals.size(); ++k)
        {
            double u = u_vals[k];
            if (u >= -1.0 && u <= 1.0)
            {
                double t;
                if (std::abs(q(0)) >= std::abs(q(1)))
                {
                    const double p = 0.5 * u * (u - 1.0) * p0(0) + 0.5 * u * (u + 1.0) * p1(0) + (1.0 - u * u) * p2(0);
                    t = (p - r(0)) / q(0);
                }
                else
                {
                    const double p = 0.5 * u * (u - 1.0) * p0(1) + 0.5 * u * (u + 1.0) * p1(1) + (1.0 - u * u) * p2(1);
                    t = (p - r(1)) / q(1);
                }
                t_vals.push_back(std::make_pair(t, libMesh::Point(u, 0.0, 0.0)));
            }
        }
        break;
    }
    default:
    {
        TBOX_ERROR("intersect_line_with_edge():"
                   << "  element type "
                   << libMesh::Utility::enum_to_string<libMesh::ElemType>(elem->type())
                   << " is not supported at this time.\n");
    }
    }
    return;
} // intersect_line_with_edge

// WARNING: This code is specialized to the case in which q is a unit vector
// aligned with the coordinate axes.
inline void
intersect_line_with_face(std::vector<std::pair<double, libMesh::Point> >& t_vals,
                         libMesh::Face* elem,
                         libMesh::Point r,
                         libMesh::VectorValue<double> q)
{
    t_vals.resize(0);
    switch (elem->type())
    {
    case libMesh::TRI3:
    {
        // Linear interpolation:
        //
        //    (1-u-v)*p0 + u*p1 + v*p2 = r + t*q
        //
        // Factor the interpolation formula:
        //
        //    (p1-p0)*u + (p2-p0)*v + p0 = r + t*q
        //
        // Solve a small linear system for u and v.
        const libMesh::Point& p0 = *elem->get_node(0);
        const libMesh::Point& p1 = *elem->get_node(1);
        const libMesh::Point& p2 = *elem->get_node(2);
        double A00, A10, A01, A11, C1, C2;
        if (q(0) != 0.0)
        {
            A00 = p1(1) - p0(1);
            A01 = p2(1) - p0(1);
            C1 = p0(1) - r(1);
            A10 = p1(2) - p0(2);
            A11 = p2(2) - p0(2);
            C2 = p0(2) - r(2);
        }
        else if (q(1) != 0.0)
        {
            A00 = p1(0) - p0(0);
            A01 = p2(0) - p0(0);
            C1 = p0(0) - r(0);
            A10 = p1(2) - p0(2);
            A11 = p2(2) - p0(2);
            C2 = p0(2) - r(2);
        }
        else
        {
            A00 = p1(0) - p0(0);
            A01 = p2(0) - p0(0);
            C1 = p0(0) - r(0);
            A10 = p1(1) - p0(1);
            A11 = p2(1) - p0(1);
            C2 = p0(1) - r(1);
        }
        const double det = A00 * A11 - A10 * A01;
        if (std::abs(det) > std::numeric_limits<double>::epsilon())
        {
            const double u = (A01 * C2 - A11 * C1) / det;
            const double v = -(A00 * C2 - A10 * C1) / det;

            // Look for intersections within the element interior.
            if (u >= 0.0 && v >= 0.0 && (u + v) <= 1.0)
            {
                double t;
                if (std::abs(q(0)) >= std::abs(q(1)) && std::abs(q(0)) >= std::abs(q(2)))
                {
                    const double p = u * p0(0) + v * p1(0) + (1.0 - u - v) * p2(0);
                    t = (p - r(0)) / q(0);
                }
                else if (std::abs(q(1)) >= std::abs(q(2)))
                {
                    const double p = u * p0(1) + v * p1(1) + (1.0 - u - v) * p2(1);
                    t = (p - r(1)) / q(1);
                }
                else
                {
                    const double p = u * p0(2) + v * p1(2) + (1.0 - u - v) * p2(2);
                    t = (p - r(2)) / q(2);
                }
                t_vals.push_back(std::make_pair(t, libMesh::Point(u, v, 0.0)));
            }
        }
        break;
    }
    case libMesh::QUAD4:
    {
        const libMesh::Point& p00 = *elem->get_node(0);
        const libMesh::Point& p10 = *elem->get_node(1);
        const libMesh::Point& p11 = *elem->get_node(2);
        const libMesh::Point& p01 = *elem->get_node(3);

        const libMesh::Point a = p11 - p10 - p01 + p00;
        const libMesh::Point b = p10 - p00;
        const libMesh::Point c = p01 - p00;
        const libMesh::Point d = p00;

        double A1, A2, B1, B2, C1, C2, D1, D2;
        if (q(0) != 0.0)
        {
            A1 = a(1);
            A2 = a(2);
            B1 = b(1);
            B2 = b(2);
            C1 = c(1);
            C2 = c(2);
            D1 = d(1) - r(1);
            D2 = d(2) - r(2);
        }
        else if (q(1) != 0.0)
        {
            A1 = a(0);
            A2 = a(2);
            B1 = b(0);
            B2 = b(2);
            C1 = c(0);
            C2 = c(2);
            D1 = d(0) - r(0);
            D2 = d(2) - r(2);
        }
        else
        {
            A1 = a(0);
            A2 = a(1);
            B1 = b(0);
            B2 = b(1);
            C1 = c(0);
            C2 = c(1);
            D1 = d(0) - r(0);
            D2 = d(1) - r(1);
        }

        // (A2*C1 - A1*C2) v^2 + (A2*D1 - A1*D2 + B2*C1 - B1*C2) v + (B2*D1 - B1*D2) = 0
        std::vector<double> v_vals;
        {
            const double a = A2 * C1 - A1 * C2;
            const double b = A2 * D1 - A1 * D2 + B2 * C1 - B1 * C2;
            const double c = B2 * D1 - B1 * D2;
            const double disc = b * b - 4.0 * a * c;
            if (disc > 0.0)
            {
                const double q = -0.5 * (b + (b > 0.0 ? 1.0 : -1.0) * sqrt(disc));
                const double v0 = q / a;
                v_vals.push_back(v0);
                const double v1 = c / q;
                if (std::abs(v0 - v1) > std::numeric_limits<double>::epsilon())
                {
                    v_vals.push_back(v1);
                }
            }
        }

        for (unsigned int k = 0; k < v_vals.size(); ++k)
        {
            double v = v_vals[k];
            if (v >= 0.0 && v <= 1.0)
            {
                double u;
                const double a = v * A2 + B2;
                const double b = v * (A2 - A1) + B2 - B1;
                if (std::abs(b) >= std::abs(a))
                {
                    u = (v * (C1 - C2) + D1 - D2) / b;
                }
                else
                {
                    u = (-v * C2 - D2) / a;
                }

                if (u >= 0.0 && u <= 1.0)
                {
                    double t;
                    if (std::abs(q(0)) >= std::abs(q(1)) && std::abs(q(0)) >= std::abs(q(2)))
                    {
                        const double p = p00(0) * (1.0 - u) * (1.0 - v) + p01(0) * (1.0 - u) * v +
                                         p10(0) * u * (1.0 - v) + p11(0) * u * v;
                        t = (p - r(0)) / q(0);
                    }
                    else if (std::abs(q(1)) >= std::abs(q(2)))
                    {
                        const double p = p00(1) * (1.0 - u) * (1.0 - v) + p01(1) * (1.0 - u) * v +
                                         p10(1) * u * (1.0 - v) + p11(1) * u * v;
                        t = (p - r(1)) / q(1);
                    }
                    else
                    {
                        const double p = p00(2) * (1.0 - u) * (1.0 - v) + p01(2) * (1.0 - u) * v +
                                         p10(2) * u * (1.0 - v) + p11(2) * u * v;
                        t = (p - r(2)) / q(2);
                    }
                    t_vals.push_back(std::make_pair(t, libMesh::Point(2.0 * u - 1.0, 2.0 * v - 1.0, 0.0)));
                }
            }
        }
        break;
    }
    default:
    {
        TBOX_ERROR("intersect_line_with_face():"
                   << "  element type "
                   << libMesh::Utility::enum_to_string<libMesh::ElemType>(elem->type())
                   << " is not supported at this time.\n");
    }
    }
    return;
} // intersect_line_with_face

struct DofObjectComp : std::binary_function<const libMesh::DofObject* const, const libMesh::DofObject* const, bool>
{
    inline bool operator()(const libMesh::DofObject* const lhs, const libMesh::DofObject* const rhs)
    {
        return lhs->id() < rhs->id();
    }
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_libmesh_utilities
