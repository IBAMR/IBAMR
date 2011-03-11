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

// C++ STDLIB INCLUDES
#include <vector>

// BLITZ++ INCLUDES
#include <blitz/array.h>

// LIBMESH INCLUDES
#define LIBMESH_REQUIRE_SEPARATE_NAMESPACE
#include <dof_object.h>
#include <numeric_vector.h>
#include <point.h>
#include <tensor_value.h>
#include <vector_value.h>

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////

namespace IBTK
{
inline double
compute_interpolation(
    const unsigned int qp,
    const libMesh::NumericVector<double>& U_vec,
    const std::vector<std::vector<double> >& phi,
    const std::vector<unsigned int>& dof_indices)
{
    double U = 0.0;
    for (unsigned int k = 0; k < phi.size(); ++k)
    {
        U += U_vec(dof_indices[k])*phi[k][qp];
    }
    return U;
}// compute_interpolation

inline double
compute_interpolation(
    const int qp,
    const libMesh::NumericVector<double>& U_vec,
    const blitz::Array<double,2>& phi,
    const std::vector<unsigned int>& dof_indices)
{
    double U = 0.0;
    for (int k = 0; k < phi.extent(blitz::secondDim); ++k)
    {
        U += U_vec(dof_indices[k])*phi(qp,k);
    }
    return U;
}// compute_interpolation

inline std::vector<double>
compute_interpolation(
    const unsigned int qp,
    const libMesh::NumericVector<double>& U_vec,
    const std::vector<std::vector<double> >& phi,
    const std::vector<std::vector<unsigned int> >& dof_indices)
{
    const unsigned int n_vars = dof_indices.size();
    std::vector<double> U(n_vars,0.0);
    for (unsigned int k = 0; k < phi.size(); ++k)
    {
        for (unsigned int i = 0; i < n_vars; ++i)
        {
            U[i] += U_vec(dof_indices[i][k])*phi[k][qp];
        }
    }
    return U;
}// compute_interpolation

inline blitz::Array<double,1>
compute_interpolation(
    const int qp,
    const libMesh::NumericVector<double>& U_vec,
    const blitz::Array<double,2>& phi,
    const blitz::Array<std::vector<unsigned int>,1>& dof_indices)
{
    const unsigned int n_vars = dof_indices.size();
    blitz::Array<double,1> U(n_vars);
    U = 0.0;
    for (int k = 0; k < phi.extent(blitz::secondDim); ++k)
    {
        for (unsigned int i = 0; i < n_vars; ++i)
        {
            U(i) += U_vec(dof_indices(i)[k])*phi(qp,k);
        }
    }
    return U;
}// compute_interpolation

inline libMesh::Point
compute_coordinate(
    const unsigned int qp,
    const libMesh::NumericVector<double>& X,
    const std::vector<std::vector<double> >& phi,
    const std::vector<std::vector<unsigned int> >& dof_indices)
{
    const unsigned int dim = dof_indices.size();
    libMesh::Point X_qp;
    for (unsigned int k = 0; k < phi.size(); ++k)
    {
        for (unsigned int i = 0; i < dim; ++i)
        {
            X_qp(i) += X(dof_indices[i][k])*phi[k][qp];
        }
    }
    return X_qp;
}// compute_coordinate

inline libMesh::Point
compute_coordinate(
    const int qp,
    const libMesh::NumericVector<double>& X,
    const blitz::Array<double,2>& phi,
    const blitz::Array<std::vector<unsigned int>,1>& dof_indices)
{
    const unsigned int dim = dof_indices.size();
    libMesh::Point X_qp;
    for (int k = 0; k < phi.extent(blitz::secondDim); ++k)
    {
        for (unsigned int i = 0; i < dim; ++i)
        {
            X_qp(i) += X(dof_indices(i)[k])*phi(qp,k);
        }
    }
    return X_qp;
}// compute_coordinate

inline libMesh::TensorValue<double>
compute_coordinate_mapping_jacobian(
    const unsigned int qp,
    const libMesh::NumericVector<double>& X,
    const std::vector<std::vector<libMesh::VectorValue<double> > >& dphi,
    const std::vector<std::vector<unsigned int> >& dof_indices)
{
    const unsigned int dim = dof_indices.size();
    libMesh::TensorValue<double> dX_ds;
    for (unsigned int k = 0; k < dphi.size(); ++k)
    {
        const libMesh::VectorValue<double>& dphi_ds = dphi[k][qp];
        for (unsigned int j = 0; j < dim; ++j)
        {
            for (unsigned int i = 0; i < dim; ++i)
            {
                dX_ds(i,j) += X(dof_indices[i][k])*dphi_ds(j);
            }
        }
    }
    if (dim == 2)
    {
        dX_ds(2,2) = 1.0;
    }
    return dX_ds;
}// compute_coordinate_mapping_jacobian

inline libMesh::TensorValue<double>
compute_coordinate_mapping_jacobian(
    const unsigned int qp,
    const libMesh::NumericVector<double>& X,
    const std::vector<std::vector<libMesh::VectorValue<double> > >& dphi,
    const std::vector<std::vector<unsigned int> >& dof_indices,
    const libMesh::NumericVector<double>* const proj_strain_J_bar,
    const std::vector<std::vector<double> >* const proj_strain_phi,
    const std::vector<unsigned int>* const proj_strain_dof_indices)
{
    libMesh::TensorValue<double> dX_ds = compute_coordinate_mapping_jacobian(qp, X, dphi, dof_indices);
    if (proj_strain_J_bar == NULL || proj_strain_phi == NULL || proj_strain_dof_indices == NULL) return dX_ds;
    const unsigned int dim = dof_indices.size();
    const double J = dX_ds.det();
    const double J_bar = compute_interpolation(qp,*proj_strain_J_bar,*proj_strain_phi,*proj_strain_dof_indices);
    const double alpha = pow(J_bar/J,1.0/double(NDIM));
    dX_ds *= alpha;
    if (dim == 2) dX_ds(2,2) = 1.0;
    return dX_ds;
}// compute_coordinate_mapping_jacobian

inline libMesh::TensorValue<double>
compute_coordinate_mapping_jacobian(
    const int qp,
    const libMesh::NumericVector<double>& X,
    const blitz::Array<libMesh::VectorValue<double>,2>& dphi,
    const blitz::Array<std::vector<unsigned int>,1>& dof_indices)
{
    const unsigned int dim = dof_indices.size();
    libMesh::TensorValue<double> dX_ds;
    for (int k = 0; k < dphi.extent(blitz::secondDim); ++k)
    {
        const libMesh::VectorValue<double>& dphi_ds = dphi(qp,k);
        for (unsigned int j = 0; j < dim; ++j)
        {
            for (unsigned int i = 0; i < dim; ++i)
            {
                dX_ds(i,j) += X(dof_indices(i)[k])*dphi_ds(j);
            }
        }
    }
    if (dim == 2)
    {
        dX_ds(2,2) = 1.0;
    }
    return dX_ds;
}// compute_coordinate_mapping_jacobian

inline libMesh::TensorValue<double>
compute_coordinate_mapping_jacobian(
    const int qp,
    const libMesh::NumericVector<double>& X,
    const blitz::Array<libMesh::VectorValue<double>,2>& dphi,
    const blitz::Array<std::vector<unsigned int>,1>& dof_indices,
    const libMesh::NumericVector<double>* const proj_strain_J_bar,
    const blitz::Array<double,2>* const proj_strain_phi,
    const std::vector<unsigned int>* const proj_strain_dof_indices)
{
    libMesh::TensorValue<double> dX_ds = compute_coordinate_mapping_jacobian(qp, X, dphi, dof_indices);
    const unsigned int dim = dof_indices.size();
    if (proj_strain_J_bar == NULL || proj_strain_phi == NULL || proj_strain_dof_indices == NULL) return dX_ds;
    const double J = dX_ds.det();
    const double J_bar = compute_interpolation(qp,*proj_strain_J_bar,*proj_strain_phi,*proj_strain_dof_indices);
    const double alpha = pow(J_bar/J,1.0/double(NDIM));
    dX_ds *= alpha;
    if (dim == 2) dX_ds(2,2) = 1.0;
    return dX_ds;
}// compute_coordinate_mapping_jacobian

inline double
compute_coordinate_mapping_jacobian_det(
    const unsigned int qp,
    const libMesh::NumericVector<double>& X,
    const std::vector<std::vector<libMesh::VectorValue<double> > >& dphi,
    const std::vector<std::vector<unsigned int> >& dof_indices)
{
    libMesh::TensorValue<double> dX_ds = compute_coordinate_mapping_jacobian(qp,X,dphi,dof_indices);
    return dX_ds.det();
}// compute_coordinate_mapping_jacobian_det

inline double
compute_coordinate_mapping_jacobian_det(
    const unsigned int qp,
    const libMesh::NumericVector<double>& X,
    const std::vector<std::vector<libMesh::VectorValue<double> > >& dphi,
    const std::vector<std::vector<unsigned int> >& dof_indices,
    const libMesh::NumericVector<double>* const proj_strain_J_bar,
    const std::vector<std::vector<double> >* const proj_strain_phi,
    const std::vector<unsigned int>* const proj_strain_dof_indices)
{
    libMesh::TensorValue<double> dX_ds = compute_coordinate_mapping_jacobian(qp, X, dphi, dof_indices, proj_strain_J_bar, proj_strain_phi, proj_strain_dof_indices);
    return dX_ds.det();
}// compute_coordinate_mapping_jacobian_det

inline double
compute_coordinate_mapping_jacobian_det(
    const int qp,
    const libMesh::NumericVector<double>& X,
    const blitz::Array<libMesh::VectorValue<double>,2>& dphi,
    const blitz::Array<std::vector<unsigned int>,1>& dof_indices)
{
    libMesh::TensorValue<double> dX_ds = compute_coordinate_mapping_jacobian(qp, X, dphi, dof_indices);
    return dX_ds.det();
}// compute_coordinate_mapping_jacobian_det

inline double
compute_coordinate_mapping_jacobian_det(
    const int qp,
    const libMesh::NumericVector<double>& X,
    const blitz::Array<libMesh::VectorValue<double>,2>& dphi,
    const blitz::Array<std::vector<unsigned int>,1>& dof_indices,
    const libMesh::NumericVector<double>* const proj_strain_J_bar,
    const blitz::Array<double,2>* const proj_strain_phi,
    const std::vector<unsigned int>* const proj_strain_dof_indices)
{
    libMesh::TensorValue<double> dX_ds = compute_coordinate_mapping_jacobian(qp, X, dphi, dof_indices, proj_strain_J_bar, proj_strain_phi, proj_strain_dof_indices);
    return dX_ds.det();
}// compute_coordinate_mapping_jacobian_det

inline libMesh::TensorValue<double>
tensor_inverse_transpose(
    const libMesh::TensorValue<double>& A,
    const int dim)
{
    const double det_A = A.det();
    libMesh::TensorValue<double> A_inv_trans = 0.0;
    if (dim == 2)
    {
        A_inv_trans(0,0) =  A(1,1)/det_A;
        A_inv_trans(0,1) = -A(1,0)/det_A;
        A_inv_trans(1,0) = -A(0,1)/det_A;
        A_inv_trans(1,1) =  A(0,0)/det_A;
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

inline libMesh::TensorValue<double>
outer_product(
    const libMesh::VectorValue<double>& u,
    const libMesh::VectorValue<double>& v)
{
    const libMesh::TensorValue<double> u_prod_v(
        u(0)*v(0) , u(0)*v(1) , u(0)*v(2) ,
        u(1)*v(0) , u(1)*v(1) , u(1)*v(2) ,
        u(2)*v(0) , u(2)*v(1) , u(2)*v(2)
                                                );
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
