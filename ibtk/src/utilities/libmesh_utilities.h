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
#include <dof_map.h>
#include <dof_object.h>
#include <edge.h>
#include <face.h>
#include <fe.h>
#include <petsc_vector.h>
#include <point.h>
#include <quadrature_gauss.h>
#include <string_to_enum.h>
#include <type_tensor.h>
#include <type_vector.h>
#include <vector_value.h>

// SAMRAI INCLUDES
#include <tbox/SAMRAI_MPI.h>
#include <tbox/Utilities.h>

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////

namespace IBTK
{
class QAdaptiveGauss
    : public libMesh::QBase
{
public:
    inline
    QAdaptiveGauss(
        const unsigned int dim,
        const double point_density)
        : libMesh::QBase(dim, FIRST),
          d_point_density(point_density)
        {
            for (unsigned int n = 1; n <= 22; ++n)
            {
                d_q1d[n] = new libMesh::QGauss(1, static_cast<Order>(2*n-1));
                d_q1d[n]->init(EDGE2);

                d_qtri2d[n] = new libMesh::QGauss(2, static_cast<Order>(2*n-1));
                d_qtri2d[n]->init(TRI3);

                d_qtet3d[n] = new libMesh::QGauss(3, static_cast<Order>(2*n-1));
                d_qtet3d[n]->init(TET4);
            }
            return;
        }

    inline
    ~QAdaptiveGauss()
        {
            for (unsigned int n = 1; n <= 22; ++n)
            {
                delete d_q1d[n];
                delete d_qtri2d[n];
                delete d_qtet3d[n];
            }
        }

    inline QuadratureType
    type() const
        { return QGAUSS; }

    inline bool
    shapes_need_reinit()
        { return true; }

    inline void
    set_elem_data(
        const ElemType type,
        const blitz::Array<double,2>& X_node,
        const double* const dx)
        {
            int n_nodes;
            switch (type)
            {
                case EDGE2:
                case EDGE3:
                case EDGE4:
                    n_nodes = 2;
                    break;
                case TRI3:
                case TRI6:
                    n_nodes = 3;
                    break;
                case QUAD4:
                case QUAD8:
                case QUAD9:
                    n_nodes = 4;
                    break;
                case TET4:
                case TET10:
                    n_nodes = 4;
                    break;
                case HEX8:
                case HEX20:
                case HEX27:
                    n_nodes = 8;
                    break;
                default:
                    n_nodes = 0;
            }

            libMesh::Point elem_X[8];
            for (int k = 0; k < n_nodes; ++k)
            {
                for (int d = 0; d < NDIM; ++d)
                {
                    elem_X[k](d) = X_node(k,d);
                }
            }

            const double dx_min = *std::min_element(dx,dx+NDIM);

            int min_points;
            switch (type)
            {
                case EDGE2:
                    min_points = 2;
                    break;
                case EDGE3:
                    min_points = 3;
                    break;
                case EDGE4:
                    min_points = 4;
                    break;
                case TRI3:
                case QUAD4:
                    min_points = 2;
                    break;
                case TRI6:
                case QUAD8:
                case QUAD9:
                    min_points = 3;
                    break;
                case TET4:
                case HEX8:
                    min_points = 2;
                    break;
                case TET10:
                case HEX20:
                case HEX27:
                    min_points = 3;
                    break;
                default:
                    min_points = 0;
            }

            switch (type)
            {
                case EDGE2:
                case EDGE3:
                case EDGE4:
                {
                    const double l_max = (elem_X[1] - elem_X[0]).size();
                    const int n = std::min(std::max(min_points,static_cast<int>(std::ceil(d_point_density*l_max/dx_min))),22);
                    const libMesh::QGauss& q = *d_q1d[n];
                    _points  = q.get_points();
                    _weights = q.get_weights();
                    break;
                }
                case TRI3:
                case TRI6:
                {
                    const double l_max = std::max((elem_X[1] - elem_X[0]).size(),
                                                  std::max((elem_X[2] - elem_X[0]).size(),
                                                           (elem_X[2] - elem_X[1]).size()));
                    const int n = std::min(std::max(min_points,static_cast<int>(std::ceil(d_point_density*l_max/dx_min))),22);
                    const libMesh::QGauss& q = *d_qtri2d[n];
                    _points  = q.get_points();
                    _weights = q.get_weights();
                    break;
                }
                case QUAD4:
                case QUAD8:
                case QUAD9:
                {
                    const double l_max0 = std::max((elem_X[1] - elem_X[0]).size(),
                                                   (elem_X[2] - elem_X[3]).size());
                    const int n0 = std::min(std::max(min_points,static_cast<int>(std::ceil(d_point_density*l_max0/dx_min))),22);
                    const libMesh::QGauss& q0 = *d_q1d[n0];

                    const double l_max1 = std::max((elem_X[3] - elem_X[0]).size(),
                                                   (elem_X[2] - elem_X[1]).size());
                    const int n1 = std::min(std::max(min_points,static_cast<int>(std::ceil(d_point_density*l_max1/dx_min))),22);
                    const libMesh::QGauss& q1 = *d_q1d[n1];

                    const unsigned int n_points0 = q0.n_points();
                    const unsigned int n_points1 = q1.n_points();

                    _points .resize(n_points0 * n_points1);
                    _weights.resize(n_points0 * n_points1);

                    const std::vector<libMesh::Point>& points0 = q0.get_points();
                    const std::vector<libMesh::Point>& points1 = q1.get_points();

                    const std::vector<libMesh::Real>& weights0 = q0.get_weights();
                    const std::vector<libMesh::Real>& weights1 = q1.get_weights();

                    unsigned int qp = 0;
                    for (unsigned int i1 = 0; i1 < n_points1; ++i1)
                    {
                        for (unsigned int i0 = 0; i0 < n_points0; ++i0, ++qp)
                        {
                            _points [qp](0) = points0[i0](0);
                            _points [qp](1) = points1[i1](0);
                            _weights[qp]    = weights0[i0] * weights1[i1];
                        }
                    }
                    break;
                }
                case TET4:
                case TET10:
                {
                    const double l_max = std::max(std::max((elem_X[1] - elem_X[0]).size(),
                                                           (elem_X[2] - elem_X[0]).size()),
                                                  std::max(std::max((elem_X[3] - elem_X[0]).size(),
                                                                    (elem_X[2] - elem_X[1]).size()),
                                                           std::max((elem_X[3] - elem_X[1]).size(),
                                                                    (elem_X[3] - elem_X[2]).size())));
                    const int n = std::min(std::max(min_points,static_cast<int>(std::ceil(d_point_density*l_max/dx_min))),22);
                    const libMesh::QGauss& q = *d_qtet3d[n];
                    _points  = q.get_points();
                    _weights = q.get_weights();
                    break;
                }
                case HEX8:
                case HEX20:
                case HEX27:
                {
                    const double l_max0 = std::max(std::max((elem_X[1] - elem_X[0]).size(),
                                                            (elem_X[2] - elem_X[3]).size()),
                                                   std::max((elem_X[5] - elem_X[4]).size(),
                                                            (elem_X[6] - elem_X[7]).size()));
                    const int n0 = std::min(std::max(min_points,static_cast<int>(std::ceil(d_point_density*l_max0/dx_min))),22);
                    const libMesh::QGauss& q0 = *d_q1d[n0];

                    const double l_max1 = std::max(std::max((elem_X[3] - elem_X[0]).size(),
                                                            (elem_X[2] - elem_X[1]).size()),
                                                   std::max((elem_X[7] - elem_X[4]).size(),
                                                            (elem_X[6] - elem_X[5]).size()));
                    const int n1 = std::min(std::max(min_points,static_cast<int>(std::ceil(d_point_density*l_max1/dx_min))),22);
                    const libMesh::QGauss& q1 = *d_q1d[n1];

                    const double l_max2 = std::max(std::max((elem_X[4] - elem_X[0]).size(),
                                                            (elem_X[5] - elem_X[1]).size()),
                                                   std::max((elem_X[6] - elem_X[2]).size(),
                                                            (elem_X[7] - elem_X[3]).size()));
                    const int n2 = std::min(std::max(min_points,static_cast<int>(std::ceil(d_point_density*l_max2/dx_min))),22);
                    const libMesh::QGauss& q2 = *d_q1d[n2];

                    const unsigned int n_points0 = q0.n_points();
                    const unsigned int n_points1 = q1.n_points();
                    const unsigned int n_points2 = q2.n_points();

                    _points .resize(n_points0 * n_points1 * n_points2);
                    _weights.resize(n_points0 * n_points1 * n_points2);

                    const std::vector<libMesh::Point>& points0 = q0.get_points();
                    const std::vector<libMesh::Point>& points1 = q1.get_points();
                    const std::vector<libMesh::Point>& points2 = q2.get_points();

                    const std::vector<libMesh::Real>& weights0 = q0.get_weights();
                    const std::vector<libMesh::Real>& weights1 = q1.get_weights();
                    const std::vector<libMesh::Real>& weights2 = q2.get_weights();

                    unsigned int qp = 0;
                    for (unsigned int i2 = 0; i2 < n_points2; ++i2)
                    {
                        for (unsigned int i1 = 0; i1 < n_points1; ++i1)
                        {
                            for (unsigned int i0 = 0; i0 < n_points0; ++i0, ++qp)
                            {
                                _points [qp](0) = points0[i0](0);
                                _points [qp](1) = points1[i1](0);
                                _points [qp](2) = points2[i2](0);
                                _weights[qp]    = weights0[i0] * weights1[i1] * weights2[i2];
                            }
                        }
                    }
                    break;
                }
                default:
                    TBOX_ERROR("unsupported\n");
            }
            return;
        }

private:
    inline void
    init_1D(
        const ElemType /*type*/,
        unsigned int /*p_level=0*/)
        {
            // intentionally blank
            return;
        }

    inline void
    init_2D(
        const ElemType /*type*/,
        unsigned int /*p_level=0*/)
        {
            // intentionally blank
            return;
        }

    inline void
    init_3D(
        const ElemType /*type*/,
        unsigned int /*p_level=0*/)
        {
            // intentionally blank
            return;
        }

    libMesh::QGauss*    d_q1d[23];
    libMesh::QGauss* d_qtri2d[23];
    libMesh::QGauss* d_qtet3d[23];
    const double d_point_density;
};

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
    VecGhostRestoreLocalForm(U_global_vec, &U_local_vec);
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
    VecGhostRestoreLocalForm(U_global_vec, &U_local_vec);
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
    const blitz::Array<double,1>& strain_J_bar_node,
    const blitz::Array<double,2>& strain_phi)
{
    jacobian(dX_ds, qp, X_node, dphi);
    const int dim = X_node.extent(blitz::secondDim);
    const double J = dX_ds.det();
    const double J_bar = interpolate(qp,strain_J_bar_node,strain_phi);
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
    const blitz::Array<double,1>& strain_J_bar_node,
    const std::vector<std::vector<double> >& strain_phi)
{
    jacobian(dX_ds, qp, X_node, dphi);
    const int dim = X_node.extent(blitz::secondDim);
    const double J = dX_ds.det();
    const double J_bar = interpolate(qp,strain_J_bar_node,strain_phi);
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

inline void
make_incompressible_tensor(
    libMesh::TensorValue<double>& A,
    const int dim)
{
    const double det_A = A.det();
    const double alpha = pow(1.0/det_A,1.0/static_cast<double>(dim));
    A *= alpha;
    if (dim == 2)
    {
        A(0,2) = 0.0;
        A(1,2) = 0.0;
        A(2,0) = 0.0;
        A(2,1) = 0.0;
        A(2,2) = 1.0;
    }
    return;
}// make_incompressible_tensor

inline libMesh::TensorValue<double>
make_incompressible_tensor(
    const libMesh::TensorValue<double>& A,
    const int dim)
{
    libMesh::TensorValue<double> A_incompressible = A;
    const double det_A = A_incompressible.det();
    const double alpha = pow(1.0/det_A,1.0/static_cast<double>(dim));
    A_incompressible *= alpha;
    if (dim == 2)
    {
        A_incompressible(0,2) = 0.0;
        A_incompressible(1,2) = 0.0;
        A_incompressible(2,0) = 0.0;
        A_incompressible(2,1) = 0.0;
        A_incompressible(2,2) = 1.0;
    }
    return A_incompressible;
}// make_incompressible_tensor

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

// WARNING: This code is specialized to the case in which q is a unit vector
// aligned with the coordinate axes.
inline std::vector<std::pair<double,libMesh::Point> >
intersect_line_with_edge(
    libMesh::Edge* elem,
    libMesh::Point r,
    libMesh::VectorValue<double> q)
{
    std::vector<std::pair<double,libMesh::Point> > t_vals;
    switch (elem->type())
    {
        case libMeshEnums::EDGE2:
        {
            // Linear interpolation:
            //
            //    0.5*(1-u)*p0 + 0.5*(1+u)*p1 = r + t * q
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
                a = 0.5*(p1(1)-p0(1));
                b = 0.5*(p1(1)+p0(1))-r(1);
            }
            else
            {
                a = 0.5*(p1(0)-p0(0));
                b = 0.5*(p1(0)+p0(0))-r(0);
            }
            const double u = -b/a;

            // Look for intersections within the element interior.
            if (u >= -1.0 && u <= 1.0)
            {
                double t;
                if (std::abs(q(0)) >= std::abs(q(1)))
                {
                    const double p = p0(0)*0.5*(1.0-u) + p1(0)*0.5*(1.0+u);
                    t = (p-r(0))/q(0);
                }
                else
                {
                    const double p = p0(1)*0.5*(1.0-u) + p1(1)*0.5*(1.0+u);
                    t = (p-r(1))/q(1);
                }
                t_vals.push_back(std::make_pair(t,libMesh::Point(u,0.0,0.0)));
            }
            break;
        }
        case libMeshEnums::EDGE3:
        {
            // Quadratic interpolation:
            //
            //    0.5*u*(u-1)*p0 + 0.5*u*(u+1)*p1 + (1-u*u)*p2 = r + t * q
            //
            // Factor the interpolation formula:
            //
            //    (0.5*p0+0.5*p1-p2)*u^2 + 0.5*(p1-p0)*u + p2 = r + t * q
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
                a = (0.5*p0(1)+0.5*p1(1)-p2(1));
                b = 0.5*(p1(1)-p0(1));
                c = p2(1)-r(1);
            }
            else
            {
                a = (0.5*p0(0)+0.5*p1(0)-p2(0));
                b = 0.5*(p1(0)-p0(0));
                c = p2(0)-r(0);
            }
            const double disc = b*b - 4.0*a*c;
            std::vector<double> u_vals;
            if (disc > 0.0)
            {
                const double q = -0.5*(b+(b>0.0 ? 1.0 : -1.0)*sqrt(disc));
                const double u0 = q/a;
                u_vals.push_back(u0);
                const double u1 = c/q;
                if (std::abs(u0-u1) > std::numeric_limits<double>::epsilon())
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
                        const double p = 0.5*u*(u-1.0)*p0(0) + 0.5*u*(u+1.0)*p1(0) + (1.0-u*u)*p2(0);
                        t = (p-r(0))/q(0);
                    }
                    else
                    {
                        const double p = 0.5*u*(u-1.0)*p0(1) + 0.5*u*(u+1.0)*p1(1) + (1.0-u*u)*p2(1);
                        t = (p-r(1))/q(1);
                    }
                    t_vals.push_back(std::make_pair(t,libMesh::Point(u,0.0,0.0)));
                }
            }
            break;
        }
        default:
        {
            TBOX_ERROR("intersect_line_with_edge():"
                       << "  element type " << libMesh::Utility::enum_to_string<libMeshEnums::ElemType>(elem->type()) << " is not supported at this time.\n");
        }
    }
    return t_vals;
}// intersect_line_with_edge

// WARNING: This code is specialized to the case in which q is a unit vector
// aligned with the coordinate axes.
inline std::vector<std::pair<double,libMesh::Point> >
intersect_line_with_face(
    libMesh::Face* elem,
    libMesh::Point r,
    libMesh::VectorValue<double> q)
{
    std::vector<std::pair<double,libMesh::Point> > t_vals;
    switch (elem->type())
    {
        case libMeshEnums::TRI3:
        {
            // Linear interpolation:
            //
            //    (1-u-v)*p0 + u*p1 + v*p2 = r + t * q
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
                A00 = p1(1)-p0(1);
                A01 = p2(1)-p0(1);
                C1  = p0(1)- r(1);
                A10 = p1(2)-p0(2);
                A11 = p2(2)-p0(2);
                C2  = p0(2)- r(2);
            }
            else if (q(1) != 0.0)
            {
                A00 = p1(0)-p0(0);
                A01 = p2(0)-p0(0);
                C1  = p0(0)- r(0);
                A10 = p1(2)-p0(2);
                A11 = p2(2)-p0(2);
                C2  = p0(2)- r(2);
            }
            else
            {
                A00 = p1(0)-p0(0);
                A01 = p2(0)-p0(0);
                C1  = p0(0)- r(0);
                A10 = p1(1)-p0(1);
                A11 = p2(1)-p0(1);
                C2  = p0(1)- r(1);
            }
            const double det = A00*A11-A10*A01;
            if (std::abs(det) > std::numeric_limits<double>::epsilon())
            {
                const double u =  (A01*C2-A11*C1)/det;
                const double v = -(A00*C2-A10*C1)/det;

                // Look for intersections within the element interior.
                if (u >= 0.0 && v >= 0.0 && (u+v) <= 1.0)
                {
                    double t;
                    if (std::abs(q(0)) >= std::abs(q(1)) && std::abs(q(0)) >= std::abs(q(2)))
                    {
                        const double p = u*p0(0) + v*p1(0) + (1.0-u-v)*p2(0);
                        t = (p-r(0))/q(0);
                    }
                    else if (std::abs(q(1)) >= std::abs(q(2)))
                    {
                        const double p = u*p0(1) + v*p1(1) + (1.0-u-v)*p2(1);
                        t = (p-r(1))/q(1);
                    }
                    else
                    {
                        const double p = u*p0(2) + v*p1(2) + (1.0-u-v)*p2(2);
                        t = (p-r(2))/q(2);
                    }
                    t_vals.push_back(std::make_pair(t,libMesh::Point(u,v,0.0)));
                }
            }
            break;
        }
        case libMeshEnums::QUAD4:
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
                D1 = d(1)-r(1);
                D2 = d(2)-r(2);
            }
            else if (q(1) != 0.0)
            {
                A1 = a(0);
                A2 = a(2);
                B1 = b(0);
                B2 = b(2);
                C1 = c(0);
                C2 = c(2);
                D1 = d(0)-r(0);
                D2 = d(2)-r(2);
            }
            else
            {
                A1 = a(0);
                A2 = a(1);
                B1 = b(0);
                B2 = b(1);
                C1 = c(0);
                C2 = c(1);
                D1 = d(0)-r(0);
                D2 = d(1)-r(1);
            }

            // (A2*C1 - A1*C2) v^2 + (A2*D1 - A1*D2 + B2*C1 - B1*C2) v + (B2*D1 - B1*D2) = 0
            std::vector<double> v_vals;
            {
                const double a = A2*C1 - A1*C2;
                const double b = A2*D1 - A1*D2 + B2*C1 - B1*C2;
                const double c = B2*D1 - B1*D2;
                const double disc = b*b - 4.0*a*c;
                if (disc > 0.0)
                {
                    const double q = -0.5*(b+(b>0.0 ? 1.0 : -1.0)*sqrt(disc));
                    const double v0 = q/a;
                    v_vals.push_back(v0);
                    const double v1 = c/q;
                    if (std::abs(v0-v1) > std::numeric_limits<double>::epsilon())
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
                    const double a = v*A2+B2;
                    const double b = v*(A2-A1)+B2-B1;
                    if (std::abs(b) >= std::abs(a))
                    {
                        u = (v*(C1-C2)+D1-D2)/b;
                    }
                    else
                    {
                        u = (-v*C2-D2)/a;
                    }

                    if (u >= 0.0 && u <= 1.0)
                    {
                        double t;
                        if (std::abs(q(0)) >= std::abs(q(1)) && std::abs(q(0)) >= std::abs(q(2)))
                        {
                            const double p = p00(0)*(1.0-u)*(1.0-v) + p01(0)*(1.0-u)*v + p10(0)*u*(1.0-v) + p11(0)*u*v;
                            t = (p-r(0))/q(0);
                        }
                        else if (std::abs(q(1)) >= std::abs(q(2)))
                        {
                            const double p = p00(1)*(1.0-u)*(1.0-v) + p01(1)*(1.0-u)*v + p10(1)*u*(1.0-v) + p11(1)*u*v;
                            t = (p-r(1))/q(1);
                        }
                        else
                        {
                            const double p = p00(2)*(1.0-u)*(1.0-v) + p01(2)*(1.0-u)*v + p10(2)*u*(1.0-v) + p11(2)*u*v;
                            t = (p-r(2))/q(2);
                        }
                        t_vals.push_back(std::make_pair(t,libMesh::Point(2.0*u-1.0,2.0*v-1.0,0.0)));
                    }
                }
            }
            break;
        }
        default:
        {
            TBOX_ERROR("intersect_line_with_face():"
                       << "  element type " << libMesh::Utility::enum_to_string<libMeshEnums::ElemType>(elem->type()) << " is not supported at this time.\n");
        }
    }
    return t_vals;
}// intersect_line_with_face

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
