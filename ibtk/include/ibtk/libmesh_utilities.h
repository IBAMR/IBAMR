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

#include "ibtk/IBTK_CHKERRQ.h"
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

#include <boost/multi_array.hpp>

#include <tuple>

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
               libMesh::NumericVector<double>* const system_vec = nullptr)
        : system_name(system_name), vars(vars), grad_vars(grad_vars), system_vec(system_vec)
    {
    }

    std::string system_name;
    std::vector<int> vars, grad_vars;
    libMesh::NumericVector<double>* system_vec;
};

using ScalarMeshFcnPtr =
    void (*)(double& F,
             const libMesh::TensorValue<double>& FF,
             const libMesh::Point& x,
             const libMesh::Point& X,
             libMesh::Elem* elem,
             const std::vector<const std::vector<double>*>& system_var_data,
             const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
             double data_time,
             void* ctx);

using VectorMeshFcnPtr =
    void (*)(libMesh::VectorValue<double>& F,
             const libMesh::TensorValue<double>& FF,
             const libMesh::Point& x,
             const libMesh::Point& X,
             libMesh::Elem* elem,
             const std::vector<const std::vector<double>*>& system_var_data,
             const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
             double data_time,
             void* ctx);

using TensorMeshFcnPtr =
    void (*)(libMesh::TensorValue<double>& F,
             const libMesh::TensorValue<double>& FF,
             const libMesh::Point& x,
             const libMesh::Point& X,
             libMesh::Elem* elem,
             const std::vector<const std::vector<double>*>& system_var_data,
             const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
             double data_time,
             void* ctx);

using ScalarSurfaceFcnPtr =
    void (*)(double& F,
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

using VectorSurfaceFcnPtr =
    void (*)(libMesh::VectorValue<double>& F,
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

using TensorSurfaceFcnPtr =
    void (*)(libMesh::TensorValue<double>& F,
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

inline void
copy_and_synch(libMesh::NumericVector<double>& v_in,
               libMesh::NumericVector<double>& v_out,
               const bool close_v_in = true,
               const bool close_v_out = true)
{
#if defined(NDEBUG)
    auto v_in_petsc = static_cast<libMesh::PetscVector<double>*>(&v_in);
    auto v_out_petsc = static_cast<libMesh::PetscVector<double>*>(&v_out);
#else
    auto v_in_petsc = dynamic_cast<libMesh::PetscVector<double>*>(&v_in);
    auto v_out_petsc = dynamic_cast<libMesh::PetscVector<double>*>(&v_out);
    TBOX_ASSERT(v_in_petsc);
    TBOX_ASSERT(v_out_petsc);
#endif
    if (close_v_in) v_in.close();
    PetscErrorCode ierr = VecCopy(v_in_petsc->vec(), v_out_petsc->vec());
    IBTK_CHKERRQ(ierr);
    if (close_v_out) v_out.close();
}

inline void
batch_vec_copy(const std::vector<libMesh::PetscVector<double>*>& x_vecs,
               const std::vector<libMesh::PetscVector<double>*>& y_vecs)
{
#if defined(NDEBUG)
    TBOX_ASSERT(x_vecs.size() == y_vecs.size());
#endif
    for (unsigned int k = 0; k < x_vecs.size(); ++k)
    {
        if (!x_vecs[k] || !y_vecs[k]) continue;
        int ierr = VecCopy(x_vecs[k]->vec(), y_vecs[k]->vec());
        IBTK_CHKERRQ(ierr);
    }
}

inline void
batch_vec_copy(const std::vector<std::vector<libMesh::PetscVector<double>*> >& x_vecs,
               const std::vector<std::vector<libMesh::PetscVector<double>*> >& y_vecs)
{
#if defined(NDEBUG)
    TBOX_ASSERT(x_vecs.size() == y_vecs.size());
#endif
    for (unsigned int n = 0; n < x_vecs.size(); ++n)
    {
        batch_vec_copy(x_vecs[n], y_vecs[n]);
    }
}

inline void
batch_vec_assembly(const std::vector<libMesh::PetscVector<double>*>& vecs)
{
    for (const auto& v : vecs)
    {
        if (!v) continue;
        int ierr = VecAssemblyBegin(v->vec());
        IBTK_CHKERRQ(ierr);
    }
    for (const auto& v : vecs)
    {
        if (!v) continue;
        int ierr = VecAssemblyEnd(v->vec());
        IBTK_CHKERRQ(ierr);
    }
}

inline void
batch_vec_assembly(const std::vector<std::vector<libMesh::PetscVector<double>*> >& vecs)
{
    for (unsigned int n = 0; n < vecs.size(); ++n)
    {
        for (const auto& v : vecs[n])
        {
            if (!v) continue;
            int ierr = VecAssemblyBegin(v->vec());
            IBTK_CHKERRQ(ierr);
        }
    }
    for (unsigned int n = 0; n < vecs.size(); ++n)
    {
        for (const auto& v : vecs[n])
        {
            if (!v) continue;
            int ierr = VecAssemblyEnd(v->vec());
            IBTK_CHKERRQ(ierr);
        }
    }
}

inline void
batch_vec_ghost_update(const std::vector<libMesh::PetscVector<double>*>& vecs,
                       const InsertMode insert_mode,
                       const ScatterMode scatter_mode)
{
    for (const auto& v : vecs)
    {
        if (!v) continue;
        int ierr = VecGhostUpdateBegin(v->vec(), insert_mode, scatter_mode);
        IBTK_CHKERRQ(ierr);
    }
    for (const auto& v : vecs)
    {
        if (!v) continue;
        int ierr = VecGhostUpdateEnd(v->vec(), insert_mode, scatter_mode);
        IBTK_CHKERRQ(ierr);
    }
}

inline void
batch_vec_ghost_update(const std::vector<std::vector<libMesh::PetscVector<double>*> >& vecs,
                       const InsertMode insert_mode,
                       const ScatterMode scatter_mode)
{
    for (unsigned int n = 0; n < vecs.size(); ++n)
    {
        for (const auto& v : vecs[n])
        {
            if (!v) continue;
            int ierr = VecGhostUpdateBegin(v->vec(), insert_mode, scatter_mode);
            IBTK_CHKERRQ(ierr);
        }
    }
    for (unsigned int n = 0; n < vecs.size(); ++n)
    {
        for (const auto& v : vecs[n])
        {
            if (!v) continue;
            int ierr = VecGhostUpdateEnd(v->vec(), insert_mode, scatter_mode);
            IBTK_CHKERRQ(ierr);
        }
    }
}

/**
 * Return the quadrature key description (see QuadratureCache, FECache, and
 * FEMapCache) of a quadrature rule.
 *
 * @param[in] quad_type The type of quadrature rule to use. At the present
 * time only QGAUSS and QGRID are supported.
 *
 * @param[in] order The order of the quadrature rule.
 *
 * @param[in] use_adaptive_quadrature Whether or not the element data should
 * be read (and @p order possibly modified) to get the correct quadrature
 * point spacing.
 *
 * @param[in] point_density Parameter for computing the number of points in
 * the quadrature rule. The number of points is computed as
 * <code>ceil(point_density * hmax / dx_min)</code>, where <code>hmax</code>
 * is the maximum edge length of the deformed element (i.e., with nodal
 * coordinates given by @p X_node).
 *
 * @param[in] elem The libMesh element. The quadrature rule generated by this
 * function will always have enough points to integrate the basis functions
 * defined on the element exactly.
 *
 * @param[in] X_node the actual nodal coordinates of the element in the
 * current configuration.
 *
 * @param[in] dx_min See @p point_density.
 *
 * @seealso FEDataManager::updateQuadratureRule.
 */
std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>
getQuadratureKey(const libMesh::QuadratureType quad_type,
                 libMesh::Order order,
                 const bool use_adaptive_quadrature,
                 const double point_density,
                 const libMesh::Elem* const elem,
                 const boost::multi_array<double, 2>& X_node,
                 const double dx_min);

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
    const double* const U_local_soln = U_petsc_vec->get_array_read();
    get_values_for_interpolation(U_node, *U_petsc_vec, U_local_soln, dof_indices);
    U_petsc_vec->restore_array();
    return;
} // get_values_for_interpolation

template <class MultiArray>
inline void
get_values_for_interpolation(MultiArray& U_node,
                             libMesh::NumericVector<double>& U_vec,
                             const std::vector<std::vector<unsigned int> >& dof_indices)
{
    libMesh::PetscVector<double>* U_petsc_vec = dynamic_cast<libMesh::PetscVector<double>*>(&U_vec);
    const double* const U_local_soln = U_petsc_vec->get_array_read();
    get_values_for_interpolation(U_node, *U_petsc_vec, U_local_soln, dof_indices);
    U_petsc_vec->restore_array();
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
inline bool
intersect_line_with_edge(std::vector<std::pair<double, libMesh::Point> >& t_vals,
                         libMesh::Edge* elem,
                         libMesh::Point r,
                         libMesh::VectorValue<double> q,
                         const double tol = 0.0)
{
    bool is_interior_intersection = false;
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
        const libMesh::Point& p0 = *elem->node_ptr(0);
        const libMesh::Point& p1 = *elem->node_ptr(1);
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
        if (u >= -1.0 - tol && u <= 1.0 + tol)
        {
            is_interior_intersection = (u >= -1.0 && u <= 1.0);
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
        const libMesh::Point& p0 = *elem->node_ptr(0);
        const libMesh::Point& p1 = *elem->node_ptr(1);
        const libMesh::Point& p2 = *elem->node_ptr(2);
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
            const double q = -0.5 * (b + (b > 0.0 ? 1.0 : -1.0) * std::sqrt(disc));
            const double u0 = q / a;
            u_vals.push_back(u0);
            const double u1 = c / q;
            if (std::abs(u0 - u1) > std::numeric_limits<double>::epsilon())
            {
                u_vals.push_back(u1);
            }
        }

        // Look for intersections within the element interior.
        for (const auto& u : u_vals)
        {
            if (u >= -1.0 - tol && u <= 1.0 + tol)
            {
                is_interior_intersection = (u >= -1.0 && u <= 1.0);
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
                   << "  element type " << libMesh::Utility::enum_to_string<libMesh::ElemType>(elem->type())
                   << " is not supported at this time.\n");
    }
    }
    return is_interior_intersection;
} // intersect_line_with_edge

// WARNING: This code is specialized to the case in which q is a unit vector
// aligned with the coordinate axes.
inline bool
intersect_line_with_face(std::vector<std::pair<double, libMesh::Point> >& t_vals,
                         libMesh::Face* elem,
                         libMesh::Point r,
                         libMesh::VectorValue<double> q,
                         const double tol = 0.0)
{
    bool is_interior_intersection = false;
    t_vals.resize(0);
    switch (elem->type())
    {
    case libMesh::TRI3:
    {
        const libMesh::VectorValue<double>& p = r;
        const libMesh::VectorValue<double>& d = q;

        // Linear interpolation:
        //
        //    (1-u-v)*p0 + u*p1 + v*p2 = p + t*d
        //
        // https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
        // http://www.lighthouse3d.com/tutorials/maths/ray-triangle-intersection/
        const libMesh::Point& p0 = *elem->node_ptr(0);
        const libMesh::Point& p1 = *elem->node_ptr(1);
        const libMesh::Point& p2 = *elem->node_ptr(2);

        const libMesh::VectorValue<double> e1 = p1 - p0;
        const libMesh::VectorValue<double> e2 = p2 - p0;
        const libMesh::VectorValue<double> h = d.cross(e2);
        double a = e1 * h;
        if (std::abs(a) > std::numeric_limits<double>::epsilon())
        {
            double f = 1.0 / a;
            const libMesh::VectorValue<double> s = p - p0;
            double u = f * (s * h);
            if (u >= -tol && u <= 1.0 + tol)
            {
                const libMesh::VectorValue<double> q = s.cross(e1);
                double v = f * (d * q);
                if (v >= tol && (u + v) <= 1.0 + tol)
                {
                    double t = f * (e2 * q);
                    t_vals.push_back(std::make_pair(t, libMesh::Point(u, v, 0.0)));
                }
            }
        }
        break;
    }
    case libMesh::QUAD4:
    {
        const libMesh::Point& p00 = *elem->node_ptr(0);
        const libMesh::Point& p10 = *elem->node_ptr(1);
        const libMesh::Point& p11 = *elem->node_ptr(2);
        const libMesh::Point& p01 = *elem->node_ptr(3);

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
                const double q = -0.5 * (b + (b > 0.0 ? 1.0 : -1.0) * std::sqrt(disc));
                const double v0 = q / a;
                v_vals.push_back(v0);
                const double v1 = c / q;
                if (std::abs(v0 - v1) > std::numeric_limits<double>::epsilon())
                {
                    v_vals.push_back(v1);
                }
            }
        }

        for (const auto& v : v_vals)
        {
            if (v >= 0.0 - tol && v <= 1.0 + tol)
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

                if (u >= 0.0 - tol && u <= 1.0 + tol)
                {
                    is_interior_intersection = (u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 0.0);
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
                   << "  element type " << libMesh::Utility::enum_to_string<libMesh::ElemType>(elem->type())
                   << " is not supported at this time.\n");
    }
    }
    return is_interior_intersection;
} // intersect_line_with_face

struct DofObjectComp
{
    inline bool operator()(const libMesh::DofObject* const lhs, const libMesh::DofObject* const rhs) const
    {
        return lhs->id() < rhs->id();
    }
};

/*!
 * Recent versions of libMesh acquired a useful function that lets us extract
 * the DoFs corresponding to basis functions with node value functionals. This
 * compatibility function either calls that function directly or uses our own
 * implementation if the present libMesh is too old.
 */
inline void
get_nodal_dof_indices(const libMesh::DofMap& dof_map,
                      const libMesh::Node* const node,
                      const unsigned int variable_n,
                      std::vector<libMesh::dof_id_type>& nodal_indices)
{
#if LIBMESH_MINOR_VERSION < 2
    // See dof_map.C, circa line 2208

    // We only call this function with variable numbers 0, 1, or 2, so skip
    // implementing some stuff
    TBOX_ASSERT(variable_n != libMesh::invalid_uint);

    nodal_indices.clear();
    const unsigned int system_n = dof_map.sys_number();

    // Get the dof numbers
    const libMesh::Variable& var = dof_map.variable(variable_n);
    if (var.type().family == libMesh::SCALAR)
    {
        dof_map.SCALAR_dof_indices(nodal_indices, variable_n);
    }
    else
    {
        const int n_components = node->n_comp(system_n, variable_n);
        for (int component_n = 0; component_n != n_components; ++component_n)
        {
            libmesh_assert_not_equal_to(node->dof_number(system_n, variable_n, component_n),
                                        libMesh::DofObject::invalid_id);
            nodal_indices.push_back(node->dof_number(system_n, variable_n, component_n));
        }
    }
#else
    dof_map.dof_indices(node, nodal_indices, variable_n);
#endif
}

/**
 * Return the maximum edge length of a given element with mapped nodes.
 *
 * @p elem The given libMesh element.
 * @p X_node Mapped coordinates of the nodes.
 */
template <class MultiArray>
inline double
get_max_edge_length(const libMesh::Elem* const elem, const MultiArray& X_node)
{
#ifndef NDEBUG
    TBOX_ASSERT(elem->n_nodes() == X_node.shape()[0]);
    TBOX_ASSERT(elem->dim() == X_node.shape()[1]);
#endif
    const libMesh::ElemType elem_type = elem->type();

    double max_edge_length_2 = 0.0;
    // TODO: implement additional element-specific versions as needed.
    switch (elem_type)
    {
    case libMesh::TET4:
    {
        // all nodes are connected to all other nodes
        constexpr unsigned int edge_pairs[6][2] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 }, { 1, 3 }, { 2, 3 } };
        for (const unsigned int(&pair)[2] : edge_pairs)
        {
            const unsigned int n1 = pair[0];
            const unsigned int n2 = pair[1];
            double diff_sq = 0.0;
            for (unsigned int d = 0; d < 3; ++d) // TET4 implies 3D
                diff_sq += (X_node[n1][d] - X_node[n2][d]) * (X_node[n1][d] - X_node[n2][d]);
            max_edge_length_2 = std::max(max_edge_length_2, diff_sq);
        }
        break;
    }
    case libMesh::HEX8:
    {
        // see the connectivity diagram in cell_hex8.h to see why these are the edges
        constexpr unsigned int edge_pairs[12][2] = { { 0, 1 }, { 0, 3 }, { 0, 4 }, { 1, 2 }, { 1, 5 }, { 2, 3 },
                                                     { 2, 6 }, { 3, 7 }, { 4, 5 }, { 4, 7 }, { 5, 6 }, { 6, 7 } };
        for (const unsigned int(&pair)[2] : edge_pairs)
        {
            const unsigned int n1 = pair[0];
            const unsigned int n2 = pair[1];
            double diff_sq = 0.0;
            for (unsigned int d = 0; d < 3; ++d) // HEX8 implies 3D
                diff_sq += (X_node[n1][d] - X_node[n2][d]) * (X_node[n1][d] - X_node[n2][d]);
            max_edge_length_2 = std::max(max_edge_length_2, diff_sq);
        }
        break;
    }
    default:
        // Use the old algorithm in all other cases
        {
            const unsigned int n_vertices = elem->n_vertices();
            if (elem->dim() == 1)
            {
                for (unsigned int n1 = 0; n1 < n_vertices; ++n1)
                {
                    for (unsigned int n2 = n1 + 1; n2 < n_vertices; ++n2)
                    {
                        double diff_sq = 0.0;
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            diff_sq += (X_node[n1][d] - X_node[n2][d]) * (X_node[n1][d] - X_node[n2][d]);
                        }
                        max_edge_length_2 = std::max(max_edge_length_2, diff_sq);
                    }
                }
            }
            else
            {
                const unsigned int n_edges = elem->n_edges();
                for (unsigned int e = 0; e < n_edges; ++e)
                {
                    for (unsigned int n1 = 0; n1 < n_vertices; ++n1)
                    {
                        if (!elem->is_node_on_edge(n1, e)) continue;
                        for (unsigned int n2 = n1 + 1; n2 < n_vertices; ++n2)
                        {
                            if (!elem->is_node_on_edge(n2, e)) continue;
                            double diff_sq = 0.0;
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                diff_sq += (X_node[n1][d] - X_node[n2][d]) * (X_node[n1][d] - X_node[n2][d]);
                            }
                            max_edge_length_2 = std::max(max_edge_length_2, diff_sq);
                        }
                    }
                }
            }
        }
    }

    return std::sqrt(max_edge_length_2);
}

/*!
 * Save, in a plain text file, the libMesh partitioning, with the format
 *
 *     x,y,z,rank
 *
 * where x, y, and z are the center of an Elem and rank is the current MPI
 * rank.
 *
 * @note this function collates the output from all MPI processors in the
 * communicator assigned to @p position_system, so it is an inherently serial
 * function.
 */
void write_elem_partitioning(const std::string& file_name, const libMesh::System& position_system);

/*!
 * Save, in a plain text file, the libMesh Node partitioning, with the format
 *
 *     x,y,z,rank
 *
 * where x, y, and z are the coordinates of a Node and rank is the current MPI
 * rank.
 *
 * @note this function collates the output from all MPI processors in the
 * communicator assigned to @p position_system, so it is an inherently serial
 * function.
 */
void write_node_partitioning(const std::string& file_name, const libMesh::System& position_system);
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_libmesh_utilities
