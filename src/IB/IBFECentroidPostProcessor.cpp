// Filename: IBFEPostProcessor.cpp
// Created on 4 Dec 2013 by Boyce Griffith
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

#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "boost/multi_array.hpp"
#include "ibamr/IBFECentroidPostProcessor.h"
#include "ibamr/IBFEMethod.h"
#include "ibamr/IBFEPostProcessor.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/FEDataManager.h"
#include "ibtk/libmesh_utilities.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/dof_map.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_type.h"
#include "libmesh/fem_context.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"
#include "libmesh/system.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_tensor.h"
#include "libmesh/type_vector.h"
#include "libmesh/variant_filter_iterator.h"
#include "libmesh/vector_value.h"
#include "tbox/Utilities.h"

namespace libMesh
{
class Elem;
} // namespace libMesh

using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFECentroidPostProcessor::IBFECentroidPostProcessor(const std::string& name, FEDataManager* fe_data_manager)
    : IBFEPostProcessor(name, fe_data_manager)
{
    // intentionally blank
    return;
} // IBFECentroidPostProcessor

IBFECentroidPostProcessor::~IBFECentroidPostProcessor()
{
    // intentionally blank
    return;
} // ~IBFECentroidPostProcessor

void
IBFECentroidPostProcessor::registerScalarVariable(const std::string& name,
                                                  libMesh::FEFamily fe_family,
                                                  libMesh::Order fe_order,
                                                  ScalarMeshFcnPtr fcn,
                                                  const std::vector<SystemData>& system_data,
                                                  void* fcn_ctx)
{
    TBOX_ASSERT(fe_family == MONOMIAL);
    TBOX_ASSERT(fe_order == CONSTANT);
    IBFEPostProcessor::registerScalarVariable(name, fe_family, fe_order, fcn, system_data, fcn_ctx);
    return;
} // registerScalarVariable

void
IBFECentroidPostProcessor::registerVectorVariable(const std::string& name,
                                                  libMesh::FEFamily fe_family,
                                                  libMesh::Order fe_order,
                                                  VectorMeshFcnPtr fcn,
                                                  const std::vector<SystemData>& system_data,
                                                  void* fcn_ctx,
                                                  unsigned int dim)
{
    TBOX_ASSERT(fe_family == MONOMIAL);
    TBOX_ASSERT(fe_order == CONSTANT);
    IBFEPostProcessor::registerVectorVariable(name, fe_family, fe_order, fcn, system_data, fcn_ctx, dim);
    return;
} // registerVectorVariable

void
IBFECentroidPostProcessor::registerTensorVariable(const std::string& name,
                                                  libMesh::FEFamily fe_family,
                                                  libMesh::Order fe_order,
                                                  TensorMeshFcnPtr fcn,
                                                  const std::vector<SystemData>& system_data,
                                                  void* fcn_ctx,
                                                  unsigned int dim)
{
    TBOX_ASSERT(fe_family == MONOMIAL);
    TBOX_ASSERT(fe_order == CONSTANT);
    IBFEPostProcessor::registerTensorVariable(name, fe_family, fe_order, fcn, system_data, fcn_ctx, dim);
    return;
} // registerTensorVariable

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBFECentroidPostProcessor::reconstructVariables(double data_time)
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const int dim = mesh.mesh_dimension();
    AutoPtr<QBase> qrule = QBase::build(QGAUSS, NDIM, CONSTANT);

    // Set up all system data required to evaluate the mesh functions.
    System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_dof_map.variable_type(0));
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
    AutoPtr<FEBase> X_fe(FEBase::build(dim, X_dof_map.variable_type(0)));
    X_fe->attach_quadrature_rule(qrule.get());
    const std::vector<libMesh::Point>& q_point = X_fe->get_xyz();
    const std::vector<std::vector<double> >& phi_X = X_fe->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& dphi_X = X_fe->get_dphi();
    X_system.solution->localize(*X_system.current_local_solution);
    NumericVector<double>& X_data = *(X_system.current_local_solution);
    X_data.close();

    TBOX_ERROR("currently broken!");

    const size_t num_scalar_vars = d_scalar_var_systems.size();
    std::vector<const DofMap*> scalar_var_dof_maps(num_scalar_vars);
    std::vector<std::vector<unsigned int> > scalar_var_dof_indices(num_scalar_vars);
    std::vector<NumericVector<double>*> scalar_var_data(num_scalar_vars);
    std::vector<unsigned int> scalar_var_system_num(num_scalar_vars);
    std::vector<std::vector<const std::vector<double>*> > scalar_var_fcn_data(num_scalar_vars);
    std::vector<std::vector<const std::vector<VectorValue<double> >*> > scalar_var_grad_fcn_data(num_scalar_vars);
    TBOX_WARNING("extra function data not treated correctly here.\n");
    for (unsigned int k = 0; k < num_scalar_vars; ++k)
    {
        scalar_var_dof_maps[k] = &d_scalar_var_systems[k]->get_dof_map();
        scalar_var_data[k] = d_scalar_var_systems[k]->solution.get();
        scalar_var_system_num[k] = d_scalar_var_systems[k]->number();
        scalar_var_fcn_data[k].reserve(d_scalar_var_system_data[k].size());
        scalar_var_grad_fcn_data[k].reserve(d_scalar_var_system_data[k].size());
    }

    const size_t num_vector_vars = d_vector_var_systems.size();
    std::vector<const DofMap*> vector_var_dof_maps(num_vector_vars);
    std::vector<std::vector<std::vector<unsigned int> > > vector_var_dof_indices(num_vector_vars);
    std::vector<NumericVector<double>*> vector_var_data(num_vector_vars);
    std::vector<unsigned int> vector_var_system_num(num_vector_vars);
    std::vector<std::vector<const std::vector<double>*> > vector_var_fcn_data(num_vector_vars);
    std::vector<std::vector<const std::vector<VectorValue<double> >*> > vector_var_grad_fcn_data(num_vector_vars);
    TBOX_WARNING("extra function data not treated correctly here.\n");
    for (unsigned int k = 0; k < num_vector_vars; ++k)
    {
        vector_var_dof_maps[k] = &d_vector_var_systems[k]->get_dof_map();
        vector_var_dof_indices[k].resize(d_vector_var_dims[k]);
        vector_var_data[k] = d_vector_var_systems[k]->solution.get();
        vector_var_system_num[k] = d_vector_var_systems[k]->number();
        vector_var_fcn_data[k].reserve(d_vector_var_system_data[k].size());
        vector_var_grad_fcn_data[k].reserve(d_vector_var_system_data[k].size());
    }

    const size_t num_tensor_vars = d_tensor_var_systems.size();
    std::vector<const DofMap*> tensor_var_dof_maps(num_tensor_vars);
    std::vector<boost::multi_array<std::vector<unsigned int>, 2> > tensor_var_dof_indices(num_tensor_vars);
    std::vector<NumericVector<double>*> tensor_var_data(num_tensor_vars);
    std::vector<unsigned int> tensor_var_system_num(num_tensor_vars);
    std::vector<std::vector<const std::vector<double>*> > tensor_var_fcn_data(num_tensor_vars);
    std::vector<std::vector<const std::vector<VectorValue<double> >*> > tensor_var_grad_fcn_data(num_tensor_vars);
    TBOX_WARNING("extra function data not treated correctly here.\n");
    for (unsigned int k = 0; k < num_tensor_vars; ++k)
    {
        tensor_var_dof_maps[k] = &d_tensor_var_systems[k]->get_dof_map();
        typedef boost::multi_array<std::vector<unsigned int>, 2> array_type;
        array_type::extent_gen extents;
        tensor_var_dof_indices[k].resize(extents[d_tensor_var_dims[k]][d_tensor_var_dims[k]]);
        tensor_var_data[k] = d_tensor_var_systems[k]->solution.get();
        tensor_var_system_num[k] = d_tensor_var_systems[k]->number();
        tensor_var_fcn_data[k].reserve(d_tensor_var_system_data[k].size());
        tensor_var_grad_fcn_data[k].reserve(d_tensor_var_system_data[k].size());
    }

    // Reconstruct the variables via simple function evaluation.
    TensorValue<double> FF_qp, VV;
    libMesh::Point X_qp;
    VectorValue<double> V;
    double v;
    boost::multi_array<double, 2> X_node;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        X_fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[d], d);
        }

        const unsigned int n_qp = qrule->n_points();
        TBOX_ASSERT(n_qp == 1);
        const unsigned int qp = 0;

        get_values_for_interpolation(X_node, X_data, X_dof_indices);
        interpolate(X_qp, qp, X_node, phi_X);
        jacobian(FF_qp, qp, X_node, dphi_X);
        const libMesh::Point& s_qp = q_point[qp];

        // Scalar-valued variables.
        for (unsigned int k = 0; k < num_scalar_vars; ++k)
        {
            scalar_var_dof_maps[k]->dof_indices(elem, scalar_var_dof_indices[k], 0);
            d_scalar_var_fcns[k](v,
                                 FF_qp,
                                 X_qp,
                                 s_qp,
                                 elem,
                                 scalar_var_fcn_data[k],
                                 scalar_var_grad_fcn_data[k],
                                 data_time,
                                 d_scalar_var_fcn_ctxs[k]);
            scalar_var_data[k]->set(scalar_var_dof_indices[k][0], v);
        }

        // Vector-valued variables.
        for (unsigned int k = 0; k < num_vector_vars; ++k)
        {
            for (unsigned int i = 0; i < d_vector_var_dims[k]; ++i)
            {
                vector_var_dof_maps[k]->dof_indices(elem, vector_var_dof_indices[k][i], i);
            }
            d_vector_var_fcns[k](V,
                                 FF_qp,
                                 X_qp,
                                 s_qp,
                                 elem,
                                 vector_var_fcn_data[k],
                                 vector_var_grad_fcn_data[k],
                                 data_time,
                                 d_vector_var_fcn_ctxs[k]);
            for (unsigned int i = 0; i < d_vector_var_dims[k]; ++i)
            {
                vector_var_data[k]->set(vector_var_dof_indices[k][i][0], V(i));
            }
        }

        // Tensor-valued variables.
        for (unsigned int k = 0; k < num_tensor_vars; ++k)
        {
            for (unsigned int i = 0; i < d_tensor_var_dims[k]; ++i)
            {
                for (unsigned int j = 0; j < d_tensor_var_dims[k]; ++j)
                {
                    tensor_var_dof_maps[k]->dof_indices(
                        elem, tensor_var_dof_indices[k][i][j], j + i * d_tensor_var_dims[k]);
                }
            }
            d_tensor_var_fcns[k](VV,
                                 FF_qp,
                                 X_qp,
                                 s_qp,
                                 elem,
                                 tensor_var_fcn_data[k],
                                 tensor_var_grad_fcn_data[k],
                                 data_time,
                                 d_tensor_var_fcn_ctxs[k]);
            for (unsigned int i = 0; i < d_tensor_var_dims[k]; ++i)
            {
                for (unsigned int j = 0; j < d_tensor_var_dims[k]; ++j)
                {
                    tensor_var_data[k]->set(tensor_var_dof_indices[k][i][j][0], VV(i, j));
                }
            }
        }
    }

    // Close all vectors.
    for (unsigned int k = 0; k < num_scalar_vars; ++k)
    {
        scalar_var_data[k]->close();
    }
    for (unsigned int k = 0; k < num_vector_vars; ++k)
    {
        vector_var_data[k]->close();
    }
    for (unsigned int k = 0; k < num_tensor_vars; ++k)
    {
        tensor_var_data[k]->close();
    }
    return;
} // reconstructVariables

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
