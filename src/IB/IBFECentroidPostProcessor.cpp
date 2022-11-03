// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2022 by the IBAMR developers
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

#include "ibamr/IBFECentroidPostProcessor.h"
#include "ibamr/IBFEDirectForcingKinematics.h"
#include "ibamr/IBFEMethod.h"
#include "ibamr/IBFEPostProcessor.h"

#include "ibtk/FEDataInterpolation.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/libmesh_utilities.h"

#include "tbox/Utilities.h"

#include "libmesh/dof_map.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/libmesh_config.h"
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

#include "ibamr/namespaces.h" // IWYU pragma: keep

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace libMesh
{
class Elem;
} // namespace libMesh

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline void
get_x_and_FF(libMesh::VectorValue<double>& x,
             libMesh::TensorValue<double>& FF,
             const std::vector<double>& x_data,
             const std::vector<VectorValue<double> >& grad_x_data,
             const unsigned int dim = NDIM)
{
    x.zero();
    FF.zero();
    for (unsigned int i = 0; i < dim; ++i)
    {
        x(i) = x_data[i];
        for (unsigned int j = 0; j < dim; ++j)
        {
            FF(i, j) = grad_x_data[i](j);
        }
    }
    for (unsigned int i = dim; i < LIBMESH_DIM; ++i)
    {
        FF(i, i) = 1.0;
    }
    return;
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFECentroidPostProcessor::IBFECentroidPostProcessor(std::string name, FEDataManager* fe_data_manager)
    : IBFEPostProcessor(std::move(name), fe_data_manager)
{
    // intentionally blank
    return;
} // IBFECentroidPostProcessor

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
    MeshBase& mesh = equation_systems->get_mesh();
    const int dim = mesh.mesh_dimension();
    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, NDIM, CONSTANT);

    // Set up all system data required to evaluate the mesh functions.
    FEDataInterpolation fe(dim, d_fe_data_manager->getFEData());
    fe.attachQuadratureRule(qrule.get());
    fe.evalQuadraturePoints();

    auto& X_system = equation_systems->get_system<System>(d_fe_data_manager->getCurrentCoordinatesSystemName());
    NumericVector<double>& X_data = *(X_system.current_local_solution);
    copy_and_synch(*X_system.solution, X_data);
    std::vector<int> vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
    const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, vars, vars, &X_data);

    const size_t num_scalar_vars = d_scalar_var_systems.size();
    std::vector<const DofMap*> scalar_var_dof_maps(num_scalar_vars);
    std::vector<std::vector<unsigned int> > scalar_var_dof_indices(num_scalar_vars);
    std::vector<std::vector<size_t> > scalar_var_fcn_system_idxs(num_scalar_vars);
    for (unsigned int k = 0; k < num_scalar_vars; ++k)
    {
        scalar_var_dof_maps[k] = &d_scalar_var_systems[k]->get_dof_map();
        fe.setupInterpolatedSystemDataIndexes(
            scalar_var_fcn_system_idxs[k], d_scalar_var_system_data[k], equation_systems);
    }

    const size_t num_vector_vars = d_vector_var_systems.size();
    std::vector<const DofMap*> vector_var_dof_maps(num_vector_vars);
    std::vector<std::vector<std::vector<unsigned int> > > vector_var_dof_indices(num_vector_vars);
    std::vector<std::vector<size_t> > vector_var_fcn_system_idxs(num_vector_vars);
    for (unsigned int k = 0; k < num_vector_vars; ++k)
    {
        vector_var_dof_maps[k] = &d_vector_var_systems[k]->get_dof_map();
        vector_var_dof_indices[k].resize(d_vector_var_dims[k]);
        fe.setupInterpolatedSystemDataIndexes(
            vector_var_fcn_system_idxs[k], d_vector_var_system_data[k], equation_systems);
    }

    const size_t num_tensor_vars = d_tensor_var_systems.size();
    std::vector<const DofMap*> tensor_var_dof_maps(num_tensor_vars);
    std::vector<boost::multi_array<std::vector<unsigned int>, 2> > tensor_var_dof_indices(num_tensor_vars);
    std::vector<std::vector<size_t> > tensor_var_fcn_system_idxs(num_tensor_vars);
    for (unsigned int k = 0; k < num_tensor_vars; ++k)
    {
        tensor_var_dof_maps[k] = &d_tensor_var_systems[k]->get_dof_map();
        using array_type = boost::multi_array<std::vector<unsigned int>, 2>;
        array_type::extent_gen extents;
        tensor_var_dof_indices[k].resize(extents[d_tensor_var_dims[k]][d_tensor_var_dims[k]]);
        fe.setupInterpolatedSystemDataIndexes(
            tensor_var_fcn_system_idxs[k], d_tensor_var_system_data[k], equation_systems);
    }

    fe.init();

    const std::vector<libMesh::Point>& q_point = fe.getQuadraturePoints();

    const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
    const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
        fe.getGradVarInterpolation();

    std::vector<std::vector<const std::vector<double>*> > scalar_var_data(num_scalar_vars),
        vector_var_data(num_vector_vars), tensor_var_data(num_tensor_vars);
    std::vector<std::vector<const std::vector<VectorValue<double> >*> > scalar_grad_var_data(num_scalar_vars),
        vector_grad_var_data(num_vector_vars), tensor_grad_var_data(num_tensor_vars);

    // Reconstruct the variables via simple function evaluation.
    TensorValue<double> FF_qp, VV;
    VectorValue<double> V, x_qp;
    double v;
    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        auto elem = *el_it;
        fe.reinit(elem);
        fe.collectDataForInterpolation(elem);
        fe.interpolate(elem);
        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            const libMesh::Point& X_qp = q_point[qp];
            const std::vector<double>& x_data = fe_interp_var_data[qp][X_sys_idx];
            const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
            get_x_and_FF(x_qp, FF_qp, x_data, grad_x_data);

            // Scalar-valued variables.
            for (unsigned int k = 0; k < num_scalar_vars; ++k)
            {
                fe.setInterpolatedDataPointers(
                    scalar_var_data[k], scalar_grad_var_data[k], scalar_var_fcn_system_idxs[k], elem, qp);
                scalar_var_dof_maps[k]->dof_indices(elem, scalar_var_dof_indices[k], 0);
                d_scalar_var_fcns[k](v,
                                     FF_qp,
                                     x_qp,
                                     X_qp,
                                     elem,
                                     scalar_var_data[k],
                                     scalar_grad_var_data[k],
                                     data_time,
                                     d_scalar_var_fcn_ctxs[k]);
                d_scalar_var_systems[k]->solution->set(scalar_var_dof_indices[k][0], v);
            }

            // Vector-valued variables.
            for (unsigned int k = 0; k < num_vector_vars; ++k)
            {
                fe.setInterpolatedDataPointers(
                    vector_var_data[k], vector_grad_var_data[k], vector_var_fcn_system_idxs[k], elem, qp);
                for (unsigned int i = 0; i < d_vector_var_dims[k]; ++i)
                {
                    vector_var_dof_maps[k]->dof_indices(elem, vector_var_dof_indices[k][i], i);
                }
                d_vector_var_fcns[k](V,
                                     FF_qp,
                                     x_qp,
                                     X_qp,
                                     elem,
                                     vector_var_data[k],
                                     vector_grad_var_data[k],
                                     data_time,
                                     d_vector_var_fcn_ctxs[k]);
                for (unsigned int i = 0; i < d_vector_var_dims[k]; ++i)
                {
                    d_vector_var_systems[k]->solution->set(vector_var_dof_indices[k][i][0], V(i));
                }
            }

            // Tensor-valued variables.
            for (unsigned int k = 0; k < num_tensor_vars; ++k)
            {
                fe.setInterpolatedDataPointers(
                    tensor_var_data[k], tensor_grad_var_data[k], tensor_var_fcn_system_idxs[k], elem, qp);
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
                                     x_qp,
                                     X_qp,
                                     elem,
                                     tensor_var_data[k],
                                     tensor_grad_var_data[k],
                                     data_time,
                                     d_tensor_var_fcn_ctxs[k]);
                for (unsigned int i = 0; i < d_tensor_var_dims[k]; ++i)
                {
                    for (unsigned int j = 0; j < d_tensor_var_dims[k]; ++j)
                    {
                        d_tensor_var_systems[k]->solution->set(tensor_var_dof_indices[k][i][j][0], VV(i, j));
                    }
                }
            }
        }
    }

    // Close all vectors.
    for (unsigned int k = 0; k < num_scalar_vars; ++k)
    {
        d_scalar_var_systems[k]->solution->close();
    }
    for (unsigned int k = 0; k < num_vector_vars; ++k)
    {
        d_vector_var_systems[k]->solution->close();
    }
    for (unsigned int k = 0; k < num_tensor_vars; ++k)
    {
        d_tensor_var_systems[k]->solution->close();
    }
    return;
} // reconstructVariables

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
