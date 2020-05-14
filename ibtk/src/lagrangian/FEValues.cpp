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

#include <ibtk/FECache.h>
#include <ibtk/FEValues.h>
#include <ibtk/JacobianCalculator.h>
#include <ibtk/namespaces.h> // IWYU pragma: keep

#include <tbox/PIO.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/Utilities.h>

#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/type_vector.h"
#include <libmesh/enum_elem_type.h>
#include <libmesh/fe.h>
#include <libmesh/point.h>
#include <libmesh/quadrature.h>

#include <map>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <int dim>
FEValues<dim>::FEValues(libMesh::QBase* qrule, const FEUpdateFlags update_flags)
    : d_qrule(qrule), d_update_flags(update_flags)
{
    // set up update flag dependencies:
    if (d_update_flags & update_dphi) d_update_flags |= update_covariants;
}

template <int dim>
void
FEValues<dim>::reinit(const libMesh::Elem* elem)
{
    // some things are not yet implemented
    TBOX_ASSERT(elem->p_level() == 0);

    const libMesh::ElemType elem_type = elem->type();
    // maybe update the quadrature rule:
    if (elem_type != d_last_elem_type)
    {
        d_qrule->init(elem_type, elem->p_level());
    }

    //
    // update mapping quantities:
    //
    auto map_iter = d_mappings.find(elem_type);
    if (map_iter == d_mappings.end())
    {
        typename decltype(d_mappings)::value_type new_entry {elem_type, nullptr};
        map_iter = d_mappings.insert(map_iter, std::move(new_entry));
        const std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order> key{ elem_type,
                                                                                          d_qrule->type(),
                                                                                          d_qrule->get_order() };
        map_iter->second = Mapping<dim>::build(key, d_update_flags);
    }
    Mapping<dim>& mapping = *map_iter->second;
    mapping.reinit(elem);

    if (d_update_flags & update_JxW)
    {
        d_JxW = mapping.getJxW();
    }
    if (d_update_flags & update_quadrature_points)
    {
        d_quadrature_points = mapping.getQuadraturePoints();
    }

    //
    // update shape function quantities:
    //
    auto ref_iter = d_reference_values.find(elem_type);
    if (ref_iter == d_reference_values.end())
    {
        ref_iter = d_reference_values.insert(ref_iter, { elem_type, *d_qrule });
    }
    const ReferenceValues& ref_values = ref_iter->second;

    if (d_last_elem_type != elem_type && d_update_flags & update_phi)
    {
        const boost::multi_array<double, 2>& ref_shape_values = ref_values.d_reference_shape_values;
        d_shape_values.resize(ref_shape_values.shape()[0]);
        for (unsigned int i = 0; i < d_shape_values.size(); ++i)
        {
            d_shape_values[i].resize(0);
            d_shape_values[i].insert(d_shape_values[i].begin(),
                                     &ref_shape_values[i][0],
                                     &ref_shape_values[i][0] + ref_shape_values.shape()[1]);
        }
    }

    if (d_update_flags & update_dphi)
    {
        const boost::multi_array<libMesh::VectorValue<double>, 2>& ref_shape_gradients =
            ref_values.d_reference_shape_gradients;
        d_shape_gradients.resize(ref_shape_gradients.shape()[0]);

        const EigenAlignedVector<Eigen::Matrix<double, dim, dim> >& covariants = mapping.getCovariants();
        for (unsigned int i = 0; i < d_shape_gradients.size(); ++i)
        {
            d_shape_gradients[i].resize(ref_shape_gradients.shape()[1]);
            for (unsigned int q = 0; q < d_shape_gradients[i].size(); ++q)
            {
                const libMesh::VectorValue<double>& ref_shape_grad = ref_shape_gradients[i][q];
                Eigen::Matrix<double, dim, 1> ref_shape_grad_;
                for (unsigned int d = 0; d < dim; ++d)
                {
                    ref_shape_grad_(d, 0) = ref_shape_grad(d);
                }

                Eigen::Matrix<double, dim, 1> shape_grad_ = covariants[q] * ref_shape_grad_;
                for (unsigned int d = 0; d < dim; ++d)
                {
                    d_shape_gradients[i][q](d) = shape_grad_(d, 0);
                }
            }
        }
    }

    d_last_elem_type = elem_type;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

template <int dim>
FEValues<dim>::ReferenceValues::ReferenceValues(const libMesh::QBase& quadrature)
    : d_elem_type(quadrature.get_elem_type())
{
    const unsigned int n_nodes = get_n_nodes(d_elem_type);
    const auto order = get_default_order(d_elem_type);

    typename decltype(d_reference_shape_values)::extent_gen extents;
    d_reference_shape_values.resize(extents[n_nodes][quadrature.n_points()]);
    d_reference_shape_gradients.resize(extents[n_nodes][quadrature.n_points()]);

    using FE = libMesh::FE<dim, libMesh::LAGRANGE>;
    // values:
    for (unsigned int node_n = 0; node_n < n_nodes; ++node_n)
    {
        for (unsigned int q = 0; q < quadrature.n_points(); ++q)
        {
            d_reference_shape_values[node_n][q] = FE::shape(d_elem_type, order, node_n, quadrature.qp(q));
        }
    }

    // gradients:
    for (unsigned int node_n = 0; node_n < n_nodes; ++node_n)
    {
        for (unsigned int q = 0; q < quadrature.n_points(); ++q)
        {
            for (unsigned int d = 0; d < dim; ++d)
            {
                d_reference_shape_gradients[node_n][q](d) =
                    FE::shape_deriv(d_elem_type, order, node_n, d, quadrature.qp(q));
            }
        }
    }
}

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

// instantiations
template class FEValues<1>;
template class FEValues<2>;
template class FEValues<3>;

} // namespace IBTK

/////////////////////////////////////////////////////////////////////////////
