// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2021 by the IBAMR developers
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

#include "ibtk/IBTK_MPI.h"
#include <ibtk/FECache.h>
#include <ibtk/FEMapping.h>
#include <ibtk/FEValues.h>

#include <tbox/PIO.h>
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

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////
namespace
{
// Unfortunately libMesh::FE is templated on the finite element type, so we
// do not have a generic way to populate shape function values based on an
// FEType argument. Get around this by implementing our own template functon
// to do exactly that.
template <int dim, FEFamily fe_family>
void
fill_reference_values(const libMesh::QBase& quadrature,
                      const Order order,
                      boost::multi_array<double, 2>& reference_shape_values,
                      boost::multi_array<libMesh::VectorValue<double>, 2>& reference_shape_gradients)
{
    using FE = libMesh::FE<dim, fe_family>;
    const ElemType elem_type = quadrature.get_elem_type();
    const unsigned int n_shape_functions = FE::n_dofs(elem_type, order);

    boost::multi_array<double, 2>::extent_gen extents;
    reference_shape_values.resize(extents[n_shape_functions][quadrature.n_points()]);
    reference_shape_gradients.resize(extents[n_shape_functions][quadrature.n_points()]);

    // values:
    for (unsigned int node_n = 0; node_n < n_shape_functions; ++node_n)
    {
        for (unsigned int q = 0; q < quadrature.n_points(); ++q)
        {
            reference_shape_values[node_n][q] = FE::shape(elem_type, order, node_n, quadrature.qp(q));
        }
    }

    // gradients:
    for (unsigned int node_n = 0; node_n < n_shape_functions; ++node_n)
    {
        for (unsigned int q = 0; q < quadrature.n_points(); ++q)
        {
            for (unsigned int d = 0; d < dim; ++d)
            {
                reference_shape_gradients[node_n][q](d) =
                    FE::shape_deriv(elem_type, order, node_n, d, quadrature.qp(q));
            }
        }
    }
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

std::unique_ptr<FEValuesBase>
FEValuesBase::build(const int dim,
                    const int spacedim,
                    libMesh::QBase* qrule,
                    const libMesh::FEType fe_type,
                    const FEUpdateFlags update_flags)
{
    TBOX_ASSERT(dim <= spacedim);
    switch (dim)
    {
    case 1:
        switch (spacedim)
        {
        case 1:
            return std::unique_ptr<FEValuesBase>(new FEValues<1, 1>(qrule, fe_type, update_flags));
        case 2:
            return std::unique_ptr<FEValuesBase>(new FEValues<1, 2>(qrule, fe_type, update_flags));
        case 3:
            return std::unique_ptr<FEValuesBase>(new FEValues<1, 3>(qrule, fe_type, update_flags));
        default:
            break;
        }
        break;
    case 2:
        switch (spacedim)
        {
        case 2:
            return std::unique_ptr<FEValuesBase>(new FEValues<2, 2>(qrule, fe_type, update_flags));
        case 3:
            return std::unique_ptr<FEValuesBase>(new FEValues<2, 3>(qrule, fe_type, update_flags));
        default:
            break;
        }
        break;
    case 3:
        TBOX_ASSERT(spacedim == dim);
        return std::unique_ptr<FEValuesBase>(new FEValues<3, 3>(qrule, fe_type, update_flags));
    default:
        break;
    }

    // we shouldn't be able to get here
    TBOX_ERROR("FEValuesBase::build():\n"
               << "This function only supports dim and spacedim equal to 1, 2, "
               << "or 3 and dim <= spacedim." << std::endl);
    return {};
}

template <int dim, int spacedim>
FEValues<dim, spacedim>::FEValues(libMesh::QBase* qrule, const FEType fe_type, const FEUpdateFlags update_flags)
    : d_qrule(qrule), d_fe_type(fe_type), d_update_flags(update_flags)
{
    // set up update flag dependencies:
    if (d_update_flags & update_dphi) d_update_flags |= update_covariants;
}

template <int dim, int spacedim>
void
FEValues<dim, spacedim>::reinit(const libMesh::Elem* elem)
{
    // some things are not yet implemented
    TBOX_ASSERT(elem->p_level() == 0);
    // dim is only available at runtime with libMesh
    TBOX_ASSERT(elem->dim() == dim);
    // TODO - find a way to assert that the spatial dimension is right

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
        typename decltype(d_mappings)::value_type new_entry{ elem_type, nullptr };
        map_iter = d_mappings.insert(map_iter, std::move(new_entry));
        const quadrature_key_type key{
            elem_type, d_qrule->type(), d_qrule->get_order(), d_qrule->allow_rules_with_negative_weights
        };
        map_iter->second = FEMapping<dim, spacedim>::build(key, d_update_flags);
    }
    FEMapping<dim, spacedim>& mapping = *map_iter->second;
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
        ReferenceValues ref_values(*d_qrule, d_fe_type);
        ref_iter = d_reference_values.emplace(elem_type, std::move(ref_values)).first;
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

        const EigenAlignedVector<Eigen::Matrix<double, spacedim, dim> >& covariants = mapping.getCovariants();
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

                Eigen::Matrix<double, spacedim, 1> shape_grad_ = covariants[q] * ref_shape_grad_;
                for (unsigned int d = 0; d < spacedim; ++d)
                {
                    d_shape_gradients[i][q](d) = shape_grad_(d, 0);
                }
            }
        }
    }

    d_last_elem_type = elem_type;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

template <int dim, int spacedim>
FEValues<dim, spacedim>::ReferenceValues::ReferenceValues(const libMesh::QBase& quadrature,
                                                          const libMesh::FEType& fe_type)
    : d_elem_type(quadrature.get_elem_type()), d_fe_type(fe_type)
{
    const auto order = d_fe_type.order;

    // See the note in fill_reference_values explaining why this is necessary
    switch (fe_type.family)
    {
    case L2_LAGRANGE:
        fill_reference_values<dim, L2_LAGRANGE>(
            quadrature, order, d_reference_shape_values, d_reference_shape_gradients);
        break;
    case LAGRANGE:
        fill_reference_values<dim, LAGRANGE>(quadrature, order, d_reference_shape_values, d_reference_shape_gradients);
        break;
    case MONOMIAL:
        fill_reference_values<dim, MONOMIAL>(quadrature, order, d_reference_shape_values, d_reference_shape_gradients);
        break;
    case SCALAR:
        fill_reference_values<dim, SCALAR>(quadrature, order, d_reference_shape_values, d_reference_shape_gradients);
        break;
    default:
        TBOX_ERROR("unsupported element type");
    }
}

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

// instantiations
template class FEValues<1, 1>;
template class FEValues<1, 2>;
template class FEValues<1, 3>;
template class FEValues<2, 2>;
template class FEValues<2, 3>;
template class FEValues<3, 3>;
} // namespace IBTK

/////////////////////////////////////////////////////////////////////////////
