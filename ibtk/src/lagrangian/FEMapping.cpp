// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "ibtk/FEMapping.h"
#include "ibtk/libmesh_utilities.h"

#include "tbox/Utilities.h"

#include "libmesh/dense_matrix.h"
#include "libmesh/elem.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/type_vector.h"
#include <libmesh/libmesh_version.h>
#include <libmesh/point.h>
#include <libmesh/quadrature.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <memory>

#include "ibtk/namespaces.h"

namespace IBTK
{
//
// Helper functions
//
template <int dim, int spacedim>
inline Eigen::Matrix<double, spacedim, dim>
getCovariant(const Eigen::Matrix<double, spacedim, dim>& contravariant)
{
    return contravariant * (contravariant.transpose() * contravariant).inverse();
}

template <int dim>
inline Eigen::Matrix<double, dim, dim>
getCovariant(const Eigen::Matrix<double, dim, dim>& contravariant)
{
    return contravariant.inverse().transpose();
}

//
// PointMap
//

template <int dim, int spacedim, int n_nodes>
PointMap<dim, spacedim, n_nodes>::PointMap(const libMesh::ElemType elem_type,
                                           const std::vector<libMesh::Point>& q_points)
    : d_reference_q_points(q_points)
{
    const int n_nodes_ = n_nodes == -1 ? static_cast<int>(get_n_nodes(elem_type)) : n_nodes;
    const auto elem_order = get_default_order(elem_type);
    d_phi.resize(n_nodes_, d_reference_q_points.size());
    for (int i = 0; i < n_nodes_; ++i)
    {
        for (unsigned int q = 0; q < d_reference_q_points.size(); ++q)
        {
            using FE = libMesh::FE<dim, libMesh::LAGRANGE>;
            d_phi(i, q) = FE::shape(elem_type, elem_order, i, d_reference_q_points[q]);
        }
    }
}

template <int dim, int spacedim, int n_nodes>
void
PointMap<dim, spacedim, n_nodes>::getMappedQuadraturePoints(const libMesh::Point* nodes,
                                                            const libMesh::Point* nodes_end,
                                                            std::vector<libMesh::Point>& physical_q_points)
{
    if (n_nodes != -1) TBOX_ASSERT(nodes_end - nodes == n_nodes);
    const int n_nodes_ = n_nodes == -1 ? nodes_end - nodes : n_nodes;
    TBOX_ASSERT(d_reference_q_points.size() == physical_q_points.size());
    TBOX_ASSERT(n_nodes_ == d_phi.rows());
    // assumes same node ordering in the input node array as is stored in d_phi
    for (unsigned int q = 0; q < d_reference_q_points.size(); ++q)
    {
        physical_q_points[q] = 0.0;
        for (int i = 0; i < n_nodes_; ++i)
        {
            for (int d = 0; d < spacedim; ++d)
            {
                physical_q_points[q](d) += d_phi(i, q) * nodes[i](d);
            }
        }
    }
}

//
// QuadratureData
//

QuadratureData::QuadratureData(const QuadratureData::key_type quad_key) : d_key(quad_key)
{
    const ElemType elem_type = std::get<0>(d_key);
    const QuadratureType quad_type = std::get<1>(d_key);
    const Order order = std::get<2>(d_key);

    const int dim = get_dim(elem_type);

    std::unique_ptr<QBase> quad_rule = QBase::build(quad_type, dim, order);
    quad_rule->init(elem_type);
    d_points = quad_rule->get_points();
    d_weights = quad_rule->get_weights();
}

//
// FEMapping
//

template <>
std::unique_ptr<FEMapping<2, 2> >
FEMapping<2, 2>::build(const key_type key, const FEUpdateFlags update_flags)
{
    switch (std::get<0>(key))
    {
    case libMesh::ElemType::TRI3:
        return std::unique_ptr<FEMapping<2, 2> >(new Tri3Mapping(key, update_flags));
    case libMesh::ElemType::TRI6:
        return std::unique_ptr<FEMapping<2, 2> >(new Tri6Mapping(key, update_flags));
    case libMesh::ElemType::QUAD4:
        return std::unique_ptr<FEMapping<2, 2> >(new Quad4Mapping(key, update_flags));
    case libMesh::ElemType::QUAD9:
        return std::unique_ptr<FEMapping<2, 2> >(new Quad9Mapping(key, update_flags));
    default:
        return std::unique_ptr<FEMapping<2, 2> >(new FELagrangeMapping<2, 2>(key, std::get<0>(key), update_flags));
    }

    return {};
}

template <>
std::unique_ptr<FEMapping<3, 3> >
FEMapping<3, 3>::build(const key_type key, const FEUpdateFlags update_flags)
{
    switch (std::get<0>(key))
    {
    case libMesh::ElemType::TET4:
        return std::unique_ptr<FEMapping<3, 3> >(new Tet4Mapping(key, update_flags));
    case libMesh::ElemType::TET10:
        return std::unique_ptr<FEMapping<3, 3> >(new Tet10Mapping(key, update_flags));
    case libMesh::ElemType::HEX8:
        return std::unique_ptr<FEMapping<3, 3> >(
            new FELagrangeMapping<3, 3, 8>(key, libMesh::ElemType::HEX8, update_flags));
    case libMesh::ElemType::HEX27:
        return std::unique_ptr<FEMapping<3, 3> >(new Hex27Mapping(key, update_flags));
    default:
        return std::unique_ptr<FEMapping<3, 3> >(new FELagrangeMapping<3, 3>(key, std::get<0>(key), update_flags));
    }

    return {};
}

template <int dim, int spacedim>
std::unique_ptr<FEMapping<dim, spacedim> >
FEMapping<dim, spacedim>::build(const key_type key, const FEUpdateFlags update_flags)
{
    return std::unique_ptr<FEMapping<dim, spacedim> >(
        new FELagrangeMapping<dim, spacedim>(key, std::get<0>(key), update_flags));
}

//
// FENodalMapping
//

template <int dim, int spacedim, int n_nodes>
FENodalMapping<dim, spacedim, n_nodes>::FENodalMapping(
    const typename FENodalMapping<dim, spacedim, n_nodes>::key_type quad_key,
    const libMesh::ElemType element_mapping_type,
    const FEUpdateFlags update_flags)
    : d_quadrature_data(quad_key), d_point_map(element_mapping_type, d_quadrature_data.d_points)
{
    d_update_flags = update_flags;

    // make sure dependencies are satisfied. These dependencies are only true
    // for Lagrange-type mappings.
    {
        if (d_update_flags & FEUpdateFlags::update_JxW) d_update_flags |= update_jacobians;

        if (d_update_flags & FEUpdateFlags::update_jacobians) d_update_flags |= update_contravariants;

        if (d_update_flags & FEUpdateFlags::update_covariants) d_update_flags |= update_contravariants;
    }

    if (d_update_flags & FEUpdateFlags::update_contravariants) d_contravariants.resize(this->d_quadrature_data.size());
    if (d_update_flags & FEUpdateFlags::update_covariants) d_covariants.resize(this->d_quadrature_data.size());

    if (d_update_flags & FEUpdateFlags::update_jacobians) d_Jacobians.resize(this->d_quadrature_data.size());

    if (d_update_flags & FEUpdateFlags::update_JxW) d_JxW.resize(this->d_quadrature_data.size());

    if (d_update_flags & FEUpdateFlags::update_quadrature_points)
        d_quadrature_points.resize(this->d_quadrature_data.size());
}

template <int dim, int spacedim, int n_nodes>
void
FENodalMapping<dim, spacedim, n_nodes>::reinit(const libMesh::Elem* elem)
{
    if (d_update_flags & FEUpdateFlags::update_contravariants || d_update_flags & FEUpdateFlags::update_covariants)
        this->fillTransforms(elem);
    if (d_update_flags & FEUpdateFlags::update_jacobians) this->fillJacobians();
    if (d_update_flags & FEUpdateFlags::update_JxW) this->fillJxW();
    if (d_update_flags & FEUpdateFlags::update_quadrature_points) this->fillQuadraturePoints(elem);
}

template <int dim, int spacedim, int n_nodes>
bool
FENodalMapping<dim, spacedim, n_nodes>::isAffine() const
{
    return false;
}

template <int dim, int spacedim, int n_nodes>
void
FENodalMapping<dim, spacedim, n_nodes>::fillJacobians()
{
    for (unsigned int q = 0; q < d_contravariants.size(); ++q)
    {
        if (dim == spacedim)
        {
            d_Jacobians[q] = d_contravariants[q].determinant();
        }
        else
        {
            Eigen::Matrix<double, dim, dim> Jac = d_contravariants[q].transpose() * d_contravariants[q];
            d_Jacobians[q] = std::sqrt(Jac.determinant());
        }
        TBOX_ASSERT(d_Jacobians[q] > 0.0);

        if (isAffine()) break;
    }

    if (isAffine()) std::fill(d_Jacobians.begin() + 1, d_Jacobians.end(), d_Jacobians[0]);

    return;
}

template <int dim, int spacedim, int n_nodes>
void
FENodalMapping<dim, spacedim, n_nodes>::fillJxW()
{
    for (unsigned int q = 0; q < d_Jacobians.size(); ++q) d_JxW[q] = d_quadrature_data.d_weights[q] * d_Jacobians[q];
}

template <int dim, int spacedim, int n_nodes>
void
FENodalMapping<dim, spacedim, n_nodes>::fillQuadraturePoints(const libMesh::Elem* elem)
{
    libMesh::Point nodes[27];

    // We occasionally (e.g., TET10 and TET4) want to call the lower-order
    // mapping from the higher-order mapping, so permit elements with more
    // nodes in that code
    if (n_nodes != -1) TBOX_ASSERT(n_nodes <= static_cast<int>(elem->n_nodes()));
    const int n_nodes_ = n_nodes == -1 ? static_cast<int>(elem->n_nodes()) : n_nodes;
    TBOX_ASSERT(n_nodes_ <= 27);
    for (int node_n = 0; node_n < n_nodes_; ++node_n)
    {
        nodes[node_n] = static_cast<const libMesh::Point&>(*elem->node_ptr(node_n));
    }

    d_point_map.getMappedQuadraturePoints(std::begin(nodes), std::begin(nodes) + n_nodes_, d_quadrature_points);
}

//
// FELagrangeMapping
//

template <int dim, int spacedim, int n_nodes>
FELagrangeMapping<dim, spacedim, n_nodes>::FELagrangeMapping(
    const typename FELagrangeMapping<dim, spacedim, n_nodes>::key_type quad_key,
    const libMesh::ElemType element_mapping_type,
    const FEUpdateFlags update_flags)
    : FENodalMapping<dim, spacedim, n_nodes>(quad_key, element_mapping_type, update_flags),
      d_n_nodes(n_nodes == -1 ? static_cast<int>(get_n_nodes(std::get<0>(quad_key))) : n_nodes)
{
    if (n_nodes != -1) TBOX_ASSERT(d_n_nodes == n_nodes);
#if LIBMESH_VERSION_LESS_THAN(1, 4, 0)
    TBOX_ASSERT(d_n_nodes <= 27);
#else
    TBOX_ASSERT(d_n_nodes <= static_cast<int>(libMesh::Elem::max_n_nodes));
#endif
    typename decltype(d_dphi)::extent_gen extents;
    d_dphi.resize(extents[d_n_nodes][this->d_quadrature_data.size()]);

    for (int node_n = 0; node_n < d_n_nodes; ++node_n)
    {
        for (unsigned int q = 0; q < this->d_quadrature_data.size(); ++q)
        {
            for (unsigned int d = 0; d < dim; ++d)
                d_dphi[node_n][q][d] =
                    libMesh::FE<dim, libMesh::LAGRANGE>::shape_deriv(element_mapping_type,
                                                                     get_default_order(element_mapping_type),
                                                                     node_n,
                                                                     d,
                                                                     this->d_quadrature_data.d_points[q]);
        }
    }
}

template <int dim, int spacedim, int n_nodes>
void
FELagrangeMapping<dim, spacedim, n_nodes>::fillTransforms(const libMesh::Elem* elem)
{
    TBOX_ASSERT(this->d_update_flags & FEUpdateFlags::update_contravariants);
    TBOX_ASSERT(d_n_nodes <= static_cast<int>(elem->n_nodes()));

    // max_n_nodes is a constant defined by libMesh - currently 27
#if LIBMESH_VERSION_LESS_THAN(1, 4, 0)
    double xs[27][spacedim];
#else
    double xs[libMesh::Elem::max_n_nodes][spacedim];
#endif

    const int n_nodes_ = n_nodes == -1 ? d_n_nodes : n_nodes;
    for (int i = 0; i < n_nodes_; ++i)
    {
        const libMesh::Point p = elem->point(i);
        for (unsigned int j = 0; j < spacedim; ++j) xs[i][j] = p(j);
    }

    for (unsigned int q = 0; q < this->d_contravariants.size(); ++q)
    {
        auto& contravariant = this->d_contravariants[q];
        contravariant.setZero();
        for (int node_n = 0; node_n < n_nodes_; ++node_n)
        {
            for (unsigned int i = 0; i < spacedim; ++i)
            {
                for (unsigned int j = 0; j < dim; ++j)
                {
                    contravariant(i, j) += xs[node_n][i] * d_dphi[node_n][q][j];
                }
            }
        }
    }

    if (this->d_update_flags & FEUpdateFlags::update_covariants)
    {
        for (unsigned int q = 0; q < this->d_contravariants.size(); ++q)
        {
            this->d_covariants[q] = getCovariant(this->d_contravariants[q]);
        }
    }

    return;
}

//
// Tri3Mapping
//

Tri3Mapping::Tri3Mapping(const key_type quad_key, const FEUpdateFlags update_flags)

    : FENodalMapping<2, 2, 3>(quad_key, libMesh::TRI3, update_flags)
{
}

void
Tri3Mapping::fillTransforms(const libMesh::Elem* elem)
{
    TBOX_ASSERT(this->d_update_flags & FEUpdateFlags::update_contravariants);
    // also permit TRI6
    const auto type = elem->type();
    TBOX_ASSERT(type == libMesh::TRI3 || type == libMesh::TRI6);

    const libMesh::Point p0 = elem->point(0);
    const libMesh::Point p1 = elem->point(1);
    const libMesh::Point p2 = elem->point(2);

    Eigen::Matrix<double, 2, 2> contravariant;
    contravariant(0, 0) = p1(0) - p0(0);
    contravariant(0, 1) = p2(0) - p0(0);
    contravariant(1, 0) = p1(1) - p0(1);
    contravariant(1, 1) = p2(1) - p0(1);
    std::fill(this->d_contravariants.begin(), this->d_contravariants.end(), contravariant);

    if (this->d_update_flags & FEUpdateFlags::update_covariants)
    {
        const Eigen::Matrix<double, 2, 2> covariant = getCovariant(contravariant);
        std::fill(this->d_covariants.begin(), this->d_covariants.end(), covariant);
    }

    return;
}

bool
Tri3Mapping::isAffine() const
{
    return true;
}

//
// Tri6Mapping
//

Tri6Mapping::Tri6Mapping(const key_type quad_key, const FEUpdateFlags update_flags)
    : FELagrangeMapping<2, 2, 6>(quad_key, libMesh::ElemType::TRI6, update_flags), tri3_mapping(quad_key, update_flags)
{
}

void
Tri6Mapping::reinit(const libMesh::Elem* elem)
{
    if (elem_is_affine(elem))
    {
        tri3_mapping.reinit(elem);
        // If we ever add more fields to the mapping classes we will need to
        // duplicate them here
        std::swap(d_contravariants, tri3_mapping.d_contravariants);
        std::swap(d_covariants, tri3_mapping.d_covariants);
        std::swap(d_Jacobians, tri3_mapping.d_Jacobians);
        std::swap(d_JxW, tri3_mapping.d_JxW);
        std::swap(d_quadrature_points, tri3_mapping.d_quadrature_points);
    }
    else
        FELagrangeMapping<2, 2, 6>::reinit(elem);
}

bool
Tri6Mapping::elem_is_affine(const libMesh::Elem* elem)
{
    std::array<libMesh::Point, 6> nodes;
    for (unsigned int n = 0; n < nodes.size(); ++n) nodes[n] = elem->node_ref(n);

    // try to determine the size of the coordinates to use as the tolerance.
    double characteristic_point_size = 0.0;
    for (int d = 0; d < LIBMESH_DIM; ++d)
    {
        characteristic_point_size += std::abs(nodes[0](d));
        characteristic_point_size += std::abs(nodes[4](d));
    }
    const double tol = 1e-16 * characteristic_point_size;

    return nodes[3].absolute_fuzzy_equals(0.5 * (nodes[0] + nodes[1]), tol) &&
           nodes[4].absolute_fuzzy_equals(0.5 * (nodes[1] + nodes[2]), tol) &&
           nodes[5].absolute_fuzzy_equals(0.5 * (nodes[0] + nodes[2]), tol);
}

//
// Quad4Mapping
//

Quad4Mapping::Quad4Mapping(const key_type quad_key, const FEUpdateFlags update_flags)

    : FENodalMapping<2, 2, 4>(quad_key, libMesh::QUAD4, update_flags)
{
}

void
Quad4Mapping::fillTransforms(const libMesh::Elem* elem)
{
    TBOX_ASSERT(this->d_update_flags & FEUpdateFlags::update_contravariants);
    TBOX_ASSERT(elem->type() == std::get<0>(this->d_quadrature_data.d_key));

    // calculate constants in Jacobians here
    const libMesh::Point p0 = elem->point(0);
    const libMesh::Point p1 = elem->point(1);
    const libMesh::Point p2 = elem->point(2);
    const libMesh::Point p3 = elem->point(3);

    const double a_1 = 0.25 * (-p0(0) + p1(0) + p2(0) - p3(0));
    const double b_1 = 0.25 * (-p0(0) - p1(0) + p2(0) + p3(0));
    const double c_1 = 0.25 * (p0(0) - p1(0) + p2(0) - p3(0));
    const double a_2 = 0.25 * (-p0(1) + p1(1) + p2(1) - p3(1));
    const double b_2 = 0.25 * (-p0(1) - p1(1) + p2(1) + p3(1));
    const double c_2 = 0.25 * (p0(1) - p1(1) + p2(1) - p3(1));

    for (unsigned int i = 0; i < this->d_contravariants.size(); i++)
    {
        // calculate Jacobians here
        const double x = d_quadrature_data.d_points[i](0);
        const double y = d_quadrature_data.d_points[i](1);

        Eigen::Matrix<double, 2, 2>& contravariant = d_contravariants[i];
        contravariant(0, 0) = a_1 + c_1 * y;
        contravariant(0, 1) = b_1 + c_1 * x;
        contravariant(1, 0) = a_2 + c_2 * y;
        contravariant(1, 1) = b_2 + c_2 * x;
    }

    if (this->d_update_flags & FEUpdateFlags::update_covariants)
    {
        for (unsigned int q = 0; q < this->d_contravariants.size(); ++q)
        {
            const auto& contravariant = this->d_contravariants[q];
            d_covariants[q] = getCovariant(contravariant);
        }
    }

    return;
}

//
// Quad9Mapping
//

Quad9Mapping::Quad9Mapping(const Quad9Mapping::key_type quad_key, FEUpdateFlags update_flags)
    : FENodalMapping<2, 2, 9>(quad_key, std::get<0>(quad_key), update_flags)
{
    // This code utilizes an implementation detail of
    // QBase::tensor_product_quad where the x coordinate increases fastest to
    // reconstruct the 1D quadrature rule
    d_n_oned_q_points = static_cast<std::size_t>(std::round(std::sqrt(d_quadrature_data.d_points.size())));
    TBOX_ASSERT(d_n_oned_q_points * d_n_oned_q_points == d_quadrature_data.d_points.size());
    std::vector<libMesh::Point> oned_points(d_n_oned_q_points);
    for (unsigned int q = 0; q < d_n_oned_q_points; ++q)
    {
        oned_points[q] = d_quadrature_data.d_points[q](0);
    }

    // verify that we really do have a tensor product rule
    unsigned int q = 0;
    for (unsigned int j = 0; j < d_n_oned_q_points; ++j)
    {
        for (unsigned int i = 0; i < d_n_oned_q_points; ++i)
        {
            TBOX_ASSERT(d_quadrature_data.d_points[q] == libMesh::Point(oned_points[i](0), oned_points[j](0)));
            ++q;
        }
    }

    d_phi.resize(3, d_n_oned_q_points);
    d_dphi.resize(3, d_n_oned_q_points);
    for (unsigned int i = 0; i < 3u; ++i)
    {
        for (unsigned int q = 0; q < oned_points.size(); ++q)
        {
            // This class orders the vertices in a different way to make
            // writing tensor products easier: we do a left-to-right ordering
            // 0 - 1 - 2 instead of 0 - 2 - 1.
            constexpr int ibamr_to_libmesh_ordering[3] = { 0, 2, 1 };
            using FE = libMesh::FE<1, libMesh::LAGRANGE>;
            d_phi(i, q) = FE::shape(libMesh::EDGE3, libMesh::SECOND, ibamr_to_libmesh_ordering[i], oned_points[q]);
            d_dphi(i, q) =
                FE::shape_deriv(libMesh::EDGE3, libMesh::SECOND, ibamr_to_libmesh_ordering[i], 0, oned_points[q]);
        }
    }
}

void
Quad9Mapping::fillTransforms(const libMesh::Elem* elem)
{
    TBOX_ASSERT(this->d_update_flags & FEUpdateFlags::update_contravariants);
    TBOX_ASSERT(elem->type() == std::get<0>(this->d_quadrature_data.d_key));

    constexpr std::size_t n_oned_shape_functions = 3;

    // We index points in the following way:
    //
    // i = 2 +--+--+
    //       |     |
    // i = 1 +  +  +
    //       |     |
    // i = 0 +--+--+
    //      j=0 1  2
    //
    // All 2D arrays created here are in row-major order.

    double xs[n_oned_shape_functions][n_oned_shape_functions];
    double ys[n_oned_shape_functions][n_oned_shape_functions];

    libMesh::Point points[n_oned_shape_functions][n_oned_shape_functions];
    points[0][0] = elem->point(0);
    points[0][1] = elem->point(4);
    points[0][2] = elem->point(1);
    points[1][0] = elem->point(7);
    points[1][1] = elem->point(8);
    points[1][2] = elem->point(5);
    points[2][0] = elem->point(3);
    points[2][1] = elem->point(6);
    points[2][2] = elem->point(2);

    // j is the x index, i is the y index
    for (unsigned int i = 0; i < n_oned_shape_functions; ++i)
    {
        for (unsigned int j = 0; j < n_oned_shape_functions; ++j)
        {
            xs[i][j] = points[i][j](0);
            ys[i][j] = points[i][j](1);
        }
    }

    for (unsigned int q = 0; q < this->d_contravariants.size(); q++)
    {
        Eigen::Matrix<double, 2, 2>& contravariant = d_contravariants[q];
        contravariant.setZero();
        // Exploit the fact that Quad9 is a tensor product element by indexing
        // the x component of each tensor product shape function with j and
        // the y component with i
        const unsigned int q_point_x = q % d_n_oned_q_points;
        const unsigned int q_point_y = q / d_n_oned_q_points;
        for (unsigned int i = 0; i < n_oned_shape_functions; ++i)
        {
            for (unsigned int j = 0; j < n_oned_shape_functions; ++j)
            {
                contravariant(0, 0) += xs[i][j] * d_dphi(j, q_point_x) * d_phi(i, q_point_y);
                contravariant(0, 1) += xs[i][j] * d_phi(j, q_point_x) * d_dphi(i, q_point_y);
                contravariant(1, 0) += ys[i][j] * d_dphi(j, q_point_x) * d_phi(i, q_point_y);
                contravariant(1, 1) += ys[i][j] * d_phi(j, q_point_x) * d_dphi(i, q_point_y);
            }
        }
    }

    if (this->d_update_flags & FEUpdateFlags::update_covariants)
    {
        for (unsigned int q = 0; q < this->d_contravariants.size(); ++q)
        {
            const auto& contravariant = this->d_contravariants[q];
            d_covariants[q] = getCovariant(contravariant);
        }
    }

    return;
}

//
// Tet4Mapping
//

Tet4Mapping::Tet4Mapping(const key_type quad_key, const FEUpdateFlags update_flags)

    : FENodalMapping<3, 3, 4>(quad_key, libMesh::TET4, update_flags)
{
}

void
Tet4Mapping::fillTransforms(const libMesh::Elem* elem)
{
    TBOX_ASSERT(this->d_update_flags & FEUpdateFlags::update_contravariants);
    // also permit TET10
    const auto type = elem->type();
    TBOX_ASSERT(type == libMesh::TET4 || type == libMesh::TET10);

    // calculate Jacobians here
    const libMesh::Point p0 = elem->point(0);
    const libMesh::Point p1 = elem->point(1);
    const libMesh::Point p2 = elem->point(2);
    const libMesh::Point p3 = elem->point(3);

    Eigen::Matrix<double, 3, 3> contravariant;
    contravariant(0, 0) = p1(0) - p0(0);
    contravariant(0, 1) = p2(0) - p0(0);
    contravariant(0, 2) = p3(0) - p0(0);
    contravariant(1, 0) = p1(1) - p0(1);
    contravariant(1, 1) = p2(1) - p0(1);
    contravariant(1, 2) = p3(1) - p0(1);
    contravariant(2, 0) = p1(2) - p0(2);
    contravariant(2, 1) = p2(2) - p0(2);
    contravariant(2, 2) = p3(2) - p0(2);
    std::fill(d_contravariants.begin(), d_contravariants.end(), contravariant);

    if (this->d_update_flags & FEUpdateFlags::update_covariants)
    {
        const Eigen::Matrix<double, 3, 3> covariant = getCovariant(contravariant);
        std::fill(this->d_covariants.begin(), this->d_covariants.end(), covariant);
    }

    return;
}

bool
Tet4Mapping::isAffine() const
{
    return true;
}

//
// Tet10Mapping
//

Tet10Mapping::Tet10Mapping(const key_type quad_key, const FEUpdateFlags update_flags)
    : FELagrangeMapping<3, 3, 10>(quad_key, libMesh::ElemType::TET10, update_flags),
      tet4_mapping(quad_key, update_flags)
{
}

void
Tet10Mapping::reinit(const libMesh::Elem* elem)
{
    if (elem_is_affine(elem))
    {
        tet4_mapping.reinit(elem);
        // If we ever add more fields to the mapping classes we will need to
        // duplicate them here
        std::swap(d_contravariants, tet4_mapping.d_contravariants);
        std::swap(d_covariants, tet4_mapping.d_covariants);
        std::swap(d_Jacobians, tet4_mapping.d_Jacobians);
        std::swap(d_JxW, tet4_mapping.d_JxW);
        std::swap(d_quadrature_points, tet4_mapping.d_quadrature_points);
    }
    else
        FELagrangeMapping<3, 3, 10>::reinit(elem);
}

bool
Tet10Mapping::elem_is_affine(const libMesh::Elem* elem)
{
    std::array<libMesh::Point, 10> nodes;
    for (unsigned int n = 0; n < nodes.size(); ++n) nodes[n] = elem->node_ref(n);

    // try to determine the size of the coordinates to use as the tolerance.
    double characteristic_point_size = 0.0;
    for (int d = 0; d < LIBMESH_DIM; ++d)
    {
        characteristic_point_size += std::abs(nodes[0](d));
        characteristic_point_size += std::abs(nodes[3](d));
    }
    const double tol = 1e-16 * characteristic_point_size;

    return nodes[4].absolute_fuzzy_equals(0.5 * (nodes[0] + nodes[1]), tol) &&
           nodes[5].absolute_fuzzy_equals(0.5 * (nodes[1] + nodes[2]), tol) &&
           nodes[6].absolute_fuzzy_equals(0.5 * (nodes[0] + nodes[2]), tol) &&
           nodes[7].absolute_fuzzy_equals(0.5 * (nodes[0] + nodes[3]), tol) &&
           nodes[8].absolute_fuzzy_equals(0.5 * (nodes[1] + nodes[3]), tol) &&
           nodes[9].absolute_fuzzy_equals(0.5 * (nodes[2] + nodes[3]), tol);
}

//
// Hex27Mapping
//

Hex27Mapping::Hex27Mapping(const key_type quad_key, const FEUpdateFlags update_flags)
    : FELagrangeMapping<3, 3, 27>(quad_key, libMesh::HEX27, update_flags),
      hex8_mapping(quad_key, libMesh::HEX8, update_flags)
{
}

void
Hex27Mapping::reinit(const libMesh::Elem* elem)
{
    if (elem_is_trilinear(elem))
    {
        hex8_mapping.reinit(elem);
        // If we ever add more fields to the mapping classes we will need to
        // duplicate them here
        std::swap(d_contravariants, hex8_mapping.d_contravariants);
        std::swap(d_covariants, hex8_mapping.d_covariants);
        std::swap(d_Jacobians, hex8_mapping.d_Jacobians);
        std::swap(d_JxW, hex8_mapping.d_JxW);
        std::swap(d_quadrature_points, hex8_mapping.d_quadrature_points);
    }
    else
        FELagrangeMapping<3, 3, 27>::reinit(elem);
}

bool
Hex27Mapping::elem_is_trilinear(const libMesh::Elem* elem)
{
    std::array<libMesh::Point, 27> nodes;
    for (unsigned int n = 0; n < nodes.size(); ++n) nodes[n] = elem->node_ref(n);

    // try to determine the size of the coordinates to use as the tolerance.
    double characteristic_point_size = 0.0;
    for (int d = 0; d < LIBMESH_DIM; ++d)
    {
        characteristic_point_size += std::abs(nodes[0](d));
        characteristic_point_size += std::abs(nodes[6](d));
    }
    const double tol = 1e-16 * characteristic_point_size;

    return
        // line midpoints
        nodes[8].absolute_fuzzy_equals(0.5 * (nodes[0] + nodes[1]), tol) &&
        nodes[9].absolute_fuzzy_equals(0.5 * (nodes[1] + nodes[2]), tol) &&
        nodes[10].absolute_fuzzy_equals(0.5 * (nodes[2] + nodes[3]), tol) &&
        nodes[11].absolute_fuzzy_equals(0.5 * (nodes[0] + nodes[3]), tol) &&
        nodes[12].absolute_fuzzy_equals(0.5 * (nodes[0] + nodes[4]), tol) &&
        nodes[13].absolute_fuzzy_equals(0.5 * (nodes[1] + nodes[5]), tol) &&
        nodes[14].absolute_fuzzy_equals(0.5 * (nodes[2] + nodes[6]), tol) &&
        nodes[15].absolute_fuzzy_equals(0.5 * (nodes[3] + nodes[7]), tol) &&
        nodes[16].absolute_fuzzy_equals(0.5 * (nodes[4] + nodes[5]), tol) &&
        nodes[17].absolute_fuzzy_equals(0.5 * (nodes[5] + nodes[6]), tol) &&
        nodes[18].absolute_fuzzy_equals(0.5 * (nodes[6] + nodes[7]), tol) &&
        nodes[19].absolute_fuzzy_equals(0.5 * (nodes[4] + nodes[7]), tol) &&
        // face midpoints
        nodes[20].absolute_fuzzy_equals(0.5 * (nodes[8] + nodes[10]), tol) &&
        nodes[21].absolute_fuzzy_equals(0.5 * (nodes[8] + nodes[16]), tol) &&
        nodes[22].absolute_fuzzy_equals(0.5 * (nodes[9] + nodes[17]), tol) &&
        nodes[23].absolute_fuzzy_equals(0.5 * (nodes[10] + nodes[18]), tol) &&
        nodes[24].absolute_fuzzy_equals(0.5 * (nodes[11] + nodes[19]), tol) &&
        nodes[25].absolute_fuzzy_equals(0.5 * (nodes[16] + nodes[18]), tol) &&
        // cell center
        nodes[26].absolute_fuzzy_equals(0.5 * (nodes[20] + nodes[25]), tol);
}

//
// Instantiations
//

template class FEMapping<1, 1>;
template class FEMapping<1, 2>;
template class FEMapping<1, 3>;
template class FEMapping<2, 2>;
template class FEMapping<2, 3>;
template class FEMapping<3, 3>;

template class FENodalMapping<1, 1>;
template class FENodalMapping<1, 2>;
template class FENodalMapping<1, 3>;
template class FENodalMapping<2, 2>;
template class FENodalMapping<2, 3>;
template class FENodalMapping<3, 3>;

template class FELagrangeMapping<1, 1>;
template class FELagrangeMapping<1, 2>;
template class FELagrangeMapping<1, 3>;
template class FELagrangeMapping<2, 2>;
template class FELagrangeMapping<2, 3>;
template class FELagrangeMapping<3, 3>;

template class FELagrangeMapping<2, 2, 6>;
template class FELagrangeMapping<3, 3, 10>;

template class FENodalMapping<2, 2, 3>;
template class FENodalMapping<2, 2, 4>;
template class FENodalMapping<2, 2, 6>;
template class FENodalMapping<2, 2, 9>;
template class FENodalMapping<3, 3, 4>;

} // namespace IBTK
