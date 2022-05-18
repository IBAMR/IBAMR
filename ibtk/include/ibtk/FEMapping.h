// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_FEMapping
#define included_IBTK_FEMapping

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#ifdef IBTK_HAVE_LIBMESH

#include <ibtk/FECache.h>
#include <ibtk/ibtk_utilities.h>

#include "tbox/Utilities.h"

#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_order.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/fe.h>
#include <libmesh/point.h>

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
IBTK_ENABLE_EXTRA_WARNINGS

#include <array>
#include <tuple>
#include <vector>

namespace IBTK
{
class Tri6Mapping;
class Tet10Mapping;
class Hex27Mapping;
} // namespace IBTK

namespace IBTK
{
/**
 * Internal class that computes mapped quadrature point locations for
 * Lagrange-type interpolatory elements.
 *
 * @tparam dim Logical dimension of the mesh.
 *
 * @tparam spacedim Spatial dimension of the mesh (i.e., nodes have spacedim
 * meaningful coordinates).
 *
 * @tparam n_nodes Number of nodes on an element. Defaults to -1, meaning a
 * run-time calculation of the number of nodes. This template parameter is
 * useful for first-order elements since the number is small and providing it
 * improves performance.
 */
template <int dim, int spacedim = dim, int n_nodes = -1>
class PointMap
{
public:
    PointMap(const libMesh::ElemType elem_type, const std::vector<libMesh::Point>& q_points);

    /**
     * Calculate mapped quadrature points.
     */
    void getMappedQuadraturePoints(const libMesh::Point* begin,
                                   const libMesh::Point* end,
                                   std::vector<libMesh::Point>& physical_q_points);

protected:
    /**
     * Quadrature points on the reference element.
     */
    std::vector<libMesh::Point> d_reference_q_points;

    /**
     * Table containing the values of 1D shape functions (which, with a tensor
     * product, define the mapping) at reference quadrature points.
     */
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> d_phi;
};

/*!
 * Abstract class defining the interface to a finite element mapping.
 */
template <int dim, int spacedim = dim>
class FEMapping
{
public:
    /*!
     * Recalculate relevant quantities for the provided element.
     */
    virtual void reinit(const libMesh::Elem* elem) = 0;

    /*!
     * Get the current jacobian times quadrature weight (JxW) values.
     */
    virtual const std::vector<double>& getJxW() const = 0;

    /*!
     * Get the positions of the quadrature points on the current element.
     */
    virtual const std::vector<libMesh::Point>& getQuadraturePoints() const = 0;

    /*!
     * Get the contravariants.
     */
    virtual const EigenAlignedVector<Eigen::Matrix<double, spacedim, dim> >& getContravariants() const = 0;

    /*!
     * Get the covariants.
     */
    virtual const EigenAlignedVector<Eigen::Matrix<double, spacedim, dim> >& getCovariants() const = 0;

    /*!
     * Standard 'quadrature key' alias - all the information to completely
     * define a libMesh quadrature rule.
     */
    using key_type = quadrature_key_type;

    /*!
     * Return a pointer to the correct mapping for a given quadrature key and
     * update flags object.
     */
    static std::unique_ptr<FEMapping<dim, spacedim> > build(const key_type key, const FEUpdateFlags update_flags);

    virtual ~FEMapping<dim, spacedim>() = default;

protected:
    /*!
     * Compute the contravariants and covariants. In general each mapping will
     * have to overload this function.
     */
    virtual void fillTransforms(const libMesh::Elem* elem) = 0;

    /*!
     * Compute determinants of contravariants (the Jacobians).
     */
    virtual void fillJacobians() = 0;

    /*!
     * Compute JxW values.
     */
    virtual void fillJxW() = 0;

    /*!
     * Compute the positions of quadrature points on the current element.
     */
    virtual void fillQuadraturePoints(const libMesh::Elem* elem) = 0;
};

/*!
 * Helper class that simply stores the relevant quadrature information extracted
 * from a quadrature key.
 */
struct QuadratureData
{
    /*!
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = quadrature_key_type;

    /*!
     * Constructor.
     */
    QuadratureData(const key_type quad_key);

    /*!
     *  Quadrature key.
     */
    const key_type d_key;

    /*!
     * Quadrature points on the reference element.
     */
    std::vector<libMesh::Point> d_points;

    /*!
     * Quadrature weights.
     */
    std::vector<double> d_weights;

    /*!
     * Get the size (the number of points) of the quadrature rule.
     */
    std::size_t size() const
    {
        return d_points.size();
    }
};

/*!
 * Base class for all nodal finite element mappings (i.e., mappings
 * corresponding to Lagrange-type finite element spaces).
 *
 * @tparam n_nodes Number of nodes of the element: defaults to runtime
 * calculation (-1).
 */
template <int dim, int spacedim = dim, int n_nodes = -1>
class FENodalMapping : public FEMapping<dim, spacedim>
{
public:
    /*!
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = quadrature_key_type;

    /*!
     * Constructor.
     *
     * @param[in] quad_key The quadrature key (i.e., a complete description of
     * the quadrature rule).
     * @param[in] mapping_element_type The element type used to compute the
     * mapping from the reference element to the physical element. This may be
     * different from the element type in the quadrature rule - for example,
     * one could provide TRI6 in the quadrature rule and TRI3 here.
     * @param[in] update_flags An enum describing which values need to be
     * computed on each element.
     */
    FENodalMapping(const key_type quad_key,
                   const libMesh::ElemType mapping_element_type,
                   const FEUpdateFlags update_flags);

    /*!
     * Recalculate relevant quantities for the provided element.
     */
    virtual void reinit(const libMesh::Elem* elem) override;

    virtual const std::vector<double>& getJxW() const override
    {
#ifndef NDEBUG
        TBOX_ASSERT(d_update_flags & FEUpdateFlags::update_JxW);
#endif
        return d_JxW;
    }

    virtual const std::vector<libMesh::Point>& getQuadraturePoints() const override
    {
#ifndef NDEBUG
        TBOX_ASSERT(d_update_flags & FEUpdateFlags::update_quadrature_points);
#endif
        return d_quadrature_points;
    }

    virtual const EigenAlignedVector<Eigen::Matrix<double, spacedim, dim> >& getContravariants() const override
    {
#ifndef NDEBUG
        TBOX_ASSERT(d_update_flags & FEUpdateFlags::update_contravariants);
#endif
        return d_contravariants;
    }

    virtual const EigenAlignedVector<Eigen::Matrix<double, spacedim, dim> >& getCovariants() const override
    {
#ifndef NDEBUG
        TBOX_ASSERT(d_update_flags & FEUpdateFlags::update_covariants);
#endif
        return d_covariants;
    }

protected:
    /*!
     * Information on the relevant quadrature rule.
     */
    const QuadratureData d_quadrature_data;

    /*!
     * Computed update flags for the mapping.
     */
    FEUpdateFlags d_update_flags;

    /*!
     * Array of contravariants.
     */
    EigenAlignedVector<Eigen::Matrix<double, spacedim, dim> > d_contravariants;

    /*!
     * Array of covariants (i.e., the transpose of the inverse of the
     * covariants when dim == spacedim)
     */
    EigenAlignedVector<Eigen::Matrix<double, spacedim, dim> > d_covariants;

    /*!
     * Array of Jacobians.
     */
    std::vector<double> d_Jacobians;

    /*!
     * Array of JxW values.
     */
    std::vector<double> d_JxW;

    /*!
     * Array of mapped quadrature points.
     */
    std::vector<libMesh::Point> d_quadrature_points;

    /*!
     * Object that computes quadrature point locations. This is sufficiently
     * different from the rest of the mapping code that it is implemented in
     * another class.
     */
    PointMap<dim, spacedim, n_nodes> d_point_map;

    /*!
     * Boolean indicating that the mapping is affine - if it is we can skip
     * some computations. The default implementation returns false. Inheriting
     * classes should overload this they represent affine mappings.
     */
    virtual bool isAffine() const;

    /*!
     * Compute determinants of contravariants (the Jacobians). The default
     * implementation given here is usually the correct one.
     */
    virtual void fillJacobians() override;

    /*!
     * Compute JxW values. The default implementation given here is usually the
     * correct one.
     */
    virtual void fillJxW() override;

    /*!
     * Compute the positions of quadrature points on the current element.
     */
    virtual void fillQuadraturePoints(const libMesh::Elem* elem) override;
};

/*!
 * A generic implementation for Lagrange-type elements: works for all elements
 * in that family but is less efficient than the specialized classes for
 * lower-order or tensor-product elements. Supports nonzero codimension.
 */
template <int dim, int spacedim = dim, int n_nodes = -1>
class FELagrangeMapping : public FENodalMapping<dim, spacedim, n_nodes>
{
public:
    /**
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = quadrature_key_type;

    /*!
     * Constructor.
     *
     * @param[in] quad_key The quadrature key (i.e., a complete description of
     * the quadrature rule).
     * @param[in] mapping_element_type The element type used to compute the
     * mapping from the reference element to the physical element. This may be
     * different from the element type in the quadrature rule - for example,
     * one could provide TRI6 in the quadrature rule and TRI3 here.
     * @param[in] update_flags An enum describing which values need to be
     * computed on each element.
     */
    FELagrangeMapping(const key_type quad_key,
                      const libMesh::ElemType mapping_element_type,
                      const FEUpdateFlags update_flags);

protected:
    virtual void fillTransforms(const libMesh::Elem* elem) override;

    /**
     * Number of nodes for the considered element type.
     */
    const int d_n_nodes;

    /**
     * Values of shape function gradients on the reference element at
     * quadrature points.
     */
    boost::multi_array<std::array<double, dim>, 2> d_dphi;

    friend class Hex27Mapping;
};

/*!
 * Specialization for TRI3 elements with codimension zero.
 */
class Tri3Mapping : public FENodalMapping<2, 2, 3>
{
public:
    /*!
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = quadrature_key_type;

    /*!
     * Constructor.
     */
    Tri3Mapping(const key_type quad_key, const FEUpdateFlags update_flags);

protected:
    virtual void fillTransforms(const libMesh::Elem* elem) override;

    virtual bool isAffine() const override;

    friend class Tri6Mapping;
};

/*!
 * Specialization for QUAD4 elements with codimension zero.
 */
class Quad4Mapping : public FENodalMapping<2, 2, 4>
{
public:
    /*!
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = quadrature_key_type;

    /*!
     * Constructor.
     */
    Quad4Mapping(const key_type quad_key, const FEUpdateFlags update_flags);

protected:
    virtual void fillTransforms(const libMesh::Elem* elem) override;
};

/*!
 * Specialization for QUAD9 elements with codimension zero.
 */
class Quad9Mapping : public FENodalMapping<2, 2, 9>
{
public:
    /**
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = quadrature_key_type;

    /**
     * Constructor.
     */
    Quad9Mapping(const key_type quad_key, const FEUpdateFlags update_flags);

protected:
    virtual void fillTransforms(const libMesh::Elem* elem) override;

    /**
     * Number of 1D quadrature points in the rule used to generate the 2D
     * tensor product rule.
     */
    std::size_t d_n_oned_q_points;

    /**
     * Table containing the values of 1D shape functions (which, with a tensor
     * product, define the mapping) at reference quadrature points.
     */
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> d_phi;

    /**
     * Table containing the derivatives of 1D shape functions (which, with a
     * tensor product, define the mapping) at reference quadrature points.
     */
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> d_dphi;
};

/*!
 * Specialization for TRI6 elements with codimension zero.
 */
class Tri6Mapping : public FELagrangeMapping<2, 2, 6>
{
public:
    /*!
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = quadrature_key_type;

    /*!
     * Constructor.
     */
    Tri6Mapping(const key_type quad_key, const FEUpdateFlags update_flags);

    virtual void reinit(const libMesh::Elem* elem) override;

protected:
    /*!
     * TRI3 mapping that is used whenever the given elem is affine.
     */
    Tri3Mapping tri3_mapping;

    /*!
     * Utility function that determines if the element is affine (i.e., all
     * nodes at edge midpoints are averages of corners)
     */
    static bool elem_is_affine(const libMesh::Elem* elem);
};

/*!
 * Specialization for TET4 elements.
 */
class Tet4Mapping : public FENodalMapping<3, 3, 4>
{
public:
    /*!
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = quadrature_key_type;

    /**
     * Constructor.
     */
    Tet4Mapping(const key_type quad_key, const FEUpdateFlags update_flags);

protected:
    virtual void fillTransforms(const libMesh::Elem* elem) override;

    virtual bool isAffine() const override;

    friend class Tet10Mapping;
};

/*!
 * Specialization for TET10 elements. Since, for most applications and in the
 * reference configuration, most TET10 elements are actually affine this class
 * tries use the TET4 mapping whenever possible.
 */
class Tet10Mapping : public FELagrangeMapping<3, 3, 10>
{
public:
    /*!
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = quadrature_key_type;

    /*!
     * Constructor.
     */
    Tet10Mapping(const key_type quad_key, const FEUpdateFlags update_flags);

    virtual void reinit(const libMesh::Elem* elem) override;

protected:
    /*!
     * TET4 mapping that is used whenever the given elem is affine.
     */
    Tet4Mapping tet4_mapping;

    /*!
     * Utility function that determines if the element is affine (i.e., all
     * nodes at edge midpoints are averages of corners)
     */
    static bool elem_is_affine(const libMesh::Elem* elem);
};

/*!
 * Specialization for HEX27 elements. Since, for most applications and in the
 * reference configuration, most HEX27 elements are actually trilinear, this
 * class tries use the lower degree mapping whenever possible.
 */
class Hex27Mapping : public FELagrangeMapping<3, 3, 27>
{
public:
    /*!
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = quadrature_key_type;

    /*!
     * Constructor.
     */
    Hex27Mapping(const key_type quad_key, const FEUpdateFlags update_flags);

    virtual void reinit(const libMesh::Elem* elem) override;

protected:
    /*!
     * HEX8 mapping that is used whenever the given elem is trilinear.
     */
    FELagrangeMapping<3, 3, 8> hex8_mapping;

    /*!
     * Utility function that determines if the element is trilinear (i.e., all
     * nodes at edge midpoints are averages of corners)
     */
    static bool elem_is_trilinear(const libMesh::Elem* elem);
};

// Specialization of build for 2D
template <>
std::unique_ptr<FEMapping<2, 2> > FEMapping<2, 2>::build(const key_type key, const FEUpdateFlags update_flags);

// Specialization of build for 3D
template <>
std::unique_ptr<FEMapping<3, 3> > FEMapping<3, 3>::build(const key_type key, const FEUpdateFlags update_flags);

} // namespace IBTK

#endif //#ifdef IBTK_HAVE_LIBMESH
#endif //#ifndef included_IBTK_FEMapping
