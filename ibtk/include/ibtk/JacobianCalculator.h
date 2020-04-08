// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_JacobianCalculator
#define included_IBTK_JacobianCalculator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <IBTK_config.h>

#include <ibtk/FECache.h>
#include <ibtk/ibtk_macros.h>
#include <ibtk/ibtk_utilities.h>

#include "tbox/Utilities.h"

#include <libmesh/dense_matrix.h>
#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_order.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/fe.h>
#include <libmesh/point.h>

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <array>
#include <iosfwd>
#include <tuple>
#include <vector>

namespace IBTK
{
class JacobianCalculator
{
public:
    /**
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

    /**
     * Constructor.
     */
    JacobianCalculator(const key_type quad_key);

    /**
     * Calculate the JxW values on the given element and return a reference to
     * the result.
     */
    virtual const std::vector<double>& get_JxW(const libMesh::Elem* elem) = 0;

    virtual ~JacobianCalculator() = default;

protected:
    const key_type d_quad_key;

    std::vector<libMesh::Point> d_quad_points;
    std::vector<double> d_quad_weights;
};

/*!
 * Class which can calculate various quantities related to the mapping from
 * the reference element to an element in a mesh.
 */
template <int dim, int spacedim = dim>
class Mapping : public JacobianCalculator
{
protected:
    /*!
     * Actual data computed on an element.
     */
    struct MappingData
    {
        EigenAlignedVector<Eigen::Matrix<double, spacedim, dim> > d_contravariants;

        std::vector<double> d_J;

        std::vector<double> d_JxW;

        std::vector<libMesh::Point> d_mapped_q_points;
    };

public:
    /*!
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

    /*!
     * Constructor.
     */
    Mapping(const key_type quad_key, const FEUpdateFlags update_flags);

    virtual void reinit(const libMesh::Elem* elem) = 0;

    /*!
     * Calculate the JxW values on the given element and return a reference to
     * the result.
     */
    const std::vector<double>& get_JxW(const libMesh::Elem* elem) override
    {
        reinit(elem);
        return d_JxW;
    }

protected:
    /*!
     *
     */
    FEUpdateFlags d_update_flags;

    /*!
     * Actual data computed on an element.
     */
    MappingData d_values;

    /*!
     * Convenience reference.
     */
    std::vector<double>& d_JxW = d_values.d_JxW;

    /*!
     * Convenience reference.
     */
    EigenAlignedVector<Eigen::Matrix<double, spacedim, dim> >& d_contravariants = d_values.d_contravariants;
};

/*!
 * A generic implementation for Lagrange-type elements: works for all elements
 * in that family but is less efficient than the specialized classes for
 * lower-order or tensor-product elements. Supports nonzero codimension.
 */
template <int dim, int spacedim = dim>
class LagrangeMapping : public Mapping<dim, spacedim>
{
public:
    /**
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

    /**
     * Constructor.
     */
    LagrangeMapping(const key_type quad_key, const FEUpdateFlags update_flags);

    virtual void reinit(const libMesh::Elem* elem) override;

protected:
    /**
     * Number of nodes for the considered element type.
     */
    const std::size_t d_n_nodes;

    /**
     * Values of shape function gradients on the reference element at
     * quadrature points.
     */
    boost::multi_array<std::array<double, dim>, 2> d_dphi;
};

/*
 * Specialization for TRI3 elements with codimension zero.
 */
class Tri3Mapping : public Mapping<2, 2>
{
public:
    /**
     * Explicitly use the base class' constructor (this class does not require
     * any additional setup).
     */
    using Mapping<2, 2>::Mapping;

    virtual void reinit(const libMesh::Elem* elem) override;
};

/*
 * Specialization for QUAD4 elements with codimension zero.
 */
class Quad4Mapping : public Mapping<2, 2>
{
public:
    /**
     * Explicitly use the base class' constructor (this class does not require
     * any additional setup).
     */
    using Mapping<2, 2>::Mapping;

    virtual void reinit(const libMesh::Elem* elem) override;
};

/*
 * Specialization for QUAD9 elements with codimension zero.
 */
class Quad9Mapping : public Mapping<2, 2>
{
public:
    /**
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

    /**
     * Constructor.
     */
    Quad9Mapping(const key_type quad_key, const FEUpdateFlags update_flags);

    virtual void reinit(const libMesh::Elem* elem) override;

protected:
    /**
     * Number of 1D quadrature points in the rule used to generate the 2D
     * tensor product rule.
     */
    std::size_t d_n_oned_q_points;

    /**
     * Table containing the values of 1D shape functions (which, with a tensor
     * product, define the mapping) at reference quadrature points.
     */
    libMesh::DenseMatrix<double> d_phi;

    /**
     * Table containing the derivatives of 1D shape functions (which, with a
     * tensor product, define the mapping) at reference quadrature points.
     */
    libMesh::DenseMatrix<double> d_dphi;
};

/*
 * Specialization for TET4 elements.
 */
class Tet4Mapping : public Mapping<3, 3>
{
public:
    /**
     * Explicitly use the base class' constructor (this class does not require
     * any additional setup).
     */
    using Mapping<3, 3>::Mapping;

    virtual void reinit(const libMesh::Elem* elem) override;
};
} // namespace IBTK

#endif //#ifndef included_IBTK_JacobianCalculator
