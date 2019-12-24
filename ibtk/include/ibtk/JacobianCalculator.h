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

#include <ibtk/ibtk_macros.h>

#include "tbox/Utilities.h"

#include <libmesh/dense_matrix.h>
#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_order.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/fe.h>
#include <libmesh/point.h>

IBTK_DISABLE_EXTRA_WARNINGS
#include "boost/multi_array.hpp"
IBTK_ENABLE_EXTRA_WARNINGS

#include <array>
#include <tuple>
#include <vector>

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

    std::vector<double> d_JxW;
};

/*
 * A generic implementation for Lagrange-type elements: works for all elements
 * in that family but is less efficient than the specialized classes for
 * lower-order or tensor-product elements. Only supports codimension zero.
 */
template <int dim>
class LagrangeJacobianCalculator : public JacobianCalculator
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
    LagrangeJacobianCalculator(const key_type quad_key);

    virtual const std::vector<double>& get_JxW(const libMesh::Elem* elem) override;

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
 * Specialization for Edge2 elements.
 */
class Edge2JacobianCalculator : public JacobianCalculator
{
public:
    /**
     * Explicitly use the base class' constructor (this class does not require
     * any additional setup).
     */
    using JacobianCalculator::JacobianCalculator;

    virtual const std::vector<double>& get_JxW(const libMesh::Elem* elem) override;
};

/*
 * Specialization for TRI3 elements.
 */
class Tri3JacobianCalculator : public JacobianCalculator
{
public:
    /**
     * Explicitly use the base class' constructor (this class does not require
     * any additional setup).
     */
    using JacobianCalculator::JacobianCalculator;

    virtual const std::vector<double>& get_JxW(const libMesh::Elem* elem) override;
};

/*
 * Specialization for QUAD4 elements.
 */
class Quad4JacobianCalculator : public JacobianCalculator
{
public:
    /**
     * Explicitly use the base class' constructor (this class does not require
     * any additional setup).
     */
    using JacobianCalculator::JacobianCalculator;

    virtual const std::vector<double>& get_JxW(const libMesh::Elem* elem) override;
};

/*
 * Specialization for QUAD9 elements.
 */
class Quad9JacobianCalculator : public JacobianCalculator
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
    Quad9JacobianCalculator(const key_type quad_key);

    virtual const std::vector<double>& get_JxW(const libMesh::Elem* elem) override;

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
class Tet4JacobianCalculator : public JacobianCalculator
{
public:
    /**
     * Explicitly use the base class' constructor (this class does not require
     * any additional setup).
     */
    using JacobianCalculator::JacobianCalculator;

    virtual const std::vector<double>& get_JxW(const libMesh::Elem* elem) override;
};

#endif //#ifndef included_IBTK_JacobianCalculator
