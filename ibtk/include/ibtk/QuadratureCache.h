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

#ifndef included_IBTK_QuadratureCache
#define included_IBTK_QuadratureCache

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#ifdef IBTK_HAVE_LIBMESH

#include <ibtk/libmesh_utilities.h>

#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_order.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/quadrature.h>

#include <map>
#include <memory>
#include <tuple>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/**
 * \brief Class storing multiple libMesh quadrature objects. We assume that
 * quadrature rules are uniquely determined by the element type, quadrature
 * type, and approximation order. There are several places in IBTK where we
 * make this assumption, e.g., we will use data from two quadrature rules
 * assumed to be equal (by this metric) to initialize FEMap objects.
 *
 * This class essentially provides a wrapper around std::map to manage
 * libMesh::QBase (and classes inheriting from it) objects.
 */
class QuadratureCache
{
public:
    /**
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = quadrature_key_type;

    /**
     * Type of values stored by this class that are accessible through
     * <code>operator[]</code>.
     */
    using value_type = libMesh::QBase;

    /**
     * Constructor. Sets up a cache of Quadrature objects.
     *
     * @param dim The topological dimension of the relevant libMesh::Mesh: see
     * libMesh::MeshBase::mesh_dimension() for more information.
     */
    QuadratureCache(unsigned int dim);

    /**
     * Return a reference to a Quadrature object that matches the specified
     * quadrature rule type and order.
     *
     * @param quad_key a tuple of enums that completely describes
     * a libMesh quadrature rule.
     */
    value_type& operator[](const key_type& quad_key);

    /**
     * Clear the cache.
     */
    void clear()
    {
        d_quadratures.clear();
    }

protected:
    /**
     * Topological dimension of the FE mesh.
     */
    unsigned int d_dim;

    /**
     * Managed libMesh::Quadrature objects.
     */
    std::map<key_type, std::unique_ptr<libMesh::QBase> > d_quadratures;
};

inline QuadratureCache::QuadratureCache(const unsigned int dim) : d_dim(dim)
{
}

inline QuadratureCache::value_type&
QuadratureCache::operator[](const QuadratureCache::key_type& quad_key)
{
    TBOX_ASSERT(static_cast<unsigned int>(get_dim(std::get<0>(quad_key))) == d_dim);
    auto it = d_quadratures.find(quad_key);
    if (it == d_quadratures.end())
    {
        const libMesh::ElemType elem_type = std::get<0>(quad_key);
        const libMesh::QuadratureType quad_type = std::get<1>(quad_key);
        const libMesh::Order order = std::get<2>(quad_key);
        const bool allow_rules_with_negative_weights = std::get<3>(quad_key);

        libMesh::QBase& new_quad =
            *(*d_quadratures.emplace(quad_key, libMesh::QBase::build(quad_type, d_dim, order)).first).second;
        new_quad.allow_rules_with_negative_weights = allow_rules_with_negative_weights;
        new_quad.init(elem_type);
        return new_quad;
    }
    else
    {
        return *(it->second);
    }
}

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifdef IBTK_HAVE_LIBMESH
#endif //#ifndef included_IBTK_QuadratureCache
