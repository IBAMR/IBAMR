// Filename: QuadratureCache.h
// Created on 30 Jan 2019 by David Wells
//
// Copyright (c) 2019, Boyce Griffith
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

#ifndef included_IBTK_QuadratureCache
#define included_IBTK_QuadratureCache

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_order.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/quadrature.h>

#include <map>
#include <memory>
#include <tuple>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

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
    using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

    /**
     * Type of values stored by this class that are accessible through
     * <code>operator[]</code>.
     */
    using value_type = libMesh::QBase;

    /**
     * Constructor. Sets up a cache of Quadrature objects.
     *
     * @param dim The dimension of the Quadrature object.
     */
    QuadratureCache(const unsigned int dim);

    /**
     * Return a reference to a Quadrature object that matches the specified
     * quadrature rule type and order.
     *
     * @param quad_key a tuple of enums that completely describes
     * a libMesh quadrature rule.
     */
    value_type &
    operator[](const key_type &quad_key);

protected:
    /**
     * Dimension of the FE mesh.
     */
    const unsigned int dim;

    /**
     * Managed libMesh::Quadrature objects.
     */
    std::map<key_type, std::unique_ptr<libMesh::QBase>> quadratures;
};

inline
QuadratureCache::QuadratureCache(const unsigned int dim)
    : dim(dim)
{}

inline
QuadratureCache::value_type &
QuadratureCache::operator[](const QuadratureCache::key_type &quad_key)
{
    auto it = quadratures.find(quad_key);
    if (it == quadratures.end())
    {
        const libMesh::ElemType elem_type = std::get<0>(quad_key);
        const libMesh::QuadratureType quad_type = std::get<1>(quad_key);
        const libMesh::Order order = std::get<2>(quad_key);

        libMesh::QBase &new_quad = *(
            *quadratures.emplace(
                quad_key, libMesh::QBase::build(quad_type, dim, order)).first).second;
        new_quad.init(elem_type);
        return new_quad;
    }
    else
    {
        return *(it->second);
    }
}

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_QuadratureCache
