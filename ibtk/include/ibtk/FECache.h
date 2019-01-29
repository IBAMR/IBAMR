// Filename: FECache.h
// Created on 25 Jan 2019 by David Wells
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <libmesh/enum_order.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/fe.h>
#include <libmesh/quadrature.h>

#include <map>
#include <memory>
#include <utility>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

/**
 * \brief Class storing multiple libMesh::FE objects, each corresponding to a
 * different quadrature rule. Each FE object is configured with a quadrature
 * rule corresponding to the provided <code>quad_key</code> parameter.
 *
 * The Lagrangian-Eulerian interaction code uses different quadrature rules on
 * different Elems to account for the change in size, over time, of each
 * corresponding grid cell. Since libMesh::FE objects cache values that are
 * independent of the current cell (such as shape function values) but
 * *dependent* upon the quadrature rule, it is much more efficient to store
 * one FE object for each quadrature rule instead of constantly recomputing,
 * e.g., shape function values.
 *
 * This class essentially provides a wrapper around std::map to manage FE
 * objects and the quadrature rules they use. The keys are descriptions of
 * quadrature rules.
 */
class FECache
{
public:
    /**
     * Constructor. Sets up a cache of FE objects calculating values for the
     * given FEType argument. All cached FE objects have the same FEType.
     *
     * @param dim The dimension of the FE object.
     *
     * @param fe_type The libMesh FEType object describing the relevant finite
     * element.
     */
    FECache(const unsigned int dim, const libMesh::FEType &fe_type);

    /**
     * Return a reference to an FE object that matches the specified
     * quadrature rule type and order.
     *
     * @param quad_key a pair of enums that completely describes
     * a libMesh quadrature rule.
     */
    libMesh::FEBase &
    operator[](const std::pair<libMesh::QuadratureType, libMesh::Order> &quad_key);

    /**
     * Return the FEType stored by the current FECache.
     */
    libMesh::FEType
    getFEType() const;

protected:
    /**
     * Dimension of the FE mesh.
     */
    const unsigned int dim;

    /**
     * Object describing the finite element type.
     */
    const libMesh::FEType fe_type;

    /**
     * Managed libMesh::Quadrature objects. These are attached to the FE
     * objects.
     */
    std::map<std::pair<libMesh::QuadratureType, libMesh::Order>,
             std::unique_ptr<libMesh::QBase>> quadratures;

    /**
     * Managed libMesh::FE objects of specified dimension and family.
     */
    std::map<std::pair<libMesh::QuadratureType, libMesh::Order>,
             std::unique_ptr<libMesh::FEBase>> fes;
};

inline
FECache::FECache(const unsigned int dim, const libMesh::FEType &fe_type)
    : dim(dim)
    , fe_type(fe_type)
{}

inline
libMesh::FEType
FECache::getFEType() const
{
    return fe_type;
}

inline
libMesh::FEBase &
FECache::operator[](const std::pair<libMesh::QuadratureType, libMesh::Order> &quad_key)
{
    auto it = fes.find(quad_key);
    if (it == fes.end())
    {
        // we should also need a new Quadrature object unless something has
        // gone wrong
#ifndef NDBEBUG
        TBOX_ASSERT(quadratures.find(quad_key) == quadratures.end());
#endif // ifndef NDEBUG
        std::unique_ptr<libMesh::QBase> &new_quad = (*quadratures.emplace(
            quad_key, libMesh::QBase::build(quad_key.first, dim, quad_key.second)).first).second;

        libMesh::FEBase &fe = *(*fes.emplace(quad_key,
                                             libMesh::FEBase::build(dim, fe_type)).first).second;
        fe.attach_quadrature_rule(new_quad.get());
        return fe;
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
