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

#ifndef included_IBTK_FECache
#define included_IBTK_FECache

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/QuadratureCache.h>

#include <tbox/Utilities.h>

#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_order.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/fe.h>

#include <map>
#include <memory>
#include <tuple>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

/**
 * Enumeration describing the various update options available for
 * libMesh::FEBase objects stored by FECache. Multiple flags can be enabled
 * with bitwise or operations, e.g.,
 *
 * @code
 * const IBTK::FEUpdateFlags update_flags = IBTK::update_phi | IBTK::update_dphi;
 * @endcode
 *
 * See IBTK::FECache for more information.
 */
enum FEUpdateFlags
{
    /**
     * Do not update anything.
     */
    update_default = 0,

    /**
     * Update phi (shape function values).
     */
    update_phi = 1,

    /**
     * Update dphi (shape function gradients).
     */
    update_dphi = 2
};

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
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

    /**
     * Type of values stored by this class that are accessible through
     * <code>operator[]</code>.
     */
    using value_type = libMesh::FEBase;

    /**
     * Constructor. Sets up a cache of FE objects calculating values for the
     * given FEType argument. All cached FE objects have the same FEType.
     *
     * @param dim The dimension of the FE object.
     *
     * @param fe_type The libMesh FEType object describing the relevant finite
     * element.
     *
     * @param flags FEUpdateFlags indicating which values should be calculated
     * by each libMesh::FEBase object.
     */
    FECache(const unsigned int dim, const libMesh::FEType& fe_type, const FEUpdateFlags flags);

    /**
     * Return a reference to an FE object that matches the specified
     * quadrature rule type and order on the given element.
     *
     * @param quad_key a tuple of enums that completely describes
     * a libMesh quadrature rule.
     *
     * @param elem Pointer to an element. This will be used to, if necessary,
     * reinitialize the returned FE object.
     */
    value_type& operator()(const key_type& quad_key, const libMesh::Elem* elem);

    /**
     * Return the FEUpdateFlags stored by the current FECache.
     */
    FEUpdateFlags getFEUpdateFlags() const;

    /**
     * Return the FEType stored by the current FECache.
     */
    libMesh::FEType getFEType() const;

protected:
    /**
     * Dimension of the FE mesh.
     */
    const unsigned int d_dim;

    /**
     * Object describing the finite element type.
     */
    const libMesh::FEType d_fe_type;

    /**
     * Update flags for the current object.
     */
    const FEUpdateFlags d_update_flags;

    /**
     * Managed libMesh::Quadrature objects. These are attached to the FE
     * objects.
     */
    QuadratureCache d_quadrature_cache;

    /**
     * Managed libMesh::FE objects of specified dimension and family.
     */
    std::map<key_type, std::unique_ptr<libMesh::FEBase> > d_fes;
};

inline FECache::FECache(const unsigned int dim, const libMesh::FEType& fe_type, const FEUpdateFlags flags)
    : d_dim(dim), d_fe_type(fe_type), d_update_flags(flags), d_quadrature_cache(d_dim)
{
}

inline FEUpdateFlags
FECache::getFEUpdateFlags() const
{
    return d_update_flags;
}

inline libMesh::FEType
FECache::getFEType() const
{
    return d_fe_type;
}

inline FECache::value_type&
FECache::operator()(const FECache::key_type& quad_key, const libMesh::Elem* elem)
{
#ifndef NDEBUG
    TBOX_ASSERT(elem->type() == std::get<0>(quad_key));
#endif
    auto it = d_fes.find(quad_key);
    if (it == d_fes.end())
    {
        libMesh::QBase& quad = d_quadrature_cache[quad_key];
        libMesh::FEBase& fe = *(*d_fes.emplace(quad_key, libMesh::FEBase::build(d_dim, d_fe_type)).first).second;
        fe.attach_quadrature_rule(&quad);

        if (d_update_flags & FEUpdateFlags::update_phi) fe.get_phi();
        if (d_update_flags & FEUpdateFlags::update_dphi) fe.get_dphi();

        fe.reinit(elem);

        return fe;
    }
    else
    {
        libMesh::FEBase& fe = *(it->second);
        // TODO: we need better reinitialization logic than hardcoding in
        // libMesh element types.
        //
        // Subdivision elements have a variable number of degrees of freedom
        // per cell (and the values at quadrature points depend on both the
        // current and neighbor element geometries) so these values must
        // always be recomputed.
        if (d_update_flags & FEUpdateFlags::update_dphi || d_fe_type.family == libMesh::FEFamily::SUBDIVISION)
            fe.reinit(elem);
        return fe;
    }
}

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_FECache
