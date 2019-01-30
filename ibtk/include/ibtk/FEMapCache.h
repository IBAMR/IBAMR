// Filename: FEMapCache.h
// Created on 29 Jan 2019 by David Wells
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

#ifndef included_IBTK_FEMapCache
#define included_IBTK_FEMapCache

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_order.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/fe_map.h>
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
 * \brief Class storing multiple libMesh::FEMap objects, each corresponding to
 * a different quadrature rule. Each FEMap object is configured with a
 * quadrature rule corresponding to the provided <code>quad_key</code>
 * parameter.
 *
 * In some cases we only need to recalculate the products of the Jacobians and
 * quadrature weights, but not the shape function values: at the present time
 * this is not possible to do in libMesh through the standard FEBase
 * interface. Hence, in IBAMR, we cache the FE (which compute shape function
 * values) and FEMap (which compute Jacobians) objects separately and only
 * call <code>reinit</code> on the appropriate object when necessary.
 *
 * This class essentially provides a wrapper around std::map to manage FEMap
 * objects and the quadrature rules they use. The keys are descriptions of
 * quadrature rules.
 *
 * @note At the present time the only values accessible through the FEMap
 * objects stored by this class are the Jacobians and JxW values: no second
 * derivative or physical quadrature point information is computed.
 */
class FEMapCache
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
    using value_type = libMesh::FEMap;

    /**
     * Constructor. Sets up a cache of FE objects calculating values for the
     * given FEType argument. All cached FE objects have the same FEType.
     *
     * @param dim The dimension of the relevant libMesh::Mesh.
     *
     * @param dim The dimension of the relevant libMesh::Mesh.
     */
    FEMapCache(const unsigned int dim);

    /**
     * Return a reference to an FEMap object that matches the specified
     * quadrature rule type and order.
     *
     * @param quad_key a tuple of enums that completely describes
     * a libMesh quadrature rule.
     */
    libMesh::FEMap &
    operator[](const key_type &quad_key);

protected:
    /**
     * Dimension of the FE mesh.
     */
    const unsigned int dim;

    /**
     * Managed libMesh::Quadrature objects. These are attached to the FE
     * objects.
     */
    std::map<key_type, std::unique_ptr<libMesh::QBase>> quadratures;

    /**
     * Managed libMesh::FEMap objects of specified dimension and family.
     */
    std::map<key_type, libMesh::FEMap> fe_maps;
};

inline
FEMapCache::FEMapCache(const unsigned int dim)
    : dim(dim)
{}

inline
libMesh::FEMap &
FEMapCache::operator[](const FEMapCache::key_type &quad_key)
{
    const libMesh::ElemType elem_type = std::get<0>(quad_key);
    const libMesh::QuadratureType quad_type = std::get<1>(quad_key);
    const libMesh::Order order = std::get<2>(quad_key);

    auto it = fe_maps.find(quad_key);
    if (it == fe_maps.end())
    {
        // we should also need a new Quadrature object unless something has
        // gone wrong
#ifndef NDBEBUG
        TBOX_ASSERT(quadratures.find(quad_key) == quadratures.end());
#endif // ifndef NDEBUG
        std::unique_ptr<libMesh::QBase> &new_quad = (
            *quadratures.emplace(
                quad_key, libMesh::QBase::build(quad_type, dim, order)).first).second;
        new_quad->init(elem_type);

        libMesh::FEMap &fe_map = fe_maps[quad_key];
        // Calling this function enables JxW calculations
        fe_map.get_JxW();

        // Doing this may not work with future (1.4.0 or up) versions of
        // libMesh. In particular; init_reference_to_physical_map is
        // undocumented (and almost surely is not intended for use by anyone
        // but libMesh developers) and *happens* to not read any geometric or
        // topological information from the Elem argument (just the default
        // order and type).
        std::unique_ptr<libMesh::Elem> exemplar_elem(libMesh::Elem::build(elem_type));

        // This is one of very few functions in libMesh that is templated on
        // the dimension (not spatial dimension) of the mesh
        switch (dim)
        {
        case 1:
            fe_map.init_reference_to_physical_map<1>(new_quad->get_points(), exemplar_elem.get());
            break;
        case 2:
            fe_map.init_reference_to_physical_map<2>(new_quad->get_points(), exemplar_elem.get());
            break;
        case 3:
            fe_map.init_reference_to_physical_map<3>(new_quad->get_points(), exemplar_elem.get());
            break;
        default:
            TBOX_ASSERT(false);
        }
        return fe_map;
    }
    else
    {
        return it->second;
    }
}



/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_FEMapCache
