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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_FEMappingCache
#define included_IBTK_FEMappingCache

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#ifdef IBTK_HAVE_LIBMESH

#include <ibtk/FECache.h>
#include <ibtk/FEMapping.h>
#include <ibtk/libmesh_utilities.h>

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
 * \brief Class storing multiple IBTK::FEMapping objects. We assume that
 * quadrature rules are uniquely determined by the element type, quadrature
 * type, and approximation order. There are several places in IBTK where we
 * make this assumption, e.g., we will use data from two quadrature rules
 * assumed to be equal (by this metric) to initialize different libMesh
 * objects.
 *
 * This class essentially provides a wrapper around std::map to manage
 * IBTK::FEMapping (and classes inheriting from it) objects.
 */
template <int dim, int spacedim = dim>
class FEMappingCache
{
public:
    /**
     * Constructor.
     */
    FEMappingCache(const FEUpdateFlags update_flags);

    /**
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = quadrature_key_type;

    /**
     * Type of values stored by this class that are accessible through
     * <code>operator[]</code>.
     */
    using value_type = FEMapping<dim, spacedim>;

    /**
     * Return a reference to a jacobian calculator object that matches the specified
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
        d_mappings.clear();
    }

protected:
    /**
     * Update Flags to use when creating mappings.
     */
    const FEUpdateFlags d_update_flags;

    /**
     * Managed IBTK::FEMapping objects.
     */
    std::map<key_type, std::unique_ptr<FEMapping<dim, spacedim> > > d_mappings;
};

template <int dim, int spacedim>
FEMappingCache<dim, spacedim>::FEMappingCache(const FEUpdateFlags update_flags) : d_update_flags(update_flags)
{
}

template <int dim, int spacedim>
inline typename FEMappingCache<dim, spacedim>::value_type&
FEMappingCache<dim, spacedim>::operator[](const FEMappingCache<dim, spacedim>::key_type& quad_key)
{
    auto it = d_mappings.find(quad_key);
    if (it == d_mappings.end())
    {
        const libMesh::ElemType elem_type = std::get<0>(quad_key);
        const int elem_dim = get_dim(elem_type);
        TBOX_ASSERT(elem_dim == dim);

        auto mapping = FEMapping<dim, spacedim>::build(quad_key, d_update_flags);
        value_type& new_mapping = *(*d_mappings.emplace(quad_key, std::move(mapping)).first).second;
        return new_mapping;
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

#endif //#ifdef IBTK_HAVE_LIBMESH
#endif //#ifndef included_IBTK_FEMappingCache
