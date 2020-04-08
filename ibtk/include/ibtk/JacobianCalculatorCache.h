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

#ifndef included_IBTK_JacobianCalculatorCache
#define included_IBTK_JacobianCalculatorCache

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/JacobianCalculator.h>
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
 * \brief Class storing multiple JacobianCalculator objects. We assume that
 * quadrature rules are uniquely determined by the element type, quadrature
 * type, and approximation order. There are several places in IBTK where we
 * make this assumption, e.g., we will use data from two quadrature rules
 * assumed to be equal (by this metric) to initialize FEMap objects.
 *
 * This class essentially provides a wrapper around std::map to manage
 * IBTK::JacobianCalculator (and classes inheriting from it) objects.
 */
class JacobianCalculatorCache
{
public:
    /**
     * Constructor. Takes, as argument, the spatial dimension of the
     * mesh. Here the spatial dimension refers to the number of relevant
     * coordinates in each node: e.g., for a 2D surface mesh (comprised of
     * QUAD4 elements) in a 3D space the spatial dimension is 3, while for a
     * 2D mesh (also comprised of QUAD4 elements) in a 2D space the spatial
     * dimension is 2.
     *
     * @seealso libMesh::MeshBase::spatial_dimension() defines the spatial
     * dimension in the same way.
     */
    JacobianCalculatorCache(const int spatial_dimension);

    /**
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

    /**
     * Type of values stored by this class that are accessible through
     * <code>operator[]</code>.
     */
    using value_type = JacobianCalculator;

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
        d_jacobian_calculators.clear();
    }

protected:
    /**
     * Spatial dimension of the mesh.
     */
    const int d_spatial_dimension;

    /**
     * Managed libMesh::Quadrature objects.
     */
    std::map<key_type, std::unique_ptr<JacobianCalculator> > d_jacobian_calculators;
};

inline JacobianCalculatorCache::JacobianCalculatorCache(const int spatial_dimension)
    : d_spatial_dimension(spatial_dimension)
{
    TBOX_ASSERT(0 < spatial_dimension && spatial_dimension <= 3);
}

inline JacobianCalculatorCache::value_type&
    JacobianCalculatorCache::operator[](const JacobianCalculatorCache::key_type& quad_key)
{
    auto it = d_jacobian_calculators.find(quad_key);
    if (it == d_jacobian_calculators.end())
    {
        const libMesh::ElemType elem_type = std::get<0>(quad_key);
        const int dim = get_dim(elem_type);

        std::unique_ptr<JacobianCalculator> jac_calc;
        switch (d_spatial_dimension)
        {
        case 1:
            jac_calc.reset(new LagrangeMapping<1>(quad_key, FEUpdateFlags::update_JxW));
            break;
        case 2:
            switch (elem_type)
            {
            case libMesh::EDGE2:
            case libMesh::EDGE3:
            case libMesh::EDGE4:
                jac_calc.reset(new LagrangeMapping<1, 2>(quad_key, FEUpdateFlags::update_JxW));
                break;
            case libMesh::TRI3:
                jac_calc.reset(new Tri3Mapping(quad_key, FEUpdateFlags::update_JxW));
                break;
            case libMesh::QUAD4:
                jac_calc.reset(new Quad4Mapping(quad_key, FEUpdateFlags::update_JxW));
                break;
            case libMesh::TRI6:
            case libMesh::QUAD8:
                jac_calc.reset(new LagrangeMapping<2, 2>(quad_key, FEUpdateFlags::update_JxW));
                break;
            case libMesh::QUAD9:
                jac_calc.reset(new Quad9Mapping(quad_key, FEUpdateFlags::update_JxW));
                break;
            default:
                TBOX_ERROR("unimplemented element type");
            }
            break;
        case 3:
            if (dim == 1)
                jac_calc.reset(new LagrangeMapping<1, 3>(quad_key, FEUpdateFlags::update_JxW));
            else if (dim == 2)
                jac_calc.reset(new LagrangeMapping<2, 3>(quad_key, FEUpdateFlags::update_JxW));
            else if (elem_type == libMesh::TET4)
                jac_calc.reset(new Tet4Mapping(quad_key, FEUpdateFlags::update_JxW));
            else
                jac_calc.reset(new LagrangeMapping<3, 3>(quad_key, FEUpdateFlags::update_JxW));
            break;
        default:
            TBOX_ERROR("unimplemented spatial dimension");
        }

        JacobianCalculator& new_jacob = *(*d_jacobian_calculators.emplace(quad_key, std::move(jac_calc)).first).second;
        return new_jacob;
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

#endif //#ifndef included_IBTK_JacobianCalculatorCache
