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
 * IBTK:JacobianCalculator (and classes inheriting from it) objects.
 */
class JacobianCalculatorCache
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
     * Managed libMesh::Quadrature objects.
     */
    std::map<key_type, std::unique_ptr<JacobianCalculator> > d_jacobian_calculators;
};

inline JacobianCalculatorCache::value_type& JacobianCalculatorCache::
operator[](const JacobianCalculatorCache::key_type& quad_key)
{
    auto it = d_jacobian_calculators.find(quad_key);
    if (it == d_jacobian_calculators.end())
    {
        const libMesh::ElemType elem_type = std::get<0>(quad_key);

        std::unique_ptr<JacobianCalculator> jac_calc;
        switch (elem_type)
        {
		case libMesh::EDGE2:
            jac_calc.reset(new Edge2JacobianCalculator(quad_key));
            break;
        case libMesh::TRI3:
            jac_calc.reset(new Tri3JacobianCalculator(quad_key));
            break;
        case libMesh::TRI6:
            jac_calc.reset(new LagrangeJacobianCalculator<2>(quad_key));
            break;
        case libMesh::QUAD4:
            jac_calc.reset(new Quad4JacobianCalculator(quad_key));
            break;
        case libMesh::QUAD9:
            jac_calc.reset(new Quad9JacobianCalculator(quad_key));
            break;
        case libMesh::TET4:
            jac_calc.reset(new Tet4JacobianCalculator(quad_key));
            break;
        case libMesh::TET10:
        case libMesh::HEX8:
        case libMesh::HEX27:
            jac_calc.reset(new LagrangeJacobianCalculator<3>(quad_key));
            break;
        default:
            TBOX_ERROR("unimplemented element type");
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
