// Filename: JacobianCalculatorCache.h
// Created on 19 Nov 2019 by David Wells & Jordan Brown
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
