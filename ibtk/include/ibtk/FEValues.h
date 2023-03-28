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

#ifndef included_IBTK_FEValues
#define included_IBTK_FEValues

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#ifdef IBTK_HAVE_LIBMESH

#include <ibtk/FECache.h>
#include <ibtk/FEMapping.h>

#include <tbox/Utilities.h>

#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/fe.h>
#include <libmesh/fe_type.h>
#include <libmesh/point.h>
#include <libmesh/quadrature.h>
#include <libmesh/type_vector.h>

#include <map>
#include <vector>

namespace IBTK
{
/**
 * Class defining the interface to FEValues in a dimension-independent way to
 * improve compatibility with libMesh.
 */
class FEValuesBase
{
public:
    virtual ~FEValuesBase() = default;

    virtual void reinit(const libMesh::Elem* elem) = 0;

    inline const std::vector<double>& getJxW() const
    {
        return d_JxW;
    }

    inline const std::vector<libMesh::Point>& getQuadraturePoints() const
    {
        return d_quadrature_points;
    }

    inline const std::vector<std::vector<double> >& getShapeValues() const
    {
        return d_shape_values;
    }

    inline const std::vector<std::vector<libMesh::VectorValue<double> > >& getShapeGradients() const
    {
        return d_shape_gradients;
    }

    static std::unique_ptr<FEValuesBase> build(const int dim,
                                               const int spacedim,
                                               libMesh::QBase* qrule,
                                               const libMesh::FEType fe_type,
                                               const FEUpdateFlags update_flags);

protected:
    std::vector<double> d_JxW;

    std::vector<libMesh::Point> d_quadrature_points;

    std::vector<std::vector<double> > d_shape_values;

    std::vector<std::vector<libMesh::VectorValue<double> > > d_shape_gradients;
};

/**
 * Class like libMesh::FE for element shape function calculations, but
 * optimized for isoparametric Lagrange finite elements.
 */
template <int dim, int spacedim = dim>
class FEValues : public FEValuesBase
{
public:
    FEValues(libMesh::QBase* qrule, const libMesh::FEType fe_type, const FEUpdateFlags update_flags);

    virtual void reinit(const libMesh::Elem* elem) override;

protected:
    libMesh::QBase* d_qrule;

    const libMesh::FEType d_fe_type;

    /**
     * Reference values, extracted from libMesh.
     */
    struct ReferenceValues
    {
        ReferenceValues(const libMesh::QBase& quadrature, const libMesh::FEType& fe_type);

        const libMesh::ElemType d_elem_type;

        const libMesh::FEType d_fe_type;

        /**
         * shape values, indexed by quadrature point and then by shape function
         * index
         */
        boost::multi_array<double, 2> d_reference_shape_values;

        /**
         * shape gradients, indexed by quadrature point and then by shape function
         * index.
         */
        boost::multi_array<libMesh::VectorValue<double>, 2> d_reference_shape_gradients;
    };

    /*
     * Mappings, indexed by element type.
     */
    std::map<libMesh::ElemType, std::unique_ptr<FEMapping<dim, spacedim> > > d_mappings;

    /*
     * Reference values, indexed by element type.
     */
    std::map<libMesh::ElemType, ReferenceValues> d_reference_values;

    /*
     * Things to actually recompute.
     */
    FEUpdateFlags d_update_flags;

    /**
     * Last element type. We can avoid reinitializing some things if this
     * matches the current element type.
     */
    libMesh::ElemType d_last_elem_type = libMesh::ElemType::INVALID_ELEM;
};
} // namespace IBTK

#endif //#ifdef IBTK_HAVE_LIBMESH
#endif //#ifndef included_IBTK_FEValues
