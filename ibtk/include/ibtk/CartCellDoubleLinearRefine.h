// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_CartCellDoubleLinearRefine
#define included_IBTK_CartCellDoubleLinearRefine

#include <ibtk/config.h>

#include <tbox/Pointer.h>

#include <Box.h>
#include <IntVector.h>
#include <RefineOperator.h>

#include <string>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

namespace IBTK
{
/*!
 * \brief Tensor-product linear interpolation of cell-centered double data.
 *
 * Each fine cell value is interpolated linearly, in each coordinate, between the
 * values at the two nearest coarse cell centers. Between the center of a coarse
 * cell that touches a physical boundary and that boundary, the coarse cell's
 * value is used unchanged (constant extension), so no coarse cells outside the
 * physical domain are read. Source and destination data must have the same depth,
 * and each must have the same ghost cell width in every direction.
 *
 * The name of the operator is "IBTK_LINEAR_REFINE". Add an instance to the grid geometry with
 * addSpatialRefineOperator() to look it up by that name, or pass it to a schedule directly.
 */
class CartCellDoubleLinearRefine : public SAMRAI::xfer::RefineOperator<NDIM>
{
public:
    CartCellDoubleLinearRefine() = default;
    ~CartCellDoubleLinearRefine() override = default;

    /*! \copydoc SAMRAI::xfer::RefineOperator::findRefineOperator */
    bool findRefineOperator(const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>>& var,
                            const std::string& op_name) const override;
    /*! \brief Return the operator's registration name. */
    const std::string& getOperatorName() const override;
    int getOperatorPriority() const override;
    /*! \brief Return the one-cell coarse stencil width. */
    SAMRAI::hier::IntVector<NDIM> getStencilWidth() const override;

    /*! \brief Refine source coarse data into the requested fine box. */
    void refine(SAMRAI::hier::Patch<NDIM>& fine,
                const SAMRAI::hier::Patch<NDIM>& coarse,
                int dst_component,
                int src_component,
                const SAMRAI::hier::Box<NDIM>& fine_box,
                const SAMRAI::hier::IntVector<NDIM>& ratio) const override;

private:
    static const std::string s_op_name;
};
} // namespace IBTK

#endif
