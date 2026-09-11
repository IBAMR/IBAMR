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

#ifndef included_IBTK_CartCellDoubleLinearCoarsen
#define included_IBTK_CartCellDoubleLinearCoarsen

#include <ibtk/config.h>

#include <tbox/Pointer.h>

#include <Box.h>
#include <CoarsenOperator.h>
#include <IntVector.h>

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
 * \brief Restrict cell-centered double data with the scaled adjoint of linear interpolation.
 *
 * Restriction is the cell-volume adjoint of CartCellDoubleLinearRefine,
 * including its constant extension at physical boundaries.
 *
 * The declared stencil width is in fine cells and must be at least
 * floor(ratio(axis)/2) in each direction, for every level using this operator.
 * Callers must allocate and fill the source ghosts before coarsening. Source and
 * destination data must have matching depths and uniform ghost widths.
 * The operator is registered as "IBTK_LINEAR_COARSEN".
 */
class CartCellDoubleLinearCoarsen : public SAMRAI::xfer::CoarsenOperator<NDIM>
{
public:
    /*! \brief Construct an operator with the declared nonnegative fine-cell stencil width. */
    explicit CartCellDoubleLinearCoarsen(SAMRAI::hier::IntVector<NDIM> gcw = SAMRAI::hier::IntVector<NDIM>(1));
    /*! \brief Destructor. */
    ~CartCellDoubleLinearCoarsen() override = default;

    /*! \copydoc SAMRAI::xfer::CoarsenOperator::findCoarsenOperator */
    bool findCoarsenOperator(const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>>& var,
                             const std::string& op_name) const override;
    /*! \brief Return the operator's registration name. */
    const std::string& getOperatorName() const override;
    /*! \brief Return priority zero. */
    int getOperatorPriority() const override;
    /*! \brief Return the declared fine-cell stencil width. */
    SAMRAI::hier::IntVector<NDIM> getStencilWidth() const override;

    /*! \brief Coarsen source fine data into the requested coarse box. */
    void coarsen(SAMRAI::hier::Patch<NDIM>& coarse,
                 const SAMRAI::hier::Patch<NDIM>& fine,
                 int dst_component,
                 int src_component,
                 const SAMRAI::hier::Box<NDIM>& coarse_box,
                 const SAMRAI::hier::IntVector<NDIM>& ratio) const override;

private:
    static const std::string s_op_name;
    const SAMRAI::hier::IntVector<NDIM> d_gcw;
};
} // namespace IBTK

#endif
