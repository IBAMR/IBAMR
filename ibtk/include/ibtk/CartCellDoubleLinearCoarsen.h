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
 * \brief Restrict cell-centered double data with the cell-volume adjoint of linear interpolation.
 *
 * With P the interpolation of CartCellDoubleLinearRefine (including its constant
 * extension at physical boundaries) and V_f and V_c the diagonal matrices of
 * fine and coarse cell volumes, restriction is V_c^{-1} P^T V_f.
 *
 * The declared stencil width is in fine cells and must be at least
 * floor(ratio(axis)/2) in each direction, for every level using this operator.
 * Callers must allocate and fill the source ghosts before coarsening. Source and
 * destination data must have the same depth, and each must have the same ghost
 * cell width in every direction.
 *
 * The result is the adjoint of the prolongation over the fine level only if the
 * source ghost cells hold the values of the neighboring fine cells of the same
 * level (including periodic images) where such cells exist, and zero elsewhere.
 * Fine cells that are not on the level do not contribute to the prolongation
 * matrix, so they must not contribute here. Nonzero values in ghost cells at a
 * coarse-fine boundary are used as given.
 *
 * At a coarse-fine interface the result is this adjoint over the fine level, not
 * normalized: with zero ghost cells, a constant restricts to less than that
 * constant in the coarse cells next to the interface (0.875 of it in 2D with
 * refinement ratio 2). The restriction from the prolongation matrix, the
 * transpose scaled by IBTK::PETScMatUtilities::constructRestrictionScalingOp(),
 * normalizes those rows, and it also has nonzero rows in coarse cells that the
 * fine level does not cover, which this operator never writes. The two agree on
 * levels without a coarse-fine interface, which are the ones that the tests cover.
 *
 * The name of the operator is "IBTK_LINEAR_COARSEN". Add an instance to the grid geometry with
 * addSpatialCoarsenOperator() to look it up by that name, or pass it to a schedule directly.
 */
class CartCellDoubleLinearCoarsen : public SAMRAI::xfer::CoarsenOperator<NDIM>
{
public:
    /*! \brief Construct an operator with the declared nonnegative fine-cell stencil width. */
    explicit CartCellDoubleLinearCoarsen(SAMRAI::hier::IntVector<NDIM> gcw = SAMRAI::hier::IntVector<NDIM>(1));
    ~CartCellDoubleLinearCoarsen() override = default;

    /*! \copydoc SAMRAI::xfer::CoarsenOperator::findCoarsenOperator */
    bool findCoarsenOperator(const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>>& var,
                             const std::string& op_name) const override;
    /*! \brief Return the operator's registration name. */
    const std::string& getOperatorName() const override;
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
