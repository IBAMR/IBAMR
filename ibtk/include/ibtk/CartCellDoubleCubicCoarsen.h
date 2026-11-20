// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2026 by the IBAMR developers
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

#ifndef included_IBTK_CartCellDoubleCubicCoarsen
#define included_IBTK_CartCellDoubleCubicCoarsen

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <tbox/Pointer.h>

#include <Box.h>
#include <CartesianCellDoubleWeightedAverage.h>
#include <CoarsenOperator.h>
#include <IntVector.h>
#include <PatchLevel.h>

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

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CartCellDoubleCubicCoarsen is a concrete
 * SAMRAI::xfer::CoarsenOperator for restricting cell-centered double precision
 * patch data via cubic interpolation.
 *
 * \note This coarsen operator requires that the refinement ratio between coarse
 * and fine index spaces is even.
 *
 * \note With use_source_ghost_layer = true and a refinement ratio of 2 along an axis, the cubic stencil uses one layer
 * of fine values from the ghost cells of the source data. These ghost cells must contain the values of the fine level
 * where its boxes cover them, and NaN elsewhere. Near the boundary of the level, the operator uses one-sided stencils
 * or the two-point average. If use_source_ghost_layer is false or a ratio is odd, and any ratio is less than 4, the
 * operator emits a warning and uses simple weighted averaging on all axes.
 */
class CartCellDoubleCubicCoarsen : public SAMRAI::xfer::CoarsenOperator<NDIM>
{
public:
    /*!
     * \brief Constructor.
     * \param use_source_ghost_layer Whether to use the first ghost layer of the source data, which must contain the
     * values of the fine level where the level covers it, and NaN elsewhere.
     */
    explicit CartCellDoubleCubicCoarsen(bool use_source_ghost_layer = false);

    /*!
     * \brief Destructor.
     */
    ~CartCellDoubleCubicCoarsen() = default;

    /*!
     * \brief Set the fine level, whose boxes determine which fine cells are covered, or release it with nullptr.
     * Required for use_source_ghost_layer = true at a refinement ratio less than 4.
     */
    void setFineLevel(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> fine_level);

    /*!
     * \name Implementation of SAMRAI::xfer::CoarsenOperator interface.
     */
    //\{

    /*!
     * Return true if the coarsening operation matches the variable and name
     * string identifier request; false, otherwise.
     */
    bool findCoarsenOperator(const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>>& var,
                             const std::string& op_name) const override;

    /*!
     * Return name string identifier of the coarsening operation.
     */
    const std::string& getOperatorName() const override;

    /*!
     * Return the priority of this operator relative to other coarsening
     * operators.  The SAMRAI transfer routines guarantee that coarsening using
     * operators with lower priority will be performed before those with higher
     * priority.
     */
    int getOperatorPriority() const override;

    /*!
     * Return the stencil width associated with the coarsening operator.  The
     * SAMRAI transfer routines guarantee that the source patch will contain
     * sufficient ghost cell data surrounding the interior to satisfy the
     * stencil width requirements for each coarsening operator.  The width is one
     * for an operator constructed with use_source_ghost_layer = true, which uses one
     * fine ghost layer at a refinement ratio of 2, and zero otherwise.
     */
    SAMRAI::hier::IntVector<NDIM> getStencilWidth() const override;

    /*!
     * Coarsen the source component on the fine patch to the destination
     * component on the coarse patch. The coarsening operation is performed on
     * the intersection of the destination patch and the coarse box.  The fine
     * patch is guaranteed to contain sufficient data for the stencil width of
     * the coarsening operator.
     */
    void coarsen(SAMRAI::hier::Patch<NDIM>& coarse,
                 const SAMRAI::hier::Patch<NDIM>& fine,
                 int dst_component,
                 int src_component,
                 const SAMRAI::hier::Box<NDIM>& coarse_box,
                 const SAMRAI::hier::IntVector<NDIM>& ratio) const override;

    //\}

protected:
private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CartCellDoubleCubicCoarsen(const CartCellDoubleCubicCoarsen& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CartCellDoubleCubicCoarsen& operator=(const CartCellDoubleCubicCoarsen& that) = delete;

    /*!
     * The operator name.
     */
    static const std::string s_op_name;

    /*!
     * "Backup" coarsen operator for refinement ratios less than 4 without the source ghost layer.
     */
    SAMRAI::geom::CartesianCellDoubleWeightedAverage<NDIM> d_weighted_average_coarsen_op;

    /*!
     * Whether to use the first ghost layer of the source data, which contains the values of the fine level where the
     * level covers it.
     */
    bool d_use_source_ghost_layer;

    /*!
     * The fine level, whose boxes determine which fine cells are covered.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> d_fine_level;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_CartCellDoubleCubicCoarsen
