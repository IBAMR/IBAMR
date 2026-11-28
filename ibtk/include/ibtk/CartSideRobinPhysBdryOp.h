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

#ifndef included_IBTK_CartSideRobinPhysBdryOp
#define included_IBTK_CartSideRobinPhysBdryOp

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/RobinPhysBdryPatchStrategy.h>

#include <tbox/Array.h>

#include <ComponentSelector.h>
#include <IntVector.h>

#include <set>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BoundaryBox;
template <int DIM>
class Patch;
} // namespace hier
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CartSideRobinPhysBdryOp is a concrete
 * SAMRAI::xfer::RefinePatchStrategy for setting Robin boundary conditions at
 * physical boundaries for side-centered vector-valued quantities.
 *
 * \note This class is intended to be used to specify boundary conditions for
 * MAC vector fields and may not work correctly for other types of data.
 *
 * \note Patch data with depth greater than one is treated as a collection of
 * MAC vector fields.  Such data requires NDIM boundary condition objects per
 * field, with the object for component \p axis of field \p k stored at index
 * NDIM * k + axis.
 *
 * \warning Presently, this class only supports pure Dirichlet or pure Neumann
 * boundary conditions for the normal component of the vector field.  Mixed
 * (Robin) boundary conditions are \em not supported in the normal direction.
 *
 * A fill to ghost width w sets the first w ghost layers normal to each boundary. The "LINEAR" interpolation uses the w
 * cells next to the boundary. The "QUADRATIC" interpolation uses w + 1 cells for the tangential components, and for the
 * normal component two cells where the condition is Dirichlet and w cells elsewhere. These cells must lie in the patch
 * or in ghost cells set by an earlier fill. The fill at the opposite boundary does not use the faces on a Dirichlet
 * boundary. The stencil width for "QUADRATIC" interpolation is the ghost width of the data. At corners and edges, the
 * "QUADRATIC" interpolation is exact for quadratic functions and uses only the ghost values at the adjacent boundaries
 * and the two cells nearest to the corner along each boundary normal.
 */
class CartSideRobinPhysBdryOp : public RobinPhysBdryPatchStrategy
{
public:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor leaves the class in an undefined state.  In
     * particular, the patch data indices requiring cell filling and the
     * boundary condition specifications must be set prior to use of method
     * setPhysicalBoundaryConditions().
     *
     * \see setPatchDataIndex
     * \see setPatchDataIndices
     * \see setPhysicalBcCoef
     * \see setPhysicalBcCoefs
     * \see setHomogeneousBc
     */
    CartSideRobinPhysBdryOp();

    /*!
     * \brief Constructor to fill boundary conditions for vector-valued
     * quantities.  Distinct boundary condition objects are provided for each
     * vector component.
     *
     * \param patch_data_index  Patch data index requiring ghost cell filling.
     * \param bc_coefs          Robin boundary conditions to use with this class.
     * \param homogeneous_bc    Whether to employ the homogeneous form of the boundary
     * conditions.
     * \param type              Type of interpolation to use. Valid options are "LINEAR"
     * or "QUADRATIC". Defaults to "LINEAR".
     */
    CartSideRobinPhysBdryOp(int patch_data_index,
                            const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                            bool homogeneous_bc = false,
                            std::string type = "LINEAR");

    /*!
     * \brief Constructor to fill boundary conditions for vector-valued
     * quantities.  Distinct boundary condition objects are provided for each
     * vector component.
     *
     * \param patch_data_indices  Collection of patch data indices which require ghost cell
     * filling.
     * \param bc_coefs            Robin boundary conditions to use with this class.
     * \param homogeneous_bc      Whether to employ the homogeneous form of the boundary
     * conditions.
     * \param type              Type of interpolation to use. Valid options are "LINEAR"
     * or "QUADRATIC". Defaults to "LINEAR".
     */
    CartSideRobinPhysBdryOp(const std::set<int>& patch_data_indices,
                            const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                            bool homogeneous_bc = false,
                            std::string type = "LINEAR");

    /*!
     * \brief Constructor to fill boundary conditions for vector-valued
     * quantities.  Distinct boundary condition objects are provided for each
     * vector component.
     *
     * \param patch_data_indices  Collection of patch data indices which require ghost cell
     * filling.
     * \param bc_coefs            Robin boundary conditions to use with this class.
     * \param homogeneous_bc      Whether to employ the homogeneous form of the boundary
     * conditions.
     * \param type              Type of interpolation to use. Valid options are "LINEAR"
     * or "QUADRATIC". Defaults to "LINEAR".
     */
    CartSideRobinPhysBdryOp(const SAMRAI::hier::ComponentSelector& patch_data_indices,
                            const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                            bool homogeneous_bc = false,
                            std::string type = "LINEAR");

    /*!
     * \brief Destructor.
     */
    ~CartSideRobinPhysBdryOp();

    /*!
     * \name Implementation of SAMRAI::xfer::RefinePatchStrategy interface.
     */
    //\{

    /*!
     * Function to set data associated with the given list of patch data indices
     * at patch boundaries that intersect the physical domain boundary.  The
     * patch data components set in this routine correspond to the "scratch"
     * components specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * \param patch                Patch on which to fill boundary data.
     * \param fill_time            Double simulation time for boundary filling.
     * \param ghost_width_to_fill  Integer vector describing maximum ghost width to fill over
     * all registered scratch components.
     */
    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch<NDIM>& patch,
                                       double fill_time,
                                       const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill) override;

    /*!
     * Function to return maximum stencil width needed over user-defined data
     * interpolation operations.  This is needed to determine the correct
     * interpolation data dependencies.
     *
     * The width is the largest ghost cell width of the data for "QUADRATIC" interpolation and one otherwise.
     */
    SAMRAI::hier::IntVector<NDIM> getRefineOpStencilWidth() const override;

    //\}

    /*!
     * Function to accumulate data near physical boundaries from values set in
     * the ghost cell region using the adjoint of the operator used to
     * extrapolate the ghost cell values.  This function can be used to
     * construct the adjoint of linear operators that use ghost cell data.
     *
     * \param patch                Patch on which to fill boundary data.
     * \param fill_time            Double simulation time for boundary filling.
     * \param ghost_width_to_fill  Integer vector describing maximum ghost width to fill over
     * all registered scratch components.
     *
     * \note Emits an error if the boundary data of a boundary condition object depend on the patch data (see
     * ExtendedRobinBcCoefStrategy::dependsOnTargetPatchData()), because the result would not be the adjoint.
     */
    void accumulateFromPhysicalBoundaryData(SAMRAI::hier::Patch<NDIM>& patch,
                                            double fill_time,
                                            const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill) override;

    /*!
     * \brief Set the boundary faces of \a patch at which the normal component of the registered data has a Dirichlet
     * condition (b = 0 in a*u + b*du/dn = g) to g/a, or to zero if the boundary conditions are homogeneous.
     */
    void setDataAtDirichletBoundaries(SAMRAI::hier::Patch<NDIM>& patch, double fill_time);

protected:
private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CartSideRobinPhysBdryOp(const CartSideRobinPhysBdryOp& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CartSideRobinPhysBdryOp& operator=(const CartSideRobinPhysBdryOp& that) = delete;

    /*!
     * \brief Set the boundary conditions for normal components along the
     * co-dimension one boundary.
     */
    void
    fillGhostCellValuesCodim1Normal(int patch_data_idx,
                                    const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM>>& physical_codim1_boxes,
                                    double fill_time,
                                    const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill,
                                    SAMRAI::hier::Patch<NDIM>& patch,
                                    bool adjoint_op);

    /*!
     * \brief Set the boundary conditions for transverse components along the
     * co-dimension one boundary.
     */
    void fillGhostCellValuesCodim1Transverse(
        int patch_data_idx,
        const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM>>& physical_codim1_boxes,
        double fill_time,
        const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill,
        SAMRAI::hier::Patch<NDIM>& patch,
        bool adjoint_op);

    /*!
     * \brief Set the boundary conditions along the co-dimension two boundary.
     */
    void fillGhostCellValuesCodim2(int patch_data_idx,
                                   const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM>>& physical_codim2_boxes,
                                   const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill,
                                   const SAMRAI::hier::Patch<NDIM>& patch,
                                   bool adjoint_op);

#if (NDIM > 2)
    /*!
     * \brief Set the boundary conditions along the co-dimension three boundary.
     */
    void fillGhostCellValuesCodim3(int patch_data_idx,
                                   const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM>>& physical_codim3_boxes,
                                   const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill,
                                   const SAMRAI::hier::Patch<NDIM>& patch,
                                   bool adjoint_op);
#endif

    std::string d_type = "LINEAR";
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_CartSideRobinPhysBdryOp
