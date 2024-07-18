// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2023 by the IBAMR developers
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

#ifndef included_IBTK_CartSideDoubleQuadraticCFInterpolation
#define included_IBTK_CartSideDoubleQuadraticCFInterpolation

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "CartesianSideDoubleConservativeLinearRefine.h"
#include "ComponentSelector.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "RefineOperator.h"
#include "SideVariable.h"
#include "tbox/Pointer.h"

#include <set>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class CoarseFineBoundary;
template <int DIM>
class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CartSideDoubleQuadraticCFInterpolation is a concrete
 * SAMRAI::xfer::RefinePatchStrategy which sets coarse-fine interface ghost cell
 * values for side-centered double precision patch data via quadratic
 * interpolation in the normal and tangential directions at coarse-fine
 * interfaces.
 *
 *
 * For the tangential component, postprocessRefine() computes a cubic approximation to A using coarse values O. As a
 * follow up, computeNormalExtension() will compute a quadratic approximation using the value A and the interior values
 * x to compute the ghost cell G.
 *
 * ******O************************
 * *           *     *     *     *
 * *           *     *     *     *
 * *           *     *     *     *
 * *           *******************
 * *           *     *     *     *
 * *           *     *     *     *
 * *           *     *     *     *
 * ******O************************
 * *     |     *     *     *     *
 * *     |     *     *     *     *
 * *     |     *     *     *     *
 * *     A--G--***x*****x*****x***
 * *           *     *     *     *
 * *           *     *     *     *
 * *           *     *     *     *
 * ******O************************
 * *           *     *     *     *
 * *           *     *     *     *
 * *           *     *     *     *
 * *           *******************
 * *           *     *     *     *
 * *           *     *     *     *
 * *           *     *     *     *
 * ******O************************
 *
 * For the normal component, postprocessRefine() computes a quadratic approximation to A using coarse values O. As a
 * follow up, computeNormalExtension() will compute a quadratic approximation using the value A and the interior values
 * x to compute the ghost cell G.
 * *******************************
 * *           *     *     *     *
 * *           *     *     *     *
 * *           *     *     *     *
 * O           *******************
 * *           *     *     *     *
 * *           *     *     *     *
 * *           *     *     *     *
 * *******************************
 * *     |     *     *     *     *
 * A     G     x     x     x     *
 * *     |     *     *     *     *
 * O     ------******************
 * *           *     *     *     *
 * *           *     *     *     *
 * *           *     *     *     *
 * *******************************
 * *           *     *     *     *
 * *           *     *     *     *
 * *           *     *     *     *
 * O           *******************
 * *           *     *     *     *
 * *           *     *     *     *
 * *           *     *     *     *
 * *******************************
 *
 * Note: If this class is used in conjunction with RefineAlgorithm, to correctly fill in ghost cells,
 * computeNormalExtension() must be called after RefineAlgorithm::refine() is completed.
 */
class CartSideDoubleQuadraticCFInterpolation : public CoarseFineBoundaryRefinePatchStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    CartSideDoubleQuadraticCFInterpolation();

    /*!
     * \brief Destructor.
     */
    ~CartSideDoubleQuadraticCFInterpolation();

    /*!
     * \name SAMRAI::xfer::RefinePatchStrategy interface.
     */
    //\{

    /*!
     * Function to set data associated with the given list of patch data indices
     * at patch boundaries that intersect the physical domain boundary.  The
     * patch data components set in this routine correspond to the "scratch"
     * components specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * Presently, the implementation does nothing.
     *
     * \param patch                Patch on which to fill boundary data.
     * \param fill_time            Double simulation time for boundary filling.
     * \param ghost_width_to_fill  Integer vector describing maximum ghost width to fill over all registered scratch
     * components.
     */
    void setPhysicalBoundaryConditions(SAMRAI::hier::PatchNd& patch,
                                       double fill_time,
                                       const SAMRAI::hier::IntVectorNd& ghost_width_to_fill) override;

    /*!
     * Function to return maximum stencil width needed over user-defined data
     * interpolation operations.  This is needed to determine the correct
     * interpolation data dependencies.
     */
    SAMRAI::hier::IntVectorNd getRefineOpStencilWidth() const override;

    /*!
     * Function to perform user-defined preprocess data refine operations.  This
     * member function is called before standard refine operations (expressed
     * using concrete subclasses of the SAMRAI::xfer::RefineOperator base
     * class).  The preprocess function must refine data from the scratch
     * components of the coarse patch into the scratch components of the fine
     * patch on the specified fine box region.  Recall that the scratch
     * components are specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * Presently, the implementation does nothing.
     *
     * \param fine      Fine patch containing destination data.
     * \param coarse    Coarse patch containing source data.
     * \param fine_box  Box region on fine patch into which data is refined.
     * \param ratio     Integer vector containing ratio relating index space between coarse and fine patches.
     */
    void preprocessRefine(SAMRAI::hier::PatchNd& fine,
                          const SAMRAI::hier::PatchNd& coarse,
                          const SAMRAI::hier::BoxNd& fine_box,
                          const SAMRAI::hier::IntVectorNd& ratio) override;

    /*!
     * Function to perform user-defined postprocess data refine operations.
     * This member function is called after standard refine operations
     * (expressed using concrete subclasses of the SAMRAI::xfer::RefineOperator
     * base class).  The postprocess function must refine data from the scratch
     * components of the coarse patch into the scratch components of the fine
     * patch on the specified fine box region.  Recall that the scratch
     * components are specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * This function computes a quadratic approximation in the tangential direction. To complete the approximation of
     * ghost cells, computeNormalExtension() must be called after postprocessRefine().
     *
     * \param fine      Fine patch containing destination data.
     * \param coarse    Coarse patch containing source data.
     * \param fine_box  Box region on fine patch into which data is refined.
     * \param ratio     Integer vector containing ratio relating index space between coarse and fine patches.
     */
    void postprocessRefine(SAMRAI::hier::PatchNd& fine,
                           const SAMRAI::hier::PatchNd& coarse,
                           const SAMRAI::hier::BoxNd& fine_box,
                           const SAMRAI::hier::IntVectorNd& ratio) override;

    //\}

    /*!
     * \name Extension of SAMRAI::xfer::RefinePatchStrategy interface to support more
     * complex coarse-fine interface discretizations.
     */
    //\{

    /*!
     * Whether or not to employ a consistent interpolation scheme at "Type 2"
     * coarse-fine interface ghost cells.
     */
    void setConsistentInterpolationScheme(bool consistent_type_2_bdry) override;

    /*!
     * \brief Reset the patch data index operated upon by this class.
     */
    void setPatchDataIndex(int patch_data_index) override;

    /*!
     * \brief Reset the patch data indices operated upon by this class.
     */
    void setPatchDataIndices(const std::set<int>& patch_data_indices) override;

    /*!
     * \brief Reset the patch data indices operated upon by this class.
     */
    void setPatchDataIndices(const SAMRAI::hier::ComponentSelector& patch_data_indices) override;

    /*!
     * Set the patch hierarchy used in constructing coarse-fine interface
     * boundary boxes.
     */
    void setPatchHierarchy(SAMRAIPointer<SAMRAI::hier::PatchHierarchyNd> hierarchy) override;

    /*!
     * Clear the patch hierarchy used in constructing coarse-fine interface
     * boundary boxes.
     */
    void clearPatchHierarchy() override;

    /*!
     * Compute the normal extension of fine data at coarse-fine interfaces.
     *
     * This function assumes that the first ghost cell is filled with a reasonable value, see the class description.
     */
    void computeNormalExtension(SAMRAI::hier::PatchNd& patch,
                                const SAMRAI::hier::IntVectorNd& ratio,
                                const SAMRAI::hier::IntVectorNd& ghost_width_to_fill) override;

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
    CartSideDoubleQuadraticCFInterpolation(const CartSideDoubleQuadraticCFInterpolation& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CartSideDoubleQuadraticCFInterpolation& operator=(const CartSideDoubleQuadraticCFInterpolation& that) = delete;

    /*!
     * The patch data indices corresponding to the "scratch" patch data that is
     * operated on by this class.
     */
    std::set<int> d_patch_data_indices;

    /*!
     * Boolean value indicating whether we are enforcing a consistent
     * interpolation scheme at "Type 2" coarse-fine interface ghost cells.
     */
    bool d_consistent_type_2_bdry = false;

    /*!
     * Refine operator employed to fill coarse grid ghost cell values.
     */
    SAMRAIPointer<SAMRAI::xfer::RefineOperatorNd> d_refine_op =
        new SAMRAI::geom::CartesianSideDoubleConservativeLinearRefineNd();

    /*!
     * Cached hierarchy-related information.
     */
    SAMRAIPointer<SAMRAI::hier::PatchHierarchyNd> d_hierarchy;
    std::vector<SAMRAI::hier::CoarseFineBoundary<NDIM> > d_cf_boundary;
    SAMRAIPointer<SAMRAI::pdat::SideVariableNd<int> > d_sc_indicator_var =
        new SAMRAI::pdat::SideVariableNd<int>("CartSideDoubleQuadraticCFInterpolation::sc_indicator_var");
    int d_sc_indicator_idx = IBTK::invalid_index;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_CartSideDoubleQuadraticCFInterpolation
