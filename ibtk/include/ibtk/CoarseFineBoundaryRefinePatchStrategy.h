// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_CoarseFineBoundaryRefinePatchStrategy
#define included_IBTK_CoarseFineBoundaryRefinePatchStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "Box.h"
#include "ComponentSelector.h"
#include "IntVector.h"
#include "RefinePatchStrategy.h"

#include <set>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace tbox
{
template <class TYPE>
class Pointer;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CoarseFineBoundaryRefinePatchStrategy is a subclass of the
 * abstract base class SAMRAI::xfer::RefinePatchStrategy that extends the
 * functionality of SAMRAI::xfer::RefinePatchStrategy to facilitate the
 * implementation of coarse-fine interface discretizations.
 */
class CoarseFineBoundaryRefinePatchStrategy : public SAMRAI::xfer::RefinePatchStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor.
     */
    CoarseFineBoundaryRefinePatchStrategy() = default;

    /*!
     * \brief Destructor.
     */
    virtual ~CoarseFineBoundaryRefinePatchStrategy() = default;

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
     * \param patch                Patch on which to fill boundary data.
     * \param fill_time            Double simulation time for boundary filling.
     * \param ghost_width_to_fill  Integer vector describing maximum ghost width to fill over
     *all
     *registered scratch components.
     */
    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch<NDIM>& patch,
                                       double fill_time,
                                       const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill) override = 0;

    /*!
     * Function to return maximum stencil width needed over user-defined data
     * interpolation operations.  This is needed to determine the correct
     * interpolation data dependencies.
     */
    SAMRAI::hier::IntVector<NDIM> getRefineOpStencilWidth() const override = 0;

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
     * \param fine      Fine patch containing destination data.
     * \param coarse    Coarse patch containing source data.
     * \param fine_box  Box region on fine patch into which data is refined.
     * \param ratio     Integer vector containing ratio relating index space between coarse and
     *fine
     *patches.
     */
    void preprocessRefine(SAMRAI::hier::Patch<NDIM>& fine,
                          const SAMRAI::hier::Patch<NDIM>& coarse,
                          const SAMRAI::hier::Box<NDIM>& fine_box,
                          const SAMRAI::hier::IntVector<NDIM>& ratio) override = 0;

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
     * \param fine      Fine patch containing destination data.
     * \param coarse    Coarse patch containing source data.
     * \param fine_box  Box region on fine patch into which data is refined.
     * \param ratio     Integer vector containing ratio relating index space between coarse and
     *fine
     *patches.
     */
    void postprocessRefine(SAMRAI::hier::Patch<NDIM>& fine,
                           const SAMRAI::hier::Patch<NDIM>& coarse,
                           const SAMRAI::hier::Box<NDIM>& fine_box,
                           const SAMRAI::hier::IntVector<NDIM>& ratio) override = 0;

    //\}

    /*!
     * \name Extension of SAMRAI::xfer::RefinePatchStrategy interface to support more
     * complex coarse-fine interface discretizations.
     */
    //\{

    /*!
     * Whether or not to employ a consistent interpolation scheme at "Type 2"
     * coarse-fine interface ghost cells.
     *
     * \note This subclasses may choose not to support a consistent "Type 2"
     * coarse-fine interface ghost cell interpolation scheme.
     */
    virtual void setConsistentInterpolationScheme(bool consistent_type_2_bdry) = 0;

    /*!
     * \brief Reset the patch data index operated upon by this class.
     */
    virtual void setPatchDataIndex(int patch_data_index) = 0;

    /*!
     * \brief Reset the patch data indices operated upon by this class.
     */
    virtual void setPatchDataIndices(const std::set<int>& patch_data_indices) = 0;

    /*!
     * \brief Reset the patch data indices operated upon by this class.
     */
    virtual void setPatchDataIndices(const SAMRAI::hier::ComponentSelector& patch_data_indices) = 0;

    /*!
     * Set the patch hierarchy used in constructing coarse-fine interface
     * boundary boxes.
     */
    virtual void setPatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy) = 0;

    /*!
     * Clear the patch hierarchy used in constructing coarse-fine interface
     * boundary boxes.
     */
    virtual void clearPatchHierarchy() = 0;

    /*!
     * Compute the normal extension of fine data at coarse-fine interfaces.
     */
    virtual void computeNormalExtension(SAMRAI::hier::Patch<NDIM>& patch,
                                        const SAMRAI::hier::IntVector<NDIM>& ratio,
                                        const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill) = 0;

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
    CoarseFineBoundaryRefinePatchStrategy(const CoarseFineBoundaryRefinePatchStrategy& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CoarseFineBoundaryRefinePatchStrategy& operator=(const CoarseFineBoundaryRefinePatchStrategy& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_CoarseFineBoundaryRefinePatchStrategy
