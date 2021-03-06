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

#ifndef included_IBTK_RefinePatchStrategySet
#define included_IBTK_RefinePatchStrategySet

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "Box.h"
#include "BoxList.h"
#include "IntVector.h"
#include "RefinePatchStrategy.h"

#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class RefinePatchStrategySet is a utility class that allows multiple
 * SAMRAI::xfer::RefinePatchStrategy objects to be employed by a single
 * SAMRAI::xfer::RefineSchedule.
 */
class RefinePatchStrategySet : public SAMRAI::xfer::RefinePatchStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor.
     */
    template <typename InputIterator>
    RefinePatchStrategySet(InputIterator first, InputIterator last, bool managed = true)
        : d_strategy_set(first, last), d_managed(managed)
    {
        // intentionally blank
        return;
    } // RefinePatchStrategySet

    /*!
     * \brief Destructor.
     *
     * \note The patch strategy objects provided to the constructor are deleted
     * by this class destructor.
     */
    ~RefinePatchStrategySet();

    /*!
     * Function to set data associated with the given list of patch data indices
     * at patch boundaries that intersect the physical domain boundary.  The
     * specific boundary conditions are determined by the user.  The patch data
     * components set in this routine correspond to the scratch components
     * specified in calls to the registerRefine() function in the
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
                                       const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill) override;

    /*!
     * Function to return maximum stencil width needed over user-defined data
     * interpolation operations.  This is needed to determine the correct
     * interpolation data dependencies.
     */
    SAMRAI::hier::IntVector<NDIM> getRefineOpStencilWidth() const override;

    /*!
     * Function to perform user-defined preprocess data refine operations.  This
     * member function is called before standard refine operations (expressed
     * using concrete subclasses of the SAMRAI::xfer::RefineOperator base
     * class).  The preprocess function refines data from the scratch components
     * of the coarse patch into the scratch components of the fine patch on the
     * specified fine box region.  Recall that the scratch components are
     * specified in calls to the registerRefine() function in the
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
                          const SAMRAI::hier::IntVector<NDIM>& ratio) override;

    /*!
     * Function to perform user-defined preprocess data refine operations.  This
     * member function is called after standard refine operations (expressed
     * using concrete subclasses of the SAMRAI::xfer::RefineOperator base
     * class).  The postprocess function refines data from the scratch
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
                           const SAMRAI::hier::IntVector<NDIM>& ratio) override;

    /*!
     * Function to perform user-defined refine operations.  This member function
     * is called before standard refining operations (expressed using concrete
     * subclasses of the SAMRAI::xfer::RefineOperator base class).  The
     * preprocess function refines data from the scratch components of the
     * coarse patch into the scratch components of the fine patch on the
     * specified fine box regions.
     *
     * \param fine        Fine patch containing destination data.
     * \param coarse      Coarse patch containing source data.
     * \param fine_boxes  List of box regions on fine patch into which data is refined.
     * \param ratio       Integer vector containing ratio relating index space between coarse
     *and
     *fine patches.
     */
    void preprocessRefineBoxes(SAMRAI::hier::Patch<NDIM>& fine,
                               const SAMRAI::hier::Patch<NDIM>& coarse,
                               const SAMRAI::hier::BoxList<NDIM>& fine_boxes,
                               const SAMRAI::hier::IntVector<NDIM>& ratio) override;

    /*!
     * Function to perform user-defined refine operations.  This member function
     * is called after standard refining operations (expressed using concrete
     * subclasses of the SAMRAI::xfer::RefineOperator base class).  The
     * postprocess function refines data from the scratch components of the
     * coarse patch into the scratch components of the fine patch on the
     * specified fine box regions.
     *
     * \param fine        Fine patch containing destination data.
     * \param coarse      Coarse patch containing source data.
     * \param fine_boxes  List of box regions on fine patch into which data is refined.
     * \param ratio       Integer vector containing ratio relating index space between coarse
     *and
     *fine patches.
     */
    void postprocessRefineBoxes(SAMRAI::hier::Patch<NDIM>& fine,
                                const SAMRAI::hier::Patch<NDIM>& coarse,
                                const SAMRAI::hier::BoxList<NDIM>& fine_boxes,
                                const SAMRAI::hier::IntVector<NDIM>& ratio) override;

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    RefinePatchStrategySet() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    RefinePatchStrategySet(const RefinePatchStrategySet& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    RefinePatchStrategySet& operator=(const RefinePatchStrategySet& that) = delete;

    /*!
     * \brief The set of SAMRAI::xfer:RefinePatchStrategy objects.
     */
    std::vector<SAMRAI::xfer::RefinePatchStrategy<NDIM>*> d_strategy_set;

    /*!
     * \brief Boolean value that indicates whether this class should provide
     * memory management for the strategy objects.
     */
    const bool d_managed;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_RefinePatchStrategySet
