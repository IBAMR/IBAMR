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

#ifndef included_IBTK_RefinePatchStrategySet
#define included_IBTK_RefinePatchStrategySet

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIBox.h"
#include "SAMRAIBoxList.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPatch.h"
#include "SAMRAIRefinePatchStrategy.h"

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
class RefinePatchStrategySet : public SAMRAIRefinePatchStrategy
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
    void setPhysicalBoundaryConditions(SAMRAIPatch& patch,
                                       double fill_time,
                                       const SAMRAIIntVector& ghost_width_to_fill) override;

    /*!
     * Function to return maximum stencil width needed over user-defined data
     * interpolation operations.  This is needed to determine the correct
     * interpolation data dependencies.
     */
    SAMRAIIntVector getRefineOpStencilWidth() const override;

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
    void preprocessRefine(SAMRAIPatch& fine,
                          const SAMRAIPatch& coarse,
                          const SAMRAIBox& fine_box,
                          const SAMRAIIntVector& ratio) override;

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
    void postprocessRefine(SAMRAIPatch& fine,
                           const SAMRAIPatch& coarse,
                           const SAMRAIBox& fine_box,
                           const SAMRAIIntVector& ratio) override;

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
    void preprocessRefineBoxes(SAMRAIPatch& fine,
                               const SAMRAIPatch& coarse,
                               const SAMRAIBoxList& fine_boxes,
                               const SAMRAIIntVector& ratio) override;

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
    void postprocessRefineBoxes(SAMRAIPatch& fine,
                                const SAMRAIPatch& coarse,
                                const SAMRAIBoxList& fine_boxes,
                                const SAMRAIIntVector& ratio) override;

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
    std::vector<SAMRAIRefinePatchStrategy*> d_strategy_set;

    /*!
     * \brief Boolean value that indicates whether this class should provide
     * memory management for the strategy objects.
     */
    const bool d_managed;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_RefinePatchStrategySet
