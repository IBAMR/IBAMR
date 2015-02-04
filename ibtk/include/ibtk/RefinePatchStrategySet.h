// Filename: RefinePatchStrategySet.h
// Created on 11 Sep 2006 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_RefinePatchStrategySet
#define included_RefinePatchStrategySet

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

#include "Box.h"
#include "BoxList.h"
#include "IntVector.h"
#include "RefinePatchStrategy.h"

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
                                       const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill);

    /*!
     * Function to return maximum stencil width needed over user-defined data
     * interpolation operations.  This is needed to determine the correct
     * interpolation data dependencies.
     */
    SAMRAI::hier::IntVector<NDIM> getRefineOpStencilWidth() const;

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
                          const SAMRAI::hier::IntVector<NDIM>& ratio);

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
                           const SAMRAI::hier::IntVector<NDIM>& ratio);

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
                               const SAMRAI::hier::IntVector<NDIM>& ratio);

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
                                const SAMRAI::hier::IntVector<NDIM>& ratio);

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    RefinePatchStrategySet();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    RefinePatchStrategySet(const RefinePatchStrategySet& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    RefinePatchStrategySet& operator=(const RefinePatchStrategySet& that);

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

#endif //#ifndef included_RefinePatchStrategySet
