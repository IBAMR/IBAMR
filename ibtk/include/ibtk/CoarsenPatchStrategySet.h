// Filename: CoarsenPatchStrategySet.h
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

#ifndef included_CoarsenPatchStrategySet
#define included_CoarsenPatchStrategySet

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

#include "CoarsenPatchStrategy.h"
#include "IntVector.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Box;
template <int DIM>
class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CoarsenPatchStrategySet is a utility class that allows multiple
 * SAMRAI::xfer::CoarsenPatchStrategy objects to be employed by a single
 * SAMRAI::xfer::CoarsenSchedule.
 */
class CoarsenPatchStrategySet : public SAMRAI::xfer::CoarsenPatchStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor.
     */
    template <typename InputIterator>
    CoarsenPatchStrategySet(InputIterator first, InputIterator last, bool managed = true)
        : d_strategy_set(first, last), d_managed(managed)
    {
        // intentionally blank
        return;
    } // CoarsenPatchStrategySet

    /*!
     * \brief Destructor.
     */
    ~CoarsenPatchStrategySet();

    /*!
     * Return maximum stencil width needed over all user-defined data coarsening
     * operations.  This is needed to determine the correct coarsening data
     * dependencies.
     */
    SAMRAI::hier::IntVector<NDIM> getCoarsenOpStencilWidth() const;

    /*!
     * Perform user-defined coarsening operations.  This member function is
     * called before standard coarsening operations (expressed using concrete
     * subclasses of the SAMRAI::xfer::CoarsenOperator base class).  The
     * preprocess function moves data from the source components on the fine
     * patch into the source components on the coarse patch in the specified
     * coarse box region.  Recall that the source components are specified in
     * calls to the registerCoarsen() function in the
     * SAMRAI::xfer::CoarsenAlgorithm class.
     *
     * \param coarse      Coarse patch containing destination data.
     * \param fine        Fine patch containing source data.
     * \param coarse_box  Box region on coarse patch into which data is coarsened.
     * \param ratio       Integer vector containing ratio relating index space between coarse
     *and
     *fine patches.
     */
    void preprocessCoarsen(SAMRAI::hier::Patch<NDIM>& coarse,
                           const SAMRAI::hier::Patch<NDIM>& fine,
                           const SAMRAI::hier::Box<NDIM>& coarse_box,
                           const SAMRAI::hier::IntVector<NDIM>& ratio);

    /*!
     * Perform user-defined coarsening operations.  This member function is
     * called after standard coarsening operations (expressed using concrete
     * subclasses of the SAMRAI::xfer::CoarsenOperator base class).  The
     * postprocess function moves data from the source components on the fine
     * patch into the source components on the coarse patch in the specified
     * coarse box region.  Recall that the source components are specified in
     * calls to the registerCoarsen() function in the
     * SAMRAI::xfer::CoarsenAlgorithm class.
     *
     * \param coarse      Coarse patch containing destination data.
     * \param fine        Fine patch containing source data.
     * \param coarse_box  Box region on coarse patch into which data is copied.
     * \param ratio       Integer vector containing ratio
     */
    void postprocessCoarsen(SAMRAI::hier::Patch<NDIM>& coarse,
                            const SAMRAI::hier::Patch<NDIM>& fine,
                            const SAMRAI::hier::Box<NDIM>& coarse_box,
                            const SAMRAI::hier::IntVector<NDIM>& ratio);

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CoarsenPatchStrategySet();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CoarsenPatchStrategySet(const CoarsenPatchStrategySet& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CoarsenPatchStrategySet& operator=(const CoarsenPatchStrategySet& that);

    /*!
     * \brief The set of SAMRAI::xfer:CoarsenPatchStrategy objects.
     */
    std::vector<SAMRAI::xfer::CoarsenPatchStrategy<NDIM>*> d_strategy_set;

    /*!
     * \brief Boolean value that indicates whether this class should provide
     * memory management for the strategy objects.
     */
    const bool d_managed;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CoarsenPatchStrategySet
