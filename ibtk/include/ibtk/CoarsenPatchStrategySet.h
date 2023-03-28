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

#ifndef included_IBTK_CoarsenPatchStrategySet
#define included_IBTK_CoarsenPatchStrategySet

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "Box.h"
#include "CoarsenPatchStrategy.h"
#include "IntVector.h"

#include <vector>

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
    SAMRAI::hier::IntVector<NDIM> getCoarsenOpStencilWidth() const override;

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
                           const SAMRAI::hier::IntVector<NDIM>& ratio) override;

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
                            const SAMRAI::hier::IntVector<NDIM>& ratio) override;

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CoarsenPatchStrategySet() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CoarsenPatchStrategySet(const CoarsenPatchStrategySet& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CoarsenPatchStrategySet& operator=(const CoarsenPatchStrategySet& that) = delete;

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

#endif //#ifndef included_IBTK_CoarsenPatchStrategySet
