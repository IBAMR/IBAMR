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

#ifndef included_IBTK_ExtendedRobinBcCoefStrategy
#define included_IBTK_ExtendedRobinBcCoefStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/ibtk_utilities.h"

#include "RobinBcCoefStrategy.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class ExtendedRobinBcCoefStrategy is a subclass of the abstract base
 * class SAMRAI::solv::RobinBcCoefStrategy that extends the functionality of
 * SAMRAI::solv::RobinBcCoefStrategy to allow for the specification of patch
 * data descriptor indices that are required for filling, and the specification
 * of whether homogeneous or inhomogeneous boundary data should be set.
 */
class ExtendedRobinBcCoefStrategy : public SAMRAI::solv::RobinBcCoefStrategy<NDIM>
{
public:
    /*!
     * \brief Empty default constructor.
     */
    ExtendedRobinBcCoefStrategy() = default;

    /*!
     * \brief Empty virtual destructor.
     */
    virtual ~ExtendedRobinBcCoefStrategy() = default;

    /*!
     * \name Extended SAMRAI::solv::RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Set the target data index.
     */
    virtual void setTargetPatchDataIndex(int target_data_idx);

    /*!
     * \brief Clear the target data index.
     */
    virtual void clearTargetPatchDataIndex();

    /*!
     * \brief Set whether the class is filling homogeneous or inhomogeneous
     * boundary conditions.
     */
    virtual void setHomogeneousBc(bool homogeneous_bc);

    //\}

protected:
    /*
     * The patch data index corresponding to the data to be filled.
     */
    int d_target_data_idx = IBTK::invalid_index;

    /*
     * Whether to use homogeneous boundary conditions.
     */
    bool d_homogeneous_bc = false;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    ExtendedRobinBcCoefStrategy(const ExtendedRobinBcCoefStrategy& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    ExtendedRobinBcCoefStrategy& operator=(const ExtendedRobinBcCoefStrategy& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_ExtendedRobinBcCoefStrategy
