// Filename: ExtendedRobinBcCoefStrategy.h
// Created on 16 May 2007 by Boyce Griffith
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

#ifndef included_ExtendedRobinBcCoefStrategy
#define included_ExtendedRobinBcCoefStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

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
    ExtendedRobinBcCoefStrategy();

    /*!
     * \brief Empty virtual destructor.
     */
    virtual ~ExtendedRobinBcCoefStrategy();

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
    int d_target_data_idx;

    /*
     * Whether to use homogeneous boundary conditions.
     */
    bool d_homogeneous_bc;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    ExtendedRobinBcCoefStrategy(const ExtendedRobinBcCoefStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    ExtendedRobinBcCoefStrategy& operator=(const ExtendedRobinBcCoefStrategy& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_ExtendedRobinBcCoefStrategy
