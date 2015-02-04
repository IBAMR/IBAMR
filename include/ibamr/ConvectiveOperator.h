// Filename: ConvectiveOperator.h
// Created on 21 Aug 2011 by Boyce Griffith
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

#ifndef included_ConvectiveOperator
#define included_ConvectiveOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "ibamr/ibamr_enums.h"
#include "ibtk/GeneralOperator.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class ConvectiveOperator is a abstract class for an implementation of
 * a convective differencing operator.
 */
class ConvectiveOperator : public IBTK::GeneralOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    ConvectiveOperator(const std::string& object_name, ConvectiveDifferencingType difference_form);

    /*!
     * \brief Destructor.
     */
    ~ConvectiveOperator();

    /*!
     * \brief Set the patch data index corresponding to the advection velocity
     * to be used when computing the convective derivative.
     */
    void setAdvectionVelocity(int u_idx);

    /*!
     * \brief Get the patch data index corresponding to the advection velocity
     * used when computing the convective derivative.
     */
    int getAdvectionVelocity() const;

    /*!
     * \brief Set the convective differencing form to be used by the operator.
     */
    void setConvectiveDifferencingType(ConvectiveDifferencingType difference_form);

    /*!
     * \brief Get the convective differencing form used by the operator.
     */
    ConvectiveDifferencingType getConvectiveDifferencingType() const;

    /*!
     * \brief Compute N = u * grad Q.
     */
    virtual void applyConvectiveOperator(int Q_idx, int N_idx) = 0;

    /*!
     * \name General operator functionality.
     */
    //\{

    /*!
     * \brief Compute \f$y=F[x]\f$.
     *
     * Before calling apply(), the form of the vectors \a x and \a y should be
     * set properly by the user on all patch interiors on the specified range of
     * levels in the patch hierarchy.  The user is responsible for all data
     * management for the quantities associated with the vectors.  In
     * particular, patch data in these vectors must be allocated prior to
     * calling this method.
     *
     * \param x input vector
     * \param y output vector, i.e., \f$y=F[x]\f$
     *
     * <b>Conditions on Parameters:</b>
     * - vectors \a x and \a y must have same hierarchy
     * - vectors \a x and \a y must have same structure, depth, etc.
     *
     * In general, the vectors \a x and \a y \em cannot be the same.
     *
     * \see initializeOperatorState
     */
    void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x, SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y);

    //\}

protected:
    /*!
     * Enumerated type that determines which form of differencing to use.
     */
    ConvectiveDifferencingType d_difference_form;

    /*!
     * The advection velocity patch data descriptor index.
     */
    int d_u_idx;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    ConvectiveOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    ConvectiveOperator(const ConvectiveOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    ConvectiveOperator& operator=(const ConvectiveOperator& that);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_ConvectiveOperator
