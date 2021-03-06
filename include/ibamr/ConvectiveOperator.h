// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_ConvectiveOperator
#define included_IBAMR_ConvectiveOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/ibamr_enums.h"

#include "ibtk/GeneralOperator.h"
#include "ibtk/ibtk_utilities.h"

#include <string>

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
    ConvectiveOperator(std::string object_name, ConvectiveDifferencingType difference_form);

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
    void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
               SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    //\}

protected:
    /*!
     * Enumerated type that determines which form of differencing to use.
     */
    ConvectiveDifferencingType d_difference_form;

    /*!
     * The advection velocity patch data descriptor index.
     */
    int d_u_idx = IBTK::invalid_index;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    ConvectiveOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    ConvectiveOperator(const ConvectiveOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    ConvectiveOperator& operator=(const ConvectiveOperator& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_ConvectiveOperator
