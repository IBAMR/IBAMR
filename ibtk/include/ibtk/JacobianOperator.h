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

#ifndef included_IBTK_JacobianOperator
#define included_IBTK_JacobianOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LinearOperator.h"

#include "tbox/Pointer.h"

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

namespace IBTK
{
/*!
 * \brief Class JacobianOperator provides an abstract interface for the
 * specification of general operators to compute Jacobian-vector products, i.e.,
 * \f$ F'[x]v \f$.
 */
class JacobianOperator : public LinearOperator
{
public:
    /*!
     * \brief Constructor.
     */
    JacobianOperator(std::string object_name);

    /*!
     * \brief Empty destructor.
     */
    ~JacobianOperator() = default;

    /*!
     * \name General Jacobian functionality.
     */
    //\{

    /*!
     * \brief Compute hierarchy dependent data required for evaluating F'[x].
     *
     * \param x value where the Jacobian is to be evaluated
     */
    virtual void formJacobian(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x) = 0;

    /*!
     * \brief Return the vector where the Jacobian is evaluated.
     *
     * \note Implementations of this member function are permitted to return a
     * NULL pointer if the operator is not initialized, or if formJacobian() has
     * not been called.
     */
    virtual SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > getBaseVector() const = 0;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    JacobianOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    JacobianOperator(const JacobianOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    JacobianOperator& operator=(const JacobianOperator& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_JacobianOperator
