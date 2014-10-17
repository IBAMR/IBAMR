// Filename: JacobianOperator.h
// Created on 20 Apr 2005 by Boyce Griffith
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

#ifndef included_JacobianOperator
#define included_JacobianOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "ibtk/LinearOperator.h"
#include "tbox/Pointer.h"

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
    JacobianOperator(const std::string& object_name);

    /*!
     * \brief Empty destructor.
     */
    ~JacobianOperator();

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
    JacobianOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    JacobianOperator(const JacobianOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    JacobianOperator& operator=(const JacobianOperator& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_JacobianOperator
