// Filename: NormOps.h
// Created on 08 Dec 2008 by Boyce Griffith
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

#ifndef included_NormOps
#define included_NormOps

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI
/////////////////////////////// INCLUDES /////////////////////////////////////

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class NormOps provides functionality for computing discrete vector
 * norms.
 */
class NormOps
{
public:
    /*!
     * \brief Compute the discrete L1 norm of the SAMRAI vector.
     */
    static double L1Norm(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>* samrai_vector, bool local_only = false);

    /*!
     * \brief Compute the discrete L2 norm of the SAMRAI vector.
     */
    static double L2Norm(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>* samrai_vector, bool local_only = false);

    /*!
     * \brief Compute the discrete max-norm of the SAMRAI vector.
     */
    static double maxNorm(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>* samrai_vector, bool local_only = false);

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    NormOps();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    NormOps(const NormOps& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    NormOps& operator=(const NormOps& that);

    /*!
     * \brief Compute the local L1 norm.
     */
    static double L1Norm_local(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>* samrai_vector);

    /*!
     * \brief Compute the local L2 norm.
     */
    static double L2Norm_local(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>* samrai_vector);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_NormOps
