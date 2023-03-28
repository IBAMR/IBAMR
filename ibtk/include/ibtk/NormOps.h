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

#ifndef included_IBTK_NormOps
#define included_IBTK_NormOps

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

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
    NormOps() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    NormOps(const NormOps& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    NormOps& operator=(const NormOps& that) = delete;

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

#endif //#ifndef included_IBTK_NormOps
