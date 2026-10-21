// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2023 by the IBAMR developers
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

#include <optional>

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
 *
 * By default, the L1 and L2 norms add the per-patch contributions in traversal order and combine the processes with a
 * single sum reduction. The option <code>-ibtk_sorted_norm_summation true</code> selects instead a summation in
 * ascending order of the per-patch and per-process contributions, at the cost of gathering every process's
 * contribution. For a fixed set of patches, the norms then do not depend on the order of the patches or processes or
 * on the order of the reduction. The contributions are rounded sums over the cells of each patch, so a different
 * decomposition of the same cells into patches can still change the last bits of the norms. The option is read from
 * the PETSc options database on the first use of a norm, and setSortedSummation() overrides it.
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

    /*!
     * \brief Set whether the L1 and L2 norms sum their per-patch and per-process contributions in ascending order,
     * overriding the option <code>-ibtk_sorted_norm_summation</code>.
     */
    static void setSortedSummation(bool sorted_summation);

    /*!
     * \brief Return whether the L1 and L2 norms use sorted summation.
     */
    static bool getSortedSummation();

    /*!
     * \brief Return the L2 norm given the sum of squares over this process's data, as L2Norm() does with plain
     * summation: the sum over the processes, unless local_only is true, followed by a square root.
     */
    static double L2NormFromSumOfSquares(double local_sum_of_squares, bool local_only = false);

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

    /*!
     * \brief Compute the local sum of the per-patch L1 norms in traversal order.
     */
    static double L1Norm_local_unsorted(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>* samrai_vector);

    /*!
     * \brief Compute the local sum of squares in traversal order.
     *
     * The sum is accumulated in a single variable, in the following order: the components of the vector in order; for
     * each component, the levels from the coarsest to the finest level of the vector; for each level, the patches in
     * the iteration order of the level; for each patch, the dot product of the component's data with itself over the
     * interior of the patch, weighted by the control volume of the component if the vector has one, as computed by the
     * patch norm operations of SAMRAI, which for side-, face-, and edge-centered data add the dot products of the
     * coordinate directions in order, skipping the directions that side-centered data do not allocate. Components of
     * every Cartesian data centering contribute.
     */
    static double L2NormSquared_local_unsorted(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>* samrai_vector);

    /*!
     * \brief Whether sorted summation is used; unset until the option is read or setSortedSummation() is called.
     */
    static std::optional<bool> s_sorted_summation;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_NormOps
