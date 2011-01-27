// Filename: FACPreconditionerStrategy.h
// Created on 10 Sep 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#ifndef included_FACPreconditionerStrategy
#define included_FACPreconditionerStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/FACPreconditioner.h>

// SAMRAI INCLUDES
#include <tbox/ConstPointer.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class FACPreconditionerStrategy provides an interface for specifying
 * the problem-specific operations needed to implement a specific FAC
 * preconditioner.
 *
 * This class is similar to the SAMRAI class SAMRAI::solv::FACOperatorStrategy,
 * except that certain methods required by the SAMRAI::solv::FACOperatorStrategy
 * interface have been modified (specifically, smoothError()) or removed
 * (specifically, restrictSolution() and postprocessOneCycle()).  This interface
 * also requires the implementation of a method, prolongError(), that is
 * optimized for the case in which the FAC preconditioner algorithm is
 * configured not to use pre-smoothing sweeps.
 *
 * \see FACPreconditioner
 */
class FACPreconditionerStrategy
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Empty default constructor.
     */
    FACPreconditionerStrategy();

    /*!
     * \brief Empty virtual desctructor.
     */
    virtual
    ~FACPreconditionerStrategy();

    /*!
     * \brief Method to allow the FACPreconditioner object to register itself
     * with the concrete FACPreconditionerStrategy.
     *
     * \note A default empty implementation is provided.
     */
    virtual void
    setFACPreconditioner(
        SAMRAI::tbox::ConstPointer<FACPreconditioner> preconditioner);

    /*!
     * \brief Restrict the residual from the source vector to the destination
     * vector on the specified level of the patch hierarchy.
     *
     * \note Implementations must support the case in which source and dest are
     * the same vector.
     */
    virtual void
    restrictResidual(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& source,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& dest,
        int dest_level_num) = 0;

    /*!
     * \brief Prolong the error from the source vector to the destination
     * vector on the specified level of the patch hierarchy.
     *
     * \note Implementations must support the case in which source and dest are
     * the same vector.
     */
    virtual void
    prolongError(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& source,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& dest,
        int dest_level_num) = 0;

    /*!
     * \brief Prolong the error from the source vector to the destination vector
     * on the specified level of the patch hierarchy and correct the fine-grid
     * error.
     *
     * \note Implementations must support the case in which source and dest are
     * the same vector.
     */
    virtual void
    prolongErrorAndCorrect(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& source,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& dest,
        int dest_level_num) = 0;

    /*!
     * \brief Smooth the error by the specified number of sweeps on the
     * specified level of the patch hierarchy.
     */
    virtual void
    smoothError(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& error,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
        int level_num,
        int num_sweeps,
        bool performing_pre_sweeps,
        bool performing_post_sweeps) = 0;

    /*!
     * \brief Solve the system of equations on the coarsest level of the patch
     * hierarchy.
     *
     * \return true if the solver converged to specified tolerance, false otherwise
     */
    virtual bool
    solveCoarsestLevel(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& error,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
        int coarsest_level_num) = 0;

    /*!
     * \brief Compute the composite-grid residual on the specified level of the
     * patch hierarchy.
     */
    virtual void
    computeResidual(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& solution,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& rhs,
        int level_num) = 0;

    /*!
     * \brief Initialize any hierarchy-dependent data.
     *
     * \note A default empty implementation is provided.
     */
    virtual void
    initializeOperatorState(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& solution,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& rhs);

    /*!
     * \brief Deallocate any hierarchy-dependent data initialized by
     * initializeOperatorState().
     *
     * \note A default empty implementation is provided.
     */
    virtual void
    deallocateOperatorState();

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FACPreconditionerStrategy(
        const FACPreconditionerStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FACPreconditionerStrategy&
    operator=(
        const FACPreconditionerStrategy& that);
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/FACPreconditionerStrategy.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_FACPreconditionerStrategy
