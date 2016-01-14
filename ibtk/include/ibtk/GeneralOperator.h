// Filename: GeneralOperator.h
// Created on 18 Nov 2003 by Boyce Griffith
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

#ifndef included_GeneralOperator
#define included_GeneralOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <iosfwd>
#include <string>
#include <utility>

#include "ibtk/HierarchyMathOps.h"
#include "tbox/DescribedClass.h"
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
 * \brief Class GeneralOperator provides an abstract interface for the
 * specification of general operators to compute \f$ y=F[x] \f$ and \f$ z=F[x]+y
 * \f$.
 */
class GeneralOperator : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     */
    GeneralOperator(const std::string& object_name, bool homogeneous_bc = false);

    /*!
     * \brief Empty virtual destructor.
     */
    virtual ~GeneralOperator();

    /*!
     * \name General operator functionality.
     */
    //\{

    /*!
     * \brief Return the object name.
     */
    const std::string& getName() const;

    /*!
     * \brief Return whether the operator is initialized.
     */
    virtual bool getIsInitialized() const;

    /*!
     * \brief Set whether the operator should use homogeneous boundary
     * conditions.
     */
    virtual void setHomogeneousBc(bool homogeneous_bc);

    /*!
     * \brief Return whether the operator is using homogeneous boundary
     * conditions.
     */
    virtual bool getHomogeneousBc() const;

    /*!
     * \brief Set the time at which the solution is to be evaluated.
     */
    virtual void setSolutionTime(double solution_time);

    /*!
     * \brief Get the time at which the solution is being evaluated.
     */
    virtual double getSolutionTime() const;

    /*!
     * \brief Set the current time interval.
     */
    virtual void setTimeInterval(double current_time, double new_time);

    /*!
     * \brief Get the current time interval.
     */
    virtual std::pair<double, double> getTimeInterval() const;

    /*!
     * \brief Get the current time step size.
     */
    virtual double getDt() const;

    /*!
     * \brief Set the HierarchyMathOps object used by the operator.
     */
    virtual void setHierarchyMathOps(SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops);

    /*!
     * \brief Get the HierarchyMathOps object used by the operator.
     */
    virtual SAMRAI::tbox::Pointer<HierarchyMathOps> getHierarchyMathOps() const;

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
     * \note Subclasses may require that the operator be initialized prior to
     * calling apply().
     *
     * \see initializeOperatorState
     */
    virtual void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                       SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) = 0;

    /*!
     * \brief Compute \f$z=F[x]+y\f$.
     *
     * Before calling applyAdd(), the form of the vectors \a x, \a y, and \a z
     * should be set properly by the user on all patch interiors on the
     * specified range of levels in the patch hierarchy.  The user is
     * responsible for all data management for the quantities associated with
     * the vectors.  In particular, patch data in these vectors must be
     * allocated prior to calling this method.
     *
     * \param x input vector
     * \param y input vector
     * \param z output vector, i.e., \f$z=F[x]+y\f$
     *
     * <b>Conditions on Parameters:</b>
     * - vectors \a x, \a y, and \a z must have same hierarchy
     * - vectors \a x, \a y, and \a z must have same structure, depth, etc.
     *
     * In general, the vectors \a x and \a y \em cannot be the same.
     *
     * \note Subclasses may require that the operator be initialized prior to
     * calling applyAdd().
     *
     * \see initializeOperatorState
     *
     * \note A default implementation is provided which employs apply() and
     * SAMRAI::solv::SAMRAIVectorReal::add().
     */
    virtual void applyAdd(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                          SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y,
                          SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& z);

    /*!
     * \brief Compute hierarchy dependent data required for computing y=F[x] and
     * z=F[x]+y.
     *
     * The vector arguments for apply(), applyAdd(), etc, need not match those
     * for initializeOperatorState().  However, there must be a certain degree
     * of similarity, including
     * - hierarchy configuration (hierarchy pointer and level range)
     * - number, type and alignment of vector component data
     * - ghost cell widths of data in the input and output vectors
     *
     * \note It is generally necessary to reinitialize the operator state when
     * the hierarchy configuration changes.
     *
     * It is safe to call initializeOperatorState() when the state is already
     * initialized.  In this case, the operator state is first deallocated and
     * then reinitialized.
     *
     * Conditions on arguments:
     * - input and output vectors must have same hierarchy
     * - input and output vectors must have same structure, depth, etc.
     *
     * Call deallocateOperatorState() to remove any data allocated by this
     * method.
     *
     * \see deallocateOperatorState
     *
     * \param in input vector
     * \param out output vector
     */
    virtual void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeOperatorState().
     *
     * \note Subclasses are required to be implemented so that it is safe to
     * call deallocateOperatorState() when the operator state is already
     * deallocated.
     *
     * \see initializeOperatorState
     */
    virtual void deallocateOperatorState();

    /*!
     * \brief Modify the RHS vector to account for boundary conditions.
     *
     * Before calling this function, the form of the vector y should be set
     * properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in this vector should be allocated.
     * The user is responsible for managing the storage for the vectors.
     *
     * \note The operator MUST be initialized prior to calling modifyRhsForBcs.
     *
     * \see initializeOperatorState
     *
     * \note A default implementation does not modify the RHS vector y.
     */
    virtual void modifyRhsForBcs(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y);

    /*!
     * \brief Impose boudary conditions in the solution vector.
     *
     * Before calling this function, the form of the vector y should be set
     * properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in this vector should be allocated.
     * The user is responsible for managing the storage for the vectors.
     *
     * \note The operator MUST be initialized prior to calling imposeSolBcs.
     *
     * \see initializeOperatorState
     *
     * \note A default implementation does not modify the sol vector u.
     */
    virtual void imposeSolBcs(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& u);

    //\}

    /*!
     * \name Logging functions.
     */
    //\{

    /*!
     * \brief Enable or disable logging.
     */
    virtual void setLoggingEnabled(bool enable_logging = true);

    /*!
     * \brief Determine whether logging is enabled or disabled.
     */
    virtual bool getLoggingEnabled() const;

    /*!
     * \brief Print class data to stream.
     */
    virtual void printClassData(std::ostream& stream);

    //\}

protected:
    // Object name.
    const std::string d_object_name;

    // Boolean value to indicate whether the preconditioner is presently
    // initialized.
    bool d_is_initialized;

    // Operator configuration.
    bool d_homogeneous_bc;
    double d_solution_time, d_current_time, d_new_time;

    // Mathematical operators.
    SAMRAI::tbox::Pointer<HierarchyMathOps> d_hier_math_ops;
    bool d_hier_math_ops_external;

    // Logging configuration.
    bool d_enable_logging;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    GeneralOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    GeneralOperator(const GeneralOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    GeneralOperator& operator=(const GeneralOperator& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_GeneralOperator
