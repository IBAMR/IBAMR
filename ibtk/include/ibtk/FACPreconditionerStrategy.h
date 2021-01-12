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

#ifndef included_IBTK_FACPreconditionerStrategy
#define included_IBTK_FACPreconditionerStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/FACPreconditioner.h"

#include "tbox/ConstPointer.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

#include <iosfwd>
#include <string>
#include <utility>

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
class FACPreconditionerStrategy : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     */
    FACPreconditionerStrategy(std::string object_name, bool homogeneous_bc = false);

    /*!
     * \brief Empty virtual desctructor.
     */
    virtual ~FACPreconditionerStrategy() = default;

    /*!
     * \brief Return the object name.
     */
    const std::string& getName() const;

    /*!
     * \brief Return whether the operator is initialized.
     */
    virtual bool getIsInitialized() const;

    /*!
     * \brief Method to allow the FACPreconditioner object to register itself
     * with the concrete FACPreconditionerStrategy.
     */
    virtual void setFACPreconditioner(SAMRAI::tbox::ConstPointer<FACPreconditioner> preconditioner);

    /*!
     * \brief Set whether the solver should use homogeneous boundary conditions.
     */
    virtual void setHomogeneousBc(bool homogeneous_bc);

    /*!
     * \brief Return whether the solver is using homogeneous boundary
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
     * \brief Zero-out the provided vector on the specified level of the patch
     * hierarchy.
     */
    virtual void setToZero(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error, int level_num) = 0;

    /*!
     * \brief Restrict the residual from the source vector to the destination
     * vector on the specified level of the patch hierarchy.
     *
     * \note Implementations must support the case in which source and dest are
     * the same vector.
     */
    virtual void restrictResidual(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& source,
                                  SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& dest,
                                  int dest_level_num) = 0;

    /*!
     * \brief Prolong the error from the source vector to the destination
     * vector on the specified level of the patch hierarchy.
     *
     * \note Implementations must support the case in which source and dest are
     * the same vector.
     */
    virtual void prolongError(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& source,
                              SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& dest,
                              int dest_level_num) = 0;

    /*!
     * \brief Prolong the error from the source vector to the destination vector
     * on the specified level of the patch hierarchy and correct the fine-grid
     * error.
     *
     * \note Implementations must support the case in which source and dest are
     * the same vector.
     */
    virtual void prolongErrorAndCorrect(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& source,
                                        SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& dest,
                                        int dest_level_num) = 0;

    /*!
     * \brief Smooth the error by the specified number of sweeps on the
     * specified level of the patch hierarchy.
     */
    virtual void smoothError(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
                             const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
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
    virtual bool solveCoarsestLevel(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
                                    const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                                    int coarsest_level_num) = 0;

    /*!
     * \brief Compute the composite-grid residual on the specified range of
     * levels of the patch hierarchy.
     */
    virtual void computeResidual(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs,
                                 int coarsest_level_num,
                                 int finest_level_num) = 0;

    /*!
     * \brief Initialize any hierarchy-dependent data.
     */
    virtual void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs);

    /*!
     * \brief Deallocate any hierarchy-dependent data initialized by
     * initializeOperatorState().
     */
    virtual void deallocateOperatorState();

    /*!
     * \brief Allocate scratch data.
     */
    virtual void allocateScratchData();

    /*!
     * \brief Deallocate scratch data.
     */
    virtual void deallocateScratchData();

    /*!
     * \name Logging functions.
     */
    //\{

    /*!
     * \brief Print class data to stream.
     */
    virtual void printClassData(std::ostream& stream);

    //\}

protected:
    /*!
     * \brief Return a SAMRAIVectorReal object that corresponds to the given
     * object but restricted to a single level of the patch hierarchy.
     */
    virtual SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >
    getLevelSAMRAIVectorReal(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& vec, int level_num) const;

    // Pointer to the FACPreconditioner that is using this operator.
    SAMRAI::tbox::ConstPointer<IBTK::FACPreconditioner> d_preconditioner;

    // Object name.
    const std::string d_object_name;

    // Boolean value to indicate whether the preconditioner is presently
    // initialized.
    bool d_is_initialized = false;

    // Solver configuration.
    bool d_homogeneous_bc;
    double d_solution_time = std::numeric_limits<double>::quiet_NaN(),
           d_current_time = std::numeric_limits<double>::quiet_NaN(),
           d_new_time = std::numeric_limits<double>::quiet_NaN();

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    FACPreconditionerStrategy() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FACPreconditionerStrategy(const FACPreconditionerStrategy& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FACPreconditionerStrategy& operator=(const FACPreconditionerStrategy& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_FACPreconditionerStrategy
