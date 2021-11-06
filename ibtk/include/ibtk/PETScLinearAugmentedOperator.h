// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2021 by the IBAMR developers
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

#ifndef included_IBTK_PETScLinearAugmentedOperator
#define included_IBTK_PETScLinearAugmentedOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LinearOperator.h"

#include "petscvec.h"

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
 * \brief Class PETScLinearAugmentedOperator provides an abstract interface for the specification of linear operators to
 * compute \f$y = Ax\f$ and \f$z=Ax+y\f$ where \f$x \f$ contains both Eulerian degrees of freedom and extra "augmented"
 * degrees of freedom. When used with PETScAugmentedKrylovLinearSolver, d_aug_x_vec is set before any call to a linear
 * operator function. Implementations should fill appropriate data in d_aug_y_vec.
 *
 * \see PETScAugmentedKyrlovLinearSolver
 */
class PETScLinearAugmentedOperator : public LinearOperator
{
public:
    /*!
     * \brief Constructor.
     */
    PETScLinearAugmentedOperator(std::string object_name, bool homogeneous_bc = false);

    /*!
     * \brief Empty destructor.
     */
    virtual ~PETScLinearAugmentedOperator();

    /*!
     * \brief Set the current augmented LHS vec.
     */
    virtual void setAugmentedVec(const Vec& vec);

    /*!
     * \brief Get the augmented RHS vec.
     */
    virtual const Vec& getAugmentedVec() const;

    /*!
     * \brief Set the augmented RHS vec.
     *
     * This function is used to modify the vector for boundary conditions. The modified vector is stored in d_aug_y_vec.
     */
    virtual void setAugmentedRhsForBcs(Vec& aug_y);

    /*!
     * \name Linear operator functionality.
     */
    //\{

    /*!
     * \brief Modify y to account for boundary conditions.
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
     * \note A default implementation evaluates y := y - A*0.
     */
    void modifyRhsForBcs(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    /*!
     * \brief compute \f$z = A[x] + y\f$.
     *
     * \note This function is not currently implemented and should not be used.
     */
    void applyAdd(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                  SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y,
                  SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& z) override;

    //\}
protected:
    Vec d_aug_x_vec = nullptr, d_aug_y_vec = nullptr, d_aug_rhs_y_vec = nullptr;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScLinearAugmentedOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScLinearAugmentedOperator(const PETScLinearAugmentedOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScLinearAugmentedOperator& operator=(const PETScLinearAugmentedOperator& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_PETScLinearAugmentedOperator
