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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/PETScLinearAugmentedOperator.h"

#include "Box.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Pointer.h"

#include <string>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScLinearAugmentedOperator::PETScLinearAugmentedOperator(std::string object_name, bool homogeneous_bc)
    : LinearOperator(std::move(object_name), homogeneous_bc)
{
    // intentionally blank
    return;
} // PETScLinearAugmentedOperator()

PETScLinearAugmentedOperator::~PETScLinearAugmentedOperator()
{
    deallocateOperatorState();
    return;
} // ~PETScLinearAugmentedOperator()

/*!
 * \brief Set the current augmented LHS vec.
 */
void
PETScLinearAugmentedOperator::setAugmentedVec(const Vec& vec)
{
    d_aug_x_vec = vec;
}

/*!
 * \brief Get the augmented RHS vec.
 */
const Vec&
PETScLinearAugmentedOperator::getAugmentedVec() const
{
    return d_aug_b_vec;
}

/*!
 * \brief Set the augmented RHS vec.
 *
 * This function is used to modify the vector for boundary conditions. The modified vector is stored in d_aug_b_vec.
 */
void
PETScLinearAugmentedOperator::setAugmentedRhsForBcs(Vec& aug_y)
{
    d_aug_x_vec = aug_y;
}

void
PETScLinearAugmentedOperator::modifyRhsForBcs(SAMRAIVectorReal<NDIM, double>& y)
{
    if (d_homogeneous_bc) return;

    // Set y := y - A*0, i.e., shift the right-hand-side vector to account for
    // inhomogeneous boundary conditions.
    // Prepare copies for Cartesian data
    Pointer<SAMRAIVectorReal<NDIM, double> > x = y.cloneVector("");
    Pointer<SAMRAIVectorReal<NDIM, double> > b = y.cloneVector("");
    x->allocateVectorData();
    b->allocateVectorData();
    x->setToScalar(0.0);
    // Prepare copies for augmented data.
    Vec temp_x_copy, temp_b_copy;
    int ierr = VecCopy(temp_x_copy, d_aug_x_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecCopy(temp_b_copy, d_aug_b_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecScale(d_aug_x_vec, 0.0);
    IBTK_CHKERRQ(ierr);
    // Apply the operator. Note that apply() also works on d_aug_x_vec and d_aug_b_vec.
    apply(*x, *b);
    // Now subtract y from the result.
    y.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&y, false), b);
    // And subtract d_aug_b_vec from temp_x_copy.
    ierr = VecAYPX(d_aug_b_vec, -1.0, temp_b_copy);
    IBTK_CHKERRQ(ierr);
    // Now reset and free the vectors
    x->freeVectorComponents();
    b->freeVectorComponents();
    ierr = VecCopy(temp_x_copy, d_aug_x_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&temp_b_copy);
    IBTK_CHKERRQ(ierr);
    return;
} // modifyRhsForBcs

void
PETScLinearAugmentedOperator::applyAdd(SAMRAIVectorReal<NDIM, double>& x,
                                       SAMRAIVectorReal<NDIM, double>& y,
                                       SAMRAIVectorReal<NDIM, double>& z)
{
    TBOX_ERROR(d_object_name + "::applyAdd() is not currently implemented.\n");
    return;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
