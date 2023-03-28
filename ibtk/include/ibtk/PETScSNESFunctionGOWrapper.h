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

#ifndef included_IBTK_PETScSNESFunctionGOWrapper
#define included_IBTK_PETScSNESFunctionGOWrapper

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/GeneralOperator.h"

#include "IntVector.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Pointer.h"

#include "petscsnes.h"
#include "petscsys.h"
#include "petscvec.h"

#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PETScSNESFunctionGOWrapper provides a GeneralOperator interface
 * for a <A HREF="http://www.mcs.anl.gov/petsc">PETSc</A> SNES nonlinear
 * function.
 */
class PETScSNESFunctionGOWrapper : public GeneralOperator
{
public:
    /*!
     * \brief Constructor.
     *
     * Construct a general operator wrapper corresponding to the provided PETSc
     * SNES function object.
     */
    PETScSNESFunctionGOWrapper(std::string object_name,
                               SNES petsc_snes,
                               PetscErrorCode (*petsc_snes_form_func)(SNES, Vec, Vec, void*),
                               void* petsc_snes_func_ctx);

    /*!
     * \brief Destructor.
     */
    ~PETScSNESFunctionGOWrapper();

    /*!
     * \name Functions to access the underlying PETSc function.
     */
    //\{

    /*!
     * \brief Get the PETSc SNES object "wrapped" by this object.
     */
    const SNES& getPETScSNES() const;

    /*!
     * \brief Get the function pointer to the PETSc SNES function "wrapped" by
     * this object.
     */
    PetscErrorCode (*getPETScSNESFormFunction())(SNES, Vec, Vec, void*);

    /*!
     * \brief Get the PETSc SNES function context object "wrapped" by this
     * object.
     */
    void* getPETScSNESFunctionContext() const;

    //\}

    /*!
     * \name General operator functionality.
     */
    //\{

    /*!
     * \brief Compute y=F[x].
     *
     * Before calling this function, the form of the vectors x and y should be
     * set properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in these vectors should be allocated.
     * The user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must have enough
     *   ghost cells for computation of F[x]).
     *
     * In general, the vectors x and y \em cannot be the same.
     *
     * \note The operator NEED NOT be initialized prior to calling apply.
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y output: y=F[x]
     */
    void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
               SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    /*!
     * \brief Compute hierarchy dependent data required for computing y=F[x].
     *
     * The vector arguments for apply() need not match those for
     * initializeOperatorState().  However, there must be a certain degree of
     * similarity, including
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
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out) override;

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeOperatorState().
     *
     * Remove all hierarchy dependent data set by initializeOperatorState().  It
     * is safe to call deallocateOperatorState() when the operator state is
     * already deallocated.
     *
     * \see initializeOperatorState
     */
    void deallocateOperatorState() override;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScSNESFunctionGOWrapper() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScSNESFunctionGOWrapper(const PETScSNESFunctionGOWrapper& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScSNESFunctionGOWrapper& operator=(const PETScSNESFunctionGOWrapper& that) = delete;

    const SNES d_petsc_snes;
    PetscErrorCode (*const d_petsc_snes_form_func)(SNES, Vec, Vec, void*);
    void* const d_petsc_snes_func_ctx;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_x, d_y;
    Vec d_petsc_x = nullptr, d_petsc_y = nullptr;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_PETScSNESFunctionGOWrapper
