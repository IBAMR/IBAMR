// Filename: PETScSNESJacobianJOWrapper.h
// Created on 27 Dec 2003 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#ifndef included_IBTK_PETScSNESJacobianJOWrapper
#define included_IBTK_PETScSNESJacobianJOWrapper

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "IntVector.h"
#include "SAMRAIVectorReal.h"
#include "ibtk/JacobianOperator.h"
#include "petscmat.h"
#include "petscsnes.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Pointer.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PETScSNESJacobianJOWrapper provides a JacobianOperator interface
 * for a <A HREF="http://www.mcs.anl.gov/petsc">PETSc</A> SNES Jacobian.
 */
class PETScSNESJacobianJOWrapper : public JacobianOperator
{
public:
    /*!
     * \brief Constructor.
     *
     * Construct a Jacobian operator wrapper corresponding to the provided PETSc
     * SNES Jacobian object.
     */
    PETScSNESJacobianJOWrapper(const std::string& object_name,
                               const SNES& petsc_snes,
                               PetscErrorCode (*petsc_snes_form_jac)(SNES, Vec, Mat, Mat, void*),
                               void* petsc_snes_jac_ctx);

    /*!
     * \brief Destructor.
     */
    ~PETScSNESJacobianJOWrapper() override;

    /*!
     * \name Functions to access the underlying PETSc function.
     */
    //\{

    /*!
     * \brief Get the PETSc SNES object "wrapped" by this object.
     */
    const SNES& getPETScSNES() const;

    /*!
     * \brief Get the function pointer to the PETSc SNES Jacobian "wrapped" by
     * this object.
     */
    PetscErrorCode (*getPETScSNESFormJacobian())(SNES, Vec, Mat, Mat, void*);

    /*!
     * \brief Get the PETSc SNES Jacobian context object "wrapped" by this
     * object.
     */
    void* getPETScSNESJacobianContext() const;

    //\}

    /*!
     * \name General Jacobian functionality.
     */
    //\{

    /*!
     * \brief Compute hierarchy dependent data required for evaluating F'[x].
     *
     * \param x value where the Jacobian is to be evaluated
     */
    void formJacobian(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x) override;

    /*!
     * \brief Return the vector where the Jacobian is evaluated.
     */
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > getBaseVector() const override;

    //\}

    /*!
     * \name Linear operator functionality.
     */
    //\{

    /*!
     * \brief Compute y=Ax.
     *
     * Before calling this function, the form of the vectors x and y should be
     * set properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in these vectors should be allocated.
     * The user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must have enough
     *   ghost cells for computation of Ax).
     *
     * In general, the vectors x and y \em cannot be the same.
     *
     * \note The operator NEED NOT be initialized prior to calling apply.
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y output: y=Ax
     */
    void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
               SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    /*!
     * \brief Compute z=Ax+y.
     *
     * Before calling this function, the form of the vectors x, y, and z should
     * be set properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in these vectors should be allocated.
     * The user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must have enough
     *   ghost cells for computation of Ax).
     *
     * In general, the vectors x and z \em cannot be the same.
     *
     * \note The operator NEED NOT be initialized prior to calling apply.
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y input
     * \param z output: z=Ax+y
     */
    void applyAdd(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                  SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y,
                  SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& z) override;

    /*!
     * \brief Compute hierarchy dependent data required for computing y=Ax and
     * z=Ax+y.
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
    PETScSNESJacobianJOWrapper() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScSNESJacobianJOWrapper(const PETScSNESJacobianJOWrapper& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScSNESJacobianJOWrapper& operator=(const PETScSNESJacobianJOWrapper& that) = delete;

    const SNES d_petsc_snes;
    Mat d_petsc_snes_jac;
    PetscErrorCode (*const d_petsc_snes_form_jac)(SNES, Vec, Mat, Mat, void*);
    void* const d_petsc_snes_jac_ctx;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_x, d_y, d_z;
    Vec d_petsc_x, d_petsc_y, d_petsc_z;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_PETScSNESJacobianJOWrapper
