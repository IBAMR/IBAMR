// Filename: PETScMatLOWrapper.h
// Created on 26 Dec 2003 by Boyce Griffith
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

#ifndef included_PETScMatLOWrapper
#define included_PETScMatLOWrapper

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "IntVector.h"
#include "SAMRAIVectorReal.h"
#include "ibtk/LinearOperator.h"
#include "petscmat.h"
#include "petscvec.h"
#include "tbox/Pointer.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PETScMatLOWrapper provides a LinearOperator interface for a <A
 * HREF="http://www.mcs.anl.gov/petsc">PETSc</A> Mat object.
 */
class PETScMatLOWrapper : public LinearOperator
{
public:
    /*!
     * \brief Constructor.
     *
     * Construct a linear operator wrapper corresponding to the provided PETSc
     * Mat object.
     *
     * \param object_name name of the linear operator
     * \param petsc_mat PETSc Mat object
     */
    PETScMatLOWrapper(const std::string& object_name, const Mat& petsc_mat);

    /*!
     * \brief Destructor.
     */
    ~PETScMatLOWrapper();

    /*!
     * \name Functions to access the underlying PETSc object.
     */
    //\{

    /*!
     * \brief Get the PETSc Mat object "wrapped" by this object.
     */
    const Mat& getPETScMat() const;

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
    void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x, SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y);

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
                  SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& z);

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
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out);

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
    void deallocateOperatorState();

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScMatLOWrapper();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScMatLOWrapper(const PETScMatLOWrapper& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScMatLOWrapper& operator=(const PETScMatLOWrapper& that);

    const Mat d_petsc_mat;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_x, d_y, d_z;
    Vec d_petsc_x, d_petsc_y, d_petsc_z;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PETScMatLOWrapper
