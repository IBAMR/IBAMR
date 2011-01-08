// Filename: PETScVecOps.h
// Created on 23 Jan 2009 by Boyce Griffith
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

#ifndef included_PETScVecOps
#define included_PETScVecOps

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petscvec.h>

// SAMRAI INCLUDES
#include <tbox/Database.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <ostream>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PETScVecOps provides utility functions for <A
 * HREF="http://www-unix.mcs.anl.gov/petsc">PETSc</A> Vec objects.
 */
class PETScVecOps
{
public:
    /*!
     * \brief Sort modes for accumulating values.
     *
     * \note Default is: NO_SORT.
     */
    enum SortMode {NO_SORT=0, SORT_INCREASING_MAGNITUDE=1, SORT_DECREASING_MAGNITUDE=2};
    static SortMode s_sort_mode;

    /*!
     * \brief Floating point precision modes for accumulating values.
     *
     * \note Default is: DOUBLE.
     */
    enum PrecisionMode {DOUBLE=0, DOUBLE_DOUBLE=1, QUAD_DOUBLE=2};
    static PrecisionMode s_precision_mode;

    /*!
     * \brief Summation modes for accumulating values.
     *
     * \note Default is: RECURSIVE_SUMMATION.
     */
    enum SummationMode {RECURSIVE_SUMMATION=0, COMPENSATED_SUMMATION=1};
    static SummationMode s_summation_mode;

    /*!
     * \brief Set static configuration options from a user-supplied database.
     */
    static void
    setFromDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * \brief Output class configuration.
     */
    static void
    printClassData(
        std::ostream& os);

    /*!
     * \brief Accurate version of the PETSc function VecSetValues.
     */
    static PetscErrorCode
    VecSetValues(
        Vec xin,
        PetscInt ni,
        const PetscInt ix[],
        const PetscScalar y[],
        InsertMode addv);

    /*!
     * \brief Accurate version of the PETSc function VecSetValuesBlocked.
     */
    static PetscErrorCode
    VecSetValuesBlocked(
        Vec xin,
        PetscInt ni,
        const PetscInt ix[],
        const PetscScalar yin[],
        InsertMode addv);

    /*!
     * \brief Accurate version of the PETSc function VecAssemblyBegin.
     */
    static PetscErrorCode
    VecAssemblyBegin(
        Vec xin);

    /*!
     * \brief Accurate version of the PETSc function VecAssemblyEnd.
     */
    static PetscErrorCode
    VecAssemblyEnd(
        Vec xin);

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScVecOps();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScVecOps(
        const PETScVecOps& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScVecOps&
    operator=(
        const PETScVecOps& that);
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "PETScVecOps.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PETScVecOps
