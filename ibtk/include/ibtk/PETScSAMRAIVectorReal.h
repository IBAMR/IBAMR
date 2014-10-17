// Filename: PETScSAMRAIVectorReal.h
// Created on 10 Nov 2004 by Boyce Griffith
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

#ifndef included_PETScSAMRAIVectorReal
#define included_PETScSAMRAIVectorReal

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "IntVector.h"
#include "SAMRAIVectorReal.h"
#include "mpi.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Pointer.h"
// IWYU pragma: no_include "petscmath.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PETScSAMRAIVectorReal is a class for wrapping
 * SAMRAI::solv::SAMRAIVectorReal objects so that they may be used with the <A
 * HREF="http://www.mcs.anl.gov/petsc">PETSc</A> solver package.
 *
 * Class PETScSAMRAIVectorReal wraps a real-valued SAMRAI vector (see
 * SAMRAI::solv::SAMRAIVectorReal class) object so that it may be used with the
 * PETSc solver package.  A SAMRAI vector is defined as a collection of patch
 * data components and associated operations living on some subset of levels in
 * a structured AMR mesh hierarchy.
 *
 * Observe that there are only three public member functions in this class.
 * They are used to create and destroy PETSc vectors (i.e., \p Vec objects) and
 * to obtain the SAMRAI vector associated with the PETSc vector.  In particular,
 * note that the constructor and destructor of this class are protected members.
 * The construction and destruction of instances of this class may occur only
 * through the static member functions that create and destroy PETSc vector
 * objects.
 *
 * Finally, we remark that PETSc allows vectors with complex-valued entries.
 * This class and the class SAMRAI::solv::SAMRAIVectorReal assume real-values
 * vectors, i.e., data of type \p double or \p float.  The (currently
 * unimplemented) class PETScSAMRAIVectorComplex must be used for complex data.
 *
 * \see SAMRAI::solv::SAMRAIVectorReal
 */
class PETScSAMRAIVectorReal
{
public:
    /*!
     * Create and return a new PETSc vector object that wraps the SAMRAI vector
     * object, so that the SAMRAI vector may be manipulated by PETSc routines.
     * It is important to note that this function does not allocate storage for
     * the vector data.  Data must be allocated through the SAMRAI vector object
     * directly.
     *
     * \note Each call to createPETScVector() should be matched with a
     * corresponding call to destroyPETScVector().
     */
    static Vec createPETScVector(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> > samrai_vec,
                                 MPI_Comm comm = PETSC_COMM_WORLD);

    /*!
     * Destroy a given PETSc vector object.  It is important to note that this
     * function does not deallocate storage for the vector data.  Vector data
     * must be deallocated through the SAMRAI vector object.
     *
     * \note Each call to createPETScVector() should be matched with a
     * corresponding call to destroyPETScVector().
     */
    static void destroyPETScVector(Vec petsc_vec);

    /*!
     * Return pointer to the SAMRAI vector object associated with the given
     * PETSc vector object.
     */
    static SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> > getSAMRAIVector(Vec petsc_vec);

    /*!
     * Replace the SAMRAI vector object associated with the given PETSc vector
     * object.
     */
    static void
    replaceSAMRAIVector(Vec petsc_vec,
                        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> > samrai_vec);

protected:
    /*
     * Constructor for PETScSAMRAIVectorReal is protected so that an object of
     * this class cannot be constructed without calling the static member
     * function createPETScVector, which is used to construct a PETSc vector and
     * associated "wrapper" that allows the PETSc vector to manipulate the
     * SAMRAI::solv::SAMRAIVectorReal data.
     *
     * The boolean argument is used to control whether the SAMRAI vector is
     * destroyed when the associated PETSc vector is destroyed.  This should
     * happen if the PETSc vector is created within PETSc via a duplicate (i.e.,
     * clone) operation, but not otherwise.
     */
    PETScSAMRAIVectorReal(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> > samrai_vector,
                          bool vector_created_via_duplicate,
                          MPI_Comm comm);

    /*
     * Destructor for PETScSAMRAIVectorReal is protected so that an object of
     * this class cannot be destroyed without calling the static member function
     * destroyPETScVector, which is used to destroy the PETSc vector and
     * associated "wrapper."
     */
    ~PETScSAMRAIVectorReal();

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScSAMRAIVectorReal();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScSAMRAIVectorReal(const PETScSAMRAIVectorReal& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScSAMRAIVectorReal& operator=(const PETScSAMRAIVectorReal& that);

    static PetscErrorCode VecDuplicate_SAMRAI(Vec v, Vec* newv);

    static PetscErrorCode VecDestroy_SAMRAI(Vec v);

    /*
     * Vector data is maintained in the SAMRAI vector structure.
     */
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> > d_samrai_vector;

    /*
     * PETSc vector object corresponding to this PETScAbstractVectorReal object.
     */
    Vec d_petsc_vector;
    bool d_vector_created_via_duplicate;
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/PETScSAMRAIVectorReal-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PETScSAMRAIVectorReal
