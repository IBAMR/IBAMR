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

#ifndef included_IBTK_PETScSAMRAIVectorReal
#define included_IBTK_PETScSAMRAIVectorReal

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "IntVector.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Pointer.h"

#include "petscsys.h"
#include "petscvec.h"

#include <mpi.h>
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
     * Get a pointer to the SAMRAI vector object associated with the given
     * PETSc vector object.
     *
     * \note The SAMRAI vector must be restored by calling restoreSAMRAIVector().
     */
    static void getSAMRAIVector(Vec petsc_vec,
                                SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> >* samrai_vec);

    /*!
     * Restore the SAMRAI vector object associated with the given PETSc vector object.
     */
    static void
    restoreSAMRAIVector(Vec petsc_vec,
                        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> >* samrai_vec);

    /*!
     * Get a pointer to the SAMRAI vector object associated with the given
     * PETSc vector object.  This vector must be treated as read only.
     *
     * \note The SAMRAI vector must be restored by calling restoreSAMRAIVectorRead().
     */
    static void
    getSAMRAIVectorRead(Vec petsc_vec,
                        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> >* samrai_vec);

    /*!
     * Restore the SAMRAI vector object associated with the given PETSc vector object.
     */
    static void
    restoreSAMRAIVectorRead(Vec petsc_vec,
                            SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> >* samrai_vec);

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
    PETScSAMRAIVectorReal() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScSAMRAIVectorReal(const PETScSAMRAIVectorReal& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScSAMRAIVectorReal& operator=(const PETScSAMRAIVectorReal& that) = delete;

    static PetscErrorCode VecDuplicate_SAMRAI(Vec v, Vec* newv);

    static PetscErrorCode VecDestroy_SAMRAI(Vec v);

    static PetscErrorCode VecDuplicateVecs_SAMRAI(Vec v, PetscInt m, Vec* V[]);

    static PetscErrorCode VecDestroyVecs_SAMRAI(PetscInt m, Vec vv[]);

    static PetscErrorCode VecDot_SAMRAI(Vec x, Vec y, PetscScalar* val);

    static PetscErrorCode VecMDot_SAMRAI(Vec x, PetscInt nv, const Vec* y, PetscScalar* val);

    static PetscErrorCode VecNorm_SAMRAI(Vec x, NormType type, PetscScalar* val);

    static PetscErrorCode VecTDot_SAMRAI(Vec x, Vec y, PetscScalar* val);

    static PetscErrorCode VecMTDot_SAMRAI(Vec x, PetscInt nv, const Vec* y, PetscScalar* val);

    static PetscErrorCode VecScale_SAMRAI(Vec x, PetscScalar alpha);

    static PetscErrorCode VecCopy_SAMRAI(Vec x, Vec y);

    static PetscErrorCode VecSet_SAMRAI(Vec x, PetscScalar alpha);

    static PetscErrorCode VecSwap_SAMRAI(Vec x, Vec y);

    static PetscErrorCode VecAXPY_SAMRAI(Vec y, PetscScalar alpha, Vec x);

    static PetscErrorCode VecAXPBY_SAMRAI(Vec y, PetscScalar alpha, PetscScalar beta, Vec x);

    static PetscErrorCode VecMAXPY_SAMRAI(Vec y, PetscInt nv, const PetscScalar* alpha, Vec* x);

    static PetscErrorCode VecAYPX_SAMRAI(Vec y, const PetscScalar alpha, Vec x);

    static PetscErrorCode VecWAXPY_SAMRAI(Vec w, PetscScalar alpha, Vec x, Vec y);

    static PetscErrorCode
    VecAXPBYPCZ_SAMRAI(Vec z, PetscScalar alpha, PetscScalar beta, PetscScalar gamma, Vec x, Vec y);

    static PetscErrorCode VecPointwiseMult_SAMRAI(Vec w, Vec x, Vec y);

    static PetscErrorCode VecPointwiseDivide_SAMRAI(Vec w, Vec x, Vec y);

    static PetscErrorCode VecGetSize_SAMRAI(Vec v, PetscInt* size);

    static PetscErrorCode VecGetLocalSize_SAMRAI(Vec v, PetscInt* size);

    static PetscErrorCode VecMax_SAMRAI(Vec x, PetscInt* p, PetscScalar* val);

    static PetscErrorCode VecMin_SAMRAI(Vec x, PetscInt* p, PetscScalar* val);

    static PetscErrorCode VecSetRandom_SAMRAI(Vec x, PetscRandom rctx);

    static PetscErrorCode VecDot_local_SAMRAI(Vec x, Vec y, PetscScalar* val);

    static PetscErrorCode VecTDot_local_SAMRAI(Vec x, Vec y, PetscScalar* val);

    static PetscErrorCode VecNorm_local_SAMRAI(Vec x, NormType type, PetscScalar* val);

    static PetscErrorCode VecMDot_local_SAMRAI(Vec x, PetscInt nv, const Vec* y, PetscScalar* val);

    static PetscErrorCode VecMTDot_local_SAMRAI(Vec x, PetscInt nv, const Vec* y, PetscScalar* val);

    static PetscErrorCode VecMaxPointwiseDivide_SAMRAI(Vec x, Vec y, PetscScalar* max);

    static PetscErrorCode VecDotNorm2_SAMRAI(Vec s, Vec t, PetscScalar* dp, PetscScalar* nm);

    /*
     * Vector data is maintained in the SAMRAI vector structure.
     */
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> > d_samrai_vector;

    /*
     * PETSc vector object corresponding to this PETScAbstractVectorReal object.
     */
    Vec d_petsc_vector;
    bool d_vector_created_via_duplicate, d_vector_checked_out_read_write = false, d_vector_checked_out_read = false;
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/PETScSAMRAIVectorReal-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_PETScSAMRAIVectorReal
