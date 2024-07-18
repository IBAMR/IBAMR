// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2023 by the IBAMR developers
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

#ifndef included_IBTK_SAMRAIScopedVectorDuplicate
#define included_IBTK_SAMRAIScopedVectorDuplicate

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <tbox/Pointer.h>

#include <SAMRAIVectorReal.h>

#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/**
 * Wrapper class around a SAMRAIVectorReal with RAII semantics (i.e., objects
 * will deallocate data and free its patch indices). Creates a new vector and
 * initializes its values to zero.
 *
 * @note The name of this class is analogous to the meaning of duplicate in
 * PETSc's VecDuplicate() function.
 */
template <typename TYPE>
class SAMRAIScopedVectorDuplicate
{
public:
    /*!
     * Constructor. Sets up a vector equivalent to @p vector but does not copy values.
     */
    SAMRAIScopedVectorDuplicate(const SAMRAIPointer<SAMRAI::solv::SAMRAIVectorRealNd<TYPE> >& vector,
                                const std::string& name = "");

    /*!
     * Constructor. Sets up a vector equivalent to @p vector but does not copy values.
     */
    SAMRAIScopedVectorDuplicate(const SAMRAI::solv::SAMRAIVectorRealNd<TYPE>& vector, const std::string& name = "");

    /*!
     * Conversion operator to a SAMRAI vector.
     */
    operator SAMRAI::solv::SAMRAIVectorRealNd<TYPE>&();

    /*!
     * Conversion operator to non-owning pointer to a SAMRAI vector.
     *
     * @note This operator should only be used with APIs which expect a pointer
     * to a vector. Since this pointer is non-owning it is not an implementation
     * of reference-counting to this object.
     */
    operator SAMRAIPointer<SAMRAI::solv::SAMRAIVectorRealNd<TYPE> >();

    /*!
     * Get the components of the vectors.
     */
    std::vector<SAMRAIPointer<SAMRAI::solv::SAMRAIVectorRealNd<TYPE> > > getComponentVectors() const;

    /*!
     * Destructor. Removes the cloned patch index and deallocates data.
     */
    ~SAMRAIScopedVectorDuplicate();

protected:
    SAMRAIPointer<SAMRAI::solv::SAMRAIVectorRealNd<TYPE> > d_vector;
};
} // namespace IBTK

#include "ibtk/private/SAMRAIScopedVectorDuplicate-inl.h"

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_SAMRAIScopedVectorDuplicate
