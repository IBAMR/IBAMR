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

#include <unordered_set>
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
    template <typename INPUT_TYPE = TYPE>
    SAMRAIScopedVectorDuplicate(const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, INPUT_TYPE> >& vector,
                                const std::string& name = "");

    /*!
     * Constructor. Sets up a vector equivalent to @p vector but does not copy values.
     */
    template <typename INPUT_TYPE = TYPE>
    SAMRAIScopedVectorDuplicate(const SAMRAI::solv::SAMRAIVectorReal<NDIM, INPUT_TYPE>& vector,
                                const std::string& name = "");

    /*!
     * Conversion operator to a SAMRAI vector.
     */
    operator SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE>&();

    /*!
     * Conversion operator to non-owning pointer to a SAMRAI vector.
     *
     * @note This operator should only be used with APIs which expect a pointer
     * to a vector. Since this pointer is non-owning it is not an implementation
     * of reference-counting to this object.
     */
    operator SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> >();

    /*!
     * Get the components of the vectors.
     */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> > > getComponentVectors() const;

    /*!
     * Transform from a different precision representation.
     */
    template <typename SRC_TYPE>
    void transformFromVector(const SAMRAIScopedVectorDuplicate<SRC_TYPE>& src_vec);

    /*!
     * Transform from a different precision representation.
     */
    template <typename SRC_TYPE>
    void transformFromVector(const SAMRAI::solv::SAMRAIVectorReal<NDIM, SRC_TYPE>& src_vec);

    /*!
     * Transform to a different precision representation.
     */
    template <typename DST_TYPE>
    void transformToVector(SAMRAIScopedVectorDuplicate<DST_TYPE>& dst_vec) const;

    /*!
     * Transform to a different precision representation.
     */
    template <typename DST_TYPE>
    void transformToVector(SAMRAI::solv::SAMRAIVectorReal<NDIM, DST_TYPE>& dst_vec) const;

    /*!
     * Destructor. Removes the cloned patch index and deallocates data.
     */
    ~SAMRAIScopedVectorDuplicate();

protected:
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> > d_vector;

    template <class SRC_TYPE, class DST_TYPE>
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > cloneVariable(SAMRAI::hier::Variable<NDIM>* input_variable,
                                                                       int input_descriptor_index) const;

    template <class DATA_TYPE>
    SAMRAI::hier::IntVector<NDIM> getGhostCellWidth(int data_id) const;

    void allocatePatchData(int data_id) const;

    void freePatchData(int data_id) const;

    template <class SRC_TYPE, class DST_TYPE>
    void transformPatchData(int src_data_id, int dst_data_id) const;

    std::unordered_set<int> d_managed_cv_data_ids;
};
} // namespace IBTK

#include "ibtk/private/SAMRAIScopedVectorDuplicate-inl.h"

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_SAMRAIScopedVectorDuplicate
