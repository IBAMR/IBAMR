// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_LIndexSetDataFactory
#define included_IBTK_LIndexSetDataFactory

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LSetDataFactory.h"

#include "Box.h"
#include "IntVector.h"
#include "tbox/Arena.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class PatchData;
template <int DIM>
class PatchDataFactory;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LIndexSetPatchDataFactory provides a
 * SAMRAI::hier::PatchDataFactory class corresponding to patch data of type
 * LIndexSetData.
 */
template <class T>
class LIndexSetDataFactory : public LSetDataFactory<T>
{
public:
    /*!
     * The default constructor for the LIndexSetDataFactory class.  The ghost
     * cell width argument gives the default width for all data objects created
     * with this factory.
     */
    LIndexSetDataFactory(SAMRAI::hier::IntVector<NDIM> ghosts);

    /*!
     * Virtual destructor for the data factory class.
     */
    virtual ~LIndexSetDataFactory() = default;

    /*!
     * Virtual factory function to allocate a concrete data object.  The default
     * information about the object (e.g., ghost cell width) is taken from the
     * factory.  If no memory pool is provided, the allocation routine assumes
     * some default memory pool.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >
    allocate(const SAMRAI::hier::Box<NDIM>& box, SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena> pool = NULL) const override;

    /*!
     * Virtual factory function to allocate a concrete data object.  The default
     * information about the object (e.g., ghost cell width) is taken from the
     * factory.  If no memory pool is provided, the allocation routine assumes
     * some default memory pool.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >
    allocate(const SAMRAI::hier::Patch<NDIM>& patch,
             SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena> pool = NULL) const override;

    /*!
     * Calculate the amount of memory needed to store the data object, including
     * object data but not dynamically allocated data.
     */
    size_t getSizeOfMemory(const SAMRAI::hier::Box<NDIM>& box) const override;

    /*!
     * Virtual function to clone the data factory.  This will return a new
     * instantiation of the factory with the same properties (e.g., same type).
     * The properties of the cloned factory can then be changed without
     * modifying the original.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchDataFactory<NDIM> >
    cloneFactory(const SAMRAI::hier::IntVector<NDIM>& ghosts) override;

    /*!
     * Return whether it is valid to copy this LIndexSetDataFactory to the
     * supplied destination patch data factory. It will return true if dst_pdf
     * is a LIndexSetDataFactory, false otherwise.
     */
    bool validCopyTo(const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchDataFactory<NDIM> >& dst_pdf) const override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LIndexSetDataFactory() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LIndexSetDataFactory(const LIndexSetDataFactory<T>& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LIndexSetDataFactory& operator=(const LIndexSetDataFactory<T>& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LIndexSetDataFactory
