// Filename: LNodeIndexDataFactory.h
// Created on 04 Jun 2007 by Boyce Griffith
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

#ifndef included_LNodeIndexDataFactory
#define included_LNodeIndexDataFactory

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/LNodeIndexSet.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CellGeometry.h>
#include <IndexDataFactory.h>
#include <PatchData.h>
#include <PatchDataFactory.h>
#include <tbox/Arena.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LNodeIndexDataFactory is a specialization of the templated class
 * SAMRAI::pdat::IndexDataFactory used to construct objects of type
 * LNodeIndexData.
 */
class LNodeIndexDataFactory
    : public SAMRAI::pdat::IndexDataFactory<NDIM,LNodeIndexSet,SAMRAI::pdat::CellGeometry<NDIM> >
{
public:
    /*!
     * The default constructor for the LNodeIndexDataFactory class.  The ghost
     * cell width argument gives the default width for all data objects created
     * with this factory.
     */
    LNodeIndexDataFactory(
        const SAMRAI::hier::IntVector<NDIM>& ghosts);

    /*!
     * Virtual destructor for the data factory class.
     */
    virtual
    ~LNodeIndexDataFactory();

    /*!
     * Virtual factory function to allocate a concrete data object.  The default
     * information about the object (e.g., ghost cell width) is taken from the
     * factory.  If no memory pool is provided, the allocation routine assumes
     * some default memory pool.
     */
    virtual SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >
    allocate(
        const SAMRAI::hier::Box<NDIM>& box,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena> pool=NULL) const;

    /*!
     * Virtual factory function to allocate a concrete data object.  The default
     * information about the object (e.g., ghost cell width) is taken from the
     * factory.  If no memory pool is provided, the allocation routine assumes
     * some default memory pool.
     */
    virtual SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >
    allocate(
        const SAMRAI::hier::Patch<NDIM>& patch,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena> pool=NULL) const;

    /*!
     * Calculate the amount of memory needed to store the data object, including
     * object data but not dynamically allocated data.
     */
    virtual size_t
    getSizeOfMemory(
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * Virtual function to clone the data factory.  This will return a new
     * instantiation of the factory with the same properties (e.g., same type).
     * The properties of the cloned factory can then be changed without
     * modifying the original.
     */
    virtual SAMRAI::tbox::Pointer<SAMRAI::hier::PatchDataFactory<NDIM> >
    cloneFactory(
        const SAMRAI::hier::IntVector<NDIM>& ghosts);

    /*!
     * Return whether it is valid to copy this LNodeIndexDataFactory to the
     * supplied destination patch data factory. It will return true if dst_pdf
     * is a LNodeIndexDataFactory, false otherwise.
     */
    virtual bool
    validCopyTo(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchDataFactory<NDIM> >& dst_pdf) const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LNodeIndexDataFactory();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LNodeIndexDataFactory(
        const LNodeIndexDataFactory& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LNodeIndexDataFactory&
    operator=(
        const LNodeIndexDataFactory& that);
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/LNodeIndexDataFactory.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeIndexDataFactory
