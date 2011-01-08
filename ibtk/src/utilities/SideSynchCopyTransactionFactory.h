// Filename: SideSynchCopyTransactionFactory.h
// Created on 17 Dec 2009 by Boyce Griffith
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

#ifndef included_SideSynchCopyTransactionFactory
#define included_SideSynchCopyTransactionFactory

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <PatchLevel.h>
#include <RefineClasses.h>
#include <RefineTransactionFactory.h>
#include <tbox/Arena.h>
#include <tbox/Pointer.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class SideSynchCopyTransactionFactory is a concrete subclass of the
 * SAMRAI::xfer::RefineTransactionFactory base class that allocates side
 * synchronization transaction objects for a SAMRAI::xfer::RefineSchedule
 * object.
 *
 * \see SideSynchCopyTransaction
 * \see SAMRAI::xfer::RefineTransactionFactory
 */
class SideSynchCopyTransactionFactory
    : public SAMRAI::xfer::RefineTransactionFactory<NDIM>
{
public:
    /*!
     * \brief Default constructor.
     */
    SideSynchCopyTransactionFactory();

    /*!
     * \brief Virtual destructor for base class.
     */
    virtual
    ~SideSynchCopyTransactionFactory();

    /*!
     * \brief Set the array of SAMRAI::xfer::RefineClass::Data items used by the
     * transactions.
     */
    void
    setRefineItems(
        const SAMRAI::xfer::RefineClasses<NDIM>::Data** refine_items,
        int num_refine_items);

    /*!
     * \brief Clear the array of SAMRAI::xfer::RefineClass<NDIM>::Data items used by the
     * transactions.
     */
    void
    unsetRefineItems();

    /*!
     * \brief Allocate an SideSynchCopyTransaction object.
     *
     * \param dst_level               SAMRAI::tbox::Pointer to destination patch level.
     * \param src_level               SAMRAI::tbox::Pointer to source patch level.
     * \param overlap                 SAMRAI::tbox::Pointer to overlap region between patches.
     * \param dst_patch_id            Integer index of destination patch in destination patch level.
     * \param src_patch_id            Integer index of source patch in source patch level.
     * \param ritem_id                Integer index of SAMRAI::xfer::RefineClass<NDIM>::Data item associated with transaction.
     * \param box                     Optional const reference to box defining region of refine transaction.  (Default is an empty box.)
     * \param use_time_interpolation  Optional boolean flag indicating whether the refine transaction involves time interpolation.  (Default is false.)
     * \param pool                    Optional pointer to memory pool from which the refine transaction may be allocated.  (Default is NULL.)
     */
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Transaction>
    allocate(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > dst_level,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > src_level,
        SAMRAI::tbox::Pointer<SAMRAI::hier::BoxOverlap<NDIM> > overlap,
        int dst_patch_id,
        int src_patch_id,
        int ritem_id,
        const SAMRAI::hier::Box<NDIM>& box=SAMRAI::hier::Box<NDIM>(),
        bool use_time_interpolation=false,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena> pool=(SAMRAI::tbox::Arena*)NULL) const;

private:

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SideSynchCopyTransactionFactory(
        const SideSynchCopyTransactionFactory& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    void
    operator=(
        const SideSynchCopyTransactionFactory& that);

    const SAMRAI::xfer::RefineClasses<NDIM>::Data** d_refine_items;
    int d_number_refine_items;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/SideSynchCopyTransactionFactory.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_SideSynchCopyTransactionFactory
