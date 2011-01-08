// Filename: VecCellDataFactory.h
// Created on 09 Apr 2010 by Boyce Griffith
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

#ifndef included_VecCellDataFactory
#define included_VecCellDataFactory

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <PatchDataFactory.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * Class VecCellDataFactory is a factory class used to allocate new instances of
 * VecCellData objects.  It is a subclass of the patch data factory class and
 * cell data is a subclass of patch data.  Both the factory and data classes are
 * templated on the type of the contained object (e.g., double or int).
 *
 * \see VecCellData
 * \see VecCellVariable
 * \see SAMRAI::pdat::PatchDataFactory
 */
template<class TYPE>
class VecCellDataFactory
    : public SAMRAI::hier::PatchDataFactory<NDIM>
{
public:
    /*!
     * The default constructor for the cell data factory class.  The ghost cell
     * width and depth (number of components) arguments give the defaults for
     * all cell data objects created with this factory.
     */
    VecCellDataFactory(
        int depth,
        const SAMRAI::hier::IntVector<NDIM>& ghosts);

    /*!
     * Virtual destructor for the cell data factory class.
     */
    virtual
    ~VecCellDataFactory<TYPE>();

    /*!
     * \brief Abstract virtual function to clone a patch data factory.
     *
     * This will return a new instantiation of the abstract factory with the
     * same properties.  The properties of the cloned factory can then be
     * changed without modifying the original.
     *
     * \param ghosts default ghost cell width for concrete classes created from
     *        the factory.
     */
    virtual SAMRAI::tbox::Pointer<SAMRAI::hier::PatchDataFactory<NDIM> >
    cloneFactory(
        const SAMRAI::hier::IntVector<NDIM>& ghosts);

    /*!
     * Virtual factory function to allocate a concrete cell data object.  The
     * default information about the object (e.g., ghost cell width) is taken
     * from the factory.  If no memory pool is provided, then the allocation
     * routine assumes some default memory pool.
     */
    virtual SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >
    allocate(
        const SAMRAI::hier::Box<NDIM>& box,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena> pool=SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena>(NULL)) const;

    /*!
     * Virtual factory function to allocate a concrete cell data object.  Same
     * as above function, except passes in a patch instead of a box.
     */
    virtual SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >
    allocate(
        const SAMRAI::hier::Patch<NDIM>& patch,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena> pool=SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena>(NULL)) const;

    /*!
     * Allocate the box geometry object associated with the patch data.  This
     * information will be used in the computation of intersections and data
     * dependencies between objects.
     */
    virtual SAMRAI::tbox::Pointer<SAMRAI::hier::BoxGeometry<NDIM> >
    getBoxGeometry(
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * Get the default depth (number of components).  This is the default depth
     * that will be used in the instantiation of cell data objects.
     */
    int
    getDefaultDepth() const;

    /*!
     * Set the default depth (number of components).  This is the default depth
     * that will be used in the instantiation of cell data objects.
     */
    void
    setDefaultDepth(
        const int depth);

    /*!
     * Calculate the amount of memory needed to store the cell data object,
     * including object data and dynamically allocated data.
     */
    virtual size_t
    getSizeOfMemory(
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * Return a boolean true value indicating that the cell data quantities will
     * always be treated as though fine values represent them on coarse-fine
     * interfaces.  See the VecCellVariable class header file for more
     * information.
     */
    virtual bool
    fineBoundaryRepresentsVariable() const;

    /*!
     * Return false since the cell data index space matches the cell-centered
     * index space for AMR patches.  Thus, cell data does not live on patch
     * borders.
     */
    virtual bool
    dataLivesOnPatchBorder() const;

    /*!
     * Return whether it is valid to copy this VecCellDataFactory to the
     * supplied destination patch data factory.  It will return true if dst_pdf
     * is a VecCellDataFactory, false otherwise.
     */
    virtual bool
    validCopyTo(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchDataFactory<NDIM> >& dst_pdf) const;

    /*!
     * Return pointer to a multiblock data translator
     */
    virtual SAMRAI::hier::MultiblockDataTranslator<NDIM>*
    getMultiblockDataTranslator();

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    VecCellDataFactory();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VecCellDataFactory(
        const VecCellDataFactory& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VecCellDataFactory&
    operator=(
        const VecCellDataFactory& that);

    int d_depth;
    SAMRAI::hier::MultiblockDataTranslator<NDIM>* d_mb_trans;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibtk/VecCellDataFactory.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_VecCellDataFactory
