// Filename: VecCellData.h
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

#ifndef included_VecCellData
#define included_VecCellData

/////////////////////////////// INCLUDES /////////////////////////////////////

// BLITZ++ INCLUDES
#include <blitz/array.h>

// SAMRAI INCLUDES
#include <CellData.h>
#include <CellIndex.h>
#include <CellIterator.h>
#include <PatchData.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class VecCellData<NDIM> provides an implementation for vector-valued
 * data defined at cell centers on AMR patches.  It is derived from the
 * SAMRAI::hier::PatchData interface common to all SAMRAI patch data types.
 * Given a CELL-centered AMR index space box, a cell data object represents data
 * of some template TYPE and depth at the centers of the cells in the box.
 * Here, depth indicates the number of data values at each cell index location.
 * The SAMRAI::pdat::CellGeometry class provides the translation between the
 * standard SAMRAI cell-centered AMR index space and cell-centered data.
 *
 * A cell-centerd data array is stored in (d,i,...,k) order, where i,...,k are
 * spatial indices and d indicates the depth at that location.  Memory
 * allocation is in column-major ordering (e.g., Fortran style) so that the
 * leftmost index runs fastest in memory.  For example, a three-dimensional cell
 * data object defined over a box [l0:u0,l1:u1,l2:u2] holds a data array
 * dimensioned as \verbatim

 [ 0 : depth-1,
 l0 : u0 ,
 l1 : u1 ,
 l2 : u2 ]

 * \endverbatim
 * Other spatial dimensions are represented similarly.
 *
 * The data type TYPE must define a default constructor (i.e., taking no
 * arguments) and also the copy assignment operator.
 *
 * \see pdat::ArrayData
 * \see hier::PatchData
 * \see pdat::VecCellDataFactory
 * \see pdat::CellIndex
 * \see pdat::CellIterator
 * \see pdat::CellGeometry
 */

template<class TYPE>
class VecCellData
    : public SAMRAI::hier::PatchData<NDIM>
{
public:
    /*!
     * \brief Calculate the amount of memory needed to represent cell- centered
     * data over a CELL-centered AMR index space box.
     *
     * This function assumes that the amount of memory needed for TYPE is
     * sizeof(TYPE).  If this is not the case, then a specialized function must
     * be defined.
     *
     * \param box const Box reference describing the interior of the standard
     *        CELL-centered index box over which the cell data object will be
     *        created.
     * \param depth gives the number of components for each spatial location in
     *        the array.
     * \param ghosts const IntVector reference indicating the width of the ghost
     *        cell region around the box over which the node data will be
     *        allocated.
     */
    static size_t
    getSizeOfData(
        const SAMRAI::hier::Box<NDIM>& box,
        const unsigned int depth,
        const SAMRAI::hier::IntVector<NDIM>& ghosts);

    /*!
     * \brief The constructor for a cell data object.
     *
     * \param box const Box reference describing the interior of the standard
     *        CELL-centered index box over which the cell data object will be
     *        created.
     * \param depth gives the number of components for each spatial location in
     *        the array.
     * \param ghosts const IntVector reference indicating the width of the ghost
     *        cell region around the box over which the node data will be
     *        allocated.
     */
    VecCellData(
        const SAMRAI::hier::Box<NDIM>& box,
        const unsigned int depth,
        const SAMRAI::hier::IntVector<NDIM>& ghosts);

    /*!
     * \brief The virtual destructor for a cell data object.
     */
    virtual
    ~VecCellData<TYPE>();

    /*!
     * \brief Return the depth (e.g., the number of components in each spatial
     * location) of the array.
     */
    unsigned int
    getDepth() const;

    /*!
     * \brief Get a pointer to the beginning of the cell centered array.
     */
    TYPE*
    getPointer();

    /*!
     * \brief Get a const pointer to the beginning of the cell centered array.
     */
    const TYPE*
    getPointer() const;

    /*!
     * \brief Return reference to cell data entry corresponding to a given cell
     * index and depth.
     */
    blitz::Array<TYPE,1>
    operator()(
        const SAMRAI::pdat::CellIndex<NDIM>& i);

    /*!
     * \brief Return reference to cell data entry corresponding to a given cell
     * index and depth.
     */
    TYPE&
    operator()(
        const SAMRAI::pdat::CellIndex<NDIM>& i,
        const unsigned int depth);

    /*!
     * \brief Return a const reference to cell data entry corresponding to a
     * given cell index and depth.
     */
    const TYPE&
    operator()(
        const SAMRAI::pdat::CellIndex<NDIM>& i,
        const unsigned int depth) const;

    /*!
     * \brief Return a reference to the Blitz++ array object for the cell
     * centered data object.
     */
    blitz::Array<TYPE,NDIM+1>&
    getArray();

    /*!
     * \brief Return a const reference to the Blitz++ array object for the cell
     * centered data object.
     */
    const blitz::Array<TYPE,NDIM+1>&
    getArray() const;

    /*!
     * \brief A fast copy from source to destination (i.e., this) patch data
     * object.
     *
     * Data is copied where there is overlap in the underlying index space.  The
     * copy is performed on the interior plus the ghost cell width (for both the
     * source and destination).  Currently, source data must be CellData or
     * VecCellData of the TYPE.  If not, then an unrecoverable error results.
     */
    void
    copy(
        const SAMRAI::hier::PatchData<NDIM>& src);

    /*!
     * \brief A fast copy from source (i.e., this) to destination patch data
     * object.
     *
     * Data is copied where there is overlap in the underlying index space.  The
     * copy is performed on the interior plus the ghost cell width (for both the
     * source and destination).  Currently, destination data must be CellData or
     * VecCellData of the same TYPE.  If not, then an unrecoverable error
     * results.
     */
    void
    copy2(
        SAMRAI::hier::PatchData<NDIM>& dst) const;

    /*!
     * \brief Copy data from source to destination (i.e., this) patch data
     * object on the given overlap.
     *
     * Source data must be VecCellData of the same TYPE and the overlap must be
     * a CellOverlap.  If not, then an unrecoverable error results.
     */
    void
    copy(
        const SAMRAI::hier::PatchData<NDIM>& src,
        const SAMRAI::hier::BoxOverlap<NDIM>& overlap);

    /*!
     * \brief Copy data from source (i.e., this) to destination patch data
     * object on the given overlap.
     *
     * Destination data must be VecCellData of the same TYPE and the overlap
     * must be a CellOverlap.  If not, then an unrecoverable error results.
     */
    void
    copy2(
        SAMRAI::hier::PatchData<NDIM>& dst,
        const SAMRAI::hier::BoxOverlap<NDIM>& overlap) const;

    /*!
     * \brief Copy data from source to destination (i.e., this) patch data
     * object on the given CELL-centered AMR index box.
     */
    void
    copyOnBox(
        const VecCellData<TYPE>& src,
        const SAMRAI::hier::Box<NDIM>& box);

    /*!
     * \brief Copy data from source to destination (i.e., this) patch data
     * object on the given CELL-centered AMR index box.
     */
    void
    copyOnBox(
        const VecCellData<TYPE>& src,
        const SAMRAI::hier::Box<NDIM>& box,
        const SAMRAI::hier::IntVector<NDIM>& src_offset);

    /*!
     * \brief Copy data from source to destination (i.e., this) patch data
     * object on the given CELL-centered AMR index box.
     */
    void
    copyOnBox(
        const SAMRAI::pdat::CellData<NDIM,TYPE>& src,
        const SAMRAI::hier::Box<NDIM>& box);

    /*!
     * \brief Copy data from source to destination (i.e., this) patch data
     * object on the given CELL-centered AMR index box.
     */
    void
    copy2OnBox(
        SAMRAI::pdat::CellData<NDIM,TYPE>& dst,
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * \brief Return true if the patch data object can estimate the stream size
     * required to fit its data using only index space information (i.e., a
     * box).
     *
     * This routine is defined for the standard types (bool, char, double,
     * float, int, and dcomplex).
     */
    bool
    canEstimateStreamSizeFromBox() const;

    /*!
     * \brief Return the number of bytes needed to stream the data in this patch
     * data object lying in the specified box overlap region.
     *
     * This routine is defined for the standard types (bool, char, double,
     * float, int, and dcomplex).
     */
    int
    getDataStreamSize(
        const SAMRAI::hier::BoxOverlap<NDIM>& overlap) const;

    /*!
     * \brief Pack data to stream from this patch data object over the specified
     * box overlap region.  The overlap must be a CellOverlap.
     */
    void
    packStream(
        SAMRAI::tbox::AbstractStream& stream,
        const SAMRAI::hier::BoxOverlap<NDIM>& overlap) const;

    /*!
     * \brief Unpack data from stream into this patch data object over the
     * specified box overlap region.  The overlap must be a CellOverlap.
     */
    void
    unpackStream(
        SAMRAI::tbox::AbstractStream& stream,
        const SAMRAI::hier::BoxOverlap<NDIM>& overlap);

    /*!
     * \brief Fill all values at depth d with the value t.
     */
    void
    fill(
        const TYPE& t,
        const unsigned int d=0);

    /*!
     * \brief Fill all values at depth d within the box with the value t.
     */
    void
    fill(
        const TYPE& t,
        const SAMRAI::hier::Box<NDIM>& box,
        const unsigned int d=0);

    /*!
     * \brief Fill all depth components with value t.
     */
    void
    fillAll(
        const TYPE& t);

    /*!
     * \brief Fill all depth components within the box with value t.
     */
    void
    fillAll(
        const TYPE& t,
        const SAMRAI::hier::Box<NDIM>& box);

    /*!
     * \brief Print all cell data values residing in the specified box.  If the
     * depth of the array is greater than one, all depths are printed.
     *
     * \param box const reference to box over whioch to print data.  Note that
     *        box is assumed to reside in standard cell-centered index space and
     *        will be converted to cell index space.
     * \param os reference to output stream.
     * \param prec integer precision for printing floating point numbers (i.e.,
     *        TYPE = float, double, or dcomplex).  The default is 12 decimal
     *        places for double and complex floating point numbers, and the
     *        default is 6 decimal places floats.  For other types, this value
     *        is ignored.
     */
    void
    print(
        const SAMRAI::hier::Box<NDIM>& box,
        std::ostream& os=SAMRAI::tbox::plog,
        const int prec=12) const;

    /*!
     * \brief Print all cell data values at the given array depth in the
     * specified box.
     *
     * \param box const reference to box over whioch to print data.  Note that
     *        box is assumed to reside in standard cell-centered index space and
     *        will be converted to cell index space.
     *
     * \param depth integer depth component, must satisfy 0 <= depth < actual
     *        depth of data array.
     * \param os reference to output stream.
     * \param prec integer precision for printing floating point numbers (i.e.,
     *        TYPE = float, double, or dcomplex).  The default is 12 decimal
     *        places for double and complex floating point numbers, and the
     *        default is 6 decimal places floats.  For other types, this value
     *        is ignored.
     */
    void
    print(
        const SAMRAI::hier::Box<NDIM>& box,
        const unsigned int depth,
        std::ostream& os=SAMRAI::tbox::plog,
        const int prec=12) const;

    /*!
     * Check that class version and restart file version are equal.  If so, read
     * data members from the database.
     */
    void
    getSpecializedFromDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database);

    /*!
     * Write out the class version number and other data members to the
     * database.
     */
    void
    putSpecializedToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database);

    /*!
     * The cell iterator iterates over the elements of a cell centered box
     * geometry.  This typedef is a convenience for using the CellIterator<NDIM>
     * class.
     */
    typedef SAMRAI::pdat::CellIterator<NDIM> Iterator;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    VecCellData();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VecCellData(
        const VecCellData& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VecCellData&
    operator=(
        const VecCellData& that);

    /*!
     * \brief Get data from database.
     */
    void
    getArrayFromDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database);

    /*!
     * \brief Put data to database.
     */
    void
    putArrayToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database);

    unsigned int d_depth;
    blitz::Array<TYPE,NDIM+1> d_data;

};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibtk/VecCellData.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_VecCellData
