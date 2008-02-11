#ifndef included_LNodeIndexDataFactory2
#define included_LNodeIndexDataFactory2

// Filename: LNodeIndexDataFactory2.h
// Last modified: <07.Feb.2008 00:30:33 griffith@box221.cims.nyu.edu>
// Created on 04 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LNodeIndexSet.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CellDataFactory.h>
#include <PatchData.h>
#include <PatchDataFactory.h>
#include <tbox/Arena.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class LNodeIndexDataFactory2 is a specialization of the templated
 * class SAMRAI::pdat::CellDataFactory used to construct objects of type
 * LNodeIndexData2.
 */
class LNodeIndexDataFactory2
    : public SAMRAI::pdat::CellDataFactory<NDIM,LNodeIndexSet>
{
public:
    /*!
     * The default constructor for the LNodeIndexDataFactory2 class.  The ghost
     * cell width argument gives the default width for all data objects created
     * with this factory.
     */
    LNodeIndexDataFactory2(
        const SAMRAI::hier::IntVector<NDIM>& ghosts);

    /*!
     * Virtual destructor for the data factory class.
     */
    virtual
    ~LNodeIndexDataFactory2();

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
     * Return whether it is valid to copy this LNodeIndexDataFactory2 to the
     * supplied destination patch data factory. It will return true if dst_pdf
     * is a LNodeIndexDataFactory2, false otherwise.
     */
    virtual bool
    validCopyTo(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchDataFactory<NDIM> >& dst_pdf) const;

    /*!
     * The remaining SAMRAI::hier::PatchDataFactory<NDIM> functionality of
     * LNodeIndexDataFactory2 is provided by the
     * SAMRAI::pdat::CellDataFactory<NDIM,LNodeIndexSet> class.
     */
    using SAMRAI::pdat::CellDataFactory<NDIM,LNodeIndexSet>::getBoxGeometry;
    using SAMRAI::pdat::CellDataFactory<NDIM,LNodeIndexSet>::getGhostCellWidth;
    using SAMRAI::pdat::CellDataFactory<NDIM,LNodeIndexSet>::fineBoundaryRepresentsVariable;
    using SAMRAI::pdat::CellDataFactory<NDIM,LNodeIndexSet>::dataLivesOnPatchBorder;
    using SAMRAI::pdat::CellDataFactory<NDIM,LNodeIndexSet>::getMultiblockDataTranslator;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LNodeIndexDataFactory2();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LNodeIndexDataFactory2(
        const LNodeIndexDataFactory2& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LNodeIndexDataFactory2&
    operator=(
        const LNodeIndexDataFactory2& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/LNodeIndexDataFactory2.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeIndexDataFactory2
