#ifndef included_LNodeIndexDataFactory
#define included_LNodeIndexDataFactory

// Filename: LNodeIndexDataFactory.h
// Last modified: <04.Jun.2007 12:32:26 griffith@box221.cims.nyu.edu>
// Created on 01 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LNodeIndexSet.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <IndexDataFactory.h>
#include <PatchData.h>
#include <PatchDataFactory.h>
#include <tbox/Arena.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class LNodeIndexDataFactory is a specialization of the templated class
 * SAMRAI::pdat::IndexDataFactory used to construct objects of type
 * LNodeIndexData.
 */
class LNodeIndexDataFactory
    : public SAMRAI::pdat::IndexDataFactory<NDIM,LNodeIndexSet>
{
public:
    /*!
     * The default constructor for the LNodeIndexDataFactory class.  The ghost
     * cell width argument gives the default width for all irregular data
     * objects created with this factory.
     */
    LNodeIndexDataFactory(
        const SAMRAI::hier::IntVector<NDIM>& ghosts);

    /*!
     * Virtual destructor for the irregular data factory class.
     */
    virtual
    ~LNodeIndexDataFactory();

    /*!
     * Virtual function to clone the irregular data factory.  This will return a
     * new instantiation of the factory with the same properties (e.g., same
     * type and ghost cell width).  The properties of the cloned factory can
     * then be changed (e.g., change the ghost cell width) without modifying the
     * original.
     */
    virtual SAMRAI::tbox::Pointer<SAMRAI::hier::PatchDataFactory<NDIM> >
    cloneFactory();

    /*!
     * Virtual factory function to allocate a concrete irregular data object.
     * The default information about the object (e.g., ghost cell width) is
     * taken from the factory.  If no memory pool is provided, the allocation
     * routine assumes some default memory pool.
     */
    virtual SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >
    allocate(
        const SAMRAI::hier::Box<NDIM>& box,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena> pool=NULL) const;

    /*!
     * Calculate the amount of memory needed to store the irregular data object,
     * including object data but not dynamically allocated data.  Because the
     * irregular data list can grow and shrink, it would be impossible to
     * estimate the necessary amount of memory.  Instead, dynamic data is
     * allocated via the standard new/free mechanisms.
     */
    virtual size_t
    getSizeOfMemory(
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * The remaining SAMRAI::hier::PatchDataFactory<NDIM> functionality of
     * LNodeIndexDataFactory is provided by the
     * SAMRAI::pdat::IndexDataFactory<NDIM,LNodeIndexSet> class.
     */
    using SAMRAI::pdat::IndexDataFactory<NDIM,LNodeIndexSet>::getBoxGeometry;
    using SAMRAI::pdat::IndexDataFactory<NDIM,LNodeIndexSet>::getDefaultGhostCellWidth;
    using SAMRAI::pdat::IndexDataFactory<NDIM,LNodeIndexSet>::setDefaultGhostCellWidth;

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
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/LNodeIndexDataFactory.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeIndexDataFactory
