#ifndef included_LNodeIndexData2
#define included_LNodeIndexData2

// Filename: LNodeIndexData2.h
// Last modified: <04.Jun.2007 16:29:12 griffith@box221.cims.nyu.edu>
// Created on 04 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/ArrayData_specialized_LNodeIndexSet.I>
#include <ibamr/LNodeIndexSet.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CellData.h>
#include <IntVector.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// FORWARD DECLARATIONS /////////////////////////

namespace IBAMR
{
class LDataManager;
}// namespace IBAMR

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class LNodeIndexData2 is a specialization of the templated class
 * SAMRAI::pdat::CellData that provides access to the Lagrangian and PETSc
 * index information and any Stashable data associated with the the Lagrangian
 * nodes in the interior and ghost cell region of a Cartesian grid patch.
 */
class LNodeIndexData2
    : public SAMRAI::pdat::CellData<NDIM,LNodeIndexSet>
{
public:
    friend class LDataManager;

    /*!
     * The cell iterator iterates over the elements of a cell centered box
     * geometry.
     */
    typedef SAMRAI::pdat::CellData<NDIM,LNodeIndexSet>::Iterator Iterator;

    /*!
     * The constructor for an SAMRAI::pdat::CellData<NDIM> object.  The box
     * describes the interior of the index space and the ghosts vector describes
     * the ghost nodes in each coordinate direction.
     */
    LNodeIndexData2(
        const SAMRAI::hier::Box<NDIM>& box,
        const SAMRAI::hier::IntVector<NDIM>& ghosts,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena> pool=SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena>(NULL));

    /*!
     * The virtual destructor for an LNodeIndexData2 object.
     */
    virtual
    ~LNodeIndexData2();

    /*!
     * Determines whether the patch data subclass can estimate the necessary
     * stream size using only index space information.
     */
    virtual bool
    canEstimateStreamSizeFromBox() const;

    /*!
     * Most of the SAMRAI::hier::PatchData<NDIM> functionality of class
     * LNodeIndexData is provided by class
     * SAMRAI::pdat::CellData<NDIM,LNodeIndexSet>.
     */
    using SAMRAI::pdat::CellData<NDIM,LNodeIndexSet>::copy;
    using SAMRAI::pdat::CellData<NDIM,LNodeIndexSet>::copy2;
    using SAMRAI::pdat::CellData<NDIM,LNodeIndexSet>::operator();
    using SAMRAI::pdat::CellData<NDIM,LNodeIndexSet>::getDataStreamSize;
    using SAMRAI::pdat::CellData<NDIM,LNodeIndexSet>::packStream;
    using SAMRAI::pdat::CellData<NDIM,LNodeIndexSet>::unpackStream;
    using SAMRAI::pdat::CellData<NDIM,LNodeIndexSet>::getSpecializedFromDatabase;
    using SAMRAI::pdat::CellData<NDIM,LNodeIndexSet>::putSpecializedToDatabase;

    /*!
     * \return A constant refrence to the set of local PETSc indices to
     * Lagrangian nodes which lie in the patch interior.
     *
     * \note This should be an ordered, contiguous set of indices.
     */
    const std::vector<int>&
    getInteriorLocalIndices() const;

    /*!
     * \return A constant refrence to the set of local PETSc indices to
     * Lagrangian nodes which lie in ghost cell region of this
     * SAMRAI::hier::PatchData<NDIM> object.
     */
    const std::vector<int>&
    getGhostLocalIndices() const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LNodeIndexData2();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LNodeIndexData2(
        const LNodeIndexData2& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LNodeIndexData2&
    operator=(
        const LNodeIndexData2& that);

    /*!
     * Each LNodeIndexData2 object maintains collections of local (PETSc)
     * indices corresponding to the local and ghost nodes which lie in the patch
     * data ghost box.
     */
    std::vector<int> d_interior_local_indices;
    std::vector<int> d_ghost_local_indices;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibamr/LNodeIndexData2.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeIndexData2
