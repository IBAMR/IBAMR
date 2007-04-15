#ifndef included_LNodeIndexData
#define included_LNodeIndexData

// Filename: LNodeIndexData.h
// Created on 01 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <14.Apr.2007 21:04:56 griffith@box221.cims.nyu.edu>

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LNodeIndexSet.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <IndexData.h>
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
 * \brief Class LNodeIndexData is a specialization of the templated
 * SAMRAI::pdat::IndexData<NDIM> class.  It extends the SAMRAI::pdat::IndexData
 * class by providing the IDs of the Lagrangian nodes in the patch interior and
 * in the ghost cell region of the patch.
 */
class LNodeIndexData
    : public SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>
{
public:
    /*!
     * \name Friend declarations.
     */
    //\{
    friend class SAMRAI::pdat::IndexIterator<NDIM,LNodeIndexSet>;
    friend class LDataManager;
    //\}

    /*!
     * The constructor for an SAMRAI::pdat::IndexData<NDIM> object.  The box
     * describes the interior of the index space and the ghosts vector describes
     * the ghost nodes in each coordinate direction.
     */
    LNodeIndexData(
        const SAMRAI::hier::Box<NDIM>& box,
        const SAMRAI::hier::IntVector<NDIM>& ghosts);

    /*!
     * The virtual destructor for an LNodeIndexData object.
     */
    virtual ~LNodeIndexData();

    /*!
     * A fast copy between the source and destination.  All data is copied from
     * the source into the destination where there is overlap in the index
     * space.
     */
    virtual void copy(
        const SAMRAI:: hier::PatchData<NDIM>& src);
    virtual void copy2(
        SAMRAI::hier::PatchData<NDIM>& dst) const;

    /*!
     * Copy data from the source into the destination using the designated
     * overlap descriptor.  The overlap description should have been computed
     * previously from computeIntersection().
     */
    virtual void copy(
        const SAMRAI::hier::PatchData<NDIM>& src,
        const SAMRAI::hier::BoxOverlap<NDIM>& overlap);
    virtual void copy2(
        SAMRAI::hier::PatchData<NDIM>& dst,
        const SAMRAI::hier::BoxOverlap<NDIM>& overlap) const;

    /*!
     * The SAMRAI::hier::PatchData<NDIM> functionality of LNodeIndexData is
     * provided by the SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet> class.
     */
    using SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>::canEstimateStreamSizeFromBox;
    using SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>::getDataStreamSize;
    using SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>::packStream;
    using SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>::unpackStream;
    using SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>::appendItem;
    using SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>::addItem;
    using SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>::removeItem;
    using SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>::removeInsideBox;
    using SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>::removeOutsideBox;
    using SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>::removeGhostItems;
    using SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>::removeAllItems;
    using SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>::isElement;
    using SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>::getItem;
    using SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>::getSpecializedFromDatabase;
    using SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>::putSpecializedToDatabase;

    /*!
     * \return A constant refrence to the set of local PETSc indices to
     * Lagrangian nodes which lie in the patch interior.
     *
     * \note This should be an ordered, contiguous set of indices.
     */
    const std::vector<int>& getInteriorLocalIndices() const;

    /*!
     * \return A constant refrence to the set of local PETSc indices to
     * Lagrangian nodes which lie in ghost cell region of this
     * SAMRAI::hier::PatchData<NDIM> object.
     */
    const std::vector<int>& getGhostLocalIndices() const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LNodeIndexData();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LNodeIndexData(
        const LNodeIndexData& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LNodeIndexData& operator=(
        const LNodeIndexData& that);

    /*!
     * Each LNodeIndexData object maintains collections of local (PETSc) indices
     * corresponding to the local and ghost nodes which lie in the patch data
     * ghost box.
     */
    std::vector<int> d_interior_local_indices;
    std::vector<int> d_ghost_local_indices;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibamr/LNodeIndexData.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeIndexData
