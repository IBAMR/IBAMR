#ifndef included_LNodeIndex
#define included_LNodeIndex

// Filename: LNodeIndex.h
// Created on 28 Feb 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <11.Apr.2007 02:43:20 boyce@trasnaform2.local>

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/Stashable.h>

// SAMRAI INCLUDES
#include <Index.h>
#include <IntVector.h>
#include <tbox/AbstractStream.h>
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class LNodeIndex provides index information about a single node of a
 * Lagrangian mesh.
 */
class LNodeIndex
    : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \name Friend declarations.
     */
    //\{
    friend bool operator<(
        const LNodeIndex&,
        const LNodeIndex&);
    //\}

    /*!
     * \brief Default constructor.
     */
    LNodeIndex(
        const int lagrangian_nidx=-1,
        const int local_petsc_nidx=-1,
        double* const X_ptr=NULL,
        const std::vector<SAMRAI::tbox::Pointer<Stashable> >& stash_data=vector<SAMRAI::tbox::Pointer<Stashable> >());

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    LNodeIndex(
        const LNodeIndex& from);

    /*!
     * \brief Destructor.
     *
     * The LNodeIndex destructor does nothing interesting.
     */
    ~LNodeIndex();

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LNodeIndex& operator=(
        const LNodeIndex& that);

    /*!
     * \return The Lagrangian index refrenced by this LNodeIndex.
     */
    int getLagrangianIndex() const;

    /*!
     * \brief Reset the Lagrangian index refrenced by this LNodeIndex.
     */
    void setLagrangianIndex(
        const int lagrangian_nidx);

    /*!
     * \return The local PETSc index refrenced by this LNodeIndex.
     */
    int getLocalPETScIndex() const;

    /*!
     * \brief Reset the local PETSc index refrenced by this LNodeIndex.
     */
    void setLocalPETScIndex(
        const int local_petsc_nidx);

    /*!
     * \return A pointer to the physical location of the node refrenced by this
     * LNodeIndex.
     */
    double* getNodeLocation() const;

    /*!
     * \brief Reset the pointer to the physical location of the node refrenced
     * by this LNodeIndex.
     */
    void setNodeLocation(
        double* const X_ptr);

    /*!
     * \return A constant refrence to any additional data associated with the
     * node refrenced by this LNodeIndex.
     */
    const std::vector<SAMRAI::tbox::Pointer<Stashable> >& getStashData() const;

    /*!
     * \return A non-constant refrence to any additional data associated with
     * the node refrenced by this LNodeIndex.
     */
    std::vector<SAMRAI::tbox::Pointer<Stashable> >& getStashData();

    /*!
     * \brief Copy data from the source.
     *
     * \note The cell index of the destination object is src_index + src_offset.
     */
    void copySourceItem(
        const SAMRAI::hier::Index<NDIM>& src_index,
        const SAMRAI::hier::IntVector<NDIM>& src_offset,
        const LNodeIndex& src_item);

    /*!
     * \brief Return an upper bound on the amount of space required to pack the
     * object to a buffer.
     */
    size_t getDataStreamSize() const;

    /*!
     * \brief Pack data into the output stream.
     */
    void packStream(
        SAMRAI::tbox::AbstractStream& stream);

    /*!
     * \brief Unpack data from the input stream.
     */
    void unpackStream(
        SAMRAI::tbox::AbstractStream& stream,
        const SAMRAI::hier::IntVector<NDIM>& offset);

private:
    /*!
     * Assign that to this.
     */
    void assignThatToThis(
        const LNodeIndex& that);

    int d_lagrangian_nidx;  // the fixed global Lagrangian index

    int d_local_petsc_nidx; // the local PETSc index

    double* d_X_ptr;        // a pointer to the physical location of
                            // the node

    // a (possibly empty) collection of objects which are associated with the
    // node
    std::vector<SAMRAI::tbox::Pointer<Stashable> > d_stash_data;
};

/*!
 * \brief Less-than comparison operator.
 *
 * \param lhs  The left-hand-side of the < operator.
 * \param rhs  The right-hand-side of the < operator.
 *
 * \return Whether lhs < rhs.
 *
 * The ordering is determined on the physical locations of the nodes.  When a
 * set of indices is sorted according to operator<(), the nodes are in the
 * "Fortan" ordering according to their physical location.
 */
bool operator<(
    const LNodeIndex& lhs,
    const LNodeIndex& rhs);

}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibamr/LNodeIndex.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeIndex
