#ifndef included_LNodeIndexSet
#define included_LNodeIndexSet

// Filename: LNodeIndexSet.h
// Last modified: <17.Apr.2007 19:32:16 griffith@box221.cims.nyu.edu>
// Created on 29 Feb 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LNodeIndex.h>

// SAMRAI INCLUDES
#include <Index.h>
#include <IntVector.h>
#include <tbox/AbstractStream.h>
#include <tbox/Database.h>
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

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
 * \brief Class LNodeIndexSet provides interprocessor communications and
 * database access funtionality to a collection of LNodexIndex objects.
 *
 * \note This class meets the required specification for use with the templated
 * class SAMRAI::pdat::IndexData.
 */
class LNodeIndexSet
    : public SAMRAI::tbox::DescribedClass
{
public:
    friend class LDataManager;

    /*!
     * \brief The type of the collection.
     */
    typedef std::vector<SAMRAI::tbox::Pointer<LNodeIndex> > IndexSet;

    /*!
     * \brief The type of object, T, stored in the collection.
     */
    typedef IndexSet::value_type value_type;

    /*!
     * \brief SAMRAI::tbox::Pointer to T.
     */
    typedef IndexSet::pointer pointer;

    /*!
     * \brief Reference to T.
     */
    typedef IndexSet::reference reference;

    /*!
     * \brief Const refrence to T.
     */
    typedef IndexSet::const_reference const_reference;

    /*!
     * \brief An unsigned integral type.
     */
    typedef IndexSet::size_type size_type;

    /*!
     * \brief A signed integral type.
     */
    typedef IndexSet::difference_type difference_type;

    /*!
     * \brief Iterator used to iterate through the set.
     */
    typedef IndexSet::iterator iterator;

    /*!
     * \brief Const iterator used to iterate through the collection.
     */
    typedef IndexSet::const_iterator const_iterator;

    /*!
     * \brief Default constructor.
     */
    LNodeIndexSet();

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    LNodeIndexSet(
        const LNodeIndexSet& from);

    /*!
     * \brief Destructor.
     */
    ~LNodeIndexSet();

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LNodeIndexSet&
    operator=(
        const LNodeIndexSet& that);

    /*!
     * \return A const_iterator pointing to the beginning of the set of indices.
     */
    const_iterator
    begin() const;

    /*!
     * \return An iterator pointing to the beginning of the set of indices.
     */
    iterator
    begin();

    /*!
     * \return A const_iterator pointing to the end of the set of indices.
     */
    const_iterator
    end() const;

    /*!
     * \return An iterator pointing to the end of the set of indices.
     */
    iterator
    end();

    /*!
     * \return The size of the set.
     */
    size_type
    size() const;

    /*!
     * \return Whether the set is empty.
     */
    bool
    empty() const;

    /*!
     * \brief Insert a new element at the end (of the set).
     */
    void
    push_back(
        const value_type& value);

    /*!
     * \return A const reference to the periodic offset.
     *
     * \note If the LNodeIndexSet lives in cell i, the index of the source
     * object is src_index = i - offset.
     */
    const SAMRAI::hier::IntVector<NDIM>&
    getPeriodicOffset() const;

    /*!
     * \brief Set the value of the periodic offset.
     *
     * \note If the LNodeIndexSet lives in cell i, the index of the source
     * object is src_index = i - offset.
     */
    void
    setPeriodicOffset(
        const SAMRAI::hier::IntVector<NDIM>& offset);

    /*!
     * \brief Copy data from the source.
     *
     * \note The index of the destination object is src_index + src_offset.
     */
    void
    copySourceItem(
        const SAMRAI::hier::Index<NDIM>& src_index,
        const SAMRAI::hier::IntVector<NDIM>& src_offset,
        const LNodeIndexSet& src_item);

    /*!
     * \brief Return an upper bound on the amount of space required to pack the
     * object to a buffer.
     */
    size_t
    getDataStreamSize() const;

    /*!
     * \brief Pack data into the output stream.
     */
    void
    packStream(
        SAMRAI::tbox::AbstractStream& stream);

    /*!
     * \brief Unpack data from the input stream.
     */
    void
    unpackStream(
        SAMRAI::tbox::AbstractStream& stream,
        const SAMRAI::hier::IntVector<NDIM>& offset);

    /*!
     * \brief Unpack data from a database.
     */
    void
    getFromDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& database);

    /*!
     * \brief Pack data into a database.
     */
    void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& database);

private:
    /*!
     * \brief Assign that to this.
     */
    void
    assignThatToThis(
        const LNodeIndexSet& that);

    /*!
     * \brief Reorder the collection of indices.
     */
    void
    reorderCollection();

    /*!
     * \brief Get rid of any excess capacity in the collection.
     */
    void
    trimToFit();

    /*!
     * \brief The collection of indices.
     */
    IndexSet d_set;

    /*!
     * \brief The periodic offset.
     */
    SAMRAI::hier::IntVector<NDIM> d_offset;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibamr/LNodeIndexSet.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeIndexSet
