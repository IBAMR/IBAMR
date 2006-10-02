//
// LNodeIndexSet.h
//
// Created on 29 Feb 2004
//         by Boyce Griffith (boyce@bigboy.speakeasy.net).
//
// Last modified: <29.Jun.2005 16:27:45 boyce@mstu1.cims.nyu.edu>
//

#ifndef included_LNodeIndexSet
#define included_LNodeIndexSet

// STL INCLUDES
//
#include <vector>

// SAMRAI-tools INCLUDES
//
#include "LNodeIndex.h"

// SAMRAI INCLUDES
//
#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#include "Index.h"
#include "IntVector.h"
#include "tbox/AbstractStream.h"
#include "tbox/Database.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

using namespace SAMRAI;
using namespace std;

// FORWARD DECLARATIONS
//
class LDataManager;

// CLASS DEFINITION
//

/*!
 * @brief Class LNodeIndexSet provides interprocessor communications
 * and database access funtionality to a collection of LNodexIndex
 * objects.  This class meets the required specification for use with
 * the templated pdat::IndexData<NDIM> class.
 */
class LNodeIndexSet
    : public tbox::DescribedClass
{
public:
    //@{ @name Class typedefs.

    /*!
     * @brief The type of the collection.
     */
    typedef vector<tbox::Pointer<LNodeIndex> > IndexSet;

    /*!
     * @brief The type of object, T, stored in the collection.
     */
    typedef IndexSet::value_type value_type;

    /*!
     * @brief tbox::Pointer to T.
     */
    typedef IndexSet::pointer pointer;

    /*!
     * @brief Reference to T.
     */
    typedef IndexSet::reference reference;
    
    /*!
     * @brief Const refrence to T.
     */
    typedef IndexSet::const_reference const_reference;

    /*!
     * @brief An unsigned integral type.
     */
    typedef IndexSet::size_type size_type;

    /*!
     * @brief A signed integral type.
     */
    typedef IndexSet::difference_type difference_type;

    /*!
     * @brief Iterator used to iterate through the set.
     */
    typedef IndexSet::iterator iterator;

    /*!
     * @brief Const iterator used to iterate through the collection.
     */
    typedef IndexSet::const_iterator const_iterator;

    //@}
    
    //@{ @name Friend declarations.
    friend class LDataManager;
    //@}
    
    /*!
     * @brief Default constructor.
     */
    LNodeIndexSet();
    
    /*!
     * @brief Copy constructor.
     *
     * @param from The value to copy to this object.
     */
    LNodeIndexSet(
        const LNodeIndexSet& from);
    
    /*!
     * @brief Destructor.
     */
    ~LNodeIndexSet();
    
    /*!
     * @brief Assignment operator.
     *
     * @param that The value to assign to this object.
     * 
     * @return A reference to this object.
     */
    LNodeIndexSet& operator=(
        const LNodeIndexSet& that);

    /*!
     * @return A const_iterator pointing to the beginning of the set
     * of indices.
     */
    const_iterator begin() const;

    /*!
     * @return An iterator pointing to the beginning of the set of
     * indices.
     */
    iterator begin();

    /*!
     * @return A const_iterator pointing to the end of the set of
     * indices.
     */
    const_iterator end() const;

    /*!
     * @return An iterator pointing to the end of the set of indices.
     */
    iterator end();

    /*!
     * @return The size of the set.
     */
    size_type size() const;

    /*!
     * @return Whether the set is empty.
     */
    bool empty() const;

    /*!
     * @brief Insert a new element at the end (of the set).
     */
    void push_back(
        const value_type& value);

    /*!
     * @return A const reference to the periodic offset.
     *
     * NOTE: If the LNodeIndexSet lives in cell i, the index of the
     * source object is src_index = i - offset.
     */
    const hier::IntVector<NDIM>& getPeriodicOffset() const;
    
    /*!
     * @brief Set the value of the periodic offset.
     *
     * NOTE: If the LNodeIndexSet lives in cell i, the index of the
     * source object is src_index = i - offset.
     */
    void setPeriodicOffset(
        const hier::IntVector<NDIM>& offset);
    
    /*!
     * @brief Copy data from the source.
     *
     * NOTE: The index of the destination object is src_index +
     * src_offset.
     */
    void copySourceItem(
        const hier::Index<NDIM>& src_index,
        const hier::IntVector<NDIM>& src_offset,
        const LNodeIndexSet& src_item);

    /*!
     * @brief Return an upper bound on the amount of space required to
     * pack the object to a buffer.
     */
    size_t getDataStreamSize() const;

    /*!
     * @brief Pack data into the output stream.
     */
    void packStream(
        tbox::AbstractStream& stream);

    /*!
     * @brief Unpack data from the input stream.
     */
    void unpackStream(
        tbox::AbstractStream& stream,
        const hier::IntVector<NDIM>& offset);

    /*!
     * @brief Unpack data from a database.
     */
    void getFromDatabase(
        tbox::Pointer<tbox::Database>& database);

    /*!
     * @brief Pack data into a database.
     */
    void putToDatabase(
        tbox::Pointer<tbox::Database>& database);

private:
    /*!
     * @brief Assign that to this.
     */
    void assignThatToThis(
        const LNodeIndexSet& that);
    
    /*!
     * @brief Reorder the collection of indices.
     */
    void reorderCollection();

    /*!
     * @brief Get rid of any excess capacity in the collection.
     */
    void trimToFit();

    /*!
     * @brief The collection of indices.
     */
    IndexSet d_set;

    /*!
     * @brief The periodic offset.
     */
    hier::IntVector<NDIM> d_offset;
};

// INLINED FUNCTION DEFINITIONS
//
#ifndef DEBUG_NO_INLINE
#include "LNodeIndexSet.I"
#endif

#endif //#ifndef included_LNodeIndexSet

//////////////////////////////////////////////////////////////////////////////
