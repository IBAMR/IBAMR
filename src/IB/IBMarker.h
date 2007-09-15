#ifndef included_IBMarker
#define included_IBMarker

// Filename: IBMarker.h
// Last modified: <13.Sep.2007 01:19:46 griffith@box221.cims.nyu.edu>
// Created on 12 Sep 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

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
 * \brief Class IBMarker provides interprocessor communications and
 * database access funtionality to a collection of LNodexIndex objects.
 *
 * \note This class meets the required specification for use with the templated
 * class SAMRAI::pdat::IndexData.
 */
class IBMarker
    : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    IBMarker();

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    IBMarker(
        const IBMarker& from);

    /*!
     * \brief Destructor.
     */
    ~IBMarker();

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBMarker&
    operator=(
        const IBMarker& that);

    /*!
     * \return The number of discrete markers associated with this object.
     */
    int
    getNumberOfMarkers() const;

    /*!
     * \return A const reference to the marker positions.
     */
    const std::vector<double>&
    getPositions() const;

    /*!
     * \return A non-const reference to the marker positions.
     */
    std::vector<double>&
    getPositions();

    /*!
     * \brief Set the marker positions.
     */
    void
    setPositions(
        const std::vector<double>& X);

    /*!
     * \return A const reference to the marker velocities.
     */
    const std::vector<double>&
    getVelocities() const;

    /*!
     * \return A non-const reference to the marker velocities.
     */
    std::vector<double>&
    getVelocities();

    /*!
     * \brief Set the marker velocities.
     */
    void
    setVelocities(
        const std::vector<double>& U);

    /*!
     * \return A const reference to the marker indices.
     */
    const std::vector<int>&
    getIndices() const;

    /*!
     * \return A non-const reference to the marker indices.
     */
    std::vector<int>&
    getIndices();

    /*!
     * \brief Set the marker indices.
     */
    void
    setIndices(
        const std::vector<int>& idxs);

    /*!
     * \return A const reference to the periodic offset.
     *
     * \note If the IBMarker lives in cell i, the index of the source object is
     * src_index = i - offset.
     */
    const SAMRAI::hier::IntVector<NDIM>&
    getPeriodicOffset() const;

    /*!
     * \brief Set the value of the periodic offset.
     *
     * \note If the IBMarker lives in cell i, the index of the source object is
     * src_index = i - offset.
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
        const IBMarker& src_item);

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
        SAMRAI::tbox::AbstractStream& stream) const;

    /*!
     * \brief Unpack data from the input stream.
     */
    void
    unpackStream(
        SAMRAI::tbox::AbstractStream& stream,
        const SAMRAI::hier::IntVector<NDIM>& offset);

    /*!
     * \brief Pack data into a database.
     */
    void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& database) const;

    /*!
     * \brief Unpack data from a database.
     */
    void
    getFromDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& database);

private:
    /*!
     * \brief The marker positions and velocities.
     */
    std::vector<double> d_X, d_U;

    /*!
     * \brief The marker indices.
     */
    std::vector<int> d_idx;

    /*!
     * \brief The periodic offset.
     */
    SAMRAI::hier::IntVector<NDIM> d_offset;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibamr/IBMarker.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBMarker
