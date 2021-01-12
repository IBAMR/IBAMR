// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_LNode
#define included_IBTK_LNode

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LNodeIndex.h"
#include "ibtk/Streamable.h"
#include "ibtk/ibtk_utilities.h"

#include "IntVector.h"
#include "tbox/Pointer.h"

#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Index;
} // namespace hier
namespace tbox
{
class AbstractStream;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LNode is the basic element of an LMesh.
 *
 * Class LNode provides Lagrangian and <A
 * HREF="http://www.mcs.anl.gov/petsc">PETSc</A> indexing information and data
 * storage for a single node of a Lagrangian mesh.
 */
class LNode : public LNodeIndex
{
public:
    /*!
     * \brief Default constructor.
     *
     * \note Any nonzero periodic offset/displacement must already be registered
     * with any provided node data items.
     */
    LNode(int lagrangian_nidx = -1,
          int global_petsc_nidx = -1,
          int local_petsc_nidx = -1,
          const SAMRAI::hier::IntVector<NDIM>& initial_periodic_offset = SAMRAI::hier::IntVector<NDIM>(0),
          const SAMRAI::hier::IntVector<NDIM>& current_periodic_offset = SAMRAI::hier::IntVector<NDIM>(0),
          const Vector& initial_periodic_displacement = Vector::Zero(),
          const Vector& current_periodic_displacement = Vector::Zero(),
          const std::vector<SAMRAI::tbox::Pointer<Streamable> >& node_data =
              std::vector<SAMRAI::tbox::Pointer<Streamable> >());

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    LNode(const LNode& from);

    /*!
     * \brief Constructor that unpacks data from an input stream.
     */
    LNode(SAMRAI::tbox::AbstractStream& stream, const SAMRAI::hier::IntVector<NDIM>& offset);

    /*!
     * \brief Destructor.
     */
    ~LNode();

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LNode& operator=(const LNode& that);

    /*!
     * \return A constant reference to any additional data items associated with
     * the node referenced by this LNode object.
     */
    const std::vector<SAMRAI::tbox::Pointer<Streamable> >& getNodeData() const;

    /*!
     * \brief Reset the collection of additional data items associated with the
     * node referenced by this LNode object.
     */
    void setNodeData(const std::vector<SAMRAI::tbox::Pointer<Streamable> >& node_data);

    /*!
     * \brief Append a data item to the collection of data items associated with
     * this node.  The appended item will appear at the end of the vector of
     * node data items associated with this node.
     */
    void appendNodeDataItem(const SAMRAI::tbox::Pointer<Streamable>& node_data_item);

    /*!
     * \brief Remove a data item to the collection of data items associated with
     * this node.  If the argument is not associated with the collection of node
     * data associated with this node, this method will have no effect.
     *
     * \note Removing items is potentially an inefficient operation.
     */
    void removeNodeDataItem(const SAMRAI::tbox::Pointer<Streamable>& node_data_item);

    /*!
     * \return A pointer to the first data item of type T associated with the
     * node referenced by the LNode object.
     *
     * If no object of the specified type is encountered, this method returns a
     * null pointer.
     *
     * \note It is possible for multiple objects of the same type to be
     * associated with each node.  This method returns only the \em first such
     * object encountered in the collection of data items associated with the
     * node.
     */
    template <typename T>
    T* getNodeDataItem() const;

    /*!
     * \return A vector of pointers to all data items of type T associated with
     * the node referenced by the LNode object.
     *
     * If no object of the specified type is encountered, this method returns an
     * empty vector.
     *
     * \note It is possible for multiple objects of the same type to be
     * associated with each node.  This method returns a vector of \em all such
     * objects encountered in the collection of data items associated with the
     * node.
     */
    template <typename T>
    std::vector<T*> getNodeDataVector() const;

    /*!
     * \brief Indicate that the LNode object has been shifted across a periodic
     * boundary.
     */
    void registerPeriodicShift(const SAMRAI::hier::IntVector<NDIM>& offset, const Vector& displacement) override;

    /*!
     * \brief Copy data from the source.
     *
     * \note The cell index of the destination object is src_index + src_offset.
     */
    void copySourceItem(const SAMRAI::hier::Index<NDIM>& src_index,
                        const SAMRAI::hier::IntVector<NDIM>& src_offset,
                        const LNodeIndex& src_item) override;

    /*!
     * \brief Return an upper bound on the amount of space required to pack the
     * object to a buffer.
     */
    size_t getDataStreamSize() const override;

    /*!
     * \brief Pack data into the output stream.
     */
    void packStream(SAMRAI::tbox::AbstractStream& stream) override;

    /*!
     * \brief Unpack data from the input stream.
     */
    virtual void unpackStream(SAMRAI::tbox::AbstractStream& stream,
                              const SAMRAI::hier::IntVector<NDIM>& offset) override;

private:
    /*!
     * Assign that to this.
     */
    void assignThatToThis(const LNode& that);

    /*!
     * Setup node data map.
     */
    void setupNodeDataTypeArray();

    // a (possibly empty) collection of data objects that are associated with
    // the node
    std::vector<SAMRAI::tbox::Pointer<Streamable> > d_node_data;
    static const short int MAX_SIZE = 8;
    Streamable* d_node_data_type_arr[MAX_SIZE];
};

} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/LNode-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LNode
