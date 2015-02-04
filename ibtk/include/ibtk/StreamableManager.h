// Filename: StreamableManager.h
// Created on 14 Jun 2004 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

#ifndef included_StreamableManager
#define included_StreamableManager

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <map>
#include <vector>

#include "ibtk/StreamableFactory.h"
#include "tbox/Pointer.h"

namespace IBTK
{
class Streamable;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class IntVector;
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
 * \brief Class StreamableManager is a singleton manager class that organizes the
 * actual packing and unpacking of concrete Streamable objects for
 * SAMRAI::tbox::AbstractStream based communication.
 *
 * \see Streamable
 * \see StreamableFactory
 */
class StreamableManager
{
public:
    /*!
     * Return a pointer to the instance of the Streamable manager.  All access to
     * the singleton StreamableManager object is through the getManager()
     * function.
     *
     * Note that when the manager is accessed for the first time, the
     * freeManager static method is registered with the ShutdownRegistry class.
     * Consequently, an allocated manager is freed at program completion.  Thus,
     * users of this class do not explicitly allocate or deallocate the manager
     * instances.
     *
     * \return A pointer to the data manager instance.
     */
    static StreamableManager* getManager();

    /*!
     * Deallocate the StreamableManager instance.
     *
     * It is not necessary to call this function at program termination, since
     * it is automatically called by the ShutdownRegistry class.
     */
    static void freeManager();

    /*!
     * \return A integer value reserved for unregistered StreamableFactory
     * objects.  A concrete StreamableFactory object must use this as its initial
     * streamable ID.
     */
    static int getUnregisteredID();

    /*!
     * Check to see if a StreamableFactory has been registered with the manager.
     *
     * \return true if the factory has been registered, false otherwise.
     *
     * \note This method simply checks to see if a StreamableFactory with the
     * same Streamable ID has been registered with the manager.  Every different
     * Streamable/StreamableFactory type \em must have a unique ID.
     */
    bool checkFactoryRegistration(SAMRAI::tbox::Pointer<StreamableFactory> factory);

    /*!
     * Register a StreamableFactory with the manager.
     *
     * Each factory object registered with the manager is provided with a unique
     * ID.
     *
     * \note To ensure that each MPI process uses the same streamable ID for each
     * streamable class registered with the manager, this method is collective on
     * all MPI processes!
     */
    int registerFactory(SAMRAI::tbox::Pointer<StreamableFactory> factory);

    /*!
     * \brief Return an upper bound on the amount of space required to pack a
     * Streamable object to a buffer.
     */
    size_t getDataStreamSize(SAMRAI::tbox::Pointer<Streamable> data_item) const;

    /*!
     * \brief Return an upper bound on the amount of space required to pack a
     * vector of Streamable objects to a buffer.
     */
    size_t getDataStreamSize(const std::vector<SAMRAI::tbox::Pointer<Streamable> >& data_items) const;

    /*!
     * \brief Pack a Streamable object into the output stream.
     */
    void packStream(SAMRAI::tbox::AbstractStream& stream, SAMRAI::tbox::Pointer<Streamable> data_item);

    /*!
     * \brief Pack a vector of Streamable objects into the output stream.
     */
    void packStream(SAMRAI::tbox::AbstractStream& stream, std::vector<SAMRAI::tbox::Pointer<Streamable> >& data_items);

    /*!
     * \brief Unpack a Streamable object from the data stream.
     */
    SAMRAI::tbox::Pointer<Streamable> unpackStream(SAMRAI::tbox::AbstractStream& stream,
                                                   const SAMRAI::hier::IntVector<NDIM>& offset);

    /*!
     * \brief Unpack a vector of Streamable objects from the data stream.
     */
    void unpackStream(SAMRAI::tbox::AbstractStream& stream,
                      const SAMRAI::hier::IntVector<NDIM>& offset,
                      std::vector<SAMRAI::tbox::Pointer<Streamable> >& data_items);

protected:
    /*!
     * \brief Constructor.
     */
    StreamableManager();

    /*!
     * \brief Destructor.
     */
    ~StreamableManager();

    /*!
     * Generate a unique ID number.
     *
     * Every call to createUniqueID() returns a different integer, simplifying
     * the task of generating ID numbers for StreamableFactory objects.
     */
    static int createUniqueID();

private:
    typedef std::map<int, SAMRAI::tbox::Pointer<StreamableFactory> > StreamableFactoryMap;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StreamableManager(const StreamableManager& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StreamableManager& operator=(const StreamableManager& that);

    /*!
     * Static data members used to control access to and destruction of
     * singleton data manager instance.
     */
    static StreamableManager* s_data_manager_instance;
    static bool s_registered_callback;
    static unsigned char s_shutdown_priority;

    /*!
     * Static data members used to simplify the process of assigning unique ID
     * numbers to each registered StreamableFactory object.
     */
    static int s_current_id_number;
    static const int s_unregistered_id_number;

    /*!
     * Map from Streamable ID to registered StreamableFactory objects.
     */
    StreamableFactoryMap d_factory_map;
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/StreamableManager-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_StreamableManager
