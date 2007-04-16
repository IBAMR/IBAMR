#ifndef included_StashableManager
#define included_StashableManager

// Filename: StashableManager.h
// Created on 14 Jun 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <16.Apr.2007 01:40:20 boyce@trasnaform2.local>

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/Stashable.h>
#include <ibamr/StashableFactory.h>

// SAMRAI INCLUDES
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <map>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class StashableManager is a singleton manager class that handles
 * object packing and unpacking for SAMRAI::tbox::AbstractStream based
 * communication via the Stashable and StashableFactory interfaces.
 *
 * \see Stashable
 * \see StashableFactory
 */
class StashableManager
{
public:
    /*!
     * Return a pointer to the instance of the Stashable manager.  All access to
     * the singleton StashableManager object is through the getManager()
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
    static StashableManager*
    getManager();

    /*!
     * Deallocate the StashableManager instance.
     *
     * It is not necessary to call this function at program termination, since
     * it is automatically called by the ShutdownRegistry class.
     */
    static void
    freeManager();

    /*!
     * \return A integer value reserved for unregistered StashableFactory
     * objects.  A concrete StashableFactory object must use this as its initial
     * stashable ID.
     */
    static int
    getUnregisteredID();

    /*!
     * Check to see if a StashableFactory has been registered with the manager.
     *
     * \return true if the factory has been registered, false otherwise.
     *
     * \note This method simply checks to see if a StashableFactory with the
     * same Stashable ID has been registered with the manager.  Every different
     * Stashable/StashableFactory type \em must have a unique ID.
     */
    bool
    checkFactoryRegistration(
        SAMRAI::tbox::Pointer<StashableFactory> factory);

    /*!
     * Register a StashableFactory with the manager.
     *
     * Each factory object registered with the manager is provided with a unique
     * ID.
     *
     * \note To ensure that each MPI process uses the same stashable ID for each
     * stashable class registered with the manager, this method is collective on
     * all MPI processes!
     */
    int
    registerFactory(
        SAMRAI::tbox::Pointer<StashableFactory> factory);

    /*!
     * \brief Return an upper bound on the amount of space required to pack a
     * Stashable object to a buffer.
     */
    size_t
    getDataStreamSize(
        const SAMRAI::tbox::Pointer<Stashable>& stash_data) const;

    /*!
     * \brief Return an upper bound on the amount of space required to pack a
     * vector of Stashable objects to a buffer.
     */
    size_t
    getDataStreamSize(
        const std::vector<SAMRAI::tbox::Pointer<Stashable> >& stash_data) const;

    /*!
     * \brief Pack a Stashable object into the output stream.
     */
    void
    packStream(
        SAMRAI::tbox::AbstractStream& stream,
        SAMRAI::tbox::Pointer<Stashable>& stash_data);

    /*!
     * \brief Pack a vector of Stashable objects into the output stream.
     */
    void
    packStream(
        SAMRAI::tbox::AbstractStream& stream,
        std::vector<SAMRAI::tbox::Pointer<Stashable> >& stash_data);

    /*!
     * \brief Unpack a Stashable object from the input stream.
     */
    void
    unpackStream(
        SAMRAI::tbox::AbstractStream& stream,
        const SAMRAI::hier::IntVector<NDIM>& offset,
        SAMRAI::tbox::Pointer<Stashable>& stash_data);

    /*!
     * \brief Unpack a vector of Stashable objects from the input stream.
     */
    void
    unpackStream(
        SAMRAI::tbox::AbstractStream& stream,
        const SAMRAI::hier::IntVector<NDIM>& offset,
        std::vector<SAMRAI::tbox::Pointer<Stashable> >& stash_data);

protected:
    /*!
     * \brief Constructor.
     */
    StashableManager();

    /*!
     * \brief Destructor.
     */
    ~StashableManager();

    /*!
     * Generate a unique ID number.
     *
     * Every call to getUniqueID() returns a different integer, simplifying the
     * task of generating ID numbers for StashableFactory objects.
     */
    static int
    getUniqueID();

private:
    typedef std::map<int,SAMRAI::tbox::Pointer<StashableFactory> > StashableFactoryMap;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StashableManager(
        const StashableManager& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StashableManager&
    operator=(
        const StashableManager& that);

    /*!
     * Static data members used to control access to and destruction of
     * singleton data manager instance.
     */
    static StashableManager* s_data_manager_instance;
    static bool s_registered_callback;
    static unsigned char s_shutdown_priority;

    /*!
     * Static data members used to simplify the process of assigning unique ID
     * numbers to each registered StashableFactory object.
     */
    static int s_current_id_number;
    static const int s_unregistered_number;

    /*!
     * Map from Stashable ID to registered StashableFactory objects.
     */
    StashableFactoryMap d_factory_map;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibamr/StashableManager.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_StashableManager
