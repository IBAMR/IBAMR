#ifndef included_SpringForceSpec
#define included_SpringForceSpec

// Filename: SpringForceSpec.h
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)
// Last modified: <04.Oct.2006 19:53:25 boyce@boyce-griffiths-powerbook-g4-15.local>

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/Stashable.h>

// SAMRAI INCLUDES
#include <tbox/AbstractStream.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * @brief Class SpringForceSpec provies a mechanism for specifying a
 * network of spring forces.
 */
class SpringForceSpec
    : public Stashable
{
public:
    /*!
     * @brief Register this class and its factory class with the
     * singleton StashableManager object.  This method must be called
     * before any SpringForceSpec objects are created.
     *
     * NOTE: This method is collective on all MPI processes.  This is
     * done to ensure that all processes employ the same stashable ID
     * for the SpringForceSpec class.
     */
    static void registerWithStashableManager();

    /*!
     * @brief Returns a boolean indicating whether the class has been
     * registered with the singleton StashableManager object.
     */
    static bool getIsRegisteredWithStashableManager();

    /*!
     * @brief Default constructor.
     */
    SpringForceSpec(
        const std::vector<int>& dst_idxs=vector<int>(),
        const std::vector<double>& stiffnesses=vector<double>(),
        const std::vector<double>& rest_lengths=vector<double>());

    /*!
     * @brief Destructor.
     */
    ~SpringForceSpec();

    /*!
     * @return The number of links attatched to the source node.
     */
    int getNumberOfLinks() const;

    /*!
     * @return A const refrence to the destination node indices for
     * each spring attached to the source node.
     */
    const std::vector<int>& getDestinationNodeIndices() const;

    /*!
     * @return A non-const refrence to the destination node indices
     * for each spring attached to the source node.
     */
    std::vector<int>& getDestinationNodeIndices();

    /*!
     * @return A const reference to the stiffnesses of each spring
     * attached to the source node.
     */
    const std::vector<double>& getStiffnesses() const;

    /*!
     * @return A non-const reference to the stiffnesses of each spring
     * attached to the source node.
     */
    std::vector<double>& getStiffnesses();

    /*!
     * @return A const reference to the resting length of each spring
     * attached to the source node.
     */
    const std::vector<double>& getRestingLengths() const;

    /*!
     * @return A non-const reference to the resting length of each
     * spring attached to the source node.
     */
    std::vector<double>& getRestingLengths();

    /*!
     * @brief Return the unique identifier used to specify the
     * StashableFactory object used by the StashableManager to extract
     * Stashable objects from data streams.
     */
    int getStashableID() const;

    /*!
     * @brief Return an upper bound on the amount of space required to
     * pack the object to a buffer.
     */
    size_t getDataStreamSize() const;

    /*!
     * @brief Pack data into the output stream.
     */
    void packStream(
        SAMRAI::tbox::AbstractStream& stream);

private:
    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    SpringForceSpec(
        const SpringForceSpec& from);

    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     *
     * @return A reference to this object.
     */
    SpringForceSpec& operator=(
        const SpringForceSpec& that);

    /*
     * Indicates whether the factory has been registered with the
     * StashableManager.
     */
    static bool s_registered_factory;

    /*
     * The stashable ID for this object type.
     */
    static int s_stashable_id;

    /*
     * Data required to define the spring force.
     */
    int d_num_links;
    std::vector<int> d_dst_idxs;
    std::vector<double> d_stiffnesses, d_rest_lengths;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/SpringForceSpec.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_SpringForceSpec
