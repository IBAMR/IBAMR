#ifndef included_IBStandardForceSpec
#define included_IBStandardForceSpec

// Filename: IBStandardForceSpec.h
// Last modified: <18.Jan.2007 16:43:01 boyce@bigboy.nyconnect.com>
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

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
 * @brief Class IBStandardForceSpec provies a mechanism for specifying
 * a network of spring forces, with support for target points.
 */
class IBStandardForceSpec
    : public Stashable
{
public:
    /*!
     * @brief Register this class and its factory class with the
     * singleton StashableManager object.  This method must be called
     * before any IBStandardForceSpec objects are created.
     *
     * NOTE: This method is collective on all MPI processes.  This is
     * done to ensure that all processes employ the same stashable ID
     * for the IBStandardForceSpec class.
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
    IBStandardForceSpec(
        const std::vector<int>& dst_idxs=vector<int>(),
        const std::vector<double>& stiffnesses=vector<double>(),
        const std::vector<double>& rest_lengths=vector<double>(),
        const std::vector<double>& X_target=vector<double>(),
        const double kappa_target=0.0);

    /*!
     * @brief Destructor.
     */
    ~IBStandardForceSpec();

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
     * \return A const pointer to the target point location.
     */
    const std::vector<double>& getTargetPosition() const;

    /*!
     * \return A const reference to the stiffnesses of the spring that
     * attaches the IB point to the target point.  Note that the
     * stiffness may be zero in the case that this node is not
     * tethered to a target point.
     */
    const double& getTargetStiffness() const;

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
     * @brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    IBStandardForceSpec();

    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    IBStandardForceSpec(
        const IBStandardForceSpec& from);

    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     *
     * @return A reference to this object.
     */
    IBStandardForceSpec& operator=(
        const IBStandardForceSpec& that);

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
     * Data required to define the spring forces.
     */
    int d_num_links;
    std::vector<int> d_dst_idxs;
    std::vector<double> d_stiffnesses, d_rest_lengths;

    /*
     * Data required to define the target penalty force.
     */
    std::vector<double> d_X_target;
    double d_kappa_target;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBStandardForceSpec.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBStandardForceSpec
