#ifndef included_TargetPointForceSpec
#define included_TargetPointForceSpec

// Filename: TargetPointForceSpec.h
// Created on 23 Oct 2006 by Boyce Griffith (boyce@bigboy.nyconnect.com)
// Last modified: <27.Oct.2006 00:22:57 boyce@bigboy.nyconnect.com>

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
 * @brief Class TargetPointForceSpec provies a mechanism for
 * specifying target points in an IB computation.
 */
class TargetPointForceSpec
    : public Stashable
{
public:
    /*!
     * \brief Register this class and its factory class with the
     * singleton StashableManager object.  This method must be called
     * before any TargetPointForceSpec objects are created.
     *
     * \note This method is collective on all MPI processes.  This is
     * done to ensure that all processes employ the same stashable ID
     * for the TargetPointForceSpec class.
     */
    static void registerWithStashableManager();

    /*!
     * \brief Returns a boolean indicating whether the class has been
     * registered with the singleton StashableManager object.
     */
    static bool getIsRegisteredWithStashableManager();

    /*!
     * \brief Constructor.
     */
    TargetPointForceSpec(
        const std::vector<double>& X,
        const double kappa);

    /*!
     * \brief Constructor.
     */
    TargetPointForceSpec(
        const double* const X,
        const double kappa);

    /*!
     * \brief Destructor.
     */
    ~TargetPointForceSpec();

    /*!
     * \return A const pointer to the target point location.
     */
    const std::vector<double>& getPosition() const;

    /*!
     * \return A const reference to the stiffnesses of the spring that
     * attaches the IB point to the target point.
     */
    const double& getStiffness() const;

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
    TargetPointForceSpec();

    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    TargetPointForceSpec(
        const TargetPointForceSpec& from);

    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     *
     * @return A reference to this object.
     */
    TargetPointForceSpec& operator=(
        const TargetPointForceSpec& that);

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
     * Data required to define the target penalty force.
     */
    const std::vector<double> d_X;
    const double d_kappa;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/TargetPointForceSpec.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_TargetPointForceSpec
