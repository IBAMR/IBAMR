#ifndef included_IBEulerianSourceSetter
#define included_IBEulerianSourceSetter

// Filename: IBEulerianSourceSetter.h
// Last modified: <17.Apr.2007 19:45:22 griffith@box221.cims.nyu.edu>
// Created on 18 Jun 2005 by Boyce Griffith (boyce@bigboy.verizon.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/SetDataStrategy.h>

// SAMRAI INCLUDES
#include <Patch.h>
#include <Variable.h>
#include <tbox/Pointer.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBEulerianSourceSetter is used to communicate the Eulerian fluid
 * source-sink distribution computed by class IBHierarchyIntegrator to the
 * incompressible Navier-Stokes solver.
 */
class IBEulerianSourceSetter
    : public IBTK::SetDataStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    IBEulerianSourceSetter(
        const std::string& object_name,
        const int Q_current_idx,
        const int Q_new_idx,
        const int Q_half_idx);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBEulerianSourceSetter();

    /*!
     * \brief Set the current and new times for the present timestep.
     */
    void
    setTimeInterval(
        const double current_time,
        const double new_time);

    /*!
     * \name Methods to set the data.
     */
    //\{

    /*!
     * \note This concrete IBTK::SetDataStrategy is time-dependent.
     */
    virtual bool
    isTimeDependent() const;

    /*!
     * Set the data on the patch interior.
     */
    virtual void
    setDataOnPatch(
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::hier::Patch<NDIM>& patch,
        const double data_time,
        const bool initial_time=false);

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBEulerianSourceSetter();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBEulerianSourceSetter(
        const IBEulerianSourceSetter& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBEulerianSourceSetter&
    operator=(
        const IBEulerianSourceSetter& that);

    /*!
     * The current and new time for the present timestep.
     */
    double d_current_time, d_new_time;

    /*!
     * Patch data descriptor indices for the current, new, and half-time
     * source/sink data.
     */
    const int d_Q_current_idx, d_Q_new_idx, d_Q_half_idx;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBEulerianSourceSetter.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBEulerianSourceSetter
