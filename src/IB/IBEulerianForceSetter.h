#ifndef included_IBEulerianForceSetter
#define included_IBEulerianForceSetter

// Filename: IBEulerianForceSetter.h
// Last modified: <03.Nov.2009 21:10:45 griffith@griffith-macbook-pro.local>
// Created on 28 Sep 2004 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)

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
 * \brief Class IBEulerianForceSetter is used to communicate the Eulerian body
 * force computed by class IBHierarchyIntegrator to the incompressible
 * Navier-Stokes solver.
 */
class IBEulerianForceSetter
    : public IBTK::SetDataStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    IBEulerianForceSetter(
        const std::string& object_name,
        const int F_current_idx,
        const int F_new_idx,
        const int F_half_idx);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBEulerianForceSetter();

    /*!
     * \brief Set the current and new times for the present timestep.
     */
    void
    setTimeInterval(
        const double current_time,
        const double new_time);

    /*!
     * \brief Register an optional additional body force specification which
     * will be added to the IB force.
     */
    void
    registerBodyForceSpecification(
        SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> F_setter);

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
    IBEulerianForceSetter();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBEulerianForceSetter(
        const IBEulerianForceSetter& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBEulerianForceSetter&
    operator=(
        const IBEulerianForceSetter& that);

    /*!
     * The current and new time for the present timestep.
     */
    double d_current_time, d_new_time;

    /*!
     * Patch data descriptor indices for the current, new, and half-time force
     * data.
     */
    const int d_F_current_idx, d_F_new_idx, d_F_half_idx;

    /*!
     * Optional body force generator.
     */
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> d_body_force_setter;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBEulerianForceSetter.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBEulerianForceSetter
