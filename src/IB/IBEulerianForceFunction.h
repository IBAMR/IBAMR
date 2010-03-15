#ifndef included_IBEulerianForceFunction
#define included_IBEulerianForceFunction

// Filename: IBEulerianForceFunction.h
// Last modified: <15.Mar.2010 00:14:23 griffith@griffith-macbook-pro.local>
// Created on 28 Sep 2004 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/CartGridFunction.h>

// SAMRAI INCLUDES
#include <Patch.h>
#include <Variable.h>
#include <tbox/Pointer.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBEulerianForceFunction is used to communicate the Eulerian body
 * force computed by class IBHierarchyIntegrator to the incompressible
 * Navier-Stokes solver.
 */
class IBEulerianForceFunction
    : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    IBEulerianForceFunction(
        const std::string& object_name,
        const int F_current_idx,
        const int F_new_idx,
        const int F_half_idx);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBEulerianForceFunction();

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
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> F_fcn);

    /*!
     * \name Methods to set the data.
     */
    //\{

    /*!
     * \note This concrete IBTK::CartGridFunction is time-dependent.
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
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
        const double data_time,
        const bool initial_time=false,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level=SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL));

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBEulerianForceFunction();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBEulerianForceFunction(
        const IBEulerianForceFunction& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBEulerianForceFunction&
    operator=(
        const IBEulerianForceFunction& that);

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
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_body_force_fcn;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBEulerianForceFunction.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBEulerianForceFunction
