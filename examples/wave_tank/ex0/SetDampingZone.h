// Filename: SetDampingZone.h
// Created on Sep 04, 2018 by Amneet Bhalla

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_SetDampingZone
#define included_SetDampingZone

///////////////////////////// INCLUDES ///////////////////////////////////

#include <Variable.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibtk/ibtk_utilities.h>
#include <tbox/Pointer.h>

namespace IBTK
{
class HierarchyMathOps;
}

/*!
 * A callback function to damp water waves near the outflow boundary.
 *
 * \param damping_coef_idx a patch data index for setting damping zone coefficent in the fluid solver.
 * \param ctx is the pointer to DampingZone class object.
 */

void callDampingZoneCallbackFunction(int damping_coef_idx,
                                     SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                     int cycle_num,
                                     double time,
                                     double current_time,
                                     double new_time,
                                     void* ctx);
class DampingZone
{
    /*!
     * \brief Class DampingZone is a class that computes damping zone coefficents to be
     * included in the fluid momentum equation.
     */
public:
    /*!
     * The only constructor of this class.
     */
    DampingZone(const std::string& object_name,
                SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> ins_hier_integrator,
                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * Destructor for this class.
     */
    ~DampingZone();

    /*!
     * Set the density based on the current level set information
     */
    void setDampingZoneCoefPatchData(int damping_coef_idx,
                                     SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                     const int cycle_num,
                                     const double time,
                                     const double current_time,
                                     const double new_time);

    //////////////// PRIVATE /////////////////////////////

private:
    /*!
     * Default constructor is not implemented and should not be used.
     */
    DampingZone();

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    DampingZone& operator=(const DampingZone& that);

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    DampingZone(const DampingZone& from);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to advection-diffusion solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> d_ins_hier_integrator;

    /*!
     * Parameters to control damping zone coefficients.
     */
    double d_x_zone_start, d_x_zone_end, d_theta;

}; // DampingZone

#endif // #ifndef included_DampingZone
