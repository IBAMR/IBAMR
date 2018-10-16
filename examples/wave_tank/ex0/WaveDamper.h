// Filename: WaveDamper.h
// Created on Sep 04, 2018 by Amneet Bhalla and Nishant Nangia

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_WaveDamper
#define included_WaveDamper

///////////////////////////// INCLUDES ///////////////////////////////////

#include <Variable.h>
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibtk/ibtk_utilities.h>
#include <tbox/Pointer.h>

namespace IBTK
{
class HierarchyMathOps;
}

/*!
 * A callback function to damp water waves in some relaxation zone at the end of the wave tank.
 *
 * \param damping_coef_idx a patch data index for setting damping zone coefficent in the fluid solver.
 * \param ctx is the pointer to WaveDamper class object.
 */

void callWaveDamperCallbackFunction(double current_time,
                                    double new_time,
                                    bool skip_synchronize_new_state_data,
                                    int num_cycles,
                                    void* ctx);
class WaveDamper
{
    /*!
     * \brief Class WaveDamper is a class that damps water waves at the end of the wave tank.
     */
public:
    /*!
     * The only constructor of this class.
     */
    WaveDamper(const std::string& object_name,
               SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> ins_hier_integrator,
               SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_hier_integrator,
               SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > phi_var,
               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * Destructor for this class.
     */
    ~WaveDamper();

    /*!
     * Relax the outlet waves
     */
    void DampOutletWaves(double current_time, double new_time, bool skip_synchronize_new_state_data, int num_cycles);

    //////////////// PRIVATE /////////////////////////////

private:
    /*!
     * Default constructor is not implemented and should not be used.
     */
    WaveDamper();

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    WaveDamper& operator=(const WaveDamper& that);

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    WaveDamper(const WaveDamper& from);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to the INS solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> d_ins_hier_integrator;

    /*!
     * Pointer to advection-diffusion solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_hier_integrator;

    /*!
     * Pointer to the level set variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_phi_var;

    /*!
     * Parameters to control damping zone coefficients.
     */
    double d_x_zone_start, d_x_zone_end, d_depth;

    /*
     * Value of alpha.
     */
    double d_alpha;

}; // WaveDamper

#endif // #ifndef included_WaveDamper
