// Filename: WaveGenerationStrategy.h
// Created on 15 Oct 2019 by Amneet Bhalla
//
// Copyright (c) 2002-2019, Amneet Bhalla
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
//    * Neither the name of The University of North Carolina nor the names of its
//      contributors may be used to endorse or promote products derived from
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

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_WaveGenerationStrategy
#define included_WaveGenerationStrategy

///////////////////////////// INCLUDES ///////////////////////////////////
#include <string>

#include "Variable.h"
#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibtk/ibtk_utilities.h"
#include "tbox/Pointer.h"

namespace IBAMR
{
/*!
 * A collection of post processing call back functions to be hooked into
 * IBAMR::HierarchyIntegrator class to employ a wave generation zone of some
 * prescribed length within the computational domain. Additional strategies
 * can be implemented as simple C functions within this collection.
 *
 * The instance of a WaveGenerationStrategy and the particular generation function
 * should be registered with an appropriate hierarchy integrator via
 * registerPostprocessIntegrateHierarchyCallback.
 *
 * \param ctx is the pointer to WaveGenerationStrategy struct.
 */

namespace WaveGenerationFunctions
{
void callStokesWaveRelaxationCallbackFunction(double current_time,
                                              double new_time,
                                              bool skip_synchronize_new_state_data,
                                              int num_cycles,
                                              void* ctx);
} // namespace WaveGenerationFunctions

/*! Struct for generating waves based on relaxation method.
 *
 */
struct WaveGenerationStrategy
{
    /*
     * Pointers to the fluid and advection-diffusion integrators.
     */
    SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> d_ins_hier_integrator;
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_hier_integrator;

    /*
     * Pointer to the level set variable representing the wave interface.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_phi_var;

    /*
     * Start and end coordinates of the damping zone, and the damping coefficient.
     */
    double d_x_zone_start, d_x_zone_end;
    double d_alpha;
    int d_sign_gas_phase = 1;

    /*!
     * Number of interface cells to represent air-water thickness.
     */
    double d_num_interface_cells;
};

/*!
 * Base class for generating Stokes wave.
 */
class StokesWaveGenerator : public WaveGenerationStrategy
{
public:
    /*!
     * \brief Constructor of the class.
     */
    StokesWaveGenerator(const std::string& object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~StokesWaveGenerator() = default;

    /*!
     * Get surface elevation at a specified horizontal position and time.
     */
    virtual double getSurfaceElevation(double x, double time) const = 0;

    /*!
     * Get velocity component at a specified position and time.
     */
    virtual double getVelocity(double x, double z_plus_d, double time, int comp_idx) const = 0;

    /*
     * \brief Get Stokes wave parameters.
     */
    double getWaterDepth() const
    {
        return d_depth;
    } // getWaterDept

    double getWaveAngularFrequency() const
    {
        return d_omega;
    } // getWaveAngularFrequency

    double getWaveNumber() const
    {
        return d_wave_number;
    } // getWaveNumber

    double getWaveAmplitude() const
    {
        return d_amplitude;
    } // getWaveAmplitude

    double getGravity() const
    {
        return d_gravity;
    } // getGravity

protected:
    /*!
     * Book-keeping.
     */
    std::string d_object_name;

    /*!
     * \brief Wave parameters.
     *
     * \param d_wave_number  : Wave number of dominant wave component [$2\pi/m$]
     * \param d_amplitude    : Amplitude of the dominant wave component [m]
     * \param d_depth        : Depth of water, from sea bed to still water level [m]
     * \param d_gravity      : Acceleration due to gravity [$m/s^2$]
     * \param d_omega        : Angular frequency [$2 \pi/s$] (optional)
     *
     * \NOTE The class calculates a more accurate value of omega from the expansion coefficients
     * and the provided value in not used.
     */
    double d_depth, d_omega, d_wave_number, d_amplitude, d_gravity;

    /*!
     * If we are calculating in deep water limit.
     */
    bool d_deep_water_limit = false;

    /*
     *  \brief Stokes wave generator data members.
     */
    /*!
     * Get wave parameters from input db.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

}; // StokesWaveGenerator

class FifthOrderStokesWaveGenerator : public StokesWaveGenerator
{
public:
    FifthOrderStokesWaveGenerator(const std::string& object_name,
                                  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * Get surface elevation at a specified horizontal position and time.
     */
    double getSurfaceElevation(double x, double time) const;

    /*!
     * Get velocity component at a specified position and time.
     */
    double getVelocity(double x, double z_plus_d, double time, int comp_idx) const;

private:
    /*!
     * Initialize Stokes coefficients.
     */
    void initStokesCoefficients();

    /*!
     * Stokes coefficients.
     */
    double d_A[6][6], d_B[6][6], d_C[5];
    double d_p[5], d_eta[5];
};

} // namespace IBAMR

#endif // #ifndef included_WaveGenerationStrategy
