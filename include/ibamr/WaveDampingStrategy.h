// Filename: WaveDampingStrategy.h
// Created on 30 Jan 2019 by Nishant Nangia and Amneet Bhalla
//
// Copyright (c) 2002-2019, Nishant Nangia and Amneet Bhalla
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

#ifndef included_WaveDampingStrategy
#define included_WaveDampingStrategy

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
 * IBAMR::HierarchyIntegrator class to employ a wave damping zone of some
 * prescribed length within the computational domain. Additional strategies
 * can be implemented as simple C functions within this collection.
 *
 * The instance of a WaveDampingStrategy and the particular damping function
 * should be registered with an appropriate hierarchy integrator via
 * registerPostprocessIntegrateHierarchyCallback.
 *
 * \param ctx is the pointer to WaveDampingStrategy struct.
 */

void callRelaxationZoneCallbackFunction(double current_time,
                                        double new_time,
                                        bool skip_synchronize_new_state_data,
                                        int num_cycles,
                                        void* ctx);

void callConservedWaveAbsorbingCallbackFunction(double current_time,
                                                double new_time,
                                                bool skip_synchronize_new_state_data,
                                                int num_cycles,
                                                void* ctx);

/*!
 * A struct holding the required information used by the wave damping function.
 */
struct WaveDampingStrategy
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
     * Start and end coordinates of the damping zone, water depth, and damping coefficient.
     */
    double d_x_zone_start, d_x_zone_end, d_depth;
    double d_alpha;
};

/*!
 * A struct holding the required information used by the variable alpha damping method
 */
struct MassConservationFunctor
{
    std::pair<double, double> operator()(const double& alpha);
    double d_dt;
    SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > d_patch_hierarchy;
    int d_u_idx;
    int d_phi_scratch_idx;
    int d_I_idx;
    int d_dI_idx;
    WaveDampingStrategy* d_ptr_wave_damper;
};
} // namespace IBAMR

#endif // #ifndef included_WaveDampingStrategy
