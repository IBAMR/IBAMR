// Filename: WaveUtilities.h
// Created on 4 Nov 2019 by Amneet Bhalla
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

#ifndef included_WaveUtilities
#define included_WaveUtilities

///////////////////////////// INCLUDES ///////////////////////////////////
#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"

#include "Variable.h"
#include "VariableContext.h"
#include "tbox/Pointer.h"

#include <string>

namespace IBAMR
{
/*!
 * Struct for generating waves at channel inlet based on relaxation method.
 */
struct WaveGenerationData
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
     *  Pointer to phi variable's new context.
     *  \note We modify the phi value after the end of each timestep.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_phi_new_ctx;

    /*
     * Start and end coordinates of the generation zone, and the damping coefficient.
     */
    double d_x_zone_start, d_x_zone_end;
    double d_alpha;
    int d_sign_gas_phase = 1;

    /*!
     * Number of interface cells to represent air-water thickness.
     */
    double d_num_interface_cells;

}; // WaveGenerationData

/*!
 * Struct for damping waves at channel outlet based on relaxation method.
 */
struct WaveDampingData
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
     *  Pointer to phi variable's new context.
     *  \note We modify the phi value after the end of each timestep.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_phi_new_ctx;

    /*
     * Start and end coordinates of the damping zone, water depth, and damping coefficient.
     */
    double d_x_zone_start, d_x_zone_end, d_depth;
    double d_alpha;
    int d_sign_gas_phase = 1;

}; // WaveDampingData

} // namespace IBAMR

#endif // #ifndef included_WaveUtilities
