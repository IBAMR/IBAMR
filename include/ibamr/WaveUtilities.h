// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_WaveUtilities
#define included_WaveUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

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
