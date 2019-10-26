// Filename: SetFluidGasSolidDensity.h
// Created on Nov 15, 2017 by Nishant Nangia
//
// Copyright (c) 2002-2018, Nishant Nangia and Amneet Bhalla
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

#ifndef included_SetFluidGasSolidDensity
#define included_SetFluidGasSolidDensity

///////////////////////////// INCLUDES ///////////////////////////////////

#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <Variable.h>

namespace IBTK
{
class HierarchyMathOps;
}

/*!
 * Pre processing call back function to be hooked into IBAMR::VCINSStaggeredHierarchyIntegratorclass.
 *
 * \param rho_idx a patch data index for the current density variable maintained by the integrator.
 * \param ctx is the pointer to SetFluidGasSolidDensity class object.
 */

void callSetFluidGasSolidDensityCallbackFunction(int rho_idx,
                                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
                                                 SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                 const int cycle_num,
                                                 const double time,
                                                 const double current_time,
                                                 const double new_time,
                                                 void* ctx);

class SetFluidGasSolidDensity
{
    /*!
     * \brief Class SetFluidGasSolidDensity is a utility class which sets the fluid and
     * solid Eulerian density based on the current level set information
     */
public:
    /*!
     * The only constructor of this class.
     */
    SetFluidGasSolidDensity(const std::string& object_name,
                            SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_solid_var,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_gas_var,
                            const double rho_fluid,
                            const double rho_gas,
                            const double rho_solid,
                            const int ls_reinit_interval,
                            const double num_solid_interface_cells,
                            const double num_gas_interface_cells);

    /*!
     * Destructor for this class.
     */
    ~SetFluidGasSolidDensity();

    /*!
     * Set the density based on the current level set information
     */
    void setDensityPatchData(int rho_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
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
    SetFluidGasSolidDensity();

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    SetFluidGasSolidDensity& operator=(const SetFluidGasSolidDensity& that);

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    SetFluidGasSolidDensity(const SetFluidGasSolidDensity& from);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to advection-diffusion solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Level set variables
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_solid_var, d_ls_gas_var;

    /*!
     * Density of the fluid and solid.
     */
    double d_rho_fluid, d_rho_gas, d_rho_solid;

    /*!
     * Level set reinitialization interval
     */
    int d_ls_reinit_interval;

    /*!
     * Number of interface cells
     */
    double d_num_solid_interface_cells, d_num_gas_interface_cells;

}; // SetFluidGasSolidDensity

#endif // #ifndef included_SetFluidGasSolidDensity
