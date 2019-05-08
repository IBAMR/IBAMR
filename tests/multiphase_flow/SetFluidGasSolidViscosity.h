// Filename: SetFluidGasSolidViscosity.h
//
// Copyright (c) 2002-2017, Amneet Bhalla and Nishant Nangia
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
//    * Neither the name of New York University nor the names of its
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
// POSSIBILITY OF SUCH DAMAGE

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic SAMRAI objects
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibtk/CartGridFunctionSet.h>

#ifndef included_IBAMR_multiphase_flow_SetFluidGasSolidViscosity
#define included_IBAMR_multiphase_flow_SetFluidGasSolidViscosity

/*!
 * \brief Class SetFluidGasSolidViscosity is a utility class which sets the fluid and
 * solid Eulerian density based on the current level set information
 */
class SetFluidGasSolidViscosity
{
public:
    /*!
     * The only constructor of this class.
     */
    SetFluidGasSolidViscosity(const std::string& object_name,
                              SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_solid_var,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_gas_var,
                              const double mu_fluid,
                              const double mu_gas,
                              const double mu_solid,
                              const int ls_reinit_interval,
                              const double num_solid_interface_cells,
                              const double num_gas_interface_cells,
                              const bool set_mu_solid)
        : d_object_name(object_name),
          d_adv_diff_solver(adv_diff_solver),
          d_ls_solid_var(ls_solid_var),
          d_ls_gas_var(ls_gas_var),
          d_mu_fluid(mu_fluid),
          d_mu_gas(mu_gas),
          d_mu_solid(mu_solid),
          d_ls_reinit_interval(ls_reinit_interval),
          d_num_solid_interface_cells(num_solid_interface_cells),
          d_num_gas_interface_cells(num_gas_interface_cells),
          d_set_mu_solid(set_mu_solid)
    {
        // intentionally left blank
        return;
    } // SetFluidGasSolidViscosity

    /*!
     * Set the density based on the current level set information
     */
    void setViscosityPatchData(int mu_idx,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
                               SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                               const int /*cycle_num*/,
                               const double time,
                               const double current_time,
                               const double new_time)
    {
        // Get the current level set information
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        int ls_solid_idx = -1;
        int ls_gas_idx = -1;

        if (SAMRAI::tbox::MathUtilities<double>::equalEps(time, current_time))
        {
            ls_solid_idx = var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getCurrentContext());
        }
        else if (SAMRAI::tbox::MathUtilities<double>::equalEps(time, new_time))
        {
            ls_solid_idx = var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getNewContext());
        }
        else
        {
            TBOX_ERROR("This statement should not be reached");
        }

        if (SAMRAI::tbox::MathUtilities<double>::equalEps(time, current_time))
        {
            ls_gas_idx = var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_solver->getCurrentContext());
        }
        else if (SAMRAI::tbox::MathUtilities<double>::equalEps(time, new_time))
        {
            ls_gas_idx = var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_solver->getNewContext());
        }
        else
        {
            TBOX_ERROR("This statement should not be reached");
        }

        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > mu_cc_var = mu_var;
        if (mu_cc_var)
        {
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
                        patch->getPatchGeometry();
                    const double* const patch_dx = patch_geom->getDx();
                    const double alpha = d_num_solid_interface_cells * patch_dx[0];
                    const double beta = d_num_gas_interface_cells * patch_dx[1];

                    const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > ls_solid_data =
                        patch->getPatchData(ls_solid_idx);
                    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > ls_gas_data =
                        patch->getPatchData(ls_gas_idx);
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);

                    // Calderer et al, 2014
                    for (SAMRAI::hier::Box<NDIM>::Iterator it(patch_box); it; it++)
                    {
                        SAMRAI::pdat::CellIndex<NDIM> ci(it());
                        const double phi_s = (*ls_solid_data)(ci);
                        const double phi_g = (*ls_gas_data)(ci);
                        double Hphi_s, Hphi_g;
                        if (phi_s < -alpha)
                        {
                            Hphi_s = 0.0;
                        }
                        else if (std::abs(phi_s) <= alpha)
                        {
                            Hphi_s = 0.5 + 0.5 * phi_s / alpha + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_s / alpha);
                        }
                        else
                        {
                            Hphi_s = 1.0;
                        }

                        if (phi_g < -beta)
                        {
                            Hphi_g = 0.0;
                        }
                        else if (std::abs(phi_g) <= beta)
                        {
                            Hphi_g = 0.5 + 0.5 * phi_g / beta + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_g / beta);
                        }
                        else
                        {
                            Hphi_g = 1.0;
                        }

                        // First, compute the viscosity of the "flowing" phases
                        const double mu_flow = (d_mu_fluid - d_mu_gas) * Hphi_g + d_mu_gas;

                        // Next, set the viscosity of the solid phase in the usual way
                        if (d_set_mu_solid)
                            (*mu_data)(ci) = (mu_flow - d_mu_solid) * Hphi_s + d_mu_solid;
                        else
                            (*mu_data)(ci) = mu_flow;
                    }
                }
            }
        }
        else
        {
            // Erroring out if any other centered is used for mu
            TBOX_ERROR("This statement should not have been reached");
        }

        return;
    } // setViscosityPatchData

    //////////////// PRIVATE /////////////////////////////
private:
    SetFluidGasSolidViscosity& operator=(const SetFluidGasSolidViscosity&) = delete;
    SetFluidGasSolidViscosity(const SetFluidGasSolidViscosity&) = delete;

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * SAMRAI::tbox::Pointer to advection-diffusion solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Level set variables
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_solid_var, d_ls_gas_var;

    /*!
     * Density of the fluid and solid.
     */
    double d_mu_fluid, d_mu_gas, d_mu_solid;

    /*!
     * Level set reinitialization interval
     */
    int d_ls_reinit_interval;

    /*!
     * Number of cells over which to transition between values
     */
    double d_num_solid_interface_cells, d_num_gas_interface_cells;

    /*!
     * Whether or not to set viscosity in the solid region
     */
    bool d_set_mu_solid;

}; // SetFluidGasSolidViscosity

inline void
callSetFluidGasSolidViscosityCallbackFunction(int mu_idx,
                                              SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
                                              SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                              const int cycle_num,
                                              const double time,
                                              const double current_time,
                                              const double new_time,
                                              void* ctx)
{
    // Set the density from the level set information
    static SetFluidGasSolidViscosity* ptr_SetFluidGasSolidViscosity = static_cast<SetFluidGasSolidViscosity*>(ctx);
    ptr_SetFluidGasSolidViscosity->setViscosityPatchData(
        mu_idx, mu_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;

} // callSetFluidGasSolidViscosityCallBackFunction

#endif
