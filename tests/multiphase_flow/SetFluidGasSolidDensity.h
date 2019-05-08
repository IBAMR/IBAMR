// Filename: SetFluidGasSolidDensity.h
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

#ifndef included_IBAMR_multiphase_flow_SetFluidGasSolidDensity
#define included_IBAMR_multiphase_flow_SetFluidGasSolidDensity

// Various options to setting side-centered densities
#define SMOOTH_SC_RHO 0
#define DESJARDINS_SC_RHO 0
#define HARMONIC_CC_TO_SC_RHO 1

/*!
 * \brief Class SetFluidGasSolidDensity is a utility class which sets the fluid and
 * solid Eulerian density based on the current level set information
 */
class SetFluidGasSolidDensity
{
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
                            const double num_gas_interface_cells)
        : d_object_name(object_name),
          d_adv_diff_solver(adv_diff_solver),
          d_ls_solid_var(ls_solid_var),
          d_ls_gas_var(ls_gas_var),
          d_rho_fluid(rho_fluid),
          d_rho_gas(rho_gas),
          d_rho_solid(rho_solid),
          d_ls_reinit_interval(ls_reinit_interval),
          d_num_solid_interface_cells(num_solid_interface_cells),
          d_num_gas_interface_cells(num_gas_interface_cells)
    {
        // intentionally left blank
        return;
    } // SetFluidGasSolidDensity

    /*!
     * Set the density based on the current level set information
     */
    void setDensityPatchData(int rho_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
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

        // Set the density based on the cell centered level set
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > rho_cc_var = rho_var;

        if (rho_cc_var)
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
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > rho_data =
                        patch->getPatchData(rho_idx);

                    // Li et al, 2015
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

                        // First, compute the density of the "flowing" phases
                        const double rho_flow = (d_rho_fluid - d_rho_gas) * Hphi_g + d_rho_gas;

                        // Next, set the density of the solid phase
                        (*rho_data)(ci) = (rho_flow - d_rho_solid) * Hphi_s + d_rho_solid;
                    }
                }
            }
        }

        // Setting side centered density directly
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > rho_sc_var = rho_var;
        if (rho_sc_var)
        {
            // Note, this method requires ghost cells to be filled for the level set variable
            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* ls_solid_bc_coef =
                d_adv_diff_solver->getPhysicalBcCoefs(d_ls_solid_var).front();
            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* ls_gas_bc_coef =
                d_adv_diff_solver->getPhysicalBcCoefs(d_ls_gas_var).front();
            SAMRAI::hier::IntVector<NDIM> cell_ghosts = 1;
            const int ls_solid_scratch_idx = var_db->registerVariableAndContext(
                d_ls_solid_var, var_db->getContext(d_object_name + "::SOLID::SCRATCH"), cell_ghosts);
            const int ls_gas_scratch_idx = var_db->registerVariableAndContext(
                d_ls_gas_var, var_db->getContext(d_object_name + "::GAS::SCRATCH"), cell_ghosts);
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                patch_hierarchy->getPatchLevel(ln)->allocatePatchData(ls_solid_scratch_idx, time);
                patch_hierarchy->getPatchLevel(ln)->allocatePatchData(ls_gas_scratch_idx, time);
            }
            typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent
                InterpolationTransactionComponent;
            std::vector<InterpolationTransactionComponent> ls_transaction_comps(2);
            ls_transaction_comps[0] = InterpolationTransactionComponent(ls_solid_scratch_idx,
                                                                        ls_solid_idx,
                                                                        "CONSERVATIVE_LINEAR_REFINE",
                                                                        false,
                                                                        "CONSERVATIVE_COARSEN",
                                                                        "LINEAR",
                                                                        false,
                                                                        ls_solid_bc_coef);
            ls_transaction_comps[1] = InterpolationTransactionComponent(ls_gas_scratch_idx,
                                                                        ls_gas_idx,
                                                                        "CONSERVATIVE_LINEAR_REFINE",
                                                                        false,
                                                                        "CONSERVATIVE_COARSEN",
                                                                        "LINEAR",
                                                                        false,
                                                                        ls_gas_bc_coef);
            SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> hier_bdry_fill =
                new IBTK::HierarchyGhostCellInterpolation();
            hier_bdry_fill->initializeOperatorState(ls_transaction_comps, patch_hierarchy);
            hier_bdry_fill->fillData(time);

            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
                        patch->getPatchGeometry();

                    const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > ls_solid_data =
                        patch->getPatchData(ls_solid_scratch_idx);
                    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > ls_gas_data =
                        patch->getPatchData(ls_gas_scratch_idx);
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > rho_data =
                        patch->getPatchData(rho_idx);

                    // Compute the indicators for both level sets
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (SAMRAI::hier::Box<NDIM>::Iterator it(
                                 SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(patch_box, axis));
                             it;
                             it++)
                        {
                            SAMRAI::pdat::SideIndex<NDIM> si(it(), axis, 0);
                            const double phi_solid_lower = (*ls_solid_data)(si.toCell(0));
                            const double phi_solid_upper = (*ls_solid_data)(si.toCell(1));
                            const double phi_gas_lower = (*ls_gas_data)(si.toCell(0));
                            const double phi_gas_upper = (*ls_gas_data)(si.toCell(1));

                            if (DESJARDINS_SC_RHO)
                            {
                                // SETTING 1: Desjardins way to set side-centered density
                                double h_solid, h_gas;
                                if (phi_solid_lower >= 0.0 && phi_solid_upper >= 0.0)
                                {
                                    h_solid = 1.0;
                                }
                                else if (phi_solid_lower <= 0.0 && phi_solid_upper <= 0.0)
                                {
                                    h_solid = 0.0;
                                }
                                else
                                {
                                    h_solid = (std::max(phi_solid_lower, 0.0) + std::max(phi_solid_upper, 0.0)) /
                                              (std::abs(phi_solid_lower) + std::abs(phi_solid_upper));
                                }
                                if (phi_gas_lower >= 0.0 && phi_gas_upper >= 0.0)
                                {
                                    h_gas = 1.0;
                                }
                                else if (phi_gas_lower <= 0.0 && phi_gas_upper <= 0.0)
                                {
                                    h_gas = 0.0;
                                }
                                else
                                {
                                    h_gas = (std::max(phi_gas_lower, 0.0) + std::max(phi_gas_upper, 0.0)) /
                                            (std::abs(phi_gas_lower) + std::abs(phi_gas_upper));
                                }

                                // First, compute the density of the "flowing" phases
                                const double rho_flow = d_rho_gas + (d_rho_fluid - d_rho_gas) * h_gas;

                                // Next, compute the density taking into account the solid phase
                                (*rho_data)(si) = d_rho_solid + (rho_flow - d_rho_solid) * h_solid;
                            }
                            else if (HARMONIC_CC_TO_SC_RHO)
                            {
                                // SETTING 2: Set rho on cell centers and harmonic average to side centers
                                const double* const patch_dx = patch_geom->getDx();
                                const double alpha = d_num_solid_interface_cells * patch_dx[0];
                                const double beta = d_num_gas_interface_cells * patch_dx[1];
                                double h_solid_lower, h_solid_upper, h_gas_lower, h_gas_upper;

                                if (phi_solid_lower < -alpha)
                                {
                                    h_solid_lower = 0.0;
                                }
                                else if (std::abs(phi_solid_lower) <= alpha)
                                {
                                    h_solid_lower = 0.5 + 0.5 * phi_solid_lower / alpha +
                                                    1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_solid_lower / alpha);
                                }
                                else
                                {
                                    h_solid_lower = 1.0;
                                }
                                if (phi_solid_upper < -alpha)
                                {
                                    h_solid_upper = 0.0;
                                }
                                else if (std::abs(phi_solid_upper) <= alpha)
                                {
                                    h_solid_upper = 0.5 + 0.5 * phi_solid_upper / alpha +
                                                    1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_solid_upper / alpha);
                                }
                                else
                                {
                                    h_solid_upper = 1.0;
                                }

                                if (phi_gas_lower < -beta)
                                {
                                    h_gas_lower = 0.0;
                                }
                                else if (std::abs(phi_gas_lower) <= beta)
                                {
                                    h_gas_lower = 0.5 + 0.5 * phi_gas_lower / beta +
                                                  1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_gas_lower / beta);
                                }
                                else
                                {
                                    h_gas_lower = 1.0;
                                }
                                if (phi_gas_upper < -beta)
                                {
                                    h_gas_upper = 0.0;
                                }
                                else if (std::abs(phi_gas_upper) <= beta)
                                {
                                    h_gas_upper = 0.5 + 0.5 * phi_gas_upper / beta +
                                                  1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_gas_upper / beta);
                                }
                                else
                                {
                                    h_gas_upper = 1.0;
                                }
                                const double rho_flow_lower = (d_rho_fluid - d_rho_gas) * h_gas_lower + d_rho_gas;
                                const double rho_flow_upper = (d_rho_fluid - d_rho_gas) * h_gas_upper + d_rho_gas;

                                const double rho_full_lower =
                                    (rho_flow_lower - d_rho_solid) * h_solid_lower + d_rho_solid;
                                const double rho_full_upper =
                                    (rho_flow_upper - d_rho_solid) * h_solid_upper + d_rho_solid;

                                (*rho_data)(si) =
                                    2.0 * rho_full_upper * rho_full_lower / (rho_full_upper + rho_full_lower);
                            }
                            else if (SMOOTH_SC_RHO)
                            {
                                // SETTING 3: Simple average of phi onto side centers and set rho_sc directly
                                double h_solid, h_gas;
                                const double* const patch_dx = patch_geom->getDx();
                                const double alpha = d_num_solid_interface_cells * patch_dx[0];
                                const double beta = d_num_gas_interface_cells * patch_dx[1];
                                const double phi_solid = 0.5 * (phi_solid_lower + phi_solid_upper);
                                const double phi_gas = 0.5 * (phi_gas_lower + phi_gas_upper);

                                if (phi_solid < -alpha)
                                {
                                    h_solid = 0.0;
                                }
                                else if (std::abs(phi_solid) <= alpha)
                                {
                                    h_solid = 0.5 + 0.5 * phi_solid / alpha +
                                              1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_solid / alpha);
                                }
                                else
                                {
                                    h_solid = 1.0;
                                }

                                if (phi_gas < -beta)
                                {
                                    h_gas = 0.0;
                                }
                                else if (std::abs(phi_gas) <= beta)
                                {
                                    h_gas = 0.5 + 0.5 * phi_gas / beta +
                                            1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_gas / beta);
                                }
                                else
                                {
                                    h_gas = 1.0;
                                }
                                const double rho_flow = (d_rho_fluid - d_rho_gas) * h_gas + d_rho_gas;
                                const double rho_full = (rho_flow - d_rho_solid) * h_solid + d_rho_solid;

                                (*rho_data)(si) = rho_full;
                            }
                            else
                            {
                                TBOX_ERROR("No side centered density setting was chosen");
                            }
                        }
                    }
                }
            }

            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(ls_solid_scratch_idx);
                patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(ls_gas_scratch_idx);
            }
            var_db->removePatchDataIndex(ls_solid_scratch_idx);
            var_db->removePatchDataIndex(ls_gas_scratch_idx);
        }

        return;
    } // setDensityPatchData

    //////////////// PRIVATE /////////////////////////////

private:
    SetFluidGasSolidDensity& operator=(const SetFluidGasSolidDensity&) = delete;

    SetFluidGasSolidDensity(const SetFluidGasSolidDensity&) = delete;

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

inline void
callSetFluidGasSolidDensityCallbackFunction(int rho_idx,
                                            SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
                                            SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                            const int cycle_num,
                                            const double time,
                                            const double current_time,
                                            const double new_time,
                                            void* ctx)
{
    // Set the density from the level set information
    static SetFluidGasSolidDensity* ptr_SetFluidGasSolidDensity = static_cast<SetFluidGasSolidDensity*>(ctx);
    ptr_SetFluidGasSolidDensity->setDensityPatchData(
        rho_idx, rho_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;

} // callSetFluidGasSolidDensityCallBackFunction

#endif
