// Filename: LSLocateGasInterface.h
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

#ifndef included_IBAMR_multiphase_flow_LSLocateGasInterface
#define included_IBAMR_multiphase_flow_LSLocateGasInterface

/*!
 * \brief class LSLocateGasInterface is a utility class which is used to
 * identify the circular interface for level set computations
 */
class LSLocateGasInterface
{
public:
    LSLocateGasInterface(const std::string& object_name,
                         SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                         const double init_height)
        : d_object_name(object_name), d_adv_diff_solver(adv_diff_solver), d_ls_var(ls_var), d_init_height(init_height)
    {
        return;
    } // LSLocateGasInterface

    /*!
     * Reinitialize the level set information
     */
    void setLevelSetPatchData(int D_idx,
                              SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                              const double /*time*/,
                              const bool initial_time)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();

        // If not the intial time, set the level set to the current value maintained by the integrator
        if (!initial_time)
        {
            SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
            const int ls_current_idx =
                var_db->mapVariableAndContextToIndex(d_ls_var, d_adv_diff_solver->getCurrentContext());
            SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(
                patch_hierarchy, coarsest_ln, finest_ln);

            hier_cc_data_ops.copyData(D_idx, ls_current_idx);

            return;
        }

        // Set the initial condition for locating the interface
        const double H = d_init_height;

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > D_data = patch->getPatchData(D_idx);
                for (SAMRAI::hier::Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    SAMRAI::pdat::CellIndex<NDIM> ci(it());

                    // Get physical coordinates
                    IBTK::Vector coord = IBTK::Vector::Zero();
                    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
                        patch->getPatchGeometry();
                    const double* patch_X_lower = patch_geom->getXLower();
                    const SAMRAI::hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
                    const double* const patch_dx = patch_geom->getDx();
                    for (int d = 0; d < NDIM; ++d)
                    {
                        coord[d] =
                            patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                    }

                    const double distance = NDIM < 3 ? coord[1] - H : coord[2] - H;

                    // Initialize the locator data to be zero on the interface,
                    // negative inside, and positive outside
                    (*D_data)(ci) = -distance;
                }
            }
        }
        return;
    } // setLevelSetPatchData

    //////////////// PRIVATE /////////////////////////////

private:
    LSLocateGasInterface& operator=(const LSLocateGasInterface&) = delete;
    LSLocateGasInterface(const LSLocateGasInterface&) = delete;

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * SAMRAI::tbox::Pointer to the advection-diffusion solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Level set variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_var;

    /*!
     * Initial level set information.
     */
    double d_init_height;
};

inline void
callLSLocateGasInterfaceCallbackFunction(int D_idx,
                                         SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                         double time,
                                         bool initial_time,
                                         void* ctx)
{
    // Set the level set information
    static LSLocateGasInterface* ptr_LSLocateGasInterface = static_cast<LSLocateGasInterface*>(ctx);
    ptr_LSLocateGasInterface->setLevelSetPatchData(D_idx, hier_math_ops, time, initial_time);

    return;
} // callLSLocateGasInterfaceCallbackFunction

#endif
