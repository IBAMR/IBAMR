// Filename: LSLocateCircularInterface.h
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
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>

#include <ibtk/HierarchyMathOps.h>

#include <cmath>

#ifndef included_IBAMR_multiphase_flow_LSLocateCircularInterface
#define included_IBAMR_multiphase_flow_LSLocateCircularInterface

// Struct to maintain the properties of the circular interface
struct CircularInterface
{
    IBTK::Vector X0;
    double R;
};

/*!
 * \brief class LSLocateCircularInterface is a utility class which is used to identify
 * the circular interface for level set computations
 */
class LSLocateCircularInterface
{
public:
    /*!
     * The only constructor of this class.
     */
    LSLocateCircularInterface(const std::string& object_name,
                              SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                              CircularInterface* init_circle)
        : d_object_name(object_name), d_adv_diff_solver(adv_diff_solver), d_ls_var(ls_var), d_circle(init_circle)
    {
        return;
    } // LSLocateCircularInterface

    /*!
     * Reinitialize the level set information
     */
    void setLevelSetPatchData(int D_idx,
                              SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                              const double /*time*/,
                              const bool /*initial_time*/)
    {
        // In this version of this class, the initial level set location is set to be
        // exact since we always know the radius of the ball

        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();

        // Set the initial condition for locating the interface
        const double& R = d_circle->R;
        const IBTK::Vector& X0 = d_circle->X0;

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

                    const double distance =
                        std::sqrt(std::pow((coord[0] - X0(0)), 2.0) + std::pow((coord[1] - X0(1)), 2.0)
#if (NDIM == 3)
                                  + std::pow((coord[2] - X0(2)), 2.0)
#endif
                        );

                    (*D_data)(ci) = distance - R;
                }
            }
        }
        return;
    } // setLevelSetPatchData

    //////////////// PRIVATE /////////////////////////////

private:
    LSLocateCircularInterface& operator=(const LSLocateCircularInterface&) = delete;

    LSLocateCircularInterface(const LSLocateCircularInterface&) = delete;

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * SAMRAI::tbox::Pointer to the advection-diffusion solver
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Level set variable
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_var;

    /*!
     * Initial level set information.
     */
    CircularInterface* d_circle;
};

inline void
callLSLocateCircularInterfaceCallbackFunction(int D_idx,
                                              SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                              double time,
                                              bool initial_time,
                                              void* ctx)
{
    // Set the level set information
    static LSLocateCircularInterface* ptr_LSLocateCircularInterface = static_cast<LSLocateCircularInterface*>(ctx);
    ptr_LSLocateCircularInterface->setLevelSetPatchData(D_idx, hier_math_ops, time, initial_time);

    return;
} // callLSLocateCircularInterfaceCallbackFunction

#endif
