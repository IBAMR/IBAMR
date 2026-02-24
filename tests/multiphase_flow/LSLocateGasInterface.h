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

// Config files

#include "ibtk/samrai_compatibility_names.h"

#include <SAMRAI_config.h>

// Headers for basic SAMRAI objects
#include "SAMRAIBox.h"
#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAICartesianPatchGeometry.h"
#include "SAMRAICellData.h"
#include "SAMRAICellIndex.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIHierarchyCellDataOpsReal.h"
#include "SAMRAIIndex.h"
#include "SAMRAILoadBalancer.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIStandardTagAndInitialize.h"
#include "SAMRAIVariableDatabase.h"

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
                         SAMRAIPointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                         SAMRAIPointer<SAMRAICellVariable<double>> ls_var,
                         const double init_height)
        : d_object_name(object_name), d_adv_diff_solver(adv_diff_solver), d_ls_var(ls_var), d_init_height(init_height)
    {
        return;
    } // LSLocateGasInterface

    /*!
     * Reinitialize the level set information
     */
    void setLevelSetPatchData(int D_idx,
                              SAMRAIPointer<IBTK::HierarchyMathOps> hier_math_ops,
                              const double /*time*/,
                              const bool initial_time)
    {
        SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy = hier_math_ops->getPatchHierarchy();
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();

        // If not the intial time, set the level set to the current value maintained by the integrator
        if (!initial_time)
        {
            SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
            const int ls_current_idx =
                var_db->mapVariableAndContextToIndex(d_ls_var, d_adv_diff_solver->getCurrentContext());
            SAMRAIHierarchyCellDataOpsReal<double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

            hier_cc_data_ops.copyData(D_idx, ls_current_idx);

            return;
        }

        // Set the initial condition for locating the interface
        const double H = d_init_height;

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                const SAMRAIBox& patch_box = patch->getBox();
                SAMRAIPointer<SAMRAICellData<double>> D_data = patch->getPatchData(D_idx);
                for (SAMRAIBox::Iterator it(patch_box); it; it++)
                {
                    SAMRAICellIndex ci(it());

                    // Get physical coordinates
                    IBTK::Vector coord = IBTK::Vector::Zero();
                    SAMRAIPointer<SAMRAICartesianPatchGeometry> patch_geom = patch->getPatchGeometry();
                    const double* patch_X_lower = patch_geom->getXLower();
                    const SAMRAIIndex& patch_lower_idx = patch_box.lower();
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
    SAMRAIPointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Level set variable.
     */
    SAMRAIPointer<SAMRAICellVariable<double>> d_ls_var;

    /*!
     * Initial level set information.
     */
    double d_init_height;
};

inline void
callLSLocateGasInterfaceCallbackFunction(int D_idx,
                                         SAMRAIPointer<IBTK::HierarchyMathOps> hier_math_ops,
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
