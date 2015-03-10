// Filename: AdvDiffPredictorCorrectorHyperbolicPatchOps.cpp
// Created on 19 Mar 2004 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <map>
#include <ostream>
#include <set>
#include <string>

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/algs/HyperbolicLevelIntegrator.h"
#include "IBAMR_config.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "ibamr/AdvDiffPredictorCorrectorHyperbolicPatchOps.h"
#include "ibamr/AdvectorExplicitPredictorPatchOps.h"
#include "ibamr/AdvectorPredictorCorrectorHyperbolicPatchOps.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartGridFunction.h"
#include "SAMRAI/tbox/Database.h"

#include "SAMRAI/tbox/Utilities.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define ADV_DIFF_CONSDIFF_FC IBAMR_FC_FUNC_(adv_diff_consdiff2d, ADV_DIFF_CONSDIFF2D)
#define ADV_DIFF_CONSDIFFWITHDIVSOURCE_FC                                                                              \
    IBAMR_FC_FUNC_(adv_diff_consdiffwithdivsource2d, ADV_DIFF_CONSDIFFWITHDIVSOURCE2D)
#endif

#if (NDIM == 3)
#define ADV_DIFF_CONSDIFF_FC IBAMR_FC_FUNC_(adv_diff_consdiff3d, ADV_DIFF_CONSDIFF3D)
#define ADV_DIFF_CONSDIFFWITHDIVSOURCE_FC                                                                              \
    IBAMR_FC_FUNC_(adv_diff_consdiffwithdivsource3d, ADV_DIFF_CONSDIFFWITHDIVSOURCE3D)
#endif

extern "C" {
void ADV_DIFF_CONSDIFF_FC(const double*,
                          const double&,
#if (NDIM == 2)
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const double*,
                          const double*,
#endif
#if (NDIM == 3)
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const double*,
                          const double*,
                          const double*,
#endif
                          double*);

void ADV_DIFF_CONSDIFFWITHDIVSOURCE_FC(const double*,
                                       const double&,
#if (NDIM == 2)
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const double*,
                                       const double*,
                                       const double*,
                                       const double*,
                                       const double*,
                                       const double*,
#endif
#if (NDIM == 3)
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const double*,
                                       const double*,
                                       const double*,
                                       const double*,
                                       const double*,
                                       const double*,
                                       const double*,
                                       const double*,
                                       const double*,
#endif
                                       double*);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffPredictorCorrectorHyperbolicPatchOps::AdvDiffPredictorCorrectorHyperbolicPatchOps(
    const std::string& object_name,
    boost::shared_ptr<Database> input_db,
    boost::shared_ptr<AdvectorExplicitPredictorPatchOps> explicit_predictor,
    boost::shared_ptr<CartesianGridGeometry> grid_geom,
    bool register_for_restart)
    : AdvectorPredictorCorrectorHyperbolicPatchOps(object_name,
                                                   input_db,
                                                   explicit_predictor,
                                                   grid_geom,
                                                   register_for_restart)
{
    d_overwrite_tags = false;
    return;
}

AdvDiffPredictorCorrectorHyperbolicPatchOps::~AdvDiffPredictorCorrectorHyperbolicPatchOps()
{
    // intentionally blank
    return;
}

void AdvDiffPredictorCorrectorHyperbolicPatchOps::conservativeDifferenceOnPatch(Patch& patch,
                                                                                const double /*time*/,
                                                                                const double dt,
                                                                                bool /*at_synchronization*/)
{
    const Box& patch_box = patch.getBox();
    const Index& ilower = patch_box.lower();
    const Index& iupper = patch_box.upper();

    auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch.getPatchGeometry());
    const double* const dx = pgeom->getDx();

    for (auto cit = d_Q_var.begin(); cit != d_Q_var.end(); ++cit)
    {
        auto Q_var = *cit;
        auto Q_data = BOOST_CAST<CellData<double> >(patch.getPatchData(Q_var, getDataContext()));
        auto u_var = d_Q_u_map[Q_var];
        if (u_var)
        {
            const bool conservation_form = d_Q_difference_form[Q_var] == CONSERVATIVE;
            const bool u_is_div_free = d_u_is_div_free[u_var];

            auto flux_integral_var = d_flux_integral_var[Q_var];
            auto q_integral_var = d_q_integral_var[Q_var];
            auto u_integral_var = d_u_integral_var[u_var];

            auto flux_integral_data = BOOST_CAST<FaceData<double> >(
                conservation_form ? patch.getPatchData(flux_integral_var, getDataContext()) : NULL);
            auto q_integral_data = BOOST_CAST<FaceData<double> >(
                !conservation_form || !u_is_div_free ? patch.getPatchData(q_integral_var, getDataContext()) : NULL);
            auto u_integral_data = BOOST_CAST<FaceData<double> >(
                !conservation_form || !u_is_div_free ? patch.getPatchData(u_integral_var, getDataContext()) : NULL);

            const IntVector& Q_data_ghost_cells = Q_data->getGhostCellWidth();
            const IntVector& flux_integral_data_ghost_cells =
                (flux_integral_data ? flux_integral_data->getGhostCellWidth() : IntVector::getZero(DIM));
            const IntVector& q_integral_data_ghost_cells =
                (q_integral_data ? q_integral_data->getGhostCellWidth() : IntVector::getZero(DIM));
            const IntVector& u_integral_data_ghost_cells =
                (u_integral_data ? u_integral_data->getGhostCellWidth() : IntVector::getZero(DIM));

            switch (d_Q_difference_form[Q_var])
            {
            case CONSERVATIVE:
            {
                for (int depth = 0; depth < Q_data->getDepth(); ++depth)
                {
                    if (u_is_div_free)
                    {
#if (NDIM == 2)
                        ADV_DIFF_CONSDIFF_FC(dx, dt, ilower(0), iupper(0), ilower(1), iupper(1),
                                             flux_integral_data_ghost_cells(0), flux_integral_data_ghost_cells(1),
                                             Q_data_ghost_cells(0), Q_data_ghost_cells(1),
                                             flux_integral_data->getPointer(0, depth),
                                             flux_integral_data->getPointer(1, depth), Q_data->getPointer(depth));
#endif
#if (NDIM == 3)
                        ADV_DIFF_CONSDIFF_FC(dx, dt, ilower(0), iupper(0), ilower(1), iupper(1), ilower(2), iupper(2),
                                             flux_integral_data_ghost_cells(0), flux_integral_data_ghost_cells(1),
                                             flux_integral_data_ghost_cells(2), Q_data_ghost_cells(0),
                                             Q_data_ghost_cells(1), Q_data_ghost_cells(2),
                                             flux_integral_data->getPointer(0, depth),
                                             flux_integral_data->getPointer(1, depth),
                                             flux_integral_data->getPointer(2, depth), Q_data->getPointer(depth));
#endif
                    }
                    else
                    {
#if (NDIM == 2)
                        ADV_DIFF_CONSDIFFWITHDIVSOURCE_FC(
                            dx, dt, ilower(0), iupper(0), ilower(1), iupper(1), flux_integral_data_ghost_cells(0),
                            flux_integral_data_ghost_cells(1), q_integral_data_ghost_cells(0),
                            q_integral_data_ghost_cells(1), u_integral_data_ghost_cells(0),
                            u_integral_data_ghost_cells(1), Q_data_ghost_cells(0), Q_data_ghost_cells(1),
                            flux_integral_data->getPointer(0, depth), flux_integral_data->getPointer(1, depth),
                            q_integral_data->getPointer(0, depth), q_integral_data->getPointer(1, depth),
                            u_integral_data->getPointer(0), u_integral_data->getPointer(1), Q_data->getPointer(depth));
#endif
#if (NDIM == 3)
                        ADV_DIFF_CONSDIFFWITHDIVSOURCE_FC(
                            dx, dt, ilower(0), iupper(0), ilower(1), iupper(1), ilower(2), iupper(2),
                            flux_integral_data_ghost_cells(0), flux_integral_data_ghost_cells(1),
                            flux_integral_data_ghost_cells(2), q_integral_data_ghost_cells(0),
                            q_integral_data_ghost_cells(1), q_integral_data_ghost_cells(2),
                            u_integral_data_ghost_cells(0), u_integral_data_ghost_cells(1),
                            u_integral_data_ghost_cells(2), Q_data_ghost_cells(0), Q_data_ghost_cells(1),
                            Q_data_ghost_cells(2), flux_integral_data->getPointer(0, depth),
                            flux_integral_data->getPointer(1, depth), flux_integral_data->getPointer(2, depth),
                            q_integral_data->getPointer(0, depth), q_integral_data->getPointer(1, depth),
                            q_integral_data->getPointer(2, depth), u_integral_data->getPointer(0),
                            u_integral_data->getPointer(1), u_integral_data->getPointer(2), Q_data->getPointer(depth));
#endif
                    }
                }
                break;
            }
            case ADVECTIVE:
            {
                CellData<double> N_data(patch_box, Q_data->getDepth(), IntVector::getZero(DIM));
                d_explicit_predictor->computeAdvectiveDerivative(N_data, *u_integral_data, *q_integral_data, patch);
                PatchCellDataOpsReal<double> patch_cc_data_ops;
                patch_cc_data_ops.scale(Q_data, -1.0 / (dt * dt),
                                        boost::shared_ptr<CellData<double> >(&N_data, NullDeleter()), patch_box);
                break;
            }
            default:
            {
                TBOX_ERROR(
                    "AdvDiffPredictorCorrectorHyperbolicPatchOps::"
                    "conservativeDifferenceOnPatch():"
                    "\n"
                    << "  unsupported differencing form: "
                    << enum_to_string<ConvectiveDifferencingType>(d_Q_difference_form[Q_var]) << " \n"
                    << "  valid choices are: ADVECTIVE, CONSERVATIVE\n");
            }
            }
        }
        else
        {
            Q_data->fillAll(0.0);
        }
    }
    return;
}

void
AdvDiffPredictorCorrectorHyperbolicPatchOps::preprocessAdvanceLevelState(const boost::shared_ptr<PatchLevel>& level,
                                                                         double current_time,
                                                                         double /*dt*/,
                                                                         bool /*first_step*/,
                                                                         bool /*last_step*/,
                                                                         bool /*regrid_advance*/)
{
    if (!d_compute_init_velocity) return;

    // Update the advection velocity (or velocities).
    for (auto cit = d_u_var.begin(); cit != d_u_var.end(); ++cit)
    {
        auto u_var = *cit;
        if (d_u_fcn[u_var] && d_u_fcn[u_var]->isTimeDependent())
        {
            VariableDatabase* var_db = VariableDatabase::getDatabase();
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, d_integrator->getScratchContext());
            d_u_fcn[u_var]->setDataOnPatchLevel(u_idx, u_var, level, current_time);
        }
    }
    return;
}

void
AdvDiffPredictorCorrectorHyperbolicPatchOps::postprocessAdvanceLevelState(const boost::shared_ptr<PatchLevel>& level,
                                                                          double current_time,
                                                                          double dt,
                                                                          bool /*first_step*/,
                                                                          bool /*last_step*/,
                                                                          bool /*regrid_advance*/)
{
    if (!d_compute_final_velocity) return;

    // Update the advection velocity (or velocities).
    for (auto cit = d_u_var.begin(); cit != d_u_var.end(); ++cit)
    {
        auto u_var = *cit;
        if (d_u_fcn[u_var] && d_u_fcn[u_var]->isTimeDependent())
        {
            VariableDatabase* var_db = VariableDatabase::getDatabase();
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, d_integrator->getScratchContext());
            d_u_fcn[u_var]->setDataOnPatchLevel(u_idx, u_var, level, current_time + dt);
        }
    }
    return;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
