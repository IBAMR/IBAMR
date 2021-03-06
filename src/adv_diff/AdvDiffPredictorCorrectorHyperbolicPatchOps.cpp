// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/AdvDiffPredictorCorrectorHyperbolicPatchOps.h"
#include "ibamr/AdvectorExplicitPredictorPatchOps.h"
#include "ibamr/AdvectorPredictorCorrectorHyperbolicPatchOps.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CartGridFunction.h"

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "HyperbolicLevelIntegrator.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchCellDataOpsReal.h"
#include "PatchData.h"
#include "PatchLevel.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <map>
#include <ostream>
#include <set>
#include <string>
#include <utility>

#include "ibamr/namespaces.h" // IWYU pragma: keep

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

extern "C"
{
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
    std::string object_name,
    Pointer<Database> input_db,
    Pointer<AdvectorExplicitPredictorPatchOps> explicit_predictor,
    Pointer<CartesianGridGeometry<NDIM> > grid_geom,
    bool register_for_restart)
    : AdvectorPredictorCorrectorHyperbolicPatchOps(std::move(object_name),
                                                   input_db,
                                                   explicit_predictor,
                                                   grid_geom,
                                                   register_for_restart)
{
    d_overwrite_tags = false;
    return;
} // AdvDiffPredictorCorrectorHyperbolicPatchOps

void
AdvDiffPredictorCorrectorHyperbolicPatchOps::conservativeDifferenceOnPatch(Patch<NDIM>& patch,
                                                                           const double /*time*/,
                                                                           const double dt,
                                                                           bool /*at_synchronization*/)
{
    const Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& ilower = patch_box.lower();
    const hier::Index<NDIM>& iupper = patch_box.upper();

    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    for (const auto& Q_var : d_Q_var)
    {
        Pointer<CellData<NDIM, double> > Q_data = patch.getPatchData(Q_var, getDataContext());
        Pointer<FaceVariable<NDIM, double> > u_var = d_Q_u_map[Q_var];
        if (u_var)
        {
            const bool conservation_form = d_Q_difference_form[Q_var] == CONSERVATIVE;
            const bool u_is_div_free = d_u_is_div_free[u_var];

            Pointer<FaceVariable<NDIM, double> > flux_integral_var = d_flux_integral_var[Q_var];
            Pointer<FaceVariable<NDIM, double> > q_integral_var = d_q_integral_var[Q_var];
            Pointer<FaceVariable<NDIM, double> > u_integral_var = d_u_integral_var[u_var];

            Pointer<FaceData<NDIM, double> > flux_integral_data =
                (conservation_form ? patch.getPatchData(flux_integral_var, getDataContext()) :
                                     Pointer<PatchData<NDIM> >(nullptr));
            Pointer<FaceData<NDIM, double> > q_integral_data =
                (!conservation_form || !u_is_div_free ? patch.getPatchData(q_integral_var, getDataContext()) :
                                                        Pointer<PatchData<NDIM> >(nullptr));
            Pointer<FaceData<NDIM, double> > u_integral_data =
                (!conservation_form || !u_is_div_free ? patch.getPatchData(u_integral_var, getDataContext()) :
                                                        Pointer<PatchData<NDIM> >(nullptr));

            const IntVector<NDIM>& Q_data_ghost_cells = Q_data->getGhostCellWidth();
            const IntVector<NDIM>& flux_integral_data_ghost_cells =
                (flux_integral_data ? flux_integral_data->getGhostCellWidth() : 0);
            const IntVector<NDIM>& q_integral_data_ghost_cells =
                (q_integral_data ? q_integral_data->getGhostCellWidth() : 0);
            const IntVector<NDIM>& u_integral_data_ghost_cells =
                (u_integral_data ? u_integral_data->getGhostCellWidth() : 0);

            switch (d_Q_difference_form[Q_var])
            {
            case CONSERVATIVE:
            {
                for (int depth = 0; depth < Q_data->getDepth(); ++depth)
                {
                    if (u_is_div_free)
                    {
#if (NDIM == 2)
                        ADV_DIFF_CONSDIFF_FC(dx,
                                             dt,
                                             ilower(0),
                                             iupper(0),
                                             ilower(1),
                                             iupper(1),
                                             flux_integral_data_ghost_cells(0),
                                             flux_integral_data_ghost_cells(1),
                                             Q_data_ghost_cells(0),
                                             Q_data_ghost_cells(1),
                                             flux_integral_data->getPointer(0, depth),
                                             flux_integral_data->getPointer(1, depth),
                                             Q_data->getPointer(depth));
#endif
#if (NDIM == 3)
                        ADV_DIFF_CONSDIFF_FC(dx,
                                             dt,
                                             ilower(0),
                                             iupper(0),
                                             ilower(1),
                                             iupper(1),
                                             ilower(2),
                                             iupper(2),
                                             flux_integral_data_ghost_cells(0),
                                             flux_integral_data_ghost_cells(1),
                                             flux_integral_data_ghost_cells(2),
                                             Q_data_ghost_cells(0),
                                             Q_data_ghost_cells(1),
                                             Q_data_ghost_cells(2),
                                             flux_integral_data->getPointer(0, depth),
                                             flux_integral_data->getPointer(1, depth),
                                             flux_integral_data->getPointer(2, depth),
                                             Q_data->getPointer(depth));
#endif
                    }
                    else
                    {
#if (NDIM == 2)
                        ADV_DIFF_CONSDIFFWITHDIVSOURCE_FC(dx,
                                                          dt,
                                                          ilower(0),
                                                          iupper(0),
                                                          ilower(1),
                                                          iupper(1),
                                                          flux_integral_data_ghost_cells(0),
                                                          flux_integral_data_ghost_cells(1),
                                                          q_integral_data_ghost_cells(0),
                                                          q_integral_data_ghost_cells(1),
                                                          u_integral_data_ghost_cells(0),
                                                          u_integral_data_ghost_cells(1),
                                                          Q_data_ghost_cells(0),
                                                          Q_data_ghost_cells(1),
                                                          flux_integral_data->getPointer(0, depth),
                                                          flux_integral_data->getPointer(1, depth),
                                                          q_integral_data->getPointer(0, depth),
                                                          q_integral_data->getPointer(1, depth),
                                                          u_integral_data->getPointer(0),
                                                          u_integral_data->getPointer(1),
                                                          Q_data->getPointer(depth));
#endif
#if (NDIM == 3)
                        ADV_DIFF_CONSDIFFWITHDIVSOURCE_FC(dx,
                                                          dt,
                                                          ilower(0),
                                                          iupper(0),
                                                          ilower(1),
                                                          iupper(1),
                                                          ilower(2),
                                                          iupper(2),
                                                          flux_integral_data_ghost_cells(0),
                                                          flux_integral_data_ghost_cells(1),
                                                          flux_integral_data_ghost_cells(2),
                                                          q_integral_data_ghost_cells(0),
                                                          q_integral_data_ghost_cells(1),
                                                          q_integral_data_ghost_cells(2),
                                                          u_integral_data_ghost_cells(0),
                                                          u_integral_data_ghost_cells(1),
                                                          u_integral_data_ghost_cells(2),
                                                          Q_data_ghost_cells(0),
                                                          Q_data_ghost_cells(1),
                                                          Q_data_ghost_cells(2),
                                                          flux_integral_data->getPointer(0, depth),
                                                          flux_integral_data->getPointer(1, depth),
                                                          flux_integral_data->getPointer(2, depth),
                                                          q_integral_data->getPointer(0, depth),
                                                          q_integral_data->getPointer(1, depth),
                                                          q_integral_data->getPointer(2, depth),
                                                          u_integral_data->getPointer(0),
                                                          u_integral_data->getPointer(1),
                                                          u_integral_data->getPointer(2),
                                                          Q_data->getPointer(depth));
#endif
                    }
                }
                break;
            }
            case ADVECTIVE:
            {
                CellData<NDIM, double> N_data(patch_box, Q_data->getDepth(), 0);
                d_explicit_predictor->computeAdvectiveDerivative(N_data, *u_integral_data, *q_integral_data, patch);
                PatchCellDataOpsReal<NDIM, double> patch_cc_data_ops;
                patch_cc_data_ops.scale(
                    Q_data, -1.0 / (dt * dt), Pointer<CellData<NDIM, double> >(&N_data, false), patch_box);
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
} // conservativeDifferenceOnPatch

void
AdvDiffPredictorCorrectorHyperbolicPatchOps::preprocessAdvanceLevelState(const Pointer<PatchLevel<NDIM> >& level,
                                                                         double current_time,
                                                                         double /*dt*/,
                                                                         bool /*first_step*/,
                                                                         bool /*last_step*/,
                                                                         bool /*regrid_advance*/)
{
    if (!d_compute_init_velocity) return;

    // Update the advection velocity (or velocities).
    for (const auto& u_var : d_u_var)
    {
        if (d_u_fcn[u_var] && d_u_fcn[u_var]->isTimeDependent())
        {
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, d_integrator->getScratchContext());
            d_u_fcn[u_var]->setDataOnPatchLevel(u_idx, u_var, level, current_time);
        }
    }
    return;
} // preprocessAdvanceLevelState

void
AdvDiffPredictorCorrectorHyperbolicPatchOps::postprocessAdvanceLevelState(const Pointer<PatchLevel<NDIM> >& level,
                                                                          double current_time,
                                                                          double dt,
                                                                          bool /*first_step*/,
                                                                          bool /*last_step*/,
                                                                          bool /*regrid_advance*/)
{
    if (!d_compute_final_velocity) return;

    // Update the advection velocity (or velocities).
    for (const auto& u_var : d_u_var)
    {
        if (d_u_fcn[u_var] && d_u_fcn[u_var]->isTimeDependent())
        {
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, d_integrator->getScratchContext());
            d_u_fcn[u_var]->setDataOnPatchLevel(u_idx, u_var, level, current_time + dt);
        }
    }
    return;
} // postprocessAdvanceLevelState

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
