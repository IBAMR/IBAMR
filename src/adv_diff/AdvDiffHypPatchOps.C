// Filename: AdvDiffHypPatchOps.C
// Created on 19 Mar 2004 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "AdvDiffHypPatchOps.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <FaceData.h>
#include <Index.h>
#include <PatchCellDataOpsReal.h>
#include <PatchData.h>
#include <VariableDatabase.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

// C++ STDLIB INCLUDES
#include <vector>

// FORTRAN ROUTINES
#if (NDIM == 1)
#define ADV_DIFF_CONSDIFF_FC FC_FUNC_(adv_diff_consdiff1d, ADV_DIFF_CONSDIFF1D)
#define ADV_DIFF_CONSDIFFWITHDIVSOURCE_FC FC_FUNC_(adv_diff_consdiffwithdivsource1d, ADV_DIFF_CONSDIFFWITHDIVSOURCE1D)
#endif

#if (NDIM == 2)
#define ADV_DIFF_CONSDIFF_FC FC_FUNC_(adv_diff_consdiff2d, ADV_DIFF_CONSDIFF2D)
#define ADV_DIFF_CONSDIFFWITHDIVSOURCE_FC FC_FUNC_(adv_diff_consdiffwithdivsource2d, ADV_DIFF_CONSDIFFWITHDIVSOURCE2D)
#endif

#if (NDIM == 3)
#define ADV_DIFF_CONSDIFF_FC FC_FUNC_(adv_diff_consdiff3d, ADV_DIFF_CONSDIFF3D)
#define ADV_DIFF_CONSDIFFWITHDIVSOURCE_FC FC_FUNC_(adv_diff_consdiffwithdivsource3d, ADV_DIFF_CONSDIFFWITHDIVSOURCE3D)
#endif

extern "C"
{
    void
    ADV_DIFF_CONSDIFF_FC(
        const double*, const double&,
#if (NDIM == 1)
        const int& , const int& ,
        const int& ,
        const int& ,
        const double* ,
#endif
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
#endif
        double*);

    void
    ADV_DIFF_CONSDIFFWITHDIVSOURCE_FC(
        const double*, const double&,
#if (NDIM == 1)
        const int& , const int& ,
        const int& ,
        const int& ,
        const int& ,
        const int& ,
        const double* ,
        const double* ,
        const double* ,
#endif
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const double* , const double* ,
        const double* , const double* ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
#endif
        double*);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Pointer<Timer> t_conservative_difference_on_patch;
static Pointer<Timer> t_preprocess_advance_level_state;
static Pointer<Timer> t_postprocess_advance_level_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffHypPatchOps::AdvDiffHypPatchOps(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<GodunovAdvector> godunov_advector,
    Pointer<CartesianGridGeometry<NDIM> > grid_geom,
    bool register_for_restart)
    : AdvectHypPatchOps(object_name, input_db, godunov_advector, grid_geom, register_for_restart)
{
    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_conservative_difference_on_patch = TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHypPatchOps::conservativeDifferenceOnPatch()");
        t_preprocess_advance_level_state = TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHypPatchOps::preprocessAdvanceLevelState()");
        t_postprocess_advance_level_state = TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHypPatchOps::postprocessAdvanceLevelState()");
        timers_need_init = false;
    }
    return;
}// AdvDiffHypPatchOps

AdvDiffHypPatchOps::~AdvDiffHypPatchOps()
{
    // intentionally blank
    return;
}// ~AdvDiffHypPatchOps

const std::string&
AdvDiffHypPatchOps::getName() const
{
    return AdvectHypPatchOps::getName();
}// getName

///
///  The following routines:
///
///      conservativeDifferenceOnPatch(),
///      preprocessAdvanceLevelState(),
///      postprocessAdvanceLevelState()
///
///  are redefined from the AdvectHypPatchOps base class.
///

void
AdvDiffHypPatchOps::conservativeDifferenceOnPatch(
    Patch<NDIM>& patch,
    const double time,
    const double dt,
    bool at_synchronization)
{
    t_conservative_difference_on_patch->start();

    (void) time;
    (void) at_synchronization;

    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const Index<NDIM>& ilower = patch.getBox().lower();
    const Index<NDIM>& iupper = patch.getBox().upper();

    const Box<NDIM>& patch_box = patch.getBox();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.size() == d_flux_integral_var.size());
    TBOX_ASSERT(d_Q_var.size() == d_q_integral_var.size());
#endif

    typedef std::vector<Pointer<CellVariable<NDIM,double> > > CellVariableVector;

    Pointer<FaceData<NDIM,double> > u_integral_data =
        (!d_u_integral_var.isNull()
         ? patch.getPatchData(d_u_integral_var, getDataContext())
         : Pointer<PatchData<NDIM> >(NULL));

    const IntVector<NDIM>& u_integral_data_ghost_cells =
        (!d_u_integral_var.isNull()
         ? u_integral_data->getGhostCellWidth()
         : 0);

    for (CellVariableVector::size_type l = 0; l < d_Q_var.size(); ++l)
    {
        Pointer<CellData<NDIM,double> > Q_data =
            patch.getPatchData(d_Q_var[l], getDataContext());
        Pointer<FaceData<NDIM,double> > flux_integral_data =
            (d_Q_in_consv_form[l]
             ? patch.getPatchData(d_flux_integral_var[l], getDataContext())
             : Pointer<PatchData<NDIM> >(NULL));
        Pointer<FaceData<NDIM,double> > q_integral_data =
            ((!d_u_is_div_free) || (!d_Q_in_consv_form[l])
             ? patch.getPatchData(d_q_integral_var[l], getDataContext())
             : Pointer<PatchData<NDIM> >(NULL));

        const IntVector<NDIM>& Q_data_ghost_cells = Q_data->getGhostCellWidth();
        const IntVector<NDIM>& flux_integral_data_ghost_cells =
            (d_Q_in_consv_form[l]
             ? flux_integral_data->getGhostCellWidth()
             : 0);
        const IntVector<NDIM>& q_integral_data_ghost_cells =
            ((!d_u_is_div_free) || (!d_Q_in_consv_form[l])
             ? q_integral_data->getGhostCellWidth()
             : 0);

        if (d_Q_in_consv_form[l])
        {
            for (int depth = 0; depth < Q_data->getDepth(); ++depth)
            {
                if (d_u_is_div_free)
                {
#if (NDIM == 1)
                    ADV_DIFF_CONSDIFF_FC(
                        dx,dt,
                        ilower(0),iupper(0),
                        flux_integral_data_ghost_cells(0),
                        Q_data_ghost_cells(0),
                        flux_integral_data->getPointer(0,depth),
                        Q_data->getPointer(depth));
#endif
#if (NDIM == 2)
                    ADV_DIFF_CONSDIFF_FC(
                        dx,dt,
                        ilower(0),iupper(0),ilower(1),iupper(1),
                        flux_integral_data_ghost_cells(0),flux_integral_data_ghost_cells(1),
                        Q_data_ghost_cells(0),Q_data_ghost_cells(1),
                        flux_integral_data->getPointer(0,depth),
                        flux_integral_data->getPointer(1,depth),
                        Q_data->getPointer(depth));
#endif
#if (NDIM == 3)
                    ADV_DIFF_CONSDIFF_FC(
                        dx,dt,
                        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                        flux_integral_data_ghost_cells(0),flux_integral_data_ghost_cells(1),flux_integral_data_ghost_cells(2),
                        Q_data_ghost_cells(0),Q_data_ghost_cells(1),Q_data_ghost_cells(2),
                        flux_integral_data->getPointer(0,depth),
                        flux_integral_data->getPointer(1,depth),
                        flux_integral_data->getPointer(2,depth),
                        Q_data->getPointer(depth));
#endif
                }
                else
                {
#if (NDIM == 1)
                    ADV_DIFF_CONSDIFFWITHDIVSOURCE_FC(
                        dx,dt,
                        ilower(0),iupper(0),
                        flux_integral_data_ghost_cells(0),
                        q_integral_data_ghost_cells(0),
                        u_integral_data_ghost_cells(0),
                        Q_data_ghost_cells(0),
                        flux_integral_data->getPointer(0,depth),
                        q_integral_data->getPointer(0,depth),
                        u_integral_data->getPointer(0),
                        Q_data->getPointer(depth));
#endif
#if (NDIM == 2)
                    ADV_DIFF_CONSDIFFWITHDIVSOURCE_FC(
                        dx,dt,
                        ilower(0),iupper(0),ilower(1),iupper(1),
                        flux_integral_data_ghost_cells(0),flux_integral_data_ghost_cells(1),
                        q_integral_data_ghost_cells(0),q_integral_data_ghost_cells(1),
                        u_integral_data_ghost_cells(0),u_integral_data_ghost_cells(1),
                        Q_data_ghost_cells(0),Q_data_ghost_cells(1),
                        flux_integral_data->getPointer(0,depth),
                        flux_integral_data->getPointer(1,depth),
                        q_integral_data->getPointer(0,depth),
                        q_integral_data->getPointer(1,depth),
                        u_integral_data->getPointer(0),
                        u_integral_data->getPointer(1),
                        Q_data->getPointer(depth));
#endif
#if (NDIM == 3)
                    ADV_DIFF_CONSDIFFWITHDIVSOURCE_FC(
                        dx,dt,
                        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                        flux_integral_data_ghost_cells(0),flux_integral_data_ghost_cells(1),flux_integral_data_ghost_cells(2),
                        q_integral_data_ghost_cells(0),q_integral_data_ghost_cells(1),q_integral_data_ghost_cells(2),
                        u_integral_data_ghost_cells(0),u_integral_data_ghost_cells(1),u_integral_data_ghost_cells(2),
                        Q_data_ghost_cells(0),Q_data_ghost_cells(1),Q_data_ghost_cells(2),
                        flux_integral_data->getPointer(0,depth),
                        flux_integral_data->getPointer(1,depth),
                        flux_integral_data->getPointer(2,depth),
                        q_integral_data->getPointer(0,depth),
                        q_integral_data->getPointer(1,depth),
                        q_integral_data->getPointer(2,depth),
                        u_integral_data->getPointer(0),
                        u_integral_data->getPointer(1),
                        u_integral_data->getPointer(2),
                        Q_data->getPointer(depth));
#endif
                }
            }
        }
        else
        {
            Pointer<PatchCellDataOpsReal<NDIM,double> > patch_cc_data_ops =
                new PatchCellDataOpsReal<NDIM,double>();
            Pointer<CellData<NDIM,double> > N_data =
                new CellData<NDIM,double>(patch_box,Q_data->getDepth(),0);

            d_godunov_advector->computeAdvectiveDerivative(
                *N_data, *u_integral_data, *q_integral_data, patch);

            patch_cc_data_ops->scale(Q_data,       // dst
                                     -1.0/(dt*dt), // alpha
                                     N_data,       // src1
                                     patch_box);
        }
    }

    t_conservative_difference_on_patch->stop();
    return;
}// conservativeDifferenceOnPatch

#if 0
// NOTE: Should the following method be removed?
void
AdvDiffHypPatchOps::preprocessAdvanceLevelState(
    const Pointer<PatchLevel<NDIM> >& level,
    double current_time,
    double dt,
    bool first_step,
    bool last_step,
    bool regrid_advance)
{
    t_preprocess_advance_level_state->start();

    (void) dt;
    (void) first_step;
    (void) last_step;
    (void) regrid_advance;

    // Update the advection velocity.
    if (!d_u_fcn.isNull() && d_u_fcn->isTimeDependent() && d_compute_init_velocity)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int u_idx = var_db->mapVariableAndContextToIndex(
            d_u_var, d_integrator->getScratchContext());
        d_u_fcn->setDataOnPatchLevel(u_idx, d_u_var, level, current_time);
    }

    t_preprocess_advance_level_state->stop();
    return;
}// preprocessAdvanceLevelState
#endif

void
AdvDiffHypPatchOps::postprocessAdvanceLevelState(
    const Pointer<PatchLevel<NDIM> >& level,
    double current_time,
    double dt,
    bool first_step,
    bool last_step,
    bool regrid_advance)
{
    t_postprocess_advance_level_state->start();

    (void) first_step;
    (void) last_step;
    (void) regrid_advance;

    // Update the advection velocity.
    if (!d_u_fcn.isNull() && d_u_fcn->isTimeDependent() && d_compute_final_velocity)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int u_idx = var_db->mapVariableAndContextToIndex(
            d_u_var, d_integrator->getNewContext());
        d_u_fcn->setDataOnPatchLevel(u_idx, d_u_var, level, current_time+dt);
    }

    t_postprocess_advance_level_state->stop();
    return;
}// postprocessAdvanceLevelState

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::AdvDiffHypPatchOps>;

//////////////////////////////////////////////////////////////////////////////
