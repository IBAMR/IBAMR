// Filename: AdvDiffHypPatchOps.C
// Created on 19 Mar 2004 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
// POSSIBILITY OF SUCH DAMAGE.

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

// C++ STDLIB INCLUDES
#include <vector>

// FORTRAN ROUTINES
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

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffHypPatchOps::AdvDiffHypPatchOps(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<GodunovAdvector> godunov_advector,
    Pointer<CartesianGridGeometry<NDIM> > grid_geom,
    bool register_for_restart)
    : AdvectHypPatchOps(object_name, input_db, godunov_advector, grid_geom, register_for_restart)
{
    // intentionally blank
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
    const double /*time*/,
    const double dt,
    bool /*at_synchronization*/)
{
    const Box<NDIM>& patch_box = patch.getBox();
    const Index<NDIM>& ilower = patch_box.lower();
    const Index<NDIM>& iupper = patch_box.upper();

    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        Pointer<FaceVariable<NDIM,double> > u_var = d_Q_u_map[Q_var];

        const bool conservation_form = d_Q_difference_form[Q_var] == CONSERVATIVE;
        const bool u_is_div_free = d_u_is_div_free[u_var];

        Pointer<FaceVariable<NDIM,double> > flux_integral_var = d_flux_integral_var[Q_var];
        Pointer<FaceVariable<NDIM,double> > q_integral_var = d_q_integral_var[Q_var];
        Pointer<FaceVariable<NDIM,double> > u_integral_var = d_u_integral_var[u_var];

        Pointer<CellData<NDIM,double> > Q_data = patch.getPatchData(Q_var, getDataContext());
        Pointer<FaceData<NDIM,double> > flux_integral_data =
            (conservation_form
             ? patch.getPatchData(flux_integral_var, getDataContext())
             : Pointer<PatchData<NDIM> >(NULL));
        Pointer<FaceData<NDIM,double> > q_integral_data =
            (!conservation_form || !u_is_div_free
             ? patch.getPatchData(q_integral_var, getDataContext())
             : Pointer<PatchData<NDIM> >(NULL));
        Pointer<FaceData<NDIM,double> > u_integral_data =
            (!conservation_form || !u_is_div_free
             ? patch.getPatchData(u_integral_var, getDataContext())
             : Pointer<PatchData<NDIM> >(NULL));

        const IntVector<NDIM>& Q_data_ghost_cells = Q_data->getGhostCellWidth();
        const IntVector<NDIM>& flux_integral_data_ghost_cells =
            (!flux_integral_data.isNull()
             ? flux_integral_data->getGhostCellWidth()
             : 0);
        const IntVector<NDIM>& q_integral_data_ghost_cells =
            (!q_integral_data.isNull()
             ? q_integral_data->getGhostCellWidth()
             : 0);
        const IntVector<NDIM>& u_integral_data_ghost_cells =
            (!u_integral_data.isNull()
             ? u_integral_data->getGhostCellWidth()
             : 0);

        switch (d_Q_difference_form[Q_var])
        {
            case CONSERVATIVE:
            {
                for (int depth = 0; depth < Q_data->getDepth(); ++depth)
                {
                    if (u_is_div_free)
                    {
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
                break;
            }
            case ADVECTIVE:
            {
                CellData<NDIM,double> N_data(patch_box,Q_data->getDepth(),0);
                d_godunov_advector->computeAdvectiveDerivative(N_data, *u_integral_data, *q_integral_data, patch);
                PatchCellDataOpsReal<NDIM,double> patch_cc_data_ops;
                patch_cc_data_ops.scale(Q_data, -1.0/(dt*dt), Pointer<CellData<NDIM,double> >(&N_data,false), patch_box);
                break;
            }
            default:
            {
                TBOX_ERROR("AdvDiffHypPatchOps::conservativeDifferenceOnPatch():\n"
                           << "  unsupported differencing form: " << enum_to_string<ConvectiveDifferencingType>(d_Q_difference_form[Q_var]) << " \n"
                           << "  valid choices are: ADVECTIVE, CONSERVATIVE\n");
            }
        }
    }
    return;
}// conservativeDifferenceOnPatch

void
AdvDiffHypPatchOps::preprocessAdvanceLevelState(
    const Pointer<PatchLevel<NDIM> >& level,
    double current_time,
    double /*dt*/,
    bool /*first_step*/,
    bool /*last_step*/,
    bool /*regrid_advance*/)
{
    if (!d_compute_init_velocity) return;

    // Update the advection velocity (or velocities).
    for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin();
         cit != d_u_var.end(); ++cit)
    {
        Pointer<FaceVariable<NDIM,double> > u_var = *cit;
        if (!d_u_fcn[u_var].isNull() && d_u_fcn[u_var]->isTimeDependent())
        {
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, d_integrator->getScratchContext());
            d_u_fcn[u_var]->setDataOnPatchLevel(u_idx, u_var, level, current_time);
        }
    }
    return;
}// preprocessAdvanceLevelState

void
AdvDiffHypPatchOps::postprocessAdvanceLevelState(
    const Pointer<PatchLevel<NDIM> >& level,
    double current_time,
    double dt,
    bool /*first_step*/,
    bool /*last_step*/,
    bool /*regrid_advance*/)
{
    if (!d_compute_final_velocity) return;

    // Update the advection velocity (or velocities).
    for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin();
         cit != d_u_var.end(); ++cit)
    {
        Pointer<FaceVariable<NDIM,double> > u_var = *cit;
        if (!d_u_fcn[u_var].isNull() && d_u_fcn[u_var]->isTimeDependent())
        {
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, d_integrator->getScratchContext());
            d_u_fcn[u_var]->setDataOnPatchLevel(u_idx, u_var, level, current_time+dt);
        }
    }
    return;
}// postprocessAdvanceLevelState

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
