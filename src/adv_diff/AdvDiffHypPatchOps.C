// Filename: AdvDiffHypPatchOps.C
// Last modified: <25.Oct.2006 17:18:55 boyce@bigboy.nyconnect.com>
// Created on 19 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

#include "AdvDiffHypPatchOps.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#define included_IBAMR_config
#include <IBAMR_config.h>
#endif

#ifndef included_SAMRAI_config
#define included_SAMRAI_config
#include <SAMRAI_config.h>
#endif

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
#include <cassert>
#include <vector>

// FORTRAN ROUTINES
#if (NDIM == 1)
#define ADV_DIFF_CONSDIFF_F77 F77_FUNC_(adv_diff_consdiff1d, ADV_DIFF_CONSDIFF1D)
#define ADV_DIFF_CONSDIFFWITHDIVSOURCE_F77 F77_FUNC_(adv_diff_consdiffwithdivsource1d, ADV_DIFF_CONSDIFFWITHDIVSOURCE1D)
#endif

#if (NDIM == 2)
#define ADV_DIFF_CONSDIFF_F77 F77_FUNC_(adv_diff_consdiff2d, ADV_DIFF_CONSDIFF2D)
#define ADV_DIFF_CONSDIFFWITHDIVSOURCE_F77 F77_FUNC_(adv_diff_consdiffwithdivsource2d, ADV_DIFF_CONSDIFFWITHDIVSOURCE2D)
#endif

#if (NDIM == 3)
#define ADV_DIFF_CONSDIFF_F77 F77_FUNC_(adv_diff_consdiff3d, ADV_DIFF_CONSDIFF3D)
#define ADV_DIFF_CONSDIFFWITHDIVSOURCE_F77 F77_FUNC_(adv_diff_consdiffwithdivsource3d, ADV_DIFF_CONSDIFFWITHDIVSOURCE3D)
#endif

extern "C"
{
    void ADV_DIFF_CONSDIFF_F77(
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

    void ADV_DIFF_CONSDIFFWITHDIVSOURCE_F77(
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
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_conservative_difference_on_patch;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_preprocess_advance_level_state;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_postprocess_advance_level_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffHypPatchOps::AdvDiffHypPatchOps(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    SAMRAI::tbox::Pointer<GodunovAdvector> godunov_advector,
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom,
    bool register_for_restart)
    : GodunovHypPatchOps(object_name, input_db, godunov_advector, grid_geom, register_for_restart)
{
    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_conservative_difference_on_patch = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHypPatchOps::conservativeDifferenceOnPatch()");
        t_preprocess_advance_level_state = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHypPatchOps::preprocessAdvanceLevelState()");
        t_postprocess_advance_level_state = SAMRAI::tbox::TimerManager::getManager()->
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

///
///  The following routines:
///
///      conservativeDifferenceOnPatch(),
///      preprocessAdvanceLevelState(),
///      postprocessAdvanceLevelState()
///
///  are redefined from the GodunovHypPatchOps base class.
///

void
AdvDiffHypPatchOps::conservativeDifferenceOnPatch(
    SAMRAI::hier::Patch<NDIM>& patch,
    const double time,
    const double dt,
    bool at_synchronization)
{
    t_conservative_difference_on_patch->start();

    (void) time;
    (void) at_synchronization;

    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const SAMRAI::hier::Index<NDIM>& ilower = patch.getBox().lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch.getBox().upper();

    const SAMRAI::hier::Box<NDIM>& patch_box = patch.getBox();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_Q_vars.size() == d_flux_integral_vars.size());
    assert(d_Q_vars.size() == d_q_integral_vars.size());
#endif

    typedef std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > CellVariableVector;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_integral_data =
        (!d_u_integral_var.isNull()
         ? patch.getPatchData(d_u_integral_var, getDataContext())
         : SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >(NULL));

    const SAMRAI::hier::IntVector<NDIM>& u_integral_data_ghost_cells =
        (!d_u_integral_var.isNull()
         ? u_integral_data->getGhostCellWidth()
         : 0);

    for (CellVariableVector::size_type l = 0; l < d_Q_vars.size(); ++l)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Q_data =
            patch.getPatchData(d_Q_vars[l], getDataContext());
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > flux_integral_data =
            (d_Q_conservation_form[l]
             ? patch.getPatchData(d_flux_integral_vars[l], getDataContext())
             : SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >(NULL));
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > q_integral_data =
            ((!d_u_is_div_free) || (!d_Q_conservation_form[l])
             ? patch.getPatchData(d_q_integral_vars[l], getDataContext())
             : SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >(NULL));

        const SAMRAI::hier::IntVector<NDIM>& Q_data_ghost_cells = Q_data->getGhostCellWidth();
        const SAMRAI::hier::IntVector<NDIM>& flux_integral_data_ghost_cells =
            (d_Q_conservation_form[l]
             ? flux_integral_data->getGhostCellWidth()
             : 0);
        const SAMRAI::hier::IntVector<NDIM>& q_integral_data_ghost_cells =
            ((!d_u_is_div_free) || (!d_Q_conservation_form[l])
             ? q_integral_data->getGhostCellWidth()
             : 0);

        if (d_Q_conservation_form[l])
        {
            for (int depth = 0; depth < Q_data->getDepth(); ++depth)
            {
                if (d_u_is_div_free)
                {
#if (NDIM == 1)
                    ADV_DIFF_CONSDIFF_F77(
                        dx,dt,
                        ilower(0),iupper(0),
                        flux_integral_data_ghost_cells(0),
                        Q_data_ghost_cells(0),
                        flux_integral_data->getPointer(0,depth),
                        Q_data->getPointer(depth));
#endif
#if (NDIM == 2)
                    ADV_DIFF_CONSDIFF_F77(
                        dx,dt,
                        ilower(0),iupper(0),ilower(1),iupper(1),
                        flux_integral_data_ghost_cells(0),flux_integral_data_ghost_cells(1),
                        Q_data_ghost_cells(0),Q_data_ghost_cells(1),
                        flux_integral_data->getPointer(0,depth),
                        flux_integral_data->getPointer(1,depth),
                        Q_data->getPointer(depth));
#endif
#if (NDIM == 3)
                    ADV_DIFF_CONSDIFF_F77(
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
                    ADV_DIFF_CONSDIFFWITHDIVSOURCE_F77(
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
                    ADV_DIFF_CONSDIFFWITHDIVSOURCE_F77(
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
                    ADV_DIFF_CONSDIFFWITHDIVSOURCE_F77(
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
            SAMRAI::tbox::Pointer<SAMRAI::math::PatchCellDataOpsReal<NDIM,double> > patch_cc_data_ops =
                new SAMRAI::math::PatchCellDataOpsReal<NDIM,double>();
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > N_data =
                new SAMRAI::pdat::CellData<NDIM,double>(patch_box,Q_data->getDepth(),0);

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

void
AdvDiffHypPatchOps::preprocessAdvanceLevelState(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& level,
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
    if (!d_u_set.isNull() && d_u_set->isTimeDependent() &&
        d_compute_init_velocity)
    {
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        const int u_idx = var_db->mapVariableAndContextToIndex(
            d_u_var, d_integrator->getScratchContext());
        d_u_set->setDataOnPatchLevel(u_idx, d_u_var, level, current_time);
    }

    t_preprocess_advance_level_state->stop();
    return;
}// preprocessAdvanceLevelState

void
AdvDiffHypPatchOps::postprocessAdvanceLevelState(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& level,
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
    if (!d_u_set.isNull() && d_u_set->isTimeDependent() &&
        d_compute_final_velocity)
    {
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        const int u_idx = var_db->mapVariableAndContextToIndex(
            d_u_var, d_integrator->getNewContext());
        d_u_set->setDataOnPatchLevel(u_idx, d_u_var, level, current_time+dt);
    }

    t_postprocess_advance_level_state->stop();
    return;
}// postprocessAdvanceLevelState

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::AdvDiffHypPatchOps>;

//////////////////////////////////////////////////////////////////////////////
