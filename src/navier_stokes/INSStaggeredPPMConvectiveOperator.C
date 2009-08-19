// Filename: INSStaggeredPPMConvectiveOperator.C
// Last modified: <19.Aug.2009 17:36:05 griffith@boyce-griffiths-mac-pro.local>
// Created on 08 May 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

#include "INSStaggeredPPMConvectiveOperator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/CartExtrapPhysBdryOp.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <FaceData.h>
#include <FaceGeometry.h>
#include <SideData.h>
#include <SideGeometry.h>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define ADVECT_DERIVATIVE_FC FC_FUNC_(advect_derivative2d, ADVECT_DERIVATIVE2D)
#define CONVECT_DERIVATIVE_FC FC_FUNC_(convect_derivative2d, CONVECT_DERIVATIVE2D)
#define GODUNOV_PREDICT_FC FC_FUNC_(godunov_predict2d, GODUNOV_PREDICT2D)
#define NAVIER_STOKES_INTERP_COMPS_FC FC_FUNC_(navier_stokes_interp_comps2d, NAVIER_STOKES_INTERP_COMPS2D)
#define NAVIER_STOKES_RESET_ADV_VELOCITY_FC FC_FUNC_(navier_stokes_reset_adv_velocity2d, NAVIER_STOKES_RESET_ADV_VELOCITY2D)
#endif

#if (NDIM == 3)
#define ADVECT_DERIVATIVE_FC FC_FUNC_(advect_derivative3d, ADVECT_DERIVATIVE3D)
#define CONVECT_DERIVATIVE_FC FC_FUNC_(convect_derivative3d, CONVECT_DERIVATIVE3D)
#define GODUNOV_PREDICT_FC FC_FUNC_(godunov_predict3d, GODUNOV_PREDICT3D)
#define NAVIER_STOKES_INTERP_COMPS_FC FC_FUNC_(navier_stokes_interp_comps3d, NAVIER_STOKES_INTERP_COMPS3D)
#define NAVIER_STOKES_RESET_ADV_VELOCITY_FC FC_FUNC_(navier_stokes_reset_adv_velocity3d, NAVIER_STOKES_RESET_ADV_VELOCITY3D)
#endif

extern "C"
{
    void
    ADVECT_DERIVATIVE_FC(
        const double*,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const double* , const double* ,
        const int& , const int& ,
        double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
        double*
#endif
                          );

    void
    CONVECT_DERIVATIVE_FC(
        const double*,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const double* , const double* ,
        const int& , const int& ,
        double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
        double*
#endif
                           );

    void
    GODUNOV_PREDICT_FC(
        const double* , const double& ,
#if (NDIM == 3)
        const unsigned int& ,
#endif
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , double* , double* , double* , double* ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        double* , double* ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , double* , double* , double* , double* , double* ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        double* , double* , double* ,
        double* , double* , double*
#endif
                                    );

    void
    NAVIER_STOKES_INTERP_COMPS_FC(
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        double* , double* ,
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        double* , double* , double* ,
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        double* , double* , double* ,
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        double* , double* , double*
#endif
                                   );

    void
    NAVIER_STOKES_RESET_ADV_VELOCITY_FC(
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        double* , double* ,
        const int& , const int& ,
        const double* , const double* ,
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        double* , double* ,
        const int& , const int& ,
        const double* , const double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        double* , double* , double* ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        double* , double* , double* ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        double* , double* , double* ,
        const int& , const int& , const int& ,
        const double* , const double* , const double*
#endif
                                         );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// NOTE: The number of ghost cells required by the Godunov advection scheme
// depends on the order of the reconstruction.  These values were chosen to work
// with xsPPM7 (the modified piecewise parabolic method of Rider, Greenough, and
// Kamm).
static const int GADVECTG = 4;

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredPPMConvectiveOperator::INSStaggeredPPMConvectiveOperator(
    const INSCoefs& problem_coefs,
    const bool conservation_form)
  : d_is_initialized(false),
    d_problem_coefs(problem_coefs),
    d_conservation_form(conservation_form),
    d_refine_alg(NULL),
    d_refine_op(NULL),
    d_refine_scheds(),
    d_hierarchy(NULL),
    d_coarsest_ln(-1),
    d_finest_ln(-1),
    d_U_var(NULL),
    d_U_scratch_idx(-1)
{
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> context = var_db->getContext("INSStaggeredPPMConvectiveOperator::CONTEXT");

    const std::string U_var_name = "INSStaggeredPPMConvectiveOperator::U";
    d_U_var = var_db->getVariable(U_var_name);
    if (d_U_var.isNull())
    {
        d_U_var = new SAMRAI::pdat::SideVariable<NDIM,double>(U_var_name);
        d_U_scratch_idx = var_db->registerVariableAndContext(d_U_var, context, SAMRAI::hier::IntVector<NDIM>(GADVECTG));
    }
    else
    {
        d_U_scratch_idx = var_db->mapVariableAndContextToIndex(d_U_var, context);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_U_scratch_idx >= 0);
#endif
    return;
}// INSStaggeredPPMConvectiveOperator

INSStaggeredPPMConvectiveOperator::~INSStaggeredPPMConvectiveOperator()
{
    deallocateOperatorState();
    return;
}// ~INSStaggeredPPMConvectiveOperator

void
INSStaggeredPPMConvectiveOperator::applyConvectiveOperator(
    const int U_idx,
    const int N_idx)
{
    if (!d_is_initialized)
    {
        TBOX_ERROR("INSStaggeredPPMConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to applyConvectiveOperator\n");
    }

    // The timestep is set to equal zero since the data is already centered in
    // time.
    static const double dt = 0.0;

    // Setup communications schedules.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > refine_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_alg->registerRefine(d_U_scratch_idx, // destination
                               U_idx,           // source
                               d_U_scratch_idx, // temporary work space
                               d_refine_op);

    // Compute the convective derivative.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        refine_alg->resetSchedule(d_refine_scheds[ln]);
        d_refine_scheds[ln]->fillData(0.0);
        d_refine_alg->resetSchedule(d_refine_scheds[ln]);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

            const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::hier::IntVector<NDIM>& patch_lower = patch_box.lower();
            const SAMRAI::hier::IntVector<NDIM>& patch_upper = patch_box.upper();

            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > N_data = patch->getPatchData(N_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_data = patch->getPatchData(d_U_scratch_idx);

            const SAMRAI::hier::IntVector<NDIM> ghosts = SAMRAI::hier::IntVector<NDIM>(GADVECTG);
            SAMRAI::hier::Box<NDIM> side_boxes[NDIM];
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> >  U_adv_data[NDIM];
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > U_half_data[NDIM];
            for (int axis = 0; axis < NDIM; ++axis)
            {
                side_boxes [axis] = SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(patch_box,axis);
                U_adv_data [axis] = new SAMRAI::pdat::FaceData<NDIM,double>(side_boxes[axis],1,ghosts);
                U_half_data[axis] = new SAMRAI::pdat::FaceData<NDIM,double>(side_boxes[axis],1,ghosts);
            }
#if (NDIM == 2)
            NAVIER_STOKES_INTERP_COMPS_FC(
                patch_lower(0), patch_upper(0),
                patch_lower(1), patch_upper(1),
                U_data->getGhostCellWidth()(0),         U_data->getGhostCellWidth()(1),
                U_data->getPointer(0),                  U_data->getPointer(1),
                side_boxes[0].lower(0),                 side_boxes[0].upper(0),
                side_boxes[0].lower(1),                 side_boxes[0].upper(1),
                U_adv_data[0]->getGhostCellWidth()(0),  U_adv_data[0]->getGhostCellWidth()(1),
                U_adv_data[0]->getPointer(0),           U_adv_data[0]->getPointer(1),
                side_boxes[1].lower(0),                 side_boxes[1].upper(0),
                side_boxes[1].lower(1),                 side_boxes[1].upper(1),
                U_adv_data[1]->getGhostCellWidth()(0),  U_adv_data[1]->getGhostCellWidth()(1),
                U_adv_data[1]->getPointer(0),           U_adv_data[1]->getPointer(1));
#endif
#if (NDIM == 3)
            NAVIER_STOKES_INTERP_COMPS_FC(
                patch_lower(0), patch_upper(0),
                patch_lower(1), patch_upper(1),
                patch_lower(2), patch_upper(2),
                U_data->getGhostCellWidth()(0),         U_data->getGhostCellWidth()(1),         U_data->getGhostCellWidth()(2),
                U_data->getPointer(0),                  U_data->getPointer(1),                  U_data->getPointer(2),
                side_boxes[0].lower(0),                 side_boxes[0].upper(0),
                side_boxes[0].lower(1),                 side_boxes[0].upper(1),
                side_boxes[0].lower(2),                 side_boxes[0].upper(2),
                U_adv_data[0]->getGhostCellWidth()(0),  U_adv_data[0]->getGhostCellWidth()(1),  U_adv_data[0]->getGhostCellWidth()(2),
                U_adv_data[0]->getPointer(0),           U_adv_data[0]->getPointer(1),           U_adv_data[0]->getPointer(2),
                side_boxes[1].lower(0),                 side_boxes[1].upper(0),
                side_boxes[1].lower(1),                 side_boxes[1].upper(1),
                side_boxes[1].lower(2),                 side_boxes[1].upper(2),
                U_adv_data[1]->getGhostCellWidth()(0),  U_adv_data[1]->getGhostCellWidth()(1),  U_adv_data[1]->getGhostCellWidth()(2),
                U_adv_data[1]->getPointer(0),           U_adv_data[1]->getPointer(1),           U_adv_data[1]->getPointer(2),
                side_boxes[2].lower(0),                 side_boxes[2].upper(0),
                side_boxes[2].lower(1),                 side_boxes[2].upper(1),
                side_boxes[2].lower(2),                 side_boxes[2].upper(2),
                U_adv_data[2]->getGhostCellWidth()(0),  U_adv_data[2]->getGhostCellWidth()(1),  U_adv_data[2]->getGhostCellWidth()(2),
                U_adv_data[2]->getPointer(0),           U_adv_data[2]->getPointer(1),           U_adv_data[2]->getPointer(2));
#endif
            for (int axis = 0; axis < NDIM; ++axis)
            {
                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > dU_data =
                    new SAMRAI::pdat::SideData<NDIM,double>(U_data->getBox(), U_data->getDepth(), U_data->getGhostCellWidth());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_L_data =
                    new SAMRAI::pdat::SideData<NDIM,double>(U_data->getBox(), U_data->getDepth(), U_data->getGhostCellWidth());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_R_data =
                    new SAMRAI::pdat::SideData<NDIM,double>(U_data->getBox(), U_data->getDepth(), U_data->getGhostCellWidth());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_scratch1_data =
                    new SAMRAI::pdat::SideData<NDIM,double>(U_data->getBox(), U_data->getDepth(), U_data->getGhostCellWidth());
#if (NDIM == 3)
                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_scratch2_data =
                    new SAMRAI::pdat::SideData<NDIM,double>(U_data->getBox(), U_data->getDepth(), U_data->getGhostCellWidth());
#endif
                SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > U_half_scratch_data =
                    new SAMRAI::pdat::FaceData<NDIM,double>( U_half_data[axis]->getBox(), U_half_data[axis]->getDepth(), U_half_data[axis]->getGhostCellWidth());
#if (NDIM == 2)
                GODUNOV_PREDICT_FC(
                    dx, dt,
                    side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                    side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                    U_data->getGhostCellWidth()(0), U_data->getGhostCellWidth()(1),
                    U_data       ->getPointer(axis),       U_scratch1_data->getPointer(axis),
                    dU_data      ->getPointer(axis),       U_L_data       ->getPointer(axis),       U_R_data->getPointer(axis),
                    U_adv_data [axis]  ->getGhostCellWidth()(0), U_adv_data [axis]    ->getGhostCellWidth()(1),
                    U_half_data[axis]  ->getGhostCellWidth()(0), U_half_data[axis]    ->getGhostCellWidth()(1),
                    U_adv_data [axis]  ->getPointer(0),          U_adv_data [axis]    ->getPointer(1),
                    U_half_scratch_data->getPointer(0),          U_half_scratch_data  ->getPointer(1),
                    U_half_data[axis]  ->getPointer(0),          U_half_data[axis]    ->getPointer(1));
#endif
#if (NDIM == 3)
                bool using_full_ctu = true;
                GODUNOV_PREDICT_FC(
                    dx, dt,
                    static_cast<unsigned int>(using_full_ctu),
                    side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                    side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                    side_boxes[axis].lower(2), side_boxes[axis].upper(2),
                    U_data->getGhostCellWidth()(0), U_data->getGhostCellWidth()(1), U_data->getGhostCellWidth()(2),
                    U_data       ->getPointer(axis),       U_scratch1_data->getPointer(axis),       U_scratch2_data->getPointer(axis),
                    dU_data      ->getPointer(axis),       U_L_data       ->getPointer(axis),       U_R_data       ->getPointer(axis),
                    U_adv_data [axis]  ->getGhostCellWidth()(0), U_adv_data [axis]    ->getGhostCellWidth()(1), U_adv_data [axis]    ->getGhostCellWidth()(2),
                    U_half_data[axis]  ->getGhostCellWidth()(0), U_half_data[axis]    ->getGhostCellWidth()(1), U_half_data[axis]    ->getGhostCellWidth()(2),
                    U_adv_data [axis]  ->getPointer(0),          U_adv_data [axis]    ->getPointer(1),          U_adv_data [axis]    ->getPointer(2),
                    U_half_scratch_data->getPointer(0),          U_half_scratch_data  ->getPointer(1),          U_half_scratch_data  ->getPointer(2),
                    U_half_data[axis]  ->getPointer(0),          U_half_data[axis]    ->getPointer(1),          U_half_data[axis]    ->getPointer(2));
#endif
            }
#if (NDIM == 2)
            NAVIER_STOKES_RESET_ADV_VELOCITY_FC(
                side_boxes[0].lower(0), side_boxes[0].upper(0),
                side_boxes[0].lower(1), side_boxes[0].upper(1),
                U_adv_data [0]->getGhostCellWidth()(0), U_adv_data [0]->getGhostCellWidth()(1),
                U_adv_data [0]->getPointer(0),          U_adv_data [0]->getPointer(1),
                U_half_data[0]->getGhostCellWidth()(0), U_half_data[0]->getGhostCellWidth()(1),
                U_half_data[0]->getPointer(0),          U_half_data[0]->getPointer(1),
                side_boxes[1].lower(0), side_boxes[1].upper(0),
                side_boxes[1].lower(1), side_boxes[1].upper(1),
                U_adv_data [1]->getGhostCellWidth()(0), U_adv_data [1]->getGhostCellWidth()(1),
                U_adv_data [1]->getPointer(0),          U_adv_data [1]->getPointer(1),
                U_half_data[1]->getGhostCellWidth()(0), U_half_data[1]->getGhostCellWidth()(1),
                U_half_data[1]->getPointer(0),          U_half_data[1]->getPointer(1));
#endif
#if (NDIM == 3)
            NAVIER_STOKES_RESET_ADV_VELOCITY_FC(
                side_boxes[0].lower(0), side_boxes[0].upper(0),
                side_boxes[0].lower(1), side_boxes[0].upper(1),
                side_boxes[0].lower(2), side_boxes[0].upper(2),
                U_adv_data [0]->getGhostCellWidth()(0), U_adv_data [0]->getGhostCellWidth()(1), U_adv_data [0]->getGhostCellWidth()(2),
                U_adv_data [0]->getPointer(0),          U_adv_data [0]->getPointer(1),          U_adv_data [0]->getPointer(2),
                U_half_data[0]->getGhostCellWidth()(0), U_half_data[0]->getGhostCellWidth()(1), U_half_data[0]->getGhostCellWidth()(2),
                U_half_data[0]->getPointer(0),          U_half_data[0]->getPointer(1),          U_half_data[0]->getPointer(2),
                side_boxes[1].lower(0), side_boxes[1].upper(0),
                side_boxes[1].lower(1), side_boxes[1].upper(1),
                side_boxes[1].lower(2), side_boxes[1].upper(2),
                U_adv_data [1]->getGhostCellWidth()(0), U_adv_data [1]->getGhostCellWidth()(1), U_adv_data [1]->getGhostCellWidth()(2),
                U_adv_data [1]->getPointer(0),          U_adv_data [1]->getPointer(1),          U_adv_data [1]->getPointer(2),
                U_half_data[1]->getGhostCellWidth()(0), U_half_data[1]->getGhostCellWidth()(1), U_half_data[1]->getGhostCellWidth()(2),
                U_half_data[1]->getPointer(0),          U_half_data[1]->getPointer(1),          U_half_data[1]->getPointer(2),
                side_boxes[2].lower(0), side_boxes[2].upper(0),
                side_boxes[2].lower(1), side_boxes[2].upper(1),
                side_boxes[2].lower(2), side_boxes[2].upper(2),
                U_adv_data [2]->getGhostCellWidth()(0), U_adv_data [2]->getGhostCellWidth()(1), U_adv_data [2]->getGhostCellWidth()(2),
                U_adv_data [2]->getPointer(0),          U_adv_data [2]->getPointer(1),          U_adv_data [2]->getPointer(2),
                U_half_data[2]->getGhostCellWidth()(0), U_half_data[2]->getGhostCellWidth()(1), U_half_data[2]->getGhostCellWidth()(2),
                U_half_data[2]->getPointer(0),          U_half_data[2]->getPointer(1),          U_half_data[2]->getPointer(2));
#endif
            for (int axis = 0; axis < NDIM; ++axis)
            {
                if (d_conservation_form)
                {
#if (NDIM == 2)
                    CONVECT_DERIVATIVE_FC(
                        dx,
                        side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                        side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                        U_adv_data [axis]->getGhostCellWidth()(0), U_adv_data [axis]->getGhostCellWidth()(1),
                        U_half_data[axis]->getGhostCellWidth()(0), U_half_data[axis]->getGhostCellWidth()(1),
                        U_adv_data [axis]->getPointer(0),          U_adv_data [axis]->getPointer(1),
                        U_half_data[axis]->getPointer(0),          U_half_data[axis]->getPointer(1),
                        N_data->getGhostCellWidth()(0), N_data->getGhostCellWidth()(1),
                        N_data->getPointer(axis));
#endif
#if (NDIM == 3)
                    CONVECT_DERIVATIVE_FC(
                        dx,
                        side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                        side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                        side_boxes[axis].lower(2), side_boxes[axis].upper(2),
                        U_adv_data [axis]->getGhostCellWidth()(0), U_adv_data [axis]->getGhostCellWidth()(1), U_adv_data [axis]->getGhostCellWidth()(2),
                        U_half_data[axis]->getGhostCellWidth()(0), U_half_data[axis]->getGhostCellWidth()(1), U_half_data[axis]->getGhostCellWidth()(2),
                        U_adv_data [axis]->getPointer(0),          U_adv_data [axis]->getPointer(1),          U_adv_data [axis]->getPointer(2),
                        U_half_data[axis]->getPointer(0),          U_half_data[axis]->getPointer(1),          U_half_data[axis]->getPointer(2),
                        N_data->getGhostCellWidth()(0), N_data->getGhostCellWidth()(1), N_data->getGhostCellWidth()(2),
                        N_data->getPointer(axis));
#endif
                }
                else
                {
#if (NDIM == 2)
                    ADVECT_DERIVATIVE_FC(
                        dx,
                        side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                        side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                        U_adv_data [axis]->getGhostCellWidth()(0), U_adv_data [axis]->getGhostCellWidth()(1),
                        U_half_data[axis]->getGhostCellWidth()(0), U_half_data[axis]->getGhostCellWidth()(1),
                        U_adv_data [axis]->getPointer(0),          U_adv_data [axis]->getPointer(1),
                        U_half_data[axis]->getPointer(0),          U_half_data[axis]->getPointer(1),
                        N_data->getGhostCellWidth()(0), N_data->getGhostCellWidth()(1),
                        N_data->getPointer(axis));
#endif
#if (NDIM == 3)
                    ADVECT_DERIVATIVE_FC(
                        dx,
                        side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                        side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                        side_boxes[axis].lower(2), side_boxes[axis].upper(2),
                        U_adv_data [axis]->getGhostCellWidth()(0), U_adv_data [axis]->getGhostCellWidth()(1), U_adv_data [axis]->getGhostCellWidth()(2),
                        U_half_data[axis]->getGhostCellWidth()(0), U_half_data[axis]->getGhostCellWidth()(1), U_half_data[axis]->getGhostCellWidth()(2),
                        U_adv_data [axis]->getPointer(0),          U_adv_data [axis]->getPointer(1),          U_adv_data [axis]->getPointer(2),
                        U_half_data[axis]->getPointer(0),          U_half_data[axis]->getPointer(1),          U_half_data[axis]->getPointer(2),
                        N_data->getGhostCellWidth()(0), N_data->getGhostCellWidth()(1), N_data->getGhostCellWidth()(2),
                        N_data->getPointer(axis));
#endif
                }
            }
        }
    }
    return;
}// applyConvectiveOperator

void
INSStaggeredPPMConvectiveOperator::apply(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y)
{
    // Initialize the operator (if necessary).
    const bool deallocate_at_completion = !d_is_initialized;
    if (!d_is_initialized) initializeOperatorState(x,y);

    // Get the vector components.
    const int U_idx = x.getComponentDescriptorIndex(0);
    const int N_idx = y.getComponentDescriptorIndex(0);

    // Compute the action of the operator.
    applyConvectiveOperator(U_idx, N_idx);

    // Deallocate the operator (if necessary).
    if (deallocate_at_completion) deallocateOperatorState();
    return;
}// apply

void
INSStaggeredPPMConvectiveOperator::initializeOperatorState(
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& out)
{
    if (d_is_initialized) deallocateOperatorState();

    // Get the hierarchy configuration.
    d_hierarchy = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln = in.getFinestLevelNumber();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
#endif
    // Setup the refine algorithm, operator, patch strategy, and schedules.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    d_refine_op = grid_geom->lookupRefineOperator(d_U_var, "CONSERVATIVE_LINEAR_REFINE");
    d_refine_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    d_refine_alg->registerRefine(d_U_scratch_idx,                   // destination
                                 in.getComponentDescriptorIndex(0), // source
                                 d_U_scratch_idx,                   // temporary work space
                                 d_refine_op);
    d_refine_strategy = new IBTK::CartExtrapPhysBdryOp(d_U_scratch_idx, BDRY_EXTRAP_TYPE);
    d_refine_scheds.resize(d_finest_ln+1);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_refine_scheds[ln] = d_refine_alg->createSchedule(level, ln-1, d_hierarchy, d_refine_strategy);
    }

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_U_scratch_idx))
        {
            level->allocatePatchData(d_U_scratch_idx);
        }
    }
    d_is_initialized = true;
    return;
}// initializeOperatorState

void
INSStaggeredPPMConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_U_scratch_idx))
        {
            level->deallocatePatchData(d_U_scratch_idx);
        }
    }

    // Deallocate the refine algorithm, operator, patch strategy, and schedules.
    d_refine_op.setNull();
    d_refine_alg.setNull();
    d_refine_strategy.setNull();
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        d_refine_scheds[ln].setNull();
    }
    d_refine_scheds.clear();

    d_is_initialized = false;
    return;
}// deallocateOperatorState

void
INSStaggeredPPMConvectiveOperator::enableLogging(
    bool enabled)
{
    // intentionally blank
    return;
}// enableLogging

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::INSStaggeredPPMConvectiveOperator>;

//////////////////////////////////////////////////////////////////////////////
