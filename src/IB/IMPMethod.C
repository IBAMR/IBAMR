// Filename: IMPMethod.C
// Created on 16 Oct 2012 by Boyce Griffith
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

#include "IMPMethod.h"

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
#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/MaterialPointSpec.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IndexUtilities.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/LNodeSetData.h>
#include <ibtk/libmesh_utilities.h>

// LIBMESH INCLUDES
using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{

#if 1
static const int kernel_width = 2;

inline double
kernel(
    double x)
{
    x += 2.;
    const double x2 = x*x;
    const double x3 = x*x2;
    if (x <= 0.)
        return 0.;
    else if (x <= 1.)
        return .1666666666666667*x3;
    else if (x <= 2.)
        return 2.*x2-.5000000000000000*x3-2.*x+.6666666666666667;
    else if (x <= 3.)
        return 10.*x-4.*x2+.5000000000000000*x3-7.333333333333333;
    else if (x <= 4.)
        return 10.66666666666667-8.*x+2.*x2-.1666666666666667*x3;
    else
        return 0.;
}// kernel

inline double
kernel_diff(
    double x)
{
    x += 2.;
    const double x2 = x*x;
    if (x <= 0.)
        return 0.;
    else if (x <= 1.)
        return .5000000000000000*x2;
    else if (x <= 2.)
        return 4.*x-1.500000000000000*x2-2.;
    else if (x <= 3.)
        return 10.-8.*x+1.500000000000000*x2;
    else if (x <= 4.)
        return -8.+4.*x-.5000000000000000*x2;
    else
        return 0.;
}// kernel_diff
#endif

#if 0
static const int kernel_width = 3;

inline double
kernel(
    double x)
{
    x += 3.;
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x*x3;
    const double x5 = x*x4;
    if (x <= 0.)
        return 0.;
    else if (x <= 1.)
        return .8333333333333333e-2*x5;
    else if (x <= 2.)
        return .2500000000000000*x4-.4166666666666667e-1*x5-.5000000000000000*x3+.5000000000000000*x2-.2500000000000000*x+.5000000000000000e-1;
    else if (x <= 3.)
        return 4.500000000000000*x3-1.*x4+.8333333333333333e-1*x5-9.500000000000000*x2+9.750000000000000*x-3.950000000000000;
    else if (x <= 4.)
        return 35.50000000000000*x2-10.50000000000000*x3+1.500000000000000*x4-.8333333333333333e-1*x5-57.75000000000000*x+36.55000000000000;
    else if (x <= 5.)
        return 102.2500000000000*x-44.50000000000000*x2+9.500000000000000*x3-1.*x4+.4166666666666667e-1*x5-91.45000000000000;
    else if (x <= 6.)
        return 64.80000000000000-54.*x+18.*x2-3.*x3+.2500000000000000*x4-.8333333333333333e-2*x5;
    else
        return 0.;
}// kernel

inline double
kernel_diff(
    double x)
{
    x += 3.;
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x*x3;
    if (x < 0.)
        return 0.;
    else if (x <= 1.)
        return .4166666666666667e-1*x4;
    else if (x <= 2.)
        return -.2500000000000000+x-1.500000000000000*x2+x3-.2083333333333333*x4;
    else if (x <= 3.)
        return 9.750000000000000-19.*x+13.50000000000000*x2-4.*x3+.4166666666666667*x4;
    else if (x <= 4.)
        return -57.75000000000000+71.*x-31.50000000000000*x2+6.*x3-.4166666666666667*x4;
    else if (x <= 5.)
        return 102.2500000000000-89.*x+28.50000000000000*x2-4.*x3+.2083333333333333*x4;
    else if (x <= 6.)
        return -54.+36.*x-9.*x2+x3-.4166666666666667e-1*x4;
    else
        return 0.;
}// kernel_diff
#endif

#if 0
static const int kernel_width = 4;

inline double
kernel(
    double x)
{
    x += 4.;
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x*x3;
    const double x5 = x*x4;
    const double x6 = x*x5;
    const double x7 = x*x6;
    if (x <= 0.)
        return 0.;
    else if (x <= 1.)
        return .1984126984126984e-3*x7;
    else if (x <= 2.)
        return .1111111111111111e-1*x6-.1388888888888889e-2*x7-.3333333333333333e-1*x5+.5555555555555556e-1*x4-.5555555555555556e-1*x3+.3333333333333333e-1*x2-.1111111111111111e-1*x+.1587301587301587e-2;
    else if (x <= 3.)
        return .4333333333333333*x5-.6666666666666667e-1*x6+.4166666666666667e-2*x7-1.500000000000000*x4+3.055555555555556*x3-3.700000000000000*x2+2.477777777777778*x-.7095238095238095;
    else if (x <= 4.)
        return 9.*x4-1.666666666666667*x5+.1666666666666667*x6-.6944444444444444e-2*x7-28.44444444444444*x3+53.*x2-54.22222222222222*x+23.59047619047619;
    else if (x <= 5.)
        return 96.*x3-22.11111111111111*x4+3.*x5-.2222222222222222*x6+.6944444444444444e-2*x7-245.6666666666667*x2+344.*x-203.9650793650794;
    else if (x <= 6.)
        return 483.5000000000000*x2-147.0555555555556*x3+26.50000000000000*x4-2.833333333333333*x5+.1666666666666667*x6-.4166666666666667e-2*x7-871.2777777777778*x+664.0904761904762;
    else if (x <= 7.)
        return 943.1222222222222*x-423.7000000000000*x2+104.9444444444444*x3-15.50000000000000*x4+1.366666666666667*x5-.6666666666666667e-1*x6+.1388888888888889e-2*x7-891.1095238095238;
    else if (x <= 8.)
        return 416.1015873015873-364.0888888888889*x+136.5333333333333*x2-28.44444444444444*x3+3.555555555555556*x4-.2666666666666667*x5+.1111111111111111e-1*x6-.1984126984126984e-3*x7;
    else
        return 0.;
}// kernel

inline double
kernel_diff(
    double x)
{
    x += 4.;
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x*x3;
    const double x5 = x*x4;
    const double x6 = x*x5;
    if (x <= 0.)
        return 0.;
    else if (x <= 1.)
        return .1388888888888889e-2*x6;
    else if (x <= 2.)
        return -.1111111111111111e-1+.6666666666666667e-1*x-.1666666666666667*x2+.2222222222222222*x3-.1666666666666667*x4+.6666666666666667e-1*x5-.9722222222222222e-2*x6;
    else if (x <= 3.)
        return 2.477777777777778-7.400000000000000*x+9.166666666666667*x2-6.*x3-.4000000000000000*x5+.2916666666666667e-1*x6+2.166666666666667*x4;
    else if (x <= 4.)
        return -54.22222222222222+106.*x-85.33333333333333*x2+36.*x3-8.333333333333333*x4+x5-.4861111111111111e-1*x6;
    else if (x <= 5.)
        return 344.-491.3333333333333*x+288.*x2-88.44444444444444*x3+15.*x4-1.333333333333333*x5+.4861111111111111e-1*x6;
    else if (x <= 6.)
        return -871.2777777777778+967.*x-441.1666666666667*x2+106.*x3-14.16666666666667*x4+x5-.2916666666666667e-1*x6;
    else if (x <= 7.)
        return 943.1222222222222-847.4000000000000*x+314.8333333333333*x2-62.*x3+6.833333333333333*x4-.4000000000000000*x5+.9722222222222222e-2*x6;
    else if (x <= 8.)
        return -364.0888888888889+273.0666666666667*x-85.33333333333333*x2+14.22222222222222*x3-1.333333333333333*x4+.6666666666666667e-1*x5-.1388888888888889e-2*x6;
    else
        return 0.;
}// kernel_diff
#endif

// Version of IMPMethod restart file data.
static const int IMP_METHOD_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IMPMethod::IMPMethod(
    const std::string& object_name,
    Pointer<Database> input_db,
    bool register_for_restart)
{
    // Set the object name and register it with the restart manager.
    d_object_name = object_name;
    d_registered_for_restart = false;
    if (register_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
        d_registered_for_restart = true;
    }

    // Ensure all pointers to helper objects are NULL.
    d_l_initializer = NULL;
    d_silo_writer = NULL;

    // Set some default values.
    d_ghosts = kernel_width+1;
    d_do_log = false;

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);

    // Get the Lagrangian Data Manager.
    LEInteractor::s_delta_fcn = &kernel;
    LEInteractor::s_delta_fcn_stencil_size = 2*kernel_width;
    LEInteractor::s_delta_fcn_C = std::numeric_limits<double>::quiet_NaN();
    d_l_data_manager = LDataManager::getManager(d_object_name+"::LDataManager", "USER_DEFINED", "USER_DEFINED", d_ghosts, d_registered_for_restart);
    d_ghosts = d_l_data_manager->getGhostCellWidth();

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time     = std::numeric_limits<double>::quiet_NaN();
    d_half_time    = std::numeric_limits<double>::quiet_NaN();

    // Indicate all Lagrangian data needs ghost values to be refilled, and that
    // all intermediate data needs to be initialized.
    d_X_current_needs_ghost_fill = true;
    d_X_new_needs_ghost_fill     = true;
    d_X_half_needs_ghost_fill    = true;
    d_X_half_needs_reinit        = true;
    d_U_half_needs_reinit        = true;
    return;
}// IMPMethod

IMPMethod::~IMPMethod()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
}// ~IMPMethod

void
IMPMethod::registerLInitStrategy(
    Pointer<LInitStrategy> l_initializer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!l_initializer.isNull());
#endif
    d_l_initializer = l_initializer;
    d_l_data_manager->registerLInitStrategy(d_l_initializer);
    return;
}// registerLInitStrategy

void
IMPMethod::freeLInitStrategy()
{
    d_l_initializer.setNull();
    d_l_data_manager->freeLInitStrategy();
    return;
}// freeLInitStrategy

LDataManager*
IMPMethod::getLDataManager() const
{
    return d_l_data_manager;
}// getLDataManager

void
IMPMethod::registerLSiloDataWriter(
    Pointer<LSiloDataWriter> silo_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!silo_writer.isNull());
#endif
    d_silo_writer = silo_writer;
    d_l_data_manager->registerLSiloDataWriter(d_silo_writer);
    return;
}// registerLSiloDataWriter

const IntVector<NDIM>&
IMPMethod::getMinimumGhostCellWidth() const
{
    return d_ghosts;
}// getMinimumGhostCellWidth

void
IMPMethod::setupTagBuffer(
    Array<int>& tag_buffer,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg) const
{
    const int finest_hier_ln = gridding_alg->getMaxLevels()-1;
    const int tsize = tag_buffer.size();
    tag_buffer.resizeArray(finest_hier_ln);
    for (int i = tsize; i < finest_hier_ln; ++i) tag_buffer[i] = 0;
    const int gcw = d_ghosts.max();
    for (int tag_ln = 0; tag_ln < finest_hier_ln; ++tag_ln)
    {
        const int data_ln = tag_ln+1;
        const int can_be_refined = data_ln < finest_hier_ln;
        if (!d_l_initializer->getLevelHasLagrangianData(data_ln, can_be_refined)) continue;
        tag_buffer[tag_ln] = std::max(tag_buffer[tag_ln], gcw);
    }
    for (int ln = finest_hier_ln-2; ln >= 0; --ln)
    {
        tag_buffer[ln] = std::max(tag_buffer[ln], tag_buffer[ln+1]/gridding_alg->getRatioToCoarserLevel(ln+1).max()+1);
    }
    return;
}// setupTagBuffer

void
IMPMethod::preprocessIntegrateData(
    double current_time,
    double new_time,
    int /*num_cycles*/)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_half_time = current_time+0.5*(new_time-current_time);

    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Look-up or allocate Lagangian data.
    d_X_current_data     .resize(finest_ln+1);
    d_X_new_data         .resize(finest_ln+1);
    d_X_half_data        .resize(finest_ln+1);
    d_X0_data            .resize(finest_ln+1);
    d_U_current_data     .resize(finest_ln+1);
    d_U_new_data         .resize(finest_ln+1);
    d_U_half_data        .resize(finest_ln+1);
    d_Grad_U_current_data.resize(finest_ln+1);
    d_Grad_U_new_data    .resize(finest_ln+1);
    d_Grad_U_half_data   .resize(finest_ln+1);
    d_F_current_data     .resize(finest_ln+1);
    d_F_new_data         .resize(finest_ln+1);
    d_F_half_data        .resize(finest_ln+1);
    d_tau_data           .resize(finest_ln+1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        d_X_current_data     [ln] = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME,ln);
        d_X_new_data         [ln] = d_l_data_manager->createLData("X_new",ln,NDIM);
        d_X0_data            [ln] = d_l_data_manager->getLData("X0",ln);
        d_U_current_data     [ln] = d_l_data_manager->getLData(LDataManager:: VEL_DATA_NAME,ln);
        d_U_new_data         [ln] = d_l_data_manager->createLData("U_new",ln,NDIM);
        d_Grad_U_current_data[ln] = d_l_data_manager->getLData("Grad_U",ln);
        d_Grad_U_new_data    [ln] = d_l_data_manager->createLData("Grad_U_new",ln,NDIM*NDIM);
        d_F_current_data     [ln] = d_l_data_manager->getLData("F",ln);
        d_F_new_data         [ln] = d_l_data_manager->createLData("F_new",ln,NDIM*NDIM);
        d_F_half_data        [ln] = d_l_data_manager->createLData("F_half",ln,NDIM*NDIM);
        d_tau_data           [ln] = d_l_data_manager->getLData("tau",ln);

        // Initialize new values to equal current values.
        ierr = VecCopy(     d_X_current_data[ln]->getVec(),      d_X_new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecCopy(     d_U_current_data[ln]->getVec(),      d_U_new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecCopy(d_Grad_U_current_data[ln]->getVec(), d_Grad_U_new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecCopy(     d_F_current_data[ln]->getVec(),      d_F_new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecCopy(     d_F_current_data[ln]->getVec(),     d_F_half_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
    }

    // Keep track of Lagrangian data objects that need to have ghost values
    // filled, or that need to be reinitialized.
    d_X_new_needs_ghost_fill = true;
    d_X_half_needs_reinit    = true;
    d_U_half_needs_reinit    = true;
    return;
}// preprocessIntegrateData

void
IMPMethod::postprocessIntegrateData(
    double /*current_time*/,
    double /*new_time*/,
    int /*num_cycles*/)
{
    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Reset time-dependent Lagrangian data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        ierr = VecSwap(     d_X_current_data[ln]->getVec(),      d_X_new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecSwap(     d_U_current_data[ln]->getVec(),      d_U_new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecSwap(d_Grad_U_current_data[ln]->getVec(), d_Grad_U_new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecSwap(     d_F_current_data[ln]->getVec(),      d_F_new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
    }
    d_X_current_needs_ghost_fill = true;

    // Deallocate Lagrangian scratch data.
    d_X_current_data     .clear();
    d_X_new_data         .clear();
    d_X_half_data        .clear();
    d_X0_data            .clear();
    d_U_current_data     .clear();
    d_U_new_data         .clear();
    d_U_half_data        .clear();
    d_Grad_U_current_data.clear();
    d_Grad_U_new_data    .clear();
    d_Grad_U_half_data   .clear();
    d_F_current_data     .clear();
    d_F_new_data         .clear();
    d_F_half_data        .clear();
    d_tau_data           .clear();

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time     = std::numeric_limits<double>::quiet_NaN();
    d_half_time    = std::numeric_limits<double>::quiet_NaN();
    return;
}// postprocessIntegrateData

void
IMPMethod::interpolateVelocity(
    const int u_data_idx,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
    const double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    std::vector<Pointer<LData> >* U_data, * Grad_U_data, * X_data;
    bool* X_needs_ghost_fill;
    getVelocityData(&U_data, &Grad_U_data, data_time);
    getPositionData(&X_data, &X_needs_ghost_fill, data_time);

    // Synchronize Eulerian and Lagrangian values.
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        if (ln < static_cast<int>(u_synch_scheds.size()) && u_synch_scheds[ln])
        {
            u_synch_scheds[ln]->coarsenData();
        }
        if (*X_needs_ghost_fill && (*X_data)[ln]) (*X_data)[ln]->beginGhostUpdate();
    }
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        if (*X_needs_ghost_fill && (*X_data)[ln]) (*X_data)[ln]->endGhostUpdate();
    }
    *X_needs_ghost_fill = false;

    // Interpolate data from the Eulerian grid to the Lagrangian mesh.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        if (ln < static_cast<int>(u_ghost_fill_scheds.size()) && u_ghost_fill_scheds[ln])
        {
            u_ghost_fill_scheds[ln]->fillData(data_time);
        }
        double*      U_array = (*     U_data)[ln]->getLocalFormVecArray()->data();
        double* Grad_U_array = (*Grad_U_data)[ln]->getLocalFormVecArray()->data();
        double*      X_array = (*     X_data)[ln]->getLocalFormVecArray()->data();
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM,double> > u_data = patch->getPatchData(u_data_idx);
            Pointer<LNodeSetData> idx_data = patch->getPatchData(d_l_data_manager->getLNodePatchDescriptorIndex());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const x_lower = patch_geom->getXLower();
            const double* const x_upper = patch_geom->getXUpper();
            const double* const dx = patch_geom->getDx();
            Box<NDIM> side_boxes[NDIM];
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                side_boxes[axis] = SideGeometry<NDIM>::toSideBox(u_data->getGhostBox()*idx_data->getGhostBox(),axis);
                if (patch_geom->getTouchesRegularBoundary(axis, /*lower*/ 0)) side_boxes[axis].lower(axis) = patch_box.lower(axis);
                if (patch_geom->getTouchesRegularBoundary(axis, /*upper*/ 1)) side_boxes[axis].upper(axis) = patch_box.upper(axis)+1;
            }
            for (LNodeSetData::CellIterator it(idx_data->getGhostBox()); it; it++)
            {
                const Index<NDIM>& i = *it;
                LNodeSet* const node_set = idx_data->getItem(i);
                if (!node_set) continue;
                for (LNodeSet::iterator it = node_set->begin(); it != node_set->end(); ++it)
                {
                    const LNode* const node_idx = *it;
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const double* const X = &X_array[NDIM*local_idx];

                    // WARNING: As written here, this implicitly imposes u = 0 in
                    // the ghost cell region at physical boundaries.
                    const Index<NDIM> i = IndexUtilities::getCellIndex(X, x_lower, x_upper, dx, patch_box.lower(), patch_box.upper());
                    VectorValue<double> U, X_cell;
                    TensorValue<double> Grad_U;
                    U.zero();
                    Grad_U.zero();
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        X_cell(d) = x_lower[d] + dx[d]*(static_cast<double>(i(d)-patch_box.lower(d))+0.5);
                    }
                    for (unsigned int component = 0; component < NDIM; ++component)
                    {
                        blitz::Array<double,1> phi[NDIM], dphi[NDIM];
                        Box<NDIM> box(i,i);
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            if (d == component)
                            {
                                box.lower(d) -= (kernel_width-1);
                                box.upper(d) += (kernel_width  );
                            }
                            else
                            {
                                box.lower(d) -= (X[d] <= X_cell(d) ? kernel_width : kernel_width-1);
                                box.upper(d) += (X[d] >= X_cell(d) ? kernel_width : kernel_width-1);
                            }
                            phi [d].resize(blitz::Range(box.lower(d),box.upper(d)));
                            dphi[d].resize(blitz::Range(box.lower(d),box.upper(d)));
                        }
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            for (int i = box.lower(d); i <= box.upper(d); ++i)
                            {
                                const double x_grid = x_lower[d] + dx[d]*(static_cast<double>(i-patch_box.lower(d))+(d == component ? 0.0 : 0.5));
                                const double del = x_grid - X[d];
                                phi [d](i) = kernel(del/dx[d]);
                                dphi[d](i) = kernel_diff(del/dx[d])/dx[d];
                            }
                        }
                        for (Box<NDIM>::Iterator b(box*side_boxes[component]); b; b++)
                        {
                            const Index<NDIM>& i = b();
                            const double u = (*u_data)(SideIndex<NDIM>(i, component, SideIndex<NDIM>::Lower));
                            double w = 1.0;
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                w *= phi[d](i(d));
                            }
                            U(component) += u*w;
                            for (unsigned int k = 0; k < NDIM; ++k)
                            {
                                double dw_dx_k = 1.0;
                                for (unsigned int d = 0; d < NDIM; ++d)
                                {
                                    if (d == k) dw_dx_k *= dphi[d](i(d));
                                    else        dw_dx_k *=  phi[d](i(d));
                                }
                                Grad_U(component,k) -= u*dw_dx_k;
                            }
                        }
                    }
                    for (int i = 0; i < NDIM; ++i)
                    {
                        U_array[NDIM*local_idx+i] = U(i);
                        for (int j = 0; j < NDIM; ++j)
                        {
                            Grad_U_array[NDIM*NDIM*local_idx+NDIM*i+j] = Grad_U(i,j);
                        }
                    }
                }
            }
        }
        (*     U_data)[ln]->restoreArrays();
        (*Grad_U_data)[ln]->restoreArrays();
        (*     X_data)[ln]->restoreArrays();
    }
    d_U_half_needs_reinit = !MathUtilities<double>::equalEps(data_time, d_half_time);
    return;
}// interpolateVelocity

void
IMPMethod::eulerStep(
    const double current_time,
    const double new_time)
{
    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;
    std::vector<Pointer<LData> >* U_data, * Grad_U_data;
    getVelocityData(&U_data, &Grad_U_data, current_time);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        // Update the positions.
        ierr = VecWAXPY(d_X_new_data[ln]->getVec(), dt, (*U_data)[ln]->getVec(), d_X_current_data[ln]->getVec());  IBTK_CHKERRQ(ierr);

        // Update the deformation gradient.
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        blitz::Array<double,2>& F_current_array = *d_F_current_data[ln]->getVecArray();
        blitz::Array<double,2>& F_new_array     = *d_F_new_data    [ln]->getVecArray();
        blitz::Array<double,2>& F_half_array    = *d_F_half_data   [ln]->getVecArray();
        blitz::Array<double,2>& Grad_U_array = *(*Grad_U_data)[ln]->getVecArray();
        TensorValue<double> F_current, F_new, F_half, Grad_U;
        TensorValue<double> I(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int idx = node_idx->getGlobalPETScIndex();
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    F_current(i,j) = F_current_array(idx,NDIM*i+j);
                    Grad_U   (i,j) = Grad_U_array   (idx,NDIM*i+j);
                }
            }
#if (NDIM == 2)
            F_current(2,2) = 1.0;
#endif
            F_new  = tensor_inverse(I-0.50*dt*Grad_U)*(I+0.50*dt*Grad_U)*F_current;
            F_half = tensor_inverse(I-0.25*dt*Grad_U)*(I+0.25*dt*Grad_U)*F_current;
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    F_new_array (idx,NDIM*i+j) = F_new (i,j);
                    F_half_array(idx,NDIM*i+j) = F_half(i,j);
                }
            }
        }
    }
    d_X_new_needs_ghost_fill = true;
    d_X_half_needs_reinit    = true;
    return;
}// eulerStep

void
IMPMethod::midpointStep(
    const double current_time,
    const double new_time)
{
    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;
    std::vector<Pointer<LData> >* U_data, * Grad_U_data;
    getVelocityData(&U_data, &Grad_U_data, current_time+0.5*dt);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        // Update the positions.
        ierr = VecWAXPY(d_X_new_data[ln]->getVec(), dt, (*U_data)[ln]->getVec(), d_X_current_data[ln]->getVec());  IBTK_CHKERRQ(ierr);

        // Update the deformation gradient.
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        blitz::Array<double,2>& F_current_array = *d_F_current_data[ln]->getVecArray();
        blitz::Array<double,2>& F_new_array     = *d_F_new_data    [ln]->getVecArray();
        blitz::Array<double,2>& Grad_U_array = *(*Grad_U_data)[ln]->getVecArray();
        TensorValue<double> F_current, F_new, Grad_U;
        TensorValue<double> I(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int idx = node_idx->getGlobalPETScIndex();
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    F_current(i,j) = F_current_array(idx,NDIM*i+j);
                    Grad_U   (i,j) = Grad_U_array   (idx,NDIM*i+j);
                }
            }
#if (NDIM == 2)
            F_current(2,2) = 1.0;
#endif
            F_new = tensor_inverse(I-0.5*dt*Grad_U)*(I+0.5*dt*Grad_U)*F_current;
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    F_new_array(idx,NDIM*i+j) = F_new(i,j);
                }
            }
        }
    }
    d_X_new_needs_ghost_fill = true;
    d_X_half_needs_reinit    = true;
    return;
}// midpointStep

void
IMPMethod::trapezoidalStep(
    const double current_time,
    const double new_time)
{
    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;
    std::vector<Pointer<LData> >* U_current_data, * U_new_data, * Grad_U_current_data, * Grad_U_new_data;
    getVelocityData(&U_current_data, &Grad_U_current_data, current_time);
    getVelocityData(&U_new_data    , &Grad_U_new_data    ,     new_time);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        // Update the positions.
        ierr = VecWAXPY(d_X_new_data[ln]->getVec(), 0.5*dt, (*U_current_data)[ln]->getVec(), d_X_current_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecAXPY( d_X_new_data[ln]->getVec(), 0.5*dt, (*U_new_data    )[ln]->getVec()                                );  IBTK_CHKERRQ(ierr);

        // Update the deformation gradient.
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        blitz::Array<double,2>& F_current_array = *d_F_current_data[ln]->getVecArray();
        blitz::Array<double,2>& F_new_array     = *d_F_new_data    [ln]->getVecArray();
        blitz::Array<double,2>& Grad_U_current_array = *(*Grad_U_current_data)[ln]->getVecArray();
        blitz::Array<double,2>& Grad_U_new_array     = *(*Grad_U_new_data    )[ln]->getVecArray();
        TensorValue<double> F_current, F_new, Grad_U_current, Grad_U_new;
        TensorValue<double> I(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int idx = node_idx->getGlobalPETScIndex();
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    F_current     (i,j) = F_current_array     (idx,NDIM*i+j);
                    Grad_U_current(i,j) = Grad_U_current_array(idx,NDIM*i+j);
                    Grad_U_new    (i,j) = Grad_U_new_array    (idx,NDIM*i+j);
                }
            }
#if (NDIM == 2)
            F_current(2,2) = 1.0;
#endif
            F_new = tensor_inverse(I-0.5*dt*Grad_U_new)*(I+0.5*dt*Grad_U_current)*F_current;
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    F_new_array(idx,NDIM*i+j) = F_new(i,j);
                }
            }
        }
    }
    d_X_new_needs_ghost_fill = true;
    d_X_half_needs_reinit    = true;
    return;
}// trapezoidalStep

void
IMPMethod::computeLagrangianForce(
    const double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    std::vector<Pointer<LData> >* X_data, * F_data;
    bool* X_needs_ghost_fill;
    getPositionData(&X_data, &X_needs_ghost_fill, data_time);
    getDeformationGradientData(&F_data, data_time);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        blitz::Array<double,2>&   x_array =  *(*X_data)[ln]->getVecArray();
        blitz::Array<double,2>&   X_array =  *d_X0_data[ln]->getVecArray();
        blitz::Array<double,2>&   F_array =  *(*F_data)[ln]->getVecArray();
        blitz::Array<double,2>& tau_array = *d_tau_data[ln]->getVecArray();
        TensorValue<double> FF, PP, tau;
        VectorValue<double> X, x;
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int idx = node_idx->getGlobalPETScIndex();
            MaterialPointSpec* mp_spec = node_idx->getNodeDataItem<MaterialPointSpec>();
            if (mp_spec)
            {
                for (int i = 0; i < NDIM; ++i)
                {
                    for (int j = 0; j < NDIM; ++j)
                    {
                        FF(i,j) = F_array(idx,NDIM*i+j);
                    }
                    x(i) = x_array(idx,i);
                    X(i) = X_array(idx,i);
                }
#if (NDIM == 2)
                FF(2,2) = 1.0;
#endif
                (*d_PK1_stress_fcn)(PP, FF, x, X, mp_spec->getSubdomainId(), mp_spec->getInternalVariables(), data_time, d_PK1_stress_fcn_ctx);
                tau = PP*FF.transpose();
            }
            else
            {
                tau.zero();
            }
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    tau_array(idx,NDIM*i+j) = tau(i,j);
                }
            }
        }
    }
    return;
}// computeLagrangianForce

void
IMPMethod::spreadForce(
    const int f_data_idx,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*f_prolongation_scheds*/,
    const double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    std::vector<Pointer<LData> >* X_data;
    bool* X_needs_ghost_fill;
    getPositionData(&X_data, &X_needs_ghost_fill, data_time);

    // Synchronize Lagrangian values.
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        if (*X_needs_ghost_fill && (*X_data)[ln]) (*X_data)[ln]->beginGhostUpdate();
        if (d_tau_data[ln]) d_tau_data[ln]->beginGhostUpdate();
    }
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        if (*X_needs_ghost_fill && (*X_data)[ln]) (*X_data)[ln]->endGhostUpdate();
        if (d_tau_data[ln]) d_tau_data[ln]->endGhostUpdate();
    }
    *X_needs_ghost_fill = false;

    // Spread data from the Lagrangian mesh to the Eulerian grid.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        blitz::Array<double,2>&   X_array = *( *X_data)[ln]->getGhostedLocalFormVecArray();
        blitz::Array<double,2>& tau_array = *d_tau_data[ln]->getGhostedLocalFormVecArray();
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM,double> > f_data = patch->getPatchData(f_data_idx);
            Pointer<LNodeSetData> idx_data = patch->getPatchData(d_l_data_manager->getLNodePatchDescriptorIndex());
            const Box<NDIM>& patch_box = patch->getBox();
            Box<NDIM> side_boxes[NDIM];
            for (unsigned int d = 0; d < NDIM; ++d) side_boxes[d] = SideGeometry<NDIM>::toSideBox(patch_box,d);
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const x_lower = patch_geom->getXLower();
            const double* const x_upper = patch_geom->getXUpper();
            const double* const dx = patch_geom->getDx();
            double dV_c = 1.0; for (unsigned int d = 0; d < NDIM; ++d) dV_c *= dx[d];
            for (LNodeSetData::CellIterator it(idx_data->getGhostBox()); it; it++)
            {
                const Index<NDIM>& i = *it;
                LNodeSet* const node_set = idx_data->getItem(i);
                if (!node_set) continue;
                for (LNodeSet::iterator it = node_set->begin(); it != node_set->end(); ++it)
                {
                    const LNode* const node_idx = *it;
                    MaterialPointSpec* mp_spec = node_idx->getNodeDataItem<MaterialPointSpec>();
                    if (!mp_spec) continue;
                    const double wgt = mp_spec->getWeight();
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const double* const X = &X_array(local_idx,0);
                    TensorValue<double> tau;
                    for (int i = 0; i < NDIM; ++i)
                    {
                        for (int j = 0; j < NDIM; ++j)
                        {
                            tau(i,j) = tau_array(local_idx,NDIM*i+j);
                        }
                    }

                    // Weight tau using a smooth kernel function evaluated about
                    // X.
                    const Index<NDIM> i = IndexUtilities::getCellIndex(X, x_lower, x_upper, dx, patch_box.lower(), patch_box.upper());
                    VectorValue<double> X_cell;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        X_cell(d) = x_lower[d] + dx[d]*(static_cast<double>(i(d)-patch_box.lower(d))+0.5);
                    }
                    for (unsigned int component = 0; component < NDIM; ++component)
                    {
                        blitz::Array<double,1> phi[NDIM], dphi[NDIM];
                        Box<NDIM> box(i,i);
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            if (d == component)
                            {
                                box.lower(d) -= (kernel_width-1);
                                box.upper(d) += (kernel_width  );
                            }
                            else
                            {
                                box.lower(d) -= (X[d] <= X_cell(d) ? kernel_width : kernel_width-1);
                                box.upper(d) += (X[d] >= X_cell(d) ? kernel_width : kernel_width-1);
                            }
                            phi [d].resize(blitz::Range(box.lower(d),box.upper(d)));
                            dphi[d].resize(blitz::Range(box.lower(d),box.upper(d)));
                        }
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            for (int i = box.lower(d); i <= box.upper(d); ++i)
                            {
                                const double x_grid = x_lower[d] + dx[d]*(static_cast<double>(i-patch_box.lower(d))+(d == component ? 0.0 : 0.5));
                                const double del = x_grid - X[d];
                                phi [d](i) = kernel(del/dx[d]);
                                dphi[d](i) = kernel_diff(del/dx[d])/dx[d];
                            }
                        }
                        for (Box<NDIM>::Iterator b(box*side_boxes[component]); b; b++)
                        {
                            const Index<NDIM>& i = b();
                            double f = 0.0;
                            for (unsigned int k = 0; k < NDIM; ++k)
                            {
                                double dw_dx_k = 1.0;
                                for (unsigned int d = 0; d < NDIM; ++d)
                                {
                                    if (d == k) dw_dx_k *= dphi[d](i(d));
                                    else        dw_dx_k *=  phi[d](i(d));
                                }
                                f += tau(component,k)*dw_dx_k;
                            }
                            (*f_data)(SideIndex<NDIM>(i, component, SideIndex<NDIM>::Lower)) += f*wgt/dV_c;
                        }
                    }
                }
            }
        }
    }
    return;
}// spreadForce

void
IMPMethod::initializePatchHierarchy(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
    int /*u_data_idx*/,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*u_ghost_fill_scheds*/,
    int /*integrator_step*/,
    double /*init_data_time*/,
    bool initial_time)
{
    // Cache pointers to the patch hierarchy and gridding algorithm.
    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;

    // Initialize various Lagrangian data objects.
    if (initial_time)
    {
        pout << "WARNING: IMPMethod implementation currently requires that the initial velocity is *zero*.\n";
    }
    return;
}// initializePatchHierarchy

void
IMPMethod::registerPK1StressTensorFunction(
    PK1StressFcnPtr PK1_stress_fcn,
    void* PK1_stress_fcn_ctx)
{
    d_PK1_stress_fcn = PK1_stress_fcn;
    d_PK1_stress_fcn_ctx = PK1_stress_fcn_ctx;
    return;
}// registerPK1StressTensorFunction

void
IMPMethod::registerLoadBalancer(
    Pointer<LoadBalancer<NDIM> > load_balancer,
    int workload_data_idx)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!load_balancer.isNull());
#endif
    d_load_balancer = load_balancer;
    d_workload_idx = workload_data_idx;
    d_l_data_manager->registerLoadBalancer(load_balancer, workload_data_idx);
    return;
}// registerLoadBalancer

void
IMPMethod::updateWorkloadEstimates(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    int /*workload_data_idx*/)
{
    d_l_data_manager->updateWorkloadEstimates();
    return;
}// updateWorkloadEstimates

void
IMPMethod::beginDataRedistribution(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    d_l_data_manager->beginDataRedistribution();
    return;
}// beginDataRedistribution

void
IMPMethod::endDataRedistribution(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    d_l_data_manager->endDataRedistribution();
    return;
}// endDataRedistribution

void
IMPMethod::initializeLevelData(
    Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    int level_number,
    double init_data_time,
    bool can_be_refined,
    bool initial_time,
    Pointer<BasePatchLevel<NDIM> > old_level,
    bool allocate_data)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    d_l_data_manager->setPatchHierarchy(hierarchy);
    d_l_data_manager->resetLevels(0, finest_hier_level);
    d_l_data_manager->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
    if (initial_time && d_l_data_manager->levelContainsLagrangianData(level_number))
    {
        Pointer<LData> Grad_U_data = d_l_data_manager->createLData("Grad_U",level_number,NDIM*NDIM,/*manage_data*/true);
        Pointer<LData>      F_data = d_l_data_manager->createLData("F"     ,level_number,NDIM*NDIM,/*manage_data*/true);
        Pointer<LData>    tau_data = d_l_data_manager->createLData("tau"   ,level_number,NDIM*NDIM,/*manage_data*/true);
        if (d_silo_writer)
        {
            d_silo_writer->registerVariableData("F0", F_data, 0*NDIM, NDIM, level_number);
            d_silo_writer->registerVariableData("F1", F_data, 1*NDIM, NDIM, level_number);
#if (NDIM == 3)
            d_silo_writer->registerVariableData("F2", F_data, 2*NDIM, NDIM, level_number);
#endif
            d_silo_writer->registerVariableData("tau0", tau_data, 0*NDIM, NDIM, level_number);
            d_silo_writer->registerVariableData("tau1", tau_data, 1*NDIM, NDIM, level_number);
#if (NDIM == 3)
            d_silo_writer->registerVariableData("tau2", tau_data, 2*NDIM, NDIM, level_number);
#endif
        }

        // Initialize the deformation gradient and Kirchhoff stress.
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(level_number);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        blitz::Array<double,2>&   F_array = *  F_data->getLocalFormVecArray();
        blitz::Array<double,2>& tau_array = *tau_data->getLocalFormVecArray();
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int idx = node_idx->getLocalPETScIndex();
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    F_array  (idx,NDIM*i+j) = (i == j ? 1.0 : 0.0);
                    tau_array(idx,NDIM*i+j) = 0.0;
                }
            }
        }
    }
    if (!d_load_balancer.isNull() && d_l_data_manager->levelContainsLagrangianData(level_number))
    {
        d_load_balancer->setWorkloadPatchDataIndex(d_workload_idx, level_number);
        d_l_data_manager->updateWorkloadEstimates(level_number, level_number);
    }
    return;
}// initializeLevelData

void
IMPMethod::resetHierarchyConfiguration(
    Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    int coarsest_level,
    int finest_level)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    d_l_data_manager->setPatchHierarchy(hierarchy);
    d_l_data_manager->resetLevels(0, finest_hier_level);
    d_l_data_manager->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    return;
}// resetHierarchyConfiguration

void
IMPMethod::applyGradientDetector(
    Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    int level_number,
    double error_data_time,
    int tag_index,
    bool initial_time,
    bool uses_richardson_extrapolation_too)
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Tag cells that contain Lagrangian nodes.
    d_l_data_manager->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    return;
}// applyGradientDetector

void
IMPMethod::putToDatabase(
    Pointer<Database> db)
{
    db->putInteger("IMP_METHOD_VERSION", IMP_METHOD_VERSION);
    db->putIntegerArray("d_ghosts", d_ghosts, NDIM);
    return;
}// putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IMPMethod::getPositionData(
    std::vector<Pointer<LData> >** X_data,
    bool** X_needs_ghost_fill,
    double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        *X_data = &d_X_current_data;
        *X_needs_ghost_fill = &d_X_current_needs_ghost_fill;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
            if (!d_X_half_data[ln])
            {
                d_X_half_data[ln] = d_l_data_manager->createLData("X_half",ln,NDIM);
                d_X_half_needs_reinit = true;
            }
        }
        if (d_X_half_needs_reinit)
        {
            reinitMidpointData(d_X_current_data, d_X_new_data, d_X_half_data);
            d_X_half_needs_reinit = false;
            d_X_half_needs_ghost_fill = true;
        }
        *X_data = &d_X_half_data;
        *X_needs_ghost_fill = &d_X_half_needs_ghost_fill;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        *X_data = &d_X_new_data;
        *X_needs_ghost_fill = &d_X_new_needs_ghost_fill;
    }
    return;
}// getPositionData

void
IMPMethod::getVelocityData(
    std::vector<Pointer<LData> >** U_data,
    std::vector<Pointer<LData> >** Grad_U_data,
    double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        *U_data      = &     d_U_current_data;
        *Grad_U_data = &d_Grad_U_current_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
            if (!d_U_half_data[ln])
            {
                d_U_half_data     [ln] = d_l_data_manager->createLData("U_half"     ,ln,NDIM     );
                d_Grad_U_half_data[ln] = d_l_data_manager->createLData("Grad_U_half",ln,NDIM*NDIM);
                d_U_half_needs_reinit = true;
            }
        }
        if (d_U_half_needs_reinit)
        {
            reinitMidpointData(     d_U_current_data,      d_U_new_data,      d_U_half_data);
            reinitMidpointData(d_Grad_U_current_data, d_Grad_U_new_data, d_Grad_U_half_data);
            d_U_half_needs_reinit = false;
        }
        *U_data      = &     d_U_half_data;
        *Grad_U_data = &d_Grad_U_half_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        *U_data      = &     d_U_new_data;
        *Grad_U_data = &d_Grad_U_new_data;
    }
    return;
}// getVelocityData

void
IMPMethod::getDeformationGradientData(
    std::vector<Pointer<LData> >** F_data,
    double data_time)
{
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        *F_data = &d_F_current_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        *F_data = &d_F_half_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        *F_data = &d_F_new_data;
    }
    return;
}// getDeformationGradientData

void
IMPMethod::reinitMidpointData(
    const std::vector<Pointer<LData> >& current_data,
    const std::vector<Pointer<LData> >& new_data,
    const std::vector<Pointer<LData> >& half_data)
{
    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        ierr = VecAXPBYPCZ(half_data[ln]->getVec(), 0.5, 0.5, 0.0, current_data[ln]->getVec(), new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
    }
    return;
}// reinitMidpointData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IMPMethod::getFromInput(
    Pointer<Database> db,
    bool is_from_restart)
{
    if (!is_from_restart)
    {
        if (db->isInteger("min_ghost_cell_width"))
        {
            d_ghosts = db->getInteger("min_ghost_cell_width");
        }
        else if (db->isDouble("min_ghost_cell_width"))
        {
            d_ghosts = static_cast<int>(std::ceil(db->getDouble("min_ghost_cell_width")));
        }
    }
    if      (db->keyExists("do_log"        )) d_do_log = db->getBool("do_log"        );
    else if (db->keyExists("enable_logging")) d_do_log = db->getBool("enable_logging");
    return;
}// getFromInput

void
IMPMethod::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to "
                   << d_object_name << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("IMP_METHOD_VERSION");
    if (ver != IMP_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    db->getIntegerArray("d_ghosts", d_ghosts, NDIM);
    return;
}// getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
