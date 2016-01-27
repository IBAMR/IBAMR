// Filename: INSStaggeredStabilizedPPMConvectiveOperator.cpp
// Created on 07 Sep 2012 by Boyce Griffith
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
#include <cmath>
#include <ostream>
#include <string>
#include <vector>

#include "ArrayData.h"
#include "Box.h"
#include "BoxArray.h"
#include "CartesianPatchGeometry.h"
#include "FaceData.h"
#include "FaceIndex.h"
#include "FaceIterator.h"
#include "GridGeometry.h"
#include "IBAMR_config.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "boost/array.hpp"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSStaggeredStabilizedPPMConvectiveOperator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/ibtk_utilities.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

// FORTRAN ROUTINES
#if (NDIM == 2)
#define ADVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(advect_derivative2d, ADVECT_DERIVATIVE2D)
#define CONVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(convect_derivative2d, CONVECT_DERIVATIVE2D)
#define GODUNOV_EXTRAPOLATE_FC IBAMR_FC_FUNC_(godunov_extrapolate2d, GODUNOV_EXTRAPOLATE2D)
#define NAVIER_STOKES_INTERP_COMPS_FC IBAMR_FC_FUNC_(navier_stokes_interp_comps2d, NAVIER_STOKES_INTERP_COMPS2D)
#define NAVIER_STOKES_RESET_ADV_VELOCITY_FC                                                                            \
    IBAMR_FC_FUNC_(navier_stokes_reset_adv_velocity2d, NAVIER_STOKES_RESET_ADV_VELOCITY2D)
#define SKEW_SYM_DERIVATIVE_FC IBAMR_FC_FUNC_(skew_sym_derivative2d, SKEW_SYM_DERIVATIVE2D)
#endif

#if (NDIM == 3)
#define ADVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(advect_derivative3d, ADVECT_DERIVATIVE3D)
#define CONVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(convect_derivative3d, CONVECT_DERIVATIVE3D)
#define GODUNOV_EXTRAPOLATE_FC IBAMR_FC_FUNC_(godunov_extrapolate3d, GODUNOV_EXTRAPOLATE3D)
#define NAVIER_STOKES_INTERP_COMPS_FC IBAMR_FC_FUNC_(navier_stokes_interp_comps3d, NAVIER_STOKES_INTERP_COMPS3D)
#define NAVIER_STOKES_RESET_ADV_VELOCITY_FC                                                                            \
    IBAMR_FC_FUNC_(navier_stokes_reset_adv_velocity3d, NAVIER_STOKES_RESET_ADV_VELOCITY3D)
#define SKEW_SYM_DERIVATIVE_FC IBAMR_FC_FUNC_(skew_sym_derivative3d, SKEW_SYM_DERIVATIVE3D)
#endif

extern "C" {
void ADVECT_DERIVATIVE_FC(const double*,
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
                          const double*,
                          const double*,
                          const int&,
                          const int&,
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
                          const double*,
                          const double*,
                          const double*,
                          const int&,
                          const int&,
                          const int&,
#endif
                          double*);

void CONVECT_DERIVATIVE_FC(const double*,
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
                           const double*,
                           const double*,
                           const int&,
                           const int&,
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
                           const double*,
                           const double*,
                           const double*,
                           const int&,
                           const int&,
                           const int&,
#endif
                           double*);

void GODUNOV_EXTRAPOLATE_FC(
#if (NDIM == 2)
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    double*,
    double*,
    double*,
    double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
    double*,
    double*
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
    const double*,
    double*,
    double*,
    double*,
    double*,
    double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
    const double*,
    double*,
    double*,
    double*
#endif
    );

void NAVIER_STOKES_INTERP_COMPS_FC(
#if (NDIM == 2)
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    double*,
    double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    double*,
    double*
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
    const double*,
    const double*,
    const double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    double*,
    double*,
    double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    double*,
    double*,
    double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    double*,
    double*,
    double*
#endif
    );

void NAVIER_STOKES_RESET_ADV_VELOCITY_FC(
#if (NDIM == 2)
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    double*,
    double*,
    const int&,
    const int&,
    const double*,
    const double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    double*,
    double*,
    const int&,
    const int&,
    const double*,
    const double*
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
    double*,
    double*,
    double*,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
    const double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    double*,
    double*,
    double*,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
    const double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    double*,
    double*,
    double*,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
    const double*
#endif
    );

void SKEW_SYM_DERIVATIVE_FC(const double*,
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
                            const double*,
                            const double*,
                            const int&,
                            const int&,
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
                            const double*,
                            const double*,
                            const double*,
                            const int&,
                            const int&,
                            const int&,
#endif
                            double*);
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

inline double
smooth_kernel(const double r)
{
    return std::abs(r) < 1.0 ? 0.5 * (cos(M_PI * r) + 1.0) : 0.0;
} // smooth_kernel

// Timers.
static Timer* t_apply_convective_operator;
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredStabilizedPPMConvectiveOperator::INSStaggeredStabilizedPPMConvectiveOperator(
    const std::string& object_name,
    Pointer<Database> input_db,
    const ConvectiveDifferencingType difference_form,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
    : ConvectiveOperator(object_name, difference_form),
      d_open_bdry(array_constant<bool, 2 * NDIM>(false)),
      d_width(array_constant<double, 2 * NDIM>(0.0)),
      d_bc_coefs(bc_coefs),
      d_bdry_extrap_type("CONSTANT"),
      d_hierarchy(NULL),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_U_var(NULL),
      d_U_scratch_idx(-1)
{
    if (d_difference_form != ADVECTIVE && d_difference_form != CONSERVATIVE && d_difference_form != SKEW_SYMMETRIC)
    {
        TBOX_ERROR(
            "INSStaggeredStabilizedPPMConvectiveOperator::"
            "INSStaggeredStabilizedPPMConvectiveOperator():\n"
            << "  unsupported differencing form: "
            << enum_to_string<ConvectiveDifferencingType>(d_difference_form)
            << " \n"
            << "  valid choices are: ADVECTIVE, CONSERVATIVE, SKEW_SYMMETRIC\n");
    }

    if (input_db)
    {
        if (input_db->keyExists("bdry_extrap_type")) d_bdry_extrap_type = input_db->getString("bdry_extrap_type");
        for (unsigned int location_index = 0; location_index < 2 * NDIM; ++location_index)
        {
            std::ostringstream width_stream;
            width_stream << "width_" << location_index;
            const std::string width_key = width_stream.str();
            if (input_db->keyExists(width_key))
            {
                d_width[location_index] = input_db->getDouble(width_key);
            }
            if (d_width[location_index] > 0.0) d_open_bdry[location_index] = true;
        }
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("INSStaggeredStabilizedPPMConvectiveOperator::CONTEXT");

    const std::string U_var_name = "INSStaggeredStabilizedPPMConvectiveOperator::U";
    d_U_var = var_db->getVariable(U_var_name);
    if (d_U_var)
    {
        d_U_scratch_idx = var_db->mapVariableAndContextToIndex(d_U_var, context);
    }
    else
    {
        d_U_var = new SideVariable<NDIM, double>(U_var_name);
        d_U_scratch_idx = var_db->registerVariableAndContext(d_U_var, context, IntVector<NDIM>(GADVECTG));
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(d_U_scratch_idx >= 0);
#endif

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_apply_convective_operator = TimerManager::getManager()->getTimer(
            "IBAMR::INSStaggeredStabilizedPPMConvectiveOperator::applyConvectiveOperator()");
        t_apply = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredStabilizedPPMConvectiveOperator::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer(
            "IBAMR::INSStaggeredStabilizedPPMConvectiveOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer(
            "IBAMR::INSStaggeredStabilizedPPMConvectiveOperator::deallocateOperatorState()"););
    return;
} // INSStaggeredStabilizedPPMConvectiveOperator

INSStaggeredStabilizedPPMConvectiveOperator::~INSStaggeredStabilizedPPMConvectiveOperator()
{
    deallocateOperatorState();
    return;
} // ~INSStaggeredStabilizedPPMConvectiveOperator

void
INSStaggeredStabilizedPPMConvectiveOperator::applyConvectiveOperator(const int U_idx, const int N_idx)
{
    IBAMR_TIMER_START(t_apply_convective_operator);
#if !defined(NDEBUG)
    if (!d_is_initialized)
    {
        TBOX_ERROR("INSStaggeredStabilizedPPMConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to applyConvectiveOperator\n");
    }
    TBOX_ASSERT(U_idx == d_u_idx);
#endif

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_U_scratch_idx);
    }

    // Fill ghost cell values for all components.
    static const bool homogeneous_bc = false;
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(1);
    transaction_comps[0] = InterpolationTransactionComponent(d_U_scratch_idx,
                                                             U_idx,
                                                             "CONSERVATIVE_LINEAR_REFINE",
                                                             false,
                                                             "CONSERVATIVE_COARSEN",
                                                             d_bdry_extrap_type,
                                                             false,
                                                             d_bc_coefs);
    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(homogeneous_bc);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(d_bc_coefs, NULL, d_U_scratch_idx, -1, homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_bc_coefs, NULL);
    //  d_bc_helper->enforceDivergenceFreeConditionAtBoundary(d_U_scratch_idx);
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);

    // Compute the convective derivative.
    Pointer<GridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const IntVector<NDIM>& ratio = level->getRatio();
        const Box<NDIM> domain_box = Box<NDIM>::refine(grid_geometry->getPhysicalDomain()[0], ratio);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();
            const double* const x_lower = patch_geom->getXLower();
            const double* const x_upper = patch_geom->getXUpper();

            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();
            const IntVector<NDIM>& patch_upper = patch_box.upper();

            Pointer<SideData<NDIM, double> > N_data = patch->getPatchData(N_idx);
            Pointer<SideData<NDIM, double> > N_upwind_data =
                new SideData<NDIM, double>(N_data->getBox(), N_data->getDepth(), N_data->getGhostCellWidth());
            Pointer<SideData<NDIM, double> > U_data = patch->getPatchData(d_U_scratch_idx);

            const IntVector<NDIM> ghosts = IntVector<NDIM>(1);
            boost::array<Box<NDIM>, NDIM> side_boxes;
            boost::array<Pointer<FaceData<NDIM, double> >, NDIM> U_adv_data;
            boost::array<Pointer<FaceData<NDIM, double> >, NDIM> U_half_data;
            boost::array<Pointer<FaceData<NDIM, double> >, NDIM> U_half_upwind_data;
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                side_boxes[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                U_adv_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                U_half_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                U_half_upwind_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
            }

// Interpolate the staggered-grid velocity field onto the faces of
// the control volumes.
#if (NDIM == 2)
            NAVIER_STOKES_INTERP_COMPS_FC(patch_lower(0),
                                          patch_upper(0),
                                          patch_lower(1),
                                          patch_upper(1),
                                          U_data->getGhostCellWidth()(0),
                                          U_data->getGhostCellWidth()(1),
                                          U_data->getPointer(0),
                                          U_data->getPointer(1),
                                          side_boxes[0].lower(0),
                                          side_boxes[0].upper(0),
                                          side_boxes[0].lower(1),
                                          side_boxes[0].upper(1),
                                          U_adv_data[0]->getGhostCellWidth()(0),
                                          U_adv_data[0]->getGhostCellWidth()(1),
                                          U_adv_data[0]->getPointer(0),
                                          U_adv_data[0]->getPointer(1),
                                          side_boxes[1].lower(0),
                                          side_boxes[1].upper(0),
                                          side_boxes[1].lower(1),
                                          side_boxes[1].upper(1),
                                          U_adv_data[1]->getGhostCellWidth()(0),
                                          U_adv_data[1]->getGhostCellWidth()(1),
                                          U_adv_data[1]->getPointer(0),
                                          U_adv_data[1]->getPointer(1));
#endif
#if (NDIM == 3)
            NAVIER_STOKES_INTERP_COMPS_FC(patch_lower(0),
                                          patch_upper(0),
                                          patch_lower(1),
                                          patch_upper(1),
                                          patch_lower(2),
                                          patch_upper(2),
                                          U_data->getGhostCellWidth()(0),
                                          U_data->getGhostCellWidth()(1),
                                          U_data->getGhostCellWidth()(2),
                                          U_data->getPointer(0),
                                          U_data->getPointer(1),
                                          U_data->getPointer(2),
                                          side_boxes[0].lower(0),
                                          side_boxes[0].upper(0),
                                          side_boxes[0].lower(1),
                                          side_boxes[0].upper(1),
                                          side_boxes[0].lower(2),
                                          side_boxes[0].upper(2),
                                          U_adv_data[0]->getGhostCellWidth()(0),
                                          U_adv_data[0]->getGhostCellWidth()(1),
                                          U_adv_data[0]->getGhostCellWidth()(2),
                                          U_adv_data[0]->getPointer(0),
                                          U_adv_data[0]->getPointer(1),
                                          U_adv_data[0]->getPointer(2),
                                          side_boxes[1].lower(0),
                                          side_boxes[1].upper(0),
                                          side_boxes[1].lower(1),
                                          side_boxes[1].upper(1),
                                          side_boxes[1].lower(2),
                                          side_boxes[1].upper(2),
                                          U_adv_data[1]->getGhostCellWidth()(0),
                                          U_adv_data[1]->getGhostCellWidth()(1),
                                          U_adv_data[1]->getGhostCellWidth()(2),
                                          U_adv_data[1]->getPointer(0),
                                          U_adv_data[1]->getPointer(1),
                                          U_adv_data[1]->getPointer(2),
                                          side_boxes[2].lower(0),
                                          side_boxes[2].upper(0),
                                          side_boxes[2].lower(1),
                                          side_boxes[2].upper(1),
                                          side_boxes[2].lower(2),
                                          side_boxes[2].upper(2),
                                          U_adv_data[2]->getGhostCellWidth()(0),
                                          U_adv_data[2]->getGhostCellWidth()(1),
                                          U_adv_data[2]->getGhostCellWidth()(2),
                                          U_adv_data[2]->getPointer(0),
                                          U_adv_data[2]->getPointer(1),
                                          U_adv_data[2]->getPointer(2));
#endif

            // Compute the first-order upwind discretization.
            if (patch_geom->getTouchesRegularBoundary())
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    const ArrayData<NDIM, double>& U_array_data = U_data->getArrayData(axis);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        for (FaceIterator<NDIM> ic(side_boxes[axis], d); ic; ic++)
                        {
                            const FaceIndex<NDIM>& i = ic();
                            const double u_ADV = (*U_adv_data[axis])(i);
                            const double U_lower = U_array_data(i.toCell(0), 0);
                            const double U_upper = U_array_data(i.toCell(1), 0);
                            (*U_half_upwind_data[axis])(i) =
                                (u_ADV > 1.0e-8) ? U_lower : (u_ADV < 1.0e-8) ? U_upper : 0.5 * (U_lower + U_upper);
                        }
                    }
                }
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    switch (d_difference_form)
                    {
                    case CONSERVATIVE:
#if (NDIM == 2)
                        CONVECT_DERIVATIVE_FC(dx,
                                              side_boxes[axis].lower(0),
                                              side_boxes[axis].upper(0),
                                              side_boxes[axis].lower(1),
                                              side_boxes[axis].upper(1),
                                              U_adv_data[axis]->getGhostCellWidth()(0),
                                              U_adv_data[axis]->getGhostCellWidth()(1),
                                              U_half_upwind_data[axis]->getGhostCellWidth()(0),
                                              U_half_upwind_data[axis]->getGhostCellWidth()(1),
                                              U_adv_data[axis]->getPointer(0),
                                              U_adv_data[axis]->getPointer(1),
                                              U_half_upwind_data[axis]->getPointer(0),
                                              U_half_upwind_data[axis]->getPointer(1),
                                              N_upwind_data->getGhostCellWidth()(0),
                                              N_upwind_data->getGhostCellWidth()(1),
                                              N_upwind_data->getPointer(axis));
#endif
#if (NDIM == 3)
                        CONVECT_DERIVATIVE_FC(dx,
                                              side_boxes[axis].lower(0),
                                              side_boxes[axis].upper(0),
                                              side_boxes[axis].lower(1),
                                              side_boxes[axis].upper(1),
                                              side_boxes[axis].lower(2),
                                              side_boxes[axis].upper(2),
                                              U_adv_data[axis]->getGhostCellWidth()(0),
                                              U_adv_data[axis]->getGhostCellWidth()(1),
                                              U_adv_data[axis]->getGhostCellWidth()(2),
                                              U_half_upwind_data[axis]->getGhostCellWidth()(0),
                                              U_half_upwind_data[axis]->getGhostCellWidth()(1),
                                              U_half_upwind_data[axis]->getGhostCellWidth()(2),
                                              U_adv_data[axis]->getPointer(0),
                                              U_adv_data[axis]->getPointer(1),
                                              U_adv_data[axis]->getPointer(2),
                                              U_half_upwind_data[axis]->getPointer(0),
                                              U_half_upwind_data[axis]->getPointer(1),
                                              U_half_upwind_data[axis]->getPointer(2),
                                              N_upwind_data->getGhostCellWidth()(0),
                                              N_upwind_data->getGhostCellWidth()(1),
                                              N_upwind_data->getGhostCellWidth()(2),
                                              N_upwind_data->getPointer(axis));
#endif
                        break;
                    case ADVECTIVE:
#if (NDIM == 2)
                        ADVECT_DERIVATIVE_FC(dx,
                                             side_boxes[axis].lower(0),
                                             side_boxes[axis].upper(0),
                                             side_boxes[axis].lower(1),
                                             side_boxes[axis].upper(1),
                                             U_adv_data[axis]->getGhostCellWidth()(0),
                                             U_adv_data[axis]->getGhostCellWidth()(1),
                                             U_half_upwind_data[axis]->getGhostCellWidth()(0),
                                             U_half_upwind_data[axis]->getGhostCellWidth()(1),
                                             U_adv_data[axis]->getPointer(0),
                                             U_adv_data[axis]->getPointer(1),
                                             U_half_upwind_data[axis]->getPointer(0),
                                             U_half_upwind_data[axis]->getPointer(1),
                                             N_upwind_data->getGhostCellWidth()(0),
                                             N_upwind_data->getGhostCellWidth()(1),
                                             N_upwind_data->getPointer(axis));
#endif
#if (NDIM == 3)
                        ADVECT_DERIVATIVE_FC(dx,
                                             side_boxes[axis].lower(0),
                                             side_boxes[axis].upper(0),
                                             side_boxes[axis].lower(1),
                                             side_boxes[axis].upper(1),
                                             side_boxes[axis].lower(2),
                                             side_boxes[axis].upper(2),
                                             U_adv_data[axis]->getGhostCellWidth()(0),
                                             U_adv_data[axis]->getGhostCellWidth()(1),
                                             U_adv_data[axis]->getGhostCellWidth()(2),
                                             U_half_upwind_data[axis]->getGhostCellWidth()(0),
                                             U_half_upwind_data[axis]->getGhostCellWidth()(1),
                                             U_half_upwind_data[axis]->getGhostCellWidth()(2),
                                             U_adv_data[axis]->getPointer(0),
                                             U_adv_data[axis]->getPointer(1),
                                             U_adv_data[axis]->getPointer(2),
                                             U_half_upwind_data[axis]->getPointer(0),
                                             U_half_upwind_data[axis]->getPointer(1),
                                             U_half_upwind_data[axis]->getPointer(2),
                                             N_upwind_data->getGhostCellWidth()(0),
                                             N_upwind_data->getGhostCellWidth()(1),
                                             N_upwind_data->getGhostCellWidth()(2),
                                             N_upwind_data->getPointer(axis));
#endif
                        break;
                    case SKEW_SYMMETRIC:
#if (NDIM == 2)
                        SKEW_SYM_DERIVATIVE_FC(dx,
                                               side_boxes[axis].lower(0),
                                               side_boxes[axis].upper(0),
                                               side_boxes[axis].lower(1),
                                               side_boxes[axis].upper(1),
                                               U_adv_data[axis]->getGhostCellWidth()(0),
                                               U_adv_data[axis]->getGhostCellWidth()(1),
                                               U_half_upwind_data[axis]->getGhostCellWidth()(0),
                                               U_half_upwind_data[axis]->getGhostCellWidth()(1),
                                               U_adv_data[axis]->getPointer(0),
                                               U_adv_data[axis]->getPointer(1),
                                               U_half_upwind_data[axis]->getPointer(0),
                                               U_half_upwind_data[axis]->getPointer(1),
                                               N_upwind_data->getGhostCellWidth()(0),
                                               N_upwind_data->getGhostCellWidth()(1),
                                               N_upwind_data->getPointer(axis));
#endif
#if (NDIM == 3)
                        SKEW_SYM_DERIVATIVE_FC(dx,
                                               side_boxes[axis].lower(0),
                                               side_boxes[axis].upper(0),
                                               side_boxes[axis].lower(1),
                                               side_boxes[axis].upper(1),
                                               side_boxes[axis].lower(2),
                                               side_boxes[axis].upper(2),
                                               U_adv_data[axis]->getGhostCellWidth()(0),
                                               U_adv_data[axis]->getGhostCellWidth()(1),
                                               U_adv_data[axis]->getGhostCellWidth()(2),
                                               U_half_upwind_data[axis]->getGhostCellWidth()(0),
                                               U_half_upwind_data[axis]->getGhostCellWidth()(1),
                                               U_half_upwind_data[axis]->getGhostCellWidth()(2),
                                               U_adv_data[axis]->getPointer(0),
                                               U_adv_data[axis]->getPointer(1),
                                               U_adv_data[axis]->getPointer(2),
                                               U_half_upwind_data[axis]->getPointer(0),
                                               U_half_upwind_data[axis]->getPointer(1),
                                               U_half_upwind_data[axis]->getPointer(2),
                                               N_upwind_data->getGhostCellWidth()(0),
                                               N_upwind_data->getGhostCellWidth()(1),
                                               N_upwind_data->getGhostCellWidth()(2),
                                               N_upwind_data->getPointer(axis));
#endif
                        break;
                    default:
                        TBOX_ERROR(
                            "INSStaggeredStabilizedConvectiveOperator::"
                            "applyConvectiveOperator():\n"
                            << "  unsupported differencing form: "
                            << enum_to_string<ConvectiveDifferencingType>(d_difference_form)
                            << " \n"
                            << "  valid choices are: ADVECTIVE, CONSERVATIVE, "
                               "SKEW_SYMMETRIC\n");
                    }
                }
            }

            // Compute the xsPPM7 discretization.
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                Pointer<SideData<NDIM, double> > dU_data =
                    new SideData<NDIM, double>(U_data->getBox(), U_data->getDepth(), U_data->getGhostCellWidth());
                Pointer<SideData<NDIM, double> > U_L_data =
                    new SideData<NDIM, double>(U_data->getBox(), U_data->getDepth(), U_data->getGhostCellWidth());
                Pointer<SideData<NDIM, double> > U_R_data =
                    new SideData<NDIM, double>(U_data->getBox(), U_data->getDepth(), U_data->getGhostCellWidth());
                Pointer<SideData<NDIM, double> > U_scratch1_data =
                    new SideData<NDIM, double>(U_data->getBox(), U_data->getDepth(), U_data->getGhostCellWidth());
#if (NDIM == 3)
                Pointer<SideData<NDIM, double> > U_scratch2_data =
                    new SideData<NDIM, double>(U_data->getBox(), U_data->getDepth(), U_data->getGhostCellWidth());
#endif
#if (NDIM == 2)
                GODUNOV_EXTRAPOLATE_FC(side_boxes[axis].lower(0),
                                       side_boxes[axis].upper(0),
                                       side_boxes[axis].lower(1),
                                       side_boxes[axis].upper(1),
                                       U_data->getGhostCellWidth()(0),
                                       U_data->getGhostCellWidth()(1),
                                       U_data->getPointer(axis),
                                       U_scratch1_data->getPointer(axis),
                                       dU_data->getPointer(axis),
                                       U_L_data->getPointer(axis),
                                       U_R_data->getPointer(axis),
                                       U_adv_data[axis]->getGhostCellWidth()(0),
                                       U_adv_data[axis]->getGhostCellWidth()(1),
                                       U_half_data[axis]->getGhostCellWidth()(0),
                                       U_half_data[axis]->getGhostCellWidth()(1),
                                       U_adv_data[axis]->getPointer(0),
                                       U_adv_data[axis]->getPointer(1),
                                       U_half_data[axis]->getPointer(0),
                                       U_half_data[axis]->getPointer(1));
#endif
#if (NDIM == 3)
                GODUNOV_EXTRAPOLATE_FC(side_boxes[axis].lower(0),
                                       side_boxes[axis].upper(0),
                                       side_boxes[axis].lower(1),
                                       side_boxes[axis].upper(1),
                                       side_boxes[axis].lower(2),
                                       side_boxes[axis].upper(2),
                                       U_data->getGhostCellWidth()(0),
                                       U_data->getGhostCellWidth()(1),
                                       U_data->getGhostCellWidth()(2),
                                       U_data->getPointer(axis),
                                       U_scratch1_data->getPointer(axis),
                                       U_scratch2_data->getPointer(axis),
                                       dU_data->getPointer(axis),
                                       U_L_data->getPointer(axis),
                                       U_R_data->getPointer(axis),
                                       U_adv_data[axis]->getGhostCellWidth()(0),
                                       U_adv_data[axis]->getGhostCellWidth()(1),
                                       U_adv_data[axis]->getGhostCellWidth()(2),
                                       U_half_data[axis]->getGhostCellWidth()(0),
                                       U_half_data[axis]->getGhostCellWidth()(1),
                                       U_half_data[axis]->getGhostCellWidth()(2),
                                       U_adv_data[axis]->getPointer(0),
                                       U_adv_data[axis]->getPointer(1),
                                       U_adv_data[axis]->getPointer(2),
                                       U_half_data[axis]->getPointer(0),
                                       U_half_data[axis]->getPointer(1),
                                       U_half_data[axis]->getPointer(2));
#endif
            }
#if (NDIM == 2)
            NAVIER_STOKES_RESET_ADV_VELOCITY_FC(side_boxes[0].lower(0),
                                                side_boxes[0].upper(0),
                                                side_boxes[0].lower(1),
                                                side_boxes[0].upper(1),
                                                U_adv_data[0]->getGhostCellWidth()(0),
                                                U_adv_data[0]->getGhostCellWidth()(1),
                                                U_adv_data[0]->getPointer(0),
                                                U_adv_data[0]->getPointer(1),
                                                U_half_data[0]->getGhostCellWidth()(0),
                                                U_half_data[0]->getGhostCellWidth()(1),
                                                U_half_data[0]->getPointer(0),
                                                U_half_data[0]->getPointer(1),
                                                side_boxes[1].lower(0),
                                                side_boxes[1].upper(0),
                                                side_boxes[1].lower(1),
                                                side_boxes[1].upper(1),
                                                U_adv_data[1]->getGhostCellWidth()(0),
                                                U_adv_data[1]->getGhostCellWidth()(1),
                                                U_adv_data[1]->getPointer(0),
                                                U_adv_data[1]->getPointer(1),
                                                U_half_data[1]->getGhostCellWidth()(0),
                                                U_half_data[1]->getGhostCellWidth()(1),
                                                U_half_data[1]->getPointer(0),
                                                U_half_data[1]->getPointer(1));
#endif
#if (NDIM == 3)
            NAVIER_STOKES_RESET_ADV_VELOCITY_FC(side_boxes[0].lower(0),
                                                side_boxes[0].upper(0),
                                                side_boxes[0].lower(1),
                                                side_boxes[0].upper(1),
                                                side_boxes[0].lower(2),
                                                side_boxes[0].upper(2),
                                                U_adv_data[0]->getGhostCellWidth()(0),
                                                U_adv_data[0]->getGhostCellWidth()(1),
                                                U_adv_data[0]->getGhostCellWidth()(2),
                                                U_adv_data[0]->getPointer(0),
                                                U_adv_data[0]->getPointer(1),
                                                U_adv_data[0]->getPointer(2),
                                                U_half_data[0]->getGhostCellWidth()(0),
                                                U_half_data[0]->getGhostCellWidth()(1),
                                                U_half_data[0]->getGhostCellWidth()(2),
                                                U_half_data[0]->getPointer(0),
                                                U_half_data[0]->getPointer(1),
                                                U_half_data[0]->getPointer(2),
                                                side_boxes[1].lower(0),
                                                side_boxes[1].upper(0),
                                                side_boxes[1].lower(1),
                                                side_boxes[1].upper(1),
                                                side_boxes[1].lower(2),
                                                side_boxes[1].upper(2),
                                                U_adv_data[1]->getGhostCellWidth()(0),
                                                U_adv_data[1]->getGhostCellWidth()(1),
                                                U_adv_data[1]->getGhostCellWidth()(2),
                                                U_adv_data[1]->getPointer(0),
                                                U_adv_data[1]->getPointer(1),
                                                U_adv_data[1]->getPointer(2),
                                                U_half_data[1]->getGhostCellWidth()(0),
                                                U_half_data[1]->getGhostCellWidth()(1),
                                                U_half_data[1]->getGhostCellWidth()(2),
                                                U_half_data[1]->getPointer(0),
                                                U_half_data[1]->getPointer(1),
                                                U_half_data[1]->getPointer(2),
                                                side_boxes[2].lower(0),
                                                side_boxes[2].upper(0),
                                                side_boxes[2].lower(1),
                                                side_boxes[2].upper(1),
                                                side_boxes[2].lower(2),
                                                side_boxes[2].upper(2),
                                                U_adv_data[2]->getGhostCellWidth()(0),
                                                U_adv_data[2]->getGhostCellWidth()(1),
                                                U_adv_data[2]->getGhostCellWidth()(2),
                                                U_adv_data[2]->getPointer(0),
                                                U_adv_data[2]->getPointer(1),
                                                U_adv_data[2]->getPointer(2),
                                                U_half_data[2]->getGhostCellWidth()(0),
                                                U_half_data[2]->getGhostCellWidth()(1),
                                                U_half_data[2]->getGhostCellWidth()(2),
                                                U_half_data[2]->getPointer(0),
                                                U_half_data[2]->getPointer(1),
                                                U_half_data[2]->getPointer(2));
#endif
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                switch (d_difference_form)
                {
                case CONSERVATIVE:
#if (NDIM == 2)
                    CONVECT_DERIVATIVE_FC(dx,
                                          side_boxes[axis].lower(0),
                                          side_boxes[axis].upper(0),
                                          side_boxes[axis].lower(1),
                                          side_boxes[axis].upper(1),
                                          U_adv_data[axis]->getGhostCellWidth()(0),
                                          U_adv_data[axis]->getGhostCellWidth()(1),
                                          U_half_data[axis]->getGhostCellWidth()(0),
                                          U_half_data[axis]->getGhostCellWidth()(1),
                                          U_adv_data[axis]->getPointer(0),
                                          U_adv_data[axis]->getPointer(1),
                                          U_half_data[axis]->getPointer(0),
                                          U_half_data[axis]->getPointer(1),
                                          N_data->getGhostCellWidth()(0),
                                          N_data->getGhostCellWidth()(1),
                                          N_data->getPointer(axis));
#endif
#if (NDIM == 3)
                    CONVECT_DERIVATIVE_FC(dx,
                                          side_boxes[axis].lower(0),
                                          side_boxes[axis].upper(0),
                                          side_boxes[axis].lower(1),
                                          side_boxes[axis].upper(1),
                                          side_boxes[axis].lower(2),
                                          side_boxes[axis].upper(2),
                                          U_adv_data[axis]->getGhostCellWidth()(0),
                                          U_adv_data[axis]->getGhostCellWidth()(1),
                                          U_adv_data[axis]->getGhostCellWidth()(2),
                                          U_half_data[axis]->getGhostCellWidth()(0),
                                          U_half_data[axis]->getGhostCellWidth()(1),
                                          U_half_data[axis]->getGhostCellWidth()(2),
                                          U_adv_data[axis]->getPointer(0),
                                          U_adv_data[axis]->getPointer(1),
                                          U_adv_data[axis]->getPointer(2),
                                          U_half_data[axis]->getPointer(0),
                                          U_half_data[axis]->getPointer(1),
                                          U_half_data[axis]->getPointer(2),
                                          N_data->getGhostCellWidth()(0),
                                          N_data->getGhostCellWidth()(1),
                                          N_data->getGhostCellWidth()(2),
                                          N_data->getPointer(axis));
#endif
                    break;
                case ADVECTIVE:
#if (NDIM == 2)
                    ADVECT_DERIVATIVE_FC(dx,
                                         side_boxes[axis].lower(0),
                                         side_boxes[axis].upper(0),
                                         side_boxes[axis].lower(1),
                                         side_boxes[axis].upper(1),
                                         U_adv_data[axis]->getGhostCellWidth()(0),
                                         U_adv_data[axis]->getGhostCellWidth()(1),
                                         U_half_data[axis]->getGhostCellWidth()(0),
                                         U_half_data[axis]->getGhostCellWidth()(1),
                                         U_adv_data[axis]->getPointer(0),
                                         U_adv_data[axis]->getPointer(1),
                                         U_half_data[axis]->getPointer(0),
                                         U_half_data[axis]->getPointer(1),
                                         N_data->getGhostCellWidth()(0),
                                         N_data->getGhostCellWidth()(1),
                                         N_data->getPointer(axis));
#endif
#if (NDIM == 3)
                    ADVECT_DERIVATIVE_FC(dx,
                                         side_boxes[axis].lower(0),
                                         side_boxes[axis].upper(0),
                                         side_boxes[axis].lower(1),
                                         side_boxes[axis].upper(1),
                                         side_boxes[axis].lower(2),
                                         side_boxes[axis].upper(2),
                                         U_adv_data[axis]->getGhostCellWidth()(0),
                                         U_adv_data[axis]->getGhostCellWidth()(1),
                                         U_adv_data[axis]->getGhostCellWidth()(2),
                                         U_half_data[axis]->getGhostCellWidth()(0),
                                         U_half_data[axis]->getGhostCellWidth()(1),
                                         U_half_data[axis]->getGhostCellWidth()(2),
                                         U_adv_data[axis]->getPointer(0),
                                         U_adv_data[axis]->getPointer(1),
                                         U_adv_data[axis]->getPointer(2),
                                         U_half_data[axis]->getPointer(0),
                                         U_half_data[axis]->getPointer(1),
                                         U_half_data[axis]->getPointer(2),
                                         N_data->getGhostCellWidth()(0),
                                         N_data->getGhostCellWidth()(1),
                                         N_data->getGhostCellWidth()(2),
                                         N_data->getPointer(axis));
#endif
                    break;
                case SKEW_SYMMETRIC:
#if (NDIM == 2)
                    SKEW_SYM_DERIVATIVE_FC(dx,
                                           side_boxes[axis].lower(0),
                                           side_boxes[axis].upper(0),
                                           side_boxes[axis].lower(1),
                                           side_boxes[axis].upper(1),
                                           U_adv_data[axis]->getGhostCellWidth()(0),
                                           U_adv_data[axis]->getGhostCellWidth()(1),
                                           U_half_data[axis]->getGhostCellWidth()(0),
                                           U_half_data[axis]->getGhostCellWidth()(1),
                                           U_adv_data[axis]->getPointer(0),
                                           U_adv_data[axis]->getPointer(1),
                                           U_half_data[axis]->getPointer(0),
                                           U_half_data[axis]->getPointer(1),
                                           N_data->getGhostCellWidth()(0),
                                           N_data->getGhostCellWidth()(1),
                                           N_data->getPointer(axis));
#endif
#if (NDIM == 3)
                    SKEW_SYM_DERIVATIVE_FC(dx,
                                           side_boxes[axis].lower(0),
                                           side_boxes[axis].upper(0),
                                           side_boxes[axis].lower(1),
                                           side_boxes[axis].upper(1),
                                           side_boxes[axis].lower(2),
                                           side_boxes[axis].upper(2),
                                           U_adv_data[axis]->getGhostCellWidth()(0),
                                           U_adv_data[axis]->getGhostCellWidth()(1),
                                           U_adv_data[axis]->getGhostCellWidth()(2),
                                           U_half_data[axis]->getGhostCellWidth()(0),
                                           U_half_data[axis]->getGhostCellWidth()(1),
                                           U_half_data[axis]->getGhostCellWidth()(2),
                                           U_adv_data[axis]->getPointer(0),
                                           U_adv_data[axis]->getPointer(1),
                                           U_adv_data[axis]->getPointer(2),
                                           U_half_data[axis]->getPointer(0),
                                           U_half_data[axis]->getPointer(1),
                                           U_half_data[axis]->getPointer(2),
                                           N_data->getGhostCellWidth()(0),
                                           N_data->getGhostCellWidth()(1),
                                           N_data->getGhostCellWidth()(2),
                                           N_data->getPointer(axis));
#endif
                    break;
                default:
                    TBOX_ERROR(
                        "INSStaggeredStabilizedPPMConvectiveOperator::applyConvectiveOperator("
                        "):\n"
                        << "  unsupported differencing form: "
                        << enum_to_string<ConvectiveDifferencingType>(d_difference_form)
                        << " \n"
                        << "  valid choices are: ADVECTIVE, CONSERVATIVE, SKEW_SYMMETRIC\n");
                }
            }

            // Blend together the low-order and high-order discretizations at
            // physical boundaries.
            if (patch_geom->getTouchesRegularBoundary())
            {
                Pointer<SideData<NDIM, double> > N_PPM_data =
                    new SideData<NDIM, double>(N_data->getBox(), N_data->getDepth(), N_data->getGhostCellWidth());
                N_PPM_data->copy(*N_data);
                for (unsigned int location_index = 0; location_index < 2 * NDIM; ++location_index)
                {
                    const unsigned int axis = location_index / 2;
                    const unsigned int side = location_index % 2;
                    const bool is_lower = side == 0;
                    if (d_open_bdry[location_index] && patch_geom->getTouchesRegularBoundary(axis, side))
                    {
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            Box<NDIM> bdry_box = domain_box;
                            const double width = d_width[location_index];
                            const int offset = static_cast<int>(width / dx[axis]);
                            if (is_lower)
                            {
                                bdry_box.upper(axis) = domain_box.lower(axis) + offset;
                            }
                            else
                            {
                                bdry_box.lower(axis) = domain_box.upper(axis) - offset;
                            }
                            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box * patch_box, d)); b; b++)
                            {
                                const Index<NDIM>& i = b();
                                const SideIndex<NDIM> i_s(i, d, SideIndex<NDIM>::Lower);
                                const double x =
                                    x_lower[axis] + dx[axis] * static_cast<double>(i(axis) - patch_box.lower(axis));
                                const double x_bdry = (is_lower ? x_lower[axis] : x_upper[axis]);
                                const double fac = smooth_kernel((x - x_bdry) / width);
                                (*N_data)(i_s) = fac * (*N_upwind_data)(i_s) + (1.0 - fac) * (*N_PPM_data)(i_s);
                            }
                        }
                    }
                }
            }
        }
    }

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_U_scratch_idx);
    }

    IBAMR_TIMER_STOP(t_apply_convective_operator);
    return;
} // applyConvectiveOperator

void
INSStaggeredStabilizedPPMConvectiveOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                                     const SAMRAIVectorReal<NDIM, double>& out)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    if (d_is_initialized) deallocateOperatorState();

    // Get the hierarchy configuration.
    d_hierarchy = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln = in.getFinestLevelNumber();
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
#else
    NULL_USE(out);
#endif

    // Setup the interpolation transaction information.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    d_transaction_comps.resize(1);
    d_transaction_comps[0] = InterpolationTransactionComponent(d_U_scratch_idx,
                                                               in.getComponentDescriptorIndex(0),
                                                               "CONSERVATIVE_LINEAR_REFINE",
                                                               false,
                                                               "CONSERVATIVE_COARSEN",
                                                               d_bdry_extrap_type,
                                                               false,
                                                               d_bc_coefs);

    // Initialize the interpolation operators.
    d_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(d_transaction_comps, d_hierarchy);

    // Initialize the BC helper.
    d_bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
    d_bc_helper->cacheBcCoefData(d_bc_coefs, d_solution_time, d_hierarchy);

    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
INSStaggeredStabilizedPPMConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

    // Deallocate the communications operators and BC helpers.
    d_hier_bdry_fill.setNull();
    d_bc_helper.setNull();

    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
