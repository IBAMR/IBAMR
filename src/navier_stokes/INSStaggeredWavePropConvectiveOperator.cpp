// Filename: INSStaggeredWavePropConvectiveOperator.cpp
// Created on 12 Jan 2018 by Aaron Barrett
//
// Copyright (c) 2002-2017, Boyce Griffith
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
#include "ibamr/INSStaggeredWavePropConvectiveOperator.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "FaceData.h"
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
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "boost/array.hpp"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

extern "C"
{
#if (NDIM == 2)
    void adv_diff_wp_convective_op2d_(const double*,
                                      const int&,
                                      const double*,
                                      const double*,
                                      const int&,
                                      const double*,
                                      const int&,
                                      const int&,
                                      const int&,
                                      const int&,
                                      const int&,
                                      const int&,
                                      const double*,
                                      const int&);
    void navier_stokes_interp_comps2d_(const int&,
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
                                       double*);
#endif
#if (NDIM == 3)
    void adv_diff_wp_convective_op3d_(const double* q_data,
                                      const int& q_gcw,
                                      const double* u_data_0,
                                      const double* u_data_1,
                                      const double* u_data_2,
                                      const int& u_gcw,
                                      const double* r_data,
                                      const int& r_gcw,
                                      const int& depth,
                                      const int& ilower0,
                                      const int& ilower1,
                                      const int& ilower2,
                                      const int& iupper0,
                                      const int& iupper1,
                                      const int& iupper2,
                                      const double* dx,
                                      const int& k);

    void navier_stokes_interp_comps3d_(const int&,
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
                                       double*);
#endif
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredWavePropConvectiveOperator::INSStaggeredWavePropConvectiveOperator(
    const std::string& object_name,
    Pointer<Database> input_db,
    const ConvectiveDifferencingType difference_form,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
    : ConvectiveOperator(object_name, difference_form),
      d_bc_coefs(bc_coefs),
      d_bdry_extrap_type("CONSTANT"),
      d_hierarchy(NULL),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_U_var(NULL),
      d_U_scratch_idx(-1)
{
    if (d_difference_form != ADVECTIVE /* && d_difference_form != CONSERVATIVE && d_difference_form != SKEW_SYMMETRIC*/)
    {
        TBOX_ERROR("INSStaggeredWavePropConvectiveOperator::INSStaggeredWavePropConvectiveOperator():\n"
                   << "  unsupported differencing form: "
                   << enum_to_string<ConvectiveDifferencingType>(d_difference_form) << " \n"
                   << "  valid choices are: ADVECTIVE\n");
    }

    if (input_db)
    {
        if (input_db->keyExists("bdry_extrap_type")) d_bdry_extrap_type = input_db->getString("bdry_extrap_type");
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("INSStaggeredWavePropConvectiveOperator::CONTEXT");

    d_k = 3;

    const std::string U_var_name = "INSStaggeredWavePropConvectiveOperator::U";
    d_U_var = var_db->getVariable(U_var_name);
    if (d_U_var)
    {
        d_U_scratch_idx = var_db->mapVariableAndContextToIndex(d_U_var, context);
    }
    else
    {
        d_U_var = new SideVariable<NDIM, double>(U_var_name);
        d_U_scratch_idx = var_db->registerVariableAndContext(d_U_var, context, IntVector<NDIM>(d_k + 1));
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(d_U_scratch_idx >= 0);
#endif
    return;
} // INSStaggeredWavePropConvectiveOperator

INSStaggeredWavePropConvectiveOperator::~INSStaggeredWavePropConvectiveOperator()
{
    deallocateOperatorState();
    return;
} // ~INSStaggeredWavePropConvectiveOperator

void
INSStaggeredWavePropConvectiveOperator::applyConvectiveOperator(const int U_idx, const int N_idx)
{
#if !defined(NDEBUG)
    if (!d_is_initialized)
    {
        TBOX_ERROR("INSStaggeredWavePropConvectiveOperator::applyConvectiveOperator():\n"
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
    HierarchyMathOps hier_math_ops("HierarchyMathOps", d_hierarchy);
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
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);

    // Compute the convective derivative.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();
            const IntVector<NDIM>& patch_upper = patch_box.upper();

            Pointer<SideData<NDIM, double> > N_data = patch->getPatchData(N_idx);
            const IntVector<NDIM> N_gcw = N_data->getGhostCellWidth();
            Pointer<SideData<NDIM, double> > U_data = patch->getPatchData(d_U_scratch_idx);
            const IntVector<NDIM> U_gcw = U_data->getGhostCellWidth();

            boost::array<Box<NDIM>, NDIM> side_boxes;
            boost::array<Pointer<FaceData<NDIM, double> >, NDIM> U_adv_data;
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                side_boxes[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                U_adv_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, IntVector<NDIM>(1));
            }
#if (NDIM == 2)
            navier_stokes_interp_comps2d_(patch_lower(0),
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
            navier_stokes_interp_comps3d_(patch_lower(0),
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
            // Do differencing
            for (int axis = 0; axis < NDIM; ++axis)
            {
#if (NDIM == 2)
                adv_diff_wp_convective_op2d_(U_data->getPointer(axis),
                                             U_gcw.max(),
                                             U_adv_data[axis]->getPointer(0),
                                             U_adv_data[axis]->getPointer(1),
                                             U_adv_data[axis]->getGhostCellWidth().max(),
                                             N_data->getPointer(axis),
                                             N_gcw.max(),
                                             1,
                                             side_boxes[axis].lower(0),
                                             side_boxes[axis].lower(1),
                                             side_boxes[axis].upper(0),
                                             side_boxes[axis].upper(1),
                                             dx,
                                             d_k);
#endif
#if (NDIM == 3)
                adv_diff_wp_convective_op3d_(U_data->getPointer(axis),
                                             U_gcw.max(),
                                             U_adv_data[axis]->getPointer(0),
                                             U_adv_data[axis]->getPointer(1),
                                             U_adv_data[axis]->getPointer(2),
                                             U_adv_data[axis]->getGhostCellWidth().max(),
                                             N_data->getPointer(axis),
                                             N_gcw.max(),
                                             1,
                                             side_boxes[axis].lower(0),
                                             side_boxes[axis].lower(1),
                                             side_boxes[axis].lower(2),
                                             side_boxes[axis].upper(0),
                                             side_boxes[axis].upper(1),
                                             side_boxes[axis].upper(2),
                                             dx,
                                             d_k);
#endif
            }
        }
    }

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_U_scratch_idx);
    }

    return;
} // applyConvectiveOperator

void
INSStaggeredWavePropConvectiveOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                                const SAMRAIVectorReal<NDIM, double>& out)
{
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

    return;
} // initializeOperatorState

void
INSStaggeredWavePropConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    // Deallocate the communications operators and BC helpers.
    d_hier_bdry_fill.setNull();
    d_bc_helper.setNull();

    d_is_initialized = false;

    return;
} // deallocateOperatorState
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
