// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
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

#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSStaggeredWavePropConvectiveOperator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "FaceData.h"
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
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <array>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

#if (NDIM == 2)
#define ADV_DIFF_WP_CONVECTIVE_OP_FC IBAMR_FC_FUNC_(adv_diff_wp_convective_op2d, ADV_DIFF_WP_CONVECTIVE_OP2D)
#define NAVIER_STOKES_INTERP_COMPS_FC IBAMR_FC_FUNC_(navier_stokes_interp_comps2d, NAVIER_STOKES_INTERP_COMPS2D)
#endif
#if (NDIM == 3)
#define ADV_DIFF_WP_CONVECTIVE_OP_FC IBAMR_FC_FUNC_(adv_diff_wp_convective_op3d, ADV_DIFF_WP_CONVECTIVE_OP3D)
#define NAVIER_STOKES_INTERP_COMPS_FC IBAMR_FC_FUNC_(navier_stokes_interp_comps3d, NAVIER_STOKES_INTERP_COMPS3D)
#endif

extern "C"
{
#if (NDIM == 2)
    // TODO mangle names
    void ADV_DIFF_WP_CONVECTIVE_OP_FC(const double*,
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
    void NAVIER_STOKES_INTERP_COMPS_FC(const int&,
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
    void ADV_DIFF_WP_CONVECTIVE_OP_FC(const double* q_data,
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

    void NAVIER_STOKES_INTERP_COMPS_FC(const int&,
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
    std::string object_name,
    Pointer<Database> input_db,
    const ConvectiveDifferencingType difference_form,
    std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs)
    : ConvectiveOperator(std::move(object_name), difference_form), d_bc_coefs(std::move(bc_coefs))
{
    if (d_difference_form != ADVECTIVE /* && d_difference_form != CONSERVATIVE && d_difference_form != SKEW_SYMMETRIC*/)
    {
        TBOX_ERROR(
            "INSStaggeredWavePropConvectiveOperator::"
            "INSStaggeredWavePropConvectiveOperator():\n"
            << "  unsupported differencing form: " << enum_to_string<ConvectiveDifferencingType>(d_difference_form)
            << " \n"
            << "  valid choices are: ADVECTIVE\n");
    }

    if (input_db)
    {
        if (input_db->keyExists("bdry_extrap_type")) d_bdry_extrap_type = input_db->getString("bdry_extrap_type");
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("INSStaggeredWavePropConvectiveOperator::CONTEXT");

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
                   << "  operator must be initialized prior to call to "
                      "applyConvectiveOperator\n");
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
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
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
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(d_bc_coefs, nullptr, d_U_scratch_idx, -1, homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_bc_coefs, nullptr);
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

            std::array<Box<NDIM>, NDIM> side_boxes;
            std::array<Pointer<FaceData<NDIM, double> >, NDIM> U_adv_data;
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                side_boxes[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                U_adv_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, IntVector<NDIM>(1));
            }
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
            // Do differencing
            for (int axis = 0; axis < NDIM; ++axis)
            {
#if (NDIM == 2)
                ADV_DIFF_WP_CONVECTIVE_OP_FC(U_data->getPointer(axis),
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
                ADV_DIFF_WP_CONVECTIVE_OP_FC(U_data->getPointer(axis),
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
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
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
