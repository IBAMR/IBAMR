// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2021 by the IBAMR developers
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
#include "ibamr/INSCollocatedWavePropConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/samrai_compatibility_names.h"

#include "MultiblockDataTranslator.h"
#include "SAMRAIBox.h"
#include "SAMRAICartesianPatchGeometry.h"
#include "SAMRAICellData.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIDatabase.h"
#include "SAMRAIFaceData.h"
#include "SAMRAIIndex.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIRobinBcCoefStrategy.h"
#include "SAMRAISAMRAIVectorReal.h"
#include "SAMRAIUtilities.h"
#include "SAMRAIVariable.h"
#include "SAMRAIVariableContext.h"
#include "SAMRAIVariableDatabase.h"

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
#endif
#if (NDIM == 3)
#define ADV_DIFF_WP_CONVECTIVE_OP_FC IBAMR_FC_FUNC_(adv_diff_wp_convective_op3d, ADV_DIFF_WP_CONVECTIVE_OP3D)
#endif

extern "C"
{
#if (NDIM == 2)
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
#endif
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

INSCollocatedWavePropConvectiveOperator::INSCollocatedWavePropConvectiveOperator(
    std::string object_name,
    SAMRAIPointer<SAMRAIDatabase> input_db,
    const ConvectiveDifferencingType difference_form,
    std::vector<SAMRAIRobinBcCoefStrategy*> bc_coefs)
    : ConvectiveOperator(std::move(object_name), difference_form), d_bc_coefs(std::move(bc_coefs))
{
    if (d_difference_form != ADVECTIVE /* && d_difference_form != CONSERVATIVE && d_difference_form != SKEW_SYMMETRIC*/)
    {
        TBOX_ERROR(
            "INSCollocatedWavePropConvectiveOperator::"
            "INSCollocatedWavePropConvectiveOperator():\n"
            << "  unsupported differencing form: " << enum_to_string<ConvectiveDifferencingType>(d_difference_form)
            << " \n"
            << "  valid choices are: ADVECTIVE\n");
    }

    if (input_db)
    {
        if (input_db->keyExists("bdry_extrap_type")) d_bdry_extrap_type = input_db->getString("bdry_extrap_type");
    }

    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    SAMRAIPointer<SAMRAIVariableContext> context =
        var_db->getContext("INSCollocatedWavePropConvectiveOperator::CONTEXT");

    const std::string U_var_name = "INSCollocatedWavePropConvectiveOperator::U";
    d_U_var = var_db->getVariable(U_var_name);
    if (d_U_var)
    {
        d_U_scratch_idx = var_db->mapVariableAndContextToIndex(d_U_var, context);
    }
    else
    {
        d_U_var = new SAMRAICellVariable<double>(U_var_name, NDIM);
        d_U_scratch_idx = var_db->registerVariableAndContext(d_U_var, context, SAMRAIIntVector(d_k + 1));
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(d_U_scratch_idx >= 0);
#endif
    return;
} // INSCollocatedWavePropConvectiveOperator

INSCollocatedWavePropConvectiveOperator::~INSCollocatedWavePropConvectiveOperator()
{
    deallocateOperatorState();
    return;
} // ~INSCollocatedWavePropConvectiveOperator

void
INSCollocatedWavePropConvectiveOperator::applyConvectiveOperator(const int U_idx, const int N_idx)
{
#if !defined(NDEBUG)
    if (!d_is_initialized)
    {
        TBOX_ERROR("INSCollocatedWavePropConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to "
                      "applyConvectiveOperator\n");
    }
#endif

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAIPointer<SAMRAIPatchLevel> level = d_hierarchy->getPatchLevel(ln);
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
                                                             true,
                                                             "CONSERVATIVE_COARSEN",
                                                             d_bdry_extrap_type,
                                                             false,
                                                             d_bc_coefs);
    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time);
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);

    // Compute the convective derivative.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAIPointer<SAMRAIPatchLevel> level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAIPatchLevel::Iterator p(level); p; p++)
        {
            SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());

            const SAMRAIPointer<SAMRAICartesianPatchGeometry> patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            const SAMRAIBox& patch_box = patch->getBox();
            const SAMRAIIntVector& patch_lower = patch_box.lower();
            const SAMRAIIntVector& patch_upper = patch_box.upper();

            SAMRAIPointer<SAMRAICellData<double>> N_data = patch->getPatchData(N_idx);
            const SAMRAIIntVector N_gcw = N_data->getGhostCellWidth();
            SAMRAIPointer<SAMRAICellData<double>> U_data = patch->getPatchData(d_U_scratch_idx);
            const SAMRAIIntVector U_gcw = U_data->getGhostCellWidth();
            SAMRAIPointer<SAMRAIFaceData<double>> U_sp_data = patch->getPatchData(d_u_idx);
            const SAMRAIIntVector U_sp_gcw = U_sp_data->getGhostCellWidth();

            // Do differencing
            for (int d = 0; d < NDIM; ++d)
            {
#if (NDIM == 2)
                ADV_DIFF_WP_CONVECTIVE_OP_FC(U_data->getPointer(d),
                                             U_gcw.max(),
                                             U_sp_data->getPointer(0),
                                             U_sp_data->getPointer(1),
                                             U_sp_gcw.max(),
                                             N_data->getPointer(d),
                                             N_gcw.max(),
                                             1,
                                             patch_lower(0),
                                             patch_lower(1),
                                             patch_upper(0),
                                             patch_upper(1),
                                             dx,
                                             d_k);
#endif
#if (NDIM == 3)
                ADV_DIFF_WP_CONVECTIVE_OP_FC(U_data->getPointer(d),
                                             U_gcw.max(),
                                             U_sp_data->getPointer(0),
                                             U_sp_data->getPointer(1),
                                             U_sp_data->getPointer(2),
                                             U_sp_gcw.max(),
                                             N_data->getPointer(d),
                                             N_gcw.max(),
                                             1,
                                             patch_lower(0),
                                             patch_lower(1),
                                             patch_lower(2),
                                             patch_upper(0),
                                             patch_upper(1),
                                             patch_upper(2),
                                             dx,
                                             d_k);
#endif
            }
        }
    }
    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAIPointer<SAMRAIPatchLevel> level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_U_scratch_idx);
    }

    return;
} // applyConvectiveOperator

void
INSCollocatedWavePropConvectiveOperator::initializeOperatorState(const SAMRAISAMRAIVectorReal<double>& in,
                                                                 const SAMRAISAMRAIVectorReal<double>& out)
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

    d_is_initialized = true;

    return;
} // initializeOperatorState

void
INSCollocatedWavePropConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    // Deallocate the communications operators and BC helpers.
    d_hier_bdry_fill.setNull();

    d_is_initialized = false;

    return;
} // deallocateOperatorState
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
