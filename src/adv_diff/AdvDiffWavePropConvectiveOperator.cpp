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

#include "ibamr/AdvDiffPhysicalBoundaryUtilities.h"
#include "ibamr/AdvDiffWavePropConvectiveOperator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/CartExtrapPhysBdryOp.h"

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenSchedule.h"
#include "FaceData.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefinePatchStrategy.h"
#include "RefineSchedule.h"
#include "SAMRAIVectorReal.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/namespaces.h"

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

namespace IBAMR
{
// Constructor
AdvDiffWavePropConvectiveOperator::AdvDiffWavePropConvectiveOperator(
    std::string object_name,
    Pointer<CellVariable<NDIM, double> > Q_var,
    Pointer<Database> input_db,
    const ConvectiveDifferencingType differencing_form,
    std::vector<RobinBcCoefStrategy<NDIM>*> conc_bc_coefs)
    : ConvectiveOperator(std::move(object_name), differencing_form),
      d_Q_var(Q_var),
      d_conc_bc_coefs(std::move(conc_bc_coefs)),
      d_difference_form(differencing_form)
{
    if (d_difference_form != ADVECTIVE /* && d_difference_form != CONSERVATIVE && d_difference_form != SKEW_SYMMETRIC*/)
    {
        TBOX_ERROR(
            "AdvDiffWavePropConvectiveOperator::"
            "AdvDiffWavePropConvectiveOperator():\n"
            << "  unsupported differencing form: " << enum_to_string<ConvectiveDifferencingType>(d_difference_form)
            << " \n"
            << "  valid choices are: ADVECTIVE\n");
    }

    if (input_db)
    {
        if (input_db->keyExists("outflow_bdry_extrap_type"))
            d_outflow_bdry_extrap_type = input_db->getString("outflow_bdry_extrap_type");
        if (input_db->keyExists("bdry_extrap_type"))
        {
            TBOX_ERROR("AdvDiffWavePropConvectiveOperator::AdvDiffWavePropConvectiveOperator():\n"
                       << "  input database key ``bdry_extrap_type'' has been changed to "
                          "``outflow_bdry_extrap_type''\n");
        }
    }

    // Register some scratch variables
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext(d_object_name + "::CONVEC_CONTEXT");
    d_Q_scratch_idx = var_db->registerVariableAndContext(d_Q_var, context, IntVector<NDIM>(d_k + 1));

} // Constructor

AdvDiffWavePropConvectiveOperator::~AdvDiffWavePropConvectiveOperator()
{
    deallocateOperatorState();
    return;
}

void
AdvDiffWavePropConvectiveOperator::applyConvectiveOperator(int Q_idx, int Y_idx)
{
    if (!d_is_initialized)
    {
        TBOX_ERROR("AdvDiffWavePropConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to "
                      "applyConvectiveOperator\n");
    }
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    // Set up refine algorithms for Q and u.
    Pointer<RefineAlgorithm<NDIM> > refine_alg_Q = new RefineAlgorithm<NDIM>();
    Pointer<RefineOperator<NDIM> > refine_op_Q = grid_geom->lookupRefineOperator(d_Q_var, "CONSERVATIVE_LINEAR_REFINE");
    refine_alg_Q->registerRefine(d_Q_scratch_idx, Q_idx, d_Q_scratch_idx, refine_op_Q);
    // Set up coarsen algorithms for Q and u.
    Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg_Q = new CoarsenAlgorithm<NDIM>();
    Pointer<CoarsenOperator<NDIM> > coarsen_op_Q = grid_geom->lookupCoarsenOperator(d_Q_var, "CONSERVATIVE_COARSEN");
    coarsen_alg_Q->registerCoarsen(d_Q_scratch_idx, d_Q_scratch_idx, coarsen_op_Q);
    // Refine the data for Q and u
    d_ghostfill_scheds_Q.resize(d_finest_ln + 1);
    for (int level_num = d_coarsest_ln; level_num <= d_finest_ln; ++level_num)
    {
        refine_alg_Q->resetSchedule(d_ghostfill_scheds_Q[level_num]);
        d_ghostfill_scheds_Q[level_num]->fillData(d_solution_time);
        d_ghostfill_alg_Q->resetSchedule(d_ghostfill_scheds_Q[level_num]);
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > Q_data = patch->getPatchData(d_Q_scratch_idx);
            Pointer<FaceData<NDIM, double> > u_adv_data = patch->getPatchData(d_u_idx);
            AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(Q_data,
                                                                            u_adv_data,
                                                                            patch,
                                                                            d_conc_bc_coefs,
                                                                            d_solution_time,
                                                                            d_outflow_bdry_extrap_type != "NONE",
                                                                            d_homogeneous_bc);
        }
    }
    for (int level_num = d_finest_ln; level_num > d_coarsest_ln; --level_num)
    {
        d_coarsen_scheds_Q[level_num]->coarsenData();
    }

    for (int level_num = d_coarsest_ln; level_num <= d_finest_ln; ++level_num)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
            const double* dx = p_geom->getDx();
            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM> patch_lower = patch_box.lower();
            const IntVector<NDIM> patch_upper = patch_box.upper();
            Pointer<CellData<NDIM, double> > Y_data = patch->getPatchData(Y_idx);
            Pointer<CellData<NDIM, double> > Q_data_scr = patch->getPatchData(d_Q_scratch_idx);
            const IntVector<NDIM> Q_data_scr_gcw = Q_data_scr->getGhostCellWidth();
            Pointer<FaceData<NDIM, double> > U_data = patch->getPatchData(d_u_idx);
            const IntVector<NDIM> U_data_gcw = U_data->getGhostCellWidth();
            const IntVector<NDIM> Y_data_gcw = Y_data->getGhostCellWidth();
#if (NDIM == 2)
            // COMPUTE CONVECTIVE OPERATOR HERE
            ADV_DIFF_WP_CONVECTIVE_OP_FC(Q_data_scr->getPointer(),
                                         Q_data_scr_gcw.max(),
                                         U_data->getPointer(0),
                                         U_data->getPointer(1),
                                         U_data_gcw.max(),
                                         Y_data->getPointer(0),
                                         Y_data_gcw.max(),
                                         Q_data_scr->getDepth(),
                                         patch_lower(0),
                                         patch_lower(1),
                                         patch_upper(0),
                                         patch_upper(1),
                                         dx,
                                         d_k);

#endif
#if (NDIM == 3)
            ADV_DIFF_WP_CONVECTIVE_OP_FC(Q_data_scr->getPointer(),
                                         Q_data_scr_gcw.max(),
                                         U_data->getPointer(0),
                                         U_data->getPointer(1),
                                         U_data->getPointer(2),
                                         U_data_gcw.max(),
                                         Y_data->getPointer(0),
                                         Y_data_gcw.max(),
                                         Q_data_scr->getDepth(),
                                         patch_lower(0),
                                         patch_lower(1),
                                         patch_lower(2),
                                         patch_upper(0),
                                         patch_upper(1),
                                         patch_upper(2),
                                         dx,
                                         d_k);

#endif
        } // end Patch loop
    }     // end Level loop
} // end applyConvectiveOperator

void
AdvDiffWavePropConvectiveOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                           const SAMRAIVectorReal<NDIM, double>& out)
{
    if (d_is_initialized) deallocateOperatorState();
    // Get Hierarchy Information
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
    /* Set up the coasen operations. These COARSEN the data (i.e. fills data at
     * coarse interfaces)
     * General process:
     * 1) Set up a coarsen algorithm
     * 2) Register a coarsen operator with the algorithm
     * 3) Fill a coarsen schedule with the coarsen algorithm
     * 4) To actually coarsen data, use coarsen schedule -> coarsen data()
     */
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<CoarsenOperator<NDIM> > coarsen_op_Q = grid_geom->lookupCoarsenOperator(d_Q_var, "CONSERVATIVE_COARSEN");
    // Step 1) and 2)
    d_coarsen_alg_Q = new CoarsenAlgorithm<NDIM>();
    d_coarsen_alg_Q->registerCoarsen(d_Q_scratch_idx, d_Q_scratch_idx, coarsen_op_Q);
    d_coarsen_scheds_Q.resize(d_finest_ln + 1);
    // Step 3)
    for (int ln = d_coarsest_ln + 1; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM> > coarser_level = d_hierarchy->getPatchLevel(ln - 1);
        d_coarsen_scheds_Q[ln] = d_coarsen_alg_Q->createSchedule(coarser_level, level);
    }
    /* Set Refine Algorithms. This interpolates data onto finer grid
     * General process:
     * 1) Set up a refine algorithm
     * 2) Register a refine operation with the algorithm
     * 3) Fill a refine schedule with the refine algorithm
     * 4) Invoke fill data() inside refine schedule
     */
    // Note we only set up refine algorithms for Q here because u has not been set
    // yet.
    Pointer<RefineOperator<NDIM> > refine_op_Q = grid_geom->lookupRefineOperator(d_Q_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ghostfill_alg_Q = new RefineAlgorithm<NDIM>();
    d_ghostfill_alg_Q->registerRefine(d_Q_scratch_idx, in.getComponentDescriptorIndex(0), d_Q_scratch_idx, refine_op_Q);
    if (d_outflow_bdry_extrap_type != "NONE")
        d_ghostfill_strategy_Q = new CartExtrapPhysBdryOp(d_Q_scratch_idx, d_outflow_bdry_extrap_type);
    d_ghostfill_scheds_Q.resize(d_finest_ln + 1);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_ghostfill_scheds_Q[ln] =
            d_ghostfill_alg_Q->createSchedule(level, ln - 1, d_hierarchy, d_ghostfill_strategy_Q);
    }
    // Allocate Patch Data
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_Q_scratch_idx)) level->allocatePatchData(d_Q_scratch_idx);
    }
    d_is_initialized = true;
    return;
}

void
AdvDiffWavePropConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;
    // Deallocate scratch data
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_Q_scratch_idx))
        {
            level->deallocatePatchData(d_Q_scratch_idx);
        }
    }
    // Deallocate the refine algorithm, operator, patch strategy, and schedules.
    d_ghostfill_alg_Q.setNull();
    d_ghostfill_strategy_Q.setNull();
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        d_ghostfill_scheds_Q[ln].setNull();
    }
    d_ghostfill_scheds_Q.clear();

    d_is_initialized = false;
    return;
}
} // namespace IBAMR
