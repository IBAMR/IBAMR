// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2019 by the IBAMR developers
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

#include "ibamr/AdvDiffBDSQConvectiveOperator.h"
#include "ibamr/namespaces.h"

#include "ibtk/DebuggingUtilities.h"

#include <utility>

extern "C"
{
#if (NDIM == 2)
    void adv_diff_bdsq_convective_op2d_(const double*,
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
                                        const double*);
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
                                      const int&);
#endif
}

namespace IBAMR
{
// Constructor
AdvDiffBDSQConvectiveOperator::AdvDiffBDSQConvectiveOperator(std::string object_name,
                                                             Pointer<CellVariable<NDIM, double> > Q_var,
                                                             Pointer<Database> /*input_db*/,
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
            "AdvDiffBDSQConvectiveOperator::"
            "AdvDiffBDSQConvectiveOperator():\n"
            << "  unsupported differencing form: " << enum_to_string<ConvectiveDifferencingType>(d_difference_form)
            << " \n"
            << "  valid choices are: ADVECTIVE\n");
    }

    // Register some scratch variables
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext(d_object_name + "::CONVEC_CONTEXT");
    d_Q_scratch_idx = var_db->registerVariableAndContext(d_Q_var, context, IntVector<NDIM>(4));

} // Constructor

AdvDiffBDSQConvectiveOperator::~AdvDiffBDSQConvectiveOperator()
{
    deallocateOperatorState();
    return;
}

void
AdvDiffBDSQConvectiveOperator::applyConvectiveOperator(int Q_idx, int Y_idx)
{
    if (!d_is_initialized)
    {
        TBOX_ERROR("AdvDiffBDSQConvectiveOperator::applyConvectiveOperator():\n"
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
            Box<NDIM> grown_box(patch_box);
            grown_box.grow(IntVector<NDIM>(1, 1));
            Pointer<CellData<NDIM, double> > Y_data = patch->getPatchData(Y_idx);
            Pointer<CellData<NDIM, double> > Q_data = patch->getPatchData(d_Q_scratch_idx);
            Pointer<FaceData<NDIM, double> > U_data = patch->getPatchData(d_u_idx);
            FaceData<NDIM, double> F_data(patch_box, 2, IntVector<NDIM>(1, 1));
            NodeData<NDIM, double> N_data(patch_box, 1, IntVector<NDIM>(1, 1));
            for (int d = 0; d < Q_data->getDepth(); ++d)
            {
                // Compute Node data
                for (NodeIterator<NDIM> ni(grown_box); ni; ni++)
                {
                    const NodeIndex<NDIM>& n_idx = ni();
                    CellIndex<NDIM> c_idx(n_idx);
                    N_data(n_idx) =
                        ((*Q_data)(c_idx + IntVector<NDIM>(-2, -2)) -
                         7.0 *
                             ((*Q_data)(c_idx + IntVector<NDIM>(-1, -2)) + (*Q_data)(c_idx + IntVector<NDIM>(0, -2))) +
                         (*Q_data)(c_idx + IntVector<NDIM>(1, -2)) - 7.0 * (*Q_data)(c_idx + IntVector<NDIM>(-2, -1)) +
                         49.0 *
                             ((*Q_data)(c_idx + IntVector<NDIM>(-1, -1)) + (*Q_data)(c_idx + IntVector<NDIM>(0, -1))) -
                         7.0 * (*Q_data)(c_idx + IntVector<NDIM>(1, -1)) -
                         7.0 * (*Q_data)(c_idx + IntVector<NDIM>(-2, 0)) +
                         49.0 * ((*Q_data)(c_idx + IntVector<NDIM>(-1, 0)) + (*Q_data)(c_idx)) -
                         7.0 * (*Q_data)(c_idx + IntVector<NDIM>(1, 0)) + (*Q_data)(c_idx + IntVector<NDIM>(-2, 1)) -
                         7.0 * ((*Q_data)(c_idx + IntVector<NDIM>(-1, 1)) + (*Q_data)(c_idx + IntVector<NDIM>(0, 1))) +
                         (*Q_data)(c_idx + IntVector<NDIM>(1, 1))) /
                        144.0;
                }

                // Reconstruct face data on patch
                std::array<double, 6> p;
                double Qij;
                for (CellIterator<NDIM> ci(grown_box); ci; ci++)
                {
                    bool done_limiting = false;
                    const CellIndex<NDIM>& idx = ci();
                    const FaceIndex<NDIM> fi_xl(idx, 0, 0), fi_xr(idx, 0, 1), fi_yl(idx, 1, 0), fi_yr(idx, 1, 1);
                    double RH = N_data(NodeIndex<NDIM>(idx, IntVector<NDIM>(1, 1))),
                           RL = N_data(NodeIndex<NDIM>(idx, IntVector<NDIM>(1, 0))),
                           LL = N_data(NodeIndex<NDIM>(idx, IntVector<NDIM>(0, 0))),
                           LH = N_data(NodeIndex<NDIM>(idx, IntVector<NDIM>(0, 1)));
                    Qij = (*Q_data)(idx);
                    p[1] = ((RH + RL) - (LH + LL)) / (2.0 * dx[0]);
                    p[2] = ((LH + RH) - (LL + RL)) / (2.0 * dx[1]);
                    p[3] = ((RH - RL) - (LH - LL)) / (dx[0] * dx[1]);
                    p[4] = (-(*Q_data)(idx - IntVector<NDIM>(2, 0)) + 12.0 * (*Q_data)(idx - IntVector<NDIM>(1, 0)) -
                            22.0 * (*Q_data)(idx) + 12.0 * (*Q_data)(idx + IntVector<NDIM>(1, 0)) -
                            (*Q_data)(idx + IntVector<NDIM>(2, 0))) /
                           (16.0 * dx[0] * dx[0]);
                    p[5] = (-(*Q_data)(idx - IntVector<NDIM>(0, 2)) + 12.0 * (*Q_data)(idx - IntVector<NDIM>(0, 1)) -
                            22.0 * (*Q_data)(idx) + 12.0 * (*Q_data)(idx + IntVector<NDIM>(0, 1)) -
                            (*Q_data)(idx + IntVector<NDIM>(0, 2))) /
                           (16.0 * dx[1] * dx[1]);
                    p[0] = (*Q_data)(idx) - (p[4] * dx[0] * dx[0] + p[5] * dx[1] * dx[1]) / 12.0;

                    // Step 1 of limiting
                    if ((Qij > RH && Qij > RL && Qij > LH && Qij > LL) ||
                        (Qij < RH && Qij < RL && Qij < LH && Qij < LL))
                    {
                        p[0] = Qij;
                        p[1] = 0.0;
                        p[2] = 0.0;
                        p[3] = 0.0;
                        p[4] = 0.0;
                        p[5] = 0.0;
                        done_limiting = true;
                    }
                    // Step 2 of limiting
                    if (!done_limiting && !limitingStepTwo(p, dx, Q_data, idx))
                    {
                        // Step 3 of limiting
                        // Reset coefficients
                        p[1] = ((RH + RL) - (LH + LL)) / (2.0 * dx[0]);
                        p[2] = ((LH + RH) - (LL + RL)) / (2.0 * dx[1]);
                        p[3] = ((RH - RL) - (LH - LL)) / (dx[0] * dx[1]);
                        p[4] =
                            (-(*Q_data)(idx - IntVector<NDIM>(2, 0)) + 12.0 * (*Q_data)(idx - IntVector<NDIM>(1, 0)) -
                             22.0 * (*Q_data)(idx) + 12.0 * (*Q_data)(idx + IntVector<NDIM>(1, 0)) -
                             (*Q_data)(idx + IntVector<NDIM>(2, 0))) /
                            (16.0 * dx[0] * dx[0]);
                        p[5] =
                            (-(*Q_data)(idx - IntVector<NDIM>(0, 2)) + 12.0 * (*Q_data)(idx - IntVector<NDIM>(0, 1)) -
                             22.0 * (*Q_data)(idx) + 12.0 * (*Q_data)(idx + IntVector<NDIM>(0, 1)) -
                             (*Q_data)(idx + IntVector<NDIM>(0, 2))) /
                            (16.0 * dx[1] * dx[1]);
                        p[0] = (*Q_data)(idx);
                        limitingStepThree(p, dx, Q_data, idx);
                    }

                    // Evaulate polynomial at cell sides
                    F_data(fi_xl, 1) = evaluateP(p, -0.5 * dx[0], 0.0);
                    F_data(fi_xr, 0) = evaluateP(p, 0.5 * dx[0], 0.0);
                    F_data(fi_yl, 1) = evaluateP(p, 0.0, -0.5 * dx[1]);
                    F_data(fi_yr, 0) = evaluateP(p, 0.0, 0.5 * dx[1]);
                }

                // Do differencing on patch
                for (CellIterator<NDIM> ci(patch_box); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    const FaceIndex<NDIM> fi_xl(idx, 0, 0), fi_xr(idx, 0, 1), fi_yl(idx, 1, 0), fi_yr(idx, 1, 1);
                    (*Y_data)(idx, d) =
                        1.0 / dx[0] *
                        (std::max((*U_data)(fi_xl), 0.0) * (F_data(fi_xl, 1) - F_data(fi_xl, 0)) +
                         std::min((*U_data)(fi_xr), 0.0) * (F_data(fi_xr, 1) - F_data(fi_xr, 0)) +
                         0.5 * ((*U_data)(fi_xr) + (*U_data)(fi_xl)) * (F_data(fi_xr, 0) - F_data(fi_xl, 1)));
                    (*Y_data)(idx, d) +=
                        1.0 / dx[1] *
                        (std::max((*U_data)(fi_yl), 0.0) * (F_data(fi_yl, 1) - F_data(fi_yl, 0)) +
                         std::min((*U_data)(fi_yr), 0.0) * (F_data(fi_yr, 1) - F_data(fi_yr, 0)) +
                         0.5 * ((*U_data)(fi_yr) + (*U_data)(fi_yl)) * (F_data(fi_yr, 0) - F_data(fi_yl, 1)));
                }
            }

        } // end Patch loop
    }     // end Level loop
} // end applyConvectiveOperator

void
AdvDiffBDSQConvectiveOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
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
AdvDiffBDSQConvectiveOperator::deallocateOperatorState()
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

bool
AdvDiffBDSQConvectiveOperator::limitingStepTwo(std::array<double, 6>& p,
                                               const double* const dx,
                                               Pointer<CellData<NDIM, double> > Q_data,
                                               const CellIndex<NDIM>& idx)
{
    bool t1xx = false, t2xx = false, t1yy = false, t2yy = false;
    double cmpx = std::min(std::abs(p[1] + p[3] * dx[1] * 0.5), std::abs(p[1] - p[3] * dx[1] * 0.5)),
           cmpy = std::min(std::abs(p[2] + p[3] * dx[0] * 0.5), std::abs(p[2] - p[3] * dx[0] * 0.5));

    if ((p[1] + p[3] * dx[1] * 0.5) * (p[1] - p[3] * dx[1] * 0.5) < 0.0)
        t1xx = true;
    else if (cmpx < (dx[0] * std::abs(p[4])))
        t2xx = true;
    if ((p[2] + p[3] * dx[0] * 0.5) * (p[2] - p[3] * dx[0] * 0.5) < 0.0)
        t1yy = true;
    else if (cmpy < (dx[1] * std::abs(p[4])))
        t2yy = true;

    if (t1xx && (t1yy || t2yy))
        p[4] = 0.0;
    else if (t2xx && (t1yy || t2yy))
        p[4] = std::copysign(1.0, p[4]) * cmpx / dx[0];

    if (t1yy && (t1xx || t2xx))
        p[5] = 0.0;
    else if (t2yy && (t1xx || t2xx))
        p[5] = std::copysign(1.0, p[5]) * cmpy / dx[1];

    p[0] = (*Q_data)(idx) - (p[4] * dx[0] * dx[0] + p[5] * dx[1] * dx[1]) / 12.0;

    if (testValues(p, dx, Q_data, idx)) return true;

    cmpx = std::min(std::abs(p[1] + p[3] * dx[1] * 0.5), std::abs(p[1] - p[3] * dx[1] * 0.5)),
    cmpy = std::min(std::abs(p[2] + p[3] * dx[0] * 0.5), std::abs(p[2] - p[3] * dx[0] * 0.5));
    if ((p[1] + p[3] * dx[1] * 0.5) * (p[1] - p[3] * dx[1] * 0.5) < 0)
        p[4] = 0.0;
    else if (cmpx < std::abs(dx[0] * p[4]))
        p[4] = std::copysign(1.0, p[4]) * cmpx / dx[0];
    if ((p[2] + p[3] * dx[0] * 0.5) * (p[2] - p[3] * dx[0] * 0.5) < 0)
        p[5] = 0.0;
    else if (cmpy < std::abs(dx[1] * p[5]))
        p[5] = std::copysign(1.0, p[5]) * cmpy / dx[1];

    p[0] = (*Q_data)(idx) - (p[4] * dx[0] * dx[0] + p[5] * dx[1] * dx[1]) / 12.0;

    if (testValues(p, dx, Q_data, idx)) return true;

    return false;
}

bool
AdvDiffBDSQConvectiveOperator::testValues(std::array<double, 6>& p,
                                          const double* const dx,
                                          Pointer<CellData<NDIM, double> > Q_data,
                                          const CellIndex<NDIM>& idx)
{
    // Form array of Q vals
    std::array<std::array<double, 3>, 3> Q;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            Q[i][j] = (*Q_data)(idx + IntVector<NDIM>(1 - i, 1 - j));
        }
    }
    // Lower left
    double val = evaluateP(p, -0.5 * dx[0], -0.5 * dx[1]);
    if ((val < min(Q[0][0], Q[0][1], Q[1][1], Q[0][1])) || (val > max(Q[0][0], Q[0][1], Q[1][1], Q[0][1])))
        return false;
    // Lower right
    val = evaluateP(p, 0.5 * dx[0], -0.5 * dx[1]);
    if ((val < min(Q[2][0], Q[2][1], Q[1][1], Q[1][0])) || (val > max(Q[2][0], Q[2][1], Q[1][1], Q[1][0])))
        return false;
    // Upper left
    val = evaluateP(p, -0.5 * dx[0], 0.5 * dx[1]);
    if ((val < min(Q[0][2], Q[0][1], Q[1][1], Q[1][2])) || (val > max(Q[0][2], Q[0][1], Q[1][1], Q[1][2])))
        return false;
    // Top right
    val = evaluateP(p, 0.5 * dx[0], 0.5 * dx[1]);
    if ((val < min(Q[2][2], Q[2][1], Q[1][1], Q[0][2])) || (val > max(Q[2][2], Q[2][1], Q[1][1], Q[0][2])))
        return false;

    // Top edge
    if (std::abs(p[1] + p[3] * dx[1] * 0.5) < std::abs(p[4] * dx[0]))
    {
        double pt = -(p[1] + p[3] * dx[1] * 0.5) / (2.0 * p[4]);
        val = evaluateP(p, pt, 0.5 * dx[1]);
        if (pt > 0)
        {
            if ((val > max(Q[2][2], Q[2][1], Q[1][1], Q[1][2])) || val < min(Q[2][2], Q[2][1], Q[1][1], Q[1][2]))
                return false;
        }
        else
        {
            if ((val > max(Q[0][2], Q[0][1], Q[1][1], Q[1][2])) || val < min(Q[0][2], Q[0][1], Q[1][1], Q[1][2]))
                return false;
        }
    }
    // Bottom edge
    if (std::abs(p[1] - p[3] * dx[1] * 0.5) < std::abs(p[4] * dx[0]))
    {
        double pt = -(p[1] - p[3] * dx[1] * 0.5) / (2.0 * p[4]);
        val = evaluateP(p, pt, -0.5 * dx[1]);
        if (pt > 0)
        {
            if ((val > max(Q[2][0], Q[2][1], Q[1][1], Q[1][0])) || val < min(Q[2][0], Q[2][1], Q[1][1], Q[1][0]))
                return false;
        }
        else
        {
            if ((val > max(Q[2][0], Q[2][1], Q[1][1], Q[1][0])) || val < min(Q[2][0], Q[2][1], Q[1][1], Q[1][0]))
                return false;
        }
    }
    // Right edge
    if (std::abs(p[2] + p[3] * dx[0] * 0.5) < std::abs(p[5] * dx[1]))
    {
        double pt = -(p[2] + p[3] * dx[0] * 0.5) / (2.0 * p[5]);
        val = evaluateP(p, 0.5 * dx[0], pt);
        if (pt > 0)
        {
            if ((val > max(Q[2][2], Q[2][1], Q[1][1], Q[1][2])) || val < min(Q[2][2], Q[2][1], Q[1][1], Q[1][2]))
                return false;
        }
        else
        {
            if ((val > max(Q[2][0], Q[2][1], Q[1][1], Q[1][0])) || val < min(Q[2][0], Q[2][1], Q[1][1], Q[1][0]))
                return false;
        }
    }
    // Left edge
    if (std::abs(p[2] - p[3] * dx[0] * 0.5) < std::abs(p[5] * dx[1]))
    {
        double pt = -(p[2] - p[3] * dx[0] * 0.5) / (2.0 * p[5]);
        val = evaluateP(p, -0.5 * dx[0], pt);
        if (pt > 0)
        {
            if ((val > max(Q[0][2], Q[0][1], Q[1][1], Q[1][2])) || val < min(Q[0][2], Q[0][1], Q[1][1], Q[1][2]))
                return false;
        }
        else
        {
            if ((val > max(Q[0][0], Q[0][1], Q[1][1], Q[1][0])) || val < min(Q[0][0], Q[0][1], Q[1][1], Q[1][0]))
                return false;
        }
    }
    return true;
}

bool
AdvDiffBDSQConvectiveOperator::limitingStepThree(std::array<double, 6>& p,
                                                 const double* const dx,
                                                 Pointer<CellData<NDIM, double> > Q_data,
                                                 const CellIndex<NDIM>& idx)
{
    double Qij = (*Q_data)(idx);
    // Limit linear and bilinear coefficients
    double LL = p[0] - 0.5 * dx[0] * p[1] - 0.5 * dx[1] * p[2] + 0.25 * dx[0] * dx[1] * p[3],
           LH = p[0] - 0.5 * dx[0] * p[1] + 0.5 * dx[1] * p[2] - 0.25 * dx[0] * dx[1] * p[3],
           RL = p[0] + 0.5 * dx[0] * p[1] - 0.5 * dx[1] * p[2] - 0.25 * dx[0] * dx[1] * p[3],
           RH = p[0] + 0.5 * dx[0] * p[1] + 0.5 * dx[1] * p[2] + 0.25 * dx[0] * dx[1] * p[3];

    std::array<double, 4> minN, maxN, Nvals = { LL, LH, RL, RH };
    minN[0] = min((*Q_data)(idx + IntVector<NDIM>(-1, -1)),
                  (*Q_data)(idx + IntVector<NDIM>(-1, 0)),
                  (*Q_data)(idx + IntVector<NDIM>(0, -1)),
                  (*Q_data)(idx));
    maxN[0] = max((*Q_data)(idx + IntVector<NDIM>(-1, -1)),
                  (*Q_data)(idx + IntVector<NDIM>(-1, 0)),
                  (*Q_data)(idx + IntVector<NDIM>(0, -1)),
                  (*Q_data)(idx));
    minN[1] = min((*Q_data)(idx + IntVector<NDIM>(-1, 1)),
                  (*Q_data)(idx + IntVector<NDIM>(-1, 0)),
                  (*Q_data)(idx + IntVector<NDIM>(0, 1)),
                  (*Q_data)(idx));
    maxN[1] = max((*Q_data)(idx + IntVector<NDIM>(-1, 1)),
                  (*Q_data)(idx + IntVector<NDIM>(-1, 0)),
                  (*Q_data)(idx + IntVector<NDIM>(0, 1)),
                  (*Q_data)(idx));
    minN[2] = min((*Q_data)(idx + IntVector<NDIM>(1, -1)),
                  (*Q_data)(idx + IntVector<NDIM>(1, 0)),
                  (*Q_data)(idx + IntVector<NDIM>(0, -1)),
                  (*Q_data)(idx));
    maxN[2] = max((*Q_data)(idx + IntVector<NDIM>(1, -1)),
                  (*Q_data)(idx + IntVector<NDIM>(1, 0)),
                  (*Q_data)(idx + IntVector<NDIM>(0, -1)),
                  (*Q_data)(idx));
    minN[3] = min((*Q_data)(idx + IntVector<NDIM>(1, 1)),
                  (*Q_data)(idx + IntVector<NDIM>(1, 0)),
                  (*Q_data)(idx + IntVector<NDIM>(0, 1)),
                  (*Q_data)(idx));
    maxN[3] = max((*Q_data)(idx + IntVector<NDIM>(1, 1)),
                  (*Q_data)(idx + IntVector<NDIM>(1, 0)),
                  (*Q_data)(idx + IntVector<NDIM>(0, 1)),
                  (*Q_data)(idx));

    bool done = true;
    for (int i = 0; i < 4; ++i)
    {
        if (Nvals[i] > maxN[i] || Nvals[i] < minN[i])
        {
            done = false;
            Nvals[i] = std::max(std::min(Nvals[i], maxN[i]), minN[i]);
        }
    }

    int iter = 0;
    while (!done && iter < 3)
    {
        double sumdiff = 0.0;
        for (const auto& N : Nvals) sumdiff += N;
        sumdiff -= 4.0 * (*Q_data)(idx);
        if (sumdiff > 0)
        {
            int kdp = 0;
            for (const auto& N : Nvals)
                if (N > Qij + 1.0e-10) kdp++;
            if (kdp == 0)
            {
                done = true;
                continue;
            }
            for (int i = 0; i < 4; ++i)
            {
                if (Nvals[i] > Qij + 1.0e-10)
                {
                    double redfac = std::min(sumdiff / kdp, Nvals[i] - minN[i]);
                    kdp--;
                    sumdiff -= redfac;
                    Nvals[i] -= redfac;
                }
            }
        }
        else
        {
            int kdp = 0;
            for (const auto& N : Nvals)
                if (N < Qij - 1.0e-10) kdp++;
            if (kdp == 0)
            {
                done = true;
                continue;
            }
            for (int i = 0; i < 4; ++i)
            {
                if (Nvals[i] < Qij - 1.0e-10)
                {
                    double redfac = std::max(sumdiff / kdp, Nvals[i] - maxN[i]);
                    kdp--;
                    sumdiff -= redfac;
                    Nvals[i] -= redfac;
                }
            }
        }
        iter++;
    }
    plog << "iter: " << iter << "\n";
    // Recompute slopes
    p[1] = ((Nvals[3] + Nvals[2]) - (Nvals[1] + Nvals[0])) / (2.0 * dx[0]);
    p[2] = ((Nvals[1] + Nvals[3]) - (Nvals[0] + Nvals[2])) / (2.0 * dx[1]);
    p[3] = ((Nvals[3] - Nvals[2]) - (Nvals[1] - Nvals[0])) / (dx[0] * dx[1]);

    double cmpx = std::min(std::abs(p[1] + p[3] * dx[1] * 0.5), std::abs(p[1] - p[3] * dx[1] * 0.5));
    double cmpy = std::min(std::abs(p[2] + p[3] * dx[0] * 0.5), std::abs(p[2] - p[3] * dx[0] * 0.5));
    if ((p[1] + p[3] * dx[1] * 0.5) * (p[1] - p[3] * dx[1] * 0.5) < 0.0)
        p[4] = 0.0;
    else if (cmpx < (dx[0] * std::abs(p[4])))
        p[4] = std::copysign(1.0, p[4]) * cmpx / dx[0];

    if ((p[2] + p[3] * dx[0] * 0.5) * (p[2] - p[3] * dx[0] * 0.5) < 0.0)
        p[5] = 0.0;
    else if (cmpy < (dx[1] * std::abs(p[5])))
        p[5] = std::copysign(1.0, p[5]) * cmpx / dx[1];

    p[0] = (*Q_data)(idx) - (p[4] * dx[0] * dx[0] + p[5] * dx[1] * dx[1]) / 12.0;

    if (!testValues(p, dx, Q_data, idx))
    {
        p[4] = p[5] = 0.0;
        p[0] = (*Q_data)(idx);
    }
    return true;
}

double
AdvDiffBDSQConvectiveOperator::evaluateP(const std::array<double, 6>& p, const double x, const double y)
{
    return p[0] + x * p[1] + y * p[2] + x * y * p[3] + x * x * p[4] + y * y * p[5];
}

double
AdvDiffBDSQConvectiveOperator::max(const double a, const double b, const double c, const double d)
{
    return std::max(a, std::max(b, std::max(c, d)));
}

double
AdvDiffBDSQConvectiveOperator::min(const double a, const double b, const double c, const double d)
{
    return std::min(a, std::min(b, std::min(c, d)));
}

} // namespace IBAMR
