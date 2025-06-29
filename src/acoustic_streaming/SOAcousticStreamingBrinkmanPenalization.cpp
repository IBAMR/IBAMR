// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
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

#include "ibamr/SOAcousticStreamingBrinkmanPenalization.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/ibtk_utilities.h"

#include <string>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{

namespace
{
inline hier::Index<NDIM>
get_shift(int dir, int shift)
{
    SAMRAI::hier::Index<NDIM> iv(0);
    iv(dir) = shift;
    return iv;
} // get_shift
} // namespace

/////////////////////////////// PUBLIC //////////////////////////////////////
SOAcousticStreamingBrinkmanPenalization::SOAcousticStreamingBrinkmanPenalization(
    std::string object_name,
    Pointer<HierarchyIntegrator> time_integrator,
    Pointer<Database> input_db,
    bool register_for_restart)
    : BrinkmanPenalizationMethod(std::move(object_name), time_integrator, input_db, register_for_restart)
{
    return;
} // SOAcousticStreamingBrinkmanPenalization

void
SOAcousticStreamingBrinkmanPenalization::setFOVelocityPatchDataIndices(int U1_real_idx, int U1_imag_idx)
{
    d_U1_real_idx = U1_real_idx;
    d_U1_imag_idx = U1_imag_idx;
    return;
} // setFOVelocityPatchDataIndices

void
SOAcousticStreamingBrinkmanPenalization::setAcousticAngularFrequency(double omega)
{
    d_acoustic_freq = omega;
    return;
} // setAcousticAngularFrequency

void
SOAcousticStreamingBrinkmanPenalization::computeBrinkmanVelocity(int b_idx, double time, int /*cycle_num*/)
{
    // Ghost fill the level set values.
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_time_integrator->getPatchHierarchy();
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent transaction_comp(d_ls_scratch_idx,
                                                       d_ls_new_idx,
                                                       "CONSERVATIVE_LINEAR_REFINE",
                                                       false,
                                                       "CONSERVATIVE_COARSEN",
                                                       "LINEAR",
                                                       false,
                                                       d_ls_bc_coef);

    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(transaction_comp, patch_hierarchy);
    hier_bdry_fill->fillData(time);

    // Set Stokes drift velocity in u_idx
    int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > finest_level = patch_hierarchy->getPatchLevel(finest_ln);
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        double vol_cell = 1.0;
        for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
        const double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

        Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(d_ls_scratch_idx);
        Pointer<SideData<NDIM, double> > b_data = patch->getPatchData(b_idx);
        Pointer<SideData<NDIM, double> > Ur = patch->getPatchData(d_U1_real_idx);
        Pointer<SideData<NDIM, double> > Ui = patch->getPatchData(d_U1_imag_idx);

        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
            {
                SideIndex<NDIM> is(it(), axis, SideIndex<NDIM>::Lower);

                const double phi_lower = (*phi_data)(is.toCell(0));
                const double phi_upper = (*phi_data)(is.toCell(1));
                const double phi = 0.5 * (phi_lower + phi_upper);
                const double Hphi = IBTK::smooth_heaviside(phi, alpha);

                if (phi <= alpha)
                {
                    const SideIndex<NDIM> is_e(it(), axis, SideIndex<NDIM>::Upper);
                    const SideIndex<NDIM> is_w(it() + get_shift(axis, -1), axis, SideIndex<NDIM>::Lower);
                    const double g_normal =
                        ((*Ui)(is) * ((*Ur)(is_e) - (*Ur)(is_w)) - (*Ur)(is) * ((*Ui)(is_e) - (*Ui)(is_w))) /
                        (2.0 * patch_dx[axis]);

                    double g_tangential = 0.0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis) continue;

                        const SideIndex<NDIM> is_n(it() + get_shift(d, 1), axis, SideIndex<NDIM>::Lower);
                        const SideIndex<NDIM> is_s(it() + get_shift(d, -1), axis, SideIndex<NDIM>::Lower);

                        const SideIndex<NDIM> is_e(it(), d, SideIndex<NDIM>::Lower);
                        const SideIndex<NDIM> is_ne(it(), d, SideIndex<NDIM>::Upper);
                        const SideIndex<NDIM> is_w(it() + get_shift(axis, -1), d, SideIndex<NDIM>::Lower);
                        const SideIndex<NDIM> is_nw(it() + get_shift(axis, -1), d, SideIndex<NDIM>::Upper);

                        g_tangential += 0.25 * ((*Ui)(is_e) + (*Ui)(is_ne) + (*Ui)(is_w) + (*Ui)(is_nw)) *
                                        ((*Ur)(is_n) - (*Ur)(is_s)) / (2 * patch_dx[d]);
                        g_tangential -= 0.25 * ((*Ur)(is_e) + (*Ur)(is_ne) + (*Ur)(is_w) + (*Ur)(is_nw)) *
                                        ((*Ui)(is_n) - (*Ui)(is_s)) / (2 * patch_dx[d]);
                    }

                    const double stokes_drift = -(g_normal + g_tangential) / (2.0 * d_acoustic_freq);
                    (*b_data)(is) = stokes_drift * d_penalty_factor * (1.0 - Hphi);
                }
            }
        }
    }

    return;
} // computeBrinkmanVelocity

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
