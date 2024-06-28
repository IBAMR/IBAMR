// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2023 by the IBAMR developers
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

#include "ibamr/AdvDiffCUIConvectiveOperator.h"

#include "ibamr/namespaces.h" // IWYU pragma: keep

// FORTRAN ROUTINES
#if (NDIM == 2)
#define CUI_EXTRAPOLATE_FC IBAMR_FC_FUNC_(cui_extrapolate2d, CUI_EXTRAPOLATE2D)
#endif

#if (NDIM == 3)
#define CUI_EXTRAPOLATE_FC IBAMR_FC_FUNC_(cui_extrapolate3d, CUI_EXTRAPOLATE3D)
#endif

extern "C"
{
    void CUI_EXTRAPOLATE_FC(
#if (NDIM == 2)
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double*,
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
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// NOTE: The number of ghost cells required by the advection scheme
// These values were chosen to work with CUI (the cubic interpolation
// upwind method of Waterson and Deconinck).
static const int Q_MIN_GCW = 2;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffCUIConvectiveOperator::AdvDiffCUIConvectiveOperator(std::string object_name,
                                                           Pointer<CellVariable<NDIM, double> > Q_var,
                                                           Pointer<Database> input_db,
                                                           const ConvectiveDifferencingType difference_form,
                                                           std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs)
    : CellConvectiveOperator(std::move(object_name), Q_var, Q_MIN_GCW, input_db, difference_form, std::move(bc_coefs))
{
    // intentionally blank
    return;
} // AdvDiffCUIConvectiveOperator

void
AdvDiffCUIConvectiveOperator::interpolateToFaceOnPatch(FaceData<NDIM, double>& q_interp_data,
                                                       const CellData<NDIM, double>& Q_cell_data,
                                                       const FaceData<NDIM, double>& u_data,
                                                       const Patch<NDIM>& patch)
{
    const auto& patch_box = patch.getBox();
    const auto& patch_lower = patch_box.lower();
    const auto& patch_upper = patch_box.upper();

    const IntVector<NDIM>& Q_cell_data_gcw = Q_cell_data.getGhostCellWidth();
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_cell_data_gcw.min() == Q_cell_data_gcw.max());
#endif
    const IntVector<NDIM>& u_data_gcw = u_data.getGhostCellWidth();
#if !defined(NDEBUG)
    TBOX_ASSERT(u_data_gcw.min() == u_data_gcw.max());
#endif
    const IntVector<NDIM>& q_interp_data_gcw = q_interp_data.getGhostCellWidth();
#if !defined(NDEBUG)
    TBOX_ASSERT(q_interp_data_gcw.min() == q_interp_data_gcw.max());
#endif
    const CellData<NDIM, double>& Q0_data = Q_cell_data;
    CellData<NDIM, double> Q1_data(patch_box, 1, Q_cell_data_gcw);
#if (NDIM == 3)
    CellData<NDIM, double> Q2_data(patch_box, 1, Q_cell_data_gcw);
#endif

    // Interpolate from cell centers to cell faces.
    for (int d = 0; d < Q_cell_data.getDepth(); ++d)
    {
        CUI_EXTRAPOLATE_FC(
#if (NDIM == 2)
            patch_lower(0),
            patch_upper(0),
            patch_lower(1),
            patch_upper(1),
            Q_cell_data_gcw(0),
            Q_cell_data_gcw(1),
            Q0_data.getPointer(d),
            Q1_data.getPointer(),
            u_data_gcw(0),
            u_data_gcw(1),
            q_interp_data_gcw(0),
            q_interp_data_gcw(1),
            u_data.getPointer(0),
            u_data.getPointer(1),
            q_interp_data.getPointer(0, d),
            q_interp_data.getPointer(1, d)
#endif
#if (NDIM == 3)
                patch_lower(0),
            patch_upper(0),
            patch_lower(1),
            patch_upper(1),
            patch_lower(2),
            patch_upper(2),
            Q_cell_data_gcw(0),
            Q_cell_data_gcw(1),
            Q_cell_data_gcw(2),
            Q0_data.getPointer(d),
            Q1_data.getPointer(),
            Q2_data.getPointer(),
            u_data_gcw(0),
            u_data_gcw(1),
            u_data_gcw(2),
            q_interp_data_gcw(0),
            q_interp_data_gcw(1),
            q_interp_data_gcw(2),
            u_data.getPointer(0),
            u_data.getPointer(1),
            u_data.getPointer(2),
            q_interp_data.getPointer(0, d),
            q_interp_data.getPointer(1, d),
            q_interp_data.getPointer(2, d)
#endif
        );
    }
    return;
} // interpolateToFaceOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
