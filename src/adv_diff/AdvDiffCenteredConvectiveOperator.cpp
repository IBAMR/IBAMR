// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#include "ibamr/AdvDiffCenteredConvectiveOperator.h"

#include "ibamr/namespaces.h" // IWYU pragma: keep

// FORTRAN ROUTINES
#if (NDIM == 2)
#define C_TO_F_CWISE_INTERP_2ND_FC IBAMR_FC_FUNC(ctofcwiseinterp2nd2d, CTOFINTERP2ND2D)
#endif

#if (NDIM == 3)
#define C_TO_F_CWISE_INTERP_2ND_FC IBAMR_FC_FUNC(ctofcwiseinterp2nd3d, CTOFINTERP2ND3D)
#endif

extern "C"
{
    void C_TO_F_CWISE_INTERP_2ND_FC(
#if (NDIM == 2)
        double*,
        double*,
        const int&,
        const double*,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&
#endif
#if (NDIM == 3)
        double*,
        double*,
        double*,
        const int&,
        const double*,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&
#endif
    );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Minimum number of ghosts cells used for each variable quantity.
static const int Q_MIN_GCW = 1;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffCenteredConvectiveOperator::AdvDiffCenteredConvectiveOperator(std::string object_name,
                                                                     Pointer<CellVariable<NDIM, double> > Q_var,
                                                                     Pointer<Database> input_db,
                                                                     const ConvectiveDifferencingType difference_form,
                                                                     std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs)
    : CellConvectiveOperator(std::move(object_name), Q_var, Q_MIN_GCW, input_db, difference_form, std::move(bc_coefs))
{
    // intentionally blank
    return;
} // AdvDiffCenteredConvectiveOperator

void
AdvDiffCenteredConvectiveOperator::interpolateToFaceOnPatch(FaceData<NDIM, double>& q_interp_data,
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

    // Interpolate from cell centers to cell faces.
    for (int d = 0; d < Q_cell_data.getDepth(); ++d)
    {
        C_TO_F_CWISE_INTERP_2ND_FC(
#if (NDIM == 2)
            q_interp_data.getPointer(0, d),
            q_interp_data.getPointer(1, d),
            q_interp_data_gcw.min(),
            Q_cell_data.getPointer(d),
            Q_cell_data_gcw.min(),
            patch_lower(0),
            patch_upper(0),
            patch_lower(1),
            patch_upper(1)
#endif
#if (NDIM == 3)
                q_interp_data.getPointer(0, d),
            q_interp_data.getPointer(1, d),
            q_interp_data.getPointer(2, d),
            q_interp_data_gcw.min(),
            Q_cell_data.getPointer(d),
            Q_cell_data_gcw.min(),
            patch_lower(0),
            patch_upper(0),
            patch_lower(1),
            patch_upper(1),
            patch_lower(2),
            patch_upper(2)
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
