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
#include "ibamr/INSVCStaggeredConservativeMassMomentumRKIntegrator.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/IndexUtilities.h"

#include "BasePatchHierarchy.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CoarseFineBoundary.h"
#include "FaceData.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchySideDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

#include <array>
#include <limits>
#include <memory>
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

// FORTRAN ROUTINES
#if (NDIM == 2)
#define CONVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(convect_derivative2d, CONVECT_DERIVATIVE2D)
#define VC_UPDATE_DENSITY_FC IBAMR_FC_FUNC_(vc_update_density2d, VC_UPDATE_DENSITY2D)
#define VC_MASS_CONSERVATION_ERROR_FC IBAMR_FC_FUNC_(vc_mass_conservation_error2d, VC_MASS_CONSERVATION_ERROR2D)
#define NAVIER_STOKES_INTERP_COMPS_FC IBAMR_FC_FUNC_(navier_stokes_interp_comps2d, NAVIER_STOKES_INTERP_COMPS2D)
#define VC_NAVIER_STOKES_UPWIND_QUANTITY_FC                                                                            \
    IBAMR_FC_FUNC_(vc_navier_stokes_upwind_quantity2d, VC_NAVIER_STOKES_UPWIND_QUANTITY2D)
#define VC_NAVIER_STOKES_CUI_QUANTITY_FC                                                                               \
    IBAMR_FC_FUNC_(vc_navier_stokes_cui_quantity2d, VC_NAVIER_STOKES_CUI_QUANTITY2D)
#define VC_NAVIER_STOKES_FBICS_QUANTITY_FC                                                                             \
    IBAMR_FC_FUNC_(vc_navier_stokes_fbics_quantity2d, VC_NAVIER_STOKES_FBICS_QUANTITY2D)
#define VC_NAVIER_STOKES_MGAMMA_QUANTITY_FC                                                                            \
    IBAMR_FC_FUNC_(vc_navier_stokes_mgamma_quantity2d, VC_NAVIER_STOKES_MGAMMA_QUANTITY2D)
#define GODUNOV_EXTRAPOLATE_FC IBAMR_FC_FUNC_(godunov_extrapolate2d, GODUNOV_EXTRAPOLATE2D)
#define VC_NAVIER_STOKES_COMPUTE_MOMENTUM_FC                                                                           \
    IBAMR_FC_FUNC_(vc_navier_stokes_compute_momentum2d, VC_NAVIER_STOKES_COMPUTE_MOMENTUM2D)
#endif

#if (NDIM == 3)
#define CONVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(convect_derivative3d, CONVECT_DERIVATIVE3D)
#define VC_UPDATE_DENSITY_FC IBAMR_FC_FUNC_(vc_update_density3d, VC_UPDATE_DENSITY3D)
#define VC_MASS_CONSERVATION_ERROR_FC IBAMR_FC_FUNC_(vc_mass_conservation_error3d, VC_MASS_CONSERVATION_ERROR3D)
#define NAVIER_STOKES_INTERP_COMPS_FC IBAMR_FC_FUNC_(navier_stokes_interp_comps3d, NAVIER_STOKES_INTERP_COMPS3D)
#define VC_NAVIER_STOKES_UPWIND_QUANTITY_FC                                                                            \
    IBAMR_FC_FUNC_(vc_navier_stokes_upwind_quantity3d, VC_NAVIER_STOKES_UPWIND_QUANTITY3D)
#define VC_NAVIER_STOKES_CUI_QUANTITY_FC                                                                               \
    IBAMR_FC_FUNC_(vc_navier_stokes_cui_quantity3d, VC_NAVIER_STOKES_CUI_QUANTITY3D)
#define VC_NAVIER_STOKES_FBICS_QUANTITY_FC                                                                             \
    IBAMR_FC_FUNC_(vc_navier_stokes_fbics_quantity3d, VC_NAVIER_STOKES_FBICS_QUANTITY3D)
#define VC_NAVIER_STOKES_MGAMMA_QUANTITY_FC                                                                            \
    IBAMR_FC_FUNC_(vc_navier_stokes_mgamma_quantity3d, VC_NAVIER_STOKES_MGAMMA_QUANTITY3D)
#define GODUNOV_EXTRAPOLATE_FC IBAMR_FC_FUNC_(godunov_extrapolate3d, GODUNOV_EXTRAPOLATE3D)
#define VC_NAVIER_STOKES_COMPUTE_MOMENTUM_FC                                                                           \
    IBAMR_FC_FUNC_(vc_navier_stokes_compute_momentum3d, VC_NAVIER_STOKES_COMPUTE_MOMENTUM3D)
#endif

extern "C"
{
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
    void VC_UPDATE_DENSITY_FC(const double*,
                              const double&,
                              const double&,
                              const double&,
                              const double&,
#if (NDIM == 2)
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const double*,
                              const int&,
                              const int&,
                              const double*,
                              const int&,
                              const int&,
                              const double*,
                              const double*,
                              const int&,
                              const int&,
                              const double*,
                              const double*,
                              const int&,
                              const int&,
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
                              const double*,
                              const int&,
                              const int&,
                              const int&,
                              const double*,
                              const int&,
                              const int&,
                              const int&,
                              const double*,
                              const double*,
                              const double*,
                              const int&,
                              const int&,
                              const int&,
                              const double*,
                              const double*,
                              const double*,
                              const int&,
                              const int&,
                              const int&,
                              const double*,
                              const int&,
                              const int&,
                              const int&,

#endif
                              double*);

    void VC_MASS_CONSERVATION_ERROR_FC(const double*,
                                       const double&,
#if (NDIM == 2)
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const double*,
                                       const int&,
                                       const int&,
                                       const double*,
                                       const int&,
                                       const int&,
                                       const double*,
                                       const double*,
                                       const int&,
                                       const int&,
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
                                       const double*,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const double*,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const double*,
                                       const double*,
                                       const double*,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const double*,
                                       const double*,
                                       const double*,
                                       const int&,
                                       const int&,
                                       const int&,

#endif
                                       double*);

    void VC_SSP_RK2_UPDATE_DENSITY_FC(const double*,
                                      const double&,
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
                                      const double*,
                                      const double*,
                                      const int&,
                                      const int&,
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
                                      const double*,
                                      const double*,
                                      const double*,
                                      const int&,
                                      const int&,
                                      const int&,
                                      const double*,
                                      const double*,
                                      const double*,
                                      const int&,
                                      const int&,
                                      const int&,
                                      const double*,
                                      const int&,
                                      const int&,
                                      const int&,
#endif
                                      double*);

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

    void VC_NAVIER_STOKES_UPWIND_QUANTITY_FC(
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
        const double*,
        const double*,
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
        const double*,
        const double*,
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
        const double*,
        const double*,
        const double*,
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
        const double*,
        const double*,
        const double*,
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
        const double*,
        const double*,
        const double*,
        const int&,
        const int&,
        const int&,
        double*,
        double*,
        double*
#endif
    );

    void VC_NAVIER_STOKES_CUI_QUANTITY_FC(
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
        const double*,
        const double*,
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
        const double*,
        const double*,
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
        const double*,
        const double*,
        const double*,
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
        const double*,
        const double*,
        const double*,
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
        const double*,
        const double*,
        const double*,
        const int&,
        const int&,
        const int&,
        double*,
        double*,
        double*
#endif
    );

    void VC_NAVIER_STOKES_FBICS_QUANTITY_FC(
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
        const double*,
        const double*,
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
        const double*,
        const double*,
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
        const double*,
        const double*,
        const double*,
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
        const double*,
        const double*,
        const double*,
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
        const double*,
        const double*,
        const double*,
        const int&,
        const int&,
        const int&,
        double*,
        double*,
        double*
#endif
    );

    void VC_NAVIER_STOKES_MGAMMA_QUANTITY_FC(
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
        const double*,
        const double*,
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
        const double*,
        const double*,
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
        const double*,
        const double*,
        const double*,
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
        const double*,
        const double*,
        const double*,
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
        const double*,
        const double*,
        const double*,
        const int&,
        const int&,
        const int&,
        double*,
        double*,
        double*
#endif
    );

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

    void VC_NAVIER_STOKES_COMPUTE_MOMENTUM_FC(
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
        const double*,
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
        const double*,
        const double*,
        const double*
#endif
    );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// NOTE: The number of ghost cells required by the convection scheme depends
// on the chosen convective limiter, which will be set via input file
static const int GUPWINDG = 2;
static const int GCUIG = 3;
static const int GFBICSG = 3;
static const int GMGAMMAG = 3;
static const int GPPMG = 4;
static const int NOGHOSTS = 0;

// Number of ghost cells to fill at coarse fine interface to enforce divergence
// free condition
static const int CF_GHOST_WIDTH = 1;

// Timers.
static Timer* t_apply_convective_operator;
static Timer* t_integrate;
static Timer* t_initialize_integrator;
static Timer* t_deallocate_integrator;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSVCStaggeredConservativeMassMomentumRKIntegrator::INSVCStaggeredConservativeMassMomentumRKIntegrator(
    std::string object_name,
    Pointer<Database> input_db)
    : STSMassFluxIntegrator(object_name, input_db)
{
    if (input_db)
    {
        if (input_db->keyExists("bdry_extrap_type"))
        {
            d_velocity_bdry_extrap_type = input_db->getString("bdry_extrap_type");
        }
        if (input_db->keyExists("velocity_bdry_extrap_type"))
        {
            d_velocity_bdry_extrap_type = input_db->getString("velocity_bdry_extrap_type");
        }
        if (input_db->keyExists("convective_limiter"))
        {
            d_velocity_convective_limiter = string_to_enum<LimiterType>(input_db->getString("convective_limiter"));
        }
        if (input_db->keyExists("velocity_convective_limiter"))
        {
            d_velocity_convective_limiter =
                string_to_enum<LimiterType>(input_db->getString("velocity_convective_limiter"));
        }
    }

    switch (d_density_convective_limiter)
    {
    case UPWIND:
        d_density_limiter_gcw = GUPWINDG;
        break;
    case CUI:
        d_density_limiter_gcw = GCUIG;
        break;
    case FBICS:
        d_density_limiter_gcw = GFBICSG;
        break;
    case MGAMMA:
        d_density_limiter_gcw = GMGAMMAG;
        break;
    case PPM:
        d_density_limiter_gcw = GPPMG;
        break;
    default:
        TBOX_ERROR(
            "INSVCStaggeredConservativeMassMomentumRKIntegrator::"
            "INSVCStaggeredConservativeMassMomentumRKIntegrator():\n"
            << "  unsupported density convective limiter: "
            << IBAMR::enum_to_string<LimiterType>(d_density_convective_limiter) << " \n"
            << "  valid choices are: UPWIND, CUI, FBICS, MGAMMA, PPM\n");
    }

    switch (d_velocity_convective_limiter)
    {
    case UPWIND:
        d_velocity_limiter_gcw = GUPWINDG;
        break;
    case CUI:
        d_velocity_limiter_gcw = GCUIG;
        break;
    case FBICS:
        d_velocity_limiter_gcw = GFBICSG;
        break;
    case MGAMMA:
        d_velocity_limiter_gcw = GMGAMMAG;
        break;
    case PPM:
        d_velocity_limiter_gcw = GPPMG;
        break;
    default:
        TBOX_ERROR(
            "INSVCStaggeredConservativeMassMomentumRKIntegrator::"
            "INSVCStaggeredConservativeMassMomentumRKIntegrator():\n"
            << "  unsupported velocity convective limiter: "
            << IBAMR::enum_to_string<LimiterType>(d_velocity_convective_limiter) << " \n"
            << "  valid choices are: UPWIND, CUI, FBICS, MGAMMA, PPM\n");
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context =
        var_db->getContext("INSVCStaggeredConservativeMassMomentumRKIntegrator::CONTEXT");

    d_V_var = new SideVariable<NDIM, double>(d_object_name + "::V");
    d_V_scratch_idx = var_db->registerVariableAndContext(
        d_V_var, var_db->getContext(d_object_name + "::V::SCRATCH"), IntVector<NDIM>(d_velocity_limiter_gcw));
    d_V_composite_idx = var_db->registerVariableAndContext(
        d_V_var, var_db->getContext(d_object_name + "::V::COMPOSITE"), IntVector<NDIM>(NOGHOSTS));
#if !defined(NDEBUG)
    TBOX_ASSERT(d_V_scratch_idx >= 0);
    TBOX_ASSERT(d_V_composite_idx >= 0);
#endif

    d_rho_var = new SideVariable<NDIM, double>(d_object_name + "::RHO_SC");
    d_rho_scratch_idx = var_db->registerVariableAndContext(
        d_rho_var, var_db->getContext(d_object_name + "::RHO_SC::SCRATCH"), IntVector<NDIM>(d_density_limiter_gcw));
    d_rho_composite_idx = var_db->registerVariableAndContext(
        d_rho_var, var_db->getContext(d_object_name + "::RHO_SC::COMPOSITE"), IntVector<NDIM>(NOGHOSTS));
    d_rho_new_idx = var_db->registerVariableAndContext(
        d_rho_var, var_db->getContext(d_object_name + "::RHO_SC::NEW"), IntVector<NDIM>(NOGHOSTS));
#if !defined(NDEBUG)
    TBOX_ASSERT(d_rho_scratch_idx >= 0);
    TBOX_ASSERT(d_rho_new_idx >= 0);
#endif

    d_S_var = new SideVariable<NDIM, double>(d_object_name + "::S");
    d_S_scratch_idx = var_db->registerVariableAndContext(d_S_var, context, IntVector<NDIM>(NOGHOSTS));

#if !defined(NDEBUG)
    TBOX_ASSERT(d_S_scratch_idx >= 0);
#endif

    d_E_var = new SideVariable<NDIM, double>(d_object_name + "::E");
    d_E_scratch_idx = var_db->registerVariableAndContext(d_E_var, context, IntVector<NDIM>(NOGHOSTS));

#if !defined(NDEBUG)
    TBOX_ASSERT(d_E_scratch_idx >= 0);
#endif

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply_convective_operator =
                      TimerManager::getManager()->getTimer("IBAMR::INSVCStaggeredConservativeMassMomentumRKIntegrator::"
                                                           "applyConvectiveOperator()");
                  t_integrate = TimerManager::getManager()->getTimer(
                      "IBAMR::INSVCStaggeredConservativeMassMomentumRKIntegrator::integrate("
                      ")");
                  t_initialize_integrator =
                      TimerManager::getManager()->getTimer("IBAMR::INSVCStaggeredConservativeMassMomentumRKIntegrator::"
                                                           "initializeSTSIntegrator()");
                  t_deallocate_integrator =
                      TimerManager::getManager()->getTimer("IBAMR::INSVCStaggeredConservativeMassMomentumRKIntegrator::"
                                                           "deallocateSTSIntegrator()"););
    return;
} // INSVCStaggeredConservativeMassMomentumRKIntegrator

INSVCStaggeredConservativeMassMomentumRKIntegrator::~INSVCStaggeredConservativeMassMomentumRKIntegrator()
{
    deallocateSTSIntegrator();
    return;
} // ~INSVCStaggeredConservativeMassMomentumRKIntegrator

void
INSVCStaggeredConservativeMassMomentumRKIntegrator::integrate(double dt)
{
    // Get hierarchy operation object
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    d_hier_sc_data_ops =
        hier_ops_manager->getOperationsDouble(new SideVariable<NDIM, double>("sc_var"), d_hierarchy, true);

    IBAMR_TIMER_START(t_integrate)
#if !defined(NDEBUG)
    if (!d_is_initialized)
    {
        TBOX_ERROR("INSVCStaggeredConservativeMassMomentumRKIntegrator::integrate():\n"
                   << "  time integrator must be initialized prior to call to "
                      "integrate()\n");
    }

    TBOX_ASSERT(d_rho_current_idx >= 0);
    TBOX_ASSERT(d_V_old_idx >= 0);
    TBOX_ASSERT(d_V_current_idx >= 0);
    TBOX_ASSERT(d_V_new_idx >= 0);
#endif

#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(dt, getTimeStepSize()));
#endif

    if (d_V_old_idx == d_V_current_idx)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_dt_prev <= 0.0);
#endif
        // Ensure that previous time step is set for initial times
        d_dt_prev = dt;
    }
#if !defined(NDEBUG)
    if (!(dt > 0.0))
    {
        TBOX_ERROR("INSVCStaggeredConservativeMassMomentumRKIntegrator::integrate():\n"
                   << " invalid time step size dt = " << dt << "\n");
    }
#endif

// Assertions for velocity interpolation and extrapolation
#if !defined(NDEBUG)
    if (d_cycle_num < 0)
    {
        TBOX_ERROR("INSVCStaggeredConservativeMassMomentumRKIntegrator::integrate():\n"
                   << "  invalid cycle number = " << d_cycle_num << "\n");
    }
    if (d_dt_prev <= 0.0 && d_density_time_stepping_type != FORWARD_EULER)
    {
        TBOX_ERROR("INSVCStaggeredConservativeMassMomentumRKIntegrator::integrate():\n"
                   << "  invalid previous time step size = " << d_dt_prev << "\n");
    }
#endif

    // Fill ghost cell values
    static const bool homogeneous_bc = false;
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;

    // Fill ghost cells for current density
    std::vector<InterpolationTransactionComponent> rho_transaction_comps({ { d_rho_scratch_idx,
                                                                             d_rho_current_idx,
                                                                             "CONSERVATIVE_LINEAR_REFINE",
                                                                             false,
                                                                             "CONSERVATIVE_COARSEN",
                                                                             d_density_bdry_extrap_type,
                                                                             false,
                                                                             d_rho_bc_coefs } });
    d_hier_rho_bdry_fill->resetTransactionComponents(rho_transaction_comps);
    d_hier_rho_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_rho_bdry_fill->fillData(d_current_time);
    d_hier_rho_bdry_fill->resetTransactionComponents(d_rho_transaction_comps);

    // Fill ghost cells for the velocity used to compute the density update.
    // Note, enforce divergence free condition on all physical boundaries to
    // ensure boundedness of density update.
    d_hier_sc_data_ops->copyData(d_V_composite_idx,
                                 d_V_current_idx,
                                 /*interior_only*/ true);
    std::vector<InterpolationTransactionComponent> v_transaction_comps({ { d_V_scratch_idx,
                                                                           d_V_composite_idx,
                                                                           "CONSERVATIVE_LINEAR_REFINE",
                                                                           false,
                                                                           "CONSERVATIVE_COARSEN",
                                                                           d_velocity_bdry_extrap_type,
                                                                           false,
                                                                           d_u_bc_coefs } });
    d_hier_v_bdry_fill->resetTransactionComponents(v_transaction_comps);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
        d_u_bc_coefs, nullptr, d_V_scratch_idx, -1, homogeneous_bc);
    d_hier_v_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_v_bdry_fill->fillData(d_current_time);
    d_bc_helper->enforceDivergenceFreeConditionAtBoundary(
        d_V_scratch_idx, d_coarsest_ln, d_finest_ln, StaggeredStokesPhysicalBoundaryHelper::ALL_BDRY);
    enforceDivergenceFreeConditionAtCoarseFineInterface(d_V_scratch_idx);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_u_bc_coefs, nullptr);
    d_hier_v_bdry_fill->resetTransactionComponents(d_v_transaction_comps);

    // Compute the old mass
    const int wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
    const double old_mass = d_hier_sc_data_ops->integral(d_rho_current_idx, wgt_sc_idx);
    if (d_enable_logging)
    {
        plog << "INSVCStaggeredConservativeMassMomentumRKIntegrator::integrate(): "
                "old mass in the domain = "
             << old_mass << "\n";
    }

    double eval_time = d_current_time;
    // When cycle_num > 0, midpoint rule (RK-2) is applied for density integration.
    if (d_cycle_num > 0)
    {
        eval_time = d_current_time + dt / 2.0;

        // Fill ghost cells for new density and velocity, if needed
        d_hier_sc_data_ops->linearSum(d_rho_composite_idx, 0.5, d_rho_current_idx, 0.5, d_rho_new_idx);
        std::vector<InterpolationTransactionComponent> update_transaction_comps({ { d_rho_scratch_idx,
                                                                                    d_rho_composite_idx,
                                                                                    "CONSERVATIVE_LINEAR_REFINE",
                                                                                    false,
                                                                                    "CONSERVATIVE_COARSEN",
                                                                                    d_density_bdry_extrap_type,
                                                                                    false,
                                                                                    d_rho_bc_coefs } });
        d_hier_rho_bdry_fill->resetTransactionComponents(update_transaction_comps);
        d_hier_rho_bdry_fill->setHomogeneousBc(homogeneous_bc);
        d_hier_rho_bdry_fill->fillData(eval_time);
        d_hier_rho_bdry_fill->resetTransactionComponents(d_rho_transaction_comps);

        // Compute an approximation to velocity at eval_time. Note, enforce
        // divergence free condition on all physical boundaries to ensure
        // boundedness of density update.
        d_hier_sc_data_ops->linearSum(
            d_V_composite_idx, 0.5, d_V_current_idx, 0.5, d_V_new_idx, /*interior_only*/ true);

        std::vector<InterpolationTransactionComponent> v_update_transaction_comps({ { d_V_scratch_idx,
                                                                                      d_V_composite_idx,
                                                                                      "CONSERVATIVE_LINEAR_REFINE",
                                                                                      false,
                                                                                      "CONSERVATIVE_COARSEN",
                                                                                      d_velocity_bdry_extrap_type,
                                                                                      false,
                                                                                      d_u_bc_coefs } });
        d_hier_v_bdry_fill->resetTransactionComponents(v_update_transaction_comps);
        StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
            d_u_bc_coefs, nullptr, d_V_scratch_idx, -1, homogeneous_bc);
        d_hier_v_bdry_fill->setHomogeneousBc(homogeneous_bc);
        d_hier_v_bdry_fill->fillData(eval_time);
        d_bc_helper->enforceDivergenceFreeConditionAtBoundary(
            d_V_scratch_idx, d_coarsest_ln, d_finest_ln, StaggeredStokesPhysicalBoundaryHelper::ALL_BDRY);
        enforceDivergenceFreeConditionAtCoarseFineInterface(d_V_scratch_idx);
        StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_u_bc_coefs, nullptr);
        d_hier_v_bdry_fill->resetTransactionComponents(d_v_transaction_comps);
    }

    // Compute the source term.
    if (d_S_fcn)
    {
        d_S_fcn->setDataOnPatchHierarchy(d_S_scratch_idx, d_S_var, d_hierarchy, eval_time);
    }
    else
    {
        d_hier_sc_data_ops->setToScalar(d_S_scratch_idx, 0.0);
    }

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

            Pointer<SideData<NDIM, double> > N_data = patch->getPatchData(d_N_idx);
            Pointer<SideData<NDIM, double> > V_data = patch->getPatchData(d_V_scratch_idx);
            Pointer<SideData<NDIM, double> > R_cur_data = patch->getPatchData(d_rho_current_idx);
            Pointer<SideData<NDIM, double> > R_pre_data = patch->getPatchData(d_rho_scratch_idx);
            Pointer<SideData<NDIM, double> > R_new_data = patch->getPatchData(d_rho_new_idx);
            Pointer<SideData<NDIM, double> > R_src_data = patch->getPatchData(d_S_scratch_idx);
            Pointer<SideData<NDIM, double> > E_data = patch->getPatchData(d_E_scratch_idx);

            // Define variables that live on the "faces" of control
            // volumes centered about side-centered staggered velocity
            // components
            const IntVector<NDIM> ghosts = IntVector<NDIM>(1);
            std::array<Box<NDIM>, NDIM> side_boxes;
            std::array<Pointer<FaceData<NDIM, double> >, NDIM> V_adv_data;
            std::array<Pointer<FaceData<NDIM, double> >, NDIM> V_half_data;
            std::array<Pointer<FaceData<NDIM, double> >, NDIM> R_half_data;
            std::array<Pointer<FaceData<NDIM, double> >, NDIM> P_half_data;
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                side_boxes[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                V_adv_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                V_half_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                R_half_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                P_half_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
            }
            // Interpolate velocity components onto "faces" using simple averages.
            computeAdvectionVelocity(V_adv_data, V_data, patch_lower, patch_upper, side_boxes);

            // Upwind side-centered densities onto faces.
            interpolateSideQuantity(R_half_data,
                                    V_adv_data,
                                    R_pre_data,
                                    patch_lower,
                                    patch_upper,
                                    side_boxes,
                                    d_density_convective_limiter);

            // Upwind side-centered velocities onto faces. To be used while computing
            // convective derivative.
            interpolateSideQuantity(
                V_half_data, V_adv_data, V_data, patch_lower, patch_upper, side_boxes, d_velocity_convective_limiter);

            IBAMR_TIMER_START(t_apply_convective_operator);

            // Compute the convective derivative with the penultimate density and
            // velocity, if necessary.
            computeConvectiveDerivative(
                N_data, P_half_data, V_adv_data, R_half_data, V_half_data, side_boxes, dx, patch);

            IBAMR_TIMER_STOP(t_apply_convective_operator);

            // Compute the updated density.
            const double a0 = 1.0;
            const double a1 = 0.0;
            const double a2 = 1.0;
            computeDensityUpdate(R_new_data,
                                 a0,
                                 R_cur_data,
                                 a1,
                                 R_pre_data,
                                 a2,
                                 V_adv_data,
                                 R_half_data,
                                 R_src_data,
                                 side_boxes,
                                 dt,
                                 dx);

            Pointer<SideData<NDIM, double> > V_cur_data = patch->getPatchData(d_V_current_idx);
            computeErrorOfMassConservationEquation(
                E_data, R_new_data, R_cur_data, V_adv_data, R_half_data, side_boxes, dt, dx);

            // subtract Error*U from the convective operator.
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);

                    (*N_data)(si) -= (*V_cur_data)(si) * (*E_data)(si);
                }
            }
        }
    }

    // Refill boundary values of newest density
    const double new_time = d_current_time + dt;
    std::vector<InterpolationTransactionComponent> new_transaction_comps({ { d_rho_scratch_idx,
                                                                             d_rho_new_idx,
                                                                             "CONSERVATIVE_LINEAR_REFINE",
                                                                             false,
                                                                             "CONSERVATIVE_COARSEN",
                                                                             d_density_bdry_extrap_type,
                                                                             false,
                                                                             d_rho_bc_coefs } });
    d_hier_rho_bdry_fill->resetTransactionComponents(new_transaction_comps);
    d_hier_rho_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_rho_bdry_fill->fillData(new_time);
    d_hier_rho_bdry_fill->resetTransactionComponents(d_rho_transaction_comps);

    d_hier_sc_data_ops->copyData(d_rho_new_idx,
                                 d_rho_scratch_idx,
                                 /*interior_only*/ true);

    // Compute the new mass
    const double new_mass = d_hier_sc_data_ops->integral(d_rho_new_idx, wgt_sc_idx);
    if (d_enable_logging)
    {
        plog << "INSVCStaggeredConservativeMassMomentumRKIntegrator::integrate(): "
                "new mass in the domain = "
             << new_mass << "\n";
        plog << "INSVCStaggeredConservativeMassMomentumRKIntegrator::integrate(): "
                "change in mass = "
             << new_mass - old_mass << "\n";
    }

    // Reset select options
    d_N_idx = -1;
    d_rho_current_idx = -1;
    d_V_old_idx = -1;
    d_V_current_idx = -1;
    d_V_new_idx = -1;
    d_cycle_num = -1;
    d_dt_prev = -1.0;

    IBAMR_TIMER_STOP(t_integrate);
    return;
} // integrate

void
INSVCStaggeredConservativeMassMomentumRKIntegrator::initializeSTSIntegrator(
    Pointer<BasePatchHierarchy<NDIM> > base_hierarchy)
{
    IBAMR_TIMER_START(t_initialize_integrator);

    // Get the hierarchy configuration.
    Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    d_hierarchy = hierarchy;

    if (d_is_initialized) deallocateSTSIntegrator();

    // Update the level numbers.
    d_coarsest_ln = 0;
    d_finest_ln = d_hierarchy->getFinestLevelNumber();

    // Setup the interpolation transaction information.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    d_rho_transaction_comps.resize(1);
    d_rho_transaction_comps[0] = InterpolationTransactionComponent(d_rho_scratch_idx,
                                                                   d_rho_composite_idx,
                                                                   "CONSERVATIVE_LINEAR_REFINE",
                                                                   false,
                                                                   "CONSERVATIVE_COARSEN",
                                                                   d_density_bdry_extrap_type,
                                                                   false,
                                                                   d_rho_bc_coefs);

    d_v_transaction_comps.resize(1);
    d_v_transaction_comps[0] = InterpolationTransactionComponent(d_V_scratch_idx,
                                                                 d_V_composite_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 false,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 d_velocity_bdry_extrap_type,
                                                                 false,
                                                                 d_u_bc_coefs);

    // Initialize the interpolation operators.
    d_hier_rho_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_rho_bdry_fill->initializeOperatorState(d_rho_transaction_comps, d_hierarchy);
    d_hier_v_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_v_bdry_fill->initializeOperatorState(d_v_transaction_comps, d_hierarchy);

    // Initialize the BC helper.
    d_bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
    d_bc_helper->cacheBcCoefData(d_u_bc_coefs, d_solution_time, d_hierarchy);

    // Create the coarse-fine boundary boxes.
    d_cf_boundary.resize(d_finest_ln + 1);
    const IntVector<NDIM>& max_ghost_width = CF_GHOST_WIDTH;
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        d_cf_boundary[ln] = CoarseFineBoundary<NDIM>(*d_hierarchy, ln, max_ghost_width);
    }

    // Allocate data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_V_scratch_idx)) level->allocatePatchData(d_V_scratch_idx);
        if (!level->checkAllocated(d_V_composite_idx)) level->allocatePatchData(d_V_composite_idx);
        if (!level->checkAllocated(d_rho_scratch_idx)) level->allocatePatchData(d_rho_scratch_idx);
        if (!level->checkAllocated(d_rho_new_idx)) level->allocatePatchData(d_rho_new_idx);
        if (!level->checkAllocated(d_rho_composite_idx)) level->allocatePatchData(d_rho_composite_idx);
        if (!level->checkAllocated(d_S_scratch_idx)) level->allocatePatchData(d_S_scratch_idx);
        if (!level->checkAllocated(d_E_scratch_idx)) level->allocatePatchData(d_E_scratch_idx);
    }

    if (!d_hier_math_ops_external)
    {
        d_hier_math_ops = new HierarchyMathOps("INSVCStaggeredConservativeMassMomentumRKIntegrator::HierarchyMathOps",
                                               d_hierarchy,
                                               d_coarsest_ln,
                                               d_finest_ln);
    }
    else
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_hier_math_ops);
#endif
    }

    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_integrator);
    return;
} // initializeSTSIntegrator

void
INSVCStaggeredConservativeMassMomentumRKIntegrator::deallocateSTSIntegrator()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_integrator);

    // Deallocate the communications operators and BC helpers.
    d_hier_rho_bdry_fill.setNull();
    d_hier_v_bdry_fill.setNull();
    d_bc_helper.setNull();

    // Deallocate data.
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_V_scratch_idx)) level->deallocatePatchData(d_V_scratch_idx);
        if (level->checkAllocated(d_V_composite_idx)) level->deallocatePatchData(d_V_composite_idx);
        if (level->checkAllocated(d_rho_scratch_idx)) level->deallocatePatchData(d_rho_scratch_idx);
        if (level->checkAllocated(d_rho_new_idx)) level->deallocatePatchData(d_rho_new_idx);
        if (level->checkAllocated(d_rho_composite_idx)) level->deallocatePatchData(d_rho_composite_idx);
        if (level->checkAllocated(d_S_scratch_idx)) level->deallocatePatchData(d_S_scratch_idx);
        if (level->checkAllocated(d_E_scratch_idx)) level->deallocatePatchData(d_E_scratch_idx);
    }

    // Deallocate coarse-fine boundary object.
    d_cf_boundary.clear();

    // Deallocate hierarchy math operations object.
    if (!d_hier_math_ops_external) d_hier_math_ops.setNull();

    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_integrator);
    return;
} // deallocateOperatorState

void
INSVCStaggeredConservativeMassMomentumRKIntegrator::setSideCenteredVelocityBoundaryConditions(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(u_bc_coefs.size() == NDIM);
#endif
    d_u_bc_coefs = u_bc_coefs;
    return;
} // setSideCenteredVelocityBoundaryConditions

int
INSVCStaggeredConservativeMassMomentumRKIntegrator::getUpdatedSideCenteredDensityPatchDataIndex()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_rho_new_idx >= 0);
#endif
    return d_rho_new_idx;
} // getUpdatedSideCenteredDensityPatchDataIndex

void
INSVCStaggeredConservativeMassMomentumRKIntegrator::setMassDensitySourceTerm(const Pointer<CartGridFunction> S_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(S_fcn);
#endif
    d_S_fcn = S_fcn;
    return;
} // setMassDensitySourceTerm
/////////////////////////////// PROTECTED ////////////////////////////////////

void
INSVCStaggeredConservativeMassMomentumRKIntegrator::computeAdvectionVelocity(
    std::array<Pointer<FaceData<NDIM, double> >, NDIM> U_adv_data,
    const Pointer<SideData<NDIM, double> > U_data,
    const IntVector<NDIM>& patch_lower,
    const IntVector<NDIM>& patch_upper,
    const std::array<Box<NDIM>, NDIM>& side_boxes)
{
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
} // computeAdvectionVelocity

void
INSVCStaggeredConservativeMassMomentumRKIntegrator::interpolateSideQuantity(
    std::array<Pointer<FaceData<NDIM, double> >, NDIM> Q_half_data,
    const std::array<Pointer<FaceData<NDIM, double> >, NDIM> U_adv_data,
    const Pointer<SideData<NDIM, double> > Q_data,
    const IntVector<NDIM>& patch_lower,
    const IntVector<NDIM>& patch_upper,
    const std::array<Box<NDIM>, NDIM>& side_boxes,
    const LimiterType& convective_limiter)
{
    switch (convective_limiter)
    {
    case UPWIND:
        // Upwind side-centered densities onto faces.
        VC_NAVIER_STOKES_UPWIND_QUANTITY_FC(patch_lower(0),
                                            patch_upper(0),
                                            patch_lower(1),
                                            patch_upper(1),
#if (NDIM == 3)
                                            patch_lower(2),
                                            patch_upper(2),
#endif
                                            Q_data->getGhostCellWidth()(0),
                                            Q_data->getGhostCellWidth()(1),
#if (NDIM == 3)
                                            Q_data->getGhostCellWidth()(2),
#endif
                                            Q_data->getPointer(0),
                                            Q_data->getPointer(1),
#if (NDIM == 3)
                                            Q_data->getPointer(2),
#endif
                                            side_boxes[0].lower(0),
                                            side_boxes[0].upper(0),
                                            side_boxes[0].lower(1),
                                            side_boxes[0].upper(1),
#if (NDIM == 3)
                                            side_boxes[0].lower(2),
                                            side_boxes[0].upper(2),
#endif
                                            U_adv_data[0]->getGhostCellWidth()(0),
                                            U_adv_data[0]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                            U_adv_data[0]->getGhostCellWidth()(2),
#endif
                                            U_adv_data[0]->getPointer(0),
                                            U_adv_data[0]->getPointer(1),
#if (NDIM == 3)
                                            U_adv_data[0]->getPointer(2),
#endif
                                            Q_half_data[0]->getGhostCellWidth()(0),
                                            Q_half_data[0]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                            Q_half_data[0]->getGhostCellWidth()(2),
#endif
                                            Q_half_data[0]->getPointer(0),
                                            Q_half_data[0]->getPointer(1),
#if (NDIM == 3)
                                            Q_half_data[0]->getPointer(2),
#endif
                                            side_boxes[1].lower(0),
                                            side_boxes[1].upper(0),
                                            side_boxes[1].lower(1),
                                            side_boxes[1].upper(1),
#if (NDIM == 3)
                                            side_boxes[1].lower(2),
                                            side_boxes[1].upper(2),
#endif
                                            U_adv_data[1]->getGhostCellWidth()(0),
                                            U_adv_data[1]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                            U_adv_data[1]->getGhostCellWidth()(2),
#endif
                                            U_adv_data[1]->getPointer(0),
                                            U_adv_data[1]->getPointer(1),
#if (NDIM == 3)
                                            U_adv_data[1]->getPointer(2),
#endif
                                            Q_half_data[1]->getGhostCellWidth()(0),
                                            Q_half_data[1]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                            Q_half_data[1]->getGhostCellWidth()(2),
#endif
                                            Q_half_data[1]->getPointer(0),
                                            Q_half_data[1]->getPointer(1)
#if (NDIM == 3)
                                                ,
                                            Q_half_data[1]->getPointer(2),
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
                                            Q_half_data[2]->getGhostCellWidth()(0),
                                            Q_half_data[2]->getGhostCellWidth()(1),
                                            Q_half_data[2]->getGhostCellWidth()(2),
                                            Q_half_data[2]->getPointer(0),
                                            Q_half_data[2]->getPointer(1),
                                            Q_half_data[2]->getPointer(2)
#endif
        );

        break;
    case CUI:
        // Upwind side-centered densities onto faces.
        VC_NAVIER_STOKES_CUI_QUANTITY_FC(patch_lower(0),
                                         patch_upper(0),
                                         patch_lower(1),
                                         patch_upper(1),
#if (NDIM == 3)
                                         patch_lower(2),
                                         patch_upper(2),
#endif
                                         Q_data->getGhostCellWidth()(0),
                                         Q_data->getGhostCellWidth()(1),
#if (NDIM == 3)
                                         Q_data->getGhostCellWidth()(2),
#endif
                                         Q_data->getPointer(0),
                                         Q_data->getPointer(1),
#if (NDIM == 3)
                                         Q_data->getPointer(2),
#endif
                                         side_boxes[0].lower(0),
                                         side_boxes[0].upper(0),
                                         side_boxes[0].lower(1),
                                         side_boxes[0].upper(1),
#if (NDIM == 3)
                                         side_boxes[0].lower(2),
                                         side_boxes[0].upper(2),
#endif
                                         U_adv_data[0]->getGhostCellWidth()(0),
                                         U_adv_data[0]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                         U_adv_data[0]->getGhostCellWidth()(2),
#endif
                                         U_adv_data[0]->getPointer(0),
                                         U_adv_data[0]->getPointer(1),
#if (NDIM == 3)
                                         U_adv_data[0]->getPointer(2),
#endif
                                         Q_half_data[0]->getGhostCellWidth()(0),
                                         Q_half_data[0]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                         Q_half_data[0]->getGhostCellWidth()(2),
#endif
                                         Q_half_data[0]->getPointer(0),
                                         Q_half_data[0]->getPointer(1),
#if (NDIM == 3)
                                         Q_half_data[0]->getPointer(2),
#endif
                                         side_boxes[1].lower(0),
                                         side_boxes[1].upper(0),
                                         side_boxes[1].lower(1),
                                         side_boxes[1].upper(1),
#if (NDIM == 3)
                                         side_boxes[1].lower(2),
                                         side_boxes[1].upper(2),
#endif
                                         U_adv_data[1]->getGhostCellWidth()(0),
                                         U_adv_data[1]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                         U_adv_data[1]->getGhostCellWidth()(2),
#endif
                                         U_adv_data[1]->getPointer(0),
                                         U_adv_data[1]->getPointer(1),
#if (NDIM == 3)
                                         U_adv_data[1]->getPointer(2),
#endif
                                         Q_half_data[1]->getGhostCellWidth()(0),
                                         Q_half_data[1]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                         Q_half_data[1]->getGhostCellWidth()(2),
#endif
                                         Q_half_data[1]->getPointer(0),
                                         Q_half_data[1]->getPointer(1)
#if (NDIM == 3)
                                             ,
                                         Q_half_data[1]->getPointer(2),
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
                                         Q_half_data[2]->getGhostCellWidth()(0),
                                         Q_half_data[2]->getGhostCellWidth()(1),
                                         Q_half_data[2]->getGhostCellWidth()(2),
                                         Q_half_data[2]->getPointer(0),
                                         Q_half_data[2]->getPointer(1),
                                         Q_half_data[2]->getPointer(2)
#endif
        );

        break;
    case FBICS:
        VC_NAVIER_STOKES_FBICS_QUANTITY_FC(patch_lower(0),
                                           patch_upper(0),
                                           patch_lower(1),
                                           patch_upper(1),
#if (NDIM == 3)
                                           patch_lower(2),
                                           patch_upper(2),
#endif
                                           Q_data->getGhostCellWidth()(0),
                                           Q_data->getGhostCellWidth()(1),
#if (NDIM == 3)
                                           Q_data->getGhostCellWidth()(2),
#endif
                                           Q_data->getPointer(0),
                                           Q_data->getPointer(1),
#if (NDIM == 3)
                                           Q_data->getPointer(2),
#endif
                                           side_boxes[0].lower(0),
                                           side_boxes[0].upper(0),
                                           side_boxes[0].lower(1),
                                           side_boxes[0].upper(1),
#if (NDIM == 3)
                                           side_boxes[0].lower(2),
                                           side_boxes[0].upper(2),
#endif
                                           U_adv_data[0]->getGhostCellWidth()(0),
                                           U_adv_data[0]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                           U_adv_data[0]->getGhostCellWidth()(2),
#endif
                                           U_adv_data[0]->getPointer(0),
                                           U_adv_data[0]->getPointer(1),
#if (NDIM == 3)
                                           U_adv_data[0]->getPointer(2),
#endif
                                           Q_half_data[0]->getGhostCellWidth()(0),
                                           Q_half_data[0]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                           Q_half_data[0]->getGhostCellWidth()(2),
#endif
                                           Q_half_data[0]->getPointer(0),
                                           Q_half_data[0]->getPointer(1),
#if (NDIM == 3)
                                           Q_half_data[0]->getPointer(2),
#endif
                                           side_boxes[1].lower(0),
                                           side_boxes[1].upper(0),
                                           side_boxes[1].lower(1),
                                           side_boxes[1].upper(1),
#if (NDIM == 3)
                                           side_boxes[1].lower(2),
                                           side_boxes[1].upper(2),
#endif
                                           U_adv_data[1]->getGhostCellWidth()(0),
                                           U_adv_data[1]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                           U_adv_data[1]->getGhostCellWidth()(2),
#endif
                                           U_adv_data[1]->getPointer(0),
                                           U_adv_data[1]->getPointer(1),
#if (NDIM == 3)
                                           U_adv_data[1]->getPointer(2),
#endif
                                           Q_half_data[1]->getGhostCellWidth()(0),
                                           Q_half_data[1]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                           Q_half_data[1]->getGhostCellWidth()(2),
#endif
                                           Q_half_data[1]->getPointer(0),
                                           Q_half_data[1]->getPointer(1)
#if (NDIM == 3)
                                               ,
                                           Q_half_data[1]->getPointer(2),
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
                                           Q_half_data[2]->getGhostCellWidth()(0),
                                           Q_half_data[2]->getGhostCellWidth()(1),
                                           Q_half_data[2]->getGhostCellWidth()(2),
                                           Q_half_data[2]->getPointer(0),
                                           Q_half_data[2]->getPointer(1),
                                           Q_half_data[2]->getPointer(2)
#endif
        );
        break;
    case MGAMMA:
        VC_NAVIER_STOKES_MGAMMA_QUANTITY_FC(patch_lower(0),
                                            patch_upper(0),
                                            patch_lower(1),
                                            patch_upper(1),
#if (NDIM == 3)
                                            patch_lower(2),
                                            patch_upper(2),
#endif
                                            Q_data->getGhostCellWidth()(0),
                                            Q_data->getGhostCellWidth()(1),
#if (NDIM == 3)
                                            Q_data->getGhostCellWidth()(2),
#endif
                                            Q_data->getPointer(0),
                                            Q_data->getPointer(1),
#if (NDIM == 3)
                                            Q_data->getPointer(2),
#endif
                                            side_boxes[0].lower(0),
                                            side_boxes[0].upper(0),
                                            side_boxes[0].lower(1),
                                            side_boxes[0].upper(1),
#if (NDIM == 3)
                                            side_boxes[0].lower(2),
                                            side_boxes[0].upper(2),
#endif
                                            U_adv_data[0]->getGhostCellWidth()(0),
                                            U_adv_data[0]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                            U_adv_data[0]->getGhostCellWidth()(2),
#endif
                                            U_adv_data[0]->getPointer(0),
                                            U_adv_data[0]->getPointer(1),
#if (NDIM == 3)
                                            U_adv_data[0]->getPointer(2),
#endif
                                            Q_half_data[0]->getGhostCellWidth()(0),
                                            Q_half_data[0]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                            Q_half_data[0]->getGhostCellWidth()(2),
#endif
                                            Q_half_data[0]->getPointer(0),
                                            Q_half_data[0]->getPointer(1),
#if (NDIM == 3)
                                            Q_half_data[0]->getPointer(2),
#endif
                                            side_boxes[1].lower(0),
                                            side_boxes[1].upper(0),
                                            side_boxes[1].lower(1),
                                            side_boxes[1].upper(1),
#if (NDIM == 3)
                                            side_boxes[1].lower(2),
                                            side_boxes[1].upper(2),
#endif
                                            U_adv_data[1]->getGhostCellWidth()(0),
                                            U_adv_data[1]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                            U_adv_data[1]->getGhostCellWidth()(2),
#endif
                                            U_adv_data[1]->getPointer(0),
                                            U_adv_data[1]->getPointer(1),
#if (NDIM == 3)
                                            U_adv_data[1]->getPointer(2),
#endif
                                            Q_half_data[1]->getGhostCellWidth()(0),
                                            Q_half_data[1]->getGhostCellWidth()(1),
#if (NDIM == 3)
                                            Q_half_data[1]->getGhostCellWidth()(2),
#endif
                                            Q_half_data[1]->getPointer(0),
                                            Q_half_data[1]->getPointer(1)
#if (NDIM == 3)
                                                ,
                                            Q_half_data[1]->getPointer(2),
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
                                            Q_half_data[2]->getGhostCellWidth()(0),
                                            Q_half_data[2]->getGhostCellWidth()(1),
                                            Q_half_data[2]->getGhostCellWidth()(2),
                                            Q_half_data[2]->getPointer(0),
                                            Q_half_data[2]->getPointer(1),
                                            Q_half_data[2]->getPointer(2)
#endif
        );
        break;
    case PPM:
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            Pointer<SideData<NDIM, double> > dQ_data =
                new SideData<NDIM, double>(Q_data->getBox(), Q_data->getDepth(), Q_data->getGhostCellWidth());
            Pointer<SideData<NDIM, double> > Q_L_data =
                new SideData<NDIM, double>(Q_data->getBox(), Q_data->getDepth(), Q_data->getGhostCellWidth());
            Pointer<SideData<NDIM, double> > Q_R_data =
                new SideData<NDIM, double>(Q_data->getBox(), Q_data->getDepth(), Q_data->getGhostCellWidth());
            Pointer<SideData<NDIM, double> > Q_scratch1_data =
                new SideData<NDIM, double>(Q_data->getBox(), Q_data->getDepth(), Q_data->getGhostCellWidth());
#if (NDIM == 3)
            Pointer<SideData<NDIM, double> > Q_scratch2_data =
                new SideData<NDIM, double>(Q_data->getBox(), Q_data->getDepth(), Q_data->getGhostCellWidth());
#endif
#if (NDIM == 2)
            GODUNOV_EXTRAPOLATE_FC(side_boxes[axis].lower(0),
                                   side_boxes[axis].upper(0),
                                   side_boxes[axis].lower(1),
                                   side_boxes[axis].upper(1),
                                   Q_data->getGhostCellWidth()(0),
                                   Q_data->getGhostCellWidth()(1),
                                   Q_data->getPointer(axis),
                                   Q_scratch1_data->getPointer(axis),
                                   dQ_data->getPointer(axis),
                                   Q_L_data->getPointer(axis),
                                   Q_R_data->getPointer(axis),
                                   U_adv_data[axis]->getGhostCellWidth()(0),
                                   U_adv_data[axis]->getGhostCellWidth()(1),
                                   Q_half_data[axis]->getGhostCellWidth()(0),
                                   Q_half_data[axis]->getGhostCellWidth()(1),
                                   U_adv_data[axis]->getPointer(0),
                                   U_adv_data[axis]->getPointer(1),
                                   Q_half_data[axis]->getPointer(0),
                                   Q_half_data[axis]->getPointer(1));
#endif
#if (NDIM == 3)
            GODUNOV_EXTRAPOLATE_FC(side_boxes[axis].lower(0),
                                   side_boxes[axis].upper(0),
                                   side_boxes[axis].lower(1),
                                   side_boxes[axis].upper(1),
                                   side_boxes[axis].lower(2),
                                   side_boxes[axis].upper(2),
                                   Q_data->getGhostCellWidth()(0),
                                   Q_data->getGhostCellWidth()(1),
                                   Q_data->getGhostCellWidth()(2),
                                   Q_data->getPointer(axis),
                                   Q_scratch1_data->getPointer(axis),
                                   Q_scratch2_data->getPointer(axis),
                                   dQ_data->getPointer(axis),
                                   Q_L_data->getPointer(axis),
                                   Q_R_data->getPointer(axis),
                                   U_adv_data[axis]->getGhostCellWidth()(0),
                                   U_adv_data[axis]->getGhostCellWidth()(1),
                                   U_adv_data[axis]->getGhostCellWidth()(2),
                                   Q_half_data[axis]->getGhostCellWidth()(0),
                                   Q_half_data[axis]->getGhostCellWidth()(1),
                                   Q_half_data[axis]->getGhostCellWidth()(2),
                                   U_adv_data[axis]->getPointer(0),
                                   U_adv_data[axis]->getPointer(1),
                                   U_adv_data[axis]->getPointer(2),
                                   Q_half_data[axis]->getPointer(0),
                                   Q_half_data[axis]->getPointer(1),
                                   Q_half_data[axis]->getPointer(2));
#endif
        }
        break;
    default:
        TBOX_ERROR(
            "INSVCStaggeredConservativeMassMomentumRKIntegrator::"
            "interpolateSideQuantity():\n"
            << "  unsupported convective limiter: " << IBAMR::enum_to_string<LimiterType>(convective_limiter) << " \n"
            << "  valid choices are: UPWIND, CUI, FBICS, MGAMMA, PPM\n");
    }
} // interpolateSideQuantity

void
INSVCStaggeredConservativeMassMomentumRKIntegrator::computeConvectiveDerivative(
    Pointer<SideData<NDIM, double> > N_data,
    std::array<Pointer<FaceData<NDIM, double> >, NDIM> P_half_data,
    const std::array<Pointer<FaceData<NDIM, double> >, NDIM> U_adv_data,
    const std::array<Pointer<FaceData<NDIM, double> >, NDIM> R_half_data,
    const std::array<Pointer<FaceData<NDIM, double> >, NDIM> U_half_data,
    const std::array<Box<NDIM>, NDIM>& side_boxes,
    const double* const dx,
    const Pointer<Patch<NDIM> >& /*patch*/)
{
// Compute the upwinded momentum P_half = R_half * U_half
#if (NDIM == 2)
    VC_NAVIER_STOKES_COMPUTE_MOMENTUM_FC(side_boxes[0].lower(0),
                                         side_boxes[0].upper(0),
                                         side_boxes[0].lower(1),
                                         side_boxes[0].upper(1),
                                         P_half_data[0]->getGhostCellWidth()(0),
                                         P_half_data[0]->getGhostCellWidth()(1),
                                         P_half_data[0]->getPointer(0),
                                         P_half_data[0]->getPointer(1),
                                         R_half_data[0]->getGhostCellWidth()(0),
                                         R_half_data[0]->getGhostCellWidth()(1),
                                         R_half_data[0]->getPointer(0),
                                         R_half_data[0]->getPointer(1),
                                         U_half_data[0]->getGhostCellWidth()(0),
                                         U_half_data[0]->getGhostCellWidth()(1),
                                         U_half_data[0]->getPointer(0),
                                         U_half_data[0]->getPointer(1),
                                         side_boxes[1].lower(0),
                                         side_boxes[1].upper(0),
                                         side_boxes[1].lower(1),
                                         side_boxes[1].upper(1),
                                         P_half_data[1]->getGhostCellWidth()(0),
                                         P_half_data[1]->getGhostCellWidth()(1),
                                         P_half_data[1]->getPointer(0),
                                         P_half_data[1]->getPointer(1),
                                         R_half_data[1]->getGhostCellWidth()(0),
                                         R_half_data[1]->getGhostCellWidth()(1),
                                         R_half_data[1]->getPointer(0),
                                         R_half_data[1]->getPointer(1),
                                         U_half_data[1]->getGhostCellWidth()(0),
                                         U_half_data[1]->getGhostCellWidth()(1),
                                         U_half_data[1]->getPointer(0),
                                         U_half_data[1]->getPointer(1));
#endif
#if (NDIM == 3)
    VC_NAVIER_STOKES_COMPUTE_MOMENTUM_FC(side_boxes[0].lower(0),
                                         side_boxes[0].upper(0),
                                         side_boxes[0].lower(1),
                                         side_boxes[0].upper(1),
                                         side_boxes[0].lower(2),
                                         side_boxes[0].upper(2),
                                         P_half_data[0]->getGhostCellWidth()(0),
                                         P_half_data[0]->getGhostCellWidth()(1),
                                         P_half_data[0]->getGhostCellWidth()(2),
                                         P_half_data[0]->getPointer(0),
                                         P_half_data[0]->getPointer(1),
                                         P_half_data[0]->getPointer(2),
                                         R_half_data[0]->getGhostCellWidth()(0),
                                         R_half_data[0]->getGhostCellWidth()(1),
                                         R_half_data[0]->getGhostCellWidth()(2),
                                         R_half_data[0]->getPointer(0),
                                         R_half_data[0]->getPointer(1),
                                         R_half_data[0]->getPointer(2),
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
                                         P_half_data[1]->getGhostCellWidth()(0),
                                         P_half_data[1]->getGhostCellWidth()(1),
                                         P_half_data[1]->getGhostCellWidth()(2),
                                         P_half_data[1]->getPointer(0),
                                         P_half_data[1]->getPointer(1),
                                         P_half_data[1]->getPointer(2),
                                         R_half_data[1]->getGhostCellWidth()(0),
                                         R_half_data[1]->getGhostCellWidth()(1),
                                         R_half_data[1]->getGhostCellWidth()(2),
                                         R_half_data[1]->getPointer(0),
                                         R_half_data[1]->getPointer(1),
                                         R_half_data[1]->getPointer(2),
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
                                         P_half_data[2]->getGhostCellWidth()(0),
                                         P_half_data[2]->getGhostCellWidth()(1),
                                         P_half_data[2]->getGhostCellWidth()(2),
                                         P_half_data[2]->getPointer(0),
                                         P_half_data[2]->getPointer(1),
                                         P_half_data[2]->getPointer(2),
                                         R_half_data[2]->getGhostCellWidth()(0),
                                         R_half_data[2]->getGhostCellWidth()(1),
                                         R_half_data[2]->getGhostCellWidth()(2),
                                         R_half_data[2]->getPointer(0),
                                         R_half_data[2]->getPointer(1),
                                         R_half_data[2]->getPointer(2),
                                         U_half_data[2]->getGhostCellWidth()(0),
                                         U_half_data[2]->getGhostCellWidth()(1),
                                         U_half_data[2]->getGhostCellWidth()(2),
                                         U_half_data[2]->getPointer(0),
                                         U_half_data[2]->getPointer(1),
                                         U_half_data[2]->getPointer(2));
#endif

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
#if (NDIM == 2)
        CONVECT_DERIVATIVE_FC(dx,
                              side_boxes[axis].lower(0),
                              side_boxes[axis].upper(0),
                              side_boxes[axis].lower(1),
                              side_boxes[axis].upper(1),
                              U_adv_data[axis]->getGhostCellWidth()(0),
                              U_adv_data[axis]->getGhostCellWidth()(1),
                              P_half_data[axis]->getGhostCellWidth()(0),
                              P_half_data[axis]->getGhostCellWidth()(1),
                              U_adv_data[axis]->getPointer(0),
                              U_adv_data[axis]->getPointer(1),
                              P_half_data[axis]->getPointer(0),
                              P_half_data[axis]->getPointer(1),
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
                              P_half_data[axis]->getGhostCellWidth()(0),
                              P_half_data[axis]->getGhostCellWidth()(1),
                              P_half_data[axis]->getGhostCellWidth()(2),
                              U_adv_data[axis]->getPointer(0),
                              U_adv_data[axis]->getPointer(1),
                              U_adv_data[axis]->getPointer(2),
                              P_half_data[axis]->getPointer(0),
                              P_half_data[axis]->getPointer(1),
                              P_half_data[axis]->getPointer(2),
                              N_data->getGhostCellWidth()(0),
                              N_data->getGhostCellWidth()(1),
                              N_data->getGhostCellWidth()(2),
                              N_data->getPointer(axis));
#endif
    }

} // computeConvectiveDerivative

void
INSVCStaggeredConservativeMassMomentumRKIntegrator::computeDensityUpdate(
    Pointer<SideData<NDIM, double> > R_data,
    const double& a0,
    const Pointer<SideData<NDIM, double> > R0_data,
    const double& a1,
    const Pointer<SideData<NDIM, double> > R1_data,
    const double& a2,
    const std::array<Pointer<FaceData<NDIM, double> >, NDIM> U_adv_data,
    const std::array<Pointer<FaceData<NDIM, double> >, NDIM> R_half_data,
    const Pointer<SideData<NDIM, double> > S_data,
    const std::array<Box<NDIM>, NDIM>& side_boxes,
    const double& dt,
    const double* const dx)
{
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
#if (NDIM == 2)
        VC_UPDATE_DENSITY_FC(dx,
                             dt,
                             a0,
                             a1,
                             a2,
                             side_boxes[axis].lower(0),
                             side_boxes[axis].upper(0),
                             side_boxes[axis].lower(1),
                             side_boxes[axis].upper(1),
                             R0_data->getGhostCellWidth()(0),
                             R0_data->getGhostCellWidth()(1),
                             R0_data->getPointer(axis),
                             R1_data->getGhostCellWidth()(0),
                             R1_data->getGhostCellWidth()(1),
                             R1_data->getPointer(axis),
                             U_adv_data[axis]->getGhostCellWidth()(0),
                             U_adv_data[axis]->getGhostCellWidth()(1),
                             U_adv_data[axis]->getPointer(0),
                             U_adv_data[axis]->getPointer(1),
                             R_half_data[axis]->getGhostCellWidth()(0),
                             R_half_data[axis]->getGhostCellWidth()(1),
                             R_half_data[axis]->getPointer(0),
                             R_half_data[axis]->getPointer(1),
                             S_data->getGhostCellWidth()(0),
                             S_data->getGhostCellWidth()(1),
                             S_data->getPointer(axis),
                             R_data->getGhostCellWidth()(0),
                             R_data->getGhostCellWidth()(1),
                             R_data->getPointer(axis));
#endif
#if (NDIM == 3)
        VC_UPDATE_DENSITY_FC(dx,
                             dt,
                             a0,
                             a1,
                             a2,
                             side_boxes[axis].lower(0),
                             side_boxes[axis].upper(0),
                             side_boxes[axis].lower(1),
                             side_boxes[axis].upper(1),
                             side_boxes[axis].lower(2),
                             side_boxes[axis].upper(2),
                             R0_data->getGhostCellWidth()(0),
                             R0_data->getGhostCellWidth()(1),
                             R0_data->getGhostCellWidth()(2),
                             R0_data->getPointer(axis),
                             R1_data->getGhostCellWidth()(0),
                             R1_data->getGhostCellWidth()(1),
                             R1_data->getGhostCellWidth()(2),
                             R1_data->getPointer(axis),
                             U_adv_data[axis]->getGhostCellWidth()(0),
                             U_adv_data[axis]->getGhostCellWidth()(1),
                             U_adv_data[axis]->getGhostCellWidth()(2),
                             U_adv_data[axis]->getPointer(0),
                             U_adv_data[axis]->getPointer(1),
                             U_adv_data[axis]->getPointer(2),
                             R_half_data[axis]->getGhostCellWidth()(0),
                             R_half_data[axis]->getGhostCellWidth()(1),
                             R_half_data[axis]->getGhostCellWidth()(2),
                             R_half_data[axis]->getPointer(0),
                             R_half_data[axis]->getPointer(1),
                             R_half_data[axis]->getPointer(2),
                             S_data->getGhostCellWidth()(0),
                             S_data->getGhostCellWidth()(1),
                             S_data->getGhostCellWidth()(2),
                             S_data->getPointer(axis),
                             R_data->getGhostCellWidth()(0),
                             R_data->getGhostCellWidth()(1),
                             R_data->getGhostCellWidth()(2),
                             R_data->getPointer(axis));
#endif
    }
} // computeDensityUpdate

void
INSVCStaggeredConservativeMassMomentumRKIntegrator::computeErrorOfMassConservationEquation(
    Pointer<SideData<NDIM, double> > E_data,
    const Pointer<SideData<NDIM, double> > Rnew_data,
    const Pointer<SideData<NDIM, double> > Rold_data,
    const std::array<Pointer<FaceData<NDIM, double> >, NDIM> U_adv_data,
    const std::array<Pointer<FaceData<NDIM, double> >, NDIM> R_half_data,
    const std::array<Box<NDIM>, NDIM>& side_boxes,
    const double& dt,
    const double* const dx)
{
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
#if (NDIM == 2)
        VC_MASS_CONSERVATION_ERROR_FC(dx,
                                      dt,
                                      side_boxes[axis].lower(0),
                                      side_boxes[axis].upper(0),
                                      side_boxes[axis].lower(1),
                                      side_boxes[axis].upper(1),
                                      Rnew_data->getGhostCellWidth()(0),
                                      Rnew_data->getGhostCellWidth()(1),
                                      Rnew_data->getPointer(axis),
                                      Rold_data->getGhostCellWidth()(0),
                                      Rold_data->getGhostCellWidth()(1),
                                      Rold_data->getPointer(axis),
                                      U_adv_data[axis]->getGhostCellWidth()(0),
                                      U_adv_data[axis]->getGhostCellWidth()(1),
                                      U_adv_data[axis]->getPointer(0),
                                      U_adv_data[axis]->getPointer(1),
                                      R_half_data[axis]->getGhostCellWidth()(0),
                                      R_half_data[axis]->getGhostCellWidth()(1),
                                      R_half_data[axis]->getPointer(0),
                                      R_half_data[axis]->getPointer(1),
                                      E_data->getGhostCellWidth()(0),
                                      E_data->getGhostCellWidth()(1),
                                      E_data->getPointer(axis));
#endif
#if (NDIM == 3)
        VC_MASS_CONSERVATION_ERROR_FC(dx,
                                      dt,
                                      side_boxes[axis].lower(0),
                                      side_boxes[axis].upper(0),
                                      side_boxes[axis].lower(1),
                                      side_boxes[axis].upper(1),
                                      side_boxes[axis].lower(2),
                                      side_boxes[axis].upper(2),
                                      Rnew_data->getGhostCellWidth()(0),
                                      Rnew_data->getGhostCellWidth()(1),
                                      Rnew_data->getGhostCellWidth()(2),
                                      Rnew_data->getPointer(axis),
                                      Rold_data->getGhostCellWidth()(0),
                                      Rold_data->getGhostCellWidth()(1),
                                      Rold_data->getGhostCellWidth()(2),
                                      Rold_data->getPointer(axis),
                                      U_adv_data[axis]->getGhostCellWidth()(0),
                                      U_adv_data[axis]->getGhostCellWidth()(1),
                                      U_adv_data[axis]->getGhostCellWidth()(2),
                                      U_adv_data[axis]->getPointer(0),
                                      U_adv_data[axis]->getPointer(1),
                                      U_adv_data[axis]->getPointer(2),
                                      R_half_data[axis]->getGhostCellWidth()(0),
                                      R_half_data[axis]->getGhostCellWidth()(1),
                                      R_half_data[axis]->getGhostCellWidth()(2),
                                      R_half_data[axis]->getPointer(0),
                                      R_half_data[axis]->getPointer(1),
                                      R_half_data[axis]->getPointer(2),
                                      E_data->getGhostCellWidth()(0),
                                      E_data->getGhostCellWidth()(1),
                                      E_data->getGhostCellWidth()(2),
                                      E_data->getPointer(axis));
#endif
    }
} // computeErrorOfMassConservationEquation

void
INSVCStaggeredConservativeMassMomentumRKIntegrator::enforceDivergenceFreeConditionAtCoarseFineInterface(int U_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
#endif
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(U_idx);
            const int patch_ln = patch->getPatchLevelNumber();
            const int patch_num = patch->getPatchNumber();
            const Array<BoundaryBox<NDIM> >& cf_bdry_codim1_boxes = d_cf_boundary[patch_ln].getBoundaries(patch_num, 1);
            const int n_cf_bdry_codim1_boxes = cf_bdry_codim1_boxes.size();
            if (n_cf_bdry_codim1_boxes == 0) continue;

            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM> ghost_width_to_fill = CF_GHOST_WIDTH;
            for (int n = 0; n < n_cf_bdry_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = cf_bdry_codim1_boxes[n];
                const unsigned int location_index = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;
                const bool is_lower = location_index % 2 == 0;
                const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
                for (Box<NDIM>::Iterator b(bc_fill_box); b; b++)
                {
                    // Place i_s on the c-f interface.
                    hier::Index<NDIM> i_s = b();

                    // Work out from the coarse-fine interface to fill the ghost cell
                    // values so that the velocity field satisfies the discrete
                    // divergence-free condition.
                    for (int k = 0; k < u_data->getGhostCellWidth()(bdry_normal_axis);
                         ++k, i_s(bdry_normal_axis) += (is_lower ? -1 : +1))
                    {
                        // Determine the ghost cell value so that the divergence of
                        // the velocity field is zero in the ghost cell.
                        SideIndex<NDIM> i_g_s(
                            i_s, bdry_normal_axis, is_lower ? SideIndex<NDIM>::Lower : SideIndex<NDIM>::Upper);
                        (*u_data)(i_g_s) = 0.0;
                        double div_u_g = 0.0;
                        for (unsigned int axis = 0; axis < NDIM; ++axis)
                        {
                            const SideIndex<NDIM> i_g_s_upper(i_s, axis, SideIndex<NDIM>::Upper);
                            const SideIndex<NDIM> i_g_s_lower(i_s, axis, SideIndex<NDIM>::Lower);
                            div_u_g +=
                                ((*u_data)(i_g_s_upper) - (*u_data)(i_g_s_lower)) * dx[bdry_normal_axis] / dx[axis];
                        }
                        (*u_data)(i_g_s) = (is_lower ? +1.0 : -1.0) * div_u_g;
                    }
                }
            }
        }
    }
    return;
} // enforceDivergenceFreeConditionAtCoarseFineInterface

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
