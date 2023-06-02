// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#include "ibtk/CartSideRobinPhysBdryOp.h"
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "ComponentSelector.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "RobinBcCoefStrategy.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <array>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

// FORTRAN ROUTINES
#if (NDIM == 2)
#define CC_ROBIN_PHYS_BDRY_OP_1_X_FC IBTK_FC_FUNC(ccrobinphysbdryop1x2d, CCROBINPHYSBDRYOP1X2D)
#define CC_ROBIN_PHYS_BDRY_OP_1_Y_FC IBTK_FC_FUNC(ccrobinphysbdryop1y2d, CCROBINPHYSBDRYOP1Y2D)
#define CC_ROBIN_PHYS_BDRY_QUAD_OP_1_X_FC IBTK_FC_FUNC(ccrobinphysbdryquadop1x2d, CCROBINPHYSBDRYQUADOP1X2D)
#define CC_ROBIN_PHYS_BDRY_QUAD_OP_1_Y_FC IBTK_FC_FUNC(ccrobinphysbdryquadop1y2d, CCROBINPHYSBDRYQUADOP1Y2D)
#define SC_ROBIN_PHYS_BDRY_OP_1_X_FC IBTK_FC_FUNC(scrobinphysbdryop1x2d, SCROBINPHYSBDRYOP1X2D)
#define SC_ROBIN_PHYS_BDRY_OP_1_Y_FC IBTK_FC_FUNC(scrobinphysbdryop1y2d, SCROBINPHYSBDRYOP1Y2D)
#define SC_ROBIN_PHYS_BDRY_QUAD_OP_1_X_FC IBTK_FC_FUNC(scrobinphysbdryquadop1x2d, SCROBINPHYSBDRYQUADOP1X2D)
#define SC_ROBIN_PHYS_BDRY_QUAD_OP_1_Y_FC IBTK_FC_FUNC(scrobinphysbdryquadop1y2d, SCROBINPHYSBDRYQUADOP1Y2D)
#define SC_ROBIN_PHYS_BDRY_OP_2_FC IBTK_FC_FUNC(scrobinphysbdryop22d, SCROBINPHYSBDRYOP22D)
#endif // if (NDIM == 2)

#if (NDIM == 3)
#define CC_ROBIN_PHYS_BDRY_OP_1_X_FC IBTK_FC_FUNC(ccrobinphysbdryop1x3d, CCROBINPHYSBDRYOP1X3D)
#define CC_ROBIN_PHYS_BDRY_OP_1_Y_FC IBTK_FC_FUNC(ccrobinphysbdryop1y3d, CCROBINPHYSBDRYOP1Y3D)
#define CC_ROBIN_PHYS_BDRY_OP_1_Z_FC IBTK_FC_FUNC(ccrobinphysbdryop1z3d, CCROBINPHYSBDRYOP1Z3D)
#define CC_ROBIN_PHYS_BDRY_QUAD_OP_1_X_FC IBTK_FC_FUNC(ccrobinphysbdryquadop1x3d, CCROBINPHYSBDRYQUADOP1X3D)
#define CC_ROBIN_PHYS_BDRY_QUAD_OP_1_Y_FC IBTK_FC_FUNC(ccrobinphysbdryquadop1y3d, CCROBINPHYSBDRYQUADOP1Y3D)
#define CC_ROBIN_PHYS_BDRY_QUAD_OP_1_Z_FC IBTK_FC_FUNC(ccrobinphysbdryquadop1z3d, CCROBINPHYSBDRYQUADOP1Z3D)
#define SC_ROBIN_PHYS_BDRY_OP_1_X_FC IBTK_FC_FUNC(scrobinphysbdryop1x3d, SCROBINPHYSBDRYOP1X3D)
#define SC_ROBIN_PHYS_BDRY_OP_1_Y_FC IBTK_FC_FUNC(scrobinphysbdryop1y3d, SCROBINPHYSBDRYOP1Y3D)
#define SC_ROBIN_PHYS_BDRY_OP_1_Z_FC IBTK_FC_FUNC(scrobinphysbdryop1z3d, SCROBINPHYSBDRYOP1Z3D)
#define SC_ROBIN_PHYS_BDRY_QUAD_OP_1_X_FC IBTK_FC_FUNC(scrobinphysbdryquadop1x3d, SCROBINPHYSBDRYQUADOP1X3D)
#define SC_ROBIN_PHYS_BDRY_QUAD_OP_1_Y_FC IBTK_FC_FUNC(scrobinphysbdryquadop1y3d, SCROBINPHYSBDRYQUADOP1Y3D)
#define SC_ROBIN_PHYS_BDRY_QUAD_OP_1_Z_FC IBTK_FC_FUNC(scrobinphysbdryquadop1z3d, SCROBINPHYSBDRYQUADOP1Z3D)
#define SC_ROBIN_PHYS_BDRY_OP_2_FC IBTK_FC_FUNC(scrobinphysbdryop23d, SCROBINPHYSBDRYOP22D)
#define CC_ROBIN_PHYS_BDRY_OP_2_FC IBTK_FC_FUNC(ccrobinphysbdryop23d, CCROBINPHYSBDRYOP23D)
#define SC_ROBIN_PHYS_BDRY_OP_3_FC IBTK_FC_FUNC(scrobinphysbdryop33d, SCROBINPHYSBDRYOP32D)
#endif // if (NDIM == 3)

extern "C"
{
    void CC_ROBIN_PHYS_BDRY_OP_1_X_FC(double* U,
                                      const int& U_gcw,
                                      const double* acoef,
                                      const double* bcoef,
                                      const double* gcoef,
                                      const int& location_index,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
#if (NDIM == 3)
                                      const int& ilower2,
                                      const int& iupper2,
#endif
                                      const int& blower1,
                                      const int& bupper1,
#if (NDIM == 3)
                                      const int& blower2,
                                      const int& bupper2,
#endif
                                      const double* dx,
                                      const int& adjoint_op);

    void CC_ROBIN_PHYS_BDRY_QUAD_OP_1_X_FC(double* U,
                                           const int& U_gcw,
                                           const double* acoef,
                                           const double* bcoef,
                                           const double* gcoef,
                                           const int& location_index,
                                           const int& ilower0,
                                           const int& iupper0,
                                           const int& ilower1,
                                           const int& iupper1,
#if (NDIM == 3)
                                           const int& ilower2,
                                           const int& iupper2,
#endif
                                           const int& blower1,
                                           const int& bupper1,
#if (NDIM == 3)
                                           const int& blower2,
                                           const int& bupper2,
#endif
                                           const double* dx);

    void CC_ROBIN_PHYS_BDRY_OP_1_Y_FC(double* U,
                                      const int& U_gcw,
                                      const double* acoef,
                                      const double* bcoef,
                                      const double* gcoef,
                                      const int& location_index,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
#if (NDIM == 3)
                                      const int& ilower2,
                                      const int& iupper2,
#endif
                                      const int& blower0,
                                      const int& bupper0,
#if (NDIM == 3)
                                      const int& blower2,
                                      const int& bupper2,
#endif
                                      const double* dx,
                                      const int& adjoint_op);

    void CC_ROBIN_PHYS_BDRY_QUAD_OP_1_Y_FC(double* U,
                                           const int& U_gcw,
                                           const double* acoef,
                                           const double* bcoef,
                                           const double* gcoef,
                                           const int& location_index,
                                           const int& ilower0,
                                           const int& iupper0,
                                           const int& ilower1,
                                           const int& iupper1,
#if (NDIM == 3)
                                           const int& ilower2,
                                           const int& iupper2,
#endif
                                           const int& blower0,
                                           const int& bupper0,
#if (NDIM == 3)
                                           const int& blower2,
                                           const int& bupper2,
#endif
                                           const double* dx);

#if (NDIM == 3)
    void CC_ROBIN_PHYS_BDRY_OP_1_Z_FC(double* U,
                                      const int& U_gcw,
                                      const double* acoef,
                                      const double* bcoef,
                                      const double* gcoef,
                                      const int& location_index,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
                                      const int& ilower2,
                                      const int& iupper2,
                                      const int& blower0,
                                      const int& bupper0,
                                      const int& blower1,
                                      const int& bupper1,
                                      const double* dx,
                                      const int& adjoint_op);

    void CC_ROBIN_PHYS_BDRY_QUAD_OP_1_Z_FC(double* U,
                                           const int& U_gcw,
                                           const double* acoef,
                                           const double* bcoef,
                                           const double* gcoef,
                                           const int& location_index,
                                           const int& ilower0,
                                           const int& iupper0,
                                           const int& ilower1,
                                           const int& iupper1,
                                           const int& ilower2,
                                           const int& iupper2,
                                           const int& blower0,
                                           const int& bupper0,
                                           const int& blower1,
                                           const int& bupper1,
                                           const double* dx);
#endif

    void SC_ROBIN_PHYS_BDRY_OP_1_X_FC(double* u0,
                                      const int& u_gcw,
                                      const double* acoef,
                                      const double* bcoef,
                                      const double* gcoef,
                                      const int& location_index,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
#if (NDIM == 3)
                                      const int& ilower2,
                                      const int& iupper2,
#endif
                                      const int& blower1,
                                      const int& bupper1,
#if (NDIM == 3)
                                      const int& blower2,
                                      const int& bupper2,
#endif
                                      const double* dx,
                                      const int& adjoint_op);

    void SC_ROBIN_PHYS_BDRY_OP_1_Y_FC(double* u1,
                                      const int& u_gcw,
                                      const double* acoef,
                                      const double* bcoef,
                                      const double* gcoef,
                                      const int& location_index,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
#if (NDIM == 3)
                                      const int& ilower2,
                                      const int& iupper2,
#endif
                                      const int& blower0,
                                      const int& bupper0,
#if (NDIM == 3)
                                      const int& blower2,
                                      const int& bupper2,
#endif
                                      const double* dx,
                                      const int& adjoint_op);

#if (NDIM == 3)
    void SC_ROBIN_PHYS_BDRY_OP_1_Z_FC(double* u2,
                                      const int& u_gcw,
                                      const double* acoef,
                                      const double* bcoef,
                                      const double* gcoef,
                                      const int& location_index,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
                                      const int& ilower2,
                                      const int& iupper2,
                                      const int& blower0,
                                      const int& bupper0,
                                      const int& blower1,
                                      const int& bupper1,
                                      const double* dx,
                                      const int& adjoint_op);
#endif

    void SC_ROBIN_PHYS_BDRY_QUAD_OP_1_X_FC(double* u1,
                                           const int& u_gcw,
                                           const double* acoef,
                                           const double* bcoef,
                                           const double* gcoef,
                                           const int& location_index,
                                           const int& ilower0,
                                           const int& iupper0,
                                           const int& ilower1,
                                           const int& iupper1,
#if (NDIM == 3)
                                           const int& ilower2,
                                           const int& iupper2,
#endif
                                           const int& blower0,
                                           const int& bupper0,
#if (NDIM == 3)
                                           const int& blower2,
                                           const int& bupper2,
#endif
                                           const double* dx);

    void SC_ROBIN_PHYS_BDRY_QUAD_OP_1_Y_FC(double* u1,
                                           const int& u_gcw,
                                           const double* acoef,
                                           const double* bcoef,
                                           const double* gcoef,
                                           const int& location_index,
                                           const int& ilower0,
                                           const int& iupper0,
                                           const int& ilower1,
                                           const int& iupper1,
#if (NDIM == 3)
                                           const int& ilower2,
                                           const int& iupper2,
#endif
                                           const int& blower0,
                                           const int& bupper0,
#if (NDIM == 3)
                                           const int& blower2,
                                           const int& bupper2,
#endif
                                           const double* dx);

#if (NDIM == 3)
    void SC_ROBIN_PHYS_BDRY_QUAD_OP_1_Z_FC(double* u2,
                                           const int& u_gcw,
                                           const double* acoef,
                                           const double* bcoef,
                                           const double* gcoef,
                                           const int& location_index,
                                           const int& ilower0,
                                           const int& iupper0,
                                           const int& ilower1,
                                           const int& iupper1,
                                           const int& ilower2,
                                           const int& iupper2,
                                           const int& blower0,
                                           const int& bupper0,
                                           const int& blower1,
                                           const int& bupper1,
                                           const double* dx);
#endif

    void SC_ROBIN_PHYS_BDRY_OP_2_FC(double* u0,
                                    double* u1,
#if (NDIM == 3)
                                    double* u2,
#endif
                                    const int& u_gcw,
                                    const int& location_index,
                                    const int& ilower0,
                                    const int& iupper0,
                                    const int& ilower1,
                                    const int& iupper1,
#if (NDIM == 3)
                                    const int& ilower2,
                                    const int& iupper2,
#endif
                                    const int& blower0,
                                    const int& bupper0,
                                    const int& blower1,
                                    const int& bupper1,
#if (NDIM == 3)
                                    const int& blower2,
                                    const int& bupper2,
#endif
                                    const int& adjoint_op);

#if (NDIM == 3)
    void CC_ROBIN_PHYS_BDRY_OP_2_FC(double* U,
                                    const int& U_gcw,
                                    const int& location_index,
                                    const int& ilower0,
                                    const int& iupper0,
                                    const int& ilower1,
                                    const int& iupper1,
                                    const int& ilower2,
                                    const int& iupper2,
                                    const int& blower0,
                                    const int& bupper0,
                                    const int& blower1,
                                    const int& bupper1,
                                    const int& blower2,
                                    const int& bupper2,
                                    const int& adjoint_op);

    void SC_ROBIN_PHYS_BDRY_OP_3_FC(double* u0,
                                    double* u1,
                                    double* u2,
                                    const int& u_gcw,
                                    const int& location_index,
                                    const int& ilower0,
                                    const int& iupper0,
                                    const int& ilower1,
                                    const int& iupper1,
                                    const int& ilower2,
                                    const int& iupper2,
                                    const int& blower0,
                                    const int& bupper0,
                                    const int& blower1,
                                    const int& bupper1,
                                    const int& blower2,
                                    const int& bupper2,
                                    const int& adjoint_op);
#endif
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int REFINE_OP_STENCIL_WIDTH = 1;

inline Box<NDIM>
compute_tangential_extension(const Box<NDIM>& box, const int data_axis)
{
    Box<NDIM> extended_box = box;
    extended_box.upper()(data_axis) += 1;
    return extended_box;
} // compute_tangential_extension
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartSideRobinPhysBdryOp::CartSideRobinPhysBdryOp() : RobinPhysBdryPatchStrategy()
{
    // intentionally blank
    return;
} // CartSideRobinPhysBdryOp

CartSideRobinPhysBdryOp::CartSideRobinPhysBdryOp(const int patch_data_index,
                                                 const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                 const bool homogeneous_bc,
                                                 std::string type)
    : RobinPhysBdryPatchStrategy(), d_type(std::move(type))
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coefs.size() == NDIM);
    TBOX_ASSERT(d_type == "LINEAR" || d_type == "QUADRATIC");
#endif
    setPatchDataIndex(patch_data_index);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
} // CartSideRobinPhysBdryOp

CartSideRobinPhysBdryOp::CartSideRobinPhysBdryOp(const std::set<int>& patch_data_indices,
                                                 const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                 const bool homogeneous_bc,
                                                 std::string type)
    : RobinPhysBdryPatchStrategy(), d_type(std::move(type))
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coefs.size() == NDIM);
    TBOX_ASSERT(d_type == "LINEAR" || d_type == "QUADRATIC");
#endif
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
} // CartSideRobinPhysBdryOp

CartSideRobinPhysBdryOp::CartSideRobinPhysBdryOp(const ComponentSelector& patch_data_indices,
                                                 const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                 const bool homogeneous_bc,
                                                 std::string type)
    : RobinPhysBdryPatchStrategy(), d_type(std::move(type))
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coefs.size() == NDIM);
    TBOX_ASSERT(d_type == "LINEAR" || d_type == "QUADRATIC");
#endif
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
} // CartSideRobinPhysBdryOp

CartSideRobinPhysBdryOp::~CartSideRobinPhysBdryOp()
{
    // intentionally blank
    return;
} // ~CartSideRobinPhysBdryOp

void
CartSideRobinPhysBdryOp::setPhysicalBoundaryConditions(Patch<NDIM>& patch,
                                                       const double fill_time,
                                                       const IntVector<NDIM>& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector<NDIM>(0)) return;

    // Ensure the target patch data corresponds to a side centered variable and
    // that the proper number of boundary condition objects have been provided.
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        Pointer<SideData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
        if (!patch_data)
        {
            TBOX_ERROR("CartSideRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                       << "  patch data index " << patch_data_idx
                       << " does not correspond to a side-centered double precision variable." << std::endl);
        }
        if (NDIM * patch_data->getDepth() != static_cast<int>(d_bc_coefs.size()))
        {
            TBOX_ERROR("CartSideRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                       << "  data depth for patch data index " << patch_data_idx << " is " << patch_data->getDepth()
                       << "\n"
                       << "  but " << d_bc_coefs.size()
                       << " boundary condition coefficient objects were provided to the class "
                          "constructor."
                       << std::endl);
        }
    }

    // Set the boundary conditions along the co-dimension one boundary boxes,
    // then extrapolate those values to the co-dimension two and three boundary
    // boxes.
    static const bool adjoint_op = false;
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim1Normal(
            patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch, adjoint_op);
    }
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim1Transverse(
            patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch, adjoint_op);
    }
    const Array<BoundaryBox<NDIM> > physical_codim2_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(patch);
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim2(patch_data_idx, physical_codim2_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
#if (NDIM > 2)
    const Array<BoundaryBox<NDIM> > physical_codim3_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim3Boxes(patch);
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim3(patch_data_idx, physical_codim3_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
#endif
    return;
} // setPhysicalBoundaryConditions

IntVector<NDIM>
CartSideRobinPhysBdryOp::getRefineOpStencilWidth() const
{
    return REFINE_OP_STENCIL_WIDTH;
} // getRefineOpStencilWidth

void
CartSideRobinPhysBdryOp::accumulateFromPhysicalBoundaryData(const Patch<NDIM>& patch,
                                                            const double fill_time,
                                                            const IntVector<NDIM>& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector<NDIM>(0)) return;
    if (d_type == "QUADRATIC")
    {
        TBOX_ERROR("Accumulating from physical boundaries not currently supported for QUADRATIC interpolation.");
    }

    // Ensure the target patch data corresponds to a side centered variable and
    // that the proper number of boundary condition objects have been provided.
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        Pointer<SideData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
        if (!patch_data)
        {
            TBOX_ERROR("CartSideRobinPhysBdryOp::accumulateFromPhysicalBoundaryData():\n"
                       << "  patch data index " << patch_data_idx
                       << " does not correspond to a side-centered double precision variable." << std::endl);
        }
        if (NDIM * patch_data->getDepth() != static_cast<int>(d_bc_coefs.size()))
        {
            TBOX_ERROR("CartSideRobinPhysBdryOp::accumulateFromPhysicalBoundaryData():\n"
                       << "  data depth for patch data index " << patch_data_idx << " is " << patch_data->getDepth()
                       << "\n"
                       << "  but " << d_bc_coefs.size()
                       << " boundary condition coefficient objects were provided to the class "
                          "constructor."
                       << std::endl);
        }
    }

    // Set the boundary conditions along the co-dimension one boundary boxes,
    // then extrapolate those values to the co-dimension two and three boundary
    // boxes.
    static const bool adjoint_op = true;
#if (NDIM > 2)
    const Array<BoundaryBox<NDIM> > physical_codim3_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim3Boxes(patch);
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim3(patch_data_idx, physical_codim3_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
#endif
    const Array<BoundaryBox<NDIM> > physical_codim2_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(patch);
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim2(patch_data_idx, physical_codim2_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim1Transverse(
            patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch, adjoint_op);
    }
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim1Normal(
            patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch, adjoint_op);
    }
    return;
} // accumulateFromPhysicalBoundaryData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CartSideRobinPhysBdryOp::fillGhostCellValuesCodim1Normal(const int patch_data_idx,
                                                         const Array<BoundaryBox<NDIM> >& physical_codim1_boxes,
                                                         const double fill_time,
                                                         const IntVector<NDIM>& ghost_width_to_fill,
                                                         const Patch<NDIM>& patch,
                                                         const bool adjoint_op)
{
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    if (n_physical_codim1_boxes == 0) return;

    const Box<NDIM>& patch_box = patch.getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    Pointer<SideData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
    const int patch_data_depth = patch_data->getDepth();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;
    var_db->mapIndexToVariable(patch_data_idx, var);
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR(
            "CartSideRobinPhysBdryOp::fillGhostCellValuesCodim1Normal():\n"
            "  patch data for patch data index "
            << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
#endif
    const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

    // Set the boundary condition coefficients and then set the ghost cell
    // values.
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const BoundaryBox<NDIM> trimmed_bdry_box(
            bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), location_index);
        const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
        Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        for (int d = 0; d < patch_data_depth; ++d)
        {
            RobinBcCoefStrategy<NDIM>* bc_coef = d_bc_coefs[NDIM * d + bdry_normal_axis];
            auto const extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coef);
            if (extended_bc_coef)
            {
                extended_bc_coef->setTargetPatchDataIndex(patch_data_idx);
                extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
            }
            bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, var, patch, trimmed_bdry_box, fill_time);
            if (d_homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
            if (extended_bc_coef) extended_bc_coef->clearTargetPatchDataIndex();
            if (location_index == 0 || location_index == 1)
            {
                if (d_type == "LINEAR")
                {
                    SC_ROBIN_PHYS_BDRY_OP_1_X_FC(patch_data->getPointer(bdry_normal_axis, d),
                                                 patch_data_gcw,
                                                 acoef_data->getPointer(),
                                                 bcoef_data->getPointer(),
                                                 gcoef_data->getPointer(),
                                                 location_index,
                                                 patch_box.lower(0),
                                                 patch_box.upper(0),
                                                 patch_box.lower(1),
                                                 patch_box.upper(1),
#if (NDIM == 3)
                                                 patch_box.lower(2),
                                                 patch_box.upper(2),
#endif
                                                 bc_coef_box.lower(1),
                                                 bc_coef_box.upper(1),
#if (NDIM == 3)
                                                 bc_coef_box.lower(2),
                                                 bc_coef_box.upper(2),
#endif
                                                 dx,
                                                 adjoint_op ? 1 : 0);
                }
                else if (d_type == "QUADRATIC")
                {
                    SC_ROBIN_PHYS_BDRY_QUAD_OP_1_X_FC(patch_data->getPointer(bdry_normal_axis, d),
                                                      patch_data_gcw,
                                                      acoef_data->getPointer(),
                                                      bcoef_data->getPointer(),
                                                      gcoef_data->getPointer(),
                                                      location_index,
                                                      patch_box.lower(0),
                                                      patch_box.upper(0),
                                                      patch_box.lower(1),
                                                      patch_box.upper(1),
#if (NDIM == 3)
                                                      patch_box.lower(2),
                                                      patch_box.upper(2),
#endif
                                                      bc_coef_box.lower(1),
                                                      bc_coef_box.upper(1),
#if (NDIM == 3)
                                                      bc_coef_box.lower(2),
                                                      bc_coef_box.upper(2),
#endif
                                                      dx);
                }
                else
                {
                    TBOX_ERROR("Unsupported interpolation type.");
                }
            }
            else if (location_index == 2 || location_index == 3)
            {
                if (d_type == "LINEAR")
                {
                    SC_ROBIN_PHYS_BDRY_OP_1_Y_FC(patch_data->getPointer(bdry_normal_axis, d),
                                                 patch_data_gcw,
                                                 acoef_data->getPointer(),
                                                 bcoef_data->getPointer(),
                                                 gcoef_data->getPointer(),
                                                 location_index,
                                                 patch_box.lower(0),
                                                 patch_box.upper(0),
                                                 patch_box.lower(1),
                                                 patch_box.upper(1),
#if (NDIM == 3)
                                                 patch_box.lower(2),
                                                 patch_box.upper(2),
#endif
                                                 bc_coef_box.lower(0),
                                                 bc_coef_box.upper(0),
#if (NDIM == 3)
                                                 bc_coef_box.lower(2),
                                                 bc_coef_box.upper(2),
#endif
                                                 dx,
                                                 adjoint_op ? 1 : 0);
                }
                else if (d_type == "QUADRATIC")
                {
                    SC_ROBIN_PHYS_BDRY_QUAD_OP_1_Y_FC(patch_data->getPointer(bdry_normal_axis, d),
                                                      patch_data_gcw,
                                                      acoef_data->getPointer(),
                                                      bcoef_data->getPointer(),
                                                      gcoef_data->getPointer(),
                                                      location_index,
                                                      patch_box.lower(0),
                                                      patch_box.upper(0),
                                                      patch_box.lower(1),
                                                      patch_box.upper(1),
#if (NDIM == 3)
                                                      patch_box.lower(2),
                                                      patch_box.upper(2),
#endif
                                                      bc_coef_box.lower(0),
                                                      bc_coef_box.upper(0),
#if (NDIM == 3)
                                                      bc_coef_box.lower(2),
                                                      bc_coef_box.upper(2),
#endif
                                                      dx);
                }
                else
                {
                    TBOX_ERROR("Unsupported interpolation type.");
                }
            }
#if (NDIM == 3)
            else if (location_index == 4 || location_index == 5)
            {
                if (d_type == "LINEAR")
                {
                    SC_ROBIN_PHYS_BDRY_OP_1_Z_FC(patch_data->getPointer(bdry_normal_axis, d),
                                                 patch_data_gcw,
                                                 acoef_data->getPointer(),
                                                 bcoef_data->getPointer(),
                                                 gcoef_data->getPointer(),
                                                 location_index,
                                                 patch_box.lower(0),
                                                 patch_box.upper(0),
                                                 patch_box.lower(1),
                                                 patch_box.upper(1),
                                                 patch_box.lower(2),
                                                 patch_box.upper(2),
                                                 bc_coef_box.lower(0),
                                                 bc_coef_box.upper(0),
                                                 bc_coef_box.lower(1),
                                                 bc_coef_box.upper(1),
                                                 dx,
                                                 adjoint_op ? 1 : 0);
                }
                else if (d_type == "QUADRATIC")
                {
                    SC_ROBIN_PHYS_BDRY_QUAD_OP_1_Z_FC(patch_data->getPointer(bdry_normal_axis, d),
                                                      patch_data_gcw,
                                                      acoef_data->getPointer(),
                                                      bcoef_data->getPointer(),
                                                      gcoef_data->getPointer(),
                                                      location_index,
                                                      patch_box.lower(0),
                                                      patch_box.upper(0),
                                                      patch_box.lower(1),
                                                      patch_box.upper(1),
                                                      patch_box.lower(2),
                                                      patch_box.upper(2),
                                                      bc_coef_box.lower(0),
                                                      bc_coef_box.upper(0),
                                                      bc_coef_box.lower(1),
                                                      bc_coef_box.upper(1),
                                                      dx);
                }
                else
                {
                    TBOX_ERROR("Unsupported interpolation type.");
                }
            }
#endif
        }
    }
    return;
} // fillGhostCellValuesCodim1Normal

void
CartSideRobinPhysBdryOp::fillGhostCellValuesCodim1Transverse(const int patch_data_idx,
                                                             const Array<BoundaryBox<NDIM> >& physical_codim1_boxes,
                                                             const double fill_time,
                                                             const IntVector<NDIM>& ghost_width_to_fill,
                                                             const Patch<NDIM>& patch,
                                                             const bool adjoint_op)
{
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    if (n_physical_codim1_boxes == 0) return;

    const Box<NDIM>& patch_box = patch.getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const patch_x_lower = pgeom->getXLower();
    const double* const patch_x_upper = pgeom->getXUpper();
    Pointer<SideData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
    const int patch_data_depth = patch_data->getDepth();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;
    var_db->mapIndexToVariable(patch_data_idx, var);
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR(
            "CartSideRobinPhysBdryOp::fillGhostCellValuesCodim1Transverse():\n"
            "  patch data for patch data index "
            << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
#endif
    const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

    std::array<Box<NDIM>, NDIM> side_box;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        side_box[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
    }

    const IntVector<NDIM>& ratio_to_level_zero = pgeom->getRatio();
    Array<Array<bool> > touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        touches_regular_bdry[axis].resizeArray(2);
        touches_periodic_bdry[axis].resizeArray(2);
        for (int upperlower = 0; upperlower < 2; ++upperlower)
        {
            touches_regular_bdry[axis][upperlower] = pgeom->getTouchesRegularBoundary(axis, upperlower);
            touches_periodic_bdry[axis][upperlower] = pgeom->getTouchesPeriodicBoundary(axis, upperlower);
        }
    }

    // Set the boundary condition coefficients and then set the ghost cell
    // values.
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const BoundaryBox<NDIM> trimmed_bdry_box(
            bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), location_index);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            if (axis != bdry_normal_axis)
            {
                const Box<NDIM> bc_coef_box = compute_tangential_extension(
                    PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), axis);
                Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

                // Temporarily reset the patch geometry object associated with
                // the patch so that boundary conditions are set at the correct
                // spatial locations.
                std::array<double, NDIM> shifted_patch_x_lower, shifted_patch_x_upper;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    shifted_patch_x_lower[d] = patch_x_lower[d];
                    shifted_patch_x_upper[d] = patch_x_upper[d];
                }
                shifted_patch_x_lower[axis] -= 0.5 * dx[axis];
                shifted_patch_x_upper[axis] -= 0.5 * dx[axis];
                auto pgeom = patch.getPatchGeometry();
                Patch<NDIM>& temp_patch = const_cast<Patch<NDIM>&>(patch);
                temp_patch.setPatchGeometry(new CartesianPatchGeometry<NDIM>(ratio_to_level_zero,
                                                                             touches_regular_bdry,
                                                                             touches_periodic_bdry,
                                                                             dx,
                                                                             shifted_patch_x_lower.data(),
                                                                             shifted_patch_x_upper.data()));

                // Set the boundary condition coefficients.
                for (int d = 0; d < patch_data_depth; ++d)
                {
                    RobinBcCoefStrategy<NDIM>* bc_coef = d_bc_coefs[NDIM * d + axis];
                    auto const extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coef);
                    if (extended_bc_coef)
                    {
                        extended_bc_coef->setTargetPatchDataIndex(patch_data_idx);
                        extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
                    }
                    bc_coef->setBcCoefs(
                        acoef_data, bcoef_data, gcoef_data, var, temp_patch, trimmed_bdry_box, fill_time);
                    if (d_homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
                    if (extended_bc_coef) extended_bc_coef->clearTargetPatchDataIndex();

                    // Set the boundary values.
                    if (location_index == 0 || location_index == 1)
                    {
                        if (d_type == "LINEAR")
                        {
                            CC_ROBIN_PHYS_BDRY_OP_1_X_FC(patch_data->getPointer(axis, d),
                                                         patch_data_gcw,
                                                         acoef_data->getPointer(),
                                                         bcoef_data->getPointer(),
                                                         gcoef_data->getPointer(),
                                                         location_index,
                                                         side_box[axis].lower(0),
                                                         side_box[axis].upper(0),
                                                         side_box[axis].lower(1),
                                                         side_box[axis].upper(1),
#if (NDIM == 3)
                                                         side_box[axis].lower(2),
                                                         side_box[axis].upper(2),
#endif
                                                         bc_coef_box.lower(1),
                                                         bc_coef_box.upper(1),
#if (NDIM == 3)
                                                         bc_coef_box.lower(2),
                                                         bc_coef_box.upper(2),
#endif
                                                         dx,
                                                         adjoint_op ? 1 : 0);
                        }
                        else if (d_type == "QUADRATIC")
                        {
                            CC_ROBIN_PHYS_BDRY_QUAD_OP_1_X_FC(patch_data->getPointer(axis, d),
                                                              patch_data_gcw,
                                                              acoef_data->getPointer(),
                                                              bcoef_data->getPointer(),
                                                              gcoef_data->getPointer(),
                                                              location_index,
                                                              side_box[axis].lower(0),
                                                              side_box[axis].upper(0),
                                                              side_box[axis].lower(1),
                                                              side_box[axis].upper(1),
#if (NDIM == 3)
                                                              side_box[axis].lower(2),
                                                              side_box[axis].upper(2),
#endif
                                                              bc_coef_box.lower(1),
                                                              bc_coef_box.upper(1),
#if (NDIM == 3)
                                                              bc_coef_box.lower(2),
                                                              bc_coef_box.upper(2),
#endif
                                                              dx);
                        }
                        else
                        {
                            TBOX_ERROR("Unsupported interpolation type.");
                        }
                    }
                    else if (location_index == 2 || location_index == 3)
                    {
                        if (d_type == "LINEAR")
                        {
                            CC_ROBIN_PHYS_BDRY_OP_1_Y_FC(patch_data->getPointer(axis, d),
                                                         patch_data_gcw,
                                                         acoef_data->getPointer(),
                                                         bcoef_data->getPointer(),
                                                         gcoef_data->getPointer(),
                                                         location_index,
                                                         side_box[axis].lower(0),
                                                         side_box[axis].upper(0),
                                                         side_box[axis].lower(1),
                                                         side_box[axis].upper(1),
#if (NDIM == 3)
                                                         side_box[axis].lower(2),
                                                         side_box[axis].upper(2),
#endif
                                                         bc_coef_box.lower(0),
                                                         bc_coef_box.upper(0),
#if (NDIM == 3)
                                                         bc_coef_box.lower(2),
                                                         bc_coef_box.upper(2),
#endif
                                                         dx,
                                                         adjoint_op ? 1 : 0);
                        }
                        else if (d_type == "QUADRATIC")
                        {
                            CC_ROBIN_PHYS_BDRY_QUAD_OP_1_Y_FC(patch_data->getPointer(axis, d),
                                                              patch_data_gcw,
                                                              acoef_data->getPointer(),
                                                              bcoef_data->getPointer(),
                                                              gcoef_data->getPointer(),
                                                              location_index,
                                                              side_box[axis].lower(0),
                                                              side_box[axis].upper(0),
                                                              side_box[axis].lower(1),
                                                              side_box[axis].upper(1),
#if (NDIM == 3)
                                                              side_box[axis].lower(2),
                                                              side_box[axis].upper(2),
#endif
                                                              bc_coef_box.lower(0),
                                                              bc_coef_box.upper(0),
#if (NDIM == 3)
                                                              bc_coef_box.lower(2),
                                                              bc_coef_box.upper(2),
#endif
                                                              dx);
                        }
                        else
                        {
                            TBOX_ERROR("Unsupported interpolation type.");
                        }
                    }
#if (NDIM == 3)
                    else if (location_index == 4 || location_index == 5)
                    {
                        if (d_type == "LINEAR")
                        {
                            CC_ROBIN_PHYS_BDRY_OP_1_Z_FC(patch_data->getPointer(axis, d),
                                                         patch_data_gcw,
                                                         acoef_data->getPointer(),
                                                         bcoef_data->getPointer(),
                                                         gcoef_data->getPointer(),
                                                         location_index,
                                                         side_box[axis].lower(0),
                                                         side_box[axis].upper(0),
                                                         side_box[axis].lower(1),
                                                         side_box[axis].upper(1),
                                                         side_box[axis].lower(2),
                                                         side_box[axis].upper(2),
                                                         bc_coef_box.lower(0),
                                                         bc_coef_box.upper(0),
                                                         bc_coef_box.lower(1),
                                                         bc_coef_box.upper(1),
                                                         dx,
                                                         adjoint_op ? 1 : 0);
                        }
                        else if (d_type == "QUADRATIC")
                        {
                            CC_ROBIN_PHYS_BDRY_QUAD_OP_1_Z_FC(patch_data->getPointer(axis, d),
                                                              patch_data_gcw,
                                                              acoef_data->getPointer(),
                                                              bcoef_data->getPointer(),
                                                              gcoef_data->getPointer(),
                                                              location_index,
                                                              side_box[axis].lower(0),
                                                              side_box[axis].upper(0),
                                                              side_box[axis].lower(1),
                                                              side_box[axis].upper(1),
                                                              side_box[axis].lower(2),
                                                              side_box[axis].upper(2),
                                                              bc_coef_box.lower(0),
                                                              bc_coef_box.upper(0),
                                                              bc_coef_box.lower(1),
                                                              bc_coef_box.upper(1),
                                                              dx);
                        }
                        else
                        {
                            TBOX_ERROR("Unsupported interpolation type.");
                        }
                    }
#endif
                }

                // Restore the patch geometry object.
                temp_patch.setPatchGeometry(pgeom);
            }
        }
    }
    return;
} // fillGhostCellValuesCodim1Transverse

void
CartSideRobinPhysBdryOp::fillGhostCellValuesCodim2(const int patch_data_idx,
                                                   const Array<BoundaryBox<NDIM> >& physical_codim2_boxes,
                                                   const IntVector<NDIM>& ghost_width_to_fill,
                                                   const Patch<NDIM>& patch,
                                                   const bool adjoint_op)
{
    const int n_physical_codim2_boxes = physical_codim2_boxes.size();
    if (n_physical_codim2_boxes == 0) return;

    const Box<NDIM>& patch_box = patch.getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    Pointer<SideData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
    const int patch_data_depth = patch_data->getDepth();
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR(
            "CartSideRobinPhysBdryOp::fillGhostCellValuesCodim2():\n"
            "  patch data for patch data index "
            << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
#endif
    const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

#if (NDIM == 3)
    std::array<Box<NDIM>, NDIM> side_box;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        side_box[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
    }
#endif

    for (int n = 0; n < n_physical_codim2_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim2_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        for (int d = 0; d < patch_data_depth; ++d)
        {
            SC_ROBIN_PHYS_BDRY_OP_2_FC(patch_data->getPointer(0, d),
                                       patch_data->getPointer(1, d),
#if (NDIM == 3)
                                       patch_data->getPointer(2, d),
#endif
                                       patch_data_gcw,
                                       location_index,
                                       patch_box.lower(0),
                                       patch_box.upper(0),
                                       patch_box.lower(1),
                                       patch_box.upper(1),
#if (NDIM == 3)
                                       patch_box.lower(2),
                                       patch_box.upper(2),
#endif
                                       bc_fill_box.lower(0),
                                       bc_fill_box.upper(0),
                                       bc_fill_box.lower(1),
                                       bc_fill_box.upper(1),
#if (NDIM == 3)
                                       bc_fill_box.lower(2),
                                       bc_fill_box.upper(2),
#endif
                                       adjoint_op ? 1 : 0);
#if (NDIM == 3)
            if (location_index < 4)
            {
                static const unsigned int axis = 0;
                Box<NDIM> side_bc_fill_box = SideGeometry<NDIM>::toSideBox(bc_fill_box, axis);
                CC_ROBIN_PHYS_BDRY_OP_2_FC(patch_data->getPointer(axis, d),
                                           patch_data_gcw,
                                           location_index,
                                           side_box[axis].lower(0),
                                           side_box[axis].upper(0),
                                           side_box[axis].lower(1),
                                           side_box[axis].upper(1),
                                           side_box[axis].lower(2),
                                           side_box[axis].upper(2),
                                           side_bc_fill_box.lower(0),
                                           side_bc_fill_box.upper(0),
                                           side_bc_fill_box.lower(1),
                                           side_bc_fill_box.upper(1),
                                           side_bc_fill_box.lower(2),
                                           side_bc_fill_box.upper(2),
                                           adjoint_op ? 1 : 0);
            }
            else if (location_index >= 4 && location_index < 8)
            {
                static const unsigned int axis = 1;
                Box<NDIM> side_bc_fill_box = SideGeometry<NDIM>::toSideBox(bc_fill_box, axis);
                CC_ROBIN_PHYS_BDRY_OP_2_FC(patch_data->getPointer(axis, d),
                                           patch_data_gcw,
                                           location_index,
                                           side_box[axis].lower(0),
                                           side_box[axis].upper(0),
                                           side_box[axis].lower(1),
                                           side_box[axis].upper(1),
                                           side_box[axis].lower(2),
                                           side_box[axis].upper(2),
                                           side_bc_fill_box.lower(0),
                                           side_bc_fill_box.upper(0),
                                           side_bc_fill_box.lower(1),
                                           side_bc_fill_box.upper(1),
                                           side_bc_fill_box.lower(2),
                                           side_bc_fill_box.upper(2),
                                           adjoint_op ? 1 : 0);
            }
            else if (location_index >= 8 && location_index < 12)
            {
                static const unsigned int axis = 2;
                Box<NDIM> side_bc_fill_box = SideGeometry<NDIM>::toSideBox(bc_fill_box, axis);
                CC_ROBIN_PHYS_BDRY_OP_2_FC(patch_data->getPointer(axis, d),
                                           patch_data_gcw,
                                           location_index,
                                           side_box[axis].lower(0),
                                           side_box[axis].upper(0),
                                           side_box[axis].lower(1),
                                           side_box[axis].upper(1),
                                           side_box[axis].lower(2),
                                           side_box[axis].upper(2),
                                           side_bc_fill_box.lower(0),
                                           side_bc_fill_box.upper(0),
                                           side_bc_fill_box.lower(1),
                                           side_bc_fill_box.upper(1),
                                           side_bc_fill_box.lower(2),
                                           side_bc_fill_box.upper(2),
                                           adjoint_op ? 1 : 0);
            }
#endif
        }
    }
    return;
} // fillGhostCellValuesCodim2

#if (NDIM > 2)
void
CartSideRobinPhysBdryOp::fillGhostCellValuesCodim3(const int patch_data_idx,
                                                   const Array<BoundaryBox<NDIM> >& physical_codim3_boxes,
                                                   const IntVector<NDIM>& ghost_width_to_fill,
                                                   const Patch<NDIM>& patch,
                                                   const bool adjoint_op)
{
    const int n_physical_codim3_boxes = physical_codim3_boxes.size();
    if (n_physical_codim3_boxes == 0) return;

    const Box<NDIM>& patch_box = patch.getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    Pointer<SideData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
    const int patch_data_depth = patch_data->getDepth();
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR(
            "CartSideRobinPhysBdryOp::fillGhostCellValuesCodim3():\n"
            "  patch data for patch data index "
            << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
#endif
    const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

    for (int n = 0; n < n_physical_codim3_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim3_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        for (int d = 0; d < patch_data_depth; ++d)
        {
            SC_ROBIN_PHYS_BDRY_OP_3_FC(patch_data->getPointer(0, d),
                                       patch_data->getPointer(1, d),
                                       patch_data->getPointer(2, d),
                                       patch_data_gcw,
                                       location_index,
                                       patch_box.lower(0),
                                       patch_box.upper(0),
                                       patch_box.lower(1),
                                       patch_box.upper(1),
                                       patch_box.lower(2),
                                       patch_box.upper(2),
                                       bc_fill_box.lower(0),
                                       bc_fill_box.upper(0),
                                       bc_fill_box.lower(1),
                                       bc_fill_box.upper(1),
                                       bc_fill_box.lower(2),
                                       bc_fill_box.upper(2),
                                       adjoint_op ? 1 : 0);
        }
    }
    return;
} // fillGhostCellValuesCodim3
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
