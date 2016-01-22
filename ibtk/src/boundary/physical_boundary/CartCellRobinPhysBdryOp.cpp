// Filename: CartCellRobinPhysBdryOp.cpp
// Created on 10 Feb 2007 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

#include <ostream>
#include <set>
#include <vector>

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "ComponentSelector.h"
#include "IBTK_config.h"
#include "IntVector.h"
#include "Patch.h"
#include "RobinBcCoefStrategy.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Array.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define CC_ROBIN_PHYS_BDRY_OP_1_X_FC IBTK_FC_FUNC(ccrobinphysbdryop1x2d, CCROBINPHYSBDRYOP1X2D)
#define CC_ROBIN_PHYS_BDRY_OP_1_Y_FC IBTK_FC_FUNC(ccrobinphysbdryop1y2d, CCROBINPHYSBDRYOP1Y2D)
#define CC_ROBIN_PHYS_BDRY_OP_2_FC IBTK_FC_FUNC(ccrobinphysbdryop22d, CCROBINPHYSBDRYOP22D)
#endif // if (NDIM == 2)

#if (NDIM == 3)
#define CC_ROBIN_PHYS_BDRY_OP_1_X_FC IBTK_FC_FUNC(ccrobinphysbdryop1x3d, CCROBINPHYSBDRYOP1X3D)
#define CC_ROBIN_PHYS_BDRY_OP_1_Y_FC IBTK_FC_FUNC(ccrobinphysbdryop1y3d, CCROBINPHYSBDRYOP1Y3D)
#define CC_ROBIN_PHYS_BDRY_OP_1_Z_FC IBTK_FC_FUNC(ccrobinphysbdryop1z3d, CCROBINPHYSBDRYOP1Z3D)
#define CC_ROBIN_PHYS_BDRY_OP_2_FC IBTK_FC_FUNC(ccrobinphysbdryop23d, CCROBINPHYSBDRYOP23D)
#define CC_ROBIN_PHYS_BDRY_OP_3_FC IBTK_FC_FUNC(ccrobinphysbdryop33d, CCROBINPHYSBDRYOP33D)
#endif // if (NDIM == 3)

extern "C" {
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
#endif

void CC_ROBIN_PHYS_BDRY_OP_2_FC(double* U,
                                const int& U_gcw,
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
void CC_ROBIN_PHYS_BDRY_OP_3_FC(double* U,
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
#endif
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int REFINE_OP_STENCIL_WIDTH = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp() : RobinPhysBdryPatchStrategy()
{
    // intentionally blank
    return;
} // CartCellRobinPhysBdryOp

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(const int patch_data_index,
                                                 RobinBcCoefStrategy<NDIM>* const bc_coef,
                                                 const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndex(patch_data_index);
    setPhysicalBcCoef(bc_coef);
    setHomogeneousBc(homogeneous_bc);
    return;
} // CartCellRobinPhysBdryOp

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(const std::set<int>& patch_data_indices,
                                                 RobinBcCoefStrategy<NDIM>* const bc_coef,
                                                 const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoef(bc_coef);
    setHomogeneousBc(homogeneous_bc);
    return;
} // CartCellRobinPhysBdryOp

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(const ComponentSelector& patch_data_indices,
                                                 RobinBcCoefStrategy<NDIM>* const bc_coef,
                                                 const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoef(bc_coef);
    setHomogeneousBc(homogeneous_bc);
    return;
} // CartCellRobinPhysBdryOp

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(const int patch_data_index,
                                                 const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                 const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndex(patch_data_index);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
} // CartCellRobinPhysBdryOp

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(const std::set<int>& patch_data_indices,
                                                 const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                 const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
} // CartCellRobinPhysBdryOp

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(const ComponentSelector& patch_data_indices,
                                                 const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                 const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
} // CartCellRobinPhysBdryOp

CartCellRobinPhysBdryOp::~CartCellRobinPhysBdryOp()
{
    // intentionally blank
    return;
} // ~CartCellRobinPhysBdryOp

void
CartCellRobinPhysBdryOp::setPhysicalBoundaryConditions(Patch<NDIM>& patch,
                                                       const double fill_time,
                                                       const IntVector<NDIM>& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector<NDIM>(0)) return;

    // Ensure the target patch data corresponds to a cell centered variable and
    // that the proper number of boundary condition objects have been provided.
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        Pointer<CellData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
        if (!patch_data)
        {
            TBOX_ERROR("CartCellRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                       << "  patch data index "
                       << patch_data_idx
                       << " does not correspond to a cell-centered double precision variable."
                       << std::endl);
        }
        if (patch_data->getDepth() != static_cast<int>(d_bc_coefs.size()))
        {
            TBOX_ERROR("CartCellRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                       << "  data depth for patch data index "
                       << patch_data_idx
                       << " is "
                       << patch_data->getDepth()
                       << "\n"
                       << "  but "
                       << d_bc_coefs.size()
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
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim1(
            patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch, adjoint_op);
    }
    const Array<BoundaryBox<NDIM> > physical_codim2_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(patch);
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim2(patch_data_idx, physical_codim2_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
#if (NDIM > 2)
    const Array<BoundaryBox<NDIM> > physical_codim3_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim3Boxes(patch);
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim3(patch_data_idx, physical_codim3_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
#endif
    return;
} // setPhysicalBoundaryConditions

IntVector<NDIM>
CartCellRobinPhysBdryOp::getRefineOpStencilWidth() const
{
    return REFINE_OP_STENCIL_WIDTH;
} // getRefineOpStencilWidth

void
CartCellRobinPhysBdryOp::accumulateFromPhysicalBoundaryData(Patch<NDIM>& patch,
                                                            const double fill_time,
                                                            const IntVector<NDIM>& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector<NDIM>(0)) return;

    // Ensure the target patch data corresponds to a cell centered variable and
    // that the proper number of boundary condition objects have been provided.
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        Pointer<CellData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
        if (!patch_data)
        {
            TBOX_ERROR("CartCellRobinPhysBdryOp::accumulateFromPhysicalBoundaryData():\n"
                       << "  patch data index "
                       << patch_data_idx
                       << " does not correspond to a cell-centered double precision variable."
                       << std::endl);
        }
        if (patch_data->getDepth() != static_cast<int>(d_bc_coefs.size()))
        {
            TBOX_ERROR("CartCellRobinPhysBdryOp::accumulateFromPhysicalBoundaryData():\n"
                       << "  data depth for patch data index "
                       << patch_data_idx
                       << " is "
                       << patch_data->getDepth()
                       << "\n"
                       << "  but "
                       << d_bc_coefs.size()
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
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim3(patch_data_idx, physical_codim3_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
#endif
    const Array<BoundaryBox<NDIM> > physical_codim2_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(patch);
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim2(patch_data_idx, physical_codim2_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim1(
            patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch, adjoint_op);
    }
    return;
} // accumulateFromPhysicalBoundaryData

/////////////////////////////// PROTECTED ////////////////////////////////////

void
CartCellRobinPhysBdryOp::fillGhostCellValuesCodim1(const int patch_data_idx,
                                                   const Array<BoundaryBox<NDIM> >& physical_codim1_boxes,
                                                   const double fill_time,
                                                   const IntVector<NDIM>& ghost_width_to_fill,
                                                   Patch<NDIM>& patch,
                                                   const bool adjoint_op)
{
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    if (n_physical_codim1_boxes == 0) return;

    const Box<NDIM>& patch_box = patch.getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    Pointer<CellData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
    const int patch_data_depth = patch_data->getDepth();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;
    var_db->mapIndexToVariable(patch_data_idx, var);
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR(
            "CartCellRobinPhysBdryOp::fillGhostCellValuesCodim1():\n"
            "  patch data for patch data index "
            << patch_data_idx
            << " does not have uniform ghost cell widths."
            << std::endl);
    }
#endif
    const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

    // Set the boundary condition coefficients and then set the ghost cell
    // values.
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const BoundaryBox<NDIM> trimmed_bdry_box(
            bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), bdry_box.getLocationIndex());
        const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
        Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        for (int d = 0; d < patch_data_depth; ++d)
        {
            RobinBcCoefStrategy<NDIM>* bc_coef = d_bc_coefs[d];
            ExtendedRobinBcCoefStrategy* const extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coef);
            if (extended_bc_coef)
            {
                extended_bc_coef->setTargetPatchDataIndex(patch_data_idx);
                extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
            }
            bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, var, patch, trimmed_bdry_box, fill_time);
            if (d_homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
            if (extended_bc_coef) extended_bc_coef->clearTargetPatchDataIndex();
            switch (location_index)
            {
            case 0: // lower x
            case 1: // upper x
                CC_ROBIN_PHYS_BDRY_OP_1_X_FC(patch_data->getPointer(d),
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
                                             bc_fill_box.lower(1),
                                             bc_fill_box.upper(1),
#if (NDIM == 3)
                                             bc_fill_box.lower(2),
                                             bc_fill_box.upper(2),
#endif
                                             dx,
                                             adjoint_op ? 1 : 0);
                break;
            case 2: // lower y
            case 3: // upper y
                CC_ROBIN_PHYS_BDRY_OP_1_Y_FC(patch_data->getPointer(d),
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
                                             bc_fill_box.lower(0),
                                             bc_fill_box.upper(0),
#if (NDIM == 3)
                                             bc_fill_box.lower(2),
                                             bc_fill_box.upper(2),
#endif
                                             dx,
                                             adjoint_op ? 1 : 0);
                break;
#if (NDIM == 3)
            case 4: // lower z
            case 5: // upper z
                CC_ROBIN_PHYS_BDRY_OP_1_Z_FC(patch_data->getPointer(d),
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
                                             bc_fill_box.lower(0),
                                             bc_fill_box.upper(0),
                                             bc_fill_box.lower(1),
                                             bc_fill_box.upper(1),
                                             dx,
                                             adjoint_op ? 1 : 0);
                break;
#endif
            default:
                TBOX_ASSERT(false);
            }
        }
    }
    return;
} // fillGhostCellValuesCodim1

void
CartCellRobinPhysBdryOp::fillGhostCellValuesCodim2(const int patch_data_idx,
                                                   const Array<BoundaryBox<NDIM> >& physical_codim2_boxes,
                                                   const IntVector<NDIM>& ghost_width_to_fill,
                                                   const Patch<NDIM>& patch,
                                                   const bool adjoint_op)
{
    const int n_physical_codim2_boxes = physical_codim2_boxes.size();
    if (n_physical_codim2_boxes == 0) return;

    const Box<NDIM>& patch_box = patch.getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    Pointer<CellData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
    const int patch_data_depth = patch_data->getDepth();
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR(
            "CartCellRobinPhysBdryOp::fillGhostCellValuesCodim2():\n"
            "  patch data for patch data index "
            << patch_data_idx
            << " does not have uniform ghost cell widths."
            << std::endl);
    }
#endif
    const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

    for (int n = 0; n < n_physical_codim2_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim2_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        for (int d = 0; d < patch_data_depth; ++d)
        {
            CC_ROBIN_PHYS_BDRY_OP_2_FC(patch_data->getPointer(d),
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
        }
    }
    return;
} // fillGhostCellValuesCodim2

#if (NDIM > 2)
void
CartCellRobinPhysBdryOp::fillGhostCellValuesCodim3(const int patch_data_idx,
                                                   const Array<BoundaryBox<NDIM> >& physical_codim3_boxes,
                                                   const IntVector<NDIM>& ghost_width_to_fill,
                                                   const Patch<NDIM>& patch,
                                                   const bool adjoint_op)
{
    const int n_physical_codim3_boxes = physical_codim3_boxes.size();
    if (n_physical_codim3_boxes == 0) return;

    const Box<NDIM>& patch_box = patch.getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    Pointer<CellData<NDIM, double> > patch_data = patch.getPatchData(patch_data_idx);
    const int patch_data_depth = patch_data->getDepth();
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR(
            "CartCellRobinPhysBdryOp::fillGhostCellValuesCodim3():\n"
            "  patch data for patch data index "
            << patch_data_idx
            << " does not have uniform ghost cell widths."
            << std::endl);
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
            CC_ROBIN_PHYS_BDRY_OP_3_FC(patch_data->getPointer(d),
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

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
