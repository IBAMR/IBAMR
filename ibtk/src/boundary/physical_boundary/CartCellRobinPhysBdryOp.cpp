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

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "IBTK_config.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/solv/RobinBcCoefStrategy.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "SAMRAI/tbox/Array.h"

#include "SAMRAI/tbox/Utilities.h"

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
}

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(const int patch_data_index,
                                                 boost::shared_ptr<RobinBcCoefStrategy> const bc_coef,
                                                 const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndex(patch_data_index);
    setPhysicalBcCoef(bc_coef);
    setHomogeneousBc(homogeneous_bc);
    return;
}

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(const std::set<int>& patch_data_indices,
                                                 boost::shared_ptr<RobinBcCoefStrategy> const bc_coef,
                                                 const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoef(bc_coef);
    setHomogeneousBc(homogeneous_bc);
    return;
}

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(const ComponentSelector& patch_data_indices,
                                                 boost::shared_ptr<RobinBcCoefStrategy> const bc_coef,
                                                 const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoef(bc_coef);
    setHomogeneousBc(homogeneous_bc);
    return;
}

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(const int patch_data_index,
                                                 const std::vector<boost::shared_ptr<RobinBcCoefStrategy>>& bc_coefs,
                                                 const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndex(patch_data_index);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(const std::set<int>& patch_data_indices,
                                                 const std::vector<boost::shared_ptr<RobinBcCoefStrategy>>& bc_coefs,
                                                 const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(const ComponentSelector& patch_data_indices,
                                                 const std::vector<boost::shared_ptr<RobinBcCoefStrategy>>& bc_coefs,
                                                 const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}

CartCellRobinPhysBdryOp::~CartCellRobinPhysBdryOp()
{
    // intentionally blank
    return;
}

void CartCellRobinPhysBdryOp::setPhysicalBoundaryConditions(Patch& patch,
                                                            const double fill_time,
                                                            const IntVector& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector::getZero(DIM)) return;

    // Ensure the target patch data corresponds to a cell centered variable and
    // that the proper number of boundary condition objects have been provided.
    for (auto cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        auto patch_data = BOOST_CAST<CellData<double> >(patch.getPatchData(patch_data_idx));
        if (!patch_data)
        {
            TBOX_ERROR("CartCellRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                       << "  patch data index " << patch_data_idx
                       << " does not correspond to a cell-centered double precision variable." << std::endl);
        }
        if (patch_data->getDepth() != static_cast<int>(d_bc_coefs.size()))
        {
            TBOX_ERROR("CartCellRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                       << "  data depth for patch data index " << patch_data_idx << " is " << patch_data->getDepth()
                       << "\n"
                       << "  but " << d_bc_coefs.size()
                       << " boundary condition coefficient objects were provided to the class "
                          "constructor." << std::endl);
        }
    }

    // Set the boundary conditions along the co-dimension one boundary boxes,
    // then extrapolate those values to the co-dimension two and three boundary
    // boxes.
    static const bool adjoint_op = false;
    const std::vector<BoundaryBox> physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    for (auto cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim1(patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch,
                                  adjoint_op);
    }
    const std::vector<BoundaryBox> physical_codim2_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(patch);
    for (auto cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim2(patch_data_idx, physical_codim2_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
#if (NDIM > 2)
    const std::vector<BoundaryBox> physical_codim3_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim3Boxes(patch);
    for (auto cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim3(patch_data_idx, physical_codim3_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
#endif
    return;
}

IntVector CartCellRobinPhysBdryOp::getRefineOpStencilWidth(const Dimension& dim) const
{
    return IntVector(dim, REFINE_OP_STENCIL_WIDTH);
}

void CartCellRobinPhysBdryOp::accumulateFromPhysicalBoundaryData(Patch& patch,
                                                                 const double fill_time,
                                                                 const IntVector& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector::getZero(DIM)) return;

    // Ensure the target patch data corresponds to a cell centered variable and
    // that the proper number of boundary condition objects have been provided.
    for (auto cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        auto patch_data = BOOST_CAST<CellData<double> >(patch.getPatchData(patch_data_idx));
        if (!patch_data)
        {
            TBOX_ERROR("CartCellRobinPhysBdryOp::accumulateFromPhysicalBoundaryData():\n"
                       << "  patch data index " << patch_data_idx
                       << " does not correspond to a cell-centered double precision variable." << std::endl);
        }
        if (patch_data->getDepth() != static_cast<int>(d_bc_coefs.size()))
        {
            TBOX_ERROR("CartCellRobinPhysBdryOp::accumulateFromPhysicalBoundaryData():\n"
                       << "  data depth for patch data index " << patch_data_idx << " is " << patch_data->getDepth()
                       << "\n"
                       << "  but " << d_bc_coefs.size()
                       << " boundary condition coefficient objects were provided to the class "
                          "constructor." << std::endl);
        }
    }

    // Set the boundary conditions along the co-dimension one boundary boxes,
    // then extrapolate those values to the co-dimension two and three boundary
    // boxes.
    static const bool adjoint_op = true;
#if (NDIM > 2)
    const std::vector<BoundaryBox> physical_codim3_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim3Boxes(patch);
    for (auto cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim3(patch_data_idx, physical_codim3_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
#endif
    const std::vector<BoundaryBox> physical_codim2_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(patch);
    for (auto cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim2(patch_data_idx, physical_codim2_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
    const std::vector<BoundaryBox> physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    for (auto cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        fillGhostCellValuesCodim1(patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch,
                                  adjoint_op);
    }
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void CartCellRobinPhysBdryOp::fillGhostCellValuesCodim1(const int patch_data_idx,
                                                        const std::vector<BoundaryBox>& physical_codim1_boxes,
                                                        const double fill_time,
                                                        const IntVector& ghost_width_to_fill,
                                                        Patch& patch,
                                                        const bool adjoint_op)
{
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    if (n_physical_codim1_boxes == 0) return;

    const Box& patch_box = patch.getBox();
    auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch.getPatchGeometry());
    const double* const dx = pgeom->getDx();
    auto patch_data = BOOST_CAST<CellData<double> >(patch.getPatchData(patch_data_idx));
    const int patch_data_depth = patch_data->getDepth();
    auto var_db = VariableDatabase::getDatabase();
    boost::shared_ptr<Variable> var;
    var_db->mapIndexToVariable(patch_data_idx, var);
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR(
            "CartCellRobinPhysBdryOp::fillGhostCellValuesCodim1():\n"
            "  patch data for patch data index "
            << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
    const IntVector gcw_to_fill = IntVector::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

    // Set the boundary condition coefficients and then set the ghost cell
    // values.
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox& bdry_box = physical_codim1_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const BoundaryBox trimmed_bdry_box(bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(),
                                           bdry_box.getLocationIndex());
        const Box bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
        auto acoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
        ;
        auto bcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
        ;
        auto gcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
        ;
        for (int d = 0; d < patch_data_depth; ++d)
        {
            boost::shared_ptr<RobinBcCoefStrategy> bc_coef = d_bc_coefs[d];
            auto extended_bc_coef = boost::dynamic_pointer_cast<ExtendedRobinBcCoefStrategy>(bc_coef);
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
                CC_ROBIN_PHYS_BDRY_OP_1_X_FC(patch_data->getPointer(d), patch_data_gcw, acoef_data->getPointer(),
                                             bcoef_data->getPointer(), gcoef_data->getPointer(), location_index,
                                             patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                                             patch_box.upper(1),
#if (NDIM == 3)
                                             patch_box.lower(2), patch_box.upper(2),
#endif
                                             bc_fill_box.lower(1), bc_fill_box.upper(1),
#if (NDIM == 3)
                                             bc_fill_box.lower(2), bc_fill_box.upper(2),
#endif
                                             dx, adjoint_op ? 1 : 0);
                break;
            case 2: // lower y
            case 3: // upper y
                CC_ROBIN_PHYS_BDRY_OP_1_Y_FC(patch_data->getPointer(d), patch_data_gcw, acoef_data->getPointer(),
                                             bcoef_data->getPointer(), gcoef_data->getPointer(), location_index,
                                             patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                                             patch_box.upper(1),
#if (NDIM == 3)
                                             patch_box.lower(2), patch_box.upper(2),
#endif
                                             bc_fill_box.lower(0), bc_fill_box.upper(0),
#if (NDIM == 3)
                                             bc_fill_box.lower(2), bc_fill_box.upper(2),
#endif
                                             dx, adjoint_op ? 1 : 0);
                break;
#if (NDIM == 3)
            case 4: // lower z
            case 5: // upper z
                CC_ROBIN_PHYS_BDRY_OP_1_Z_FC(patch_data->getPointer(d), patch_data_gcw, acoef_data->getPointer(),
                                             bcoef_data->getPointer(), gcoef_data->getPointer(), location_index,
                                             patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                                             patch_box.upper(1), patch_box.lower(2), patch_box.upper(2),
                                             bc_fill_box.lower(0), bc_fill_box.upper(0), bc_fill_box.lower(1),
                                             bc_fill_box.upper(1), dx, adjoint_op ? 1 : 0);
                break;
#endif
            default:
                TBOX_ASSERT(false);
            }
        }
    }
    return;
}

void CartCellRobinPhysBdryOp::fillGhostCellValuesCodim2(const int patch_data_idx,
                                                        const std::vector<BoundaryBox>& physical_codim2_boxes,
                                                        const IntVector& ghost_width_to_fill,
                                                        const Patch& patch,
                                                        const bool adjoint_op)
{
    const int n_physical_codim2_boxes = physical_codim2_boxes.size();
    if (n_physical_codim2_boxes == 0) return;

    const Box& patch_box = patch.getBox();
    auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch.getPatchGeometry());
    auto patch_data = BOOST_CAST<CellData<double> >(patch.getPatchData(patch_data_idx));
    const int patch_data_depth = patch_data->getDepth();
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR(
            "CartCellRobinPhysBdryOp::fillGhostCellValuesCodim2():\n"
            "  patch data for patch data index "
            << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
    const IntVector gcw_to_fill = IntVector::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

    for (int n = 0; n < n_physical_codim2_boxes; ++n)
    {
        const BoundaryBox& bdry_box = physical_codim2_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        for (int d = 0; d < patch_data_depth; ++d)
        {
            CC_ROBIN_PHYS_BDRY_OP_2_FC(patch_data->getPointer(d), patch_data_gcw, location_index, patch_box.lower(0),
                                       patch_box.upper(0), patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                                       patch_box.lower(2), patch_box.upper(2),
#endif
                                       bc_fill_box.lower(0), bc_fill_box.upper(0), bc_fill_box.lower(1),
                                       bc_fill_box.upper(1),
#if (NDIM == 3)
                                       bc_fill_box.lower(2), bc_fill_box.upper(2),
#endif
                                       adjoint_op ? 1 : 0);
        }
    }
    return;
}

#if (NDIM > 2)
void CartCellRobinPhysBdryOp::fillGhostCellValuesCodim3(const int patch_data_idx,
                                                        const std::vector<BoundaryBox>& physical_codim3_boxes,
                                                        const IntVector& ghost_width_to_fill,
                                                        const Patch& patch,
                                                        const bool adjoint_op)
{
    const int n_physical_codim3_boxes = physical_codim3_boxes.size();
    if (n_physical_codim3_boxes == 0) return;

    const Box& patch_box = patch.getBox();
    auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch.getPatchGeometry());
    auto patch_data = BOOST_CAST<CellData<double> >(patch.getPatchData(patch_data_idx));
    const int patch_data_depth = patch_data->getDepth();
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR(
            "CartCellRobinPhysBdryOp::fillGhostCellValuesCodim3():\n"
            "  patch data for patch data index "
            << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
    const IntVector gcw_to_fill = IntVector::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

    for (int n = 0; n < n_physical_codim3_boxes; ++n)
    {
        const BoundaryBox& bdry_box = physical_codim3_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        for (int d = 0; d < patch_data_depth; ++d)
        {
            CC_ROBIN_PHYS_BDRY_OP_3_FC(patch_data->getPointer(d), patch_data_gcw, location_index, patch_box.lower(0),
                                       patch_box.upper(0), patch_box.lower(1), patch_box.upper(1), patch_box.lower(2),
                                       patch_box.upper(2), bc_fill_box.lower(0), bc_fill_box.upper(0),
                                       bc_fill_box.lower(1), bc_fill_box.upper(1), bc_fill_box.lower(2),
                                       bc_fill_box.upper(2), adjoint_op ? 1 : 0);
        }
    }
    return;
}
#endif

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
