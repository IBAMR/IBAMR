// Filename: CartCellRobinPhysBdryOp.cpp
// Created on 10 Feb 2007 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include <stddef.h>
#include <map>
#include <ostream>
#include <string>

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartCellRobinPhysBdryOp.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "IBTK_config.h"
#include "Index.h"
#include "Patch.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAI_config.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
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

extern "C"
{
    void
    CC_ROBIN_PHYS_BDRY_OP_1_X_FC(
        double* U, const int& U_gcw,
        const double* acoef, const double* bcoef, const double* gcoef,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
#if (NDIM == 3)
        const int& ilower2, const int& iupper2,
#endif
        const int& blower1, const int& bupper1,
#if (NDIM == 3)
        const int& blower2, const int& bupper2,
#endif
        const double* dx);

    void
    CC_ROBIN_PHYS_BDRY_OP_1_Y_FC(
        double* U, const int& U_gcw,
        const double* acoef, const double* bcoef, const double* gcoef,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
#if (NDIM == 3)
        const int& ilower2, const int& iupper2,
#endif
        const int& blower0, const int& bupper0,
#if (NDIM == 3)
        const int& blower2, const int& bupper2,
#endif
        const double* dx);

#if (NDIM == 3)
    void
    CC_ROBIN_PHYS_BDRY_OP_1_Z_FC(
        double* U, const int& U_gcw,
        const double* acoef, const double* bcoef, const double* gcoef,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& ilower2, const int& iupper2,
        const int& blower0, const int& bupper0,
        const int& blower1, const int& bupper1,
        const double* dx);
#endif

    void
    CC_ROBIN_PHYS_BDRY_OP_2_FC(
        double* U, const int& U_gcw,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
#if (NDIM == 3)
        const int& ilower2, const int& iupper2,
#endif
        const int& blower0, const int& bupper0,
        const int& blower1, const int& bupper1
#if (NDIM == 3)
        ,const int& blower2,const int& bupper2
#endif
                               );

#if (NDIM == 3)
    void
    CC_ROBIN_PHYS_BDRY_OP_3_FC(
        double* U, const int& U_gcw,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& ilower2, const int& iupper2,
        const int& blower0, const int& bupper0,
        const int& blower1, const int& bupper1,
        const int& blower2, const int& bupper2);
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

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp()
    : RobinPhysBdryPatchStrategy()
{
    // intentionally blank
    return;
}// CartCellRobinPhysBdryOp

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(
    const int patch_data_index,
    RobinBcCoefStrategy<NDIM>* const bc_coef,
    const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndex(patch_data_index);
    setPhysicalBcCoef(bc_coef);
    setHomogeneousBc(homogeneous_bc);
    return;
}// CartCellRobinPhysBdryOp

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(
    const std::set<int>& patch_data_indices,
    RobinBcCoefStrategy<NDIM>* const bc_coef,
    const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoef(bc_coef);
    setHomogeneousBc(homogeneous_bc);
    return;
}// CartCellRobinPhysBdryOp

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(
    const ComponentSelector& patch_data_indices,
    RobinBcCoefStrategy<NDIM>* const bc_coef,
    const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoef(bc_coef);
    setHomogeneousBc(homogeneous_bc);
    return;
}// CartCellRobinPhysBdryOp

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(
    const int patch_data_index,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndex(patch_data_index);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// CartCellRobinPhysBdryOp

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(
    const std::set<int>& patch_data_indices,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// CartCellRobinPhysBdryOp

CartCellRobinPhysBdryOp::CartCellRobinPhysBdryOp(
    const ComponentSelector& patch_data_indices,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// CartCellRobinPhysBdryOp

CartCellRobinPhysBdryOp::~CartCellRobinPhysBdryOp()
{
    // intentionally blank
    return;
}// ~CartCellRobinPhysBdryOp

void
CartCellRobinPhysBdryOp::setPhysicalBoundaryConditions(
    Patch<NDIM>& patch,
    const double fill_time,
    const IntVector<NDIM>& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector<NDIM>(0)) return;

    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const Box<NDIM>& patch_box = patch.getBox();
    const Index<NDIM>& patch_lower = patch_box.lower();
    const Index<NDIM>& patch_upper = patch_box.upper();
    const double* const dx = pgeom->getDx();

    // Compute the boundary fill boxes.
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
#if (NDIM > 1)
    const Array<BoundaryBox<NDIM> > physical_codim2_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(patch);
    const int n_physical_codim2_boxes = physical_codim2_boxes.size();
#if (NDIM > 2)
    const Array<BoundaryBox<NDIM> > physical_codim3_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim3Boxes(patch);
    const int n_physical_codim3_boxes = physical_codim3_boxes.size();
#endif
#endif

    // Set the physical boundary condition coefficients for the specified
    // scratch patch data indices before actually filling the ghost cell values.
    std::map<int,std::vector<std::vector<Pointer<ArrayData<NDIM,double> > > > > acoefs, bcoefs, gcoefs;
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<Variable<NDIM> > var;
        var_db->mapIndexToVariable(patch_data_idx, var);
        Pointer<CellVariable<NDIM,double> > cc_var = var;
        if (!cc_var)
        {
            TBOX_ERROR("CartCellRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                       << "  patch data index " << patch_data_idx << " does not correspond to a cell-centered double precision variable." << std::endl);
        }

        Pointer<CellData<NDIM,double> > patch_data = patch.getPatchData(patch_data_idx);
        const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);
        if (patch_data->getDepth() != static_cast<int>(d_bc_coefs.size()))
        {
            TBOX_ERROR("CartCellRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                       << "  data depth for patch data index " << patch_data_idx << " is " << patch_data->getDepth() << "\n"
                       << "  but " << d_bc_coefs.size() << " boundary condition coefficient objects were provided to the class constructor." << std::endl);
        }

        acoefs[patch_data_idx].resize(patch_data->getDepth(),std::vector<Pointer<ArrayData<NDIM,double> > >(n_physical_codim1_boxes,Pointer<ArrayData<NDIM,double> >(NULL)));
        bcoefs[patch_data_idx].resize(patch_data->getDepth(),std::vector<Pointer<ArrayData<NDIM,double> > >(n_physical_codim1_boxes,Pointer<ArrayData<NDIM,double> >(NULL)));
        gcoefs[patch_data_idx].resize(patch_data->getDepth(),std::vector<Pointer<ArrayData<NDIM,double> > >(n_physical_codim1_boxes,Pointer<ArrayData<NDIM,double> >(NULL)));
        for (int depth = 0; depth < patch_data->getDepth(); ++depth)
        {
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
                const BoundaryBox<NDIM> trimmed_bdry_box(bdry_box.getBox()*bc_fill_box, bdry_box.getBoundaryType(), bdry_box.getLocationIndex());
                const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
                Pointer<ArrayData<NDIM,double> > acoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM,double> > bcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM,double> > gcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                ExtendedRobinBcCoefStrategy* const extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[depth]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->setTargetPatchDataIndex(patch_data_idx);
                    extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
                }
                d_bc_coefs[depth]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, var, patch, trimmed_bdry_box, fill_time);
                if (d_homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
                if (extended_bc_coef) extended_bc_coef->clearTargetPatchDataIndex();
                acoefs[patch_data_idx][depth][n] = acoef_data;
                bcoefs[patch_data_idx][depth][n] = bcoef_data;
                gcoefs[patch_data_idx][depth][n] = gcoef_data;
            }
        }
    }

    // Set the boundary conditions along the co-dimension one boundary boxes,
    // then extrapolate those values to the co-dimension two and three boundary
    // boxes.
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        Pointer<CellData<NDIM,double> > patch_data = patch.getPatchData(patch_data_idx);
        const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);
        for (int depth = 0; depth < patch_data->getDepth(); ++depth)
        {
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                Pointer<ArrayData<NDIM,double> > acoef_data = acoefs[patch_data_idx][depth][n];
                Pointer<ArrayData<NDIM,double> > bcoef_data = bcoefs[patch_data_idx][depth][n];
                Pointer<ArrayData<NDIM,double> > gcoef_data = gcoefs[patch_data_idx][depth][n];

                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);

                double* const U = patch_data->getPointer(depth);
                const int U_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
                if (U_gcw != (patch_data->getGhostCellWidth()).min())
                {
                    TBOX_ERROR("CartCellRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                               "  patch data for patch data index " << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
                }
#endif
                const unsigned int location_index = bdry_box.getLocationIndex();
                switch (location_index)
                {
                    case 0:  // lower x
                    case 1:  // upper x
                        CC_ROBIN_PHYS_BDRY_OP_1_X_FC(
                            U, U_gcw,
                            acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                            location_index,
                            patch_lower(0), patch_upper(0),
                            patch_lower(1), patch_upper(1),
#if (NDIM == 3)
                            patch_lower(2), patch_upper(2),
#endif
                            bc_fill_box.lower()(1), bc_fill_box.upper()(1),
#if (NDIM == 3)
                            bc_fill_box.lower()(2), bc_fill_box.upper()(2),
#endif
                            dx);
                        break;
                    case 2:  // lower y
                    case 3:  // upper y
                        CC_ROBIN_PHYS_BDRY_OP_1_Y_FC(
                            U, U_gcw,
                            acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                            location_index,
                            patch_lower(0), patch_upper(0),
                            patch_lower(1), patch_upper(1),
#if (NDIM == 3)
                            patch_lower(2), patch_upper(2),
#endif
                            bc_fill_box.lower()(0), bc_fill_box.upper()(0),
#if (NDIM == 3)
                            bc_fill_box.lower()(2), bc_fill_box.upper()(2),
#endif
                            dx);
                        break;
#if (NDIM == 3)
                    case 4:  // lower z
                    case 5:  // upper z
                        CC_ROBIN_PHYS_BDRY_OP_1_Z_FC(
                            U, U_gcw,
                            acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                            location_index,
                            patch_lower(0), patch_upper(0),
                            patch_lower(1), patch_upper(1),
                            patch_lower(2), patch_upper(2),
                            bc_fill_box.lower()(0), bc_fill_box.upper()(0),
                            bc_fill_box.lower()(1), bc_fill_box.upper()(1),
                            dx);
                        break;
#endif
                    default:
                        TBOX_ASSERT(false);
                }
            }

#if (NDIM > 1)
            for (int n = 0; n < n_physical_codim2_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim2_boxes[n];
                const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);

                double* const U = patch_data->getPointer(depth);
                const int U_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
                if (U_gcw != (patch_data->getGhostCellWidth()).min())
                {
                    TBOX_ERROR("CartCellRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                               "  patch data for patch data index " << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
                }
#endif
                const unsigned int location_index = bdry_box.getLocationIndex();
                CC_ROBIN_PHYS_BDRY_OP_2_FC(
                    U, U_gcw,
                    location_index,
                    patch_lower(0), patch_upper(0),
                    patch_lower(1), patch_upper(1),
#if (NDIM == 3)
                    patch_lower(2), patch_upper(2),
#endif
                    bc_fill_box.lower()(0), bc_fill_box.upper()(0),
                    bc_fill_box.lower()(1), bc_fill_box.upper()(1)
#if (NDIM == 3)
                    ,bc_fill_box.lower()(2), bc_fill_box.upper()(2)
#endif
                                           );
            }
#endif

#if (NDIM > 2)
            for (int n = 0; n < n_physical_codim3_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim3_boxes[n];
                const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);

                double* const U = patch_data->getPointer(depth);
                const int U_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
                if (U_gcw != (patch_data->getGhostCellWidth()).min())
                {
                    TBOX_ERROR("CartCellRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                               "  patch data for patch data index " << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
                }
#endif
                const unsigned int location_index = bdry_box.getLocationIndex();
                CC_ROBIN_PHYS_BDRY_OP_3_FC(
                    U, U_gcw,
                    location_index,
                    patch_lower(0), patch_upper(0),
                    patch_lower(1), patch_upper(1),
                    patch_lower(2), patch_upper(2),
                    bc_fill_box.lower()(0), bc_fill_box.upper()(0),
                    bc_fill_box.lower()(1), bc_fill_box.upper()(1),
                    bc_fill_box.lower()(2), bc_fill_box.upper()(2));
            }
#endif
        }
    }
    return;
}// setPhysicalBoundaryConditions

IntVector<NDIM>
CartCellRobinPhysBdryOp::getRefineOpStencilWidth() const
{
    return REFINE_OP_STENCIL_WIDTH;
}// getRefineOpStencilWidth

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
