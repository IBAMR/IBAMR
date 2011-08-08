// Filename: CartSideRobinPhysBdryOp.C
// Created on 21 May 2008 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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

#include "CartSideRobinPhysBdryOp.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <ArrayDataBasicOps.h>
#include <CartesianPatchGeometry.h>
#include <SideData.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableDatabase.h>
#include <tbox/Array.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

// BLITZ++ INCLUDES
#include <blitz/tinyvec.h>

// C++ STDLIB INCLUDES
#include <map>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define CC_ROBIN_PHYS_BDRY_OP_1_X_FC FC_FUNC(ccrobinphysbdryop1x2d, CCROBINPHYSBDRYOP1X2D)
#define CC_ROBIN_PHYS_BDRY_OP_1_Y_FC FC_FUNC(ccrobinphysbdryop1y2d, CCROBINPHYSBDRYOP1Y2D)
#define SC_ROBIN_PHYS_BDRY_OP_1_X_FC FC_FUNC(scrobinphysbdryop1x2d, SCROBINPHYSBDRYOP1X2D)
#define SC_ROBIN_PHYS_BDRY_OP_1_Y_FC FC_FUNC(scrobinphysbdryop1y2d, SCROBINPHYSBDRYOP1Y2D)
#define SC_ROBIN_PHYS_BDRY_OP_2_FC FC_FUNC(scrobinphysbdryop22d, SCROBINPHYSBDRYOP22D)
#endif // if (NDIM == 2)

#if (NDIM == 3)
#define CC_ROBIN_PHYS_BDRY_OP_1_X_FC FC_FUNC(ccrobinphysbdryop1x3d, CCROBINPHYSBDRYOP1X3D)
#define CC_ROBIN_PHYS_BDRY_OP_1_Y_FC FC_FUNC(ccrobinphysbdryop1y3d, CCROBINPHYSBDRYOP1Y3D)
#define CC_ROBIN_PHYS_BDRY_OP_1_Z_FC FC_FUNC(ccrobinphysbdryop1z3d, CCROBINPHYSBDRYOP1Z3D)
#define SC_ROBIN_PHYS_BDRY_OP_1_X_FC FC_FUNC(scrobinphysbdryop1x3d, SCROBINPHYSBDRYOP1X3D)
#define SC_ROBIN_PHYS_BDRY_OP_1_Y_FC FC_FUNC(scrobinphysbdryop1y3d, SCROBINPHYSBDRYOP1Y3D)
#define SC_ROBIN_PHYS_BDRY_OP_1_Z_FC FC_FUNC(scrobinphysbdryop1z3d, SCROBINPHYSBDRYOP1Z3D)
#define SC_ROBIN_PHYS_BDRY_OP_2_FC FC_FUNC(scrobinphysbdryop23d, SCROBINPHYSBDRYOP22D)
#define CC_ROBIN_PHYS_BDRY_OP_2_FC FC_FUNC(ccrobinphysbdryop23d, CCROBINPHYSBDRYOP23D)
#define SC_ROBIN_PHYS_BDRY_OP_3_FC FC_FUNC(scrobinphysbdryop33d, SCROBINPHYSBDRYOP32D)
#endif // if (NDIM == 3)

extern "C"
{
    void
    CC_ROBIN_PHYS_BDRY_OP_1_X_FC(
        double* u, const int& u_gcw,
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
        double* u, const int& u_gcw,
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
    SC_ROBIN_PHYS_BDRY_OP_1_X_FC(
        double* u0, const int& u_gcw,
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
    SC_ROBIN_PHYS_BDRY_OP_1_Y_FC(
        double* u1, const int& u_gcw,
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
    SC_ROBIN_PHYS_BDRY_OP_1_Z_FC(
        double* u2, const int& u_gcw,
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
    SC_ROBIN_PHYS_BDRY_OP_2_FC(
        double* u0, double* u1,
#if (NDIM == 3)
        double* u2,
#endif
        const int& u_gcw,
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
    CC_ROBIN_PHYS_BDRY_OP_2_FC(
        double* U, const int& U_gcw,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& ilower2, const int& iupper2,
        const int& blower0, const int& bupper0,
        const int& blower1, const int& bupper1,
        const int& blower2, const int& bupper2);

    void
    SC_ROBIN_PHYS_BDRY_OP_3_FC(
        double* u0, double* u1, double* u2,
        const int& u_gcw,
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
static const int REFINE_OP_STENCIL_WIDTH = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

inline Box<NDIM>
compute_tangential_extension(
    const Box<NDIM>& box,
    const int data_axis)
{
    Box<NDIM> extended_box = box;
    extended_box.upper()(data_axis) += 1;
    return extended_box;
}// compute_tangential_extension
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartSideRobinPhysBdryOp::CartSideRobinPhysBdryOp()
    : d_patch_data_indices(),
      d_bc_coefs(),
      d_homogeneous_bc(false)
{
    // intentionally blank
    return;
}// CartSideRobinPhysBdryOp

CartSideRobinPhysBdryOp::CartSideRobinPhysBdryOp(
    const int patch_data_index,
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs,
    const bool homogeneous_bc)
    : d_patch_data_indices(),
      d_bc_coefs(),
      d_homogeneous_bc(false)
{
    setPatchDataIndex(patch_data_index);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// CartSideRobinPhysBdryOp

CartSideRobinPhysBdryOp::CartSideRobinPhysBdryOp(
    const std::set<int>& patch_data_indices,
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs,
    const bool homogeneous_bc)
    : d_patch_data_indices(),
      d_bc_coefs(),
      d_homogeneous_bc(false)
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// CartSideRobinPhysBdryOp

CartSideRobinPhysBdryOp::CartSideRobinPhysBdryOp(
    const ComponentSelector& patch_data_indices,
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs,
    const bool homogeneous_bc)
    : d_patch_data_indices(),
      d_bc_coefs(),
      d_homogeneous_bc(false)
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// CartSideRobinPhysBdryOp

CartSideRobinPhysBdryOp::~CartSideRobinPhysBdryOp()
{
    // intentionally blank
    return;
}// ~CartSideRobinPhysBdryOp

void
CartSideRobinPhysBdryOp::setPatchDataIndex(
    const int patch_data_index)
{
    std::set<int> patch_data_indices;
    patch_data_indices.insert(patch_data_index);
    setPatchDataIndices(patch_data_indices);
    return;
}// setPatchDataIndex

void
CartSideRobinPhysBdryOp::setPatchDataIndices(
    const std::set<int>& patch_data_indices)
{
    d_patch_data_indices.clear();
    d_patch_data_indices = patch_data_indices;
    return;
}// setPatchDataIndices

void
CartSideRobinPhysBdryOp::setPatchDataIndices(
    const ComponentSelector& patch_data_indices)
{
    std::set<int> patch_data_index_set;
    for (int l = 0; l < patch_data_indices.getSize(); ++l)
    {
        if (patch_data_indices.isSet(l))
        {
            const int patch_data_index = l;
            patch_data_index_set.insert(patch_data_index);
        }
    }
    setPatchDataIndices(patch_data_index_set);
    return;
}// setPatchDataIndices

void
CartSideRobinPhysBdryOp::setPhysicalBcCoefs(
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(bc_coefs[d] != NULL);
    }
#endif
    d_bc_coefs = bc_coefs;
    return;
}// setPhysicalBcCoefs

void
CartSideRobinPhysBdryOp::setHomogeneousBc(
    bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
CartSideRobinPhysBdryOp::setPhysicalBoundaryConditions(
    Patch<NDIM>& patch,
    const double fill_time,
    const IntVector<NDIM>& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector<NDIM>(0)) return;

    // Ensure the target patch data corresponds to a side centered variable and
    // that the proper number of boundary condition objects have been provided.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        var_db->mapIndexToVariable(patch_data_idx, var);
        Pointer<SideVariable<NDIM,double> > sc_var = var;
        if (sc_var.isNull())
        {
            TBOX_ERROR("CartSideRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                       << "  patch data index " << patch_data_idx << " does not correspond to a side-centered double precision variable." << std::endl);
        }
    }

    // Indicate whether we are employing homogeneous or inhomogeneous boundary
    // conditions for all extended Robin BC coef strategy objects employed by
    // this object.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (extended_bc_coef != NULL)
        {
            extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
        }
    }

    // Compute the boundary boxes.
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
#if (NDIM > 1)
    const Array<BoundaryBox<NDIM> > physical_codim2_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(patch);
#if (NDIM > 2)
    const Array<BoundaryBox<NDIM> > physical_codim3_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim3Boxes(patch);
#endif
#endif
#if (NDIM > 1)
    // Compute the portions of the co-dimension one boundary which abut the
    // co-dimension 2 boundary.
    //
    // NOTE: There is no need to include contributions from the co-dimension 3
    // boundary.  The intersection of the grown co-dimension 3 boundary boxes
    // with the co-dimension 1 boundary is covered by the intersection of the
    // grown co-dimension 2 boundary boxes with the co-dimension 1 boundary.
    Array<BoundaryBox<NDIM> > physical_codim1_reset_boxes(physical_codim1_boxes.size()*physical_codim2_boxes.size());
    int reset_box_counter = 0;
    for (int n2 = 0; n2 < physical_codim2_boxes.size(); ++n2)
    {
        const BoundaryBox<NDIM>& codim2_bdry_box = physical_codim2_boxes[n2];
        static const int GROWTH_FACTOR = 2;
        const Box<NDIM>& grown_codim2_box = Box<NDIM>::grow(codim2_bdry_box.getBox(),GROWTH_FACTOR);
        for (int n1 = 0; n1 < physical_codim1_boxes.size(); ++n1)
        {
            const BoundaryBox<NDIM>& codim1_bdry_box = physical_codim1_boxes[n1];
            const Box<NDIM>& codim1_box = codim1_bdry_box.getBox();
            const Box<NDIM> intersection_box = codim1_box*grown_codim2_box;
            if (!intersection_box.empty())
            {
                physical_codim1_reset_boxes[reset_box_counter] = BoundaryBox<NDIM>(intersection_box, codim1_bdry_box.getBoundaryType(), codim1_bdry_box.getLocationIndex());
                ++reset_box_counter;
            }
        }
    }
    physical_codim1_reset_boxes.resizeArray(reset_box_counter);
#endif
    // To set the boundary condition values:
    // (1) Compute the boundary condition coefficients along the co-dimension one
    //     boundary,
    // (2) Set the boundary conditions along the co-dimension one boundary boxes
    //     (thereby ensuring that Dirichlet boundary conditions set during step
    //     (1) are used properly in setting tangential boundary conditions),
    // (3) Recompute the boundary condition coefficients if the patch touches a
    //     co-dimension two or co-dimension three boundary box,
    // (4) Reset the boundary conditions along the co-dimension one boundary
    //     where it abuts the co-dimension two boundary,
    // (5) Extrapolate boundary values to the co-dimension two and three boundary
    //     boxes.
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        setCodimension1BdryValues(
            patch_data_idx,
            physical_codim1_boxes, fill_time, ghost_width_to_fill, patch);
    }
#if (NDIM > 1)
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        setCodimension1BdryValues(patch_data_idx, physical_codim1_reset_boxes, fill_time, ghost_width_to_fill, patch);
    }
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        setCodimension2BdryValues(patch_data_idx, physical_codim2_boxes, ghost_width_to_fill, patch);
    }
#if (NDIM > 2)
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        setCodimension3BdryValues(patch_data_idx, physical_codim3_boxes, ghost_width_to_fill, patch);
    }
#endif
#endif
    return;
}// setPhysicalBoundaryConditions

IntVector<NDIM>
CartSideRobinPhysBdryOp::getRefineOpStencilWidth() const
{
    return REFINE_OP_STENCIL_WIDTH;
}// getRefineOpStencilWidth

void
CartSideRobinPhysBdryOp::preprocessRefine(
    Patch<NDIM>& /*fine*/,
    const Patch<NDIM>& /*coarse*/,
    const Box<NDIM>& /*fine_box*/,
    const IntVector<NDIM>& /*ratio*/)
{
    // intentionally blank
    return;
}// preprocessRefine

void
CartSideRobinPhysBdryOp::postprocessRefine(
    Patch<NDIM>& /*fine*/,
    const Patch<NDIM>& /*coarse*/,
    const Box<NDIM>& /*fine_box*/,
    const IntVector<NDIM>& /*ratio*/)
{
    // intentionally blank
    return;
}// postprocessRefine

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CartSideRobinPhysBdryOp::setCodimension1BdryValues(
    const int patch_data_idx,
    const Array<BoundaryBox<NDIM> >& physical_codim1_boxes,
    const double fill_time,
    const IntVector<NDIM>& ghost_width_to_fill,
    Patch<NDIM>& patch)
{
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    if (n_physical_codim1_boxes == 0) return;

    Pointer<SideData<NDIM,double> > patch_data = patch.getPatchData(patch_data_idx);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;
    var_db->mapIndexToVariable(patch_data_idx, var);
    const int U_gcw = (patch_data->getGhostCellWidth()).max();
#ifdef DEBUG_CHECK_ASSERTIONS
    if (U_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR("CartSideRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                   "  patch data for patch data index " << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
#endif
    const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);
    blitz::TinyVector<double*,NDIM> U;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        U[axis] = patch_data->getPointer(axis);
    }

    const Box<NDIM>& patch_box = patch.getBox();
    const Index<NDIM>& patch_lower = patch_box.lower();
    const Index<NDIM>& patch_upper = patch_box.upper();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const patch_x_lower = pgeom->getXLower();
    const double* const patch_x_upper = pgeom->getXUpper();

    blitz::TinyVector<Box<NDIM>,NDIM> side_box;
    blitz::TinyVector<Index<NDIM>,NDIM> side_box_lower, side_box_upper;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        side_box[axis] = SideGeometry<NDIM>::toSideBox(patch_box,axis);
        side_box_lower[axis] = side_box[axis].lower();
        side_box_upper[axis] = side_box[axis].upper();
    }

    // Set the normal components.
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const unsigned int location_index   = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const unsigned int axis = bdry_normal_axis;
        const BoundaryBox<NDIM> trimmed_bdry_box(bdry_box.getBox()*bc_fill_box, bdry_box.getBoundaryType(), location_index);
        const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
        Pointer<ArrayData<NDIM,double> > acoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM,double> > bcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM,double> > gcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);

        // Set the boundary condition coefficients.
        RobinBcCoefStrategy<NDIM>* bc_coef = d_bc_coefs[axis];
        ExtendedRobinBcCoefStrategy* const extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coef);
        const bool using_extended_robin_bc_coef = extended_bc_coef != NULL;
        if (using_extended_robin_bc_coef)
        {
            extended_bc_coef->setTargetPatchDataIndex(patch_data_idx);
        }
        bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, var, patch, trimmed_bdry_box, fill_time);
        if (d_homogeneous_bc && !using_extended_robin_bc_coef) gcoef_data->fillAll(0.0);

        // Set the boundary values.
        if (location_index == 0 || location_index == 1)
        {
            SC_ROBIN_PHYS_BDRY_OP_1_X_FC(
                U[axis], U_gcw,
                acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                location_index,
                patch_lower(0), patch_upper(0),
                patch_lower(1), patch_upper(1),
#if (NDIM == 3)
                patch_lower(2), patch_upper(2),
#endif
                bc_coef_box.lower()(1), bc_coef_box.upper()(1),
#if (NDIM == 3)
                bc_coef_box.lower()(2), bc_coef_box.upper()(2),
#endif
                dx);
        }
        else if (location_index == 2 || location_index == 3)
        {
            SC_ROBIN_PHYS_BDRY_OP_1_Y_FC(
                U[axis], U_gcw,
                acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                location_index,
                patch_lower(0), patch_upper(0),
                patch_lower(1), patch_upper(1),
#if (NDIM == 3)
                patch_lower(2), patch_upper(2),
#endif
                bc_coef_box.lower()(0), bc_coef_box.upper()(0),
#if (NDIM == 3)
                bc_coef_box.lower()(2), bc_coef_box.upper()(2),
#endif
                dx);
        }
#if (NDIM == 3)
        else if (location_index == 4 || location_index == 5)
        {
            SC_ROBIN_PHYS_BDRY_OP_1_Z_FC(
                U[axis], U_gcw,
                acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                location_index,
                patch_lower(0), patch_upper(0),
                patch_lower(1), patch_upper(1),
                patch_lower(2), patch_upper(2),
                bc_coef_box.lower()(0), bc_coef_box.upper()(0),
                bc_coef_box.lower()(1), bc_coef_box.upper()(1),
                dx);
        }
#endif
    }

    // Set the transverse components.
    const IntVector<NDIM>& ratio_to_level_zero = pgeom->getRatio();
    Array<Array<bool> > touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        touches_regular_bdry [axis].resizeArray(2);
        touches_periodic_bdry[axis].resizeArray(2);
        for (int upperlower = 0; upperlower < 2; ++upperlower)
        {
            touches_regular_bdry [axis][upperlower] = pgeom->getTouchesRegularBoundary( axis,upperlower);
            touches_periodic_bdry[axis][upperlower] = pgeom->getTouchesPeriodicBoundary(axis,upperlower);
        }
    }
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const unsigned int location_index   = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const BoundaryBox<NDIM> trimmed_bdry_box(bdry_box.getBox()*bc_fill_box, bdry_box.getBoundaryType(), location_index);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            if (axis != bdry_normal_axis)
            {
                const Box<NDIM> bc_coef_box = compute_tangential_extension(PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), axis);
                Pointer<ArrayData<NDIM,double> > acoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM,double> > bcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM,double> > gcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);

                // Temporarily reset the patch geometry object associated with
                // the patch so that boundary conditions are set at the correct
                // spatial locations.
                blitz::TinyVector<double,NDIM> shifted_patch_x_lower, shifted_patch_x_upper;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    shifted_patch_x_lower[d] = patch_x_lower[d];
                    shifted_patch_x_upper[d] = patch_x_upper[d];
                }
                shifted_patch_x_lower[axis] -= 0.5*dx[axis];
                shifted_patch_x_upper[axis] -= 0.5*dx[axis];
                patch.setPatchGeometry(new CartesianPatchGeometry<NDIM>(ratio_to_level_zero, touches_regular_bdry, touches_periodic_bdry, dx, shifted_patch_x_lower.data(), shifted_patch_x_upper.data()));

                // Set the boundary condition coefficients.
                RobinBcCoefStrategy<NDIM>* bc_coef = d_bc_coefs[axis];
                ExtendedRobinBcCoefStrategy* const extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coef);
                const bool using_extended_robin_bc_coef = extended_bc_coef != NULL;
                if (using_extended_robin_bc_coef)
                {
                    extended_bc_coef->setTargetPatchDataIndex(patch_data_idx);
                }
                bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, var, patch, trimmed_bdry_box, fill_time);
                if (d_homogeneous_bc && !using_extended_robin_bc_coef) gcoef_data->fillAll(0.0);

                // Restore the original patch geometry object.
                patch.setPatchGeometry(pgeom);

                // Set the boundary values.
                if (location_index == 0 || location_index == 1)
                {
                    CC_ROBIN_PHYS_BDRY_OP_1_X_FC(
                        U[axis], U_gcw,
                        acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                        location_index,
                        side_box_lower[axis](0), side_box_upper[axis](0),
                        side_box_lower[axis](1), side_box_upper[axis](1),
#if (NDIM == 3)
                        side_box_lower[axis](2), side_box_upper[axis](2),
#endif
                        bc_coef_box.lower()(1), bc_coef_box.upper()(1),
#if (NDIM == 3)
                        bc_coef_box.lower()(2), bc_coef_box.upper()(2),
#endif
                        dx);
                }
                else if (location_index == 2 || location_index == 3)
                {
                    CC_ROBIN_PHYS_BDRY_OP_1_Y_FC(
                        U[axis], U_gcw,
                        acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                        location_index,
                        side_box_lower[axis](0), side_box_upper[axis](0),
                        side_box_lower[axis](1), side_box_upper[axis](1),
#if (NDIM == 3)
                        side_box_lower[axis](2), side_box_upper[axis](2),
#endif
                        bc_coef_box.lower()(0), bc_coef_box.upper()(0),
#if (NDIM == 3)
                        bc_coef_box.lower()(2), bc_coef_box.upper()(2),
#endif
                        dx);
                }

#if (NDIM == 3)
                else if (location_index == 4 || location_index == 5)
                {
                    CC_ROBIN_PHYS_BDRY_OP_1_Z_FC(
                        U[axis], U_gcw,
                        acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                        location_index,
                        side_box_lower[axis](0), side_box_upper[axis](0),
                        side_box_lower[axis](1), side_box_upper[axis](1),
                        side_box_lower[axis](2), side_box_upper[axis](2),
                        bc_coef_box.lower()(0), bc_coef_box.upper()(0),
                        bc_coef_box.lower()(1), bc_coef_box.upper()(1),
                        dx);
                }
#endif
            }
        }
    }
    return;
}// setCodimension1BdryValues

#if (NDIM > 1)

void
CartSideRobinPhysBdryOp::setCodimension2BdryValues(
    const int patch_data_idx,
    const Array<BoundaryBox<NDIM> >& physical_codim2_boxes,
    const IntVector<NDIM>& ghost_width_to_fill,
    const Patch<NDIM>& patch)
{
    const int n_physical_codim2_boxes = physical_codim2_boxes.size();
    if (n_physical_codim2_boxes == 0) return;

    Pointer<SideData<NDIM,double> > patch_data = patch.getPatchData(patch_data_idx);
    const int U_gcw = (patch_data->getGhostCellWidth()).max();
#ifdef DEBUG_CHECK_ASSERTIONS
    if (U_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR("CartSideRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                   "  patch data for patch data index " << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
#endif
    const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);
    blitz::TinyVector<double*,NDIM> U;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        U[axis] = patch_data->getPointer(axis);
    }

    const Box<NDIM>& patch_box = patch.getBox();
    const Index<NDIM>& patch_lower = patch_box.lower();
    const Index<NDIM>& patch_upper = patch_box.upper();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

#if (NDIM == 3)
    blitz::TinyVector<Box<NDIM>,NDIM> side_box;
    blitz::TinyVector<Index<NDIM>,NDIM> side_box_lower, side_box_upper;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        side_box[axis] = SideGeometry<NDIM>::toSideBox(patch_box,axis);
        side_box_lower[axis] = side_box[axis].lower();
        side_box_upper[axis] = side_box[axis].upper();
    }
#endif

    for (int n = 0; n < n_physical_codim2_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim2_boxes[n];
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const unsigned int location_index = bdry_box.getLocationIndex();
        SC_ROBIN_PHYS_BDRY_OP_2_FC(
            U[0], U[1],
#if (NDIM == 3)
            U[2],
#endif
            U_gcw,
            location_index,
            patch_lower(0), patch_upper(0),
            patch_lower(1), patch_upper(1),
#if (NDIM == 3)
            patch_lower(2), patch_upper(2),
#endif
            bc_fill_box.lower()(0), bc_fill_box.upper()(0),
            bc_fill_box.lower()(1), bc_fill_box.upper()(1)
#if (NDIM == 3)
            ,bc_fill_box.lower()(2),bc_fill_box.upper()(2)
#endif
                                   );
#if (NDIM == 3)
        if (location_index < 4)
        {
            static const unsigned int axis = 0;
            CC_ROBIN_PHYS_BDRY_OP_2_FC(
                U[axis], U_gcw,
                location_index,
                side_box_lower[axis](0), side_box_upper[axis](0),
                side_box_lower[axis](1), side_box_upper[axis](1),
                side_box_lower[axis](2), side_box_upper[axis](2),
                bc_fill_box.lower()(0), bc_fill_box.upper()(0),
                bc_fill_box.lower()(1), bc_fill_box.upper()(1),
                bc_fill_box.lower()(2), bc_fill_box.upper()(2));
        }
        else if (location_index >= 4 && location_index < 8)
        {
            static const unsigned int axis = 1;
            CC_ROBIN_PHYS_BDRY_OP_2_FC(
                U[axis], U_gcw,
                location_index,
                side_box_lower[axis](0), side_box_upper[axis](0),
                side_box_lower[axis](1), side_box_upper[axis](1),
                side_box_lower[axis](2), side_box_upper[axis](2),
                bc_fill_box.lower()(0), bc_fill_box.upper()(0),
                bc_fill_box.lower()(1), bc_fill_box.upper()(1),
                bc_fill_box.lower()(2), bc_fill_box.upper()(2));
        }
        else if (location_index >= 8 && location_index < 12)
        {
            static const unsigned int axis = 2;
            CC_ROBIN_PHYS_BDRY_OP_2_FC(
                U[axis], U_gcw,
                location_index,
                side_box_lower[axis](0), side_box_upper[axis](0),
                side_box_lower[axis](1), side_box_upper[axis](1),
                side_box_lower[axis](2), side_box_upper[axis](2),
                bc_fill_box.lower()(0), bc_fill_box.upper()(0),
                bc_fill_box.lower()(1), bc_fill_box.upper()(1),
                bc_fill_box.lower()(2), bc_fill_box.upper()(2));
        }
#endif
    }
    return;
}// setCodimension2BdryValues

#if (NDIM > 2)

void
CartSideRobinPhysBdryOp::setCodimension3BdryValues(
    const int patch_data_idx,
    const Array<BoundaryBox<NDIM> >& physical_codim3_boxes,
    const IntVector<NDIM>& ghost_width_to_fill,
    const Patch<NDIM>& patch)
{
    const int n_physical_codim3_boxes = physical_codim3_boxes.size();
    if (n_physical_codim3_boxes == 0) return;

    Pointer<SideData<NDIM,double> > patch_data = patch.getPatchData(patch_data_idx);
    const int U_gcw = (patch_data->getGhostCellWidth()).max();
#ifdef DEBUG_CHECK_ASSERTIONS
    if (U_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR("CartSideRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                   "  patch data for patch data index " << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
#endif
    const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);
    blitz::TinyVector<double*,NDIM> U;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        U[axis] = patch_data->getPointer(axis);
    }

    const Box<NDIM>& patch_box = patch.getBox();
    const Index<NDIM>& patch_lower = patch_box.lower();
    const Index<NDIM>& patch_upper = patch_box.upper();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

    for (int n = 0; n < n_physical_codim3_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim3_boxes[n];
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const unsigned int location_index = bdry_box.getLocationIndex();
        SC_ROBIN_PHYS_BDRY_OP_3_FC(
            U[0], U[1], U[2], U_gcw,
            location_index,
            patch_lower(0), patch_upper(0),
            patch_lower(1), patch_upper(1),
            patch_lower(2), patch_upper(2),
            bc_fill_box.lower()(0), bc_fill_box.upper()(0),
            bc_fill_box.lower()(1), bc_fill_box.upper()(1),
            bc_fill_box.lower()(2), bc_fill_box.upper()(2));
    }
    return;
}// setCodimension3BdryValues

#endif
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::CartSideRobinPhysBdryOp>;

//////////////////////////////////////////////////////////////////////////////
