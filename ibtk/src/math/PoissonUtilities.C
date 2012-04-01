// Filename: PoissonUtilities.C
// Created on 31 Mar 2012 by Boyce Griffith
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

#include "PoissonUtilities.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <ArrayDataBasicOps.h>
#include <CartesianPatchGeometry.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

void
PoissonUtilities::computeCCMatrixCoefficients(
    Pointer<Patch<NDIM> > patch,
    CellData<NDIM,double>& diagonal,
    SideData<NDIM,double>& off_diagonal,
    const PoissonSpecifications& poisson_spec,
    RobinBcCoefStrategy<NDIM>* bc_coef,
    double data_time)
{
    computeCCMatrixCoefficients(patch, diagonal, off_diagonal, poisson_spec, std::vector<RobinBcCoefStrategy<NDIM>*>(1,bc_coef), data_time);
    return;
}// computeCCMatrixCoefficients

void
PoissonUtilities::computeCCMatrixCoefficients(
    Pointer<Patch<NDIM> > patch,
    CellData<NDIM,double>& diagonal,
    SideData<NDIM,double>& off_diagonal,
    const PoissonSpecifications& poisson_spec,
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs,
    double data_time)
{
    computeCCMatrixCoefficients(patch, diagonal, off_diagonal, poisson_spec, std::vector<RobinBcCoefStrategy<NDIM>*>(&bc_coefs[0],&bc_coefs[0]+NDIM), data_time);
    return;
}// computeCCMatrixCoefficients

void
PoissonUtilities::computeCCMatrixCoefficients(
    Pointer<Patch<NDIM> > patch,
    CellData<NDIM,double>& diagonal,
    SideData<NDIM,double>& off_diagonal,
    const PoissonSpecifications& poisson_spec,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    double data_time)
{
    const int depth = bc_coefs.size();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(    diagonal.getDepth() == depth);
    TBOX_ASSERT(off_diagonal.getDepth() == depth);
#endif

    ArrayDataBasicOps<NDIM,double> array_ops;
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    // Compute all off-diagonal matrix coefficients for all cell sides,
    // including those that touch the physical boundary, ignoring physical
    // boundary conditions for now.
    if (!poisson_spec.dIsConstant())
    {
        off_diagonal.copy(*patch->getPatchData(poisson_spec.getDPatchDataId()));
    }
    else
    {
        off_diagonal.fill(poisson_spec.getDConstant());
    }

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        array_ops.scale(off_diagonal.getArrayData(axis),
                        1.0/(dx[axis]*dx[axis]),
                        off_diagonal.getArrayData(axis),
                        side_box);
    }

    // Compute the diagonal matrix coefficients.
    if (!poisson_spec.cIsZero() && !poisson_spec.cIsConstant())
    {
        diagonal.copy(*patch->getPatchData(poisson_spec.getCPatchDataId()));
    }
    else
    {
        if (poisson_spec.cIsZero()) diagonal.fill(0.0);
        else diagonal.fill(poisson_spec.getCConstant());
    }

    for (int d = 0; d < depth; ++d)
    {
        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            Index<NDIM> i = b();
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const SideIndex<NDIM> ilower(i, axis, SideIndex<NDIM>::Lower);
                diagonal(i,d) -= off_diagonal(ilower,d);
                const SideIndex<NDIM> iupper(i, axis, SideIndex<NDIM>::Upper);
                diagonal(i,d) -= off_diagonal(iupper,d);
            }
        }
    }

    // Modify the diagonal and off-diagonal entries to account for
    // homogeneous boundary conditions.
    //
    // Here, we follow the same linear extrapolation approach implemented in
    // class CartesianRobinBcHelper.  Namely, with u_i
    // denoting the interior cell, u_o denoting the ghost cell, and u_b and
    // u_n denoting the value and normal derivative of u at the boundary,
    //
    //     u_b = (u_i + u_o)/2   and   u_n = (u_o - u_i)/h
    //
    // Now, if
    //
    //     a*u_b + b*u_n = 0
    //
    // then
    //
    //     u_o = -((a*h - 2*b)/(a*h + 2*b))*u_i
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
        const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

        ArrayData<NDIM,double> acoef_data(bc_coef_box, 1);
        ArrayData<NDIM,double> bcoef_data(bc_coef_box, 1);

        Pointer<ArrayData<NDIM,double> > acoef_data_ptr(&acoef_data, false);
        Pointer<ArrayData<NDIM,double> > bcoef_data_ptr(&bcoef_data, false);
        Pointer<ArrayData<NDIM,double> > gcoef_data_ptr(NULL);

        for (int d = 0; d < depth; ++d)
        {
            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[d]);
             if (extended_bc_coef != NULL) extended_bc_coef->setHomogeneousBc(true);
             bc_coefs[d]->setBcCoefs(
                 acoef_data_ptr, bcoef_data_ptr, gcoef_data_ptr, NULL,
                 *patch, trimmed_bdry_box, data_time);

             const unsigned int location_index = bdry_box.getLocationIndex();
             const unsigned int bdry_normal_axis =  location_index / 2;
             const bool bdry_lower_side = (location_index % 2) == 0;
             const bool bdry_upper_side = (location_index % 2) != 0;

             for (Box<NDIM>::Iterator b(bc_coef_box); b; b++)
             {
                 const Index<NDIM>& i_s_bdry = b();
                 const double& a = acoef_data(i_s_bdry,0);
                 const double& b = bcoef_data(i_s_bdry,0);
                 const double& h = dx[bdry_normal_axis];

                 // i_s_bdry: side index located on physical boundary
                 //
                 // i_c_intr: cell index located adjacent to physical boundary
                 // in the patch interior
                 Index<NDIM> i_c_intr = i_s_bdry;
                 if (bdry_upper_side)
                 {
                     i_c_intr(bdry_normal_axis) -= 1;
                 }

                 if (bdry_lower_side)
                 {
                     const SideIndex<NDIM> ilower(i_c_intr, bdry_normal_axis, SideIndex<NDIM>::Lower);
                     diagonal(i_c_intr,d) += off_diagonal(ilower,d)*(-(a*h-2.0*b)/(a*h+2.0*b));
                     off_diagonal(ilower,d) = 0.0;
                 }

                 if (bdry_upper_side)
                 {
                     const SideIndex<NDIM> iupper(i_c_intr, bdry_normal_axis, SideIndex<NDIM>::Upper);
                     diagonal(i_c_intr,d) += off_diagonal(iupper,d)*(-(a*h-2.0*b)/(a*h+2.0*b));
                     off_diagonal(iupper,d) = 0.0;
                 }
             }
         }
    }
    return;
}// computeCCMatrixCoefficients

void
PoissonUtilities::adjustCCBoundaryRhsEntries(
    Pointer<Patch<NDIM> > patch,
    CellData<NDIM,double>& rhs_data,
    const PoissonSpecifications& poisson_spec,
    RobinBcCoefStrategy<NDIM>* bc_coef,
    double data_time)
{
    adjustCCBoundaryRhsEntries(patch, rhs_data, poisson_spec, std::vector<RobinBcCoefStrategy<NDIM>*>(1,bc_coef), data_time);
    return;
}// adjustCCBoundaryRhsEntries

void
PoissonUtilities::adjustCCBoundaryRhsEntries(
    Pointer<Patch<NDIM> > patch,
    CellData<NDIM,double>& rhs_data,
    const PoissonSpecifications& poisson_spec,
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs,
    double data_time)
{
    adjustCCBoundaryRhsEntries(patch, rhs_data, poisson_spec, std::vector<RobinBcCoefStrategy<NDIM>*>(&bc_coefs[0],&bc_coefs[0]+NDIM), data_time);
    return;
}// adjustCCBoundaryRhsEntries

void
PoissonUtilities::adjustCCBoundaryRhsEntries(
    Pointer<Patch<NDIM> > patch,
    CellData<NDIM,double>& rhs_data,
    const PoissonSpecifications& poisson_spec,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    double data_time)
{
    const int depth = bc_coefs.size();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(rhs_data.getDepth() == depth);
#endif

    const Box<NDIM>& patch_box = patch->getBox();
    OutersideData<NDIM,double> D_data(patch_box, depth);
    if (!poisson_spec.dIsConstant())
    {
        D_data.copy(*patch->getPatchData(poisson_spec.getDPatchDataId()));
    }
    else
    {
        D_data.fillAll(poisson_spec.getDConstant());
    }

    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const Array<BoundaryBox<NDIM> > codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);

    // Modify the rhs entries to account for inhomogeneous boundary conditions.
    //
    // Here, we follow the same linear extrapolation approach implemented in
    // class CartesianRobinBcHelper.  Namely, with u_i denoting
    // the interior cell, u_o denoting the ghost cell, and u_b and u_n denoting
    // the value and normal derivative of u at the boundary,
    //
    //     u_b = (u_i + u_o)/2   and   u_n = (u_o - u_i)/h
    //
    // Now, if
    //
    //     a*u_b + b*u_n = g
    //
    // then, with u_i = 0,
    //
    //     u_o = 2*h*g/(2*b + a*h)
    //
    // so that the boundary flux is
    //
    //     (u_i - u_o)/h = -2*g/(2*b + h*a)
    //
    // In this loop, we modify the rhs entries appropriately.
    const int n_bdry_boxes = codim1_boxes.size();
    for (int n = 0; n < n_bdry_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = codim1_boxes[n];
        const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
        const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

        ArrayData<NDIM,double> acoef_data(bc_coef_box, 1);
        ArrayData<NDIM,double> bcoef_data(bc_coef_box, 1);
        ArrayData<NDIM,double> gcoef_data(bc_coef_box, 1);

        Pointer<ArrayData<NDIM,double> > acoef_data_ptr(&acoef_data, false);
        Pointer<ArrayData<NDIM,double> > bcoef_data_ptr(&bcoef_data, false);
        Pointer<ArrayData<NDIM,double> > gcoef_data_ptr(&gcoef_data, false);

        for (int d = 0; d < depth; ++d)
        {
            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[d]);
             if (extended_bc_coef != NULL) extended_bc_coef->setHomogeneousBc(false);
             bc_coefs[d]->setBcCoefs(
                 acoef_data_ptr, bcoef_data_ptr, gcoef_data_ptr, NULL,
                 *patch, trimmed_bdry_box, data_time);

             const unsigned int location_index = bdry_box.getLocationIndex();
             const unsigned int bdry_normal_axis =  location_index / 2;
             const bool bdry_upper_side = (location_index % 2) != 0;
             const int bdry_side = (bdry_upper_side ? 1 : 0);

             // i_s_bdry: side index located on physical boundary
             //
             // i_c_intr: cell index located adjacent to physical boundary in
             // the patch interior
             for (Box<NDIM>::Iterator b(bc_coef_box); b; b++)
             {
                 const Index<NDIM>& i_s_bdry = b();
                 const double& a = acoef_data(i_s_bdry,0);
                 const double& b = bcoef_data(i_s_bdry,0);
                 const double& g = gcoef_data(i_s_bdry,0);
                 const double& h = dx[bdry_normal_axis];

                 Index<NDIM> i_c_intr = i_s_bdry;
                 if (bdry_upper_side)
                 {
                     i_c_intr(bdry_normal_axis) -= 1;
                 }
                 rhs_data(i_c_intr,d) += (D_data.getArrayData(bdry_normal_axis,bdry_side)(i_s_bdry,d)/h)*(-2.0*g)/(2.0*b+h*a);
             }
        }
    }
    return;
}// adjustCCBoundaryRhsEntries

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
