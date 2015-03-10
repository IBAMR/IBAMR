// Filename: PoissonUtilities.cpp
// Created on 31 Mar 2012 by Boyce Griffith
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

#include <stddef.h>
#include <functional>
#include <map>
#include <ostream>
#include <vector>

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/math/ArrayDataBasicOps.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/solv/RobinBcCoefStrategy.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/hier/Variable.h"
#include "boost/array.hpp"
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/PoissonUtilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "ibtk/ibtk_utilities.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include "SAMRAI/tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline Box compute_tangential_extension(const Box& box, const int data_axis)
{
    Box extended_box = box;
    extended_box.setUpper(data_axis, box.upper(data_axis) + 1);
    return extended_box;
}
}

void PoissonUtilities::computeCCMatrixCoefficients(boost::shared_ptr<Patch> patch,
                                                   CellData<double>& matrix_coefficients,
                                                   const std::vector<Index>& stencil,
                                                   const PoissonSpecifications& poisson_spec,
                                                   RobinBcCoefStrategy* bc_coef,
                                                   double data_time)
{
    computeCCMatrixCoefficients(patch, matrix_coefficients, stencil, poisson_spec,
                                std::vector<RobinBcCoefStrategy*>(1, bc_coef), data_time);
    return;
}

void PoissonUtilities::computeCCComplexMatrixCoefficients(boost::shared_ptr<Patch> patch,
                                                          CellData<double>& matrix_coefficients,
                                                          const std::vector<Index>& stencil,
                                                          const PoissonSpecifications& poisson_spec_real,
                                                          const PoissonSpecifications& poisson_spec_imag,
                                                          RobinBcCoefStrategy* bc_coef,
                                                          double data_time)
{
    computeCCComplexMatrixCoefficients(patch, matrix_coefficients, stencil, poisson_spec_real, poisson_spec_imag,
                                       std::vector<RobinBcCoefStrategy*>(2, bc_coef), data_time);
    return;
}

void PoissonUtilities::computeCCMatrixCoefficients(boost::shared_ptr<Patch> patch,
                                                   CellData<double>& matrix_coefficients,
                                                   const std::vector<Index>& stencil,
                                                   const PoissonSpecifications& poisson_spec,
                                                   const std::vector<RobinBcCoefStrategy*>& bc_coefs,
                                                   double data_time)
{
    const int stencil_sz = static_cast<int>(stencil.size());
    TBOX_ASSERT(stencil_sz == 2 * NDIM + 1);
    std::map<Index, int> stencil_map;
    for (int k = 0; k < stencil_sz; ++k)
    {
        stencil_map[stencil[k]] = k;
    }
    const int stencil_index_diag = stencil_map[Index::getZeroIndex(DIM)];
    boost::array<int, NDIM> stencil_index_lower, stencil_index_upper;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        Index ilower(DIM, 0), iupper(DIM, 0);
        ilower(axis) = -1;
        iupper(axis) = +1;
        stencil_index_lower[axis] = stencil_map[ilower];
        stencil_index_upper[axis] = stencil_map[iupper];
    }
    const int depth = static_cast<int>(bc_coefs.size());
    TBOX_ASSERT(matrix_coefficients.getDepth() == depth * stencil_sz);
    const Box& patch_box = patch->getBox();
    CellData<double> diagonal(patch_box, depth, IntVector::getZero(DIM));
    SideData<double> off_diagonal(patch_box, depth, IntVector::getZero(DIM));

    ArrayDataBasicOps<double> array_ops;
    auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch->getPatchGeometry());
    const double* const dx = pgeom->getDx();

    // Compute all off-diagonal matrix coefficients for all cell sides,
    // including those that touch the physical boundary; however, do not yet
    // take physical boundary conditions into account.  Boundary conditions are
    // handled subsequently.
    if (!poisson_spec.dIsConstant())
    {
        off_diagonal.copy(*patch->getPatchData(poisson_spec.getDPatchDataId()));
    }
    else
    {
        off_diagonal.fillAll(poisson_spec.getDConstant());
    }

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        Box side_box = SideGeometry::toSideBox(patch_box, axis);
        array_ops.scale(off_diagonal.getArrayData(axis), 1.0 / (dx[axis] * dx[axis]), off_diagonal.getArrayData(axis),
                        side_box);
    }

    // Compute all diagonal matrix coefficients for all cells, including those
    // that are adjacent to the physical boundary; however, do not yet take
    // physical boundary conditions into account.  Boundary conditions are
    // handled subsequently.
    if (!poisson_spec.cIsZero() && !poisson_spec.cIsConstant())
    {
        diagonal.copy(*patch->getPatchData(poisson_spec.getCPatchDataId()));
    }
    else
    {
        if (poisson_spec.cIsZero())
            diagonal.fillAll(0.0);
        else
            diagonal.fillAll(poisson_spec.getCConstant());
    }

    for (int d = 0; d < depth; ++d)
    {
        for (CellIterator b = CellGeometry::begin(patch_box); b != CellGeometry::end(patch_box); ++b)
        {
            const CellIndex& i = *b;
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const SideIndex ilower(i, axis, SideIndex::Lower);
                diagonal(i, d) -= off_diagonal(ilower, d);
                const SideIndex iupper(i, axis, SideIndex::Upper);
                diagonal(i, d) -= off_diagonal(iupper, d);
            }
        }
    }

    // Modify matrix coefficients to account for physical boundary conditions.
    const std::vector<BoundaryBox> physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    const int n_physical_codim1_boxes = static_cast<int>(physical_codim1_boxes.size());
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox& bdry_box = physical_codim1_boxes[n];
        const BoundaryBox trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
        const Box bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

        auto acoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
        auto bcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
        auto gcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);

        for (int d = 0; d < depth; ++d)
        {
            static const bool homogeneous_bc = true;
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[d]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[d]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool bdry_lower_side = (location_index % 2) == 0;
            const bool bdry_upper_side = (location_index % 2) != 0;

            // Modify the diagonal and off-diagonal entries to account for
            // homogeneous boundary conditions.
            //
            // Here, we follow the same linear extrapolation approach
            // implemented in class CartesianRobinBcHelper.  Namely, with u_i
            // denoting the interior cell, u_o denoting the ghost cell, and u_b
            // and u_n denoting the value and normal derivative of u at the
            // boundary,
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
            for (CellIterator bc = CellGeometry::begin(bc_coef_box); bc != CellGeometry::end(bc_coef_box); ++bc)
            {
                const CellIndex& i_s_bdry = *bc;
                const double& a = (*acoef_data)(i_s_bdry, 0);
                const double& b = (*bcoef_data)(i_s_bdry, 0);
                const double& h = dx[bdry_normal_axis];

                // i_s_bdry: side index located on physical boundary
                //
                // i_c_intr: cell index located adjacent to physical boundary
                // in the patch interior
                CellIndex i_c_intr = i_s_bdry;
                if (bdry_upper_side)
                {
                    i_c_intr(bdry_normal_axis) -= 1;
                }

                if (bdry_lower_side)
                {
                    const SideIndex ilower(i_c_intr, bdry_normal_axis, SideIndex::Lower);
                    diagonal(i_c_intr, d) += off_diagonal(ilower, d) * (-(a * h - 2.0 * b) / (a * h + 2.0 * b));
                    off_diagonal(ilower, d) = 0.0;
                }

                if (bdry_upper_side)
                {
                    const SideIndex iupper(i_c_intr, bdry_normal_axis, SideIndex::Upper);
                    diagonal(i_c_intr, d) += off_diagonal(iupper, d) * (-(a * h - 2.0 * b) / (a * h + 2.0 * b));
                    off_diagonal(iupper, d) = 0.0;
                }
            }
        }
    }

    // Setup the matrix coefficients.
    for (int d = 0; d < depth; ++d)
    {
        const unsigned int offset = static_cast<unsigned int>(d * stencil_sz);
        for (CellIterator b = CellGeometry::begin(patch_box); b != CellGeometry::end(patch_box); ++b)
        {
            const CellIndex& i = *b;
            matrix_coefficients(i, offset + stencil_index_diag) = diagonal(i, d);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const SideIndex ilower(i, axis, SideIndex::Lower);
                matrix_coefficients(i, offset + stencil_index_lower[axis]) = off_diagonal(ilower, d);
                const SideIndex iupper(i, axis, SideIndex::Upper);
                matrix_coefficients(i, offset + stencil_index_upper[axis]) = off_diagonal(iupper, d);
            }
        }
    }
    return;
}

void PoissonUtilities::computeCCComplexMatrixCoefficients(boost::shared_ptr<Patch> patch,
                                                          CellData<double>& matrix_coefficients,
                                                          const std::vector<Index>& stencil,
                                                          const PoissonSpecifications& poisson_spec_real,
                                                          const PoissonSpecifications& poisson_spec_imag,
                                                          const std::vector<RobinBcCoefStrategy*>& bc_coefs,
                                                          double data_time)
{
    const int stencil_sz = static_cast<int>(stencil.size());
    TBOX_ASSERT(stencil_sz == 2 * NDIM + 1);
    std::map<Index, int> stencil_map;
    for (int k = 0; k < stencil_sz; ++k)
    {
        stencil_map[stencil[k]] = k;
    }
    const int stencil_index_diag = stencil_map[Index::getZeroIndex(DIM)];
    boost::array<int, NDIM> stencil_index_lower, stencil_index_upper;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        Index ilower(DIM, 0), iupper(DIM, 0);
        ilower(axis) = -1;
        iupper(axis) = +1;
        stencil_index_lower[axis] = stencil_map[ilower];
        stencil_index_upper[axis] = stencil_map[iupper];
    }
    const int depth = static_cast<int>(bc_coefs.size());
    TBOX_ASSERT(depth >= 2 && depth % 2 == 0);
    TBOX_ASSERT(matrix_coefficients.getDepth() == depth * stencil_sz * 2);
    const Box& patch_box = patch->getBox();
    CellData<double> diagonal(patch_box, 2 * depth, IntVector::getZero(DIM));
    SideData<double> off_diagonal(patch_box, 2 * depth, IntVector::getZero(DIM));

    CellData<double> negC_imag(patch_box, 1, IntVector::getZero(DIM));
    SideData<double> negD_imag(patch_box, 1, IntVector::getZero(DIM));

    if (!poisson_spec_imag.cIsZero() && !poisson_spec_imag.cIsConstant())
    {
        negC_imag.copy(*patch->getPatchData(poisson_spec_imag.getCPatchDataId()));
    }
    else
    {
        if (poisson_spec_imag.cIsZero())
            negC_imag.fillAll(0.0);
        else
            negC_imag.fillAll(poisson_spec_imag.getCConstant());
    }
    if (!poisson_spec_imag.dIsConstant())
    {
        negD_imag.copy(*patch->getPatchData(poisson_spec_imag.getDPatchDataId()));
    }
    else
    {
        negD_imag.fill(poisson_spec_imag.getDConstant());
    }

    ArrayDataBasicOps<double> array_ops;
    auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch->getPatchGeometry());
    const double* const dx = pgeom->getDx();

    array_ops.scale(negC_imag.getArrayData(), -1.0, negC_imag.getArrayData(), patch_box);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        Box side_box = SideGeometry::toSideBox(patch_box, axis);
        array_ops.scale(negD_imag.getArrayData(axis), -1.0, negD_imag.getArrayData(axis), side_box);
    }

    // Compute all off-diagonal matrix coefficients for all cell sides,
    // including those that touch the physical boundary; however, do not yet
    // take physical boundary conditions into account.  Boundary conditions are
    // handled subsequently.

    for (int d = 0; d < depth; d = d + 2)
    {
        if (!poisson_spec_real.dIsConstant())
        {
            off_diagonal.copyDepth(
                2 * d, CPP_CAST<SideData<double>&>(*patch->getPatchData(poisson_spec_real.getDPatchDataId())), 0);
            off_diagonal.copyDepth(
                2 * d + 3, CPP_CAST<SideData<double>&>(*patch->getPatchData(poisson_spec_real.getDPatchDataId())), 0);
        }
        else
        {
            off_diagonal.fill(poisson_spec_real.getDConstant(), 2 * d);
            off_diagonal.fill(poisson_spec_real.getDConstant(), 2 * d + 3);
        }

        off_diagonal.copyDepth(2 * d + 1, negD_imag, 0);
        if (!poisson_spec_imag.dIsConstant())
        {
            off_diagonal.copyDepth(
                2 * d + 2, CPP_CAST<SideData<double>&>(*patch->getPatchData(poisson_spec_imag.getDPatchDataId())), 0);
        }
        else
        {
            off_diagonal.fill(poisson_spec_imag.getDConstant(), 2 * d + 2);
        }
    }

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        Box side_box = SideGeometry::toSideBox(patch_box, axis);
        array_ops.scale(off_diagonal.getArrayData(axis), 1.0 / (dx[axis] * dx[axis]), off_diagonal.getArrayData(axis),
                        side_box);
    }

    // Compute all diagonal matrix coefficients for all cells, including those
    // that are adjacent to the physical boundary; however, do not yet take
    // physical boundary conditions into account.  Boundary conditions are
    // handled subsequently.
    for (int d = 0; d < depth; d = d + 2)
    {
        if (!poisson_spec_real.cIsZero() && !poisson_spec_real.cIsConstant())
        {
            diagonal.copyDepth(
                2 * d, CPP_CAST<CellData<double>&>(*patch->getPatchData(poisson_spec_real.getCPatchDataId())), 0);
            diagonal.copyDepth(
                2 * d + 3, CPP_CAST<CellData<double>&>(*patch->getPatchData(poisson_spec_real.getCPatchDataId())), 0);
        }
        else
        {
            if (poisson_spec_real.cIsZero())
            {
                diagonal.fill(0.0, 2 * d);
                diagonal.fill(0.0, 2 * d + 3);
            }
            else
            {
                diagonal.fill(poisson_spec_real.getCConstant(), 2 * d);
                diagonal.fill(poisson_spec_real.getCConstant(), 2 * d + 3);
            }
        }

        diagonal.copyDepth(2 * d + 1, negC_imag, 0);
        if (!poisson_spec_imag.cIsZero() && !poisson_spec_imag.cIsConstant())
        {
            diagonal.copyDepth(
                2 * d + 2, CPP_CAST<CellData<double>&>(*patch->getPatchData(poisson_spec_imag.getCPatchDataId())), 0);
        }
        else
        {
            if (poisson_spec_imag.cIsZero())
                diagonal.fill(0.0, 2 * d + 2);
            else
                diagonal.fill(poisson_spec_real.getCConstant(), 2 * d + 2);
        }
    }

    for (int d = 0; d < 2 * depth; ++d)
    {
        for (CellIterator b = CellGeometry::begin(patch_box); b != CellGeometry::end(patch_box); ++b)
        {
            const CellIndex& i = *b;
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const SideIndex ilower(i, axis, SideIndex::Lower);
                diagonal(i, d) -= off_diagonal(ilower, d);
                const SideIndex iupper(i, axis, SideIndex::Upper);
                diagonal(i, d) -= off_diagonal(iupper, d);
            }
        }
    }

    // Modify matrix coefficients to account for physical boundary conditions.
    const std::vector<BoundaryBox> physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    const int n_physical_codim1_boxes = static_cast<int>(physical_codim1_boxes.size());
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox& bdry_box = physical_codim1_boxes[n];
        const BoundaryBox trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
        const Box bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

        std::vector<boost::shared_ptr<ArrayData<double> > > acoef_data(depth), bcoef_data(depth), gcoef_data(depth);
        for (int d = 0; d < depth; ++d)
        {
            acoef_data[d] = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
            bcoef_data[d] = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
            gcoef_data[d] = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
            static const bool homogeneous_bc = true;
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[d]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[d]->setBcCoefs(acoef_data[d], bcoef_data[d], gcoef_data[d], NULL, *patch, trimmed_bdry_box,
                                    data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data[d]->fillAll(0.0);
        }

        for (int d = 0; d < depth; d = d + 2)
        {
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool bdry_lower_side = (location_index % 2) == 0;
            const bool bdry_upper_side = (location_index % 2) != 0;

            // Modify the diagonal and off-diagonal entries to account for
            // homogeneous boundary conditions.
            //
            // Here, we follow the same linear extrapolation approach
            // implemented in class CartesianRobinBcHelper.  Namely, with u_i
            // denoting the interior cell, u_o denoting the ghost cell, and u_b
            // and u_n denoting the value and normal derivative of u at the
            // boundary,
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
            for (CellIterator bc = CellGeometry::begin(bc_coef_box); bc != CellGeometry::end(bc_coef_box); ++bc)
            {
                const CellIndex& i_s_bdry = *bc;
                const double& areal = (*acoef_data[d])(i_s_bdry, 0);
                const double& breal = (*bcoef_data[d])(i_s_bdry, 0);
                const double& aimag = (*acoef_data[d + 1])(i_s_bdry, 0);
                const double& bimag = (*bcoef_data[d + 1])(i_s_bdry, 0);

                const double& h = dx[bdry_normal_axis];

                // i_s_bdry: side index located on physical boundary
                //
                // i_c_intr: cell index located adjacent to physical boundary
                // in the patch interior
                CellIndex i_c_intr = i_s_bdry;
                if (bdry_upper_side)
                {
                    i_c_intr(bdry_normal_axis) -= 1;
                }

                if (bdry_lower_side)
                {
                    const SideIndex ilower(i_c_intr, bdry_normal_axis, SideIndex::Lower);

                    diagonal(i_c_intr, 2 * d) +=
                        off_diagonal(ilower, 2 * d) * (-(areal * h - 2.0 * breal) / (areal * h + 2.0 * breal));
                    off_diagonal(ilower, 2 * d) = 0.0;

                    diagonal(i_c_intr, 2 * d + 1) +=
                        off_diagonal(ilower, 2 * d + 1) * (-(aimag * h - 2.0 * bimag) / (aimag * h + 2.0 * bimag));
                    off_diagonal(ilower, 2 * d + 1) = 0.0;

                    diagonal(i_c_intr, 2 * d + 2) +=
                        off_diagonal(ilower, 2 * d + 2) * (-(areal * h - 2.0 * breal) / (areal * h + 2.0 * breal));
                    off_diagonal(ilower, 2 * d + 2) = 0.0;

                    diagonal(i_c_intr, 2 * d + 3) +=
                        off_diagonal(ilower, 2 * d + 3) * (-(aimag * h - 2.0 * bimag) / (aimag * h + 2.0 * bimag));
                    off_diagonal(ilower, 2 * d + 3) = 0.0;
                }

                if (bdry_upper_side)
                {
                    const SideIndex iupper(i_c_intr, bdry_normal_axis, SideIndex::Upper);

                    diagonal(i_c_intr, 2 * d) +=
                        off_diagonal(iupper, 2 * d) * (-(areal * h - 2.0 * breal) / (areal * h + 2.0 * breal));
                    off_diagonal(iupper, 2 * d) = 0.0;

                    diagonal(i_c_intr, 2 * d + 1) +=
                        off_diagonal(iupper, 2 * d + 1) * (-(aimag * h - 2.0 * bimag) / (aimag * h + 2.0 * bimag));
                    off_diagonal(iupper, 2 * d + 1) = 0.0;

                    diagonal(i_c_intr, 2 * d + 2) +=
                        off_diagonal(iupper, 2 * d + 2) * (-(areal * h - 2.0 * breal) / (areal * h + 2.0 * breal));
                    off_diagonal(iupper, 2 * d + 2) = 0.0;

                    diagonal(i_c_intr, 2 * d + 3) +=
                        off_diagonal(iupper, 2 * d + 3) * (-(aimag * h - 2.0 * bimag) / (aimag * h + 2.0 * bimag));
                    off_diagonal(iupper, 2 * d + 3) = 0.0;
                }
            }
        }
    }

    // Setup the matrix coefficients.
    for (int d = 0; d < depth; d = d + 2)
    {
        const unsigned int offset = static_cast<unsigned int>(d * stencil_sz * 2);
        for (CellIterator b = CellGeometry::begin(patch_box); b != CellGeometry::end(patch_box); ++b)
        {
            const CellIndex& i = *b;
            matrix_coefficients(i, offset + stencil_index_diag) = diagonal(i, 2 * d);
            matrix_coefficients(i, offset + stencil_index_diag + stencil_sz) = diagonal(i, 2 * d + 1);
            matrix_coefficients(i, offset + stencil_index_diag + 2 * stencil_sz) = diagonal(i, 2 * d + 2);
            matrix_coefficients(i, offset + stencil_index_diag + 3 * stencil_sz) = diagonal(i, 2 * d + 3);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const SideIndex ilower(i, axis, SideIndex::Lower);
                matrix_coefficients(i, offset + stencil_index_lower[axis]) = off_diagonal(ilower, 2 * d);
                matrix_coefficients(i, offset + stencil_index_lower[axis] + stencil_sz) =
                    off_diagonal(ilower, 2 * d + 1);
                matrix_coefficients(i, offset + stencil_index_lower[axis] + 2 * stencil_sz) =
                    off_diagonal(ilower, 2 * d + 2);
                matrix_coefficients(i, offset + stencil_index_lower[axis] + 3 * stencil_sz) =
                    off_diagonal(ilower, 2 * d + 3);

                const SideIndex iupper(i, axis, SideIndex::Upper);
                matrix_coefficients(i, offset + stencil_index_upper[axis]) = off_diagonal(iupper, 2 * d);
                matrix_coefficients(i, offset + stencil_index_upper[axis] + stencil_sz) =
                    off_diagonal(iupper, 2 * d + 1);
                matrix_coefficients(i, offset + stencil_index_upper[axis] + 2 * stencil_sz) =
                    off_diagonal(iupper, 2 * d + 2);
                matrix_coefficients(i, offset + stencil_index_upper[axis] + 3 * stencil_sz) =
                    off_diagonal(iupper, 2 * d + 3);
            }
        }
    }
    return;
}

void PoissonUtilities::computeSCMatrixCoefficients(boost::shared_ptr<Patch> patch,
                                                   SideData<double>& matrix_coefficients,
                                                   const std::vector<Index>& stencil,
                                                   const PoissonSpecifications& poisson_spec,
                                                   const std::vector<RobinBcCoefStrategy*>& bc_coefs,
                                                   double data_time)
{
    TBOX_ASSERT(bc_coefs.size() == NDIM);
    const int stencil_sz = static_cast<int>(stencil.size());
    TBOX_ASSERT(stencil_sz == 2 * NDIM + 1);
    std::map<Index, int> stencil_map;
    for (int k = 0; k < stencil_sz; ++k)
    {
        stencil_map[stencil[k]] = k;
    }
    const int stencil_index_diag = stencil_map[Index::getZeroIndex(DIM)];
    boost::array<int, NDIM> stencil_index_lower, stencil_index_upper;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        Index ilower(DIM, 0), iupper(DIM, 0);
        ilower(axis) = -1;
        iupper(axis) = +1;
        stencil_index_lower[axis] = stencil_map[ilower];
        stencil_index_upper[axis] = stencil_map[iupper];
    }
    TBOX_ASSERT(matrix_coefficients.getDepth() == stencil_sz);
    if (!(poisson_spec.cIsZero() || poisson_spec.cIsConstant()) || !poisson_spec.dIsConstant())
    {
        TBOX_ERROR(
            "PoissonUtilities::computeSCMatrixCoefficients() does not support non-constant "
            "coefficient problems\n");
    }
    const double C = (poisson_spec.cIsZero() ? 0.0 : poisson_spec.getCConstant());
    const double D = poisson_spec.getDConstant();

    const Box& patch_box = patch->getBox();
    const std::vector<BoundaryBox> physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    const int n_physical_codim1_boxes = static_cast<int>(physical_codim1_boxes.size());
    auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch->getPatchGeometry());
    const double* const dx = pgeom->getDx();
    const double* const patch_x_lower = pgeom->getXLower();
    const double* const patch_x_upper = pgeom->getXUpper();
    const IntVector& ratio_to_level_zero = pgeom->getRatio();
    CartesianPatchGeometry::TwoDimBool touches_regular_bdry(DIM);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int side = 0; side < 2; ++side)
        {
            touches_regular_bdry(axis, side) = pgeom->getTouchesRegularBoundary(axis, side);
        }
    }

    // Compute all matrix coefficients, including those on the physical
    // boundary; however, do not yet take physical boundary conditions into
    // account.  Boundary conditions are handled subsequently.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        std::vector<double> mat_vals(stencil_sz, 0.0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const double dx_sq = dx[d] * dx[d];
            mat_vals[stencil_index_lower[d]] += D / dx_sq; // lower off-diagonal
            mat_vals[stencil_index_upper[d]] += D / dx_sq; // upper off-diagonal
            mat_vals[stencil_index_diag] -= D / dx_sq;     // diagonal
            mat_vals[stencil_index_diag] -= D / dx_sq;     // diagonal
        }
        mat_vals[stencil_index_diag] += C; // diagonal

        for (int stencil_index = 0; stencil_index < stencil_sz; ++stencil_index)
        {
            matrix_coefficients.fill(mat_vals[stencil_index], stencil_index);
        }
    }

    // Modify matrix coefficients to account for physical boundary conditions
    // along boundaries which ARE NOT aligned with the data axis.
    //
    // NOTE: It important to set these values first to avoid problems at corners
    // in the physical domain.  In particular, since Dirichlet boundary
    // conditions for values located on the physical boundary override all other
    // boundary conditions, we set those values last.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const BoundaryBox& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = location_index % 2 == 0;

            if (bdry_normal_axis == axis) continue;

            const Box bc_fill_box =
                pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector::getOne(DIM));
            const BoundaryBox trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box bc_coef_box = compute_tangential_extension(
                PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), axis);

            auto acoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
            auto bcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
            auto gcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);

            // Temporarily reset the patch geometry object associated with the
            // patch so that boundary conditions are set at the correct spatial
            // locations.
            boost::array<double, NDIM> shifted_patch_x_lower, shifted_patch_x_upper;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                shifted_patch_x_lower[d] = patch_x_lower[d];
                shifted_patch_x_upper[d] = patch_x_upper[d];
            }
            shifted_patch_x_lower[axis] -= 0.5 * dx[axis];
            shifted_patch_x_upper[axis] -= 0.5 * dx[axis];
            patch->setPatchGeometry(
                boost::make_shared<CartesianPatchGeometry>(ratio_to_level_zero, touches_regular_bdry, dx,
                                                           shifted_patch_x_lower.data(), shifted_patch_x_upper.data()));

            // Set the boundary condition coefficients.
            static const bool homogeneous_bc = true;
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            // Restore the original patch geometry object.
            patch->setPatchGeometry(pgeom);

            // Modify the matrix coefficients to account for homogeneous
            // boundary conditions.
            //
            // Here, we follow the same linear extrapolation approach
            // implemented in class CartesianRobinBcHelper.  Namely, with u_i
            // denoting the interior cell, u_o denoting the ghost cell, and u_b
            // and u_n denoting the value and normal derivative of u at the
            // boundary,
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
            for (CellIterator bc = CellGeometry::begin(bc_coef_box); bc != CellGeometry::end(bc_coef_box); ++bc)
            {
                const Index& i = *bc;
                const double& a = (*acoef_data)(i, 0);
                const double& b = (*bcoef_data)(i, 0);
                const double& h = dx[bdry_normal_axis];

                Index i_intr = i;
                if (is_lower)
                {
                    i_intr(bdry_normal_axis) += 0;
                }
                else
                {
                    i_intr(bdry_normal_axis) -= 1;
                }
                const SideIndex i_s(i_intr, axis, SideIndex::Lower);

                if (is_lower)
                {
                    matrix_coefficients(i_s, stencil_index_diag) +=
                        matrix_coefficients(i_s, stencil_index_lower[bdry_normal_axis]) *
                        (-(a * h - 2.0 * b) / (a * h + 2.0 * b));
                    matrix_coefficients(i_s, stencil_index_lower[bdry_normal_axis]) = 0.0;
                }
                else
                {
                    matrix_coefficients(i_s, stencil_index_diag) +=
                        matrix_coefficients(i_s, stencil_index_upper[bdry_normal_axis]) *
                        (-(a * h - 2.0 * b) / (a * h + 2.0 * b));
                    matrix_coefficients(i_s, stencil_index_upper[bdry_normal_axis]) = 0.0;
                }
            }
        }
    }

    // Modify matrix coefficients to account for physical boundary
    // conditions along boundaries which ARE aligned with the data axis.
    //
    // NOTE: It important to set these values last to avoid problems at corners
    // in the physical domain.  In particular, since Dirichlet boundary
    // conditions for values located on the physical boundary override all other
    // boundary conditions, we set those values last.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const BoundaryBox& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = location_index % 2 == 0;

            if (bdry_normal_axis != axis) continue;

            const Box bc_fill_box =
                pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector::getOne(DIM));
            const BoundaryBox trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            auto acoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
            auto bcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
            auto gcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);

            // Set the boundary condition coefficients.
            static const bool homogeneous_bc = true;
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            // Modify the matrix coefficients to account for homogeneous
            // boundary conditions.
            //
            // With u_i denoting the interior cell, u_o denoting the ghost cell,
            // and u_b and u_n denoting the value and normal derivative of u at
            // the boundary,
            //
            //     u_n = (u_o - u_i)/(2*h)
            //
            // Now, if
            //
            //     a*u_b + b*u_n = 0
            //
            // then if b != 0,
            //
            //     u_o = u_i - 2*h*(a/b)*u_b
            //
            // and
            //
            //     -(D/h^2)*u_o = (D*2*(a/b)/h)*u_b - (D/h^2)*u_i
            //
            // If b == 0, then u_b = 0, which we enforce directly.
            for (CellIterator bc = CellGeometry::begin(bc_coef_box); bc != CellGeometry::end(bc_coef_box); ++bc)
            {
                const Index& i = *bc;
                const SideIndex i_s(i, axis, SideIndex::Lower);
                const double& a = (*acoef_data)(i, 0);
                const double& b = (*bcoef_data)(i, 0);
                if (b == 0.0)
                {
                    for (int k = 0; k < stencil_sz; ++k)
                    {
                        matrix_coefficients(i_s, k) = 0.0;
                    }
                    matrix_coefficients(i_s, stencil_index_diag) = 1.0;
                }
                else
                {
                    TBOX_ASSERT(!MathUtilities<double>::equalEps(b, 0.0));
                    if (is_lower)
                    {
                        matrix_coefficients(i_s, stencil_index_diag) -=
                            matrix_coefficients(i_s, stencil_index_lower[bdry_normal_axis]) * 2 * dx[bdry_normal_axis] *
                            a / b;
                        matrix_coefficients(i_s, stencil_index_upper[bdry_normal_axis]) +=
                            matrix_coefficients(i_s, stencil_index_lower[bdry_normal_axis]);
                        matrix_coefficients(i_s, stencil_index_lower[bdry_normal_axis]) = 0.0;
                    }
                    else
                    {
                        matrix_coefficients(i_s, stencil_index_diag) -=
                            matrix_coefficients(i_s, stencil_index_upper[bdry_normal_axis]) * 2 * dx[bdry_normal_axis] *
                            a / b;
                        matrix_coefficients(i_s, stencil_index_lower[bdry_normal_axis]) +=
                            matrix_coefficients(i_s, stencil_index_upper[bdry_normal_axis]);
                        matrix_coefficients(i_s, stencil_index_upper[bdry_normal_axis]) = 0.0;
                    }
                }
            }
        }
    }
    return;
}

void PoissonUtilities::adjustCCBoundaryRhsEntries(boost::shared_ptr<Patch> patch,
                                                  CellData<double>& rhs_data,
                                                  const PoissonSpecifications& poisson_spec,
                                                  RobinBcCoefStrategy* bc_coef,
                                                  double data_time,
                                                  bool homogeneous_bc)
{
    adjustCCBoundaryRhsEntries(patch, rhs_data, poisson_spec, std::vector<RobinBcCoefStrategy*>(1, bc_coef), data_time,
                               homogeneous_bc);
    return;
}

void PoissonUtilities::adjustCCBoundaryRhsEntries(boost::shared_ptr<Patch> patch,
                                                  CellData<double>& rhs_data,
                                                  const PoissonSpecifications& poisson_spec,
                                                  const std::vector<RobinBcCoefStrategy*>& bc_coefs,
                                                  double data_time,
                                                  bool homogeneous_bc)
{
    const int depth = static_cast<int>(bc_coefs.size());
    TBOX_ASSERT(rhs_data.getDepth() == depth);

    const Box& patch_box = patch->getBox();
    OutersideData<double> D_data(patch_box, depth);
    if (!poisson_spec.dIsConstant())
    {
        D_data.copy(*patch->getPatchData(poisson_spec.getDPatchDataId()));
    }
    else
    {
        D_data.fillAll(poisson_spec.getDConstant());
    }

    auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch->getPatchGeometry());
    const double* const dx = pgeom->getDx();
    const std::vector<BoundaryBox> codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);

    // Modify the rhs entries to account for inhomogeneous boundary conditions.
    const int n_bdry_boxes = static_cast<int>(codim1_boxes.size());
    for (int n = 0; n < n_bdry_boxes; ++n)
    {
        const BoundaryBox& bdry_box = codim1_boxes[n];
        const BoundaryBox trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
        const Box bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

        auto acoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
        auto bcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
        auto gcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);

        for (int d = 0; d < depth; ++d)
        {
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[d]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[d]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool bdry_upper_side = (location_index % 2) != 0;
            const int bdry_side = (bdry_upper_side ? 1 : 0);

            // Here, we follow the same linear extrapolation approach
            // implemented in class CartesianRobinBcHelper.  Namely, with u_i
            // denoting the interior cell, u_o denoting the ghost cell, and u_b
            // and u_n denoting the value and normal derivative of u at the
            // boundary,
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
            //
            // NOTE: i_s_bdry: side index located on physical boundary
            //       i_c_intr: cell index located adjacent to physical boundary
            //                 in the patch interior
            for (CellIterator bc = CellGeometry::begin(bc_coef_box); bc != CellGeometry::end(bc_coef_box); ++bc)
            {
                const CellIndex& i_s_bdry = *bc;
                const double& a = (*acoef_data)(i_s_bdry, 0);
                const double& b = (*bcoef_data)(i_s_bdry, 0);
                const double& g = (*gcoef_data)(i_s_bdry, 0);
                const double& h = dx[bdry_normal_axis];
                CellIndex i_c_intr = i_s_bdry;
                if (bdry_upper_side)
                {
                    i_c_intr(bdry_normal_axis) -= 1;
                }
                rhs_data(i_c_intr, d) += (D_data.getArrayData(bdry_normal_axis, bdry_side)(i_s_bdry, d) / h) *
                                         (-2.0 * g) / (2.0 * b + h * a);
            }
        }
    }
    return;
}

void PoissonUtilities::adjustCCComplexBoundaryRhsEntries(boost::shared_ptr<Patch> patch,
                                                         CellData<double>& rhs_data,
                                                         const PoissonSpecifications& poisson_spec_real,
                                                         const PoissonSpecifications& poisson_spec_imag,
                                                         RobinBcCoefStrategy* bc_coef,
                                                         double data_time,
                                                         bool homogeneous_bc)
{
    adjustCCComplexBoundaryRhsEntries(patch, rhs_data, poisson_spec_real, poisson_spec_imag,
                                      std::vector<RobinBcCoefStrategy*>(2, bc_coef), data_time, homogeneous_bc);
    return;
}

void PoissonUtilities::adjustCCComplexBoundaryRhsEntries(boost::shared_ptr<Patch> patch,
                                                         CellData<double>& rhs_data,
                                                         const PoissonSpecifications& poisson_spec_real,
                                                         const PoissonSpecifications& poisson_spec_imag,
                                                         const std::vector<RobinBcCoefStrategy*>& bc_coefs,
                                                         double data_time,
                                                         bool homogeneous_bc)
{
    const int depth = static_cast<int>(bc_coefs.size());
    TBOX_ASSERT(depth >= 2 && depth % 2 == 0);
    TBOX_ASSERT(rhs_data.getDepth() == depth);

    const Box& patch_box = patch->getBox();
    OutersideData<double> D_data_real(patch_box, 1);
    OutersideData<double> D_data_imag(patch_box, 1);
    if (!poisson_spec_real.dIsConstant())
    {
        D_data_real.copy(*patch->getPatchData(poisson_spec_real.getDPatchDataId()));
    }
    else
    {
        D_data_real.fillAll(poisson_spec_real.getDConstant());
    }

    if (!poisson_spec_imag.dIsConstant())
    {
        D_data_imag.copy(*patch->getPatchData(poisson_spec_imag.getDPatchDataId()));
    }
    else
    {
        D_data_imag.fillAll(poisson_spec_imag.getDConstant());
    }

    auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch->getPatchGeometry());
    const double* const dx = pgeom->getDx();
    const std::vector<BoundaryBox> codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);

    // Modify the rhs entries to account for inhomogeneous boundary conditions.
    const int n_bdry_boxes = static_cast<int>(codim1_boxes.size());
    for (int n = 0; n < n_bdry_boxes; ++n)
    {
        const BoundaryBox& bdry_box = codim1_boxes[n];
        const BoundaryBox trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
        const Box bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

        std::vector<boost::shared_ptr<ArrayData<double> > > acoef_data(depth), bcoef_data(depth), gcoef_data(depth);
        for (int d = 0; d < depth; ++d)
        {
            acoef_data[d] = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
            bcoef_data[d] = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
            gcoef_data[d] = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[d]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[d]->setBcCoefs(acoef_data[d], bcoef_data[d], gcoef_data[d], NULL, *patch, trimmed_bdry_box,
                                    data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data[d]->fillAll(0.0);
        }

        for (int d = 0; d < depth; d = d + 2)
        {
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool bdry_upper_side = (location_index % 2) != 0;
            const int bdry_side = (bdry_upper_side ? 1 : 0);

            // Here, we follow the same linear extrapolation approach
            // implemented in class CartesianRobinBcHelper.  Namely, with u_i
            // denoting the interior cell, u_o denoting the ghost cell, and u_b
            // and u_n denoting the value and normal derivative of u at the
            // boundary,
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
            //
            // NOTE: i_s_bdry: side index located on physical boundary
            //       i_c_intr: cell index located adjacent to physical boundary
            //                 in the patch interior
            for (CellIterator b = CellGeometry::begin(bc_coef_box); b != CellGeometry::end(bc_coef_box); ++b)
            {
                const CellIndex& i_s_bdry = *b;

                const double& areal = (*acoef_data[d])(i_s_bdry, 0);
                const double& breal = (*bcoef_data[d])(i_s_bdry, 0);
                const double& greal = (*gcoef_data[d])(i_s_bdry, 0);
                const double& aimag = (*acoef_data[d + 1])(i_s_bdry, 0);
                const double& bimag = (*bcoef_data[d + 1])(i_s_bdry, 0);
                const double& gimag = (*gcoef_data[d + 1])(i_s_bdry, 0);

                const double& h = dx[bdry_normal_axis];
                CellIndex i_c_intr = i_s_bdry;
                if (bdry_upper_side)
                {
                    i_c_intr(bdry_normal_axis) -= 1;
                }
                rhs_data(i_c_intr, d) += (D_data_real.getArrayData(bdry_normal_axis, bdry_side)(i_s_bdry, 0) / h) *
                                             (-2.0 * greal) / (2.0 * breal + h * areal) -
                                         (D_data_imag.getArrayData(bdry_normal_axis, bdry_side)(i_s_bdry, 0) / h) *
                                             (-2.0 * gimag) / (2.0 * bimag + h * aimag);

                rhs_data(i_c_intr, d + 1) += (D_data_imag.getArrayData(bdry_normal_axis, bdry_side)(i_s_bdry, 0) / h) *
                                                 (-2.0 * greal) / (2.0 * breal + h * areal) +
                                             (D_data_real.getArrayData(bdry_normal_axis, bdry_side)(i_s_bdry, 0) / h) *
                                                 (-2.0 * gimag) / (2.0 * bimag + h * aimag);
            }
        }
    }
    return;
}

void PoissonUtilities::adjustSCBoundaryRhsEntries(boost::shared_ptr<Patch> patch,
                                                  SideData<double>& rhs_data,
                                                  const PoissonSpecifications& poisson_spec,
                                                  const std::vector<RobinBcCoefStrategy*>& bc_coefs,
                                                  double data_time,
                                                  bool homogeneous_bc)
{
    TBOX_ASSERT(bc_coefs.size() == NDIM);
    if (!(poisson_spec.cIsZero() || poisson_spec.cIsConstant()) || !poisson_spec.dIsConstant())
    {
        TBOX_ERROR(
            "PoissonUtilities::adjustSCBoundaryRhsEntries() does not support non-constant "
            "coefficient problems\n");
    }
    //  const double C = (poisson_spec.cIsZero() ? 0.0 : poisson_spec.getCConstant());
    const double D = poisson_spec.getDConstant();

    const Box& patch_box = patch->getBox();
    const std::vector<BoundaryBox> physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    const int n_physical_codim1_boxes = static_cast<int>(physical_codim1_boxes.size());
    auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch->getPatchGeometry());
    const double* const dx = pgeom->getDx();
    const double* const patch_x_lower = pgeom->getXLower();
    const double* const patch_x_upper = pgeom->getXUpper();
    const IntVector& ratio_to_level_zero = pgeom->getRatio();
    PatchGeometry::TwoDimBool touches_regular_bdry(DIM);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int side = 0; side < 2; ++side)
        {
            touches_regular_bdry(axis, side) = pgeom->getTouchesRegularBoundary(axis, side);
        }
    }

    // Modify the rhs entries to account for inhomogeneous boundary conditions
    // along boundaries which ARE NOT aligned with the data axis.
    //
    // NOTE: It important to set these values first to avoid problems at corners
    // in the physical domain.  In particular, since Dirichlet boundary
    // conditions for values located on the physical boundary override all other
    // boundary conditions, we set those values last.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const BoundaryBox& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = location_index % 2 == 0;

            if (bdry_normal_axis == axis) continue;

            const Box bc_fill_box =
                pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector::getOne(DIM));
            const BoundaryBox trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box bc_coef_box = compute_tangential_extension(
                PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), axis);

            auto acoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
            auto bcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
            auto gcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);

            // Temporarily reset the patch geometry object associated with the
            // patch so that boundary conditions are set at the correct spatial
            // locations.
            boost::array<double, NDIM> shifted_patch_x_lower, shifted_patch_x_upper;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                shifted_patch_x_lower[d] = patch_x_lower[d];
                shifted_patch_x_upper[d] = patch_x_upper[d];
            }
            shifted_patch_x_lower[axis] -= 0.5 * dx[axis];
            shifted_patch_x_upper[axis] -= 0.5 * dx[axis];
            patch->setPatchGeometry(
                boost::make_shared<CartesianPatchGeometry>(ratio_to_level_zero, touches_regular_bdry, dx,
                                                           shifted_patch_x_lower.data(), shifted_patch_x_upper.data()));

            // Set the boundary condition coefficients.
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            // Restore the original patch geometry object.
            patch->setPatchGeometry(pgeom);

            // Here, we follow the same linear extrapolation approach
            // implemented in class CartesianRobinBcHelper.  Namely, with u_i
            // denoting the interior cell, u_o denoting the ghost cell, and u_b
            // and u_n denoting the value and normal derivative of u at the
            // boundary,
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
            for (CellIterator bc = CellGeometry::begin(bc_coef_box); bc != CellGeometry::end(bc_coef_box); ++bc)
            {
                const Index& i = *bc;
                const double& a = (*acoef_data)(i, 0);
                const double& b = (*bcoef_data)(i, 0);
                const double& g = (*gcoef_data)(i, 0);
                const double& h = dx[bdry_normal_axis];
                Index i_intr = i;
                if (is_lower)
                {
                    i_intr(bdry_normal_axis) += 0;
                }
                else
                {
                    i_intr(bdry_normal_axis) -= 1;
                }
                const SideIndex i_s(i_intr, axis, SideIndex::Lower);
                rhs_data(i_s) += (D / h) * (-2.0 * g) / (2.0 * b + h * a);
            }
        }
    }

    // Modify the rhs entries to account for inhomogeneous boundary conditions
    // along boundaries which ARE aligned with the data axis.
    //
    // NOTE: It important to set these values last to avoid problems at corners
    // in the physical domain.  In particular, since Dirichlet boundary
    // conditions for values located on the physical boundary override all other
    // boundary conditions, we set those values last.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const BoundaryBox& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;

            if (bdry_normal_axis != axis) continue;

            const Box bc_fill_box =
                pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector::getOne(DIM));
            const BoundaryBox trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            auto acoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
            auto bcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
            auto gcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);

            // Set the boundary condition coefficients.
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            // For the non-symmetric boundary treatment,
            //
            //     u_n = (u_o - u_i)/(2*h)
            //
            // Now, if
            //
            //     a*u_b + b*u_n = g
            //
            // then, with u_b = u_i = 0,
            //
            //     u_o = 2*h*g/b
            //
            // so that the boundary flux is
            //
            //     (u_b - u_o)/h = -2*g/b
            //
            // In this loop, we modify the rhs entries appropriately.
            //
            // NOTE: At Dirichlet boundaries, boundary values are provided by
            // the right-hand side vector.
            for (CellIterator bc = CellGeometry::begin(bc_coef_box); bc != CellGeometry::end(bc_coef_box); ++bc)
            {
                const Index& i = *bc;
                //              const double& a = (*acoef_data)(i,0);
                const double& b = (*bcoef_data)(i, 0);
                const double& g = (*gcoef_data)(i, 0);
                const double& h = dx[bdry_normal_axis];
                const SideIndex i_s_bdry(i, bdry_normal_axis, SideIndex::Lower);
                if (b != 0.0)
                {
                    TBOX_ASSERT(!MathUtilities<double>::equalEps(b, 0.0));
                    rhs_data(i_s_bdry) += (D / h) * (-2.0 * g) / b;
                }
            }
        }
    }
    return;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
