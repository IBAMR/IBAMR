// Filename: CartExtrapPhysBdryOp.cpp
// Created on 30 Sep 2006 by Boyce Griffith
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

#include <stdlib.h>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/FaceIterator.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/pdat/NodeIterator.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/pdat/SideIterator.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "boost/array.hpp"
#include "ibtk/CartExtrapPhysBdryOp.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "SAMRAI/tbox/Array.h"

#include "SAMRAI/tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int REFINE_OP_STENCIL_WIDTH = 0;

template <typename D, typename I>
inline double
compute_linear_extrap(D& patch_data, const I& i, const I& i_intr, const IntVector& i_shft, const int depth)
{
    double ret_val = patch_data(i_intr, depth);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (i_shft(d) != 0)
        {
            const I& i_intr0 = i_intr;
            I i_intr1 = i_intr;
            i_intr1(d) += i_shft(d);

            const double& f0 = patch_data(i_intr0, depth);
            const double& f1 = patch_data(i_intr1, depth);

            const double du = f0 - f1;
            const double delta = std::abs(i(d) - i_intr(d));

            ret_val += du * delta;
        }
    }
    return ret_val;
}

template <typename D, typename I>
inline double compute_quadratic_extrap(D& patch_data,
                                       const I& i,
                                       const I& i_intr,
                                       const IntVector& i_shft,
                                       const int depth,
                                       const int codim)
{
    if (codim == 1)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (i_shft(d) != 0)
            {
#if 1
                const I& i_intr0 = i_intr;

                I i_intr1 = i_intr;
                i_intr1(d) += i_shft(d);

                I i_intr2 = i_intr1;
                i_intr2(d) += i_shft(d);

                const double& f0 = patch_data(i_intr0, depth);
                const double& f1 = patch_data(i_intr1, depth);
                const double& f2 = patch_data(i_intr2, depth);

                const double x = std::abs(i(d) - i_intr(d));

                return (1.0 / 2.0 * f2 - f1 + 1.0 / 2.0 * f0) * x * x +
                       (-1.0 / 2.0 * f2 + 2.0 * f1 - 3.0 / 2.0 * f0) * x + f0;
#endif

#if 0
                // NOTE: The following only works in general for the case that
                // the ghost cell width is >= 3.
                const I& i_intr0 = i_intr;

                I i_intr2 = i_intr;
                i_intr2(d) += 2*i_shft(d);

                I i_intr3 = i_intr2;
                i_intr3(d) += i_shft(d);

                const double& f0 = patch_data(i_intr0,depth);
                const double& f2 = patch_data(i_intr2,depth);
                const double& f3 = patch_data(i_intr3,depth);

                const double x = std::abs(i(d)-i_intr(d));

                return (1.0/3.0*f3-1.0/2.0*f2+1.0/6.0*f0)*x*x+(-2.0/3.0*f3+3.0/2.0*f2-5.0/6.0*f0)*x+f0;
#endif
            }
        }
    }
    else
    {
        return compute_linear_extrap(patch_data, i, i_intr, i_shft, depth);
    }
    return 0.0; // this statement should not be reached
}
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartExtrapPhysBdryOp::CartExtrapPhysBdryOp() : d_patch_data_indices(), d_extrap_type("NULL")
{
    // intentionally blank
    return;
}

CartExtrapPhysBdryOp::CartExtrapPhysBdryOp(const int patch_data_index, const std::string& extrap_type)
    : d_patch_data_indices(), d_extrap_type("NULL")
{
    setPatchDataIndex(patch_data_index);
    setExtrapolationType(extrap_type);
    return;
}

CartExtrapPhysBdryOp::CartExtrapPhysBdryOp(const std::set<int>& patch_data_indices, const std::string& extrap_type)
    : d_patch_data_indices(), d_extrap_type("NULL")
{
    setPatchDataIndices(patch_data_indices);
    setExtrapolationType(extrap_type);
    return;
}

CartExtrapPhysBdryOp::CartExtrapPhysBdryOp(const ComponentSelector& patch_data_indices, const std::string& extrap_type)
    : d_patch_data_indices(), d_extrap_type("NULL")
{
    setPatchDataIndices(patch_data_indices);
    setExtrapolationType(extrap_type);
    return;
}

CartExtrapPhysBdryOp::~CartExtrapPhysBdryOp()
{
    // intentionally blank
    return;
}

void CartExtrapPhysBdryOp::setPatchDataIndex(const int patch_data_index)
{
    std::set<int> patch_data_indices;
    patch_data_indices.insert(patch_data_index);
    setPatchDataIndices(patch_data_indices);
    return;
}

void CartExtrapPhysBdryOp::setPatchDataIndices(const std::set<int>& patch_data_indices)
{
    d_patch_data_indices.clear();
    d_patch_data_indices = patch_data_indices;
    return;
}

void CartExtrapPhysBdryOp::setPatchDataIndices(const ComponentSelector& patch_data_indices)
{
    std::set<int> patch_data_index_set;
    for (auto l = 0; l < patch_data_indices.getSize(); ++l)
    {
        if (patch_data_indices.isSet(l))
        {
            const int patch_data_index = l;
            patch_data_index_set.insert(patch_data_index);
        }
    }
    setPatchDataIndices(patch_data_index_set);
    return;
}

void CartExtrapPhysBdryOp::setExtrapolationType(const std::string& extrap_type)
{
    // Ensure that the extrapolation type is supported by this class.
    if (extrap_type != "CONSTANT" && extrap_type != "LINEAR" && extrap_type != "QUADRATIC")
    {
        TBOX_ERROR("CartExtrapPhysBdryOp::setExtrapolationType():\n"
                   << "  unknown extrapolation type: " << extrap_type << "\n"
                   << "  valid selections are: CONSTANT, LINEAR, or QUADRATIC" << std::endl);
    }

    if (extrap_type == "QUADRATIC")
    {
        IBTK_DO_ONCE(TBOX_WARNING("CartExtrapPhysBdryOp::setExtrapolationType():\n"
                                  << "  extrapolation type " << extrap_type
                                  << " generally requires large ghost cell widths" << std::endl););
    }

    d_extrap_type = extrap_type;
    return;
}

void CartExtrapPhysBdryOp::setPhysicalBoundaryConditions(Patch& patch,
                                                         const double /*fill_time*/,
                                                         const IntVector& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector::getZero(DIM)) return;

    auto pgeom = patch.getPatchGeometry();
    const Box& patch_box = patch.getBox();

    std::vector<std::pair<Box, std::pair<int, int> > > bdry_fill_boxes;

#if (NDIM > 1)
#if (NDIM > 2)
    // Compute the co-dimension three boundary fill boxes.
    const std::vector<BoundaryBox> physical_codim3_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim3Boxes(patch);
    const int n_physical_codim3_boxes = physical_codim3_boxes.size();
    for (int n = 0; n < n_physical_codim3_boxes; ++n)
    {
        const BoundaryBox& bdry_box = physical_codim3_boxes[n];
        const Box bdry_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
        const unsigned int location_index = bdry_box.getLocationIndex();
        const int codim = 3;
        bdry_fill_boxes.push_back(std::make_pair(bdry_fill_box, std::make_pair(location_index, codim)));
    }
#endif
    // Compute the co-dimension two boundary fill boxes.
    const std::vector<BoundaryBox> physical_codim2_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(patch);
    const auto n_physical_codim2_boxes = physical_codim2_boxes.size();
    for (auto n = 0; n < n_physical_codim2_boxes; ++n)
    {
        const BoundaryBox& bdry_box = physical_codim2_boxes[n];
        const Box bdry_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
        const unsigned int location_index = bdry_box.getLocationIndex();
        const int codim = 2;
        bdry_fill_boxes.push_back(std::make_pair(bdry_fill_box, std::make_pair(location_index, codim)));
    }
#endif
    // Compute the co-dimension one boundary fill boxes.
    const std::vector<BoundaryBox> physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    const auto n_physical_codim1_boxes = physical_codim1_boxes.size();
    for (auto n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox& bdry_box = physical_codim1_boxes[n];
        const Box bdry_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
        const unsigned int location_index = bdry_box.getLocationIndex();
        const int codim = 1;
        bdry_fill_boxes.push_back(std::make_pair(bdry_fill_box, std::make_pair(location_index, codim)));
    }

    // Set the boundary values.
    setPhysicalBoundaryConditions_cell(patch, bdry_fill_boxes);
    setPhysicalBoundaryConditions_face(patch, bdry_fill_boxes);
    setPhysicalBoundaryConditions_node(patch, bdry_fill_boxes);
    setPhysicalBoundaryConditions_side(patch, bdry_fill_boxes);
    return;
}

IntVector CartExtrapPhysBdryOp::getRefineOpStencilWidth(const Dimension& dim) const
{
    return IntVector(dim, REFINE_OP_STENCIL_WIDTH);
}

void CartExtrapPhysBdryOp::preprocessRefine(Patch& /*fine*/,
                                            const Patch& /*coarse*/,
                                            const Box& /*fine_box*/,
                                            const IntVector& /*ratio*/)
{
    // intentionally blank
    return;
}

void CartExtrapPhysBdryOp::postprocessRefine(Patch& /*fine*/,
                                             const Patch& /*coarse*/,
                                             const Box& /*fine_box*/,
                                             const IntVector& /*ratio*/)
{
    // intentionally blank
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void CartExtrapPhysBdryOp::setPhysicalBoundaryConditions_cell(
    Patch& patch,
    const std::vector<std::pair<Box, std::pair<int, int> > >& bdry_fill_boxes)
{
    const Box& patch_box = patch.getBox();
    const Index& patch_lower = patch_box.lower();
    const Index& patch_upper = patch_box.upper();

    const int extrap_type =
        (d_extrap_type == "CONSTANT" ? 0 : (d_extrap_type == "LINEAR" ? 1 : (d_extrap_type == "QUADRATIC" ? 2 : -1)));
    if (extrap_type != 0 && extrap_type != 1 && extrap_type != 2)
    {
        TBOX_ERROR("CartExtrapPhysBdryOp::setPhysicalBoundaryConditions():\n"
                   << "  unknown extrapolation type: " << d_extrap_type << "\n"
                   << "  valid selections are: CONSTANT, LINEAR, or QUADRATIC" << std::endl);
    }

    // Set the physical boundary conditions for the specified patch data
    // indices.
    for (auto cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        auto var_db = VariableDatabase::getDatabase();
        boost::shared_ptr<Variable> var;
        var_db->mapIndexToVariable(patch_data_idx, var);
        auto cc_var = boost::dynamic_pointer_cast<CellVariable<double> >(var);
        if (!cc_var) continue;

        auto patch_data = BOOST_CAST<CellData<double> >(patch.getPatchData(patch_data_idx));
        const Box& ghost_box = patch_data->getGhostBox();

        // Loop over the boundary fill boxes and extrapolate the data.
        for (auto it = bdry_fill_boxes.begin(); it != bdry_fill_boxes.end(); ++it)
        {
            const Box& bdry_fill_box = it->first;
            const unsigned int location_index = it->second.first;
            const int codim = it->second.second;
#if (NDIM == 2)
            const boost::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isLower(location_index, codim,
                                                                                             1) } };
            const boost::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isUpper(location_index, codim,
                                                                                             1) } };
#endif
#if (NDIM == 3)
            const boost::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isLower(location_index, codim, 1),
                                                          PhysicalBoundaryUtilities::isLower(location_index, codim,
                                                                                             2) } };
            const boost::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isUpper(location_index, codim, 1),
                                                          PhysicalBoundaryUtilities::isUpper(location_index, codim,
                                                                                             2) } };
#endif
            // Loop over the boundary box indices and compute the nearest
            // interior index.
            for (int depth = 0; depth < patch_data->getDepth(); ++depth)
            {
                const Box it_box = bdry_fill_box * ghost_box;
                for (auto b = CellGeometry::begin(it_box), e = CellGeometry::end(it_box); b != e; ++b)
                {
                    const CellIndex& i = *b;
                    CellIndex i_intr = i;
                    IntVector i_shft = IntVector::getZero(DIM);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (is_lower[d])
                        {
                            i_intr(d) = patch_lower(d);
                            i_shft(d) = +1; // use interior data for extrapolation
                        }
                        else if (is_upper[d])
                        {
                            i_intr(d) = patch_upper(d);
                            i_shft(d) = -1; // use interior data for extrapolation
                        }
                    }

                    // Perform constant, linear, or quadratic extrapolation.
                    switch (extrap_type)
                    {
                    case 0:
                        (*patch_data)(i, depth) = (*patch_data)(i_intr, depth);
                        break;
                    case 1:
                        (*patch_data)(i, depth) = compute_linear_extrap(*patch_data, i, i_intr, i_shft, depth);
                        break;
                    case 2:
                        (*patch_data)(i, depth) =
                            compute_quadratic_extrap(*patch_data, i, i_intr, i_shft, depth, codim);
                        break;
                    }
                }
            }
        }
    }
    return;
}

void CartExtrapPhysBdryOp::setPhysicalBoundaryConditions_face(
    Patch& patch,
    const std::vector<std::pair<Box, std::pair<int, int> > >& bdry_fill_boxes)
{
    const Box& patch_box = patch.getBox();
    const Index& patch_lower = patch_box.lower();
    const Index& patch_upper = patch_box.upper();

    const int extrap_type =
        (d_extrap_type == "CONSTANT" ? 0 : (d_extrap_type == "LINEAR" ? 1 : (d_extrap_type == "QUADRATIC" ? 2 : -1)));
    if (extrap_type != 0 && extrap_type != 1 && extrap_type != 2)
    {
        TBOX_ERROR("CartExtrapPhysBdryOp::setPhysicalBoundaryConditions():\n"
                   << "  unknown extrapolation type: " << d_extrap_type << "\n"
                   << "  valid selections are: CONSTANT, LINEAR, or QUADRATIC" << std::endl);
    }

    // Set the physical boundary conditions for the specified patch data
    // indices.
    for (auto cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        auto var_db = VariableDatabase::getDatabase();
        boost::shared_ptr<Variable> var;
        var_db->mapIndexToVariable(patch_data_idx, var);
        auto fc_var = boost::dynamic_pointer_cast<FaceVariable<double> >(var);
        if (!fc_var) continue;

        auto patch_data = BOOST_CAST<FaceData<double> >(patch.getPatchData(patch_data_idx));
        const Box& ghost_box = patch_data->getGhostBox();

        // Loop over the boundary fill boxes and extrapolate the data.
        for (auto it = bdry_fill_boxes.begin(); it != bdry_fill_boxes.end(); ++it)
        {
            const Box& bdry_fill_box = it->first;
            const unsigned int location_index = it->second.first;
            const int codim = it->second.second;
#if (NDIM == 2)
            const boost::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isLower(location_index, codim,
                                                                                             1) } };
            const boost::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isUpper(location_index, codim,
                                                                                             1) } };
#endif
#if (NDIM == 3)
            const boost::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isLower(location_index, codim, 1),
                                                          PhysicalBoundaryUtilities::isLower(location_index, codim,
                                                                                             2) } };
            const boost::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isUpper(location_index, codim, 1),
                                                          PhysicalBoundaryUtilities::isUpper(location_index, codim,
                                                                                             2) } };
#endif
            for (int depth = 0; depth < patch_data->getDepth(); ++depth)
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    const Box it_box = bdry_fill_box * ghost_box;
                    for (auto b = FaceGeometry::begin(it_box, axis), e = FaceGeometry::end(it_box, axis); b != e; ++b)
                    {
                        const FaceIndex i = *b;
                        FaceIndex i_bdry = i;
                        IntVector i_shft = IntVector::getZero(DIM);
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            if (is_lower[d])
                            {
                                i_bdry((NDIM - axis + d) % NDIM) = patch_lower(d);
                                i_shft((NDIM - axis + d) % NDIM) = +1; // use interior data for extrapolation
                            }
                            else if (is_upper[d])
                            {
                                if (axis != d)
                                {
                                    i_bdry((NDIM - axis + d) % NDIM) = patch_upper(d);
                                }
                                else
                                {
                                    i_bdry((NDIM - axis + d) % NDIM) = patch_upper(d) + 1;
                                }
                                i_shft((NDIM - axis + d) % NDIM) = -1; // use interior data for extrapolation
                            }
                        }

                        // Perform constant, linear, or quadratic extrapolation.
                        switch (extrap_type)
                        {
                        case 0:
                            (*patch_data)(i, depth) = (*patch_data)(i_bdry, depth);
                            break;
                        case 1:
                            (*patch_data)(i, depth) = compute_linear_extrap(*patch_data, i, i_bdry, i_shft, depth);
                            break;
                        case 2:
                            (*patch_data)(i, depth) =
                                compute_quadratic_extrap(*patch_data, i, i_bdry, i_shft, depth, codim);
                            break;
                        }
                    }
                }
            }
        }
    }
    return;
}

void CartExtrapPhysBdryOp::setPhysicalBoundaryConditions_node(
    Patch& patch,
    const std::vector<std::pair<Box, std::pair<int, int> > >& bdry_fill_boxes)
{
    const Box& patch_box = patch.getBox();
    const Index& patch_lower = patch_box.lower();
    const Index& patch_upper = patch_box.upper();

    const int extrap_type =
        (d_extrap_type == "CONSTANT" ? 0 : (d_extrap_type == "LINEAR" ? 1 : (d_extrap_type == "QUADRATIC" ? 2 : -1)));
    if (extrap_type != 0 && extrap_type != 1 && extrap_type != 2)
    {
        TBOX_ERROR("CartExtrapPhysBdryOp::setPhysicalBoundaryConditions():\n"
                   << "  unknown extrapolation type: " << d_extrap_type << "\n"
                   << "  valid selections are: CONSTANT, LINEAR, or QUADRATIC" << std::endl);
    }

    // Set the physical boundary conditions for the specified patch data
    // indices.
    for (auto cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        auto var_db = VariableDatabase::getDatabase();
        boost::shared_ptr<Variable> var;
        var_db->mapIndexToVariable(patch_data_idx, var);
        auto nc_var = boost::dynamic_pointer_cast<NodeVariable<double> >(var);
        if (!nc_var) continue;

        auto patch_data = BOOST_CAST<NodeData<double> >(patch.getPatchData(patch_data_idx));
        const Box& ghost_box = patch_data->getGhostBox();

        // Loop over the boundary fill boxes and extrapolate the data.
        for (auto it = bdry_fill_boxes.begin(); it != bdry_fill_boxes.end(); ++it)
        {
            const Box& bdry_fill_box = it->first;
            const unsigned int location_index = it->second.first;
            const int codim = it->second.second;
#if (NDIM == 2)
            const boost::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isLower(location_index, codim,
                                                                                             1) } };
            const boost::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isUpper(location_index, codim,
                                                                                             1) } };
#endif
#if (NDIM == 3)
            const boost::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isLower(location_index, codim, 1),
                                                          PhysicalBoundaryUtilities::isLower(location_index, codim,
                                                                                             2) } };
            const boost::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isUpper(location_index, codim, 1),
                                                          PhysicalBoundaryUtilities::isUpper(location_index, codim,
                                                                                             2) } };
#endif
            // Loop over the boundary box indices and compute the
            // nearest interior index.
            for (int depth = 0; depth < patch_data->getDepth(); ++depth)
            {
                const Box it_box = bdry_fill_box * ghost_box;
                for (auto b = NodeGeometry::begin(it_box), e = NodeGeometry::end(it_box); b != e; ++b)
                {
                    const NodeIndex& i = *b;
                    NodeIndex i_bdry = i;
                    IntVector i_shft = IntVector::getZero(DIM);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (is_lower[d])
                        {
                            i_bdry(d) = patch_lower(d);
                            i_shft(d) = +1; // use interior data for extrapolation
                        }
                        else if (is_upper[d])
                        {
                            i_bdry(d) = patch_upper(d) + 1;
                            i_shft(d) = -1; // use interior data for extrapolation
                        }
                    }

                    // Perform constant, linear, or quadratic extrapolation.
                    switch (extrap_type)
                    {
                    case 0:
                        (*patch_data)(i, depth) = (*patch_data)(i_bdry, depth);
                        break;
                    case 1:
                        (*patch_data)(i, depth) = compute_linear_extrap(*patch_data, i, i_bdry, i_shft, depth);
                        break;
                    case 2:
                        (*patch_data)(i, depth) =
                            compute_quadratic_extrap(*patch_data, i, i_bdry, i_shft, depth, codim);
                        break;
                    }
                }
            }
        }
    }
    return;
}

void CartExtrapPhysBdryOp::setPhysicalBoundaryConditions_side(
    Patch& patch,
    const std::vector<std::pair<Box, std::pair<int, int> > >& bdry_fill_boxes)
{
    const Box& patch_box = patch.getBox();
    const Index& patch_lower = patch_box.lower();
    const Index& patch_upper = patch_box.upper();

    const int extrap_type =
        (d_extrap_type == "CONSTANT" ? 0 : (d_extrap_type == "LINEAR" ? 1 : (d_extrap_type == "QUADRATIC" ? 2 : -1)));
    if (extrap_type != 0 && extrap_type != 1 && extrap_type != 2)
    {
        TBOX_ERROR("CartExtrapPhysBdryOp::setPhysicalBoundaryConditions():\n"
                   << "  unknown extrapolation type: " << d_extrap_type << "\n"
                   << "  valid selections are: CONSTANT, LINEAR, or QUADRATIC" << std::endl);
    }

    // Set the physical boundary conditions for the specified patch data
    // indices.
    for (auto cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int patch_data_idx = (*cit);
        auto var_db = VariableDatabase::getDatabase();
        boost::shared_ptr<Variable> var;
        var_db->mapIndexToVariable(patch_data_idx, var);
        auto sc_var = boost::dynamic_pointer_cast<SideVariable<double> >(var);
        if (!sc_var) continue;
        
        auto patch_data = BOOST_CAST<SideData<double> >(patch.getPatchData(patch_data_idx));
        const Box& ghost_box = patch_data->getGhostBox();

        // Loop over the boundary fill boxes and extrapolate the data.
        for (auto it = bdry_fill_boxes.begin(); it != bdry_fill_boxes.end(); ++it)
        {
            const Box& bdry_fill_box = it->first;
            const unsigned int location_index = it->second.first;
            const int codim = it->second.second;
#if (NDIM == 2)
            const boost::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isLower(location_index, codim,
                                                                                             1) } };
            const boost::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isUpper(location_index, codim,
                                                                                             1) } };
#endif
#if (NDIM == 3)
            const boost::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isLower(location_index, codim, 1),
                                                          PhysicalBoundaryUtilities::isLower(location_index, codim,
                                                                                             2) } };
            const boost::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                          PhysicalBoundaryUtilities::isUpper(location_index, codim, 1),
                                                          PhysicalBoundaryUtilities::isUpper(location_index, codim,
                                                                                             2) } };
#endif
            for (int depth = 0; depth < patch_data->getDepth(); ++depth)
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    const Box it_box = bdry_fill_box * ghost_box;
                    for (auto b = SideGeometry::begin(it_box, axis), e = SideGeometry::end(it_box, axis); b != e; ++b)
                    {
                        const SideIndex i = *b;
                        SideIndex i_bdry = i;
                        IntVector i_shft = IntVector::getZero(DIM);
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            if (is_lower[d])
                            {
                                i_bdry(d) = patch_lower(d);
                                i_shft(d) = +1; // use interior data for extrapolation
                            }
                            else if (is_upper[d])
                            {
                                if (axis != d)
                                {
                                    i_bdry(d) = patch_upper(d);
                                }
                                else
                                {
                                    i_bdry(d) = patch_upper(d) + 1;
                                }
                                i_shft(d) = -1; // use interior data for extrapolation
                            }
                        }

                        // Perform constant, linear, or quadratic extrapolation.
                        switch (extrap_type)
                        {
                        case 0:
                            (*patch_data)(i, depth) = (*patch_data)(i_bdry, depth);
                            break;
                        case 1:
                            (*patch_data)(i, depth) = compute_linear_extrap(*patch_data, i, i_bdry, i_shft, depth);
                            break;
                        case 2:
                            (*patch_data)(i, depth) =
                                compute_quadratic_extrap(*patch_data, i, i_bdry, i_shft, depth, codim);
                            break;
                        }
                    }
                }
            }
        }
    }
    return;
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
