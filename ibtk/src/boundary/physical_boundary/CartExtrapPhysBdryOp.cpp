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

#include "ibtk/CartExtrapPhysBdryOp.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIArray.h"
#include "SAMRAIBoundaryBox.h"
#include "SAMRAIBox.h"
#include "SAMRAICellData.h"
#include "SAMRAICellIndex.h"
#include "SAMRAICellIterator.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIComponentSelector.h"
#include "SAMRAIFaceData.h"
#include "SAMRAIFaceIndex.h"
#include "SAMRAIFaceIterator.h"
#include "SAMRAIFaceVariable.h"
#include "SAMRAIIndex.h"
#include "SAMRAIIntVector.h"
#include "SAMRAINodeData.h"
#include "SAMRAINodeIndex.h"
#include "SAMRAINodeIterator.h"
#include "SAMRAINodeVariable.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchGeometry.h"
#include "SAMRAIPointer.h"
#include "SAMRAISideData.h"
#include "SAMRAISideIndex.h"
#include "SAMRAISideIterator.h"
#include "SAMRAISideVariable.h"
#include "SAMRAIUtilities.h"
#include "SAMRAIVariable.h"
#include "SAMRAIVariableDatabase.h"

#include <array>
#include <ostream>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int REFINE_OP_STENCIL_WIDTH = 0;

template <typename D, typename I>
inline double
compute_linear_extrap(D& patch_data, const I& i, const I& i_intr, const SAMRAIIntVector& i_shft, const int depth)
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
} // compute_linear_extrap

template <typename D, typename I>
inline double
compute_quadratic_extrap(D& patch_data,
                         const I& i,
                         const I& i_intr,
                         const SAMRAIIntVector& i_shft,
                         const int depth,
                         const int codim)
{
    if (codim == 1)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (i_shft(d) != 0)
            {
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
                       (1.0 / 2.0 * f2 - 2.0 * f1 + 3.0 / 2.0 * f0) * x + f0;
            }
        }
    }
    else
    {
        return compute_linear_extrap(patch_data, i, i_intr, i_shft, depth);
    }
    return 0.0; // this statement should not be reached
} // compute_quadratic_extrap
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartExtrapPhysBdryOp::CartExtrapPhysBdryOp(const int patch_data_index, const std::string& extrap_type)
{
    setPatchDataIndex(patch_data_index);
    setExtrapolationType(extrap_type);
    return;
} // CartExtrapPhysBdryOp

CartExtrapPhysBdryOp::CartExtrapPhysBdryOp(const std::set<int>& patch_data_indices, const std::string& extrap_type)
{
    setPatchDataIndices(patch_data_indices);
    setExtrapolationType(extrap_type);
    return;
} // CartExtrapPhysBdryOp

CartExtrapPhysBdryOp::CartExtrapPhysBdryOp(const SAMRAIComponentSelector& patch_data_indices,
                                           const std::string& extrap_type)
{
    setPatchDataIndices(patch_data_indices);
    setExtrapolationType(extrap_type);
    return;
} // CartExtrapPhysBdryOp

void
CartExtrapPhysBdryOp::setPatchDataIndex(const int patch_data_index)
{
    std::set<int> patch_data_indices;
    patch_data_indices.insert(patch_data_index);
    setPatchDataIndices(patch_data_indices);
    return;
} // setPatchDataIndex

void
CartExtrapPhysBdryOp::setPatchDataIndices(const std::set<int>& patch_data_indices)
{
    d_patch_data_indices.clear();
    d_patch_data_indices = patch_data_indices;
    return;
} // setPatchDataIndices

void
CartExtrapPhysBdryOp::setPatchDataIndices(const SAMRAIComponentSelector& patch_data_indices)
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
} // setPatchDataIndices

void
CartExtrapPhysBdryOp::setExtrapolationType(const std::string& extrap_type)
{
    // Ensure that the extrapolation type is supported by this class.
    if (extrap_type != "CONSTANT" && extrap_type != "LINEAR" && extrap_type != "QUADRATIC")
    {
        TBOX_ERROR("CartExtrapPhysBdryOp::setExtrapolationType():\n"
                   << "  unknown extrapolation type: " << extrap_type << "\n"
                   << "  valid selections are: CONSTANT, LINEAR, or QUADRATIC" << std::endl);
    }
    d_extrap_type = extrap_type;
    return;
} // setExtrapolationType

void
CartExtrapPhysBdryOp::setPhysicalBoundaryConditions(SAMRAIPatch& patch,
                                                    const double /*fill_time*/,
                                                    const SAMRAIIntVector& ghost_width_to_fill)
{
    if (ghost_width_to_fill == SAMRAIIntVector(0)) return;

    SAMRAIPointer<SAMRAIPatchGeometry> pgeom = patch.getPatchGeometry();
    const SAMRAIBox& patch_box = patch.getBox();

    std::vector<std::pair<SAMRAIBox, std::pair<int, int> > > bdry_fill_boxes;

#if (NDIM > 1)
#if (NDIM > 2)
    // Compute the co-dimension three boundary fill boxes.
    const SAMRAIArray<SAMRAIBoundaryBox> physical_codim3_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim3Boxes(patch);
    const int n_physical_codim3_boxes = physical_codim3_boxes.size();
    for (int n = 0; n < n_physical_codim3_boxes; ++n)
    {
        const SAMRAIBoundaryBox& bdry_box = physical_codim3_boxes[n];
        const SAMRAIBox bdry_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
        const unsigned int location_index = bdry_box.getLocationIndex();
        const int codim = 3;
        bdry_fill_boxes.push_back(std::make_pair(bdry_fill_box, std::make_pair(location_index, codim)));
    }
#endif
    // Compute the co-dimension two boundary fill boxes.
    const SAMRAIArray<SAMRAIBoundaryBox> physical_codim2_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(patch);
    const int n_physical_codim2_boxes = physical_codim2_boxes.size();
    for (int n = 0; n < n_physical_codim2_boxes; ++n)
    {
        const SAMRAIBoundaryBox& bdry_box = physical_codim2_boxes[n];
        const SAMRAIBox bdry_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
        const unsigned int location_index = bdry_box.getLocationIndex();
        const int codim = 2;
        bdry_fill_boxes.push_back(std::make_pair(bdry_fill_box, std::make_pair(location_index, codim)));
    }
#endif
    // Compute the co-dimension one boundary fill boxes.
    const SAMRAIArray<SAMRAIBoundaryBox> physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const SAMRAIBoundaryBox& bdry_box = physical_codim1_boxes[n];
        const SAMRAIBox bdry_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
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
} // setPhysicalBoundaryConditions

SAMRAIIntVector
CartExtrapPhysBdryOp::getRefineOpStencilWidth() const
{
    return REFINE_OP_STENCIL_WIDTH;
} // getRefineOpStencilWidth

void
CartExtrapPhysBdryOp::preprocessRefine(SAMRAIPatch& /*fine*/,
                                       const SAMRAIPatch& /*coarse*/,
                                       const SAMRAIBox& /*fine_box*/,
                                       const SAMRAIIntVector& /*ratio*/)
{
    // intentionally blank
    return;
} // preprocessRefine

void
CartExtrapPhysBdryOp::postprocessRefine(SAMRAIPatch& /*fine*/,
                                        const SAMRAIPatch& /*coarse*/,
                                        const SAMRAIBox& /*fine_box*/,
                                        const SAMRAIIntVector& /*ratio*/)
{
    // intentionally blank
    return;
} // postprocessRefine

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CartExtrapPhysBdryOp::setPhysicalBoundaryConditions_cell(
    SAMRAIPatch& patch,
    const std::vector<std::pair<SAMRAIBox, std::pair<int, int> > >& bdry_fill_boxes)
{
    const SAMRAIBox& patch_box = patch.getBox();
    const SAMRAIIndex& patch_lower = patch_box.lower();
    const SAMRAIIndex& patch_upper = patch_box.upper();

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
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
        SAMRAIPointer<SAMRAIVariable> var;
        var_db->mapIndexToVariable(patch_data_idx, var);
        SAMRAIPointer<SAMRAICellVariable<double> > cc_var = var;
        if (!cc_var) continue;

        SAMRAIPointer<SAMRAICellData<double> > patch_data = patch.getPatchData(patch_data_idx);
        const SAMRAIBox& ghost_box = patch_data->getGhostBox();

        // Loop over the boundary fill boxes and extrapolate the data.
        for (const auto& bdry_fill_box_pair : bdry_fill_boxes)
        {
            const SAMRAIBox& bdry_fill_box = bdry_fill_box_pair.first;
            const unsigned int location_index = bdry_fill_box_pair.second.first;
            const int codim = bdry_fill_box_pair.second.second;
#if (NDIM == 2)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 1) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 1) } };
#endif
#if (NDIM == 3)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 2) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 2) } };
#endif
            // Loop over the boundary box indices and compute the nearest
            // interior index.
            for (int depth = 0; depth < patch_data->getDepth(); ++depth)
            {
                for (SAMRAICellIterator b(bdry_fill_box * ghost_box); b; b++)
                {
                    const SAMRAICellIndex& i = b();
                    SAMRAICellIndex i_intr = i;
                    SAMRAIIntVector i_shft = 0;
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
} // setPhysicalBoundaryConditions_cell

void
CartExtrapPhysBdryOp::setPhysicalBoundaryConditions_face(
    SAMRAIPatch& patch,
    const std::vector<std::pair<SAMRAIBox, std::pair<int, int> > >& bdry_fill_boxes)
{
    const SAMRAIBox& patch_box = patch.getBox();
    const SAMRAIIndex& patch_lower = patch_box.lower();
    const SAMRAIIndex& patch_upper = patch_box.upper();

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
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
        SAMRAIPointer<SAMRAIVariable> var;
        var_db->mapIndexToVariable(patch_data_idx, var);
        SAMRAIPointer<SAMRAIFaceVariable<double> > fc_var = var;
        if (!fc_var) continue;
        SAMRAIPointer<SAMRAIFaceData<double> > patch_data = patch.getPatchData(patch_data_idx);
        const SAMRAIBox& ghost_box = patch_data->getGhostBox();

        // Loop over the boundary fill boxes and extrapolate the data.
        for (const auto& bdry_fill_box_pair : bdry_fill_boxes)
        {
            const SAMRAIBox& bdry_fill_box = bdry_fill_box_pair.first;
            const unsigned int location_index = bdry_fill_box_pair.second.first;
            const int codim = bdry_fill_box_pair.second.second;
#if (NDIM == 2)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 1) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 1) } };
#endif
#if (NDIM == 3)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 2) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 2) } };
#endif
            for (int depth = 0; depth < patch_data->getDepth(); ++depth)
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    for (SAMRAIFaceIterator b(bdry_fill_box * ghost_box, axis); b; b++)
                    {
                        const SAMRAIFaceIndex i = b();
                        SAMRAIFaceIndex i_bdry = i;
                        SAMRAIIntVector i_shft = 0;
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
} // setPhysicalBoundaryConditions_face

void
CartExtrapPhysBdryOp::setPhysicalBoundaryConditions_node(
    SAMRAIPatch& patch,
    const std::vector<std::pair<SAMRAIBox, std::pair<int, int> > >& bdry_fill_boxes)
{
    const SAMRAIBox& patch_box = patch.getBox();
    const SAMRAIIndex& patch_lower = patch_box.lower();
    const SAMRAIIndex& patch_upper = patch_box.upper();

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
    for (int patch_data_idx : d_patch_data_indices)
    {
        SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
        SAMRAIPointer<SAMRAIVariable> var;
        var_db->mapIndexToVariable(patch_data_idx, var);
        SAMRAIPointer<SAMRAINodeVariable<double> > nc_var = var;
        if (!nc_var) continue;
        SAMRAIPointer<SAMRAINodeData<double> > patch_data = patch.getPatchData(patch_data_idx);
        const SAMRAIBox& ghost_box = patch_data->getGhostBox();

        // Loop over the boundary fill boxes and extrapolate the data.
        for (const auto& bdry_fill_box_pair : bdry_fill_boxes)
        {
            const SAMRAIBox& bdry_fill_box = bdry_fill_box_pair.first;
            const unsigned int location_index = bdry_fill_box_pair.second.first;
            const int codim = bdry_fill_box_pair.second.second;
#if (NDIM == 2)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 1) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 1) } };
#endif
#if (NDIM == 3)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 2) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 2) } };
#endif
            // Loop over the boundary box indices and compute the
            // nearest interior index.
            for (int depth = 0; depth < patch_data->getDepth(); ++depth)
            {
                for (SAMRAINodeIterator b(bdry_fill_box * ghost_box); b; b++)
                {
                    const SAMRAINodeIndex& i = b();
                    SAMRAINodeIndex i_bdry = i;
                    SAMRAIIntVector i_shft = 0;
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
} // setPhysicalBoundaryConditions_node

void
CartExtrapPhysBdryOp::setPhysicalBoundaryConditions_side(
    SAMRAIPatch& patch,
    const std::vector<std::pair<SAMRAIBox, std::pair<int, int> > >& bdry_fill_boxes)
{
    const SAMRAIBox& patch_box = patch.getBox();
    const SAMRAIIndex& patch_lower = patch_box.lower();
    const SAMRAIIndex& patch_upper = patch_box.upper();

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
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
        SAMRAIPointer<SAMRAIVariable> var;
        var_db->mapIndexToVariable(patch_data_idx, var);
        SAMRAIPointer<SAMRAISideVariable<double> > sc_var = var;
        if (!sc_var) continue;
        SAMRAIPointer<SAMRAISideData<double> > patch_data = patch.getPatchData(patch_data_idx);
        const SAMRAIBox& ghost_box = patch_data->getGhostBox();

        // Loop over the boundary fill boxes and extrapolate the data.
        for (const auto& bdry_fill_box_map : bdry_fill_boxes)
        {
            const SAMRAIBox& bdry_fill_box = bdry_fill_box_map.first;
            const unsigned int location_index = bdry_fill_box_map.second.first;
            const int codim = bdry_fill_box_map.second.second;
#if (NDIM == 2)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 1) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 1) } };
#endif
#if (NDIM == 3)
            const std::array<bool, NDIM> is_lower = { { PhysicalBoundaryUtilities::isLower(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isLower(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isLower(
                                                            location_index, codim, 2) } };
            const std::array<bool, NDIM> is_upper = { { PhysicalBoundaryUtilities::isUpper(location_index, codim, 0),
                                                        PhysicalBoundaryUtilities::isUpper(location_index, codim, 1),
                                                        PhysicalBoundaryUtilities::isUpper(
                                                            location_index, codim, 2) } };
#endif
            for (int depth = 0; depth < patch_data->getDepth(); ++depth)
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    for (SAMRAISideIterator b(bdry_fill_box * ghost_box, axis); b; b++)
                    {
                        const SAMRAISideIndex i = b();
                        SAMRAISideIndex i_bdry = i;
                        SAMRAIIntVector i_shft = 0;
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
} // setPhysicalBoundaryConditions_side

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
