// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
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

#include <ibtk/private/VectorSignature.h>

#include <tbox/Pointer.h>

#include <Box.h>
#include <Index.h>
#include <PatchDataFactory.h>
#include <PatchDescriptor.h>
#include <VariableDatabase.h>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
namespace detail
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

VectorSignature::VectorSignature(const SAMRAIVectorReal<NDIM, double>& vec)
    : d_hierarchy(vec.getPatchHierarchy()),
      d_coarsest_ln(vec.getCoarsestLevelNumber()),
      d_finest_ln(vec.getFinestLevelNumber())
{
    const Pointer<PatchDescriptor<NDIM>> patch_descriptor = VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor();
    for (int comp = 0; comp < vec.getNumberOfComponents(); ++comp)
    {
        Component c;
        c.data_idx = vec.getComponentDescriptorIndex(comp);
        const Pointer<PatchDataFactory<NDIM>> factory = patch_descriptor->getPatchDataFactory(c.data_idx);
        const PatchDataFactory<NDIM>& factory_ref = *factory;
        c.variable = vec.getComponentVariable(comp);
        c.factory_type = typeid(factory_ref);
        c.memory_size = factory->getSizeOfMemory(Box<NDIM>(hier::Index<NDIM>(0), hier::Index<NDIM>(0)));
        c.ghost_width = factory->getGhostCellWidth();
        c.control_volume_idx = vec.getControlVolumeIndex(comp);
        d_components.push_back(c);
    }
    return;
} // VectorSignature

bool
VectorSignature::operator==(const VectorSignature& other) const
{
    return d_coarsest_ln == other.d_coarsest_ln && d_finest_ln == other.d_finest_ln &&
           hasSameStructure(other, /*compare_data_indices*/ true);
} // operator==

bool
VectorSignature::hasSameStructure(const VectorSignature& other) const
{
    return hasSameStructure(other, /*compare_data_indices*/ false);
} // hasSameStructure

/////////////////////////////// PRIVATE //////////////////////////////////////

bool
VectorSignature::hasSameStructure(const VectorSignature& other, const bool compare_data_indices) const
{
    if (d_hierarchy.isNull() || d_hierarchy.getPointer() != other.d_hierarchy.getPointer() ||
        d_components.size() != other.d_components.size())
    {
        return false;
    }
    for (std::size_t k = 0; k < d_components.size(); ++k)
    {
        const Component& a = d_components[k];
        const Component& b = other.d_components[k];
        if (a.variable.getPointer() != b.variable.getPointer() || a.factory_type != b.factory_type ||
            a.memory_size != b.memory_size || !(a.ghost_width == b.ghost_width) ||
            a.control_volume_idx != b.control_volume_idx || (compare_data_indices && a.data_idx != b.data_idx))
        {
            return false;
        }
    }
    return true;
} // hasSameStructure

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
