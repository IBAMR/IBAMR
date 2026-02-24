// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_SAMRAIScopedVectorDuplicate_inl_h
#define included_IBTK_SAMRAIScopedVectorDuplicate_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/samrai_compatibility_names.h"
#include <ibtk/SAMRAIScopedVectorDuplicate.h>
#include <ibtk/ibtk_utilities.h>

#include "SAMRAIPointer.h"
#include "SAMRAISAMRAIVectorReal.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
template <typename TYPE>
SAMRAIScopedVectorDuplicate<TYPE>::SAMRAIScopedVectorDuplicate(
    const SAMRAIPointer<SAMRAISAMRAIVectorReal<TYPE>>& vector,
    const std::string& name)
    : SAMRAIScopedVectorDuplicate(checked_dereference(vector), name)
{
}

template <typename TYPE>
SAMRAIScopedVectorDuplicate<TYPE>::SAMRAIScopedVectorDuplicate(const SAMRAISAMRAIVectorReal<TYPE>& vector,
                                                               const std::string& name)
{
    d_vector = vector.cloneVector(name);
    d_vector->allocateVectorData();
    d_vector->setToScalar(TYPE(0.0), /*interior_only*/ false);
}

template <typename TYPE>
SAMRAIScopedVectorDuplicate<TYPE>::operator SAMRAISAMRAIVectorReal<TYPE>&()
{
    return *d_vector;
}

template <typename TYPE>
SAMRAIScopedVectorDuplicate<TYPE>::operator SAMRAIPointer<SAMRAISAMRAIVectorReal<TYPE>>()
{
    return SAMRAIPointer<SAMRAISAMRAIVectorReal<double>>(&*d_vector, false);
}

template <typename TYPE>
std::vector<SAMRAIPointer<SAMRAISAMRAIVectorReal<TYPE>>>
SAMRAIScopedVectorDuplicate<TYPE>::getComponentVectors() const
{
    // Setup SAMRAIVectorReal objects to correspond to the individual vector
    // components.
    std::vector<SAMRAIPointer<SAMRAISAMRAIVectorReal<TYPE>>> comps;
    for (int comp = 0; comp < d_vector->getNumberOfComponents(); ++comp)
    {
        comps.emplace_back(new SAMRAISAMRAIVectorReal<TYPE>(d_vector->getName() + "_component_" + std::to_string(comp),
                                                            d_vector->getPatchHierarchy(),
                                                            d_vector->getCoarsestLevelNumber(),
                                                            d_vector->getFinestLevelNumber()));
        comps.back()->addComponent(d_vector->getComponentVariable(comp),
                                   d_vector->getComponentDescriptorIndex(comp),
                                   d_vector->getControlVolumeIndex(comp));
    }
    return comps;
}

template <typename TYPE>
SAMRAIScopedVectorDuplicate<TYPE>::~SAMRAIScopedVectorDuplicate()
{
    if (d_vector)
    {
        deallocate_vector_data(*d_vector);
        free_vector_components(*d_vector);
    }
}
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_SAMRAIScopedVectorDuplicate_inl_h
