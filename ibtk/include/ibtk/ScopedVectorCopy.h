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

#ifndef included_IBTK_ScopedVectorCopy
#define included_IBTK_ScopedVectorCopy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <SAMRAIVectorReal.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/**
 * Wrapper class around a SAMRAIVectorReal with RAII semantics.
 */
template <typename TYPE>
class ScopedVectorCopy
{
public:
    ScopedVectorCopy(const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> >& vector,
                     const std::string& name = "")
        : ScopedVectorCopy(*vector, name)
    {
    }

    ScopedVectorCopy(const SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE>& vector, const std::string& name = "")
    {
        d_vector = vector.cloneVector(name);
        d_vector->allocateVectorData();
        d_vector->copyVector(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> >(
                                 const_cast<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE>*>(&vector), false),
                             false);
    }

    operator SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE>()
    {
        return *d_vector;
    }

    operator SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> >()
    {
        return d_vector;
    }

    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> > > getComponentVectors() const
    {
        // Setup SAMRAIVectorReal objects to correspond to the individual vector
        // components.
        std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> > > comps;
        for (int comp = 0; comp < d_vector->getNumberOfComponents(); ++comp)
        {
            comps.emplace_back(new SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE>(d_vector->getName() + "_component_" +
                                                                                  std::to_string(comp),
                                                                              d_vector->getPatchHierarchy(),
                                                                              d_vector->getCoarsestLevelNumber(),
                                                                              d_vector->getFinestLevelNumber()));
            comps.back()->addComponent(d_vector->getComponentVariable(comp),
                                       d_vector->getComponentDescriptorIndex(comp),
                                       d_vector->getControlVolumeIndex(comp));
        }
        return comps;
    }

    ~ScopedVectorCopy()
    {
        if (d_vector)
        {
            deallocate_vector_data(*d_vector);
            free_vector_components(*d_vector);
        }
    }

private:
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, TYPE> > d_vector;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_ScopedVectorCopy
