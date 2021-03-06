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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_MaterialPointSpec_inl
#define included_IBAMR_MaterialPointSpec_inl

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/MaterialPointSpec.h"

#include "ibtk/StreamableManager.h"

#include "tbox/PIO.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline bool
MaterialPointSpec::getIsRegisteredWithStreamableManager()
{
    return (STREAMABLE_CLASS_ID != IBTK::StreamableManager::getUnregisteredID());
} // getIsRegisteredWithStreamableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline MaterialPointSpec::MaterialPointSpec(const int point_idx,
                                            const double weight,
                                            const libMesh::subdomain_id_type subdomain_id,
                                            const std::vector<double>& internal_vars)
    : d_point_idx(point_idx), d_weight(weight), d_subdomain_id(subdomain_id), d_internal_vars(internal_vars)
{
#if !defined(NDEBUG)
    if (!getIsRegisteredWithStreamableManager())
    {
        TBOX_ERROR("MaterialPointSpec::MaterialPointSpec():\n"
                   << "  must call MaterialPointSpec::registerWithStreamableManager() before\n"
                   << "  creating any MaterialPointSpec objects.\n");
    }
#endif
    return;
} // MaterialPointSpec

inline MaterialPointSpec::~MaterialPointSpec()
{
    // intentionally blank
    return;
} // ~MaterialPointSpec

inline const int&
MaterialPointSpec::getPointIndex() const
{
    return d_point_idx;
} // getPointIndex

inline int&
MaterialPointSpec::getPointIndex()
{
    return d_point_idx;
} // getPointIndex

inline const double&
MaterialPointSpec::getWeight() const
{
    return d_weight;
} // getWeight

inline double&
MaterialPointSpec::getWeight()
{
    return d_weight;
} // getWeight

inline const libMesh::subdomain_id_type&
MaterialPointSpec::getSubdomainId() const
{
    return d_subdomain_id;
} // getSubdomainId

inline libMesh::subdomain_id_type&
MaterialPointSpec::getSubdomainId()
{
    return d_subdomain_id;
} // getSubdomainId

inline const std::vector<double>&
MaterialPointSpec::getInternalVariables() const
{
    return d_internal_vars;
} // getInternalVariables

inline std::vector<double>&
MaterialPointSpec::getInternalVariables()
{
    return d_internal_vars;
} // getInternalVariables

inline int
MaterialPointSpec::getStreamableClassID() const
{
    return STREAMABLE_CLASS_ID;
} // getStreamableClassID

inline size_t
MaterialPointSpec::getDataStreamSize() const
{
    return (3 * SAMRAI::tbox::AbstractStream::sizeofInt() +
            (1 + d_internal_vars.size()) * SAMRAI::tbox::AbstractStream::sizeofDouble());
} // getDataStreamSize

inline void
MaterialPointSpec::packStream(SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_point_idx, 1);
    stream.pack(&d_weight, 1);
    const int subdomain_id = d_subdomain_id;
    stream.pack(&subdomain_id, 1);
    const int n_internal_vars = static_cast<int>(d_internal_vars.size());
    stream.pack(&n_internal_vars);
    if (n_internal_vars) stream.pack(&d_internal_vars[0], n_internal_vars);
    return;
} // packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_MaterialPointSpec_inl
