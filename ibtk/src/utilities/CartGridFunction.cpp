// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2022 by the IBAMR developers
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

#include "ibtk/CartGridFunction.h"

#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "Variable.h"

#include <string>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartGridFunction::CartGridFunction(std::string object_name) : d_object_name(std::move(object_name))
{
    // intentionally blank
    return;
} // CartGridFunction

void
CartGridFunction::setDataOnPatchHierarchy(const int data_idx,
                                          Pointer<Variable<NDIM> > var,
                                          Pointer<PatchHierarchy<NDIM> > hierarchy,
                                          const double data_time,
                                          const bool initial_time,
                                          const int coarsest_ln_in,
                                          const int finest_ln_in)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    const int coarsest_ln = (coarsest_ln_in == invalid_level_number ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == invalid_level_number ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }
    return;
} // setDataOnPatchHierarchy

void
CartGridFunction::setDataOnPatchLevel(const int data_idx,
                                      Pointer<Variable<NDIM> > var,
                                      Pointer<PatchLevel<NDIM> > level,
                                      const double data_time,
                                      const bool initial_time)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(level);
#endif
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        setDataOnPatch(data_idx, var, level->getPatch(p()), data_time, initial_time, level);
    }
    return;
} // setDataOnPatchLevel

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
