// Filename: CartGridFunction.cpp
// Created on 15 Mar 2004 by Boyce Griffith
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

#include <ostream>
#include <string>

#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "Variable.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartGridFunction::CartGridFunction(const std::string& object_name) : d_object_name(object_name)
{
    // intentionally blank
    return;
} // CartGridFunction

CartGridFunction::~CartGridFunction()
{
    // intentionally blank
    return;
} // ~CartGridFunction

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
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
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
