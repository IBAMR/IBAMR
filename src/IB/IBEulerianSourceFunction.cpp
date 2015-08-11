// Filename: IBEulerianSourceFunction.cpp
// Created on 18 Jun 2005 by Boyce Griffith
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

#include "SAMRAI/pdat/CellData.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/math/PatchCellDataBasicOps.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartGridFunction.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{

class Variable;

class PatchLevel;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

////////////////////////////// PUBLIC ///////////////////////////////////////

IBHierarchyIntegrator::IBEulerianSourceFunction::IBEulerianSourceFunction(const IBHierarchyIntegrator* const ib_solver)
    : CartGridFunction(ib_solver->getName() + "::IBEulerianSourceFunction"), d_ib_solver(ib_solver)
{
    // intentionally blank
    return;
}

IBHierarchyIntegrator::IBEulerianSourceFunction::~IBEulerianSourceFunction()
{
    // intentionally blank
    return;
}

bool IBHierarchyIntegrator::IBEulerianSourceFunction::isTimeDependent() const
{
    return true;
}

void IBHierarchyIntegrator::IBEulerianSourceFunction::setDataOnPatch(const int data_idx,
                                                                     const boost::shared_ptr<Variable>& /*var*/,
                                                                     const boost::shared_ptr<Patch>& patch,
                                                                     const double /*data_time*/,
                                                                     const bool initial_time,
                                                                     const boost::shared_ptr<PatchLevel>& /*level*/)
{
    auto q_cc_data = BOOST_CAST<CellData<double>>(patch->getPatchData(data_idx));
    q_cc_data->fillAll(0.0);
    if (initial_time) return;
    auto q_ib_cc_data = BOOST_CAST<CellData<double>>(patch->getPatchData(d_ib_solver->d_q_idx));
    PatchCellDataBasicOps<double> patch_ops;
    patch_ops.add(q_cc_data, q_cc_data, q_ib_cc_data, patch->getBox());
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
