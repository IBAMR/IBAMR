// Filename: IBLagrangianForceStrategy.cpp
// Created on 03 May 2005 by Boyce Griffith
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
#include <vector>

#include "ibamr/IBLagrangianForceStrategy.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "petscmat.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

namespace IBTK
{
class LData;
} // namespace IBTK

namespace IBTK
{
class LDataManager;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBLagrangianForceStrategy::IBLagrangianForceStrategy()
{
    // intentionally blank
    return;
} // IBLagrangianForceStrategy

IBLagrangianForceStrategy::~IBLagrangianForceStrategy()
{
    // intentionally blank
    return;
} // ~IBLagrangianForceStrategy

void
IBLagrangianForceStrategy::setTimeInterval(const double /*current_time*/, const double /*new_time*/)
{
    // intentionally blank
    return;
} // setTimeInterval

void
IBLagrangianForceStrategy::initializeLevelData(const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                               const int /*level_number*/,
                                               const double /*init_data_time*/,
                                               const bool /*initial_time*/,
                                               LDataManager* const /*l_data_manager*/)
{
    // intentionally blank
    return;
} // initializeLevelData

void
IBLagrangianForceStrategy::computeLagrangianForce(Pointer<LData> /*F_data*/,
                                                  Pointer<LData> /*X_data*/,
                                                  Pointer<LData> /*U_data*/,
                                                  const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                  const int /*level_number*/,
                                                  const double /*data_time*/,
                                                  LDataManager* const /*l_data_manager*/)
{
    TBOX_ERROR("IBLagrangianForceStrategy::computeLagrangianForce():\n"
               << "  not implemented for this IBLagrangianForceStrategy."
               << std::endl);
    return;
} // computeLagrangianForce

void
IBLagrangianForceStrategy::computeLagrangianForceJacobianNonzeroStructure(
    std::vector<int>& /*d_nnz*/,
    std::vector<int>& /*o_nnz*/,
    const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    const int /*level_number*/,
    LDataManager* const /*l_data_manager*/)
{
    TBOX_ERROR("IBLagrangianForceStrategy::computeLagrangianForceJacobianNonzeroStructure():\n"
               << "  Jacobian functionality not implemented for this IBLagrangianForceStrategy."
               << std::endl);
    return;
} // computeLagrangianForceJacobianNonzeroStructure

void
IBLagrangianForceStrategy::computeLagrangianForceJacobian(Mat& /*J_mat*/,
                                                          MatAssemblyType /*assembly_type*/,
                                                          const double /*X_coef*/,
                                                          Pointer<LData> /*X_data*/,
                                                          const double /*U_coef*/,
                                                          Pointer<LData> /*U_data*/,
                                                          const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                          const int /*level_number*/,
                                                          const double /*data_time*/,
                                                          LDataManager* const /*l_data_manager*/)
{
    TBOX_ERROR("IBLagrangianForceStrategy::computeLagrangianForceJacobian():\n"
               << "  Jacobian functionality not implemented for this IBLagrangianForceStrategy."
               << std::endl);
    return;
} // computeLagrangianForceJacobian

double
IBLagrangianForceStrategy::computeLagrangianEnergy(Pointer<LData> /*X_data*/,
                                                   Pointer<LData> /*U_data*/,
                                                   const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                   const int /*level_number*/,
                                                   const double /*data_time*/,
                                                   LDataManager* const /*l_data_manager*/)
{
    TBOX_ERROR("IBLagrangianForceStrategy::computeLagrangianEnergy():\n"
               << "  potential energy functionality not implemented for this "
                  "IBLagrangianForceStrategy."
               << std::endl);
    return 0.0;
} // computeLagrangianEnergy

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
