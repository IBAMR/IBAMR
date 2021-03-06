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

#include "ibamr/IBLagrangianForceStrategy.h"

#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include "petscmat.h"

#include <ostream>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

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
               << "  not implemented for this IBLagrangianForceStrategy." << std::endl);
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
               << "  Jacobian functionality not implemented for this IBLagrangianForceStrategy." << std::endl);
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
               << "  Jacobian functionality not implemented for this IBLagrangianForceStrategy." << std::endl);
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
