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

#include "ibamr/IBLagrangianForceStrategySet.h"

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LData.h"

#include "IntVector.h"
#include "PatchHierarchy.h"

#include <string>

#include "ibamr/namespaces.h" // IWYU pragma: keep

namespace IBTK
{
class LDataManager;
} // namespace IBTK

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
IBLagrangianForceStrategySet::setTimeInterval(const double current_time, const double new_time)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->setTimeInterval(current_time, new_time);
    }
    return;
} // setTimeInterval

void
IBLagrangianForceStrategySet::initializeLevelData(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                  const int level_number,
                                                  const double init_data_time,
                                                  const bool initial_time,
                                                  LDataManager* const l_data_manager)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->initializeLevelData(hierarchy, level_number, init_data_time, initial_time, l_data_manager);
    }
    return;
} // initializeLevelData

void
IBLagrangianForceStrategySet::computeLagrangianForce(Pointer<LData> F_data,
                                                     Pointer<LData> X_data,
                                                     Pointer<LData> U_data,
                                                     const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                     const int level_number,
                                                     const double data_time,
                                                     LDataManager* const l_data_manager)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->computeLagrangianForce(F_data, X_data, U_data, hierarchy, level_number, data_time, l_data_manager);
    }
    return;
} // computeLagrangianForce

void
IBLagrangianForceStrategySet::computeLagrangianForceJacobianNonzeroStructure(
    std::vector<int>& d_nnz,
    std::vector<int>& o_nnz,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    LDataManager* const l_data_manager)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->computeLagrangianForceJacobianNonzeroStructure(d_nnz, o_nnz, hierarchy, level_number, l_data_manager);
    }
    return;
} // computeLagrangianForceJacobianNonzeroStructure

void
IBLagrangianForceStrategySet::computeLagrangianForceJacobian(Mat& J_mat,
                                                             MatAssemblyType assembly_type,
                                                             const double X_coef,
                                                             Pointer<LData> X_data,
                                                             const double U_coef,
                                                             Pointer<LData> U_data,
                                                             const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                             const int level_number,
                                                             const double data_time,
                                                             LDataManager* const l_data_manager)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->computeLagrangianForceJacobian(J_mat,
                                                 MAT_FLUSH_ASSEMBLY,
                                                 X_coef,
                                                 X_data,
                                                 U_coef,
                                                 U_data,
                                                 hierarchy,
                                                 level_number,
                                                 data_time,
                                                 l_data_manager);
    }
    if (assembly_type != MAT_FLUSH_ASSEMBLY)
    {
        int ierr;
        ierr = MatAssemblyBegin(J_mat, assembly_type);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(J_mat, assembly_type);
        IBTK_CHKERRQ(ierr);
    }
    return;
} // computeLagrangianForceJacobian

double
IBLagrangianForceStrategySet::computeLagrangianEnergy(Pointer<LData> X_data,
                                                      Pointer<LData> U_data,
                                                      const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                      const int level_number,
                                                      const double data_time,
                                                      LDataManager* const l_data_manager)
{
    double ret_val = 0.0;
    for (const auto& strategy : d_strategy_set)
    {
        ret_val +=
            strategy->computeLagrangianEnergy(X_data, U_data, hierarchy, level_number, data_time, l_data_manager);
    }
    return ret_val;
} // computeLagrangianEnergy

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
