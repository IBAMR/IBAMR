// Filename: IBLagrangianForceStrategySet.cpp
// Created on 04 April 2007 by Boyce Griffith
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

#include "ibamr/IBLagrangianForceStrategySet.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LData.h"

namespace IBTK
{
class LDataManager;
} // namespace IBTK

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBLagrangianForceStrategySet::~IBLagrangianForceStrategySet()
{
    // intentionally blank
    return;
}

void IBLagrangianForceStrategySet::setTimeInterval(const double current_time, const double new_time)
{
    for (auto cit = d_strategy_set.begin();
         cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->setTimeInterval(current_time, new_time);
    }
    return;
}

void IBLagrangianForceStrategySet::initializeLevelData(const boost::shared_ptr<PatchHierarchy > hierarchy,
                                                       const int level_number,
                                                       const double init_data_time,
                                                       const bool initial_time,
                                                       LDataManager* const l_data_manager)
{
    for (auto cit = d_strategy_set.begin();
         cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->initializeLevelData(hierarchy, level_number, init_data_time, initial_time, l_data_manager);
    }
    return;
}

void IBLagrangianForceStrategySet::computeLagrangianForce(boost::shared_ptr<LData> F_data,
                                                          boost::shared_ptr<LData> X_data,
                                                          boost::shared_ptr<LData> U_data,
                                                          const boost::shared_ptr<PatchHierarchy > hierarchy,
                                                          const int level_number,
                                                          const double data_time,
                                                          LDataManager* const l_data_manager)
{
    for (auto cit = d_strategy_set.begin();
         cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->computeLagrangianForce(F_data, X_data, U_data, hierarchy, level_number, data_time, l_data_manager);
    }
    return;
}

void IBLagrangianForceStrategySet::computeLagrangianForceJacobianNonzeroStructure(
    std::vector<int>& d_nnz,
    std::vector<int>& o_nnz,
    const boost::shared_ptr<PatchHierarchy > hierarchy,
    const int level_number,
    LDataManager* const l_data_manager)
{
    for (auto cit = d_strategy_set.begin();
         cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->computeLagrangianForceJacobianNonzeroStructure(d_nnz, o_nnz, hierarchy, level_number, l_data_manager);
    }
    return;
}

void IBLagrangianForceStrategySet::computeLagrangianForceJacobian(Mat& J_mat,
                                                                  MatAssemblyType assembly_type,
                                                                  const double X_coef,
                                                                  boost::shared_ptr<LData> X_data,
                                                                  const double U_coef,
                                                                  boost::shared_ptr<LData> U_data,
                                                                  const boost::shared_ptr<PatchHierarchy > hierarchy,
                                                                  const int level_number,
                                                                  const double data_time,
                                                                  LDataManager* const l_data_manager)
{
    for (auto cit = d_strategy_set.begin();
         cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->computeLagrangianForceJacobian(J_mat,
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
}

double IBLagrangianForceStrategySet::computeLagrangianEnergy(boost::shared_ptr<LData> X_data,
                                                             boost::shared_ptr<LData> U_data,
                                                             const boost::shared_ptr<PatchHierarchy > hierarchy,
                                                             const int level_number,
                                                             const double data_time,
                                                             LDataManager* const l_data_manager)
{
    double ret_val = 0.0;
    for (auto cit = d_strategy_set.begin();
         cit != d_strategy_set.end();
         ++cit)
    {
        ret_val += (*cit)->computeLagrangianEnergy(X_data, U_data, hierarchy, level_number, data_time, l_data_manager);
    }
    return ret_val;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
