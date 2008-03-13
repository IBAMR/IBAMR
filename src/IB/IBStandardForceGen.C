// Filename: IBStandardForceGen.C
// Last modified: <12.Mar.2008 23:00:40 griffith@box221.cims.nyu.edu>
// Created on 03 May 2005 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)

#include "IBStandardForceGen.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStandardForceGen::IBStandardForceGen(
    SAMRAI::tbox::Pointer<IBSpringForceGen> spring_force_gen,
    SAMRAI::tbox::Pointer<IBBeamForceGen> beam_force_gen,
    SAMRAI::tbox::Pointer<IBTargetPointForceGen> target_point_force_gen)
    : d_force_strategy_set(NULL)
{
    std::vector<SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> > strategy_set;

    if (spring_force_gen.isNull())
    {
        TBOX_WARNING("IBStandardForceGen::initializeLevelData():\n"
                     << "  spring forces disabled." << std::endl);
    }
    else
    {
        strategy_set.push_back(spring_force_gen);
    }

    if (beam_force_gen.isNull())
    {
        TBOX_WARNING("IBStandardForceGen::initializeLevelData():\n"
                     << "  beam forces disabled." << std::endl);
    }
    else
    {
        strategy_set.push_back(beam_force_gen);
    }

    if (target_point_force_gen.isNull())
    {
        TBOX_WARNING("IBStandardForceGen::initializeLevelData():\n"
                     << "  target point forces disabled." << std::endl);
    }
    else
    {
        strategy_set.push_back(target_point_force_gen);
    }

    d_force_strategy_set = new IBLagrangianForceStrategySet(
        strategy_set.begin(), strategy_set.end());
    return;
}// IBStandardForceGen

IBStandardForceGen::~IBStandardForceGen()
{
    // intentionally blank
    return;
}// ~IBStandardForceGen

void
IBStandardForceGen::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    IBTK::LDataManager* const lag_manager)
{
    d_force_strategy_set->initializeLevelData(
        hierarchy, level_number, init_data_time, initial_time, lag_manager);
    return;
}// initializeLevelData

void
IBStandardForceGen::computeLagrangianForce(
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> F_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    d_force_strategy_set->computeLagrangianForce(
        F_data, X_data, hierarchy, level_number, data_time, lag_manager);
    return;
}// computeLagrangianForce

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBStandardForceGen>;

//////////////////////////////////////////////////////////////////////////////
