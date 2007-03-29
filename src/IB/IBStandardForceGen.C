// Filename: IBStandardForceGen.C
// Created on 03 May 2005 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)
// Last modified: <29.Mar.2007 16:17:33 griffith@box221.cims.nyu.edu>

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
    : d_spring_force_gen(spring_force_gen),
      d_beam_force_gen(beam_force_gen),
      d_target_point_force_gen(target_point_force_gen)
{
    if (d_spring_force_gen.isNull())
    {
        TBOX_WARNING("IBStandardForceGen::initializeLevelData():\n"
                     << "  spring forces disabled." << endl);
    }
    if (d_beam_force_gen.isNull())
    {
        TBOX_WARNING("IBStandardForceGen::initializeLevelData():\n"
                     << "  beam forces disabled." << endl);
    }
    if (d_target_point_force_gen.isNull())
    {
        TBOX_WARNING("IBStandardForceGen::initializeLevelData():\n"
                     << "  target point forces disabled." << endl);
    }
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
    const LDataManager* const lag_manager)
{
    // Initialize the spring force generator.
    if (!d_spring_force_gen.isNull())
    {
        d_spring_force_gen->initializeLevelData(
            hierarchy, level_number, init_data_time, initial_time, lag_manager);
    }

    // Initialize the beam force generator.
    if (!d_beam_force_gen.isNull())
    {
        d_beam_force_gen->initializeLevelData(
            hierarchy, level_number, init_data_time, initial_time, lag_manager);
    }

    // Initialize the target point force generator.
    if (!d_target_point_force_gen.isNull())
    {
        d_target_point_force_gen->initializeLevelData(
            hierarchy, level_number, init_data_time, initial_time, lag_manager);
    }
    return;
}// initializeLevelData

void
IBStandardForceGen::computeLagrangianForce(
    SAMRAI::tbox::Pointer<LNodeLevelData> F_data,
    SAMRAI::tbox::Pointer<LNodeLevelData> X_data,
    SAMRAI::tbox::Pointer<LNodeLevelData> U_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    const LDataManager* const lag_manager)
{
    // Compute the spring forces.
    if (!d_spring_force_gen.isNull())
    {
        d_spring_force_gen->computeLagrangianForce(
            F_data, X_data, U_data, hierarchy, level_number, data_time, lag_manager);
    }

    // Compute the beam forces.
    if (!d_beam_force_gen.isNull())
    {
        d_beam_force_gen->computeLagrangianForce(
            F_data, X_data, U_data, hierarchy, level_number, data_time, lag_manager);
    }

    // Compute the target point forces.
    if (!d_target_point_force_gen.isNull())
    {
        d_target_point_force_gen->computeLagrangianForce(
            F_data, X_data, U_data, hierarchy, level_number, data_time, lag_manager);
    }
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
