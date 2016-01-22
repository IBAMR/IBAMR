// Filename: ConstraintIBKinematics.cpp
// Created on 10 Dec 2011 by Amneet Bhalla
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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

/////////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include "tbox/Utilities.h"
#include "ibamr/namespaces.h"
#include "ibamr/ConstraintIBKinematics.h"

namespace IBAMR
{
namespace
{
} // namespace anonymous

ConstraintIBKinematics::StructureParameters::StructureParameters(Pointer<Database> input_db,
                                                                 LDataManager* l_data_manager)
    : d_total_nodes(0), d_tagged_pt_idx(-1), d_struct_is_self_translating(false), d_struct_is_self_rotating(false)
{
    Array<std::string> struct_names = input_db->getStringArray("structure_names");
    Array<int> struct_levels = input_db->getIntegerArray("structure_levels");

#if !defined(NDEBUG)
    TBOX_ASSERT(!struct_names.isNull());
    TBOX_ASSERT(!struct_levels.isNull());
#endif

    d_coarsest_ln = struct_levels[0];
    d_finest_ln = struct_levels[struct_levels.getSize() - 1];

    // Get the index range for the structures.
    for (int i = 0; i < struct_names.size(); ++i)
    {
        const int level = struct_levels[i];
        const int struct_id = l_data_manager->getLagrangianStructureID(struct_names[i], level);
        if (struct_id != -1)
        {
            std::pair<int, int> idx_range = l_data_manager->getLagrangianStructureIndexRange(struct_id, level);
            d_idx_range.push_back(idx_range);
            d_total_nodes += idx_range.second - idx_range.first;
        }
        else
        {
            TBOX_ERROR("StructureParameters::StructureParameters() Structure " << struct_names[i]
                                                                               << " does not exist on level "
                                                                               << level
                                                                               << std::endl);
        }
    }

    // Get options from database
    d_calculate_trans_mom = input_db->getIntegerArray("calculate_translational_momentum");
    d_calculate_rot_mom = input_db->getIntegerArray("calculate_rotational_momentum");
    for (int i = 0; i < 3; ++i)
    {
        if (d_calculate_trans_mom[i]) d_struct_is_self_translating = true;
        if (d_calculate_rot_mom[i]) d_struct_is_self_rotating = true;
    }
    // Only self translating bodies can be self rotating.
    if (d_struct_is_self_rotating && !d_struct_is_self_translating)
    {
        TBOX_ERROR("ERROR:: StructureParameters::StructureParameters( ) "
                   << "\n"
                   << "Only self translating/propelling bodies can be self rotating."
                   << std::endl);
    }

    d_lag_position_update_method = input_db->getString("lag_position_update_method");
    if (d_lag_position_update_method != "CONSTRAINT_VELOCITY" &&
        d_lag_position_update_method != "CONSTRAINT_POSITION" &&
        d_lag_position_update_method != "CONSTRAINT_EXPT_POSITION")
    {
        TBOX_ERROR(
            "ERROR:: StructureParameters::StructureParameters( ) "
            << "\n"
            << "Update methods supported are CONSTRAINT_VELOCITY, CONSTRAINT_POSITION AND  CONSTRAINT_EXPT_POSITION\n\n"
            << std::endl);
    }

    Array<int> tagged_pt_identifier = input_db->getIntegerArray("tagged_pt_identifier");
    const int level_tagged = tagged_pt_identifier[0];
    const int relative_idx_tagged = tagged_pt_identifier[1];
    for (int i = 0; i < struct_levels.size(); ++i)
    {
        if (struct_levels[i] == level_tagged)
        {
            d_tagged_pt_idx = d_idx_range[i].first + relative_idx_tagged;
            break;
        }
    }
    if (d_tagged_pt_idx == -1)
    {
        TBOX_ERROR("ERROR:: StructureParameters::StructureParameters( ) "
                   << "\n"
                   << "Could not tag Lagrangian point for this structure \n\n"
                   << std::endl);
    }

    return;

} // StructureParameters

ConstraintIBKinematics::ConstraintIBKinematics(const std::string& object_name,
                                               Pointer<Database> input_db,
                                               LDataManager* l_data_manager,
                                               bool register_for_restart)
    : d_struct_param(input_db, l_data_manager)
{
    // Set the object name and register it with the restart manager.
    d_object_name = object_name;
    d_registered_for_restart = false;
    if (register_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
        d_registered_for_restart = true;
    }

    return;

} // ConstraintIBKinematics

ConstraintIBKinematics::~ConstraintIBKinematics()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;

} // ~ConstraintIBKinematics

void ConstraintIBKinematics::putToDatabase(Pointer<Database> /*db*/)
{
    // intentionally left blank
    return;
}

} // IBAMR

//////////////////////////////////////////////////////////////////////////////
