// Filename: LDataManager.I
// Created on 10 Dec 2009 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/compiler_hints.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline std::vector<std::string>
LDataManager::getLagrangianStructureNames(
    const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif
    std::vector<std::string> ret_val;
    for (std::map<int,std::string>::const_iterator cit(d_strct_id_to_strct_name_map[level_number].begin());
         cit != d_strct_id_to_strct_name_map[level_number].end(); ++cit)
    {
        ret_val.push_back(cit->second);
    }
    return ret_val;
}// getLagrangianStructureNames

inline std::vector<int>
LDataManager::getLagrangianStructureIDs(
    const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif
    std::vector<int> ret_val;
    for (std::map<std::string,int>::const_iterator cit(d_strct_name_to_strct_id_map[level_number].begin());
         cit != d_strct_name_to_strct_id_map[level_number].end(); ++cit)
    {
        ret_val.push_back(cit->second);
    }
    return ret_val;
}// getLagrangianStructureIDs

inline int
LDataManager::getLagrangianStructureID(
    const int lagrangian_index,
    const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif
    std::map<int,int>::const_iterator cit = d_last_lag_idx_to_strct_id_map[level_number].lower_bound(lagrangian_index);
    if (UNLIKELY(cit == d_last_lag_idx_to_strct_id_map[level_number].end())) return -1;
    const int strct_id = cit->second;
#ifdef DEBUG_CHECK_ASSERTIONS
    const std::pair<int,int>& idx_range = getLagrangianStructureIndexRange(strct_id, level_number);
    TBOX_ASSERT(idx_range.first <= lagrangian_index && lagrangian_index < idx_range.second);
#endif
    return strct_id;
}// getLagrangianStructureID

inline int
LDataManager::getLagrangianStructureID(
    const std::string& structure_name,
    const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif
    std::map<std::string,int>::const_iterator cit = d_strct_name_to_strct_id_map[level_number].find(structure_name);
    if (UNLIKELY(cit == d_strct_name_to_strct_id_map[level_number].end())) return -1;
    return cit->second;
}// getLagrangianStructureID

inline std::string
LDataManager::getLagrangianStructureName(
    const int structure_id,
    const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif
    std::map<int,std::string>::const_iterator cit = d_strct_id_to_strct_name_map[level_number].find(structure_id);
    if (UNLIKELY(cit == d_strct_id_to_strct_name_map[level_number].end())) return std::string("UNKNOWN");
    return cit->second;
}// getLagrangianStructureName

inline std::pair<int,int>
LDataManager::getLagrangianStructureIndexRange(
    const int structure_id,
    const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif
    std::map<int,std::pair<int,int> >::const_iterator cit = d_strct_id_to_lag_idx_range_map[level_number].find(structure_id);
    if (UNLIKELY(cit == d_strct_id_to_lag_idx_range_map[level_number].end())) return std::make_pair(-1,-1);
    return cit->second;
}// getLagrangianStructureIndexRange

inline bool
LDataManager::getLagrangianStructureIsActivated(
    const int structure_id,
    const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif
    std::set<int>::const_iterator cit = d_inactive_strcts[level_number].getSet().find(structure_id);
    return (cit == d_inactive_strcts[level_number].getSet().end());
}// getLagrangianStructureIsActivated

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
