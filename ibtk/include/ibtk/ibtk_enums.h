// Filename: ibtk_enums.h
// Created on 13 Aug 2011 by Boyce Griffith
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

#ifndef included_ibtk_enums
#define included_ibtk_enums

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <cstring>

#include "tbox/Utilities.h"

/////////////////////////////// ENUM DEFINITIONS /////////////////////////////

namespace IBTK
{
/*!
 * \brief Routine for converting strings to enums.
 */
template <typename T>
inline T
string_to_enum(const std::string& /*val*/)
{
    TBOX_ERROR("UNSUPPORTED ENUM TYPE\n");
    return -1;
} // string_to_enum

/*!
 * \brief Routine for converting enums to strings.
 */
template <typename T>
inline std::string enum_to_string(T /*val*/)
{
    TBOX_ERROR("UNSUPPORTED ENUM TYPE\n");
    return "UNKNOWN";
} // enum_to_string

/*!
 * \brief Enumerated type for different multigrid cycle types.
 */
enum MGCycleType
{
    V_CYCLE,
    W_CYCLE,
    F_CYCLE,
    UNKNOWN_MG_CYCLE_TYPE = -1
};

template <>
inline MGCycleType
string_to_enum<MGCycleType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "V") == 0) return V_CYCLE;
    if (strcasecmp(val.c_str(), "V_CYCLE") == 0) return V_CYCLE;
    if (strcasecmp(val.c_str(), "V-CYCLE") == 0) return V_CYCLE;
    if (strcasecmp(val.c_str(), "W") == 0) return W_CYCLE;
    if (strcasecmp(val.c_str(), "W_CYCLE") == 0) return W_CYCLE;
    if (strcasecmp(val.c_str(), "W-CYCLE") == 0) return W_CYCLE;
    if (strcasecmp(val.c_str(), "F") == 0) return F_CYCLE;
    if (strcasecmp(val.c_str(), "F_CYCLE") == 0) return F_CYCLE;
    if (strcasecmp(val.c_str(), "F-CYCLE") == 0) return F_CYCLE;
    return UNKNOWN_MG_CYCLE_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<MGCycleType>(MGCycleType val)
{
    if (val == V_CYCLE) return "V_CYCLE";
    if (val == W_CYCLE) return "W_CYCLE";
    if (val == F_CYCLE) return "F_CYCLE";
    return "UNKNOWN_MG_CYCLE_TYPE";
} // enum_to_string

/*!
 * \brief Enumerated type for different regridding modes.
 */
enum RegridMode
{
    STANDARD,
    AGGRESSIVE,
    UNKNOWN_REGRID_MODE = -1
};

template <>
inline RegridMode
string_to_enum<RegridMode>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "STANDARD") == 0) return STANDARD;
    if (strcasecmp(val.c_str(), "AGGRESSIVE") == 0) return AGGRESSIVE;
    return UNKNOWN_REGRID_MODE;
} // string_to_enum

template <>
inline std::string
enum_to_string<RegridMode>(RegridMode val)
{
    if (val == STANDARD) return "STANDARD";
    if (val == AGGRESSIVE) return "AGGRESSIVE";
    return "UNKNOWN_REGRID_MODE";
} // enum_to_string

/*!
 * \brief Enumerated type for different standard data contexts.
 */
enum VariableContextType
{
    CURRENT_DATA,
    NEW_DATA,
    SCRATCH_DATA,
    UNKNOWN_VARIABLE_CONTEXT_TYPE = -1
};

template <>
inline VariableContextType
string_to_enum<VariableContextType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "CURRENT_DATA") == 0) return CURRENT_DATA;
    if (strcasecmp(val.c_str(), "NEW_DATA") == 0) return NEW_DATA;
    if (strcasecmp(val.c_str(), "SCRATCH_DATA") == 0) return SCRATCH_DATA;
    return UNKNOWN_VARIABLE_CONTEXT_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<VariableContextType>(VariableContextType val)
{
    if (val == CURRENT_DATA) return "CURRENT_DATA";
    if (val == NEW_DATA) return "NEW_DATA";
    if (val == SCRATCH_DATA) return "SCRATCH_DATA";
    return "UNKNOWN_VARIABLE_CONTEXT_TYPE";
} // enum_to_string

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_ibtk_enums
