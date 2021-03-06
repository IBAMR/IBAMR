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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_ibtk_enums
#define included_IBTK_ibtk_enums

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "tbox/Utilities.h"

#include <cstring>

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
    F_CYCLE,
    FMG_CYCLE,
    V_CYCLE,
    W_CYCLE,
    UNKNOWN_MG_CYCLE_TYPE = -1
};

template <>
inline MGCycleType
string_to_enum<MGCycleType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "F") == 0) return F_CYCLE;
    if (strcasecmp(val.c_str(), "F_CYCLE") == 0) return F_CYCLE;
    if (strcasecmp(val.c_str(), "F-CYCLE") == 0) return F_CYCLE;
    if (strcasecmp(val.c_str(), "FMG") == 0) return FMG_CYCLE;
    if (strcasecmp(val.c_str(), "FMG_CYCLE") == 0) return FMG_CYCLE;
    if (strcasecmp(val.c_str(), "FMG-CYCLE") == 0) return FMG_CYCLE;
    if (strcasecmp(val.c_str(), "V") == 0) return V_CYCLE;
    if (strcasecmp(val.c_str(), "V_CYCLE") == 0) return V_CYCLE;
    if (strcasecmp(val.c_str(), "V-CYCLE") == 0) return V_CYCLE;
    if (strcasecmp(val.c_str(), "W") == 0) return W_CYCLE;
    if (strcasecmp(val.c_str(), "W_CYCLE") == 0) return W_CYCLE;
    if (strcasecmp(val.c_str(), "W-CYCLE") == 0) return W_CYCLE;
    return UNKNOWN_MG_CYCLE_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<MGCycleType>(MGCycleType val)
{
    if (val == F_CYCLE) return "F_CYCLE";
    if (val == FMG_CYCLE) return "FMG_CYCLE";
    if (val == V_CYCLE) return "V_CYCLE";
    if (val == W_CYCLE) return "W_CYCLE";
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

/*!
 * \brief Enumerated type for different interpolation types for
 * the material properties of the viscous solver.
 */
enum VCInterpType
{
    VC_AVERAGE_INTERP = 1,
    VC_HARMONIC_INTERP = 2,
    UNKNOWN_VC_INTERP_TYPE = -1
};

template <>
inline VCInterpType
string_to_enum<VCInterpType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "VC_AVERAGE_INTERP") == 0) return VC_AVERAGE_INTERP;
    if (strcasecmp(val.c_str(), "VC_HARMONIC_INTERP") == 0) return VC_HARMONIC_INTERP;
    return UNKNOWN_VC_INTERP_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<VCInterpType>(VCInterpType val)
{
    if (val == VC_AVERAGE_INTERP) return "VC_AVERAGE_INTERP";
    if (val == VC_HARMONIC_INTERP) return "VC_HARMONIC_INTERP";
    return "UNKNOWN_VC_INTERP_TYPE";
} // enum_to_string

enum NodeOutsidePatchCheckType
{
    NODE_OUTSIDE_PERMIT = 1,
    NODE_OUTSIDE_WARN = 2,
    NODE_OUTSIDE_ERROR = 3,
    UNKNOWN_NODE_OUTSIDE_PATCH_CHECK_TYPE = -1
};

template <>
inline NodeOutsidePatchCheckType
string_to_enum<NodeOutsidePatchCheckType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "NODE_OUTSIDE_PERMIT") == 0) return NODE_OUTSIDE_PERMIT;
    if (strcasecmp(val.c_str(), "NODE_OUTSIDE_WARN") == 0) return NODE_OUTSIDE_WARN;
    if (strcasecmp(val.c_str(), "NODE_OUTSIDE_ERROR") == 0) return NODE_OUTSIDE_ERROR;
    return UNKNOWN_NODE_OUTSIDE_PATCH_CHECK_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<NodeOutsidePatchCheckType>(NodeOutsidePatchCheckType val)
{
    if (val == NODE_OUTSIDE_PERMIT) return "NODE_OUTSIDE_PERMIT";
    if (val == NODE_OUTSIDE_WARN) return "NODE_OUTSIDE_WARN";
    if (val == NODE_OUTSIDE_ERROR) return "NODE_OUTSIDE_ERROR";
    return "UNKNOWN_NODE_OUTSIDE_PATCH_CHECK_TYPE";
} // enum_to_string

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_ibtk_enums
