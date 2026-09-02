// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#include <tbox/Utilities.h>

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
inline std::string
enum_to_string(T /*val*/)
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

enum class TimePoint
{
    CURRENT_TIME,
    HALF_TIME,
    NEW_TIME,
    UNKNOWN_TIME
};

template <>
inline TimePoint
string_to_enum<TimePoint>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "CURRENT_TIME") == 0) return TimePoint::CURRENT_TIME;
    if (strcasecmp(val.c_str(), "HALF_TIME") == 0) return TimePoint::HALF_TIME;
    if (strcasecmp(val.c_str(), "NEW_TIME") == 0) return TimePoint::NEW_TIME;
    return TimePoint::UNKNOWN_TIME;
}

template <>
inline std::string
enum_to_string<TimePoint>(TimePoint val)
{
    if (val == TimePoint::CURRENT_TIME) return "CURRENT_TIME";
    if (val == TimePoint::HALF_TIME) return "HALF_TIME";
    if (val == TimePoint::NEW_TIME) return "NEW_TIME";
    return "UNKNOWN_TIME_POINT";
}

/*!
 * \brief Semantic identities for IB interpolation and spreading kernels.
 *
 * For COMPOSITE_BSPLINE_XY, X is the B-spline order in the face-normal
 * direction and Y is the order in the face-tangential direction(s). This
 * catalog identifies kernels independently of the capabilities and storage
 * requirements of a particular interpolation or spreading backend.
 *
 * PIECEWISE_CONSTANT, PIECEWISE_LINEAR, and DISCONTINUOUS_LINEAR are accepted
 * as input aliases for BSPLINE_1, BSPLINE_2, and COMPOSITE_BSPLINE_21,
 * respectively. String conversion always returns the canonical name.
 */
enum class IBKernelType
{
    BSPLINE_1,
    BSPLINE_2,
    BSPLINE_3,
    BSPLINE_4,
    BSPLINE_5,
    BSPLINE_6,
    COMPOSITE_BSPLINE_12,
    COMPOSITE_BSPLINE_21,
    COMPOSITE_BSPLINE_23,
    COMPOSITE_BSPLINE_32,
    COMPOSITE_BSPLINE_34,
    COMPOSITE_BSPLINE_43,
    COMPOSITE_BSPLINE_45,
    COMPOSITE_BSPLINE_54,
    COMPOSITE_BSPLINE_56,
    COMPOSITE_BSPLINE_65,
    PIECEWISE_CUBIC,
    IB_3,
    IB_4,
    IB_4_W8,
    IB_5,
    IB_6,
    USER_DEFINED,
    UNKNOWN
};

template <>
inline IBKernelType
string_to_enum<IBKernelType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "BSPLINE_1") == 0 || strcasecmp(val.c_str(), "PIECEWISE_CONSTANT") == 0)
        return IBKernelType::BSPLINE_1;
    if (strcasecmp(val.c_str(), "BSPLINE_2") == 0 || strcasecmp(val.c_str(), "PIECEWISE_LINEAR") == 0)
        return IBKernelType::BSPLINE_2;
    if (strcasecmp(val.c_str(), "BSPLINE_3") == 0) return IBKernelType::BSPLINE_3;
    if (strcasecmp(val.c_str(), "BSPLINE_4") == 0) return IBKernelType::BSPLINE_4;
    if (strcasecmp(val.c_str(), "BSPLINE_5") == 0) return IBKernelType::BSPLINE_5;
    if (strcasecmp(val.c_str(), "BSPLINE_6") == 0) return IBKernelType::BSPLINE_6;
    if (strcasecmp(val.c_str(), "COMPOSITE_BSPLINE_12") == 0) return IBKernelType::COMPOSITE_BSPLINE_12;
    if (strcasecmp(val.c_str(), "COMPOSITE_BSPLINE_21") == 0 || strcasecmp(val.c_str(), "DISCONTINUOUS_LINEAR") == 0)
        return IBKernelType::COMPOSITE_BSPLINE_21;
    if (strcasecmp(val.c_str(), "COMPOSITE_BSPLINE_23") == 0) return IBKernelType::COMPOSITE_BSPLINE_23;
    if (strcasecmp(val.c_str(), "COMPOSITE_BSPLINE_32") == 0) return IBKernelType::COMPOSITE_BSPLINE_32;
    if (strcasecmp(val.c_str(), "COMPOSITE_BSPLINE_34") == 0) return IBKernelType::COMPOSITE_BSPLINE_34;
    if (strcasecmp(val.c_str(), "COMPOSITE_BSPLINE_43") == 0) return IBKernelType::COMPOSITE_BSPLINE_43;
    if (strcasecmp(val.c_str(), "COMPOSITE_BSPLINE_45") == 0) return IBKernelType::COMPOSITE_BSPLINE_45;
    if (strcasecmp(val.c_str(), "COMPOSITE_BSPLINE_54") == 0) return IBKernelType::COMPOSITE_BSPLINE_54;
    if (strcasecmp(val.c_str(), "COMPOSITE_BSPLINE_56") == 0) return IBKernelType::COMPOSITE_BSPLINE_56;
    if (strcasecmp(val.c_str(), "COMPOSITE_BSPLINE_65") == 0) return IBKernelType::COMPOSITE_BSPLINE_65;
    if (strcasecmp(val.c_str(), "PIECEWISE_CUBIC") == 0) return IBKernelType::PIECEWISE_CUBIC;
    if (strcasecmp(val.c_str(), "IB_3") == 0) return IBKernelType::IB_3;
    if (strcasecmp(val.c_str(), "IB_4") == 0) return IBKernelType::IB_4;
    if (strcasecmp(val.c_str(), "IB_4_W8") == 0) return IBKernelType::IB_4_W8;
    if (strcasecmp(val.c_str(), "IB_5") == 0) return IBKernelType::IB_5;
    if (strcasecmp(val.c_str(), "IB_6") == 0) return IBKernelType::IB_6;
    if (strcasecmp(val.c_str(), "USER_DEFINED") == 0) return IBKernelType::USER_DEFINED;
    return IBKernelType::UNKNOWN;
}

template <>
inline std::string
enum_to_string<IBKernelType>(IBKernelType val)
{
    if (val == IBKernelType::BSPLINE_1) return "BSPLINE_1";
    if (val == IBKernelType::BSPLINE_2) return "BSPLINE_2";
    if (val == IBKernelType::BSPLINE_3) return "BSPLINE_3";
    if (val == IBKernelType::BSPLINE_4) return "BSPLINE_4";
    if (val == IBKernelType::BSPLINE_5) return "BSPLINE_5";
    if (val == IBKernelType::BSPLINE_6) return "BSPLINE_6";
    if (val == IBKernelType::COMPOSITE_BSPLINE_12) return "COMPOSITE_BSPLINE_12";
    if (val == IBKernelType::COMPOSITE_BSPLINE_21) return "COMPOSITE_BSPLINE_21";
    if (val == IBKernelType::COMPOSITE_BSPLINE_23) return "COMPOSITE_BSPLINE_23";
    if (val == IBKernelType::COMPOSITE_BSPLINE_32) return "COMPOSITE_BSPLINE_32";
    if (val == IBKernelType::COMPOSITE_BSPLINE_34) return "COMPOSITE_BSPLINE_34";
    if (val == IBKernelType::COMPOSITE_BSPLINE_43) return "COMPOSITE_BSPLINE_43";
    if (val == IBKernelType::COMPOSITE_BSPLINE_45) return "COMPOSITE_BSPLINE_45";
    if (val == IBKernelType::COMPOSITE_BSPLINE_54) return "COMPOSITE_BSPLINE_54";
    if (val == IBKernelType::COMPOSITE_BSPLINE_56) return "COMPOSITE_BSPLINE_56";
    if (val == IBKernelType::COMPOSITE_BSPLINE_65) return "COMPOSITE_BSPLINE_65";
    if (val == IBKernelType::PIECEWISE_CUBIC) return "PIECEWISE_CUBIC";
    if (val == IBKernelType::IB_3) return "IB_3";
    if (val == IBKernelType::IB_4) return "IB_4";
    if (val == IBKernelType::IB_4_W8) return "IB_4_W8";
    if (val == IBKernelType::IB_5) return "IB_5";
    if (val == IBKernelType::IB_6) return "IB_6";
    if (val == IBKernelType::USER_DEFINED) return "USER_DEFINED";
    return "UNKNOWN_IB_KERNEL_TYPE";
}

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_ibtk_enums
