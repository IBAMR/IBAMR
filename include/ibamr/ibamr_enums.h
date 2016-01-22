// Filename: ibamr_enums.h
// Created on 03 Mar 2011 by Boyce Griffith
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

#ifndef included_ibamr_enums
#define included_ibamr_enums

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <cstring>

#include "tbox/Utilities.h"

/////////////////////////////// ENUM DEFINITIONS /////////////////////////////

namespace IBAMR
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
 * \brief Enumerated type for different convective differencing schemes.
 */
enum ConvectiveDifferencingType
{
    ADVECTIVE,
    CONSERVATIVE,
    SKEW_SYMMETRIC,
    UNKNOWN_CONVECTIVE_DIFFERENCING_TYPE = -1
};

template <>
inline ConvectiveDifferencingType
string_to_enum<ConvectiveDifferencingType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "ADVECTION") == 0) return ADVECTIVE;
    if (strcasecmp(val.c_str(), "ADVECTIVE") == 0) return ADVECTIVE;
    if (strcasecmp(val.c_str(), "CONSERVATION") == 0) return CONSERVATIVE;
    if (strcasecmp(val.c_str(), "CONSERVATIVE") == 0) return CONSERVATIVE;
    if (strcasecmp(val.c_str(), "DIVERGENCE") == 0) return CONSERVATIVE;
    if (strcasecmp(val.c_str(), "SKEW_SYMMETRIC") == 0) return SKEW_SYMMETRIC;
    return UNKNOWN_CONVECTIVE_DIFFERENCING_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<ConvectiveDifferencingType>(ConvectiveDifferencingType val)
{
    if (val == ADVECTIVE) return "ADVECTIVE";
    if (val == CONSERVATIVE) return "CONSERVATIVE";
    if (val == SKEW_SYMMETRIC) return "SKEW_SYMMETRIC";
    return "UNKNOWN_CONVECTIVE_DIFFERENCING_TYPE";
} // enum_to_string

/*!
 * \brief Enumerated type for different limiter types
 */
enum LimiterType
{
    CTU_ONLY = 1, // the assigned integer values are used in src/advect/fortran/limitertypes.i
    MINMOD_LIMITED = 2,
    MC_LIMITED = 3,
    SUPERBEE_LIMITED = 4,
    MUSCL_LIMITED = 5,
    SECOND_ORDER = 6,
    FOURTH_ORDER = 7,
    PPM = 8,
    XSPPM7 = 9,
    UNKNOWN_LIMITER_TYPE = -1
};

template <>
inline LimiterType
string_to_enum<LimiterType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "CTU_ONLY") == 0) return CTU_ONLY;
    if (strcasecmp(val.c_str(), "MINMOD_LIMITED") == 0) return MINMOD_LIMITED;
    if (strcasecmp(val.c_str(), "MINMOD") == 0) return MINMOD_LIMITED;
    if (strcasecmp(val.c_str(), "MC_LIMITED") == 0) return MC_LIMITED;
    if (strcasecmp(val.c_str(), "MC") == 0) return MC_LIMITED;
    if (strcasecmp(val.c_str(), "SUPERBEE_LIMITED") == 0) return SUPERBEE_LIMITED;
    if (strcasecmp(val.c_str(), "SUPERBEE") == 0) return SUPERBEE_LIMITED;
    if (strcasecmp(val.c_str(), "MUSCL_LIMITED") == 0) return MUSCL_LIMITED;
    if (strcasecmp(val.c_str(), "MUSCL") == 0) return MUSCL_LIMITED;
    if (strcasecmp(val.c_str(), "SECOND_ORDER") == 0) return SECOND_ORDER;
    if (strcasecmp(val.c_str(), "FOURTH_ORDER") == 0) return FOURTH_ORDER;
    if (strcasecmp(val.c_str(), "PPM") == 0) return PPM;
    if (strcasecmp(val.c_str(), "XSPPM7") == 0) return XSPPM7;
    return UNKNOWN_LIMITER_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<LimiterType>(LimiterType val)
{
    if (val == CTU_ONLY) return "CTU_ONLY";
    if (val == MINMOD_LIMITED) return "MINMOD_LIMITED";
    if (val == MC_LIMITED) return "MC_LIMITED";
    if (val == SUPERBEE_LIMITED) return "SUPERBEE_LIMITED";
    if (val == MUSCL_LIMITED) return "MUSCL_LIMITED";
    if (val == SECOND_ORDER) return "SECOND_ORDER";
    if (val == FOURTH_ORDER) return "FOURTH_ORDER";
    if (val == PPM) return "PPM";
    if (val == XSPPM7) return "XSPPM7";
    return "UNKNOWN_LIMITER_TYPE";
} // enum_to_string

/*!
 * \brief Enumerated type for different basic time stepping schemes.
 */
enum TimeSteppingType
{
    ADAMS_BASHFORTH,
    BACKWARD_EULER,
    FORWARD_EULER,
    MIDPOINT_RULE,
    TRAPEZOIDAL_RULE,
    UNKNOWN_TIME_STEPPING_TYPE = -1
};

template <>
inline TimeSteppingType
string_to_enum<TimeSteppingType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "ADAMS_BASHFORTH") == 0) return ADAMS_BASHFORTH;
    if (strcasecmp(val.c_str(), "BACKWARD_EULER") == 0) return BACKWARD_EULER;
    if (strcasecmp(val.c_str(), "FORWARD_EULER") == 0) return FORWARD_EULER;
    if (strcasecmp(val.c_str(), "MIDPOINT_RULE") == 0) return MIDPOINT_RULE;
    if (strcasecmp(val.c_str(), "TRAPEZOIDAL_RULE") == 0) return TRAPEZOIDAL_RULE;
    if (strcasecmp(val.c_str(), "CRANK_NICOLSON") == 0) return TRAPEZOIDAL_RULE;
    return UNKNOWN_TIME_STEPPING_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<TimeSteppingType>(TimeSteppingType val)
{
    if (val == ADAMS_BASHFORTH) return "ADAMS_BASHFORTH";
    if (val == BACKWARD_EULER) return "BACKWARD_EULER";
    if (val == FORWARD_EULER) return "FORWARD_EULER";
    if (val == MIDPOINT_RULE) return "MIDPOINT_RULE";
    if (val == TRAPEZOIDAL_RULE) return "TRAPEZOIDAL_RULE";
    return "UNKNOWN_TIME_STEPPING_TYPE";
} // enum_to_string

inline bool
is_multistep_time_stepping_type(TimeSteppingType val)
{
    switch (val)
    {
    case ADAMS_BASHFORTH:
        return true;
    case BACKWARD_EULER:
    case FORWARD_EULER:
    case MIDPOINT_RULE:
    case TRAPEZOIDAL_RULE:
        return false;
    default:
        TBOX_ERROR("is_multistep_time_stepping_type(): unknown time stepping type\n");
        return false;
    }
} // is_multistep_time_stepping_type

/*!
 * \brief Enumerated type for different types of traction boundary conditions.
 */
enum TractionBcType
{
    TRACTION,
    PSEUDO_TRACTION,
    UNKNOWN_TRACTION_BC_TYPE = -1
};

template <>
inline TractionBcType
string_to_enum<TractionBcType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "TRACTION") == 0) return TRACTION;
    if (strcasecmp(val.c_str(), "PSEUDO_TRACTION") == 0) return PSEUDO_TRACTION;
    return UNKNOWN_TRACTION_BC_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<TractionBcType>(TractionBcType val)
{
    if (val == TRACTION) return "TRACTION";
    if (val == PSEUDO_TRACTION) return "PSEUDO_TRACTION";
    return "UNKNOWN_TRACTION_BC_TYPE";
} // enum_to_string

/*!
 * \brief Enumerated type for different pressure update schemes for the
 * projection method.
 */
enum ProjectionMethodType
{
    PRESSURE_UPDATE,
    PRESSURE_INCREMENT,
    UNKNOWN_PROJECTION_METHOD_TYPE = -1
};

template <>
inline ProjectionMethodType
string_to_enum<ProjectionMethodType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "PRESSURE_UPDATE") == 0) return PRESSURE_UPDATE;
    if (strcasecmp(val.c_str(), "KIM_MOIN") == 0) return PRESSURE_UPDATE;
    if (strcasecmp(val.c_str(), "PRESSURE_INCREMENT") == 0) return PRESSURE_INCREMENT;
    if (strcasecmp(val.c_str(), "BCG") == 0) return PRESSURE_INCREMENT;
    if (strcasecmp(val.c_str(), "BELL_COLELLA_GLAZ") == 0) return PRESSURE_INCREMENT;
    return UNKNOWN_PROJECTION_METHOD_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<ProjectionMethodType>(ProjectionMethodType val)
{
    if (val == PRESSURE_UPDATE) return "PRESSURE_UPDATE";
    if (val == PRESSURE_INCREMENT) return "PRESSURE_INCREMENT";
    return "UNKNOWN_PROJECTION_METHOD_TYPE";
} // enum_to_string

/*!
 * \brief Enumerated type for different forms of the stochastic stress tensor.
 */
enum StochasticStressTensorType
{
    UNCORRELATED,
    SYMMETRIC,
    SYMMETRIC_TRACELESS,
    UNKNOWN_STOCHASTIC_STRESS_TENSOR_TYPE = -1
};

template <>
inline StochasticStressTensorType
string_to_enum<StochasticStressTensorType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "UNCORRELATED") == 0) return UNCORRELATED;
    if (strcasecmp(val.c_str(), "SYMMETRIC") == 0) return SYMMETRIC;
    if (strcasecmp(val.c_str(), "SYMMETRIC_TRACELESS") == 0) return SYMMETRIC_TRACELESS;
    return UNKNOWN_STOCHASTIC_STRESS_TENSOR_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<StochasticStressTensorType>(StochasticStressTensorType val)
{
    if (val == UNCORRELATED) return "UNCORRELATED";
    if (val == SYMMETRIC) return "SYMMETRIC";
    if (val == SYMMETRIC_TRACELESS) return "SYMMETRIC_TRACELESS";
    return "UNKNOWN_STOCHASTIC_STRESS_TENSOR_TYPE";
} // enum_to_string

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_ibamr_enums
