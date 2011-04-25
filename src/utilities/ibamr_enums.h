// Filename: ibamr_enums.h
// Created on 03 Mar 2011 by Boyce Griffith
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

#ifndef included_ibamr_enums
#define included_ibamr_enums

/////////////////////////////// ENUM DEFINITIONS /////////////////////////////

#include <cstring>
#include <tbox/Utilities.h>

/*!
 * \brief Routine for converting strings to enums.
 */
template<typename T>
inline T
string_to_enum(
    const std::string& val)
{
    TBOX_ERROR("UNSUPPORTED ENUM TYPE\n");
    return -1;
}// string_to_enum

/*!
 * \brief Routine for converting enums to strings.
 */
template<typename T>
inline std::string
enum_to_string(
    const T& val)
{
    TBOX_ERROR("UNSUPPORTED ENUM TYPE\n");
    return "UNKNOWN";
}// enum_to_string

/*!
 * \brief Enumerated type for different convective differencing schemes.
 */
enum ConvectiveDifferencingType
{
    ADVECTIVE,
    CONSERVATIVE,
    SKEW_SYMMETRIC,
    UNKNOWN_CONVECTIVE_DIFFERENCING_TYPE=-1
};

template<>
inline ConvectiveDifferencingType
string_to_enum<ConvectiveDifferencingType>(
    const std::string& val)
{
    if (strcasecmp(val.c_str(), "ADVECTIVE"     ) == 0) return ADVECTIVE;
    if (strcasecmp(val.c_str(), "CONSERVATIVE"  ) == 0) return CONSERVATIVE;
    if (strcasecmp(val.c_str(), "DIVERGENCE"    ) == 0) return CONSERVATIVE;
    if (strcasecmp(val.c_str(), "SKEW_SYMMETRIC") == 0) return SKEW_SYMMETRIC;
    return UNKNOWN_CONVECTIVE_DIFFERENCING_TYPE;
}// string_to_enum

template<>
inline std::string
enum_to_string<ConvectiveDifferencingType>(
    const ConvectiveDifferencingType& val)
{
    if (val == ADVECTIVE     ) return "ADVECTIVE";
    if (val == CONSERVATIVE  ) return "CONSERVATIVE";
    if (val == SKEW_SYMMETRIC) return "SKEW_SYMMETRIC";
    return "UNKNOWN_CONVECTIVE_DIFFERENCING_TYPE";
}// enum_to_string

/*!
 * \brief Enumerated type for different pressure update schemes for the
 * projection method.
 */
enum ProjectionMethodType
{
    PRESSURE_UPDATE,
    PRESSURE_INCREMENT,
    UNKNOWN_PROJECTION_METHOD_TYPE=-1
};

template<>
inline ProjectionMethodType
string_to_enum<ProjectionMethodType>(
    const std::string& val)
{
    if (strcasecmp(val.c_str(), "PRESSURE_UPDATE"   ) == 0) return PRESSURE_UPDATE;
    if (strcasecmp(val.c_str(), "PRESSURE_INCREMENT") == 0) return PRESSURE_INCREMENT;
    return UNKNOWN_PROJECTION_METHOD_TYPE;
}// string_to_enum

template<>
inline std::string
enum_to_string<ProjectionMethodType>(
    const ProjectionMethodType& val)
{
    if (val == PRESSURE_UPDATE   ) return "PRESSURE_UPDATE";
    if (val == PRESSURE_INCREMENT) return "PRESSURE_INCREMENT";
    return "UNKNOWN_PROJECTION_METHOD_TYPE";
}// enum_to_string

/*!
 * \brief Enumerated type for different regridding modes.
 */
enum RegridMode
{
    STANDARD,
    AGGRESSIVE,
    UNKNOWN_REGRID_MODE=-1
};

template<>
inline RegridMode
string_to_enum<RegridMode>(
    const std::string& val)
{
    if (strcasecmp(val.c_str(), "STANDARD"  ) == 0) return STANDARD;
    if (strcasecmp(val.c_str(), "AGGRESSIVE") == 0) return AGGRESSIVE;
    return UNKNOWN_REGRID_MODE;
}// string_to_enum

template<>
inline std::string
enum_to_string<RegridMode>(
    const RegridMode& val)
{
    if (val == STANDARD  ) return "STANDARD";
    if (val == AGGRESSIVE) return "AGGRESSIVE";
    return "UNKNOWN_REGRID_MODE";
}// enum_to_string

/*!
 * \brief Enumerated type for different Stokes preconditioners.
 */
enum StokesPreconditionerType
{
    NONE,
    PROJECTION_METHOD,
    VANKA,
    BLOCK_FACTORIZATION,
    UNKNOWN_STOKES_PRECONDITIONER_TYPE=-1
};

template<>
inline StokesPreconditionerType
string_to_enum<StokesPreconditionerType>(
    const std::string& val)
{
    if (strcasecmp(val.c_str(), "NONE"               ) == 0) return NONE;
    if (strcasecmp(val.c_str(), "PROJECTION_METHOD"  ) == 0) return PROJECTION_METHOD;
    if (strcasecmp(val.c_str(), "VANKA"              ) == 0) return VANKA;
    if (strcasecmp(val.c_str(), "BLOCK_FACTORIZATION") == 0) return BLOCK_FACTORIZATION;
    return UNKNOWN_STOKES_PRECONDITIONER_TYPE;
}// string_to_enum

template<>
inline std::string
enum_to_string<StokesPreconditionerType>(
    const StokesPreconditionerType& val)
{
    if (val == NONE               ) return "NONE";
    if (val == PROJECTION_METHOD  ) return "PROJECTION_METHOD";
    if (val == VANKA              ) return "VANKA";
    if (val == BLOCK_FACTORIZATION) return "BLOCK_FACTORIZATION";
    return "UNKNOWN_STOKES_PRECONDITIONER_TYPE";
}// enum_to_string

/*!
 * \brief Enumerated type for different viscous timestepping methods.
 */
enum ViscousTimesteppingType
{
    BACKWARD_EULER,
    CRANK_NICOLSON,
    TGA,
    UNKNOWN_VISCOUS_TIMESTEPPING_METHOD=-1
};

template<>
inline ViscousTimesteppingType
string_to_enum<ViscousTimesteppingType>(
    const std::string& val)
{
    if (strcasecmp(val.c_str(), "BACKWARD_EULER") == 0) return BACKWARD_EULER;
    if (strcasecmp(val.c_str(), "CRANK_NICOLSON") == 0) return CRANK_NICOLSON;
    if (strcasecmp(val.c_str(), "TGA"           ) == 0) return TGA;
    return UNKNOWN_VISCOUS_TIMESTEPPING_METHOD;
}// string_to_enum

template<>
inline std::string
enum_to_string<ViscousTimesteppingType>(
    const ViscousTimesteppingType& val)
{
    if (val == BACKWARD_EULER) return "BACKWARD_EULER";
    if (val == CRANK_NICOLSON) return "CRANK_NICOLSON";
    if (val == TGA           ) return "TGA";
    return "UNKNOWN_VISCOUS_TIMESTEPPING_METHOD";
}// enum_to_string

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_ibamr_enums
