// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2021 by the IBAMR developers
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

#ifndef included_IBAMR_ibamr_enums
#define included_IBAMR_ibamr_enums

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "tbox/Utilities.h"

#include <cstring>

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
    UPWIND = 10,
    CUI = 11,
    FBICS = 12,
    MGAMMA = 13,
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
    if (strcasecmp(val.c_str(), "UPWIND") == 0) return UPWIND;
    if (strcasecmp(val.c_str(), "CUI") == 0) return CUI;
    if (strcasecmp(val.c_str(), "FBICS") == 0) return FBICS;
    if (strcasecmp(val.c_str(), "MGAMMA") == 0) return MGAMMA;
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
    if (val == UPWIND) return "UPWIND";
    if (val == CUI) return "CUI";
    if (val == FBICS) return "FBICS";
    if (val == MGAMMA) return "MGAMMA";
    return "UNKNOWN_LIMITER_TYPE";
} // enum_to_string

/*!
 * \brief Enumerated type for different level set methods
 */
enum LevelSetType
{
    FAST_SWEEPING = 1,
    RELAXATION = 2,
    UNKNOWN_LEVEL_SET_TYPE = -1
};

template <>
inline LevelSetType
string_to_enum<LevelSetType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "FAST_SWEEPING") == 0) return FAST_SWEEPING;
    if (strcasecmp(val.c_str(), "RELAXATION") == 0) return RELAXATION;
    return UNKNOWN_LEVEL_SET_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<LevelSetType>(LevelSetType val)
{
    if (val == FAST_SWEEPING) return "FAST_SWEEPING";
    if (val == RELAXATION) return "RELAXATION";
    return "UNKNOWN_LEVEL_SET_TYPE";
} // enum_to_string

/*!
 * \brief Enumerated type for different level set order
 */
enum LevelSetOrder
{
    FIRST_ORDER_LS = 1,
    SECOND_ORDER_LS = 2,
    THIRD_ORDER_ENO_LS = 3,
    THIRD_ORDER_WENO_LS = 4,
    FIFTH_ORDER_WENO_LS = 5,
    UNKNOWN_LEVEL_SET_ORDER = -1
};

template <>
inline LevelSetOrder
string_to_enum<LevelSetOrder>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "FIRST_ORDER") == 0) return FIRST_ORDER_LS;
    if (strcasecmp(val.c_str(), "SECOND_ORDER") == 0) return SECOND_ORDER_LS;
    if (strcasecmp(val.c_str(), "THIRD_ORDER_ENO") == 0) return THIRD_ORDER_ENO_LS;
    if (strcasecmp(val.c_str(), "THIRD_ORDER_WENO") == 0) return THIRD_ORDER_WENO_LS;
    if (strcasecmp(val.c_str(), "FIFTH_ORDER_WENO") == 0) return FIFTH_ORDER_WENO_LS;
    return UNKNOWN_LEVEL_SET_ORDER;
} // string_to_enum

template <>
inline std::string
enum_to_string<LevelSetOrder>(LevelSetOrder val)
{
    if (val == FIRST_ORDER_LS) return "FIRST_ORDER";
    if (val == SECOND_ORDER_LS) return "SECOND_ORDER";
    if (val == THIRD_ORDER_ENO_LS) return "THIRD_ORDER_ENO";
    if (val == THIRD_ORDER_WENO_LS) return "THIRD_ORDER_WENO";
    if (val == FIFTH_ORDER_WENO_LS) return "FIFTH_ORDER_WENO";
    return "UNKNOWN_LEVEL_SET_ORDER";
} // enum_to_string

/*!
 * \brief Enumerated type for pressure formulations.
 */
enum PressureProjectionType
{
    CONSISTENT_PROJECTION,
    LUMPED_PROJECTION,
    STABILIZED_PROJECTION,
    UNKNOWN_PRESSURE_TYPE = -1
};

template <>
inline PressureProjectionType
string_to_enum<PressureProjectionType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "CONSISTENT_PROJECTION") == 0) return CONSISTENT_PROJECTION;
    if (strcasecmp(val.c_str(), "LUMPED_PROJECTION") == 0) return LUMPED_PROJECTION;
    if (strcasecmp(val.c_str(), "STABILIZED_PROJECTION") == 0) return STABILIZED_PROJECTION;
    return UNKNOWN_PRESSURE_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<PressureProjectionType>(PressureProjectionType val)
{
    if (val == CONSISTENT_PROJECTION) return "CONSISTENT_PROJECTION";
    if (val == LUMPED_PROJECTION) return "LUMPED_PROJECTION";
    if (val == STABILIZED_PROJECTION) return "STABILIZED_PROJECTION";
    return "UNKNOWN_PRESSURE_TYPE";
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
    SSPRK2,
    SSPRK3,
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
    if (strcasecmp(val.c_str(), "SSPRK1") == 0) return FORWARD_EULER;
    if (strcasecmp(val.c_str(), "SSPRK2") == 0) return SSPRK2;
    if (strcasecmp(val.c_str(), "SSPRK3") == 0) return SSPRK3;
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
    if (val == SSPRK2) return "SSPRK2";
    if (val == SSPRK3) return "SSPRK3";
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

/*!
 * \brief Enumerated type for different forms of dense mobility matrix.
 */
enum MobilityMatrixType
{
    READ_FROM_FILE,
    RPY,
    EMPIRICAL,
    UNKNOWN_MOBILITY_MATRIX_TYPE = -1
};

template <>
inline MobilityMatrixType
string_to_enum<MobilityMatrixType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "READ_FROM_FILE") == 0) return READ_FROM_FILE;
    if (strcasecmp(val.c_str(), "RPY") == 0) return RPY;
    if (strcasecmp(val.c_str(), "EMPIRICAL") == 0) return EMPIRICAL;
    return UNKNOWN_MOBILITY_MATRIX_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<MobilityMatrixType>(MobilityMatrixType val)
{
    if (val == READ_FROM_FILE) return "READ_FROM_FILE";
    if (val == RPY) return "RPY";
    if (val == EMPIRICAL) return "EMPIRICAL";
    return "UNKNOWN_MOBILITY_MATRIX_TYPE";
} // enum_to_string

/*!
 * \brief Enumerated type for different direct methods for dense mobility
 *  matrix inversion.
 */
enum MobilityMatrixInverseType
{
    LAPACK_CHOLESKY,
    LAPACK_LU,
    LAPACK_SVD,
    UNKNOWN_MOBILITY_MATRIX_INVERSE_TYPE = -1
};

template <>
inline MobilityMatrixInverseType
string_to_enum<MobilityMatrixInverseType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "LAPACK_CHOLESKY") == 0) return LAPACK_CHOLESKY;
    if (strcasecmp(val.c_str(), "LAPACK_LU") == 0) return LAPACK_LU;
    if (strcasecmp(val.c_str(), "LAPACK_SVD") == 0) return LAPACK_SVD;
    return UNKNOWN_MOBILITY_MATRIX_INVERSE_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<MobilityMatrixInverseType>(MobilityMatrixInverseType val)
{
    if (val == LAPACK_CHOLESKY) return "LAPACK_CHOLESKY";
    if (val == LAPACK_LU) return "LAPACK_LU";
    if (val == LAPACK_SVD) return "LAPACK_SVD";
    return "UNKNOWN_MOBILITY_MATRIX_INVERSE_TYPE";
} // enum_to_string

/*!
 * \brief Enumerated type for different possible libMesh partitioners.
 */
enum LibmeshPartitionerType
{
    LIBMESH_DEFAULT,
    SAMRAI_BOX,
    UNKNOWN_LIBMESH_PARTITIONER_TYPE = -1
};

template <>
inline LibmeshPartitionerType
string_to_enum<LibmeshPartitionerType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "LIBMESH_DEFAULT") == 0) return LIBMESH_DEFAULT;
    if (strcasecmp(val.c_str(), "SAMRAI_BOX") == 0) return SAMRAI_BOX;
    return UNKNOWN_LIBMESH_PARTITIONER_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<LibmeshPartitionerType>(LibmeshPartitionerType val)
{
    if (val == LIBMESH_DEFAULT) return "LIBMESH_DEFAULT";
    if (val == SAMRAI_BOX) return "SAMRAI_BOX";
    return "UNKNOWN_LIBMESH_PARTITIONER_TYPE";
} // enum_to_string

/*!
 * \brief Enumerated type for different evolution of SPD tensors.
 */
enum TensorEvolutionType
{
    STANDARD,
    SQUARE_ROOT,
    LOGARITHM,
    UNKNOWN_TENSOR_EVOLUTION_TYPE = -1
};

template <>
inline TensorEvolutionType
string_to_enum<TensorEvolutionType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "STANDARD") == 0) return STANDARD;
    if (strcasecmp(val.c_str(), "SQUARE_ROOT") == 0) return SQUARE_ROOT;
    if (strcasecmp(val.c_str(), "LOGARITHM") == 0) return LOGARITHM;
    return UNKNOWN_TENSOR_EVOLUTION_TYPE;
}

template <>
inline std::string
enum_to_string<TensorEvolutionType>(TensorEvolutionType val)
{
    if (val == STANDARD) return "STANDARD";
    if (val == SQUARE_ROOT) return "SQUARE_ROOT";
    if (val == LOGARITHM) return "LOGARITHM";
    return "UNKNOWN_TENSOR_EVOLUTION_TYPE";
}

/*!
 * \brief Enumerated type for different Advection-Diffusion Brinkman penalization boundary conditions
 */
enum AdvDiffBrinkmanPenalizationBcType
{
    DIRICHLET = 1,
    NEUMANN = 2,
    ROBIN = 3,
    UNKNOWN_BRINKMAN_BC_TYPE = -1
};

template <>
inline AdvDiffBrinkmanPenalizationBcType
string_to_enum<AdvDiffBrinkmanPenalizationBcType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "DIRICHLET") == 0) return DIRICHLET;
    if (strcasecmp(val.c_str(), "NEUMANN") == 0) return NEUMANN;
    if (strcasecmp(val.c_str(), "ROBIN") == 0) return ROBIN;
    return UNKNOWN_BRINKMAN_BC_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<AdvDiffBrinkmanPenalizationBcType>(AdvDiffBrinkmanPenalizationBcType val)
{
    if (val == DIRICHLET) return "DIRICHLET";
    if (val == NEUMANN) return "NEUMANN";
    if (val == ROBIN) return "ROBIN";
    return "UNKNOWN_BRINKMAN_BC_TYPE";
} // enum_to_string

/*!
 * \brief Enumerated type for different indicator functions.
 */
enum IndicatorFunctionType
{
    SMOOTH = 1,
    DISCONTINUOUS = 2,
    UNKNOWN_INDICATOR_FUNC_TYPE = -1
};

template <>
inline IndicatorFunctionType
string_to_enum<IndicatorFunctionType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "SMOOTH") == 0) return SMOOTH;
    if (strcasecmp(val.c_str(), "DISCONTINUOUS") == 0) return DISCONTINUOUS;
    return UNKNOWN_INDICATOR_FUNC_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<IndicatorFunctionType>(IndicatorFunctionType val)
{
    if (val == SMOOTH) return "SMOOTH";
    if (val == DISCONTINUOUS) return "DISCONTINUOUS";
    return "UNKNOWN_INDICATOR_FUNC_TYPE";
} // enum_to_string
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_ibamr_enums
