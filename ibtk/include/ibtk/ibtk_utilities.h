// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2024 by the IBAMR developers
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

#ifndef included_IBTK_ibtk_utilities
#define included_IBTK_ibtk_utilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "PatchHierarchy.h"
#include "SAMRAIVectorReal.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include "Eigen/Core" // IWYU pragma: export
#include "Eigen/StdVector"
IBTK_ENABLE_EXTRA_WARNINGS

#include <algorithm>
#include <array>
#include <utility>

/////////////////////////////// MACRO DEFINITIONS ////////////////////////////

#define IBTK_BIT_SET(bitfield, b) ((bitfield) |= (1 << (b)))
#define IBTK_BIT_CLEAR(bitfield, b) ((bitfield) &= ~(1 << (b)))
#define IBTK_BIT_TOGGLE(bitfield, b) ((bitfield) ^= (1 << (b)))
#define IBTK_BIT_CHECK(bitfield, b) ((bitfield) & (1 << (b)))

#define IBTK_DO_ONCE(task)                                                                                             \
    do                                                                                                                 \
    {                                                                                                                  \
        static bool done = false;                                                                                      \
        if (done == false)                                                                                             \
        {                                                                                                              \
            task;                                                                                                      \
            done = true;                                                                                               \
        }                                                                                                              \
    } while (0);

/**
 * Define the deprecation macro when the compiler specifies it is available
 */
#if defined(__has_cpp_attribute)
#if __has_cpp_attribute(deprecated)
#define IBTK_DEPRECATED(msg) [[deprecated(msg)]]
#else
#define IBTK_DEPRECATED(msg)
#endif
#else
#define IBTK_DEPRECATED(msg)
#endif

#define IBTK_DEPRECATED_CLASS1(deprecated_class_name)                                                                  \
    IBTK_DO_ONCE({                                                                                                     \
        SAMRAI::tbox::pout << "WARNING: class " << deprecated_class_name                                               \
                           << " is deprecated and may be removed in the future." << std::endl;                         \
    });

#define IBTK_DEPRECATED_CLASS2(deprecated_class_name, new_class_name)                                                  \
    IBTK_DO_ONCE({                                                                                                     \
        SAMRAI::tbox::pout << "WARNING: class " << deprecated_class_name                                               \
                           << " is deprecated and may be removed in the future.\n"                                     \
                           << "Please update your code to use class " << new_class_name << "." << std::endl;           \
    });

#define IBTK_DEPRECATED_FUNCTION1(deprecated_function_name)                                                            \
    IBTK_DO_ONCE({                                                                                                     \
        SAMRAI::tbox::pout << "WARNING: function " << deprecated_function_name                                         \
                           << " is deprecated and may be removed in the future." << std::endl;                         \
    });

#define IBTK_DEPRECATED_FUNCTION2(deprecated_function_name, new_function_name)                                         \
    IBTK_DO_ONCE({                                                                                                     \
        SAMRAI::tbox::pout << "WARNING: function " << deprecated_function_name                                         \
                           << " is deprecated and may be removed in the future.\n"                                     \
                           << "Please update your code to use function " << new_function_name << "." << std::endl;     \
    });

#define IBTK_DEPRECATED_MEMBER_FUNCTION1(class_name, deprecated_method_name)                                           \
    IBTK_DO_ONCE({                                                                                                     \
        SAMRAI::tbox::pout << "WARNING: member function " << deprecated_method_name << " of class " << class_name      \
                           << " is deprecated and may be removed in the future." << std::endl;                         \
    });

#define IBTK_DEPRECATED_MEMBER_FUNCTION2(class_name, deprecated_method_name, new_method_name)                          \
    IBTK_DO_ONCE({                                                                                                     \
        SAMRAI::tbox::pout << "WARNING: member function " << deprecated_method_name << " of class " << class_name      \
                           << " is deprecated and may be removed in the future.\n"                                     \
                           << "Please update your code to use member function " << new_method_name << "."              \
                           << std::endl;                                                                               \
    });

#define IBTK_DEPRECATED_FUNCTIONALITY(message)                                                                         \
    IBTK_DO_ONCE({ SAMRAI::tbox::pout << "WARNING: " << message << "." << std::endl; });

namespace IBTK
{
static const bool ENABLE_TIMERS = true;
}

#define IBTK_TIMER_START(timer)                                                                                        \
    do                                                                                                                 \
    {                                                                                                                  \
        if (IBTK::ENABLE_TIMERS) timer->start();                                                                       \
    } while (0);

#define IBTK_TIMER_STOP(timer)                                                                                         \
    do                                                                                                                 \
    {                                                                                                                  \
        if (IBTK::ENABLE_TIMERS) timer->stop();                                                                        \
    } while (0);

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////

namespace IBTK
{
/*!
 * Check whether the relative difference between a and b are within the threshold eps.
 *
 * \note This function should be used with caution to check numbers close to zero. In this case, consider using the
 * abs_equal_eps function.
 */
inline bool
rel_equal_eps(double a, double b, double eps = std::sqrt(std::numeric_limits<double>::epsilon()))
{
    return (a == b) || (std::abs(a - b) / std::max(std::abs(a), std::abs(b))) < eps;
}

/*!
 * \brief Check whether the absolute difference between a and b are within the threshold eps.
 *
 * \note This function should be used with caution to check numbers that have large magnitudes. In these cases, consider
 * using the rel_equal_eps function.
 */
inline bool
abs_equal_eps(double a, double b, double eps = std::sqrt(std::numeric_limits<double>::epsilon()))
{
    return std::abs(a - b) < eps;
}

inline std::string
get_data_time_str(const double data_time, const double current_time, const double new_time)
{
    const double half_time = 0.5 * (current_time + new_time);
    if (rel_equal_eps(data_time, current_time))
    {
        return "current";
    }
    else if (rel_equal_eps(data_time, half_time))
    {
        return "half";
    }
    else if (rel_equal_eps(data_time, new_time))
    {
        return "new";
    }
    else
    {
        return "unknown";
    }
}

/*!
 * Get the smallest cell width on the specified level. This operation is
 * collective.
 */
double get_min_patch_dx(const SAMRAI::hier::PatchLevel<NDIM>& patch_level);

template <class T, unsigned N>
inline std::array<T, N>
array_constant(const T& v)
{
    std::array<T, N> arr;
    std::fill(arr.begin(), arr.end(), v);
    return arr;
} // array_constant

template <class T, unsigned N>
inline std::array<T, N>
array_one()
{
    std::array<T, N> arr;
    std::fill(arr.begin(), arr.end(), 1);
    return arr;
} // array_one

template <class T, unsigned N>
inline std::array<T, N>
array_zero()
{
    std::array<T, N> arr;
    std::fill(arr.begin(), arr.end(), 0);
    return arr;
} // array_zero

inline bool
level_can_be_refined(int level_number, int max_levels)
{
    const int finest_level_number = max_levels - 1;
    return level_number < finest_level_number;
}

/*!
 * Convert a Voigt notation index to the corresponding symmetric tensor index. This function only returns the upper
 * triangular index.
 */
inline std::pair<int, int>
voigt_to_tensor_idx(const int k)
{
    if (k < NDIM)
        return std::make_pair(k, k);
    else
#if (NDIM == 2)
        return std::make_pair(0, 1);
#endif
#if (NDIM == 3)
    return std::make_pair(k > 3 ? 0 : 1, k > 4 ? 1 : 2);
#endif
}

/*!
 * Convert a symmetric tensor index to the corresponding Voigt notation index.
 */
inline int
tensor_idx_to_voigt(const std::pair<int, int>& idx)
{
    if (idx.first == idx.second)
        return idx.first;
    else
        return 3 * NDIM - 3 - idx.first - idx.second;
}

/*!
 * Smooth heaviside function.
 */
inline double
smooth_heaviside(const double& phi, const double& alpha)
{
    double Hphi = 1.0;
    if (phi < -alpha)
    {
        Hphi = 0.0;
    }
    else if (std::abs(phi) <= alpha)
    {
        Hphi = 0.5 + 0.5 * phi / alpha + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi / alpha);
    }
    return Hphi;
}

/*!
 * Smooth delta function.
 */
inline double
smooth_delta(const double& phi, const double& alpha)
{
    double delta = 0.0;
    if (std::abs(phi) <= alpha)
    {
        delta = 0.5 / alpha + 1.0 / (2.0 * alpha) * std::cos(M_PI * phi / alpha);
    }
    return delta;
}

/*!
 * Discontinuous heaviside function.
 */
inline double
discontinuous_heaviside(const double& phi)
{
    double Hphi = 1.0;
    if (phi < 0.0)
    {
        Hphi = 0.0;
    }
    else if (phi == 0.0)
    {
        Hphi = 0.5;
    }
    return Hphi;
}

/*!
 * Eigen types have special alignment requirements and require a specific
 * memory allocator. This is a convenience type alias for a
 * <code>std::vector</code> with the correct allocator.
 */
template <typename T>
using EigenAlignedVector = std::vector<T, Eigen::aligned_allocator<T> >;

using Matrix2d = Eigen::Matrix<double, 2, 2>;
using Vector2d = Eigen::Matrix<double, 2, 1>;
using ColumnVector2d = Eigen::Matrix<double, 2, 1>;
using RowVector2d = Eigen::Matrix<double, 1, 2>;

using Matrix3d = Eigen::Matrix<double, 3, 3>;
using Vector3d = Eigen::Matrix<double, 3, 1>;
using ColumnVector3d = Eigen::Matrix<double, 3, 1>;
using RowVector3d = Eigen::Matrix<double, 1, 3>;

using MatrixNd = Eigen::Matrix<double, NDIM, NDIM>;
using VectorNd = Eigen::Matrix<double, NDIM, 1>;
using ColumnVectorNd = Eigen::Matrix<double, NDIM, 1>;
using RowVectorNd = Eigen::Matrix<double, 1, NDIM>;

using MatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using VectorXd = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using ColumnVectorXd = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using RowVectorXd = Eigen::Matrix<double, 1, Eigen::Dynamic>;

using Matrix = MatrixNd;
using Point = VectorNd;
using Vector = VectorNd;

static const int s_max_free_dofs = NDIM * (NDIM + 1) / 2;
using RigidDOFVector = Eigen::Matrix<double, s_max_free_dofs, 1>;
using FreeRigidDOFVector = Eigen::Matrix<int, s_max_free_dofs, 1>;
using RDV = RigidDOFVector;
using FRDV = FreeRigidDOFVector;

static const int invalid_level_number = -1;
static const int invalid_index = -1;

/*!
 * Deallocate a SAMRAIVectorReal.
 */
inline void
deallocate_vector_data(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                       int coarsest_ln = invalid_level_number,
                       int finest_ln = invalid_level_number)
{
    const int coarsest_ln_in = x.getCoarsestLevelNumber();
    const int finest_ln_in = x.getFinestLevelNumber();

    // By default, deallocate data over the full range of levels set in the vector.
    if (coarsest_ln == invalid_level_number) coarsest_ln = coarsest_ln_in;
    if (finest_ln == invalid_level_number) finest_ln = finest_ln_in;

    // If the coarsest level in the vector is finer than the finest level in the patch hierarchy, then there (should
    // be!) no data to deallocate.
    const auto& hierarchy = x.getPatchHierarchy();
    if (coarsest_ln > hierarchy->getFinestLevelNumber()) return;

    // Ensure that we only deallocate vector data on a range of levels that is actually in the hierarchy:
    finest_ln = std::min(finest_ln, hierarchy->getFinestLevelNumber());
    x.resetLevels(coarsest_ln, finest_ln);
    x.deallocateVectorData();

    // Restore the range of level numbers (even if that range of levels no longer makes sense for the hierarchy
    // configuration):
    x.resetLevels(coarsest_ln_in, finest_ln_in);
    return;
}

/*!
 * Free the components of a SAMRAIVectorReal.
 */
inline void
free_vector_components(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                       int coarsest_ln = invalid_level_number,
                       int finest_ln = invalid_level_number)
{
    const int coarsest_ln_in = x.getCoarsestLevelNumber();
    const int finest_ln_in = x.getFinestLevelNumber();

    // By default, deallocate data over the full range of levels set in the vector.
    if (coarsest_ln == invalid_level_number) coarsest_ln = coarsest_ln_in;
    if (finest_ln == invalid_level_number) finest_ln = finest_ln_in;

    // If the coarsest level in the vector is finer than the finest level in the patch hierarchy, then there (should
    // be!) no components to free on those levels.
    const auto& hierarchy = x.getPatchHierarchy();
    if (coarsest_ln > hierarchy->getFinestLevelNumber()) return;

    // Ensure that we only free vector components on a range of levels that is actually in the hierarchy:
    finest_ln = std::min(finest_ln, hierarchy->getFinestLevelNumber());
    x.resetLevels(coarsest_ln, finest_ln);
    x.freeVectorComponents();

    // Restore the range of level numbers (even if that range of levels no longer makes sense for the hierarchy
    // configuration):
    x.resetLevels(coarsest_ln_in, finest_ln_in);
    return;
}

/*!
 * Utility function which asserts that the SAMRAI pointer is not null.
 *
 * This is useful for writing generic code in which we might want to assert that
 * a pointer is not null in the initialization list where we cannot use the
 * normal assertion macro.
 */
template <class T>
const T&
checked_dereference(const SAMRAI::tbox::Pointer<T>& p)
{
#ifndef NDEBUG
    TBOX_ASSERT(p);
#endif
    return *p;
}

/*!
 * Same idea, but for a non-const pointer.
 */
template <class T>
T&
checked_dereference(SAMRAI::tbox::Pointer<T>& p)
{
#ifndef NDEBUG
    TBOX_ASSERT(p);
#endif
    return *p;
}

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_ibtk_utilities
