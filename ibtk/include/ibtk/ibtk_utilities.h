// Filename: ibtk_utilities.h
// Created on 27 Jan 2011 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#ifndef included_IBTK_ibtk_utilities
#define included_IBTK_ibtk_utilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <algorithm>
#include <array>

#include "Eigen/Core" // IWYU pragma: export
#include "Eigen/StdVector"
#include "tbox/PIO.h"
#include "tbox/Utilities.h"

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
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_ibtk_utilities
