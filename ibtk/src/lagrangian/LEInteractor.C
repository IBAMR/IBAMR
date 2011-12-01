// Filename: LEInteractor.C
// Created on 14 Jul 2004 by Boyce Griffith
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

#include "LEInteractor.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/LNodeIndex.h>
#include <ibtk/LNodeIndexSet.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/compiler_hints.h>
#include <ibtk/namespaces.h>

// BLITZ++ INCLUDES
#include <blitz/array.h>

// SAMRAI INCLUDES
#include <CellIndex.h>
#include <CartesianPatchGeometry.h>
#include <Index.h>
#include <IntVector.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <functional>
#include <vector>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define LAGRANGIAN_PIECEWISE_CONSTANT_INTERP_FC FC_FUNC_(lagrangian_piecewise_constant_interp2d, LAGRANGIAN_PIECEWISE_CONSTANT_INTERP2D)
#define LAGRANGIAN_PIECEWISE_CONSTANT_SPREAD_FC FC_FUNC_(lagrangian_piecewise_constant_spread2d, LAGRANGIAN_PIECEWISE_CONSTANT_SPREAD2D)

#define LAGRANGIAN_PIECEWISE_LINEAR_INTERP_FC FC_FUNC_(lagrangian_piecewise_linear_interp2d, LAGRANGIAN_PIECEWISE_LINEAR_INTERP2D)
#define LAGRANGIAN_PIECEWISE_LINEAR_SPREAD_FC FC_FUNC_(lagrangian_piecewise_linear_spread2d, LAGRANGIAN_PIECEWISE_LINEAR_SPREAD2D)

#define LAGRANGIAN_WIDE_PIECEWISE_LINEAR_INTERP_FC FC_FUNC_(lagrangian_wide_piecewise_linear_interp2d, LAGRANGIAN_WIDE_PIECEWISE_LINEAR_INTERP2D)
#define LAGRANGIAN_WIDE_PIECEWISE_LINEAR_SPREAD_FC FC_FUNC_(lagrangian_wide_piecewise_linear_spread2d, LAGRANGIAN_WIDE_PIECEWISE_LINEAR_SPREAD2D)

#define LAGRANGIAN_PIECEWISE_CUBIC_INTERP_FC FC_FUNC_(lagrangian_piecewise_cubic_interp2d, LAGRANGIAN_PIECEWISE_CUBIC_INTERP2D)
#define LAGRANGIAN_PIECEWISE_CUBIC_SPREAD_FC FC_FUNC_(lagrangian_piecewise_cubic_spread2d, LAGRANGIAN_PIECEWISE_CUBIC_SPREAD2D)

#define LAGRANGIAN_WIDE_PIECEWISE_CUBIC_INTERP_FC FC_FUNC_(lagrangian_wide_piecewise_cubic_interp2d, LAGRANGIAN_WIDE_PIECEWISE_CUBIC_INTERP2D)
#define LAGRANGIAN_WIDE_PIECEWISE_CUBIC_SPREAD_FC FC_FUNC_(lagrangian_wide_piecewise_cubic_spread2d, LAGRANGIAN_WIDE_PIECEWISE_CUBIC_SPREAD2D)

#define LAGRANGIAN_IB_3_INTERP_FC FC_FUNC_(lagrangian_ib_3_interp2d, LAGRANGIAN_IB_3_INTERP2D)
#define LAGRANGIAN_IB_3_SPREAD_FC FC_FUNC_(lagrangian_ib_3_spread2d, LAGRANGIAN_IB_3_SPREAD2D)

#define LAGRANGIAN_WIDE_IB_3_INTERP_FC FC_FUNC_(lagrangian_wide_ib_3_interp2d, LAGRANGIAN_WIDE_IB_3_INTERP2D)
#define LAGRANGIAN_WIDE_IB_3_SPREAD_FC FC_FUNC_(lagrangian_wide_ib_3_spread2d, LAGRANGIAN_WIDE_IB_3_SPREAD2D)

#define LAGRANGIAN_IB_4_INTERP_FC FC_FUNC_(lagrangian_ib_4_interp2d, LAGRANGIAN_IB_4_INTERP2D)
#define LAGRANGIAN_IB_4_SPREAD_FC FC_FUNC_(lagrangian_ib_4_spread2d, LAGRANGIAN_IB_4_SPREAD2D)
#define LAGRANGIAN_IB_4_SPREAD_XP_FC FC_FUNC_(lagrangian_ib_4_spread_xp2d, LAGRANGIAN_IB_4_SPREAD_XP2D)

#define LAGRANGIAN_WIDE_IB_4_INTERP_FC FC_FUNC_(lagrangian_wide_ib_4_interp2d, LAGRANGIAN_WIDE_IB_4_INTERP2D)
#define LAGRANGIAN_WIDE_IB_4_SPREAD_FC FC_FUNC_(lagrangian_wide_ib_4_spread2d, LAGRANGIAN_WIDE_IB_4_SPREAD2D)

#define LAGRANGIAN_IB_6_INTERP_FC FC_FUNC_(lagrangian_ib_6_interp2d, LAGRANGIAN_IB_6_INTERP2D)
#define LAGRANGIAN_IB_6_SPREAD_FC FC_FUNC_(lagrangian_ib_6_spread2d, LAGRANGIAN_IB_6_SPREAD2D)
#define LAGRANGIAN_IB_6_SPREAD_XP_FC FC_FUNC_(lagrangian_ib_6_spread_xp2d, LAGRANGIAN_IB_6_SPREAD_XP2D)
#endif

#if (NDIM == 3)
#define LAGRANGIAN_PIECEWISE_CONSTANT_INTERP_FC FC_FUNC_(lagrangian_piecewise_constant_interp3d, LAGRANGIAN_PIECEWISE_CONSTANT_INTERP3D)
#define LAGRANGIAN_PIECEWISE_CONSTANT_SPREAD_FC FC_FUNC_(lagrangian_piecewise_constant_spread3d, LAGRANGIAN_PIECEWISE_CONSTANT_SPREAD3D)

#define LAGRANGIAN_IB_3_INTERP_FC FC_FUNC_(lagrangian_ib_3_interp3d, LAGRANGIAN_IB_3_INTERP3D)
#define LAGRANGIAN_IB_3_SPREAD_FC FC_FUNC_(lagrangian_ib_3_spread3d, LAGRANGIAN_IB_3_SPREAD3D)

#define LAGRANGIAN_IB_4_INTERP_FC FC_FUNC_(lagrangian_ib_4_interp3d, LAGRANGIAN_IB_4_INTERP3D)
#define LAGRANGIAN_IB_4_SPREAD_FC FC_FUNC_(lagrangian_ib_4_spread3d, LAGRANGIAN_IB_4_SPREAD3D)
#define LAGRANGIAN_IB_4_SPREAD_XP_FC FC_FUNC_(lagrangian_ib_4_spread_xp3d, LAGRANGIAN_IB_4_SPREAD_XP3D)

#define LAGRANGIAN_IB_6_INTERP_FC FC_FUNC_(lagrangian_ib_6_interp3d, LAGRANGIAN_IB_6_INTERP3D)
#define LAGRANGIAN_IB_6_SPREAD_FC FC_FUNC_(lagrangian_ib_6_spread3d, LAGRANGIAN_IB_6_SPREAD3D)
#define LAGRANGIAN_IB_6_SPREAD_XP_FC FC_FUNC_(lagrangian_ib_6_spread_xp3d, LAGRANGIAN_IB_6_SPREAD_XP3D)
#endif

extern "C"
{
    void
    LAGRANGIAN_PIECEWISE_CONSTANT_INTERP_FC(
        const double* , const double* , const double* , const int& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        const double* ,
        const int* , const double* , const int& ,
        const double* , double*
                                            );

    void
    LAGRANGIAN_PIECEWISE_CONSTANT_SPREAD_FC(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        double*
                                            );

#if (NDIM == 2)
    void
    LAGRANGIAN_PIECEWISE_LINEAR_INTERP_FC(
        const double* , const double* , const double* , const int& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        const double* ,
        const int* , const double* , const int& ,
        const double* , double*
                                          );

    void
    LAGRANGIAN_PIECEWISE_LINEAR_SPREAD_FC(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        double*
                                          );

    void
    LAGRANGIAN_WIDE_PIECEWISE_LINEAR_INTERP_FC(
        const double* , const double* , const double* , const int& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        const double* ,
        const int* , const double* , const int& ,
        const double* , double*
                                               );

    void
    LAGRANGIAN_WIDE_PIECEWISE_LINEAR_SPREAD_FC(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        double*
                                               );
    void
    LAGRANGIAN_PIECEWISE_CUBIC_INTERP_FC(
        const double* , const double* , const double* , const int& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        const double* ,
        const int* , const double* , const int& ,
        const double* , double*
                                         );

    void
    LAGRANGIAN_PIECEWISE_CUBIC_SPREAD_FC(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        double*
                                         );

    void
    LAGRANGIAN_WIDE_PIECEWISE_CUBIC_INTERP_FC(
        const double* , const double* , const double* , const int& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        const double* ,
        const int* , const double* , const int& ,
        const double* , double*
                                              );

    void
    LAGRANGIAN_WIDE_PIECEWISE_CUBIC_SPREAD_FC(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        double*
                                              );
#endif
    void
    LAGRANGIAN_IB_3_INTERP_FC(
        const double* , const double* , const double* , const int& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        const double* ,
        const int* , const double* , const int& ,
        const double* , double*
                              );

    void
    LAGRANGIAN_IB_3_SPREAD_FC(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        double*
                              );
#if (NDIM == 2)
    void
    LAGRANGIAN_WIDE_IB_3_INTERP_FC(
        const double* , const double* , const double* , const int& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        const double* ,
        const int* , const double* , const int& ,
        const double* , double*
                                   );

    void
    LAGRANGIAN_WIDE_IB_3_SPREAD_FC(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        double*
                                   );
#endif
    void
    LAGRANGIAN_IB_4_INTERP_FC(
        const double* , const double* , const double* , const int& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int* , const int* ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int* , const int* ,
        const int& , const int& , const int& ,
#endif
        const double* ,
        const int* , const double* , const int& ,
        const double* , double*
                              );

    void
    LAGRANGIAN_IB_4_SPREAD_FC(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int* , const int* ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int* , const int* ,
        const int& , const int& , const int& ,
#endif
        double*
                              );

    void
    LAGRANGIAN_IB_4_SPREAD_XP_FC(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int* , const int* ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int* , const int* ,
        const int& , const int& , const int& ,
#endif
        double*
                                 );
#if (NDIM == 2)
    void
    LAGRANGIAN_WIDE_IB_4_INTERP_FC(
        const double* , const double* , const double* , const int& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        const double* ,
        const int* , const double* , const int& ,
        const double* , double*
                                   );

    void
    LAGRANGIAN_WIDE_IB_4_SPREAD_FC(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        double*
                                   );
#endif
    void
    LAGRANGIAN_IB_6_INTERP_FC(
        const double* , const double* , const double* , const int& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int* , const int* ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int* , const int* ,
        const int& , const int& , const int& ,
#endif
        const double* ,
        const int* , const double* , const int& ,
        const double* , double*
                              );

    void
    LAGRANGIAN_IB_6_SPREAD_FC(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int* , const int* ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int* , const int* ,
        const int& , const int& , const int& ,
#endif
        double*
                              );

    void
    LAGRANGIAN_IB_6_SPREAD_XP_FC(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int* , const int* ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int* , const int* ,
        const int& , const int& , const int& ,
#endif
        double*
                                 );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline double
ib4_delta_fcn(
    double r)
{
    r = std::abs(r);
    if (r < 1.0)
    {
        const double t2 = r * r;
        const double t6 = sqrt(-0.4e1 * t2 + 0.4e1 * r + 0.1e1);
        return -r / 0.4e1 + 0.3e1 / 0.8e1 + t6 / 0.8e1;
    }
    else if (r < 2.0)
    {
        const double t2 = r * r;
        const double t6 = sqrt(0.12e2 * r - 0.7e1 - 0.4e1 * t2);
        return -r / 0.4e1 + 0.5e1 / 0.8e1 - t6 / 0.8e1;
    }
    else
    {
        return 0.0;
    }
}// ib4_delta_fcn

struct SortModeComp
    : std::binary_function<std::pair<const LNodeIndex*,blitz::TinyVector<double,NDIM> >,std::pair<const LNodeIndex*,blitz::TinyVector<double,NDIM> >,bool>
{
    inline bool
    operator()(
        const std::pair<const LNodeIndex*,blitz::TinyVector<double,NDIM> >& lhs,
        const std::pair<const LNodeIndex*,blitz::TinyVector<double,NDIM> >& rhs) const
        {
            if (LEInteractor::s_sort_mode == LEInteractor::SORT_INCREASING_LAG_IDX) return lhs.first->getLagrangianIndex() < rhs.first->getLagrangianIndex();
            if (LEInteractor::s_sort_mode == LEInteractor::SORT_DECREASING_LAG_IDX) return lhs.first->getLagrangianIndex() > rhs.first->getLagrangianIndex();
            return true;
        }
};
}

double (*LEInteractor::s_delta_fcn)(double r) = &ib4_delta_fcn;
int LEInteractor::s_delta_fcn_stencil_size = 4;
double LEInteractor::s_delta_fcn_C = 3.0/8.0;
LEInteractor::SortMode LEInteractor::s_sort_mode = NO_SORT;
LEInteractor::PrecisionMode LEInteractor::s_precision_mode  = DOUBLE;

void
LEInteractor::setFromDatabase(
    Pointer<Database> db)
{
    if (db.isNull()) return;
    const std::string debug_sort_mode_str = db->getStringWithDefault("debug_sort_mode", "NO_SORT");
    const std::string precision_mode_str = db->getStringWithDefault("precision_mode", "DOUBLE");

    if (debug_sort_mode_str == "NO_SORT")
    {
        s_sort_mode = NO_SORT;
    }
    else if (debug_sort_mode_str == "SORT_INCREASING_LAG_IDX")
    {
        s_sort_mode = SORT_INCREASING_LAG_IDX;
    }
    else if (debug_sort_mode_str == "SORT_DECREASING_LAG_IDX")
    {
        s_sort_mode = SORT_DECREASING_LAG_IDX;
    }
    else
    {
        TBOX_ERROR("LEInteractor::setFromDatabase():\n"
                   << ":  invalid debug_sort_mode: " << debug_sort_mode_str << ".\n"
                   << "   Choices are: NO_SORT, SORT_INCREASING_LAG_IDX, SORT_DECREASING_LAG_IDX.\n");
    }

    if (precision_mode_str == "DOUBLE")
    {
        s_precision_mode = DOUBLE;
    }
    else if (precision_mode_str == "DOUBLE_DOUBLE")
    {
        s_precision_mode = DOUBLE_DOUBLE;
    }
    else
    {
        TBOX_ERROR("LEInteractor::setFromDatabase():\n"
                   << ":  invalid precision_mode: " << precision_mode_str << ".\n"
                   << "   Choices are: DOUBLE, DOUBLE_DOUBLE.\n");
    }
    return;
}// setFromDatabase

void
LEInteractor::printClassData(
    std::ostream& os)
{
    os << "LEInteractor::printClassData():\n";
    os << "  s_sort_mode = ";
    if (s_sort_mode == NO_SORT)
    {
        os << "NO_SORT";
    }
    else if (s_sort_mode == SORT_INCREASING_LAG_IDX)
    {
        os << "SORT_INCREASING_LAG_IDX";
    }
    else if (s_sort_mode == SORT_DECREASING_LAG_IDX)
    {
        os << "SORT_DECREASING_LAG_IDX";
    }
    else
    {
        os << "UNKNOWN";
    }
    os << "\n";
    os << "  s_precision_mode = ";
    if (s_precision_mode == DOUBLE)
    {
        os << "DOUBLE";
    }
    else if (s_precision_mode == DOUBLE_DOUBLE)
    {
        os << "DOUBLE_DOUBLE";
    }
    else
    {
        os << "UNKNOWN";
    }
    os << "\n";
    return;
}// printClassData

/////////////////////////////// PUBLIC ///////////////////////////////////////

int
LEInteractor::getStencilSize(
    const std::string& weighting_fcn)
{
    if (weighting_fcn == "PIECEWISE_CONSTANT") return 1;
#if (NDIM == 2)
    if (weighting_fcn == "PIECEWISE_LINEAR") return 2;
    if (weighting_fcn == "WIDE_PIECEWISE_LINEAR") return 4;
    if (weighting_fcn == "PIECEWISE_CUBIC") return 4;
    if (weighting_fcn == "WIDE_PIECEWISE_CUBIC") return 8;
#endif
    if (weighting_fcn == "IB_3") return 4;
#if (NDIM == 2)
    if (weighting_fcn == "WIDE_IB_3") return 6;
#endif
    if (weighting_fcn == "IB_4") return 4;
#if (NDIM == 2)
    if (weighting_fcn == "WIDE_IB_4") return 8;
#endif
    if (weighting_fcn == "IB_6") return 6;
    if (weighting_fcn == "USER_DEFINED") return s_delta_fcn_stencil_size;
    TBOX_ERROR("LEInteractor::getStencilSize()\n"
               << "  Unknown weighting function "
               << weighting_fcn << std::endl);
    return -1;
}// getStencilSize

double
LEInteractor::getC(
    const std::string& weighting_fcn)
{
    if (weighting_fcn == "PIECEWISE_CONSTANT") return 1.0;
#if (NDIM == 2)
    if (weighting_fcn == "PIECEWISE_LINEAR" ||
        weighting_fcn == "WIDE_PIECEWISE_LINEAR" ||
        weighting_fcn == "PIECEWISE_CUBIC" ||
        weighting_fcn == "WIDE_PIECEWISE_CUBIC")
    {
        TBOX_ERROR("LEInteractor::getC()\n"
                   << "  Weighting function " << weighting_fcn << " does not satisfy the quadratic (sum-of-squares) condition.\n"
                   << "  Consequently, the value of C = sum_{j} phi(r-j)^2 is dependent on the value of r." << std::endl);
    }
#endif
    if (weighting_fcn == "IB_3") return 0.5;
#if (NDIM == 2)
    if (weighting_fcn == "WIDE_IB_3") return 0.25*getC("IB_3");
#endif
    if (weighting_fcn == "IB_4") return 3.0/8.0;
#if (NDIM == 2)
    if (weighting_fcn == "WIDE_IB_4") return 0.25*getC("IB_4");
#endif
    if (weighting_fcn == "IB_6") return 67.0/128.0;
    if (weighting_fcn == "USER_DEFINED") return s_delta_fcn_C;
    TBOX_ERROR("LEInteractor::getC()\n"
               << "  Unknown weighting function "
               << weighting_fcn << std::endl);
    return -1;
}// getC

template<class T>
void
LEInteractor::interpolate(
    Pointer<LData> Q_data,
    const Pointer<LData> X_data,
    const Pointer<LIndexSetData<T> > idx_data,
    const Pointer<CellData<NDIM,double> > q_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& interp_box,
    const IntVector<NDIM>& periodic_shift,
    const std::string& interp_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!Q_data.isNull());
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!X_data.isNull());
    TBOX_ASSERT(!idx_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_data->getDepth() == static_cast<unsigned int>(q_data->getDepth()));
    TBOX_ASSERT(X_data->getDepth() == NDIM);
#endif
    interpolate(Q_data->getGhostedLocalFormVecArray()->data(), Q_data->getDepth(),
                X_data->getGhostedLocalFormVecArray()->data(), X_data->getDepth(),
                idx_data,
                q_data,
                patch, interp_box, periodic_shift, interp_fcn);
    Q_data->restoreArrays();
    X_data->restoreArrays();
    return;
}// interpolate

template<class T>
void
LEInteractor::interpolate(
    Pointer<LData> Q_data,
    const Pointer<LData> X_data,
    const Pointer<LIndexSetData<T> > idx_data,
    const Pointer<NodeData<NDIM,double> > q_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& interp_box,
    const IntVector<NDIM>& periodic_shift,
    const std::string& interp_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!Q_data.isNull());
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!X_data.isNull());
    TBOX_ASSERT(!idx_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_data->getDepth() == static_cast<unsigned int>(q_data->getDepth()));
    TBOX_ASSERT(X_data->getDepth() == NDIM);
#endif
    interpolate(Q_data->getGhostedLocalFormVecArray()->data(), Q_data->getDepth(),
                X_data->getGhostedLocalFormVecArray()->data(), X_data->getDepth(),
                idx_data,
                q_data,
                patch, interp_box, periodic_shift, interp_fcn);
    Q_data->restoreArrays();
    X_data->restoreArrays();
    return;
}// interpolate

template<class T>
void
LEInteractor::interpolate(
    Pointer<LData> Q_data,
    const Pointer<LData> X_data,
    const Pointer<LIndexSetData<T> > idx_data,
    const Pointer<SideData<NDIM,double> > q_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& interp_box,
    const IntVector<NDIM>& periodic_shift,
    const std::string& interp_fcn)
{
    if (Q_data->getDepth() != NDIM || q_data->getDepth() != 1)
    {
        TBOX_ERROR("LEInteractor::interpolate():\n"
                   << "  side-centered interpolation requires vector-valued data.\n");
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!Q_data.isNull());
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!X_data.isNull());
    TBOX_ASSERT(!idx_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_data->getDepth() == NDIM);
    TBOX_ASSERT(X_data->getDepth() == NDIM);
    TBOX_ASSERT(q_data->getDepth() == 1);
#endif
    interpolate(Q_data->getGhostedLocalFormVecArray()->data(), Q_data->getDepth(),
                X_data->getGhostedLocalFormVecArray()->data(), X_data->getDepth(),
                idx_data,
                q_data,
                patch, interp_box, periodic_shift, interp_fcn);
    Q_data->restoreArrays();
    X_data->restoreArrays();
    return;
}// interpolate

template<class T>
void
LEInteractor::interpolate(
    double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const Pointer<LIndexSetData<T> > idx_data,
    const Pointer<CellData<NDIM,double> > q_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& interp_box,
    const IntVector<NDIM>& periodic_shift,
    const std::string& interp_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!idx_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_depth == q_data->getDepth());
    TBOX_ASSERT(X_depth == NDIM);
#else
    NULL_USE(X_depth);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    blitz::TinyVector<int,NDIM> patch_touches_lower_physical_bdry(0);
    blitz::TinyVector<int,NDIM> patch_touches_upper_physical_bdry(0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,upper);
    }

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    std::vector<double> periodic_offsets;
    buildLocalIndices(local_indices, periodic_offsets, interp_box, patch, periodic_shift, idx_data);

    // Interpolate.
    if (!local_indices.empty())
    {
        interpolate(Q_data, Q_depth, X_data,
                    q_data->getPointer(), q_data->getBox(), q_data->getGhostCellWidth(), q_data->getDepth(),
                    x_lower, x_upper, dx,
                    patch_touches_lower_physical_bdry, patch_touches_upper_physical_bdry,
                    local_indices, periodic_offsets,
                    interp_fcn);
    }
    return;
}// interpolate

template<class T>
void
LEInteractor::interpolate(
    double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const Pointer<LIndexSetData<T> > idx_data,
    const Pointer<NodeData<NDIM,double> > q_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& interp_box,
    const IntVector<NDIM>& periodic_shift,
    const std::string& interp_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!idx_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_depth == q_data->getDepth());
    TBOX_ASSERT(X_depth == NDIM);
#else
    NULL_USE(X_depth);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    blitz::TinyVector<int,NDIM> patch_touches_lower_physical_bdry(0);
    blitz::TinyVector<int,NDIM> patch_touches_upper_physical_bdry(0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,upper);
    }

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    std::vector<double> periodic_offsets;
    buildLocalIndices(local_indices, periodic_offsets, interp_box, patch, periodic_shift, idx_data);

    // Interpolate.
    if (!local_indices.empty())
    {
        blitz::TinyVector<double,NDIM> x_lower_node, x_upper_node;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            x_lower_node[d] = x_lower[d]-0.5*dx[d];
            x_upper_node[d] = x_upper[d]+0.5*dx[d];
        }
        interpolate(Q_data, Q_depth, X_data,
                    q_data->getPointer(), NodeGeometry<NDIM>::toNodeBox(q_data->getBox()), q_data->getGhostCellWidth(), q_data->getDepth(),
                    x_lower_node.data(), x_upper_node.data(), dx,
                    patch_touches_lower_physical_bdry, patch_touches_upper_physical_bdry,
                    local_indices, periodic_offsets,
                    interp_fcn);
    }
    return;
}// interpolate

template<class T>
void
LEInteractor::interpolate(
    double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const Pointer<LIndexSetData<T> > idx_data,
    const Pointer<SideData<NDIM,double> > q_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& interp_box,
    const IntVector<NDIM>& periodic_shift,
    const std::string& interp_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!idx_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_depth == NDIM);
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(q_data->getDepth() == 1);
#else
    NULL_USE(X_depth);
#endif
    if (Q_depth != NDIM || q_data->getDepth() != 1)
    {
        TBOX_ERROR("LEInteractor::interpolate():\n"
                   << "  side-centered interpolation requires vector-valued data.\n");
    }

    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    blitz::TinyVector<int,NDIM> patch_touches_lower_physical_bdry(0);
    blitz::TinyVector<int,NDIM> patch_touches_upper_physical_bdry(0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,upper);
    }

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    std::vector<double> periodic_offsets;
    buildLocalIndices(local_indices, periodic_offsets, interp_box, patch, periodic_shift, idx_data);

    // Interpolate.
    if (!local_indices.empty())
    {
        blitz::TinyVector<double,NDIM> x_lower_axis, x_upper_axis;
        const int local_sz = (*std::max_element(local_indices.begin(),local_indices.end()))+1;
        std::vector<double> Q_data_axis(local_sz);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                x_lower_axis[d] = x_lower[d];
                x_upper_axis[d] = x_upper[d];
            }
            x_lower_axis[axis] -= 0.5*dx[axis];
            x_upper_axis[axis] += 0.5*dx[axis];
            interpolate(&Q_data_axis[0], 1, X_data,
                        q_data->getPointer(axis), SideGeometry<NDIM>::toSideBox(q_data->getBox(), axis), q_data->getGhostCellWidth(), 1,
                        x_lower_axis.data(), x_upper_axis.data(), dx,
                        patch_touches_lower_physical_bdry, patch_touches_upper_physical_bdry,
                        local_indices, periodic_offsets,
                        interp_fcn);
            for (unsigned int k = 0; k < local_indices.size(); ++k)
            {
                Q_data[NDIM*local_indices[k]+axis] = Q_data_axis[local_indices[k]];
            }
        }
    }
    return;
}// interpolate

void
LEInteractor::interpolate(
    std::vector<double>& Q_data,
    const int Q_depth,
    const std::vector<double>& X_data,
    const int X_depth,
    const Pointer<CellData<NDIM,double> > q_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& interp_box,
    const std::string& interp_fcn)
{
    if (Q_data.size() == 0) return;
    interpolate(&Q_data[0], Q_data.size(), Q_depth,
                &X_data[0], X_data.size(), X_depth,
                q_data, patch, interp_box, interp_fcn);
}// interpolate

void
LEInteractor::interpolate(
    std::vector<double>& Q_data,
    const int Q_depth,
    const std::vector<double>& X_data,
    const int X_depth,
    const Pointer<NodeData<NDIM,double> > q_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& interp_box,
    const std::string& interp_fcn)
{
    if (Q_data.size() == 0) return;
    interpolate(&Q_data[0], Q_data.size(), Q_depth,
                &X_data[0], X_data.size(), X_depth,
                q_data, patch, interp_box, interp_fcn);
}// interpolate

void
LEInteractor::interpolate(
    std::vector<double>& Q_data,
    const int Q_depth,
    const std::vector<double>& X_data,
    const int X_depth,
    const Pointer<SideData<NDIM,double> > q_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& interp_box,
    const std::string& interp_fcn)
{
    if (Q_data.size() == 0) return;
    interpolate(&Q_data[0], Q_data.size(), Q_depth,
                &X_data[0], X_data.size(), X_depth,
                q_data, patch, interp_box, interp_fcn);
}// interpolate

void
LEInteractor::interpolate(
    double* const Q_data,
    const int Q_size,
    const int Q_depth,
    const double* const X_data,
    const int X_size,
    const int X_depth,
    const Pointer<CellData<NDIM,double> > q_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& interp_box,
    const std::string& interp_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_depth == q_data->getDepth());
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(Q_size/Q_depth == X_size/X_depth);
#else
    NULL_USE(Q_size);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    blitz::TinyVector<int,NDIM> patch_touches_lower_physical_bdry(0);
    blitz::TinyVector<int,NDIM> patch_touches_upper_physical_bdry(0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, interp_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_offsets(NDIM*local_indices.size());

    // Interpolate.
    if (!local_indices.empty())
    {
        interpolate(Q_data, Q_depth, X_data,
                    q_data->getPointer(), q_data->getBox(), q_data->getGhostCellWidth(), q_data->getDepth(),
                    x_lower, x_upper, dx,
                    patch_touches_lower_physical_bdry, patch_touches_upper_physical_bdry,
                    local_indices, periodic_offsets,
                    interp_fcn);
    }
    return;
}// interpolate

void
LEInteractor::interpolate(
    double* const Q_data,
    const int Q_size,
    const int Q_depth,
    const double* const X_data,
    const int X_size,
    const int X_depth,
    const Pointer<NodeData<NDIM,double> > q_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& interp_box,
    const std::string& interp_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_depth == q_data->getDepth());
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(Q_size/Q_depth == X_size/X_depth);
#else
    NULL_USE(Q_size);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    blitz::TinyVector<int,NDIM> patch_touches_lower_physical_bdry(0);
    blitz::TinyVector<int,NDIM> patch_touches_upper_physical_bdry(0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, interp_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_offsets(NDIM*local_indices.size());

    // Interpolate.
    if (!local_indices.empty())
    {
        blitz::TinyVector<double,NDIM> x_lower_node, x_upper_node;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            x_lower_node[d] = x_lower[d]-0.5*dx[d];
            x_upper_node[d] = x_upper[d]+0.5*dx[d];
        }
        interpolate(Q_data, Q_depth, X_data,
                    q_data->getPointer(), NodeGeometry<NDIM>::toNodeBox(q_data->getBox()), q_data->getGhostCellWidth(), q_data->getDepth(),
                    x_lower_node.data(), x_upper_node.data(), dx,
                    patch_touches_lower_physical_bdry, patch_touches_upper_physical_bdry,
                    local_indices, periodic_offsets,
                    interp_fcn);
    }
    return;
}// interpolate

void
LEInteractor::interpolate(
    double* const Q_data,
    const int Q_size,
    const int Q_depth,
    const double* const X_data,
    const int X_size,
    const int X_depth,
    const Pointer<SideData<NDIM,double> > q_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& interp_box,
    const std::string& interp_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_depth == NDIM);
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(Q_size/Q_depth == X_size/X_depth);
    TBOX_ASSERT(q_data->getDepth() == 1);
#else
    NULL_USE(Q_size);
#endif
    if (Q_depth != NDIM || q_data->getDepth() != 1)
    {
        TBOX_ERROR("LEInteractor::interpolate():\n"
                   << "  side-centered interpolation requires vector-valued data.\n");
    }

    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    blitz::TinyVector<int,NDIM> patch_touches_lower_physical_bdry(0);
    blitz::TinyVector<int,NDIM> patch_touches_upper_physical_bdry(0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, interp_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_offsets(NDIM*local_indices.size());

    // Interpolate.
    if (!local_indices.empty())
    {
        blitz::TinyVector<double,NDIM> x_lower_axis, x_upper_axis;
        const int local_sz = (*std::max_element(local_indices.begin(),local_indices.end()))+1;
        std::vector<double> Q_data_axis(local_sz);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                x_lower_axis[d] = x_lower[d];
                x_upper_axis[d] = x_upper[d];
            }
            x_lower_axis[axis] -= 0.5*dx[axis];
            x_upper_axis[axis] += 0.5*dx[axis];
            interpolate(&Q_data_axis[0], 1, X_data,
                        q_data->getPointer(axis), SideGeometry<NDIM>::toSideBox(q_data->getBox(), axis), q_data->getGhostCellWidth(), 1,
                        x_lower_axis.data(), x_upper_axis.data(), dx,
                        patch_touches_lower_physical_bdry, patch_touches_upper_physical_bdry,
                        local_indices, periodic_offsets,
                        interp_fcn);
            for (unsigned int k = 0; k < local_indices.size(); ++k)
            {
                Q_data[NDIM*local_indices[k]+axis] = Q_data_axis[local_indices[k]];
            }
        }
    }
    return;
}// interpolate

template<class T>
void
LEInteractor::spread(
    Pointer<CellData<NDIM,double> > q_data,
    const Pointer<LData> Q_data,
    const Pointer<LData> X_data,
    const Pointer<LIndexSetData<T> > idx_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& spread_box,
    const IntVector<NDIM>& periodic_shift,
    const std::string& spread_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!Q_data.isNull());
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!X_data.isNull());
    TBOX_ASSERT(!idx_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_data->getDepth() == static_cast<unsigned int>(q_data->getDepth()));
    TBOX_ASSERT(X_data->getDepth() == NDIM);
#endif
    spread(q_data,
           Q_data->getGhostedLocalFormVecArray()->data(), Q_data->getDepth(),
           X_data->getGhostedLocalFormVecArray()->data(), X_data->getDepth(),
           idx_data,
           patch, spread_box, periodic_shift, spread_fcn);
    Q_data->restoreArrays();
    X_data->restoreArrays();
    return;
}// spread

template<class T>
void
LEInteractor::spread(
    Pointer<NodeData<NDIM,double> > q_data,
    const Pointer<LData> Q_data,
    const Pointer<LData> X_data,
    const Pointer<LIndexSetData<T> > idx_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& spread_box,
    const IntVector<NDIM>& periodic_shift,
    const std::string& spread_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!Q_data.isNull());
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!X_data.isNull());
    TBOX_ASSERT(!idx_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_data->getDepth() == static_cast<unsigned int>(q_data->getDepth()));
    TBOX_ASSERT(X_data->getDepth() == NDIM);
#endif
    spread(q_data,
           Q_data->getGhostedLocalFormVecArray()->data(), Q_data->getDepth(),
           X_data->getGhostedLocalFormVecArray()->data(), X_data->getDepth(),
           idx_data,
           patch, spread_box, periodic_shift, spread_fcn);
    Q_data->restoreArrays();
    X_data->restoreArrays();
    return;
}// spread

template<class T>
void
LEInteractor::spread(
    Pointer<SideData<NDIM,double> > q_data,
    const Pointer<LData> Q_data,
    const Pointer<LData> X_data,
    const Pointer<LIndexSetData<T> > idx_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& spread_box,
    const IntVector<NDIM>& periodic_shift,
    const std::string& spread_fcn)
{
    if (Q_data->getDepth() != NDIM || q_data->getDepth() != 1)
    {
        TBOX_ERROR("LEInteractor::spread():\n"
                   << "  side-centered spreading requires vector-valued data.\n");
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!Q_data.isNull());
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!X_data.isNull());
    TBOX_ASSERT(!idx_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(q_data->getDepth() == 1);
    TBOX_ASSERT(Q_data->getDepth() == NDIM);
    TBOX_ASSERT(X_data->getDepth() == NDIM);
#endif
    spread(q_data,
           Q_data->getGhostedLocalFormVecArray()->data(), Q_data->getDepth(),
           X_data->getGhostedLocalFormVecArray()->data(), X_data->getDepth(),
           idx_data,
           patch, spread_box, periodic_shift, spread_fcn);
    Q_data->restoreArrays();
    X_data->restoreArrays();
    return;
}// spread

template<class T>
void
LEInteractor::spread(
    Pointer<CellData<NDIM,double> > q_data,
    const double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const Pointer<LIndexSetData<T> > idx_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& spread_box,
    const IntVector<NDIM>& periodic_shift,
    const std::string& spread_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!idx_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_depth == q_data->getDepth());
    TBOX_ASSERT(X_depth == NDIM);
#else
    NULL_USE(X_depth);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    blitz::TinyVector<int,NDIM> patch_touches_lower_physical_bdry(0);
    blitz::TinyVector<int,NDIM> patch_touches_upper_physical_bdry(0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,upper);
    }

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    std::vector<double> periodic_offsets;
    buildLocalIndices(local_indices, periodic_offsets, spread_box, patch, periodic_shift, idx_data);

    // Spread.
    if (!local_indices.empty())
    {
        spread(q_data->getPointer(), q_data->getBox(), q_data->getGhostCellWidth(), q_data->getDepth(),
               Q_data, Q_depth, X_data,
               x_lower, x_upper, dx,
               patch_touches_lower_physical_bdry, patch_touches_upper_physical_bdry,
               local_indices, periodic_offsets,
               spread_fcn);
    }
    return;
}// spread

template<class T>
void
LEInteractor::spread(
    Pointer<NodeData<NDIM,double> > q_data,
    const double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const Pointer<LIndexSetData<T> > idx_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& spread_box,
    const IntVector<NDIM>& periodic_shift,
    const std::string& spread_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!idx_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_depth == q_data->getDepth());
    TBOX_ASSERT(X_depth == NDIM);
#else
    NULL_USE(X_depth);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    blitz::TinyVector<int,NDIM> patch_touches_lower_physical_bdry(0);
    blitz::TinyVector<int,NDIM> patch_touches_upper_physical_bdry(0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,upper);
    }

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    std::vector<double> periodic_offsets;
    buildLocalIndices(local_indices, periodic_offsets, spread_box, patch, periodic_shift, idx_data);

    // Spread.
    if (!local_indices.empty())
    {
        blitz::TinyVector<double,NDIM> x_lower_node, x_upper_node;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            x_lower_node[d] = x_lower[d]-0.5*dx[d];
            x_upper_node[d] = x_upper[d]+0.5*dx[d];
        }
        spread(q_data->getPointer(), NodeGeometry<NDIM>::toNodeBox(q_data->getBox()), q_data->getGhostCellWidth(), q_data->getDepth(),
               Q_data, Q_depth, X_data,
               x_lower_node.data(), x_upper_node.data(), dx,
               patch_touches_lower_physical_bdry, patch_touches_upper_physical_bdry,
               local_indices, periodic_offsets,
               spread_fcn);
    }
    return;
}// spread

template<class T>
void
LEInteractor::spread(
    Pointer<SideData<NDIM,double> > q_data,
    const double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const Pointer<LIndexSetData<T> > idx_data,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& spread_box,
    const IntVector<NDIM>& periodic_shift,
    const std::string& spread_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!idx_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(q_data->getDepth() == 1);
    TBOX_ASSERT(Q_depth == NDIM);
    TBOX_ASSERT(X_depth == NDIM);
#else
    NULL_USE(X_depth);
#endif
    if (Q_depth != NDIM || q_data->getDepth() != 1)
    {
        TBOX_ERROR("LEInteractor::spread():\n"
                   << "  side-centered spreading requires vector-valued data.\n");
    }

    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    blitz::TinyVector<int,NDIM> patch_touches_lower_physical_bdry(0);
    blitz::TinyVector<int,NDIM> patch_touches_upper_physical_bdry(0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,upper);
    }

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    std::vector<double> periodic_offsets;
    buildLocalIndices(local_indices, periodic_offsets, spread_box, patch, periodic_shift, idx_data);

    // Spread.
    if (!local_indices.empty())
    {
        blitz::TinyVector<double,NDIM> x_lower_axis, x_upper_axis;
        const int local_sz = (*std::max_element(local_indices.begin(),local_indices.end()))+1;
        std::vector<double> Q_data_axis(local_sz);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                x_lower_axis[d] = x_lower[d];
                x_upper_axis[d] = x_upper[d];
            }
            x_lower_axis[axis] -= 0.5*dx[axis];
            x_upper_axis[axis] += 0.5*dx[axis];
            for (unsigned int k = 0; k < local_indices.size(); ++k)
            {
                Q_data_axis[local_indices[k]] = Q_data[NDIM*local_indices[k]+axis];
            }
            spread(q_data->getPointer(axis), SideGeometry<NDIM>::toSideBox(q_data->getBox(),axis), q_data->getGhostCellWidth(), 1,
                   &Q_data_axis[0], 1, X_data,
                   x_lower_axis.data(), x_upper_axis.data(), dx,
                   patch_touches_lower_physical_bdry, patch_touches_upper_physical_bdry,
                   local_indices, periodic_offsets,
                   spread_fcn);
        }
    }
    return;
}// spread

void
LEInteractor::spread(
    Pointer<CellData<NDIM,double> > q_data,
    const std::vector<double>& Q_data,
    const int Q_depth,
    const std::vector<double>& X_data,
    const int X_depth,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& spread_box,
    const std::string& interp_fcn)
{
    if (Q_data.size() == 0) return;
    spread(q_data,
           &Q_data[0], Q_data.size(), Q_depth,
           &X_data[0], X_data.size(), X_depth,
           patch, spread_box, interp_fcn);
}// spread

void
LEInteractor::spread(
    Pointer<NodeData<NDIM,double> > q_data,
    const std::vector<double>& Q_data,
    const int Q_depth,
    const std::vector<double>& X_data,
    const int X_depth,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& spread_box,
    const std::string& interp_fcn)
{
    if (Q_data.size() == 0) return;
    spread(q_data,
           &Q_data[0], Q_data.size(), Q_depth,
           &X_data[0], X_data.size(), X_depth,
           patch, spread_box, interp_fcn);
}// spread

void
LEInteractor::spread(
    Pointer<SideData<NDIM,double> > q_data,
    const std::vector<double>& Q_data,
    const int Q_depth,
    const std::vector<double>& X_data,
    const int X_depth,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& spread_box,
    const std::string& interp_fcn)
{
    if (Q_data.size() == 0) return;
    spread(q_data,
           &Q_data[0], Q_data.size(), Q_depth,
           &X_data[0], X_data.size(), X_depth,
           patch, spread_box, interp_fcn);
}// spread

void
LEInteractor::spread(
    Pointer<CellData<NDIM,double> > q_data,
    const double* const Q_data,
    const int Q_size,
    const int Q_depth,
    const double* const X_data,
    const int X_size,
    const int X_depth,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& spread_box,
    const std::string& spread_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_depth == q_data->getDepth());
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(Q_size/Q_depth == X_size/X_depth);
#else
    NULL_USE(Q_size);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    blitz::TinyVector<int,NDIM> patch_touches_lower_physical_bdry(0);
    blitz::TinyVector<int,NDIM> patch_touches_upper_physical_bdry(0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, spread_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_offsets(NDIM*local_indices.size());

    // Spread.
    if (!local_indices.empty())
    {
        spread(q_data->getPointer(), q_data->getBox(), q_data->getGhostCellWidth(), q_data->getDepth(),
               Q_data, Q_depth, X_data,
               x_lower, x_upper, dx,
               patch_touches_lower_physical_bdry, patch_touches_upper_physical_bdry,
               local_indices, periodic_offsets,
               spread_fcn);
    }
    return;
}// spread

void
LEInteractor::spread(
    Pointer<NodeData<NDIM,double> > q_data,
    const double* const Q_data,
    const int Q_size,
    const int Q_depth,
    const double* const X_data,
    const int X_size,
    const int X_depth,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& spread_box,
    const std::string& spread_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_depth == q_data->getDepth());
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(Q_size/Q_depth == X_size/X_depth);
#else
    NULL_USE(Q_size);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    blitz::TinyVector<int,NDIM> patch_touches_lower_physical_bdry(0);
    blitz::TinyVector<int,NDIM> patch_touches_upper_physical_bdry(0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, spread_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_offsets(NDIM*local_indices.size());

    // Spread.
    if (!local_indices.empty())
    {
        blitz::TinyVector<double,NDIM> x_lower_node, x_upper_node;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            x_lower_node[d] = x_lower[d]-0.5*dx[d];
            x_upper_node[d] = x_upper[d]+0.5*dx[d];
        }
        spread(q_data->getPointer(), NodeGeometry<NDIM>::toNodeBox(q_data->getBox()), q_data->getGhostCellWidth(), q_data->getDepth(),
               Q_data, Q_depth, X_data,
               x_lower_node.data(), x_upper_node.data(), dx,
               patch_touches_lower_physical_bdry, patch_touches_upper_physical_bdry,
               local_indices, periodic_offsets,
               spread_fcn);
    }
    return;
}// spread

void
LEInteractor::spread(
    Pointer<SideData<NDIM,double> > q_data,
    const double* const Q_data,
    const int /*Q_size*/,
    const int Q_depth,
    const double* const X_data,
    const int X_size,
    const int X_depth,
    const Pointer<Patch<NDIM> > patch,
    const Box<NDIM>& spread_box,
    const std::string& spread_fcn)
{
    if (Q_depth != NDIM || q_data->getDepth() != 1)
    {
        TBOX_ERROR("LEInteractor::spread():\n"
                   << "  side-centered spreading requires vector-valued data.\n");
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!q_data.isNull());
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(Q_depth == NDIM);
    TBOX_ASSERT(X_depth == NDIM);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    blitz::TinyVector<int,NDIM> patch_touches_lower_physical_bdry(0);
    blitz::TinyVector<int,NDIM> patch_touches_upper_physical_bdry(0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis,upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, spread_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_offsets(NDIM*local_indices.size());

    // Spread.
    if (!local_indices.empty())
    {
        blitz::TinyVector<double,NDIM> x_lower_axis, x_upper_axis;
        const int local_sz = (*std::max_element(local_indices.begin(),local_indices.end()))+1;
        std::vector<double> Q_data_axis(local_sz);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                x_lower_axis[d] = x_lower[d];
                x_upper_axis[d] = x_upper[d];
            }
            x_lower_axis[axis] -= 0.5*dx[axis];
            x_upper_axis[axis] += 0.5*dx[axis];
            for (unsigned int k = 0; k < local_indices.size(); ++k)
            {
                Q_data_axis[local_indices[k]] = Q_data[NDIM*local_indices[k]+axis];
            }
            spread(q_data->getPointer(axis), SideGeometry<NDIM>::toSideBox(q_data->getBox(),axis), q_data->getGhostCellWidth(), 1,
                   &Q_data_axis[0], 1, X_data,
                   x_lower_axis.data(), x_upper_axis.data(), dx,
                   patch_touches_lower_physical_bdry, patch_touches_upper_physical_bdry,
                   local_indices, periodic_offsets,
                   spread_fcn);
        }
    }
    return;
}// spread

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
LEInteractor::interpolate(
    double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const double* const q_data,
    const Box<NDIM>& q_data_box,
    const IntVector<NDIM>& q_gcw,
    const int q_depth,
    const double* const x_lower,
    const double* const x_upper,
    const double* const dx,
    const blitz::TinyVector<int,NDIM>& patch_touches_lower_physical_bdry,
    const blitz::TinyVector<int,NDIM>& patch_touches_upper_physical_bdry,
    const std::vector<int>& local_indices,
    const std::vector<double>& periodic_offsets,
    const std::string& interp_fcn)
{
    if (local_indices.empty()) return;
    const int local_indices_size = local_indices.size();
    const IntVector<NDIM>& ilower = q_data_box.lower();
    const IntVector<NDIM>& iupper = q_data_box.upper();
    if (interp_fcn == "PIECEWISE_CONSTANT")
    {
        LAGRANGIAN_PIECEWISE_CONSTANT_INTERP_FC(
            dx,x_lower,x_upper,q_depth,
#if (NDIM == 2)
            ilower(0),iupper(0),ilower(1),iupper(1),
            q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            q_gcw(0),q_gcw(1),q_gcw(2),
#endif
            q_data,
            &local_indices[0], &periodic_offsets[0], local_indices_size,
            X_data, Q_data);
    }
    else if (interp_fcn == "PIECEWISE_LINEAR")
    {
#if (NDIM == 2)
        LAGRANGIAN_PIECEWISE_LINEAR_INTERP_FC(
            dx,x_lower,x_upper,q_depth,
#if (NDIM == 2)
            ilower(0),iupper(0),ilower(1),iupper(1),
            q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            q_gcw(0),q_gcw(1),q_gcw(2),
#endif
            q_data,
            &local_indices[0], &periodic_offsets[0], local_indices_size,
            X_data,Q_data);
#else
        TBOX_ERROR("LEInteractor::interpolate()\n" <<
                   "  Weighting function " << interp_fcn << " is only supported in 2D" << std::endl);
#endif
    }
    else if (interp_fcn == "WIDE_PIECEWISE_LINEAR")
    {
#if (NDIM == 2)
        LAGRANGIAN_WIDE_PIECEWISE_LINEAR_INTERP_FC(
            dx,x_lower,x_upper,q_depth,
#if (NDIM == 2)
            ilower(0),iupper(0),ilower(1),iupper(1),
            q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            q_gcw(0),q_gcw(1),q_gcw(2),
#endif
            q_data,
            &local_indices[0], &periodic_offsets[0], local_indices_size,
            X_data,Q_data);
#else
        TBOX_ERROR("LEInteractor::interpolate()\n" <<
                   "  Weighting function " << interp_fcn << " is only supported in 2D" << std::endl);
#endif
    }
    else if (interp_fcn == "PIECEWISE_CUBIC")
    {
#if (NDIM == 2)
        LAGRANGIAN_PIECEWISE_CUBIC_INTERP_FC(
            dx,x_lower,x_upper,q_depth,
#if (NDIM == 2)
            ilower(0),iupper(0),ilower(1),iupper(1),
            q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            q_gcw(0),q_gcw(1),q_gcw(2),
#endif
            q_data,
            &local_indices[0], &periodic_offsets[0], local_indices_size,
            X_data,Q_data);
#else
        TBOX_ERROR("LEInteractor::interpolate()\n" <<
                   "  Weighting function " << interp_fcn << " is only supported in 2D" << std::endl);
#endif
    }
    else if (interp_fcn == "WIDE_PIECEWISE_CUBIC")
    {
#if (NDIM == 2)
        LAGRANGIAN_WIDE_PIECEWISE_CUBIC_INTERP_FC(
            dx,x_lower,x_upper,q_depth,
#if (NDIM == 2)
            ilower(0),iupper(0),ilower(1),iupper(1),
            q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            q_gcw(0),q_gcw(1),q_gcw(2),
#endif
            q_data,
            &local_indices[0], &periodic_offsets[0], local_indices_size,
            X_data,Q_data);
#else
        TBOX_ERROR("LEInteractor::interpolate()\n" <<
                   "  Weighting function " << interp_fcn << " is only supported in 2D" << std::endl);
#endif
    }
    else if (interp_fcn == "IB_3")
    {
        LAGRANGIAN_IB_3_INTERP_FC(
            dx,x_lower,x_upper,q_depth,
#if (NDIM == 2)
            ilower(0),iupper(0),ilower(1),iupper(1),
            q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            q_gcw(0),q_gcw(1),q_gcw(2),
#endif
            q_data,
            &local_indices[0], &periodic_offsets[0], local_indices_size,
            X_data,Q_data);
    }
    else if (interp_fcn == "WIDE_IB_3")
    {
#if (NDIM == 2)
        LAGRANGIAN_WIDE_IB_3_INTERP_FC(
            dx,x_lower,x_upper,q_depth,
#if (NDIM == 2)
            ilower(0),iupper(0),ilower(1),iupper(1),
            q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            q_gcw(0),q_gcw(1),q_gcw(2),
#endif
            q_data,
            &local_indices[0], &periodic_offsets[0], local_indices_size,
            X_data,Q_data);
#else
        TBOX_ERROR("LEInteractor::interpolate()\n" <<
                   "  Weighting function " << interp_fcn << " is only supported in 2D" << std::endl);
#endif
    }
    else if (interp_fcn == "IB_4")
    {
        LAGRANGIAN_IB_4_INTERP_FC(
            dx,x_lower,x_upper,q_depth,
#if (NDIM == 2)
            ilower(0),iupper(0),ilower(1),iupper(1),
            &patch_touches_lower_physical_bdry[0],
            &patch_touches_upper_physical_bdry[0],
            q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            &patch_touches_lower_physical_bdry[0],
            &patch_touches_upper_physical_bdry[0],
            q_gcw(0),q_gcw(1),q_gcw(2),
#endif
            q_data,
            &local_indices[0], &periodic_offsets[0], local_indices_size,
            X_data,Q_data);
    }
    else if (interp_fcn == "WIDE_IB_4")
    {
#if (NDIM == 2)
        LAGRANGIAN_WIDE_IB_4_INTERP_FC(
            dx,x_lower,x_upper,q_depth,
#if (NDIM == 2)
            ilower(0),iupper(0),ilower(1),iupper(1),
            q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            q_gcw(0),q_gcw(1),q_gcw(2),
#endif
            q_data,
            &local_indices[0], &periodic_offsets[0], local_indices_size,
            X_data,Q_data);
#else
        TBOX_ERROR("LEInteractor::interpolate()\n" <<
                   "  Weighting function " << interp_fcn << " is only supported in 2D" << std::endl);
#endif
    }
    else if (interp_fcn == "IB_6")
    {
        LAGRANGIAN_IB_6_INTERP_FC(
            dx,x_lower,x_upper,q_depth,
#if (NDIM == 2)
            ilower(0),iupper(0),ilower(1),iupper(1),
            &patch_touches_lower_physical_bdry[0],
            &patch_touches_upper_physical_bdry[0],
            q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            &patch_touches_lower_physical_bdry[0],
            &patch_touches_upper_physical_bdry[0],
            q_gcw(0),q_gcw(1),q_gcw(2),
#endif
            q_data,
            &local_indices[0], &periodic_offsets[0], local_indices_size,
            X_data,Q_data);
    }
    else if (interp_fcn == "USER_DEFINED")
    {
        userDefinedInterpolate(
            Q_data, Q_depth, X_data,
            q_data, q_data_box, q_gcw, q_depth,
            x_lower, x_upper, dx,
            &local_indices[0], &periodic_offsets[0], local_indices_size);
    }
    else
    {
        TBOX_ERROR("LEInteractor::interpolate()\n" <<
                   "  Unknown interpolation weighting function "
                   << interp_fcn << std::endl);
    }
    return;
}// interpolate

void
LEInteractor::spread(
    double* const q_data,
    const Box<NDIM>& q_data_box,
    const IntVector<NDIM>& q_gcw,
    const int q_depth,
    const double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const double* const x_lower,
    const double* const x_upper,
    const double* const dx,
    const blitz::TinyVector<int,NDIM>& patch_touches_lower_physical_bdry,
    const blitz::TinyVector<int,NDIM>& patch_touches_upper_physical_bdry,
    const std::vector<int>& local_indices,
    const std::vector<double>& periodic_offsets,
    const std::string& spread_fcn)
{
    if (local_indices.empty()) return;
    const int local_indices_size = local_indices.size();
    const IntVector<NDIM>& ilower = q_data_box.lower();
    const IntVector<NDIM>& iupper = q_data_box.upper();
    if (spread_fcn == "PIECEWISE_CONSTANT")
    {
        if (s_precision_mode == DOUBLE)
        {
            LAGRANGIAN_PIECEWISE_CONSTANT_SPREAD_FC(
                dx,x_lower,x_upper,q_depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data, Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data);
        }
        else if (s_precision_mode == DOUBLE_DOUBLE)
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  extended precision not currently supported for PIECEWISE_CONSTANT delta function.\n");
        }
        else
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  invalid precision mode; s_precision_mode = " << s_precision_mode << ".\n");
        }
    }
    else if (spread_fcn == "PIECEWISE_LINEAR")
    {
#if (NDIM == 2)
        if (s_precision_mode == DOUBLE)
        {
            LAGRANGIAN_PIECEWISE_LINEAR_SPREAD_FC(
                dx,x_lower,x_upper,q_depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data, Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data);
        }
        else if (s_precision_mode == DOUBLE_DOUBLE)
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  extended precision not currently supported for PIECEWISE_LINEAR delta function.\n");
        }
        else
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  invalid precision mode; s_precision_mode = " << s_precision_mode << ".\n");
        }
#else
        TBOX_ERROR("LEInteractor::spread()\n" <<
                   "  Weighting function " << spread_fcn << " is only supported in 2D" << std::endl);
#endif
    }
    else if (spread_fcn == "WIDE_PIECEWISE_LINEAR")
    {
#if (NDIM == 2)
        if (s_precision_mode == DOUBLE)
        {
            LAGRANGIAN_WIDE_PIECEWISE_LINEAR_SPREAD_FC(
                dx,x_lower,x_upper,q_depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data, Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data);
        }
        else if (s_precision_mode == DOUBLE_DOUBLE)
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  extended precision not currently supported for WIDE_PIECEWISE_LINEAR delta function.\n");
        }
        else
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  invalid precision mode; s_precision_mode = " << s_precision_mode << ".\n");
        }
#else
        TBOX_ERROR("LEInteractor::spread()\n" <<
                   "  Weighting function " << spread_fcn << " is only supported in 2D" << std::endl);
#endif
    }
    else if (spread_fcn == "PIECEWISE_CUBIC")
    {
#if (NDIM == 2)
        if (s_precision_mode == DOUBLE)
        {
            LAGRANGIAN_PIECEWISE_CUBIC_SPREAD_FC(
                dx,x_lower,x_upper,q_depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data, Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data);
        }
        else if (s_precision_mode == DOUBLE_DOUBLE)
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  extended precision not currently supported for PIECEWISE_CUBIC delta function.\n");
        }
        else
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  invalid precision mode; s_precision_mode = " << s_precision_mode << ".\n");
        }
#else
        TBOX_ERROR("LEInteractor::spread()\n" <<
                   "  Weighting function " << spread_fcn << " is only supported in 2D" << std::endl);
#endif
    }
    else if (spread_fcn == "WIDE_PIECEWISE_CUBIC")
    {
#if (NDIM == 2)
        if (s_precision_mode == DOUBLE)
        {
            LAGRANGIAN_WIDE_PIECEWISE_CUBIC_SPREAD_FC(
                dx,x_lower,x_upper,q_depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data, Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data);
        }
        else if (s_precision_mode == DOUBLE_DOUBLE)
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  extended precision not currently supported for WIDE_PIECEWISE_CUBIC delta function.\n");
        }
        else
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  invalid precision mode; s_precision_mode = " << s_precision_mode << ".\n");
        }
#else
        TBOX_ERROR("LEInteractor::spread()\n" <<
                   "  Weighting function " << spread_fcn << " is only supported in 2D" << std::endl);
#endif
    }
    else if (spread_fcn == "IB_3")
    {
        if (s_precision_mode == DOUBLE)
        {
            LAGRANGIAN_IB_3_SPREAD_FC(
                dx,x_lower,x_upper,q_depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data, Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data);
        }
        else if (s_precision_mode == DOUBLE_DOUBLE)
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  extended precision not currently supported for IB_3 delta function.\n");
        }
        else
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  invalid precision mode; s_precision_mode = " << s_precision_mode << ".\n");
        }
    }
    else if (spread_fcn == "WIDE_IB_3")
    {
#if (NDIM == 2)
        if (s_precision_mode == DOUBLE)
        {
            LAGRANGIAN_WIDE_IB_3_SPREAD_FC(
                dx,x_lower,x_upper,q_depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data, Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data);
        }
        else if (s_precision_mode == DOUBLE_DOUBLE)
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  extended precision not currently supported for WIDE_IB_3 delta function.\n");
        }
        else
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  invalid precision mode; s_precision_mode = " << s_precision_mode << ".\n");
        }
#else
        TBOX_ERROR("LEInteractor::spread()\n" <<
                   "  Weighting function " << spread_fcn << " is only supported in 2D" << std::endl);
#endif
    }
    else if (spread_fcn == "IB_4")
    {
        if (s_precision_mode == DOUBLE)
        {
            LAGRANGIAN_IB_4_SPREAD_FC(
                dx,x_lower,x_upper,q_depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data, Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                &patch_touches_lower_physical_bdry[0],
                &patch_touches_upper_physical_bdry[0],
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                &patch_touches_lower_physical_bdry[0],
                &patch_touches_upper_physical_bdry[0],
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data);
        }
        else if (s_precision_mode == DOUBLE_DOUBLE)
        {
            LAGRANGIAN_IB_4_SPREAD_XP_FC(
                dx,x_lower,x_upper,q_depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data, Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                &patch_touches_lower_physical_bdry[0],
                &patch_touches_upper_physical_bdry[0],
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                &patch_touches_lower_physical_bdry[0],
                &patch_touches_upper_physical_bdry[0],
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data);
        }
        else
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  invalid precision mode; s_precision_mode = " << s_precision_mode << ".\n");
        }
    }
    else if (spread_fcn == "WIDE_IB_4")
    {
#if (NDIM == 2)
        if (s_precision_mode == DOUBLE)
        {
            LAGRANGIAN_WIDE_IB_4_SPREAD_FC(
                dx,x_lower,x_upper,q_depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data, Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data);
        }
        else if (s_precision_mode == DOUBLE_DOUBLE)
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  extended precision not currently supported for WIDE_IB_4 delta function.\n");
        }
        else
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  invalid precision mode; s_precision_mode = " << s_precision_mode << ".\n");
        }
#else
        TBOX_ERROR("LEInteractor::spread()\n" <<
                   "  Weighting function " << spread_fcn << " is only supported in 2D" << std::endl);
#endif
    }
    else if (spread_fcn == "IB_6")
    {
        if (s_precision_mode == DOUBLE)
        {
            LAGRANGIAN_IB_6_SPREAD_FC(
                dx,x_lower,x_upper,q_depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data,Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                &patch_touches_lower_physical_bdry[0],
                &patch_touches_upper_physical_bdry[0],
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                &patch_touches_lower_physical_bdry[0],
                &patch_touches_upper_physical_bdry[0],
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data);
        }
        else if (s_precision_mode == DOUBLE_DOUBLE)
        {
            LAGRANGIAN_IB_6_SPREAD_XP_FC(
                dx,x_lower,x_upper,q_depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data,Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                &patch_touches_lower_physical_bdry[0],
                &patch_touches_upper_physical_bdry[0],
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                &patch_touches_lower_physical_bdry[0],
                &patch_touches_upper_physical_bdry[0],
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data);
        }
        else
        {
            TBOX_ERROR("LEInteractor::spread():\n"
                       << "  invalid precision mode; s_precision_mode = " << s_precision_mode << ".\n");
        }
    }
    else if (spread_fcn == "USER_DEFINED")
    {
        userDefinedSpread(
            q_data, q_data_box, q_gcw, q_depth,
            x_lower, x_upper, dx,
            Q_data, Q_depth, X_data,
            &local_indices[0], &periodic_offsets[0], local_indices_size);
    }
    else
    {
        TBOX_ERROR("LEInteractor::spread()\n" <<
                   "  Unknown spreading weighting function "
                   << spread_fcn << std::endl);
    }
    return;
}// spread

template<class T>
void
LEInteractor::buildLocalIndices(
    std::vector<int>& local_indices,
    std::vector<double>& periodic_offsets,
    const Box<NDIM>& box,
    const Pointer<Patch<NDIM> > patch,
    const IntVector<NDIM>& periodic_shift,
    const Pointer<LIndexSetData<T> > idx_data)
{
    local_indices   .clear();
    periodic_offsets.clear();
    const unsigned int upper_bound = idx_data->getInteriorLocalPETScIndices().size() + idx_data->getGhostLocalPETScIndices().size();
    if (upper_bound == 0) return;
    local_indices   .reserve(     upper_bound);
    periodic_offsets.reserve(NDIM*upper_bound);

    const Box<NDIM>& patch_box = patch->getBox();
    const Index<NDIM>& ilower = patch_box.lower();
    const Index<NDIM>& iupper = patch_box.upper();
    const Box<NDIM>& ghost_box = idx_data->getGhostBox();

    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    blitz::TinyVector<bool,NDIM> patch_touches_lower_periodic_bdry, patch_touches_upper_periodic_bdry;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        patch_touches_lower_periodic_bdry[axis] = pgeom->getTouchesPeriodicBoundary(axis,0);
        patch_touches_upper_periodic_bdry[axis] = pgeom->getTouchesPeriodicBoundary(axis,1);
    }

    blitz::TinyVector<int   ,NDIM>      offset;
    blitz::TinyVector<double,NDIM> node_offset;
    if (s_sort_mode == NO_SORT)
    {
        if (box == patch_box)
        {
            local_indices = idx_data->getInteriorLocalPETScIndices();
            periodic_offsets.resize(NDIM*local_indices.size(),0.0);
        }
        else if (box == ghost_box)
        {
            local_indices = idx_data->getInteriorLocalPETScIndices();
            periodic_offsets.resize(NDIM*local_indices.size(),0.0);
            local_indices   .insert(local_indices   .end(), idx_data->getGhostLocalPETScIndices().begin(), idx_data->getGhostLocalPETScIndices().end());
            periodic_offsets.insert(periodic_offsets.end(), idx_data->getGhostPeriodicOffsets  ().begin(), idx_data->getGhostPeriodicOffsets  ().end());
        }
        else
        {
            for (typename LIndexSetData<T>::SetIterator it(*idx_data); it; it++)
            {
                const Index<NDIM>& i = it.getIndex();
                if (!box.contains(i)) continue;

                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    if      (patch_touches_lower_periodic_bdry[d] && i(d) < ilower(d))
                    {
                        offset[d] = -periodic_shift(d);  // X is ABOVE the top    of the patch --- need to shift DOWN
                    }
                    else if (patch_touches_upper_periodic_bdry[d] && i(d) > iupper(d))
                    {
                        offset[d] = +periodic_shift(d);  // X is BELOW the bottom of the patch --- need to shift UP
                    }
                    else
                    {
                        offset[d] = 0;
                    }
                }

                const LSet<T>& node_set = it.getItem();
                const typename LSet<T>::size_type num_ids = node_set.size();
                const unsigned int old_size = local_indices.size();
                const unsigned int new_size = old_size+num_ids;
                local_indices   .resize(     new_size);
                periodic_offsets.resize(NDIM*new_size);
                for (typename LSet<T>::size_type n = 0; n < num_ids; ++n)
                {
                    local_indices[old_size+n] = node_set[n]->getLocalPETScIndex();
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        periodic_offsets[NDIM*(old_size+n)+d] = static_cast<double>(offset[d])*dx[d];
                    }
                }
            }
        }
    }
    else
    {
        if (s_sort_mode != SORT_INCREASING_LAG_IDX && s_sort_mode != SORT_DECREASING_LAG_IDX)
        {
            TBOX_ERROR("LEInteractor::buildLocalIndices():\n"
                       << "  invalid debug sort mode; s_sort_mode = " << s_sort_mode << ".\n");
        }

        std::vector<std::pair<const LNodeIndex*,blitz::TinyVector<double,NDIM> > > box_idxs_and_offsets;
        box_idxs_and_offsets.reserve(upper_bound);
        for (typename LIndexSetData<T> ::SetIterator it(*idx_data); it; it++)
        {
            const Index<NDIM>& i = it.getIndex();
            if (!box.contains(i)) continue;

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if      (patch_touches_lower_periodic_bdry[d] && i(d) < ilower(d))
                {
                    offset[d] = -periodic_shift(d);  // X is ABOVE the top    of the patch --- need to shift DOWN
                }
                else if (patch_touches_upper_periodic_bdry[d] && i(d) > iupper(d))
                {
                    offset[d] = +periodic_shift(d);  // X is BELOW the bottom of the patch --- need to shift UP
                }
                else
                {
                    offset[d] = 0;
                }
                node_offset[d] = static_cast<double>(offset[d])*dx[d];
            }

            const LSet<T>& node_set = it.getItem();
            for (typename LSet<T>::const_iterator cit = node_set.begin(); cit != node_set.end(); ++cit)
            {
                box_idxs_and_offsets.push_back(std::make_pair(cit->getPointer(),node_offset));
            }

        }

        std::sort(box_idxs_and_offsets.begin(),box_idxs_and_offsets.end(),SortModeComp());
        for (unsigned int n = 0; n < box_idxs_and_offsets.size(); ++n)
        {
            local_indices.push_back(box_idxs_and_offsets[n].first->getLocalPETScIndex());
            periodic_offsets.insert(periodic_offsets.end(),box_idxs_and_offsets[n].second.data(),box_idxs_and_offsets[n].second.data()+NDIM);
        }
    }
    return;
}// buildLocalIndices

void
LEInteractor::buildLocalIndices(
    std::vector<int>& local_indices,
    const Box<NDIM>& box,
    const Pointer<Patch<NDIM> > patch,
    const double* const X_data,
    const int X_size,
    const int X_depth)
{
    local_indices.clear();
    const int upper_bound = X_size/X_depth;
    if (upper_bound == 0) return;

    const Box<NDIM>& patch_box = patch->getBox();
    const CellIndex<NDIM>& patch_lower = patch_box.lower();
    const CellIndex<NDIM>& patch_upper = patch_box.upper();

    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    const double* const patch_x_lower = patch_geom->getXLower();
    const double* const patch_x_upper = patch_geom->getXUpper();
    const double* const patch_dx = patch_geom->getDx();

    local_indices.reserve(upper_bound);
    for (int k = 0; k < X_size/X_depth; ++k)
    {
        const double* const X = &X_data[NDIM*k];
        const Index<NDIM> i = IndexUtilities::getCellIndex(X, patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper);
        if (box.contains(i)) local_indices.push_back(k);
    }
    return;
}// buildLocalIndices

void
LEInteractor::userDefinedInterpolate(
    double* Q,
    const int Q_depth,
    const double* const X,
    const double* const q,
    const Box<NDIM>& q_data_box,
    const int* const q_gcw,
    const int q_depth,
    const double* const x_lower,
    const double* const /*x_upper*/,
    const double* const dx,
    const int* const local_indices,
    const double* const X_shift,
    const int num_local_indices)
{
    const int* const ilower = q_data_box.lower();
    const int* const iupper = q_data_box.upper();
    const blitz::TinyVector<int,NDIM+1> shape(q_data_box.numberCells(0)+2*q_gcw[0],
                                              q_data_box.numberCells(1)+2*q_gcw[1],
#if (NDIM == 3)
                                              q_data_box.numberCells(2)+2*q_gcw[2],
#endif
                                              q_depth);
    blitz::GeneralArrayStorage<NDIM+1> storage_order = blitz::ColumnMajorArray<NDIM+1>();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        storage_order.base()[d] = ilower[d]-q_gcw[d];
    }
    const blitz::Array<double,NDIM+1> q_data(const_cast<double*>(q), shape, blitz::neverDeleteData, storage_order);

    blitz::TinyVector<double,NDIM> X_cell;
    blitz::TinyVector<int,NDIM> stencil_center, stencil_lower, stencil_upper;
    for (int l = 0; l < num_local_indices; ++l)
    {
        const int s = local_indices[l];

        // Determine the Cartesian cell in which X(s) is located.
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            stencil_center[d] = static_cast<int>(std::floor((X[d+s*NDIM]+X_shift[d+l*NDIM]-x_lower[d])/dx[d])) + ilower[d];
            X_cell[d] = x_lower[d]+(static_cast<double>(stencil_center[d]-ilower[d])+0.5)*dx[d];
        }

        // Determine the interpolation stencil corresponding to the position of
        // X(s) within the cell.
        if (s_delta_fcn_stencil_size%2 == 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (X[d+s*NDIM] < X_cell[d])
                {
                    stencil_lower[d] = stencil_center[d]-s_delta_fcn_stencil_size/2;
                    stencil_upper[d] = stencil_center[d]+s_delta_fcn_stencil_size/2-1;
                }
                else
                {
                    stencil_lower[d] = stencil_center[d]-s_delta_fcn_stencil_size/2+1;
                    stencil_upper[d] = stencil_center[d]+s_delta_fcn_stencil_size/2;
                }
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                stencil_lower[d] = stencil_center[d]-s_delta_fcn_stencil_size/2;
                stencil_upper[d] = stencil_center[d]+s_delta_fcn_stencil_size/2;
            }
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            stencil_lower[d] = std::min(std::max(stencil_lower[d],ilower[d]-q_gcw[d]),iupper[d]+q_gcw[d]);
            stencil_upper[d] = std::min(std::max(stencil_upper[d],ilower[d]-q_gcw[d]),iupper[d]+q_gcw[d]);
        }

        // Compute the delta function weights.
        blitz::Array<double,1> w0(blitz::Range(stencil_lower[0],stencil_upper[0]));
        for (int ic0 = stencil_lower[0]; ic0 <= stencil_upper[0]; ++ic0)
        {
            w0(ic0) = s_delta_fcn((X[0+s*NDIM]+X_shift[0+l*NDIM]-(X_cell[0]+static_cast<double>(ic0-stencil_center[0])*dx[0]))/dx[0]);
        }

        blitz::Array<double,1> w1(blitz::Range(stencil_lower[1],stencil_upper[1]));
        for (int ic1 = stencil_lower[1]; ic1 <= stencil_upper[1]; ++ic1)
        {
            w1(ic1) = s_delta_fcn((X[1+s*NDIM]+X_shift[1+l*NDIM]-(X_cell[1]+static_cast<double>(ic1-stencil_center[1])*dx[1]))/dx[1]);
        }
#if (NDIM == 3)
        blitz::Array<double,1> w2(blitz::Range(stencil_lower[2],stencil_upper[2]));
        for (int ic2 = stencil_lower[2]; ic2 <= stencil_upper[2]; ++ic2)
        {
            w2(ic2) = s_delta_fcn((X[2+s*NDIM]+X_shift[2+l*NDIM]-(X_cell[2]+static_cast<double>(ic2-stencil_center[2])*dx[2]))/dx[2]);
        }
#endif
        // Interpolate u onto V.
        for (int d = 0; d < Q_depth; ++d)
        {
            Q[d+s*Q_depth] = 0.0;
#if (NDIM == 3)
            for (int ic2 = stencil_lower[2]; ic2 <= stencil_upper[2]; ++ic2)
            {
#endif
                for (int ic1 = stencil_lower[1]; ic1 <= stencil_upper[1]; ++ic1)
                {
                    for (int ic0 = stencil_lower[0]; ic0 <= stencil_upper[0]; ++ic0)
                    {
#if (NDIM == 2)
                        Q[d+s*Q_depth] += w0(ic0)*w1(ic1)*q_data(ic0,ic1,d);
#endif
#if (NDIM == 3)
                        Q[d+s*Q_depth] += w0(ic0)*w1(ic1)*w2(ic2)*q_data(ic0,ic1,ic2,d);
#endif
                    }
                }
#if (NDIM == 3)
            }
#endif
        }
    }
    return;
}// userDefinedInterpolate

void
LEInteractor::userDefinedSpread(
    double* q,
    const Box<NDIM>& q_data_box,
    const int* const q_gcw,
    const int q_depth,
    const double* const x_lower,
    const double* const /*x_upper*/,
    const double* const dx,
    const double* const Q,
    const int Q_depth,
    const double* const X,
    const int* const local_indices,
    const double* const X_shift,
    const int num_local_indices)
{
    const int* const ilower = q_data_box.lower();
    const int* const iupper = q_data_box.upper();
    const blitz::TinyVector<int,NDIM+1> shape(q_data_box.numberCells(0)+2*q_gcw[0],
                                              q_data_box.numberCells(1)+2*q_gcw[1],
#if (NDIM == 3)
                                              q_data_box.numberCells(2)+2*q_gcw[2],
#endif
                                              q_depth);
    blitz::GeneralArrayStorage<NDIM+1> storage_order = blitz::ColumnMajorArray<NDIM+1>();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        storage_order.base()[d] = ilower[d]-q_gcw[d];
    }
    blitz::Array<double,NDIM+1> q_data(q, shape, blitz::neverDeleteData, storage_order);

    blitz::TinyVector<double,NDIM> X_cell;
    blitz::TinyVector<int,NDIM> stencil_center, stencil_lower, stencil_upper;
    for (int l = 0; l < num_local_indices; ++l)
    {
        const int s = local_indices[l];

        // Determine the Cartesian cell in which X(s) is located.
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            stencil_center[d] = static_cast<int>(std::floor((X[d+s*NDIM]+X_shift[d+l*NDIM]-x_lower[d])/dx[d])) + ilower[d];
            X_cell[d] = x_lower[d]+(static_cast<double>(stencil_center[d]-ilower[d])+0.5)*dx[d];
        }

        // Determine the interpolation stencil corresponding to the position of
        // X(s) within the cell.
        if (s_delta_fcn_stencil_size%2 == 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (X[d+s*NDIM] < X_cell[d])
                {
                    stencil_lower[d] = stencil_center[d]-s_delta_fcn_stencil_size/2;
                    stencil_upper[d] = stencil_center[d]+s_delta_fcn_stencil_size/2-1;
                }
                else
                {
                    stencil_lower[d] = stencil_center[d]-s_delta_fcn_stencil_size/2+1;
                    stencil_upper[d] = stencil_center[d]+s_delta_fcn_stencil_size/2;
                }
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                stencil_lower[d] = stencil_center[d]-s_delta_fcn_stencil_size/2;
                stencil_upper[d] = stencil_center[d]+s_delta_fcn_stencil_size/2;
            }
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            stencil_lower[d] = std::min(std::max(stencil_lower[d],ilower[d]-q_gcw[d]),iupper[d]+q_gcw[d]);
            stencil_upper[d] = std::min(std::max(stencil_upper[d],ilower[d]-q_gcw[d]),iupper[d]+q_gcw[d]);
        }


        // Compute the delta function weights.
        blitz::Array<double,1> w0(blitz::Range(stencil_lower[0],stencil_upper[0]));
        for (int ic0 = stencil_lower[0]; ic0 <= stencil_upper[0]; ++ic0)
        {
            w0(ic0) = s_delta_fcn((X[0+s*NDIM]+X_shift[0+l*NDIM]-(X_cell[0]+static_cast<double>(ic0-stencil_center[0])*dx[0]))/dx[0]);
        }

        blitz::Array<double,1> w1(blitz::Range(stencil_lower[1],stencil_upper[1]));
        for (int ic1 = stencil_lower[1]; ic1 <= stencil_upper[1]; ++ic1)
        {
            w1(ic1) = s_delta_fcn((X[1+s*NDIM]+X_shift[1+l*NDIM]-(X_cell[1]+static_cast<double>(ic1-stencil_center[1])*dx[1]))/dx[1]);
        }
#if (NDIM == 3)
        blitz::Array<double,1> w2(blitz::Range(stencil_lower[2],stencil_upper[2]));
        for (int ic2 = stencil_lower[2]; ic2 <= stencil_upper[2]; ++ic2)
        {
            w2(ic2) = s_delta_fcn((X[2+s*NDIM]+X_shift[2+l*NDIM]-(X_cell[2]+static_cast<double>(ic2-stencil_center[2])*dx[2]))/dx[2]);
        }
#endif
        // Spread V onto u.
        for (int d = 0; d < Q_depth; ++d)
        {
#if (NDIM == 3)
            for (int ic2 = stencil_lower[2]; ic2 <= stencil_upper[2]; ++ic2)
            {
#endif
                for (int ic1 = stencil_lower[1]; ic1 <= stencil_upper[1]; ++ic1)
                {
                    for (int ic0 = stencil_lower[0]; ic0 <= stencil_upper[0]; ++ic0)
                    {
#if (NDIM == 2)
                        q_data(ic0,ic1,d) += w0(ic0)*w1(ic1)*Q[d+s*Q_depth]/(dx[0]*dx[1]);
#endif
#if (NDIM == 3)
                        q_data(ic0,ic1,ic2,d) += w0(ic0)*w1(ic1)*w2(ic2)*Q[d+s*Q_depth]/(dx[0]*dx[1]*dx[2]);
#endif
                    }
                }
#if (NDIM == 3)
            }
#endif
        }
    }
    return;
}// userDefinedSpread

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <ibtk/LNode.h>

template void IBTK::LEInteractor::interpolate(
    SAMRAI::tbox::Pointer<LData> Q_data,
    const SAMRAI::tbox::Pointer<LData> X_data,
    const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const SAMRAI::hier::Box<NDIM>& interp_box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& interp_fcn);

template void IBTK::LEInteractor::interpolate(
    SAMRAI::tbox::Pointer<LData> Q_data,
    const SAMRAI::tbox::Pointer<LData> X_data,
    const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM,double> > q_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const SAMRAI::hier::Box<NDIM>& interp_box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& interp_fcn);

template void IBTK::LEInteractor::interpolate(
    SAMRAI::tbox::Pointer<LData> Q_data,
    const SAMRAI::tbox::Pointer<LData> X_data,
    const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > q_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const SAMRAI::hier::Box<NDIM>& interp_box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& interp_fcn);

template void IBTK::LEInteractor::interpolate(
    double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const SAMRAI::hier::Box<NDIM>& interp_box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& interp_fcn);

template void IBTK::LEInteractor::interpolate(
    double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM,double> > q_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const SAMRAI::hier::Box<NDIM>& interp_box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& interp_fcn);

template void IBTK::LEInteractor::interpolate(
    double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > q_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const SAMRAI::hier::Box<NDIM>& interp_box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& interp_fcn);

template void IBTK::LEInteractor::spread(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
    const SAMRAI::tbox::Pointer<LData> Q_data,
    const SAMRAI::tbox::Pointer<LData> X_data,
    const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const SAMRAI::hier::Box<NDIM>& spread_box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& spread_fcn);

template void IBTK::LEInteractor::spread(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM,double> > q_data,
    const SAMRAI::tbox::Pointer<LData> Q_data,
    const SAMRAI::tbox::Pointer<LData> X_data,
    const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const SAMRAI::hier::Box<NDIM>& spread_box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& spread_fcn);

template void IBTK::LEInteractor::spread(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > q_data,
    const SAMRAI::tbox::Pointer<LData> Q_data,
    const SAMRAI::tbox::Pointer<LData> X_data,
    const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const SAMRAI::hier::Box<NDIM>& spread_box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& spread_fcn);

template void IBTK::LEInteractor::spread(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
    const double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const SAMRAI::hier::Box<NDIM>& spread_box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& spread_fcn);

template void IBTK::LEInteractor::spread(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM,double> > q_data,
    const double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const SAMRAI::hier::Box<NDIM>& spread_box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& spread_fcn);

template void IBTK::LEInteractor::spread(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > q_data,
    const double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const SAMRAI::hier::Box<NDIM>& spread_box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& spread_fcn);

template void IBTK::LEInteractor::buildLocalIndices(
    std::vector<int>& local_indices,
    std::vector<double>& periodic_offsets,
    const SAMRAI::hier::Box<NDIM>& box,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data);

//////////////////////////////////////////////////////////////////////////////
