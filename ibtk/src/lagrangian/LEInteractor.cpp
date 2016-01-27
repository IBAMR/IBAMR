// Filename: LEInteractor.cpp
// Created on 14 Jul 2004 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <algorithm>
#include <cmath>
#include <ostream>
#include <string>
#include <vector>

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "EdgeData.h"
#include "EdgeGeometry.h"
#include "IBTK_config.h"
#include "Index.h"
#include "IntVector.h"
#include "NodeData.h"
#include "NodeGeometry.h"
#include "Patch.h"
#include "SideData.h"
#include "SideGeometry.h"
#include <Eigen/Dense>
#include "boost/array.hpp"
#include "boost/multi_array.hpp"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LData.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/LIndexSetData.h"
#include "ibtk/LSet.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define LAGRANGIAN_PIECEWISE_CONSTANT_INTERP_FC                                                                        \
    IBTK_FC_FUNC_(lagrangian_piecewise_constant_interp2d, LAGRANGIAN_PIECEWISE_CONSTANT_INTERP2D)
#define LAGRANGIAN_PIECEWISE_CONSTANT_SPREAD_FC                                                                        \
    IBTK_FC_FUNC_(lagrangian_piecewise_constant_spread2d, LAGRANGIAN_PIECEWISE_CONSTANT_SPREAD2D)

#define LAGRANGIAN_DISCONTINUOUS_LINEAR_INTERP_FC                                                                      \
    IBTK_FC_FUNC_(lagrangian_discontinuous_linear_interp2d, LAGRANGIAN_DISCONTINUOUS_LINEAR_INTERP2D)
#define LAGRANGIAN_DISCONTINUOUS_LINEAR_SPREAD_FC                                                                      \
    IBTK_FC_FUNC_(lagrangian_discontinuous_linear_spread2d, LAGRANGIAN_DISCONTINUOUS_LINEAR_SPREAD2D)

#define LAGRANGIAN_PIECEWISE_LINEAR_INTERP_FC                                                                          \
    IBTK_FC_FUNC_(lagrangian_piecewise_linear_interp2d, LAGRANGIAN_PIECEWISE_LINEAR_INTERP2D)
#define LAGRANGIAN_PIECEWISE_LINEAR_SPREAD_FC                                                                          \
    IBTK_FC_FUNC_(lagrangian_piecewise_linear_spread2d, LAGRANGIAN_PIECEWISE_LINEAR_SPREAD2D)

#define LAGRANGIAN_PIECEWISE_CUBIC_INTERP_FC                                                                           \
    IBTK_FC_FUNC_(lagrangian_piecewise_cubic_interp2d, LAGRANGIAN_PIECEWISE_CUBIC_INTERP2D)
#define LAGRANGIAN_PIECEWISE_CUBIC_SPREAD_FC                                                                           \
    IBTK_FC_FUNC_(lagrangian_piecewise_cubic_spread2d, LAGRANGIAN_PIECEWISE_CUBIC_SPREAD2D)

#define LAGRANGIAN_IB_3_INTERP_FC IBTK_FC_FUNC_(lagrangian_ib_3_interp2d, LAGRANGIAN_IB_3_INTERP2D)
#define LAGRANGIAN_IB_3_SPREAD_FC IBTK_FC_FUNC_(lagrangian_ib_3_spread2d, LAGRANGIAN_IB_3_SPREAD2D)

#define LAGRANGIAN_IB_4_INTERP_FC IBTK_FC_FUNC_(lagrangian_ib_4_interp2d, LAGRANGIAN_IB_4_INTERP2D)
#define LAGRANGIAN_IB_4_SPREAD_FC IBTK_FC_FUNC_(lagrangian_ib_4_spread2d, LAGRANGIAN_IB_4_SPREAD2D)

#define LAGRANGIAN_IB_4_W8_INTERP_FC IBTK_FC_FUNC_(lagrangian_ib_4_w8_interp2d, LAGRANGIAN_IB_4_W8_INTERP2D)
#define LAGRANGIAN_IB_4_W8_SPREAD_FC IBTK_FC_FUNC_(lagrangian_ib_4_w8_spread2d, LAGRANGIAN_IB_4_W8_SPREAD2D)

#define LAGRANGIAN_IB_6_INTERP_FC IBTK_FC_FUNC_(lagrangian_ib_6_interp2d, LAGRANGIAN_IB_6_INTERP2D)
#define LAGRANGIAN_IB_6_SPREAD_FC IBTK_FC_FUNC_(lagrangian_ib_6_spread2d, LAGRANGIAN_IB_6_SPREAD2D)
#endif

#if (NDIM == 3)
#define LAGRANGIAN_PIECEWISE_CONSTANT_INTERP_FC                                                                        \
    IBTK_FC_FUNC_(lagrangian_piecewise_constant_interp3d, LAGRANGIAN_PIECEWISE_CONSTANT_INTERP3D)
#define LAGRANGIAN_PIECEWISE_CONSTANT_SPREAD_FC                                                                        \
    IBTK_FC_FUNC_(lagrangian_piecewise_constant_spread3d, LAGRANGIAN_PIECEWISE_CONSTANT_SPREAD3D)

#define LAGRANGIAN_DISCONTINUOUS_LINEAR_INTERP_FC                                                                      \
    IBTK_FC_FUNC_(lagrangian_discontinuous_linear_interp3d, LAGRANGIAN_DISCONTINUOUS_LINEAR_INTERP3D)
#define LAGRANGIAN_DISCONTINUOUS_LINEAR_SPREAD_FC                                                                      \
    IBTK_FC_FUNC_(lagrangian_discontinuous_linear_spread3d, LAGRANGIAN_DISCONTINUOUS_LINEAR_SPREAD3D)

#define LAGRANGIAN_PIECEWISE_LINEAR_INTERP_FC                                                                          \
    IBTK_FC_FUNC_(lagrangian_piecewise_linear_interp3d, LAGRANGIAN_PIECEWISE_LINEAR_INTERP3D)
#define LAGRANGIAN_PIECEWISE_LINEAR_SPREAD_FC                                                                          \
    IBTK_FC_FUNC_(lagrangian_piecewise_linear_spread3d, LAGRANGIAN_PIECEWISE_LINEAR_SPREAD3D)

#define LAGRANGIAN_PIECEWISE_CUBIC_INTERP_FC                                                                           \
    IBTK_FC_FUNC_(lagrangian_piecewise_cubic_interp3d, LAGRANGIAN_PIECEWISE_CUBIC_INTERP3D)
#define LAGRANGIAN_PIECEWISE_CUBIC_SPREAD_FC                                                                           \
    IBTK_FC_FUNC_(lagrangian_piecewise_cubic_spread3d, LAGRANGIAN_PIECEWISE_CUBIC_SPREAD3D)

#define LAGRANGIAN_IB_3_INTERP_FC IBTK_FC_FUNC_(lagrangian_ib_3_interp3d, LAGRANGIAN_IB_3_INTERP3D)
#define LAGRANGIAN_IB_3_SPREAD_FC IBTK_FC_FUNC_(lagrangian_ib_3_spread3d, LAGRANGIAN_IB_3_SPREAD3D)

#define LAGRANGIAN_IB_4_INTERP_FC IBTK_FC_FUNC_(lagrangian_ib_4_interp3d, LAGRANGIAN_IB_4_INTERP3D)
#define LAGRANGIAN_IB_4_SPREAD_FC IBTK_FC_FUNC_(lagrangian_ib_4_spread3d, LAGRANGIAN_IB_4_SPREAD3D)

#define LAGRANGIAN_IB_4_W8_INTERP_FC IBTK_FC_FUNC_(lagrangian_ib_4_w8_interp3d, LAGRANGIAN_IB_4_W8_INTERP3D)
#define LAGRANGIAN_IB_4_W8_SPREAD_FC IBTK_FC_FUNC_(lagrangian_ib_4_w8_spread3d, LAGRANGIAN_IB_4_W8_SPREAD3D)

#define LAGRANGIAN_IB_6_INTERP_FC IBTK_FC_FUNC_(lagrangian_ib_6_interp3d, LAGRANGIAN_IB_6_INTERP3D)
#define LAGRANGIAN_IB_6_SPREAD_FC IBTK_FC_FUNC_(lagrangian_ib_6_spread3d, LAGRANGIAN_IB_6_SPREAD3D)
#endif

extern "C" {
void LAGRANGIAN_PIECEWISE_CONSTANT_INTERP_FC(const double*,
                                             const double*,
                                             const double*,
                                             const int&,
#if (NDIM == 2)
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
#endif
#if (NDIM == 3)
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
#endif
                                             const double*,
                                             const int*,
                                             const double*,
                                             const int&,
                                             const double*,
                                             double*);

void LAGRANGIAN_PIECEWISE_CONSTANT_SPREAD_FC(const double*,
                                             const double*,
                                             const double*,
                                             const int&,
                                             const int*,
                                             const double*,
                                             const int&,
                                             const double*,
                                             const double*,
#if (NDIM == 2)
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
#endif
#if (NDIM == 3)
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
#endif
                                             double*);

void LAGRANGIAN_DISCONTINUOUS_LINEAR_INTERP_FC(const double*,
                                               const double*,
                                               const double*,
                                               const int&,
                                               const int&,
#if (NDIM == 2)
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
#endif
#if (NDIM == 3)
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
#endif
                                               const double*,
                                               const int*,
                                               const double*,
                                               const int&,
                                               const double*,
                                               double*);

void LAGRANGIAN_DISCONTINUOUS_LINEAR_SPREAD_FC(const double*,
                                               const double*,
                                               const double*,
                                               const int&,
                                               const int&,
                                               const int*,
                                               const double*,
                                               const int&,
                                               const double*,
                                               const double*,
#if (NDIM == 2)
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
#endif
#if (NDIM == 3)
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
#endif
                                               double*);

void LAGRANGIAN_PIECEWISE_LINEAR_INTERP_FC(const double*,
                                           const double*,
                                           const double*,
                                           const int&,
#if (NDIM == 2)
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
#endif
#if (NDIM == 3)
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
#endif
                                           const double*,
                                           const int*,
                                           const double*,
                                           const int&,
                                           const double*,
                                           double*);

void LAGRANGIAN_PIECEWISE_LINEAR_SPREAD_FC(const double*,
                                           const double*,
                                           const double*,
                                           const int&,
                                           const int*,
                                           const double*,
                                           const int&,
                                           const double*,
                                           const double*,
#if (NDIM == 2)
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
#endif
#if (NDIM == 3)
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
                                           const int&,
#endif
                                           double*);

void LAGRANGIAN_PIECEWISE_CUBIC_INTERP_FC(const double*,
                                          const double*,
                                          const double*,
                                          const int&,
#if (NDIM == 2)
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
#endif
#if (NDIM == 3)
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
#endif
                                          const double*,
                                          const int*,
                                          const double*,
                                          const int&,
                                          const double*,
                                          double*);

void LAGRANGIAN_PIECEWISE_CUBIC_SPREAD_FC(const double*,
                                          const double*,
                                          const double*,
                                          const int&,
                                          const int*,
                                          const double*,
                                          const int&,
                                          const double*,
                                          const double*,
#if (NDIM == 2)
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
#endif
#if (NDIM == 3)
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
#endif
                                          double*);

void LAGRANGIAN_IB_3_INTERP_FC(const double*,
                               const double*,
                               const double*,
                               const int&,
#if (NDIM == 2)
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
#endif
#if (NDIM == 3)
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
#endif
                               const double*,
                               const int*,
                               const double*,
                               const int&,
                               const double*,
                               double*);

void LAGRANGIAN_IB_3_SPREAD_FC(const double*,
                               const double*,
                               const double*,
                               const int&,
                               const int*,
                               const double*,
                               const int&,
                               const double*,
                               const double*,
#if (NDIM == 2)
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
#endif
#if (NDIM == 3)
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
#endif
                               double*);

void LAGRANGIAN_IB_4_INTERP_FC(const double*,
                               const double*,
                               const double*,
                               const int&,
#if (NDIM == 2)
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
#endif
#if (NDIM == 3)
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
#endif
                               const double*,
                               const int*,
                               const double*,
                               const int&,
                               const double*,
                               double*);

void LAGRANGIAN_IB_4_SPREAD_FC(const double*,
                               const double*,
                               const double*,
                               const int&,
                               const int*,
                               const double*,
                               const int&,
                               const double*,
                               const double*,
#if (NDIM == 2)
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
#endif
#if (NDIM == 3)
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
#endif
                               double*);

void LAGRANGIAN_IB_4_W8_INTERP_FC(const double*,
                                  const double*,
                                  const double*,
                                  const int&,
#if (NDIM == 2)
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
#endif
#if (NDIM == 3)
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
#endif
                                  const double*,
                                  const int*,
                                  const double*,
                                  const int&,
                                  const double*,
                                  double*);

void LAGRANGIAN_IB_4_W8_SPREAD_FC(const double*,
                                  const double*,
                                  const double*,
                                  const int&,
                                  const int*,
                                  const double*,
                                  const int&,
                                  const double*,
                                  const double*,
#if (NDIM == 2)
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
#endif
#if (NDIM == 3)
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
#endif
                                  double*);

void LAGRANGIAN_IB_6_INTERP_FC(const double*,
                               const double*,
                               const double*,
                               const int&,
#if (NDIM == 2)
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
#endif
#if (NDIM == 3)
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
#endif
                               const double*,
                               const int*,
                               const double*,
                               const int&,
                               const double*,
                               double*);

void LAGRANGIAN_IB_6_SPREAD_FC(const double*,
                               const double*,
                               const double*,
                               const int&,
                               const int*,
                               const double*,
                               const int&,
                               const double*,
                               const double*,
#if (NDIM == 2)
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
#endif
#if (NDIM == 3)
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
#endif
                               double*);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline double
ib4_kernel_fcn(double r)
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
}

inline int
NINT(double a)
{
    return (a >= 0.0 ? static_cast<int>(a + 0.5) : static_cast<int>(a - 0.5));
}

typedef boost::multi_array<double, 1> Weight;
typedef boost::array<Weight, NDIM> TensorProductWeights;
typedef boost::multi_array<double, NDIM> MLSWeight;

void
perform_mls(const int stencil_sz,
            const double* const X,
            const int* const stencil_lower,
            const int* const /*stencil_upper*/,
            const double* const p_start,
            const double* const dx,
            const ArrayData<NDIM, double>& mask_data,
            const TensorProductWeights& D,
            MLSWeight& Psi)
{
    MLSWeight::extent_gen extents;
    MLSWeight T;

#if (NDIM == 2)
    T.resize(extents[stencil_sz][stencil_sz]);
    Psi.resize(extents[stencil_sz][stencil_sz]);
#elif(NDIM == 3)
    T.resize(extents[stencil_sz][stencil_sz][stencil_sz]);
    Psi.resize(extents[stencil_sz][stencil_sz][stencil_sz]);
#endif

    // Compute the tensor product of the weights.
    double x[NDIM], p_j, p_k;
#if (NDIM == 3)
    for (int i2 = 0; i2 < stencil_sz; ++i2)
    {
        const int ic2 = stencil_lower[2] + i2;
#endif
        for (int i1 = 0; i1 < stencil_sz; ++i1)
        {
            const int ic1 = stencil_lower[1] + i1;
            for (int i0 = 0; i0 < stencil_sz; ++i0)
            {
                const int ic0 = stencil_lower[0] + i0;
#if (NDIM == 2)
                const Index<NDIM> idx(ic0, ic1);
                T[i1][i0] = D[0][i0] * D[1][i1] * mask_data(idx, /*depth*/ 0);
#elif(NDIM == 3)
            const Index<NDIM> idx(ic0, ic1, ic2);
            T[i2][i1][i0] = D[0][i0] * D[1][i1] * D[2][i2] * mask_data(idx, /*depth*/ 0);
#endif
            }
        }
#if (NDIM == 3)
    }
#endif

    // Set the Gram matrix and the RHS.
    // Here we are solving the equation of the type G L = p, in which p
    // is the vector of basis functions that we want to reproduce, G is Gram
    // matrix and L is Lagrange muliplier which imposes the reproducibilty constraint.
    Eigen::Matrix<double, NDIM + 1, NDIM + 1> G;
    G.setZero();
    Eigen::Matrix<double, NDIM + 1, 1> p, L;
    p[0] = 1.0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        p[d + 1] = X[d];
    }

    for (int j = 0; j <= NDIM; ++j)
    {
        for (int k = 0; k <= NDIM; ++k)
        {
#if (NDIM == 3)
            for (int i2 = 0; i2 < stencil_sz; ++i2)
            {
                x[2] = p_start[2] + i2 * dx[2];
#endif
                for (int i1 = 0; i1 < stencil_sz; ++i1)
                {
                    x[1] = p_start[1] + i1 * dx[1];
                    for (int i0 = 0; i0 < stencil_sz; ++i0)
                    {
                        x[0] = p_start[0] + i0 * dx[0];

#if (NDIM == 2)
                        p_j = j == 0 ? 1.0 : (j == 1 ? x[0] : x[1]);
                        p_k = k == 0 ? 1.0 : (k == 1 ? x[0] : x[1]);
                        G(j, k) += p_j * p_k * T[i1][i0];
#elif(NDIM == 3)
                    p_j = j == 0 ? 1.0 : (j == 1 ? x[0] : j == 2 ? x[1] : x[2]);
                    p_k = k == 0 ? 1.0 : (k == 1 ? x[0] : k == 2 ? x[1] : x[2]);
                    G(j, k) += p_j * p_k * T[i2][i1][i0];
#endif
                    }
                }
#if (NDIM == 3)
            }
#endif
        }
    }

    // Solve the system for L
    L = G.ldlt().solve(p);

    // Find the modified weights using the Lagrange multiplier and to-be-reproduced
    // polynomial basis.
    std::fill(Psi.origin(), Psi.origin() + Psi.num_elements(), 0.0);
#if (NDIM == 3)
    for (int i2 = 0; i2 < stencil_sz; ++i2)
    {
        x[2] = p_start[2] + i2 * dx[2];
#endif
        for (int i1 = 0; i1 < stencil_sz; ++i1)
        {
            x[1] = p_start[1] + i1 * dx[1];
            for (int i0 = 0; i0 < stencil_sz; ++i0)
            {
                x[0] = p_start[0] + i0 * dx[0];
#if (NDIM == 2)
                for (int j = 0; j <= 2; ++j)
                {
                    p_j = j == 0 ? 1.0 : (j == 1 ? x[0] : x[1]);
                    Psi[i1][i0] += L[j] * p_j;
                }
                Psi[i1][i0] *= T[i1][i0];
#elif(NDIM == 3)
            for (int j = 0; j <= 3; ++j)
            {
                p_j = j == 0 ? 1.0 : (j == 1 ? x[0] : j == 2 ? x[1] : x[2]);
                Psi[i2][i1][i0] += L[j] * p_j;
            }
            Psi[i2][i1][i0] *= T[i2][i1][i0];
#endif
            }
        }
#if (NDIM == 3)
    }
#endif

    return;

} // perform_mls

void
get_mls_weights(const std::string& kernel_fcn,
                const double* const X,
                const double* const X_shift,
                const double* const dx,
                const double* const x_lower,
                const int* const ilower,
                const ArrayData<NDIM, double>& mask_data,
                int* stencil_lower,
                int* stencil_upper,
                MLSWeight& Psi)
{
    Weight::extent_gen extents;

    if (kernel_fcn == "IB_4")
    {
        // Resize some arrays.
        const int stencil_sz = LEInteractor::getStencilSize("IB_4");
        TensorProductWeights D;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            D[d].resize(extents[stencil_sz]);
        }

        // Determine the interpolation stencil corresponding to the position
        // of X within the cell and compute the regular IB weights.
        double X_dx, q, r, p_start[NDIM];
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dx = (X[d] + X_shift[d] - x_lower[d]) / dx[d];
            stencil_lower[d] = NINT(X_dx) + ilower[d] - 2;
            stencil_upper[d] = stencil_lower[d] + 3;
            r = X_dx - ((stencil_lower[d] + 1 - ilower[d]) + 0.5);
            p_start[d] = X[d] - (r + 1) * dx[d];
            q = std::sqrt(1.0 + 4.0 * r * (1.0 - r));
            D[d][0] = 0.125 * (3.0 - 2.0 * r - q);
            D[d][1] = 0.125 * (3.0 - 2.0 * r + q);
            D[d][2] = 0.125 * (1.0 + 2.0 * r + q);
            D[d][3] = 0.125 * (1.0 + 2.0 * r - q);
        }
        perform_mls(stencil_sz, X, stencil_lower, stencil_upper, p_start, dx, mask_data, D, Psi);
    }
    else if (kernel_fcn == "USER_DEFINED")
    {
        boost::array<double, NDIM> X_cell;
        boost::array<int, NDIM> stencil_center;

        // Determine the Cartesian cell in which X is located.
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            stencil_center[d] = static_cast<int>(std::floor((X[d] + X_shift[d] - x_lower[d]) / dx[d])) + ilower[d];
            X_cell[d] = x_lower[d] + (static_cast<double>(stencil_center[d] - ilower[d]) + 0.5) * dx[d];
        }

        // Determine the interpolation stencil corresponding to the position of
        // X within the cell.
        if (LEInteractor::s_kernel_fcn_stencil_size % 2 == 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (X[d] < X_cell[d])
                {
                    stencil_lower[d] = stencil_center[d] - LEInteractor::s_kernel_fcn_stencil_size / 2;
                    stencil_upper[d] = stencil_center[d] + LEInteractor::s_kernel_fcn_stencil_size / 2 - 1;
                }
                else
                {
                    stencil_lower[d] = stencil_center[d] - LEInteractor::s_kernel_fcn_stencil_size / 2 + 1;
                    stencil_upper[d] = stencil_center[d] + LEInteractor::s_kernel_fcn_stencil_size / 2;
                }
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                stencil_lower[d] = stencil_center[d] - LEInteractor::s_kernel_fcn_stencil_size / 2;
                stencil_upper[d] = stencil_center[d] + LEInteractor::s_kernel_fcn_stencil_size / 2;
            }
        }

        // Compute the kernel function weights.
        TensorProductWeights D;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            D[d].resize(extents[LEInteractor::s_kernel_fcn_stencil_size]);
        }

        double p_start[NDIM];
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            p_start[d] = X_cell[d] + static_cast<double>(stencil_lower[d] - stencil_center[d]) * dx[d];
            for (int k = 0, j = stencil_lower[d]; j <= stencil_upper[d]; ++j, ++k)
            {
                D[d][k] = LEInteractor::s_kernel_fcn(
                    (X[d] + X_shift[d] - (X_cell[d] + static_cast<double>(j - stencil_center[d]) * dx[d])) / dx[d]);
            }
        }
        perform_mls(
            LEInteractor::s_kernel_fcn_stencil_size, X, stencil_lower, stencil_upper, p_start, dx, mask_data, D, Psi);
    }

    return;
} // get_mls_weights

void
interpolate_data(const int stencil_sz,
                 const int* const ig_lower,
                 const int* const ig_upper,
                 const int* const stencil_lower,
                 const int* const stencil_upper,
                 const ArrayData<NDIM, double>& q_data,
                 const int q_comp,
                 const MLSWeight& Psi,
                 double& Q)
{
    int istart[NDIM], istop[NDIM];
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        istart[d] = std::max(ig_lower[d] - stencil_lower[d], 0);
        istop[d] = (stencil_sz - 1) - std::max(stencil_upper[d] - ig_upper[d], 0);
    }

    // Interpolate q onto Q using the modified weights.
    Q = 0.0;
#if (NDIM == 3)
    for (int i2 = istart[2]; i2 <= istop[2]; ++i2)
    {
        const int ic2 = stencil_lower[2] + i2;
#endif
        for (int i1 = istart[1]; i1 <= istop[1]; ++i1)
        {
            const int ic1 = stencil_lower[1] + i1;
            for (int i0 = istart[0]; i0 <= istop[0]; ++i0)
            {
                const int ic0 = stencil_lower[0] + i0;
#if (NDIM == 2)
                const Index<NDIM> idx(ic0, ic1);
                Q += q_data(idx, q_comp) * Psi[i1][i0];
#elif(NDIM == 3)
            const Index<NDIM> idx(ic0, ic1, ic2);
            Q += q_data(idx, q_comp) * Psi[i2][i1][i0];
#endif
            }
        }
#if (NDIM == 3)
    }
#endif

    return;

} // interpolate_data

void
spread_data(const int stencil_sz,
            const int* const ig_lower,
            const int* const ig_upper,
            const int* const stencil_lower,
            const int* const stencil_upper,
            const double* const dx,
            ArrayData<NDIM, double>& q_data,
            const int q_comp,
            const MLSWeight& Psi,
            const double& Q)
{
    int istart[NDIM], istop[NDIM];
    double fac = 1.0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        istart[d] = std::max(ig_lower[d] - stencil_lower[d], 0);
        istop[d] = (stencil_sz - 1) - std::max(stencil_upper[d] - ig_upper[d], 0);
        fac /= dx[d];
    }

#if (NDIM == 3)
    for (int i2 = istart[2]; i2 <= istop[2]; ++i2)
    {
        const int ic2 = stencil_lower[2] + i2;
#endif
        for (int i1 = istart[1]; i1 <= istop[1]; ++i1)
        {
            const int ic1 = stencil_lower[1] + i1;
            for (int i0 = istart[0]; i0 <= istop[0]; ++i0)
            {
                const int ic0 = stencil_lower[0] + i0;
#if (NDIM == 2)
                const Index<NDIM> idx(ic0, ic1);
                q_data(idx, q_comp) += Q * Psi[i1][i0] * fac;
#elif(NDIM == 3)
            const Index<NDIM> idx(ic0, ic1, ic2);
            q_data(idx, q_comp) += Q * Psi[i2][i1][i0] * fac;
#endif
            }
        }
#if (NDIM == 3)
    }
#endif
} // spread_data
}

double (*LEInteractor::s_kernel_fcn)(double r) = &ib4_kernel_fcn;
int LEInteractor::s_kernel_fcn_stencil_size = 4;

void LEInteractor::setFromDatabase(Pointer<Database> /*db*/)
{
    // intentionally blank
    return;
}

void
LEInteractor::printClassData(std::ostream& os)
{
    os << "LEInteractor::printClassData():\n";
    return;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

int
LEInteractor::getStencilSize(const std::string& kernel_fcn)
{
    if (kernel_fcn == "PIECEWISE_CONSTANT") return 1;
    if (kernel_fcn == "DISCONTINUOUS_LINEAR") return 2;
    if (kernel_fcn == "PIECEWISE_LINEAR") return 2;
    if (kernel_fcn == "PIECEWISE_CUBIC") return 4;
    if (kernel_fcn == "IB_3") return 4;
    if (kernel_fcn == "IB_4") return 4;
    if (kernel_fcn == "IB_4_W8") return 8;
    if (kernel_fcn == "IB_6") return 6;
    if (kernel_fcn == "USER_DEFINED") return s_kernel_fcn_stencil_size;
    TBOX_ERROR("LEInteractor::getStencilSize()\n"
               << "  Unknown kernel function "
               << kernel_fcn
               << std::endl);
    return -1;
}

int
LEInteractor::getMinimumGhostWidth(const std::string& kernel_fcn)
{
    return static_cast<int>(floor(0.5 * getStencilSize(kernel_fcn))) + 1;
}

template <class T>
void
LEInteractor::interpolate(Pointer<LData> Q_data,
                          const Pointer<LData> X_data,
                          const Pointer<LIndexSetData<T> > idx_data,
                          const Pointer<CellData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const IntVector<NDIM>& periodic_shift,
                          const std::string& interp_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_data);
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(X_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_data->getDepth() == static_cast<unsigned int>(q_data->getDepth()));
    TBOX_ASSERT(X_data->getDepth() == NDIM);
#endif
    interpolate(Q_data->getGhostedLocalFormVecArray()->data(),
                Q_data->getDepth(),
                X_data->getGhostedLocalFormVecArray()->data(),
                X_data->getDepth(),
                idx_data,
                q_data,
                patch,
                interp_box,
                periodic_shift,
                interp_fcn);
    Q_data->restoreArrays();
    X_data->restoreArrays();
    return;
}

template <class T>
void
LEInteractor::interpolate(Pointer<LData> Q_data,
                          const Pointer<LData> X_data,
                          const Pointer<LIndexSetData<T> > idx_data,
                          const Pointer<NodeData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const IntVector<NDIM>& periodic_shift,
                          const std::string& interp_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_data);
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(X_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_data->getDepth() == static_cast<unsigned int>(q_data->getDepth()));
    TBOX_ASSERT(X_data->getDepth() == NDIM);
#endif
    interpolate(Q_data->getGhostedLocalFormVecArray()->data(),
                Q_data->getDepth(),
                X_data->getGhostedLocalFormVecArray()->data(),
                X_data->getDepth(),
                idx_data,
                q_data,
                patch,
                interp_box,
                periodic_shift,
                interp_fcn);
    Q_data->restoreArrays();
    X_data->restoreArrays();
    return;
}

template <class T>
void
LEInteractor::interpolate(Pointer<LData> Q_data,
                          const Pointer<LData> X_data,
                          const Pointer<LIndexSetData<T> > idx_data,
                          const Pointer<SideData<NDIM, double> > q_data,
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
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_data);
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(X_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_data->getDepth() == NDIM);
    TBOX_ASSERT(X_data->getDepth() == NDIM);
    TBOX_ASSERT(q_data->getDepth() == 1);
#endif
    interpolate(Q_data->getGhostedLocalFormVecArray()->data(),
                Q_data->getDepth(),
                X_data->getGhostedLocalFormVecArray()->data(),
                X_data->getDepth(),
                idx_data,
                q_data,
                patch,
                interp_box,
                periodic_shift,
                interp_fcn);
    Q_data->restoreArrays();
    X_data->restoreArrays();
    return;
}

template <class T>
void
LEInteractor::interpolate(Pointer<LData> Q_data,
                          const Pointer<LData> X_data,
                          const Pointer<LIndexSetData<T> > idx_data,
                          const Pointer<EdgeData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const IntVector<NDIM>& periodic_shift,
                          const std::string& interp_fcn)
{
    if (NDIM != 3 || Q_data->getDepth() != NDIM || q_data->getDepth() != 1)
    {
        TBOX_ERROR("LEInteractor::interpolate():\n"
                   << "  edge-centered interpolation requires 3D vector-valued data.\n");
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_data);
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(X_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_data->getDepth() == NDIM);
    TBOX_ASSERT(X_data->getDepth() == NDIM);
    TBOX_ASSERT(q_data->getDepth() == 1);
#endif
    interpolate(Q_data->getGhostedLocalFormVecArray()->data(),
                Q_data->getDepth(),
                X_data->getGhostedLocalFormVecArray()->data(),
                X_data->getDepth(),
                idx_data,
                q_data,
                patch,
                interp_box,
                periodic_shift,
                interp_fcn);
    Q_data->restoreArrays();
    X_data->restoreArrays();
    return;
}

template <class T>
void
LEInteractor::interpolate(double* const Q_data,
                          const int Q_depth,
                          const double* const X_data,
                          const int X_depth,
                          const Pointer<LIndexSetData<T> > idx_data,
                          const Pointer<CellData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const IntVector<NDIM>& periodic_shift,
                          const std::string& interp_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
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
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    std::vector<double> periodic_shifts;
    buildLocalIndices(local_indices, periodic_shifts, interp_box, patch, periodic_shift, idx_data);

    // Interpolate.
    if (!local_indices.empty())
    {
        interpolate(Q_data,
                    Q_depth,
                    X_data,
                    q_data->getPointer(),
                    q_data->getBox(),
                    q_data->getGhostCellWidth(),
                    q_data->getDepth(),
                    x_lower,
                    x_upper,
                    dx,
                    patch_touches_lower_physical_bdry,
                    patch_touches_upper_physical_bdry,
                    local_indices,
                    periodic_shifts,
                    interp_fcn);
    }
    return;
}

template <class T>
void
LEInteractor::interpolate(double* const Q_data,
                          const int Q_depth,
                          const double* const X_data,
                          const int X_depth,
                          const Pointer<LIndexSetData<T> > idx_data,
                          const Pointer<NodeData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const IntVector<NDIM>& periodic_shift,
                          const std::string& interp_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
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
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    std::vector<double> periodic_shifts;
    buildLocalIndices(local_indices, periodic_shifts, interp_box, patch, periodic_shift, idx_data);

    // Interpolate.
    if (!local_indices.empty())
    {
        boost::array<double, NDIM> x_lower_node, x_upper_node;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            x_lower_node[d] = x_lower[d] - 0.5 * dx[d];
            x_upper_node[d] = x_upper[d] + 0.5 * dx[d];
        }
        interpolate(Q_data,
                    Q_depth,
                    X_data,
                    q_data->getPointer(),
                    NodeGeometry<NDIM>::toNodeBox(q_data->getBox()),
                    q_data->getGhostCellWidth(),
                    q_data->getDepth(),
                    x_lower_node.data(),
                    x_upper_node.data(),
                    dx,
                    patch_touches_lower_physical_bdry,
                    patch_touches_upper_physical_bdry,
                    local_indices,
                    periodic_shifts,
                    interp_fcn);
    }
    return;
}

template <class T>
void
LEInteractor::interpolate(double* const Q_data,
                          const int Q_depth,
                          const double* const X_data,
                          const int X_depth,
                          const Pointer<LIndexSetData<T> > idx_data,
                          const Pointer<SideData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const IntVector<NDIM>& periodic_shift,
                          const std::string& interp_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
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
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    std::vector<double> periodic_shifts;
    buildLocalIndices(local_indices, periodic_shifts, interp_box, patch, periodic_shift, idx_data);

    // Interpolate.
    if (!local_indices.empty())
    {
        boost::array<double, NDIM> x_lower_axis, x_upper_axis;
        const int local_sz = (*std::max_element(local_indices.begin(), local_indices.end())) + 1;
        std::vector<double> Q_data_axis(local_sz);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                x_lower_axis[d] = x_lower[d];
                x_upper_axis[d] = x_upper[d];
            }
            x_lower_axis[axis] -= 0.5 * dx[axis];
            x_upper_axis[axis] += 0.5 * dx[axis];
            interpolate(&Q_data_axis[0],
                        /*Q_depth*/ 1,
                        X_data,
                        q_data->getPointer(axis),
                        SideGeometry<NDIM>::toSideBox(q_data->getBox(), axis),
                        q_data->getGhostCellWidth(),
                        /*q_depth*/ 1,
                        x_lower_axis.data(),
                        x_upper_axis.data(),
                        dx,
                        patch_touches_lower_physical_bdry,
                        patch_touches_upper_physical_bdry,
                        local_indices,
                        periodic_shifts,
                        interp_fcn,
                        axis);
            for (unsigned int k = 0; k < local_indices.size(); ++k)
            {
                Q_data[NDIM * local_indices[k] + axis] = Q_data_axis[local_indices[k]];
            }
        }
    }
    return;
}

template <class T>
void
LEInteractor::interpolate(double* const Q_data,
                          const int Q_depth,
                          const double* const X_data,
                          const int X_depth,
                          const Pointer<LIndexSetData<T> > idx_data,
                          const Pointer<EdgeData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const IntVector<NDIM>& periodic_shift,
                          const std::string& interp_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_depth == NDIM);
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(q_data->getDepth() == 1);
#else
    NULL_USE(X_depth);
#endif
    if (NDIM != 3 || Q_depth != NDIM || q_data->getDepth() != 1)
    {
        TBOX_ERROR("LEInteractor::interpolate():\n"
                   << "  edge-centered interpolation requires 3D vector-valued data.\n");
    }

    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    std::vector<double> periodic_shifts;
    buildLocalIndices(local_indices, periodic_shifts, interp_box, patch, periodic_shift, idx_data);

    // Interpolate.
    if (!local_indices.empty())
    {
        boost::array<double, NDIM> x_lower_axis, x_upper_axis;
        const int local_sz = (*std::max_element(local_indices.begin(), local_indices.end())) + 1;
        std::vector<double> Q_data_axis(local_sz);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                x_lower_axis[d] = x_lower[d];
                x_upper_axis[d] = x_upper[d];
                if (d != axis)
                {
                    x_lower_axis[d] -= 0.5 * dx[d];
                    x_upper_axis[d] += 0.5 * dx[d];
                }
            }
            interpolate(&Q_data_axis[0],
                        /*Q_depth*/ 1,
                        X_data,
                        q_data->getPointer(axis),
                        EdgeGeometry<NDIM>::toEdgeBox(q_data->getBox(), axis),
                        q_data->getGhostCellWidth(),
                        /*q_depth*/ 1,
                        x_lower_axis.data(),
                        x_upper_axis.data(),
                        dx,
                        patch_touches_lower_physical_bdry,
                        patch_touches_upper_physical_bdry,
                        local_indices,
                        periodic_shifts,
                        interp_fcn,
                        axis);
            for (unsigned int k = 0; k < local_indices.size(); ++k)
            {
                Q_data[NDIM * local_indices[k] + axis] = Q_data_axis[local_indices[k]];
            }
        }
    }
    return;
}

void
LEInteractor::interpolate(std::vector<double>& Q_data,
                          const int Q_depth,
                          const std::vector<double>& X_data,
                          const int X_depth,
                          const Pointer<CellData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const std::string& interp_fcn)
{
    if (Q_data.empty()) return;
    interpolate(&Q_data[0],
                static_cast<int>(Q_data.size()),
                Q_depth,
                &X_data[0],
                static_cast<int>(X_data.size()),
                X_depth,
                q_data,
                patch,
                interp_box,
                interp_fcn);
}

void
LEInteractor::interpolate(std::vector<double>& Q_data,
                          const int Q_depth,
                          const std::vector<double>& X_data,
                          const int X_depth,
                          const Pointer<CellData<NDIM, double> > mask_data,
                          const Pointer<CellData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const std::string& interp_fcn)
{
    if (Q_data.empty()) return;
    interpolate(&Q_data[0],
                static_cast<int>(Q_data.size()),
                Q_depth,
                &X_data[0],
                static_cast<int>(X_data.size()),
                X_depth,
                mask_data,
                q_data,
                patch,
                interp_box,
                interp_fcn);
}

void
LEInteractor::interpolate(std::vector<double>& Q_data,
                          const int Q_depth,
                          const std::vector<double>& X_data,
                          const int X_depth,
                          const Pointer<NodeData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const std::string& interp_fcn)
{
    if (Q_data.empty()) return;
    interpolate(&Q_data[0],
                static_cast<int>(Q_data.size()),
                Q_depth,
                &X_data[0],
                static_cast<int>(X_data.size()),
                X_depth,
                q_data,
                patch,
                interp_box,
                interp_fcn);
}

void
LEInteractor::interpolate(std::vector<double>& Q_data,
                          const int Q_depth,
                          const std::vector<double>& X_data,
                          const int X_depth,
                          const Pointer<SideData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const std::string& interp_fcn)
{
    if (Q_data.empty()) return;
    interpolate(&Q_data[0],
                static_cast<int>(Q_data.size()),
                Q_depth,
                &X_data[0],
                static_cast<int>(X_data.size()),
                X_depth,
                q_data,
                patch,
                interp_box,
                interp_fcn);
}

void
LEInteractor::interpolate(std::vector<double>& Q_data,
                          const int Q_depth,
                          const std::vector<double>& X_data,
                          const int X_depth,
                          const Pointer<SideData<NDIM, double> > mask_data,
                          const Pointer<SideData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const std::string& interp_fcn)
{
    if (Q_data.empty()) return;
    interpolate(&Q_data[0],
                static_cast<int>(Q_data.size()),
                Q_depth,
                &X_data[0],
                static_cast<int>(X_data.size()),
                X_depth,
                mask_data,
                q_data,
                patch,
                interp_box,
                interp_fcn);
}

void
LEInteractor::interpolate(std::vector<double>& Q_data,
                          const int Q_depth,
                          const std::vector<double>& X_data,
                          const int X_depth,
                          const Pointer<EdgeData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const std::string& interp_fcn)
{
    if (Q_data.empty()) return;
    interpolate(&Q_data[0],
                static_cast<int>(Q_data.size()),
                Q_depth,
                &X_data[0],
                static_cast<int>(X_data.size()),
                X_depth,
                q_data,
                patch,
                interp_box,
                interp_fcn);
}

void
LEInteractor::interpolate(double* const Q_data,
                          const int Q_size,
                          const int Q_depth,
                          const double* const X_data,
                          const int X_size,
                          const int X_depth,
                          const Pointer<CellData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const std::string& interp_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_depth == q_data->getDepth());
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(Q_size / Q_depth == X_size / X_depth);
#else
    NULL_USE(Q_size);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, interp_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_shifts(NDIM * local_indices.size());

    // Interpolate.
    if (!local_indices.empty())
    {
        interpolate(Q_data,
                    Q_depth,
                    X_data,
                    q_data->getPointer(),
                    q_data->getBox(),
                    q_data->getGhostCellWidth(),
                    q_data->getDepth(),
                    x_lower,
                    x_upper,
                    dx,
                    patch_touches_lower_physical_bdry,
                    patch_touches_upper_physical_bdry,
                    local_indices,
                    periodic_shifts,
                    interp_fcn);
    }
    return;
}

void
LEInteractor::interpolate(double* const Q_data,
                          const int Q_size,
                          const int Q_depth,
                          const double* const X_data,
                          const int X_size,
                          const int X_depth,
                          const Pointer<CellData<NDIM, double> > mask_data,
                          const Pointer<CellData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const std::string& interp_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_depth == q_data->getDepth());
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(Q_size / Q_depth == X_size / X_depth);
    TBOX_ASSERT(mask_data);
    TBOX_ASSERT(mask_data->getDepth() == 1);
#else
    NULL_USE(Q_size);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();
    const Box<NDIM>& patch_box = patch->getBox();
    const IntVector<NDIM>& ilower = patch_box.lower();
    const IntVector<NDIM>& iupper = patch_box.upper();

    // Get ghost cell width info.
    const IntVector<NDIM>& q_gcw = q_data->getGhostCellWidth();
    const IntVector<NDIM>& mask_gcw = q_data->getGhostCellWidth();
    const int stencil_size = getStencilSize(interp_fcn);
    const int min_ghosts = getMinimumGhostWidth(interp_fcn);
    const int q_gcw_min = q_gcw.min();
    const int mask_gcw_min = mask_gcw.min();
    if (q_gcw_min < min_ghosts)
    {
        TBOX_ERROR("LEInteractor::interpolate(): insufficient ghost cells for Eulerian field data:\n"
                   << "  kernel function          = "
                   << interp_fcn
                   << "\n"
                   << "  kernel stencil size      = "
                   << stencil_size
                   << "\n"
                   << "  minimum ghost cell width = "
                   << min_ghosts
                   << "\n"
                   << "  ghost cell width         = "
                   << q_gcw_min
                   << "\n");
    }
    if (mask_gcw_min < stencil_size)
    {
        TBOX_ERROR("LEInteractor::interpolate(): insufficient ghost cells for Eulerian mask data:\n"
                   << "  kernel function          = "
                   << interp_fcn
                   << "\n"
                   << "  kernel stencil size      = "
                   << stencil_size
                   << "\n"
                   << "  minimum ghost cell width = "
                   << stencil_size
                   << "\n"
                   << "  ghost cell width         = "
                   << mask_gcw_min
                   << "\n");
    }
    const IntVector<NDIM> ig_lower = ilower - q_gcw;
    const IntVector<NDIM> ig_upper = iupper + q_gcw;

    // Get boundary info.
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, interp_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_shifts(NDIM * local_indices.size());

    // Interpolate.
    const int nindices = static_cast<int>(local_indices.size());
    if (nindices)
    {
        IntVector<NDIM> stencil_lower, stencil_upper;
        for (int k = 0; k < nindices; ++k)
        {
            int s = local_indices[k];
            MLSWeight Psi;
            const int stencil_sz = LEInteractor::getStencilSize(interp_fcn);
            get_mls_weights(interp_fcn,
                            &X_data[s * NDIM],
                            &periodic_shifts[k * NDIM],
                            dx,
                            x_lower,
                            ilower,
                            mask_data->getArrayData(),
                            stencil_lower,
                            stencil_upper,
                            Psi);

            for (int comp = 0; comp < Q_depth; ++comp)
            {
                interpolate_data(stencil_sz,
                                 ig_lower,
                                 ig_upper,
                                 stencil_lower,
                                 stencil_upper,
                                 q_data->getArrayData(),
                                 comp,
                                 Psi,
                                 Q_data[s * Q_depth + comp]);
            }
        }
    }

    return;
}

void
LEInteractor::interpolate(double* const Q_data,
                          const int Q_size,
                          const int Q_depth,
                          const double* const X_data,
                          const int X_size,
                          const int X_depth,
                          const Pointer<NodeData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const std::string& interp_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_depth == q_data->getDepth());
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(Q_size / Q_depth == X_size / X_depth);
#else
    NULL_USE(Q_size);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, interp_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_shifts(NDIM * local_indices.size());

    // Interpolate.
    if (!local_indices.empty())
    {
        boost::array<double, NDIM> x_lower_node, x_upper_node;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            x_lower_node[d] = x_lower[d] - 0.5 * dx[d];
            x_upper_node[d] = x_upper[d] + 0.5 * dx[d];
        }
        interpolate(Q_data,
                    Q_depth,
                    X_data,
                    q_data->getPointer(),
                    NodeGeometry<NDIM>::toNodeBox(q_data->getBox()),
                    q_data->getGhostCellWidth(),
                    q_data->getDepth(),
                    x_lower_node.data(),
                    x_upper_node.data(),
                    dx,
                    patch_touches_lower_physical_bdry,
                    patch_touches_upper_physical_bdry,
                    local_indices,
                    periodic_shifts,
                    interp_fcn);
    }
    return;
}

void
LEInteractor::interpolate(double* const Q_data,
                          const int Q_size,
                          const int Q_depth,
                          const double* const X_data,
                          const int X_size,
                          const int X_depth,
                          const Pointer<SideData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const std::string& interp_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_depth == NDIM);
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(Q_size / Q_depth == X_size / X_depth);
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
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, interp_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_shifts(NDIM * local_indices.size());

    // Interpolate.
    if (!local_indices.empty())
    {
        boost::array<double, NDIM> x_lower_axis, x_upper_axis;
        const int local_sz = (*std::max_element(local_indices.begin(), local_indices.end())) + 1;
        std::vector<double> Q_data_axis(local_sz);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                x_lower_axis[d] = x_lower[d];
                x_upper_axis[d] = x_upper[d];
            }
            x_lower_axis[axis] -= 0.5 * dx[axis];
            x_upper_axis[axis] += 0.5 * dx[axis];
            interpolate(&Q_data_axis[0],
                        /*Q_depth*/ 1,
                        X_data,
                        q_data->getPointer(axis),
                        SideGeometry<NDIM>::toSideBox(q_data->getBox(), axis),
                        q_data->getGhostCellWidth(),
                        /*q_depth*/ 1,
                        x_lower_axis.data(),
                        x_upper_axis.data(),
                        dx,
                        patch_touches_lower_physical_bdry,
                        patch_touches_upper_physical_bdry,
                        local_indices,
                        periodic_shifts,
                        interp_fcn,
                        axis);
            for (unsigned int k = 0; k < local_indices.size(); ++k)
            {
                Q_data[NDIM * local_indices[k] + axis] = Q_data_axis[local_indices[k]];
            }
        }
    }
    return;
}

void
LEInteractor::interpolate(double* const Q_data,
                          const int Q_size,
                          const int Q_depth,
                          const double* const X_data,
                          const int X_size,
                          const int X_depth,
                          const Pointer<SideData<NDIM, double> > mask_data,
                          const Pointer<SideData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const std::string& interp_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_depth == NDIM);
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(Q_size / Q_depth == X_size / X_depth);
    TBOX_ASSERT(q_data->getDepth() == 1);
    TBOX_ASSERT(mask_data);
    TBOX_ASSERT(mask_data->getDepth() == 1);
#else
    NULL_USE(Q_size);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    const Box<NDIM>& patch_box = patch->getBox();
    const IntVector<NDIM>& ilower = patch_box.lower();

    // Get ghost cell width info.
    const IntVector<NDIM>& q_gcw = q_data->getGhostCellWidth();
    const IntVector<NDIM>& mask_gcw = mask_data->getGhostCellWidth();
    const int stencil_size = getStencilSize(interp_fcn);
    const int min_ghosts = getMinimumGhostWidth(interp_fcn);
    const int q_gcw_min = q_gcw.min();
    const int mask_gcw_min = mask_gcw.min();
    if (q_gcw_min < min_ghosts)
    {
        TBOX_ERROR("LEInteractor::interpolate(): insufficient ghost cells for Eulerian field data:\n"
                   << "  kernel function          = "
                   << interp_fcn
                   << "\n"
                   << "  kernel stencil size      = "
                   << stencil_size
                   << "\n"
                   << "  minimum ghost cell width = "
                   << min_ghosts
                   << "\n"
                   << "  ghost cell width         = "
                   << q_gcw_min
                   << "\n");
    }
    if (mask_gcw_min < stencil_size)
    {
        TBOX_ERROR("LEInteractor::interpolate(): insufficient ghost cells for Eulerian mask data:\n"
                   << "  kernel function          = "
                   << interp_fcn
                   << "\n"
                   << "  kernel stencil size      = "
                   << stencil_size
                   << "\n"
                   << "  minimum ghost cell width = "
                   << stencil_size
                   << "\n"
                   << "  ghost cell width         = "
                   << mask_gcw_min
                   << "\n");
    }

    // Get boundary info.
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, interp_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_shifts(NDIM * local_indices.size());

    // Interpolate.
    const int nindices = static_cast<int>(local_indices.size());
    if (nindices)
    {
        boost::array<double, NDIM> x_lower_axis, x_upper_axis;
        IntVector<NDIM> stencil_lower, stencil_upper;

        for (int axis = 0; axis < NDIM; ++axis)
        {
            Box<NDIM> data_box = SideGeometry<NDIM>::toSideBox(q_data->getBox(), axis);
            const IntVector<NDIM> ig_lower = data_box.lower() - q_gcw;
            const IntVector<NDIM> ig_upper = data_box.upper() + q_gcw;

            for (int d = 0; d < NDIM; ++d)
            {
                x_lower_axis[d] = x_lower[d];
                x_upper_axis[d] = x_upper[d];
            }
            x_lower_axis[axis] -= 0.5 * dx[axis];
            x_upper_axis[axis] += 0.5 * dx[axis];

            for (int k = 0; k < nindices; ++k)
            {
                int s = local_indices[k];
                MLSWeight Psi;
                const int stencil_sz = LEInteractor::getStencilSize(interp_fcn);
                get_mls_weights(interp_fcn,
                                &X_data[s * NDIM],
                                &periodic_shifts[k * NDIM],
                                dx,
                                x_lower_axis.data(),
                                ilower,
                                mask_data->getArrayData(axis),
                                stencil_lower,
                                stencil_upper,
                                Psi);
                interpolate_data(stencil_sz,
                                 ig_lower,
                                 ig_upper,
                                 stencil_lower,
                                 stencil_upper,
                                 q_data->getArrayData(axis),
                                 0,
                                 Psi,
                                 Q_data[s * Q_depth + axis]);
            }
        }
    }

    return;
}

void
LEInteractor::interpolate(double* const Q_data,
                          const int Q_size,
                          const int Q_depth,
                          const double* const X_data,
                          const int X_size,
                          const int X_depth,
                          const Pointer<EdgeData<NDIM, double> > q_data,
                          const Pointer<Patch<NDIM> > patch,
                          const Box<NDIM>& interp_box,
                          const std::string& interp_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_depth == NDIM);
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(Q_size / Q_depth == X_size / X_depth);
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
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, interp_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_shifts(NDIM * local_indices.size());

    // Interpolate.
    if (!local_indices.empty())
    {
        boost::array<double, NDIM> x_lower_axis, x_upper_axis;
        const int local_sz = (*std::max_element(local_indices.begin(), local_indices.end())) + 1;
        std::vector<double> Q_data_axis(local_sz);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                x_lower_axis[d] = x_lower[d];
                x_upper_axis[d] = x_upper[d];
                if (d != axis)
                {
                    x_lower_axis[d] -= 0.5 * dx[d];
                    x_upper_axis[d] += 0.5 * dx[d];
                }
            }
            interpolate(&Q_data_axis[0],
                        /*Q_depth*/ 1,
                        X_data,
                        q_data->getPointer(axis),
                        EdgeGeometry<NDIM>::toEdgeBox(q_data->getBox(), axis),
                        q_data->getGhostCellWidth(),
                        /*q_depth*/ 1,
                        x_lower_axis.data(),
                        x_upper_axis.data(),
                        dx,
                        patch_touches_lower_physical_bdry,
                        patch_touches_upper_physical_bdry,
                        local_indices,
                        periodic_shifts,
                        interp_fcn,
                        axis);
            for (unsigned int k = 0; k < local_indices.size(); ++k)
            {
                Q_data[NDIM * local_indices[k] + axis] = Q_data_axis[local_indices[k]];
            }
        }
    }
    return;
}

template <class T>
void
LEInteractor::spread(Pointer<CellData<NDIM, double> > q_data,
                     const Pointer<LData> Q_data,
                     const Pointer<LData> X_data,
                     const Pointer<LIndexSetData<T> > idx_data,
                     const Pointer<Patch<NDIM> > patch,
                     const Box<NDIM>& spread_box,
                     const IntVector<NDIM>& periodic_shift,
                     const std::string& spread_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_data);
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(X_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_data->getDepth() == static_cast<unsigned int>(q_data->getDepth()));
    TBOX_ASSERT(X_data->getDepth() == NDIM);
#endif
    spread(q_data,
           Q_data->getGhostedLocalFormVecArray()->data(),
           Q_data->getDepth(),
           X_data->getGhostedLocalFormVecArray()->data(),
           X_data->getDepth(),
           idx_data,
           patch,
           spread_box,
           periodic_shift,
           spread_fcn);
    Q_data->restoreArrays();
    X_data->restoreArrays();
    return;
}

template <class T>
void
LEInteractor::spread(Pointer<NodeData<NDIM, double> > q_data,
                     const Pointer<LData> Q_data,
                     const Pointer<LData> X_data,
                     const Pointer<LIndexSetData<T> > idx_data,
                     const Pointer<Patch<NDIM> > patch,
                     const Box<NDIM>& spread_box,
                     const IntVector<NDIM>& periodic_shift,
                     const std::string& spread_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_data);
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(X_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_data->getDepth() == static_cast<unsigned int>(q_data->getDepth()));
    TBOX_ASSERT(X_data->getDepth() == NDIM);
#endif
    spread(q_data,
           Q_data->getGhostedLocalFormVecArray()->data(),
           Q_data->getDepth(),
           X_data->getGhostedLocalFormVecArray()->data(),
           X_data->getDepth(),
           idx_data,
           patch,
           spread_box,
           periodic_shift,
           spread_fcn);
    Q_data->restoreArrays();
    X_data->restoreArrays();
    return;
}

template <class T>
void
LEInteractor::spread(Pointer<SideData<NDIM, double> > q_data,
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
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_data);
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(X_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(q_data->getDepth() == 1);
    TBOX_ASSERT(Q_data->getDepth() == NDIM);
    TBOX_ASSERT(X_data->getDepth() == NDIM);
#endif
    spread(q_data,
           Q_data->getGhostedLocalFormVecArray()->data(),
           Q_data->getDepth(),
           X_data->getGhostedLocalFormVecArray()->data(),
           X_data->getDepth(),
           idx_data,
           patch,
           spread_box,
           periodic_shift,
           spread_fcn);
    Q_data->restoreArrays();
    X_data->restoreArrays();
    return;
}

template <class T>
void
LEInteractor::spread(Pointer<EdgeData<NDIM, double> > q_data,
                     const Pointer<LData> Q_data,
                     const Pointer<LData> X_data,
                     const Pointer<LIndexSetData<T> > idx_data,
                     const Pointer<Patch<NDIM> > patch,
                     const Box<NDIM>& spread_box,
                     const IntVector<NDIM>& periodic_shift,
                     const std::string& spread_fcn)
{
    if (NDIM != 3 || Q_data->getDepth() != NDIM || q_data->getDepth() != 1)
    {
        TBOX_ERROR("LEInteractor::spread():\n"
                   << "  edge-centered interpolation requires 3D vector-valued data.\n");
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_data);
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(X_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(q_data->getDepth() == 1);
    TBOX_ASSERT(Q_data->getDepth() == NDIM);
    TBOX_ASSERT(X_data->getDepth() == NDIM);
#endif
    spread(q_data,
           Q_data->getGhostedLocalFormVecArray()->data(),
           Q_data->getDepth(),
           X_data->getGhostedLocalFormVecArray()->data(),
           X_data->getDepth(),
           idx_data,
           patch,
           spread_box,
           periodic_shift,
           spread_fcn);
    Q_data->restoreArrays();
    X_data->restoreArrays();
    return;
}

template <class T>
void
LEInteractor::spread(Pointer<CellData<NDIM, double> > q_data,
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
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
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
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    std::vector<double> periodic_shifts;
    buildLocalIndices(local_indices, periodic_shifts, spread_box, patch, periodic_shift, idx_data);

    // Spread.
    if (!local_indices.empty())
    {
        spread(q_data->getPointer(),
               q_data->getBox(),
               q_data->getGhostCellWidth(),
               q_data->getDepth(),
               Q_data,
               Q_depth,
               X_data,
               x_lower,
               x_upper,
               dx,
               patch_touches_lower_physical_bdry,
               patch_touches_upper_physical_bdry,
               local_indices,
               periodic_shifts,
               spread_fcn);
    }
    return;
}

template <class T>
void
LEInteractor::spread(Pointer<NodeData<NDIM, double> > q_data,
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
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
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
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    std::vector<double> periodic_shifts;
    buildLocalIndices(local_indices, periodic_shifts, spread_box, patch, periodic_shift, idx_data);

    // Spread.
    if (!local_indices.empty())
    {
        boost::array<double, NDIM> x_lower_node, x_upper_node;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            x_lower_node[d] = x_lower[d] - 0.5 * dx[d];
            x_upper_node[d] = x_upper[d] + 0.5 * dx[d];
        }
        spread(q_data->getPointer(),
               NodeGeometry<NDIM>::toNodeBox(q_data->getBox()),
               q_data->getGhostCellWidth(),
               q_data->getDepth(),
               Q_data,
               Q_depth,
               X_data,
               x_lower_node.data(),
               x_upper_node.data(),
               dx,
               patch_touches_lower_physical_bdry,
               patch_touches_upper_physical_bdry,
               local_indices,
               periodic_shifts,
               spread_fcn);
    }
    return;
}

template <class T>
void
LEInteractor::spread(Pointer<SideData<NDIM, double> > q_data,
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
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
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
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    std::vector<double> periodic_shifts;
    buildLocalIndices(local_indices, periodic_shifts, spread_box, patch, periodic_shift, idx_data);

    // Spread.
    if (!local_indices.empty())
    {
        boost::array<double, NDIM> x_lower_axis, x_upper_axis;
        const int local_sz = (*std::max_element(local_indices.begin(), local_indices.end())) + 1;
        std::vector<double> Q_data_axis(local_sz);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                x_lower_axis[d] = x_lower[d];
                x_upper_axis[d] = x_upper[d];
            }
            x_lower_axis[axis] -= 0.5 * dx[axis];
            x_upper_axis[axis] += 0.5 * dx[axis];
            for (unsigned int k = 0; k < local_indices.size(); ++k)
            {
                Q_data_axis[local_indices[k]] = Q_data[NDIM * local_indices[k] + axis];
            }
            spread(q_data->getPointer(axis),
                   SideGeometry<NDIM>::toSideBox(q_data->getBox(), axis),
                   q_data->getGhostCellWidth(),
                   /*q_depth*/ 1,
                   &Q_data_axis[0],
                   /*Q_depth*/ 1,
                   X_data,
                   x_lower_axis.data(),
                   x_upper_axis.data(),
                   dx,
                   patch_touches_lower_physical_bdry,
                   patch_touches_upper_physical_bdry,
                   local_indices,
                   periodic_shifts,
                   spread_fcn,
                   axis);
        }
    }
    return;
}

template <class T>
void
LEInteractor::spread(Pointer<EdgeData<NDIM, double> > q_data,
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
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(idx_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(q_data->getDepth() == 1);
    TBOX_ASSERT(Q_depth == NDIM);
    TBOX_ASSERT(X_depth == NDIM);
#else
    NULL_USE(X_depth);
#endif
    if (NDIM != 3 || Q_depth != NDIM || q_data->getDepth() != 1)
    {
        TBOX_ERROR("LEInteractor::spread():\n"
                   << "  edge-centered interpolation requires 3D vector-valued data.\n");
    }

    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    std::vector<double> periodic_shifts;
    buildLocalIndices(local_indices, periodic_shifts, spread_box, patch, periodic_shift, idx_data);

    // Spread.
    if (!local_indices.empty())
    {
        boost::array<double, NDIM> x_lower_axis, x_upper_axis;
        const int local_sz = (*std::max_element(local_indices.begin(), local_indices.end())) + 1;
        std::vector<double> Q_data_axis(local_sz);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                x_lower_axis[d] = x_lower[d];
                x_upper_axis[d] = x_upper[d];
                if (d != axis)
                {
                    x_lower_axis[d] -= 0.5 * dx[d];
                    x_upper_axis[d] += 0.5 * dx[d];
                }
            }
            for (unsigned int k = 0; k < local_indices.size(); ++k)
            {
                Q_data_axis[local_indices[k]] = Q_data[NDIM * local_indices[k] + axis];
            }
            spread(q_data->getPointer(axis),
                   EdgeGeometry<NDIM>::toEdgeBox(q_data->getBox(), axis),
                   q_data->getGhostCellWidth(),
                   /*q_depth*/ 1,
                   &Q_data_axis[0],
                   /*Q_depth*/ 1,
                   X_data,
                   x_lower_axis.data(),
                   x_upper_axis.data(),
                   dx,
                   patch_touches_lower_physical_bdry,
                   patch_touches_upper_physical_bdry,
                   local_indices,
                   periodic_shifts,
                   spread_fcn,
                   axis);
        }
    }
    return;
}

void
LEInteractor::spread(Pointer<CellData<NDIM, double> > q_data,
                     const std::vector<double>& Q_data,
                     const int Q_depth,
                     const std::vector<double>& X_data,
                     const int X_depth,
                     const Pointer<Patch<NDIM> > patch,
                     const Box<NDIM>& spread_box,
                     const std::string& spread_fcn)
{
    if (Q_data.empty()) return;
    spread(q_data,
           &Q_data[0],
           static_cast<int>(Q_data.size()),
           Q_depth,
           &X_data[0],
           static_cast<int>(X_data.size()),
           X_depth,
           patch,
           spread_box,
           spread_fcn);
}

void
LEInteractor::spread(Pointer<CellData<NDIM, double> > mask_data,
                     Pointer<CellData<NDIM, double> > q_data,
                     const std::vector<double>& Q_data,
                     const int Q_depth,
                     const std::vector<double>& X_data,
                     const int X_depth,
                     const Pointer<Patch<NDIM> > patch,
                     const Box<NDIM>& spread_box,
                     const std::string& spread_fcn)
{
    if (Q_data.empty()) return;
    spread(mask_data,
           q_data,
           &Q_data[0],
           static_cast<int>(Q_data.size()),
           Q_depth,
           &X_data[0],
           static_cast<int>(X_data.size()),
           X_depth,
           patch,
           spread_box,
           spread_fcn);
}

void
LEInteractor::spread(Pointer<NodeData<NDIM, double> > q_data,
                     const std::vector<double>& Q_data,
                     const int Q_depth,
                     const std::vector<double>& X_data,
                     const int X_depth,
                     const Pointer<Patch<NDIM> > patch,
                     const Box<NDIM>& spread_box,
                     const std::string& spread_fcn)
{
    if (Q_data.empty()) return;
    spread(q_data,
           &Q_data[0],
           static_cast<int>(Q_data.size()),
           Q_depth,
           &X_data[0],
           static_cast<int>(X_data.size()),
           X_depth,
           patch,
           spread_box,
           spread_fcn);
}

void
LEInteractor::spread(Pointer<SideData<NDIM, double> > q_data,
                     const std::vector<double>& Q_data,
                     const int Q_depth,
                     const std::vector<double>& X_data,
                     const int X_depth,
                     const Pointer<Patch<NDIM> > patch,
                     const Box<NDIM>& spread_box,
                     const std::string& spread_fcn)
{
    if (Q_data.empty()) return;
    spread(q_data,
           &Q_data[0],
           static_cast<int>(Q_data.size()),
           Q_depth,
           &X_data[0],
           static_cast<int>(X_data.size()),
           X_depth,
           patch,
           spread_box,
           spread_fcn);
}

void
LEInteractor::spread(Pointer<SideData<NDIM, double> > mask_data,
                     Pointer<SideData<NDIM, double> > q_data,
                     const std::vector<double>& Q_data,
                     const int Q_depth,
                     const std::vector<double>& X_data,
                     const int X_depth,
                     const Pointer<Patch<NDIM> > patch,
                     const Box<NDIM>& spread_box,
                     const std::string& spread_fcn)
{
    if (Q_data.empty()) return;
    spread(mask_data,
           q_data,
           &Q_data[0],
           static_cast<int>(Q_data.size()),
           Q_depth,
           &X_data[0],
           static_cast<int>(X_data.size()),
           X_depth,
           patch,
           spread_box,
           spread_fcn);
}

void
LEInteractor::spread(Pointer<EdgeData<NDIM, double> > q_data,
                     const std::vector<double>& Q_data,
                     const int Q_depth,
                     const std::vector<double>& X_data,
                     const int X_depth,
                     const Pointer<Patch<NDIM> > patch,
                     const Box<NDIM>& spread_box,
                     const std::string& spread_fcn)
{
    if (Q_data.empty()) return;
    spread(q_data,
           &Q_data[0],
           static_cast<int>(Q_data.size()),
           Q_depth,
           &X_data[0],
           static_cast<int>(X_data.size()),
           X_depth,
           patch,
           spread_box,
           spread_fcn);
}

void
LEInteractor::spread(Pointer<CellData<NDIM, double> > q_data,
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
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_depth == q_data->getDepth());
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(Q_size / Q_depth == X_size / X_depth);
#else
    NULL_USE(Q_size);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, spread_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_shifts(NDIM * local_indices.size());

    // Spread.
    if (!local_indices.empty())
    {
        spread(q_data->getPointer(),
               q_data->getBox(),
               q_data->getGhostCellWidth(),
               q_data->getDepth(),
               Q_data,
               Q_depth,
               X_data,
               x_lower,
               x_upper,
               dx,
               patch_touches_lower_physical_bdry,
               patch_touches_upper_physical_bdry,
               local_indices,
               periodic_shifts,
               spread_fcn);
    }
    return;
}

void
LEInteractor::spread(Pointer<CellData<NDIM, double> > mask_data,
                     Pointer<CellData<NDIM, double> > q_data,
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
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_depth == q_data->getDepth());
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(Q_size / Q_depth == X_size / X_depth);
    TBOX_ASSERT(mask_data);
    TBOX_ASSERT(mask_data->getDepth() == 1);
#else
    NULL_USE(Q_size);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();
    const Box<NDIM>& patch_box = patch->getBox();
    const IntVector<NDIM>& ilower = patch_box.lower();
    const IntVector<NDIM>& iupper = patch_box.upper();

    // Get ghost cell width info.
    const IntVector<NDIM>& q_gcw = q_data->getGhostCellWidth();
    const IntVector<NDIM>& mask_gcw = mask_data->getGhostCellWidth();
    const int stencil_size = getStencilSize(spread_fcn);
    const int min_ghosts = getMinimumGhostWidth(spread_fcn);
    const int q_gcw_min = q_gcw.min();
    const int mask_gcw_min = mask_gcw.min();
    if (q_gcw_min < min_ghosts)
    {
        TBOX_ERROR("LEInteractor::spread(): insufficient ghost cells for Eulerian field data:\n"
                   << "  kernel function          = "
                   << spread_fcn
                   << "\n"
                   << "  kernel stencil size      = "
                   << stencil_size
                   << "\n"
                   << "  minimum ghost cell width = "
                   << min_ghosts
                   << "\n"
                   << "  ghost cell width         = "
                   << q_gcw_min
                   << "\n");
    }
    if (mask_gcw_min < stencil_size)
    {
        TBOX_ERROR("LEInteractor::spread(): insufficient ghost cells for Eulerian mask data:\n"
                   << "  kernel function          = "
                   << spread_fcn
                   << "\n"
                   << "  kernel stencil size      = "
                   << stencil_size
                   << "\n"
                   << "  minimum ghost cell width = "
                   << stencil_size
                   << "\n"
                   << "  ghost cell width         = "
                   << mask_gcw_min
                   << "\n");
    }
    const IntVector<NDIM> ig_lower = ilower - q_gcw;
    const IntVector<NDIM> ig_upper = iupper + q_gcw;

    // Get boundary info.
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, spread_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_shifts(NDIM * local_indices.size());

    // Spread.
    const int nindices = static_cast<int>(local_indices.size());
    if (nindices)
    {
        IntVector<NDIM> stencil_lower, stencil_upper;
        for (int k = 0; k < nindices; ++k)
        {
            int s = local_indices[k];
            MLSWeight Psi;
            const int stencil_sz = LEInteractor::getStencilSize(spread_fcn);
            get_mls_weights(spread_fcn,
                            &X_data[s * NDIM],
                            &periodic_shifts[k * NDIM],
                            dx,
                            x_lower,
                            ilower,
                            mask_data->getArrayData(),
                            stencil_lower,
                            stencil_upper,
                            Psi);

            for (int comp = 0; comp < Q_depth; ++comp)
            {
                spread_data(stencil_sz,
                            ig_lower,
                            ig_upper,
                            stencil_lower,
                            stencil_upper,
                            dx,
                            q_data->getArrayData(),
                            comp,
                            Psi,
                            Q_data[s * Q_depth + comp]);
            }
        }
    }

    return;
}

void
LEInteractor::spread(Pointer<NodeData<NDIM, double> > q_data,
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
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_depth == q_data->getDepth());
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(Q_size / Q_depth == X_size / X_depth);
#else
    NULL_USE(Q_size);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, spread_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_shifts(NDIM * local_indices.size());

    // Spread.
    if (!local_indices.empty())
    {
        boost::array<double, NDIM> x_lower_node, x_upper_node;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            x_lower_node[d] = x_lower[d] - 0.5 * dx[d];
            x_upper_node[d] = x_upper[d] + 0.5 * dx[d];
        }
        spread(q_data->getPointer(),
               NodeGeometry<NDIM>::toNodeBox(q_data->getBox()),
               q_data->getGhostCellWidth(),
               q_data->getDepth(),
               Q_data,
               Q_depth,
               X_data,
               x_lower_node.data(),
               x_upper_node.data(),
               dx,
               patch_touches_lower_physical_bdry,
               patch_touches_upper_physical_bdry,
               local_indices,
               periodic_shifts,
               spread_fcn);
    }
    return;
}

void
LEInteractor::spread(Pointer<SideData<NDIM, double> > q_data,
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
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_depth == NDIM);
    TBOX_ASSERT(X_depth == NDIM);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, spread_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_shifts(NDIM * local_indices.size());

    // Spread.
    if (!local_indices.empty())
    {
        boost::array<double, NDIM> x_lower_axis, x_upper_axis;
        const int local_sz = (*std::max_element(local_indices.begin(), local_indices.end())) + 1;
        std::vector<double> Q_data_axis(local_sz);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                x_lower_axis[d] = x_lower[d];
                x_upper_axis[d] = x_upper[d];
            }
            x_lower_axis[axis] -= 0.5 * dx[axis];
            x_upper_axis[axis] += 0.5 * dx[axis];
            for (unsigned int k = 0; k < local_indices.size(); ++k)
            {
                Q_data_axis[local_indices[k]] = Q_data[NDIM * local_indices[k] + axis];
            }
            spread(q_data->getPointer(axis),
                   SideGeometry<NDIM>::toSideBox(q_data->getBox(), axis),
                   q_data->getGhostCellWidth(),
                   /*q_depth*/ 1,
                   &Q_data_axis[0],
                   /*Q_depth*/ 1,
                   X_data,
                   x_lower_axis.data(),
                   x_upper_axis.data(),
                   dx,
                   patch_touches_lower_physical_bdry,
                   patch_touches_upper_physical_bdry,
                   local_indices,
                   periodic_shifts,
                   spread_fcn,
                   axis);
        }
    }
    return;
}

void
LEInteractor::spread(Pointer<SideData<NDIM, double> > mask_data,
                     Pointer<SideData<NDIM, double> > q_data,
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
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_depth == NDIM);
    TBOX_ASSERT(X_depth == NDIM);
    TBOX_ASSERT(Q_size / Q_depth == X_size / X_depth);
    TBOX_ASSERT(q_data->getDepth() == 1);
    TBOX_ASSERT(mask_data);
#else
    NULL_USE(Q_size);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    const Box<NDIM>& patch_box = patch->getBox();
    const IntVector<NDIM>& ilower = patch_box.lower();

    // Get ghost cell width info.
    const IntVector<NDIM>& q_gcw = q_data->getGhostCellWidth();
    const IntVector<NDIM>& mask_gcw = mask_data->getGhostCellWidth();
    const int stencil_size = getStencilSize(spread_fcn);
    const int min_ghosts = getMinimumGhostWidth(spread_fcn);
    const int q_gcw_min = q_gcw.min();
    const int mask_gcw_min = mask_gcw.min();
    if (q_gcw_min < min_ghosts)
    {
        TBOX_ERROR("LEInteractor::interpolate(): insufficient ghost cells for Eulerian field data:"
                   << "  kernel function          = "
                   << spread_fcn
                   << "\n"
                   << "  kernel stencil size      = "
                   << stencil_size
                   << "\n"
                   << "  minimum ghost cell width = "
                   << min_ghosts
                   << "\n"
                   << "  ghost cell width         = "
                   << q_gcw_min
                   << "\n");
    }
    if (mask_gcw_min < stencil_size)
    {
        TBOX_ERROR("LEInteractor::interpolate(): insufficient ghost cells for Eulerian mask data:"
                   << "  kernel function          = "
                   << spread_fcn
                   << "\n"
                   << "  kernel stencil size      = "
                   << stencil_size
                   << "\n"
                   << "  minimum ghost cell width = "
                   << stencil_size
                   << "\n"
                   << "  ghost cell width         = "
                   << mask_gcw_min
                   << "\n");
    }

    // Determine the boundary info.
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, spread_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_shifts(NDIM * local_indices.size());

    // Spread.
    const int nindices = static_cast<int>(local_indices.size());
    if (nindices)
    {
        boost::array<double, NDIM> x_lower_axis, x_upper_axis;
        IntVector<NDIM> stencil_lower, stencil_upper;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            Box<NDIM> data_box = SideGeometry<NDIM>::toSideBox(q_data->getBox(), axis);
            const IntVector<NDIM> ig_lower = data_box.lower() - q_gcw;
            const IntVector<NDIM> ig_upper = data_box.upper() + q_gcw;

            for (int d = 0; d < NDIM; ++d)
            {
                x_lower_axis[d] = x_lower[d];
                x_upper_axis[d] = x_upper[d];
            }
            x_lower_axis[axis] -= 0.5 * dx[axis];
            x_upper_axis[axis] += 0.5 * dx[axis];

            for (int k = 0; k < nindices; ++k)
            {
                int s = local_indices[k];
                MLSWeight Psi;
                const int stencil_sz = LEInteractor::getStencilSize(spread_fcn);
                get_mls_weights(spread_fcn,
                                &X_data[s * NDIM],
                                &periodic_shifts[k * NDIM],
                                dx,
                                x_lower_axis.data(),
                                ilower,
                                mask_data->getArrayData(axis),
                                stencil_lower,
                                stencil_upper,
                                Psi);
                spread_data(stencil_sz,
                            ig_lower,
                            ig_upper,
                            stencil_lower,
                            stencil_upper,
                            dx,
                            q_data->getArrayData(axis),
                            0,
                            Psi,
                            Q_data[s * Q_depth + axis]);
            }
        }
    }
    return;
}

void
LEInteractor::spread(Pointer<EdgeData<NDIM, double> > q_data,
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
    if (NDIM != 3 || Q_depth != NDIM || q_data->getDepth() != 1)
    {
        TBOX_ERROR("LEInteractor::spread():\n"
                   << "  edge-centered interpolation requires 3D vector-valued data.\n");
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(q_data);
    TBOX_ASSERT(patch);
    TBOX_ASSERT(Q_depth == NDIM);
    TBOX_ASSERT(X_depth == NDIM);
#endif
    // Determine the patch geometry.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();
    boost::array<int, NDIM> patch_touches_lower_physical_bdry(array_zero<int, NDIM>());
    boost::array<int, NDIM> patch_touches_upper_physical_bdry(array_zero<int, NDIM>());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        static const int lower = 0;
        patch_touches_lower_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, lower);
        static const int upper = 1;
        patch_touches_upper_physical_bdry[axis] = pgeom->getTouchesRegularBoundary(axis, upper);
    }

    // Generate a list of local indices which lie in the specified box and set
    // all periodic offsets to zero.
    std::vector<int> local_indices;
    buildLocalIndices(local_indices, spread_box, patch, X_data, X_size, X_depth);
    std::vector<double> periodic_shifts(NDIM * local_indices.size());

    // Spread.
    if (!local_indices.empty())
    {
        boost::array<double, NDIM> x_lower_axis, x_upper_axis;
        const int local_sz = (*std::max_element(local_indices.begin(), local_indices.end())) + 1;
        std::vector<double> Q_data_axis(local_sz);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                x_lower_axis[d] = x_lower[d];
                x_upper_axis[d] = x_upper[d];
                if (d != axis)
                {
                    x_lower_axis[axis] -= 0.5 * dx[axis];
                    x_upper_axis[axis] += 0.5 * dx[axis];
                }
            }
            for (unsigned int k = 0; k < local_indices.size(); ++k)
            {
                Q_data_axis[local_indices[k]] = Q_data[NDIM * local_indices[k] + axis];
            }
            spread(q_data->getPointer(axis),
                   EdgeGeometry<NDIM>::toEdgeBox(q_data->getBox(), axis),
                   q_data->getGhostCellWidth(),
                   /*q_depth*/ 1,
                   &Q_data_axis[0],
                   /*Q_depth*/ 1,
                   X_data,
                   x_lower_axis.data(),
                   x_upper_axis.data(),
                   dx,
                   patch_touches_lower_physical_bdry,
                   patch_touches_upper_physical_bdry,
                   local_indices,
                   periodic_shifts,
                   spread_fcn,
                   axis);
        }
    }
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
LEInteractor::interpolate(double* const Q_data,
                          const int Q_depth,
                          const double* const X_data,
                          const double* const q_data,
                          const Box<NDIM>& q_data_box,
                          const IntVector<NDIM>& q_gcw,
                          const int q_depth,
                          const double* const x_lower,
                          const double* const x_upper,
                          const double* const dx,
                          const boost::array<int, NDIM>& /*patch_touches_lower_physical_bdry*/,
                          const boost::array<int, NDIM>& /*patch_touches_upper_physical_bdry*/,
                          const std::vector<int>& local_indices,
                          const std::vector<double>& periodic_shifts,
                          const std::string& interp_fcn,
                          const int axis)
{
    const int stencil_size = getStencilSize(interp_fcn);
    const int min_ghosts = getMinimumGhostWidth(interp_fcn);
    const int q_gcw_min = q_gcw.min();
    if (q_gcw_min < min_ghosts)
    {
        TBOX_ERROR("LEInteractor::interpolate(): insufficient ghost cells:"
                   << "  kernel function          = "
                   << interp_fcn
                   << "\n"
                   << "  kernel stencil size      = "
                   << stencil_size
                   << "\n"
                   << "  minimum ghost cell width = "
                   << min_ghosts
                   << "\n"
                   << "  ghost cell width         = "
                   << q_gcw_min
                   << "\n");
    }
    if (local_indices.empty()) return;
    const int local_indices_size = static_cast<int>(local_indices.size());
    const IntVector<NDIM>& ilower = q_data_box.lower();
    const IntVector<NDIM>& iupper = q_data_box.upper();
    if (interp_fcn == "PIECEWISE_CONSTANT")
    {
        LAGRANGIAN_PIECEWISE_CONSTANT_INTERP_FC(dx,
                                                x_lower,
                                                x_upper,
                                                q_depth,
#if (NDIM == 2)
                                                ilower(0),
                                                iupper(0),
                                                ilower(1),
                                                iupper(1),
                                                q_gcw(0),
                                                q_gcw(1),
#endif
#if (NDIM == 3)
                                                ilower(0),
                                                iupper(0),
                                                ilower(1),
                                                iupper(1),
                                                ilower(2),
                                                iupper(2),
                                                q_gcw(0),
                                                q_gcw(1),
                                                q_gcw(2),
#endif
                                                q_data,
                                                &local_indices[0],
                                                &periodic_shifts[0],
                                                local_indices_size,
                                                X_data,
                                                Q_data);
    }
    else if (interp_fcn == "DISCONTINUOUS_LINEAR")
    {
        LAGRANGIAN_DISCONTINUOUS_LINEAR_INTERP_FC(dx,
                                                  x_lower,
                                                  x_upper,
                                                  q_depth,
                                                  axis,
#if (NDIM == 2)
                                                  ilower(0),
                                                  iupper(0),
                                                  ilower(1),
                                                  iupper(1),
                                                  q_gcw(0),
                                                  q_gcw(1),
#endif
#if (NDIM == 3)
                                                  ilower(0),
                                                  iupper(0),
                                                  ilower(1),
                                                  iupper(1),
                                                  ilower(2),
                                                  iupper(2),
                                                  q_gcw(0),
                                                  q_gcw(1),
                                                  q_gcw(2),
#endif
                                                  q_data,
                                                  &local_indices[0],
                                                  &periodic_shifts[0],
                                                  local_indices_size,
                                                  X_data,
                                                  Q_data);
    }
    else if (interp_fcn == "PIECEWISE_LINEAR")
    {
        LAGRANGIAN_PIECEWISE_LINEAR_INTERP_FC(dx,
                                              x_lower,
                                              x_upper,
                                              q_depth,
#if (NDIM == 2)
                                              ilower(0),
                                              iupper(0),
                                              ilower(1),
                                              iupper(1),
                                              q_gcw(0),
                                              q_gcw(1),
#endif
#if (NDIM == 3)
                                              ilower(0),
                                              iupper(0),
                                              ilower(1),
                                              iupper(1),
                                              ilower(2),
                                              iupper(2),
                                              q_gcw(0),
                                              q_gcw(1),
                                              q_gcw(2),
#endif
                                              q_data,
                                              &local_indices[0],
                                              &periodic_shifts[0],
                                              local_indices_size,
                                              X_data,
                                              Q_data);
    }
    else if (interp_fcn == "PIECEWISE_CUBIC")
    {
        LAGRANGIAN_PIECEWISE_CUBIC_INTERP_FC(dx,
                                             x_lower,
                                             x_upper,
                                             q_depth,
#if (NDIM == 2)
                                             ilower(0),
                                             iupper(0),
                                             ilower(1),
                                             iupper(1),
                                             q_gcw(0),
                                             q_gcw(1),
#endif
#if (NDIM == 3)
                                             ilower(0),
                                             iupper(0),
                                             ilower(1),
                                             iupper(1),
                                             ilower(2),
                                             iupper(2),
                                             q_gcw(0),
                                             q_gcw(1),
                                             q_gcw(2),
#endif
                                             q_data,
                                             &local_indices[0],
                                             &periodic_shifts[0],
                                             local_indices_size,
                                             X_data,
                                             Q_data);
    }
    else if (interp_fcn == "IB_3")
    {
        LAGRANGIAN_IB_3_INTERP_FC(dx,
                                  x_lower,
                                  x_upper,
                                  q_depth,
#if (NDIM == 2)
                                  ilower(0),
                                  iupper(0),
                                  ilower(1),
                                  iupper(1),
                                  q_gcw(0),
                                  q_gcw(1),
#endif
#if (NDIM == 3)
                                  ilower(0),
                                  iupper(0),
                                  ilower(1),
                                  iupper(1),
                                  ilower(2),
                                  iupper(2),
                                  q_gcw(0),
                                  q_gcw(1),
                                  q_gcw(2),
#endif
                                  q_data,
                                  &local_indices[0],
                                  &periodic_shifts[0],
                                  local_indices_size,
                                  X_data,
                                  Q_data);
    }
    else if (interp_fcn == "IB_4")
    {
        LAGRANGIAN_IB_4_INTERP_FC(dx,
                                  x_lower,
                                  x_upper,
                                  q_depth,
#if (NDIM == 2)
                                  ilower(0),
                                  iupper(0),
                                  ilower(1),
                                  iupper(1),
                                  q_gcw(0),
                                  q_gcw(1),
#endif
#if (NDIM == 3)
                                  ilower(0),
                                  iupper(0),
                                  ilower(1),
                                  iupper(1),
                                  ilower(2),
                                  iupper(2),
                                  q_gcw(0),
                                  q_gcw(1),
                                  q_gcw(2),
#endif
                                  q_data,
                                  &local_indices[0],
                                  &periodic_shifts[0],
                                  local_indices_size,
                                  X_data,
                                  Q_data);
    }
    else if (interp_fcn == "IB_4_W8")
    {
        LAGRANGIAN_IB_4_W8_INTERP_FC(dx,
                                     x_lower,
                                     x_upper,
                                     q_depth,
#if (NDIM == 2)
                                     ilower(0),
                                     iupper(0),
                                     ilower(1),
                                     iupper(1),
                                     q_gcw(0),
                                     q_gcw(1),
#endif
#if (NDIM == 3)
                                     ilower(0),
                                     iupper(0),
                                     ilower(1),
                                     iupper(1),
                                     ilower(2),
                                     iupper(2),
                                     q_gcw(0),
                                     q_gcw(1),
                                     q_gcw(2),
#endif
                                     q_data,
                                     &local_indices[0],
                                     &periodic_shifts[0],
                                     local_indices_size,
                                     X_data,
                                     Q_data);
    }
    else if (interp_fcn == "IB_6")
    {
        LAGRANGIAN_IB_6_INTERP_FC(dx,
                                  x_lower,
                                  x_upper,
                                  q_depth,
#if (NDIM == 2)
                                  ilower(0),
                                  iupper(0),
                                  ilower(1),
                                  iupper(1),
                                  q_gcw(0),
                                  q_gcw(1),
#endif
#if (NDIM == 3)
                                  ilower(0),
                                  iupper(0),
                                  ilower(1),
                                  iupper(1),
                                  ilower(2),
                                  iupper(2),
                                  q_gcw(0),
                                  q_gcw(1),
                                  q_gcw(2),
#endif
                                  q_data,
                                  &local_indices[0],
                                  &periodic_shifts[0],
                                  local_indices_size,
                                  X_data,
                                  Q_data);
    }
    else if (interp_fcn == "USER_DEFINED")
    {
        userDefinedInterpolate(Q_data,
                               Q_depth,
                               X_data,
                               q_data,
                               q_data_box,
                               q_gcw,
                               q_depth,
                               x_lower,
                               x_upper,
                               dx,
                               &local_indices[0],
                               &periodic_shifts[0],
                               local_indices_size);
    }
    else
    {
        TBOX_ERROR("LEInteractor::interpolate()\n"
                   << "  Unknown interpolation kernel function "
                   << interp_fcn
                   << std::endl);
    }
    return;
}

void
LEInteractor::spread(double* const q_data,
                     const Box<NDIM>& q_data_box,
                     const IntVector<NDIM>& q_gcw,
                     const int q_depth,
                     const double* const Q_data,
                     const int Q_depth,
                     const double* const X_data,
                     const double* const x_lower,
                     const double* const x_upper,
                     const double* const dx,
                     const boost::array<int, NDIM>& patch_touches_lower_physical_bdry,
                     const boost::array<int, NDIM>& patch_touches_upper_physical_bdry,
                     const std::vector<int>& local_indices,
                     const std::vector<double>& periodic_shifts,
                     const std::string& spread_fcn,
                     const int axis)
{
    const int stencil_size = getStencilSize(spread_fcn);
    const int min_ghosts = getMinimumGhostWidth(spread_fcn);
    const int q_gcw_min = q_gcw.min();
    bool patch_touches_physical_bdry = false;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        patch_touches_physical_bdry = patch_touches_physical_bdry || patch_touches_lower_physical_bdry[d];
        patch_touches_physical_bdry = patch_touches_physical_bdry || patch_touches_upper_physical_bdry[d];
    }
    if (patch_touches_physical_bdry && q_gcw_min < min_ghosts)
    {
        TBOX_ERROR("LEInteractor::spread(): insufficient ghost cells at physical boundary:"
                   << "  kernel function          = "
                   << spread_fcn
                   << "\n"
                   << "  kernel stencil size      = "
                   << stencil_size
                   << "\n"
                   << "  minimum ghost cell width = "
                   << min_ghosts
                   << "\n"
                   << "  ghost cell width         = "
                   << q_gcw_min
                   << "\n");
    }
    if (local_indices.empty()) return;
    const int local_indices_size = static_cast<int>(local_indices.size());
    const IntVector<NDIM>& ilower = q_data_box.lower();
    const IntVector<NDIM>& iupper = q_data_box.upper();
    if (spread_fcn == "PIECEWISE_CONSTANT")
    {
        LAGRANGIAN_PIECEWISE_CONSTANT_SPREAD_FC(dx,
                                                x_lower,
                                                x_upper,
                                                q_depth,
                                                &local_indices[0],
                                                &periodic_shifts[0],
                                                local_indices_size,
                                                X_data,
                                                Q_data,
#if (NDIM == 2)
                                                ilower(0),
                                                iupper(0),
                                                ilower(1),
                                                iupper(1),
                                                q_gcw(0),
                                                q_gcw(1),
#endif
#if (NDIM == 3)
                                                ilower(0),
                                                iupper(0),
                                                ilower(1),
                                                iupper(1),
                                                ilower(2),
                                                iupper(2),
                                                q_gcw(0),
                                                q_gcw(1),
                                                q_gcw(2),
#endif
                                                q_data);
    }
    else if (spread_fcn == "DISCONTINUOUS_LINEAR")
    {
        LAGRANGIAN_DISCONTINUOUS_LINEAR_SPREAD_FC(dx,
                                                  x_lower,
                                                  x_upper,
                                                  q_depth,
                                                  axis,
                                                  &local_indices[0],
                                                  &periodic_shifts[0],
                                                  local_indices_size,
                                                  X_data,
                                                  Q_data,
#if (NDIM == 2)
                                                  ilower(0),
                                                  iupper(0),
                                                  ilower(1),
                                                  iupper(1),
                                                  q_gcw(0),
                                                  q_gcw(1),
#endif
#if (NDIM == 3)
                                                  ilower(0),
                                                  iupper(0),
                                                  ilower(1),
                                                  iupper(1),
                                                  ilower(2),
                                                  iupper(2),
                                                  q_gcw(0),
                                                  q_gcw(1),
                                                  q_gcw(2),
#endif
                                                  q_data);
    }
    else if (spread_fcn == "PIECEWISE_LINEAR")
    {
        LAGRANGIAN_PIECEWISE_LINEAR_SPREAD_FC(dx,
                                              x_lower,
                                              x_upper,
                                              q_depth,
                                              &local_indices[0],
                                              &periodic_shifts[0],
                                              local_indices_size,
                                              X_data,
                                              Q_data,
#if (NDIM == 2)
                                              ilower(0),
                                              iupper(0),
                                              ilower(1),
                                              iupper(1),
                                              q_gcw(0),
                                              q_gcw(1),
#endif
#if (NDIM == 3)
                                              ilower(0),
                                              iupper(0),
                                              ilower(1),
                                              iupper(1),
                                              ilower(2),
                                              iupper(2),
                                              q_gcw(0),
                                              q_gcw(1),
                                              q_gcw(2),
#endif
                                              q_data);
    }
    else if (spread_fcn == "PIECEWISE_CUBIC")
    {
        LAGRANGIAN_PIECEWISE_CUBIC_SPREAD_FC(dx,
                                             x_lower,
                                             x_upper,
                                             q_depth,
                                             &local_indices[0],
                                             &periodic_shifts[0],
                                             local_indices_size,
                                             X_data,
                                             Q_data,
#if (NDIM == 2)
                                             ilower(0),
                                             iupper(0),
                                             ilower(1),
                                             iupper(1),
                                             q_gcw(0),
                                             q_gcw(1),
#endif
#if (NDIM == 3)
                                             ilower(0),
                                             iupper(0),
                                             ilower(1),
                                             iupper(1),
                                             ilower(2),
                                             iupper(2),
                                             q_gcw(0),
                                             q_gcw(1),
                                             q_gcw(2),
#endif
                                             q_data);
    }
    else if (spread_fcn == "IB_3")
    {
        LAGRANGIAN_IB_3_SPREAD_FC(dx,
                                  x_lower,
                                  x_upper,
                                  q_depth,
                                  &local_indices[0],
                                  &periodic_shifts[0],
                                  local_indices_size,
                                  X_data,
                                  Q_data,
#if (NDIM == 2)
                                  ilower(0),
                                  iupper(0),
                                  ilower(1),
                                  iupper(1),
                                  q_gcw(0),
                                  q_gcw(1),
#endif
#if (NDIM == 3)
                                  ilower(0),
                                  iupper(0),
                                  ilower(1),
                                  iupper(1),
                                  ilower(2),
                                  iupper(2),
                                  q_gcw(0),
                                  q_gcw(1),
                                  q_gcw(2),
#endif
                                  q_data);
    }
    else if (spread_fcn == "IB_4")
    {
        LAGRANGIAN_IB_4_SPREAD_FC(dx,
                                  x_lower,
                                  x_upper,
                                  q_depth,
                                  &local_indices[0],
                                  &periodic_shifts[0],
                                  local_indices_size,
                                  X_data,
                                  Q_data,
#if (NDIM == 2)
                                  ilower(0),
                                  iupper(0),
                                  ilower(1),
                                  iupper(1),
                                  q_gcw(0),
                                  q_gcw(1),
#endif
#if (NDIM == 3)
                                  ilower(0),
                                  iupper(0),
                                  ilower(1),
                                  iupper(1),
                                  ilower(2),
                                  iupper(2),
                                  q_gcw(0),
                                  q_gcw(1),
                                  q_gcw(2),
#endif
                                  q_data);
    }
    else if (spread_fcn == "IB_4_W8")
    {
        LAGRANGIAN_IB_4_W8_SPREAD_FC(dx,
                                     x_lower,
                                     x_upper,
                                     q_depth,
                                     &local_indices[0],
                                     &periodic_shifts[0],
                                     local_indices_size,
                                     X_data,
                                     Q_data,
#if (NDIM == 2)
                                     ilower(0),
                                     iupper(0),
                                     ilower(1),
                                     iupper(1),
                                     q_gcw(0),
                                     q_gcw(1),
#endif
#if (NDIM == 3)
                                     ilower(0),
                                     iupper(0),
                                     ilower(1),
                                     iupper(1),
                                     ilower(2),
                                     iupper(2),
                                     q_gcw(0),
                                     q_gcw(1),
                                     q_gcw(2),
#endif
                                     q_data);
    }
    else if (spread_fcn == "IB_6")
    {
        LAGRANGIAN_IB_6_SPREAD_FC(dx,
                                  x_lower,
                                  x_upper,
                                  q_depth,
                                  &local_indices[0],
                                  &periodic_shifts[0],
                                  local_indices_size,
                                  X_data,
                                  Q_data,
#if (NDIM == 2)
                                  ilower(0),
                                  iupper(0),
                                  ilower(1),
                                  iupper(1),
                                  q_gcw(0),
                                  q_gcw(1),
#endif
#if (NDIM == 3)
                                  ilower(0),
                                  iupper(0),
                                  ilower(1),
                                  iupper(1),
                                  ilower(2),
                                  iupper(2),
                                  q_gcw(0),
                                  q_gcw(1),
                                  q_gcw(2),
#endif
                                  q_data);
    }
    else if (spread_fcn == "USER_DEFINED")
    {
        userDefinedSpread(q_data,
                          q_data_box,
                          q_gcw,
                          q_depth,
                          x_lower,
                          x_upper,
                          dx,
                          Q_data,
                          Q_depth,
                          X_data,
                          &local_indices[0],
                          &periodic_shifts[0],
                          local_indices_size);
    }
    else
    {
        TBOX_ERROR("LEInteractor::spread()\n"
                   << "  Unknown spreading kernel function "
                   << spread_fcn
                   << std::endl);
    }
    return;
}

template <class T>
void
LEInteractor::buildLocalIndices(std::vector<int>& local_indices,
                                std::vector<double>& periodic_shifts,
                                const Box<NDIM>& box,
                                const Pointer<Patch<NDIM> > patch,
                                const IntVector<NDIM>& periodic_shift,
                                const Pointer<LIndexSetData<T> > idx_data)
{
    local_indices.clear();
    periodic_shifts.clear();
    const size_t upper_bound = idx_data->getLocalPETScIndices().size();
    if (upper_bound == 0) return;
    local_indices.reserve(upper_bound);
    periodic_shifts.reserve(NDIM * upper_bound);

    const Box<NDIM>& patch_box = patch->getBox();
    const Index<NDIM>& ilower = patch_box.lower();
    const Index<NDIM>& iupper = patch_box.upper();
    const Box<NDIM>& ghost_box = idx_data->getGhostBox();

    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    boost::array<bool, NDIM> patch_touches_lower_periodic_bdry, patch_touches_upper_periodic_bdry;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        patch_touches_lower_periodic_bdry[axis] = pgeom->getTouchesPeriodicBoundary(axis, 0);
        patch_touches_upper_periodic_bdry[axis] = pgeom->getTouchesPeriodicBoundary(axis, 1);
    }

    if (box == patch_box)
    {
        local_indices = idx_data->getInteriorLocalPETScIndices();
        periodic_shifts = idx_data->getInteriorPeriodicShifts();
    }
    else if (box == ghost_box)
    {
        local_indices = idx_data->getLocalPETScIndices();
        periodic_shifts = idx_data->getPeriodicShifts();
    }
    else
    {
        for (typename LIndexSetData<T>::SetIterator it(*idx_data); it; it++)
        {
            const Index<NDIM>& i = it.getIndex();
            if (!box.contains(i)) continue;

            boost::array<int, NDIM> offset;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (patch_touches_lower_periodic_bdry[d] && i(d) < ilower(d))
                {
                    offset[d] = -periodic_shift(d); // X is ABOVE the top    of the patch --- need
                                                    // to shift DOWN
                }
                else if (patch_touches_upper_periodic_bdry[d] && i(d) > iupper(d))
                {
                    offset[d] = +periodic_shift(d); // X is BELOW the bottom of the patch ---
                                                    // need to shift UP
                }
                else
                {
                    offset[d] = 0;
                }
            }
            const LSet<T>& idx_set = it.getItem();
            for (typename LSet<T>::const_iterator n = idx_set.begin(); n != idx_set.end(); ++n)
            {
                const typename LSet<T>::value_type& idx = *n;
                local_indices.push_back(idx->getLocalPETScIndex());
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    periodic_shifts.push_back(static_cast<double>(offset[d]) * dx[d]);
                }
            }
        }
    }
    return;
}

void
LEInteractor::buildLocalIndices(std::vector<int>& local_indices,
                                const Box<NDIM>& box,
                                const Pointer<Patch<NDIM> > patch,
                                const double* const X_data,
                                const int X_size,
                                const int X_depth)
{
    local_indices.clear();
    const int upper_bound = X_size / X_depth;
    if (upper_bound == 0) return;

    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    local_indices.reserve(upper_bound);
    for (int k = 0; k < X_size / X_depth; ++k)
    {
        const double* const X = &X_data[NDIM * k];
        const Index<NDIM> i = IndexUtilities::getCellIndex(X, patch_geom, patch_box);
        if (box.contains(i)) local_indices.push_back(k);
    }
    return;
}

void
LEInteractor::userDefinedInterpolate(double* Q,
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
    typedef boost::multi_array_types::extent_range range;
    boost::const_multi_array_ref<double, NDIM + 1> q_data(
        q,
        (boost::extents[range(ilower[0] - q_gcw[0], iupper[0] + q_gcw[0] + 1)][range(ilower[1] - q_gcw[1],
                                                                                     iupper[1] + q_gcw[1] + 1)]
#if (NDIM == 3)
                       [range(ilower[2] - q_gcw[2], iupper[2] + q_gcw[2] + 1)]
#endif
                       [range(0, q_depth)]),
        boost::fortran_storage_order());
    boost::array<double, NDIM> X_cell;
    boost::array<int, NDIM> stencil_center, stencil_lower, stencil_upper;
    for (int l = 0; l < num_local_indices; ++l)
    {
        const int s = local_indices[l];

        // Determine the Cartesian cell in which X(s) is located.
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            stencil_center[d] =
                static_cast<int>(std::floor((X[d + s * NDIM] + X_shift[d + l * NDIM] - x_lower[d]) / dx[d])) +
                ilower[d];
            X_cell[d] = x_lower[d] + (static_cast<double>(stencil_center[d] - ilower[d]) + 0.5) * dx[d];
        }

        // Determine the interpolation stencil corresponding to the position of
        // X(s) within the cell.
        if (s_kernel_fcn_stencil_size % 2 == 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (X[d + s * NDIM] < X_cell[d])
                {
                    stencil_lower[d] = stencil_center[d] - s_kernel_fcn_stencil_size / 2;
                    stencil_upper[d] = stencil_center[d] + s_kernel_fcn_stencil_size / 2 - 1;
                }
                else
                {
                    stencil_lower[d] = stencil_center[d] - s_kernel_fcn_stencil_size / 2 + 1;
                    stencil_upper[d] = stencil_center[d] + s_kernel_fcn_stencil_size / 2;
                }
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                stencil_lower[d] = stencil_center[d] - s_kernel_fcn_stencil_size / 2;
                stencil_upper[d] = stencil_center[d] + s_kernel_fcn_stencil_size / 2;
            }
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            stencil_lower[d] = std::min(std::max(stencil_lower[d], ilower[d] - q_gcw[d]), iupper[d] + q_gcw[d]);
            stencil_upper[d] = std::min(std::max(stencil_upper[d], ilower[d] - q_gcw[d]), iupper[d] + q_gcw[d]);
        }

        // Compute the kernel function weights.
        boost::multi_array<double, 1> w0(boost::extents[range(stencil_lower[0], stencil_upper[0] + 1)]);
        for (int ic0 = stencil_lower[0]; ic0 <= stencil_upper[0]; ++ic0)
        {
            w0[ic0] = s_kernel_fcn((X[0 + s * NDIM] + X_shift[0 + l * NDIM] -
                                    (X_cell[0] + static_cast<double>(ic0 - stencil_center[0]) * dx[0])) /
                                   dx[0]);
        }

        boost::multi_array<double, 1> w1(boost::extents[range(stencil_lower[1], stencil_upper[1] + 1)]);
        for (int ic1 = stencil_lower[1]; ic1 <= stencil_upper[1]; ++ic1)
        {
            w1[ic1] = s_kernel_fcn((X[1 + s * NDIM] + X_shift[1 + l * NDIM] -
                                    (X_cell[1] + static_cast<double>(ic1 - stencil_center[1]) * dx[1])) /
                                   dx[1]);
        }
#if (NDIM == 3)
        boost::multi_array<double, 1> w2(boost::extents[range(stencil_lower[2], stencil_upper[2] + 1)]);
        for (int ic2 = stencil_lower[2]; ic2 <= stencil_upper[2]; ++ic2)
        {
            w2[ic2] = s_kernel_fcn((X[2 + s * NDIM] + X_shift[2 + l * NDIM] -
                                    (X_cell[2] + static_cast<double>(ic2 - stencil_center[2]) * dx[2])) /
                                   dx[2]);
        }
#endif
        // Interpolate u onto V.
        for (int d = 0; d < Q_depth; ++d)
        {
            Q[d + s * Q_depth] = 0.0;
#if (NDIM == 3)
            for (int ic2 = stencil_lower[2]; ic2 <= stencil_upper[2]; ++ic2)
            {
#endif
                for (int ic1 = stencil_lower[1]; ic1 <= stencil_upper[1]; ++ic1)
                {
                    for (int ic0 = stencil_lower[0]; ic0 <= stencil_upper[0]; ++ic0)
                    {
#if (NDIM == 2)
                        Q[d + s * Q_depth] += w0[ic0] * w1[ic1] * q_data[ic0][ic1][d];
#endif
#if (NDIM == 3)
                        Q[d + s * Q_depth] += w0[ic0] * w1[ic1] * w2[ic2] * q_data[ic0][ic1][ic2][d];
#endif
                    }
                }
#if (NDIM == 3)
            }
#endif
        }
    }
    return;
}

void
LEInteractor::userDefinedSpread(double* q,
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
    typedef boost::multi_array_types::extent_range range;
    boost::multi_array_ref<double, NDIM + 1> q_data(
        q,
        (boost::extents[range(ilower[0] - q_gcw[0], iupper[0] + q_gcw[0] + 1)][range(ilower[1] - q_gcw[1],
                                                                                     iupper[1] + q_gcw[1] + 1)]
#if (NDIM == 3)
                       [range(ilower[2] - q_gcw[2], iupper[2] + q_gcw[2] + 1)]
#endif
                       [range(0, q_depth)]),
        boost::fortran_storage_order());
    boost::array<double, NDIM> X_cell;
    boost::array<int, NDIM> stencil_center, stencil_lower, stencil_upper;
    for (int l = 0; l < num_local_indices; ++l)
    {
        const int s = local_indices[l];

        // Determine the Cartesian cell in which X(s) is located.
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            stencil_center[d] =
                static_cast<int>(std::floor((X[d + s * NDIM] + X_shift[d + l * NDIM] - x_lower[d]) / dx[d])) +
                ilower[d];
            X_cell[d] = x_lower[d] + (static_cast<double>(stencil_center[d] - ilower[d]) + 0.5) * dx[d];
        }

        // Determine the interpolation stencil corresponding to the position of
        // X(s) within the cell.
        if (s_kernel_fcn_stencil_size % 2 == 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (X[d + s * NDIM] < X_cell[d])
                {
                    stencil_lower[d] = stencil_center[d] - s_kernel_fcn_stencil_size / 2;
                    stencil_upper[d] = stencil_center[d] + s_kernel_fcn_stencil_size / 2 - 1;
                }
                else
                {
                    stencil_lower[d] = stencil_center[d] - s_kernel_fcn_stencil_size / 2 + 1;
                    stencil_upper[d] = stencil_center[d] + s_kernel_fcn_stencil_size / 2;
                }
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                stencil_lower[d] = stencil_center[d] - s_kernel_fcn_stencil_size / 2;
                stencil_upper[d] = stencil_center[d] + s_kernel_fcn_stencil_size / 2;
            }
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            stencil_lower[d] = std::min(std::max(stencil_lower[d], ilower[d] - q_gcw[d]), iupper[d] + q_gcw[d]);
            stencil_upper[d] = std::min(std::max(stencil_upper[d], ilower[d] - q_gcw[d]), iupper[d] + q_gcw[d]);
        }

        // Compute the kernel function weights.
        boost::multi_array<double, 1> w0(boost::extents[range(stencil_lower[0], stencil_upper[0] + 1)]);
        for (int ic0 = stencil_lower[0]; ic0 <= stencil_upper[0]; ++ic0)
        {
            w0[ic0] = s_kernel_fcn((X[0 + s * NDIM] + X_shift[0 + l * NDIM] -
                                    (X_cell[0] + static_cast<double>(ic0 - stencil_center[0]) * dx[0])) /
                                   dx[0]);
        }

        boost::multi_array<double, 1> w1(boost::extents[range(stencil_lower[1], stencil_upper[1] + 1)]);
        for (int ic1 = stencil_lower[1]; ic1 <= stencil_upper[1]; ++ic1)
        {
            w1[ic1] = s_kernel_fcn((X[1 + s * NDIM] + X_shift[1 + l * NDIM] -
                                    (X_cell[1] + static_cast<double>(ic1 - stencil_center[1]) * dx[1])) /
                                   dx[1]);
        }
#if (NDIM == 3)
        boost::multi_array<double, 1> w2(boost::extents[range(stencil_lower[2], stencil_upper[2] + 1)]);
        for (int ic2 = stencil_lower[2]; ic2 <= stencil_upper[2]; ++ic2)
        {
            w2[ic2] = s_kernel_fcn((X[2 + s * NDIM] + X_shift[2 + l * NDIM] -
                                    (X_cell[2] + static_cast<double>(ic2 - stencil_center[2]) * dx[2])) /
                                   dx[2]);
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
                        q_data[ic0][ic1][d] += w0[ic0] * w1[ic1] * Q[d + s * Q_depth] / (dx[0] * dx[1]);
#endif
#if (NDIM == 3)
                        q_data[ic0][ic1][ic2][d] +=
                            w0[ic0] * w1[ic1] * w2[ic2] * Q[d + s * Q_depth] / (dx[0] * dx[1] * dx[2]);
#endif
                    }
                }
#if (NDIM == 3)
            }
#endif
        }
    }
    return;
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

template void IBTK::LEInteractor::interpolate(SAMRAI::tbox::Pointer<LData> Q_data,
                                              const SAMRAI::tbox::Pointer<LData> X_data,
                                              const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                              const SAMRAI::hier::Box<NDIM>& interp_box,
                                              const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                              const std::string& interp_fcn);

template void IBTK::LEInteractor::interpolate(SAMRAI::tbox::Pointer<LData> Q_data,
                                              const SAMRAI::tbox::Pointer<LData> X_data,
                                              const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > q_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                              const SAMRAI::hier::Box<NDIM>& interp_box,
                                              const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                              const std::string& interp_fcn);

template void IBTK::LEInteractor::interpolate(SAMRAI::tbox::Pointer<LData> Q_data,
                                              const SAMRAI::tbox::Pointer<LData> X_data,
                                              const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                              const SAMRAI::hier::Box<NDIM>& interp_box,
                                              const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                              const std::string& interp_fcn);

template void IBTK::LEInteractor::interpolate(SAMRAI::tbox::Pointer<LData> Q_data,
                                              const SAMRAI::tbox::Pointer<LData> X_data,
                                              const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > q_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                              const SAMRAI::hier::Box<NDIM>& interp_box,
                                              const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                              const std::string& interp_fcn);

template void IBTK::LEInteractor::interpolate(double* const Q_data,
                                              const int Q_depth,
                                              const double* const X_data,
                                              const int X_depth,
                                              const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                              const SAMRAI::hier::Box<NDIM>& interp_box,
                                              const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                              const std::string& interp_fcn);

template void IBTK::LEInteractor::interpolate(double* const Q_data,
                                              const int Q_depth,
                                              const double* const X_data,
                                              const int X_depth,
                                              const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > q_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                              const SAMRAI::hier::Box<NDIM>& interp_box,
                                              const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                              const std::string& interp_fcn);

template void IBTK::LEInteractor::interpolate(double* const Q_data,
                                              const int Q_depth,
                                              const double* const X_data,
                                              const int X_depth,
                                              const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                              const SAMRAI::hier::Box<NDIM>& interp_box,
                                              const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                              const std::string& interp_fcn);

template void IBTK::LEInteractor::interpolate(double* const Q_data,
                                              const int Q_depth,
                                              const double* const X_data,
                                              const int X_depth,
                                              const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > q_data,
                                              const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                              const SAMRAI::hier::Box<NDIM>& interp_box,
                                              const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                              const std::string& interp_fcn);

template void IBTK::LEInteractor::spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                                         const SAMRAI::tbox::Pointer<LData> Q_data,
                                         const SAMRAI::tbox::Pointer<LData> X_data,
                                         const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                         const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                         const SAMRAI::hier::Box<NDIM>& spread_box,
                                         const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                         const std::string& spread_fcn);

template void IBTK::LEInteractor::spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > q_data,
                                         const SAMRAI::tbox::Pointer<LData> Q_data,
                                         const SAMRAI::tbox::Pointer<LData> X_data,
                                         const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                         const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                         const SAMRAI::hier::Box<NDIM>& spread_box,
                                         const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                         const std::string& spread_fcn);

template void IBTK::LEInteractor::spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                                         const SAMRAI::tbox::Pointer<LData> Q_data,
                                         const SAMRAI::tbox::Pointer<LData> X_data,
                                         const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                         const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                         const SAMRAI::hier::Box<NDIM>& spread_box,
                                         const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                         const std::string& spread_fcn);

template void IBTK::LEInteractor::spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > q_data,
                                         const SAMRAI::tbox::Pointer<LData> Q_data,
                                         const SAMRAI::tbox::Pointer<LData> X_data,
                                         const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                         const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                         const SAMRAI::hier::Box<NDIM>& spread_box,
                                         const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                         const std::string& spread_fcn);

template void IBTK::LEInteractor::spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                                         const double* const Q_data,
                                         const int Q_depth,
                                         const double* const X_data,
                                         const int X_depth,
                                         const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                         const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                         const SAMRAI::hier::Box<NDIM>& spread_box,
                                         const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                         const std::string& spread_fcn);

template void IBTK::LEInteractor::spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > q_data,
                                         const double* const Q_data,
                                         const int Q_depth,
                                         const double* const X_data,
                                         const int X_depth,
                                         const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                         const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                         const SAMRAI::hier::Box<NDIM>& spread_box,
                                         const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                         const std::string& spread_fcn);

template void IBTK::LEInteractor::spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                                         const double* const Q_data,
                                         const int Q_depth,
                                         const double* const X_data,
                                         const int X_depth,
                                         const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                         const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                         const SAMRAI::hier::Box<NDIM>& spread_box,
                                         const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                         const std::string& spread_fcn);

template void IBTK::LEInteractor::spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > q_data,
                                         const double* const Q_data,
                                         const int Q_depth,
                                         const double* const X_data,
                                         const int X_depth,
                                         const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data,
                                         const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                         const SAMRAI::hier::Box<NDIM>& spread_box,
                                         const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                         const std::string& spread_fcn);

template void IBTK::LEInteractor::buildLocalIndices(std::vector<int>& local_indices,
                                                    std::vector<double>& periodic_shifts,
                                                    const SAMRAI::hier::Box<NDIM>& box,
                                                    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                                    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                                    const SAMRAI::tbox::Pointer<LIndexSetData<LNode> > idx_data);

//////////////////////////////////////////////////////////////////////////////
