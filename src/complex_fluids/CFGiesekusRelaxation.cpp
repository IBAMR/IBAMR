// Filename: CFGiesekusRelaxation.cpp
// Created on 23 Oct 2019 by Aaron Barrett
//
// Copyright (c) 2002-2019, Boyce Griffith
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

#include "ibamr/CFGiesekusRelaxation.h"
#include "ibamr/namespaces.h"

// Namespace
namespace IBAMR
{
CFGiesekusRelaxation::CFGiesekusRelaxation(const std::string& object_name, Pointer<Database> input_db)
    : CFRelaxationOperator(object_name, input_db)
{
    d_lambda = input_db->getDouble("relaxation_time");
    d_alpha = input_db->getDouble("alpha");
    return;
} // Constructor

void
CFGiesekusRelaxation::setDataOnPatch(const int data_idx,
                                     Pointer<Variable<NDIM> > /*var*/,
                                     Pointer<Patch<NDIM> > patch,
                                     const double /*data_time*/,
                                     const bool initial_time,
                                     Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CellData<NDIM, double> > ret_data = patch->getPatchData(data_idx);
    Pointer<CellData<NDIM, double> > in_data = patch->getPatchData(d_W_cc_idx);
    ret_data->fillAll(0.0);
    if (initial_time) return;
    const double l_inv = 1.0 / d_lambda;
    for (CellIterator<NDIM> i(patch_box); i; i++)
    {
        const CellIndex<NDIM>& idx = i();
        MatrixNd mat;
#if (NDIM == 2)
        mat(0, 0) = (*in_data)(idx, 0);
        mat(1, 1) = (*in_data)(idx, 1);
        mat(0, 1) = mat(1, 0) = (*in_data)(idx, 2);
#endif
#if (NDIM == 3)
        mat(0, 0) = (*in_data)(idx, 0);
        mat(1, 1) = (*in_data)(idx, 1);
        mat(2, 2) = (*in_data)(idx, 2);
        mat(1, 2) = mat(2, 1) = (*in_data)(idx, 3);
        mat(0, 2) = mat(2, 0) = (*in_data)(idx, 4);
        mat(0, 1) = mat(1, 0) = (*in_data)(idx, 5);
#endif
        mat = convertToConformation(mat);
#if (NDIM == 2)
        double Qxx = mat(0, 0);
        double Qyy = mat(1, 1);
        double Qxy = mat(0, 1);
        (*ret_data)(idx, 0) =
            l_inv * (-1.0 * (d_alpha * (Qxx * Qxx + Qxy * Qxy) + (1.0 - 2.0 * d_alpha) * Qxx + (d_alpha - 1.0)));
        (*ret_data)(idx, 1) =
            l_inv * (-1.0 * (d_alpha * (Qyy * Qyy + Qxy * Qxy) + (1.0 - 2.0 * d_alpha) * Qyy + (d_alpha - 1.0)));
        (*ret_data)(idx, 2) = l_inv * (-1.0 * (d_alpha * (Qxx * Qxy + Qxy * Qyy) + (1.0 - 2.0 * d_alpha) * Qxy));
#endif
#if (NDIM == 3)
        double Qxx = mat(0, 0);
        double Qyy = mat(1, 1);
        double Qzz = mat(2, 2);
        double Qxy = mat(0, 1);
        double Qxz = mat(0, 2);
        double Qyz = mat(1, 2);
        (*ret_data)(idx, 0) = l_inv * (1.0 - Qxx - d_alpha * ((-1.0 + Qxx) * (-1.0 + Qxx) + Qxy * Qxy + Qxz * Qxz));
        (*ret_data)(idx, 1) = l_inv * (1.0 - Qyy - d_alpha * ((-1.0 + Qyy) * (-1.0 + Qyy) + Qxy * Qxy + Qyz * Qyz));
        (*ret_data)(idx, 2) = l_inv * (1.0 - Qzz - d_alpha * ((-1.0 + Qzz) * (-1.0 + Qzz) + Qxz * Qxz + Qyz * Qyz));
        (*ret_data)(idx, 3) = l_inv * (-Qyz - d_alpha * (Qxy * Qxz + (-1.0 + Qyy) * Qyz + Qyz * (-1.0 + Qzz)));
        (*ret_data)(idx, 4) = l_inv * (-Qxz - d_alpha * ((-1.0 + Qxx) * Qxz + Qxz * Qyz + Qxz * (-1.0 + Qzz)));
        (*ret_data)(idx, 5) = l_inv * (-Qxy - d_alpha * ((-1.0 + Qxx) * Qxy + Qxy * (-1.0 + Qyy) + Qxz * Qyz));
#endif
    }
} // setDataOnPatch

} // namespace IBAMR
