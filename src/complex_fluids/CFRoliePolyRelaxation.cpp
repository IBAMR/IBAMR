// Filename: CFRoliePolyRelaxation.cpp
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

#include "ibamr/CFRoliePolyRelaxation.h"
#include "ibamr/namespaces.h"

// Namespace
namespace IBAMR
{
CFRoliePolyRelaxation::CFRoliePolyRelaxation(const std::string& object_name, Pointer<Database> input_db)
    : CFRelaxationOperator(object_name, input_db)
{
    d_lambda_d = input_db->getDouble("lambda_d");
    d_lambda_R = input_db->getDouble("lambda_R");
    d_beta = input_db->getDouble("beta");
    d_delta = input_db->getDouble("delta");
    return;
} // Constructor

void
CFRoliePolyRelaxation::setDataOnPatch(const int data_idx,
                                      Pointer<Variable<NDIM> > /*var*/,
                                      Pointer<Patch<NDIM> > patch,
                                      const double /*data_time*/,
                                      const bool initial_time,
                                      Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
    Pointer<CellData<NDIM, double> > ret_data = patch->getPatchData(data_idx);
    Pointer<CellData<NDIM, double> > in_data = patch->getPatchData(d_W_cc_idx);
    ret_data->fillAll(0.0);
    double tr = 0.0;
    if (initial_time) return;
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
        tr = Qxx + Qyy;
        (*ret_data)(idx, 0) =
            -1.0 / d_lambda_d * (Qxx - 1.0) -
            2.0 * (1.0 - sqrt(2.0 / tr)) / d_lambda_R * (Qxx + d_beta * pow(tr / 2.0, d_delta) * (Qxx - 1.0));
        (*ret_data)(idx, 1) =
            -1.0 / d_lambda_d * (Qyy - 1.0) -
            2.0 * (1.0 - sqrt(2.0 / tr)) / d_lambda_R * (Qyy + d_beta * pow(tr / 2.0, d_delta) * (Qyy - 1.0));
        (*ret_data)(idx, 2) = -1.0 / d_lambda_d * (Qxy)-2.0 * (1.0 - sqrt(2.0 / tr)) / d_lambda_R *
                              (Qxy + d_beta * pow(tr / 2.0, d_delta) * (Qxy));
#endif
#if (NDIM == 3)
        double Qxx = mat(0, 0);
        double Qyy = mat(1, 1);
        double Qzz = mat(2, 2);
        double Qxy = mat(0, 1);
        double Qxz = mat(0, 2);
        double Qyz = mat(1, 2);
        tr = Qxx + Qyy + Qzz;
        (*ret_data)(idx, 0) =
            -1.0 / d_lambda_d * (Qxx - 1.0) -
            2.0 * (1.0 - sqrt(3.0 / tr)) / d_lambda_R * (Qxx + d_beta * pow(tr / 3.0, d_delta) * (Qxx - 1.0));
        (*ret_data)(idx, 1) =
            -1.0 / d_lambda_d * (Qyy - 1.0) -
            2.0 * (1.0 - sqrt(3.0 / tr)) / d_lambda_R * (Qyy + d_beta * pow(tr / 3.0, d_delta) * (Qyy - 1.0));
        (*ret_data)(idx, 2) =
            -1.0 / d_lambda_d * (Qzz - 1.0) -
            2.0 * (1.0 - sqrt(3.0 / tr)) / d_lambda_R * (Qzz + d_beta * pow(tr / 3.0, d_delta) * (Qzz - 1.0));
        (*ret_data)(idx, 3) = -1.0 / d_lambda_d * (Qyz)-2.0 * (1.0 - sqrt(3.0 / tr)) / d_lambda_R *
                              (Qyz + d_beta * pow(tr / 3.0, d_delta) * (Qyz));
        (*ret_data)(idx, 4) = -1.0 / d_lambda_d * (Qxz)-2.0 * (1.0 - sqrt(3.0 / tr)) / d_lambda_R *
                              (Qxz + d_beta * pow(tr / 3.0, d_delta) * (Qxz));
        (*ret_data)(idx, 5) = -1.0 / d_lambda_d * (Qxy)-2.0 * (1.0 - sqrt(3.0 / tr)) / d_lambda_R *
                              (Qxy + d_beta * pow(tr / 3.0, d_delta) * (Qxy));
#endif
    }
} // setDataOnPatch

} // namespace IBAMR
