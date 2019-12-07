// Filename: CFRelaxationOperator.cpp
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

#include "ibamr/CFRelaxationOperator.h"
#include "ibamr/namespaces.h"

// Namespace
namespace IBAMR
{
CFRelaxationOperator::CFRelaxationOperator(const std::string& object_name, Pointer<Database> input_db)
    : CartGridFunction(object_name)
{
    if (input_db)
    {
        d_evolve_type = string_to_enum<TensorEvolutionType>(input_db->getStringWithDefault("evolve_type", "STANDARD"));
    }
    return;
} // Constructor

void
CFRelaxationOperator::setPatchDataIndex(int data_index)
{
    d_W_cc_idx = data_index;
    return;
} // setPatchDataIndex

bool
CFRelaxationOperator::isTimeDependent() const
{
    return true;
} // isTimeDependent

MatrixNd
CFRelaxationOperator::convertToConformation(const MatrixNd& mat)
{
    switch (d_evolve_type)
    {
    case SQUARE_ROOT:
        return mat * mat;
        break;
    case LOGARITHM:
        return mat.exp();
        break;
    case STANDARD:
        return mat;
        break;
    case UNKNOWN_TENSOR_EVOLUTION_TYPE:
        TBOX_ERROR(d_object_name << ":\n"
                                 << "  Uknown tensor evolution type.");
        return mat;
        break;
    default:
        TBOX_ERROR("Should not reach this statement.");
        break;
    }
    return mat;
}

} // namespace IBAMR
