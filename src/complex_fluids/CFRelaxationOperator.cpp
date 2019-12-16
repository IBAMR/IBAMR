// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

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
