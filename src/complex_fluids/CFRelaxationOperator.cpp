// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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

#include "ibtk/ibtk_utilities.h"

#include "tbox/Utilities.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include <unsupported/Eigen/MatrixFunctions>
IBTK_ENABLE_EXTRA_WARNINGS

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

// Namespace
namespace IBAMR
{
CFRelaxationOperator::CFRelaxationOperator(const std::string& object_name, Pointer<Database> input_db)
    : CartGridFunction(object_name)
{
    if (input_db)
    {
        if (input_db->keyExists("evolution_type"))
            d_evolve_type = string_to_enum<TensorEvolutionType>(input_db->getString("evolution_type"));
        if (input_db->keyExists("evolve_type"))
            d_evolve_type = string_to_enum<TensorEvolutionType>(input_db->getString("evolve_type"));
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
    case LOGARITHM:
        return mat.exp();
    case STANDARD:
        return mat;
    case UNKNOWN_TENSOR_EVOLUTION_TYPE:
        TBOX_ERROR(d_object_name << ":\n"
                                 << "  Uknown tensor evolution type.");
        return mat;
    default:
        TBOX_ERROR("Should not reach this statement.");
        break;
    }
    return mat;
}

} // namespace IBAMR
