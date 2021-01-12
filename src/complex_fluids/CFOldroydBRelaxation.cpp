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

#include "ibamr/CFOldroydBRelaxation.h"

#include "ibtk/ibtk_utilities.h"

#include "CellData.h"
#include "CellIterator.h"
#include "Patch.h"
#include "tbox/Database.h"

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Box;
template <int DIM>
class Variable;
} // namespace hier
namespace pdat
{
template <int DIM>
class CellIndex;
} // namespace pdat
} // namespace SAMRAI

// Namespace
namespace IBAMR
{
CFOldroydBRelaxation::CFOldroydBRelaxation(const std::string& object_name, Pointer<Database> input_db)
    : CFRelaxationOperator(object_name, input_db)
{
    d_lambda = input_db->getDouble("relaxation_time");
    return;
} // Constructor

void
CFOldroydBRelaxation::setDataOnPatch(const int data_idx,
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
        (*ret_data)(idx, 0) = l_inv * (1.0 - mat(0, 0));
        (*ret_data)(idx, 1) = l_inv * (1.0 - mat(1, 1));
        (*ret_data)(idx, 2) = l_inv * (-mat(1, 0));
#endif
#if (NDIM == 3)
        (*ret_data)(idx, 0) = l_inv * (1.0 - mat(0, 0));
        (*ret_data)(idx, 1) = l_inv * (1.0 - mat(1, 1));
        (*ret_data)(idx, 2) = l_inv * (1.0 - mat(2, 2));
        (*ret_data)(idx, 3) = l_inv * (-mat(1, 2));
        (*ret_data)(idx, 4) = l_inv * (-mat(0, 2));
        (*ret_data)(idx, 5) = l_inv * (-mat(0, 1));
#endif
    }
} // setDataOnPatch

} // namespace IBAMR
