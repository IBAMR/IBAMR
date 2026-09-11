// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIndex.h>
#include <Patch.h>

#include <cmath>

#include "PointwiseLevelSet.h"

#include <ibamr/app_namespaces.h>

namespace MultiphaseExamples
{
PointwiseLevelSet::PointwiseLevelSet(const std::string& object_name) : CartGridFunction(object_name)
{
}

void
PointwiseLevelSet::setDataOnPatch(const int data_idx,
                                  Pointer<Variable<NDIM>> /*var*/,
                                  Pointer<Patch<NDIM>> patch,
                                  const double data_time,
                                  bool /*initial_time*/,
                                  Pointer<PatchLevel<NDIM>> /*patch_level*/)
{
    const Box<NDIM>& box = patch->getBox();
    Pointer<CellData<NDIM, double>> data = patch->getPatchData(data_idx);
    TBOX_ASSERT(data && data->getDepth() == 1);
    Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
    const double* X_lower = geometry->getXLower();
    const double* dx = geometry->getDx();
    for (Box<NDIM>::Iterator i(box); i; i++)
    {
        IBTK::Vector X = IBTK::Vector::Zero();
        for (int d = 0; d < NDIM; ++d)
        {
            X[d] = X_lower[d] + dx[d] * (static_cast<double>(i()(d) - box.lower()(d)) + 0.5);
        }
        (*data)(CellIndex<NDIM>(i())) = evaluateSignedDistance(X, data_time);
    }
}

SphereLevelSet::SphereLevelSet(const std::string& object_name, const IBTK::Vector& center, const double radius)
    : PointwiseLevelSet(object_name), d_center(center), d_radius(radius)
{
}

bool
SphereLevelSet::isTimeDependent() const
{
    return false;
}

double
SphereLevelSet::evaluateSignedDistance(const IBTK::Vector& X, double /*time*/) const
{
    return std::sqrt(std::pow(X[0] - d_center[0], 2.0) + std::pow(X[1] - d_center[1], 2.0)
#if (NDIM == 3)
                     + std::pow(X[2] - d_center[2], 2.0)
#endif
                         ) -
           d_radius;
}

PlaneLevelSet::PlaneLevelSet(const std::string& object_name, const int axis, const double height)
    : PointwiseLevelSet(object_name), d_axis(axis), d_height(height)
{
    TBOX_ASSERT(0 <= axis && axis < NDIM);
}

bool
PlaneLevelSet::isTimeDependent() const
{
    return false;
}

double
PlaneLevelSet::evaluateSignedDistance(const IBTK::Vector& X, double /*time*/) const
{
    return X[d_axis] - d_height;
}
} // namespace MultiphaseExamples
