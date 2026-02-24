// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "ibtk/samrai_compatibility_names.h"
// SAMRAI INCLUDES
#include "SAMRAIBox.h"
#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAICartesianPatchGeometry.h"
#include "SAMRAIIndex.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAISideData.h"
#include "SAMRAISideGeometry.h"
#include "SAMRAISideIndex.h"
#include "SAMRAIUtilities.h"
#include "SAMRAIVariable.h"

namespace
{
inline double
smooth_kernel(const double r)
{
    return std::abs(r) < 1.0 ? 0.5 * (cos(M_PI * r) + 1.0) : 0.0;
} // smooth_kernel
} // namespace

////////////////////////////// PUBLIC ///////////////////////////////////////

hagen_poiseuille_FeedbackForcer::hagen_poiseuille_FeedbackForcer(
    const double height,
    const double diameter,
    const INSHierarchyIntegrator* fluid_solver,
    const SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy)
    : d_H(height), d_D(diameter), d_fluid_solver(fluid_solver), d_patch_hierarchy(patch_hierarchy)
{
    // intentionally blank
    return;
} // hagen_poiseuille_FeedbackForcer

hagen_poiseuille_FeedbackForcer::~hagen_poiseuille_FeedbackForcer()
{
    // intentionally blank
    return;
} // ~hagen_poiseuille_FeedbackForcer

bool
hagen_poiseuille_FeedbackForcer::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
hagen_poiseuille_FeedbackForcer::setDataOnPatch(const int data_idx,
                                                SAMRAIPointer<SAMRAIVariable> /*var*/,
                                                SAMRAIPointer<SAMRAIPatch> patch,
                                                const double /*data_time*/,
                                                const bool initial_time,
                                                SAMRAIPointer<SAMRAIPatchLevel> /*patch_level*/)
{
    SAMRAIPointer<SAMRAISideData<double>> F_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(F_data);
#endif
    F_data->fillAll(0.0);
    if (initial_time) return;
    const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
    const double dt = d_fluid_solver->getCurrentTimeStepSize();
    const double rho = d_fluid_solver->getStokesSpecifications()->getRho();
    const double kappa = cycle_num >= 0 ? 0.25 * rho / dt : 0.0;
    SAMRAIPointer<SAMRAISideData<double>> U_current_data =
        patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getCurrentContext());
    SAMRAIPointer<SAMRAISideData<double>> U_new_data =
        patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getNewContext());
#if !defined(NDEBUG)
    TBOX_ASSERT(U_current_data);
#endif
    const SAMRAIBox& patch_box = patch->getBox();
    SAMRAIPointer<SAMRAICartesianPatchGeometry> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geometry = d_patch_hierarchy->getGridGeometry();
    const double* const dx_coarsest = grid_geometry->getDx();
    double dx_finest[NDIM];
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    const SAMRAIIntVector& finest_ratio = d_patch_hierarchy->getPatchLevel(finest_ln)->getRatio();
    for (int d = 0; d < NDIM; ++d)
    {
        dx_finest[d] = dx_coarsest[d] / static_cast<double>(finest_ratio(d));
    }
    const SAMRAIBox domain_box = SAMRAIBox::refine(grid_geometry->getPhysicalDomain()[0], pgeom->getRatio());

    // Clamp the velocity along the top/bottom of the domain (but not in the
    // interior of the structure).
    static const int axis = 0;
    const double L = 4.0 * dx_finest[axis];
    const int offset = static_cast<int>(L / dx[axis]);
    for (int side = 0; side <= 1; ++side)
    {
        const bool is_lower = side == 0;
        if (pgeom->getTouchesRegularBoundary(axis, side))
        {
            SAMRAIBox bdry_box = domain_box;
            if (side == 0)
            {
                bdry_box.upper(axis) = domain_box.lower(axis) + offset;
            }
            else
            {
                bdry_box.lower(axis) = domain_box.upper(axis) - offset;
            }
            bdry_box = bdry_box * patch_box;
            for (int component = 0; component < NDIM; ++component)
            {
                for (SAMRAIBox::Iterator b(SAMRAISideGeometry::toSideBox(bdry_box, component)); b; b++)
                {
                    const SAMRAIIndex& i = b();
                    const SAMRAISideIndex i_s(i, component, SAMRAISideIndex::Lower);
                    const double U_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                    const double U_new = U_new_data ? (*U_new_data)(i_s) : 0.0;
                    const double U = (cycle_num > 0) ? 0.5 * (U_new + U_current) : U_current;

                    double X[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == component ? 0.0 : 0.5));
                    }
                    double fac;
                    if (component == axis)
                    {
                        const double H = d_H;
                        const double D = d_D;
                        if ((X[1] > 0.5 * H - 0.5 * D) && (X[1] < 0.5 * H + 0.5 * D))
                        {
                            fac = 0.0;
                        }
                        else
                        {
                            fac = 1.0;
                        }
                    }
                    else
                    {
                        fac = 1.0;
                    }
                    const double x_bdry = (is_lower ? x_lower[axis] : x_upper[axis]);
                    fac *= smooth_kernel((X[axis] - x_bdry) / L);
                    (*F_data)(i_s) = -fac * kappa * U;
                }
            }
        }
    }
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
