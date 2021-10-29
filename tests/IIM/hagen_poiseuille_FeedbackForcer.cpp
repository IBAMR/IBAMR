// Filename: hagen_poiseuille_FeedbackForcer.cpp
// Created by IBAMR developers - 2015-2021

#include "hagen_poiseuille_FeedbackForcer.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <SideData.h>
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline double smooth_kernel(const double r)
{
    return std::abs(r) < 1.0 ? 0.5 * (cos(M_PI * r) + 1.0) : 0.0;
} // smooth_kernel
}

////////////////////////////// PUBLIC ///////////////////////////////////////

hagen_poiseuille_FeedbackForcer::hagen_poiseuille_FeedbackForcer(const double height,
                               const double diameter,
                               const INSHierarchyIntegrator* fluid_solver,
                               const Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
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

bool hagen_poiseuille_FeedbackForcer::isTimeDependent() const
{
    return true;
} // isTimeDependent

void hagen_poiseuille_FeedbackForcer::setDataOnPatch(const int data_idx,
                                    Pointer<Variable<NDIM> > /*var*/,
                                    Pointer<Patch<NDIM> > patch,
                                    const double /*data_time*/,
                                    const bool initial_time,
                                    Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    Pointer<SideData<NDIM, double> > F_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(F_data);
#endif
    F_data->fillAll(0.0);
    if (initial_time) return;
    const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
    const double dt = d_fluid_solver->getCurrentTimeStepSize();
    const double rho = d_fluid_solver->getStokesSpecifications()->getRho();
    const double kappa = cycle_num >= 0 ? 0.25 * rho / dt : 0.0;
    Pointer<SideData<NDIM, double> > U_current_data =
        patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getCurrentContext());
    Pointer<SideData<NDIM, double> > U_new_data =
        patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getNewContext());
#if !defined(NDEBUG)
    TBOX_ASSERT(U_current_data);
#endif
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = d_patch_hierarchy->getGridGeometry();
    const double* const dx_coarsest = grid_geometry->getDx();
    double dx_finest[NDIM];
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    const IntVector<NDIM>& finest_ratio = d_patch_hierarchy->getPatchLevel(finest_ln)->getRatio();
    for (int d = 0; d < NDIM; ++d)
    {
        dx_finest[d] = dx_coarsest[d] / static_cast<double>(finest_ratio(d));
    }
    const Box<NDIM> domain_box = Box<NDIM>::refine(grid_geometry->getPhysicalDomain()[0], pgeom->getRatio());

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
            Box<NDIM> bdry_box = domain_box;
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
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box, component)); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);
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
                        if ((X[1] > 0.5*H - 0.5 * D) && (X[1] < 0.5 * H + 0.5 * D))
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
