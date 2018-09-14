// Filename: SetDampingZone.cpp
// Created on Sep 04, 2018 by Amneet Bhalla

// APPLICATION INCLUDES
#include "SetDampingZone.h"

#include <CartesianGridGeometry.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/HierarchyMathOps.h>

// C++ INCLUDES

/////////////////////////////// STATIC ///////////////////////////////////////

void
callDampingZoneCallbackFunction(int damping_coef_idx,
                                Pointer<HierarchyMathOps> hier_math_ops,
                                int cycle_num,
                                double time,
                                double current_time,
                                double new_time,
                                void* ctx)
{
    static DampingZone* ptr_DampingZone = static_cast<DampingZone*>(ctx);
    ptr_DampingZone->setDampingZoneCoefPatchData(
        damping_coef_idx, hier_math_ops, cycle_num, time, current_time, new_time);
    return;

} // callDampingZoneCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////

DampingZone::DampingZone(const std::string& object_name,
                         Pointer<INSVCStaggeredHierarchyIntegrator> ins_hier_integrator,
                         Pointer<Database> input_db)
    : d_object_name(object_name), d_ins_hier_integrator(ins_hier_integrator)
{
    d_x_zone_start = input_db->getDouble("x_zone_start");
    d_x_zone_end = input_db->getDouble("x_zone_end");
    d_theta = input_db->getDouble("theta");
    return;
} // DampingZone

DampingZone::~DampingZone()
{
    // intentionally left blank
    return;

} //~DampingZone

void
DampingZone::setDampingZoneCoefPatchData(int damping_coef_idx,
                                         Pointer<HierarchyMathOps> hier_math_ops,
                                         const int cycle_num,
                                         const double time,
                                         const double current_time,
                                         const double new_time)
{
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();

    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // Get interpolated density variable
    const int rho_ins_idx = d_ins_hier_integrator->getLinearOperatorRhoPatchDataIndex();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double* const patch_x_lower = patch_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();

            Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_ins_idx);
            Pointer<SideData<NDIM, double> > L_data = patch->getPatchData(damping_coef_idx);

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    Index<NDIM> i = it();
                    SideIndex<NDIM> i_side(i, axis, 0);
                    const double x_posn =
                        patch_x_lower[0] +
                        patch_dx[0] * (static_cast<double>(i(0) - patch_lower(0)) + axis == 0 ? 0.0 : 0.5);

                    if (x_posn >= d_x_zone_start)
                    {
                        (*L_data)(i_side, 0) = (*rho_data)(i_side, 0) * d_theta * (x_posn - d_x_zone_start) /
                                               (d_x_zone_end - d_x_zone_start);
                    }
                    else
                    {
                        (*L_data)(i_side, 0) = 0.0;
                    }
                }
            }
        }
    }
    return;
} // setDampingZoneCoefPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
