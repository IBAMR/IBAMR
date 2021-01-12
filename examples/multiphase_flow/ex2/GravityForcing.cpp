// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "GravityForcing.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

GravityForcing::GravityForcing(const std::string& object_name,
                               Pointer<INSVCStaggeredHierarchyIntegrator> ins_hierarchy_integrator,
                               std::vector<double> grav_const)
    : d_object_name(object_name), d_ins_hierarchy_integrator(ins_hierarchy_integrator), d_grav_const(grav_const)
{
    // intentionally blank
    return;
} // GravityForcing

GravityForcing::~GravityForcing()
{
    // intentionally blank
    return;
} // ~GravityForcing

bool
GravityForcing::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
GravityForcing::setDataOnPatchHierarchy(const int data_idx,
                                        Pointer<hier::Variable<NDIM> > /*var*/,
                                        Pointer<PatchHierarchy<NDIM> > hierarchy,
                                        const double /*data_time*/,
                                        const bool /*initial_time*/,
                                        const int coarsest_ln_in,
                                        const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);

    // Get interpolated density variable
    const int rho_ins_idx = d_ins_hierarchy_integrator->getLinearOperatorRhoPatchDataIndex();

#if !defined(NDEBUG)
    TBOX_ASSERT(rho_ins_idx >= 0);
#endif
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& box = patch->getBox();
            Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(data_idx);
            const Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_ins_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(box, axis)); it; it++)
                {
                    SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);
                    (*f_data)(s_i) = ((*rho_data)(s_i)) * d_grav_const[axis];
                }
            }
        }
    }
    return;
} // setDataOnPatchHierarchy

void
GravityForcing::setDataOnPatch(const int data_idx,
                               Pointer<hier::Variable<NDIM> > /*var*/,
                               Pointer<Patch<NDIM> > patch,
                               const double /*data_time*/,
                               const bool initial_time,
                               Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    if (initial_time)
    {
        Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(data_idx);
        f_data->fillAll(0.0);
    }
    // Intentionally left blank

} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
