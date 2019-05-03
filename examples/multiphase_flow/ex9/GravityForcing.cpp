// Filename: GravityForcing.cpp
// Created on 13 Nov 2017 by Nishant Nangia
//
// Copyright (c) 2002-2019, Amneet Bhalla and Nishant Nangia
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

#include "GravityForcing.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <IBAMR_config.h>
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
                                        Pointer<Variable<NDIM> > /*var*/,
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
                               Pointer<Variable<NDIM> > /*var*/,
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
