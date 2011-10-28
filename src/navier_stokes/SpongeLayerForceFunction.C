// Filename: SpongeLayerForceFunction.C
// Created on 28 Oct 2011 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "SpongeLayerForceFunction.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// SAMRAI INCLUDES
#include <CellData.h>
#include <SideData.h>

// C++ STDLIB INCLUDES

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

////////////////////////////// PUBLIC ///////////////////////////////////////

SpongeLayerForceFunction::SpongeLayerForceFunction(
    const std::string& object_name,
    const Pointer<Database> input_db,
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry)
    : CartGridFunction(object_name),
      d_grid_geometry(grid_geometry),
      d_u_var(NULL),
      d_u_ctx(NULL)
{
    if (input_db.isNull()) return;
    d_kappa = input_db->getDoubleWithDefault("kappa",0.0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            d_forcing_enabled[0][axis][d] = false;
            d_forcing_enabled[1][axis][d] = false;
        }

        std::ostringstream lower_forcing_stream;
        lower_forcing_stream << "lower_forcing_" << axis;
        const std::string lower_forcing_key = lower_forcing_stream.str();
        if (input_db->keyExists(lower_forcing_key))
        {
            tbox::Array<int> lower_forcing_comps = input_db->getIntegerArray(lower_forcing_key);
            for (int k = 0; k < lower_forcing_comps.size(); ++k)
            {
                d_forcing_enabled[0][axis][lower_forcing_comps[k]] = true;
            }
        }
        std::ostringstream lower_width_stream;
        lower_width_stream << "lower_width_" << axis;
        const std::string lower_width_key = lower_width_stream.str();
        d_width[0][axis] = input_db->getDoubleWithDefault(lower_width_key, 0.0);

        std::ostringstream upper_forcing_stream;
        upper_forcing_stream << "upper_forcing_" << axis;
        const std::string upper_forcing_key = upper_forcing_stream.str();
        if (input_db->keyExists(upper_forcing_key))
        {
            tbox::Array<int> upper_forcing_comps = input_db->getIntegerArray(upper_forcing_key);
            for (int k = 0; k < upper_forcing_comps.size(); ++k)
            {
                d_forcing_enabled[1][axis][upper_forcing_comps[k]] = true;
            }
        }
        std::ostringstream upper_width_stream;
        upper_width_stream << "upper_width_" << axis;
        const std::string upper_width_key = upper_width_stream.str();
        d_width[1][axis] = input_db->getDoubleWithDefault(upper_width_key, 0.0);
    }
    return;
}// SpongeLayerForceFunction

SpongeLayerForceFunction::~SpongeLayerForceFunction()
{
    // intentionally blank
    return;
}// ~SpongeLayerForceFunction

void
SpongeLayerForceFunction::setVelocityVariableAndContext(
    Pointer<Variable<NDIM> > u_var,
    Pointer<VariableContext> u_ctx)
{
    d_u_var = u_var;
    d_u_ctx = u_ctx;
    return;
}// setVelocityVariableAndContext

bool
SpongeLayerForceFunction::isTimeDependent() const
{
    return true;
}// isTimeDependent

void
SpongeLayerForceFunction::setDataOnPatch(
    const int data_idx,
    Pointer<Variable<NDIM> > /*var*/,
    Pointer<Patch<NDIM> > patch,
    const double /*data_time*/,
    const bool initial_time,
    Pointer<PatchLevel<NDIM> > /*level*/)
{
    Pointer<PatchData<NDIM> > f_data = patch->getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!f_data.isNull());
#endif
    Pointer<CellData<NDIM,double> > f_cc_data = f_data;
    Pointer<SideData<NDIM,double> > f_sc_data = f_data;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!f_cc_data.isNull() || !f_sc_data.isNull());
#endif
    if (!f_cc_data.isNull()) f_cc_data->fillAll(0.0);
    if (!f_sc_data.isNull()) f_sc_data->fillAll(0.0);
    if (initial_time) return;
    if (!f_cc_data.isNull()) setDataOnPatchCell(data_idx, patch);
    if (!f_sc_data.isNull()) setDataOnPatchSide(data_idx, patch);
    return;
}// setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
SpongeLayerForceFunction::setDataOnPatchCell(
    const int data_idx,
    Pointer<Patch<NDIM> > patch)
{
    Pointer<CellData<NDIM,double> > U_data = patch->getPatchData(d_u_var, d_u_ctx);
    Pointer<CellData<NDIM,double> > F_data = patch->getPatchData(data_idx);
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    blitz::TinyVector<int,NDIM> offset[2];
    for (unsigned int side = 0; side <= 1; ++side)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            offset[side][d] = static_cast<int>(d_width[side][d]/dx[d])-1;
        }
    }
    const IntVector<NDIM>& ratio = pgeom->getRatio();
    const Box<NDIM> domain_box = Box<NDIM>::refine(d_grid_geometry->getPhysicalDomain()[0],ratio);
    const Index<NDIM>& domain_lower = domain_box.lower();
    const Index<NDIM>& domain_upper = domain_box.upper();
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (unsigned int side = 0; side <= 1; ++side)
        {
            if (!pgeom->getTouchesRegularBoundary(axis,side)) continue;
            Box<NDIM> bdry_box = domain_box;
            if (side == 0)
            {
                bdry_box.upper()(axis) = domain_lower(axis)+offset[side][axis];
            }
            else
            {
                bdry_box.lower()(axis) = domain_upper(axis)-offset[side][axis];
            }
            for (unsigned int component = 0; component < NDIM; ++component)
            {
                if (!d_forcing_enabled[side][axis][component]) continue;
                for (Box<NDIM>::Iterator b(bdry_box*patch_box); b; b++)
                {
                    const Index<NDIM>& i = b();
                    (*F_data)(i,component) = d_kappa*(0.0 - (*U_data)(i,component));
                }
            }
        }
    }
    return;
}// setDataOnPatchCell

void
SpongeLayerForceFunction::setDataOnPatchSide(
    const int data_idx,
    Pointer<Patch<NDIM> > patch)
{
    Pointer<SideData<NDIM,double> > U_data = patch->getPatchData(d_u_var, d_u_ctx);
    Pointer<SideData<NDIM,double> > F_data = patch->getPatchData(data_idx);
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    blitz::TinyVector<int,NDIM> offset[2];
    for (unsigned int side = 0; side <= 1; ++side)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            offset[side][d] = static_cast<int>(d_width[side][d]/dx[d])-1;
        }
    }
    const IntVector<NDIM>& ratio = pgeom->getRatio();
    const Box<NDIM> domain_box = Box<NDIM>::refine(d_grid_geometry->getPhysicalDomain()[0],ratio);
    const Index<NDIM>& domain_lower = domain_box.lower();
    const Index<NDIM>& domain_upper = domain_box.upper();
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (unsigned int side = 0; side <= 1; ++side)
        {
            if (!pgeom->getTouchesRegularBoundary(axis,side)) continue;
            for (unsigned int component = 0; component < NDIM; ++component)
            {
                if (!d_forcing_enabled[side][axis][component]) continue;
                Box<NDIM> bdry_box = domain_box;
                if (side == 0)
                {
                    bdry_box.upper()(axis) = domain_lower(axis)+(offset[side][axis] + (axis==component ? 1 : 0));
                }
                else
                {
                    bdry_box.lower()(axis) = domain_upper(axis)-(offset[side][axis] + (axis==component ? 1 : 0));
                }
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box*patch_box,component)); b; b++)
                {
                    const Index<NDIM>& i = b();
                    const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);
                    (*F_data)(i_s) = d_kappa*(0.0 - (*U_data)(i_s));
                }
            }
        }
    }
    return;
}// setDataOnPatchSide

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
