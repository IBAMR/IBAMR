// Filename: SpongeLayerForceFunction.cpp
// Created on 28 Oct 2011 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <cmath>
#include <iosfwd>
#include <ostream>
#include <string>

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxArray.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "boost/array.hpp"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/SpongeLayerForceFunction.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartGridFunction.h"
#include "ibtk/ibtk_utilities.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Database.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{

class PatchLevel;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline double smooth_kernel(const double r)
{
    return std::abs(r) < 1.0 ? 0.5 * (cos(M_PI * r) + 1.0) : 0.0;
} // smooth_kernel
}

////////////////////////////// PUBLIC ///////////////////////////////////////

SpongeLayerForceFunction::SpongeLayerForceFunction(const std::string& object_name,
                                                   const boost::shared_ptr<Database> input_db,
                                                   const INSHierarchyIntegrator* fluid_solver,
                                                   boost::shared_ptr<CartesianGridGeometry > grid_geometry)
    : CartGridFunction(object_name), d_forcing_enabled(array_constant<std::vector<bool>, 2 * NDIM>(std::vector<bool>(NDIM))),
      d_width(array_constant<double, 2 * NDIM>(0.0)), d_fluid_solver(fluid_solver), d_grid_geometry(grid_geometry)
{
    if (input_db)
    {
        for (unsigned int location_index = 0; location_index < 2 * NDIM; ++location_index)
        {
            for (unsigned int d = 0; d < NDIM; ++d) d_forcing_enabled[location_index][d] = false;
            std::ostringstream forcing_enabled_stream;
            forcing_enabled_stream << "forcing_enabled_" << location_index;
            const std::string forcing_enabled_key = forcing_enabled_stream.str();
            if (input_db->keyExists(forcing_enabled_key))
            {
                d_forcing_enabled[location_index] = input_db->getBoolArray(forcing_enabled_key);
            }
            std::ostringstream width_stream;
            width_stream << "width_" << location_index;
            const std::string width_key = width_stream.str();
            if (input_db->keyExists(width_key))
            {
                d_width[location_index] = input_db->getDouble(width_key);
            }
        }
    }
    return;
} // SpongeLayerForceFunction

SpongeLayerForceFunction::~SpongeLayerForceFunction()
{
    // intentionally blank
    return;
} // ~SpongeLayerForceFunction

bool SpongeLayerForceFunction::isTimeDependent() const
{
    return true;
} // isTimeDependent

void SpongeLayerForceFunction::setDataOnPatch(const int data_idx,
                                              boost::shared_ptr<Variable > /*var*/,
                                              boost::shared_ptr<Patch > patch,
                                              const double /*data_time*/,
                                              const bool initial_time,
                                              boost::shared_ptr<PatchLevel > /*level*/)
{
    boost::shared_ptr<PatchData > f_data = patch->getPatchData(data_idx);
    TBOX_ASSERT(f_data);
    boost::shared_ptr<CellData<double> > f_cc_data = f_data;
    boost::shared_ptr<SideData<double> > f_sc_data = f_data;
    TBOX_ASSERT(f_cc_data || f_sc_data);
    if (f_cc_data) f_cc_data->fillAll(0.0);
    if (f_sc_data) f_sc_data->fillAll(0.0);
    if (initial_time) return;
    const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
    const double dt = d_fluid_solver->getCurrentTimeStepSize();
    const double rho = d_fluid_solver->getStokesSpecifications()->getRho();
    const double kappa = cycle_num >= 0 ? 0.5 * rho / dt : 0.0;
    boost::shared_ptr<PatchData > u_current_data =
        patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getCurrentContext());
    boost::shared_ptr<PatchData > u_new_data =
        patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getNewContext());
    TBOX_ASSERT(u_current_data);
    if (f_cc_data) setDataOnPatchCell(f_data, u_current_data, u_new_data, kappa, patch);
    if (f_sc_data) setDataOnPatchSide(f_data, u_current_data, u_new_data, kappa, patch);
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void SpongeLayerForceFunction::setDataOnPatchCell(boost::shared_ptr<CellData<double> > F_data,
                                                  boost::shared_ptr<CellData<double> > U_current_data,
                                                  boost::shared_ptr<CellData<double> > U_new_data,
                                                  const double kappa,
                                                  boost::shared_ptr<Patch > patch)
{
    TBOX_ASSERT(F_data && U_current_data);
    const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
    const Box& patch_box = patch->getBox();
    boost::shared_ptr<CartesianPatchGeometry > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const IntVector& ratio = pgeom->getRatio();
    const Box domain_box = Box::refine(d_grid_geometry->getPhysicalDomain()[0], ratio);
    for (unsigned int location_index = 0; location_index < 2 * NDIM; ++location_index)
    {
        const unsigned int axis = location_index / 2;
        const unsigned int side = location_index % 2;
        const bool is_lower = side == 0;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (d_forcing_enabled[location_index][d] && pgeom->getTouchesRegularBoundary(axis, side))
            {
                Box bdry_box = domain_box;
                const int offset = static_cast<int>(d_width[location_index] / dx[axis]);
                if (is_lower)
                {
                    bdry_box.upper(axis) = domain_box.lower(axis) + offset;
                }
                else
                {
                    bdry_box.lower(axis) = domain_box.upper(axis) - offset;
                }
                for (CellIterator b(bdry_box * patch_box); b; b++)
                {
                    const CellIndex& i = *b;
                    const double U_current = U_current_data ? (*U_current_data)(i, d) : 0.0;
                    const double U_new = U_new_data ? (*U_new_data)(i, d) : 0.0;
                    const double U = (cycle_num > 0) ? 0.5 * (U_new + U_current) : U_current;
                    const double x =
                        x_lower[axis] + dx[axis] * (static_cast<double>(i(axis) - patch_box.lower(axis)) + 0.5);
                    const double x_bdry = (is_lower ? x_lower[axis] : x_upper[axis]);
                    (*F_data)(i, d) = smooth_kernel((x - x_bdry) / d_width[location_index]) * kappa * (0.0 - U);
                }
            }
        }
    }
    return;
} // setDataOnPatchCell

void SpongeLayerForceFunction::setDataOnPatchSide(boost::shared_ptr<SideData<double> > F_data,
                                                  boost::shared_ptr<SideData<double> > U_current_data,
                                                  boost::shared_ptr<SideData<double> > U_new_data,
                                                  const double kappa,
                                                  boost::shared_ptr<Patch > patch)
{
    TBOX_ASSERT(F_data && U_current_data);
    const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
    const Box& patch_box = patch->getBox();
    boost::shared_ptr<CartesianPatchGeometry > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    const IntVector& ratio = pgeom->getRatio();
    const Box domain_box = Box::refine(d_grid_geometry->getPhysicalDomain()[0], ratio);
    for (unsigned int location_index = 0; location_index < 2 * NDIM; ++location_index)
    {
        const unsigned int axis = location_index / 2;
        const unsigned int side = location_index % 2;
        const bool is_lower = side == 0;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (d_forcing_enabled[location_index][d] && pgeom->getTouchesRegularBoundary(axis, side))
            {
                Box bdry_box = domain_box;
                const int offset = static_cast<int>(d_width[location_index] / dx[axis]);
                if (is_lower)
                {
                    bdry_box.upper(axis) = domain_box.lower(axis) + offset;
                }
                else
                {
                    bdry_box.lower(axis) = domain_box.upper(axis) - offset;
                }
                for (Box::iterator b(SideGeometry::toSideBox(bdry_box * patch_box, d)); b; b++)
                {
                    const Index& i = *b;
                    const SideIndex i_s(i, d, SideIndex::Lower);
                    const double U_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                    const double U_new = U_new_data ? (*U_new_data)(i_s) : 0.0;
                    const double U = (cycle_num > 0) ? 0.5 * (U_new + U_current) : U_current;
                    const double x = x_lower[axis] + dx[axis] * static_cast<double>(i(axis) - patch_box.lower(axis));
                    const double x_bdry = (is_lower ? x_lower[axis] : x_upper[axis]);
                    (*F_data)(i_s) = smooth_kernel((x - x_bdry) / d_width[location_index]) * kappa * (0.0 - U);
                }
            }
        }
    }
    return;
} // setDataOnPatchSide

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
