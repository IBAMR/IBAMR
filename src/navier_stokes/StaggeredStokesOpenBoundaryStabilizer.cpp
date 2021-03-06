// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/StaggeredStokesOpenBoundaryStabilizer.h"
#include "ibamr/StokesSpecifications.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "BoxArray.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "Variable.h"
#include "VariableContext.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <cmath>
#include <ostream>
#include <string>

#include "ibamr/namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline double
smooth_kernel(const double r)
{
    return std::abs(r) < 1.0 ? 0.5 * (std::cos(M_PI * r) + 1.0) : 0.0;
} // smooth_kernel
} // namespace

////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesOpenBoundaryStabilizer::StaggeredStokesOpenBoundaryStabilizer(
    const std::string& object_name,
    Pointer<Database> input_db,
    const INSHierarchyIntegrator* fluid_solver,
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry)
    : CartGridFunction(object_name),
      d_open_bdry(array_constant<bool, 2 * NDIM>(false)),
      d_inflow_bdry(array_constant<bool, 2 * NDIM>(false)),
      d_outflow_bdry(array_constant<bool, 2 * NDIM>(false)),
      d_width(array_constant<double, 2 * NDIM>(0.0)),
      d_fluid_solver(fluid_solver),
      d_grid_geometry(grid_geometry)
{
    if (input_db)
    {
        for (unsigned int location_index = 0; location_index < 2 * NDIM; ++location_index)
        {
            const std::string stabilization_type_key = "stabilization_type_" + std::to_string(location_index);
            if (input_db->keyExists(stabilization_type_key))
            {
                const std::string stabilization_type = input_db->getString(stabilization_type_key);
                if (stabilization_type == "INFLOW")
                {
                    d_open_bdry[location_index] = true;
                    d_inflow_bdry[location_index] = true;
                    d_outflow_bdry[location_index] = false;
                }
                else if (stabilization_type == "OUTFLOW")
                {
                    d_open_bdry[location_index] = true;
                    d_inflow_bdry[location_index] = false;
                    d_outflow_bdry[location_index] = true;
                }
                else if (stabilization_type != "NONE")
                {
                    TBOX_ERROR(
                        "StaggeredStokesOpenBoundaryStabilizer::"
                        "StaggeredStokesOpenBoundaryStabilizer():\n"
                        << "  unsupported stabilization type: ``" << stabilization_type << "''\n"
                        << "  supported values are: ``INFLOW'', ``OUTFLOW'', or "
                           "``NONE''\n");
                }
            }
            const std::string width_key = "width_" + std::to_string(location_index);
            if (input_db->keyExists(width_key))
            {
                d_width[location_index] = input_db->getDouble(width_key);
            }
        }
    }
    return;
} // StaggeredStokesOpenBoundaryStabilizer

bool
StaggeredStokesOpenBoundaryStabilizer::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
StaggeredStokesOpenBoundaryStabilizer::setDataOnPatch(const int data_idx,
                                                      Pointer<Variable<NDIM> > /*var*/,
                                                      Pointer<Patch<NDIM> > patch,
                                                      const double /*data_time*/,
                                                      const bool initial_time,
                                                      Pointer<PatchLevel<NDIM> > /*level*/)
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
    const double kappa = cycle_num >= 0 ? 0.5 * rho / dt : 0.0;
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
    const IntVector<NDIM>& ratio = pgeom->getRatio();
    const Box<NDIM> domain_box = Box<NDIM>::refine(d_grid_geometry->getPhysicalDomain()[0], ratio);
    for (unsigned int location_index = 0; location_index < 2 * NDIM; ++location_index)
    {
        const unsigned int axis = location_index / 2;
        const unsigned int side = location_index % 2;
        const bool is_lower = side == 0;
        if (d_open_bdry[location_index] && pgeom->getTouchesRegularBoundary(axis, side))
        {
            Box<NDIM> bdry_box = domain_box;
            const int offset = static_cast<int>(d_width[location_index] / dx[axis]);
            if (is_lower)
            {
                bdry_box.upper(axis) = domain_box.lower(axis) + offset;
            }
            else
            {
                bdry_box.lower(axis) = domain_box.upper(axis) - offset;
            }
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box * patch_box, axis)); b; b++)
            {
                const hier::Index<NDIM>& i = b();
                const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                const double U_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                const double U_new = U_new_data ? (*U_new_data)(i_s) : 0.0;
                const double U = (cycle_num > 0) ? 0.5 * (U_new + U_current) : U_current;
                const double n = is_lower ? -1.0 : +1.0;
                if ((d_inflow_bdry[location_index] && U * n > 0.0) || (d_outflow_bdry[location_index] && U * n < 0.0))
                {
                    const double x = x_lower[axis] + dx[axis] * static_cast<double>(i(axis) - patch_box.lower(axis));
                    const double x_bdry = (is_lower ? x_lower[axis] : x_upper[axis]);
                    (*F_data)(i_s) = smooth_kernel((x - x_bdry) / d_width[location_index]) * kappa * (0.0 - U);
                }
            }
        }
    }
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
