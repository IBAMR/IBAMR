// Filename: SurfaceTensionForceFunction.cpp
// Created on 16 Dec 2017 by Amneet Bhalla
//
// Copyright (c) 2002-2017, Amneet Bhalla
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

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "IBAMR_config.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchData.h"
#include "SideData.h"
#include "Variable.h"
#include "VariableContext.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/SurfaceTensionForceFunction.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartGridFunction.h"
#include "ibtk/ibtk_utilities.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
} // namespace hier
} // namespace SAMRAI

// FORTRAN ROUTINES
#if (NDIM == 2)
#define MOLLIFY_IB_4_FC IBAMR_FC_FUNC(mollify_ib_4_2d, MOLLIFY_IB_4_2D)
#define SC_NORMAL_FC IBAMR_FC_FUNC(sc_normal_2d, SC_NORMAL_2D)
#define CC_CURVATURE_FC IBAMR_FC_FUNC(cc_curvature_2d, CC_CURVATURE_2D)
#define SC_SURFACE_TENSION_FORCE_FC IBAMR_FC_FUNC(sc_surface_tension_force_2d, SC_SURFACE_TENSION_FORCE_2D)
#endif

#if (NDIM == 3)
#define MOLLIFY_IB_4_FC IBAMR_FC_FUNC(mollify_ib_4_3d, MOLLIFY_IB_4_3D)
#define SC_NORMAL_FC IBAMR_FC_FUNC(sc_normal_3d, SC_NORMAL_3D)
#define CC_CURVATURE_FC IBAMR_FC_FUNC(cc_curvature_3d, CC_CURVATURE_3D)
#define SC_SURFACE_TENSION_FORCE_FC IBAMR_FC_FUNC(sc_surface_tension_force_3d, SC_SURFACE_TENSION_FORCE_3D)
#endif

extern "C" {
void MOLLIFY_IB_4_FC(double* V,
                     const int& V_gcw,
                     const double* U,
                     const int& U_gcw,
                     const int& ilower0,
                     const int& iupper0,
                     const int& ilower1,
                     const int& iupper1
#if (NDIM == 3)
                     ,
                     const int& ilower2,
                     const int& iupper2
#endif
                     );

void SC_NORMAL_FC(double* N00,
                  double* N01,
#if (NDIM == 3)
                  double* N02,
#endif
                  double* N10,
                  double* N11,
#if (NDIM == 3)
                  double* N12,
                  double* N20,
                  double* N21,
                  double* N22,
#endif
                  const int& N_gcw,
                  const double* U,
                  const int& U_gcw,
                  const int& ilower0,
                  const int& iupper0,
                  const int& ilower1,
                  const int& iupper1,
#if (NDIM == 3)
                  const int& ilower2,
                  const int& iupper2,
#endif
                  const double* dx);

void CC_CURVATURE_FC(double* K,
                     const int& K_gcw,
                     const double* N00,
                     const double* N01,
#if (NDIM == 3)
                     const double* N02,
#endif
                     const double* N10,
                     const double* N11,
#if (NDIM == 3)
                     const double* N12,
                     const double* N20,
                     const double* N21,
                     const double* N22,
#endif
                     const int& N_gcw,
                     const int& ilower0,
                     const int& iupper0,
                     const int& ilower1,
                     const int& iupper1,
#if (NDIM == 3)
                     const int& ilower2,
                     const int& iupper2,
#endif
                     const double* dx);

void SC_SURFACE_TENSION_FORCE_FC(double* F0,
                                 double* F1,
#if (NDIM == 3)
                                 double* F2,
#endif
                                 const int& F_gcw,
                                 const double* K,
                                 const int& K_gcw,
                                 const double* N00,
                                 const double* N11,
#if (NDIM == 3)
                                 const double* N22,
#endif
                                 const int& N_gcw,
                                 const int& ilower0,
                                 const int& iupper0,
                                 const int& ilower1,
                                 const int& iupper1,
#if (NDIM == 3)
                                 const int& ilower2,
                                 const int& iupper2,
#endif
                                 const double& sigma);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline double
smooth_kernel(const double r)
{
    return std::abs(r) < 1.0 ? 0.5 * (cos(M_PI * r) + 1.0) : 0.0;
} // smooth_kernel
}

////////////////////////////// PUBLIC ///////////////////////////////////////

SurfaceTensionForceFunction::SurfaceTensionForceFunction(const std::string& object_name,
                                                         const Pointer<Database> input_db,
                                                         const INSHierarchyIntegrator* fluid_solver,
                                                         Pointer<CartesianGridGeometry<NDIM> > grid_geometry)
    : CartGridFunction(object_name), d_fluid_solver(fluid_solver), d_grid_geometry(grid_geometry)
{
    // Set some default values
    d_phi_idx = -1;
    d_kernel_fcn = "none";
    d_sigma = 1.0;

    if (input_db)
    {
        d_kernel_fcn = input_db->getStringWithDefault("kernel", d_kernel_fcn);
        d_kernel_fcn = input_db->getStringWithDefault("smoother", d_kernel_fcn);
        d_kernel_fcn = input_db->getStringWithDefault("kernel_fcn", d_kernel_fcn);
        d_kernel_fcn = input_db->getStringWithDefault("smoother_fcn", d_kernel_fcn);

        d_sigma = input_db->getDoubleWithDefault("sigma", d_sigma);
        d_sigma = input_db->getDoubleWithDefault("surface_tension_coef", d_sigma);
    }
    return;
} // SurfaceTensionForceFunction

SurfaceTensionForceFunction::~SurfaceTensionForceFunction()
{
    // intentionally blank
    return;
} // ~SurfaceTensionForceFunction

void
SurfaceTensionForceFunction::setIndicatorPatchDataIndex(int phi_idx)
{
    d_phi_idx = phi_idx;
    return;
} // setIndicatorPatchDataIndex

void
SurfaceTensionForceFunction::setSmoother(const std::string& kernel_fcn)
{
    d_kernel_fcn = kernel_fcn;
    return;
} // setSmoother

void
SurfaceTensionForceFunction::setSurfaceTensionCoef(double sigma)
{
    d_sigma = sigma;
    return;
} // setSurfaceTensionCoef

bool
SurfaceTensionForceFunction::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
SurfaceTensionForceFunction::setDataOnPatchHierarchy(const int data_idx,
                                                     Pointer<Variable<NDIM> > var,
                                                     Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                     const double data_time,
                                                     const bool initial_time,
                                                     const int coarsest_ln_in,
                                                     const int finest_ln_in)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif

    // Allocate data for smooth indicator function.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > phi_var;
    var_db->mapIndexToVariable(d_phi_idx, phi_var);
    Pointer<CellVariable<NDIM, double> > phi_cc_var = phi_var;
#if !defined(NDEBUG)
    TBOX_ASSERT(!phi_var.isNull());
#endif

    IntVector<NDIM> cell_ghosts = getMinimumGhostWidth(d_kernel_fcn);
    d_smooth_phi_idx =
        var_db->registerVariableAndContext(phi_cc_var, var_db->getContext(d_object_name + "::SMOOTH"), cell_ghosts);
    const int scratch_phi_idx =
        var_db->registerVariableAndContext(phi_cc_var, var_db->getContext(d_object_name + "::SCRATCH"), cell_ghosts);

    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_smooth_phi_idx, data_time);
        hierarchy->getPatchLevel(ln)->allocatePatchData(scratch_phi_idx, data_time);
    }

    // Fill smooth phi with given phi.
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);
    hier_cc_data_ops.copyData(d_phi_idx, scratch_phi_idx, /*interior_only*/ true);

    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent scratch_phi_transaction(
        scratch_phi_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "QUADRATIC", false, NULL);
    Pointer<HierarchyGhostCellInterpolation> fill_op = new HierarchyGhostCellInterpolation();
    fill_op->initializeOperatorState(scratch_phi_transaction, hierarchy, coarsest_ln, finest_ln);
    fill_op->fillData(data_time);

    // Mollify phi.
    if (d_kernel_fcn == "none")
    {
        hier_cc_data_ops.copyData(scratch_phi_idx, d_smooth_phi_idx, /*interior_only*/ false);
    }
    else
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());

                Pointer<CellData<NDIM, double> > smooth_phi_data = patch->getPatchData(d_smooth_phi_idx);
                Pointer<CellData<NDIM, double> > scratch_phi_data = patch->getPatchData(scratch_phi_idx);
                double* const V = smooth_phi_data->getPointer(0);
                const double* const U = smooth_phi_data->getPointer(0);
                const int V_gcw = (smooth_phi_data->getGhostCellWidth()).max();
                const int U_gcw = (scratch_phi_data->getGhostCellWidth()).max();

                const Box<NDIM>& patch_box = patch->getBox();

                if (d_kernel_fcn == "IB_4")
                {
                    MOLLIFY_IB_4_FC(V,
                                    V_gcw,
                                    U,
                                    U_gcw,
                                    patch_box.lower(0),
                                    patch_box.upper(0),
                                    patch_box.lower(1),
                                    patch_box.upper(1)
#if (NDIM == 3)
                                        ,
                                    patch_box.lower(2),
                                    patch_box.upper(2)
#endif
                                        );
                }
            }
        }
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        InterpolationTransactionComponent smooth_phi_transaction(
            d_smooth_phi_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "QUADRATIC", false, NULL);
        fill_op->resetTransactionComponent(smooth_phi_transaction);
        fill_op->fillData(data_time);
    }

    // Fill data on each patch level
    CartGridFunction::setDataOnPatchHierarchy(
        data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);

    // Deallocate and remove scratch/smooth phi.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(scratch_phi_idx);
        hierarchy->getPatchLevel(ln)->deallocatePatchData(d_smooth_phi_idx);
    }
    var_db->removePatchDataIndex(scratch_phi_idx);
    var_db->removePatchDataIndex(d_smooth_phi_idx);

    return;
} // setDataOnPatchHierarchy

void
SurfaceTensionForceFunction::setDataOnPatch(const int data_idx,
                                            Pointer<Variable<NDIM> > /*var*/,
                                            Pointer<Patch<NDIM> > patch,
                                            const double data_time,
                                            const bool initial_time,
                                            Pointer<PatchLevel<NDIM> > level)
{
    Pointer<PatchData<NDIM> > f_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(f_data);
#endif
    Pointer<CellData<NDIM, double> > f_cc_data = f_data;
    Pointer<SideData<NDIM, double> > f_sc_data = f_data;
#if !defined(NDEBUG)
    TBOX_ASSERT(f_cc_data || f_sc_data);
#endif
    if (f_cc_data) f_cc_data->fillAll(0.0);
    if (f_sc_data) f_sc_data->fillAll(0.0);

    if (f_cc_data) setDataOnPatchCell(f_cc_data, patch, data_time, initial_time, level);
    if (f_sc_data) setDataOnPatchSide(f_sc_data, patch, data_time, initial_time, level);
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
SurfaceTensionForceFunction::setDataOnPatchCell(Pointer<CellData<NDIM, double> > /*F_data*/,
                                                Pointer<Patch<NDIM> > /*patch*/,
                                                const double /*data_time*/,
                                                const bool /*initial_time*/,
                                                Pointer<PatchLevel<NDIM> > /*level*/)
{
    TBOX_ERROR(
        "SurfaceTensionForceFunction::setDataOnPatchCell() Cell centered surface force tension is not implemented yet."
        << std::endl);

    return;
} // setDataOnPatchCell

void
SurfaceTensionForceFunction::setDataOnPatchSide(Pointer<SideData<NDIM, double> > F_data,
                                                Pointer<Patch<NDIM> > patch,
                                                const double /*data_time*/,
                                                const bool /*initial_time*/,
                                                Pointer<PatchLevel<NDIM> > /*level*/)
{
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    // First find normal in terms of gradient of phi.
    // n = grad(phi)
    SideData<NDIM, double> N(patch_box, /*depth*/ NDIM, /*gcw*/ IntVector<NDIM>(2));
    Pointer<CellData<NDIM, double> > U = patch->getPatchData(d_smooth_phi_idx);

    SC_NORMAL_FC(N.getPointer(0, 0),
                 N.getPointer(0, 1),
#if (NDIM == 3)
                 N.getPointer(0, 2),
#endif
                 N.getPointer(1, 0),
                 N.getPointer(1, 1),
#if (NDIM == 3)
                 N.getPointer(1, 2),
                 N.getPointer(2, 0),
                 N.getPointer(2, 1),
                 N.getPointer(2, 2),
#endif
                 N.getGhostCellWidth().max(),
                 U->getPointer(),
                 U->getGhostCellWidth().max(),
                 patch_box.lower(0),
                 patch_box.upper(0),
                 patch_box.lower(1),
                 patch_box.upper(1),
#if (NDIM == 3)
                 patch_box.lower(2),
                 patch_box.upper(2),
#endif
                 dx);

    // Next find the cell centered curvature.
    // K = -div (n/|n|)
    CellData<NDIM, double> K(patch_box, /*depth*/ 1, /*gcw*/ IntVector<NDIM>(1));
    CC_CURVATURE_FC(K.getPointer(),
                    K.getGhostCellWidth().max(),
                    N.getPointer(0, 0),
                    N.getPointer(0, 1),
#if (NDIM == 3)
                    N.getPointer(0, 2),
#endif
                    N.getPointer(1, 0),
                    N.getPointer(1, 1),
#if (NDIM == 3)
                    N.getPointer(1, 2),
                    N.getPointer(2, 0),
                    N.getPointer(2, 1),
                    N.getPointer(2, 2),
#endif
                    N.getGhostCellWidth().max(),
                    patch_box.lower(0),
                    patch_box.upper(0),
                    patch_box.lower(1),
                    patch_box.upper(1),
#if (NDIM == 3)
                    patch_box.lower(2),
                    patch_box.upper(2),
#endif
                    dx);

    // Compute the surface tension force
    // F = sigma * K * grad(phi)
    SC_SURFACE_TENSION_FORCE_FC(F_data->getPointer(0),
                                F_data->getPointer(1),
#if (NDIM == 3)
                                F_data->getPointer(2),
#endif
                                F_data->getGhostCellWidth().max(),
                                K.getPointer(),
                                K.getGhostCellWidth().max(),
                                N.getPointer(0, 0),
                                N.getPointer(1, 1),
#if (NDIM == 3)
                                N.getPointer(2, 2),
#endif
                                N.getGhostCellWidth().max(),
                                patch_box.lower(0),
                                patch_box.upper(0),
                                patch_box.lower(1),
                                patch_box.upper(1),
#if (NDIM == 3)
                                patch_box.lower(2),
                                patch_box.upper(2),
#endif
                                d_sigma);

    return;
} // setDataOnPatchSide

int
SurfaceTensionForceFunction::getStencilSize(const std::string& kernel_fcn)
{
    if (kernel_fcn == "IB_4") return 4;
    TBOX_ERROR("SurfaceTensionForceFunction::getStencilSize()\n"
               << "  Unknown kernel function "
               << kernel_fcn
               << std::endl);
    return -1;

} // getStencilSize

int
SurfaceTensionForceFunction::getMinimumGhostWidth(const std::string& kernel_fcn)
{
    if (kernel_fcn == "none")
    {
        return 2;
    }
    else
    {
        return std::max(2, static_cast<int>(floor(0.5 * getStencilSize(kernel_fcn))) + 1);
    }
} // getMinimumGhostWidth

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
