// Filename: TurbulenceSSTKOmegaSourceFunction.cpp
// Created on 29 Sep 2019 by Ramakrishnan Thirumalaisamy
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
#include "IBAMR_config.h"

#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibamr/TurbulenceSSTKOmegaSourceFunction.h"
#include "ibamr/TwoEquationTurbulenceHierarchyIntegrator.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/CartGridFunction.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchData.h"
#include "SideData.h"
#include "Variable.h"
#include "VariableContext.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <cmath>
#include <iosfwd>
#include <ostream>
#include <string>

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
#define SST_K_EQN_BUOYANCY IBAMR_FC_FUNC(sst_k_eqn_buoyancy_2d, SST_K_EQN_BUOYANCY_2D)
#define SST_K_EQN_PRODUCTION IBAMR_FC_FUNC(sst_k_eqn_production_2d, SST_K_EQN_PRODUCTION_2D)
#define NEAR_WALL_PRODUCTION IBAMR_FC_FUNC(near_wall_production_2d, NEAR_WALL_PRODUCTION_2D)
#define SST_W_EQN_PRODUCTION IBAMR_FC_FUNC(sst_w_eqn_production_2d, SST_W_EQN_PRODUCTION_2D)
#define SST_W_EQN_CROSSDIFFUSION IBAMR_FC_FUNC(sst_w_eqn_crossdiffusion_2d, SST_W_EQN_CROSSDIFFUSION_2D)
#endif

#if (NDIM == 3)
#define SST_K_EQN_BUOYANCY IBAMR_FC_FUNC(sst_k_eqn_buoyancy_3d, SST_K_EQN_BUOYANCY_3D)
#define SST_K_EQN_PRODUCTION IBAMR_FC_FUNC(sst_k_eqn_production_3d, SST_K_EQN_PRODUCTION_3D)
#define NEAR_WALL_PRODUCTION IBAMR_FC_FUNC(near_wall_production_3d, NEAR_WALL_PRODUCTION_3D)
#define SST_W_EQN_PRODUCTION IBAMR_FC_FUNC(sst_w_eqn_production_3d, SST_W_EQN_PRODUCTION_3D)
#define SST_W_EQN_CROSSDIFFUSION IBAMR_FC_FUNC(sst_w_eqn_crossdiffusion_3d, SST_W_EQN_CROSSDIFFUSION_3D)
#endif

extern "C"
{
    void SST_K_EQN_PRODUCTION(const double*,
                              const int&,
                              const double*,
                              const int&,
                              const double*,
                              const int&,
                              const double*,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
#if (NDIM == 3)
                              const int&,
                              const int&,
#endif
                              const double&);

    void NEAR_WALL_PRODUCTION(const double*,
                              const int&,
                              const double*,
                              const int&,
                              const double*,
                              const int&,
                              const double*,
                              const int&,
#if (NDIM == 3)
                              const double*,
                              const int&,
#endif
                              const double&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
#if (NDIM == 3)
                              const int&,
                              const int&,
#endif
                              const int&,
                              const int&,
                              const int&,
                              const int&,
#if (NDIM == 3)
                              const int&,
                              const int&,
#endif
                              const double*);

    void SST_K_EQN_BUOYANCY(const double*,
                            const int&,
                            const double*,
                            const int&,
                            const double*,
                            const int&,
                            const double&,
                            const double&,
#if (NDIM == 3)
                            const double&,
#endif
                            const double&,
                            const int&,
                            const int&,
                            const int&,
                            const int&,
#if (NDIM == 3)
                            const int&,
                            const int&,
#endif
                            const double*);

    void SST_W_EQN_PRODUCTION(const double*,
                              const int&,
                              const double*,
                              const int&,
                              const double*,
                              const int&,
                              const double*,
                              const int&,
                              const double*,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
#if (NDIM == 3)
                              const int&,
                              const int&,
#endif
                              const double&,
                              const double&);

    void SST_W_EQN_CROSSDIFFUSION(const double*,
                                  const int&,
                                  const double*,
                                  const int&,
                                  const double*,
                                  const int&,
                                  const double*,
                                  const int&,
                                  const double*,
                                  const int&,
                                  const double&,
                                  const int&,
                                  const int&,
                                  const int&,
                                  const int&,
#if (NDIM == 3)
                                  const int&,
                                  const int&,
#endif
                                  const double*);
}

namespace IBAMR
{
namespace
{
// turbulent Prandtl number
static const double SIGMA_T = 0.85;
static const double BETA_STAR = 0.09;
static const double ALPHA_1 = 0.5532;
static const double ALPHA_2 = 0.4403;
static const double SIGMA_W2 = 0.856;
static const double KAPPA = 0.4187;
static const double B = 5.25;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;
} // namespace

TurbulenceSSTKOmegaSourceFunction::TurbulenceSSTKOmegaSourceFunction(
    const std::string& object_name,
    Pointer<SAMRAI::tbox::Database> input_db,
    TwoEquationTurbulenceHierarchyIntegrator* turb_hier_integrator,
    INSVCStaggeredHierarchyIntegrator* ins_hier_integrator)
    : CartGridFunction(object_name),
      d_turb_hier_integrator(turb_hier_integrator),
      d_ins_hier_integrator(ins_hier_integrator)
{
    double gravity[NDIM];
    if (input_db)
    {
        input_db->getDoubleArray("gravity", gravity, NDIM);
        d_gravity.x() = gravity[0];
        d_gravity.y() = gravity[1];
        if (NDIM == 3) d_gravity.z() = gravity[2];
    }
}

bool
TurbulenceSSTKOmegaSourceFunction::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
TurbulenceSSTKOmegaSourceFunction::setDataOnPatchHierarchy(const int data_idx,
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
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    d_hierarchy = hierarchy;

    // Hierarchy data operations object
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);
    HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(hierarchy, coarsest_ln, finest_ln);
    HierarchyNodeDataOpsReal<NDIM, double> hier_nc_data_ops(hierarchy, coarsest_ln, finest_ln);

    // Get the data index of the turbulent viscosity associated with new index.
    Pointer<CellVariable<NDIM, double> > mu_t_var = d_ins_hier_integrator->getTurbulentViscosityVariable();
    d_mu_t_new_idx = var_db->mapVariableAndContextToIndex(mu_t_var, d_ins_hier_integrator->getNewContext());

    // Get the data index of the molecular viscosity associated with new index.
    Pointer<CellVariable<NDIM, double> > mu_var = d_ins_hier_integrator->getViscosityVariable();
    d_mu_new_idx = var_db->mapVariableAndContextToIndex(mu_var, d_ins_hier_integrator->getNewContext());

    // Get the data index of the density associated with new index.
    Pointer<CellVariable<NDIM, double> > rho_cc_var = d_turb_hier_integrator->getCellCenteredMassDensityVariable();
    d_rho_new_idx = var_db->mapVariableAndContextToIndex(rho_cc_var, d_turb_hier_integrator->getNewContext());
    d_rho_scratch_idx = var_db->mapVariableAndContextToIndex(rho_cc_var, d_turb_hier_integrator->getScratchContext());

    // copying the new index data into scratch index data.
    hier_cc_data_ops.copyData(d_rho_scratch_idx, d_rho_new_idx);

    // filling ghost cells for density. INSVCStaggeredConservativeHierarchyIntegrator
    // class works with side-centered density which has NDIM components whereas
    // TwoEquationTurbulenceHierarchyIntegrator class works with cell-centered density.
    // Therefore, consider any one element of the bc_coef vector for density bc_coef object.
    std::vector<RobinBcCoefStrategy<NDIM>*> rho_bc_coefs = d_ins_hier_integrator->getMassDensityBoundaryConditions();
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent rho_bc_component(d_rho_scratch_idx,
                                                       DATA_REFINE_TYPE,
                                                       USE_CF_INTERPOLATION,
                                                       DATA_COARSEN_TYPE,
                                                       BDRY_EXTRAP_TYPE,
                                                       CONSISTENT_TYPE_2_BDRY,
                                                       rho_bc_coefs[0]);
    HierarchyGhostCellInterpolation rho_ghost_cell_fill_op;
    rho_ghost_cell_fill_op.initializeOperatorState(rho_bc_component, hierarchy, coarsest_ln, finest_ln);
    rho_ghost_cell_fill_op.fillData(data_time);

    // filling ghost cells for molecular viscosity
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* mu_bc_coef = d_ins_hier_integrator->getViscosityBoundaryConditions();
    d_mu_scratch_idx = var_db->mapVariableAndContextToIndex(mu_var, d_ins_hier_integrator->getScratchContext());
    hier_cc_data_ops.copyData(d_mu_scratch_idx, d_mu_new_idx);
    InterpolationTransactionComponent mu_bc_component(d_mu_scratch_idx,
                                                      DATA_REFINE_TYPE,
                                                      USE_CF_INTERPOLATION,
                                                      DATA_COARSEN_TYPE,
                                                      BDRY_EXTRAP_TYPE,
                                                      CONSISTENT_TYPE_2_BDRY,
                                                      mu_bc_coef);
    HierarchyGhostCellInterpolation mu_ghost_cell_fill_op;
    mu_ghost_cell_fill_op.initializeOperatorState(mu_bc_component, d_hierarchy);
    mu_ghost_cell_fill_op.fillData(data_time);

    // Get the data index of the k variable associated with the new index.
    Pointer<CellVariable<NDIM, double> > k_var = d_turb_hier_integrator->getKVariable();
    d_k_new_idx = var_db->mapVariableAndContextToIndex(k_var, d_turb_hier_integrator->getNewContext());
    d_k_scratch_idx = var_db->mapVariableAndContextToIndex(k_var, d_turb_hier_integrator->getScratchContext());
    // copying the new index data into scratch index data
    hier_cc_data_ops.copyData(d_k_scratch_idx, d_k_new_idx);

    // filling ghost cells for k variable.
    RobinBcCoefStrategy<NDIM>* k_bc_coef = d_turb_hier_integrator->getPhysicalBcCoefKEquation();
    InterpolationTransactionComponent k_ghost_cc_interpolation(d_k_scratch_idx,
                                                               DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               k_bc_coef);
    HierarchyGhostCellInterpolation k_ghost_cell_fill_op;
    k_ghost_cell_fill_op.initializeOperatorState(k_ghost_cc_interpolation, hierarchy, coarsest_ln, finest_ln);
    k_ghost_cell_fill_op.fillData(data_time);

    // Get the data index of the w variable, associated with the new index.
    Pointer<CellVariable<NDIM, double> > w_var = d_turb_hier_integrator->getWVariable();
    d_w_new_idx = var_db->mapVariableAndContextToIndex(w_var, d_turb_hier_integrator->getNewContext());
    d_w_scratch_idx = var_db->mapVariableAndContextToIndex(w_var, d_turb_hier_integrator->getScratchContext());
    // copying the new index data into scratch index data
    hier_cc_data_ops.copyData(d_w_scratch_idx, d_w_new_idx);

    // filling ghost cells for  w variable.
    RobinBcCoefStrategy<NDIM>* w_bc_coef = d_turb_hier_integrator->getPhysicalBcCoefWEquation();
    InterpolationTransactionComponent w_ghost_cc_interpolation(d_w_scratch_idx,
                                                               DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               w_bc_coef);
    HierarchyGhostCellInterpolation w_ghost_cell_fill_op;
    w_ghost_cell_fill_op.initializeOperatorState(w_ghost_cc_interpolation, hierarchy, coarsest_ln, finest_ln);
    w_ghost_cell_fill_op.fillData(data_time);

    Pointer<NodeVariable<NDIM, double> > tau_w_var = d_ins_hier_integrator->getTauwVariable();
    d_tau_w_idx = var_db->mapVariableAndContextToIndex(tau_w_var, d_ins_hier_integrator->getScratchContext());

    Pointer<SideVariable<NDIM, double> > U_tau_var = d_ins_hier_integrator->getUtauVariable();
    d_U_tau_idx = var_db->mapVariableAndContextToIndex(U_tau_var, d_ins_hier_integrator->getScratchContext());

    // Get the data index of the blending function,F1, associated with scratch index.
    d_f1_scratch_idx = getBlendingFunctionVariableIndex(d_turb_hier_integrator);

    // Get the data index of the production variable, associated with scratch index.
    d_p_scratch_idx = getProductionVariableIndex(d_turb_hier_integrator);

    // Get the wall location indices.
    d_wall_location_index = getWallLocationIndices(d_turb_hier_integrator);

    // Fill data on each patch level.
    CartGridFunction::setDataOnPatchHierarchy(
        data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);

    return;
}

// setting the data on patches
void
TurbulenceSSTKOmegaSourceFunction::setDataOnPatch(int data_idx,
                                                  Pointer<Variable<NDIM> > var,
                                                  Pointer<Patch<NDIM> > patch,
                                                  double data_time,
                                                  bool initial_time,
                                                  Pointer<PatchLevel<NDIM> > level)
{
    Pointer<CellData<NDIM, double> > f_data = patch->getPatchData(data_idx);
    if (var->getName() == "turbulent_kinetic_energy::F")
    {
        setDataOnPatchCellForK(f_data, patch, data_time, initial_time, level);
    }
    else if (var->getName() == "turbulent_specific_dissipation_rate::F")
    {
        setDataOnPatchCellForOmega(f_data, patch, data_time, initial_time, level);
    }

    return;
}
/////////////////////////////// PRIVATE ////////////////////////////
int
TurbulenceSSTKOmegaSourceFunction::getBlendingFunctionVariableIndex(
    Pointer<TwoEquationTurbulenceHierarchyIntegrator> turb_hier_integrator)
{
    return turb_hier_integrator->d_F1_scratch_idx;
} // getBlendingFunctionVariableIndex

int
TurbulenceSSTKOmegaSourceFunction::getProductionVariableIndex(
    Pointer<TwoEquationTurbulenceHierarchyIntegrator> turb_hier_integrator)
{
    return turb_hier_integrator->d_p_scratch_idx;
} // getProductionVariableIndex

SAMRAI::tbox::Array<int>
TurbulenceSSTKOmegaSourceFunction::getWallLocationIndices(
    Pointer<TwoEquationTurbulenceHierarchyIntegrator> turb_hier_integrator)
{
    return turb_hier_integrator->d_wall_location_index;
} // getWallLocationIndices

void
TurbulenceSSTKOmegaSourceFunction::setDataOnPatchCellForK(Pointer<CellData<NDIM, double> > k_f_data,
                                                          Pointer<Patch<NDIM> > patch,
                                                          const double data_time,
                                                          const bool initial_time,
                                                          Pointer<PatchLevel<NDIM> > level)
{
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    Pointer<CellData<NDIM, double> > mu_t_data = patch->getPatchData(d_mu_t_new_idx);
    Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(d_mu_scratch_idx);
    Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(d_rho_scratch_idx);
    Pointer<NodeData<NDIM, double> > tau_w_data = patch->getPatchData(d_tau_w_idx);
    Pointer<SideData<NDIM, double> > U_tau_data = patch->getPatchData(d_U_tau_idx);
    const IntVector<NDIM>& U_tau_ghost_cells = U_tau_data->getGhostCellWidth();

    // find P = min (G, 10*Beta_star*k*w)
    Pointer<CellData<NDIM, double> > k_data = patch->getPatchData(d_k_scratch_idx);
    Pointer<CellData<NDIM, double> > w_data = patch->getPatchData(d_w_new_idx);
    Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(d_p_scratch_idx);
    SST_K_EQN_PRODUCTION(k_f_data->getPointer(),
                         k_f_data->getGhostCellWidth().max(),
                         p_data->getPointer(),
                         p_data->getGhostCellWidth().max(),
                         k_data->getPointer(),
                         k_data->getGhostCellWidth().max(),
                         w_data->getPointer(),
                         w_data->getGhostCellWidth().max(),
                         patch_box.lower(0),
                         patch_box.upper(0),
                         patch_box.lower(1),
                         patch_box.upper(1),
#if (NDIM == 3)
                         patch_box.lower(2),
                         patch_box.upper(2),
#endif
                         BETA_STAR);

    // Add production term calculated from wall law here.
    Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    const double* const patch_dx = patch_geom->getDx();

    for (unsigned int i = 0; i < d_wall_location_index.size(); i++)
    {
        const unsigned int axis = d_wall_location_index[i] / 2;
        const unsigned int side = d_wall_location_index[i] % 2;
        const bool is_lower = side == 0;
        if (patch_geom->getTouchesRegularBoundary(axis, side))
        {
            Box<NDIM> bdry_box = patch_box;
            if (is_lower)
            {
                bdry_box.upper(axis) = patch_box.lower(axis);
            }
            else
            {
                bdry_box.lower(axis) = patch_box.upper(axis);
            }

            Box<NDIM> trim_box = patch_box * bdry_box;
            NEAR_WALL_PRODUCTION(k_f_data->getPointer(),
                                 k_f_data->getGhostCellWidth().max(),
                                 tau_w_data->getPointer(),
                                 tau_w_data->getGhostCellWidth().max(),
                                 U_tau_data->getPointer(0),
                                 U_tau_ghost_cells(0),
                                 U_tau_data->getPointer(1),
                                 U_tau_ghost_cells(1),
#if (NDIM == 3)
                                 U_tau_data->getPointer(2),
                                 U_tau_ghost_cells(2),
#endif
                                 KAPPA,
                                 d_wall_location_index[i],
                                 trim_box.lower(0),
                                 trim_box.upper(0),
                                 trim_box.lower(1),
                                 trim_box.upper(1),
#if (NDIM == 3)
                                 trim_box.lower(2),
                                 trim_box.upper(2),
#endif
                                 patch_box.lower(0),
                                 patch_box.upper(0),
                                 patch_box.lower(1),
                                 patch_box.upper(1),
#if (NDIM == 3)
                                 patch_box.lower(2),
                                 patch_box.upper(2),
#endif

                                 patch_dx);
        }
    }

    // routine to calculate buoyancy term
    SST_K_EQN_BUOYANCY(k_f_data->getPointer(),
                       k_f_data->getGhostCellWidth().max(),
                       mu_t_data->getPointer(),
                       mu_t_data->getGhostCellWidth().max(),
                       rho_data->getPointer(),
                       rho_data->getGhostCellWidth().max(),
                       d_gravity[0],
                       d_gravity[1],
#if (NDIM == 3)
                       d_gravity[2],
#endif
                       SIGMA_T,
                       patch_box.lower(0),
                       patch_box.upper(0),
                       patch_box.lower(1),
                       patch_box.upper(1),
#if (NDIM == 3)
                       patch_box.lower(2),
                       patch_box.upper(2),
#endif
                       dx);
    return;
}

void
TurbulenceSSTKOmegaSourceFunction::setDataOnPatchCellForOmega(Pointer<CellData<NDIM, double> > w_f_data,
                                                              Pointer<Patch<NDIM> > patch,
                                                              const double data_time,
                                                              const bool initial_time,
                                                              Pointer<PatchLevel<NDIM> > level)
{
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    Pointer<CellData<NDIM, double> > mu_t_data = patch->getPatchData(d_mu_t_new_idx);

    Pointer<CellData<NDIM, double> > k_data = patch->getPatchData(d_k_scratch_idx);
    Pointer<CellData<NDIM, double> > w_data = patch->getPatchData(d_w_scratch_idx);
    Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(d_rho_new_idx);
    Pointer<CellData<NDIM, double> > f1_data = patch->getPatchData(d_f1_scratch_idx);
    Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(d_p_scratch_idx);

    // routine to calculate (alpha/nut)*G
    SST_W_EQN_PRODUCTION(w_f_data->getPointer(),
                         w_f_data->getGhostCellWidth().max(),
                         p_data->getPointer(),
                         p_data->getGhostCellWidth().max(),
                         mu_t_data->getPointer(),
                         mu_t_data->getGhostCellWidth().max(),
                         rho_data->getPointer(),
                         rho_data->getGhostCellWidth().max(),
                         f1_data->getPointer(),
                         f1_data->getGhostCellWidth().max(),
                         patch_box.lower(0),
                         patch_box.upper(0),
                         patch_box.lower(1),
                         patch_box.upper(1),
#if (NDIM == 3)
                         patch_box.lower(2),
                         patch_box.upper(2),
#endif
                         ALPHA_1,
                         ALPHA_2);

    SST_W_EQN_CROSSDIFFUSION(w_f_data->getPointer(),
                             w_f_data->getGhostCellWidth().max(),
                             rho_data->getPointer(),
                             rho_data->getGhostCellWidth().max(),
                             f1_data->getPointer(),
                             f1_data->getGhostCellWidth().max(),
                             k_data->getPointer(),
                             k_data->getGhostCellWidth().max(),
                             w_data->getPointer(),
                             w_data->getGhostCellWidth().max(),
                             SIGMA_W2,
                             patch_box.lower(0),
                             patch_box.upper(0),
                             patch_box.lower(1),
                             patch_box.upper(1),
#if (NDIM == 3)
                             patch_box.lower(2),
                             patch_box.upper(2),
#endif
                             dx);
}

} // namespace IBAMR
