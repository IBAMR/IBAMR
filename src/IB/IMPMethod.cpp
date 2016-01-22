// Filename: IMPMethod.cpp
// Created on 16 Oct 2012 by Boyce Griffith
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

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <algorithm>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CoarsenSchedule.h"
#include "GriddingAlgorithm.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineSchedule.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "boost/array.hpp"
#include "boost/math/special_functions/round.hpp"
#include "boost/multi_array.hpp"
#include "ibamr/IMPMethod.h"
#include "ibamr/MaterialPointSpec.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/LIndexSetData.h"
#include "ibtk/LInitStrategy.h"
#include "ibtk/LMesh.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeSet.h"
#include "ibtk/LNodeSetData.h"
#include "ibtk/LSet.h"
#include "ibtk/LSetData.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/libmesh_utilities.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_tensor.h"
#include "libmesh/type_vector.h"
#include "libmesh/vector_value.h"
#include "petscvec.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

namespace IBTK
{
} // namespace IBTK

using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const std::string KERNEL_FCN = "IB_6";
static const int kernel_width = 3;
void
kernel(const double X,
       const double patch_x_lower,
       const double /*patch_x_upper*/,
       const double dx,
       const int patch_box_lower,
       const int /*patch_box_upper*/,
       int& stencil_box_lower,
       int& stencil_box_upper,
       boost::multi_array<double, 1>& phi,
       boost::multi_array<double, 1>& dphi)
{
    const double X_o_dx = (X - patch_x_lower) / dx;
    stencil_box_lower = boost::math::round(X_o_dx) + patch_box_lower - kernel_width;
    stencil_box_upper = stencil_box_lower + 2 * kernel_width - 1;
    const double r = 1.0 - X_o_dx + ((stencil_box_lower + kernel_width - 1 - patch_box_lower) + 0.5);

    const double r2 = r * r;
    const double r3 = r * r2;
    const double r4 = r * r3;
    const double r5 = r * r4;
    const double r6 = r * r5;

    static const double K = (59.0 / 60.0) * (1.0 - sqrt(1.0 - (3220.0 / 3481.0)));
    static const double K2 = K * K;

    static const double alpha = 28.0;

    const double beta = (9.0 / 4.0) - (3.0 / 2.0) * (K + r2) + ((22.0 / 3.0) - 7.0 * K) * r - (7.0 / 3.0) * r3;
    const double dbeta = ((22.0 / 3.0) - 7.0 * K) - 3.0 * r - 7.0 * r2;

    const double gamma = (1.0 / 4.0) * (((161.0 / 36.0) - (59.0 / 6.0) * K + 5.0 * K2) * (1.0 / 2.0) * r2 +
                                        (-(109.0 / 24.0) + 5.0 * K) * (1.0 / 3.0) * r4 + (5.0 / 18.0) * r6);
    const double dgamma = (1.0 / 4.0) * (((161.0 / 36.0) - (59.0 / 6.0) * K + 5.0 * K2) * r +
                                         (-(109.0 / 24.0) + 5.0 * K) * (4.0 / 3.0) * r3 + (5.0 / 3.0) * r5);

    const double discr = beta * beta - 4.0 * alpha * gamma;

    phi[0] = (-beta + copysign(1.0, (3.0 / 2.0) - K) * sqrt(discr)) / (2.0 * alpha);
    phi[1] =
        -3.0 * phi[0] - (1.0 / 16.0) + (1.0 / 8.0) * (K + r2) + (1.0 / 12.0) * (3.0 * K - 1.0) * r + (1.0 / 12.0) * r3;
    phi[2] = 2.0 * phi[0] + (1.0 / 4.0) + (1.0 / 6.0) * (4.0 - 3.0 * K) * r - (1.0 / 6.0) * r3;
    phi[3] = 2.0 * phi[0] + (5.0 / 8.0) - (1.0 / 4.0) * (K + r2);
    phi[4] = -3.0 * phi[0] + (1.0 / 4.0) - (1.0 / 6.0) * (4.0 - 3.0 * K) * r + (1.0 / 6.0) * r3;
    phi[5] = phi[0] - (1.0 / 16.0) + (1.0 / 8.0) * (K + r2) - (1.0 / 12.0) * (3.0 * K - 1.0) * r - (1.0 / 12.0) * r3;

    dphi[0] = -(dbeta * phi[0] + dgamma) / (2.0 * alpha * phi[0] + beta);
    dphi[1] = -3.0 * dphi[0] + (1.0 / 12.0) * (3.0 * K - 1.0) + (1.0 / 4.0) * r + (1.0 / 4.0) * r2;
    dphi[2] = 2.0 * dphi[0] + (1.0 / 6.0) * (4.0 - 3.0 * K) - (1.0 / 2.0) * r2;
    dphi[3] = 2.0 * dphi[0] - (1.0 / 2.0) * r;
    dphi[4] = -3.0 * dphi[0] - (1.0 / 6.0) * (4.0 - 3.0 * K) + (1.0 / 2.0) * r2;
    dphi[5] = dphi[0] - (1.0 / 12.0) * (3.0 * K - 1.0) + (1.0 / 4.0) * r - (1.0 / 4.0) * r2;
    return;
}

// Version of IMPMethod restart file data.
static const int IMP_METHOD_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IMPMethod::IMPMethod(const std::string& object_name, Pointer<Database> input_db, bool register_for_restart)
{
    // Set the object name and register it with the restart manager.
    d_object_name = object_name;
    d_registered_for_restart = false;
    if (register_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
        d_registered_for_restart = true;
    }

    // Ensure all pointers to helper objects are NULL.
    d_l_initializer = NULL;
    d_silo_writer = NULL;

    // Set some default values.
    d_ghosts = LEInteractor::getMinimumGhostWidth(KERNEL_FCN);
    d_do_log = false;

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Get the Lagrangian Data Manager.
    d_l_data_manager = LDataManager::getManager(
        d_object_name + "::LDataManager", KERNEL_FCN, KERNEL_FCN, d_ghosts, d_registered_for_restart);
    d_ghosts = d_l_data_manager->getGhostCellWidth();

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_half_time = std::numeric_limits<double>::quiet_NaN();

    // Indicate all Lagrangian data needs ghost values to be refilled, and that
    // all intermediate data needs to be initialized.
    d_X_current_needs_ghost_fill = true;
    d_X_new_needs_ghost_fill = true;
    d_X_half_needs_ghost_fill = true;
    d_X_half_needs_reinit = true;
    d_U_half_needs_reinit = true;
    return;
} // IMPMethod

IMPMethod::~IMPMethod()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
} // ~IMPMethod

void
IMPMethod::registerLInitStrategy(Pointer<LInitStrategy> l_initializer)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(l_initializer);
#endif
    d_l_initializer = l_initializer;
    d_l_data_manager->registerLInitStrategy(d_l_initializer);
    return;
} // registerLInitStrategy

void
IMPMethod::freeLInitStrategy()
{
    d_l_initializer.setNull();
    d_l_data_manager->freeLInitStrategy();
    return;
} // freeLInitStrategy

LDataManager*
IMPMethod::getLDataManager() const
{
    return d_l_data_manager;
} // getLDataManager

void
IMPMethod::registerLSiloDataWriter(Pointer<LSiloDataWriter> silo_writer)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(silo_writer);
#endif
    d_silo_writer = silo_writer;
    d_l_data_manager->registerLSiloDataWriter(d_silo_writer);
    return;
} // registerLSiloDataWriter

const IntVector<NDIM>&
IMPMethod::getMinimumGhostCellWidth() const
{
    return d_ghosts;
} // getMinimumGhostCellWidth

void
IMPMethod::setupTagBuffer(Array<int>& tag_buffer, Pointer<GriddingAlgorithm<NDIM> > gridding_alg) const
{
    const int finest_hier_ln = gridding_alg->getMaxLevels() - 1;
    const int tsize = tag_buffer.size();
    tag_buffer.resizeArray(finest_hier_ln);
    for (int i = tsize; i < finest_hier_ln; ++i) tag_buffer[i] = 0;
    const int gcw = d_ghosts.max();
    for (int tag_ln = 0; tag_ln < finest_hier_ln; ++tag_ln)
    {
        const int data_ln = tag_ln + 1;
        const int can_be_refined = data_ln < finest_hier_ln;
        if (!d_l_initializer->getLevelHasLagrangianData(data_ln, can_be_refined)) continue;
        tag_buffer[tag_ln] = std::max(tag_buffer[tag_ln], gcw);
    }
    for (int ln = finest_hier_ln - 2; ln >= 0; --ln)
    {
        tag_buffer[ln] =
            std::max(tag_buffer[ln], tag_buffer[ln + 1] / gridding_alg->getRatioToCoarserLevel(ln + 1).max() + 1);
    }
    return;
} // setupTagBuffer

void
IMPMethod::preprocessIntegrateData(double current_time, double new_time, int /*num_cycles*/)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_half_time = current_time + 0.5 * (new_time - current_time);

    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Look-up or allocate Lagangian data.
    d_X_current_data.resize(finest_ln + 1);
    d_X_new_data.resize(finest_ln + 1);
    d_X_half_data.resize(finest_ln + 1);
    d_X0_data.resize(finest_ln + 1);
    d_U_current_data.resize(finest_ln + 1);
    d_U_new_data.resize(finest_ln + 1);
    d_U_half_data.resize(finest_ln + 1);
    d_Grad_U_current_data.resize(finest_ln + 1);
    d_Grad_U_new_data.resize(finest_ln + 1);
    d_Grad_U_half_data.resize(finest_ln + 1);
    d_F_current_data.resize(finest_ln + 1);
    d_F_new_data.resize(finest_ln + 1);
    d_F_half_data.resize(finest_ln + 1);
    d_tau_data.resize(finest_ln + 1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        d_X_current_data[ln] = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME, ln);
        d_X_new_data[ln] = d_l_data_manager->createLData("X_new", ln, NDIM);
        d_X0_data[ln] = d_l_data_manager->getLData("X0", ln);
        d_U_current_data[ln] = d_l_data_manager->getLData(LDataManager::VEL_DATA_NAME, ln);
        d_U_new_data[ln] = d_l_data_manager->createLData("U_new", ln, NDIM);
        d_Grad_U_current_data[ln] = d_l_data_manager->getLData("Grad_U", ln);
        d_Grad_U_new_data[ln] = d_l_data_manager->createLData("Grad_U_new", ln, NDIM * NDIM);
        d_F_current_data[ln] = d_l_data_manager->getLData("F", ln);
        d_F_new_data[ln] = d_l_data_manager->createLData("F_new", ln, NDIM * NDIM);
        d_F_half_data[ln] = d_l_data_manager->createLData("F_half", ln, NDIM * NDIM);
        d_tau_data[ln] = d_l_data_manager->getLData("tau", ln);

        // Initialize new values to equal current values.
        ierr = VecCopy(d_X_current_data[ln]->getVec(), d_X_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
        ierr = VecCopy(d_U_current_data[ln]->getVec(), d_U_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
        ierr = VecCopy(d_Grad_U_current_data[ln]->getVec(), d_Grad_U_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
        ierr = VecCopy(d_F_current_data[ln]->getVec(), d_F_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
        ierr = VecCopy(d_F_current_data[ln]->getVec(), d_F_half_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
    }

    // Keep track of Lagrangian data objects that need to have ghost values
    // filled, or that need to be reinitialized.
    d_X_new_needs_ghost_fill = true;
    d_X_half_needs_reinit = true;
    d_U_half_needs_reinit = true;
    return;
} // preprocessIntegrateData

void
IMPMethod::postprocessIntegrateData(double /*current_time*/, double /*new_time*/, int /*num_cycles*/)
{
    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Reset time-dependent Lagrangian data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        ierr = VecSwap(d_X_current_data[ln]->getVec(), d_X_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
        ierr = VecSwap(d_U_current_data[ln]->getVec(), d_U_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
        ierr = VecSwap(d_Grad_U_current_data[ln]->getVec(), d_Grad_U_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
        ierr = VecSwap(d_F_current_data[ln]->getVec(), d_F_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
    }
    d_X_current_needs_ghost_fill = true;

    // Deallocate Lagrangian scratch data.
    d_X_current_data.clear();
    d_X_new_data.clear();
    d_X_half_data.clear();
    d_X0_data.clear();
    d_U_current_data.clear();
    d_U_new_data.clear();
    d_U_half_data.clear();
    d_Grad_U_current_data.clear();
    d_Grad_U_new_data.clear();
    d_Grad_U_half_data.clear();
    d_F_current_data.clear();
    d_F_new_data.clear();
    d_F_half_data.clear();
    d_tau_data.clear();

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_half_time = std::numeric_limits<double>::quiet_NaN();
    return;
} // postprocessIntegrateData

void
IMPMethod::interpolateVelocity(const int u_data_idx,
                               const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                               const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                               const double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Determine the type of data centering.
    Pointer<hier::Variable<NDIM> > u_var;
    var_db->mapIndexToVariable(u_data_idx, u_var);
    Pointer<SideVariable<NDIM, double> > u_sc_var = u_var;
    const bool sc_data = u_sc_var;
    TBOX_ASSERT(sc_data);

    // Synchronize Eulerian and Lagrangian values.
    std::vector<Pointer<LData> > *U_data, *Grad_U_data, *X_data;
    bool* X_needs_ghost_fill;
    getVelocityData(&U_data, &Grad_U_data, data_time);
    getPositionData(&X_data, &X_needs_ghost_fill, data_time);
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        if (ln < static_cast<int>(u_synch_scheds.size()) && u_synch_scheds[ln])
        {
            u_synch_scheds[ln]->coarsenData();
        }
        if (*X_needs_ghost_fill && (*X_data)[ln]) (*X_data)[ln]->beginGhostUpdate();
    }
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        if (*X_needs_ghost_fill && (*X_data)[ln]) (*X_data)[ln]->endGhostUpdate();
    }
    *X_needs_ghost_fill = false;

    // Interpolate data from the Eulerian grid to the Lagrangian mesh.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        if (ln < static_cast<int>(u_ghost_fill_scheds.size()) && u_ghost_fill_scheds[ln])
        {
            u_ghost_fill_scheds[ln]->fillData(data_time);
        }
        boost::multi_array_ref<double, 2>& U_array = *(*U_data)[ln]->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& Grad_U_array = *(*Grad_U_data)[ln]->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& X_array = *(*X_data)[ln]->getLocalFormVecArray();
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_data_idx);
            Pointer<LNodeSetData> idx_data = patch->getPatchData(d_l_data_manager->getLNodePatchDescriptorIndex());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const x_lower = patch_geom->getXLower();
            const double* const x_upper = patch_geom->getXUpper();
            const double* const dx = patch_geom->getDx();
            Box<NDIM> side_boxes[NDIM];
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                side_boxes[axis] = SideGeometry<NDIM>::toSideBox(u_data->getGhostBox() * idx_data->getGhostBox(), axis);
            }
            for (LNodeSetData::CellIterator it(idx_data->getGhostBox()); it; it++)
            {
                const Index<NDIM>& i = *it;
                LNodeSet* const node_set = idx_data->getItem(i);
                if (!node_set) continue;
                for (LNodeSet::iterator it = node_set->begin(); it != node_set->end(); ++it)
                {
                    const LNode* const node_idx = *it;
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const double* const X = &X_array[local_idx][0];
                    VectorValue<double> U;
                    TensorValue<double> Grad_U;

                    // Interpolate U and Grad U using a smoothed kernel
                    // function evaluated about X.
                    Box<NDIM> stencil_box;
                    const int stencil_size = LEInteractor::getStencilSize(KERNEL_FCN);
                    boost::array<boost::multi_array<double, 1>, NDIM> phi, dphi;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        phi[d].resize(boost::extents[stencil_size]);
                        dphi[d].resize(boost::extents[stencil_size]);
                    }
                    for (unsigned int component = 0; component < NDIM; ++component)
                    {
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            kernel(X[d],
                                   x_lower[d] + (d == component ? -0.5 * dx[d] : 0.0),
                                   x_upper[d] + (d == component ? +0.5 * dx[d] : 0.0),
                                   dx[d],
                                   patch_box.lower(d),
                                   patch_box.upper(d) + (d == component ? 1 : 0),
                                   stencil_box.lower(d),
                                   stencil_box.upper(d),
                                   phi[d],
                                   dphi[d]);
                        }
                        for (Box<NDIM>::Iterator b(stencil_box * side_boxes[component]); b; b++)
                        {
                            const Index<NDIM>& i = b();
                            const Index<NDIM> i_shift = i - stencil_box.lower();
                            const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);
                            const double u = (*u_data)(i_s);
                            double w = 1.0;
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                w *= phi[d][i_shift(d)];
                            }
                            U(component) += u * w;
                            for (unsigned int k = 0; k < NDIM; ++k)
                            {
                                double dw_dx_k = 1.0;
                                for (unsigned int d = 0; d < NDIM; ++d)
                                {
                                    if (d == k)
                                    {
                                        dw_dx_k *= dphi[d][i_shift(d)] / dx[d];
                                    }
                                    else
                                    {
                                        dw_dx_k *= phi[d][i_shift(d)];
                                    }
                                }
                                Grad_U(component, k) -= u * dw_dx_k;
                            }
                        }
                    }
                    for (int i = 0; i < NDIM; ++i)
                    {
                        U_array[local_idx][i] = U(i);
                        for (int j = 0; j < NDIM; ++j)
                        {
                            Grad_U_array[local_idx][NDIM * i + j] = Grad_U(i, j);
                        }
                    }
                }
            }
        }
        (*U_data)[ln]->restoreArrays();
        (*Grad_U_data)[ln]->restoreArrays();
        (*X_data)[ln]->restoreArrays();
    }
    d_U_half_needs_reinit = !MathUtilities<double>::equalEps(data_time, d_half_time);
    return;
} // interpolateVelocity

void
IMPMethod::eulerStep(const double current_time, const double new_time)
{
    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    std::vector<Pointer<LData> > *U_data, *Grad_U_data;
    getVelocityData(&U_data, &Grad_U_data, current_time);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        // Update the positions.
        ierr = VecWAXPY(d_X_new_data[ln]->getVec(), dt, (*U_data)[ln]->getVec(), d_X_current_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);

        // Update the deformation gradient.
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        boost::multi_array_ref<double, 2>& F_current_array = *d_F_current_data[ln]->getVecArray();
        boost::multi_array_ref<double, 2>& F_new_array = *d_F_new_data[ln]->getVecArray();
        boost::multi_array_ref<double, 2>& F_half_array = *d_F_half_data[ln]->getVecArray();
        boost::multi_array_ref<double, 2>& Grad_U_array = *(*Grad_U_data)[ln]->getVecArray();
        TensorValue<double> F_current, F_new, F_half, Grad_U;
        TensorValue<double> I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int idx = node_idx->getGlobalPETScIndex();
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    F_current(i, j) = F_current_array[idx][NDIM * i + j];
                    Grad_U(i, j) = Grad_U_array[idx][NDIM * i + j];
                }
            }
#if (NDIM == 2)
            F_current(2, 2) = 1.0;
#endif
            F_new = tensor_inverse(I - 0.50 * dt * Grad_U) * (I + 0.50 * dt * Grad_U) * F_current;
            F_half = tensor_inverse(I - 0.25 * dt * Grad_U) * (I + 0.25 * dt * Grad_U) * F_current;
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    F_new_array[idx][NDIM * i + j] = F_new(i, j);
                    F_half_array[idx][NDIM * i + j] = F_half(i, j);
                }
            }
        }
    }
    d_X_new_needs_ghost_fill = true;
    d_X_half_needs_reinit = true;
    return;
} // eulerStep

void
IMPMethod::midpointStep(const double current_time, const double new_time)
{
    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    std::vector<Pointer<LData> > *U_data, *Grad_U_data;
    getVelocityData(&U_data, &Grad_U_data, current_time + 0.5 * dt);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        // Update the positions.
        ierr = VecWAXPY(d_X_new_data[ln]->getVec(), dt, (*U_data)[ln]->getVec(), d_X_current_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);

        // Update the deformation gradient.
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        boost::multi_array_ref<double, 2>& F_current_array = *d_F_current_data[ln]->getVecArray();
        boost::multi_array_ref<double, 2>& F_new_array = *d_F_new_data[ln]->getVecArray();
        boost::multi_array_ref<double, 2>& Grad_U_array = *(*Grad_U_data)[ln]->getVecArray();
        TensorValue<double> F_current, F_new, Grad_U;
        TensorValue<double> I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int idx = node_idx->getGlobalPETScIndex();
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    F_current(i, j) = F_current_array[idx][NDIM * i + j];
                    Grad_U(i, j) = Grad_U_array[idx][NDIM * i + j];
                }
            }
#if (NDIM == 2)
            F_current(2, 2) = 1.0;
#endif
            F_new = tensor_inverse(I - 0.5 * dt * Grad_U) * (I + 0.5 * dt * Grad_U) * F_current;
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    F_new_array[idx][NDIM * i + j] = F_new(i, j);
                }
            }
        }
    }
    d_X_new_needs_ghost_fill = true;
    d_X_half_needs_reinit = true;
    return;
} // midpointStep

void
IMPMethod::trapezoidalStep(const double current_time, const double new_time)
{
    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    std::vector<Pointer<LData> > *U_current_data, *U_new_data, *Grad_U_current_data, *Grad_U_new_data;
    getVelocityData(&U_current_data, &Grad_U_current_data, current_time);
    getVelocityData(&U_new_data, &Grad_U_new_data, new_time);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        // Update the positions.
        ierr = VecWAXPY(
            d_X_new_data[ln]->getVec(), 0.5 * dt, (*U_current_data)[ln]->getVec(), d_X_current_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPY(d_X_new_data[ln]->getVec(), 0.5 * dt, (*U_new_data)[ln]->getVec());
        IBTK_CHKERRQ(ierr);

        // Update the deformation gradient.
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        boost::multi_array_ref<double, 2>& F_current_array = *d_F_current_data[ln]->getVecArray();
        boost::multi_array_ref<double, 2>& F_new_array = *d_F_new_data[ln]->getVecArray();
        boost::multi_array_ref<double, 2>& Grad_U_current_array = *(*Grad_U_current_data)[ln]->getVecArray();
        boost::multi_array_ref<double, 2>& Grad_U_new_array = *(*Grad_U_new_data)[ln]->getVecArray();
        TensorValue<double> F_current, F_new, Grad_U_current, Grad_U_new;
        TensorValue<double> I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int idx = node_idx->getGlobalPETScIndex();
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    F_current(i, j) = F_current_array[idx][NDIM * i + j];
                    Grad_U_current(i, j) = Grad_U_current_array[idx][NDIM * i + j];
                    Grad_U_new(i, j) = Grad_U_new_array[idx][NDIM * i + j];
                }
            }
#if (NDIM == 2)
            F_current(2, 2) = 1.0;
#endif
            F_new = tensor_inverse(I - 0.5 * dt * Grad_U_new) * (I + 0.5 * dt * Grad_U_current) * F_current;
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    F_new_array[idx][NDIM * i + j] = F_new(i, j);
                }
            }
        }
    }
    d_X_new_needs_ghost_fill = true;
    d_X_half_needs_reinit = true;
    return;
} // trapezoidalStep

void
IMPMethod::computeLagrangianForce(const double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    std::vector<Pointer<LData> > *X_data, *F_data;
    bool* X_needs_ghost_fill;
    getPositionData(&X_data, &X_needs_ghost_fill, data_time);
    getDeformationGradientData(&F_data, data_time);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        boost::multi_array_ref<double, 2>& x_array = *(*X_data)[ln]->getVecArray();
        boost::multi_array_ref<double, 2>& X_array = *d_X0_data[ln]->getVecArray();
        boost::multi_array_ref<double, 2>& F_array = *(*F_data)[ln]->getVecArray();
        boost::multi_array_ref<double, 2>& tau_array = *d_tau_data[ln]->getVecArray();
        TensorValue<double> FF, PP, tau;
        VectorValue<double> X, x;
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int idx = node_idx->getGlobalPETScIndex();
            MaterialPointSpec* mp_spec = node_idx->getNodeDataItem<MaterialPointSpec>();
            if (mp_spec && d_PK1_stress_fcn)
            {
                for (int i = 0; i < NDIM; ++i)
                {
                    for (int j = 0; j < NDIM; ++j)
                    {
                        FF(i, j) = F_array[idx][NDIM * i + j];
                    }
                    x(i) = x_array[idx][i];
                    X(i) = X_array[idx][i];
                }
#if (NDIM == 2)
                FF(2, 2) = 1.0;
#endif
                (*d_PK1_stress_fcn)(PP,
                                    FF,
                                    x,
                                    X,
                                    mp_spec->getSubdomainId(),
                                    mp_spec->getInternalVariables(),
                                    data_time,
                                    d_PK1_stress_fcn_ctx);
                tau = PP * FF.transpose();
            }
            else
            {
                tau.zero();
            }
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    tau_array[idx][NDIM * i + j] = tau(i, j);
                }
            }
        }
    }
    return;
} // computeLagrangianForce

void
IMPMethod::spreadForce(const int f_data_idx,
                       RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                       const std::vector<Pointer<RefineSchedule<NDIM> > >& /*f_prolongation_scheds*/,
                       const double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Determine the type of data centering.
    Pointer<hier::Variable<NDIM> > f_var;
    var_db->mapIndexToVariable(f_data_idx, f_var);
    Pointer<SideVariable<NDIM, double> > f_sc_var = f_var;
    const bool sc_data = f_sc_var;
    TBOX_ASSERT(sc_data);

    // Make a copy of the Eulerian data.
    const int f_copy_data_idx = var_db->registerClonedPatchDataIndex(f_var, f_data_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(f_copy_data_idx);
    }
    Pointer<HierarchyDataOpsReal<NDIM, double> > f_data_ops =
        HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(f_var, d_hierarchy, true);
    f_data_ops->swapData(f_copy_data_idx, f_data_idx);
    f_data_ops->setToScalar(f_data_idx, 0.0, /*interior_only*/ false);

    // Synchronize Lagrangian values.
    std::vector<Pointer<LData> >* X_data;
    bool* X_needs_ghost_fill;
    getPositionData(&X_data, &X_needs_ghost_fill, data_time);
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        if (*X_needs_ghost_fill && (*X_data)[ln]) (*X_data)[ln]->beginGhostUpdate();
        if (d_tau_data[ln]) d_tau_data[ln]->beginGhostUpdate();
    }
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        if (*X_needs_ghost_fill && (*X_data)[ln]) (*X_data)[ln]->endGhostUpdate();
        if (d_tau_data[ln]) d_tau_data[ln]->endGhostUpdate();
    }
    *X_needs_ghost_fill = false;

    // Spread data from the Lagrangian mesh to the Eulerian grid.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        boost::multi_array_ref<double, 2>& X_array = *(*X_data)[ln]->getGhostedLocalFormVecArray();
        boost::multi_array_ref<double, 2>& tau_array = *d_tau_data[ln]->getGhostedLocalFormVecArray();
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_data_idx);
            Pointer<LNodeSetData> idx_data = patch->getPatchData(d_l_data_manager->getLNodePatchDescriptorIndex());
            const Box<NDIM>& patch_box = patch->getBox();
            Box<NDIM> side_boxes[NDIM];
            for (unsigned int d = 0; d < NDIM; ++d) side_boxes[d] = SideGeometry<NDIM>::toSideBox(patch_box, d);
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const x_lower = patch_geom->getXLower();
            const double* const x_upper = patch_geom->getXUpper();
            const double* const dx = patch_geom->getDx();
            double dV_c = 1.0;
            for (unsigned int d = 0; d < NDIM; ++d) dV_c *= dx[d];
            for (LNodeSetData::CellIterator it(idx_data->getGhostBox()); it; it++)
            {
                const Index<NDIM>& i = *it;
                LNodeSet* const node_set = idx_data->getItem(i);
                if (!node_set) continue;
                for (LNodeSet::iterator it = node_set->begin(); it != node_set->end(); ++it)
                {
                    const LNode* const node_idx = *it;
                    MaterialPointSpec* mp_spec = node_idx->getNodeDataItem<MaterialPointSpec>();
                    if (!mp_spec) continue;
                    const double wgt = mp_spec->getWeight();
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const double* const X = &X_array[local_idx][0];
                    TensorValue<double> tau;
                    for (int i = 0; i < NDIM; ++i)
                    {
                        for (int j = 0; j < NDIM; ++j)
                        {
                            tau(i, j) = tau_array[local_idx][NDIM * i + j];
                        }
                    }

                    // Weight tau using a smooth kernel function evaluated about
                    // X.
                    Box<NDIM> stencil_box;
                    const int stencil_size = LEInteractor::getStencilSize(KERNEL_FCN);
                    boost::array<boost::multi_array<double, 1>, NDIM> phi, dphi;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        phi[d].resize(boost::extents[stencil_size]);
                        dphi[d].resize(boost::extents[stencil_size]);
                    }
                    for (unsigned int component = 0; component < NDIM; ++component)
                    {
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            kernel(X[d],
                                   x_lower[d] + (d == component ? -0.5 * dx[d] : 0.0),
                                   x_upper[d] + (d == component ? +0.5 * dx[d] : 0.0),
                                   dx[d],
                                   patch_box.lower(d),
                                   patch_box.upper(d) + (d == component ? 1 : 0),
                                   stencil_box.lower(d),
                                   stencil_box.upper(d),
                                   phi[d],
                                   dphi[d]);
                        }
                        for (Box<NDIM>::Iterator b(stencil_box * side_boxes[component]); b; b++)
                        {
                            const Index<NDIM>& i = b();
                            const Index<NDIM> i_shift = i - stencil_box.lower();
                            const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);
                            double f = 0.0;
                            for (unsigned int k = 0; k < NDIM; ++k)
                            {
                                double dw_dx_k = 1.0;
                                for (unsigned int d = 0; d < NDIM; ++d)
                                {
                                    if (d == k)
                                    {
                                        dw_dx_k *= dphi[d][i_shift(d)] / dx[d];
                                    }
                                    else
                                    {
                                        dw_dx_k *= phi[d][i_shift(d)];
                                    }
                                }
                                f += tau(component, k) * dw_dx_k;
                            }
                            (*f_data)(i_s) += f * wgt / dV_c;
                        }
                    }
                }
            }
            if (f_phys_bdry_op)
            {
                f_phys_bdry_op->setPatchDataIndex(f_data_idx);
                f_phys_bdry_op->accumulateFromPhysicalBoundaryData(*patch, data_time, f_data->getGhostCellWidth());
            }
        }
    }

    // Accumulate data.
    f_data_ops->swapData(f_copy_data_idx, f_data_idx);
    f_data_ops->add(f_data_idx, f_data_idx, f_copy_data_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(f_copy_data_idx);
    }
    var_db->removePatchDataIndex(f_copy_data_idx);
    return;
} // spreadForce

void
IMPMethod::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                    Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
                                    int /*u_data_idx*/,
                                    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
                                    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*u_ghost_fill_scheds*/,
                                    int /*integrator_step*/,
                                    double /*init_data_time*/,
                                    bool initial_time)
{
    // Cache pointers to the patch hierarchy and gridding algorithm.
    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;

    // Initialize various Lagrangian data objects.
    if (initial_time)
    {
        pout << "WARNING: IMPMethod implementation currently requires that the initial "
                "velocity is "
                "*zero*.\n";
    }
    return;
} // initializePatchHierarchy

void
IMPMethod::registerPK1StressTensorFunction(PK1StressFcnPtr PK1_stress_fcn, void* PK1_stress_fcn_ctx)
{
    d_PK1_stress_fcn = PK1_stress_fcn;
    d_PK1_stress_fcn_ctx = PK1_stress_fcn_ctx;
    return;
} // registerPK1StressTensorFunction

void
IMPMethod::registerLoadBalancer(Pointer<LoadBalancer<NDIM> > load_balancer, int workload_data_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(load_balancer);
#endif
    d_load_balancer = load_balancer;
    d_workload_idx = workload_data_idx;
    d_l_data_manager->registerLoadBalancer(load_balancer, workload_data_idx);
    return;
} // registerLoadBalancer

void
IMPMethod::updateWorkloadEstimates(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/, int /*workload_data_idx*/)
{
    d_l_data_manager->updateWorkloadEstimates();
    return;
} // updateWorkloadEstimates

void IMPMethod::beginDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                        Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    d_l_data_manager->beginDataRedistribution();
    return;
} // beginDataRedistribution

void IMPMethod::endDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                      Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    d_l_data_manager->endDataRedistribution();
    return;
} // endDataRedistribution

void
IMPMethod::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                               int level_number,
                               double init_data_time,
                               bool can_be_refined,
                               bool initial_time,
                               Pointer<BasePatchLevel<NDIM> > old_level,
                               bool allocate_data)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    d_l_data_manager->setPatchHierarchy(hierarchy);
    d_l_data_manager->setPatchLevels(0, finest_hier_level);
    d_l_data_manager->initializeLevelData(
        hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
    if (initial_time && d_l_data_manager->levelContainsLagrangianData(level_number))
    {
        Pointer<LData> Grad_U_data =
            d_l_data_manager->createLData("Grad_U", level_number, NDIM * NDIM, /*manage_data*/ true);
        Pointer<LData> F_data = d_l_data_manager->createLData("F",
                                                              level_number,
                                                              NDIM * NDIM,
                                                              /*manage_data*/ true);
        Pointer<LData> tau_data = d_l_data_manager->createLData("tau", level_number, NDIM * NDIM, /*manage_data*/ true);
        if (d_silo_writer)
        {
            d_silo_writer->registerVariableData("F0", F_data, 0 * NDIM, NDIM, level_number);
            d_silo_writer->registerVariableData("F1", F_data, 1 * NDIM, NDIM, level_number);
#if (NDIM == 3)
            d_silo_writer->registerVariableData("F2", F_data, 2 * NDIM, NDIM, level_number);
#endif
            d_silo_writer->registerVariableData("tau0", tau_data, 0 * NDIM, NDIM, level_number);
            d_silo_writer->registerVariableData("tau1", tau_data, 1 * NDIM, NDIM, level_number);
#if (NDIM == 3)
            d_silo_writer->registerVariableData("tau2", tau_data, 2 * NDIM, NDIM, level_number);
#endif
        }

        // Initialize the deformation gradient and Kirchhoff stress.
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(level_number);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        boost::multi_array_ref<double, 2>& F_array = *F_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& tau_array = *tau_data->getLocalFormVecArray();
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int idx = node_idx->getLocalPETScIndex();
            for (int i = 0; i < NDIM; ++i)
            {
                for (int j = 0; j < NDIM; ++j)
                {
                    F_array[idx][NDIM * i + j] = (i == j ? 1.0 : 0.0);
                    tau_array[idx][NDIM * i + j] = 0.0;
                }
            }
        }
    }
    if (d_load_balancer && d_l_data_manager->levelContainsLagrangianData(level_number))
    {
        d_load_balancer->setWorkloadPatchDataIndex(d_workload_idx, level_number);
        d_l_data_manager->updateWorkloadEstimates(level_number, level_number);
    }
    return;
} // initializeLevelData

void
IMPMethod::resetHierarchyConfiguration(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                       int coarsest_level,
                                       int finest_level)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    d_l_data_manager->setPatchHierarchy(hierarchy);
    d_l_data_manager->setPatchLevels(0, finest_hier_level);
    d_l_data_manager->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    return;
} // resetHierarchyConfiguration

void
IMPMethod::applyGradientDetector(Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
                                 int level_number,
                                 double error_data_time,
                                 int tag_index,
                                 bool initial_time,
                                 bool uses_richardson_extrapolation_too)
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
#endif
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Tag cells that contain Lagrangian nodes.
    d_l_data_manager->applyGradientDetector(
        hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    return;
} // applyGradientDetector

void
IMPMethod::putToDatabase(Pointer<Database> db)
{
    db->putInteger("IMP_METHOD_VERSION", IMP_METHOD_VERSION);
    db->putIntegerArray("d_ghosts", d_ghosts, NDIM);
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IMPMethod::getPositionData(std::vector<Pointer<LData> >** X_data, bool** X_needs_ghost_fill, double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        *X_data = &d_X_current_data;
        *X_needs_ghost_fill = &d_X_current_needs_ghost_fill;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
            if (!d_X_half_data[ln])
            {
                d_X_half_data[ln] = d_l_data_manager->createLData("X_half", ln, NDIM);
                d_X_half_needs_reinit = true;
            }
        }
        if (d_X_half_needs_reinit)
        {
            reinitMidpointData(d_X_current_data, d_X_new_data, d_X_half_data);
            d_X_half_needs_reinit = false;
            d_X_half_needs_ghost_fill = true;
        }
        *X_data = &d_X_half_data;
        *X_needs_ghost_fill = &d_X_half_needs_ghost_fill;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        *X_data = &d_X_new_data;
        *X_needs_ghost_fill = &d_X_new_needs_ghost_fill;
    }
    return;
} // getPositionData

void
IMPMethod::getVelocityData(std::vector<Pointer<LData> >** U_data,
                           std::vector<Pointer<LData> >** Grad_U_data,
                           double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        *U_data = &d_U_current_data;
        *Grad_U_data = &d_Grad_U_current_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
            if (!d_U_half_data[ln])
            {
                d_U_half_data[ln] = d_l_data_manager->createLData("U_half", ln, NDIM);
                d_Grad_U_half_data[ln] = d_l_data_manager->createLData("Grad_U_half", ln, NDIM * NDIM);
                d_U_half_needs_reinit = true;
            }
        }
        if (d_U_half_needs_reinit)
        {
            reinitMidpointData(d_U_current_data, d_U_new_data, d_U_half_data);
            reinitMidpointData(d_Grad_U_current_data, d_Grad_U_new_data, d_Grad_U_half_data);
            d_U_half_needs_reinit = false;
        }
        *U_data = &d_U_half_data;
        *Grad_U_data = &d_Grad_U_half_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        *U_data = &d_U_new_data;
        *Grad_U_data = &d_Grad_U_new_data;
    }
    return;
} // getVelocityData

void
IMPMethod::getDeformationGradientData(std::vector<Pointer<LData> >** F_data, double data_time)
{
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        *F_data = &d_F_current_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        *F_data = &d_F_half_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        *F_data = &d_F_new_data;
    }
    return;
} // getDeformationGradientData

void
IMPMethod::reinitMidpointData(const std::vector<Pointer<LData> >& current_data,
                              const std::vector<Pointer<LData> >& new_data,
                              const std::vector<Pointer<LData> >& half_data)
{
    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        ierr = VecAXPBYPCZ(half_data[ln]->getVec(), 0.5, 0.5, 0.0, current_data[ln]->getVec(), new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
    }
    return;
} // reinitMidpointData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IMPMethod::getFromInput(Pointer<Database> db, bool is_from_restart)
{
    if (!is_from_restart)
    {
        if (db->isInteger("min_ghost_cell_width"))
        {
            d_ghosts = db->getInteger("min_ghost_cell_width");
        }
        else if (db->isDouble("min_ghost_cell_width"))
        {
            d_ghosts = static_cast<int>(std::ceil(db->getDouble("min_ghost_cell_width")));
        }
    }
    if (db->keyExists("do_log"))
        d_do_log = db->getBool("do_log");
    else if (db->keyExists("enable_logging"))
        d_do_log = db->getBool("enable_logging");
    return;
} // getFromInput

void
IMPMethod::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file."
                                 << std::endl);
    }
    int ver = db->getInteger("IMP_METHOD_VERSION");
    if (ver != IMP_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    db->getIntegerArray("d_ghosts", d_ghosts, NDIM);
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
