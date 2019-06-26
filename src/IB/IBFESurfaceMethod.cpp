// Filename: IBFESurfaceMethod.cpp
// Created on 19 May 2018 by Boyce Griffith
// Copyright (c) 2002-2017, Boyce Griffith
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

#include "ibamr/IBFESurfaceMethod.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/FEDataInterpolation.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_macros.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
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
#include "SideData.h"
#include "SideIndex.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

#include "libmesh/boundary_info.h"
#include "libmesh/compare_types.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/edge.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_type.h"
#include "libmesh/fem_context.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/system.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_tensor.h"
#include "libmesh/type_vector.h"
#include "libmesh/variant_filter_iterator.h"
#include "libmesh/vector_value.h"

#include "petscvec.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include "boost/math/special_functions/round.hpp"
#include "boost/multi_array.hpp"
IBTK_ENABLE_EXTRA_WARNINGS

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace SAMRAI
{
namespace xfer
{
template <int DIM>
class RefineSchedule;
template <int DIM>
class CoarsenSchedule;
} // namespace xfer
} // namespace SAMRAI

using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of IBFESurfaceMethod restart file data.
static const int IBFE_METHOD_VERSION = 1;

std::string
libmesh_restart_file_name(const std::string& restart_dump_dirname,
                          unsigned int time_step_number,
                          unsigned int part,
                          const std::string& extension)
{
    std::ostringstream file_name_prefix;
    file_name_prefix << restart_dump_dirname << "/libmesh_data_part_" << part << "." << std::setw(6)
                     << std::setfill('0') << std::right << time_step_number << "." << extension;
    return file_name_prefix.str();
}
} // namespace

const std::string IBFESurfaceMethod::COORDS_SYSTEM_NAME = "IB coordinates system";
const std::string IBFESurfaceMethod::COORD_MAPPING_SYSTEM_NAME = "IB coordinate mapping system";
const std::string IBFESurfaceMethod::FORCE_SYSTEM_NAME = "IB force system";
const std::string IBFESurfaceMethod::VELOCITY_SYSTEM_NAME = "IB velocity system";
const std::string IBFESurfaceMethod::NORMAL_VELOCITY_SYSTEM_NAME = "IB normal velocity system";
const std::string IBFESurfaceMethod::TANGENTIAL_VELOCITY_SYSTEM_NAME = "IB tangential velocity system";
const std::string IBFESurfaceMethod::PRESSURE_JUMP_SYSTEM_NAME = "IB [[p]] system";
const std::string IBFESurfaceMethod::WSS_SYSTEM_NAME = "One sided interior/exterior shear stress system";
const std::string IBFESurfaceMethod::P_SYSTEM_NAME = "One sided interior/exterior pressure system";
const std::string IBFESurfaceMethod::TAU_SYSTEM_NAME = "traction system";
const std::array<std::string, NDIM> IBFESurfaceMethod::VELOCITY_JUMP_SYSTEM_NAME = {
    { "IB velocity [[du]] jump system",
      "IB velocity [[dv]] jump system"
#if (NDIM == 3)
      ,
      "IB velocity [[dw]] jump system"
#endif
    }
};

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFESurfaceMethod::IBFESurfaceMethod(const std::string& object_name,
                                     Pointer<Database> input_db,
                                     MeshBase* mesh,
                                     int max_level_number,
                                     bool register_for_restart,
                                     const std::string& restart_read_dirname,
                                     unsigned int restart_restore_number)
{
    commonConstructor(object_name,
                      input_db,
                      std::vector<MeshBase*>(d_num_parts, mesh),
                      max_level_number,
                      register_for_restart,
                      restart_read_dirname,
                      restart_restore_number);
    return;
} // IBFESurfaceMethod

IBFESurfaceMethod::IBFESurfaceMethod(const std::string& object_name,
                                     Pointer<Database> input_db,
                                     const std::vector<MeshBase*>& meshes,
                                     int max_level_number,
                                     bool register_for_restart,
                                     const std::string& restart_read_dirname,
                                     unsigned int restart_restore_number)
    : d_num_parts(static_cast<int>(meshes.size()))
{
    commonConstructor(object_name,
                      input_db,
                      meshes,
                      max_level_number,
                      register_for_restart,
                      restart_read_dirname,
                      restart_restore_number);
    return;
} // IBFESurfaceMethod

IBFESurfaceMethod::~IBFESurfaceMethod()
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        delete d_equation_systems[part];
    }
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
} // ~IBFESurfaceMethod

FEDataManager*
IBFESurfaceMethod::getFEDataManager(const unsigned int part) const
{
    TBOX_ASSERT(d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_num_parts);
    return d_fe_data_managers[part];
} // getFEDataManager

void
IBFESurfaceMethod::registerInitialCoordinateMappingFunction(const CoordinateMappingFcnData& data,
                                                            const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_coordinate_mapping_fcn_data[part] = data;
    return;
} // registerInitialCoordinateMappingFunction

void
IBFESurfaceMethod::registerInitialVelocityFunction(const InitialVelocityFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_initial_velocity_fcn_data[part] = data;
    return;
} // registerInitialVelocityFunction

void
IBFESurfaceMethod::registerLagSurfacePressureFunction(const LagSurfacePressureFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_surface_pressure_fcn_data[part] = data;
    return;
} // registerLagSurfacePressureFunction

void
IBFESurfaceMethod::registerLagSurfaceForceFunction(const LagSurfaceForceFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_surface_force_fcn_data[part] = data;
    return;
} // registerLagSurfaceForceFunction

const VectorValue<double>&
IBFESurfaceMethod::getSurfaceForceIntegral(const unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_lag_surface_force_integral[part];
}

const IntVector<NDIM>&
IBFESurfaceMethod::getMinimumGhostCellWidth() const
{
    return d_ghosts;
} // getMinimumGhostCellWidth

void
IBFESurfaceMethod::setupTagBuffer(Array<int>& tag_buffer, Pointer<GriddingAlgorithm<NDIM> > gridding_alg) const
{
    const int finest_hier_ln = gridding_alg->getMaxLevels() - 1;
    const int tsize = tag_buffer.size();
    tag_buffer.resizeArray(finest_hier_ln);
    for (int i = tsize; i < finest_hier_ln; ++i) tag_buffer[i] = 0;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        const int gcw = d_fe_data_managers[part]->getGhostCellWidth().max();
        const int tag_ln = d_fe_data_managers[part]->getLevelNumber() - 1;
        if (tag_ln >= 0 && tag_ln < finest_hier_ln)
        {
            tag_buffer[tag_ln] = std::max(tag_buffer[tag_ln], gcw);
        }
    }
    for (int ln = finest_hier_ln - 2; ln >= 0; --ln)
    {
        tag_buffer[ln] =
            std::max(tag_buffer[ln], tag_buffer[ln + 1] / gridding_alg->getRatioToCoarserLevel(ln + 1).max() + 1);
    }
    return;
} // setupTagBuffer

void
IBFESurfaceMethod::preprocessIntegrateData(double current_time, double new_time, int /*num_cycles*/)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_half_time = current_time + 0.5 * (new_time - current_time);

    // Extract the FE data.
    d_X_systems.resize(d_num_parts);
    d_X_current_vecs.resize(d_num_parts);
    d_X_new_vecs.resize(d_num_parts);
    d_X_half_vecs.resize(d_num_parts);
    d_X_IB_ghost_vecs.resize(d_num_parts);

    d_U_systems.resize(d_num_parts);
    d_U_current_vecs.resize(d_num_parts);
    d_U_new_vecs.resize(d_num_parts);
    d_U_half_vecs.resize(d_num_parts);

    d_U_n_systems.resize(d_num_parts);
    d_U_n_current_vecs.resize(d_num_parts);
    d_U_n_new_vecs.resize(d_num_parts);
    d_U_n_half_vecs.resize(d_num_parts);

    d_U_t_systems.resize(d_num_parts);
    d_U_t_current_vecs.resize(d_num_parts);
    d_U_t_new_vecs.resize(d_num_parts);
    d_U_t_half_vecs.resize(d_num_parts);

    d_F_systems.resize(d_num_parts);
    d_F_half_vecs.resize(d_num_parts);
    d_F_IB_ghost_vecs.resize(d_num_parts);

    d_P_jump_systems.resize(d_num_parts);
    d_P_jump_half_vecs.resize(d_num_parts);
    d_P_jump_IB_ghost_vecs.resize(d_num_parts);

    d_P_systems.resize(d_num_parts);
    d_P_half_vecs.resize(d_num_parts);
    d_P_IB_ghost_vecs.resize(d_num_parts);

    d_DU_jump_systems.resize(d_num_parts);
    d_DU_jump_half_vecs.resize(d_num_parts);
    d_DU_jump_IB_ghost_vecs.resize(d_num_parts);

    d_WSS_systems.resize(d_num_parts);
    d_WSS_half_vecs.resize(d_num_parts);
    d_WSS_IB_ghost_vecs.resize(d_num_parts);

    d_TAU_systems.resize(d_num_parts);
    d_TAU_half_vecs.resize(d_num_parts);
    d_TAU_IB_ghost_vecs.resize(d_num_parts);

    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_X_systems[part] = &d_equation_systems[part]->get_system(COORDS_SYSTEM_NAME);
        d_X_current_vecs[part] = dynamic_cast<PetscVector<double>*>(d_X_systems[part]->current_local_solution.get());
        d_X_new_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_X_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
        d_X_half_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_X_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
        d_X_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedCoordsVector(/*localize_data*/ false));

        d_U_systems[part] = &d_equation_systems[part]->get_system(VELOCITY_SYSTEM_NAME);
        d_U_current_vecs[part] = dynamic_cast<PetscVector<double>*>(d_U_systems[part]->current_local_solution.get());
        d_U_new_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
        d_U_half_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_current_vecs[part]->clone().release()); // WARNING: must be manually deleted

        d_U_n_systems[part] = &d_equation_systems[part]->get_system(NORMAL_VELOCITY_SYSTEM_NAME);
        d_U_n_current_vecs[part] =
            dynamic_cast<PetscVector<double>*>(d_U_n_systems[part]->current_local_solution.get());
        d_U_n_new_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_n_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
        d_U_n_half_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_n_current_vecs[part]->clone().release()); // WARNING: must be manually deleted

        d_U_t_systems[part] = &d_equation_systems[part]->get_system(TANGENTIAL_VELOCITY_SYSTEM_NAME);
        d_U_t_current_vecs[part] =
            dynamic_cast<PetscVector<double>*>(d_U_t_systems[part]->current_local_solution.get());
        d_U_t_new_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_t_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
        d_U_t_half_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_t_current_vecs[part]->clone().release()); // WARNING: must be manually deleted

        d_F_systems[part] = &d_equation_systems[part]->get_system(FORCE_SYSTEM_NAME);
        d_F_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_F_systems[part]->current_local_solution.get());
        d_F_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(FORCE_SYSTEM_NAME, /*localize_data*/ false));

        if (d_use_pressure_jump_conditions)
        {
            d_P_jump_systems[part] = &d_equation_systems[part]->get_system(PRESSURE_JUMP_SYSTEM_NAME);
            d_P_jump_half_vecs[part] =
                dynamic_cast<PetscVector<double>*>(d_P_jump_systems[part]->current_local_solution.get());
            d_P_jump_IB_ghost_vecs[part] =
                dynamic_cast<PetscVector<double>*>(d_fe_data_managers[part]->buildGhostedSolutionVector(
                    PRESSURE_JUMP_SYSTEM_NAME, /*localize_data*/ false));

            d_P_systems[part] = &d_equation_systems[part]->get_system(P_SYSTEM_NAME);
            d_P_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_P_systems[part]->current_local_solution.get());
            d_P_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
                d_fe_data_managers[part]->buildGhostedSolutionVector(P_SYSTEM_NAME, /*localize_data*/ false));
        }

        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                d_DU_jump_systems[part][d] = &d_equation_systems[part]->get_system(VELOCITY_JUMP_SYSTEM_NAME[d]);
                d_DU_jump_half_vecs[part][d] =
                    dynamic_cast<PetscVector<double>*>(d_DU_jump_systems[part][d]->current_local_solution.get());
                d_DU_jump_IB_ghost_vecs[part][d] =
                    dynamic_cast<PetscVector<double>*>(d_fe_data_managers[part]->buildGhostedSolutionVector(
                        VELOCITY_JUMP_SYSTEM_NAME[d], /*localize_data*/ false));
            }

            d_WSS_systems[part] = &d_equation_systems[part]->get_system(WSS_SYSTEM_NAME);
            d_WSS_half_vecs[part] =
                dynamic_cast<PetscVector<double>*>(d_WSS_systems[part]->current_local_solution.get());
            d_WSS_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
                d_fe_data_managers[part]->buildGhostedSolutionVector(WSS_SYSTEM_NAME, /*localize_data*/ false));
        }
        if (d_use_velocity_jump_conditions && d_use_pressure_jump_conditions)
        {
            d_TAU_systems[part] = &d_equation_systems[part]->get_system(TAU_SYSTEM_NAME);
            d_TAU_half_vecs[part] =
                dynamic_cast<PetscVector<double>*>(d_TAU_systems[part]->current_local_solution.get());
            d_TAU_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
                d_fe_data_managers[part]->buildGhostedSolutionVector(TAU_SYSTEM_NAME, /*localize_data*/ false));
        }

        // Initialize X^{n+1/2} and X^{n+1} to equal X^{n}, and initialize
        // U^{n+1/2} and U^{n+1} to equal U^{n}.
        *d_X_current_vecs[part] = *d_X_systems[part]->solution;
        *d_X_new_vecs[part] = *d_X_current_vecs[part];
        *d_X_half_vecs[part] = *d_X_current_vecs[part];

        *d_U_current_vecs[part] = *d_U_systems[part]->solution;
        *d_U_new_vecs[part] = *d_U_current_vecs[part];
        *d_U_half_vecs[part] = *d_U_current_vecs[part];

        *d_U_n_current_vecs[part] = *d_U_n_systems[part]->solution;
        *d_U_n_new_vecs[part] = *d_U_n_current_vecs[part];
        *d_U_n_half_vecs[part] = *d_U_n_current_vecs[part];

        *d_U_t_current_vecs[part] = *d_U_t_systems[part]->solution;
        *d_U_t_new_vecs[part] = *d_U_t_current_vecs[part];
        *d_U_t_half_vecs[part] = *d_U_t_current_vecs[part];

        *d_F_half_vecs[part] = *d_F_systems[part]->solution;

        if (d_use_pressure_jump_conditions)
        {
            *d_P_jump_half_vecs[part] = *d_P_jump_systems[part]->solution;
        }

        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                *d_DU_jump_half_vecs[part][d] = *d_DU_jump_systems[part][d]->solution;
            }

            *d_WSS_half_vecs[part] = *d_WSS_systems[part]->solution;
        }
        if (d_use_velocity_jump_conditions && d_use_pressure_jump_conditions)
        {
            *d_TAU_half_vecs[part] = *d_TAU_systems[part]->solution;
        }
    }
    return;
} // preprocessIntegrateData

void
IBFESurfaceMethod::calculateInterfacialFluidForces(double data_time)
{
    batch_vec_ghost_update(
        { d_X_new_vecs, d_U_new_vecs, d_U_n_new_vecs, d_U_t_new_vecs, d_F_half_vecs }, INSERT_VALUES, SCATTER_FORWARD);
    batch_vec_ghost_update(
        { d_WSS_half_vecs, d_P_half_vecs, d_P_jump_half_vecs, d_TAU_half_vecs }, INSERT_VALUES, SCATTER_FORWARD);

    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        if ((!d_use_pressure_jump_conditions || !d_use_velocity_jump_conditions) && d_compute_fluid_traction)
        {
            TBOX_ERROR(
                d_object_name << ": To compute the traction both velocity and preussure jumps need to be turned on!"
                              << std::endl);
        }

        if (d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
        {
            const int coarsest_ln = 0;
            const int finest_ln = d_hierarchy->getFinestLevelNumber();
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            const int p_data_idx = var_db->mapVariableAndContextToIndex(
                getINSHierarchyIntegrator()->getPressureVariable(), getINSHierarchyIntegrator()->getScratchContext());
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                const Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                if (!level->checkAllocated(mask_scratch_idx)) level->allocatePatchData(mask_scratch_idx);
            }

            HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(d_hierarchy, coarsest_ln, finest_ln);
            hier_cc_data_ops.copyData(mask_scratch_idx, p_data_idx, /*interior only*/ false);

            RefineAlgorithm<NDIM> ghost_fill_alg_p;
            ghost_fill_alg_p.registerRefine(mask_scratch_idx, mask_scratch_idx, mask_scratch_idx, NULL);
            Pointer<RefineSchedule<NDIM> > ghost_fill_schd_p =
                ghost_fill_alg_p.createSchedule(d_hierarchy->getPatchLevel(finest_ln));
            ghost_fill_schd_p->fillData(data_time);
            interpolatePressureForTraction(mask_scratch_idx, data_time, part);

            if (d_compute_fluid_traction && d_traction_activation_time <= data_time)
            {
                computeFluidTraction(data_time, part);
            }
        }
    }

} // calculateInterfacialFluidForces

void
IBFESurfaceMethod::postprocessIntegrateData(double /*current_time*/, double /*new_time*/, int /*num_cycles*/)
{
    batch_vec_ghost_update(
        { d_X_new_vecs, d_U_new_vecs, d_U_n_new_vecs, d_U_t_new_vecs, d_F_half_vecs }, INSERT_VALUES, SCATTER_FORWARD);
    batch_vec_ghost_update(
        { d_WSS_half_vecs, d_P_half_vecs, d_P_jump_half_vecs, d_TAU_half_vecs }, INSERT_VALUES, SCATTER_FORWARD);

    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        if ((!d_use_pressure_jump_conditions || !d_use_velocity_jump_conditions) && d_compute_fluid_traction)
        {
            TBOX_ERROR(
                d_object_name << ": To compute the traction both velocity and preussure jumps need to be turned on!"
                              << std::endl);
        }

        if (d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
        {
            const int coarsest_ln = 0;
            const int finest_ln = d_hierarchy->getFinestLevelNumber();
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            const int p_data_idx = var_db->mapVariableAndContextToIndex(
                getINSHierarchyIntegrator()->getPressureVariable(), getINSHierarchyIntegrator()->getScratchContext());
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                const Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                if (!level->checkAllocated(mask_scratch_idx)) level->allocatePatchData(mask_scratch_idx);
            }

            HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(d_hierarchy, coarsest_ln, finest_ln);
            hier_cc_data_ops.copyData(mask_scratch_idx, p_data_idx, /*interior only*/ false);

            RefineAlgorithm<NDIM> ghost_fill_alg_p;
            ghost_fill_alg_p.registerRefine(mask_scratch_idx, mask_scratch_idx, mask_scratch_idx, NULL);
            Pointer<RefineSchedule<NDIM> > ghost_fill_schd_p =
                ghost_fill_alg_p.createSchedule(d_hierarchy->getPatchLevel(finest_ln));
            ghost_fill_schd_p->fillData(d_new_time);
            interpolatePressureForTraction(mask_scratch_idx, d_new_time, part);

            if (d_compute_fluid_traction && d_traction_activation_time <= d_current_time)
            {
                computeFluidTraction(d_new_time, part);
            }
        }

        // Reset time-dependent Lagrangian data.
        *d_X_systems[part]->solution = *d_X_new_vecs[part];
        *d_X_systems[part]->current_local_solution = *d_X_new_vecs[part];
        delete d_X_new_vecs[part];
        delete d_X_half_vecs[part];

        *d_U_systems[part]->solution = *d_U_new_vecs[part];
        *d_U_systems[part]->current_local_solution = *d_U_new_vecs[part];
        delete d_U_new_vecs[part];
        delete d_U_half_vecs[part];

        *d_U_n_systems[part]->solution = *d_U_n_new_vecs[part];
        *d_U_n_systems[part]->current_local_solution = *d_U_n_new_vecs[part];
        delete d_U_n_new_vecs[part];
        delete d_U_n_half_vecs[part];

        *d_U_t_systems[part]->solution = *d_U_t_new_vecs[part];
        *d_U_t_systems[part]->current_local_solution = *d_U_t_new_vecs[part];
        delete d_U_t_new_vecs[part];
        delete d_U_t_half_vecs[part];

        *d_F_systems[part]->solution = *d_F_half_vecs[part];
        *d_F_systems[part]->current_local_solution = *d_F_half_vecs[part];

        if (d_use_pressure_jump_conditions)
        {
            *d_P_jump_systems[part]->solution = *d_P_jump_half_vecs[part];
            *d_P_jump_systems[part]->current_local_solution = *d_P_jump_half_vecs[part];

            *d_P_systems[part]->solution = *d_P_half_vecs[part];
            *d_P_systems[part]->current_local_solution = *d_P_half_vecs[part];
        }

        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                *d_DU_jump_systems[part][d]->solution = *d_DU_jump_half_vecs[part][d];
                *d_DU_jump_systems[part][d]->current_local_solution = *d_DU_jump_half_vecs[part][d];
            }
            *d_WSS_systems[part]->solution = *d_WSS_half_vecs[part];
            *d_WSS_systems[part]->current_local_solution = *d_WSS_half_vecs[part];
        }

        if (d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
        {
            *d_TAU_systems[part]->solution = *d_TAU_half_vecs[part];
            *d_TAU_systems[part]->current_local_solution = *d_TAU_half_vecs[part];
        }

        // Update the coordinate mapping dX = X - s.
        updateCoordinateMapping(part);
    }

    d_X_systems.clear();
    d_X_current_vecs.clear();
    d_X_new_vecs.clear();
    d_X_half_vecs.clear();
    d_X_IB_ghost_vecs.clear();

    d_U_systems.clear();
    d_U_current_vecs.clear();
    d_U_new_vecs.clear();
    d_U_half_vecs.clear();

    d_U_n_systems.clear();
    d_U_n_current_vecs.clear();
    d_U_n_new_vecs.clear();
    d_U_n_half_vecs.clear();

    d_U_t_systems.clear();
    d_U_t_current_vecs.clear();
    d_U_t_new_vecs.clear();
    d_U_t_half_vecs.clear();

    d_F_systems.clear();
    d_F_half_vecs.clear();
    d_F_IB_ghost_vecs.clear();

    d_P_jump_systems.clear();
    d_P_jump_half_vecs.clear();
    d_P_jump_IB_ghost_vecs.clear();

    d_DU_jump_systems.clear();
    d_DU_jump_half_vecs.clear();
    d_DU_jump_IB_ghost_vecs.clear();

    d_WSS_systems.clear();
    d_WSS_half_vecs.clear();
    d_WSS_IB_ghost_vecs.clear();

    d_P_systems.clear();
    d_P_half_vecs.clear();
    d_P_IB_ghost_vecs.clear();

    d_TAU_systems.clear();
    d_TAU_half_vecs.clear();
    d_TAU_IB_ghost_vecs.clear();

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_half_time = std::numeric_limits<double>::quiet_NaN();
    return;
} // postprocessIntegrateData

void
IBFESurfaceMethod::interpolateVelocity(const int u_data_idx,
                                       const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
                                       const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                       const double data_time)
{
    const double mu = getINSHierarchyIntegrator()->getStokesSpecifications()->getMu();

    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        NumericVector<double>* U_vec = nullptr;
        NumericVector<double>* U_n_vec = nullptr;
        NumericVector<double>* U_t_vec = nullptr;
        NumericVector<double>* X_vec = nullptr;
        NumericVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
        const std::array<PetscVector<double>*, NDIM>& DU_jump_ghost_vec = d_DU_jump_IB_ghost_vecs[part];
        if (MathUtilities<double>::equalEps(data_time, d_current_time))
        {
            U_vec = d_U_current_vecs[part];
            U_n_vec = d_U_n_current_vecs[part];
            U_t_vec = d_U_t_current_vecs[part];
            X_vec = d_X_current_vecs[part];
        }
        else if (MathUtilities<double>::equalEps(data_time, d_half_time))
        {
            U_vec = d_U_half_vecs[part];
            U_n_vec = d_U_n_half_vecs[part];
            U_t_vec = d_U_t_half_vecs[part];
            X_vec = d_X_half_vecs[part];
        }
        else if (MathUtilities<double>::equalEps(data_time, d_new_time))
        {
            U_vec = d_U_new_vecs[part];
            U_n_vec = d_U_n_new_vecs[part];
            U_t_vec = d_U_t_new_vecs[part];
            X_vec = d_X_new_vecs[part];
        }
        copy_and_synch(*X_vec, *X_ghost_vec);

        NumericVector<double>* WSS_vec = d_WSS_half_vecs[part];

        // Extract the mesh.
        EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
        const MeshBase& mesh = equation_systems->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();
        std::unique_ptr<QBase> qrule;

        // Extract the FE systems and DOF maps, and setup the FE object.
        System& U_system = *d_U_systems[part];
        const DofMap& U_dof_map = U_system.get_dof_map();
        FEDataManager::SystemDofMapCache& U_dof_map_cache =
            *d_fe_data_managers[part]->getDofMapCache(VELOCITY_SYSTEM_NAME);
        FEType U_fe_type = U_dof_map.variable_type(0);
        for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(U_dof_map.variable_type(d) == U_fe_type);

        System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
        const DofMap& X_dof_map = X_system.get_dof_map();
        FEDataManager::SystemDofMapCache& X_dof_map_cache =
            *d_fe_data_managers[part]->getDofMapCache(COORDS_SYSTEM_NAME);
        FEType X_fe_type = X_dof_map.variable_type(0);
        for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
        TBOX_ASSERT(U_fe_type == X_fe_type);

        std::array<System*, NDIM> DU_jump_system;
        std::array<const DofMap*, NDIM> DU_jump_dof_map;
        std::array<FEDataManager::SystemDofMapCache*, NDIM> DU_jump_dof_map_cache;
        FEType DU_jump_fe_type = INVALID_FE;
        if (d_use_velocity_jump_conditions)
        {
            for (unsigned d = 0; d < NDIM; ++d)
            {
                DU_jump_system[d] = &equation_systems->get_system(VELOCITY_JUMP_SYSTEM_NAME[d]);
                DU_jump_dof_map[d] = &DU_jump_system[d]->get_dof_map();
                DU_jump_dof_map_cache[d] = d_fe_data_managers[part]->getDofMapCache(VELOCITY_JUMP_SYSTEM_NAME[d]);
            }
            DU_jump_fe_type = DU_jump_system[0]->variable_type(0);
            for (unsigned i = 0; i < NDIM; ++i)
            {
                for (unsigned j = 0; j < NDIM; ++j)
                {
                    TBOX_ASSERT(DU_jump_fe_type == DU_jump_system[i]->variable_type(j));
                }
            }
        }

        System* WSS_system;
        const DofMap* WSS_dof_map;
        FEDataManager::SystemDofMapCache* WSS_dof_map_cache;
        FEType WSS_fe_type = INVALID_FE;
        if (d_use_velocity_jump_conditions)
        {
            WSS_system = &equation_systems->get_system(WSS_SYSTEM_NAME);
            WSS_dof_map = &WSS_system->get_dof_map();
            WSS_dof_map_cache = d_fe_data_managers[part]->getDofMapCache(WSS_SYSTEM_NAME);
            WSS_fe_type = WSS_dof_map->variable_type(0);
            for (unsigned i = 0; i < NDIM; ++i)
            {
                TBOX_ASSERT(WSS_dof_map->variable_type(i) == WSS_fe_type);
            }
        }
        TBOX_ASSERT(WSS_fe_type == DU_jump_fe_type);

        std::unique_ptr<FEBase> fe_X = FEBase::build(dim, X_fe_type);
        const std::vector<double>& JxW = fe_X->get_JxW();
        const std::vector<std::vector<double> >& phi_X = fe_X->get_phi();
        std::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi_X;
        dphi_dxi_X[0] = &fe_X->get_dphidxi();
        if (NDIM > 2) dphi_dxi_X[1] = &fe_X->get_dphideta();

        std::unique_ptr<FEBase> fe_DU_jump = FEBase::build(dim, DU_jump_fe_type);
        const std::vector<std::vector<double> >& phi_DU_jump = fe_DU_jump->get_phi();

        // Communicate any unsynchronized ghost data and extract the underlying
        // solution data.
        for (const auto& u_ghost_fill_sched : u_ghost_fill_scheds)
        {
            if (u_ghost_fill_sched) u_ghost_fill_sched->fillData(data_time);
        }
        X_ghost_vec->close();
        auto X_petsc_vec = static_cast<PetscVector<double>*>(X_ghost_vec);
        Vec X_global_vec = X_petsc_vec->vec();
        Vec X_local_vec;
        VecGhostGetLocalForm(X_global_vec, &X_local_vec);
        double* X_local_soln;
        VecGetArray(X_local_vec, &X_local_soln);
        std::unique_ptr<NumericVector<double> > X0_vec = X_petsc_vec->clone();
        copy_and_synch(X_system.get_vector("INITIAL_COORDINATES"), *X0_vec);

        // Loop over the patches to interpolate values to the element quadrature
        // points from the grid, then use these values to compute the projection
        // of the interpolated velocity field onto the FE basis functions.
        std::unique_ptr<NumericVector<double> > U_rhs_vec = U_vec->zero_clone();
        std::array<DenseVector<double>, NDIM> U_rhs_e;
        std::unique_ptr<NumericVector<double> > U_n_rhs_vec = U_n_vec->zero_clone();
        std::array<DenseVector<double>, NDIM> U_n_rhs_e;
        std::unique_ptr<NumericVector<double> > U_t_rhs_vec = U_t_vec->zero_clone();
        std::array<DenseVector<double>, NDIM> U_t_rhs_e;

        std::unique_ptr<NumericVector<double> > WSS_rhs_vec;
        std::array<DenseVector<double>, NDIM> WSS_rhs_e;
        if (d_use_velocity_jump_conditions)
        {
            WSS_rhs_vec = WSS_vec->zero_clone();
        }

        boost::multi_array<double, 2> x_node;
        std::array<boost::multi_array<double, 2>, NDIM> DU_jump_node;
        std::vector<double> U_qp, U_io_qp, WSS_qp, n_qp, x_qp, x_io_qp;
        std::array<std::vector<double>, NDIM> DU_jump_qp;
        VectorValue<double> U, WSS, U_n, U_t, n;
        std::array<VectorValue<double>, 2> dx_dxi;

        std::vector<libMesh::dof_id_type> dof_id_scratch;
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_fe_data_managers[part]->getLevelNumber());
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            // The relevant collection of elements.
            const std::vector<Elem*>& patch_elems =
                d_fe_data_managers[part]->getActivePatchElementMap()[local_patch_num];
            const size_t num_active_patch_elems = patch_elems.size();
            if (!num_active_patch_elems) continue;
            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();

            double diag_dis = 0.0;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                diag_dis += patch_dx[d] * patch_dx[d];
            }
            const double dh = d_wss_calc_width * sqrt(diag_dis);
            const int u_ghost_num = static_cast<int>(ceil(2.0 * dh / patch_dx_min));
            std::array<double, NDIM> x_lower_gh, x_upper_gh;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                x_lower_gh[d] = patch_x_lower[d] - (static_cast<double>(u_ghost_num)) * patch_dx[d];
                x_upper_gh[d] = patch_x_upper[d] + (static_cast<double>(u_ghost_num)) * patch_dx[d];
            }

            double* x_upper_ghost = &x_upper_gh[0];
            double* x_lower_ghost = &x_lower_gh[0];

            // Setup vectors to store the values of U, DU_j, x, and n at the
            // quadrature points.
            unsigned int n_qpoints_patch = 0;
            for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems[e_idx];
                const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                get_values_for_interpolation(x_node, *X_ghost_vec, X_dof_indices);
                FEDataManager::updateInterpQuadratureRule(qrule, d_default_interp_spec, elem, x_node, patch_dx_min);
                n_qpoints_patch += qrule->n_points();
            }

            if (!n_qpoints_patch) continue;
            U_qp.resize(NDIM * n_qpoints_patch);
            WSS_qp.resize(NDIM * n_qpoints_patch);
            U_io_qp.resize(NDIM * n_qpoints_patch);
            x_qp.resize(NDIM * n_qpoints_patch);
            x_io_qp.resize(NDIM * n_qpoints_patch);
            n_qp.resize(NDIM * n_qpoints_patch);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                DU_jump_qp[axis].resize(NDIM * n_qpoints_patch);
            }
            std::fill(U_qp.begin(), U_qp.end(), 0.0);
            std::fill(U_io_qp.begin(), U_io_qp.end(), 0.0);
            std::fill(WSS_qp.begin(), WSS_qp.end(), 0.0);

            // Loop over the elements and compute the positions of the quadrature points.
            qrule.reset();
            unsigned int qp_offset = 0;
            for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems[e_idx];
                const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                get_values_for_interpolation(x_node, *X_petsc_vec, X_local_soln, X_dof_indices);
                if (d_use_velocity_jump_conditions)
                {
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        const auto& DU_jump_dof_indices = DU_jump_dof_map_cache[axis]->dof_indices(elem);
                        get_values_for_interpolation(DU_jump_node[axis], *DU_jump_ghost_vec[axis], DU_jump_dof_indices);
                    }
                }
                const bool qrule_changed =
                    FEDataManager::updateInterpQuadratureRule(qrule, d_default_interp_spec, elem, x_node, patch_dx_min);
                if (qrule_changed) fe_X->attach_quadrature_rule(qrule.get());
                fe_X->reinit(elem);

                if (qrule_changed) fe_DU_jump->attach_quadrature_rule(qrule.get());
                fe_DU_jump->reinit(elem);

                const unsigned int n_nodes = elem->n_nodes();
                const unsigned int n_qpoints = qrule->n_points();

                // Zero out the values prior to accumulation.
                double* x_begin = &x_qp[NDIM * qp_offset];
                std::fill(x_begin, x_begin + NDIM * n_qpoints, 0.0);

                double* x_io_begin = &x_io_qp[NDIM * qp_offset];
                std::fill(x_io_begin, x_io_begin + NDIM * n_qpoints, 0.0);

                double* n_begin = &n_qp[NDIM * qp_offset];
                std::fill(n_begin, n_begin + NDIM * n_qpoints, 0.0);

                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    double* DU_jump_begin = &DU_jump_qp[axis][NDIM * qp_offset];
                    std::fill(DU_jump_begin, DU_jump_begin + NDIM * n_qpoints, 0.0);
                }

                // Interpolate x, du, and dv at the quadrature points via
                // accumulation, e.g., x(qp) = sum_k x_k * phi_k(qp) for each
                // qp.
                //
                // TODO: Switch to using interpolate() helper functions instead of manually accumulating.
                for (unsigned int k = 0; k < n_nodes; ++k)
                {
                    for (unsigned int qp = 0; qp < n_qpoints; ++qp)
                    {
                        const double& p = phi_X[k][qp];
                        const double& p2 = phi_DU_jump[k][qp];
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            x_qp[NDIM * (qp_offset + qp) + d] += x_node[k][d] * p;
                        }
                        if (d_use_velocity_jump_conditions)
                        {
                            for (unsigned int axis = 0; axis < NDIM; ++axis)
                            {
                                for (unsigned int d = 0; d < NDIM; ++d)
                                {
                                    DU_jump_qp[axis][NDIM * (qp_offset + qp) + d] += DU_jump_node[axis][k][d] * p2;
                                }
                            }
                        }
                    }
                }
                for (unsigned int qp = 0; qp < n_qpoints; ++qp)
                {
                    for (unsigned int l = 0; l < NDIM - 1; ++l)
                    {
                        interpolate(dx_dxi[l], qp, x_node, *dphi_dxi_X[l]);
                    }
                    if (NDIM == 2)
                    {
                        dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                    }
                    n = (dx_dxi[0].cross(dx_dxi[1])).unit();
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        n_qp[NDIM * (qp_offset + qp) + d] = n(d);
                        x_io_qp[NDIM * (qp_offset + qp) + d] =
                            x_qp[NDIM * (qp_offset + qp) + d] + (d_traction_interior_side ? -1.0 : 1.0) * n(d) * dh;
                    }
                }
                qp_offset += n_qpoints;
            }

            // Interpolate values from the Cartesian grid patch to the
            // quadrature points.
            //
            // NOTE: Values are interpolated only to those quadrature points
            // that are within the patch interior.
            const Box<NDIM>& interp_box = patch->getBox();
            Pointer<PatchData<NDIM> > u_data = patch->getPatchData(u_data_idx);

            const Box<NDIM> ghost_box = Box<NDIM>::grow(patch->getBox(), IntVector<NDIM>(u_ghost_num));

            Pointer<CellData<NDIM, double> > u_cc_data = u_data;
            if (u_cc_data)
            {
                LEInteractor::interpolate(
                    U_qp, NDIM, x_qp, NDIM, u_cc_data, patch, interp_box, d_default_interp_spec.kernel_fcn);
            }
            Pointer<SideData<NDIM, double> > u_sc_data = u_data;
            if (u_sc_data && !d_use_velocity_jump_conditions)
            {
                LEInteractor::interpolate(
                    U_qp, NDIM, x_qp, NDIM, u_sc_data, patch, interp_box, d_default_interp_spec.kernel_fcn);
            }
            else if (u_sc_data && d_use_velocity_jump_conditions)
            {
                LEInteractor::interpolate(
                    U_io_qp, NDIM, x_io_qp, NDIM, u_sc_data, patch, ghost_box, d_default_interp_spec.kernel_fcn);

                const IntVector<NDIM>& u_gcw = u_sc_data->getGhostCellWidth();
                const int u_depth = u_sc_data->getDepth();
                TBOX_ASSERT(u_depth == 1);

                // Keep the quadrature points that are inside the interpolation box.
                std::vector<int> local_indices;
                for (unsigned int k = 0; k < n_qpoints_patch; ++k)
                {
                    const double* const x = &x_qp[NDIM * k];
                    const Index<NDIM> i = IndexUtilities::getCellIndex(x, patch_geom, patch_box);
                    if (interp_box.contains(i)) local_indices.push_back(k);

                    const double* const x_io = &x_io_qp[NDIM * k];
                    const Index<NDIM> io = IndexUtilities::getCellIndex(
                        x_io, x_lower_ghost, x_upper_ghost, patch_geom->getDx(), ghost_box.lower(), ghost_box.upper());

                    // Some kind of assertation can be applied here using the indices of the cells away from the
                    // interfce
                }
                if (local_indices.empty()) continue;
                Index<NDIM> ic_lower, ic_upper, ic_center;
                std::array<std::array<double, 2>, NDIM> w, wr;
                std::vector<double> U_axis(n_qpoints_patch, 0.0);
                std::vector<double> U_axis_o(n_qpoints_patch, 0.0);
                Box<NDIM> side_boxes[NDIM];

                for (int axis = 0; axis < NDIM; ++axis)
                {
                    side_boxes[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                }
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    IBTK::Point x_lower_axis, x_upper_axis;

                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        x_lower_axis[d] = patch_x_lower[d];
                        x_upper_axis[d] = patch_x_upper[d];
                    }
                    x_lower_axis[axis] -= 0.5 * patch_dx[axis];
                    x_upper_axis[axis] += 0.5 * patch_dx[axis];

                    const Index<NDIM>& ilower = side_boxes[axis].lower();
                    const Index<NDIM>& iupper = side_boxes[axis].upper();

                    typedef boost::multi_array_types::extent_range range;
                    boost::const_multi_array_ref<double, NDIM> u_sc_data_array(
                        u_sc_data->getPointer(axis),
                        (boost::extents[range(ilower[0] - u_gcw[0], iupper[0] + u_gcw[0] + 1)]
                                       [range(ilower[1] - u_gcw[1], iupper[1] + u_gcw[1] + 1)]
#if (NDIM == 3)
                                       [range(ilower[2] - u_gcw[2], iupper[2] + u_gcw[2] + 1)]
#endif
                         ),
                        boost::fortran_storage_order());

                    for (unsigned int k = 0; k < local_indices.size(); ++k)
                    {
                        const int s = local_indices[k];
                        IBTK::Point x, x_cell, xo, x_cell_o;
                        const double* const dx = patch_dx;
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            x[d] = x_qp[s * NDIM + d];
                            ic_center[d] = ilower[d] + boost::math::iround((x[d] - x_lower_axis[d]) / dx[d] - 0.5);
                            x_cell[d] = x_lower_axis[d] + ((ic_center[d] - ilower[d]) + 0.5) * dx[d];
                            if (x[d] <= x_cell[d])
                            {
                                ic_lower[d] = ic_center[d] - 1;
                                ic_upper[d] = ic_center[d];
                            }
                            else
                            {
                                ic_lower[d] = ic_center[d];
                                ic_upper[d] = ic_center[d] + 1;
                            }

                            if (x[d] <= x_cell[d])
                            {
                                w[d][0] = (x_cell[d] - x[d]) / dx[d];
                            }
                            else
                            {
                                w[d][0] = 1.0 + (x_cell[d] - x[d]) / dx[d];
                            }
                            w[d][1] = 1.0 - w[d][0];
                            wr[d][0] = +w[d][0];
                            wr[d][1] = -w[d][1];
                        }

                        boost::multi_array<double, NDIM + 1> Ujump(
                            boost::extents[range(ic_lower[0], ic_upper[0] + 1)][range(ic_lower[1], ic_upper[1] + 1)]
#if (NDIM == 3)
                                          [range(ic_lower[2], ic_upper[2] + 1)]
#endif
                                          [range(0, NDIM)]);

                        boost::multi_array<double, NDIM + 1> interpCoeff(
                            boost::extents[range(ic_lower[0], ic_upper[0] + 1)][range(ic_lower[1], ic_upper[1] + 1)]
#if (NDIM == 3)
                                          [range(ic_lower[2], ic_upper[2] + 1)]
#endif
                                          [range(0, NDIM)]);

                        VectorValue<double> norm_vec, du_jump, coeff_vec, wrc;
                        // Loop over indices to calculate the interp coefficients (Lower=0, Upper=1)

                        for (int d = 0; d < NDIM; ++d) norm_vec(d) = n_qp[s * NDIM + d];

                        Box<NDIM> stencil_box(ic_lower, ic_upper);

                        for (int d = 0; d < NDIM; ++d)
                        {
                            for (BoxIterator<NDIM> b(stencil_box); b; b++)
                            {
                                const Index<NDIM>& ic = b();
                                for (int j = 0; j < NDIM; ++j) wrc(j) = wr[j][ic_upper[j] - ic[j]];
#if (NDIM == 2)
                                interpCoeff[ic[0]][ic[1]][d] = (norm_vec * wrc) * norm_vec(d);
#endif
#if (NDIM == 3)
                                interpCoeff[ic[0]][ic[1]][ic[2]][d] = (norm_vec * wrc) * norm_vec(d);
#endif
                            }
                        }

                        for (int d = 0; d < NDIM; ++d)
                        {
                            for (BoxIterator<NDIM> b(stencil_box); b; b++)
                            {
                                const Index<NDIM>& ic = b();
                                for (int j = 0; j < NDIM; ++j) du_jump(j) = DU_jump_qp[d][s * NDIM + j];
#if (NDIM == 2)
                                coeff_vec =
                                    VectorValue<double>(interpCoeff[ic[0]][ic[1]][0], interpCoeff[ic[0]][ic[1]][1]);
                                Ujump[ic[0]][ic[1]][d] = dx[0] * w[0][ic[0] - ic_lower[0]] * w[1][ic[1] - ic_lower[1]] *
                                                         (coeff_vec * du_jump);
#endif

#if (NDIM == 3)
                                coeff_vec = VectorValue<double>(interpCoeff[ic[0]][ic[1]][ic[2]][0],
                                                                interpCoeff[ic[0]][ic[1]][ic[2]][1],
                                                                interpCoeff[ic[0]][ic[1]][ic[2]][2]);
                                Ujump[ic[0]][ic[1]][ic[2]][d] = dx[0] * w[0][ic[0] - ic_lower[0]] *
                                                                w[1][ic[1] - ic_lower[1]] * w[2][ic[2] - ic_lower[2]] *
                                                                (norm_vec * du_jump) * norm_vec(d);
#endif
                            }
                        }
                        // Accumulate the value of U at the current location.
                        U_axis[s] = 0.0;

                        for (BoxIterator<NDIM> b(stencil_box); b; b++)
                        {
                            const Index<NDIM>& ic = b();
#if (NDIM == 2)

                            U_axis[s] +=
                                w[0][ic[0] - ic_lower[0]] * w[1][ic[1] - ic_lower[1]] * u_sc_data_array[ic[0]][ic[1]];
                            const double nproj = n_qp[s * NDIM + 0] * wr[0][ic_upper[0] - ic[0]] +
                                                 n_qp[s * NDIM + 1] * wr[1][ic_upper[1] - ic[1]];

                            const double CC = (nproj > 0.0) ? Ujump[ic[0]][ic[1]][axis] : 0.0;
                            U_axis[s] -= CC / mu;
#endif
#if (NDIM == 3)

                            U_axis[s] += w[0][ic[0] - ic_lower[0]] * w[1][ic[1] - ic_lower[1]] *
                                         w[2][ic[2] - ic_lower[2]] * u_sc_data_array[ic[0]][ic[1]][ic[2]];
                            const double nproj = n_qp[s * NDIM + 0] * wr[0][ic_upper[0] - ic[0]] +
                                                 n_qp[s * NDIM + 1] * wr[1][ic_upper[1] - ic[1]] +
                                                 n_qp[s * NDIM + 2] * wr[2][ic_upper[2] - ic[2]];

                            const double CC = (nproj > 0.0) ? Ujump[ic[0]][ic[1]][ic[2]][axis] : 0.0;
                            U_axis[s] -= CC / mu;
#endif
                        }
                    }

                    for (unsigned int k = 0; k < local_indices.size(); ++k)
                    {
                        U_qp[NDIM * local_indices[k] + axis] = U_axis[local_indices[k]];

                        if (dh != 0.0)
                        {
                            WSS_qp[NDIM * local_indices[k] + axis] =
                                mu * (d_traction_interior_side ? -1.0 : 1.0) * (1.0 / dh) *
                                (U_io_qp[NDIM * local_indices[k] + axis] - U_qp[NDIM * local_indices[k] + axis]);
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name << ": The width for the wall shear stress hasn't been setup!"
                                                     << std::endl);
                        }
                    }
                }
            }
            // Loop over the elements and accumulate the right-hand-side values.
            qrule.reset();
            qp_offset = 0;
            for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems[e_idx];
                const auto& U_dof_indices = U_dof_map_cache.dof_indices(elem);
                const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    U_rhs_e[d].resize(static_cast<int>(U_dof_indices[d].size()));
                    U_n_rhs_e[d].resize(static_cast<int>(U_dof_indices[d].size()));
                    U_t_rhs_e[d].resize(static_cast<int>(U_dof_indices[d].size()));
                    if (d_use_velocity_jump_conditions)
                    {
                        const auto& WSS_dof_indices = WSS_dof_map_cache->dof_indices(elem);
                        WSS_rhs_e[d].resize(static_cast<int>(WSS_dof_indices.size()));
                    }
                }
                get_values_for_interpolation(x_node, *X_ghost_vec, X_dof_indices);
                const bool qrule_changed =
                    FEDataManager::updateInterpQuadratureRule(qrule, d_default_interp_spec, elem, x_node, patch_dx_min);
                if (qrule_changed) fe_X->attach_quadrature_rule(qrule.get());

                fe_X->reinit(elem);

                if (qrule_changed)
                {
                    fe_DU_jump->attach_quadrature_rule(qrule.get());
                }
                fe_DU_jump->reinit(elem);

                const unsigned int n_qpoints = qrule->n_points();
                const size_t n_basis = U_rhs_e[0].size();
                const size_t n_basis_DU_jump = WSS_rhs_e[0].size();
                for (unsigned int qp = 0; qp < n_qpoints; ++qp)
                {
                    for (unsigned int k = 0; k < NDIM - 1; ++k)
                    {
                        interpolate(dx_dxi[k], qp, x_node, *dphi_dxi_X[k]);
                    }
                    if (NDIM == 2)
                    {
                        dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                    }
                    n = (dx_dxi[0].cross(dx_dxi[1])).unit();
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        U(d) = U_qp[NDIM * (qp_offset + qp) + d];
                        if (d_use_velocity_jump_conditions)
                        {
                            WSS(d) = WSS_qp[NDIM * (qp_offset + qp) + d];
                        }
                    }
                    U_n = (U * n) * n;
                    U_t = U - U_n;
                    for (unsigned int k = 0; k < n_basis; ++k)
                    {
                        const double p_JxW = phi_X[k][qp] * JxW[qp];
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            U_rhs_e[d](k) += U(d) * p_JxW;
                        }
                    }

                    for (unsigned int k = 0; k < n_basis_DU_jump; ++k)
                    {
                        const double p_JxW = phi_DU_jump[k][qp] * JxW[qp];
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            U_n_rhs_e[d](k) += U_n(d) * p_JxW;
                            U_t_rhs_e[d](k) += U_t(d) * p_JxW;
                            if (d_use_velocity_jump_conditions)
                            {
                                WSS_rhs_e[d](k) += WSS(d) * p_JxW;
                            }
                        }
                    }
                }
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    dof_id_scratch = U_dof_indices[d];
                    U_dof_map.constrain_element_vector(U_rhs_e[d], dof_id_scratch);
                    U_rhs_vec->add_vector(U_rhs_e[d], dof_id_scratch);

                    dof_id_scratch = U_dof_indices[d];
                    U_dof_map.constrain_element_vector(U_n_rhs_e[d], dof_id_scratch);
                    U_n_rhs_vec->add_vector(U_n_rhs_e[d], dof_id_scratch);

                    dof_id_scratch = U_dof_indices[d];
                    U_dof_map.constrain_element_vector(U_t_rhs_e[d], dof_id_scratch);
                    U_t_rhs_vec->add_vector(U_t_rhs_e[d], dof_id_scratch);

                    if (d_use_velocity_jump_conditions)
                    {
                        const auto& WSS_dof_indices = WSS_dof_map_cache->dof_indices(elem);
                        dof_id_scratch = WSS_dof_indices[d];
                        WSS_dof_map->constrain_element_vector(WSS_rhs_e[d], dof_id_scratch);
                        WSS_rhs_vec->add_vector(WSS_rhs_e[d], dof_id_scratch);
                    }
                }
                qp_offset += n_qpoints;
            }
        }
        U_rhs_vec->close();
        U_n_rhs_vec->close();
        U_t_rhs_vec->close();
        if (d_use_velocity_jump_conditions)
        {
            WSS_rhs_vec->close();
            d_fe_data_managers[part]->computeL2Projection(
                *WSS_vec, *WSS_rhs_vec, WSS_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
        }
        // Solve for the nodal values.
        d_fe_data_managers[part]->computeL2Projection(
            *U_vec, *U_rhs_vec, VELOCITY_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
        d_fe_data_managers[part]->computeL2Projection(
            *U_n_vec, *U_n_rhs_vec, VELOCITY_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
        d_fe_data_managers[part]->computeL2Projection(
            *U_t_vec, *U_t_rhs_vec, VELOCITY_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
    }
    return;
} // interpolateVelocity

void
IBFESurfaceMethod::forwardEulerStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        if (d_use_direct_forcing)
        {
            ierr = VecCopy(d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            ierr =
                VecWAXPY(d_X_new_vecs[part]->vec(), dt, d_U_current_vecs[part]->vec(), d_X_current_vecs[part]->vec());
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecAXPBYPCZ(
            d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_X_new_vecs[part]->close();
        d_X_half_vecs[part]->close();
    }
    return;
} // eulerStep

void
IBFESurfaceMethod::midpointStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        if (d_use_direct_forcing)
        {
            ierr = VecCopy(d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            ierr = VecWAXPY(d_X_new_vecs[part]->vec(), dt, d_U_half_vecs[part]->vec(), d_X_current_vecs[part]->vec());
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecAXPBYPCZ(
            d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_X_new_vecs[part]->close();
        d_X_half_vecs[part]->close();
    }
    return;
} // midpointStep

void
IBFESurfaceMethod::trapezoidalStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        if (d_use_direct_forcing)
        {
            ierr = VecCopy(d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            ierr = VecWAXPY(
                d_X_new_vecs[part]->vec(), 0.5 * dt, d_U_current_vecs[part]->vec(), d_X_current_vecs[part]->vec());
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPY(d_X_new_vecs[part]->vec(), 0.5 * dt, d_U_new_vecs[part]->vec());
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecAXPBYPCZ(
            d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_X_new_vecs[part]->close();
        d_X_half_vecs[part]->close();
    }
    return;
} // trapezoidalStep

void
IBFESurfaceMethod::computeLagrangianForce(const double data_time)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
        const MeshBase& mesh = equation_systems->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        // Setup global and elemental right-hand-side vectors.
        NumericVector<double>* F_vec = d_F_half_vecs[part];
        std::unique_ptr<NumericVector<double> > F_rhs_vec = F_vec->zero_clone();
        std::array<DenseVector<double>, NDIM> F_rhs_e;
        VectorValue<double>& F_integral = d_lag_surface_force_integral[part];
        F_integral.zero();

        NumericVector<double>* X_vec = d_X_half_vecs[part];
        double surface_area = 0.0;

        NumericVector<double>* P_jump_vec = d_P_jump_half_vecs[part];
        std::unique_ptr<NumericVector<double> > P_jump_rhs_vec;
        DenseVector<double> P_jump_rhs_e;
        if (d_use_pressure_jump_conditions)
        {
            P_jump_rhs_vec = P_jump_vec->zero_clone();
        }
        double P_jump_rhs_integral = 0.0;

        std::array<NumericVector<double>*, NDIM> DU_jump_vec;
        std::array<std::unique_ptr<NumericVector<double> >, NDIM> DU_jump_rhs_vec;
        std::array<std::vector<DenseVector<double> >, NDIM> DU_jump_rhs_e;
        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                DU_jump_vec[d] = d_DU_jump_half_vecs[part][d];
                DU_jump_rhs_vec[d] = DU_jump_vec[d]->zero_clone();
            }
        }

        // Extract the FE systems and DOF maps, and setup the FE objects.
        System& F_system = equation_systems->get_system(FORCE_SYSTEM_NAME);
        const DofMap& F_dof_map = F_system.get_dof_map();
        FEDataManager::SystemDofMapCache& F_dof_map_cache =
            *d_fe_data_managers[part]->getDofMapCache(FORCE_SYSTEM_NAME);
        FEType F_fe_type = F_dof_map.variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(F_dof_map.variable_type(d) == F_fe_type);
        }

        System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
        const DofMap& X_dof_map = X_system.get_dof_map();
        FEDataManager::SystemDofMapCache& X_dof_map_cache =
            *d_fe_data_managers[part]->getDofMapCache(COORDS_SYSTEM_NAME);
        FEType X_fe_type = X_dof_map.variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
        }
        TBOX_ASSERT(X_fe_type == F_fe_type);
        NumericVector<double>& X0_vec = X_system.get_vector("INITIAL_COORDINATES");

        System* P_jump_system;
        const DofMap* P_jump_dof_map;
        FEDataManager::SystemDofMapCache* P_jump_dof_map_cache;
        FEType P_jump_fe_type = INVALID_FE;
        if (d_use_pressure_jump_conditions)
        {
            P_jump_system = &equation_systems->get_system(PRESSURE_JUMP_SYSTEM_NAME);
            P_jump_dof_map = &P_jump_system->get_dof_map();
            P_jump_dof_map_cache = d_fe_data_managers[part]->getDofMapCache(PRESSURE_JUMP_SYSTEM_NAME);
            P_jump_fe_type = P_jump_dof_map->variable_type(0);
        }

        std::array<System*, NDIM> DU_jump_system;
        std::array<DofMap*, NDIM> DU_jump_dof_map;
        std::array<FEDataManager::SystemDofMapCache*, NDIM> DU_jump_dof_map_cache;
        FEType DU_jump_fe_type = INVALID_FE;
        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                DU_jump_system[i] = &equation_systems->get_system(VELOCITY_JUMP_SYSTEM_NAME[i]);
                DU_jump_dof_map[i] = &DU_jump_system[i]->get_dof_map();
                DU_jump_dof_map_cache[i] = d_fe_data_managers[part]->getDofMapCache(VELOCITY_JUMP_SYSTEM_NAME[i]);
            }
            DU_jump_fe_type = DU_jump_dof_map[0]->variable_type(0);
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                for (unsigned int j = 0; j < NDIM; ++j)
                {
                    TBOX_ASSERT(DU_jump_dof_map[i]->variable_type(j) == DU_jump_fe_type);
                }
            }
        }
        TBOX_ASSERT(P_jump_fe_type == DU_jump_fe_type);

        std::unique_ptr<QBase> qrule = QBase::build(d_default_quad_type[part], dim, d_default_quad_order[part]);

        std::unique_ptr<FEBase> fe_X = FEBase::build(dim, X_fe_type);
        fe_X->attach_quadrature_rule(qrule.get());
        const std::vector<double>& JxW = fe_X->get_JxW();
        const std::vector<std::vector<double> >& phi_X = fe_X->get_phi();
        std::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi_X;
        dphi_dxi_X[0] = &fe_X->get_dphidxi();
        if (NDIM > 2) dphi_dxi_X[1] = &fe_X->get_dphideta();

        std::unique_ptr<FEBase> fe_P_jump = FEBase::build(dim, P_jump_fe_type);
        fe_P_jump->attach_quadrature_rule(qrule.get());
        const std::vector<std::vector<double> >& phi_P_jump = fe_P_jump->get_phi();

        FEDataInterpolation fe_interpolator(dim, d_fe_data_managers[part]);
        fe_interpolator.attachQuadratureRule(qrule.get());

        std::vector<size_t> surface_force_fcn_system_idxs;
        fe_interpolator.setupInterpolatedSystemDataIndexes(
            surface_force_fcn_system_idxs, d_lag_surface_force_fcn_data[part].system_data, equation_systems);
        std::vector<size_t> surface_pressure_fcn_system_idxs;
        fe_interpolator.setupInterpolatedSystemDataIndexes(
            surface_pressure_fcn_system_idxs, d_lag_surface_pressure_fcn_data[part].system_data, equation_systems);
        fe_interpolator.init(/*use_IB_ghosted_vecs*/ false);

        std::vector<const std::vector<double>*> surface_force_var_data, surface_pressure_var_data;
        std::vector<const std::vector<VectorValue<double> >*> surface_force_grad_var_data,
            surface_pressure_grad_var_data;

        // Loop over the elements to compute the right-hand side vector.
        boost::multi_array<double, 2> X_node, x_node;
        double DU[NDIM][NDIM];
        TensorValue<double> FF;
        VectorValue<double> F, F_b, F_s, F_qp, N, X, n, x;
        std::array<VectorValue<double>, 2> dX_dxi, dx_dxi;
        std::vector<libMesh::dof_id_type> dof_id_scratch;
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            const auto& F_dof_indices = F_dof_map_cache.dof_indices(elem);
            const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                F_rhs_e[d].resize(static_cast<int>(F_dof_indices[d].size()));
            }
            if (d_use_pressure_jump_conditions)
            {
                const auto& P_jump_dof_indices = P_jump_dof_map_cache->dof_indices(elem);
                P_jump_rhs_e.resize(static_cast<int>(P_jump_dof_indices[0].size()));
            }
            if (d_use_velocity_jump_conditions)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    const auto& DU_jump_dof_indices = DU_jump_dof_map_cache[d]->dof_indices(elem);
                    for (unsigned int k = 0; k < NDIM; ++k)
                    {
                        DU_jump_rhs_e[d][k].resize(static_cast<int>(DU_jump_dof_indices[k].size()));
                    }
                }
            }

            fe_X->reinit(elem);

            fe_interpolator.reinit(elem);
            fe_interpolator.collectDataForInterpolation(elem);
            fe_interpolator.interpolate(elem);

            fe_P_jump->reinit(elem);

            get_values_for_interpolation(x_node, *X_vec, X_dof_indices);
            get_values_for_interpolation(X_node, X0_vec, X_dof_indices);
            const unsigned int n_qpoints = qrule->n_points();
            const size_t n_basis = phi_X.size();
            const size_t n_basis2 = phi_P_jump.size();
            for (unsigned int qp = 0; qp < n_qpoints; ++qp)
            {
                interpolate(X, qp, X_node, phi_X);
                interpolate(x, qp, x_node, phi_X);
                for (unsigned int k = 0; k < NDIM - 1; ++k)
                {
                    interpolate(dX_dxi[k], qp, X_node, *dphi_dxi_X[k]);
                    interpolate(dx_dxi[k], qp, x_node, *dphi_dxi_X[k]);
                }
                if (NDIM == 2)
                {
                    dX_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                    dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                }

                // Construct unit vectors in the reference and current
                // configurations.
                N = dX_dxi[0].cross(dX_dxi[1]);
                const double dA = N.norm();
                N = N.unit();
                n = dx_dxi[0].cross(dx_dxi[1]);
                const double da = n.norm();
                n = n.unit();

                F.zero();

                if (d_lag_surface_pressure_fcn_data[part].fcn)
                {
                    // Compute the value of the pressure at the quadrature point
                    // and add the corresponding force to the right-hand-side
                    // vector.
                    double P = 0;
                    fe_interpolator.setInterpolatedDataPointers(surface_pressure_var_data,
                                                                surface_pressure_grad_var_data,
                                                                surface_pressure_fcn_system_idxs,
                                                                elem,
                                                                qp);
                    d_lag_surface_pressure_fcn_data[part].fcn(P,
                                                              n,
                                                              N,
                                                              FF,
                                                              x,
                                                              X,
                                                              elem,
                                                              /*side*/ 0,
                                                              surface_pressure_var_data,
                                                              surface_pressure_grad_var_data,
                                                              data_time,
                                                              d_lag_surface_pressure_fcn_data[part].ctx);
                    F -= P * n * da / dA;
                }

                if (d_lag_surface_force_fcn_data[part].fcn)
                {
                    // Compute the value of the surface force at the quadrature
                    // point and add the corresponding force to the
                    // right-hand-side vector.
                    fe_interpolator.setInterpolatedDataPointers(
                        surface_force_var_data, surface_force_grad_var_data, surface_force_fcn_system_idxs, elem, qp);
                    d_lag_surface_force_fcn_data[part].fcn(F_s,
                                                           n,
                                                           N,
                                                           FF,
                                                           x,
                                                           X,
                                                           elem,
                                                           /*side*/ 0,
                                                           surface_force_var_data,
                                                           surface_force_grad_var_data,
                                                           data_time,
                                                           d_lag_surface_force_fcn_data[part].ctx);
                    F += F_s;
                }

                const double P_j = F * n * dA / da;
                for (unsigned int i = 0; i < NDIM; ++i)
                    for (unsigned int k = 0; k < NDIM; ++k)
                        DU[i][k] = -(dA / da) * (F(i) - F * n * n(i)) * n(k); // [Ux] , [Uy], [Uz]

                for (unsigned int d = 0; d < NDIM; ++d) F_integral(d) += F(d) * JxW[qp];

                // Remote the part of the force that has been already included in the jump

                if (d_use_pressure_jump_conditions && !d_use_velocity_jump_conditions)
                {
                    F -= (F * n) * n;
                }
                if (!d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
                {
                    F = (F * n) * n;
                }
                if (d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
                {
                    F = 0.0;
                }

                // Add the boundary forces to the right-hand-side vector.
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    F_qp = F * phi_X[k][qp] * JxW[qp];
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        F_rhs_e[i](k) += F_qp(i);
                    }
                }
                for (unsigned int k = 0; k < n_basis2; ++k)
                {
                    if (d_use_pressure_jump_conditions)
                    {
                        P_jump_rhs_e(k) += P_j * phi_P_jump[k][qp] * JxW[qp];
                    }
                    if (d_use_velocity_jump_conditions)
                    {
                        for (unsigned int i = 0; i < NDIM; ++i)
                        {
                            for (unsigned int j = 0; j < NDIM; ++j)
                            {
                                DU_jump_rhs_e[i][j](k) += DU[i][j] * phi_P_jump[k][qp] * JxW[qp];
                            }
                        }
                    }
                }
                if (d_use_pressure_jump_conditions)
                {
                    P_jump_rhs_integral += P_j * JxW[qp];
                    surface_area += JxW[qp];
                }
            }

            // Apply constraints (e.g., enforce periodic boundary conditions)
            // and add the elemental contributions to the global vector.
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                dof_id_scratch = F_dof_indices[i];
                F_dof_map.constrain_element_vector(F_rhs_e[i], dof_id_scratch);
                F_rhs_vec->add_vector(F_rhs_e[i], dof_id_scratch);
                if (d_use_velocity_jump_conditions)
                {
                    const auto& DU_jump_dof_indices = DU_jump_dof_map_cache[i]->dof_indices(elem);
                    for (unsigned int k = 0; k < NDIM; ++k)
                    {
                        dof_id_scratch = DU_jump_dof_indices[k];
                        DU_jump_dof_map[i]->constrain_element_vector(DU_jump_rhs_e[i][k], dof_id_scratch);
                        DU_jump_rhs_vec[i]->add_vector(DU_jump_rhs_e[i][k], dof_id_scratch);
                    }
                }
            }
            if (d_use_pressure_jump_conditions)
            {
                const auto& P_jump_dof_indices = P_jump_dof_map_cache->dof_indices(elem);
                dof_id_scratch = P_jump_dof_indices[0];
                P_jump_dof_map->constrain_element_vector(P_jump_rhs_e, dof_id_scratch);
                P_jump_rhs_vec->add_vector(P_jump_rhs_e, dof_id_scratch);
            }
        }

        SAMRAI_MPI::sumReduction(&F_integral(0), NDIM);

        // Solve for F.
        d_fe_data_managers[part]->computeL2Projection(
            *F_vec, *F_rhs_vec, FORCE_SYSTEM_NAME, d_use_consistent_mass_matrix);
        if (d_use_pressure_jump_conditions)
        {
            d_fe_data_managers[part]->computeL2Projection(
                *P_jump_vec, *P_jump_rhs_vec, PRESSURE_JUMP_SYSTEM_NAME, d_use_consistent_mass_matrix);
            P_jump_rhs_integral = SAMRAI_MPI::sumReduction(P_jump_rhs_integral);
            surface_area = SAMRAI_MPI::sumReduction(surface_area);
            if (d_normalize_pressure_jump) P_jump_vec->add(-P_jump_rhs_integral / surface_area);
            P_jump_vec->close();
        }
        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                d_fe_data_managers[part]->computeL2Projection(
                    *DU_jump_vec[d], *DU_jump_rhs_vec[d], VELOCITY_JUMP_SYSTEM_NAME[d], d_use_consistent_mass_matrix);
                DU_jump_vec[d]->close();
            }
        }
    }
    return;
} // computeLagrangianForce

void
IBFESurfaceMethod::spreadForce(const int f_data_idx,
                               RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                               const std::vector<Pointer<RefineSchedule<NDIM> > >& /*f_prolongation_scheds*/,
                               const double data_time)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        PetscVector<double>* X_vec = d_X_half_vecs[part];
        PetscVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
        PetscVector<double>* F_vec = d_F_half_vecs[part];
        PetscVector<double>* F_ghost_vec = d_F_IB_ghost_vecs[part];
        copy_and_synch(*X_vec, *X_ghost_vec);
        copy_and_synch(*F_vec, *F_ghost_vec);
        d_fe_data_managers[part]->spread(
            f_data_idx, *F_ghost_vec, *X_ghost_vec, FORCE_SYSTEM_NAME, f_phys_bdry_op, data_time);
        PetscVector<double>* P_jump_ghost_vec;
        std::array<PetscVector<double>*, NDIM> DU_jump_ghost_vec;
        std::array<PetscVector<double>*, NDIM> DU_jump_vec;
        if (d_use_pressure_jump_conditions)
        {
            PetscVector<double>* P_jump_vec = d_P_jump_half_vecs[part];
            P_jump_ghost_vec = d_P_jump_IB_ghost_vecs[part];
            copy_and_synch(*P_jump_vec, *P_jump_ghost_vec);
        }
        if (d_use_velocity_jump_conditions)
        {
            DU_jump_ghost_vec = d_DU_jump_IB_ghost_vecs[part];
            DU_jump_vec = d_DU_jump_half_vecs[part];
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                copy_and_synch(*DU_jump_vec[d], *DU_jump_vec[d]);
            }
        }

        if (d_use_pressure_jump_conditions || d_use_velocity_jump_conditions)
        {
            imposeJumpConditions(f_data_idx, *P_jump_ghost_vec, DU_jump_ghost_vec, *X_ghost_vec, data_time, part);
        }
    }
    return;
} // spreadForce

FEDataManager::InterpSpec
IBFESurfaceMethod::getDefaultInterpSpec() const
{
    return d_default_interp_spec;
}

FEDataManager::SpreadSpec
IBFESurfaceMethod::getDefaultSpreadSpec() const
{
    return d_default_spread_spec;
}

void
IBFESurfaceMethod::setInterpSpec(const FEDataManager::InterpSpec& interp_spec, const unsigned int part)
{
    TBOX_ASSERT(!d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_num_parts);
    d_interp_spec[part] = interp_spec;
    return;
}

void
IBFESurfaceMethod::setSpreadSpec(const FEDataManager::SpreadSpec& spread_spec, const unsigned int part)
{
    TBOX_ASSERT(!d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_num_parts);
    d_spread_spec[part] = spread_spec;
    return;
}

void
IBFESurfaceMethod::initializeFEEquationSystems()
{
    if (d_fe_equation_systems_initialized) return;

    const bool from_restart = RestartManager::getManager()->isFromRestart();

    // Create the FE data managers that manage mappings between the FE mesh
    // parts and the Cartesian grid.
    d_equation_systems.resize(d_num_parts, nullptr);
    d_fe_data_managers.resize(d_num_parts, nullptr);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Create FE data managers.
        const std::string manager_name = "IBFESurfaceMethod FEDataManager::" + std::to_string(part);
        d_fe_data_managers[part] = FEDataManager::getManager(manager_name, d_interp_spec[part], d_spread_spec[part]);
        d_ghosts = IntVector<NDIM>::max(d_ghosts, d_fe_data_managers[part]->getGhostCellWidth());

        // Create FE equation systems objects and corresponding variables.
        d_equation_systems[part] = new EquationSystems(*d_meshes[part]);
        EquationSystems* equation_systems = d_equation_systems[part];
        d_fe_data_managers[part]->setEquationSystems(equation_systems, d_max_level_number - 1);
        d_fe_data_managers[part]->COORDINATES_SYSTEM_NAME = COORDS_SYSTEM_NAME;
        if (from_restart)
        {
            const std::string& file_name = libmesh_restart_file_name(
                d_libmesh_restart_read_dir, d_libmesh_restart_restore_number, part, d_libmesh_restart_file_extension);
            const XdrMODE xdr_mode = (d_libmesh_restart_file_extension == "xdr" ? DECODE : READ);
            const int read_mode =
                EquationSystems::READ_HEADER | EquationSystems::READ_DATA | EquationSystems::READ_ADDITIONAL_DATA;
            equation_systems->read(file_name, xdr_mode, read_mode, /*partition_agnostic*/ true);
        }
        else
        {
            auto& X_system = equation_systems->add_system<System>(COORDS_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_system.add_variable("X_" + std::to_string(d), d_fe_order[part], d_fe_family[part]);
            }
            X_system.add_vector("INITIAL_COORDINATES", /*projections*/ true, GHOSTED);

            auto& dX_system = equation_systems->add_system<System>(COORD_MAPPING_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dX_system.add_variable("dX_" + std::to_string(d), d_fe_order[part], d_fe_family[part]);
            }

            auto& U_system = equation_systems->add_system<System>(VELOCITY_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_system.add_variable("U_" + std::to_string(d), d_fe_order[part], d_fe_family[part]);
            }

            auto& U_n_system = equation_systems->add_system<System>(NORMAL_VELOCITY_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_n_system.add_variable("U_n_" + std::to_string(d), d_fe_order[part], d_fe_family[part]);
            }

            auto& U_t_system = equation_systems->add_system<System>(TANGENTIAL_VELOCITY_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_t_system.add_variable("U_t_" + std::to_string(d), d_fe_order[part], d_fe_family[part]);
            }

            auto& F_system = equation_systems->add_system<System>(FORCE_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                F_system.add_variable("F_" + std::to_string(d), d_fe_order[part], d_fe_family[part]);
            }

            if (d_use_pressure_jump_conditions)
            {
                System& P_jump_system = equation_systems->add_system<System>(PRESSURE_JUMP_SYSTEM_NAME);
                System& P_system = equation_systems->add_system<System>(P_SYSTEM_NAME);
                if (d_use_l2_lagrange_family)
                {
                    P_jump_system.add_variable("P_jump_", d_fe_order[part], L2_LAGRANGE);
                    P_system.add_variable("P_", d_fe_order[part], L2_LAGRANGE);
                }
                else
                {
                    P_jump_system.add_variable("P_jump_", d_fe_order[part], d_fe_family[part]);
                    P_system.add_variable("P_", d_fe_order[part], d_fe_family[part]);
                }
            }

            if (d_use_velocity_jump_conditions)
            {
                std::array<System*, NDIM> DU_jump_system;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    DU_jump_system[d] = &equation_systems->add_system<System>(VELOCITY_JUMP_SYSTEM_NAME[d]);
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        const std::string system_name = "DU_jump_" + std::to_string(d) + "_" + std::to_string(i);
                        // TODO: The value of the FEFamily should just be read in from the input file directly instead
                        // of using a boolean.
                        if (d_use_l2_lagrange_family)
                        {
                            DU_jump_system[d]->add_variable(system_name, d_fe_order[part], L2_LAGRANGE);
                        }
                        else
                        {
                            DU_jump_system[d]->add_variable(system_name, d_fe_order[part], d_fe_family[part]);
                        }
                    }
                }

                System& WSS_system = equation_systems->add_system<System>(WSS_SYSTEM_NAME);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    const std::string system_name = "WSS_" + std::to_string(d);
                    // TODO: The value of the FEFamily should just be read in from the input file directly instead of
                    // using a boolean.
                    if (d_use_l2_lagrange_family)
                    {
                        WSS_system.add_variable(system_name, d_fe_order[part], L2_LAGRANGE);
                    }
                    else
                    {
                        WSS_system.add_variable(system_name, d_fe_order[part], d_fe_family[part]);
                    }
                }
            }

            if (d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
            {
                auto& TAU_system = equation_systems->add_system<System>(TAU_SYSTEM_NAME);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    std::string system_name = "TAU_" + std::to_string(d);
                    // TODO: The value of the FEFamily should just be read in from the input file directly instead of
                    // using a boolean.
                    if (d_use_l2_lagrange_family)
                    {
                        TAU_system.add_variable(system_name, d_fe_order[part], L2_LAGRANGE);
                    }
                    else
                    {
                        TAU_system.add_variable(system_name, d_fe_order[part], d_fe_family[part]);
                    }
                }
            }
        }
    }
    d_fe_equation_systems_initialized = true;
    return;
}

void
IBFESurfaceMethod::initializeFEData()
{
    if (d_fe_data_initialized) return;
    initializeFEEquationSystems();
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Initialize FE equation systems.
        EquationSystems* equation_systems = d_equation_systems[part];
        if (from_restart)
        {
            equation_systems->reinit();
        }
        else
        {
            equation_systems->init();
            initializeCoordinates(part);
            initializeVelocity(part);
        }
        updateCoordinateMapping(part);

        // Assemble systems.
        auto& X_system = equation_systems->get_system<System>(COORDS_SYSTEM_NAME);
        auto& dX_system = equation_systems->get_system<System>(COORD_MAPPING_SYSTEM_NAME);
        auto& U_system = equation_systems->get_system<System>(VELOCITY_SYSTEM_NAME);
        auto& U_n_system = equation_systems->get_system<System>(NORMAL_VELOCITY_SYSTEM_NAME);
        auto& U_t_system = equation_systems->get_system<System>(TANGENTIAL_VELOCITY_SYSTEM_NAME);
        auto& F_system = equation_systems->get_system<System>(FORCE_SYSTEM_NAME);

        X_system.assemble_before_solve = false;
        X_system.assemble();

        dX_system.assemble_before_solve = false;
        dX_system.assemble();

        U_system.assemble_before_solve = false;
        U_system.assemble();

        U_n_system.assemble_before_solve = false;
        U_n_system.assemble();

        U_t_system.assemble_before_solve = false;
        U_t_system.assemble();

        F_system.assemble_before_solve = false;
        F_system.assemble();

        if (d_use_pressure_jump_conditions)
        {
            System& P_jump_system = equation_systems->get_system<System>(PRESSURE_JUMP_SYSTEM_NAME);
            P_jump_system.assemble_before_solve = false;
            P_jump_system.assemble();

            System& P_system = equation_systems->get_system<System>(P_SYSTEM_NAME);
            P_system.assemble_before_solve = false;
            P_system.assemble();
        }
        if (d_use_velocity_jump_conditions)
        {
            std::array<System*, NDIM> DU_jump_system;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                DU_jump_system[d] = &equation_systems->get_system<System>(VELOCITY_JUMP_SYSTEM_NAME[d]);
                DU_jump_system[d]->assemble_before_solve = false;
                DU_jump_system[d]->assemble();
            }
            System& WSS_system = equation_systems->get_system<System>(WSS_SYSTEM_NAME);
            WSS_system.assemble_before_solve = false;
            WSS_system.assemble();
        }
        if (d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
        {
            System& TAU_system = equation_systems->get_system<System>(TAU_SYSTEM_NAME);
            TAU_system.assemble_before_solve = false;
            TAU_system.assemble();
        }
    }
    d_fe_data_initialized = true;
    return;
} // initializeFEData

void
IBFESurfaceMethod::registerEulerianVariables()
{
    mask_var = new CellVariable<NDIM, double>(d_object_name + "::mask");
    registerVariable(mask_current_idx,
                     mask_new_idx,
                     mask_scratch_idx,
                     mask_var,
                     d_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");
    return;
} // registerEulerianVariables

void
IBFESurfaceMethod::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                            Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
                                            int /*u_data_idx*/,
                                            const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
                                            const std::vector<Pointer<RefineSchedule<NDIM> > >& /*u_ghost_fill_scheds*/,
                                            int /*integrator_step*/,
                                            double /*init_data_time*/,
                                            bool /*initial_time*/)
{
    // Cache pointers to the patch hierarchy and gridding algorithm.
    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;

    // Initialize the FE data manager.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->reinitElementMappings();
    }

    d_is_initialized = true;
    return;
} // initializePatchHierarchy

void
IBFESurfaceMethod::registerLoadBalancer(Pointer<LoadBalancer<NDIM> > load_balancer, int workload_data_idx)
{
    IBAMR_DEPRECATED_MEMBER_FUNCTION1("IBFESurfaceMethod", "registerLoadBalancer");
    TBOX_ASSERT(load_balancer);
    d_load_balancer = load_balancer;
    d_workload_idx = workload_data_idx;

    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->registerLoadBalancer(load_balancer, workload_data_idx);
    }
    return;
} // registerLoadBalancer

void
IBFESurfaceMethod::addWorkloadEstimate(Pointer<PatchHierarchy<NDIM> > hierarchy, const int workload_data_idx)
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->addWorkloadEstimate(hierarchy, workload_data_idx);
    }
    return;
} // addWorkloadEstimate

void IBFESurfaceMethod::beginDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    // intentionally blank
    return;
} // beginDataRedistribution

void IBFESurfaceMethod::endDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                              Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    if (d_is_initialized)
    {
        for (unsigned int part = 0; part < d_num_parts; ++part)
        {
            d_fe_data_managers[part]->reinitElementMappings();
        }
    }
    return;
} // endDataRedistribution

void
IBFESurfaceMethod::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                       int level_number,
                                       double init_data_time,
                                       bool can_be_refined,
                                       bool initial_time,
                                       Pointer<BasePatchLevel<NDIM> > old_level,
                                       bool allocate_data)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->setPatchHierarchy(hierarchy);
        d_fe_data_managers[part]->setPatchLevels(0, finest_hier_level);
        d_fe_data_managers[part]->initializeLevelData(
            hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
    }
    return;
} // initializeLevelData

void
IBFESurfaceMethod::resetHierarchyConfiguration(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                               int coarsest_level,
                                               int /*finest_level*/)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->setPatchHierarchy(hierarchy);
        d_fe_data_managers[part]->setPatchLevels(0, hierarchy->getFinestLevelNumber());
        d_fe_data_managers[part]->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_hier_level);
    }
    return;
} // resetHierarchyConfiguration

void
IBFESurfaceMethod::applyGradientDetector(Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
                                         int level_number,
                                         double error_data_time,
                                         int tag_index,
                                         bool initial_time,
                                         bool uses_richardson_extrapolation_too)
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->applyGradientDetector(
            hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    }
    return;
} // applyGradientDetector

void
IBFESurfaceMethod::putToDatabase(Pointer<Database> db)
{
    db->putInteger("IBFE_METHOD_VERSION", IBFE_METHOD_VERSION);
    db->putInteger("d_num_parts", d_num_parts);
    db->putIntegerArray("d_ghosts", d_ghosts, NDIM);
    db->putBool("d_use_velocity_jump_conditions", d_use_velocity_jump_conditions);
    db->putBool("d_use_l2_lagrange_family", d_use_l2_lagrange_family);
    db->putBool("d_use_pressure_jump_conditions", d_use_pressure_jump_conditions);
    db->putBool("d_compute_fluid_traction", d_compute_fluid_traction);
    db->putBool("d_traction_interior_side", d_traction_interior_side);
    db->putDouble("d_traction_activation_time", d_traction_activation_time);
    db->putBool("d_use_consistent_mass_matrix", d_use_consistent_mass_matrix);
    db->putBool("d_use_direct_forcing", d_use_direct_forcing);
    return;
} // putToDatabase

void
IBFESurfaceMethod::writeFEDataToRestartFile(const std::string& restart_dump_dirname, unsigned int time_step_number)
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        const std::string& file_name =
            libmesh_restart_file_name(restart_dump_dirname, time_step_number, part, d_libmesh_restart_file_extension);
        const XdrMODE xdr_mode = (d_libmesh_restart_file_extension == "xdr" ? ENCODE : WRITE);
        const int write_mode = EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA;
        d_equation_systems[part]->write(file_name, xdr_mode, write_mode, /*partition_agnostic*/ true);
    }
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

namespace
{
struct IndexOrder
{
    inline bool operator()(const SAMRAI::hier::Index<NDIM>& lhs, const SAMRAI::hier::Index<NDIM>& rhs) const
    {
        return (lhs(0) < rhs(0)
#if (NDIM > 1)
                || (lhs(0) == rhs(0) && lhs(1) < rhs(1))
#if (NDIM > 2)
                || (lhs(0) == rhs(0) && lhs(1) == rhs(1) && lhs(2) < rhs(2))
#endif
#endif
        );
    }
};
} // namespace

void
IBFESurfaceMethod::interpolatePressureForTraction(const int p_data_idx, const double data_time, unsigned int part)
{
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_fe_data_managers[part]->getPatchHierarchy();

    NumericVector<double>* P_vec = d_P_half_vecs[part];

    NumericVector<double>* P_jump_ghost_vec = d_P_jump_IB_ghost_vecs[part];

    NumericVector<double>* X_vec = NULL;
    NumericVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];

    std::unique_ptr<NumericVector<double> > P_rhs_vec = (*P_vec).zero_clone();
    (*P_rhs_vec).zero();
    DenseVector<double> P_rhs_e;

    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        X_vec = d_X_current_vecs[part];
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        X_vec = d_X_half_vecs[part];
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        X_vec = d_X_new_vecs[part];
    }
    copy_and_synch(*X_vec, *X_ghost_vec);

    // Extract the FE systems and DOF maps, and setup the FE object.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    std::unique_ptr<QBase> qrule;

    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    DofMap& X_dof_map = X_system.get_dof_map();
    FEDataManager::SystemDofMapCache& X_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(COORDS_SYSTEM_NAME);
    FEType X_fe_type = X_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
    }

    System& P_system = equation_systems->get_system(P_SYSTEM_NAME);
    const DofMap& P_dof_map = P_system.get_dof_map();
    FEDataManager::SystemDofMapCache& P_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(P_SYSTEM_NAME);
    FEType P_fe_type = P_dof_map.variable_type(0);

    System& P_jump_system = equation_systems->get_system(PRESSURE_JUMP_SYSTEM_NAME);
    const DofMap& P_jump_dof_map = P_jump_system.get_dof_map();
    FEDataManager::SystemDofMapCache& P_jump_dof_map_cache =
        *d_fe_data_managers[part]->getDofMapCache(PRESSURE_JUMP_SYSTEM_NAME);
    FEType P_jump_fe_type = P_jump_dof_map.variable_type(0);
    TBOX_ASSERT(P_fe_type == P_jump_fe_type);

    std::unique_ptr<FEBase> fe_X = FEBase::build(dim, X_fe_type);
    const std::vector<double>& JxW = fe_X->get_JxW();
    const std::vector<std::vector<double> >& phi_X = fe_X->get_phi();
    std::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi_X;
    dphi_dxi_X[0] = &fe_X->get_dphidxi();
    if (NDIM > 2) dphi_dxi_X[1] = &fe_X->get_dphideta();

    std::unique_ptr<FEBase> fe_P = FEBase::build(dim, P_fe_type);
    const std::vector<std::vector<double> >& phi_P = fe_P->get_phi();

    const std::vector<std::vector<Elem*> >& active_patch_element_map =
        d_fe_data_managers[part]->getActivePatchElementMap();

    boost::multi_array<double, 2> x_node;
    boost::multi_array<double, 1> P_jump_node;
    std::vector<double> x_qp, x_io_qp;
    std::vector<double> P_qp, Q_io_qp, P_jump_qp, N_qp;
    std::array<VectorValue<double>, 2> dx_dxi;

    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_fe_data_managers[part]->getLevelNumber());
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    VectorValue<double> tau1, tau2, n;
    X_ghost_vec->close();
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const std::vector<Elem*>& patch_elems = active_patch_element_map[local_patch_num];
        const size_t num_active_patch_elems = patch_elems.size();
        if (!num_active_patch_elems) continue;
        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);

        const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const x_lower = pgeom->getXLower();
        const double* const x_upper = pgeom->getXUpper();

        const double* const dx = pgeom->getDx();

        double diag_dis = 0.0;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            diag_dis += dx[d] * dx[d];
        }
        const double dh = d_p_calc_width * sqrt(diag_dis);

        const int p_ghost_num = static_cast<int>(ceil(2.0 * dh / patch_dx_min));

        std::array<double, NDIM> x_lower_gh, x_upper_gh;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            x_lower_gh[d] = x_lower[d] - (static_cast<double>(p_ghost_num)) * dx[d];
            x_upper_gh[d] = x_upper[d] + (static_cast<double>(p_ghost_num)) * dx[d];
        }

        double* x_upper_ghost = x_upper_gh.data();
        double* x_lower_ghost = x_lower_gh.data();

        // Setup vectors to store the values of U and X at the quadrature
        // points.
        //
        // All this loop is doing is computing the total number of quadraturee
        // points associated with all of the elements we are currently
        // processing.  That number is n_qp_patch.
        unsigned int n_qp_patch = 0;
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
            get_values_for_interpolation(x_node, *X_ghost_vec, X_dof_indices);
            FEDataManager::updateInterpQuadratureRule(qrule, d_default_interp_spec, elem, x_node, patch_dx_min);
            n_qp_patch += qrule->n_points();
        }

        if (!n_qp_patch) continue;
        P_qp.resize(n_qp_patch);
        Q_io_qp.resize(n_qp_patch);
        P_jump_qp.resize(n_qp_patch);
        x_io_qp.resize(NDIM * n_qp_patch);
        x_qp.resize(NDIM * n_qp_patch);
        N_qp.resize(NDIM * n_qp_patch);
        std::fill(P_qp.begin(), P_qp.end(), 0.0);
        std::fill(Q_io_qp.begin(), Q_io_qp.end(), 0.0);
        std::fill(P_jump_qp.begin(), P_jump_qp.end(), 0.0);
        std::fill(N_qp.begin(), N_qp.end(), 0.0);

        // Loop over the elements and compute the positions of the quadrature points.
        qrule.reset();
        unsigned int qp_offset = 0;
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
            get_values_for_interpolation(x_node, *X_ghost_vec, X_dof_indices);
            const auto& P_jump_dof_indices = P_jump_dof_map_cache.dof_indices(elem);
            get_values_for_interpolation(P_jump_node, *P_jump_ghost_vec, P_jump_dof_indices[0]);
            const bool qrule_changed =
                FEDataManager::updateInterpQuadratureRule(qrule, d_default_interp_spec, elem, x_node, patch_dx_min);
            if (qrule_changed)
            {
                fe_X->attach_quadrature_rule(qrule.get());
                fe_P->attach_quadrature_rule(qrule.get());
            }
            fe_X->reinit(elem);
            fe_P->reinit(elem);

            const unsigned int n_node = elem->n_nodes();
            const unsigned int n_qp = qrule->n_points();

            // Zero out the values of X, du, and dv prior to accumulation.
            double* x_begin = &x_qp[NDIM * qp_offset];
            std::fill(x_begin, x_begin + NDIM * n_qp, 0.0);

            double* x_io_begin = &x_io_qp[NDIM * qp_offset];
            std::fill(x_io_begin, x_io_begin + NDIM * n_qp, 0.0);

            double* N_begin = &N_qp[NDIM * qp_offset];
            std::fill(N_begin, N_begin + NDIM * n_qp, 0.0);

            double* P_jump_begin = &P_jump_qp[qp_offset];
            std::fill(P_jump_begin, P_jump_begin + n_qp, 0.0);

            // Interpolate X, du, and dv at all of the quadrature points
            // via accumulation, i.e., X(qp) = sum_k X_k * phi_k(qp) for
            // each qp.

            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                for (unsigned int k = 0; k < NDIM - 1; ++k)
                {
                    interpolate(dx_dxi[k], qp, x_node, *dphi_dxi_X[k]);
                }
                if (NDIM == 2)
                {
                    dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                }
                n = (dx_dxi[0].cross(dx_dxi[1])).unit();

                for (unsigned int i = 0; i < NDIM; ++i)
                {
                    for (unsigned int k = 0; k < n_node; ++k)
                    {
                        const double& p_X = phi_X[k][qp];
                        x_qp[NDIM * (qp_offset + qp) + i] += x_node[k][i] * p_X;
                    }
                    N_qp[NDIM * (qp_offset + qp) + i] = n(i);

                    // Note that here we calculate the pressure on one side as the jump plus the pressure on the other
                    // side.
                    x_io_qp[NDIM * (qp_offset + qp) + i] =
                        x_qp[NDIM * (qp_offset + qp) + i] + (d_traction_interior_side ? -1.0 : 1.0) * n(i) * dh;
                }

                for (unsigned int k = 0; k < n_node; ++k)
                {
                    const double& p_P = phi_P[k][qp];
                    P_jump_qp[qp_offset + qp] += P_jump_node[k] * p_P;
                }
            }
            qp_offset += n_qp;
        }

        // Interpolate values from the Cartesian grid patch to the quadrature
        // points.
        //
        // NOTE: Values are interpolated only to those quadrature points that
        // are within the patch interior
        const Box<NDIM>& interp_box = patch->getBox();
        Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(p_data_idx);
        const Box<NDIM> ghost_box = Box<NDIM>::grow(patch->getBox(), IntVector<NDIM>(p_ghost_num));

        LEInteractor::interpolate(
            Q_io_qp, 1, x_io_qp, NDIM, p_data, patch, ghost_box, d_default_interp_spec.kernel_fcn);

        std::vector<int> local_indices;
        local_indices.clear();
        const int upper_bound = n_qp_patch;
        if (upper_bound == 0) continue;

        local_indices.reserve(upper_bound);
        for (unsigned int k = 0; k < n_qp_patch; ++k)
        {
            const double* const xx = &x_qp[NDIM * k];
            const Index<NDIM> i = IndexUtilities::getCellIndex(xx, patch_geom, interp_box);
            if (interp_box.contains(i)) local_indices.push_back(k);

            const double* const x_io = &x_io_qp[NDIM * k];
            const Index<NDIM> ip = IndexUtilities::getCellIndex(
                x_io, x_lower_ghost, x_upper_ghost, patch_geom->getDx(), ghost_box.lower(), ghost_box.upper());
            if (!ghost_box.contains(ip) && interp_box.contains(i))
            {
                TBOX_ERROR(d_object_name << "::IBFESurfaceMethod():\n"
                                         << " the pressure interpolation ghost width hasn't beeen properly set"
                                         << std::endl);
            }
        }

        const unsigned int nindices = static_cast<int>(local_indices.size());

        if (!local_indices.empty())
        {
            for (unsigned int k = 0; k < nindices; ++k)
            {
                P_qp[local_indices[k]] = Q_io_qp[local_indices[k]];
            }
        }

        // Loop over the elements and accumulate the right-hand-side values.
        qrule.reset();
        qp_offset = 0;
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];

            const auto& P_dof_indices = P_dof_map_cache.dof_indices(elem);
            P_rhs_e.resize(static_cast<int>(P_dof_indices[0].size()));

            const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
            get_values_for_interpolation(x_node, *X_ghost_vec, X_dof_indices);
            const bool qrule_changed =
                FEDataManager::updateInterpQuadratureRule(qrule, d_default_interp_spec, elem, x_node, patch_dx_min);
            if (qrule_changed)
            {
                fe_X->attach_quadrature_rule(qrule.get());
                fe_P->attach_quadrature_rule(qrule.get());
            }
            fe_X->reinit(elem);
            fe_P->reinit(elem);
            const unsigned int n_qp = qrule->n_points();
            const size_t n_basis2 = P_dof_indices.size();
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const int idx = qp_offset + qp;
                for (unsigned int k = 0; k < n_basis2; ++k)
                {
                    P_rhs_e(k) += P_qp[idx] * JxW[qp];
                }
            }

            auto dof_id_scratch = P_dof_indices[0];
            P_dof_map.constrain_element_vector(P_rhs_e, dof_id_scratch);
            P_rhs_vec->add_vector(P_rhs_e, dof_id_scratch);

            qp_offset += n_qp;
        }
    }
    P_rhs_vec->close();

    d_fe_data_managers[part]->computeL2Projection(
        *P_vec, *P_rhs_vec, P_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);

    d_X_half_vecs[part]->close();
    d_X_current_vecs[part]->close();
    d_X_new_vecs[part]->close();

    d_P_half_vecs[part]->close();
    d_P_jump_IB_ghost_vecs[part]->close();
    d_X_IB_ghost_vecs[part]->close();

    return;
} // interpolatePressureForTraction

void
IBFESurfaceMethod::computeFluidTraction(const double data_time, unsigned int part)
{
    NumericVector<double>* WSS_vec = NULL;
    NumericVector<double>* WSS_ghost_vec = d_WSS_IB_ghost_vecs[part];

    NumericVector<double>* P_vec = NULL;
    NumericVector<double>* P_ghost_vec = d_P_IB_ghost_vecs[part];

    NumericVector<double>* P_jump_vec = NULL;
    NumericVector<double>* P_jump_ghost_vec = d_P_jump_IB_ghost_vecs[part];

    std::array<PetscVector<double>*, NDIM> DU_jump_ghost_vec;
    std::array<PetscVector<double>*, NDIM> DU_jump_vec;

    NumericVector<double>* TAU_vec = d_TAU_half_vecs[part];

    NumericVector<double>* X_vec = NULL;

    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        X_vec = d_X_current_vecs[part];
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        X_vec = d_X_half_vecs[part];
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        X_vec = d_X_new_vecs[part];
    }
    NumericVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
    copy_and_synch(*X_vec, *X_ghost_vec);

    WSS_vec = d_WSS_half_vecs[part];
    copy_and_synch(*WSS_vec, *WSS_ghost_vec);

    P_vec = d_P_half_vecs[part];
    copy_and_synch(*P_vec, *P_ghost_vec);

    P_jump_vec = d_P_jump_half_vecs[part];
    copy_and_synch(*P_jump_vec, *P_jump_ghost_vec);

    DU_jump_ghost_vec = d_DU_jump_IB_ghost_vecs[part];
    DU_jump_vec = d_DU_jump_half_vecs[part];
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        copy_and_synch(*DU_jump_vec[d], *DU_jump_vec[d]);
    }
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        copy_and_synch(*DU_jump_vec[d], *DU_jump_vec[d]);
    }

    std::unique_ptr<NumericVector<double> > TAU_rhs_vec = TAU_vec->zero_clone();
    std::array<DenseVector<double>, NDIM> TAU_rhs_e;

    // Extract the FE systems and DOF maps, and setup the FE objects.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    std::unique_ptr<QBase> qrule;
    std::array<std::vector<double>, NDIM> DU_jump_qp;

    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    FEDataManager::SystemDofMapCache& X_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(COORDS_SYSTEM_NAME);
    FEType X_fe_type = X_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
    }

    System& P_system = equation_systems->get_system(P_SYSTEM_NAME);
    const DofMap& P_dof_map = P_system.get_dof_map();
    FEDataManager::SystemDofMapCache& P_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(P_SYSTEM_NAME);
    FEType P_fe_type = P_dof_map.variable_type(0);

    System& P_jump_system = equation_systems->get_system(PRESSURE_JUMP_SYSTEM_NAME);
    const DofMap& P_jump_dof_map = P_jump_system.get_dof_map();
    FEDataManager::SystemDofMapCache& P_jump_dof_map_cache =
        *d_fe_data_managers[part]->getDofMapCache(PRESSURE_JUMP_SYSTEM_NAME);
    FEType P_jump_fe_type = P_jump_dof_map.variable_type(0);
    TBOX_ASSERT(P_jump_fe_type == P_fe_type);

    System& WSS_system = equation_systems->get_system(WSS_SYSTEM_NAME);
    const DofMap& WSS_dof_map = WSS_system.get_dof_map();
    FEDataManager::SystemDofMapCache& WSS_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(WSS_SYSTEM_NAME);
    FEType WSS_fe_type = WSS_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(WSS_dof_map.variable_type(d) == WSS_fe_type);
    }
    TBOX_ASSERT(WSS_fe_type == P_fe_type);

    std::array<DofMap*, NDIM> DU_jump_dof_map;
    std::array<FEDataManager::SystemDofMapCache*, NDIM> DU_jump_dof_map_cache;
    std::array<System*, NDIM> DU_jump_system;
    for (unsigned int i = 0; i < NDIM; ++i)
    {
        DU_jump_system[i] = &equation_systems->get_system(VELOCITY_JUMP_SYSTEM_NAME[i]);
        DU_jump_dof_map_cache[i] = d_fe_data_managers[part]->getDofMapCache(VELOCITY_JUMP_SYSTEM_NAME[i]);
        DU_jump_dof_map[i] = &DU_jump_system[i]->get_dof_map();
    }
    FEType DU_jump_fe_type = DU_jump_dof_map[0]->variable_type(0);
    for (unsigned int i = 0; i < NDIM; ++i)
    {
        for (unsigned int j = 0; j < NDIM; ++j)
        {
            TBOX_ASSERT(DU_jump_dof_map[i]->variable_type(j) == DU_jump_fe_type);
        }
    }
    TBOX_ASSERT(DU_jump_fe_type == P_fe_type);

    System& TAU_system = equation_systems->get_system(TAU_SYSTEM_NAME);
    const DofMap& TAU_dof_map = TAU_system.get_dof_map();
    FEDataManager::SystemDofMapCache& TAU_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(TAU_SYSTEM_NAME);
    FEType TAU_fe_type = TAU_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(TAU_dof_map.variable_type(d) == TAU_fe_type);
    }
    TBOX_ASSERT(TAU_fe_type == P_fe_type);

    std::unique_ptr<FEBase> fe_X = FEBase::build(dim, X_fe_type);
    const std::vector<double>& JxW = fe_X->get_JxW();
    const std::vector<std::vector<double> >& phi_X = fe_X->get_phi();
    std::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi_X;
    dphi_dxi_X[0] = &fe_X->get_dphidxi();
    if (NDIM > 2) dphi_dxi_X[1] = &fe_X->get_dphideta();

    std::unique_ptr<FEBase> fe_P = FEBase::build(dim, P_fe_type);
    const std::vector<std::vector<double> >& phi_P = fe_P->get_phi();

    X_ghost_vec->close();
    PetscVector<double>* X_petsc_vec = static_cast<PetscVector<double>*>(X_ghost_vec);
    Vec X_global_vec = X_petsc_vec->vec();
    Vec X_local_vec;
    VecGhostGetLocalForm(X_global_vec, &X_local_vec);
    double* X_local_soln;
    VecGetArray(X_local_vec, &X_local_soln);
    std::unique_ptr<NumericVector<double> > X0_vec = X_petsc_vec->clone();
    copy_and_synch(X_system.get_vector("INITIAL_COORDINATES"), *X0_vec);
    X0_vec->close();

    const std::vector<std::vector<Elem*> >& active_patch_element_map =
        d_fe_data_managers[part]->getActivePatchElementMap();

    boost::multi_array<double, 2> x_node, X_node, WSS_node, n_qp_node;
    boost::multi_array<double, 1> P_node, P_jump_node;
    std::array<boost::multi_array<double, 2>, NDIM> DU_jump_node;
    std::vector<double> x_qp, X0_qp;
    std::vector<double> P_qp, P_jump_qp, N_qp, WSS_qp, TAU_qp;
    std::array<VectorValue<double>, 2> dX_dxi, dx_dxi;
    VectorValue<double> n, N, x, X;

    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_fe_data_managers[part]->getLevelNumber());
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const std::vector<Elem*>& patch_elems = active_patch_element_map[local_patch_num];
        const size_t num_active_patch_elems = patch_elems.size();
        if (!num_active_patch_elems) continue;
        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);

        const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

        const double* const dx = pgeom->getDx();

        double diag_dis = 0.0;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            diag_dis += dx[d] * dx[d];
        }

        unsigned int n_qp_patch = 0;
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
            get_values_for_interpolation(x_node, *X_petsc_vec, X_local_soln, X_dof_indices);
            FEDataManager::updateInterpQuadratureRule(qrule, d_default_interp_spec, elem, x_node, patch_dx_min);
            n_qp_patch += qrule->n_points();
        }

        if (!n_qp_patch) continue;
        P_qp.resize(n_qp_patch);
        P_jump_qp.resize(n_qp_patch);

        x_qp.resize(NDIM * n_qp_patch);
        X0_qp.resize(NDIM * n_qp_patch);
        WSS_qp.resize(NDIM * n_qp_patch);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            DU_jump_qp[axis].resize(NDIM * n_qp_patch);
        }

        TAU_qp.resize(NDIM * n_qp_patch);
        N_qp.resize(NDIM * n_qp_patch);
        std::fill(N_qp.begin(), N_qp.end(), 0.0);
        std::fill(x_qp.begin(), x_qp.end(), 0.0);
        std::fill(WSS_qp.begin(), WSS_qp.end(), 0.0);
        std::fill(P_qp.begin(), P_qp.end(), 0.0);
        std::fill(P_jump_qp.begin(), P_jump_qp.end(), 0.0);
        std::fill(TAU_qp.begin(), TAU_qp.end(), 0.0);

        // Loop over the elements and compute the positions of the quadrature points.
        qrule.reset();
        unsigned int qp_offset = 0;
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
            const auto& WSS_dof_indices = WSS_dof_map_cache.dof_indices(elem);
            const auto& P_dof_indices = P_dof_map_cache.dof_indices(elem);
            const auto& P_jump_dof_indices = P_jump_dof_map_cache.dof_indices(elem);
            get_values_for_interpolation(P_jump_node, *P_jump_ghost_vec, P_jump_dof_indices[0]);
            get_values_for_interpolation(x_node, *X_petsc_vec, X_local_soln, X_dof_indices);
            get_values_for_interpolation(WSS_node, *WSS_ghost_vec, WSS_dof_indices);
            get_values_for_interpolation(P_node, *P_ghost_vec, P_dof_indices[0]);
            get_values_for_interpolation(X_node, *X0_vec, X_dof_indices);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const auto& DU_jump_dof_indices = DU_jump_dof_map_cache[d]->dof_indices(elem);
                get_values_for_interpolation(DU_jump_node[d], *DU_jump_ghost_vec[d], DU_jump_dof_indices);
            }

            const bool qrule_changed =
                FEDataManager::updateInterpQuadratureRule(qrule, d_default_interp_spec, elem, x_node, patch_dx_min);
            if (qrule_changed)
            {
                fe_X->attach_quadrature_rule(qrule.get());
                fe_P->attach_quadrature_rule(qrule.get());
            }
            fe_X->reinit(elem);
            fe_P->reinit(elem);

            const unsigned int n_node = elem->n_nodes();
            const unsigned int n_qp = qrule->n_points();

            // Zero out the values of X, du, and dv prior to accumulation.
            double* x_begin = &x_qp[NDIM * qp_offset];
            std::fill(x_begin, x_begin + NDIM * n_qp, 0.0);

            double* N_begin = &N_qp[NDIM * qp_offset];
            std::fill(N_begin, N_begin + NDIM * n_qp, 0.0);

            double* WSS_begin = &WSS_qp[NDIM * qp_offset];
            std::fill(WSS_begin, WSS_begin + NDIM * n_qp, 0.0);

            double* TAU_begin = &TAU_qp[NDIM * qp_offset];
            std::fill(TAU_begin, TAU_begin + NDIM * n_qp, 0.0);

            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                double* DU_jump_begin = &DU_jump_qp[axis][NDIM * qp_offset];
                std::fill(DU_jump_begin, DU_jump_begin + NDIM * n_qp, 0.0);
            }

            double* P_jump_begin = &P_jump_qp[qp_offset];
            std::fill(P_jump_begin, P_jump_begin + n_qp, 0.0);

            double* P_begin = &P_qp[qp_offset];
            std::fill(P_begin, P_begin + n_qp, 0.0);

            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                interpolate(X, qp, X_node, phi_X);
                interpolate(x, qp, x_node, phi_X);
                for (unsigned int k = 0; k < NDIM - 1; ++k)
                {
                    interpolate(dX_dxi[k], qp, X_node, *dphi_dxi_X[k]);
                    interpolate(dx_dxi[k], qp, x_node, *dphi_dxi_X[k]);
                }
                if (NDIM == 2)
                {
                    dX_dxi[1] = dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                }

                // Construct unit vectors in the reference and current
                // configurations.
                N = dX_dxi[0].cross(dX_dxi[1]);
                const double dA = N.norm();
                N = N.unit();
                n = dx_dxi[0].cross(dx_dxi[1]);
                const double da = n.norm();
                n = n.unit();

                for (unsigned int i = 0; i < NDIM; ++i)
                {
                    for (unsigned int k = 0; k < n_node; ++k)
                    {
                        const double& p_X = phi_X[k][qp];
                        x_qp[NDIM * (qp_offset + qp) + i] += x_node[k][i] * p_X;
                        const double& p_P = phi_P[k][qp];
                        WSS_qp[NDIM * (qp_offset + qp) + i] += (da / dA) * WSS_node[k][i] * p_P;
                    }
                    // TODO: This is confusing --- N_qp is set to equal n --- shouldn't this be n_qp?
                    N_qp[NDIM * (qp_offset + qp) + i] = n(i);
                }

                for (unsigned int k = 0; k < n_node; ++k)
                {
                    const double& p_P = phi_P[k][qp];
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            DU_jump_qp[axis][NDIM * (qp_offset + qp) + d] += DU_jump_node[axis][k][d] * p_P;
                        }
                    }
                }

                for (unsigned int k = 0; k < n_node; ++k)
                {
                    const double& p_P = phi_P[k][qp];
                    P_qp[qp_offset + qp] += (da / dA) * P_node[k] * p_P;
                    P_jump_qp[qp_offset + qp] += P_jump_node[k] * p_P;
                }
            }
            qp_offset += n_qp;
        }

        const Box<NDIM>& interp_box = patch->getBox();
        std::vector<int> local_indices;
        local_indices.clear();
        const int upper_bound = n_qp_patch;
        if (upper_bound == 0) continue;

        local_indices.reserve(upper_bound);
        for (unsigned int k = 0; k < n_qp_patch; ++k)
        {
            const double* const XX = &x_qp[NDIM * k];
            const Index<NDIM> i = IndexUtilities::getCellIndex(XX, patch_geom, interp_box);
            if (interp_box.contains(i)) local_indices.push_back(k);
        }

        const unsigned int nindices = static_cast<int>(local_indices.size());

        if (!local_indices.empty())
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (unsigned int k = 0; k < nindices; ++k)
                {
                    // Using the exterior traciton tau_e

                    TAU_qp[NDIM * local_indices[k] + axis] =
                        (WSS_qp[NDIM * local_indices[k] + axis] -
                         P_qp[local_indices[k]] * N_qp[NDIM * local_indices[k] + axis]);
                }
            }
        }

        // Loop over the elements and accumulate the right-hand-side values.
        qrule.reset();
        qp_offset = 0;
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
            const auto& WSS_dof_indices = WSS_dof_map_cache.dof_indices(elem);
            const auto& TAU_dof_indices = TAU_dof_map_cache.dof_indices(elem);
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                TAU_rhs_e[i].resize(static_cast<int>(TAU_dof_indices[i].size()));
            }

            get_values_for_interpolation(x_node, *X_petsc_vec, X_local_soln, X_dof_indices);
            const bool qrule_changed =
                FEDataManager::updateInterpQuadratureRule(qrule, d_default_interp_spec, elem, x_node, patch_dx_min);
            if (qrule_changed)
            {
                fe_X->attach_quadrature_rule(qrule.get());
                fe_P->attach_quadrature_rule(qrule.get());
            }
            fe_X->reinit(elem);
            fe_P->reinit(elem);
            const unsigned int n_qp = qrule->n_points();
            const size_t n_basis2 = WSS_dof_indices[0].size();
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const int idx = NDIM * (qp_offset + qp);
                for (unsigned int k = 0; k < n_basis2; ++k)
                {
                    const double p_JxW = phi_P[k][qp] * JxW[qp];
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        TAU_rhs_e[i](k) += TAU_qp[idx + i] * p_JxW;
                    }
                }
            }
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                auto dof_id_scratch = TAU_dof_indices[d];
                TAU_dof_map.constrain_element_vector(TAU_rhs_e[d], dof_id_scratch);
                TAU_rhs_vec->add_vector(TAU_rhs_e[d], dof_id_scratch);
            }
            qp_offset += n_qp;
        }
    }

    TAU_rhs_vec->close();

    d_fe_data_managers[part]->computeL2Projection(
        *TAU_vec, *TAU_rhs_vec, TAU_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);

    d_X_half_vecs[part]->close();
    d_X_current_vecs[part]->close();
    d_X_new_vecs[part]->close();
    d_TAU_half_vecs[part]->close();
    d_WSS_half_vecs[part]->close();
    d_P_jump_half_vecs[part]->close();
    d_P_half_vecs[part]->close();

    VecRestoreArray(X_local_vec, &X_local_soln);
    VecGhostRestoreLocalForm(X_global_vec, &X_local_vec);

    d_WSS_IB_ghost_vecs[part]->close();

    d_P_IB_ghost_vecs[part]->close();
    d_P_jump_IB_ghost_vecs[part]->close();
    d_X_IB_ghost_vecs[part]->close();

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_DU_jump_half_vecs[part][d]->close();
        d_DU_jump_IB_ghost_vecs[part][d]->close();
    }

    return;
} // computeFluidTraction

void
IBFESurfaceMethod::imposeJumpConditions(const int f_data_idx,
                                        PetscVector<double>& P_jump_ghost_vec,
                                        std::array<PetscVector<double>*, NDIM>& DU_jump_ghost_vec,
                                        PetscVector<double>& X_ghost_vec,
                                        const double /*data_time*/,
                                        const unsigned int part)
{
    // Extract the mesh.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // Extract the FE systems and DOF maps, and setup the FE object
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    FEDataManager::SystemDofMapCache& X_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(COORDS_SYSTEM_NAME);
    FEType X_fe_type = X_dof_map.variable_type(0);
    for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);

    System* P_jump_system;
    const DofMap* P_jump_dof_map;
    FEDataManager::SystemDofMapCache* P_jump_dof_map_cache;
    FEType P_jump_fe_type;
    if (d_use_pressure_jump_conditions)
    {
        P_jump_system = &equation_systems->get_system(PRESSURE_JUMP_SYSTEM_NAME);
        P_jump_dof_map = &P_jump_system->get_dof_map();
        P_jump_dof_map_cache = d_fe_data_managers[part]->getDofMapCache(PRESSURE_JUMP_SYSTEM_NAME);
        P_jump_fe_type = P_jump_dof_map->variable_type(0);
    }

    std::array<DofMap*, NDIM> DU_jump_dof_map;
    std::array<FEDataManager::SystemDofMapCache*, NDIM> DU_jump_dof_map_cache;
    std::array<System*, NDIM> DU_jump_system;
    FEType DU_jump_fe_type;
    if (d_use_velocity_jump_conditions)
    {
        for (unsigned int i = 0; i < NDIM; ++i)
        {
            DU_jump_system[i] = &equation_systems->get_system(VELOCITY_JUMP_SYSTEM_NAME[i]);
            DU_jump_dof_map_cache[i] = d_fe_data_managers[part]->getDofMapCache(VELOCITY_JUMP_SYSTEM_NAME[i]);
            DU_jump_dof_map[i] = &DU_jump_system[i]->get_dof_map();
            DU_jump_fe_type = DU_jump_dof_map[i]->variable_type(0);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                TBOX_ASSERT(DU_jump_dof_map[i]->variable_type(d) == DU_jump_fe_type);
            }
        }
    }

    FEType fe_type = X_fe_type;
    std::unique_ptr<FEBase> fe_X = FEBase::build(dim, fe_type);
    std::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi;
    dphi_dxi[0] = &fe_X->get_dphidxi();
    if (NDIM > 2) dphi_dxi[1] = &fe_X->get_dphideta();

    std::unique_ptr<FEBase> fe_P_jump = FEBase::build(dim, P_jump_fe_type);
    const std::vector<std::vector<double> >& phi_P_jump = fe_P_jump->get_phi();

    // Loop over the patches to impose jump conditions on the Eulerian grid.
    const std::vector<std::vector<Elem*> >& active_patch_element_map =
        d_fe_data_managers[part]->getActivePatchElementMap();
    const int level_num = d_fe_data_managers[part]->getLevelNumber();
    boost::multi_array<double, 1> P_jump_node;
    boost::multi_array<double, 2> x_node;
    std::array<boost::multi_array<double, 2>, NDIM> DU_jump_node;
    std::array<VectorValue<double>, 2> dx_dxi;
    VectorValue<double> n, jn;
    std::vector<libMesh::Point> X_node_cache, x_node_cache;
    IBTK::Point x_min, x_max;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    const IntVector<NDIM>& ratio = level->getRatio();
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const std::vector<Elem*>& patch_elems = active_patch_element_map[local_patch_num];
        const size_t num_active_patch_elems = patch_elems.size();
        if (num_active_patch_elems == 0) continue;

        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_data_idx);
        const Box<NDIM>& patch_box = patch->getBox();
        const CellIndex<NDIM>& patch_lower = patch_box.lower();
        std::array<Box<NDIM>, NDIM> side_ghost_boxes;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            side_ghost_boxes[d] = SideGeometry<NDIM>::toSideBox(f_data->getGhostBox(), d);
        }
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const x_lower = patch_geom->getXLower();
        const double* const dx = patch_geom->getDx();

        std::array<std::map<hier::Index<NDIM>, std::vector<libMesh::Point>, IndexOrder>, NDIM> intersection_points,
            intersection_ref_coords;
        std::array<std::map<hier::Index<NDIM>, std::vector<VectorValue<double> >, IndexOrder>, NDIM>
            intersection_normals;

        std::array<std::map<hier::Index<NDIM>, std::vector<libMesh::Point>, IndexOrder>, NDIM> intersection_u_points,
            intersection_u_ref_coords;
        std::array<std::map<hier::Index<NDIM>, std::vector<VectorValue<double> >, IndexOrder>, NDIM>
            intersection_u_normals;

        std::array<std::array<std::map<hier::Index<NDIM>, std::vector<libMesh::Point>, IndexOrder>, NDIM>, NDIM>
            intersectionSide_u_points, intersectionSide_u_ref_coords;
        std::array<std::array<std::map<hier::Index<NDIM>, std::vector<VectorValue<double> >, IndexOrder>, NDIM>, NDIM>
            intersectionSide_u_normals;

        // Loop over the elements.
        for (size_t e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
            get_values_for_interpolation(x_node, X_ghost_vec, X_dof_indices);
            if (d_use_pressure_jump_conditions)
            {
                const auto& P_jump_dof_indices = P_jump_dof_map_cache->dof_indices(elem);
                get_values_for_interpolation(P_jump_node, P_jump_ghost_vec, P_jump_dof_indices[0]);
            }
            if (d_use_velocity_jump_conditions)
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    const auto& DU_jump_dof_indices = DU_jump_dof_map_cache[axis]->dof_indices(elem);
                    get_values_for_interpolation(DU_jump_node[axis], *DU_jump_ghost_vec[axis], DU_jump_dof_indices);
                }
            }

            // Cache the nodal and physical coordinates of the side element,
            // determine the bounding box of the current configuration of the
            // element, and set the nodal coordinates to correspond to the
            // physical coordinates.
            const unsigned int n_nodes = elem->n_nodes();
            X_node_cache.resize(n_nodes);
            x_node_cache.resize(n_nodes);
            x_min = IBTK::Point::Constant(std::numeric_limits<double>::max());
            x_max = IBTK::Point::Constant(-std::numeric_limits<double>::max());
            for (unsigned int k = 0; k < n_nodes; ++k)
            {
                X_node_cache[k] = elem->point(k);
                libMesh::Point& x = x_node_cache[k];
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    x(d) = x_node[k][d];
                }
                if (d_perturb_fe_mesh_nodes)
                {
                    // Perturb the mesh configuration to keep the FE mesh nodes
                    // away from cell edges, nodes, and centers.
                    //
                    // This implies that we only have to deal with multiple
                    // intersections along element edges, and not at element
                    // nodes.
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        const int i_s = std::floor((x(d) - x_lower[d]) / dx[d]) + patch_lower[d];
                        for (int shift = 0; shift <= 2; ++shift)
                        {
                            const double x_s =
                                x_lower[d] + dx[d] * (static_cast<double>(i_s - patch_lower[d]) + 0.5 * shift);
                            const double tol = 1.0e-4 * dx[d];
                            if (x(d) <= x_s) x(d) = std::min(x_s - tol, x(d));
                            if (x(d) >= x_s) x(d) = std::max(x_s + tol, x(d));
                        }
                    }
                }
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    x_min[d] = std::min(x_min[d], x(d));
                    x_max[d] = std::max(x_max[d], x(d));
                }
                elem->point(k) = x;
            }
            Box<NDIM> box(IndexUtilities::getCellIndex(&x_min[0], grid_geom, ratio),
                          IndexUtilities::getCellIndex(&x_max[0], grid_geom, ratio));
            box.grow(IntVector<NDIM>(1));
            box = box * patch_box;

            // Loop over coordinate directions and look for intersections with
            // the background fluid grid.
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                Box<NDIM> extended_box = patch_box;
                extended_box.grow(IntVector<NDIM>(1));
                if (patch_geom->getTouchesRegularBoundary(axis, 1)) extended_box.upper(axis) += 1;

                Box<NDIM> side_boxes[NDIM];
                for (int d = 0; d < NDIM; ++d) side_boxes[d] = SideGeometry<NDIM>::toSideBox(extended_box, d);

                Box<NDIM> side_u_boxes[NDIM];
                for (int d = 0; d < NDIM; ++d) side_u_boxes[d] = SideGeometry<NDIM>::toSideBox(patch_box, d);

                // Setup a unit vector pointing in the coordinate direction of
                // interest.
                VectorValue<double> q;
                q(axis) = 1.0;

                // Loop over the relevant range of indices.
                Box<NDIM> axis_box = box;
                axis_box.lower(axis) = 0;
                axis_box.upper(axis) = 0;

                unsigned int SideDim[NDIM][NDIM - 1];
                for (unsigned int d = 0; d < NDIM; ++d)
                    for (unsigned int l = 0; l < NDIM - 1; ++l) SideDim[d][l] = (d + l + 1) % NDIM;

                for (BoxIterator<NDIM> b(axis_box); b; b++)
                {
                    const Index<NDIM>& i_c = b();
                    libMesh::Point r;
                    std::array<libMesh::Point, NDIM - 1> rs;

                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        r(d) = (d == axis ? 0.0 :
                                            x_lower[d] + dx[d] * (static_cast<double>(i_c(d) - patch_lower[d]) + 0.5));

                        for (unsigned int l = 0; l < NDIM - 1; ++l)
                        {
                            rs[l](d) =
                                (d == axis ? 0.0 :
                                             d == SideDim[axis][l] ?
                                             x_lower[d] + dx[d] * (static_cast<double>(i_c(d) - patch_lower[d])) :
                                             x_lower[d] + dx[d] * (static_cast<double>(i_c(d) - patch_lower[d]) +
                                                                   0.5)); // In 2D this
                        }
                    }

                    std::vector<std::pair<double, libMesh::Point> > intersections;
                    std::array<std::vector<std::pair<double, libMesh::Point> >, NDIM - 1> intersectionsSide;

                    static const double tolerance = sqrt(std::numeric_limits<double>::epsilon());

#if (NDIM == 2)
                    intersect_line_with_edge(intersections, static_cast<Edge*>(elem), r, q, tolerance);
#endif
#if (NDIM == 3)
                    intersect_line_with_face(intersections, static_cast<Face*>(elem), r, q, tolerance);
#endif
                    for (unsigned int l = 0; l < NDIM - 1; ++l)
                    {
#if (NDIM == 2)
                        intersect_line_with_edge(intersectionsSide[l], static_cast<Edge*>(elem), rs[l], q, tolerance);
#endif
#if (NDIM == 3)
                        intersect_line_with_face(intersectionsSide[l], static_cast<Face*>(elem), rs[l], q, tolerance);
#endif
                    }

                    if (d_use_pressure_jump_conditions)
                    {
                        for (unsigned int k = 0; k < intersections.size(); ++k)
                        {
                            const libMesh::Point x = r + intersections[k].first * q;
                            const libMesh::Point& xi = intersections[k].second;
                            SideIndex<NDIM> i_s(i_c, axis, 0);
                            i_s(axis) = boost::math::iround((x(axis) - x_lower[axis]) / dx[axis]) + patch_lower[axis];
                            if (extended_box.contains(i_s))
                            {
                                std::vector<libMesh::Point> ref_coords(1, xi);
                                fe_X->reinit(elem, &ref_coords);
                                fe_P_jump->reinit(elem, &ref_coords);
                                for (unsigned int l = 0; l < NDIM - 1; ++l)
                                {
                                    interpolate(dx_dxi[l], 0, x_node, *dphi_dxi[l]);
                                }
                                if (NDIM == 2)
                                {
                                    dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                                }
                                n = (dx_dxi[0].cross(dx_dxi[1])).unit();

                                // Make sure we haven't already found this
                                // intersection.
                                //
                                // (Because we are doing this in floating point
                                // arithmetic, we can't even count on the
                                // intersection being assigned to the same index!)
                                bool found_same_intersection_point = false;
                                for (int shift = -1; shift <= 1; ++shift)
                                {
                                    SideIndex<NDIM> i_s_prime = i_s;
                                    i_s_prime(axis) += shift;
                                    const std::vector<libMesh::Point>& candidate_coords =
                                        intersection_points[axis][i_s_prime];
                                    const std::vector<libMesh::Point>& candidate_ref_coords =
                                        intersection_ref_coords[axis][i_s_prime];
                                    const std::vector<VectorValue<double> >& candidate_normals =
                                        intersection_normals[axis][i_s_prime];

                                    checkDoubleCountingIntersection(axis,
                                                                    dx,
                                                                    n,
                                                                    x,
                                                                    xi,
                                                                    i_s,
                                                                    i_s_prime,
                                                                    candidate_coords,
                                                                    candidate_ref_coords,
                                                                    candidate_normals,
                                                                    found_same_intersection_point);
                                    if (found_same_intersection_point) break;
                                }

                                if (!found_same_intersection_point)
                                {
                                    // Evaluate the jump conditions and apply them
                                    // to the Eulerian grid.
                                    if (side_ghost_boxes[axis].contains(i_s))
                                    {
                                        const double C_p = interpolate(0, P_jump_node, phi_P_jump);
                                        const double sgn = n(axis) > 0.0 ? 1.0 : n(axis) < 0.0 ? -1.0 : 0.0;
                                        (*f_data)(i_s) += sgn * (C_p / dx[axis]);
                                    }

                                    // Keep track of the positions where we have
                                    // imposed jump conditions.
                                    intersection_points[axis][i_s].push_back(x);
                                    intersection_ref_coords[axis][i_s].push_back(xi);
                                    intersection_normals[axis][i_s].push_back(n);
                                }
                            }
                        }
                    }

                    if (d_use_velocity_jump_conditions)
                    {
                        for (unsigned int k = 0; k < intersections.size(); ++k)
                        {
                            libMesh::Point xu = r + intersections[k].first * q;
                            const libMesh::Point& xui = intersections[k].second;
                            SideIndex<NDIM> i_s_um(i_c, axis, 0);
                            Index<NDIM> i_c_neighbor = i_c;
                            i_c_neighbor(axis) += 1;

                            SideIndex<NDIM> i_s_up(i_c_neighbor, axis, 0);
                            i_s_up(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis] + 1.0)) +
                                           patch_lower[axis];
                            i_s_um(axis) =
                                static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis])) + patch_lower[axis];

                            if (side_u_boxes[axis].contains(i_s_up) && side_u_boxes[axis].contains(i_s_um))
                            {
                                std::vector<libMesh::Point> ref_coords(1, xui);
                                fe_X->reinit(elem, &ref_coords);
                                fe_P_jump->reinit(elem, &ref_coords);
                                for (unsigned int l = 0; l < NDIM - 1; ++l)
                                {
                                    interpolate(dx_dxi[l], 0, x_node, *dphi_dxi[l]);
                                }
                                if (NDIM == 2)
                                {
                                    dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                                }
                                n = (dx_dxi[0].cross(dx_dxi[1])).unit();

                                bool found_same_intersection_point = false;

                                for (int shift = -1; shift <= 1; ++shift)
                                {
                                    SideIndex<NDIM> i_s_prime = i_s_um;
                                    i_s_prime(axis) += shift;
                                    const std::vector<libMesh::Point>& candidate_coords =
                                        intersection_u_points[axis][i_s_prime];
                                    const std::vector<libMesh::Point>& candidate_ref_coords =
                                        intersection_u_ref_coords[axis][i_s_prime];
                                    const std::vector<VectorValue<double> >& candidate_normals =
                                        intersection_u_normals[axis][i_s_prime];

                                    checkDoubleCountingIntersection(axis,
                                                                    dx,
                                                                    n,
                                                                    xu,
                                                                    xui,
                                                                    i_s_um,
                                                                    i_s_prime,
                                                                    candidate_coords,
                                                                    candidate_ref_coords,
                                                                    candidate_normals,
                                                                    found_same_intersection_point);
                                    if (found_same_intersection_point) break;
                                }

                                if (!found_same_intersection_point)
                                {
                                    // imposed jump conditions.

                                    TBOX_ASSERT(i_s_um.getAxis() == i_s_up.getAxis());
                                    TBOX_ASSERT(i_s_up(axis) - i_s_um(axis) == 1);
                                    const double x_cell_bdry_um =
                                        x_lower[axis] +
                                        static_cast<double>(i_s_um(axis) - patch_lower[axis]) * dx[axis];

                                    const double x_cell_bdry_up =
                                        x_lower[axis] +
                                        static_cast<double>(i_s_up(axis) - patch_lower[axis]) * dx[axis];
                                    const double sdh_um = ((xu(axis) - x_cell_bdry_um)); // Signed Distance h

                                    const double sdh_up = ((xu(axis) - x_cell_bdry_up)); // Signed Distance h

                                    TBOX_ASSERT((sdh_um) < dx[axis] && sdh_um > 0);
                                    TBOX_ASSERT(fabs(sdh_up) < dx[axis] && sdh_up < 0);

                                    double C_u_um = 0;
                                    double C_u_up = 0;

                                    interpolate(&jn(0), 0, DU_jump_node[axis], phi_P_jump);
                                    C_u_up = sdh_up * jn(axis);
                                    C_u_um = sdh_um * jn(axis);

                                    const double sgn = n(axis) > 0.0 ? 1.0 : n(axis) < 0.0 ? -1.0 : 0.0;
                                    // Note that the corrections are applied to opposite sides
                                    (*f_data)(i_s_up) -= sgn * (C_u_um / (dx[axis] * dx[axis]));
                                    (*f_data)(i_s_um) += sgn * (C_u_up / (dx[axis] * dx[axis]));
                                }
                            }
                            // Keep track of the positions where we have
                            intersection_u_points[axis][i_s_um].push_back(xu);
                            intersection_u_ref_coords[axis][i_s_um].push_back(xui);
                            intersection_u_normals[axis][i_s_um].push_back(n);
                        }

                        for (unsigned int j = 0; j < NDIM - 1; ++j)
                        {
                            for (unsigned int k = 0; k < intersectionsSide[j].size(); ++k)
                            {
                                libMesh::Point xu = rs[j] + intersectionsSide[j][k].first * q;
                                const libMesh::Point& xui = intersectionsSide[j][k].second;
                                SideIndex<NDIM> i_s_up;
                                SideIndex<NDIM> i_s_um;

                                if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) >= 0.5 * dx[axis])
                                {
                                    SideIndex<NDIM> i_side_um(i_c, SideDim[axis][j], 0);
                                    Index<NDIM> i_c_neighbor = i_c;
                                    i_c_neighbor(axis) += 1;

                                    SideIndex<NDIM> i_side_up(i_c_neighbor, SideDim[axis][j], 0);

                                    i_side_up(axis) =
                                        static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis] + 0.5)) +
                                        patch_lower[axis];
                                    i_side_um(axis) =
                                        static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis])) +
                                        patch_lower[axis];
                                    i_s_up = i_side_up;
                                    i_s_um = i_side_um;
                                }
                                else if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) < 0.5 * dx[axis] &&
                                         fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) >= 0.0)
                                {
                                    SideIndex<NDIM> i_side_up(i_c, SideDim[axis][j], 0);
                                    Index<NDIM> i_c_neighbor = i_c;
                                    i_c_neighbor(axis) -= 1;
                                    SideIndex<NDIM> i_side_um(i_c_neighbor, SideDim[axis][j], 0);
                                    i_side_up(axis) =
                                        static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis])) +
                                        patch_lower[axis];
                                    i_side_um(axis) =
                                        static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis] - 0.5)) +
                                        patch_lower[axis];
                                    i_s_up = i_side_up;
                                    i_s_um = i_side_um;
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name << ":  Restart file version different than class version."
                                                             << std::endl);
                                }

                                if (side_u_boxes[axis].contains(i_s_up) && side_u_boxes[axis].contains(i_s_um))
                                {
                                    std::vector<libMesh::Point> ref_coords(1, xui);
                                    fe_X->reinit(elem, &ref_coords);
                                    fe_P_jump->reinit(elem, &ref_coords);
                                    for (unsigned int l = 0; l < NDIM - 1; ++l)
                                    {
                                        interpolate(dx_dxi[l], 0, x_node, *dphi_dxi[l]);
                                    }
                                    if (NDIM == 2)
                                    {
                                        dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                                    }
                                    n = (dx_dxi[0].cross(dx_dxi[1])).unit();

                                    bool found_same_intersection_point = false;

                                    for (int shift = -1; shift <= 1; ++shift)
                                    {
                                        SideIndex<NDIM> i_s_prime = i_s_um;
                                        i_s_prime(SideDim[axis][j]) += shift;
                                        const std::vector<libMesh::Point>& candidate_coords =
                                            intersectionSide_u_points[j][axis][i_s_prime];
                                        const std::vector<libMesh::Point>& candidate_ref_coords =
                                            intersectionSide_u_ref_coords[j][axis][i_s_prime];
                                        const std::vector<VectorValue<double> >& candidate_normals =
                                            intersectionSide_u_normals[j][axis][i_s_prime];

                                        checkDoubleCountingIntersection(axis,
                                                                        dx,
                                                                        n,
                                                                        xu,
                                                                        xui,
                                                                        i_s_um,
                                                                        i_s_prime,
                                                                        candidate_coords,
                                                                        candidate_ref_coords,
                                                                        candidate_normals,
                                                                        found_same_intersection_point);
                                        if (found_same_intersection_point) break;
                                    }

                                    if (!found_same_intersection_point)
                                    {
                                        // Evaluate the jump conditions and apply them
                                        // to the Eulerian grid.

                                        // Evaluate the jump conditions and apply them
                                        // to the Eulerian grid.
                                        if (side_u_boxes[SideDim[axis][j]].contains(i_s_um) &&
                                            side_u_boxes[SideDim[axis][j]].contains(i_s_up))
                                        {
                                            TBOX_ASSERT(i_s_up(axis) - i_s_um(axis) == 1);

                                            const double x_mid_side_up =
                                                x_lower[axis] +
                                                static_cast<double>(i_s_up(axis) - patch_lower[axis] + 0.5) * dx[axis];

                                            const double x_mid_side_um =
                                                x_lower[axis] +
                                                static_cast<double>(i_s_um(axis) - patch_lower[axis] + 0.5) * dx[axis];

                                            TBOX_ASSERT(xu(axis) <= x_mid_side_up);
                                            TBOX_ASSERT(xu(axis) > x_mid_side_um);

                                            const double sdh_up = xu(axis) - x_mid_side_up; // Signed Distance h
                                            const double sdh_um = xu(axis) - x_mid_side_um;

                                            double C_u_um = 0;
                                            double C_u_up = 0;

                                            interpolate(&jn(0), 0, DU_jump_node[SideDim[axis][j]], phi_P_jump);
                                            C_u_um = sdh_um * jn(axis);
                                            C_u_up = sdh_up * jn(axis);

                                            const double sgn = n(axis) > 0.0 ? 1.0 : n(axis) < 0.0 ? -1.0 : 0.0;

                                            (*f_data)(i_s_um) += sgn * (C_u_up / (dx[axis] * dx[axis]));

                                            (*f_data)(i_s_up) -= sgn * (C_u_um / (dx[axis] * dx[axis]));
                                        }

                                        intersectionSide_u_points[j][axis][i_s_um].push_back(xu);
                                        intersectionSide_u_ref_coords[j][axis][i_s_um].push_back(xui);
                                        intersectionSide_u_normals[j][axis][i_s_um].push_back(n);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Restore the element coordinates.
            for (unsigned int k = 0; k < n_nodes; ++k)
            {
                elem->point(k) = X_node_cache[k];
            }
        }
    }

    return;
} // imposeJumpConditions

void
IBFESurfaceMethod::checkDoubleCountingIntersection(int axis,
                                                   const double* dx,
                                                   libMesh::VectorValue<double> n,
                                                   const libMesh::Point x,
                                                   const libMesh::Point& xi,
                                                   const SideIndex<NDIM> i_s,
                                                   const SideIndex<NDIM> i_s_prime,
                                                   const std::vector<libMesh::Point> candidate_coords,
                                                   const std::vector<libMesh::Point> candidate_ref_coords,
                                                   const std::vector<libMesh::VectorValue<double> > candidate_normals,
                                                   bool found_same_intersection_point)
{
    std::vector<libMesh::Point>::const_iterator x_prime_it = candidate_coords.begin();
    std::vector<libMesh::Point>::const_iterator xi_prime_it = candidate_ref_coords.begin();
    std::vector<VectorValue<double> >::const_iterator n_prime_it = candidate_normals.begin();
    for (; x_prime_it != candidate_coords.end(); ++x_prime_it, ++xi_prime_it, ++n_prime_it)
    {
        const libMesh::Point& x_prime = *x_prime_it;
        const libMesh::Point& xi_prime = *xi_prime_it;
        const libMesh::Point& n_prime = *n_prime_it;
        // TODO: Do not use a hard-coded magic number?
        if (x.absolute_fuzzy_equals(x_prime, 1.0e-5 * dx[axis]))
        {
            // WARNING: This check is ONLY
            // guaranteed to work at edges (where
            // only two elements meet).  To avoid FE
            // mesh nodes, set
            // d_perturb_fe_mesh_nodes to true.
            found_same_intersection_point = n(axis) * n_prime(axis) > 0.0;
            if (d_do_log)
            {
                plog << "==========\n";
                plog << "multiple intersections detected:\n";
                plog << "  x    = " << x << "\n";
                plog << "  x'   = " << x_prime << "\n";
                plog << "  xi   = " << xi << "\n";
                plog << "  xi'  = " << xi_prime << "\n";
                plog << "  n    = " << n << "\n";
                plog << "  n'   = " << n_prime << "\n";
                plog << "  i_s  = " << i_s << "\n";
                plog << "  i_s' = " << i_s_prime << "\n";
                plog << "  axis = " << axis << "\n";
            }
        }
        if (found_same_intersection_point) break;
    }
    return;
}

void
IBFESurfaceMethod::initializeCoordinates(const unsigned int part)
{
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    const bool identity_mapping = !d_coordinate_mapping_fcn_data[part].fcn;

    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num))
        {
            TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
            const libMesh::Point& X = *n;
            libMesh::Point x = X;
            if (!identity_mapping)
            {
                d_coordinate_mapping_fcn_data[part].fcn(x, X, d_coordinate_mapping_fcn_data[part].ctx);
            }
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int dof_index = n->dof_number(X_sys_num, d, 0);
                X_coords.set(dof_index, x(d));
            }
        }
    }
    X_coords.close();
    X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_coords);
    copy_and_synch(X_coords, *X_system.current_local_solution);
    copy_and_synch(X_coords, X_system.get_vector("INITIAL_COORDINATES"));
    return;
} // initializeCoordinates

void
IBFESurfaceMethod::updateCoordinateMapping(const unsigned int part)
{
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    System& dX_system = equation_systems->get_system(COORD_MAPPING_SYSTEM_NAME);
    const unsigned int dX_sys_num = dX_system.number();
    NumericVector<double>& dX_coords = *dX_system.solution;

    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num))
        {
            TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
            TBOX_ASSERT(n->n_vars(dX_sys_num) == NDIM);
            const libMesh::Point& X = *n;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int X_dof_index = n->dof_number(X_sys_num, d, 0);
                const int dX_dof_index = n->dof_number(dX_sys_num, d, 0);
                dX_coords.set(dX_dof_index, X_coords(X_dof_index) - X(d));
            }
        }
    }
    copy_and_synch(dX_coords, *dX_system.current_local_solution);
    return;
} // updateCoordinateMapping

void
IBFESurfaceMethod::initializeVelocity(const unsigned int part)
{
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& U_system = equation_systems->get_system(VELOCITY_SYSTEM_NAME);
    const unsigned int U_sys_num = U_system.number();
    NumericVector<double>& U_vec = *U_system.solution;
    VectorValue<double> U;
    if (!d_initial_velocity_fcn_data[part].fcn)
    {
        U.zero();
    }
    else
    {
        for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
        {
            Node* n = *it;
            if (n->n_vars(U_sys_num))
            {
                TBOX_ASSERT(n->n_vars(U_sys_num) == NDIM);
                const libMesh::Point& X = *n;
                d_initial_velocity_fcn_data[part].fcn(U, X, d_initial_velocity_fcn_data[part].ctx);
            }
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int dof_index = n->dof_number(U_sys_num, d, 0);
                U_vec.set(dof_index, U(d));
            }
        }
    }
    U_vec.close();
    U_system.get_dof_map().enforce_constraints_exactly(U_system, &U_vec);
    copy_and_synch(U_vec, *U_system.current_local_solution);
    return;
} // initializeVelocity

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBFESurfaceMethod::commonConstructor(const std::string& object_name,
                                     Pointer<Database> input_db,
                                     const std::vector<libMesh::MeshBase*>& meshes,
                                     int max_level_number,
                                     bool register_for_restart,
                                     const std::string& restart_read_dirname,
                                     unsigned int restart_restore_number)
{
    // Set the object name and register it with the restart manager.
    d_object_name = object_name;
    d_registered_for_restart = false;
    if (register_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
        d_registered_for_restart = true;
    }
    d_libmesh_restart_read_dir = restart_read_dirname;
    d_libmesh_restart_restore_number = restart_restore_number;

    // Store the mesh pointers.
    d_meshes = meshes;
    d_max_level_number = max_level_number;

    // Set some default values.
    const bool use_adaptive_quadrature = true;
    const int point_density = 2.0;
    d_wss_calc_width = 0.0;
    d_p_calc_width = 0.0;
    d_traction_activation_time = 0.0;
    const bool use_nodal_quadrature = false;
    const bool interp_use_consistent_mass_matrix = true;
    d_default_interp_spec = FEDataManager::InterpSpec("IB_4",
                                                      QGAUSS,
                                                      INVALID_ORDER,
                                                      use_adaptive_quadrature,
                                                      point_density,
                                                      interp_use_consistent_mass_matrix,
                                                      use_nodal_quadrature);
    d_default_spread_spec = FEDataManager::SpreadSpec(
        "IB_4", QGAUSS, INVALID_ORDER, use_adaptive_quadrature, point_density, use_nodal_quadrature);
    d_ghosts = 0;
    d_use_velocity_jump_conditions = false;
    d_use_pressure_jump_conditions = false;
    d_use_l2_lagrange_family = false;
    d_compute_fluid_traction = false;

    d_fe_family.resize(d_num_parts, INVALID_FE);
    d_fe_order.resize(d_num_parts, INVALID_ORDER);
    d_default_quad_type.resize(d_num_parts, INVALID_Q_RULE);
    d_default_quad_order.resize(d_num_parts, INVALID_ORDER);

    // Initialize function data to NULL.
    d_coordinate_mapping_fcn_data.resize(d_num_parts);
    d_initial_velocity_fcn_data.resize(d_num_parts);
    d_lag_surface_pressure_fcn_data.resize(d_num_parts);
    d_lag_surface_force_fcn_data.resize(d_num_parts);
    d_lag_surface_force_integral.resize(d_num_parts);

    // Determine whether we should use first-order or second-order shape
    // functions for each part of the structure.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        const MeshBase& mesh = *meshes[part];
        bool mesh_has_first_order_elems = false;
        bool mesh_has_second_order_elems = false;
        MeshBase::const_element_iterator el_it = mesh.elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.elements_end();
        for (; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            mesh_has_first_order_elems = mesh_has_first_order_elems || elem->default_order() == FIRST;
            mesh_has_second_order_elems = mesh_has_second_order_elems || elem->default_order() == SECOND;
        }
        mesh_has_first_order_elems = SAMRAI_MPI::maxReduction(mesh_has_first_order_elems);
        mesh_has_second_order_elems = SAMRAI_MPI::maxReduction(mesh_has_second_order_elems);
        if ((mesh_has_first_order_elems && mesh_has_second_order_elems) ||
            (!mesh_has_first_order_elems && !mesh_has_second_order_elems))
        {
            TBOX_ERROR(d_object_name
                       << "::IBFESurfaceMethod():\n"
                       << "  each FE mesh part must contain only FIRST order elements or only SECOND order elements"
                       << std::endl);
        }
        d_fe_family[part] = LAGRANGE;
        d_default_quad_type[part] = QGAUSS;
        if (mesh_has_first_order_elems)
        {
            d_fe_order[part] = FIRST;
            d_default_quad_order[part] = THIRD;
        }
        if (mesh_has_second_order_elems)
        {
            d_fe_order[part] = SECOND;
            d_default_quad_order[part] = FIFTH;
        }

        // Report configuration.
        pout << "\n";
        pout << d_object_name << ": mesh part " << part << " is using "
             << Utility::enum_to_string<Order>(d_fe_order[part]) << " order "
             << Utility::enum_to_string<FEFamily>(d_fe_family[part]) << " finite elements.\n";
        pout << "\n";
    }

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Set up the interaction spec objects.
    d_interp_spec.resize(d_num_parts, d_default_interp_spec);
    d_spread_spec.resize(d_num_parts, d_default_spread_spec);

    return;
} // commonConstructor

void
IBFESurfaceMethod::getFromInput(Pointer<Database> db, bool /*is_from_restart*/)
{
    // Interpolation settings.
    if (db->isString("interp_delta_fcn"))
        d_default_interp_spec.kernel_fcn = db->getString("interp_delta_fcn");
    else if (db->isString("IB_delta_fcn"))
        d_default_interp_spec.kernel_fcn = db->getString("IB_delta_fcn");
    else if (db->isString("interp_kernel_fcn"))
        d_default_interp_spec.kernel_fcn = db->getString("interp_kernel_fcn");
    else if (db->isString("IB_kernel_fcn"))
        d_default_interp_spec.kernel_fcn = db->getString("IB_kernel_fcn");

    if (db->isString("interp_quad_type"))
        d_default_interp_spec.quad_type = Utility::string_to_enum<QuadratureType>(db->getString("interp_quad_type"));
    else if (db->isString("IB_quad_type"))
        d_default_interp_spec.quad_type = Utility::string_to_enum<QuadratureType>(db->getString("IB_quad_type"));

    if (db->isString("interp_quad_order"))
        d_default_interp_spec.quad_order = Utility::string_to_enum<Order>(db->getString("interp_quad_order"));
    else if (db->isString("IB_quad_order"))
        d_default_interp_spec.quad_order = Utility::string_to_enum<Order>(db->getString("IB_quad_order"));

    if (db->isBool("interp_use_adaptive_quadrature"))
        d_default_interp_spec.use_adaptive_quadrature = db->getBool("interp_use_adaptive_quadrature");
    else if (db->isBool("IB_use_adaptive_quadrature"))
        d_default_interp_spec.use_adaptive_quadrature = db->getBool("IB_use_adaptive_quadrature");

    if (db->isDouble("interp_point_density"))
        d_default_interp_spec.point_density = db->getDouble("interp_point_density");
    else if (db->isDouble("IB_point_density"))
        d_default_interp_spec.point_density = db->getDouble("IB_point_density");

    if (db->isBool("interp_use_consistent_mass_matrix"))
        d_default_interp_spec.use_consistent_mass_matrix = db->getBool("interp_use_consistent_mass_matrix");
    else if (db->isBool("IB_use_consistent_mass_matrix"))
        d_default_interp_spec.use_consistent_mass_matrix = db->getBool("IB_use_consistent_mass_matrix");

    if (db->isBool("interp_use_nodal_quadrature"))
        d_default_interp_spec.use_nodal_quadrature = db->getBool("interp_use_nodal_quadrature");
    else if (db->isBool("IB_use_nodal_quadrature"))
        d_default_interp_spec.use_nodal_quadrature = db->getBool("IB_use_nodal_quadrature");

    // Spreading settings.
    if (db->isString("spread_delta_fcn"))
        d_default_spread_spec.kernel_fcn = db->getString("spread_delta_fcn");
    else if (db->isString("IB_delta_fcn"))
        d_default_spread_spec.kernel_fcn = db->getString("IB_delta_fcn");
    else if (db->isString("spread_kernel_fcn"))
        d_default_spread_spec.kernel_fcn = db->getString("spread_kernel_fcn");
    else if (db->isString("IB_kernel_fcn"))
        d_default_spread_spec.kernel_fcn = db->getString("IB_kernel_fcn");

    if (db->isString("spread_quad_type"))
        d_default_spread_spec.quad_type = Utility::string_to_enum<QuadratureType>(db->getString("spread_quad_type"));
    else if (db->isString("IB_quad_type"))
        d_default_spread_spec.quad_type = Utility::string_to_enum<QuadratureType>(db->getString("IB_quad_type"));

    if (db->isString("spread_quad_order"))
        d_default_spread_spec.quad_order = Utility::string_to_enum<Order>(db->getString("spread_quad_order"));
    else if (db->isString("IB_quad_order"))
        d_default_spread_spec.quad_order = Utility::string_to_enum<Order>(db->getString("IB_quad_order"));

    if (db->isBool("spread_use_adaptive_quadrature"))
        d_default_spread_spec.use_adaptive_quadrature = db->getBool("spread_use_adaptive_quadrature");
    else if (db->isBool("IB_use_adaptive_quadrature"))
        d_default_spread_spec.use_adaptive_quadrature = db->getBool("IB_use_adaptive_quadrature");

    if (db->isDouble("spread_point_density"))
        d_default_spread_spec.point_density = db->getDouble("spread_point_density");
    else if (db->isDouble("IB_point_density"))
        d_default_spread_spec.point_density = db->getDouble("IB_point_density");

    if (db->isBool("spread_use_nodal_quadrature"))
        d_default_spread_spec.use_nodal_quadrature = db->getBool("spread_use_nodal_quadrature");
    else if (db->isBool("IB_use_nodal_quadrature"))
        d_default_spread_spec.use_nodal_quadrature = db->getBool("IB_use_nodal_quadrature");

    // Force computation settings.
    if (db->isBool("use_pressure_jump_conditions"))
        d_use_pressure_jump_conditions = db->getBool("use_pressure_jump_conditions");

    if (db->isBool("use_l2_lagrange_family")) d_use_l2_lagrange_family = db->getBool("use_l2_lagrange_family");

    if (db->isBool("use_velocity_jump_conditions"))
        d_use_velocity_jump_conditions = db->getBool("use_velocity_jump_conditions");
    if (d_use_pressure_jump_conditions || d_use_velocity_jump_conditions)
    {
        if (db->isBool("perturb_fe_mesh_nodes")) d_perturb_fe_mesh_nodes = db->getBool("perturb_fe_mesh_nodes");
        if (db->isBool("normalize_pressure_jump")) d_normalize_pressure_jump = db->getBool("normalize_pressure_jump");
    }
    if (d_use_velocity_jump_conditions && d_use_pressure_jump_conditions)
    {
        if (db->isDouble("p_calc_width")) d_p_calc_width = db->getDouble("p_calc_width");
        if (db->isDouble("wss_calc_width")) d_wss_calc_width = db->getDouble("wss_calc_width");
        if (db->isDouble("traction_activation_time"))
            d_traction_activation_time = db->getDouble("traction_activation_time");
        if (db->isBool("compute_fluid_traction")) d_compute_fluid_traction = db->getBool("compute_fluid_traction");
        if (db->isBool("traction_interior_side")) d_traction_interior_side = db->getBool("traction_interior_side");
    }
    if (db->isBool("use_consistent_mass_matrix"))
        d_use_consistent_mass_matrix = db->getBool("use_consistent_mass_matrix");
    if (db->isBool("use_direct_forcing")) d_use_direct_forcing = db->getBool("use_direct_forcing");

    // Restart settings.
    if (db->isString("libmesh_restart_file_extension"))
    {
        d_libmesh_restart_file_extension = db->getString("libmesh_restart_file_extension");
    }
    else
    {
        d_libmesh_restart_file_extension = "xdr";
    }

    // Other settings.
    if (db->isInteger("min_ghost_cell_width"))
    {
        d_ghosts = db->getInteger("min_ghost_cell_width");
    }
    else if (db->isDouble("min_ghost_cell_width"))
    {
        d_ghosts = static_cast<int>(std::ceil(db->getDouble("min_ghost_cell_width")));
    }
    if (db->keyExists("do_log"))
        d_do_log = db->getBool("do_log");
    else if (db->keyExists("enable_logging"))
        d_do_log = db->getBool("enable_logging");
    return;
} // getFromInput

void
IBFESurfaceMethod::getFromRestart()
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
                                 << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("IBFE_METHOD_VERSION");
    if (ver != IBFE_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    db->getIntegerArray("d_ghosts", d_ghosts, NDIM);
    d_use_l2_lagrange_family = db->getBool("d_use_l2_lagrange_family");
    d_use_pressure_jump_conditions = db->getBool("d_use_pressure_jump_conditions");
    d_use_velocity_jump_conditions = db->getBool("d_use_velocity_jump_conditions");
    d_compute_fluid_traction = db->getBool("d_compute_fluid_traction");
    d_traction_interior_side = db->getBool("d_d_traction_interior_side");
    d_traction_activation_time = db->getDouble("d_traction_activation_time");
    d_use_consistent_mass_matrix = db->getBool("d_use_consistent_mass_matrix");
    d_use_direct_forcing = db->getBool("d_use_direct_forcing");
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
