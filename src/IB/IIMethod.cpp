// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2019 by the IBAMR developers
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

#include "ibamr/IIMethod.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/FEDataInterpolation.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/SAMRAIDataCache.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "GriddingAlgorithm.h"
#include "Index.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "Patch.h"
#include "PatchData.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineSchedule.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
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
#include "libmesh/enum_parallel_type.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/equation_systems.h"
#include "libmesh/face.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/fem_context.h"
#include "libmesh/id_types.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_version.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/system.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_vector.h"
#include "libmesh/variant_filter_iterator.h"
#include "libmesh/vector_value.h"

#include "petscvec.h"

#include "ibamr/namespaces.h" // IWYU pragma: keep

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/math/special_functions/round.hpp>
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

namespace SAMRAI
{
namespace xfer
{
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
// Version of IIMethod restart file data.
static const int IBFE_METHOD_VERSION = 2;

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

const std::string IIMethod::COORDS_SYSTEM_NAME = "coordinates system";
const std::string IIMethod::COORD_MAPPING_SYSTEM_NAME = "coordinate mapping system";
const std::string IIMethod::FORCE_SYSTEM_NAME = "IB force system";
const std::string IIMethod::VELOCITY_SYSTEM_NAME = "velocity system";
const std::string IIMethod::NORMAL_VELOCITY_SYSTEM_NAME = "normal velocity system";
const std::string IIMethod::TANGENTIAL_VELOCITY_SYSTEM_NAME = "tangential velocity system";
const std::string IIMethod::PRESSURE_JUMP_SYSTEM_NAME = "[[p]] system";
const std::string IIMethod::WSS_IN_SYSTEM_NAME = "One sided interior wall shear stress system";
const std::string IIMethod::WSS_OUT_SYSTEM_NAME = "One sided exterior wall shear stress system";
const std::string IIMethod::PRESSURE_IN_SYSTEM_NAME = "One sided interior pressure system";
const std::string IIMethod::PRESSURE_OUT_SYSTEM_NAME = "One sided exterior pressure system";
const std::string IIMethod::TAU_IN_SYSTEM_NAME = "Interior traction system";
const std::string IIMethod::TAU_OUT_SYSTEM_NAME = "Exterior traction system";
const std::array<std::string, NDIM> IIMethod::VELOCITY_JUMP_SYSTEM_NAME = { { "velocity [[du]] jump system",
                                                                              "velocity [[dv]] jump system"
#if (NDIM == 3)
                                                                              ,
                                                                              "velocity [[dw]] jump system"
#endif
} };

/////////////////////////////// PUBLIC ///////////////////////////////////////

IIMethod::IIMethod(const std::string& object_name,
                   Pointer<Database> input_db,
                   MeshBase* mesh,
                   int max_level_number,
                   bool register_for_restart,
                   const std::string& restart_read_dirname,
                   unsigned int restart_restore_number)
{
    commonConstructor(object_name,
                      input_db,
                      std::vector<MeshBase*>(1, mesh),
                      max_level_number,
                      register_for_restart,
                      restart_read_dirname,
                      restart_restore_number);
    return;
} // IIMethod

IIMethod::IIMethod(const std::string& object_name,
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
} // IIMethod

IIMethod::~IIMethod()
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
} // ~IIMethod

FEDataManager*
IIMethod::getFEDataManager(const unsigned int part) const
{
    TBOX_ASSERT(d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_num_parts);
    return d_fe_data_managers[part];
} // getFEDataManager

void
IIMethod::registerDisconElemFamilyForJumps(const unsigned int part)
{
    TBOX_ASSERT(!d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_num_parts);
    d_use_discon_elem_for_jumps[part] = true;
    return;
} // registerDisconElemFamilyForJumps

void
IIMethod::registerTangentialVelocityMotion(const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_use_tangential_velocity[part] = true;
    return;
} // registerTangentialVelocityMotion

void
IIMethod::registerPressureJumpNormalization(const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    TBOX_ASSERT(d_use_pressure_jump_conditions);
    d_normalize_pressure_jump[part] = true;
    return;
} // registerPressureJumpNormalization

void
IIMethod::registerInitialCoordinateMappingFunction(const CoordinateMappingFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_coordinate_mapping_fcn_data[part] = data;
    return;
} // registerInitialCoordinateMappingFunction

void
IIMethod::registerInitialVelocityFunction(const InitialVelocityFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_initial_velocity_fcn_data[part] = data;
    return;
} // registerInitialVelocityFunction

void
IIMethod::registerLagSurfacePressureFunction(const LagSurfacePressureFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_surface_pressure_fcn_data[part] = data;
    return;
} // registerLagSurfacePressureFunction

void
IIMethod::registerLagSurfaceForceFunction(const LagSurfaceForceFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_surface_force_fcn_data[part] = data;
    return;
} // registerLagSurfaceForceFunction

const VectorValue<double>&
IIMethod::getSurfaceForceIntegral(const unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_lag_surface_force_integral[part];
} // getSurfaceForceIntegral

const IntVector<NDIM>&
IIMethod::getMinimumGhostCellWidth() const
{
    return d_ghosts;
} // getMinimumGhostCellWidth

void
IIMethod::setupTagBuffer(Array<int>& tag_buffer, Pointer<GriddingAlgorithm<NDIM> > gridding_alg) const
{
    const int finest_hier_ln = gridding_alg->getMaxLevels() - 1;
    const int tsize = tag_buffer.size();
    tag_buffer.resizeArray(finest_hier_ln);
    for (int i = tsize; i < finest_hier_ln; ++i) tag_buffer[i] = 0;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        const int gcw = d_fe_data_managers[part]->getGhostCellWidth().max();
        const int tag_ln = d_fe_data_managers[part]->getFinestPatchLevelNumber() - 1;
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
IIMethod::preprocessIntegrateData(double current_time, double new_time, int /*num_cycles*/)
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

    d_P_in_systems.resize(d_num_parts);
    d_P_in_half_vecs.resize(d_num_parts);
    d_P_in_IB_ghost_vecs.resize(d_num_parts);

    d_P_out_systems.resize(d_num_parts);
    d_P_out_half_vecs.resize(d_num_parts);
    d_P_out_IB_ghost_vecs.resize(d_num_parts);

    d_DU_jump_systems.resize(d_num_parts);
    d_DU_jump_half_vecs.resize(d_num_parts);
    d_DU_jump_IB_ghost_vecs.resize(d_num_parts);

    d_WSS_in_systems.resize(d_num_parts);
    d_WSS_in_half_vecs.resize(d_num_parts);
    d_WSS_in_IB_ghost_vecs.resize(d_num_parts);

    d_WSS_out_systems.resize(d_num_parts);
    d_WSS_out_half_vecs.resize(d_num_parts);
    d_WSS_out_IB_ghost_vecs.resize(d_num_parts);

    d_TAU_in_systems.resize(d_num_parts);
    d_TAU_in_half_vecs.resize(d_num_parts);
    d_TAU_in_IB_ghost_vecs.resize(d_num_parts);

    d_TAU_out_systems.resize(d_num_parts);
    d_TAU_out_half_vecs.resize(d_num_parts);
    d_TAU_out_IB_ghost_vecs.resize(d_num_parts);

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

            d_P_in_systems[part] = &d_equation_systems[part]->get_system(PRESSURE_IN_SYSTEM_NAME);
            d_P_in_half_vecs[part] =
                dynamic_cast<PetscVector<double>*>(d_P_in_systems[part]->current_local_solution.get());
            d_P_in_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
                d_fe_data_managers[part]->buildGhostedSolutionVector(PRESSURE_IN_SYSTEM_NAME, /*localize_data*/ false));

            d_P_out_systems[part] = &d_equation_systems[part]->get_system(PRESSURE_OUT_SYSTEM_NAME);
            d_P_out_half_vecs[part] =
                dynamic_cast<PetscVector<double>*>(d_P_out_systems[part]->current_local_solution.get());
            d_P_out_IB_ghost_vecs[part] =
                dynamic_cast<PetscVector<double>*>(d_fe_data_managers[part]->buildGhostedSolutionVector(
                    PRESSURE_OUT_SYSTEM_NAME, /*localize_data*/ false));
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
            d_WSS_in_systems[part] = &d_equation_systems[part]->get_system(WSS_IN_SYSTEM_NAME);
            d_WSS_in_half_vecs[part] =
                dynamic_cast<PetscVector<double>*>(d_WSS_in_systems[part]->current_local_solution.get());
            d_WSS_in_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
                d_fe_data_managers[part]->buildGhostedSolutionVector(WSS_IN_SYSTEM_NAME, /*localize_data*/ false));

            d_WSS_out_systems[part] = &d_equation_systems[part]->get_system(WSS_OUT_SYSTEM_NAME);
            d_WSS_out_half_vecs[part] =
                dynamic_cast<PetscVector<double>*>(d_WSS_out_systems[part]->current_local_solution.get());
            d_WSS_out_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
                d_fe_data_managers[part]->buildGhostedSolutionVector(WSS_OUT_SYSTEM_NAME, /*localize_data*/ false));
        }
        if (d_use_velocity_jump_conditions && d_use_pressure_jump_conditions)
        {
            d_TAU_in_systems[part] = &d_equation_systems[part]->get_system(TAU_IN_SYSTEM_NAME);
            d_TAU_in_half_vecs[part] =
                dynamic_cast<PetscVector<double>*>(d_TAU_in_systems[part]->current_local_solution.get());
            d_TAU_in_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
                d_fe_data_managers[part]->buildGhostedSolutionVector(TAU_IN_SYSTEM_NAME, /*localize_data*/ false));

            d_TAU_out_systems[part] = &d_equation_systems[part]->get_system(TAU_OUT_SYSTEM_NAME);
            d_TAU_out_half_vecs[part] =
                dynamic_cast<PetscVector<double>*>(d_TAU_out_systems[part]->current_local_solution.get());
            d_TAU_out_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
                d_fe_data_managers[part]->buildGhostedSolutionVector(TAU_OUT_SYSTEM_NAME, /*localize_data*/ false));
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
            *d_P_in_half_vecs[part] = *d_P_in_systems[part]->solution;
            *d_P_out_half_vecs[part] = *d_P_out_systems[part]->solution;
        }

        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                *d_DU_jump_half_vecs[part][d] = *d_DU_jump_systems[part][d]->solution;
            }

            *d_WSS_in_half_vecs[part] = *d_WSS_in_systems[part]->solution;
            *d_WSS_out_half_vecs[part] = *d_WSS_out_systems[part]->solution;
        }
        if (d_use_velocity_jump_conditions && d_use_pressure_jump_conditions)
        {
            *d_TAU_in_half_vecs[part] = *d_TAU_in_systems[part]->solution;
            *d_TAU_out_half_vecs[part] = *d_TAU_out_systems[part]->solution;
        }
    }
    return;
} // preprocessIntegrateData

void
IIMethod::postprocessIntegrateData(double /*current_time*/, double /*new_time*/, int /*num_cycles*/)
{
    std::vector<std::vector<libMesh::PetscVector<double>*> > vec_collection_update = {
        d_U_new_vecs, d_X_new_vecs, d_U_n_new_vecs, d_U_t_new_vecs, d_F_half_vecs
    };

    if (d_use_pressure_jump_conditions)
    {
        vec_collection_update.push_back(d_P_jump_half_vecs);
        vec_collection_update.push_back(d_P_in_half_vecs);
        vec_collection_update.push_back(d_P_out_half_vecs);
    }
    if (d_use_velocity_jump_conditions)
    {
        vec_collection_update.push_back(d_WSS_in_half_vecs);
        vec_collection_update.push_back(d_WSS_out_half_vecs);
    }
    if (d_compute_fluid_traction)
    {
        vec_collection_update.push_back(d_TAU_in_half_vecs);
        vec_collection_update.push_back(d_TAU_out_half_vecs);
    }
    batch_vec_ghost_update(vec_collection_update, INSERT_VALUES, SCATTER_FORWARD);

    if (d_compute_fluid_traction)
    {
        // Evaluate the fluid forces on the interface.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int p_data_idx = var_db->mapVariableAndContextToIndex(getINSHierarchyIntegrator()->getPressureVariable(),
                                                                    getINSHierarchyIntegrator()->getScratchContext());
        calculateInterfacialFluidForces(p_data_idx, d_new_time); // TODO: Should this be half_time?
    }

    // Update the coordinate mapping dX = X - s.
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
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
            *d_P_in_systems[part]->solution = *d_P_in_half_vecs[part];
            *d_P_in_systems[part]->current_local_solution = *d_P_in_half_vecs[part];
            *d_P_out_systems[part]->solution = *d_P_out_half_vecs[part];
            *d_P_out_systems[part]->current_local_solution = *d_P_out_half_vecs[part];
        }

        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                *d_DU_jump_systems[part][d]->solution = *d_DU_jump_half_vecs[part][d];
                *d_DU_jump_systems[part][d]->current_local_solution = *d_DU_jump_half_vecs[part][d];
            }
            *d_WSS_in_systems[part]->solution = *d_WSS_in_half_vecs[part];
            *d_WSS_in_systems[part]->current_local_solution = *d_WSS_in_half_vecs[part];
            *d_WSS_out_systems[part]->solution = *d_WSS_out_half_vecs[part];
            *d_WSS_out_systems[part]->current_local_solution = *d_WSS_out_half_vecs[part];
        }
        if (d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
        {
            *d_TAU_in_systems[part]->solution = *d_TAU_in_half_vecs[part];
            *d_TAU_in_systems[part]->current_local_solution = *d_TAU_in_half_vecs[part];
            *d_TAU_out_systems[part]->solution = *d_TAU_out_half_vecs[part];
            *d_TAU_out_systems[part]->current_local_solution = *d_TAU_out_half_vecs[part];
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

    d_WSS_in_systems.clear();
    d_WSS_in_half_vecs.clear();
    d_WSS_in_IB_ghost_vecs.clear();

    d_WSS_out_systems.clear();
    d_WSS_out_half_vecs.clear();
    d_WSS_out_IB_ghost_vecs.clear();

    d_P_in_systems.clear();
    d_P_in_half_vecs.clear();
    d_P_in_IB_ghost_vecs.clear();

    d_P_out_systems.clear();
    d_P_out_half_vecs.clear();
    d_P_out_IB_ghost_vecs.clear();

    d_TAU_in_systems.clear();
    d_TAU_in_half_vecs.clear();
    d_TAU_in_IB_ghost_vecs.clear();

    d_TAU_out_systems.clear();
    d_TAU_out_half_vecs.clear();
    d_TAU_out_IB_ghost_vecs.clear();

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_half_time = std::numeric_limits<double>::quiet_NaN();
    return;
} // postprocessIntegrateData

void
IIMethod::interpolateVelocity(const int u_data_idx,
                              const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                              const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                              const double data_time)
{
    const double mu = getINSHierarchyIntegrator()->getStokesSpecifications()->getMu();
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const int coarsest_ln = 0;
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        if (ln < static_cast<int>(u_synch_scheds.size()) && u_synch_scheds[ln])
        {
            u_synch_scheds[ln]->coarsenData();
        }
    }

    // Communicate ghost data.
    for (const auto& u_ghost_fill_sched : u_ghost_fill_scheds)
    {
        if (u_ghost_fill_sched) u_ghost_fill_sched->fillData(data_time);
    }

    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        NumericVector<double>* U_vec = nullptr;
        NumericVector<double>* U_n_vec = nullptr;
        NumericVector<double>* U_t_vec = nullptr;
        NumericVector<double>* X_vec = nullptr;
        NumericVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
        const std::array<PetscVector<double>*, NDIM> DU_jump_ghost_vec = {
            d_use_velocity_jump_conditions ? d_DU_jump_IB_ghost_vecs[part][0] : nullptr,
            d_use_velocity_jump_conditions ? d_DU_jump_IB_ghost_vecs[part][1] : nullptr,
#if (NDIM == 3)
            d_use_velocity_jump_conditions ? d_DU_jump_IB_ghost_vecs[part][2] : nullptr,
#endif
        };
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

        NumericVector<double>* WSS_in_vec = d_use_velocity_jump_conditions ? d_WSS_in_half_vecs[part] : nullptr;
        NumericVector<double>* WSS_out_vec = d_use_velocity_jump_conditions ? d_WSS_out_half_vecs[part] : nullptr;

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
        System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
        const DofMap& X_dof_map = X_system.get_dof_map();
        FEDataManager::SystemDofMapCache& X_dof_map_cache =
            *d_fe_data_managers[part]->getDofMapCache(COORDS_SYSTEM_NAME);
        std::vector<std::vector<unsigned int> > U_dof_indices(NDIM);
        std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
        FEType U_fe_type = U_dof_map.variable_type(0);
        for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(U_dof_map.variable_type(d) == U_fe_type);
        FEType X_fe_type = X_dof_map.variable_type(0);
        for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
        TBOX_ASSERT(U_fe_type == X_fe_type);

        std::array<System*, NDIM> DU_jump_system;
        std::array<const DofMap*, NDIM> DU_jump_dof_map;
        std::array<FEDataManager::SystemDofMapCache*, NDIM> DU_jump_dof_map_cache;
        std::array<std::vector<std::vector<unsigned int> >, NDIM> DU_jump_dof_indices;
        FEType DU_jump_fe_type;
        std::vector<std::vector<unsigned int> > WSS_out_dof_indices(NDIM);
        System* WSS_out_system;
        const DofMap* WSS_out_dof_map = NULL;
        FEDataManager::SystemDofMapCache* WSS_out_dof_map_cache = NULL;

        std::vector<std::vector<unsigned int> > WSS_in_dof_indices(NDIM);
        System* WSS_in_system;
        const DofMap* WSS_in_dof_map = NULL;
        FEDataManager::SystemDofMapCache* WSS_in_dof_map_cache = NULL;

        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                DU_jump_system[i] = &equation_systems->get_system(VELOCITY_JUMP_SYSTEM_NAME[i]);
                DU_jump_dof_map[i] = &DU_jump_system[i]->get_dof_map();
                DU_jump_dof_map_cache[i] = d_fe_data_managers[part]->getDofMapCache(VELOCITY_JUMP_SYSTEM_NAME[i]);
                DU_jump_fe_type = DU_jump_dof_map[i]->variable_type(0);
                for (unsigned int j = 0; j < NDIM; ++j)
                {
                    TBOX_ASSERT(DU_jump_dof_map[i]->variable_type(j) == DU_jump_fe_type);
                }
                DU_jump_dof_indices[i].resize(NDIM);
            }

            WSS_out_system = &equation_systems->get_system(WSS_OUT_SYSTEM_NAME);
            WSS_out_dof_map = &WSS_out_system->get_dof_map();
            WSS_out_dof_map_cache = d_fe_data_managers[part]->getDofMapCache(WSS_OUT_SYSTEM_NAME);
            FEType WSS_out_fe_type = WSS_out_dof_map->variable_type(0);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                TBOX_ASSERT(WSS_out_dof_map->variable_type(d) == WSS_out_fe_type);
            }

            WSS_in_system = &equation_systems->get_system(WSS_IN_SYSTEM_NAME);
            WSS_in_dof_map = &WSS_in_system->get_dof_map();
            WSS_in_dof_map_cache = d_fe_data_managers[part]->getDofMapCache(WSS_IN_SYSTEM_NAME);
            FEType WSS_in_fe_type = WSS_in_dof_map->variable_type(0);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                TBOX_ASSERT(WSS_in_dof_map->variable_type(d) == WSS_in_fe_type);
            }
        }
        FEType fe_type = U_fe_type;
        std::unique_ptr<FEBase> fe_X = FEBase::build(dim, fe_type);
        const std::vector<double>& JxW = fe_X->get_JxW();
        const std::vector<std::vector<double> >& phi_X = fe_X->get_phi();
        std::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi;
        dphi_dxi[0] = &fe_X->get_dphidxi();
        if (NDIM > 2) dphi_dxi[1] = &fe_X->get_dphideta();

        FEType fe_DU_jump_type = DU_jump_fe_type;
        std::unique_ptr<FEBase> fe_DU_jump = FEBase::build(dim, fe_DU_jump_type);
        const std::vector<double>& JxW_jump = fe_DU_jump->get_JxW();
        const std::vector<std::vector<double> >& phi_DU_jump = fe_DU_jump->get_phi();

        X_ghost_vec->close();

        // Loop over the patches to interpolate values to the element quadrature
        // points from the grid, then use these values to compute the projection
        // of the interpolated velocity field onto the FE basis functions.
        std::unique_ptr<NumericVector<double> > U_rhs_vec = U_vec->zero_clone();
        std::vector<DenseVector<double> > U_rhs_e(NDIM);
        std::unique_ptr<NumericVector<double> > U_n_rhs_vec = U_n_vec->zero_clone();
        std::vector<DenseVector<double> > U_n_rhs_e(NDIM);
        std::unique_ptr<NumericVector<double> > U_t_rhs_vec = U_t_vec->zero_clone();
        std::vector<DenseVector<double> > U_t_rhs_e(NDIM);

        std::unique_ptr<NumericVector<double> > WSS_out_rhs_vec =
            (d_use_velocity_jump_conditions ? WSS_out_vec->zero_clone() : std::unique_ptr<NumericVector<double> >());
        DenseVector<double> WSS_out_rhs_e[NDIM];

        std::unique_ptr<NumericVector<double> > WSS_in_rhs_vec =
            (d_use_velocity_jump_conditions ? WSS_in_vec->zero_clone() : std::unique_ptr<NumericVector<double> >());
        DenseVector<double> WSS_in_rhs_e[NDIM];

        boost::multi_array<double, 2> x_node;
        std::array<boost::multi_array<double, 2>, NDIM> DU_jump_node;
        std::vector<double> U_qp, U_in_qp, U_out_qp, WSS_in_qp, WSS_out_qp, n_qp, x_qp, x_in_qp, x_out_qp;
        std::array<std::vector<double>, NDIM> DU_jump_qp;
        VectorValue<double> U, WSS_in, WSS_out, U_n, U_t, n;
        std::array<VectorValue<double>, 2> dx_dxi;

        Pointer<PatchLevel<NDIM> > level =
            d_hierarchy->getPatchLevel(d_fe_data_managers[part]->getFinestPatchLevelNumber());
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
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    X_dof_map_cache.dof_indices(elem, X_dof_indices[axis], axis);
                }
                get_values_for_interpolation(x_node, *X_ghost_vec, X_dof_indices);
                FEDataManager::updateInterpQuadratureRule(qrule, d_default_interp_spec, elem, x_node, patch_dx_min);
                n_qpoints_patch += qrule->n_points();
            }

            if (!n_qpoints_patch) continue;
            U_qp.resize(NDIM * n_qpoints_patch);
            WSS_in_qp.resize(NDIM * n_qpoints_patch);
            WSS_out_qp.resize(NDIM * n_qpoints_patch);
            U_in_qp.resize(NDIM * n_qpoints_patch);
            U_out_qp.resize(NDIM * n_qpoints_patch);
            x_qp.resize(NDIM * n_qpoints_patch);
            x_in_qp.resize(NDIM * n_qpoints_patch);
            x_out_qp.resize(NDIM * n_qpoints_patch);
            n_qp.resize(NDIM * n_qpoints_patch);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                DU_jump_qp[axis].resize(NDIM * n_qpoints_patch);
            }
            std::fill(U_qp.begin(), U_qp.end(), 0.0);
            std::fill(U_in_qp.begin(), U_in_qp.end(), 0.0);
            std::fill(U_out_qp.begin(), U_out_qp.end(), 0.0);
            std::fill(WSS_in_qp.begin(), WSS_in_qp.end(), 0.0);
            std::fill(WSS_out_qp.begin(), WSS_out_qp.end(), 0.0);

            // Loop over the elements and compute the positions of the quadrature points.
            qrule.reset();
            unsigned int qp_offset = 0;
            for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems[e_idx];
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X_dof_map_cache.dof_indices(elem, X_dof_indices[d], d);
                }
                get_values_for_interpolation(x_node, *X_ghost_vec, X_dof_indices);
                if (d_use_velocity_jump_conditions)
                {
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            DU_jump_dof_map_cache[axis]->dof_indices(elem, DU_jump_dof_indices[axis][d], d);
                        }
                        get_values_for_interpolation(
                            DU_jump_node[axis], *DU_jump_ghost_vec[axis], DU_jump_dof_indices[axis]);
                    }
                }
                const bool qrule_changed =
                    FEDataManager::updateInterpQuadratureRule(qrule, d_default_interp_spec, elem, x_node, patch_dx_min);
                if (qrule_changed) fe_X->attach_quadrature_rule(qrule.get());
                fe_X->reinit(elem);
                if (d_use_velocity_jump_conditions)
                {
                    if (qrule_changed) fe_DU_jump->attach_quadrature_rule(qrule.get());
                    fe_DU_jump->reinit(elem);
                }
                const unsigned int n_nodes = elem->n_nodes();
                const unsigned int n_qpoints = qrule->n_points();

                // Zero out the values prior to accumulation.
                double* x_begin = &x_qp[NDIM * qp_offset];
                std::fill(x_begin, x_begin + NDIM * n_qpoints, 0.0);

                double* x_in_begin = &x_in_qp[NDIM * qp_offset];
                std::fill(x_in_begin, x_in_begin + NDIM * n_qpoints, 0.0);

                double* x_out_begin = &x_out_qp[NDIM * qp_offset];
                std::fill(x_out_begin, x_out_begin + NDIM * n_qpoints, 0.0);

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
                for (unsigned int k = 0; k < n_nodes; ++k)
                {
                    for (unsigned int qp = 0; qp < n_qpoints; ++qp)
                    {
                        const double& p = phi_X[k][qp];
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            x_qp[NDIM * (qp_offset + qp) + d] += x_node[k][d] * p;
                        }
                        if (d_use_velocity_jump_conditions)
                        {
                            const double& p2 = phi_DU_jump[k][qp];
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
                        interpolate(dx_dxi[l], qp, x_node, *dphi_dxi[l]);
                    }
                    if (NDIM == 2)
                    {
                        dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                    }
                    n = (dx_dxi[0].cross(dx_dxi[1])).unit();
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        n_qp[NDIM * (qp_offset + qp) + d] = n(d);
                        x_in_qp[NDIM * (qp_offset + qp) + d] = x_qp[NDIM * (qp_offset + qp) + d] - n(d) * dh;
                        x_out_qp[NDIM * (qp_offset + qp) + d] = x_qp[NDIM * (qp_offset + qp) + d] + n(d) * dh;
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
                    U_in_qp, NDIM, x_in_qp, NDIM, u_sc_data, patch, ghost_box, d_default_interp_spec.kernel_fcn);

                LEInteractor::interpolate(
                    U_out_qp, NDIM, x_out_qp, NDIM, u_sc_data, patch, ghost_box, d_default_interp_spec.kernel_fcn);

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

                    const double* const x_in = &x_in_qp[NDIM * k];
                    const Index<NDIM> in = IndexUtilities::getCellIndex(
                        x_in, x_lower_ghost, x_upper_ghost, patch_geom->getDx(), ghost_box.lower(), ghost_box.upper());

                    const double* const x_out = &x_out_qp[NDIM * k];
                    const Index<NDIM> out = IndexUtilities::getCellIndex(
                        x_out, x_lower_ghost, x_upper_ghost, patch_geom->getDx(), ghost_box.lower(), ghost_box.upper());

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
                                                                (coeff_vec * du_jump);
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
                            if (d_use_velocity_jump_conditions)
                            {
                                const double CC = (nproj > 0.0) ? Ujump[ic[0]][ic[1]][axis] : 0.0;
                                U_axis[s] -= CC / mu;
                            }
#endif
#if (NDIM == 3)

                            U_axis[s] += w[0][ic[0] - ic_lower[0]] * w[1][ic[1] - ic_lower[1]] *
                                         w[2][ic[2] - ic_lower[2]] * u_sc_data_array[ic[0]][ic[1]][ic[2]];
                            const double nproj = n_qp[s * NDIM + 0] * wr[0][ic_upper[0] - ic[0]] +
                                                 n_qp[s * NDIM + 1] * wr[1][ic_upper[1] - ic[1]] +
                                                 n_qp[s * NDIM + 2] * wr[2][ic_upper[2] - ic[2]];
                            if (d_use_velocity_jump_conditions)
                            {
                                const double CC = (nproj > 0.0) ? Ujump[ic[0]][ic[1]][ic[2]][axis] : 0.0;
                                U_axis[s] -= CC / mu;
                            }
#endif
                        }
                    }
                    if (d_use_velocity_jump_conditions)
                    {
                        for (unsigned int k = 0; k < local_indices.size(); ++k)
                        {
                            U_qp[NDIM * local_indices[k] + axis] = U_axis[local_indices[k]];
                            if (dh != 0.0)
                            {
                                WSS_in_qp[NDIM * local_indices[k] + axis] =
                                    mu * (1.0 / dh) *
                                    (U_in_qp[NDIM * local_indices[k] + axis] - U_qp[NDIM * local_indices[k] + axis]);

                                double du_dn_jump = 0.0;
                                for (int dd = 0; dd < NDIM; ++dd)
                                {
                                    du_dn_jump += DU_jump_qp[axis][NDIM * local_indices[k] + dd] *
                                                  n_qp[NDIM * local_indices[k] + dd];
                                }
                                WSS_out_qp[NDIM * local_indices[k] + axis] =
                                    (1.0 - d_exterior_calc_coef) *
                                        (du_dn_jump - WSS_in_qp[NDIM * local_indices[k] + axis]) +
                                    d_exterior_calc_coef * mu * (1.0 / dh) *
                                        (U_out_qp[NDIM * local_indices[k] + axis] -
                                         U_qp[NDIM * local_indices[k] + axis]);
                            }
                            else
                            {
                                TBOX_ERROR(d_object_name << ": The width for the wall shear stress hasn't been set up!"
                                                         << std::endl);
                            }
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
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    U_dof_map_cache.dof_indices(elem, U_dof_indices[d], d);
                    U_rhs_e[d].resize(static_cast<int>(U_dof_indices[d].size()));
                    U_n_rhs_e[d].resize(static_cast<int>(U_dof_indices[d].size()));
                    U_t_rhs_e[d].resize(static_cast<int>(U_dof_indices[d].size()));
                    X_dof_map_cache.dof_indices(elem, X_dof_indices[d], d);
                    if (d_use_velocity_jump_conditions)
                    {
                        WSS_out_dof_map_cache->dof_indices(elem, WSS_out_dof_indices[d], d);
                        WSS_out_rhs_e[d].resize(static_cast<int>(WSS_out_dof_indices[d].size()));

                        WSS_in_dof_map_cache->dof_indices(elem, WSS_in_dof_indices[d], d);
                        WSS_in_rhs_e[d].resize(static_cast<int>(WSS_in_dof_indices[d].size()));
                    }
                }
                get_values_for_interpolation(x_node, *X_ghost_vec, X_dof_indices);
                const bool qrule_changed =
                    FEDataManager::updateInterpQuadratureRule(qrule, d_default_interp_spec, elem, x_node, patch_dx_min);
                if (qrule_changed) fe_X->attach_quadrature_rule(qrule.get());

                fe_X->reinit(elem);
                if (d_use_velocity_jump_conditions)
                {
                    if (qrule_changed)
                    {
                        fe_DU_jump->attach_quadrature_rule(qrule.get());
                    }
                    fe_DU_jump->reinit(elem);
                }
                const unsigned int n_qpoints = qrule->n_points();
                const size_t n_basis = U_dof_indices[0].size();
                const size_t n_basis_jump = WSS_in_dof_indices[0].size();

                for (unsigned int qp = 0; qp < n_qpoints; ++qp)
                {
                    for (unsigned int k = 0; k < NDIM - 1; ++k)
                    {
                        interpolate(dx_dxi[k], qp, x_node, *dphi_dxi[k]);
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
                            WSS_in(d) = WSS_in_qp[NDIM * (qp_offset + qp) + d];
                            WSS_out(d) = WSS_out_qp[NDIM * (qp_offset + qp) + d];
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
                            U_n_rhs_e[d](k) += U_n(d) * p_JxW;
                            U_t_rhs_e[d](k) += U_t(d) * p_JxW;
                        }
                    }
                    if (d_use_velocity_jump_conditions)
                    {
                        for (unsigned int k = 0; k < n_basis_jump; ++k)
                        {
                            const double p_JxW_jump = phi_DU_jump[k][qp] * JxW_jump[qp];
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                WSS_in_rhs_e[d](k) += WSS_in(d) * p_JxW_jump;
                                WSS_out_rhs_e[d](k) += WSS_out(d) * p_JxW_jump;
                            }
                        }
                    }
                }
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    U_dof_map.constrain_element_vector(U_rhs_e[d], U_dof_indices[d]);
                    U_dof_map.constrain_element_vector(U_n_rhs_e[d], U_dof_indices[d]);
                    U_dof_map.constrain_element_vector(U_t_rhs_e[d], U_dof_indices[d]);
                    U_rhs_vec->add_vector(U_rhs_e[d], U_dof_indices[d]);
                    U_n_rhs_vec->add_vector(U_n_rhs_e[d], U_dof_indices[d]);
                    U_t_rhs_vec->add_vector(U_t_rhs_e[d], U_dof_indices[d]);
                    if (d_use_velocity_jump_conditions)
                    {
                        WSS_out_dof_map->constrain_element_vector(WSS_out_rhs_e[d], WSS_out_dof_indices[d]);
                        WSS_out_rhs_vec->add_vector(WSS_out_rhs_e[d], WSS_out_dof_indices[d]);
                        WSS_in_dof_map->constrain_element_vector(WSS_in_rhs_e[d], WSS_in_dof_indices[d]);
                        WSS_in_rhs_vec->add_vector(WSS_in_rhs_e[d], WSS_in_dof_indices[d]);
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
            WSS_in_rhs_vec->close();
            d_fe_data_managers[part]->computeL2Projection(
                *WSS_in_vec, *WSS_in_rhs_vec, WSS_IN_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
            WSS_in_vec->close();

            WSS_out_rhs_vec->close();
            d_fe_data_managers[part]->computeL2Projection(
                *WSS_out_vec, *WSS_out_rhs_vec, WSS_OUT_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
            WSS_out_vec->close();
        }

        // Solve for the nodal values.
        d_fe_data_managers[part]->computeL2Projection(
            *U_vec, *U_rhs_vec, VELOCITY_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
        U_vec->close();
        d_fe_data_managers[part]->computeL2Projection(
            *U_n_vec, *U_n_rhs_vec, VELOCITY_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
        U_n_vec->close();
        d_fe_data_managers[part]->computeL2Projection(
            *U_t_vec, *U_t_rhs_vec, VELOCITY_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
        U_t_vec->close();
        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int d = 0; d < NDIM; ++d) d_DU_jump_IB_ghost_vecs[part][d]->close();
        }
        d_X_IB_ghost_vecs[part]->close();
    }
    return;
} // interpolateVelocity

void
IIMethod::computeFluidTraction(const double data_time, unsigned int part)
{
    batch_vec_ghost_update({ d_WSS_in_half_vecs[part],
                             d_WSS_out_half_vecs[part],
                             d_P_in_half_vecs[part],
                             d_P_out_half_vecs[part],
                             d_TAU_in_half_vecs[part],
                             d_TAU_out_half_vecs[part],
                             d_X_new_vecs[part] },
                           INSERT_VALUES,
                           SCATTER_FORWARD);
    NumericVector<double>* WSS_in_vec = NULL;
    NumericVector<double>* WSS_in_ghost_vec = d_WSS_in_IB_ghost_vecs[part];

    NumericVector<double>* WSS_out_vec = NULL;
    NumericVector<double>* WSS_out_ghost_vec = d_WSS_out_IB_ghost_vecs[part];

    NumericVector<double>* P_in_vec = NULL;
    NumericVector<double>* P_in_ghost_vec = d_P_in_IB_ghost_vecs[part];

    NumericVector<double>* P_out_vec = NULL;
    NumericVector<double>* P_out_ghost_vec = d_P_out_IB_ghost_vecs[part];

    NumericVector<double>* TAU_in_vec = d_TAU_in_half_vecs[part];
    NumericVector<double>* TAU_out_vec = d_TAU_out_half_vecs[part];

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

    WSS_in_vec = d_WSS_in_half_vecs[part];
    copy_and_synch(*WSS_in_vec, *WSS_in_ghost_vec);

    WSS_out_vec = d_WSS_out_half_vecs[part];
    copy_and_synch(*WSS_out_vec, *WSS_out_ghost_vec);

    P_in_vec = d_P_in_half_vecs[part];
    copy_and_synch(*P_in_vec, *P_in_ghost_vec);

    P_out_vec = d_P_out_half_vecs[part];
    copy_and_synch(*P_out_vec, *P_out_ghost_vec);

    std::unique_ptr<NumericVector<double> > TAU_in_rhs_vec = TAU_in_vec->zero_clone();
    std::array<DenseVector<double>, NDIM> TAU_in_rhs_e;

    std::unique_ptr<NumericVector<double> > TAU_out_rhs_vec = TAU_out_vec->zero_clone();
    std::array<DenseVector<double>, NDIM> TAU_out_rhs_e;

    // Extract the FE systems and DOF maps, and setup the FE objects.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    std::unique_ptr<QBase> qrule;

    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    FEDataManager::SystemDofMapCache& X_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(COORDS_SYSTEM_NAME);
    FEType X_fe_type = X_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
    }

    System& P_in_system = equation_systems->get_system(PRESSURE_IN_SYSTEM_NAME);
    const DofMap& P_in_dof_map = P_in_system.get_dof_map();
    FEDataManager::SystemDofMapCache& P_in_dof_map_cache =
        *d_fe_data_managers[part]->getDofMapCache(PRESSURE_IN_SYSTEM_NAME);
    FEType P_in_fe_type = P_in_dof_map.variable_type(0);

    System& P_out_system = equation_systems->get_system(PRESSURE_OUT_SYSTEM_NAME);
    const DofMap& P_out_dof_map = P_out_system.get_dof_map();
    FEDataManager::SystemDofMapCache& P_out_dof_map_cache =
        *d_fe_data_managers[part]->getDofMapCache(PRESSURE_OUT_SYSTEM_NAME);
    FEType P_out_fe_type = P_out_dof_map.variable_type(0);

    System& WSS_in_system = equation_systems->get_system(WSS_IN_SYSTEM_NAME);
    const DofMap& WSS_in_dof_map = WSS_in_system.get_dof_map();
    FEDataManager::SystemDofMapCache& WSS_in_dof_map_cache =
        *d_fe_data_managers[part]->getDofMapCache(WSS_IN_SYSTEM_NAME);
    FEType WSS_in_fe_type = WSS_in_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(WSS_in_dof_map.variable_type(d) == WSS_in_fe_type);
    }

    System& WSS_out_system = equation_systems->get_system(WSS_OUT_SYSTEM_NAME);
    const DofMap& WSS_out_dof_map = WSS_out_system.get_dof_map();
    FEDataManager::SystemDofMapCache& WSS_out_dof_map_cache =
        *d_fe_data_managers[part]->getDofMapCache(WSS_OUT_SYSTEM_NAME);
    FEType WSS_out_fe_type = WSS_out_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(WSS_out_dof_map.variable_type(d) == WSS_out_fe_type);
    }

    System& TAU_in_system = equation_systems->get_system(TAU_IN_SYSTEM_NAME);
    const DofMap& TAU_in_dof_map = TAU_in_system.get_dof_map();
    FEDataManager::SystemDofMapCache& TAU_in_dof_map_cache =
        *d_fe_data_managers[part]->getDofMapCache(TAU_IN_SYSTEM_NAME);
    FEType TAU_in_fe_type = TAU_in_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(TAU_in_dof_map.variable_type(d) == TAU_in_fe_type);
    }

    System& TAU_out_system = equation_systems->get_system(TAU_OUT_SYSTEM_NAME);
    const DofMap& TAU_out_dof_map = TAU_out_system.get_dof_map();
    FEDataManager::SystemDofMapCache& TAU_out_dof_map_cache =
        *d_fe_data_managers[part]->getDofMapCache(TAU_OUT_SYSTEM_NAME);
    FEType TAU_out_fe_type = TAU_out_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(TAU_out_dof_map.variable_type(d) == TAU_out_fe_type);
    }

    TBOX_ASSERT(P_in_fe_type == TAU_in_fe_type);
    TBOX_ASSERT(P_out_fe_type == P_in_fe_type);

    std::unique_ptr<FEBase> fe_X = FEBase::build(dim, X_fe_type);
    const std::vector<double>& JxW = fe_X->get_JxW();
    const std::vector<std::vector<double> >& phi_X = fe_X->get_phi();
    std::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi_X;
    dphi_dxi_X[0] = &fe_X->get_dphidxi();
    if (NDIM > 2) dphi_dxi_X[1] = &fe_X->get_dphideta();

    std::unique_ptr<FEBase> fe_P = FEBase::build(dim, P_out_fe_type);
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

    boost::multi_array<double, 2> x_node, X_node, WSS_in_node, WSS_out_node, n_qp_node;
    boost::multi_array<double, 1> P_in_node, P_out_node;
    std::vector<double> x_in_qp, x_out_qp, x_qp;
    std::vector<double> P_in_qp, P_out_qp, Normal_qp, WSS_in_qp, WSS_out_qp, TAU_in_qp, TAU_out_qp;
    std::array<VectorValue<double>, 2> dX_dxi, dx_dxi;
    VectorValue<double> n, N, x, X;

    Pointer<PatchLevel<NDIM> > level =
        d_hierarchy->getPatchLevel(d_fe_data_managers[part]->getFinestPatchLevelNumber());
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
        P_in_qp.resize(n_qp_patch);
        P_out_qp.resize(n_qp_patch);
        x_qp.resize(NDIM * n_qp_patch);
        WSS_in_qp.resize(NDIM * n_qp_patch);
        WSS_out_qp.resize(NDIM * n_qp_patch);

        TAU_in_qp.resize(NDIM * n_qp_patch);
        TAU_out_qp.resize(NDIM * n_qp_patch);
        Normal_qp.resize(NDIM * n_qp_patch);
        std::fill(Normal_qp.begin(), Normal_qp.end(), 0.0);
        std::fill(x_qp.begin(), x_qp.end(), 0.0);
        std::fill(WSS_in_qp.begin(), WSS_in_qp.end(), 0.0);
        std::fill(WSS_out_qp.begin(), WSS_out_qp.end(), 0.0);
        std::fill(P_in_qp.begin(), P_in_qp.end(), 0.0);
        std::fill(P_out_qp.begin(), P_out_qp.end(), 0.0);
        std::fill(TAU_in_qp.begin(), TAU_in_qp.end(), 0.0);
        std::fill(TAU_out_qp.begin(), TAU_out_qp.end(), 0.0);

        // Loop over the elements and compute the positions of the quadrature points.
        qrule.reset();
        unsigned int qp_offset = 0;
        std::vector<libMesh::dof_id_type> dof_id_scratch;
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
            const auto& WSS_in_dof_indices = WSS_in_dof_map_cache.dof_indices(elem);
            const auto& WSS_out_dof_indices = WSS_out_dof_map_cache.dof_indices(elem);
            const auto& P_in_dof_indices = P_in_dof_map_cache.dof_indices(elem);
            const auto& P_out_dof_indices = P_out_dof_map_cache.dof_indices(elem);
            get_values_for_interpolation(x_node, *X_petsc_vec, X_local_soln, X_dof_indices);
            get_values_for_interpolation(WSS_in_node, *WSS_in_ghost_vec, WSS_in_dof_indices);
            get_values_for_interpolation(WSS_out_node, *WSS_out_ghost_vec, WSS_out_dof_indices);
            copy_dof_ids_to_vector(0, P_in_dof_indices, dof_id_scratch);
            get_values_for_interpolation(P_in_node, *P_in_ghost_vec, dof_id_scratch);
            copy_dof_ids_to_vector(0, P_out_dof_indices, dof_id_scratch);
            get_values_for_interpolation(P_out_node, *P_out_ghost_vec, dof_id_scratch);
            get_values_for_interpolation(X_node, *X0_vec, X_dof_indices);

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

            double* Normal_begin = &Normal_qp[NDIM * qp_offset];
            std::fill(Normal_begin, Normal_begin + NDIM * n_qp, 0.0);

            double* WSS_in_begin = &WSS_in_qp[NDIM * qp_offset];
            std::fill(WSS_in_begin, WSS_in_begin + NDIM * n_qp, 0.0);

            double* WSS_out_begin = &WSS_out_qp[NDIM * qp_offset];
            std::fill(WSS_out_begin, WSS_out_begin + NDIM * n_qp, 0.0);

            double* TAU_in_begin = &TAU_in_qp[NDIM * qp_offset];
            std::fill(TAU_in_begin, TAU_in_begin + NDIM * n_qp, 0.0);

            double* TAU_out_begin = &TAU_out_qp[NDIM * qp_offset];
            std::fill(TAU_out_begin, TAU_out_begin + NDIM * n_qp, 0.0);

            double* P_in_begin = &P_in_qp[qp_offset];
            std::fill(P_in_begin, P_in_begin + n_qp, 0.0);

            double* P_out_begin = &P_out_qp[qp_offset];
            std::fill(P_out_begin, P_out_begin + n_qp, 0.0);

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
                        WSS_in_qp[NDIM * (qp_offset + qp) + i] += (da / dA) * WSS_in_node[k][i] * p_P;
                        WSS_out_qp[NDIM * (qp_offset + qp) + i] += (da / dA) * WSS_out_node[k][i] * p_P;
                    }
                    Normal_qp[NDIM * (qp_offset + qp) + i] = n(i);
                }

                for (unsigned int k = 0; k < n_node; ++k)
                {
                    const double& p_P = phi_P[k][qp];
                    P_in_qp[qp_offset + qp] += (da / dA) * P_in_node[k] * p_P;
                    P_out_qp[qp_offset + qp] += (da / dA) * P_out_node[k] * p_P;
                }
            }
            qp_offset += n_qp;
        }

        const Box<NDIM>& interp_box = patch->getBox();
        std::vector<int> local_indices;
        local_indices.clear();
        const int upper_bound = n_qp_patch;
        if (upper_bound == 0) return;

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
                    // calculate both the interior and exterior fluid tracitons (tau)
                    TAU_in_qp[NDIM * local_indices[k] + axis] =
                        WSS_in_qp[NDIM * local_indices[k] + axis] -
                        P_in_qp[local_indices[k]] * Normal_qp[NDIM * local_indices[k] + axis];

                    TAU_out_qp[NDIM * local_indices[k] + axis] =
                        WSS_out_qp[NDIM * local_indices[k] + axis] -
                        P_out_qp[local_indices[k]] * Normal_qp[NDIM * local_indices[k] + axis];
                }
            }
        }

        // Loop over the elements and accumulate the right-hand-side values.
        qrule.reset();
        qp_offset = 0;
        std::vector<libMesh::dof_id_type> dof_id_in_scratch;
        std::vector<libMesh::dof_id_type> dof_id_out_scratch;
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
            const auto& TAU_in_dof_indices = TAU_in_dof_map_cache.dof_indices(elem);
            const auto& TAU_out_dof_indices = TAU_out_dof_map_cache.dof_indices(elem);

            for (unsigned int i = 0; i < NDIM; ++i)
            {
                TAU_in_rhs_e[i].resize(static_cast<int>(TAU_in_dof_indices[i].size()));
                TAU_out_rhs_e[i].resize(static_cast<int>(TAU_out_dof_indices[i].size()));
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
            const size_t n_basis2 = TAU_out_dof_indices[0].size();
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const int idx = NDIM * (qp_offset + qp);
                for (unsigned int k = 0; k < n_basis2; ++k)
                {
                    const double p_JxW = phi_P[k][qp] * JxW[qp];
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        TAU_in_rhs_e[i](k) += TAU_in_qp[idx + i] * p_JxW;
                        TAU_out_rhs_e[i](k) += TAU_out_qp[idx + i] * p_JxW;
                    }
                }
            }
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                copy_dof_ids_to_vector(d, TAU_in_dof_indices, dof_id_in_scratch);
                TAU_in_dof_map.constrain_element_vector(TAU_in_rhs_e[d], dof_id_in_scratch);
                TAU_in_rhs_vec->add_vector(TAU_in_rhs_e[d], dof_id_in_scratch);

                copy_dof_ids_to_vector(d, TAU_out_dof_indices, dof_id_out_scratch);
                TAU_out_dof_map.constrain_element_vector(TAU_out_rhs_e[d], dof_id_out_scratch);
                TAU_out_rhs_vec->add_vector(TAU_out_rhs_e[d], dof_id_out_scratch);
            }
            qp_offset += n_qp;
        }
    }

    TAU_in_rhs_vec->close();
    TAU_out_rhs_vec->close();

    d_fe_data_managers[part]->computeL2Projection(
        *TAU_in_vec, *TAU_in_rhs_vec, TAU_IN_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
    TAU_in_vec->close();

    d_fe_data_managers[part]->computeL2Projection(
        *TAU_out_vec, *TAU_out_rhs_vec, TAU_OUT_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
    TAU_out_vec->close();

    d_X_half_vecs[part]->close();
    d_X_current_vecs[part]->close();
    d_X_new_vecs[part]->close();
    d_TAU_in_half_vecs[part]->close();
    d_WSS_in_half_vecs[part]->close();
    d_P_in_half_vecs[part]->close();
    d_TAU_out_half_vecs[part]->close();
    d_WSS_out_half_vecs[part]->close();
    d_P_out_half_vecs[part]->close();

    VecRestoreArray(X_local_vec, &X_local_soln);
    VecGhostRestoreLocalForm(X_global_vec, &X_local_vec);

    d_WSS_in_IB_ghost_vecs[part]->close();
    d_WSS_out_IB_ghost_vecs[part]->close();

    d_P_in_IB_ghost_vecs[part]->close();
    d_P_out_IB_ghost_vecs[part]->close();
    d_X_IB_ghost_vecs[part]->close();

    return;
} // computeFluidTraction

void
IIMethod::extrapolatePressureForTraction(const int p_data_idx, const double data_time, unsigned int part)
{
    batch_vec_ghost_update(
        { d_P_out_half_vecs[part], d_P_jump_half_vecs[part], d_P_in_half_vecs[part], d_X_new_vecs[part] },
        INSERT_VALUES,
        SCATTER_FORWARD);

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_fe_data_managers[part]->getPatchHierarchy();

    NumericVector<double>* P_in_vec = d_P_in_half_vecs[part];
    NumericVector<double>* P_out_vec = d_P_out_half_vecs[part];
    NumericVector<double>* P_jump_ghost_vec = d_P_jump_IB_ghost_vecs[part];
    NumericVector<double>* X_vec = NULL;
    NumericVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];

    std::unique_ptr<NumericVector<double> > P_in_rhs_vec = (*P_in_vec).zero_clone();
    P_in_rhs_vec->zero();
    DenseVector<double> P_in_rhs_e;

    std::unique_ptr<NumericVector<double> > P_out_rhs_vec = (*P_out_vec).zero_clone();
    P_out_rhs_vec->zero();
    DenseVector<double> P_out_rhs_e;

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
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
    FEType X_fe_type = X_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
    }

    FEDataManager::SystemDofMapCache& P_jump_dof_map_cache =
        *d_fe_data_managers[part]->getDofMapCache(PRESSURE_JUMP_SYSTEM_NAME);
    std::vector<unsigned int> P_jump_dof_indices;
    std::unique_ptr<FEBase> fe_X = FEBase::build(dim, X_fe_type);
    const std::vector<double>& JxW = fe_X->get_JxW();
    const std::vector<std::vector<double> >& phi_X = fe_X->get_phi();
    std::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi_X;
    dphi_dxi_X[0] = &fe_X->get_dphidxi();
    if (NDIM > 2) dphi_dxi_X[1] = &fe_X->get_dphideta();

    System& P_in_system = equation_systems->get_system(PRESSURE_IN_SYSTEM_NAME);
    const DofMap& P_in_dof_map = P_in_system.get_dof_map();
    FEDataManager::SystemDofMapCache& P_in_dof_map_cache =
        *d_fe_data_managers[part]->getDofMapCache(PRESSURE_IN_SYSTEM_NAME);
    std::vector<unsigned int> P_in_dof_indices;

    System& P_out_system = equation_systems->get_system(PRESSURE_OUT_SYSTEM_NAME);
    const DofMap& P_out_dof_map = P_out_system.get_dof_map();
    FEDataManager::SystemDofMapCache& P_out_dof_map_cache =
        *d_fe_data_managers[part]->getDofMapCache(PRESSURE_OUT_SYSTEM_NAME);
    FEType P_out_fe_type = P_out_dof_map.variable_type(0);
    std::vector<unsigned int> P_out_dof_indices;

    std::unique_ptr<FEBase> fe_P = FEBase::build(dim, P_out_fe_type);
    const std::vector<std::vector<double> >& phi_P = fe_P->get_phi();

    const std::vector<std::vector<Elem*> >& active_patch_element_map =
        d_fe_data_managers[part]->getActivePatchElementMap();

    boost::multi_array<double, 2> x_node;
    std::vector<double> x_qp, x_in_qp, x_out_qp;
    boost::multi_array<double, 1> P_jump_node;
    std::vector<double> P_i_qp, P_o_qp, P_in_qp, P_out_qp, P_jump_qp, N_qp;
    std::array<VectorValue<double>, 2> dx_dxi;

    Pointer<PatchLevel<NDIM> > level =
        d_hierarchy->getPatchLevel(d_fe_data_managers[part]->getFinestPatchLevelNumber());
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

        if (d_p_calc_width == 0)
        {
            TBOX_ERROR(d_object_name << ": The width for the interfacial pressure calc hasn't been set up!"
                                     << std::endl);
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
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(elem, X_dof_indices[d], d);
            }
            get_values_for_interpolation(x_node, *X_ghost_vec, X_dof_indices);

            FEDataManager::updateInterpQuadratureRule(qrule, d_default_interp_spec, elem, x_node, patch_dx_min);

            n_qp_patch += qrule->n_points();
        }

        if (!n_qp_patch) continue;
        P_jump_qp.resize(n_qp_patch);
        P_i_qp.resize(n_qp_patch);
        P_o_qp.resize(n_qp_patch);
        P_in_qp.resize(n_qp_patch);
        P_out_qp.resize(n_qp_patch);
        x_in_qp.resize(NDIM * n_qp_patch);
        x_out_qp.resize(NDIM * n_qp_patch);
        x_qp.resize(NDIM * n_qp_patch);
        N_qp.resize(NDIM * n_qp_patch);
        std::fill(P_i_qp.begin(), P_i_qp.end(), 0.0);
        std::fill(P_o_qp.begin(), P_o_qp.end(), 0.0);
        std::fill(P_in_qp.begin(), P_in_qp.end(), 0.0);
        std::fill(P_out_qp.begin(), P_out_qp.end(), 0.0);
        std::fill(P_jump_qp.begin(), P_jump_qp.end(), 0.0);
        std::fill(N_qp.begin(), N_qp.end(), 0.0);

        // Loop over the elements and compute the positions of the quadrature points.
        qrule.reset();
        unsigned int qp_offset = 0;
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(elem, X_dof_indices[d], d);
            }
            get_values_for_interpolation(x_node, *X_ghost_vec, X_dof_indices);

            P_jump_dof_map_cache.dof_indices(elem, P_jump_dof_indices);
            get_values_for_interpolation(P_jump_node, *P_jump_ghost_vec, P_jump_dof_indices);

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

            double* x_in_begin = &x_in_qp[NDIM * qp_offset];
            std::fill(x_in_begin, x_in_begin + NDIM * n_qp, 0.0);

            double* x_out_begin = &x_out_qp[NDIM * qp_offset];
            std::fill(x_out_begin, x_out_begin + NDIM * n_qp, 0.0);

            double* N_begin = &N_qp[NDIM * qp_offset];
            std::fill(N_begin, N_begin + NDIM * n_qp, 0.0);

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
                    x_in_qp[NDIM * (qp_offset + qp) + i] = x_qp[NDIM * (qp_offset + qp) + i] - n(i) * dh;
                    x_out_qp[NDIM * (qp_offset + qp) + i] = x_qp[NDIM * (qp_offset + qp) + i] + n(i) * dh;
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
        // Note: Values are interpolated only to those quadrature points that
        // are within the patch interior

        const Box<NDIM>& interp_box = patch->getBox();

        Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(p_data_idx);

        const Box<NDIM> ghost_box = Box<NDIM>::grow(patch->getBox(), IntVector<NDIM>(p_ghost_num));

        LEInteractor::interpolate(P_i_qp, 1, x_in_qp, NDIM, p_data, patch, ghost_box, d_default_interp_spec.kernel_fcn);

        LEInteractor::interpolate(
            P_o_qp, 1, x_out_qp, NDIM, p_data, patch, ghost_box, d_default_interp_spec.kernel_fcn);

        std::vector<int> local_indices;
        local_indices.clear();
        const int upper_bound = n_qp_patch;
        if (upper_bound == 0) return;

        local_indices.reserve(upper_bound);
        for (unsigned int k = 0; k < n_qp_patch; ++k)
        {
            const double* const xx = &x_qp[NDIM * k];
            const Index<NDIM> i = IndexUtilities::getCellIndex(xx, patch_geom, interp_box);
            if (interp_box.contains(i)) local_indices.push_back(k);

            const double* const x_i = &x_in_qp[NDIM * k];
            const Index<NDIM> ip = IndexUtilities::getCellIndex(
                x_i, x_lower_ghost, x_upper_ghost, patch_geom->getDx(), ghost_box.lower(), ghost_box.upper());
            if (!ghost_box.contains(ip) && interp_box.contains(i))
                TBOX_ERROR(d_object_name << "::IIMethod():\n"
                                         << " the pressure interpolation ghost width hasn't beeen properly set"
                                         << std::endl);
            const double* const x_o = &x_out_qp[NDIM * k];
            const Index<NDIM> op = IndexUtilities::getCellIndex(
                x_o, x_lower_ghost, x_upper_ghost, patch_geom->getDx(), ghost_box.lower(), ghost_box.upper());
            if (!ghost_box.contains(op) && interp_box.contains(i))
                TBOX_ERROR(d_object_name << "::IIMethod():\n"
                                         << " the pressure interpolation ghost width hasn't beeen properly set"
                                         << std::endl);
        }

        const unsigned int nindices = static_cast<int>(local_indices.size());

        if (!local_indices.empty())
        {
            for (unsigned int k = 0; k < nindices; ++k)
            {
                P_out_qp[local_indices[k]] =
                    (1.0 - d_exterior_calc_coef) * (P_jump_qp[local_indices[k]] + P_i_qp[local_indices[k]]) +
                    d_exterior_calc_coef * P_o_qp[local_indices[k]];
                P_in_qp[local_indices[k]] = P_i_qp[local_indices[k]];
            }
        }

        // Loop over the elements and accumulate the right-hand-side values.
        qrule.reset();
        qp_offset = 0;
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];

            P_in_dof_map_cache.dof_indices(elem, P_in_dof_indices);
            P_in_rhs_e.resize(static_cast<int>(P_in_dof_indices.size()));

            P_out_dof_map_cache.dof_indices(elem, P_out_dof_indices);
            P_out_rhs_e.resize(static_cast<int>(P_out_dof_indices.size()));

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(elem, X_dof_indices[d], d);
            }
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
            const size_t n_basis2 = P_out_dof_indices.size();
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const int idx = qp_offset + qp;
                for (unsigned int k = 0; k < n_basis2; ++k)
                {
                    const double p_JxW = phi_P[k][qp] * JxW[qp];
                    P_in_rhs_e(k) += P_in_qp[idx] * p_JxW;
                    P_out_rhs_e(k) += P_out_qp[idx] * p_JxW;
                }
            }

            P_in_dof_map.constrain_element_vector(P_in_rhs_e, P_in_dof_indices);
            P_in_rhs_vec->add_vector(P_in_rhs_e, P_in_dof_indices);

            P_out_dof_map.constrain_element_vector(P_out_rhs_e, P_out_dof_indices);
            P_out_rhs_vec->add_vector(P_out_rhs_e, P_out_dof_indices);

            qp_offset += n_qp;
        }
    }
    P_in_rhs_vec->close();
    P_out_rhs_vec->close();

    d_fe_data_managers[part]->computeL2Projection(
        *P_in_vec, *P_in_rhs_vec, PRESSURE_IN_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
    P_in_vec->close();

    d_fe_data_managers[part]->computeL2Projection(
        *P_out_vec, *P_out_rhs_vec, PRESSURE_OUT_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
    P_out_vec->close();

    d_X_half_vecs[part]->close();
    d_X_current_vecs[part]->close();
    d_X_new_vecs[part]->close();

    d_P_in_half_vecs[part]->close();
    d_P_out_half_vecs[part]->close();
    d_X_IB_ghost_vecs[part]->close();
    d_P_jump_IB_ghost_vecs[part]->close();

    return;

} // extrapolatePressureForTraction

void
IIMethod::calculateInterfacialFluidForces(const int p_data_idx, double data_time)
{
    if (d_compute_fluid_traction && (!d_use_pressure_jump_conditions || !d_use_velocity_jump_conditions))
    {
        TBOX_ERROR(d_object_name << ": To compute the traction both velocity and pressure jumps need to be turned on!"
                                 << std::endl);
    }

    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const auto p_scratch_data_idx = d_eulerian_data_cache->getCachedPatchDataIndex(d_p_scratch_idx);
    RefineAlgorithm<NDIM> ghost_fill_alg_p;
    // TODO: Can we cache this algorithm/schedule?
    ghost_fill_alg_p.registerRefine(p_scratch_data_idx, p_data_idx, p_scratch_data_idx, NULL);
    Pointer<RefineSchedule<NDIM> > ghost_fill_schd_p =
        ghost_fill_alg_p.createSchedule(d_hierarchy->getPatchLevel(finest_ln));

    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        // TODO: Does this need to be re-filled for each part?
        ghost_fill_schd_p->fillData(data_time);
        extrapolatePressureForTraction(p_scratch_data_idx, data_time, part);
        computeFluidTraction(data_time, part);
    }

} // calculateInterfacialFluidForces

void
IIMethod::forwardEulerStep(const double current_time, const double new_time)
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
        else if (d_use_tangential_velocity[part])
        {
            ierr =
                VecWAXPY(d_X_new_vecs[part]->vec(), dt, d_U_t_current_vecs[part]->vec(), d_X_current_vecs[part]->vec());
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
} // forwardEulerStep

void
IIMethod::midpointStep(const double current_time, const double new_time)
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
        else if (d_use_tangential_velocity[part])
        {
            ierr = VecWAXPY(d_X_new_vecs[part]->vec(), dt, d_U_t_half_vecs[part]->vec(), d_X_current_vecs[part]->vec());
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
IIMethod::trapezoidalStep(const double current_time, const double new_time)
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
        else if (d_use_tangential_velocity[part])
        {
            ierr = VecWAXPY(
                d_X_new_vecs[part]->vec(), 0.5 * dt, d_U_t_current_vecs[part]->vec(), d_X_current_vecs[part]->vec());
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPY(d_X_new_vecs[part]->vec(), 0.5 * dt, d_U_t_new_vecs[part]->vec());
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
IIMethod::computeLagrangianForce(const double data_time)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
    batch_vec_ghost_update(d_X_half_vecs, INSERT_VALUES, SCATTER_FORWARD);
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

        NumericVector<double>* P_jump_vec = d_use_pressure_jump_conditions ? d_P_jump_half_vecs[part] : nullptr;
        std::unique_ptr<NumericVector<double> > P_jump_rhs_vec;
        DenseVector<double> P_jump_rhs_e;
        if (d_use_pressure_jump_conditions)
        {
            P_jump_rhs_vec = P_jump_vec->zero_clone();
        }
        double P_jump_rhs_integral = 0.0;

        std::array<NumericVector<double>*, NDIM> DU_jump_vec;
        std::array<std::unique_ptr<NumericVector<double> >, NDIM> DU_jump_rhs_vec;
        std::array<std::array<DenseVector<double>, NDIM>, NDIM> DU_jump_rhs_e;
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
        const DofMap* P_jump_dof_map = NULL;
        FEDataManager::SystemDofMapCache* P_jump_dof_map_cache = NULL;
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

        // The P_jump_fe_type and DU_jump_fe_type are equal only if we are applying both jumps
        // or none of the jumps.
        if (d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
        {
            TBOX_ASSERT(P_jump_fe_type == DU_jump_fe_type);
        }

        std::unique_ptr<QBase> qrule = QBase::build(d_default_quad_type[part], dim, d_default_quad_order[part]);

        std::unique_ptr<FEBase> fe_X = FEBase::build(dim, X_fe_type);
        fe_X->attach_quadrature_rule(qrule.get());
        const std::vector<double>& JxW = fe_X->get_JxW();
        const std::vector<std::vector<double> >& phi_X = fe_X->get_phi();
        std::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi_X;
        dphi_dxi_X[0] = &fe_X->get_dphidxi();
        if (NDIM > 2) dphi_dxi_X[1] = &fe_X->get_dphideta();

        FEType fe_jump_type = INVALID_FE;
        if (d_use_pressure_jump_conditions)
        {
            fe_jump_type = P_jump_fe_type;
        }
        else
        {
            fe_jump_type = DU_jump_fe_type;
        }

        std::unique_ptr<FEBase> fe_jump = FEBase::build(dim, fe_jump_type);
        fe_jump->attach_quadrature_rule(qrule.get());
        const std::vector<std::vector<double> >& phi_jump = fe_jump->get_phi();

        FEDataInterpolation fe_interpolator(dim, d_fe_data_managers[part]->getFEData());
        fe_interpolator.attachQuadratureRule(qrule.get());

        std::vector<size_t> surface_force_fcn_system_idxs;
        fe_interpolator.setupInterpolatedSystemDataIndexes(
            surface_force_fcn_system_idxs, d_lag_surface_force_fcn_data[part].system_data, equation_systems);
        std::vector<size_t> surface_pressure_fcn_system_idxs;
        fe_interpolator.setupInterpolatedSystemDataIndexes(
            surface_pressure_fcn_system_idxs, d_lag_surface_pressure_fcn_data[part].system_data, equation_systems);
        fe_interpolator.init();

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

            if (d_use_pressure_jump_conditions || d_use_velocity_jump_conditions)
            {
                fe_jump->reinit(elem);
            }

            get_values_for_interpolation(x_node, *X_vec, X_dof_indices);
            get_values_for_interpolation(X_node, X0_vec, X_dof_indices);
            const unsigned int n_qpoints = qrule->n_points();
            const size_t n_basis = phi_X.size();
            const size_t n_basis2 = phi_jump.size();
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
                        P_jump_rhs_e(k) += P_j * phi_jump[k][qp] * JxW[qp];
                    }
                    if (d_use_velocity_jump_conditions)
                    {
                        for (unsigned int i = 0; i < NDIM; ++i)
                        {
                            for (unsigned int j = 0; j < NDIM; ++j)
                            {
                                DU_jump_rhs_e[i][j](k) += DU[i][j] * phi_jump[k][qp] * JxW[qp];
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
                copy_dof_ids_to_vector(i, F_dof_indices, dof_id_scratch);
                F_dof_map.constrain_element_vector(F_rhs_e[i], dof_id_scratch);
                F_rhs_vec->add_vector(F_rhs_e[i], dof_id_scratch);
                if (d_use_velocity_jump_conditions)
                {
                    const auto& DU_jump_dof_indices = DU_jump_dof_map_cache[i]->dof_indices(elem);
                    for (unsigned int k = 0; k < NDIM; ++k)
                    {
                        copy_dof_ids_to_vector(k, DU_jump_dof_indices, dof_id_scratch);
                        DU_jump_dof_map[i]->constrain_element_vector(DU_jump_rhs_e[i][k], dof_id_scratch);
                        DU_jump_rhs_vec[i]->add_vector(DU_jump_rhs_e[i][k], dof_id_scratch);
                    }
                }
            }
            if (d_use_pressure_jump_conditions)
            {
                const auto& P_jump_dof_indices = P_jump_dof_map_cache->dof_indices(elem);
                copy_dof_ids_to_vector(0, P_jump_dof_indices, dof_id_scratch);
                P_jump_dof_map->constrain_element_vector(P_jump_rhs_e, dof_id_scratch);
                P_jump_rhs_vec->add_vector(P_jump_rhs_e, dof_id_scratch);
            }
        }

        SAMRAI_MPI::sumReduction(&F_integral(0), NDIM);

        // Solve for F.
        F_rhs_vec->close();
        d_fe_data_managers[part]->computeL2Projection(
            *F_vec, *F_rhs_vec, FORCE_SYSTEM_NAME, d_default_interp_spec.use_consistent_mass_matrix);
        F_vec->close();
        if (d_use_pressure_jump_conditions)
        {
            P_jump_rhs_vec->close();
            d_fe_data_managers[part]->computeL2Projection(*P_jump_vec,
                                                          *P_jump_rhs_vec,
                                                          PRESSURE_JUMP_SYSTEM_NAME,
                                                          d_default_interp_spec.use_consistent_mass_matrix);
            P_jump_rhs_integral = SAMRAI_MPI::sumReduction(P_jump_rhs_integral);
            surface_area = SAMRAI_MPI::sumReduction(surface_area);
            if (d_normalize_pressure_jump[part]) P_jump_vec->add(-P_jump_rhs_integral / surface_area);
            P_jump_vec->close();
        }
        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                DU_jump_rhs_vec[d]->close();
                d_fe_data_managers[part]->computeL2Projection(*DU_jump_vec[d],
                                                              *DU_jump_rhs_vec[d],
                                                              VELOCITY_JUMP_SYSTEM_NAME[d],
                                                              d_default_interp_spec.use_consistent_mass_matrix);
                DU_jump_vec[d]->close();
            }
        }
    }
    return;
} // computeLagrangianForce

void
IIMethod::spreadForce(const int f_data_idx,
                      RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                      const std::vector<Pointer<RefineSchedule<NDIM> > >& /*f_prolongation_scheds*/,
                      const double data_time)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));

    std::vector<std::vector<libMesh::PetscVector<double>*> > vec_collection_update = {
        d_X_IB_ghost_vecs, d_X_half_vecs, d_F_IB_ghost_vecs, d_F_half_vecs
    };

    if (d_use_pressure_jump_conditions)
    {
        vec_collection_update.push_back(d_P_jump_IB_ghost_vecs);
        vec_collection_update.push_back(d_P_jump_half_vecs);
    }

    if (d_use_velocity_jump_conditions)
    {
        for (unsigned part = 0; part < d_num_parts; ++part)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                vec_collection_update.push_back({ d_DU_jump_half_vecs[part][d], d_DU_jump_IB_ghost_vecs[part][d] });
            }
        }
    }

    batch_vec_ghost_update(vec_collection_update, INSERT_VALUES, SCATTER_FORWARD);

    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        PetscVector<double>* X_vec = d_X_half_vecs[part];
        PetscVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
        PetscVector<double>* F_vec = d_F_half_vecs[part];
        PetscVector<double>* F_ghost_vec = d_F_IB_ghost_vecs[part];
        X_vec->localize(*X_ghost_vec);
        F_vec->localize(*F_ghost_vec);
        d_fe_data_managers[part]->spread(
            f_data_idx, *F_ghost_vec, *X_ghost_vec, FORCE_SYSTEM_NAME, f_phys_bdry_op, data_time);
        PetscVector<double>* P_jump_vec;
        PetscVector<double>* P_jump_ghost_vec = NULL;
        std::array<PetscVector<double>*, NDIM> DU_jump_ghost_vec;
        std::array<PetscVector<double>*, NDIM> DU_jump_vec;
        if (d_use_pressure_jump_conditions)
        {
            P_jump_vec = d_P_jump_half_vecs[part];
            P_jump_ghost_vec = d_P_jump_IB_ghost_vecs[part];
            P_jump_vec->localize(*P_jump_ghost_vec);
        }
        if (d_use_velocity_jump_conditions)
        {
            for (auto d = 0; d < NDIM; ++d)
            {
                DU_jump_ghost_vec[d] = d_DU_jump_IB_ghost_vecs[part][d];
                DU_jump_vec[d] = d_DU_jump_half_vecs[part][d];
                DU_jump_vec[d]->localize(*DU_jump_ghost_vec[d]);
            }
        }

        if (d_use_pressure_jump_conditions || d_use_velocity_jump_conditions)
        {
            imposeJumpConditions(f_data_idx, *P_jump_ghost_vec, DU_jump_ghost_vec, *X_ghost_vec, data_time, part);
        }
        if (d_use_velocity_jump_conditions)
        {
            for (unsigned int d = 0; d < NDIM; ++d) d_DU_jump_IB_ghost_vecs[part][d]->close();
        }
        if (d_use_pressure_jump_conditions)
        {
            d_P_jump_IB_ghost_vecs[part]->close();
        }

        d_F_IB_ghost_vecs[part]->close();
        d_X_IB_ghost_vecs[part]->close();
    }
    return;
} // spreadForce

FEDataManager::InterpSpec
IIMethod::getDefaultInterpSpec() const
{
    return d_default_interp_spec;
} // getDefaultInterpSpec

FEDataManager::SpreadSpec
IIMethod::getDefaultSpreadSpec() const
{
    return d_default_spread_spec;
} // getDefaultSpreadSpec

void
IIMethod::setInterpSpec(const FEDataManager::InterpSpec& interp_spec, const unsigned int part)
{
    TBOX_ASSERT(!d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_num_parts);
    d_interp_spec[part] = interp_spec;
    return;
} // setInterpSpec

void
IIMethod::setSpreadSpec(const FEDataManager::SpreadSpec& spread_spec, const unsigned int part)
{
    TBOX_ASSERT(!d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_num_parts);
    d_spread_spec[part] = spread_spec;
    return;
} // setSpreadSpec

void
IIMethod::initializeFEEquationSystems()
{
    if (d_fe_equation_systems_initialized) return;

    const bool from_restart = RestartManager::getManager()->isFromRestart();

    // Create the FE data managers that manage mappings between the FE mesh
    // parts and the Cartesian grid.
    d_equation_systems.resize(d_num_parts, nullptr);
    d_fe_data_managers.resize(d_num_parts, nullptr);
    IntVector<NDIM> min_ghost_width(0);
    if (!d_eulerian_data_cache) d_eulerian_data_cache.reset(new SAMRAIDataCache());
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Create FE equation systems objects and corresponding variables.
        d_equation_systems[part] = new EquationSystems(*d_meshes[part]);
        EquationSystems* equation_systems = d_equation_systems[part];
        auto fe_data = std::make_shared<FEData>(
            d_object_name + "::FEdata::" + std::to_string(part), *equation_systems, d_registered_for_restart);

        // Create FE data managers.
        const std::string manager_name = "IIMethod FEDataManager::" + std::to_string(part);
        Pointer<InputDatabase> fe_data_manager_db(new InputDatabase(manager_name + "::input_db"));

        d_fe_data_managers[part] = FEDataManager::getManager(fe_data,
                                                             manager_name,
                                                             fe_data_manager_db,
                                                             d_max_level_number + 1,
                                                             d_interp_spec[part],
                                                             d_spread_spec[part],
                                                             d_default_workload_spec,
                                                             min_ghost_width,
                                                             d_eulerian_data_cache);
        d_ghosts = IntVector<NDIM>::max(d_ghosts, d_fe_data_managers[part]->getGhostCellWidth());
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
                System& P_in_system = equation_systems->add_system<System>(PRESSURE_IN_SYSTEM_NAME);
                System& P_out_system = equation_systems->add_system<System>(PRESSURE_OUT_SYSTEM_NAME);
                if (d_use_discon_elem_for_jumps[part])
                {
                    P_jump_system.add_variable("P_jump_", d_fe_order[part], d_pressure_jump_fe_family);
                    P_in_system.add_variable("P_in_", d_fe_order[part], d_pressure_jump_fe_family);
                    P_out_system.add_variable("P_out_", d_fe_order[part], d_pressure_jump_fe_family);
                }
                else
                {
                    P_jump_system.add_variable("P_jump_", d_fe_order[part], d_fe_family[part]);
                    P_in_system.add_variable("P_in_", d_fe_order[part], d_fe_family[part]);
                    P_out_system.add_variable("P_out_", d_fe_order[part], d_fe_family[part]);
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
                        if (d_use_discon_elem_for_jumps[part])
                        {
                            DU_jump_system[d]->add_variable(system_name, d_fe_order[part], d_velocity_jump_fe_family);
                        }
                        else
                        {
                            DU_jump_system[d]->add_variable(system_name, d_fe_order[part], d_fe_family[part]);
                        }
                    }
                }

                System& WSS_in_system = equation_systems->add_system<System>(WSS_IN_SYSTEM_NAME);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    const std::string system_name = "WSS_in_" + std::to_string(d);
                    if (d_use_discon_elem_for_jumps[part])
                    {
                        WSS_in_system.add_variable(system_name, d_fe_order[part], d_wss_fe_family);
                    }
                    else
                    {
                        WSS_in_system.add_variable(system_name, d_fe_order[part], d_fe_family[part]);
                    }
                }

                System& WSS_out_system = equation_systems->add_system<System>(WSS_OUT_SYSTEM_NAME);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    const std::string system_name = "WSS_out_" + std::to_string(d);
                    if (d_use_discon_elem_for_jumps[part])
                    {
                        WSS_out_system.add_variable(system_name, d_fe_order[part], d_wss_fe_family);
                    }
                    else
                    {
                        WSS_out_system.add_variable(system_name, d_fe_order[part], d_fe_family[part]);
                    }
                }
            }

            if (d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
            {
                auto& TAU_in_system = equation_systems->add_system<System>(TAU_IN_SYSTEM_NAME);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    std::string system_name = "TAU_IN_" + std::to_string(d);
                    if (d_use_discon_elem_for_jumps[part])
                    {
                        TAU_in_system.add_variable(system_name, d_fe_order[part], d_tau_fe_family);
                    }
                    else
                    {
                        TAU_in_system.add_variable(system_name, d_fe_order[part], d_fe_family[part]);
                    }
                }

                auto& TAU_out_system = equation_systems->add_system<System>(TAU_OUT_SYSTEM_NAME);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    std::string system_name = "TAU_OUT_" + std::to_string(d);
                    if (d_use_discon_elem_for_jumps[part])
                    {
                        TAU_out_system.add_variable(system_name, d_fe_order[part], d_tau_fe_family);
                    }
                    else
                    {
                        TAU_out_system.add_variable(system_name, d_fe_order[part], d_fe_family[part]);
                    }
                }
            }

            auto insert_parallel_into_ghosted = [](const PetscVector<Number>& parallel_vector,
                                                   PetscVector<Number>& ghosted_vector) {
                TBOX_ASSERT(parallel_vector.size() == ghosted_vector.size());
                TBOX_ASSERT(parallel_vector.local_size() == ghosted_vector.local_size());
                ghosted_vector = parallel_vector;
                ghosted_vector.close();
            };

            const std::array<std::string, 1> system_names{ { COORDS_SYSTEM_NAME } };
            const std::array<std::string, 1> vector_names{ { "INITIAL_COORDINATES" } };
            for (const std::string& system_name : system_names)
            {
                auto& system = equation_systems->get_system(system_name);
                for (const std::string& vector_name : vector_names)
                {
                    std::unique_ptr<NumericVector<double> > clone_vector;
                    if (from_restart)
                    {
                        NumericVector<double>* current = system.request_vector(vector_name);
                        if (current != nullptr)
                        {
                            clone_vector = current->clone();
                        }
                    }
                    system.remove_vector(vector_name);
                    system.add_vector(vector_name, /*projections*/ true, /*type*/ GHOSTED);

                    if (clone_vector != nullptr)
                    {
                        const auto& parallel_vector = dynamic_cast<const PetscVector<Number>&>(*clone_vector);
                        auto& ghosted_vector = dynamic_cast<PetscVector<Number>&>(system.get_vector(vector_name));
                        insert_parallel_into_ghosted(parallel_vector, ghosted_vector);
                    }
                }
            }
        }
    }
    d_fe_equation_systems_initialized = true;
    return;
} // initializeFEEquationSystems

void
IIMethod::initializeFEData()
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

            System& P_in_system = equation_systems->get_system<System>(PRESSURE_IN_SYSTEM_NAME);
            P_in_system.assemble_before_solve = false;
            P_in_system.assemble();

            System& P_out_system = equation_systems->get_system<System>(PRESSURE_OUT_SYSTEM_NAME);
            P_out_system.assemble_before_solve = false;
            P_out_system.assemble();
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
            System& WSS_in_system = equation_systems->get_system<System>(WSS_IN_SYSTEM_NAME);
            WSS_in_system.assemble_before_solve = false;
            WSS_in_system.assemble();

            System& WSS_out_system = equation_systems->get_system<System>(WSS_OUT_SYSTEM_NAME);
            WSS_out_system.assemble_before_solve = false;
            WSS_out_system.assemble();
        }

        if (d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
        {
            System& TAU_in_system = equation_systems->get_system<System>(TAU_IN_SYSTEM_NAME);
            TAU_in_system.assemble_before_solve = false;
            TAU_in_system.assemble();

            System& TAU_out_system = equation_systems->get_system<System>(TAU_OUT_SYSTEM_NAME);
            TAU_out_system.assemble_before_solve = false;
            TAU_out_system.assemble();
        }
    }
    d_fe_data_initialized = true;
    return;
} // initializeFEData

void
IIMethod::registerEulerianVariables()
{
    d_p_var = new CellVariable<NDIM, double>(d_object_name + "::p");
    registerVariable(d_p_scratch_idx, d_p_var, d_ghosts);
    return;
} // registerEulerianVariables

void
IIMethod::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
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
    d_eulerian_data_cache->setPatchHierarchy(hierarchy);
    d_eulerian_data_cache->resetLevels(0, hierarchy->getFinestLevelNumber()); // TODO: implement
                                                                              // this->getFinestPatchLevelNumber()

    d_is_initialized = true;
    return;
} // initializePatchHierarchy

void
IIMethod::registerLoadBalancer(Pointer<LoadBalancer<NDIM> > load_balancer, int workload_data_idx)
{
    IBAMR_DEPRECATED_MEMBER_FUNCTION1("IIMethod", "registerLoadBalancer");
    TBOX_ASSERT(load_balancer);
    d_load_balancer = load_balancer;
    d_workload_idx = workload_data_idx;

    return;
} // registerLoadBalancer

void
IIMethod::addWorkloadEstimate(Pointer<PatchHierarchy<NDIM> > hierarchy, const int workload_data_idx)
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->addWorkloadEstimate(hierarchy, workload_data_idx);
    }
    return;
} // addWorkloadEstimate

void IIMethod::beginDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                       Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    // intentionally blank
    return;
} // beginDataRedistribution

void IIMethod::endDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
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
IIMethod::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                              int /*level_number*/,
                              double /*init_data_time*/,
                              bool /*can_be_refined*/,
                              bool /*initial_time*/,
                              Pointer<BasePatchLevel<NDIM> > /*old_level*/,
                              bool /*allocate_data*/)
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->setPatchHierarchy(hierarchy);
    }
    return;
} // initializeLevelData

void
IIMethod::resetHierarchyConfiguration(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                      int /*coarsest_level*/,
                                      int /*finest_level*/)
{
    // const int finest_hier_level = hierarchy->getFinestLevelNumber();
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->setPatchHierarchy(hierarchy);
    }
    return;
} // resetHierarchyConfiguration

void
IIMethod::applyGradientDetector(Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
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
IIMethod::putToDatabase(Pointer<Database> db)
{
    db->putInteger("IBFE_METHOD_VERSION", IBFE_METHOD_VERSION);
    db->putInteger("d_num_parts", d_num_parts);
    db->putIntegerArray("d_ghosts", d_ghosts, NDIM);
    db->putBool("d_use_velocity_jump_conditions", d_use_velocity_jump_conditions);
    db->putBool("d_use_pressure_jump_conditions", d_use_pressure_jump_conditions);
    db->putBool("d_compute_fluid_traction", d_compute_fluid_traction);
    db->putBool("d_use_consistent_mass_matrix", d_use_consistent_mass_matrix);
    db->putBool("d_use_direct_forcing", d_use_direct_forcing);
    return;
} // putToDatabase

void
IIMethod::writeFEDataToRestartFile(const std::string& restart_dump_dirname, unsigned int time_step_number)
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
} // writeFEDataToRestartFile

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
IIMethod::imposeJumpConditions(const int f_data_idx,
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
    for (unsigned int d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);

    System* P_jump_system;
    const DofMap* P_jump_dof_map;
    FEDataManager::SystemDofMapCache* P_jump_dof_map_cache = NULL;
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
    const int level_num = d_fe_data_managers[part]->getFinestPatchLevelNumber();
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

        Box<NDIM> side_boxes[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            side_boxes[d] = SideGeometry<NDIM>::toSideBox(patch_box, d);
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
        std::vector<libMesh::dof_id_type> dof_id_scratch;
        for (size_t e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
            get_values_for_interpolation(x_node, X_ghost_vec, X_dof_indices);
            if (d_use_pressure_jump_conditions)
            {
                const auto& P_jump_dof_indices = P_jump_dof_map_cache->dof_indices(elem);
                copy_dof_ids_to_vector(0, P_jump_dof_indices, dof_id_scratch);
                get_values_for_interpolation(P_jump_node, P_jump_ghost_vec, dof_id_scratch);
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
                        const int i_s = boost::math::iround(((x(d) - x_lower[d]) / dx[d]) - 0.5) + patch_lower[d];
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
                Box<NDIM> extended_side_box = patch_box;
                extended_side_box.grow(IntVector<NDIM>(2));
                if (patch_geom->getTouchesRegularBoundary(axis, 1)) extended_box.upper(axis) += 1;

                Box<NDIM> side_u_boxes[NDIM];
                for (int d = 0; d < NDIM; ++d)
                {
                    side_u_boxes[d] = SideGeometry<NDIM>::toSideBox(extended_side_box, d);
                }

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
                    const hier::Index<NDIM>& i_c = b();
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
                                             x_lower[d] + dx[d] * (static_cast<double>(i_c(d) - patch_lower[d]) + 0.5));
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

                                    found_same_intersection_point =
                                        checkDoubleCountingIntersection(axis,
                                                                        dx,
                                                                        n,
                                                                        x,
                                                                        xi,
                                                                        i_s,
                                                                        i_s_prime,
                                                                        candidate_coords,
                                                                        candidate_ref_coords,
                                                                        candidate_normals);
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
                            i_s_up(axis) =
                                boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] + 0.5) + patch_lower[axis];
                            i_s_um(axis) =
                                boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] - 0.5) + patch_lower[axis];

                            if (extended_box.contains(i_s_up) && extended_box.contains(i_s_um))
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

                                    found_same_intersection_point =
                                        checkDoubleCountingIntersection(axis,
                                                                        dx,
                                                                        n,
                                                                        xu,
                                                                        xui,
                                                                        i_s_um,
                                                                        i_s_prime,
                                                                        candidate_coords,
                                                                        candidate_ref_coords,
                                                                        candidate_normals);
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
                                    if (side_ghost_boxes[axis].contains(i_s_up) &&
                                        side_ghost_boxes[axis].contains(i_s_um))
                                    {
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

                                    // Keep track of the positions where we have
                                    // imposed jump conditions.
                                    intersection_u_points[axis][i_s_um].push_back(xu);
                                    intersection_u_ref_coords[axis][i_s_um].push_back(xui);
                                    intersection_u_normals[axis][i_s_um].push_back(n);
                                }
                            }
                        }

                        for (unsigned int j = 0; j < NDIM - 1; ++j)
                        {
                            for (unsigned int k = 0; k < intersectionsSide[j].size(); ++k)
                            {
                                libMesh::Point xu = rs[j] + intersectionsSide[j][k].first * q;
                                const libMesh::Point& xui = intersectionsSide[j][k].second;
                                SideIndex<NDIM> i_s_up;
                                SideIndex<NDIM> i_s_um;

                                if (xu(axis) - x_lower[axis] > 0.0)
                                {
                                    if (fmod(xu(axis) - x_lower[axis], dx[axis]) >= 0.5 * dx[axis])
                                    {
                                        SideIndex<NDIM> i_side_um(i_c, SideDim[axis][j], 0);
                                        Index<NDIM> i_c_neighbor = i_c;
                                        i_c_neighbor(axis) += 1;

                                        SideIndex<NDIM> i_side_up(i_c_neighbor, SideDim[axis][j], 0);

                                        i_side_up(axis) = boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis]) +
                                                          patch_lower[axis];
                                        i_side_um(axis) =
                                            boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] - 0.5) +
                                            patch_lower[axis];
                                        i_s_up = i_side_up;
                                        i_s_um = i_side_um;
                                    }
                                    else if (fmod((xu(axis) - x_lower[axis]), dx[axis]) < 0.5 * dx[axis])
                                    {
                                        SideIndex<NDIM> i_side_up(i_c, SideDim[axis][j], 0);
                                        Index<NDIM> i_c_neighbor = i_c;
                                        i_c_neighbor(axis) -= 1;
                                        SideIndex<NDIM> i_side_um(i_c_neighbor, SideDim[axis][j], 0);
                                        i_side_up(axis) =
                                            boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] - 0.5) +
                                            patch_lower[axis];
                                        i_side_um(axis) =
                                            boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] - 1.0) +
                                            patch_lower[axis];
                                        i_s_up = i_side_up;
                                        i_s_um = i_side_um;
                                    }
                                    else
                                    {
                                        continue;
                                    }
                                }
                                else if (xu(axis) - x_lower[axis] < 0.0)
                                {
                                    if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) < 0.5 * dx[axis])
                                    {
                                        SideIndex<NDIM> i_side_um(i_c, SideDim[axis][j], 0);
                                        Index<NDIM> i_c_neighbor = i_c;
                                        i_c_neighbor(axis) += 1;

                                        SideIndex<NDIM> i_side_up(i_c_neighbor, SideDim[axis][j], 0);

                                        i_side_up(axis) = boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis]) +
                                                          patch_lower[axis];
                                        i_side_um(axis) =
                                            boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] - 0.5) +
                                            patch_lower[axis];
                                        i_s_up = i_side_up;
                                        i_s_um = i_side_um;
                                    }
                                    else
                                    {
                                        SideIndex<NDIM> i_side_up(i_c, SideDim[axis][j], 0);
                                        Index<NDIM> i_c_neighbor = i_c;
                                        i_c_neighbor(axis) -= 1;
                                        SideIndex<NDIM> i_side_um(i_c_neighbor, SideDim[axis][j], 0);
                                        i_side_up(axis) =
                                            boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] - 0.5) +
                                            patch_lower[axis];
                                        i_side_um(axis) =
                                            boost::math::iround((xu(axis) - x_lower[axis]) / dx[axis] - 1.0) +
                                            patch_lower[axis];
                                        i_s_up = i_side_up;
                                        i_s_um = i_side_um;
                                    }
                                }
                                else
                                {
                                    TBOX_ERROR(d_object_name << ":  Restart file version different than class version."
                                                             << std::endl);
                                }

                                if (extended_side_box.contains(i_s_up) && extended_side_box.contains(i_s_um))
                                {
                                    TBOX_ASSERT(i_s_up(axis) - i_s_um(axis) == 1);
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

                                        found_same_intersection_point =
                                            checkDoubleCountingIntersection(axis,
                                                                            dx,
                                                                            n,
                                                                            xu,
                                                                            xui,
                                                                            i_s_um,
                                                                            i_s_prime,
                                                                            candidate_coords,
                                                                            candidate_ref_coords,
                                                                            candidate_normals);
                                        if (found_same_intersection_point) break;
                                    }

                                    if (!found_same_intersection_point)
                                    {
                                        // Evaluate the jump conditions and apply them
                                        // to the Eulerian grid.

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
                                        if (side_ghost_boxes[SideDim[axis][j]].contains(i_s_up) &&
                                            side_ghost_boxes[SideDim[axis][j]].contains(i_s_um))
                                        {
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

bool
IIMethod::checkDoubleCountingIntersection(const int axis,
                                          const double* const dx,
                                          const libMesh::VectorValue<double>& n,
                                          const libMesh::Point& x,
                                          const libMesh::Point& xi,
                                          const SideIndex<NDIM>& i_s,
                                          const SideIndex<NDIM>& i_s_prime,
                                          const std::vector<libMesh::Point>& candidate_coords,
                                          const std::vector<libMesh::Point>& candidate_ref_coords,
                                          const std::vector<libMesh::VectorValue<double> >& candidate_normals)
{
    bool found_same_intersection_point = false;
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
    return found_same_intersection_point;
} // checkDoubleCountingIntersection

void
IIMethod::initializeCoordinates(const unsigned int part)
{
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    const bool identity_mapping = !d_coordinate_mapping_fcn_data[part].fcn;
    auto it = mesh.local_nodes_begin();
    const auto end_it = mesh.local_nodes_end();
    for (; it != end_it; ++it)
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
IIMethod::updateCoordinateMapping(const unsigned int part)
{
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    System& dX_system = equation_systems->get_system(COORD_MAPPING_SYSTEM_NAME);
    const unsigned int dX_sys_num = dX_system.number();
    NumericVector<double>& dX_coords = *dX_system.solution;
    auto it = mesh.local_nodes_begin();
    const auto end_it = mesh.local_nodes_end();
    for (; it != end_it; ++it)
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
IIMethod::initializeVelocity(const unsigned int part)
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
        auto it = mesh.local_nodes_begin();
        const auto end_it = mesh.local_nodes_end();
        for (; it != end_it; ++it)
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
IIMethod::commonConstructor(const std::string& object_name,
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
    d_max_level_number = max_level_number - 1;

    // Set some default values.
    const bool use_adaptive_quadrature = true;
    const int point_density = 2.0;
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

    d_use_tangential_velocity.resize(d_num_parts);
    d_normalize_pressure_jump.resize(d_num_parts);
    d_use_discon_elem_for_jumps.resize(d_num_parts);

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
                       << "::IIMethod():\n"
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
IIMethod::getFromInput(Pointer<Database> db, bool /*is_from_restart*/)
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
    if (d_use_pressure_jump_conditions)
    {
        if (db->isString("pressure_jump_fe_family"))
            d_pressure_jump_fe_family = Utility::string_to_enum<FEFamily>(db->getString("pressure_jump_fe_family"));
    }

    if (db->isBool("use_velocity_jump_conditions"))
        d_use_velocity_jump_conditions = db->getBool("use_velocity_jump_conditions");
    if (d_use_velocity_jump_conditions)
    {
        if (db->isDouble("wss_calc_width")) d_wss_calc_width = db->getDouble("wss_calc_width");
        if (db->isString("velocity_jump_fe_family"))
            d_velocity_jump_fe_family = Utility::string_to_enum<FEFamily>(db->getString("velocity_jump_fe_family"));
        if (db->isString("wss_fe_family"))
            d_wss_fe_family = Utility::string_to_enum<FEFamily>(db->getString("wss_fe_family"));
    }
    if (d_use_pressure_jump_conditions && d_use_velocity_jump_conditions)
    {
        if (db->isDouble("wss_calc_width")) d_wss_calc_width = db->getDouble("wss_calc_width");
        if (db->isString("tau_fe_family"))
            d_tau_fe_family = Utility::string_to_enum<FEFamily>(db->getString("tau_fe_family"));
    }
    if (d_use_pressure_jump_conditions || d_use_velocity_jump_conditions)
    {
        if (db->isBool("perturb_fe_mesh_nodes")) d_perturb_fe_mesh_nodes = db->getBool("perturb_fe_mesh_nodes");
    }
    if (db->isBool("compute_fluid_traction")) d_compute_fluid_traction = db->getBool("compute_fluid_traction");
    if (d_compute_fluid_traction)
    {
        if (db->isDouble("p_calc_width")) d_p_calc_width = db->getDouble("p_calc_width");
        if (db->isDouble("exterior_calc_coef")) d_exterior_calc_coef = db->getDouble("exterior_calc_coef");
    }
    if (db->isBool("use_consistent_mass_matrix"))
        d_use_consistent_mass_matrix = db->getBool("use_consistent_mass_matrix");
    if (db->isBool("use_direct_forcing")) d_use_direct_forcing = db->getBool("use_direct_forcing");

    // Restart settings.
    if (db->isString("libmesh_restart_file_extension"))
        d_libmesh_restart_file_extension = db->getString("libmesh_restart_file_extension");

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
IIMethod::getFromRestart()
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
    d_use_pressure_jump_conditions = db->getBool("d_use_pressure_jump_conditions");
    d_use_velocity_jump_conditions = db->getBool("d_use_velocity_jump_conditions");
    d_compute_fluid_traction = db->getBool("d_compute_fluid_traction");
    d_use_consistent_mass_matrix = db->getBool("d_use_consistent_mass_matrix");
    d_use_direct_forcing = db->getBool("d_use_direct_forcing");
    d_exterior_calc_coef = db->getDouble("exterior_calc_coef");
    d_p_calc_width = db->getDouble("p_calc_width");
    d_wss_calc_width = db->getDouble("wss_calc_width");
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
