// Filename: IBFEMethod.cpp
// Created on 5 Oct 2011 by Boyce Griffith
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

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <set>
#include <stdbool.h>
#include <stddef.h>
#include <string>
#include <utility>
#include <vector>

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
#include "SideGeometry.h"
#include "SideIndex.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "boost/array.hpp"
#include "boost/multi_array.hpp"
#include "ibamr/IBFEMethod.h"
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
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/boundary_info.h"
#include "libmesh/compare_types.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/edge.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_type.h"
#include "libmesh/fem_context.h"
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
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

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
// Version of IBFEMethod restart file data.
static const int IBFE_METHOD_VERSION = 1;

inline short int
get_dirichlet_bdry_ids(const std::vector<short int>& bdry_ids)
{
    short int dirichlet_bdry_ids = 0;
    for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end(); ++cit)
    {
        const short int bdry_id = *cit;
        if (bdry_id == FEDataManager::ZERO_DISPLACEMENT_X_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_Y_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_Z_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_XY_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_XZ_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_YZ_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_XYZ_BDRY_ID)
        {
            dirichlet_bdry_ids |= bdry_id;
        }
    }
    return dirichlet_bdry_ids;
}

inline bool
is_physical_bdry(const Elem* elem,
                 const unsigned short int side,
                 const BoundaryInfo& boundary_info,
                 const DofMap& dof_map)
{
    const std::vector<short int>& bdry_ids = boundary_info.boundary_ids(elem, side);
    bool at_physical_bdry = !elem->neighbor(side);
    for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end(); ++cit)
    {
        if (dof_map.is_periodic_boundary(*cit)) at_physical_bdry = false;
    }
    return at_physical_bdry;
}

inline bool
is_dirichlet_bdry(const Elem* elem,
                  const unsigned short int side,
                  const BoundaryInfo& boundary_info,
                  const DofMap& dof_map)
{
    if (!is_physical_bdry(elem, side, boundary_info, dof_map)) return false;
    const std::vector<short int>& bdry_ids = boundary_info.boundary_ids(elem, side);
    return get_dirichlet_bdry_ids(bdry_ids) != 0;
}

inline bool
has_physical_bdry(const Elem* elem, const BoundaryInfo& boundary_info, const DofMap& dof_map)
{
    bool has_physical_bdry = false;
    for (unsigned short int side = 0; side < elem->n_sides() && !has_physical_bdry; ++side)
    {
        has_physical_bdry = has_physical_bdry || is_physical_bdry(elem, side, boundary_info, dof_map);
    }
    return has_physical_bdry;
}

inline void
get_x_and_FF(libMesh::VectorValue<double>& x,
             libMesh::TensorValue<double>& FF,
             const std::vector<double>& x_data,
             const std::vector<VectorValue<double> >& grad_x_data,
             const unsigned int dim = NDIM)
{
    x.zero();
    FF.zero();
    for (unsigned int i = 0; i < dim; ++i)
    {
        x(i) = x_data[i];
        for (unsigned int j = 0; j < dim; ++j)
        {
            FF(i, j) = grad_x_data[i](j);
        }
    }
    for (unsigned int i = dim; i < LIBMESH_DIM; ++i)
    {
        FF(i, i) = 1.0;
    }
    return;
}

static const Real PENALTY = 1.e10;

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

inline int
NINT(double a)
{
    return (a >= 0.0 ? static_cast<int>(a + 0.5) : static_cast<int>(a - 0.5));
}
}

const std::string IBFEMethod::COORDS_SYSTEM_NAME = "IB coordinates system";
const std::string IBFEMethod::COORDS0_SYSTEM_NAME = "IB initial coordinates system";
const std::string IBFEMethod::COORD_MAPPING_SYSTEM_NAME = "IB coordinate mapping system";
const std::string IBFEMethod::P_J_SYSTEM_NAME = "IB p jump system";
const std::string IBFEMethod::H_SYSTEM_NAME = "Curvature system";
const std::string IBFEMethod::DP_J_SYSTEM_NAME = "IB dp jump system";
const std::string IBFEMethod::DU_J_SYSTEM_NAME = "IB velocity du jump system";
const std::string IBFEMethod::DV_J_SYSTEM_NAME = "IB velocity dv jump system";
const std::string IBFEMethod::DW_J_SYSTEM_NAME = "IB velocity dw jump system";
const std::string IBFEMethod::D2U_J_SYSTEM_NAME = "IB velocity d2u jump system";
const std::string IBFEMethod::D2V_J_SYSTEM_NAME = "IB velocity d2v jump system";
const std::string IBFEMethod::D2W_J_SYSTEM_NAME = "IB velocity d2w jump system";
const std::string IBFEMethod::FORCE_SYSTEM_NAME = "IB force system";
const std::string IBFEMethod::FORCE_T_SYSTEM_NAME = "IB tangential t force system";
const std::string IBFEMethod::FORCE_B_SYSTEM_NAME = "IB tangential b force system";
const std::string IBFEMethod::FORCE_N_SYSTEM_NAME = "IB normal force system";
const std::string IBFEMethod::VELOCITY_SYSTEM_NAME = "IB velocity system";
const std::string IBFEMethod::WSS_I_SYSTEM_NAME = "One sided inside shear stress system";
const std::string IBFEMethod::WSS_O_SYSTEM_NAME = "One sided outside shear stress system";

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEMethod::IBFEMethod(const std::string& object_name,
                       Pointer<Database> input_db,
                       Mesh* mesh,
                       int max_level_number,
                       bool register_for_restart,
                       const std::string& restart_read_dirname,
                       unsigned int restart_restore_number)
    : d_num_parts(1)
{
    commonConstructor(object_name,
                      input_db,
                      std::vector<Mesh*>(1, mesh),
                      max_level_number,
                      register_for_restart,
                      restart_read_dirname,
                      restart_restore_number);
    return;
} // IBFEMethod

IBFEMethod::IBFEMethod(const std::string& object_name,
                       Pointer<Database> input_db,
                       const std::vector<Mesh*>& meshes,
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
} // IBFEMethod

IBFEMethod::~IBFEMethod()
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
} // ~IBFEMethod

FEDataManager*
IBFEMethod::getFEDataManager(const unsigned int part) const
{
    TBOX_ASSERT(d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_num_parts);
    return d_fe_data_managers[part];
} // getFEDataManager

void
IBFEMethod::registerInitialCoordinateMappingFunction(const CoordinateMappingFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_coordinate_mapping_fcn_data[part] = data;
    return;
} // registerInitialCoordinateMappingFunction

void
IBFEMethod::registerLagForceFunction(const LagForceFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_force_fcn_data[part] = data;
    return;
} // registerLagForceFunction

const IntVector<NDIM>&
IBFEMethod::getMinimumGhostCellWidth() const
{
    return d_ghosts;
} // getMinimumGhostCellWidth

void
IBFEMethod::setupTagBuffer(Array<int>& tag_buffer, Pointer<GriddingAlgorithm<NDIM> > gridding_alg) const
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
IBFEMethod::preprocessIntegrateData(double current_time, double new_time, int /*num_cycles*/)
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

    d_X0_systems.resize(d_num_parts);
    d_X0_vecs.resize(d_num_parts);

    d_U_systems.resize(d_num_parts);
    d_U_current_vecs.resize(d_num_parts);
    d_U_new_vecs.resize(d_num_parts);
    d_U_half_vecs.resize(d_num_parts);

    d_F_systems.resize(d_num_parts);
    d_F_half_vecs.resize(d_num_parts);
    d_F_IB_ghost_vecs.resize(d_num_parts);

    d_F_t_systems.resize(d_num_parts);
    d_F_t_half_vecs.resize(d_num_parts);
    d_F_t_IB_ghost_vecs.resize(d_num_parts);

    d_F_b_systems.resize(d_num_parts);
    d_F_b_half_vecs.resize(d_num_parts);
    d_F_b_IB_ghost_vecs.resize(d_num_parts);
    
    d_F_n_systems.resize(d_num_parts);
    d_F_n_half_vecs.resize(d_num_parts);
    d_F_n_IB_ghost_vecs.resize(d_num_parts);

    d_P_j_systems.resize(d_num_parts);
    d_P_j_half_vecs.resize(d_num_parts);
    d_P_j_IB_ghost_vecs.resize(d_num_parts);
    
    d_H_systems.resize(d_num_parts);
    d_H_half_vecs.resize(d_num_parts);
    d_H_IB_ghost_vecs.resize(d_num_parts);

    d_dP_j_systems.resize(d_num_parts);
    d_dP_j_half_vecs.resize(d_num_parts);
    d_dP_j_IB_ghost_vecs.resize(d_num_parts);

    d_du_j_systems.resize(d_num_parts);
    d_du_j_half_vecs.resize(d_num_parts);
    d_du_j_IB_ghost_vecs.resize(d_num_parts);

    d_dv_j_systems.resize(d_num_parts);
    d_dv_j_half_vecs.resize(d_num_parts);
    d_dv_j_IB_ghost_vecs.resize(d_num_parts);

    d_dw_j_systems.resize(d_num_parts);
    d_dw_j_half_vecs.resize(d_num_parts);
    d_dw_j_IB_ghost_vecs.resize(d_num_parts);

    d_d2u_j_systems.resize(d_num_parts);
    d_d2u_j_half_vecs.resize(d_num_parts);
    d_d2u_j_IB_ghost_vecs.resize(d_num_parts);

    d_d2v_j_systems.resize(d_num_parts);
    d_d2v_j_half_vecs.resize(d_num_parts);
    d_d2v_j_IB_ghost_vecs.resize(d_num_parts);

    d_d2w_j_systems.resize(d_num_parts);
    d_d2w_j_half_vecs.resize(d_num_parts);
    d_d2w_j_IB_ghost_vecs.resize(d_num_parts);

    d_WSS_i_systems.resize(d_num_parts);
    d_WSS_i_half_vecs.resize(d_num_parts);
    d_WSS_i_IB_ghost_vecs.resize(d_num_parts);

    d_WSS_o_systems.resize(d_num_parts);
    d_WSS_o_half_vecs.resize(d_num_parts);
    d_WSS_o_IB_ghost_vecs.resize(d_num_parts);
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

        d_X0_systems[part] = &d_equation_systems[part]->get_system(COORDS0_SYSTEM_NAME);
        d_X0_vecs[part] = dynamic_cast<PetscVector<double>*>(d_X0_systems[part]->current_local_solution.get());

        d_U_systems[part] = &d_equation_systems[part]->get_system(VELOCITY_SYSTEM_NAME);
        d_U_current_vecs[part] = dynamic_cast<PetscVector<double>*>(d_U_systems[part]->current_local_solution.get());
        d_U_new_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
        d_U_half_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_current_vecs[part]->clone().release()); // WARNING: must be manually deleted

        d_F_systems[part] = &d_equation_systems[part]->get_system(FORCE_SYSTEM_NAME);
        d_F_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_F_systems[part]->current_local_solution.get());
        d_F_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(FORCE_SYSTEM_NAME, /*localize_data*/ false));

        d_F_t_systems[part] = &d_equation_systems[part]->get_system(FORCE_T_SYSTEM_NAME);
        d_F_t_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_F_t_systems[part]->current_local_solution.get());
        d_F_t_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(FORCE_T_SYSTEM_NAME, /*localize_data*/ false));
            
            
            
        d_F_b_systems[part] = &d_equation_systems[part]->get_system(FORCE_B_SYSTEM_NAME);
        d_F_b_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_F_b_systems[part]->current_local_solution.get());
        d_F_b_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(FORCE_B_SYSTEM_NAME, /*localize_data*/ false));

        d_F_n_systems[part] = &d_equation_systems[part]->get_system(FORCE_N_SYSTEM_NAME);
        d_F_n_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_F_n_systems[part]->current_local_solution.get());
        d_F_n_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(FORCE_N_SYSTEM_NAME, /*localize_data*/ false));
            
            
        d_H_systems[part] = &d_equation_systems[part]->get_system(H_SYSTEM_NAME);
        d_H_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_H_systems[part]->current_local_solution.get());
        d_H_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(H_SYSTEM_NAME, /*localize_data*/ false));


        d_du_j_systems[part] = &d_equation_systems[part]->get_system(DU_J_SYSTEM_NAME);
        d_du_j_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_du_j_systems[part]->current_local_solution.get());
        d_du_j_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(DU_J_SYSTEM_NAME, /*localize_data*/ false));

        d_dv_j_systems[part] = &d_equation_systems[part]->get_system(DV_J_SYSTEM_NAME);
        d_dv_j_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_dv_j_systems[part]->current_local_solution.get());
        d_dv_j_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(DV_J_SYSTEM_NAME, /*localize_data*/ false));

        d_dw_j_systems[part] = &d_equation_systems[part]->get_system(DW_J_SYSTEM_NAME);
        d_dw_j_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_dw_j_systems[part]->current_local_solution.get());
        d_dw_j_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(DW_J_SYSTEM_NAME, /*localize_data*/ false));

        d_d2u_j_systems[part] = &d_equation_systems[part]->get_system(D2U_J_SYSTEM_NAME);
        d_d2u_j_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_d2u_j_systems[part]->current_local_solution.get());
        d_d2u_j_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(D2U_J_SYSTEM_NAME, /*localize_data*/ false));

        d_d2v_j_systems[part] = &d_equation_systems[part]->get_system(D2V_J_SYSTEM_NAME);
        d_d2v_j_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_d2v_j_systems[part]->current_local_solution.get());
        d_d2v_j_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(D2V_J_SYSTEM_NAME, /*localize_data*/ false));
            
        d_d2w_j_systems[part] = &d_equation_systems[part]->get_system(D2W_J_SYSTEM_NAME);
        d_d2w_j_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_d2w_j_systems[part]->current_local_solution.get());
        d_d2w_j_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(D2W_J_SYSTEM_NAME, /*localize_data*/ false));
            
        d_dP_j_systems[part] = &d_equation_systems[part]->get_system(DP_J_SYSTEM_NAME);
        d_dP_j_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_dP_j_systems[part]->current_local_solution.get());
        d_dP_j_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(DP_J_SYSTEM_NAME, /*localize_data*/ false));
            
        d_P_j_systems[part] = &d_equation_systems[part]->get_system(P_J_SYSTEM_NAME);
        d_P_j_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_P_j_systems[part]->current_local_solution.get());
        d_P_j_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(P_J_SYSTEM_NAME, /*localize_data*/ false));

        d_WSS_i_systems[part] = &d_equation_systems[part]->get_system(WSS_I_SYSTEM_NAME);
        d_WSS_i_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_WSS_i_systems[part]->current_local_solution.get());
        d_WSS_i_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(WSS_I_SYSTEM_NAME, /*localize_data*/ false));

        d_WSS_o_systems[part] = &d_equation_systems[part]->get_system(WSS_O_SYSTEM_NAME);
        d_WSS_o_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_WSS_o_systems[part]->current_local_solution.get());
        d_WSS_o_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(WSS_O_SYSTEM_NAME, /*localize_data*/ false));

        // Initialize X^{n+1/2} and X^{n+1} to equal X^{n}, and initialize
        // U^{n+1/2} and U^{n+1} to equal U^{n}.
        d_X_systems[part]->solution->close();
        d_X_systems[part]->solution->localize(*d_X_current_vecs[part]);
        d_X_systems[part]->solution->localize(*d_X_new_vecs[part]);
        d_X_systems[part]->solution->localize(*d_X_half_vecs[part]);

        d_X0_systems[part]->solution->close();
        d_X0_systems[part]->solution->localize(*d_X0_vecs[part]);

        d_U_systems[part]->solution->close();
        d_U_systems[part]->solution->localize(*d_U_current_vecs[part]);
        d_U_systems[part]->solution->localize(*d_U_new_vecs[part]);
        d_U_systems[part]->solution->localize(*d_U_half_vecs[part]);

        d_F_systems[part]->solution->close();
        d_F_systems[part]->solution->localize(*d_F_half_vecs[part]);

        d_F_t_systems[part]->solution->close();
        d_F_t_systems[part]->solution->localize(*d_F_t_half_vecs[part]);
        
        
        d_F_n_systems[part]->solution->close();
        d_F_n_systems[part]->solution->localize(*d_F_n_half_vecs[part]);

        d_F_b_systems[part]->solution->close();
        d_F_b_systems[part]->solution->localize(*d_F_b_half_vecs[part]);
        
        d_H_systems[part]->solution->close();
        d_H_systems[part]->solution->localize(*d_H_half_vecs[part]);
        d_WSS_i_systems[part]->solution->close();
        d_WSS_i_systems[part]->solution->localize(*d_WSS_i_half_vecs[part]);

        d_WSS_o_systems[part]->solution->close();
        d_WSS_o_systems[part]->solution->localize(*d_WSS_o_half_vecs[part]);

        d_dP_j_systems[part]->solution->close();
        d_dP_j_systems[part]->solution->localize(*d_dP_j_half_vecs[part]);

        d_du_j_systems[part]->solution->close();
        d_du_j_systems[part]->solution->localize(*d_du_j_half_vecs[part]);

        d_dv_j_systems[part]->solution->close();
        d_dv_j_systems[part]->solution->localize(*d_dv_j_half_vecs[part]);

        d_dw_j_systems[part]->solution->close();
        d_dw_j_systems[part]->solution->localize(*d_dw_j_half_vecs[part]);

        d_d2u_j_systems[part]->solution->close();
        d_d2u_j_systems[part]->solution->localize(*d_d2u_j_half_vecs[part]);

        d_d2v_j_systems[part]->solution->close();
        d_d2v_j_systems[part]->solution->localize(*d_d2v_j_half_vecs[part]);

        d_d2w_j_systems[part]->solution->close();
        d_d2w_j_systems[part]->solution->localize(*d_d2w_j_half_vecs[part]);
        
    }

    // Update the mask data.
    getVelocityHierarchyDataOps()->copyData(mask_new_idx, mask_current_idx);
    return;
} // preprocessIntegrateData

void
IBFEMethod::postprocessIntegrateData(double /*current_time*/, double /*new_time*/, int /*num_cycles*/)
{
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        // Reset time-dependent Lagrangian data.
        d_X_new_vecs[part]->close();
        *d_X_systems[part]->solution = *d_X_new_vecs[part];
        d_X_systems[part]->solution->close();
        d_X_systems[part]->solution->localize(*d_X_systems[part]->current_local_solution);
        delete d_X_new_vecs[part];
        delete d_X_half_vecs[part];

        d_U_new_vecs[part]->close();
        *d_U_systems[part]->solution = *d_U_new_vecs[part];
        d_U_systems[part]->solution->close();
        d_U_systems[part]->solution->localize(*d_U_systems[part]->current_local_solution);
        delete d_U_new_vecs[part];
        delete d_U_half_vecs[part];

        d_F_half_vecs[part]->close();
        *d_F_systems[part]->solution = *d_F_half_vecs[part];
        d_F_systems[part]->solution->close();
        d_F_systems[part]->solution->localize(*d_F_systems[part]->current_local_solution);

        d_F_n_half_vecs[part]->close();
        *d_F_n_systems[part]->solution = *d_F_n_half_vecs[part];
        d_F_n_systems[part]->solution->close();
        d_F_n_systems[part]->solution->localize(*d_F_n_systems[part]->current_local_solution);

        d_F_t_half_vecs[part]->close();
        *d_F_t_systems[part]->solution = *d_F_t_half_vecs[part];
        d_F_t_systems[part]->solution->close();
        d_F_t_systems[part]->solution->localize(*d_F_t_systems[part]->current_local_solution);

        d_F_b_half_vecs[part]->close();
        *d_F_b_systems[part]->solution = *d_F_b_half_vecs[part];
        d_F_b_systems[part]->solution->close();
        d_F_b_systems[part]->solution->localize(*d_F_b_systems[part]->current_local_solution);

        d_P_j_half_vecs[part]->close();
        *d_P_j_systems[part]->solution = *d_P_j_half_vecs[part];
        d_P_j_systems[part]->solution->close();
        d_P_j_systems[part]->solution->localize(*d_P_j_systems[part]->current_local_solution);
        
        d_H_half_vecs[part]->close();
        *d_H_systems[part]->solution = *d_H_half_vecs[part];
        d_H_systems[part]->solution->close();
        d_H_systems[part]->solution->localize(*d_H_systems[part]->current_local_solution);


        d_dP_j_half_vecs[part]->close();
        *d_dP_j_systems[part]->solution = *d_dP_j_half_vecs[part];
        d_dP_j_systems[part]->solution->close();
        d_dP_j_systems[part]->solution->localize(*d_dP_j_systems[part]->current_local_solution);

        d_du_j_half_vecs[part]->close();
        *d_du_j_systems[part]->solution = *d_du_j_half_vecs[part];
        d_du_j_systems[part]->solution->close();
        d_du_j_systems[part]->solution->localize(*d_du_j_systems[part]->current_local_solution);

        d_dv_j_half_vecs[part]->close();
        *d_dv_j_systems[part]->solution = *d_dv_j_half_vecs[part];
        d_dv_j_systems[part]->solution->close();
        d_dv_j_systems[part]->solution->localize(*d_dv_j_systems[part]->current_local_solution);

        d_WSS_i_half_vecs[part]->close();
        *d_WSS_i_systems[part]->solution = *d_WSS_i_half_vecs[part];
        d_WSS_i_systems[part]->solution->close();
        d_WSS_i_systems[part]->solution->localize(*d_WSS_i_systems[part]->current_local_solution);

        d_WSS_o_half_vecs[part]->close();
        *d_WSS_o_systems[part]->solution = *d_WSS_o_half_vecs[part];
        d_WSS_o_systems[part]->solution->close();
        d_WSS_o_systems[part]->solution->localize(*d_WSS_o_systems[part]->current_local_solution);
        
        d_dw_j_half_vecs[part]->close();
        *d_dw_j_systems[part]->solution = *d_dw_j_half_vecs[part];
        d_dw_j_systems[part]->solution->close();
        d_dw_j_systems[part]->solution->localize(*d_dw_j_systems[part]->current_local_solution);


        d_d2u_j_half_vecs[part]->close();
        *d_d2u_j_systems[part]->solution = *d_d2u_j_half_vecs[part];
        d_d2u_j_systems[part]->solution->close();
        d_d2u_j_systems[part]->solution->localize(*d_d2u_j_systems[part]->current_local_solution);

        d_d2v_j_half_vecs[part]->close();
        *d_d2v_j_systems[part]->solution = *d_d2v_j_half_vecs[part];
        d_d2v_j_systems[part]->solution->close();
        d_d2v_j_systems[part]->solution->localize(*d_d2v_j_systems[part]->current_local_solution);
        

        d_d2w_j_half_vecs[part]->close();
        *d_d2w_j_systems[part]->solution = *d_d2w_j_half_vecs[part];
        d_d2w_j_systems[part]->solution->close();
        d_d2w_j_systems[part]->solution->localize(*d_d2w_j_systems[part]->current_local_solution);

        // Update the coordinate mapping dX = X - s.
        updateCoordinateMapping(part);
    }

    d_X_systems.clear();
    d_X_current_vecs.clear();
    d_X_new_vecs.clear();
    d_X_half_vecs.clear();
    d_X_IB_ghost_vecs.clear();

    d_X0_systems.clear();
    d_X0_vecs.clear();

    d_U_systems.clear();
    d_U_current_vecs.clear();
    d_U_new_vecs.clear();
    d_U_half_vecs.clear();

    d_F_systems.clear();
    d_F_half_vecs.clear();
    d_F_IB_ghost_vecs.clear();

    d_F_t_systems.clear();
    d_F_t_half_vecs.clear();
    d_F_t_IB_ghost_vecs.clear();

    d_F_n_systems.clear();
    d_F_n_half_vecs.clear();
    d_F_n_IB_ghost_vecs.clear();

    d_F_b_systems.clear();
    d_F_b_half_vecs.clear();
    d_F_b_IB_ghost_vecs.clear();
    
    d_H_systems.clear();
    d_H_half_vecs.clear();
    d_H_IB_ghost_vecs.clear();

    d_dP_j_systems.clear();
    d_dP_j_half_vecs.clear();
    d_dP_j_IB_ghost_vecs.clear();

   d_WSS_i_systems.clear();
    d_WSS_i_half_vecs.clear();
    d_WSS_i_IB_ghost_vecs.clear();

    d_WSS_o_systems.clear();
    d_WSS_o_half_vecs.clear();
    d_WSS_o_IB_ghost_vecs.clear();
    d_du_j_systems.clear();
    d_du_j_half_vecs.clear();
    d_du_j_IB_ghost_vecs.clear();

    d_dv_j_systems.clear();
    d_dv_j_half_vecs.clear();
    d_dv_j_IB_ghost_vecs.clear();

    d_dw_j_systems.clear();
    d_dw_j_half_vecs.clear();
    d_dw_j_IB_ghost_vecs.clear();

    d_d2u_j_systems.clear();
    d_d2u_j_half_vecs.clear();
    d_d2u_j_IB_ghost_vecs.clear();

    d_d2v_j_systems.clear();
    d_d2v_j_half_vecs.clear();
    d_d2v_j_IB_ghost_vecs.clear();

    d_d2w_j_systems.clear();
    d_d2w_j_half_vecs.clear();
    d_d2w_j_IB_ghost_vecs.clear();

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_half_time = std::numeric_limits<double>::quiet_NaN();
    return;
} // postprocessIntegrateData

//  interpolate velocity by including the jump condition 

//~ void
//~ IBFEMethod::interpolateVelocity(const int u_data_idx,
                                //~ const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
                                //~ const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                //~ const double data_time)
//~ {
//~ 
	//~ const int ndim = 2;
	//~ 
    //~ for (unsigned int part = 0; part < d_num_parts; ++part)
    //~ {
			//~ 
		//~ NumericVector<double>* du_j_ghost_vec = d_du_j_IB_ghost_vecs[part];
        //~ NumericVector<double>* dv_j_ghost_vec = d_dv_j_IB_ghost_vecs[part];
        //~ NumericVector<double>* X_vec = NULL;
        //~ NumericVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
        //~ NumericVector<double>* U_vec = NULL;
        //~ if (MathUtilities<double>::equalEps(data_time, d_current_time))
        //~ {
            //~ X_vec = d_X_current_vecs[part];
            //~ U_vec = d_U_current_vecs[part];
        //~ }
        //~ else if (MathUtilities<double>::equalEps(data_time, d_half_time))
        //~ {
            //~ X_vec = d_X_half_vecs[part];
            //~ U_vec = d_U_half_vecs[part];
        //~ }
        //~ else if (MathUtilities<double>::equalEps(data_time, d_new_time))
        //~ {
            //~ X_vec = d_X_new_vecs[part];
            //~ U_vec = d_U_new_vecs[part];
        //~ }
        //~ X_vec->localize(*X_ghost_vec);
    //~ 
		//~ // Extract the FE systems and DOF maps, and setup the FE object.
		//~ EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();	
		//~ const MeshBase& mesh = equation_systems->get_mesh();
		//~ const unsigned int dim = mesh.mesh_dimension();
		//~ AutoPtr<QBase> qrule;	
		//~ System& U_system = equation_systems->get_system(VELOCITY_SYSTEM_NAME);
		//~ const unsigned int n_vars = U_system.n_vars();
		//~ const DofMap& U_dof_map = U_system.get_dof_map();
		//~ FEDataManager::SystemDofMapCache& U_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(VELOCITY_SYSTEM_NAME);
		//~ FEType U_fe_type = U_dof_map.variable_type(0);
		//~ for (unsigned int d = 0; d < n_vars; ++d)
		//~ {
			//~ TBOX_ASSERT(U_dof_map.variable_type(d) == U_fe_type);
		//~ }
//~ 
		//~ std::vector<std::vector<unsigned int> > U_dof_indices(n_vars);
		//~ std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
		//~ 
		//~ System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
		//~ const DofMap& X_dof_map = X_system.get_dof_map();
		//~ FEDataManager::SystemDofMapCache& X_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(COORDS_SYSTEM_NAME);
//~ 
		//~ FEType X_fe_type = X_dof_map.variable_type(0);
//~ 
		//~ for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
		//~ 
		//~ AutoPtr<FEBase> U_fe_autoptr(FEBase::build(dim, U_fe_type)), X_fe_autoptr(NULL);
		//~ if (U_fe_type != X_fe_type)
		//~ {
			//~ X_fe_autoptr = AutoPtr<FEBase>(FEBase::build(dim, X_fe_type));
		//~ }
		//~ FEBase* U_fe = U_fe_autoptr.get();
		//~ FEBase* X_fe = X_fe_autoptr.get() ? X_fe_autoptr.get() : U_fe_autoptr.get();
		//~ const std::vector<double>& JxW_U = U_fe->get_JxW();
		//~ const std::vector<std::vector<double> >& phi_U = U_fe->get_phi();
		//~ const std::vector<std::vector<double> >& phi_X = X_fe->get_phi();
		//~ 
		//~ for (unsigned int k = 0; k < u_ghost_fill_scheds.size(); ++k)
		//~ {
			//~ if (u_ghost_fill_scheds[k]) u_ghost_fill_scheds[k]->fillData(data_time);
		//~ }
		//~ 
		//~ 
		//~ System& du_j_system = equation_systems->get_system(DU_J_SYSTEM_NAME);
		//~ const DofMap& du_j_dof_map = du_j_system.get_dof_map();
		//~ FEDataManager::SystemDofMapCache& du_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DU_J_SYSTEM_NAME);
		//~ FEType du_j_fe_type = du_j_dof_map.variable_type(0);
		//~ for (unsigned int d = 0; d < NDIM; ++d)
		//~ {
			 //~ TBOX_ASSERT(du_j_dof_map.variable_type(d) == du_j_fe_type);
		//~ }
		//~ std::vector<std::vector<unsigned int> > du_j_dof_indices(NDIM);
//~ 
		//~ System& dv_j_system = equation_systems->get_system(DV_J_SYSTEM_NAME);
		//~ const DofMap& dv_j_dof_map = dv_j_system.get_dof_map();
		//~ FEDataManager::SystemDofMapCache& dv_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DV_J_SYSTEM_NAME);    
		//~ FEType dv_j_fe_type = dv_j_dof_map.variable_type(0);
		//~ for (unsigned int d = 0; d < NDIM; ++d)
		//~ {
			 //~ TBOX_ASSERT(dv_j_dof_map.variable_type(d) == dv_j_fe_type);
		//~ }
		//~ std::vector<std::vector<unsigned int> > dv_j_dof_indices(NDIM);
//~ 
		//~ const std::vector<std::vector<Elem*> >& active_patch_element_map =
        //~ d_fe_data_managers[part]->getActivePatchElementMap();
		//~ const int level_num = d_fe_data_managers[part]->getLevelNumber();
		//~ //if (!X_vec.closed())/ 
		//~ (*X_ghost_vec).close();
		//~ 
		//~ PetscVector<double>* X_petsc_vec = static_cast<PetscVector<double>*>(X_ghost_vec);
		//~ Vec X_global_vec = X_petsc_vec->vec();
		//~ Vec X_local_vec;
		//~ VecGhostGetLocalForm(X_global_vec, &X_local_vec);
		//~ double* X_local_soln;
		//~ VecGetArray(X_local_vec, &X_local_soln);
		//~ AutoPtr<NumericVector<double> > U_rhs_vec = (*U_vec).zero_clone();
		//~ (*U_rhs_vec).zero();
		//~ DenseVector<double> U_rhs_e[n_vars];
		//~ boost::multi_array<double, 2> X_node;
		//~ boost::multi_array<double, 2> du_j_node, dv_j_node;
		//~ std::vector<double> U_qp, X_qp;
		//~ std::vector<double> du_j_qp, dv_j_qp;
		//~ du_j_ghost_vec->close();
		//~ PetscVector<double>* du_j_petsc_vec = static_cast<PetscVector<double>*>(du_j_ghost_vec);
		//~ Vec du_j_global_vec = du_j_petsc_vec->vec();
		//~ Vec du_j_local_vec;
		//~ VecGhostGetLocalForm(du_j_global_vec, &du_j_local_vec);
		//~ double* du_j_local_soln;
		//~ VecGetArray(du_j_local_vec, &du_j_local_soln);
		//~ dv_j_ghost_vec->close();
		//~ PetscVector<double>* dv_j_petsc_vec = static_cast<PetscVector<double>*>(dv_j_ghost_vec);
		//~ Vec dv_j_global_vec = dv_j_petsc_vec->vec();
		//~ Vec dv_j_local_vec;
		//~ VecGhostGetLocalForm(dv_j_global_vec, &dv_j_local_vec);
		//~ double* dv_j_local_soln;
		//~ VecGetArray(dv_j_local_vec, &dv_j_local_soln);
		//~ VectorValue<double> ju, jv;
		//~ Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_fe_data_managers[part]->getLevelNumber());
		//~ const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
		//~ const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift();
//~ 
        //~ 
        //~ int local_patch_num = 0;
		//~ for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
		//~ {
			//~ // The relevant collection of elements.
			//~ const std::vector<Elem*>& patch_elems = active_patch_element_map[local_patch_num];
			//~ const size_t num_active_patch_elems = patch_elems.size();
			//~ if (!num_active_patch_elems) continue;
			//~ const Pointer<Patch<NDIM> > patch = level->getPatch(p());
			//~ const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
			//~ const double* const patch_dx = patch_geom->getDx();
			//~ const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);
//~ 
			//~ // Setup vectors to store the values of U and X at the quadrature
			//~ // points.
			//~ //
			//~ // All this loop is doing is computing the total number of quadraturee
			//~ // points associated with all of the elements we are currently
			//~ // processing.  That number is n_qp_patch.
			//~ unsigned int n_qp_patch = 0;
			//~ for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
			//~ {
				//~ Elem* const elem = patch_elems[e_idx];
				//~ for (unsigned int d = 0; d < NDIM; ++d)
				//~ {
					//~ X_dof_map_cache.dof_indices(elem, X_dof_indices[d], d);
				//~ }
				//~ get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);
				//~ FEDataManager::updateInterpQuadratureRule(qrule, IBFEMethod::getDefaultInterpSpec(), elem, X_node, patch_dx_min);
				//~ n_qp_patch += qrule->n_points();
			//~ }			
			//~ 
			//~ if (!n_qp_patch) continue;
			//~ U_qp.resize(n_vars * n_qp_patch);
			//~ X_qp.resize(NDIM * n_qp_patch);
			//~ du_j_qp.resize(NDIM * n_qp_patch);
			//~ dv_j_qp.resize(NDIM * n_qp_patch);
			//~ std::fill(U_qp.begin(), U_qp.end(), 0.0);
//~ 
			//~ // Loop over the elements and compute the positions of the quadrature points.
			//~ qrule.reset();
			//~ unsigned int qp_offset = 0;
			//~ for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
			//~ {
				//~ Elem* const elem = patch_elems[e_idx];
				//~ for (unsigned int d = 0; d < NDIM; ++d)
				//~ {
					//~ X_dof_map_cache.dof_indices(elem, X_dof_indices[d], d);
					//~ du_j_dof_map_cache.dof_indices(elem, du_j_dof_indices[d], d);
					//~ dv_j_dof_map_cache.dof_indices(elem, dv_j_dof_indices[d], d);
				//~ }
				//~ get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);
				//~ get_values_for_interpolation(du_j_node, *du_j_petsc_vec, du_j_local_soln, du_j_dof_indices);
				//~ get_values_for_interpolation(dv_j_node, *dv_j_petsc_vec, dv_j_local_soln, dv_j_dof_indices);    
				//~ const bool qrule_changed = FEDataManager::updateInterpQuadratureRule(qrule, IBFEMethod::getDefaultInterpSpec(), elem, X_node, patch_dx_min);
				//~ if (qrule_changed)
				//~ {
					//~ // NOTE: Because we are only using the shape function values for
					//~ // the FE object associated with X, we only need to reinitialize
					//~ // X_fe whenever the quadrature rule changes.  In particular,
					//~ // notice that the shape function values depend only on the
					//~ // element type and quadrature rule, not on the element
					//~ // geometry.
					//~ X_fe->attach_quadrature_rule(qrule.get());
					//~ X_fe->reinit(elem);
				//~ }
				//~ const unsigned int n_node = elem->n_nodes();
				//~ const unsigned int n_qp = qrule->n_points();
				//~ 
				//~ // Zero out the values of X, du, and dv prior to accumulation.
				//~ double* X_begin = &X_qp[NDIM * qp_offset];
				//~ std::fill(X_begin, X_begin + NDIM * n_qp, 0.0);
				//~ double* du_j_begin = &du_j_qp[NDIM * qp_offset];
				//~ std::fill(du_j_begin, du_j_begin + NDIM * n_qp, 0.0);
				//~ double* dv_j_begin = &dv_j_qp[NDIM * qp_offset];
				//~ std::fill(dv_j_begin, dv_j_begin + NDIM * n_qp, 0.0);
							
				//~ // Interpolate X, du, and dv at all of the quadrature points
				//~ // via accumulation, i.e., X(qp) = sum_k X_k * phi_k(qp) for
				//~ // each qp.
				//~ for (unsigned int k = 0; k < n_node; ++k)
				//~ {
					//~ for (unsigned int qp = 0; qp < n_qp; ++qp)
					//~ {
						//~ const double& p_X = phi_X[k][qp];
						//~ for (unsigned int i = 0; i < NDIM; ++i)
						//~ {
							//~ X_qp [NDIM * (qp_offset + qp) + i] += X_node[k][i] * p_X;
							//~ du_j_qp[NDIM * (qp_offset + qp) + i] += du_j_node[k][i] * p_X;
							//~ dv_j_qp[NDIM * (qp_offset + qp) + i] += dv_j_node[k][i] * p_X;
						//~ }
					//~ }
				//~ }
				//~ qp_offset += n_qp;				
			//~ }
			//~ // Interpolate values from the Cartesian grid patch to the quadrature 
			//~ // points. 
			//~ // Note: Values are interpolated only to those quadrature points that 
			//~ // are within the patch interior
			//~ 
			//~ const Box<NDIM>& interp_box = patch->getBox();
			//~ Pointer<SideData<NDIM, double> >  u_sc_data = patch->getPatchData(u_data_idx);
			//~ Pointer<SideData<NDIM, double> >  mask_sc_data = patch->getPatchData(mask_current_idx);
//~ 
			//~ const IntVector<NDIM>& u_gcw = u_sc_data->getGhostCellWidth();
			//~ 
			//~ const int u_depth = u_sc_data->getDepth();
//~ 
			//~ const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
			//~ const double* const x_lower = pgeom->getXLower();
			//~ const double* const x_upper = pgeom->getXUpper();
			//~ const double* const dx = pgeom->getDx();
//~ 
		    //~ std::vector<int> local_indices;
			//~ local_indices.clear();
			//~ const int upper_bound = n_qp_patch;
			//~ if (upper_bound == 0) return;
//~ 
			//~ local_indices.reserve(upper_bound);
			//~ for (unsigned int k = 0; k < n_qp_patch; ++k)
			//~ {
				//~ const double* const XX = &X_qp[NDIM * k];
				//~ const Index<NDIM> i = IndexUtilities::getCellIndex(XX, patch_geom, interp_box);
				//~ if (interp_box.contains(i)) local_indices.push_back(k);
			//~ }
//~ 
//~ 
			//~ 
			//~ std::vector<double> periodic_shifts(NDIM * local_indices.size());
	//~ 
//~ 
		    //~ const int nindices = static_cast<int>(local_indices.size());
		    //~ 
		   //~ 
		    //~ typedef boost::multi_array_types::extent_range range;
//~ 
			//~ if (!local_indices.empty())
			//~ {
					//~ boost::array<int, 2> ic_trimmed_lower, ic_trimmed_upper, ic_lower, ic_upper, ic_center;
					//~ boost::array<double, 2> X_shifted, X_cell;
					//~ boost::array<double, 2> w0, w1;
					//~ boost::array<double, NDIM> x_lower_axis, x_upper_axis;
					//~ const int local_sz = (*std::max_element(local_indices.begin(), local_indices.end())) + 1;
					//~ std::vector<double> Q_data_axis(local_sz);
					//~ 
					//~ x_lower_axis[0] = x_lower_axis[1] = x_upper_axis[0] = x_upper_axis[1] = 0.0;
					//~ 
					//~ Box<NDIM> side_boxes[NDIM];
					//~ 
					//~ for (unsigned int axis = 0; axis < NDIM; ++axis)
					//~ {
						//~ 
						//~ for (unsigned int d = 0; d < NDIM; ++d)
						//~ {
							//~ x_lower_axis[d] = x_lower[d];
							//~ x_upper_axis[d] = x_upper[d];
						//~ }
						//~ x_lower_axis[axis] -= 0.5*dx[axis];
						//~ x_upper_axis[axis] += 0.5*dx[axis];
	//~ 
						//~ for (int d = 0; d < NDIM; ++d)
						//~ {
							//~ side_boxes[d] = SideGeometry<NDIM>::toSideBox(interp_box, d);
						//~ }
					//~ 
						//~ const IntVector<NDIM>& ilower = side_boxes[axis].lower();
						//~ const IntVector<NDIM>& iupper = side_boxes[axis].upper();	
#if (NDIM == 2)				
						//~ boost::const_multi_array_ref<double, ndim + 1 > u_sc_data_array(u_sc_data->getPointer(axis),(boost::extents[range(ilower[0] - u_gcw[0], iupper[0] + u_gcw[0] + 1)][range(ilower[1] - u_gcw[1],
                                                                                     //~ iupper[1] + u_gcw[1] + 1)][range(0, u_depth)]),boost::fortran_storage_order());
                                                                                     //~ 
                                                                                     //~ 
  						//~ boost::const_multi_array_ref<double, ndim + 1 > mask_sc_data_array(mask_sc_data->getPointer(axis),(boost::extents[range(ilower[0] - u_gcw[0], iupper[0] + u_gcw[0] + 1)][range(ilower[1] - u_gcw[1],
                                                                                     //~ iupper[1] + u_gcw[1] + 1)][range(0, u_depth)]),boost::fortran_storage_order());                                                                        
		//~ 
#endif		
					//~ 
					//~ 
						//~ for (unsigned int k = 0; k < nindices; ++k)
						//~ {
							//~ 
									//~ const int s = local_indices[k];
					//~ 
									//~ 
									//~ for (int d = 0; d < ndim; ++d) 
									//~ {
										//~ X_shifted[d] = X_qp[d + s * ndim] + periodic_shifts[d + k * ndim];
									//~ }
									//~ 
//~ 
								//~ 
									//~ for (unsigned int d = 0; d < ndim; ++d) 
									//~ {
										//~ ic_center[d] = ilower[d] + NINT((X_shifted[d]-x_lower_axis[d])/dx[d]-0.5);
										//~ X_cell[d] = x_lower_axis[d] + (static_cast<double>(ic_center[d]-ilower[d]) + 0.5 )*dx[d];
									//~ 
									//~ 
										//~ 
										//~ if ( X_shifted[d] <= X_cell[d] )
										//~ {
										   //~ ic_lower[d] = ic_center[d]-1;
										   //~ ic_upper[d] = ic_center[d];
//~ 
										//~ }
										//~ else
										//~ {
										   //~ ic_lower[d] = ic_center[d];
										   //~ ic_upper[d] = ic_center[d]+1;
//~ 
										//~ }
										//~ ic_trimmed_lower[d] = std::max(ic_lower[d],ilower[d]- u_gcw[d]);
										//~ ic_trimmed_upper[d] = std::min(ic_upper[d],iupper[d]+ u_gcw[d]);
									//~ }
								//~ 
//~ 
								//~ 
								//~ 
								//~ 
									//~ if ( X_shifted[0] <= X_cell[0] )
									//~ {
									   //~ w0[0] = (X_cell[0]-X_shifted[0])/dx[0];
									   //~ 
									   //~ 
									   //~ w0[1] = 1.0 - w0[0];
									//~ }
									//~ else
									//~ {
									   //~ w0[0] = 1.0 + (X_cell[0]-X_shifted[0])/dx[0];
									   //~ w0[1] = 1.0 - w0[0];
									//~ }
									//~ 
									//~ 
									//~ if ( X_shifted[1] <= X_cell[1] )
									//~ {
									   //~ w1[0] = (X_cell[1]-X_shifted[1])/dx[1];
									   //~ w1[1] = 1.0 - w1[0];
									//~ }
									//~ else
									//~ {
									   //~ w1[0] = 1.0 + (X_cell[1]-X_shifted[1])/dx[1];
									   //~ w1[1] = 1.0 - w1[0];
									//~ }
									//~ 
//~ 
								//~ 
									//~ boost::multi_array<double, 2> ujump (boost::extents[range(ic_trimmed_lower[0], ic_trimmed_upper[0] + 1)][range(ic_trimmed_lower[1], ic_trimmed_upper[1] + 1)]);
							//~ 
							//~ 
									//~ boost::multi_array<double, 2> vjump (boost::extents[range(ic_trimmed_lower[0], ic_trimmed_upper[0] + 1)][range(ic_trimmed_lower[1], ic_trimmed_upper[1] + 1)]);
									//~ 
									//~ 
									//~ 
									//~ 
									 //~ ujump[ic_trimmed_lower[0]][ic_trimmed_lower[1]] = dx[0]*w0[0]*w1[0]*(w0[1]*du_j_qp[s * ndim] + w1[1]*du_j_qp[1 + s * ndim]);
									 //~ vjump[ic_trimmed_lower[0]][ic_trimmed_lower[1]] = dx[1]*w0[0]*w1[0]*(w0[1]*dv_j_qp[s * ndim] + w1[1]*dv_j_qp[1 + s * ndim]);
									 //~ 
									 //~ ujump[ic_trimmed_upper[0]][ic_trimmed_lower[1]] = -dx[0]*w0[1]*w1[0]*(w0[0]*du_j_qp[s * ndim] - w1[1]*du_j_qp[1 + s * ndim]);
									 //~ vjump[ic_trimmed_upper[0]][ic_trimmed_lower[1]] = -dx[1]*w0[1]*w1[0]*(w0[0]*dv_j_qp[s * ndim] - w1[1]*dv_j_qp[1 + s * ndim]);
//~ 
									 //~ ujump[ic_trimmed_upper[0]][ic_trimmed_upper[1]] = -dx[0]*w0[1]*w1[1]*(w0[0]*du_j_qp[s * ndim] + w1[0]*du_j_qp[1 + s * ndim]);
									 //~ vjump[ic_trimmed_upper[0]][ic_trimmed_upper[1]] = -dx[1]*w0[1]*w1[1]*(w0[0]*dv_j_qp[s * ndim] + w1[0]*dv_j_qp[1 + s * ndim]);
//~ 
									 //~ ujump[ic_trimmed_lower[0]][ic_trimmed_upper[1]] =  dx[0]*w0[0]*w1[1]*(w0[1]*du_j_qp[s * ndim] - w1[0]*du_j_qp[1 + s * ndim]);
									 //~ vjump[ic_trimmed_lower[0]][ic_trimmed_upper[1]] =  dx[1]*w0[0]*w1[1]*(w0[1]*dv_j_qp[s * ndim] - w1[0]*dv_j_qp[1 + s * ndim]);
//~ 
									//~ 
									//~ 
									//~ // ujump[ic_trimmed_lower[0]][ic_trimmed_lower[1]] = dx[0]*w0[0]*w1[0]*(w0[1]*du_qp[s * ndim] + w1[1]*du_qp[1 + s * ndim]);
									//~ // vjump[ic_trimmed_lower[0]][ic_trimmed_lower[1]] = dx[1]*w0[0]*w1[0]*(w0[1]*dv_qp[s * ndim] + w1[1]*dv_qp[1 + s * ndim]);
									//~ // 
									//~ // ujump[ic_trimmed_upper[0]][ic_trimmed_lower[1]] = -dx[0]*w0[1]*w1[0]*(w0[0]*du_qp[s * ndim] - w1[1]*du_qp[1 + s * ndim]);
									//~ // vjump[ic_trimmed_upper[0]][ic_trimmed_lower[1]] = -dx[1]*w0[1]*w1[0]*(w0[0]*dv_qp[s * ndim] - w1[1]*dv_qp[1 + s * ndim]);
																		
//~ 
									//~ // ujump[ic_trimmed_lower[0]][ic_trimmed_upper[1]] =  dx[0]*w0[0]*w1[1]*(w0[1]*du_qp[s * ndim] - w1[0]*du_qp[1 + s * ndim]);
									//~ // vjump[ic_trimmed_lower[0]][ic_trimmed_upper[1]] = dx[1]*w0[0]*w1[1]*(w0[1]*dv_qp[s * ndim] - w1[0]*dv_qp[1 + s * ndim]);
																		
									//~ // ujump[ic_trimmed_upper[0]][ic_trimmed_upper[1]] = -dx[0]*w0[1]*w1[1]*(w0[0]*du_qp[s * ndim] + w1[0]*du_qp[1 + s * ndim]);
									//~ // vjump[ic_trimmed_upper[0]][ic_trimmed_upper[1]] = -dx[1]*w0[1]*w1[1]*(w0[0]*dv_qp[s * ndim] + w1[0]*dv_qp[1 + s * ndim]);
//~ 
								//~ 
//~ 
								//~ 
									//~ double CC= 0.0;
									//~ 
									//~ for (int d = 0; d < u_depth; ++d) 
									//~ {
											//~ Q_data_axis[s ] = 0.0;
											//~ for (int ic1 = ic_trimmed_lower[1]; ic1 <= ic_trimmed_upper[1]; ++ic1)
											//~ {
											   //~ for (int ic0 = ic_trimmed_lower[0]; ic0 <= ic_trimmed_upper[0]; ++ic0)
											   //~ {
												   //~ if (d_modify_vel_interp_jumps)
												   //~ {
														 //~ if (axis == 0  && fabs(mask_sc_data_array[ic0][ic1][d] - 1.0)< 0.0001)
														 //~ {
																//~ CC = ujump[ic0][ic1];
//~ 
														 //~ }
														 //~ else if (axis == 1 && fabs(mask_sc_data_array[ic0][ic1][d] - 1.0)< 0.0001)
														 //~ {
															  //~ CC = vjump[ic0][ic1];	
														 //~ }
												  //~ }
													//~ 
												  //~ Q_data_axis[s] = Q_data_axis[ s ] + w0[ic0 - ic_lower[0]]*w1[ ic1 - ic_lower[1]]*u_sc_data_array[ic0][ic1][d] + CC/d_mu ;
												//~ }
											//~ }
								   //~ } //depth
							//~ 
							//~ 
							//~ 
							//~ 
							//~ 
							//~ 
							//~ 
							//~ }
						//~ 
						//~ 
						//~ 
						//~ 
						//~ 
						//~ 
						//~ for (unsigned int k = 0; k < nindices; ++k)
						//~ {
							//~ U_qp[n_vars*local_indices[k] + axis] = Q_data_axis[local_indices[k]];
						//~ }
//~ 
					//~ }		 
				//~ 
			//~ }
			//~ 
		//~ 
				//~ 
		//~ // Loop over the elements and accumulate the right-hand-side values.
			//~ qrule.reset();
			//~ qp_offset = 0;
			//~ for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
			//~ {
				//~ Elem* const elem = patch_elems[e_idx];
				//~ for (unsigned int i = 0; i < n_vars; ++i)
				//~ {
					//~ U_dof_map_cache.dof_indices(elem, U_dof_indices[i], i);
					//~ U_rhs_e[i].resize(static_cast<int>(U_dof_indices[i].size()));
				//~ }
				//~ for (unsigned int d = 0; d < ndim; ++d)
				//~ {
					//~ X_dof_map_cache.dof_indices(elem, X_dof_indices[d], d);
				//~ }
				//~ get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);
				//~ const bool qrule_changed = FEDataManager::updateInterpQuadratureRule(qrule, IBFEMethod::getDefaultInterpSpec(), elem, X_node, patch_dx_min);
				//~ if (qrule_changed)
				//~ {
					//~ // NOTE: Because we are only using the shape function values for
					//~ // the FE object associated with X, we only need to reinitialize
					//~ // X_fe whenever the quadrature rule changes.  In particular,
					//~ // notice that the shape function values depend only on the
					//~ // element type and quadrature rule, not on the element
					//~ // geometry.
					//~ U_fe->attach_quadrature_rule(qrule.get());
					//~ X_fe->attach_quadrature_rule(qrule.get());
					//~ if (X_fe != U_fe) X_fe->reinit(elem);
				//~ }
				//~ U_fe->reinit(elem);
				//~ const unsigned int n_qp = qrule->n_points();
				//~ const size_t n_basis = U_dof_indices[0].size();
				//~ for (unsigned int qp = 0; qp < n_qp; ++qp)
				//~ {
					//~ const int idx = n_vars * (qp_offset + qp);
					//~ for (unsigned int k = 0; k < n_basis; ++k)
					//~ {
						//~ const double p_JxW_U = phi_U[k][qp] * JxW_U[qp];
						//~ for (unsigned int i = 0; i < n_vars; ++i)
						//~ {
							//~ U_rhs_e[i](k) += U_qp[idx + i] * p_JxW_U;
						//~ }
					//~ }
				//~ }
				//~ for (unsigned int i = 0; i < n_vars; ++i)
				//~ {
					//~ U_dof_map.constrain_element_vector(U_rhs_e[i], U_dof_indices[i]);
					//~ U_rhs_vec->add_vector(U_rhs_e[i], U_dof_indices[i]);
				//~ }
				//~ qp_offset += n_qp;
           //~ }
			//~ 
		//~ 
		//~ }
	//~ 
        //~ 
		//~ U_rhs_vec->close();
//~ 
		//~ VecRestoreArray(X_local_vec, &X_local_soln);
		//~ VecGhostRestoreLocalForm(X_global_vec, &X_local_vec);
		//~ VecRestoreArray(du_j_local_vec, &du_j_local_soln);
		//~ VecGhostRestoreLocalForm(du_j_global_vec, &du_j_local_vec);
		//~ VecRestoreArray(dv_j_local_vec, &dv_j_local_soln);
		//~ VecGhostRestoreLocalForm(dv_j_global_vec, &dv_j_local_vec);
		//~ d_fe_data_managers[part]->computeL2Projection(*U_vec, *U_rhs_vec, VELOCITY_SYSTEM_NAME,  IBFEMethod::getDefaultInterpSpec().use_consistent_mass_matrix);
//~ 
//~ 
	    //~ d_du_j_IB_ghost_vecs[part]->close();
        //~ d_dv_j_IB_ghost_vecs[part]->close();
        //~ d_X_half_vecs[part]->close();
        //~ d_X_current_vecs[part]->close();
        //~ d_X_new_vecs[part]->close();
        //~ d_U_new_vecs[part]->close();
	//~ }
        //~ 
   //~ 
 //~ 
    //~ return;
//~ 
//~ } // interpolateVelocityWithJump







void
IBFEMethod::interpolateVelocity(const int u_data_idx,
                                const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
                                const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                const double data_time)
{

    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
		NumericVector<double>* WSS_i_vec = d_WSS_i_half_vecs[part];
        NumericVector<double>* WSS_o_vec = d_WSS_o_half_vecs[part];
        NumericVector<double>* du_j_ghost_vec = d_du_j_IB_ghost_vecs[part];
        NumericVector<double>* dv_j_ghost_vec = d_dv_j_IB_ghost_vecs[part];
        NumericVector<double>* dw_j_ghost_vec = d_dw_j_IB_ghost_vecs[part];
        NumericVector<double>* X_vec = NULL;
        NumericVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
        NumericVector<double>* U_vec = NULL;
        AutoPtr<NumericVector<double> > WSS_i_rhs_vec = (*WSS_i_vec).zero_clone();
        (*WSS_i_rhs_vec).zero();
		DenseVector<double> WSS_i_rhs_e[NDIM];
		
	    AutoPtr<NumericVector<double> > WSS_o_rhs_vec = (*WSS_o_vec).zero_clone();
	    (*WSS_o_rhs_vec).zero();
		DenseVector<double> WSS_o_rhs_e[NDIM];

        if (MathUtilities<double>::equalEps(data_time, d_current_time))
        {
            X_vec = d_X_current_vecs[part];
            U_vec = d_U_current_vecs[part];
        }
        else if (MathUtilities<double>::equalEps(data_time, d_half_time))
        {
            X_vec = d_X_half_vecs[part];
            U_vec = d_U_half_vecs[part];
        }
        else if (MathUtilities<double>::equalEps(data_time, d_new_time))
        {
            X_vec = d_X_new_vecs[part];
            U_vec = d_U_new_vecs[part];
        }
        X_vec->localize(*X_ghost_vec);

        // Extract the FE systems and DOF maps, and setup the FE object.
        EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
        const MeshBase& mesh = equation_systems->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();
        AutoPtr<QBase> qrule;
        System& U_system = equation_systems->get_system(VELOCITY_SYSTEM_NAME);
        const unsigned int n_vars = U_system.n_vars();
        const DofMap& U_dof_map = U_system.get_dof_map();
        FEDataManager::SystemDofMapCache& U_dof_map_cache =
            *d_fe_data_managers[part]->getDofMapCache(VELOCITY_SYSTEM_NAME);
        FEType U_fe_type = U_dof_map.variable_type(0);
        for (unsigned int d = 0; d < n_vars; ++d)
        {
            TBOX_ASSERT(U_dof_map.variable_type(d) == U_fe_type);
        }

        std::vector<std::vector<unsigned int> > U_dof_indices(n_vars);
        std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);

        System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
        const DofMap& X_dof_map = X_system.get_dof_map();
        FEDataManager::SystemDofMapCache& X_dof_map_cache =
            *d_fe_data_managers[part]->getDofMapCache(COORDS_SYSTEM_NAME);

        FEType X_fe_type = X_dof_map.variable_type(0);

        for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);

        AutoPtr<FEBase> U_fe_autoptr(FEBase::build(dim, U_fe_type)), X_fe_autoptr(NULL);
        if (U_fe_type != X_fe_type)
        {
            X_fe_autoptr = AutoPtr<FEBase>(FEBase::build(dim, X_fe_type));
        }
        FEBase* U_fe = U_fe_autoptr.get();
        FEBase* X_fe = X_fe_autoptr.get() ? X_fe_autoptr.get() : U_fe_autoptr.get();
        const std::vector<double>& JxW_U = U_fe->get_JxW();
        const std::vector<std::vector<double> >& phi_U = U_fe->get_phi();
        const std::vector<std::vector<double> >& phi_X = X_fe->get_phi();
        const std::vector<std::vector<double> >& X_dphi_dxi = X_fe->get_dphidxi();
        const std::vector<std::vector<double> >& X_dphi_deta = X_fe->get_dphideta();

        VectorValue<double> tau1, tau2, n;

        for (unsigned int k = 0; k < u_ghost_fill_scheds.size(); ++k)
        {
            if (u_ghost_fill_scheds[k]) u_ghost_fill_scheds[k]->fillData(data_time);
        }

        System& WSS_i_system = equation_systems->get_system(WSS_I_SYSTEM_NAME);
        const DofMap& WSS_i_dof_map = WSS_i_system.get_dof_map();
        FEDataManager::SystemDofMapCache& WSS_i_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(WSS_I_SYSTEM_NAME);
        FEType WSS_i_fe_type = WSS_i_dof_map.variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(WSS_i_dof_map.variable_type(d) == WSS_i_fe_type);
        }
        std::vector<std::vector<unsigned int> > WSS_i_dof_indices(NDIM);

        System& WSS_o_system = equation_systems->get_system(WSS_O_SYSTEM_NAME);
        const DofMap& WSS_o_dof_map = WSS_o_system.get_dof_map();
        FEDataManager::SystemDofMapCache& WSS_o_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(WSS_O_SYSTEM_NAME);
        FEType WSS_o_fe_type = WSS_o_dof_map.variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(WSS_o_dof_map.variable_type(d) == WSS_o_fe_type);
        }
        std::vector<std::vector<unsigned int> > WSS_o_dof_indices(NDIM);

        System& du_j_system = equation_systems->get_system(DU_J_SYSTEM_NAME);
        const DofMap& du_j_dof_map = du_j_system.get_dof_map();
        FEDataManager::SystemDofMapCache& du_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DU_J_SYSTEM_NAME);
        FEType du_j_fe_type = du_j_dof_map.variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(du_j_dof_map.variable_type(d) == du_j_fe_type);
        }
        std::vector<std::vector<unsigned int> > du_j_dof_indices(NDIM);

        System& dv_j_system = equation_systems->get_system(DV_J_SYSTEM_NAME);
        const DofMap& dv_j_dof_map = dv_j_system.get_dof_map();
        FEDataManager::SystemDofMapCache& dv_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DV_J_SYSTEM_NAME);
        FEType dv_j_fe_type = dv_j_dof_map.variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(dv_j_dof_map.variable_type(d) == dv_j_fe_type);
        }
        std::vector<std::vector<unsigned int> > dv_j_dof_indices(NDIM);

#if (NDIM == 3)
        System& dw_j_system = equation_systems->get_system(DW_J_SYSTEM_NAME);
        const DofMap& dw_j_dof_map = dw_j_system.get_dof_map();
        FEDataManager::SystemDofMapCache& dw_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DW_J_SYSTEM_NAME);
        FEType dw_j_fe_type = dw_j_dof_map.variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(dw_j_dof_map.variable_type(d) == dw_j_fe_type);
        }
        std::vector<std::vector<unsigned int> > dw_j_dof_indices(NDIM);
#endif

        const std::vector<std::vector<Elem*> >& active_patch_element_map =
            d_fe_data_managers[part]->getActivePatchElementMap();
        const int level_num = d_fe_data_managers[part]->getLevelNumber();
        // if (!X_vec.closed())/
        (*X_ghost_vec).close();

        PetscVector<double>* X_petsc_vec = static_cast<PetscVector<double>*>(X_ghost_vec);
        Vec X_global_vec = X_petsc_vec->vec();
        Vec X_local_vec;
        VecGhostGetLocalForm(X_global_vec, &X_local_vec);
        double* X_local_soln;
        VecGetArray(X_local_vec, &X_local_soln);
        AutoPtr<NumericVector<double> > U_rhs_vec = (*U_vec).zero_clone();
        (*U_rhs_vec).zero();
        std::vector<DenseVector<double> > U_rhs_e(n_vars);
        boost::multi_array<double, 2> X_node;
        boost::multi_array<double, 2> du_j_node, dv_j_node, dw_j_node;
        std::vector<double> WSS_i_qp, WSS_o_qp, U_qp, N_qp, X_qp_p, X_qp_m;
        std::vector<double> du_j_qp, dv_j_qp, dw_j_qp;

        du_j_ghost_vec->close();
        PetscVector<double>* du_j_petsc_vec = static_cast<PetscVector<double>*>(du_j_ghost_vec);
        Vec du_j_global_vec = du_j_petsc_vec->vec();
        Vec du_j_local_vec;
        VecGhostGetLocalForm(du_j_global_vec, &du_j_local_vec);
        double* du_j_local_soln;
        VecGetArray(du_j_local_vec, &du_j_local_soln);

        dv_j_ghost_vec->close();
        PetscVector<double>* dv_j_petsc_vec = static_cast<PetscVector<double>*>(dv_j_ghost_vec);
        Vec dv_j_global_vec = dv_j_petsc_vec->vec();
        Vec dv_j_local_vec;
        VecGhostGetLocalForm(dv_j_global_vec, &dv_j_local_vec);
        double* dv_j_local_soln;
        VecGetArray(dv_j_local_vec, &dv_j_local_soln);

        dw_j_ghost_vec->close();
        PetscVector<double>* dw_j_petsc_vec = static_cast<PetscVector<double>*>(dw_j_ghost_vec);
        Vec dw_j_global_vec = dw_j_petsc_vec->vec();
        Vec dw_j_local_vec;
        VecGhostGetLocalForm(dw_j_global_vec, &dw_j_local_vec);
        double* dw_j_local_soln;
        VecGetArray(dw_j_local_vec, &dw_j_local_soln);

        VectorValue<double> ju, jv, jw;
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_fe_data_managers[part]->getLevelNumber());
        const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
        const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift();

        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            // The relevant collection of elements.
            const std::vector<Elem*>& patch_elems = active_patch_element_map[local_patch_num];
            const size_t num_active_patch_elems = patch_elems.size();
            if (!num_active_patch_elems) continue;
            const Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();
            const double dx_min = *std::min_element(dx, dx + NDIM);

            const Box<NDIM>& interp_box = patch->getBox();
            Pointer<SideData<NDIM, double> > u_sc_data = patch->getPatchData(u_data_idx);
            const IntVector<NDIM>& u_gcw = u_sc_data->getGhostCellWidth();
            const int u_depth = u_sc_data->getDepth();
            const double* const x_lower = patch_geom->getXLower();
            const double* const x_upper = patch_geom->getXUpper();

            double diag_dis = 0.0;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                diag_dis += dx[d] * dx[d];
            }
            const double dh = d_vel_interp_width * sqrt(diag_dis);

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
                    X_dof_map_cache.dof_indices(elem, X_dof_indices[d], d);
                }
                get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);
                FEDataManager::updateInterpQuadratureRule(
                    qrule, IBFEMethod::getDefaultInterpSpec(), elem, X_node, dx_min);
                n_qp_patch += qrule->n_points();
            }

            if (!n_qp_patch) continue;
            WSS_i_qp.resize(n_vars * n_qp_patch);
            WSS_o_qp.resize(n_vars * n_qp_patch);
            U_qp.resize(n_vars * n_qp_patch);
            X_qp_m.resize(NDIM * n_qp_patch);
            X_qp_p.resize(NDIM * n_qp_patch);
            du_j_qp.resize(NDIM * n_qp_patch);
            dv_j_qp.resize(NDIM * n_qp_patch);
            dw_j_qp.resize(NDIM * n_qp_patch);
            N_qp.resize(NDIM * n_qp_patch);
            std::fill(U_qp.begin(), U_qp.end(), 0.0);
            std::fill(N_qp.begin(), N_qp.end(), 0.0);
        	std::fill(WSS_i_qp.begin(), WSS_i_qp.end(), 0.0);
            std::fill(WSS_o_qp.begin(), WSS_o_qp.end(), 0.0);

            // Loop over the elements and compute the positions of the quadrature points.
            qrule.reset();
            unsigned int qp_offset = 0;

            for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems[e_idx];
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X_dof_map_cache.dof_indices(elem, X_dof_indices[d], d);
                    du_j_dof_map_cache.dof_indices(elem, du_j_dof_indices[d], d);
                    dv_j_dof_map_cache.dof_indices(elem, dv_j_dof_indices[d], d);
#if (NDIM == 3)
                    dw_j_dof_map_cache.dof_indices(elem, dw_j_dof_indices[d], d);
#endif
                }
                get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);
                get_values_for_interpolation(du_j_node, *du_j_petsc_vec, du_j_local_soln, du_j_dof_indices);
                get_values_for_interpolation(dv_j_node, *dv_j_petsc_vec, dv_j_local_soln, dv_j_dof_indices);
#if (NDIM == 3)
                get_values_for_interpolation(dw_j_node, *dw_j_petsc_vec, dw_j_local_soln, dw_j_dof_indices);
#endif

                const bool qrule_changed = FEDataManager::updateInterpQuadratureRule(
                    qrule, IBFEMethod::getDefaultInterpSpec(), elem, X_node, dx_min);
                if (qrule_changed)
                {
                    // NOTE: Because we are only using the shape function values for
                    // the FE object associated with X, we only need to reinitialize
                    // X_fe whenever the quadrature rule changes.  In particular,
                    // notice that the shape function values depend only on the
                    // element type and quadrature rule, not on the element
                    // geometry.
                    X_fe->attach_quadrature_rule(qrule.get());
                    X_fe->reinit(elem);
                }
                const unsigned int n_node = elem->n_nodes();
                const unsigned int n_qp = qrule->n_points();
                //~
                // Zero out the values of X, du, and dv prior to accumulation.
                double* X_begin_m = &X_qp_m[NDIM * qp_offset];
                std::fill(X_begin_m, X_begin_m + NDIM * n_qp, 0.0);

                double* X_begin_p = &X_qp_p[NDIM * qp_offset];
                std::fill(X_begin_p, X_begin_p + NDIM * n_qp, 0.0);

                double* du_j_begin = &du_j_qp[NDIM * qp_offset];
                std::fill(du_j_begin, du_j_begin + NDIM * n_qp, 0.0);
                double* dv_j_begin = &dv_j_qp[NDIM * qp_offset];
                std::fill(dv_j_begin, dv_j_begin + NDIM * n_qp, 0.0);
#if (NDIM == 3)
                double* dw_j_begin = &dw_j_qp[NDIM * qp_offset];
                std::fill(dw_j_begin, dw_j_begin + NDIM * n_qp, 0.0);
#endif
                double* N_begin = &N_qp[NDIM * qp_offset];
                std::fill(N_begin, N_begin + NDIM * n_qp, 0.0);

                // Interpolate X, du, and dv at all of the quadrature points
                // via accumulation, i.e., X(qp) = sum_k X_k * phi_k(qp) for
                // each qp.

                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
                    if (dim == 1)
                        tau2 = VectorValue<double>(0.0, 0.0, 1.0);
                    else
                        interpolate(&tau2(0), qp, X_node, X_dphi_deta);
                    n = tau1.cross(tau2);
                    n = n.unit();

                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        for (unsigned int k = 0; k < n_node; ++k)
                        {
                            const double& p_X = phi_X[k][qp];

                            X_qp_p[NDIM * (qp_offset + qp) + i] += X_node[k][i] * p_X;
                            X_qp_m[NDIM * (qp_offset + qp) + i] += X_node[k][i] * p_X;
                            du_j_qp[NDIM * (qp_offset + qp) + i] += du_j_node[k][i] * p_X;
                            dv_j_qp[NDIM * (qp_offset + qp) + i] += dv_j_node[k][i] * p_X;
#if (NDIM == 3)
                            dw_j_qp[NDIM * (qp_offset + qp) + i] += dw_j_node[k][i] * p_X;
#endif
                        }
                        N_qp[NDIM * (qp_offset + qp) + i] = n(i);
                        X_qp_p[NDIM * (qp_offset + qp) + i] += -n(i) * dh;
                        X_qp_m[NDIM * (qp_offset + qp) + i] += n(i) * dh;
                    }
                }
                qp_offset += n_qp;
            }
            // Interpolate values from the Cartesian grid patch to the quadrature
            // points.
            // Note: Values are interpolated only to those quadrature points that
            // are within the patch interior
            //

            std::vector<int> local_indices;
            local_indices.clear();
            const int upper_bound = n_qp_patch;
            if (upper_bound == 0) return;

            local_indices.reserve(upper_bound);
            for (unsigned int k = 0; k < n_qp_patch; ++k)
            {
                const double* const XX_p = &X_qp_p[NDIM * k];
                const double* const XX_m = &X_qp_m[NDIM * k];

                const Index<NDIM> i_p = IndexUtilities::getCellIndex(XX_p, patch_geom, interp_box);
                const Index<NDIM> i_m = IndexUtilities::getCellIndex(XX_m, patch_geom, interp_box);
                if (interp_box.contains(i_p) && interp_box.contains(i_m)) local_indices.push_back(k);
            }

            std::vector<double> periodic_shifts(NDIM * local_indices.size());

            const int nindices = static_cast<int>(local_indices.size());

            typedef boost::multi_array_types::extent_range range;

            if (!local_indices.empty())
            {
                boost::array<int, NDIM> ic_trimmed_lower_p, ic_trimmed_upper_p, ic_lower_p, ic_upper_p, ic_center_p;
                boost::array<int, NDIM> ic_trimmed_lower_m, ic_trimmed_upper_m, ic_lower_m, ic_upper_m, ic_center_m;
                boost::array<double, NDIM> X_shifted_p, X_cell_p;
                boost::array<double, NDIM> X_shifted_m, X_cell_m;
                boost::array<double, 2> w0_p, w1_p;
                boost::array<double, 2> w0_m, w1_m;
#if (NDIM == 3)
                boost::array<double, 2> w2_m, w2_p;
#endif
                boost::array<double, NDIM> x_lower_axis, x_upper_axis;
                const int local_sz = (*std::max_element(local_indices.begin(), local_indices.end())) + 1;
                std::vector<double> Q_data_axis_p(local_sz), Q_data_axis_m(local_sz);

                x_lower_axis[0] = x_lower_axis[1] = x_upper_axis[0] = x_upper_axis[1] = 0.0;
#if (NDIM == 3)
                x_lower_axis[2] = x_lower_axis[2] = 0.0;
#endif

                Box<NDIM> side_boxes[NDIM];

                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        x_lower_axis[d] = x_lower[d];
                        x_upper_axis[d] = x_upper[d];
                    }
                    x_lower_axis[axis] -= 0.5 * dx[axis];
                    x_upper_axis[axis] += 0.5 * dx[axis];

                    for (int d = 0; d < NDIM; ++d)
                    {
                        side_boxes[d] = SideGeometry<NDIM>::toSideBox(interp_box, d);
                    }

                    const IntVector<NDIM>& ilower = side_boxes[axis].lower();
                    const IntVector<NDIM>& iupper = side_boxes[axis].upper();

                    boost::const_multi_array_ref<double, NDIM + 1> u_sc_data_array(
                        u_sc_data->getPointer(axis),
                        (boost::extents[range(ilower[0] - u_gcw[0], iupper[0] + u_gcw[0] + 1)][range(
                            ilower[1] - u_gcw[1], iupper[1] + u_gcw[1] + 1)]
#if (NDIM == 3)
                                       [range(ilower[2] - u_gcw[2], iupper[2] + u_gcw[2] + 1)]
#endif
                                       [range(0, u_depth)]),
                        boost::fortran_storage_order());

                    for (unsigned int k = 0; k < nindices; ++k)
                    {
                        const int s = local_indices[k];

                        for (int d = 0; d < NDIM; ++d)
                        {
                            X_shifted_p[d] = X_qp_p[d + s * NDIM] + periodic_shifts[d + k * NDIM];
                            X_shifted_m[d] = X_qp_m[d + s * NDIM] + periodic_shifts[d + k * NDIM];
                        }

                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            ic_center_p[d] = ilower[d] + NINT((X_shifted_p[d] - x_lower_axis[d]) / dx[d] - 0.5);
                            X_cell_p[d] =
                                x_lower_axis[d] + (static_cast<double>(ic_center_p[d] - ilower[d]) + 0.5) * dx[d];

                            if (X_shifted_p[d] <= X_cell_p[d])
                            {
                                ic_lower_p[d] = ic_center_p[d] - 1;
                                ic_upper_p[d] = ic_center_p[d];
                            }
                            else
                            {
                                ic_lower_p[d] = ic_center_p[d];
                                ic_upper_p[d] = ic_center_p[d] + 1;
                            }
                            ic_trimmed_lower_p[d] = std::max(ic_lower_p[d], ilower[d] - u_gcw[d]);
                            ic_trimmed_upper_p[d] = std::min(ic_upper_p[d], iupper[d] + u_gcw[d]);

                            ic_center_m[d] = ilower[d] + NINT((X_shifted_m[d] - x_lower_axis[d]) / dx[d] - 0.5);
                            X_cell_m[d] =
                                x_lower_axis[d] + (static_cast<double>(ic_center_m[d] - ilower[d]) + 0.5) * dx[d];

                            if (X_shifted_m[d] <= X_cell_m[d])
                            {
                                ic_lower_m[d] = ic_center_m[d] - 1;
                                ic_upper_m[d] = ic_center_m[d];
                            }
                            else
                            {
                                ic_lower_m[d] = ic_center_m[d];
                                ic_upper_m[d] = ic_center_m[d] + 1;
                            }
                            ic_trimmed_lower_m[d] = std::max(ic_lower_m[d], ilower[d] - u_gcw[d]);
                            ic_trimmed_upper_m[d] = std::min(ic_upper_m[d], iupper[d] + u_gcw[d]);
                        }

                        if (X_shifted_p[0] <= X_cell_p[0])
                        {
                            w0_p[0] = (X_cell_p[0] - X_shifted_p[0]) / dx[0];

                            w0_p[1] = 1.0 - w0_p[0];
                        }
                        else
                        {
                            w0_p[0] = 1.0 + (X_cell_p[0] - X_shifted_p[0]) / dx[0];
                            w0_p[1] = 1.0 - w0_p[0];
                        }

                        if (X_shifted_m[0] <= X_cell_m[0])
                        {
                            w0_m[0] = (X_cell_m[0] - X_shifted_m[0]) / dx[0];

                            w0_m[1] = 1.0 - w0_m[0];
                        }
                        else
                        {
                            w0_m[0] = 1.0 + (X_cell_m[0] - X_shifted_m[0]) / dx[0];
                            w0_m[1] = 1.0 - w0_m[0];
                        }

                        if (X_shifted_p[1] <= X_cell_p[1])
                        {
                            w1_p[0] = (X_cell_p[1] - X_shifted_p[1]) / dx[1];
                            w1_p[1] = 1.0 - w1_p[0];
                        }
                        else
                        {
                            w1_p[0] = 1.0 + (X_cell_p[1] - X_shifted_p[1]) / dx[1];
                            w1_p[1] = 1.0 - w1_p[0];
                        }

                        if (X_shifted_m[1] <= X_cell_m[1])
                        {
                            w1_m[0] = (X_cell_m[1] - X_shifted_m[1]) / dx[1];
                            w1_m[1] = 1.0 - w1_m[0];
                        }
                        else
                        {
                            w1_m[0] = 1.0 + (X_cell_m[1] - X_shifted_m[1]) / dx[1];
                            w1_m[1] = 1.0 - w1_m[0];
                        }
#if (NDIM == 3)

                        if (X_shifted_p[2] <= X_cell_p[2])
                        {
                            w2_p[0] = (X_cell_p[2] - X_shifted_p[2]) / dx[2];
                            w2_p[1] = 1.0 - w2_p[0];
                        }
                        else
                        {
                            w2_p[0] = 1.0 + (X_cell_p[2] - X_shifted_p[2]) / dx[2];
                            w2_p[1] = 1.0 - w2_p[0];
                        }

                        if (X_shifted_m[2] <= X_cell_m[2])
                        {
                            w2_m[0] = (X_cell_m[2] - X_shifted_m[2]) / dx[2];
                            w2_m[1] = 1.0 - w2_m[0];
                        }
                        else
                        {
                            w2_m[0] = 1.0 + (X_cell_m[2] - X_shifted_m[2]) / dx[2];
                            w2_m[1] = 1.0 - w2_m[0];
                        }

#endif

                        for (int d = 0; d < u_depth; ++d)
                        {
                            Q_data_axis_p[s] = 0.0;
                            Q_data_axis_m[s] = 0.0;

#if (NDIM == 2)

                            for (int ic1 = ic_trimmed_lower_p[1]; ic1 <= ic_trimmed_upper_p[1]; ++ic1)
                            {
                                for (int ic0 = ic_trimmed_lower_p[0]; ic0 <= ic_trimmed_upper_p[0]; ++ic0)
                                {
                                    Q_data_axis_p[s] += w0_p[ic0 - ic_lower_p[0]] * w1_p[ic1 - ic_lower_p[1]] *
                                                        u_sc_data_array[ic0][ic1][d];
                                }
                            }

                            for (int ic1 = ic_trimmed_lower_m[1]; ic1 <= ic_trimmed_upper_m[1]; ++ic1)
                            {
                                for (int ic0 = ic_trimmed_lower_m[0]; ic0 <= ic_trimmed_upper_m[0]; ++ic0)
                                {
                                    Q_data_axis_m[s] += w0_m[ic0 - ic_lower_m[0]] * w1_m[ic1 - ic_lower_m[1]] *
                                                        u_sc_data_array[ic0][ic1][d];
                                }
                            }
#endif
#if (NDIM == 3)

                            for (int ic2 = ic_trimmed_lower_p[2]; ic2 <= ic_trimmed_upper_p[2]; ++ic2)
                            {
                                for (int ic1 = ic_trimmed_lower_p[1]; ic1 <= ic_trimmed_upper_p[1]; ++ic1)
                                {
                                    for (int ic0 = ic_trimmed_lower_p[0]; ic0 <= ic_trimmed_upper_p[0]; ++ic0)
                                    {
                                        Q_data_axis_p[s] += w0_p[ic0 - ic_lower_p[0]] * w1_p[ic1 - ic_lower_p[1]] *
                                                            w2_p[ic2 - ic_lower_p[2]] *
                                                            u_sc_data_array[ic0][ic1][ic2][d];
                                    }
                                }
                            }

                            for (int ic2 = ic_trimmed_lower_m[2]; ic2 <= ic_trimmed_upper_m[2]; ++ic2)
                            {
                                for (int ic1 = ic_trimmed_lower_m[1]; ic1 <= ic_trimmed_upper_m[1]; ++ic1)
                                {
                                    for (int ic0 = ic_trimmed_lower_m[0]; ic0 <= ic_trimmed_upper_m[0]; ++ic0)
                                    {
                                        Q_data_axis_m[s] += w0_m[ic0 - ic_lower_m[0]] * w1_m[ic1 - ic_lower_m[1]] *
                                                            w2_m[ic2 - ic_lower_m[2]] *
                                                            u_sc_data_array[ic0][ic1][ic2][d];
                                    }
                                }
                            }

#endif

                        } // depth
                    }

                    double CC = 0.0;

                    for (unsigned int k = 0; k < nindices; ++k)
                    {
                        if (d_modify_vel_interp_jumps)
                        {
#if (NDIM == 2)
                            if (axis == 0)
                                CC = 0.5 * dh *
                                     (du_j_qp[n_vars * local_indices[k]] * N_qp[n_vars * local_indices[k]] +
                                      du_j_qp[n_vars * local_indices[k] + 1] * N_qp[n_vars * local_indices[k] + 1]);
                            else if (axis == 1)
                                CC = 0.5 * dh *
                                     (dv_j_qp[n_vars * local_indices[k]] * N_qp[n_vars * local_indices[k]] +
                                      dv_j_qp[n_vars * local_indices[k] + 1] * N_qp[n_vars * local_indices[k] + 1]);
#endif
#if (NDIM == 3)
                            if (axis == 0)
                                CC = 0.5 * dh *
                                     (du_j_qp[n_vars * local_indices[k]] * N_qp[n_vars * local_indices[k]] +
                                      du_j_qp[n_vars * local_indices[k] + 1] * N_qp[n_vars * local_indices[k] + 1] +
                                      du_j_qp[n_vars * local_indices[k] + 2] * N_qp[n_vars * local_indices[k] + 2]);
                            else if (axis == 1)
                                CC = 0.5 * dh *
                                     (dv_j_qp[n_vars * local_indices[k]] * N_qp[n_vars * local_indices[k]] +
                                      dv_j_qp[n_vars * local_indices[k] + 1] * N_qp[n_vars * local_indices[k] + 1] +
                                      dv_j_qp[n_vars * local_indices[k] + 2] * N_qp[n_vars * local_indices[k] + 2]);
                            if (axis == 2)
                                CC = 0.5 * dh *
                                     (dw_j_qp[n_vars * local_indices[k]] * N_qp[n_vars * local_indices[k]] +
                                      dw_j_qp[n_vars * local_indices[k] + 1] * N_qp[n_vars * local_indices[k] + 1] +
                                      dw_j_qp[n_vars * local_indices[k] + 2] * N_qp[n_vars * local_indices[k] + 2]);
#endif
							U_qp[n_vars * local_indices[k] + axis] =
								0.5 * (Q_data_axis_p[local_indices[k]] + Q_data_axis_m[local_indices[k]]) - CC / d_mu;
							WSS_i_qp[n_vars * local_indices[k] + axis] = (U_qp[n_vars * local_indices[k] + axis] - Q_data_axis_m[local_indices[k]])/dh;
							WSS_o_qp[n_vars * local_indices[k] + axis] = (Q_data_axis_p[local_indices[k]] - U_qp[n_vars * local_indices[k] + axis])/dh;

                        }
                        else
                        {
							U_qp[n_vars * local_indices[k] + axis] = 0.5 * (Q_data_axis_p[local_indices[k]] + Q_data_axis_m[local_indices[k]]);
							WSS_i_qp[n_vars * local_indices[k] + axis] = 0.0;
							WSS_o_qp[n_vars * local_indices[k] + axis] = 0.0;
							
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
                for (unsigned int i = 0; i < n_vars; ++i)
                {
                    U_dof_map_cache.dof_indices(elem, U_dof_indices[i], i);
                    U_rhs_e[i].resize(static_cast<int>(U_dof_indices[i].size()));
                    WSS_i_dof_map_cache.dof_indices(elem, WSS_i_dof_indices[i], i);
                    WSS_i_rhs_e[i].resize(static_cast<int>(WSS_i_dof_indices[i].size()));
                    
                    WSS_o_dof_map_cache.dof_indices(elem, WSS_o_dof_indices[i], i);
                    WSS_o_rhs_e[i].resize(static_cast<int>(WSS_o_dof_indices[i].size()));
                }
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X_dof_map_cache.dof_indices(elem, X_dof_indices[d], d);
                }
                get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);
                const bool qrule_changed = FEDataManager::updateInterpQuadratureRule(
                    qrule, IBFEMethod::getDefaultInterpSpec(), elem, X_node, dx_min);
                if (qrule_changed)
                {
                    // NOTE: Because we are only using the shape function values for
                    // the FE object associated with X, we only need to reinitialize
                    // X_fe whenever the quadrature rule changes.  In particular,
                    // notice that the shape function values depend only on the
                    // element type and quadrature rule, not on the element
                    // geometry.
                    U_fe->attach_quadrature_rule(qrule.get());
                    X_fe->attach_quadrature_rule(qrule.get());
                    if (X_fe != U_fe) X_fe->reinit(elem);
                }
                U_fe->reinit(elem);
                const unsigned int n_qp = qrule->n_points();
                const size_t n_basis = U_dof_indices[0].size();
                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    const int idx = n_vars * (qp_offset + qp);
                    for (unsigned int k = 0; k < n_basis; ++k)
                    {
                        const double p_JxW_U = phi_U[k][qp] * JxW_U[qp];
                        for (unsigned int i = 0; i < n_vars; ++i)
                        {
                            U_rhs_e[i](k) += U_qp[idx + i] * p_JxW_U;
                            WSS_i_rhs_e[i](k) += WSS_i_qp[idx + i] * p_JxW_U;
                            WSS_o_rhs_e[i](k) += WSS_o_qp[idx + i] * p_JxW_U;
                        }
                    }
                }
                for (unsigned int i = 0; i < n_vars; ++i)
                {
                    U_dof_map.constrain_element_vector(U_rhs_e[i], U_dof_indices[i]);
                    U_rhs_vec->add_vector(U_rhs_e[i], U_dof_indices[i]);
                    WSS_i_dof_map.constrain_element_vector(WSS_i_rhs_e[i], WSS_i_dof_indices[i]);
                    WSS_i_rhs_vec->add_vector(WSS_i_rhs_e[i], WSS_i_dof_indices[i]);
                    
                    WSS_o_dof_map.constrain_element_vector(WSS_o_rhs_e[i], WSS_o_dof_indices[i]);
                    WSS_o_rhs_vec->add_vector(WSS_o_rhs_e[i], WSS_o_dof_indices[i]);
                }
                qp_offset += n_qp;
            }
        }

        U_rhs_vec->close();
        WSS_i_rhs_vec->close();
        WSS_o_rhs_vec->close();

        VecRestoreArray(X_local_vec, &X_local_soln);
        VecGhostRestoreLocalForm(X_global_vec, &X_local_vec);
        VecRestoreArray(du_j_local_vec, &du_j_local_soln);
        VecGhostRestoreLocalForm(du_j_global_vec, &du_j_local_vec);
        VecRestoreArray(dv_j_local_vec, &dv_j_local_soln);
        VecGhostRestoreLocalForm(dv_j_global_vec, &dv_j_local_vec);
        VecRestoreArray(dw_j_local_vec, &dw_j_local_soln);
        VecGhostRestoreLocalForm(dw_j_global_vec, &dw_j_local_vec);

        d_fe_data_managers[part]->computeL2Projection(
            *U_vec, *U_rhs_vec, VELOCITY_SYSTEM_NAME, IBFEMethod::getDefaultInterpSpec().use_consistent_mass_matrix);
        d_fe_data_managers[part]->computeL2Projection(
            *WSS_i_vec, *WSS_i_rhs_vec, WSS_I_SYSTEM_NAME, IBFEMethod::getDefaultInterpSpec().use_consistent_mass_matrix);
            
        d_fe_data_managers[part]->computeL2Projection(
            *WSS_o_vec, *WSS_o_rhs_vec, WSS_O_SYSTEM_NAME, IBFEMethod::getDefaultInterpSpec().use_consistent_mass_matrix);

        d_du_j_IB_ghost_vecs[part]->close();
        d_dv_j_IB_ghost_vecs[part]->close();
        d_dw_j_IB_ghost_vecs[part]->close();
        d_X_half_vecs[part]->close();
        d_X_current_vecs[part]->close();
        d_X_new_vecs[part]->close();
        d_U_new_vecs[part]->close();
        d_WSS_i_half_vecs[part]->close();
        d_WSS_o_half_vecs[part]->close();
    }

    return;
    //~
} // interpolateVelocityWithJump
 






void
IBFEMethod::forwardEulerStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr = VecWAXPY(d_X_new_vecs[part]->vec(), dt, d_U_current_vecs[part]->vec(), d_X_current_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(
            d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_X_new_vecs[part]->close();
        d_X_half_vecs[part]->close();
    }
    return;
} // eulerStep

void
IBFEMethod::midpointStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr = VecWAXPY(d_X_new_vecs[part]->vec(), dt, d_U_half_vecs[part]->vec(), d_X_current_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(
            d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_X_new_vecs[part]->close();
        d_X_half_vecs[part]->close();
    }
    return;
} // midpointStep

void
IBFEMethod::trapezoidalStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr =
            VecWAXPY(d_X_new_vecs[part]->vec(), 0.5 * dt, d_U_current_vecs[part]->vec(), d_X_current_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPY(d_X_new_vecs[part]->vec(), 0.5 * dt, d_U_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(
            d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_X_new_vecs[part]->close();
        d_X_half_vecs[part]->close();
    }
    return;
} // trapezoidalStep

void
IBFEMethod::computeLagrangianForce(const double data_time)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        computeInteriorForceDensity(*d_F_half_vecs[part],
									*d_F_n_half_vecs[part],	
									*d_F_t_half_vecs[part],
									*d_F_b_half_vecs[part],
									*d_H_half_vecs[part],
                                    *d_X_half_vecs[part],
									*d_P_j_half_vecs[part],
                                    *d_dP_j_half_vecs[part],
                                    *d_du_j_half_vecs[part],
                                    *d_dv_j_half_vecs[part],
                                    *d_dw_j_half_vecs[part],
                                    *d_d2u_j_half_vecs[part],
                                    *d_d2v_j_half_vecs[part],
                                    *d_d2w_j_half_vecs[part],
                                    data_time,
                                    part);
    }
    return;
} // computeLagrangianForce

void
IBFEMethod::spreadForce(const int f_data_idx,
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
        PetscVector<double>* P_j_vec = d_P_j_half_vecs[part];
        PetscVector<double>* P_j_ghost_vec = d_P_j_IB_ghost_vecs[part];
        PetscVector<double>* du_j_ghost_vec = d_du_j_IB_ghost_vecs[part];
        PetscVector<double>* dv_j_ghost_vec = d_dv_j_IB_ghost_vecs[part];
        PetscVector<double>* dw_j_ghost_vec = d_dw_j_IB_ghost_vecs[part];
        PetscVector<double>* du_j_vec = d_du_j_half_vecs[part];
        PetscVector<double>* dv_j_vec = d_dv_j_half_vecs[part];
        PetscVector<double>* dw_j_vec = d_dw_j_half_vecs[part];
        PetscVector<double>* d2u_j_ghost_vec = d_d2u_j_IB_ghost_vecs[part];
        PetscVector<double>* d2v_j_ghost_vec = d_d2v_j_IB_ghost_vecs[part];
        PetscVector<double>* d2w_j_ghost_vec = d_d2w_j_IB_ghost_vecs[part];
        PetscVector<double>* d2u_j_vec = d_d2u_j_half_vecs[part];
        PetscVector<double>* d2v_j_vec = d_d2v_j_half_vecs[part];
        PetscVector<double>* d2w_j_vec = d_d2w_j_half_vecs[part];


        PetscVector<double>* dP_j_ghost_vec = d_dP_j_IB_ghost_vecs[part];
        PetscVector<double>* dP_j_vec = d_dP_j_half_vecs[part];

        X_vec->localize(*X_ghost_vec);
        F_vec->localize(*F_ghost_vec);
        P_j_vec->localize(*P_j_ghost_vec);
		dP_j_vec->localize(*dP_j_ghost_vec);
        du_j_vec->localize(*du_j_ghost_vec);
        dv_j_vec->localize(*dv_j_ghost_vec);
        dw_j_vec->localize(*dw_j_ghost_vec);
        d2u_j_vec->localize(*d2u_j_ghost_vec);
        d2v_j_vec->localize(*d2v_j_ghost_vec);
		d2w_j_vec->localize(*d2w_j_ghost_vec);

        d_fe_data_managers[part]->spread(
            f_data_idx, *F_ghost_vec, *X_ghost_vec, FORCE_SYSTEM_NAME, f_phys_bdry_op, data_time);
        if (d_split_normal_force || d_split_tangential_force)
        {
            if (d_use_jump_conditions)
            {
                imposeJumpConditionsWeak(f_data_idx,
                                     mask_current_idx,
                                     *F_ghost_vec,
                                     *X_ghost_vec,
                                     *P_j_ghost_vec,
									 *dP_j_ghost_vec,
                                     *du_j_ghost_vec,
                                     *dv_j_ghost_vec,
                                     *dw_j_ghost_vec,
                                     *d2u_j_ghost_vec,
                                     *d2v_j_ghost_vec,
                                     *d2w_j_ghost_vec,
                                     data_time,
                                     part);
            }
        }
    }
    return;
} // spreadForce

FEDataManager::InterpSpec
IBFEMethod::getDefaultInterpSpec() const
{
    return d_default_interp_spec;
}

FEDataManager::SpreadSpec
IBFEMethod::getDefaultSpreadSpec() const
{
    return d_default_spread_spec;
}

void
IBFEMethod::setInterpSpec(const FEDataManager::InterpSpec& interp_spec, const unsigned int part)
{
    TBOX_ASSERT(!d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_num_parts);
    d_interp_spec[part] = interp_spec;
    return;
}

void
IBFEMethod::setSpreadSpec(const FEDataManager::SpreadSpec& spread_spec, const unsigned int part)
{
    TBOX_ASSERT(!d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_num_parts);
    d_spread_spec[part] = spread_spec;
    return;
}

void
IBFEMethod::initializeFEEquationSystems()
{
    if (d_fe_equation_systems_initialized) return;

    const bool from_restart = RestartManager::getManager()->isFromRestart();

    // Create the FE data managers that manage mappings between the FE mesh
    // parts and the Cartesian grid.
    d_equation_systems.resize(d_num_parts, NULL);
    d_fe_data_managers.resize(d_num_parts, NULL);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Create FE data managers.
        std::ostringstream manager_stream;
        manager_stream << "IBFEMethod FEDataManager::" << part;
        const std::string& manager_name = manager_stream.str();
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
            System& X_system = equation_systems->add_system<System>(COORDS_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "X_" << d;
                X_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }

            System& X0_system = equation_systems->add_system<System>(COORDS0_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "X0_" << d;
                X0_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }

            System& dX_system = equation_systems->add_system<System>(COORD_MAPPING_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "dX_" << d;
                dX_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }

            System& U_system = equation_systems->add_system<System>(VELOCITY_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "U_" << d;
                U_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }
            System& WSS_i_system = equation_systems->add_system<System>(WSS_I_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "WSS_i" << d;
                WSS_i_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }
            
            
            
            System& WSS_o_system = equation_systems->add_system<System>(WSS_O_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "WSS_o" << d;
                WSS_o_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }
            
            System& F_system = equation_systems->add_system<System>(FORCE_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "F_" << d;
                F_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }

            System& du_j_system = equation_systems->add_system<System>(DU_J_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "du_" << d;
                du_j_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }

            System& dv_j_system = equation_systems->add_system<System>(DV_J_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "dv_" << d;
                dv_j_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }

            System& dw_j_system = equation_systems->add_system<System>(DW_J_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "dw_" << d;
                dw_j_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }

            System& dP_j_system = equation_systems->add_system<System>(DP_J_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "dP_" << d;
                dP_j_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }


            System& P_j_system = equation_systems->add_system<System>(P_J_SYSTEM_NAME);
            P_j_system.add_variable("[[p]]", d_fe_order[part], d_fe_family[part]);
            
            System& F_n_system = equation_systems->add_system<System>(FORCE_N_SYSTEM_NAME);
            F_n_system.add_variable("[[Fn]]", d_fe_order[part], d_fe_family[part]);
            
            System& F_t_system = equation_systems->add_system<System>(FORCE_T_SYSTEM_NAME);
            F_t_system.add_variable("[[Ft]]", d_fe_order[part], d_fe_family[part]);
            
            
            System& F_b_system = equation_systems->add_system<System>(FORCE_B_SYSTEM_NAME);
            F_b_system.add_variable("[[Fb]]", d_fe_order[part], d_fe_family[part]);
            
            System& H_system = equation_systems->add_system<System>(H_SYSTEM_NAME);
            H_system.add_variable("[[H]]", d_fe_order[part], d_fe_family[part]);
            
            
            
            System& d2u_j_system = equation_systems->add_system<System>(D2U_J_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "d2u_" << d;
                d2u_j_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }


            System& d2v_j_system = equation_systems->add_system<System>(D2V_J_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "d2v_" << d;
                d2v_j_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }
            

            System& d2w_j_system = equation_systems->add_system<System>(D2W_J_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "d2w_" << d;
                d2w_j_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }

            
            
            
            
            
        }
    }
    d_fe_equation_systems_initialized = true;
    return;
}

void
IBFEMethod::initializeFEData()
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
        }
        updateCoordinateMapping(part);

        // Assemble systems.
        System& X_system = equation_systems->get_system<System>(COORDS_SYSTEM_NAME);
        System& X0_system = equation_systems->get_system<System>(COORDS0_SYSTEM_NAME);
        System& dX_system = equation_systems->get_system<System>(COORD_MAPPING_SYSTEM_NAME);
        System& U_system = equation_systems->get_system<System>(VELOCITY_SYSTEM_NAME);
        System& WSS_i_system = equation_systems->get_system<System>(WSS_I_SYSTEM_NAME);
        System& WSS_o_system = equation_systems->get_system<System>(WSS_O_SYSTEM_NAME);
        System& F_system = equation_systems->get_system<System>(FORCE_SYSTEM_NAME);
        System& du_j_system = equation_systems->get_system<System>(DU_J_SYSTEM_NAME);
        System& dv_j_system = equation_systems->get_system<System>(DV_J_SYSTEM_NAME);
        System& dw_j_system = equation_systems->get_system<System>(DW_J_SYSTEM_NAME);
        System& d2u_j_system = equation_systems->get_system<System>(D2U_J_SYSTEM_NAME);
        System& d2v_j_system = equation_systems->get_system<System>(D2V_J_SYSTEM_NAME);
		System& d2w_j_system = equation_systems->get_system<System>(D2W_J_SYSTEM_NAME);
		System& dP_j_system = equation_systems->get_system<System>(DP_J_SYSTEM_NAME);
        System& P_j_system = equation_systems->get_system<System>(P_J_SYSTEM_NAME);
		System& F_t_system = equation_systems->get_system<System>(FORCE_T_SYSTEM_NAME);
        System& F_n_system = equation_systems->get_system<System>(FORCE_N_SYSTEM_NAME);
        System& F_b_system = equation_systems->get_system<System>(FORCE_B_SYSTEM_NAME);
        X_system.assemble_before_solve = false;
        X_system.assemble();

        X0_system.assemble_before_solve = false;
        X0_system.assemble();

        dX_system.assemble_before_solve = false;
        dX_system.assemble();

        U_system.assemble_before_solve = false;
        U_system.assemble();
        WSS_i_system.assemble_before_solve = false;
        WSS_i_system.assemble();
        
        WSS_o_system.assemble_before_solve = false;
        WSS_o_system.assemble();

        F_system.assemble_before_solve = false;
        F_system.assemble();

        F_t_system.assemble_before_solve = false;
        F_t_system.assemble();

        F_n_system.assemble_before_solve = false;
        F_n_system.assemble();

        F_b_system.assemble_before_solve = false;
        F_b_system.assemble();

        du_j_system.assemble_before_solve = false;
        du_j_system.assemble();

        dv_j_system.assemble_before_solve = false;
        dv_j_system.assemble();

        dw_j_system.assemble_before_solve = false;
        dw_j_system.assemble();

        d2u_j_system.assemble_before_solve = false;
        d2u_j_system.assemble();

        d2v_j_system.assemble_before_solve = false;
        d2v_j_system.assemble();

        d2w_j_system.assemble_before_solve = false;
        d2w_j_system.assemble();

        dP_j_system.assemble_before_solve = false;
        dP_j_system.assemble();
        P_j_system.assemble_before_solve = false;
        P_j_system.assemble();

  // Set up boundary conditions.  Specifically, add appropriate boundary
        // IDs to the BoundaryInfo object associated with the mesh, and add DOF
        // constraints for the nodal forces and velocities.
        const MeshBase& mesh = equation_systems->get_mesh();
        DofMap& F_dof_map = F_system.get_dof_map();
        DofMap& U_dof_map = U_system.get_dof_map();
        const unsigned int F_sys_num = F_system.number();
        const unsigned int U_sys_num = U_system.number();
        DofMap& WSS_i_dof_map = WSS_i_system.get_dof_map();
        const unsigned int WSS_i_sys_num = WSS_i_system.number();
        DofMap& WSS_o_dof_map = WSS_o_system.get_dof_map();
        const unsigned int WSS_o_sys_num = WSS_o_system.number();
        MeshBase::const_element_iterator el_it = mesh.elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.elements_end();
        for (; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor(side);
                if (!at_mesh_bdry) continue;

                static const short int dirichlet_bdry_id_set[3] = { FEDataManager::ZERO_DISPLACEMENT_X_BDRY_ID,
                                                                    FEDataManager::ZERO_DISPLACEMENT_Y_BDRY_ID,
                                                                    FEDataManager::ZERO_DISPLACEMENT_Z_BDRY_ID };
                const short int dirichlet_bdry_ids =
                    get_dirichlet_bdry_ids(mesh.boundary_info->boundary_ids(elem, side));
                if (!dirichlet_bdry_ids) continue;
                for (unsigned int n = 0; n < elem->n_nodes(); ++n)
                {
                    if (!elem->is_node_on_side(n, side)) continue;

                    Node* node = elem->get_node(n);
                    mesh.boundary_info->add_node(node, dirichlet_bdry_ids);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (!(dirichlet_bdry_ids & dirichlet_bdry_id_set[d])) continue;
                        if (node->n_dofs(F_sys_num))
                        {
                            const int F_dof_index = node->dof_number(F_sys_num, d, 0);
                            DofConstraintRow F_constraint_row;
                            F_constraint_row[F_dof_index] = 1.0;
                            F_dof_map.add_constraint_row(F_dof_index, F_constraint_row, 0.0, false);
                        }
                        if (node->n_dofs(U_sys_num))
                        {
                            const int U_dof_index = node->dof_number(U_sys_num, d, 0);
                            DofConstraintRow U_constraint_row;
                            U_constraint_row[U_dof_index] = 1.0;
                            U_dof_map.add_constraint_row(U_dof_index, U_constraint_row, 0.0, false);
                        }
                         if (node->n_dofs(WSS_i_sys_num))
                        {
                            const int WSS_i_dof_index = node->dof_number(WSS_i_sys_num, d, 0);
                            DofConstraintRow WSS_i_constraint_row;
                            WSS_i_constraint_row[WSS_i_dof_index] = 1.0;
                            WSS_i_dof_map.add_constraint_row(WSS_i_dof_index, WSS_i_constraint_row, 0.0, false);
                        }
                        
                        if (node->n_dofs(WSS_o_sys_num))
                        {
                            const int WSS_o_dof_index = node->dof_number(WSS_o_sys_num, d, 0);
                            DofConstraintRow WSS_o_constraint_row;
                            WSS_o_constraint_row[WSS_o_dof_index] = 1.0;
                            WSS_o_dof_map.add_constraint_row(WSS_o_dof_index, WSS_o_constraint_row, 0.0, false);
                        }
                    }
                }
            }
        }
    }
    d_fe_data_initialized = true;
    return;
} // initializeFEData

void
IBFEMethod::registerEulerianVariables()
{
    const IntVector<NDIM> ghosts = 1;
    mask_var = new SideVariable<NDIM, double>(d_object_name + "::mask");
    registerVariable(mask_current_idx,
                     mask_new_idx,
                     mask_scratch_idx,
                     mask_var,
                     ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");
    return;
} // registerEulerianVariables

void
IBFEMethod::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
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

    d_mu = getINSHierarchyIntegrator()->getStokesSpecifications()->getMu();
    d_is_initialized = true;
    return;
} // initializePatchHierarchy

void
IBFEMethod::registerLoadBalancer(Pointer<LoadBalancer<NDIM> > load_balancer, int workload_data_idx)
{
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
IBFEMethod::updateWorkloadEstimates(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/, int /*workload_data_idx*/)
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->updateWorkloadEstimates();
    }
    return;
} // updateWorkloadEstimates

void IBFEMethod::beginDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                         Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    // intentionally blank
    return;
} // beginDataRedistribution

void IBFEMethod::endDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
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
IBFEMethod::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
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
        if (d_load_balancer && level_number == d_fe_data_managers[part]->getLevelNumber())
        {
            d_load_balancer->setWorkloadPatchDataIndex(d_workload_idx, level_number);
            d_fe_data_managers[part]->updateWorkloadEstimates(level_number, level_number);
        }
    }
    return;
} // initializeLevelData

void
IBFEMethod::resetHierarchyConfiguration(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
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
IBFEMethod::applyGradientDetector(Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
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
IBFEMethod::putToDatabase(Pointer<Database> db)
{
    db->putInteger("IBFE_METHOD_VERSION", IBFE_METHOD_VERSION);
    db->putInteger("d_num_parts", d_num_parts);
    db->putIntegerArray("d_ghosts", d_ghosts, NDIM);
    db->putBool("d_split_normal_force", d_split_normal_force);
    db->putBool("d_split_tangential_force", d_split_tangential_force);
    db->putBool("d_use_jump_conditions", d_use_jump_conditions);
    db->putDouble("d_vel_interp_width", d_vel_interp_width);
    db->putBool("d_modify_vel_interp_jumps", d_modify_vel_interp_jumps);
	db->putBool("d_use_higher_order_jump", d_use_higher_order_jump);
    //db->putString("d_fe_family", Utility::enum_to_string<FEFamily>(d_fe_family));
    //db->putString("d_fe_order", Utility::enum_to_string<Order>(d_fe_order));
    //db->putString("d_quad_type", Utility::enum_to_string<QuadratureType>(d_quad_type));
    //db->putString("d_quad_order", Utility::enum_to_string<Order>(d_quad_order));
    db->putBool("d_use_consistent_mass_matrix", d_use_consistent_mass_matrix);
    return;
} // putToDatabase

void
IBFEMethod::writeFEDataToRestartFile(const std::string& restart_dump_dirname, unsigned int time_step_number)
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

void
IBFEMethod::computeInteriorForceDensity(PetscVector<double>& F_vec,
										PetscVector<double>& F_n_vec,
										PetscVector<double>& F_t_vec,
										PetscVector<double>& F_b_vec,
										PetscVector<double>& H_vec,
                                        PetscVector<double>& X_vec,
										PetscVector<double>& P_j_vec,
                                        PetscVector<double>& dP_j_vec,
                                        PetscVector<double>& du_j_vec,
                                        PetscVector<double>& dv_j_vec,
                                        PetscVector<double>& dw_j_vec,
                                        PetscVector<double>& d2u_j_vec,
                                        PetscVector<double>& d2v_j_vec,
                                        PetscVector<double>& d2w_j_vec,
                                        const double data_time,
                                        const unsigned int part)
{
	
    // Extract the mesh.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // Setup global and elemental right-hand-side vectors.
    AutoPtr<NumericVector<double> > F_rhs_vec = F_vec.zero_clone();
    DenseVector<double> F_rhs_e[NDIM];

    AutoPtr<NumericVector<double> > F_t_rhs_vec = F_t_vec.zero_clone();
    DenseVector<double> F_t_rhs_e;
    
    AutoPtr<NumericVector<double> > H_rhs_vec = H_vec.zero_clone();
    DenseVector<double> H_rhs_e;
    
    AutoPtr<NumericVector<double> > F_n_rhs_vec = F_n_vec.zero_clone();
    DenseVector<double> F_n_rhs_e;
    
    AutoPtr<NumericVector<double> > F_b_rhs_vec = F_b_vec.zero_clone();
    DenseVector<double> F_b_rhs_e;
    
    AutoPtr<NumericVector<double> > dP_j_rhs_vec = dP_j_vec.zero_clone();
    DenseVector<double> dP_j_rhs_e[NDIM];
    
    AutoPtr<NumericVector<double> > P_j_rhs_vec = P_j_vec.zero_clone();
    DenseVector<double> P_j_rhs_e;

    AutoPtr<NumericVector<double> > du_j_rhs_vec = du_j_vec.zero_clone();
    DenseVector<double> du_j_rhs_e[NDIM];

    AutoPtr<NumericVector<double> > dv_j_rhs_vec = dv_j_vec.zero_clone();
    DenseVector<double> dv_j_rhs_e[NDIM];

    AutoPtr<NumericVector<double> > dw_j_rhs_vec = dw_j_vec.zero_clone();
    DenseVector<double> dw_j_rhs_e[NDIM];

    AutoPtr<NumericVector<double> > d2u_j_rhs_vec = d2u_j_vec.zero_clone();
    DenseVector<double> d2u_j_rhs_e[NDIM];

    AutoPtr<NumericVector<double> > d2v_j_rhs_vec = d2v_j_vec.zero_clone();
    DenseVector<double> d2v_j_rhs_e[NDIM];

    AutoPtr<NumericVector<double> > d2w_j_rhs_vec = d2w_j_vec.zero_clone();
    DenseVector<double> d2w_j_rhs_e[NDIM];

    // Extract the FE systems and DOF maps, and setup the FE objects.
    System& F_system = equation_systems->get_system(FORCE_SYSTEM_NAME);
    const DofMap& F_dof_map = F_system.get_dof_map();
    FEDataManager::SystemDofMapCache& F_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(FORCE_SYSTEM_NAME);
    FEType F_fe_type = F_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(F_dof_map.variable_type(d) == F_fe_type);
    }
    std::vector<std::vector<unsigned int> > F_dof_indices(NDIM);

    System& F_t_system = equation_systems->get_system(FORCE_T_SYSTEM_NAME);
    const DofMap& F_t_dof_map = F_t_system.get_dof_map();
    FEDataManager::SystemDofMapCache& F_t_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(FORCE_T_SYSTEM_NAME);
    TBOX_ASSERT(F_t_dof_map.variable_type(0) == F_fe_type);
    std::vector<unsigned int>  F_t_dof_indices;
    
    
    System& F_n_system = equation_systems->get_system(FORCE_N_SYSTEM_NAME);
    const DofMap& F_n_dof_map = F_n_system.get_dof_map();
    FEDataManager::SystemDofMapCache& F_n_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(FORCE_N_SYSTEM_NAME);
    TBOX_ASSERT(F_n_dof_map.variable_type(0) == F_fe_type);
    std::vector<unsigned int>  F_n_dof_indices;

    System& F_b_system = equation_systems->get_system(FORCE_B_SYSTEM_NAME);
    const DofMap& F_b_dof_map = F_b_system.get_dof_map();
    FEDataManager::SystemDofMapCache& F_b_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(FORCE_B_SYSTEM_NAME);
    TBOX_ASSERT(F_b_dof_map.variable_type(0) == F_fe_type);
    std::vector<unsigned int>  F_b_dof_indices;
    
    System& H_system = equation_systems->get_system(H_SYSTEM_NAME);
    const DofMap& H_dof_map = H_system.get_dof_map();
    FEDataManager::SystemDofMapCache& H_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(H_SYSTEM_NAME);
    TBOX_ASSERT(H_dof_map.variable_type(0) == F_fe_type);
    std::vector<unsigned int>  H_dof_indices;
    

    System& du_j_system = equation_systems->get_system(DU_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& du_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DU_J_SYSTEM_NAME);
    const DofMap& du_j_dof_map = du_j_system.get_dof_map();
    FEType du_j_fe_type = du_j_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(du_j_dof_map.variable_type(d) == du_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > du_j_dof_indices(NDIM);

    System& dv_j_system = equation_systems->get_system(DV_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& dv_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DV_J_SYSTEM_NAME);
    const DofMap& dv_j_dof_map = dv_j_system.get_dof_map();
    FEType dv_j_fe_type = dv_j_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(dv_j_dof_map.variable_type(d) == dv_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > dv_j_dof_indices(NDIM);

    System& dw_j_system = equation_systems->get_system(DW_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& dw_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DW_J_SYSTEM_NAME);
    const DofMap& dw_j_dof_map = dw_j_system.get_dof_map();
    FEType dw_j_fe_type = dw_j_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(dw_j_dof_map.variable_type(d) == dw_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > dw_j_dof_indices(NDIM);

    System& d2u_j_system = equation_systems->get_system(D2U_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& d2u_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(D2U_J_SYSTEM_NAME);
    const DofMap& d2u_j_dof_map = d2u_j_system.get_dof_map();
    FEType d2u_j_fe_type = d2u_j_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d2u_j_dof_map.variable_type(d) == d2u_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > d2u_j_dof_indices(NDIM);

    System& d2v_j_system = equation_systems->get_system(D2V_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& d2v_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(D2V_J_SYSTEM_NAME);
    const DofMap& d2v_j_dof_map = d2v_j_system.get_dof_map();
    FEType d2v_j_fe_type = d2v_j_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d2v_j_dof_map.variable_type(d) == d2v_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > d2v_j_dof_indices(NDIM);

    System& d2w_j_system = equation_systems->get_system(D2W_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& d2w_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(D2W_J_SYSTEM_NAME);
    const DofMap& d2w_j_dof_map = d2w_j_system.get_dof_map();
    FEType d2w_j_fe_type = d2w_j_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d2w_j_dof_map.variable_type(d) == d2w_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > d2w_j_dof_indices(NDIM);

    System& dP_j_system = equation_systems->get_system(DP_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& dP_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DP_J_SYSTEM_NAME);
    const DofMap& dP_j_dof_map = dP_j_system.get_dof_map();
    FEType dP_j_fe_type = dP_j_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(dP_j_dof_map.variable_type(d) == dP_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > dP_j_dof_indices(NDIM);

    System& P_j_system = equation_systems->get_system(P_J_SYSTEM_NAME);
    const DofMap& P_j_dof_map = P_j_system.get_dof_map();
    FEDataManager::SystemDofMapCache& P_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(P_J_SYSTEM_NAME);
    TBOX_ASSERT(P_j_dof_map.variable_type(0) == F_fe_type);
    std::vector<unsigned int> P_j_dof_indices;

    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& X_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    FEType X_fe_type = X_dof_map.variable_type(0);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
    }
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);

    System& X0_system = equation_systems->get_system(COORDS0_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& X0_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(COORDS0_SYSTEM_NAME);
    const DofMap& X0_dof_map = X0_system.get_dof_map();
    FEType X0_fe_type = X0_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X0_dof_map.variable_type(d) == X0_fe_type);
    }
    TBOX_ASSERT(X0_fe_type == X_fe_type);
    std::vector<std::vector<unsigned int> > X0_dof_indices(NDIM);

    std::vector<int> vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
    std::vector<int> no_vars;

    FEDataInterpolation fe(dim, d_fe_data_managers[part]);
    AutoPtr<QBase> qrule = QBase::build(d_default_quad_type[part], dim, d_default_quad_order[part]);
    fe.attachQuadratureRule(qrule.get());
    fe.evalQuadraturePoints();
    fe.evalQuadratureWeights();
    fe.registerSystem(F_system, vars, vars); // compute phi and dphi for the force system
    std::vector<size_t> force_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        force_fcn_system_idxs, d_lag_force_fcn_data[part].system_data, equation_systems);
    fe.init(/*use_IB_ghosted_vecs*/ false);

    AutoPtr<FEBase> X_fe_base(FEBase::build(dim, X_fe_type));
    X_fe_base->attach_quadrature_rule(qrule.get());
    const std::vector<std::vector<double> >& X_phi = X_fe_base->get_phi();
    const std::vector<std::vector<double> >& X_dphi_dxi = X_fe_base->get_dphidxi();
    const std::vector<std::vector<double> >& X_dphi_deta = X_fe_base->get_dphideta();
    const std::vector<std::vector<double> >& X_d2phi_dxi2 = X_fe_base->get_d2phidxi2();
  //~ const std::vector<double>& X_curvature = X_fe_base->get_curvatures();
    const std::vector<libMesh::Point>& q_point = fe.getQuadraturePoints();
    const std::vector<double>& JxW = fe.getQuadratureWeights();
    const std::vector<std::vector<double> >& phi = fe.getPhi(F_fe_type);

    std::vector<const std::vector<double>*> force_var_data;
    std::vector<const std::vector<VectorValue<double> >*> force_grad_var_data;

    PetscVector<double>* X_petsc_vec = static_cast<PetscVector<double>*>(&X_vec);
    Vec X_global_vec = X_petsc_vec->vec();
    Vec X_local_vec;
    VecGhostGetLocalForm(X_global_vec, &X_local_vec);
    double* X_local_soln;
    VecGetArray(X_local_vec, &X_local_soln);

    PetscVector<double>* X0_petsc_vec = static_cast<PetscVector<double>*>(d_X0_vecs[part]);
    Vec X0_global_vec = X0_petsc_vec->vec();
    Vec X0_local_vec;
    VecGhostGetLocalForm(X0_global_vec, &X0_local_vec);
    double* X0_local_soln;
    VecGetArray(X0_local_vec, &X0_local_soln);

    // Loop over the elements to compute the right-hand side vector.
    TensorValue<double> FF;
    VectorValue<double> F, Ft, Fb, F_qp, N, Tau1, kk, Tau2, X, n, s, tau1, tau2, x, du, dv, dw, du_j_qp, dv_j_qp, dw_j_qp, dP_j, d2u_j, d2v_j, d2w_j;
    VectorValue<double> d2u, d2v, d2u_j_qp, d2v_j_qp, d2w_j_qp, dP_j_qp, tau1unit, tau2unit;  
    double dA_da, F_t, F_b, F_n, HK;
    boost::multi_array<double, 2> X_node, X0_node;
    boost::multi_array<double, 1> F_t_node, F_b_node, F_n_node, H_node;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        P_j_dof_map_cache.dof_indices(elem, P_j_dof_indices);
        P_j_rhs_e.resize(static_cast<int>(P_j_dof_indices.size()));
        
        H_dof_map_cache.dof_indices(elem, H_dof_indices);
        H_rhs_e.resize(static_cast<int>(H_dof_indices.size()));
        
        F_t_dof_map_cache.dof_indices(elem, F_t_dof_indices);
        F_t_rhs_e.resize(static_cast<int>(F_t_dof_indices.size()));

        F_b_dof_map_cache.dof_indices(elem, F_b_dof_indices);
        F_b_rhs_e.resize(static_cast<int>(F_b_dof_indices.size()));
                
        F_n_dof_map_cache.dof_indices(elem, F_n_dof_indices);
        F_n_rhs_e.resize(static_cast<int>(F_n_dof_indices.size()));

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            F_dof_map_cache.dof_indices(elem, F_dof_indices[d], d);
            F_rhs_e[d].resize(static_cast<int>(F_dof_indices[d].size()));
            X_dof_map_cache.dof_indices(elem, X_dof_indices[d], d);
            X0_dof_map_cache.dof_indices(elem, X0_dof_indices[d], d);

            du_j_dof_map_cache.dof_indices(elem, du_j_dof_indices[d], d);
            du_j_rhs_e[d].resize(static_cast<int>(du_j_dof_indices[d].size()));

            dv_j_dof_map_cache.dof_indices(elem, dv_j_dof_indices[d], d);
            dv_j_rhs_e[d].resize(static_cast<int>(dv_j_dof_indices[d].size()));

            //~ #if (NDIM == 3)
            dw_j_dof_map_cache.dof_indices(elem, dw_j_dof_indices[d], d);
            dw_j_rhs_e[d].resize(static_cast<int>(dw_j_dof_indices[d].size()));
            //~ #endif

        }
        fe.reinit(elem);
        fe.collectDataForInterpolation(elem);
        fe.interpolate(elem);
        X_fe_base->reinit(elem);
        get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);
        get_values_for_interpolation(X0_node, *X0_petsc_vec, X0_local_soln, X0_dof_indices);
        const unsigned int n_qp = qrule->n_points();
        const size_t n_basis = phi.size();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            s = q_point[qp];

            interpolate(&X(0), qp, X0_node, X_phi);
            interpolate(&Tau1(0), qp, X0_node, X_dphi_dxi);
            if (dim == 1)
                Tau2 = VectorValue<double>(0.0, 0.0, 1.0);
            else
                interpolate(&Tau2(0), qp, X0_node, X_dphi_deta);
            N = Tau1.cross(Tau2);

            interpolate(&x(0), qp, X_node, X_phi);
            interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
			tau1unit = tau1.unit();

            if (dim == 1)
                tau2 = VectorValue<double>(0.0, 0.0, 1.0);
            else
				interpolate(&tau2(0), qp, X_node, X_dphi_deta);
				

			
			
            n = tau1.cross(tau2);
           
	        tau2unit = tau2.unit();
            dA_da = sqrt(N * N) / sqrt(n * n);

            N = N.unit();
            n = n.unit();
            
            interpolate(&kk(0), qp, X_node,  X_d2phi_dxi2);

            HK = - (n(0)*kk(0) + n(1)*kk(1))/(tau1(0)*tau1(0) + tau1(1)*tau1(1));

            if (d_lag_force_fcn_data[part].fcn)
            {
                // Compute the value of the body force at the quadrature
                // point and add the corresponding forces to the
                // right-hand-side vector.
                fe.setInterpolatedDataPointers(force_var_data, force_grad_var_data, force_fcn_system_idxs, elem, qp);
                d_lag_force_fcn_data[part].fcn(F,
                                               FF,
                                               x,
                                               s,
                                               elem,
                                               force_var_data,
                                               force_grad_var_data,
                                               data_time,
                                               d_lag_force_fcn_data[part].ctx);

                // extract the pressure jump
                double P_j = F * n * dA_da;
                du = - dA_da * (F(0) - F * n * n(0)) * n; // [ux] , [uy], [uz]
                
                dv = - dA_da * (F(1) - F * n * n(1)) * n; // [vx] , [vy], [vz]
#if (NDIM == 3)
                dw = - dA_da * (F(2) - F * n * n(2)) * n; // [wx] , [wy], [wz]
#endif
                if (d_split_normal_force) 
		{
					
			F_t = F * tau1unit;
			F_n = F * n;
#if (NDIM == 3)
			F_b = F * tau2unit;
#endif
	 	}

                // Split off the normal and/or tangential parts of the force (if requested).
                if (d_split_normal_force && !d_split_tangential_force)
                {
                    F = F - (F * n) * n;
                }
                else if (d_split_tangential_force && !d_split_normal_force)
                {
                    F = (F * n) * n;
                }
                else if (d_split_tangential_force && d_split_normal_force)
                {
                    F = 0.0;
                }

                // Accumulate the RHS values.
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    F_qp = F * phi[k][qp] * JxW[qp];
                    double F_t_qp = F_t * phi[k][qp] * JxW[qp];
                    double F_n_qp = F_n * phi[k][qp] * JxW[qp];
                    double H_qp = HK * phi[k][qp] * JxW[qp];
#if (NDIM == 3)
                    dw_j_qp = dw * phi[k][qp] * JxW[qp];
		    double F_b_qp = F_b * phi[k][qp] * JxW[qp];
#endif
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        F_rhs_e[i](k) += F_qp(i);
                        du_j_rhs_e[i](k) += du(i) * phi[k][qp] * JxW[qp];
                        dv_j_rhs_e[i](k) += dv(i) * phi[k][qp] * JxW[qp];
#if (NDIM == 3)
                        dw_j_rhs_e[i](k) += dw(i) * phi[k][qp] * JxW[qp];
#endif
                    }
                        
                    double P_j_qp = P_j * phi[k][qp] * JxW[qp];
                    P_j_rhs_e(k) += P_j_qp;
                    F_t_rhs_e(k) += F_t_qp;
					F_n_rhs_e(k) += F_n_qp;
					H_rhs_e(k) += H_qp;
#if (NDIM == 3)
					F_b_rhs_e(k) += F_b_qp;
#endif
                }
            }
        }

        // Apply constraints (e.g., enforce periodic boundary conditions)
        // and add the elemental contributions to the global vector.
        for (unsigned int i = 0; i < NDIM; ++i)
        {
            F_dof_map.constrain_element_vector(F_rhs_e[i], F_dof_indices[i]);
            F_rhs_vec->add_vector(F_rhs_e[i], F_dof_indices[i]);
	        F_vec.add_vector(F_rhs_e[i], F_dof_indices[i]);

            du_j_dof_map.constrain_element_vector(du_j_rhs_e[i], du_j_dof_indices[i]);
            du_j_rhs_vec->add_vector(du_j_rhs_e[i], du_j_dof_indices[i]);

            dv_j_dof_map.constrain_element_vector(dv_j_rhs_e[i], dv_j_dof_indices[i]);
            dv_j_rhs_vec->add_vector(dv_j_rhs_e[i], dv_j_dof_indices[i]);

#if (NDIM == 3)
            dw_j_dof_map.constrain_element_vector(dw_j_rhs_e[i], dw_j_dof_indices[i]);
            dw_j_rhs_vec->add_vector(dw_j_rhs_e[i], dw_j_dof_indices[i]);
#endif
        }
        P_j_dof_map.constrain_element_vector(P_j_rhs_e, P_j_dof_indices);
        P_j_rhs_vec->add_vector(P_j_rhs_e, P_j_dof_indices);

	    F_t_dof_map.constrain_element_vector(F_t_rhs_e, F_t_dof_indices);
        F_t_rhs_vec->add_vector(F_t_rhs_e, F_t_dof_indices);
        F_t_vec.add_vector(F_t_rhs_e, F_t_dof_indices);
        
        H_dof_map.constrain_element_vector(H_rhs_e, H_dof_indices);
        H_rhs_vec->add_vector(H_rhs_e, H_dof_indices);
        H_vec.add_vector(H_rhs_e, H_dof_indices);
        
#if (NDIM == 3)
	    F_b_dof_map.constrain_element_vector(F_b_rhs_e, F_b_dof_indices);
        F_b_rhs_vec->add_vector(F_b_rhs_e, F_b_dof_indices);
        F_b_vec.add_vector(F_b_rhs_e, F_b_dof_indices);

#endif
	    F_n_dof_map.constrain_element_vector(F_n_rhs_e, F_n_dof_indices);
        F_n_rhs_vec->add_vector(F_n_rhs_e, F_n_dof_indices);
        F_n_vec.add_vector(F_n_rhs_e, F_n_dof_indices);
    }

    // Solve for F, [[p]] and [[dU/dn]].
    
        d_fe_data_managers[part]->computeL2Projection(F_vec, *F_rhs_vec, FORCE_SYSTEM_NAME, d_use_consistent_mass_matrix);

        d_fe_data_managers[part]->computeL2Projection(F_t_vec, *F_t_rhs_vec, FORCE_T_SYSTEM_NAME, d_use_consistent_mass_matrix);
        
        d_fe_data_managers[part]->computeL2Projection(H_vec, *H_rhs_vec, H_SYSTEM_NAME, d_use_consistent_mass_matrix);
        
        d_fe_data_managers[part]->computeL2Projection(F_n_vec, *F_n_rhs_vec, FORCE_N_SYSTEM_NAME, d_use_consistent_mass_matrix);

        //~ if (d_split_normal_force)
        //~ {
        d_fe_data_managers[part]->computeL2Projection(
            P_j_vec, *P_j_rhs_vec, P_J_SYSTEM_NAME, d_use_consistent_mass_matrix);
        //~ }
        //~ if (d_split_tangential_force)
        //~ {
        d_fe_data_managers[part]->computeL2Projection(
            du_j_vec, *du_j_rhs_vec, DU_J_SYSTEM_NAME, d_use_consistent_mass_matrix);
        d_fe_data_managers[part]->computeL2Projection(
            dv_j_vec, *dv_j_rhs_vec, DV_J_SYSTEM_NAME, d_use_consistent_mass_matrix);
        //~ }

#if (NDIM == 3)
    d_fe_data_managers[part]->computeL2Projection(dw_j_vec, *dw_j_rhs_vec, DW_J_SYSTEM_NAME, d_use_consistent_mass_matrix);
    d_fe_data_managers[part]->computeL2Projection(F_b_vec, *F_b_rhs_vec, FORCE_B_SYSTEM_NAME, d_use_consistent_mass_matrix);
        
#endif

    PetscVector<double>* F_t_petsc_vec = static_cast<PetscVector<double>*>(&F_t_vec);
    Vec F_t_global_vec = F_t_petsc_vec->vec();
    Vec F_t_local_vec;
    VecGhostGetLocalForm(F_t_global_vec, &F_t_local_vec);
    double* F_t_local_soln;
    VecGetArray(F_t_local_vec, &F_t_local_soln);


    PetscVector<double>* F_n_petsc_vec = static_cast<PetscVector<double>*>(&F_n_vec);
    Vec F_n_global_vec = F_n_petsc_vec->vec();
    Vec F_n_local_vec;
    VecGhostGetLocalForm(F_n_global_vec, &F_n_local_vec);
    double* F_n_local_soln;
    VecGetArray(F_n_local_vec, &F_n_local_soln);
    
    PetscVector<double>* H_petsc_vec = static_cast<PetscVector<double>*>(&H_vec);
    Vec H_global_vec = H_petsc_vec->vec();
    Vec H_local_vec;
    VecGhostGetLocalForm(H_global_vec, &H_local_vec);
    double* H_local_soln;
    VecGetArray(H_local_vec, &H_local_soln);

#if (NDIM == 3)
    PetscVector<double>* F_b_petsc_vec = static_cast<PetscVector<double>*>(&F_b_vec);
    Vec F_b_global_vec = F_b_petsc_vec->vec();
    Vec F_b_local_vec;
    VecGhostGetLocalForm(F_b_global_vec, &F_b_local_vec);
    double* F_b_local_soln;
    VecGetArray(F_b_local_vec, &F_b_local_soln);
#endif



	//~ // Calculate the second order jumps [[dp/dx]], [[dp/dy]], [[du/x^]], [[du/dx^2]],[[dv/x^]], [[dv/dx^2]]
	
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;


        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map_cache.dof_indices(elem, X_dof_indices[d], d);
            X0_dof_map_cache.dof_indices(elem, X0_dof_indices[d], d);
            dP_j_dof_map_cache.dof_indices(elem, dP_j_dof_indices[d], d);
            dP_j_rhs_e[d].resize(static_cast<int>(dP_j_dof_indices[d].size()));
            
            d2u_j_dof_map_cache.dof_indices(elem, d2u_j_dof_indices[d], d);
            d2u_j_rhs_e[d].resize(static_cast<int>(d2u_j_dof_indices[d].size()));
            
            
            d2v_j_dof_map_cache.dof_indices(elem, d2v_j_dof_indices[d], d);
            d2v_j_rhs_e[d].resize(static_cast<int>(d2v_j_dof_indices[d].size()));  
#if (NDIM == 3)
            d2w_j_dof_map_cache.dof_indices(elem, d2w_j_dof_indices[d], d);
            d2w_j_rhs_e[d].resize(static_cast<int>(d2w_j_dof_indices[d].size()));
#endif
        }
        F_t_dof_map_cache.dof_indices(elem, F_t_dof_indices);
        F_n_dof_map_cache.dof_indices(elem, F_n_dof_indices);
#if (NDIM == 3)
	   F_b_dof_map_cache.dof_indices(elem, F_b_dof_indices);
#endif

       
        
        
        fe.reinit(elem);
        fe.collectDataForInterpolation(elem);
        fe.interpolate(elem);
        X_fe_base->reinit(elem);
        get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);
        get_values_for_interpolation(X0_node, *X0_petsc_vec, X0_local_soln, X0_dof_indices);
        
        
        
        get_values_for_interpolation(F_t_node, *F_t_petsc_vec, F_t_local_soln, F_t_dof_indices);
        
        get_values_for_interpolation(F_n_node, *F_n_petsc_vec, F_n_local_soln, F_n_dof_indices);
#if (NDIM == 3)
        get_values_for_interpolation(F_b_node, *F_b_petsc_vec, F_b_local_soln, F_b_dof_indices);
#endif         
        get_values_for_interpolation(H_node, *H_petsc_vec, H_local_soln, H_dof_indices);
        
        const unsigned int n_qp = qrule->n_points();
        const size_t n_basis = phi.size();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            s = q_point[qp];

            interpolate(&X(0), qp, X0_node, X_phi);
            interpolate(&Tau1(0), qp, X0_node, X_dphi_dxi);
            if (dim == 1)
                Tau2 = VectorValue<double>(0.0, 0.0, 1.0);
            else
                interpolate(&Tau2(0), qp, X0_node, X_dphi_deta);
            N = Tau1.cross(Tau2);

            interpolate(&x(0), qp, X_node, X_phi);
            interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
            
			tau1unit = tau1.unit();
			
            if (dim == 1)
                tau2 = VectorValue<double>(0.0, 0.0, 1.0);
            else
				interpolate(&tau2(0), qp, X_node, X_dphi_deta);
            n = tau1.cross(tau2);
           

            dA_da = sqrt(N * N) / sqrt(n * n);

            N = N.unit();
            n = n.unit();

            double* nn;

            nn = &n(0);

            if (d_lag_force_fcn_data[part].fcn)
            {
                //~ // Compute the value of the body force at the quadrature
                //~ // point and add the corresponding forces to the
                //~ // right-hand-side vector.


			//~ //	interpolate(&kappa(0), qp, X_node , X_curvature);
				
				// 0.5 * (kappa(0) + kappa(1));
				// pout << " H = "<<  X_curvature[qp] <<"\n\n";
				
				double H;
				H = interpolate(qp, H_node, X_phi);
				 // H = interpolate(qp, X_node , X_curvature);
			
			
               
                const double dP_dt_j = dA_da * interpolate(qp, F_t_node, X_dphi_dxi);
                
                
                const double Ft = dA_da * interpolate(qp, F_t_node, X_phi);
#if (NDIM == 3)
				const double Fb = dA_da * interpolate(qp, F_b_node, X_phi);
				const double dP_db_j = dA_da * interpolate(qp, F_b_node,  X_dphi_deta);
#endif                
                //~ const double dP_dn_j = dA_da * interpolate(qp, F_n_node, X_dphi_dxi);
                
                dP_j  =  dP_dt_j * n + dA_da * interpolate(qp, F_n_node, X_dphi_dxi)* tau1unit;
                //dP_dn_j*n + dP_dt_j*tau1unit;
#if (NDIM == 3)
				 dP_j  =  (dP_dt_j + dP_db_j) * n + dA_da * interpolate(qp, F_n_node, X_dphi_dxi)* tau1unit + dA_da * interpolate(qp, F_n_node, X_dphi_deta)* tau2unit;
				//dP_dn_j*n + dP_dt_j*tau1unit + dP_db_j*tau2unit;		
#endif 
                
                
                d2u_j = (dP_dt_j * n + dA_da * interpolate(qp, F_n_node, X_dphi_dxi)* tau1unit - H * Ft * tau1unit )* n(0) * n(0) + 2.0 * (-dP_dt_j * tau1unit - H * Ft * n ) * n(0) * tau1unit(0) + ( H * Ft * tau1unit ) * tau1unit(0) * tau1unit(0);   
                d2v_j = (dP_dt_j * n + dA_da * interpolate(qp, F_n_node, X_dphi_dxi)* tau1unit - H * Ft * tau1unit )* n(1) * n(1) + 2.0 * (-dP_dt_j * tau1unit - H * Ft * n ) * n(1) * tau1unit(1) + ( H * Ft * tau1unit ) * tau1unit(1) * tau1unit(1);
 #if (NDIM == 3)
                 d2w_j = (dP_j - H * Fb * tau1unit )* n(2) * n(2) + 2.0 * (-dP_dt_j * tau2unit - H * Fb * n ) * n(2) * tau1unit(2) + ( H * Fb * tau2unit ) * tau1unit(2) * tau1unit(2);

 #endif               
                
            
                //~ // Accumulate the RHS values.
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    dP_j_qp = dP_j * phi[k][qp] * JxW[qp];
                    
                    d2u_j_qp = d2u_j * phi[k][qp] * JxW[qp];
                    
                    d2v_j_qp = d2v_j * phi[k][qp] * JxW[qp];
#if (NDIM == 3)
					d2w_j_qp = d2w_j * phi[k][qp] * JxW[qp];
#endif
                     
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        dP_j_rhs_e[i](k) += dP_j_qp(i)* n(i);
                        d2u_j_rhs_e[i](k) += d2u_j_qp(i)* n(i);
                        d2v_j_rhs_e[i](k) += d2v_j_qp(i)* n(i);
#if (NDIM == 3)
						d2w_j_rhs_e[i](k) += d2w_j_qp(i)* n(i);
#endif
                    }

                }
            }
        }

        // Apply constraints (e.g., enforce periodic boundary conditions)
        // and add the elemental contributions to the global vector.
        for (unsigned int i = 0; i < NDIM; ++i)
        {
            dP_j_dof_map.constrain_element_vector(dP_j_rhs_e[i], dP_j_dof_indices[i]);
            dP_j_rhs_vec->add_vector(dP_j_rhs_e[i], dP_j_dof_indices[i]);
            
            d2u_j_dof_map.constrain_element_vector(d2u_j_rhs_e[i], d2u_j_dof_indices[i]);
            d2u_j_rhs_vec->add_vector(d2u_j_rhs_e[i], d2u_j_dof_indices[i]);
            
            
            d2v_j_dof_map.constrain_element_vector(d2v_j_rhs_e[i], d2v_j_dof_indices[i]);
            d2v_j_rhs_vec->add_vector(d2v_j_rhs_e[i], d2v_j_dof_indices[i]);
            
#if (NDIM == 3)
            d2w_j_dof_map.constrain_element_vector(d2w_j_rhs_e[i], d2w_j_dof_indices[i]);
            d2w_j_rhs_vec->add_vector(d2w_j_rhs_e[i], d2w_j_dof_indices[i]);
#endif


        }
    }





        d_fe_data_managers[part]->computeL2Projection(
            dP_j_vec, *dP_j_rhs_vec, DP_J_SYSTEM_NAME, d_use_consistent_mass_matrix);
            
            
       d_fe_data_managers[part]->computeL2Projection(
            d2u_j_vec, *d2u_j_rhs_vec, D2U_J_SYSTEM_NAME, d_use_consistent_mass_matrix);
            
            
      d_fe_data_managers[part]->computeL2Projection(
            d2v_j_vec, *d2v_j_rhs_vec, D2V_J_SYSTEM_NAME, d_use_consistent_mass_matrix);

#if (NDIM == 3)
      d_fe_data_managers[part]->computeL2Projection(
            d2w_j_vec, *d2w_j_rhs_vec, D2W_J_SYSTEM_NAME, d_use_consistent_mass_matrix);

#endif








    //~ VecRestoreArray(F_n_local_vec, &F_n_local_soln);
    //~ VecGhostRestoreLocalForm(F_n_global_vec, &F_n_local_vec);
    //~ VecRestoreArray(F_n_local_vec, &F_n_local_soln);
    //~ VecGhostRestoreLocalForm(F_n_global_vec, &F_n_local_vec);

    //~ VecRestoreArray(F_t_local_vec, &F_t_local_soln);
    //~ VecGhostRestoreLocalForm(F_t_global_vec, &F_t_local_vec);
    //~ VecRestoreArray(F_t_local_vec, &F_t_local_soln);
    //~ VecGhostRestoreLocalForm(F_t_global_vec, &F_t_local_vec);



    VecRestoreArray(X_local_vec, &X_local_soln);
    VecGhostRestoreLocalForm(X_global_vec, &X_local_vec);
    VecRestoreArray(X0_local_vec, &X0_local_soln);
    VecGhostRestoreLocalForm(X0_global_vec, &X0_local_vec);
    return;
} // computeInteriorForceDensity

void
IBFEMethod::imposeJumpConditionsWeak(const int f_data_idx,
                                 const int mask_current_idx,
                                 PetscVector<double>& /*F_ghost_vec*/,
                                 PetscVector<double>& X_ghost_vec,
                                 PetscVector<double>& P_j_ghost_vec,
                                 PetscVector<double>& dP_j_ghost_vec,
                                 PetscVector<double>& du_j_ghost_vec,
                                 PetscVector<double>& dv_j_ghost_vec,
                                 PetscVector<double>& dw_j_ghost_vec,  
                                 PetscVector<double>& d2u_j_ghost_vec,
                                 PetscVector<double>& d2v_j_ghost_vec,
								 PetscVector<double>& d2w_j_ghost_vec,
                                 const double data_time,
                                 const unsigned int part)

{
    if (!d_split_normal_force && !d_split_tangential_force) return;

    // Check to see if we need to integrate the normal surface force.
    const bool integrate_normal_force = d_split_normal_force && d_use_jump_conditions;
    const bool integrate_tangential_force = d_split_tangential_force && d_use_jump_conditions;
    // if (!integrate_normal_force) return;

    // Extract the mesh.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    TBOX_ASSERT(dim == NDIM - 1);

    std::vector<const std::vector<double>*> force_var_data;
    std::vector<const std::vector<VectorValue<double> >*> force_grad_var_data;

    System& F_system = equation_systems->get_system(FORCE_SYSTEM_NAME);
    const DofMap& F_dof_map = F_system.get_dof_map();
    FEDataManager::SystemDofMapCache& F_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(FORCE_SYSTEM_NAME);
    FEType F_fe_type = F_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(F_dof_map.variable_type(d) == F_fe_type);
    }
    std::vector<std::vector<unsigned int> > F_dof_indices(NDIM);

    std::vector<int> vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
    FEDataInterpolation fe(dim, d_fe_data_managers[part]);
    AutoPtr<QBase> qrule = QBase::build(d_default_quad_type[part], dim, d_default_quad_order[part]);
    fe.attachQuadratureRule(qrule.get());
    fe.evalQuadraturePoints();
    fe.evalQuadratureWeights();
    fe.registerSystem(F_system, vars, vars); // compute phi and dphi for the force system
    std::vector<size_t> force_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        force_fcn_system_idxs, d_lag_force_fcn_data[part].system_data, equation_systems);
    fe.init(/*use_IB_ghosted_vecs*/ false);

    // Extract the FE systems and DOF maps, and setup the FE object.
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& X_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();

    FEType X_fe_type = X_dof_map.variable_type(0);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
    }
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);


    // Setup global and elemental right-hand-side vectors.

    //~ if (integrate_normal_force)
    //~ {
    System& P_j_system = equation_systems->get_system(P_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& P_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(P_J_SYSTEM_NAME);
    DofMap& P_j_dof_map = P_j_system.get_dof_map();
    TBOX_ASSERT(P_j_dof_map.variable_type(0) == X_fe_type);
    std::vector<unsigned int> P_j_dof_indices;
    //~ }

    System& du_j_system = equation_systems->get_system(DU_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& du_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DU_J_SYSTEM_NAME);
    const DofMap& du_j_dof_map = du_j_system.get_dof_map();

    FEType du_j_fe_type = du_j_dof_map.variable_type(0);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(du_j_dof_map.variable_type(d) == du_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > du_j_dof_indices(NDIM);

    System& dv_j_system = equation_systems->get_system(DV_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& dv_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DV_J_SYSTEM_NAME);
    const DofMap& dv_j_dof_map = dv_j_system.get_dof_map();
    FEType dv_j_fe_type = dv_j_dof_map.variable_type(0);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(dv_j_dof_map.variable_type(d) == dv_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > dv_j_dof_indices(NDIM);


    System& dw_j_system = equation_systems->get_system(DW_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& dw_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DW_J_SYSTEM_NAME);
    const DofMap& dw_j_dof_map = dw_j_system.get_dof_map();
    FEType dw_j_fe_type = dw_j_dof_map.variable_type(0);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(dw_j_dof_map.variable_type(d) == dw_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > dw_j_dof_indices(NDIM);


    System& d2u_j_system = equation_systems->get_system(D2U_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& d2u_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(D2U_J_SYSTEM_NAME);
    const DofMap& d2u_j_dof_map = d2u_j_system.get_dof_map();

    FEType d2u_j_fe_type = d2u_j_dof_map.variable_type(0);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d2u_j_dof_map.variable_type(d) == d2u_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > d2u_j_dof_indices(NDIM);

    System& d2v_j_system = equation_systems->get_system(D2V_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& d2v_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(D2V_J_SYSTEM_NAME);
    const DofMap& d2v_j_dof_map = d2v_j_system.get_dof_map();
    FEType d2v_j_fe_type = d2v_j_dof_map.variable_type(0);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d2v_j_dof_map.variable_type(d) == d2v_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > d2v_j_dof_indices(NDIM);

    System& d2w_j_system = equation_systems->get_system(D2W_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& d2w_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(D2W_J_SYSTEM_NAME);
    const DofMap& d2w_j_dof_map = d2w_j_system.get_dof_map();
    FEType d2w_j_fe_type = d2w_j_dof_map.variable_type(0);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d2w_j_dof_map.variable_type(d) == d2w_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > d2w_j_dof_indices(NDIM);
    
    
    
    System& dP_j_system = equation_systems->get_system(DP_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& dP_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DP_J_SYSTEM_NAME);
    const DofMap& dP_j_dof_map = dP_j_system.get_dof_map();
    FEType dP_j_fe_type = dP_j_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(dP_j_dof_map.variable_type(d) == dP_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > dP_j_dof_indices(NDIM);


    AutoPtr<FEBase> X_fe_base(FEBase::build(dim, X_fe_type));
    const std::vector<std::vector<double> >& X_phi = X_fe_base->get_phi();
    const std::vector<std::vector<double> >& X_dphi_dxi = X_fe_base->get_dphidxi();
    const std::vector<std::vector<double> >& X_dphi_deta = X_fe_base->get_dphideta();

    PetscVector<double>* X_petsc_vec = static_cast<PetscVector<double>*>(&X_ghost_vec);
    Vec X_global_vec = X_petsc_vec->vec();
    Vec X_local_vec;
    VecGhostGetLocalForm(X_global_vec, &X_local_vec);
    double* X_local_soln;
    VecGetArray(X_local_vec, &X_local_soln);


    //~ if (integrate_normal_force)
    //~ {
    PetscVector<double>* P_j_petsc_vec = static_cast<PetscVector<double>*>(&P_j_ghost_vec);
    Vec P_j_global_vec = P_j_petsc_vec->vec();
    Vec P_j_local_vec;
    VecGhostGetLocalForm(P_j_global_vec, &P_j_local_vec);
    double* P_j_local_soln;
    VecGetArray(P_j_local_vec, &P_j_local_soln);
    
    PetscVector<double>* dP_j_petsc_vec = static_cast<PetscVector<double>*>(&dP_j_ghost_vec);
    Vec dP_j_global_vec = dP_j_petsc_vec->vec();
    Vec dP_j_local_vec;
    VecGhostGetLocalForm(dP_j_global_vec, &dP_j_local_vec);
    double* dP_j_local_soln;
    VecGetArray(dP_j_local_vec, &dP_j_local_soln);
    //~ }
    PetscVector<double>* du_j_petsc_vec = static_cast<PetscVector<double>*>(&du_j_ghost_vec);
    Vec du_j_global_vec = du_j_petsc_vec->vec();
    Vec du_j_local_vec;
    VecGhostGetLocalForm(du_j_global_vec, &du_j_local_vec);
    double* du_j_local_soln;
    VecGetArray(du_j_local_vec, &du_j_local_soln);

    PetscVector<double>* dv_j_petsc_vec = static_cast<PetscVector<double>*>(&dv_j_ghost_vec);
    Vec dv_j_global_vec = dv_j_petsc_vec->vec();
    Vec dv_j_local_vec;
    VecGhostGetLocalForm(dv_j_global_vec, &dv_j_local_vec);
    double* dv_j_local_soln;
    VecGetArray(dv_j_local_vec, &dv_j_local_soln);

    PetscVector<double>* dw_j_petsc_vec = static_cast<PetscVector<double>*>(&dw_j_ghost_vec);
    Vec dw_j_global_vec = dw_j_petsc_vec->vec();
    Vec dw_j_local_vec;
    VecGhostGetLocalForm(dw_j_global_vec, &dw_j_local_vec);
    double* dw_j_local_soln;
    VecGetArray(dw_j_local_vec, &dw_j_local_soln);

    PetscVector<double>* d2u_j_petsc_vec = static_cast<PetscVector<double>*>(&d2u_j_ghost_vec);
    Vec d2u_j_global_vec = d2u_j_petsc_vec->vec();
    Vec d2u_j_local_vec;
    VecGhostGetLocalForm(d2u_j_global_vec, &d2u_j_local_vec);
    double* d2u_j_local_soln;
    VecGetArray(d2u_j_local_vec, &d2u_j_local_soln);

    PetscVector<double>* d2v_j_petsc_vec = static_cast<PetscVector<double>*>(&d2v_j_ghost_vec);
    Vec d2v_j_global_vec = d2v_j_petsc_vec->vec();
    Vec d2v_j_local_vec;
    VecGhostGetLocalForm(d2v_j_global_vec, &d2v_j_local_vec);
    double* d2v_j_local_soln;
    VecGetArray(d2v_j_local_vec, &d2v_j_local_soln);

    PetscVector<double>* d2w_j_petsc_vec = static_cast<PetscVector<double>*>(&d2w_j_ghost_vec);
    Vec d2w_j_global_vec = d2w_j_petsc_vec->vec();
    Vec d2w_j_local_vec;
    VecGhostGetLocalForm(d2w_j_global_vec, &d2w_j_local_vec);
    double* d2w_j_local_soln;
    VecGetArray(d2w_j_local_vec, &d2w_j_local_soln);
    
    
    const std::vector<std::vector<double> >& phi = fe.getPhi(F_fe_type);

    // Loop over the patches to impose jump conditions on the Eulerian grid that
    // are determined from the interior and transmission elastic force
    // densities.
    const std::vector<std::vector<Elem*> >& active_patch_element_map =
        d_fe_data_managers[part]->getActivePatchElementMap();
    const int level_num = d_fe_data_managers[part]->getLevelNumber();
    std::vector<libMesh::Point> s_node_cache, x_node_cache;
    TensorValue<double> FF;
    VectorValue<double> F, n, s, N, tau1, tau2, Tau2, Tau1, x, X, jn, jnn;
    IBTK::Point x_min, x_max;
    //~ if (integrate_normal_force)
    boost::multi_array<double, 1> P_j_node;
    boost::multi_array<double, 2> du_j_node, dv_j_node, dw_j_node, dP_j_node;

    boost::multi_array<double, 2> d2u_j_node, d2v_j_node, d2w_j_node;

    std::vector<std::vector<unsigned int> > side_dof_indices(NDIM);

    //~ if (integrate_normal_force)
    //~ {
    std::vector<libMesh::Point> intersection_ref_coords_p;
    std::vector<libMesh::Point> intersectionSide_ref_coords_u;
    std::vector<libMesh::Point> intersectionSide2_ref_coords_u;
    std::vector<libMesh::Point> intersection_ref_coords_u;
    std::vector<libMesh::Point> intersection_coords_u;
    std::vector<libMesh::Point> intersection_coords_p;
    std::vector<libMesh::Point> intersectionSide_coords_u;
    std::vector<libMesh::Point> intersectionSide2_coords_u;

    std::vector<SideIndex<NDIM> > intersection_indices_p;
    //~ }
    //~ if (integrate_tangential_force)
    //~ {
    std::vector<SideIndex<NDIM> > intersection_indices_um;

    std::vector<SideIndex<NDIM> > intersection_indices_up;

    std::vector<SideIndex<NDIM> > intersectionSide_indices_up;

    std::vector<SideIndex<NDIM> > intersectionSide_indices_um;

    //~ #if(NDIM == 3)
    std::vector<SideIndex<NDIM> > intersectionSide2_indices_up;

    std::vector<SideIndex<NDIM> > intersectionSide2_indices_um;
    //~ #endif
    //~ }

    std::vector<std::pair<double, libMesh::Point> > intersections;
    std::vector<std::pair<double, libMesh::Point> > intersectionsSide;
    std::vector<std::pair<double, libMesh::Point> > intersectionsSide2;

    boost::multi_array<double, 2> X_node;

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
        Pointer<SideData<NDIM, double> > mask_current = patch->getPatchData(mask_current_idx);

        const Box<NDIM>& patch_box = patch->getBox();
        Box<NDIM> side_boxes[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            side_boxes[d] = SideGeometry<NDIM>::toSideBox(patch_box, d);
        }
        const CellIndex<NDIM>& patch_lower = patch_box.lower();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const x_lower = patch_geom->getXLower();
        const double* const x_upper = patch_geom->getXUpper();
        const double* const dx = patch_geom->getDx();

        //~ if (integrate_normal_force)
        //~ {
        SideData<NDIM, int> num_intersections_p(patch_box, 1, IntVector<NDIM>(0));
        num_intersections_p.fillAll(0);
        //~ }

        //~ if (integrate_tangential_force)
        //~ {
        SideData<NDIM, int> num_intersections_um(patch_box, 1, IntVector<NDIM>(0));
        num_intersections_um.fillAll(0);

        SideData<NDIM, int> num_intersections_up(patch_box, 1, IntVector<NDIM>(0));
        num_intersections_up.fillAll(0);

        SideData<NDIM, int> num_intersectionsSide_um(patch_box, 1, IntVector<NDIM>(1));
        num_intersectionsSide_um.fillAll(0);

        SideData<NDIM, int> num_intersectionsSide_up(patch_box, 1, IntVector<NDIM>(1));
        num_intersectionsSide_up.fillAll(0);

        SideData<NDIM, int> num_intersectionsSide2_um(patch_box, 1, IntVector<NDIM>(1));
        num_intersectionsSide2_um.fillAll(0);

        SideData<NDIM, int> num_intersectionsSide2_up(patch_box, 1, IntVector<NDIM>(1));
        num_intersectionsSide2_up.fillAll(0);
        //~ }

        // Loop over the elements.
        for (size_t e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const unsigned int n_node = elem->n_nodes();

            for (int d = 0; d < NDIM; ++d)
            {
                F_dof_map_cache.dof_indices(elem, F_dof_indices[d], d);
                X_dof_map_cache.dof_indices(elem, X_dof_indices[d], d);
                du_j_dof_map_cache.dof_indices(elem, du_j_dof_indices[d], d);
                dv_j_dof_map_cache.dof_indices(elem, dv_j_dof_indices[d], d);
				dw_j_dof_map_cache.dof_indices(elem, dw_j_dof_indices[d], d);
                d2u_j_dof_map_cache.dof_indices(elem, d2u_j_dof_indices[d], d);
                d2v_j_dof_map_cache.dof_indices(elem, d2v_j_dof_indices[d], d);
				d2w_j_dof_map_cache.dof_indices(elem, d2w_j_dof_indices[d], d);
                dP_j_dof_map_cache.dof_indices(elem, dP_j_dof_indices[d], d);
            }

            get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);

            if (integrate_normal_force)
            {
				P_j_dof_map_cache.dof_indices(elem, P_j_dof_indices);
                get_values_for_interpolation(P_j_node, *P_j_petsc_vec, P_j_local_soln, P_j_dof_indices);
            }
            get_values_for_interpolation(du_j_node, *du_j_petsc_vec, du_j_local_soln, du_j_dof_indices);
            get_values_for_interpolation(dv_j_node, *dv_j_petsc_vec, dv_j_local_soln, dv_j_dof_indices);
            get_values_for_interpolation(dw_j_node, *dw_j_petsc_vec, dw_j_local_soln, dw_j_dof_indices);

            get_values_for_interpolation(d2u_j_node, *d2u_j_petsc_vec, d2u_j_local_soln, d2u_j_dof_indices);

            get_values_for_interpolation(d2v_j_node, *d2v_j_petsc_vec, d2v_j_local_soln, d2v_j_dof_indices);
            
            get_values_for_interpolation(d2w_j_node, *d2w_j_petsc_vec, d2w_j_local_soln, d2w_j_dof_indices);

            get_values_for_interpolation(dP_j_node, *dP_j_petsc_vec, dP_j_local_soln, dP_j_dof_indices);

            // Cache the nodal and physical coordinates of the side element,
            // determine the bounding box of the current configuration of
            // the side element, and set the nodal coordinates to correspond
            // to the physical coordinates.
            s_node_cache.resize(n_node);
            x_node_cache.resize(n_node);
            x_min = IBTK::Point::Constant(std::numeric_limits<double>::max());
            x_max = IBTK::Point::Constant(-std::numeric_limits<double>::max());
            for (unsigned int k = 0; k < n_node; ++k)
            {
                s_node_cache[k] = elem->point(k);
                libMesh::Point& x = x_node_cache[k];
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    x(d) = X_ghost_vec(X_dof_indices[d][k]);
                    x_min[d] = std::min(x_min[d], x(d));
                    x_max[d] = std::max(x_max[d], x(d));
                }
                elem->point(k) = x;
            }
            Box<NDIM> box(IndexUtilities::getCellIndex(&x_min[0], grid_geom, ratio),
                          IndexUtilities::getCellIndex(&x_max[0], grid_geom, ratio));
            box.grow(IntVector<NDIM>(1));
            box = box * patch_box;

            //~ Box<NDIM> side_boxes[NDIM];
            //~ for (int d = 0; d < NDIM; ++d)
            //~ {
            //~ side_boxes[d] = SideGeometry<NDIM>::toSideBox(box, d);
            //~ }

            // Loop over coordinate directions and look for intersections
            // with the background fluid grid.
            if (integrate_normal_force)
            {
                intersection_ref_coords_p.clear();
                intersection_coords_p.clear();
                intersection_indices_p.clear();
            }

            if (integrate_tangential_force)
            {
                intersection_ref_coords_u.clear();
                intersection_coords_u.clear();
                intersectionSide_ref_coords_u.clear();
                intersectionSide2_ref_coords_u.clear();
                intersectionSide_coords_u.clear();
                intersectionSide2_coords_u.clear();

                intersection_indices_um.clear();
                intersection_indices_up.clear();
                intersectionSide_indices_um.clear();
                intersectionSide_indices_up.clear();

                intersectionSide2_indices_um.clear();
                intersectionSide2_indices_up.clear();
            }
            VectorValue<double> DX;
            
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
				 DX(axis) = dx[axis];
			}

            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                // Setup a unit vector pointing in the coordinate direction of interest.
                VectorValue<double> q;
                q(axis) = 1.0;
                

                // Loop over the relevant range of indices.
                Box<NDIM> axis_box = box;
                axis_box.lower(axis) = 0;
                axis_box.upper(axis) = 0;
                for (BoxIterator<NDIM> b(axis_box); b; b++)
                {
                    const Index<NDIM>& i_c = b();
                    libMesh::Point r;
                    libMesh::Point rs;
                    libMesh::Point rss;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        r(d) = (d == axis ? 0.0 : x_lower[d] +
                                                      dx[d] * (static_cast<double>(i_c(d) - patch_lower[d]) +
                                                               0.5)); // In 2D this corresponds to applying jumps for
                                                                      // [p],[ux] and [vy]
#if (NDIM == 2)
                        rs(d) = (d == axis ?
                                     0.0 :
                                     x_lower[d] + dx[d] * (static_cast<double>(i_c(d) - patch_lower[d]))); // In 2D this
                        // corresponds to
                        // applying jumps
                        // for [uy] and
                        // [vx]
#endif
                    }

#if (NDIM == 3)
                    if (axis == 0)
                    {
                        rs(0) = 0.0;
                        rs(1) = x_lower[1] + dx[1] * (static_cast<double>(i_c(1) - patch_lower[1]));
                        rs(2) = x_lower[2] + dx[2] * (static_cast<double>(i_c(2) - patch_lower[2]) + 0.5);
                        rss(0) = 0.0;
                        rss(1) = x_lower[1] + dx[1] * (static_cast<double>(i_c(1) - patch_lower[1]) + 0.5);
                        rss(2) = x_lower[2] + dx[2] * (static_cast<double>(i_c(2) - patch_lower[2]));
                    }
                    else if (axis == 1)
                    {
                        rs(0) = x_lower[0] + dx[0] * (static_cast<double>(i_c(0) - patch_lower[0]) + 0.5);
                        rs(1) = 0.0;
                        rs(2) = x_lower[2] + dx[2] * (static_cast<double>(i_c(2) - patch_lower[2]));
                        rss(0) = x_lower[0] + dx[0] * (static_cast<double>(i_c(0) - patch_lower[0]));
                        rss(1) = 0.0;
                        rss(2) = x_lower[2] + dx[2] * (static_cast<double>(i_c(2) - patch_lower[2]) + 0.5);
                    }
                    else if (axis == 2)
                    {
                        rs(0) = x_lower[0] + dx[0] * (static_cast<double>(i_c(0) - patch_lower[0]));
                        rs(1) = x_lower[1] + dx[1] * (static_cast<double>(i_c(1) - patch_lower[1]) + 0.5);
                        rs(2) = 0.0;
                        rss(0) = x_lower[0] + dx[0] * (static_cast<double>(i_c(0) - patch_lower[0]) + 0.5);
                        rss(1) = x_lower[1] + dx[1] * (static_cast<double>(i_c(1) - patch_lower[1]));
                        rss(2) = 0.0;
                    }

#endif

#if (NDIM == 2)
                    intersect_line_with_edge(intersections, static_cast<Edge*>(elem), DX, r, q);
                    intersect_line_with_edge(intersectionsSide, static_cast<Edge*>(elem), DX, rs, q);
#endif
#if (NDIM == 3)
                    intersect_line_with_face(intersections, static_cast<Face*>(elem), r, q);
                    intersect_line_with_face(intersectionsSide, static_cast<Face*>(elem), rs, q);
                    intersect_line_with_face(intersectionsSide2, static_cast<Face*>(elem), rss, q);
#endif
                    if (integrate_normal_force)
                    {
                        for (unsigned int k = 0; k < intersections.size(); ++k)
                        {
                            libMesh::Point x = r + intersections[k].first * q;
                            SideIndex<NDIM> i_s(i_c, axis, 0);
                            i_s(axis) = static_cast<int>(std::floor((x(axis) - x_lower[axis]) / dx[axis] + 0.5)) +
                                        patch_lower[axis];

                            if (side_boxes[axis].contains(i_s))
                            {
                                intersection_ref_coords_p.push_back(intersections[k].second);
                                intersection_coords_p.push_back(x);
                                intersection_indices_p.push_back(i_s);
                                num_intersections_p(i_s) += 1;
                            }
                        }
                    }

                    //~ if (integrate_tangential_force && side_boxes[axis].contains(i_c))
                    if (integrate_tangential_force)
                    {
                        for (unsigned int k = 0; k < intersections.size(); ++k)
                        {
                            libMesh::Point xs = r + intersections[k].first * q;
                            SideIndex<NDIM> i_ss(i_c, axis, 0);
                            Index<NDIM> i_c_neighbor = i_c;
                            i_c_neighbor(axis) += 1;
                            
                            SideIndex<NDIM> i_sss(i_c_neighbor, axis, 0);
                            i_sss(axis) = static_cast<int>(std::floor((xs(axis) - x_lower[axis]) / dx[axis] + 1.0)) +
                                          patch_lower[axis];
                            i_ss(axis) =
                                static_cast<int>(std::floor((xs(axis) - x_lower[axis]) / dx[axis])) + patch_lower[axis];

                            const double x_cell_bdry_um =
                                x_lower[axis] + static_cast<double>(i_ss(axis) - patch_lower[axis]) * dx[axis];

                            const double x_cell_bdry_up =
                                x_lower[axis] + static_cast<double>(i_sss(axis) - patch_lower[axis]) * dx[axis];

                            const double sdh_um = ((xs(axis) - x_cell_bdry_um)); // Signed Distance h

                            const double sdh_up = ((xs(axis) - x_cell_bdry_up)); // Signed Distance h

                            TBOX_ASSERT((sdh_um) < dx[axis] && sdh_um > 0);

                            if (side_boxes[axis].contains(i_ss) && side_boxes[axis].contains(i_sss))
                            {
                                intersection_ref_coords_u.push_back(intersections[k].second);
                                intersection_coords_u.push_back(xs);
                                intersection_indices_um.push_back(i_ss);
                                num_intersections_um(i_ss) += 1;
                               
                                intersection_indices_up.push_back(i_sss);
                                num_intersections_up(i_sss) += 1;
                                
                            }
                        }


//// This will cache the point to the right(or top) of the intersection along
// the side cell

#if (NDIM == 2)

                        for (unsigned int k = 0; k < intersectionsSide.size(); ++k)
                        {
                            libMesh::Point xu = rs + intersectionsSide[k].first * q;
                            int dd = (axis == 0 ? 1 : 0);
                            if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) >= 0.5 * dx[axis])
                            {
                                SideIndex<NDIM> i_s_um(i_c, dd, 0);
                                Index<NDIM> i_c_neighbor = i_c;
                                i_c_neighbor(axis) += 1;
                             
                                SideIndex<NDIM> i_s_up(i_c_neighbor, dd, 0);
                          
                                i_s_up(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis] + 0.5)) +
                                             patch_lower[axis];             
                                i_s_um(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis])) +
                                              patch_lower[axis];
                                              
                                if (side_boxes[axis].contains(i_s_up) && side_boxes[axis].contains(i_s_um))
                                {
                                    intersectionSide_ref_coords_u.push_back(intersectionsSide[k].second);
                                    intersectionSide_coords_u.push_back(xu);
                                    intersectionSide_indices_up.push_back(i_s_up);
                                    num_intersectionsSide_up(i_s_up) += 1;
                                    intersectionSide_indices_um.push_back(i_s_um);
                                    num_intersectionsSide_um(i_s_um) += 1;
                                 
                                }
                            }
                            else if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) < 0.5 * dx[axis] &&
                                     fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) >= 0.0)
                            {
								
								SideIndex<NDIM> i_s_up(i_c, dd, 0);
								Index<NDIM> i_c_neighbor = i_c;
								i_c_neighbor(axis) -= 1;
                                SideIndex<NDIM> i_s_um(i_c_neighbor, dd, 0);
                                i_s_up(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis])) +
                                             patch_lower[axis];
                                i_s_um(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis] - 0.5)) +
											 patch_lower[axis];

                                if (side_boxes[axis].contains(i_s_up) &&  side_boxes[axis].contains(i_s_um))
                                {
                                    intersectionSide_indices_up.push_back(i_s_up);
                                    num_intersectionsSide_up(i_s_up) += 1;
                                    intersectionSide_indices_um.push_back(i_s_um);
                                    num_intersectionsSide_up(i_s_um) += 1;
                                    
                                    intersectionSide_ref_coords_u.push_back(intersectionsSide[k].second);
                                    intersectionSide_coords_u.push_back(xu);
                                }
                            }
                            else
                            {
								pout<< "improper side index found!"<<"\n\n";
						    }
                        }

       

#endif
                        int SideDim[2];
                        SideDim[0] = SideDim[1] = 0;
                        if (axis == 0)
                        {
                            SideDim[0] = 1;
                            SideDim[1] = 2;
                        }
                        else if (axis == 1)
                        {
                            SideDim[0] = 2;
                            SideDim[1] = 0;
                        }
                        else if (axis == 2)
                        {
                            SideDim[0] = 0;
                            SideDim[1] = 1;
                        }

#if (NDIM == 3)

                        for (unsigned int k = 0; k < intersectionsSide.size(); ++k)
                        {
                            libMesh::Point xu = rs + intersectionsSide[k].first * q;
                            if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) >= 0.5 * dx[axis])
                            {
                                SideIndex<NDIM> i_s_um(i_c, SideDim[0], 0);
                                Index<NDIM> i_c_neighbor = i_c;
                                i_c_neighbor(axis) += 1;

                                SideIndex<NDIM> i_s_up(i_c_neighbor, SideDim[0], 0);

                                i_s_up(axis) =
                                    static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis] + 0.5)) +
                                    patch_lower[axis];
                                i_s_um(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis])) +
                                               patch_lower[axis];

                                if (side_boxes[axis].contains(i_s_up) && side_boxes[axis].contains(i_s_um))
                                {
                                    TBOX_ASSERT(i_s_up(axis) - i_s_um(axis) == 1);
                                    intersectionSide_ref_coords_u.push_back(intersectionsSide[k].second);
                                    intersectionSide_coords_u.push_back(xu);
                                    intersectionSide_indices_up.push_back(i_s_up);
                                    num_intersectionsSide_up(i_s_up) += 1;
                                    intersectionSide_indices_um.push_back(i_s_um);
                                    num_intersectionsSide_um(i_s_um) += 1;
                                }
                            }
                            else if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) < 0.5 * dx[axis] &&
                                     fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) >= 0.0)
                            {
                                SideIndex<NDIM> i_s_up(i_c, SideDim[0], 0);
                                Index<NDIM> i_c_neighbor = i_c;
                                i_c_neighbor(axis) -= 1;
                                SideIndex<NDIM> i_s_um(i_c_neighbor, SideDim[0], 0);

                                i_s_up(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis])) +
                                               patch_lower[axis];
                                i_s_um(axis) =
                                    static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis] - 0.5)) +
                                    patch_lower[axis];

                                if (side_boxes[axis].contains(i_s_up) && side_boxes[axis].contains(i_s_um))
                                {
                                    intersectionSide_indices_up.push_back(i_s_up);
                                    num_intersectionsSide_up(i_s_up) += 1;
                                    intersectionSide_indices_um.push_back(i_s_um);
                                    num_intersectionsSide_up(i_s_um) += 1;
                                    intersectionSide_coords_u.push_back(xu);
                                    intersectionSide_ref_coords_u.push_back(intersectionsSide[k].second);
                                }
                            }
                            else
                            {
                                pout << "improper side index found!"
                                     << "\n\n";
                            }
                        }

                        for (unsigned int k = 0; k < intersectionsSide2.size(); ++k)
                        {
                            libMesh::Point xu = rss + intersectionsSide2[k].first * q;
                            if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) >= 0.5 * dx[axis])
                            {
                                SideIndex<NDIM> i_s_um(i_c, SideDim[1], 0);
                                Index<NDIM> i_c_neighbor = i_c;
                                i_c_neighbor(axis) += 1;

                                SideIndex<NDIM> i_s_up(i_c_neighbor, SideDim[1], 0);

                                i_s_up(axis) =
                                    static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis] + 0.5)) +
                                    patch_lower[axis];
                                i_s_um(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis])) +
                                               patch_lower[axis];

                                if (side_boxes[axis].contains(i_s_up) && side_boxes[axis].contains(i_s_um))
                                {
                                    intersectionSide2_ref_coords_u.push_back(intersectionsSide2[k].second);
                                    intersectionSide2_coords_u.push_back(xu);
                                    intersectionSide2_indices_up.push_back(i_s_up);
                                    num_intersectionsSide2_up(i_s_up) += 1;
                                    intersectionSide2_indices_um.push_back(i_s_um);
                                    num_intersectionsSide2_um(i_s_um) += 1;
                                }
                            }
                            else if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) < 0.5 * dx[axis] &&
                                     fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) >= 0.0)
                            {
                                SideIndex<NDIM> i_s_up(i_c, SideDim[1], 0);
                                Index<NDIM> i_c_neighbor = i_c;
                                i_c_neighbor(axis) -= 1;
                                SideIndex<NDIM> i_s_um(i_c_neighbor, SideDim[1], 0);

                                i_s_up(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis])) +
                                               patch_lower[axis];
                                i_s_um(axis) =
                                    static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis] - 0.5)) +
                                    patch_lower[axis];

                                if (side_boxes[axis].contains(i_s_up) && side_boxes[axis].contains(i_s_um))
                                {
                                    intersectionSide2_indices_up.push_back(i_s_up);
                                    num_intersectionsSide2_up(i_s_up) += 1;
                                    intersectionSide2_indices_um.push_back(i_s_um);
                                    num_intersectionsSide2_up(i_s_um) += 1;

                                    intersectionSide2_ref_coords_u.push_back(intersectionsSide2[k].second);
                                    intersectionSide2_coords_u.push_back(xu);
                                }
                            }
                            else
                            {
                                pout << "improper side index found!"
                                     << "\n\n";
                            }
                        }

#endif
                    }
                }
            }

            // Restore the element coordinates.
            for (unsigned int k = 0; k < n_node; ++k)
            {
                elem->point(k) = s_node_cache[k];
            }

            // If there are no intersection points, then continue to the next side.
            //~ if (intersection_ref_coords.empty()) continue;



            // Evaluate the jump conditions and apply them to the Eulerian grid.
            if (integrate_normal_force && !intersection_ref_coords_p.empty())
            {
                X_fe_base->attach_quadrature_rule(qrule.get());
                X_fe_base->reinit(elem, &intersection_ref_coords_p);
                //~ static const double TOL = sqrt(std::numeric_limits<double>::epsilon());
                //~ fe_face->reinit(elem, side, TOL, &intersection_ref_coords);
                const size_t n_qp = intersection_ref_coords_p.size();
                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    libMesh::Point xx = intersection_coords_p[qp];
                    interpolate(&x(0), qp, X_node, X_phi);
                    interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
                    if (dim == 1)
                        tau2 = VectorValue<double>(0.0, 0.0, 1.0);
                    else
                    interpolate(&tau2(0), qp, X_node, X_dphi_deta);
                    n = tau1.cross(tau2);
                    n = n.unit();
		
                    const SideIndex<NDIM>& i_s = intersection_indices_p[qp];
                    const unsigned int axis = i_s.getAxis();
#if !defined(NDEBUG)
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            const double x_lower_bound = x_lower[d] +
                                                         (static_cast<double>(i_s(d) - patch_lower[d]) - 0.5) * dx[d] -
                                                         sqrt(std::numeric_limits<double>::epsilon());
                            const double x_upper_bound = x_lower[d] +
                                                         (static_cast<double>(i_s(d) - patch_lower[d]) + 0.5) * dx[d] +
                                                         sqrt(std::numeric_limits<double>::epsilon());
                            TBOX_ASSERT(x_lower_bound <= xx(d) && xx(d) <= x_upper_bound);
                        }
                        else
                        {
                            const double x_intersection =
                                x_lower[d] + (static_cast<double>(i_s(d) - patch_lower[d]) + 0.5) * dx[d];
                            const double x_interp = xx(d);
                            const double rel_diff =
                                std::abs(x_intersection - x_interp) /
                                std::max(1.0, std::max(std::abs(x_intersection), std::abs(x_interp)));
                            TBOX_ASSERT(rel_diff <= sqrt(std::numeric_limits<double>::epsilon()));
                        }
                    }
#endif
                    // Impose the jump conditions.
                    
                    interpolate(&jn(0), qp, dP_j_node, X_phi);
                    
                    const double x_cell_bdry =
                                            x_lower[axis] + static_cast<double>(i_s(axis) - patch_lower[axis]) *
                                            dx[axis];
											const double h = x_cell_bdry + (x(axis) > x_cell_bdry ? +0.5 : -0.5) * dx[axis] -
											x(axis);
                    double C_p = interpolate(qp, P_j_node, X_phi);

					if (d_use_higher_order_jump) C_p -= h* jn(axis);

                    (*f_data)(i_s) += (n(axis) > 0.0 ? +1.0 : -1.0) * (C_p / dx[axis]);
                }
            }

            // Applying the jump for viscous term (in 2D u_xx, u_yy,
            // v_xx and v_yy)

            // Now use the computed weak solution and the intersection info to actually apply the jump

            if (integrate_tangential_force && !intersection_ref_coords_u.empty())
            {

                // Evaluate the jump conditions and apply them to the Eulerian
                // grid.
                // This set of intersection points will take care of the jump in the pressure
                // as well as the jump in the viscous term for the  following components
                // 1) u_xx (in the x-momentum equation)
                // 2) v_yy (in the y-momentum equation)
                // 3) w_zz (in the z-momentum equation)
				TBOX_ASSERT(intersection_indices_um.size() == intersection_indices_up.size());

                                X_fe_base->attach_quadrature_rule(qrule.get());
                X_fe_base->reinit(elem, &intersection_ref_coords_u);

                const size_t n_qp = intersection_ref_coords_u.size();

                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    interpolate(&x(0), qp, X_node, X_phi);

                    libMesh::Point xx = intersection_coords_u[qp];

                    interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
                    if (dim == 1)
                        tau2 = VectorValue<double>(0.0, 0.0, 1.0);
                    else
                        interpolate(&tau2(0), qp, X_node, X_dphi_deta);
                    n = tau1.cross(tau2);
                    n = n.unit();
					
                    const SideIndex<NDIM>& i_s_um = intersection_indices_um[qp];
                    const unsigned int axis = i_s_um.getAxis();
                    const SideIndex<NDIM>& i_s_up = intersection_indices_up[qp];

                    if (side_boxes[axis].contains(i_s_um) && side_boxes[axis].contains(i_s_up))
                    {
                        TBOX_ASSERT(i_s_um.getAxis() == i_s_up.getAxis());
                        TBOX_ASSERT(i_s_up(axis) - i_s_um(axis) == 1);
                        const double x_cell_bdry_um =
                            x_lower[axis] + static_cast<double>(i_s_um(axis) - patch_lower[axis]) * dx[axis];
                            
                        const double x_cell_bdry_up =
                            x_lower[axis] + static_cast<double>(i_s_up(axis) - patch_lower[axis]) * dx[axis];
                        const double sdh_um = ((xx(axis) - x_cell_bdry_um)); // Signed Distance h

                        const double sdh_up = ((xx(axis) - x_cell_bdry_up)); // Signed Distance h

                        TBOX_ASSERT((sdh_um) < dx[axis] && sdh_um > 0);
                        TBOX_ASSERT(fabs(sdh_up) < dx[axis] && sdh_up < 0);

                        double C_u_um = 0;
                        double C_u_up = 0;
#if (NDIM == 2)

                        if (axis == 0)
                        {
							interpolate(&jn(0), qp, du_j_node, X_phi);
							interpolate(&jnn(0), qp, d2u_j_node, X_phi);
                            C_u_up = sdh_up * jn(0);
                            C_u_um = sdh_um * jn(0);
                            if (d_use_higher_order_jump)
                            {
								C_u_up += 0.5 * sdh_up * sdh_up * jnn(0);
								C_u_um += 0.5 * sdh_um * sdh_um * jnn(0);
							}
                        }
                        else if (axis == 1)
                        {
                            interpolate(&jn(0), qp, dv_j_node, X_phi);
							interpolate(&jnn(0), qp, d2v_j_node, X_phi);
                            C_u_up = sdh_up * jn(1);
                            C_u_um = sdh_um * jn(1);
                            if (d_use_higher_order_jump)
                            {
							    C_u_up += 0.5 * sdh_up * sdh_up * jnn(1);
								C_u_um += 0.5 * sdh_um * sdh_um * jnn(1);
							}
                        }
#endif

#if (NDIM == 3)
                        if (axis == 0)
                        {
                            interpolate(&jn(0), qp, du_j_node, X_phi);
							interpolate(&jnn(0), qp, d2u_j_node, X_phi);
                            C_u_up = sdh_up * jn(0);
                            C_u_um = sdh_um * jn(0);
                            if (d_use_higher_order_jump)
                            {
							    C_u_up += 0.5 * sdh_up * sdh_up * jnn(0);
								C_u_um += 0.5 * sdh_um * sdh_um * jnn(0);
							}
                        }
                        else if (axis == 1)
                        {
                            interpolate(&jn(0), qp, dv_j_node, X_phi);
							interpolate(&jnn(0), qp, d2v_j_node, X_phi);
                            C_u_up = sdh_up * jn(1) ;
                            C_u_um = sdh_um * jn(1) ;
                            if (d_use_higher_order_jump)
                            {
							    C_u_up += 0.5 * sdh_up * sdh_up * jnn(1);
								C_u_um += 0.5 * sdh_um * sdh_um * jnn(1);
							}
                        }
                        else if (axis == 2)
                        {
                            interpolate(&jn(0), qp, dw_j_node, X_phi);
							interpolate(&jnn(0), qp, d2w_j_node, X_phi);
                            C_u_up = sdh_up * jn(2);
                            C_u_um = sdh_um * jn(2);
                            if (d_use_higher_order_jump)
                            {
							    C_u_up += 0.5 * sdh_up * sdh_up * jnn(2);
								C_u_um += 0.5 * sdh_um * sdh_um * jnn(2);
							}
                        }
#endif
                        (*f_data)(i_s_up) += - (n(axis) > 0.0 ? 1.0 : -1.0) * (C_u_um / (dx[axis] * dx[axis]));
                        (*f_data)(i_s_um) += - (n(axis) > 0.0 ? -1.0 : 1.0) * (C_u_up / (dx[axis] * dx[axis]));
                        
                    }
                }
            }

 

#if (NDIM == 2)
            // Apply the jumps of the viscous term on the direction of cell sides
            // For the velocity field U=(u,v) this will take care of corrections
            // in the following components done in a direction by direction manner:
            // 1) u_yy, u_zz  (in the x-momentum equation)
            // 2) v_xx, v_zz  (in the y-momentum equation)
           

            if (integrate_tangential_force && !intersectionSide_ref_coords_u.empty())
            {
				 TBOX_ASSERT(intersectionSide_indices_um.size() == intersectionSide_indices_up.size());
                X_fe_base->reinit(elem, &intersectionSide_ref_coords_u);
                const size_t n_qp = intersectionSide_indices_up.size();

                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    const SideIndex<NDIM>& i_side_up = intersectionSide_indices_up[qp];
                    
                    const unsigned int dd = i_side_up.getAxis();
                    
                    const unsigned int axis = (dd == 0 ? 1 : 0);

                    libMesh::Point xx = intersectionSide_coords_u[qp];

                    const SideIndex<NDIM>& i_side_um = intersectionSide_indices_um[qp];
                    interpolate(&x(0), qp, X_node, X_phi);
                    interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
                    if (dim == 1)
                        tau2 = VectorValue<double>(0.0, 0.0, 1.0);
                    else
                        interpolate(&tau2(0), qp, X_node, X_dphi_deta);
                    n = tau1.cross(tau2);
                    n = n.unit();

                    if (side_boxes[dd].contains(i_side_up) && side_boxes[dd].contains(i_side_um) )
                    {
                        // Impose the jump conditions.
                        
                        TBOX_ASSERT( i_side_up(axis) - i_side_um(axis) == 1);
                        
                        const double x_mid_side_up =
                            x_lower[axis] + static_cast<double>(i_side_up(axis) - patch_lower[axis] + 0.5) * dx[axis];
                            
                        const double x_mid_side_um =
                            x_lower[axis] + static_cast<double>(i_side_um(axis) - patch_lower[axis] + 0.5) * dx[axis];

                        TBOX_ASSERT(xx(axis) <= x_mid_side_up);
                        TBOX_ASSERT(xx(axis) > x_mid_side_um);

                        const double sdh_up = xx(axis) - x_mid_side_up; // Signed Distance h
                        const double sdh_um = xx(axis) - x_mid_side_um;

                        double C_u_um = 0;
                        double C_u_up = 0;

                        if (dd == 0)
                        {
                            interpolate(&jn(0), qp, du_j_node, X_phi);
							interpolate(&jnn(0), qp, d2u_j_node, X_phi);
                            C_u_um = sdh_um * jn(1) ;
                            C_u_up = sdh_up * jn(1) ;
                            if (d_use_higher_order_jump)
                            {
							    C_u_up += 0.5 * sdh_up * sdh_up * jnn(1);
								C_u_um += 0.5 * sdh_um * sdh_um * jnn(1);
							}
                        }
                        else
                        {
                            interpolate(&jn(0), qp, dv_j_node, X_phi);
							interpolate(&jnn(0), qp, d2v_j_node, X_phi);
                            C_u_um = sdh_um * jn(0);
                            C_u_up = sdh_up * jn(0);
                            if (d_use_higher_order_jump)
                            {
							    C_u_up += 0.5 * sdh_up * sdh_up * jnn(0);
								C_u_um += 0.5 * sdh_um * sdh_um * jnn(0);
							}
                            
                        }
						 
                       (*f_data)(i_side_um) += - (n(axis) > 0.0 ? -1.0 : 1.0) * (C_u_up / (dx[axis] * dx[axis]));
                       
                       (*f_data)(i_side_up) += - (n(axis) > 0.0 ? 1.0 : -1.0) * (C_u_um / (dx[axis] * dx[axis]));


                    }
                }
            }

#endif

#if (NDIM == 3)

            if (integrate_tangential_force && !intersectionSide_ref_coords_u.empty())
            {
                TBOX_ASSERT(intersectionSide_indices_um.size() == intersectionSide_indices_up.size());
                X_fe_base->reinit(elem, &intersectionSide_ref_coords_u);
                const size_t n_qp = intersectionSide_indices_up.size();

                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    libMesh::Point xx = intersectionSide_coords_u[qp];
                    const SideIndex<NDIM>& i_side_up = intersectionSide_indices_up[qp];
                    const SideIndex<NDIM>& i_side_um = intersectionSide_indices_um[qp];
                    const unsigned int dd = i_side_up.getAxis();
                    const unsigned int ddd = i_side_um.getAxis();
                    TBOX_ASSERT(dd == ddd);

                    interpolate(&x(0), qp, X_node, X_phi);
                    interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
                    if (dim == 1)
                        tau2 = VectorValue<double>(0.0, 0.0, 1.0);
                    else
                        interpolate(&tau2(0), qp, X_node, X_dphi_deta);
                    n = tau1.cross(tau2);
                    n = n.unit();

                    if (side_boxes[dd].contains(i_side_up) && side_boxes[dd].contains(i_side_um))
                    {
                        int axis;

                        if (dd == 0)
                        {
                            axis = 2;
                        }
                        else if (dd == 1)
                        {
                            axis = 0;
                        }
                        else
                        {
                            axis = 1;
                        }

                        // Impose the jump conditions.
                        // *********************************************************************************************************//

                        TBOX_ASSERT(i_side_up(axis) - i_side_um(axis) == 1);

                        double x_mid_side_up =
                            x_lower[axis] + static_cast<double>(i_side_up(axis) - patch_lower[axis] + 0.5) * dx[axis];

                        double x_mid_side_um =
                            x_lower[axis] + static_cast<double>(i_side_um(axis) - patch_lower[axis] + 0.5) * dx[axis];

                        TBOX_ASSERT(xx(axis) <= x_mid_side_up);
                        TBOX_ASSERT(xx(axis) > x_mid_side_um);

                        double sdh_up = xx(axis) - x_mid_side_up; // Signed Distance h
                        double sdh_um = xx(axis) - x_mid_side_um;

                        double C_u_um, C_u_up;

                        if (dd == 0)
                        {
                            interpolate(&jn(0), qp, du_j_node, X_phi);
							interpolate(&jnn(0), qp, d2u_j_node, X_phi);
                            C_u_um = sdh_um * jn(axis) ;
                            C_u_up = sdh_up * jn(axis) ;
                            if (d_use_higher_order_jump)
                            {
							    C_u_up += 0.5 * sdh_up * sdh_up * jnn(axis);
								C_u_um += 0.5 * sdh_um * sdh_um * jnn(axis);
							}
                        }
                        else if (dd == 1)
                        {
                            interpolate(&jn(0), qp, dv_j_node, X_phi);
							interpolate(&jnn(0), qp, d2v_j_node, X_phi);
                            C_u_um = sdh_um * jn(axis);
                            C_u_up = sdh_up * jn(axis);
                            if (d_use_higher_order_jump)
                            {
							    C_u_up +=  0.5 * sdh_up * sdh_up * jnn(axis);
								C_u_um +=  0.5 * sdh_um * sdh_um * jnn(axis);
							}
                        }
                        else
                        {
                            interpolate(&jn(0), qp, dw_j_node, X_phi);
							interpolate(&jn(0), qp, d2w_j_node, X_phi);
                            C_u_um = sdh_um * jn(axis);
                            C_u_up = sdh_up * jn(axis);
                            if (d_use_higher_order_jump)
                            {
							    C_u_um +=  0.5 * sdh_um * sdh_um * jnn(axis);
								C_u_up +=  0.5 * sdh_up * sdh_up * jnn(axis);
							}
                        }

                        (*f_data)(i_side_um) += - (n(axis) > 0.0 ? -1.0 : 1.0) * (C_u_up / (dx[axis] * dx[axis]));

                        (*f_data)(i_side_up) += - (n(axis) > 0.0 ? 1.0 : -1.0) * (C_u_um / (dx[axis] * dx[axis]));
                    }
                }
            }

            if (integrate_tangential_force && !intersectionSide2_ref_coords_u.empty())
            {
                TBOX_ASSERT(intersectionSide2_indices_um.size() == intersectionSide2_indices_up.size());
                X_fe_base->reinit(elem, &intersectionSide2_ref_coords_u);
                const size_t n_qp = intersectionSide2_indices_up.size();

                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    const SideIndex<NDIM>& i_side_up = intersectionSide2_indices_up[qp];

                    const SideIndex<NDIM>& i_side_um = intersectionSide2_indices_um[qp];

                    const unsigned int dd = i_side_up.getAxis();
                    const unsigned int ddd = i_side_um.getAxis();
                    TBOX_ASSERT(dd == ddd);

                    libMesh::Point xx = intersectionSide2_coords_u[qp];

                    interpolate(&x(0), qp, X_node, X_phi);
                    interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
                    if (dim == 1)
                        tau2 = VectorValue<double>(0.0, 0.0, 1.0);
                    else
                        interpolate(&tau2(0), qp, X_node, X_dphi_deta);
                    n = tau1.cross(tau2);
                    n = n.unit();

                    if (side_boxes[dd].contains(i_side_up) && side_boxes[dd].contains(i_side_um))
                    {
                        int axis;

                        if (dd == 0)
                        {
                            axis = 1;
                        }
                        else if (dd == 1)
                        {
                            axis = 2;
                        }
                        else
                        {
                            axis = 0;
                        }

                        // Impose the jump conditions.
                        // *********************************************************************************************************//
                        TBOX_ASSERT(i_side_up(axis) - i_side_um(axis) == 1);

                        double x_mid_side_up =
                            x_lower[axis] + static_cast<double>(i_side_up(axis) - patch_lower[axis] + 0.5) * dx[axis];

                        double x_mid_side_um =
                            x_lower[axis] + static_cast<double>(i_side_um(axis) - patch_lower[axis] + 0.5) * dx[axis];

                        TBOX_ASSERT(xx(axis) <= x_mid_side_up);
                        TBOX_ASSERT(xx(axis) > x_mid_side_um);

                        double sdh_up = xx(axis) - x_mid_side_up; // Signed Distance h
                        double sdh_um = xx(axis) - x_mid_side_um;

                        double C_u_um, C_u_up;

                        if (dd == 0)
                        {
                            interpolate(&jn(0), qp, du_j_node, X_phi);
							interpolate(&jnn(0), qp, du_j_node, X_phi);
                            C_u_um = sdh_um * jn(axis); 
                            C_u_up = sdh_up * jn(axis); 
                            if (d_use_higher_order_jump)
                            {
							    C_u_um +=  0.5 * sdh_um * sdh_um * jnn(axis);
								C_u_up +=  0.5 * sdh_up * sdh_up * jnn(axis);
							}
                        }
                        else if (dd == 1)
                        {
                            interpolate(&jn(0), qp, dv_j_node, X_phi);
							interpolate(&jnn(0), qp, d2v_j_node, X_phi);
                            C_u_um = sdh_um * jn(axis) ;
                            C_u_up = sdh_up * jn(axis) ;
                            if (d_use_higher_order_jump)
                            {
							    C_u_um +=  0.5 * sdh_um * sdh_um * jnn(axis);
								C_u_up +=  0.5 * sdh_up * sdh_up * jnn(axis);
							}
                        }
                        else
                        {
                            interpolate(&jn(0), qp, dw_j_node, X_phi);
							interpolate(&jnn(0), qp, d2w_j_node, X_phi);
                            C_u_um = sdh_um * jn(axis);
                            C_u_up = sdh_up * jn(axis);
                            if (d_use_higher_order_jump)
                            {
							    C_u_um +=  0.5 * sdh_um * sdh_um * jnn(axis);
								C_u_up +=  0.5 * sdh_up * sdh_up * jnn(axis);
							}
                        }

                        (*f_data)(i_side_um) += - (n(axis) > 0.0 ? -1.0 : 1.0) * (C_u_up / (dx[axis] * dx[axis]));

                        (*f_data)(i_side_up) += - (n(axis) > 0.0 ? 1.0 : -1.0) * (C_u_um / (dx[axis] * dx[axis]));

                        // **********************************************************************************************************//

                        // **********************************************************************************************************//
                    }
                }
            }

#endif
        }
    }

    VecRestoreArray(X_local_vec, &X_local_soln);
    VecGhostRestoreLocalForm(X_global_vec, &X_local_vec);


    VecRestoreArray(P_j_local_vec, &P_j_local_soln);
    VecGhostRestoreLocalForm(P_j_global_vec, &P_j_local_vec);

    VecRestoreArray(dP_j_local_vec, &dP_j_local_soln);
    VecGhostRestoreLocalForm(dP_j_global_vec, &dP_j_local_vec);

    VecRestoreArray(du_j_local_vec, &du_j_local_soln);
    VecGhostRestoreLocalForm(du_j_global_vec, &du_j_local_vec);

    VecRestoreArray(dv_j_local_vec, &dv_j_local_soln);
    VecGhostRestoreLocalForm(dv_j_global_vec, &dv_j_local_vec);

    VecRestoreArray(dw_j_local_vec, &dw_j_local_soln);
    VecGhostRestoreLocalForm(dw_j_global_vec, &dw_j_local_vec);


    VecRestoreArray(d2u_j_local_vec, &d2u_j_local_soln);
    VecGhostRestoreLocalForm(d2u_j_global_vec, &d2u_j_local_vec);

    VecRestoreArray(d2v_j_local_vec, &d2v_j_local_soln);
    VecGhostRestoreLocalForm(d2v_j_global_vec, &d2v_j_local_vec);

    VecRestoreArray(d2w_j_local_vec, &d2w_j_local_soln);
    VecGhostRestoreLocalForm(d2w_j_global_vec, &d2w_j_local_vec);
    return;



    return;
} // imposeJumpConditionsWeak




void
IBFEMethod::imposeJumpConditionsPointWise(const int f_data_idx,
                                 const int mask_current_idx,
                                 PetscVector<double>& /*F_ghost_vec*/,
                                 PetscVector<double>& X_ghost_vec,
                                 PetscVector<double>& P_j_ghost_vec,
                                 PetscVector<double>& dP_j_ghost_vec,
                                 PetscVector<double>& du_j_ghost_vec,
                                 PetscVector<double>& dv_j_ghost_vec,
                                 PetscVector<double>& dw_j_ghost_vec,  
                                 PetscVector<double>& d2u_j_ghost_vec,
                                 PetscVector<double>& d2v_j_ghost_vec,
								 PetscVector<double>& d2w_j_ghost_vec,
                                 const double data_time,
                                 const unsigned int part)

{
    if (!d_split_normal_force && !d_split_tangential_force) return;

    // Check to see if we need to integrate the normal surface force.
    const bool integrate_normal_force = d_split_normal_force && d_use_jump_conditions;
    const bool integrate_tangential_force = d_split_tangential_force && d_use_jump_conditions;
    // if (!integrate_normal_force) return;

    // Extract the mesh.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    TBOX_ASSERT(dim == NDIM - 1);

    std::vector<const std::vector<double>*> force_var_data;
    std::vector<const std::vector<VectorValue<double> >*> force_grad_var_data;


    System& F_system = equation_systems->get_system(FORCE_SYSTEM_NAME);
    const DofMap& F_dof_map = F_system.get_dof_map();
    FEDataManager::SystemDofMapCache& F_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(FORCE_SYSTEM_NAME);
    FEType F_fe_type = F_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(F_dof_map.variable_type(d) == F_fe_type);
    }
    std::vector<std::vector<unsigned int> > F_dof_indices(NDIM);

    std::vector<int> vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
    FEDataInterpolation fe(dim, d_fe_data_managers[part]);
    AutoPtr<QBase> qrule = QBase::build(d_default_quad_type[part], dim, d_default_quad_order[part]);
    fe.attachQuadratureRule(qrule.get());
    fe.evalQuadraturePoints();
    fe.evalQuadratureWeights();
    fe.registerSystem(F_system, vars, vars); // compute phi and dphi for the force system
    std::vector<size_t> force_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        force_fcn_system_idxs, d_lag_force_fcn_data[part].system_data, equation_systems);
    fe.init(/*use_IB_ghosted_vecs*/ true);
    

     const std::vector<libMesh::Point>& q_point = fe.getQuadraturePoints();
    // Extract the FE systems and DOF maps, and setup the FE object.
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& X_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();

    FEType X_fe_type = X_dof_map.variable_type(0);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == X_fe_type);
    }
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);

   
    System& X0_system = equation_systems->get_system(COORDS0_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& X0_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(COORDS0_SYSTEM_NAME);
    const DofMap& X0_dof_map = X0_system.get_dof_map();
    FEType X0_fe_type = X0_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X0_dof_map.variable_type(d) == X0_fe_type);
    }
    TBOX_ASSERT(X0_fe_type == X_fe_type);
    std::vector<std::vector<unsigned int> > X0_dof_indices(NDIM);

    // Setup global and elemental right-hand-side vectors.

    //~ if (integrate_normal_force)
    //~ {
    System& P_j_system = equation_systems->get_system(P_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& P_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(P_J_SYSTEM_NAME);
    DofMap& P_j_dof_map = P_j_system.get_dof_map();
    TBOX_ASSERT(P_j_dof_map.variable_type(0) == X_fe_type);
    std::vector<unsigned int> P_j_dof_indices;
    //~ }

    System& du_j_system = equation_systems->get_system(DU_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& du_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DU_J_SYSTEM_NAME);
    const DofMap& du_j_dof_map = du_j_system.get_dof_map();

    FEType du_j_fe_type = du_j_dof_map.variable_type(0);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(du_j_dof_map.variable_type(d) == du_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > du_j_dof_indices(NDIM);

    System& dv_j_system = equation_systems->get_system(DV_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& dv_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DV_J_SYSTEM_NAME);
    const DofMap& dv_j_dof_map = dv_j_system.get_dof_map();
    FEType dv_j_fe_type = dv_j_dof_map.variable_type(0);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(dv_j_dof_map.variable_type(d) == dv_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > dv_j_dof_indices(NDIM);


    System& dw_j_system = equation_systems->get_system(DW_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& dw_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DW_J_SYSTEM_NAME);
    const DofMap& dw_j_dof_map = dw_j_system.get_dof_map();
    FEType dw_j_fe_type = dw_j_dof_map.variable_type(0);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(dw_j_dof_map.variable_type(d) == dw_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > dw_j_dof_indices(NDIM);


    System& d2u_j_system = equation_systems->get_system(D2U_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& d2u_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(D2U_J_SYSTEM_NAME);
    const DofMap& d2u_j_dof_map = d2u_j_system.get_dof_map();

    FEType d2u_j_fe_type = d2u_j_dof_map.variable_type(0);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d2u_j_dof_map.variable_type(d) == d2u_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > d2u_j_dof_indices(NDIM);

    System& d2v_j_system = equation_systems->get_system(D2V_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& d2v_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(D2V_J_SYSTEM_NAME);
    const DofMap& d2v_j_dof_map = d2v_j_system.get_dof_map();
    FEType d2v_j_fe_type = d2v_j_dof_map.variable_type(0);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d2v_j_dof_map.variable_type(d) == d2v_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > d2v_j_dof_indices(NDIM);

    System& d2w_j_system = equation_systems->get_system(D2W_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& d2w_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(D2W_J_SYSTEM_NAME);
    const DofMap& d2w_j_dof_map = d2w_j_system.get_dof_map();
    FEType d2w_j_fe_type = d2w_j_dof_map.variable_type(0);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d2w_j_dof_map.variable_type(d) == d2w_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > d2w_j_dof_indices(NDIM);
    
    
    
    System& dP_j_system = equation_systems->get_system(DP_J_SYSTEM_NAME);
    FEDataManager::SystemDofMapCache& dP_j_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(DP_J_SYSTEM_NAME);
    const DofMap& dP_j_dof_map = dP_j_system.get_dof_map();
    FEType dP_j_fe_type = dP_j_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(dP_j_dof_map.variable_type(d) == dP_j_fe_type);
    }
    std::vector<std::vector<unsigned int> > dP_j_dof_indices(NDIM);


    AutoPtr<FEBase> X_fe_base(FEBase::build(dim, X_fe_type));
    const std::vector<std::vector<double> >& X_phi = X_fe_base->get_phi();
    const std::vector<std::vector<double> >& X_dphi_dxi = X_fe_base->get_dphidxi();
    const std::vector<std::vector<double> >& X_dphi_deta = X_fe_base->get_dphideta();

    PetscVector<double>* X_petsc_vec = static_cast<PetscVector<double>*>(&X_ghost_vec);
    Vec X_global_vec = X_petsc_vec->vec();
    Vec X_local_vec;
    VecGhostGetLocalForm(X_global_vec, &X_local_vec);
    double* X_local_soln;
    VecGetArray(X_local_vec, &X_local_soln);


    PetscVector<double>* X0_petsc_vec = static_cast<PetscVector<double>*>(d_X0_vecs[part]);
    Vec X0_global_vec = X0_petsc_vec->vec();
    Vec X0_local_vec;
    VecGhostGetLocalForm(X0_global_vec, &X0_local_vec);
    double* X0_local_soln;
    VecGetArray(X0_local_vec, &X0_local_soln);

    //~ if (integrate_normal_force)
    //~ {
    PetscVector<double>* P_j_petsc_vec = static_cast<PetscVector<double>*>(&P_j_ghost_vec);
    Vec P_j_global_vec = P_j_petsc_vec->vec();
    Vec P_j_local_vec;
    VecGhostGetLocalForm(P_j_global_vec, &P_j_local_vec);
    double* P_j_local_soln;
    VecGetArray(P_j_local_vec, &P_j_local_soln);
    
    PetscVector<double>* dP_j_petsc_vec = static_cast<PetscVector<double>*>(&dP_j_ghost_vec);
    Vec dP_j_global_vec = dP_j_petsc_vec->vec();
    Vec dP_j_local_vec;
    VecGhostGetLocalForm(dP_j_global_vec, &dP_j_local_vec);
    double* dP_j_local_soln;
    VecGetArray(dP_j_local_vec, &dP_j_local_soln);
    //~ }
    PetscVector<double>* du_j_petsc_vec = static_cast<PetscVector<double>*>(&du_j_ghost_vec);
    Vec du_j_global_vec = du_j_petsc_vec->vec();
    Vec du_j_local_vec;
    VecGhostGetLocalForm(du_j_global_vec, &du_j_local_vec);
    double* du_j_local_soln;
    VecGetArray(du_j_local_vec, &du_j_local_soln);

    PetscVector<double>* dv_j_petsc_vec = static_cast<PetscVector<double>*>(&dv_j_ghost_vec);
    Vec dv_j_global_vec = dv_j_petsc_vec->vec();
    Vec dv_j_local_vec;
    VecGhostGetLocalForm(dv_j_global_vec, &dv_j_local_vec);
    double* dv_j_local_soln;
    VecGetArray(dv_j_local_vec, &dv_j_local_soln);

    PetscVector<double>* dw_j_petsc_vec = static_cast<PetscVector<double>*>(&dw_j_ghost_vec);
    Vec dw_j_global_vec = dw_j_petsc_vec->vec();
    Vec dw_j_local_vec;
    VecGhostGetLocalForm(dw_j_global_vec, &dw_j_local_vec);
    double* dw_j_local_soln;
    VecGetArray(dw_j_local_vec, &dw_j_local_soln);

    PetscVector<double>* d2u_j_petsc_vec = static_cast<PetscVector<double>*>(&d2u_j_ghost_vec);
    Vec d2u_j_global_vec = d2u_j_petsc_vec->vec();
    Vec d2u_j_local_vec;
    VecGhostGetLocalForm(d2u_j_global_vec, &d2u_j_local_vec);
    double* d2u_j_local_soln;
    VecGetArray(d2u_j_local_vec, &d2u_j_local_soln);

    PetscVector<double>* d2v_j_petsc_vec = static_cast<PetscVector<double>*>(&d2v_j_ghost_vec);
    Vec d2v_j_global_vec = d2v_j_petsc_vec->vec();
    Vec d2v_j_local_vec;
    VecGhostGetLocalForm(d2v_j_global_vec, &d2v_j_local_vec);
    double* d2v_j_local_soln;
    VecGetArray(d2v_j_local_vec, &d2v_j_local_soln);

    PetscVector<double>* d2w_j_petsc_vec = static_cast<PetscVector<double>*>(&d2w_j_ghost_vec);
    Vec d2w_j_global_vec = d2w_j_petsc_vec->vec();
    Vec d2w_j_local_vec;
    VecGhostGetLocalForm(d2w_j_global_vec, &d2w_j_local_vec);
    double* d2w_j_local_soln;
    VecGetArray(d2w_j_local_vec, &d2w_j_local_soln);
    
    
    const std::vector<std::vector<double> >& phi = fe.getPhi(F_fe_type);

    // Loop over the patches to impose jump conditions on the Eulerian grid that
    // are determined from the interior and transmission elastic force
    // densities.
    const std::vector<std::vector<Elem*> >& active_patch_element_map =
        d_fe_data_managers[part]->getActivePatchElementMap();
    const int level_num = d_fe_data_managers[part]->getLevelNumber();
    std::vector<libMesh::Point> s_node_cache, x_node_cache;
    TensorValue<double> FF;
    VectorValue<double> F, F_s, n, s, N, tau1, tau2, Tau2, Tau1, x, X, jn, jnn, tau1unit, tau2unit;
    double dA_da;
    IBTK::Point x_min, x_max;
    //~ if (integrate_normal_force)
    boost::multi_array<double, 1> P_j_node;
    boost::multi_array<double, 2> du_j_node, dv_j_node, dw_j_node, dP_j_node;

    boost::multi_array<double, 2> d2u_j_node, d2v_j_node, d2w_j_node;

    std::vector<std::vector<unsigned int> > side_dof_indices(NDIM);

    //~ if (integrate_normal_force)
    //~ {
    std::vector<libMesh::Point> intersection_ref_coords_p;
    std::vector<libMesh::Point> intersectionSide_ref_coords_u;
    std::vector<libMesh::Point> intersectionSide2_ref_coords_u;
    std::vector<libMesh::Point> intersection_ref_coords_u;
    std::vector<libMesh::Point> intersection_coords_u;
    std::vector<libMesh::Point> intersection_coords_p;
    std::vector<libMesh::Point> intersectionSide_coords_u;
    std::vector<libMesh::Point> intersectionSide2_coords_u;

    std::vector<SideIndex<NDIM> > intersection_indices_p;
    //~ }
    //~ if (integrate_tangential_force)
    //~ {
    std::vector<SideIndex<NDIM> > intersection_indices_um;

    std::vector<SideIndex<NDIM> > intersection_indices_up;

    std::vector<SideIndex<NDIM> > intersectionSide_indices_up;

    std::vector<SideIndex<NDIM> > intersectionSide_indices_um;

    //~ #if(NDIM == 3)
    std::vector<SideIndex<NDIM> > intersectionSide2_indices_up;

    std::vector<SideIndex<NDIM> > intersectionSide2_indices_um;
    //~ #endif
    //~ }

    std::vector<std::pair<double, libMesh::Point> > intersections;
    std::vector<std::pair<double, libMesh::Point> > intersectionsSide;
    std::vector<std::pair<double, libMesh::Point> > intersectionsSide2;

    boost::multi_array<double, 2> X_node, X0_node;

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
        Pointer<SideData<NDIM, double> > mask_current = patch->getPatchData(mask_current_idx);

        const Box<NDIM>& patch_box = patch->getBox();
        Box<NDIM> side_boxes[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            side_boxes[d] = SideGeometry<NDIM>::toSideBox(patch_box, d);
        }
        const CellIndex<NDIM>& patch_lower = patch_box.lower();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const x_lower = patch_geom->getXLower();
        const double* const x_upper = patch_geom->getXUpper();
        const double* const dx = patch_geom->getDx();

        SideData<NDIM, int> num_intersections_p(patch_box, 1, IntVector<NDIM>(0));
        num_intersections_p.fillAll(0);
    
        SideData<NDIM, int> num_intersections_um(patch_box, 1, IntVector<NDIM>(0));
        num_intersections_um.fillAll(0);

        SideData<NDIM, int> num_intersections_up(patch_box, 1, IntVector<NDIM>(0));
        num_intersections_up.fillAll(0);

        SideData<NDIM, int> num_intersectionsSide_um(patch_box, 1, IntVector<NDIM>(1));
        num_intersectionsSide_um.fillAll(0);

        SideData<NDIM, int> num_intersectionsSide_up(patch_box, 1, IntVector<NDIM>(1));
        num_intersectionsSide_up.fillAll(0);

        SideData<NDIM, int> num_intersectionsSide2_um(patch_box, 1, IntVector<NDIM>(1));
        num_intersectionsSide2_um.fillAll(0);

        SideData<NDIM, int> num_intersectionsSide2_up(patch_box, 1, IntVector<NDIM>(1));
        num_intersectionsSide2_up.fillAll(0);
     
 
        // Loop over the elements.
        for (size_t e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const unsigned int n_node = elem->n_nodes();

            for (int d = 0; d < NDIM; ++d)
            {
                F_dof_map_cache.dof_indices(elem, F_dof_indices[d], d);
                X_dof_map_cache.dof_indices(elem, X_dof_indices[d], d);
                X0_dof_map_cache.dof_indices(elem, X0_dof_indices[d], d);
                du_j_dof_map_cache.dof_indices(elem, du_j_dof_indices[d], d);
                dv_j_dof_map_cache.dof_indices(elem, dv_j_dof_indices[d], d);
				dw_j_dof_map_cache.dof_indices(elem, dw_j_dof_indices[d], d);
                d2u_j_dof_map_cache.dof_indices(elem, d2u_j_dof_indices[d], d);
                d2v_j_dof_map_cache.dof_indices(elem, d2v_j_dof_indices[d], d);
				d2w_j_dof_map_cache.dof_indices(elem, d2w_j_dof_indices[d], d);
                dP_j_dof_map_cache.dof_indices(elem, dP_j_dof_indices[d], d);
            }
            
            
           fe.reinit(elem);
           fe.collectDataForInterpolation(elem);
           fe.interpolate(elem);
        // X_fe_base->reinit(elem);
//~ 
            get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, X_dof_indices);
            get_values_for_interpolation(X0_node, *X0_petsc_vec, X0_local_soln, X0_dof_indices);

            if (integrate_normal_force)
            {
				P_j_dof_map_cache.dof_indices(elem, P_j_dof_indices);
                get_values_for_interpolation(P_j_node, *P_j_petsc_vec, P_j_local_soln, P_j_dof_indices);
            }
            get_values_for_interpolation(du_j_node, *du_j_petsc_vec, du_j_local_soln, du_j_dof_indices);
            get_values_for_interpolation(dv_j_node, *dv_j_petsc_vec, dv_j_local_soln, dv_j_dof_indices);
            get_values_for_interpolation(dw_j_node, *dw_j_petsc_vec, dw_j_local_soln, dw_j_dof_indices);

            get_values_for_interpolation(d2u_j_node, *d2u_j_petsc_vec, d2u_j_local_soln, d2u_j_dof_indices);
            get_values_for_interpolation(d2v_j_node, *d2v_j_petsc_vec, d2v_j_local_soln, d2v_j_dof_indices);
            get_values_for_interpolation(d2w_j_node, *d2w_j_petsc_vec, d2w_j_local_soln, d2w_j_dof_indices);
            get_values_for_interpolation(dP_j_node, *dP_j_petsc_vec, dP_j_local_soln, dP_j_dof_indices);

            // Cache the nodal and physical coordinates of the side element,
            // determine the bounding box of the current configuration of
            // the side element, and set the nodal coordinates to correspond
            // to the physical coordinates.
            s_node_cache.resize(n_node);
            x_node_cache.resize(n_node);
            x_min = IBTK::Point::Constant(std::numeric_limits<double>::max());
            x_max = IBTK::Point::Constant(-std::numeric_limits<double>::max());
            for (unsigned int k = 0; k < n_node; ++k)
            {
                s_node_cache[k] = elem->point(k);
                libMesh::Point& x = x_node_cache[k];
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    x(d) = X_ghost_vec(X_dof_indices[d][k]);
                    x_min[d] = std::min(x_min[d], x(d));
                    x_max[d] = std::max(x_max[d], x(d));
                }
                elem->point(k) = x;
            }
            Box<NDIM> box(IndexUtilities::getCellIndex(&x_min[0], grid_geom, ratio),
                          IndexUtilities::getCellIndex(&x_max[0], grid_geom, ratio));
            box.grow(IntVector<NDIM>(1));
            box = box * patch_box;


            // Loop over coordinate directions and look for intersections
            // with the background fluid grid.
            if (integrate_normal_force)
            {
                intersection_ref_coords_p.clear();
                intersection_coords_p.clear();
                intersection_indices_p.clear();
            }

            if (integrate_tangential_force)
            {
                intersection_ref_coords_u.clear();
                intersection_coords_u.clear();
                intersectionSide_ref_coords_u.clear();
                intersectionSide2_ref_coords_u.clear();
                intersectionSide_coords_u.clear();
                intersectionSide2_coords_u.clear();

                intersection_indices_um.clear();
                intersection_indices_up.clear();
                intersectionSide_indices_um.clear();
                intersectionSide_indices_up.clear();

                intersectionSide2_indices_um.clear();
                intersectionSide2_indices_up.clear();
            }
            VectorValue<double> DX;
            
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
				 DX(axis) = dx[axis];
			}

            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                // Setup a unit vector pointing in the coordinate direction of interest.
                VectorValue<double> q;
                q(axis) = 1.0;
                

                // Loop over the relevant range of indices.
                Box<NDIM> axis_box = box;
                axis_box.lower(axis) = 0;
                axis_box.upper(axis) = 0;
                for (BoxIterator<NDIM> b(axis_box); b; b++)
                {
                    const Index<NDIM>& i_c = b();
                    libMesh::Point r;
                    libMesh::Point rs;
                    libMesh::Point rss;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        r(d) = (d == axis ? 0.0 : x_lower[d] +
                                                      dx[d] * (static_cast<double>(i_c(d) - patch_lower[d]) +
                                                               0.5)); // In 2D this corresponds to applying jumps for
                                                                      // [p],[ux] and [vy]
#if (NDIM == 2)
                        rs(d) = (d == axis ?
                                     0.0 :
                                     x_lower[d] + dx[d] * (static_cast<double>(i_c(d) - patch_lower[d]))); // In 2D this
                        // corresponds to
                        // applying jumps
                        // for [uy] and
                        // [vx]
#endif
                    }

#if (NDIM == 3)
                    if (axis == 0)
                    {
                        rs(0) = 0.0;
                        rs(1) = x_lower[1] + dx[1] * (static_cast<double>(i_c(1) - patch_lower[1]));
                        rs(2) = x_lower[2] + dx[2] * (static_cast<double>(i_c(2) - patch_lower[2]) + 0.5);
                        rss(0) = 0.0;
                        rss(1) = x_lower[1] + dx[1] * (static_cast<double>(i_c(1) - patch_lower[1]) + 0.5);
                        rss(2) = x_lower[2] + dx[2] * (static_cast<double>(i_c(2) - patch_lower[2]));
                    }
                    else if (axis == 1)
                    {
                        rs(0) = x_lower[0] + dx[0] * (static_cast<double>(i_c(0) - patch_lower[0]) + 0.5);
                        rs(1) = 0.0;
                        rs(2) = x_lower[2] + dx[2] * (static_cast<double>(i_c(2) - patch_lower[2]));
                        rss(0) = x_lower[0] + dx[0] * (static_cast<double>(i_c(0) - patch_lower[0]));
                        rss(1) = 0.0;
                        rss(2) = x_lower[2] + dx[2] * (static_cast<double>(i_c(2) - patch_lower[2]) + 0.5);
                    }
                    else if (axis == 2)
                    {
                        rs(0) = x_lower[0] + dx[0] * (static_cast<double>(i_c(0) - patch_lower[0]));
                        rs(1) = x_lower[1] + dx[1] * (static_cast<double>(i_c(1) - patch_lower[1]) + 0.5);
                        rs(2) = 0.0;
                        rss(0) = x_lower[0] + dx[0] * (static_cast<double>(i_c(0) - patch_lower[0]) + 0.5);
                        rss(1) = x_lower[1] + dx[1] * (static_cast<double>(i_c(1) - patch_lower[1]));
                        rss(2) = 0.0;
                    }

#endif

#if (NDIM == 2)
                    intersect_line_with_edge(intersections, static_cast<Edge*>(elem), DX, r, q);
                    intersect_line_with_edge(intersectionsSide, static_cast<Edge*>(elem), DX, rs, q);
#endif
#if (NDIM == 3)
                    intersect_line_with_face(intersections, static_cast<Face*>(elem), r, q);
                    intersect_line_with_face(intersectionsSide, static_cast<Face*>(elem), rs, q);
                    intersect_line_with_face(intersectionsSide2, static_cast<Face*>(elem), rss, q);
#endif
                    if (integrate_normal_force)
                    {
                        for (unsigned int k = 0; k < intersections.size(); ++k)
                        {
                            libMesh::Point x = r + intersections[k].first * q;
                            SideIndex<NDIM> i_s(i_c, axis, 0);
                            i_s(axis) = static_cast<int>(std::floor((x(axis) - x_lower[axis]) / dx[axis] + 0.5)) +
                                        patch_lower[axis];

                            if (side_boxes[axis].contains(i_s))
                            {
                                intersection_ref_coords_p.push_back(intersections[k].second);
                                intersection_coords_p.push_back(x);
                                intersection_indices_p.push_back(i_s);
                                num_intersections_p(i_s) += 1;
                            }
                        }
                    }

  
                    if (integrate_tangential_force)
                    {
                        for (unsigned int k = 0; k < intersections.size(); ++k)
                        {
                            libMesh::Point xs = r + intersections[k].first * q;
                            SideIndex<NDIM> i_ss(i_c, axis, 0);
                            Index<NDIM> i_c_neighbor = i_c;
                            i_c_neighbor(axis) += 1;
                            
                            SideIndex<NDIM> i_sss(i_c_neighbor, axis, 0);
                            i_sss(axis) = static_cast<int>(std::floor((xs(axis) - x_lower[axis]) / dx[axis] + 1.0)) +
                                          patch_lower[axis];
                            i_ss(axis) =
                                static_cast<int>(std::floor((xs(axis) - x_lower[axis]) / dx[axis])) + patch_lower[axis];

                            const double x_cell_bdry_um =
                                x_lower[axis] + static_cast<double>(i_ss(axis) - patch_lower[axis]) * dx[axis];

                            const double x_cell_bdry_up =
                                x_lower[axis] + static_cast<double>(i_sss(axis) - patch_lower[axis]) * dx[axis];

                            const double sdh_um = ((xs(axis) - x_cell_bdry_um)); // Signed Distance h

                            const double sdh_up = ((xs(axis) - x_cell_bdry_up)); // Signed Distance h

                            TBOX_ASSERT((sdh_um) < dx[axis] && sdh_um > 0);

                            if (side_boxes[axis].contains(i_ss) && side_boxes[axis].contains(i_sss))
                            {
                                intersection_ref_coords_u.push_back(intersections[k].second);
                                intersection_coords_u.push_back(xs);
                                intersection_indices_um.push_back(i_ss);
                                num_intersections_um(i_ss) += 1;
                               
                                intersection_indices_up.push_back(i_sss);
                                num_intersections_up(i_sss) += 1;
                                
                            }
                        }


//// This will cache the point to the right(or top) of the intersection along
// the side cell

#if (NDIM == 2)

                        for (unsigned int k = 0; k < intersectionsSide.size(); ++k)
                        {
                            libMesh::Point xu = rs + intersectionsSide[k].first * q;
                            int dd = (axis == 0 ? 1 : 0);
                            if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) >= 0.5 * dx[axis])
                            {
                                SideIndex<NDIM> i_s_um(i_c, dd, 0);
                                Index<NDIM> i_c_neighbor = i_c;
                                i_c_neighbor(axis) += 1;
                             
                                SideIndex<NDIM> i_s_up(i_c_neighbor, dd, 0);
                          
                                i_s_up(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis] + 0.5)) +
                                             patch_lower[axis];             
                                i_s_um(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis])) +
                                              patch_lower[axis];
                                              
                                if (side_boxes[axis].contains(i_s_up) && side_boxes[axis].contains(i_s_um))
                                {
                                    intersectionSide_ref_coords_u.push_back(intersectionsSide[k].second);
                                    intersectionSide_coords_u.push_back(xu);
                                    intersectionSide_indices_up.push_back(i_s_up);
                                    num_intersectionsSide_up(i_s_up) += 1;
                                    intersectionSide_indices_um.push_back(i_s_um);
                                    num_intersectionsSide_um(i_s_um) += 1;
                                 
                                }
                            }
                            else if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) < 0.5 * dx[axis] &&
                                     fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) >= 0.0)
                            {
								
								SideIndex<NDIM> i_s_up(i_c, dd, 0);
								Index<NDIM> i_c_neighbor = i_c;
								i_c_neighbor(axis) -= 1;
                                SideIndex<NDIM> i_s_um(i_c_neighbor, dd, 0);
                                i_s_up(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis])) +
                                             patch_lower[axis];
                                i_s_um(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis] - 0.5)) +
											 patch_lower[axis];

                                if (side_boxes[axis].contains(i_s_up) &&  side_boxes[axis].contains(i_s_um))
                                {
                                    intersectionSide_indices_up.push_back(i_s_up);
                                    num_intersectionsSide_up(i_s_up) += 1;
                                    intersectionSide_indices_um.push_back(i_s_um);
                                    num_intersectionsSide_up(i_s_um) += 1;
                                    
                                    intersectionSide_ref_coords_u.push_back(intersectionsSide[k].second);
                                    intersectionSide_coords_u.push_back(xu);
                                }
                            }
                            else
                            {
								pout<< "improper side index found!"<<"\n\n";
						    }
                        }

       

#endif
                        int SideDim[2];
                        SideDim[0] = SideDim[1] = 0;
                        if (axis == 0)
                        {
                            SideDim[0] = 1;
                            SideDim[1] = 2;
                        }
                        else if (axis == 1)
                        {
                            SideDim[0] = 2;
                            SideDim[1] = 0;
                        }
                        else if (axis == 2)
                        {
                            SideDim[0] = 0;
                            SideDim[1] = 1;
                        }

#if (NDIM == 3)

                        for (unsigned int k = 0; k < intersectionsSide.size(); ++k)
                        {
                            libMesh::Point xu = rs + intersectionsSide[k].first * q;
                            if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) >= 0.5 * dx[axis])
                            {
                                SideIndex<NDIM> i_s_um(i_c, SideDim[0], 0);
                                Index<NDIM> i_c_neighbor = i_c;
                                i_c_neighbor(axis) += 1;

                                SideIndex<NDIM> i_s_up(i_c_neighbor, SideDim[0], 0);

                                i_s_up(axis) =
                                    static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis] + 0.5)) +
                                    patch_lower[axis];
                                i_s_um(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis])) +
                                               patch_lower[axis];

                                if (side_boxes[axis].contains(i_s_up) && side_boxes[axis].contains(i_s_um))
                                {
                                    TBOX_ASSERT(i_s_up(axis) - i_s_um(axis) == 1);
                                    intersectionSide_ref_coords_u.push_back(intersectionsSide[k].second);
                                    intersectionSide_coords_u.push_back(xu);
                                    intersectionSide_indices_up.push_back(i_s_up);
                                    num_intersectionsSide_up(i_s_up) += 1;
                                    intersectionSide_indices_um.push_back(i_s_um);
                                    num_intersectionsSide_um(i_s_um) += 1;
                                }
                            }
                            else if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) < 0.5 * dx[axis] &&
                                     fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) >= 0.0)
                            {
                                SideIndex<NDIM> i_s_up(i_c, SideDim[0], 0);
                                Index<NDIM> i_c_neighbor = i_c;
                                i_c_neighbor(axis) -= 1;
                                SideIndex<NDIM> i_s_um(i_c_neighbor, SideDim[0], 0);

                                i_s_up(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis])) +
                                               patch_lower[axis];
                                i_s_um(axis) =
                                    static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis] - 0.5)) +
                                    patch_lower[axis];

                                if (side_boxes[axis].contains(i_s_up) && side_boxes[axis].contains(i_s_um))
                                {
                                    intersectionSide_indices_up.push_back(i_s_up);
                                    num_intersectionsSide_up(i_s_up) += 1;
                                    intersectionSide_indices_um.push_back(i_s_um);
                                    num_intersectionsSide_up(i_s_um) += 1;
                                    intersectionSide_coords_u.push_back(xu);
                                    intersectionSide_ref_coords_u.push_back(intersectionsSide[k].second);
                                }
                            }
                            else
                            {
                                pout << "improper side index found!"
                                     << "\n\n";
                            }
                        }

                        for (unsigned int k = 0; k < intersectionsSide2.size(); ++k)
                        {
                            libMesh::Point xu = rss + intersectionsSide2[k].first * q;
                            if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) >= 0.5 * dx[axis])
                            {
                                SideIndex<NDIM> i_s_um(i_c, SideDim[1], 0);
                                Index<NDIM> i_c_neighbor = i_c;
                                i_c_neighbor(axis) += 1;

                                SideIndex<NDIM> i_s_up(i_c_neighbor, SideDim[1], 0);

                                i_s_up(axis) =
                                    static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis] + 0.5)) +
                                    patch_lower[axis];
                                i_s_um(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis])) +
                                               patch_lower[axis];

                                if (side_boxes[axis].contains(i_s_up) && side_boxes[axis].contains(i_s_um))
                                {
                                    intersectionSide2_ref_coords_u.push_back(intersectionsSide2[k].second);
                                    intersectionSide2_coords_u.push_back(xu);
                                    intersectionSide2_indices_up.push_back(i_s_up);
                                    num_intersectionsSide2_up(i_s_up) += 1;
                                    intersectionSide2_indices_um.push_back(i_s_um);
                                    num_intersectionsSide2_um(i_s_um) += 1;
                                }
                            }
                            else if (fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) < 0.5 * dx[axis] &&
                                     fmod(fabs(xu(axis) - x_lower[axis]), dx[axis]) >= 0.0)
                            {
                                SideIndex<NDIM> i_s_up(i_c, SideDim[1], 0);
                                Index<NDIM> i_c_neighbor = i_c;
                                i_c_neighbor(axis) -= 1;
                                SideIndex<NDIM> i_s_um(i_c_neighbor, SideDim[1], 0);

                                i_s_up(axis) = static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis])) +
                                               patch_lower[axis];
                                i_s_um(axis) =
                                    static_cast<int>(std::floor((xu(axis) - x_lower[axis]) / dx[axis] - 0.5)) +
                                    patch_lower[axis];

                                if (side_boxes[axis].contains(i_s_up) && side_boxes[axis].contains(i_s_um))
                                {
                                    intersectionSide2_indices_up.push_back(i_s_up);
                                    num_intersectionsSide2_up(i_s_up) += 1;
                                    intersectionSide2_indices_um.push_back(i_s_um);
                                    num_intersectionsSide2_up(i_s_um) += 1;

                                    intersectionSide2_ref_coords_u.push_back(intersectionsSide2[k].second);
                                    intersectionSide2_coords_u.push_back(xu);
                                }
                            }
                            else
                            {
                                pout << "improper side index found!"
                                     << "\n\n";
                            }
                        }

#endif
                    }
                }
            }

            // Restore the element coordinates.
            for (unsigned int k = 0; k < n_node; ++k)
            {
                elem->point(k) = s_node_cache[k];
            }

            // If there are no intersection points, then continue to the next side.



            // Evaluate the jump conditions and apply them to the Eulerian grid.
            if (integrate_normal_force && !intersection_ref_coords_p.empty())
            {
                X_fe_base->attach_quadrature_rule(qrule.get());
                X_fe_base->reinit(elem, &intersection_ref_coords_p);
                const size_t n_qp = intersection_ref_coords_p.size();
                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    libMesh::Point& xx = intersection_coords_p[qp];
                    libMesh::Point& XX = intersection_ref_coords_p[qp];
		
                    const SideIndex<NDIM>& i_s = intersection_indices_p[qp];
                    const unsigned int axis = i_s.getAxis();
#if !defined(NDEBUG)
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            const double x_lower_bound = x_lower[d] +
                                                         (static_cast<double>(i_s(d) - patch_lower[d]) - 0.5) * dx[d] -
                                                         sqrt(std::numeric_limits<double>::epsilon());
                            const double x_upper_bound = x_lower[d] +
                                                         (static_cast<double>(i_s(d) - patch_lower[d]) + 0.5) * dx[d] +
                                                         sqrt(std::numeric_limits<double>::epsilon());
                            TBOX_ASSERT(x_lower_bound <= xx(d) && xx(d) <= x_upper_bound);
                        }
                        else
                        {
                            const double x_intersection =
                                x_lower[d] + (static_cast<double>(i_s(d) - patch_lower[d]) + 0.5) * dx[d];
                            const double x_interp = xx(d);
                            const double rel_diff =
                                std::abs(x_intersection - x_interp) /
                                std::max(1.0, std::max(std::abs(x_intersection), std::abs(x_interp)));
                            TBOX_ASSERT(rel_diff <= sqrt(std::numeric_limits<double>::epsilon()));
                        }
                    }
#endif
                    //~ // Impose the jump conditions.
//~ 
//~ 
//~ 
                    interpolate(&X(0), qp, X0_node, X_phi);
					interpolate(&Tau1(0), qp, X0_node, X_dphi_dxi);
					if (dim == 1)
						Tau2 = VectorValue<double>(0.0, 0.0, 1.0);
					else
						interpolate(&Tau2(0), qp, X0_node, X_dphi_deta);
					N = Tau1.cross(Tau2);

					interpolate(&x(0), qp, X_node, X_phi);
					interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
					tau1unit = tau1.unit();

				if (dim == 1)
					tau2 = VectorValue<double>(0.0, 0.0, 1.0);
				else
					interpolate(&tau2(0), qp, X_node, X_dphi_deta);
	
					n = tau1.cross(tau2);
					tau2unit = tau2.unit();
					dA_da = sqrt(N * N) / sqrt(n * n);
					N = N.unit();
					n = n.unit();

				
					double* nn;

					nn = &n(0);
					
					F.zero();
//~ 
					//~ if (d_lag_force_fcn_data[part].fcn)
					//~ {
						// Compute the value of the body force at the quadrature
						//~ // point and add the corresponding forces to the
						//~ // right-hand-side vector.
						//~ fe.setInterpolatedDataPointers(force_var_data, force_grad_var_data, force_fcn_system_idxs, elem, qp);
						//~ d_lag_force_fcn_data[part].fcn(F_s,
													   //~ FF,
													   //~ xx,
													   //~ XX,
													   //~ elem,
													   //~ nn,
													   //~ force_var_data,
													   //~ force_grad_var_data,
													   //~ data_time,
					for (unsigned int dod = 0; dod < NDIM; ++dod)
					{                                 //~ d_lag_force_fcn_data[part].ctx);
					     F(dod) +=  1000 * (X(dod) - x(dod));
					                             
					}
					
					F *=dA_da;
					

                    
                    
                    interpolate(&jn(0), qp, dP_j_node, X_phi);
                    
                    const double x_cell_bdry =
                                            x_lower[axis] + static_cast<double>(i_s(axis) - patch_lower[axis]) *
                                            dx[axis];
											const double h = x_cell_bdry + (x(axis) > x_cell_bdry ? +0.5 : -0.5) * dx[axis] -
											x(axis);
                       const double C_p = F * n;

                    (*f_data)(i_s) += (n(axis) > 0.0 ? +1.0 : -1.0) * (C_p / dx[axis]);
                }
            }

            // Applying the jump for viscous term (in 2D u_xx, u_yy,
            // v_xx and v_yy)

            // Now use the computed weak solution and the intersection info to actually apply the jump

            if (integrate_tangential_force && !intersection_ref_coords_u.empty())
            {

                // Evaluate the jump conditions and apply them to the Eulerian
                // grid.
                // This set of intersection points will take care of the jump in the pressure
                // as well as the jump in the viscous term for the  following components
                // 1) u_xx (in the x-momentum equation)
                // 2) v_yy (in the y-momentum equation)
                // 3) w_zz (in the z-momentum equation)
				TBOX_ASSERT(intersection_indices_um.size() == intersection_indices_up.size());

                X_fe_base->attach_quadrature_rule(qrule.get());
                X_fe_base->reinit(elem, &intersection_ref_coords_u);

                const size_t n_qp = intersection_ref_coords_u.size();

                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {

                    					
					libMesh::Point& xx = intersection_coords_u[qp];
					libMesh::Point& XX = intersection_ref_coords_u[qp];

					
                    const SideIndex<NDIM>& i_s_um = intersection_indices_um[qp];
                    const unsigned int axis = i_s_um.getAxis();
                    const SideIndex<NDIM>& i_s_up = intersection_indices_up[qp];

                    if (side_boxes[axis].contains(i_s_um) && side_boxes[axis].contains(i_s_up))
                    {
                        TBOX_ASSERT(i_s_um.getAxis() == i_s_up.getAxis());
                        TBOX_ASSERT(i_s_up(axis) - i_s_um(axis) == 1);
                        const double x_cell_bdry_um =
                            x_lower[axis] + static_cast<double>(i_s_um(axis) - patch_lower[axis]) * dx[axis];
                            
                        const double x_cell_bdry_up =
                            x_lower[axis] + static_cast<double>(i_s_up(axis) - patch_lower[axis]) * dx[axis];
                        const double sdh_um = ((xx(axis) - x_cell_bdry_um)); // Signed Distance h

                        const double sdh_up = ((xx(axis) - x_cell_bdry_up)); // Signed Distance h

                        TBOX_ASSERT((sdh_um) < dx[axis] && sdh_um > 0);
                        TBOX_ASSERT(fabs(sdh_up) < dx[axis] && sdh_up < 0);
                        
                        
                        
                        
                        
                        
                    interpolate(&X(0), qp, X0_node, X_phi);
					interpolate(&Tau1(0), qp, X0_node, X_dphi_dxi);
					if (dim == 1)
						Tau2 = VectorValue<double>(0.0, 0.0, 1.0);
					else
						interpolate(&Tau2(0), qp, X0_node, X_dphi_deta);
					N = Tau1.cross(Tau2);

					interpolate(&x(0), qp, X_node, X_phi);
					interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
					tau1unit = tau1.unit();

				if (dim == 1)
					tau2 = VectorValue<double>(0.0, 0.0, 1.0);
				else
					interpolate(&tau2(0), qp, X_node, X_dphi_deta);
	
					n = tau1.cross(tau2);
					tau2unit = tau2.unit();
					dA_da = sqrt(N * N) / sqrt(n * n);
					N = N.unit();
					n = n.unit();


					double* nn;

					nn = &n(0);
					
					F.zero();

					
					for (unsigned int dod = 0; dod < NDIM; ++dod)
					{                            
					     F(dod) +=  1000 * (X(dod) - x(dod));
					                           
					}
					
					F *=dA_da;
					

                        
                        
                        
                        
                        
                        
                        
                        
                        
                        

                        double C_u_um = 0;
                        double C_u_up = 0;
#if (NDIM == 2)

                        if (axis == 0)
                        {
							interpolate(&jn(0), qp, du_j_node, X_phi);
							interpolate(&jnn(0), qp, d2u_j_node, X_phi);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(0);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(0);
                        }
                        else if (axis == 1)
                        {
                            interpolate(&jn(0), qp, dv_j_node, X_phi);
							interpolate(&jnn(0), qp, d2v_j_node, X_phi);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(1);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);//  + 0.5 * sdh_um * sdh_um * jnn(1);
                        }
#endif

#if (NDIM == 3)
                        if (axis == 0)
                        {
                            interpolate(&jn(0), qp, du_j_node, X_phi);
							interpolate(&jnn(0), qp, d2u_j_node, X_phi);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(0);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(0);
                        }
                        else if (axis == 1)
                        {
                            interpolate(&jn(0), qp, dv_j_node, X_phi);
			    interpolate(&jnn(0), qp, d2v_j_node, X_phi);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(1);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(1);
                        }
                        else if (axis == 2)
                        {
                            interpolate(&jn(0), qp, dw_j_node, X_phi);
							interpolate(&jnn(0), qp, d2w_j_node, X_phi);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(2);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(2);
                        }
#endif
                        (*f_data)(i_s_up) +=  (n(axis) > 0.0 ? 1.0 : -1.0) * (C_u_um / (dx[axis] * dx[axis]));
                        (*f_data)(i_s_um) +=  (n(axis) > 0.0 ? -1.0 : 1.0) * (C_u_up / (dx[axis] * dx[axis]));
                        
                    }
                }
            }

 

#if (NDIM == 2)
            // Apply the jumps of the viscous term on the direction of cell sides
            // For the velocity field U=(u,v) this will take care of corrections
            // in the following components done in a direction by direction manner:
            // 1) u_yy, u_zz  (in the x-momentum equation)
            // 2) v_xx, v_zz  (in the y-momentum equation)
           

            if (integrate_tangential_force && !intersectionSide_ref_coords_u.empty())
            {
				 TBOX_ASSERT(intersectionSide_indices_um.size() == intersectionSide_indices_up.size());
                X_fe_base->reinit(elem, &intersectionSide_ref_coords_u);
                const size_t n_qp = intersectionSide_indices_up.size();

                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    const SideIndex<NDIM>& i_side_up = intersectionSide_indices_up[qp];
                    
                    const unsigned int dd = i_side_up.getAxis();
                    
                    const unsigned int axis = (dd == 0 ? 1 : 0);

                    libMesh::Point& xx = intersectionSide_coords_u[qp];
					libMesh::Point& XX = intersectionSide_ref_coords_u[qp];

                    const SideIndex<NDIM>& i_side_um = intersectionSide_indices_um[qp];
               

                    if (side_boxes[dd].contains(i_side_up) && side_boxes[dd].contains(i_side_um) )
                    {
                        // Impose the jump conditions.
                        
                        TBOX_ASSERT( i_side_up(axis) - i_side_um(axis) == 1);
                        
                        const double x_mid_side_up =
                            x_lower[axis] + static_cast<double>(i_side_up(axis) - patch_lower[axis] + 0.5) * dx[axis];
                            
                        const double x_mid_side_um =
                            x_lower[axis] + static_cast<double>(i_side_um(axis) - patch_lower[axis] + 0.5) * dx[axis];

                        TBOX_ASSERT(xx(axis) <= x_mid_side_up);
                        TBOX_ASSERT(xx(axis) > x_mid_side_um);

                        const double sdh_up = xx(axis) - x_mid_side_up; // Signed Distance h
                        const double sdh_um = xx(axis) - x_mid_side_um;

                        double C_u_um = 0;
                        double C_u_up = 0;
                        
                        
                        
                        
                        
                        
                    interpolate(&X(0), qp, X0_node, X_phi);
					interpolate(&Tau1(0), qp, X0_node, X_dphi_dxi);
					if (dim == 1)
						Tau2 = VectorValue<double>(0.0, 0.0, 1.0);
					else
						interpolate(&Tau2(0), qp, X0_node, X_dphi_deta);
					N = Tau1.cross(Tau2);

					interpolate(&x(0), qp, X_node, X_phi);
					interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
					tau1unit = tau1.unit();

					if (dim == 1)
						tau2 = VectorValue<double>(0.0, 0.0, 1.0);
					else
						interpolate(&tau2(0), qp, X_node, X_dphi_deta);
	
					n = tau1.cross(tau2);
					tau2unit = tau2.unit();
					dA_da = sqrt(N * N) / sqrt(n * n);
					N = N.unit();
					n = n.unit();


					double* nn;

					nn = &n(0);
					
					F.zero();
					

					
					for (unsigned int dod = 0; dod < NDIM; ++dod)
					{                            
					     F(dod) +=  1000 * (X(dod) - x(dod));
					                           
					}
					
					F *=dA_da;
					
                        
                        
                        
                        
                        
                        
                        
                        

                        if (dd == 0)
                        {
                            interpolate(&jn(0), qp, du_j_node, X_phi);
							interpolate(&jnn(0), qp, d2u_j_node, X_phi);
                            C_u_um = sdh_um * (F(1) - F*n*n(1))*n(1);// + 0.5 * sdh_um * sdh_um * jnn(1) ;
                            C_u_up = sdh_up * (F(1) - F*n*n(1))*n(1);// + 0.5 * sdh_up * sdh_up * jnn(1) ;
                        }
                        else
                        {
                            interpolate(&jn(0), qp, dv_j_node, X_phi);
							interpolate(&jnn(0), qp, d2v_j_node, X_phi);
                            C_u_um = sdh_um * (F(0) - F*n*n(0))*n(0);// + 0.5 * sdh_um * sdh_um * jnn(0);
                            C_u_up = sdh_up * (F(0) - F*n*n(0))*n(0);// + 0.5 * sdh_up * sdh_up * jnn(0);
                            
                        }
						 
                       (*f_data)(i_side_um) +=  (n(axis) > 0.0 ? -1.0 : 1.0) * (C_u_up / (dx[axis] * dx[axis]));
                       
                       (*f_data)(i_side_up) +=  (n(axis) > 0.0 ? 1.0 : -1.0) * (C_u_um / (dx[axis] * dx[axis]));


                    }
                }
            }

#endif

#if (NDIM == 3)

            if (integrate_tangential_force && !intersectionSide_ref_coords_u.empty())
            {
                TBOX_ASSERT(intersectionSide_indices_um.size() == intersectionSide_indices_up.size());
                X_fe_base->reinit(elem, &intersectionSide_ref_coords_u);
                const size_t n_qp = intersectionSide_indices_up.size();

                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    libMesh::Point& xx = intersectionSide_coords_u[qp];
                    libMesh::Point& XX = intersectionSide_ref_coords_u[qp];
                    const SideIndex<NDIM>& i_side_up = intersectionSide_indices_up[qp];
                    const SideIndex<NDIM>& i_side_um = intersectionSide_indices_um[qp];
                    const unsigned int dd = i_side_up.getAxis();
                    const unsigned int ddd = i_side_um.getAxis();
                    TBOX_ASSERT(dd == ddd);

                    
                    if (side_boxes[dd].contains(i_side_up) && side_boxes[dd].contains(i_side_um))
                    {
                        int axis;

                        if (dd == 0)
                        {
                            axis = 2;
                        }
                        else if (dd == 1)
                        {
                            axis = 0;
                        }
                        else
                        {
                            axis = 1;
                        }

                        // Impose the jump conditions.

                        TBOX_ASSERT(i_side_up(axis) - i_side_um(axis) == 1);

                        double x_mid_side_up =
                            x_lower[axis] + static_cast<double>(i_side_up(axis) - patch_lower[axis] + 0.5) * dx[axis];

                        double x_mid_side_um =
                            x_lower[axis] + static_cast<double>(i_side_um(axis) - patch_lower[axis] + 0.5) * dx[axis];

                        TBOX_ASSERT(xx(axis) <= x_mid_side_up);
                        TBOX_ASSERT(xx(axis) > x_mid_side_um);

                        double sdh_up = xx(axis) - x_mid_side_up; // Signed Distance h
                        double sdh_um = xx(axis) - x_mid_side_um;

                        double C_u_um, C_u_up;
                        
                        
                                                
                    interpolate(&X(0), qp, X0_node, X_phi);
					interpolate(&Tau1(0), qp, X0_node, X_dphi_dxi);
					if (dim == 1)
						Tau2 = VectorValue<double>(0.0, 0.0, 1.0);
					else
						interpolate(&Tau2(0), qp, X0_node, X_dphi_deta);
					N = Tau1.cross(Tau2);

					interpolate(&x(0), qp, X_node, X_phi);
					interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
					tau1unit = tau1.unit();

					if (dim == 1)
						tau2 = VectorValue<double>(0.0, 0.0, 1.0);
					else
						interpolate(&tau2(0), qp, X_node, X_dphi_deta);
	
					n = tau1.cross(tau2);
					tau2unit = tau2.unit();
					dA_da = sqrt(N * N) / sqrt(n * n);
					N = N.unit();
					n = n.unit();


					double* nn;

					nn = &n(0);
					

					F.zero();


					for (unsigned int dod = 0; dod < NDIM; ++dod)
					{                            
					     F(dod) +=  1000 * (X(dod) - x(dod));
					                           
					}
					
					F *=dA_da;
					
					
					
                        

                        if (dd == 0)
                        {
                            interpolate(&jn(0), qp, du_j_node, X_phi);
							interpolate(&jnn(0), qp, d2u_j_node, X_phi);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(axis);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(axis);
                        }
                        else if (dd == 1)
                        {
                            interpolate(&jn(0), qp, dv_j_node, X_phi);
							interpolate(&jnn(0), qp, d2v_j_node, X_phi);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(axis);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(axis);
                        }
                        else
                        {
                            interpolate(&jn(0), qp, dw_j_node, X_phi);
							interpolate(&jn(0), qp, d2w_j_node, X_phi);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(axis);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(axis);
                        }

                        (*f_data)(i_side_um) += - (n(axis) > 0.0 ? -1.0 : 1.0) * (C_u_up / (dx[axis] * dx[axis]));

                        (*f_data)(i_side_up) += - (n(axis) > 0.0 ? 1.0 : -1.0) * (C_u_um / (dx[axis] * dx[axis]));
                    }
                }
            }

            if (integrate_tangential_force && !intersectionSide2_ref_coords_u.empty())
            {
                TBOX_ASSERT(intersectionSide2_indices_um.size() == intersectionSide2_indices_up.size());
                X_fe_base->reinit(elem, &intersectionSide2_ref_coords_u);
                const size_t n_qp = intersectionSide2_indices_up.size();

                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    const SideIndex<NDIM>& i_side_up = intersectionSide2_indices_up[qp];

                    const SideIndex<NDIM>& i_side_um = intersectionSide2_indices_um[qp];

                    const unsigned int dd = i_side_up.getAxis();
                    const unsigned int ddd = i_side_um.getAxis();
                    TBOX_ASSERT(dd == ddd);

                    libMesh::Point xx = intersectionSide2_coords_u[qp];
                    libMesh::Point& XX = intersectionSide2_ref_coords_u[qp];

                    if (side_boxes[dd].contains(i_side_up) && side_boxes[dd].contains(i_side_um))
                    {
                        int axis;

                        if (dd == 0)
                        {
                            axis = 1;
                        }
                        else if (dd == 1)
                        {
                            axis = 2;
                        }
                        else
                        {
                            axis = 0;
                        }

                        // Impose the jump conditions.

                        TBOX_ASSERT(i_side_up(axis) - i_side_um(axis) == 1);

                        double x_mid_side_up =
                            x_lower[axis] + static_cast<double>(i_side_up(axis) - patch_lower[axis] + 0.5) * dx[axis];

                        double x_mid_side_um =
                            x_lower[axis] + static_cast<double>(i_side_um(axis) - patch_lower[axis] + 0.5) * dx[axis];

                        TBOX_ASSERT(xx(axis) <= x_mid_side_up);
                        TBOX_ASSERT(xx(axis) > x_mid_side_um);

                        double sdh_up = xx(axis) - x_mid_side_up; // Signed Distance h
                        double sdh_um = xx(axis) - x_mid_side_um;

                        double C_u_um, C_u_up;
                        
                        
                        
                        
                        
                        
                        
                        
                        
                    interpolate(&X(0), qp, X0_node, X_phi);
					interpolate(&Tau1(0), qp, X0_node, X_dphi_dxi);
					if (dim == 1)
						Tau2 = VectorValue<double>(0.0, 0.0, 1.0);
					else
						interpolate(&Tau2(0), qp, X0_node, X_dphi_deta);
					N = Tau1.cross(Tau2);

					interpolate(&x(0), qp, X_node, X_phi);
					interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
					tau1unit = tau1.unit();

					if (dim == 1)
						tau2 = VectorValue<double>(0.0, 0.0, 1.0);
					else
						interpolate(&tau2(0), qp, X_node, X_dphi_deta);
	
					n = tau1.cross(tau2);
					tau2unit = tau2.unit();
					dA_da = sqrt(N * N) / sqrt(n * n);
					N = N.unit();
					n = n.unit();


					double* nn;

					nn = &n(0);
					F.zero();
					



					for (unsigned int dod = 0; dod < NDIM; ++dod)
					{                            
					     F(dod) +=  1000 * (X(dod) - x(dod));
					                           
					}
					
					F *=dA_da;
                        
                        
                        
                        
                        
                        
                        
                        
                        

                        if (dd == 0)
                        {
                            interpolate(&jn(0), qp, du_j_node, X_phi);
							interpolate(&jnn(0), qp, du_j_node, X_phi);
                            C_u_um = sdh_um *  (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(axis);
                            C_u_up = sdh_up *  (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(axis);
                        }
                        else if (dd == 1)
                        {
                            interpolate(&jn(0), qp, dv_j_node, X_phi);
							interpolate(&jnn(0), qp, d2v_j_node, X_phi);
                            C_u_um = sdh_um *  (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(axis);
                            C_u_up = sdh_up *  (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(axis);
                        }
                        else
                        {
                            interpolate(&jn(0), qp, dw_j_node, X_phi);
							interpolate(&jnn(0), qp, d2w_j_node, X_phi);
                            C_u_um = sdh_um *  (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(axis);
                            C_u_up = sdh_up *  (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(axis);
                        }

                        (*f_data)(i_side_um) +=  (n(axis) > 0.0 ? -1.0 : 1.0) * (C_u_up / (dx[axis] * dx[axis]));

                        (*f_data)(i_side_up) +=  (n(axis) > 0.0 ? 1.0 : -1.0) * (C_u_um / (dx[axis] * dx[axis]));


                    }
                }
            }

#endif
            // Applying the jump for viscous term (in 2D u_xx, u_yy,
            // v_xx and v_yy)

            // Now use the computed weak solution and the intersection info to actually apply the jump

            if (integrate_tangential_force && !intersection_ref_coords_u.empty())
            {

                // Evaluate the jump conditions and apply them to the Eulerian
                // grid.
                // This set of intersection points will take care of the jump in the pressure
                // as well as the jump in the viscous term for the  following components
                // 1) u_xx (in the x-momentum equation)
                // 2) v_yy (in the y-momentum equation)
                // 3) w_zz (in the z-momentum equation)
				TBOX_ASSERT(intersection_indices_um.size() == intersection_indices_up.size());

                X_fe_base->attach_quadrature_rule(qrule.get());
                X_fe_base->reinit(elem, &intersection_ref_coords_u);

                const size_t n_qp = intersection_ref_coords_u.size();

                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {

                    					
					libMesh::Point& xx = intersection_coords_u[qp];
					libMesh::Point& XX = intersection_ref_coords_u[qp];

					
                    const SideIndex<NDIM>& i_s_um = intersection_indices_um[qp];
                    const unsigned int axis = i_s_um.getAxis();
                    const SideIndex<NDIM>& i_s_up = intersection_indices_up[qp];

                    if (side_boxes[axis].contains(i_s_um) && side_boxes[axis].contains(i_s_up))
                    {
                        TBOX_ASSERT(i_s_um.getAxis() == i_s_up.getAxis());
                        TBOX_ASSERT(i_s_up(axis) - i_s_um(axis) == 1);
                        const double x_cell_bdry_um =
                            x_lower[axis] + static_cast<double>(i_s_um(axis) - patch_lower[axis]) * dx[axis];
                            
                        const double x_cell_bdry_up =
                            x_lower[axis] + static_cast<double>(i_s_up(axis) - patch_lower[axis]) * dx[axis];
                        const double sdh_um = ((xx(axis) - x_cell_bdry_um)); // Signed Distance h

                        const double sdh_up = ((xx(axis) - x_cell_bdry_up)); // Signed Distance h

                        TBOX_ASSERT((sdh_um) < dx[axis] && sdh_um > 0);
                        TBOX_ASSERT(fabs(sdh_up) < dx[axis] && sdh_up < 0);
                        
                        
                        
                        
                        
                        
                    interpolate(&X(0), qp, X0_node, X_phi);
					interpolate(&Tau1(0), qp, X0_node, X_dphi_dxi);
					if (dim == 1)
						Tau2 = VectorValue<double>(0.0, 0.0, 1.0);
					else
						interpolate(&Tau2(0), qp, X0_node, X_dphi_deta);
					N = Tau1.cross(Tau2);

					interpolate(&x(0), qp, X_node, X_phi);
					interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
					tau1unit = tau1.unit();

				if (dim == 1)
					tau2 = VectorValue<double>(0.0, 0.0, 1.0);
				else
					interpolate(&tau2(0), qp, X_node, X_dphi_deta);
	
					n = tau1.cross(tau2);
					tau2unit = tau2.unit();
					dA_da = sqrt(N * N) / sqrt(n * n);
					N = N.unit();
					n = n.unit();


					double* nn;

					nn = &n(0);
					
					F.zero();


					for (unsigned int dod = 0; dod < NDIM; ++dod)
					{                            
					     F(dod) +=  1000 * (X(dod) - x(dod));
					                           
					}
					
					//~ if (d_lag_force_fcn_data[part].fcn)
					//~ {
						//~ 
						//~ // Compute the value of the body force at the quadrature
						//~ // point and add the corresponding forces to the
						//~ // right-hand-side vector.
						//~ fe.setInterpolatedDataPointers(force_var_data, force_grad_var_data, force_fcn_system_idxs, elem, qp);
						//~ d_lag_force_fcn_data[part].fcn(F_s,
													   //~ FF,
													   //~ xx,
													   //~ XX,
													   //~ elem,
													   //~ nn,
													   //~ force_var_data,
													   //~ force_grad_var_data,
													   //~ data_time,
					                                   //~ d_lag_force_fcn_data[part].ctx);
					     //~ F += F_s;
					                             //~ 
					//~ }
					
					F *=dA_da;
					

                        
                        
                        
                        
                        
                        
                        
                        
                        
                        

                        double C_u_um = 0;
                        double C_u_up = 0;
#if (NDIM == 2)

                        if (axis == 0)
                        {
							interpolate(&jn(0), qp, du_j_node, X_phi);
							interpolate(&jnn(0), qp, d2u_j_node, X_phi);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(0);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(0);
                        }
                        else if (axis == 1)
                        {
                            interpolate(&jn(0), qp, dv_j_node, X_phi);
							interpolate(&jnn(0), qp, d2v_j_node, X_phi);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(1);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);//  + 0.5 * sdh_um * sdh_um * jnn(1);
                        }
#endif

#if (NDIM == 3)
                        if (axis == 0)
                        {
                            interpolate(&jn(0), qp, du_j_node, X_phi);
							interpolate(&jnn(0), qp, d2u_j_node, X_phi);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(0);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(0);
                        }
                        else if (axis == 1)
                        {
                            interpolate(&jn(0), qp, dv_j_node, X_phi);
			    interpolate(&jnn(0), qp, d2v_j_node, X_phi);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(1);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(1);
                        }
                        else if (axis == 2)
                        {
                            interpolate(&jn(0), qp, dw_j_node, X_phi);
							interpolate(&jnn(0), qp, d2w_j_node, X_phi);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(2);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(2);
                        }
#endif
                        (*f_data)(i_s_up) +=  (n(axis) > 0.0 ? 1.0 : -1.0) * (C_u_um / (dx[axis] * dx[axis]));
                        (*f_data)(i_s_um) +=  (n(axis) > 0.0 ? -1.0 : 1.0) * (C_u_up / (dx[axis] * dx[axis]));
                        
                    }
                }
            }

 

#if (NDIM == 2)
            // Apply the jumps of the viscous term on the direction of cell sides
            // For the velocity field U=(u,v) this will take care of corrections
            // in the following components done in a direction by direction manner:
            // 1) u_yy, u_zz  (in the x-momentum equation)
            // 2) v_xx, v_zz  (in the y-momentum equation)
           

            if (integrate_tangential_force && !intersectionSide_ref_coords_u.empty())
            {
				 TBOX_ASSERT(intersectionSide_indices_um.size() == intersectionSide_indices_up.size());
                X_fe_base->reinit(elem, &intersectionSide_ref_coords_u);
                const size_t n_qp = intersectionSide_indices_up.size();

                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    const SideIndex<NDIM>& i_side_up = intersectionSide_indices_up[qp];
                    
                    const unsigned int dd = i_side_up.getAxis();
                    
                    const unsigned int axis = (dd == 0 ? 1 : 0);

                    libMesh::Point& xx = intersectionSide_coords_u[qp];
					libMesh::Point& XX = intersectionSide_ref_coords_u[qp];

                    const SideIndex<NDIM>& i_side_um = intersectionSide_indices_um[qp];
               

                    if (side_boxes[dd].contains(i_side_up) && side_boxes[dd].contains(i_side_um) )
                    {
                        // Impose the jump conditions.
                        
                        TBOX_ASSERT( i_side_up(axis) - i_side_um(axis) == 1);
                        
                        const double x_mid_side_up =
                            x_lower[axis] + static_cast<double>(i_side_up(axis) - patch_lower[axis] + 0.5) * dx[axis];
                            
                        const double x_mid_side_um =
                            x_lower[axis] + static_cast<double>(i_side_um(axis) - patch_lower[axis] + 0.5) * dx[axis];

                        TBOX_ASSERT(xx(axis) <= x_mid_side_up);
                        TBOX_ASSERT(xx(axis) > x_mid_side_um);

                        const double sdh_up = xx(axis) - x_mid_side_up; // Signed Distance h
                        const double sdh_um = xx(axis) - x_mid_side_um;

                        double C_u_um = 0;
                        double C_u_up = 0;
                        
                        
                        
                        
                        
                        
                    interpolate(&X(0), qp, X0_node, X_phi);
					interpolate(&Tau1(0), qp, X0_node, X_dphi_dxi);
					if (dim == 1)
						Tau2 = VectorValue<double>(0.0, 0.0, 1.0);
					else
						interpolate(&Tau2(0), qp, X0_node, X_dphi_deta);
					N = Tau1.cross(Tau2);

					interpolate(&x(0), qp, X_node, X_phi);
					interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
					tau1unit = tau1.unit();

					if (dim == 1)
						tau2 = VectorValue<double>(0.0, 0.0, 1.0);
					else
						interpolate(&tau2(0), qp, X_node, X_dphi_deta);
	
					n = tau1.cross(tau2);
					tau2unit = tau2.unit();
					dA_da = sqrt(N * N) / sqrt(n * n);
					N = N.unit();
					n = n.unit();


					double* nn;

					nn = &n(0);
					
					F.zero();
					


					for (unsigned int dod = 0; dod < NDIM; ++dod)
					{                            
					     F(dod) +=  1000 * (XX(dod) - xx(dod));
					                           
					}
					
					//~ if (d_lag_force_fcn_data[part].fcn)
					//~ {
						//~ // Compute the value of the body force at the quadrature
						//~ // point and add the corresponding forces to the
						//~ // right-hand-side vector.
						//~ fe.setInterpolatedDataPointers(force_var_data, force_grad_var_data, force_fcn_system_idxs, elem, qp);
						//~ d_lag_force_fcn_data[part].fcn(F_s,
													   //~ FF,
													   //~ xx,
													   //~ XX,
													   //~ elem,
													   //~ nn,
													   //~ force_var_data,
													   //~ force_grad_var_data,
													   //~ data_time,
					                                   //~ d_lag_force_fcn_data[part].ctx);
					     //~ F += F_s;
					                             //~ 
					//~ }
					
					F *=dA_da;
					
                        
                        
                        
                        
                        
                        
                        
                        

                        if (dd == 0)
                        {
                            interpolate(&jn(0), qp, du_j_node, X_phi);
							interpolate(&jnn(0), qp, d2u_j_node, X_phi);
                            C_u_um = sdh_um * (F(1) - F*n*n(1))*n(1);// + 0.5 * sdh_um * sdh_um * jnn(1) ;
                            C_u_up = sdh_up * (F(1) - F*n*n(1))*n(1);// + 0.5 * sdh_up * sdh_up * jnn(1) ;
                        }
                        else
                        {
                            interpolate(&jn(0), qp, dv_j_node, X_phi);
							interpolate(&jnn(0), qp, d2v_j_node, X_phi);
                            C_u_um = sdh_um * (F(0) - F*n*n(0))*n(0);// + 0.5 * sdh_um * sdh_um * jnn(0);
                            C_u_up = sdh_up * (F(0) - F*n*n(0))*n(0);// + 0.5 * sdh_up * sdh_up * jnn(0);
                            
                        }
						 
                       (*f_data)(i_side_um) +=  (n(axis) > 0.0 ? -1.0 : 1.0) * (C_u_up / (dx[axis] * dx[axis]));
                       
                       (*f_data)(i_side_up) +=  (n(axis) > 0.0 ? 1.0 : -1.0) * (C_u_um / (dx[axis] * dx[axis]));


                    }
                }
            }

#endif

#if (NDIM == 3)

            if (integrate_tangential_force && !intersectionSide_ref_coords_u.empty())
            {
                TBOX_ASSERT(intersectionSide_indices_um.size() == intersectionSide_indices_up.size());
                X_fe_base->reinit(elem, &intersectionSide_ref_coords_u);
                const size_t n_qp = intersectionSide_indices_up.size();

                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    libMesh::Point& xx = intersectionSide_coords_u[qp];
                    libMesh::Point& XX = intersectionSide_ref_coords_u[qp];
                    const SideIndex<NDIM>& i_side_up = intersectionSide_indices_up[qp];
                    const SideIndex<NDIM>& i_side_um = intersectionSide_indices_um[qp];
                    const unsigned int dd = i_side_up.getAxis();
                    const unsigned int ddd = i_side_um.getAxis();
                    TBOX_ASSERT(dd == ddd);

                    
                    if (side_boxes[dd].contains(i_side_up) && side_boxes[dd].contains(i_side_um))
                    {
                        int axis;

                        if (dd == 0)
                        {
                            axis = 2;
                        }
                        else if (dd == 1)
                        {
                            axis = 0;
                        }
                        else
                        {
                            axis = 1;
                        }

                        // Impose the jump conditions.

                        TBOX_ASSERT(i_side_up(axis) - i_side_um(axis) == 1);

                        double x_mid_side_up =
                            x_lower[axis] + static_cast<double>(i_side_up(axis) - patch_lower[axis] + 0.5) * dx[axis];

                        double x_mid_side_um =
                            x_lower[axis] + static_cast<double>(i_side_um(axis) - patch_lower[axis] + 0.5) * dx[axis];

                        TBOX_ASSERT(xx(axis) <= x_mid_side_up);
                        TBOX_ASSERT(xx(axis) > x_mid_side_um);

                        double sdh_up = xx(axis) - x_mid_side_up; // Signed Distance h
                        double sdh_um = xx(axis) - x_mid_side_um;

                        double C_u_um, C_u_up;
                        
                        
                                                
                    interpolate(&X(0), qp, X0_node, X_phi);
					interpolate(&Tau1(0), qp, X0_node, X_dphi_dxi);
					if (dim == 1)
						Tau2 = VectorValue<double>(0.0, 0.0, 1.0);
					else
						interpolate(&Tau2(0), qp, X0_node, X_dphi_deta);
					N = Tau1.cross(Tau2);

					interpolate(&x(0), qp, X_node, X_phi);
					interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
					tau1unit = tau1.unit();

					if (dim == 1)
						tau2 = VectorValue<double>(0.0, 0.0, 1.0);
					else
						interpolate(&tau2(0), qp, X_node, X_dphi_deta);
	
					n = tau1.cross(tau2);
					tau2unit = tau2.unit();
					dA_da = sqrt(N * N) / sqrt(n * n);
					N = N.unit();
					n = n.unit();


					double* nn;

					nn = &n(0);
					

					F.zero();
					
					
					for (unsigned int dod = 0; dod < NDIM; ++dod)
					{                            
					     F(dod) +=  1000 * (XX(dod) - xx(dod));
					                           
					}
					//~ if (d_lag_force_fcn_data[part].fcn)
					//~ {
						//~ // Compute the value of the body force at the quadrature
						//~ // point and add the corresponding forces to the
						//~ // right-hand-side vector.
						//~ fe.setInterpolatedDataPointers(force_var_data, force_grad_var_data, force_fcn_system_idxs, elem, qp);
						//~ d_lag_force_fcn_data[part].fcn(F_s,
													   //~ FF,
													   //~ xx,
													   //~ XX,
													   //~ elem,
													   //~ nn,
													   //~ force_var_data,
													   //~ force_grad_var_data,
													   //~ data_time,
					                                   //~ d_lag_force_fcn_data[part].ctx);
					     //~ F += F_s;
					                             //~ 
					//~ }
					
					F *=dA_da;
					
					
					
                        

                        if (dd == 0)
                        {
                            interpolate(&jn(0), qp, du_j_node, X_phi);
							interpolate(&jnn(0), qp, d2u_j_node, X_phi);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(axis);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(axis);
                        }
                        else if (dd == 1)
                        {
                            interpolate(&jn(0), qp, dv_j_node, X_phi);
							interpolate(&jnn(0), qp, d2v_j_node, X_phi);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(axis);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(axis);
                        }
                        else
                        {
                            interpolate(&jn(0), qp, dw_j_node, X_phi);
							interpolate(&jn(0), qp, d2w_j_node, X_phi);
                            C_u_um = sdh_um * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(axis);
                            C_u_up = sdh_up * (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(axis);
                        }

                        (*f_data)(i_side_um) += - (n(axis) > 0.0 ? -1.0 : 1.0) * (C_u_up / (dx[axis] * dx[axis]));

                        (*f_data)(i_side_up) += - (n(axis) > 0.0 ? 1.0 : -1.0) * (C_u_um / (dx[axis] * dx[axis]));
                    }
                }
            }

            if (integrate_tangential_force && !intersectionSide2_ref_coords_u.empty())
            {
                TBOX_ASSERT(intersectionSide2_indices_um.size() == intersectionSide2_indices_up.size());
                X_fe_base->reinit(elem, &intersectionSide2_ref_coords_u);
                const size_t n_qp = intersectionSide2_indices_up.size();

                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    const SideIndex<NDIM>& i_side_up = intersectionSide2_indices_up[qp];

                    const SideIndex<NDIM>& i_side_um = intersectionSide2_indices_um[qp];

                    const unsigned int dd = i_side_up.getAxis();
                    const unsigned int ddd = i_side_um.getAxis();
                    TBOX_ASSERT(dd == ddd);

                    libMesh::Point xx = intersectionSide2_coords_u[qp];
                    libMesh::Point& XX = intersectionSide2_ref_coords_u[qp];

                    if (side_boxes[dd].contains(i_side_up) && side_boxes[dd].contains(i_side_um))
                    {
                        int axis;

                        if (dd == 0)
                        {
                            axis = 1;
                        }
                        else if (dd == 1)
                        {
                            axis = 2;
                        }
                        else
                        {
                            axis = 0;
                        }

                        // Impose the jump conditions.

                        TBOX_ASSERT(i_side_up(axis) - i_side_um(axis) == 1);

                        double x_mid_side_up =
                            x_lower[axis] + static_cast<double>(i_side_up(axis) - patch_lower[axis] + 0.5) * dx[axis];

                        double x_mid_side_um =
                            x_lower[axis] + static_cast<double>(i_side_um(axis) - patch_lower[axis] + 0.5) * dx[axis];

                        TBOX_ASSERT(xx(axis) <= x_mid_side_up);
                        TBOX_ASSERT(xx(axis) > x_mid_side_um);

                        double sdh_up = xx(axis) - x_mid_side_up; // Signed Distance h
                        double sdh_um = xx(axis) - x_mid_side_um;

                        double C_u_um, C_u_up;
                        
                        
                        
                        
                        
                        
                        
                        
                        
                    interpolate(&X(0), qp, X0_node, X_phi);
					interpolate(&Tau1(0), qp, X0_node, X_dphi_dxi);
					if (dim == 1)
						Tau2 = VectorValue<double>(0.0, 0.0, 1.0);
					else
						interpolate(&Tau2(0), qp, X0_node, X_dphi_deta);
					N = Tau1.cross(Tau2);

					interpolate(&x(0), qp, X_node, X_phi);
					interpolate(&tau1(0), qp, X_node, X_dphi_dxi);
					tau1unit = tau1.unit();

					if (dim == 1)
						tau2 = VectorValue<double>(0.0, 0.0, 1.0);
					else
						interpolate(&tau2(0), qp, X_node, X_dphi_deta);
	
					n = tau1.cross(tau2);
					tau2unit = tau2.unit();
					dA_da = sqrt(N * N) / sqrt(n * n);
					N = N.unit();
					n = n.unit();


					double* nn;

					nn = &n(0);
					F.zero();
					
					
					for (unsigned int dod = 0; dod < NDIM; ++dod)
					{                            
					     F(dod) +=  1000 * (XX(dod) - xx(dod));
					                           
					}

					//~ if (d_lag_force_fcn_data[part].fcn)
					//~ {
						//~ // Compute the value of the body force at the quadrature
						//~ // point and add the corresponding forces to the
						//~ // right-hand-side vector.
						//~ fe.setInterpolatedDataPointers(force_var_data, force_grad_var_data, force_fcn_system_idxs, elem, qp);
						//~ d_lag_force_fcn_data[part].fcn(F_s,
													   //~ FF,
													   //~ xx,
													   //~ XX,
													   //~ elem,
													   //~ nn,
													   //~ force_var_data,
													   //~ force_grad_var_data,
													   //~ data_time,
					                                   //~ d_lag_force_fcn_data[part].ctx);
					     //~ F += F_s;
					                             //~ 
					//~ }
					
					F *=dA_da;
                        
                        
                        
                        
                        
                        
                        
                        
                        

                        if (dd == 0)
                        {
                            interpolate(&jn(0), qp, du_j_node, X_phi);
							interpolate(&jnn(0), qp, du_j_node, X_phi);
                            C_u_um = sdh_um *  (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(axis);
                            C_u_up = sdh_up *  (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(axis);
                        }
                        else if (dd == 1)
                        {
                            interpolate(&jn(0), qp, dv_j_node, X_phi);
							interpolate(&jnn(0), qp, d2v_j_node, X_phi);
                            C_u_um = sdh_um *  (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(axis);
                            C_u_up = sdh_up *  (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(axis);
                        }
                        else
                        {
                            interpolate(&jn(0), qp, dw_j_node, X_phi);
							interpolate(&jnn(0), qp, d2w_j_node, X_phi);
                            C_u_um = sdh_um *  (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_um * sdh_um * jnn(axis);
                            C_u_up = sdh_up *  (F(axis) - F*n*n(axis))*n(axis);// + 0.5 * sdh_up * sdh_up * jnn(axis);
                        }

                        (*f_data)(i_side_um) +=  (n(axis) > 0.0 ? -1.0 : 1.0) * (C_u_up / (dx[axis] * dx[axis]));

                        (*f_data)(i_side_up) +=  (n(axis) > 0.0 ? 1.0 : -1.0) * (C_u_um / (dx[axis] * dx[axis]));


                    }
                }
            }

#endif
        }
  
  
    }



    VecRestoreArray(X_local_vec, &X_local_soln);
    VecGhostRestoreLocalForm(X_global_vec, &X_local_vec);


    VecRestoreArray(X0_local_vec, &X0_local_soln);
    VecGhostRestoreLocalForm(X0_global_vec, &X0_local_vec);

    VecRestoreArray(P_j_local_vec, &P_j_local_soln);
    VecGhostRestoreLocalForm(P_j_global_vec, &P_j_local_vec);

    VecRestoreArray(dP_j_local_vec, &dP_j_local_soln);
    VecGhostRestoreLocalForm(dP_j_global_vec, &dP_j_local_vec);

    VecRestoreArray(du_j_local_vec, &du_j_local_soln);
    VecGhostRestoreLocalForm(du_j_global_vec, &du_j_local_vec);

    VecRestoreArray(dv_j_local_vec, &dv_j_local_soln);
    VecGhostRestoreLocalForm(dv_j_global_vec, &dv_j_local_vec);

    VecRestoreArray(dw_j_local_vec, &dw_j_local_soln);
    VecGhostRestoreLocalForm(dw_j_global_vec, &dw_j_local_vec);


    VecRestoreArray(d2u_j_local_vec, &d2u_j_local_soln);
    VecGhostRestoreLocalForm(d2u_j_global_vec, &d2u_j_local_vec);

    VecRestoreArray(d2v_j_local_vec, &d2v_j_local_soln);
    VecGhostRestoreLocalForm(d2v_j_global_vec, &d2v_j_local_vec);

    VecRestoreArray(d2w_j_local_vec, &d2w_j_local_soln);
    VecGhostRestoreLocalForm(d2w_j_global_vec, &d2w_j_local_vec);
    return;



    return;
} // imposeJumpConditionsPointWise

void
IBFEMethod::initializeCoordinates(const unsigned int part)
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
    X_system.solution->localize(*X_system.current_local_solution);

    // Keep track of the initial coordinates.
    System& X0_system = equation_systems->get_system(COORDS0_SYSTEM_NAME);
    NumericVector<double>& X0_coords = *X0_system.solution;
    X0_coords = X_coords;
    X0_system.solution->localize(*X0_system.current_local_solution);
    return;
} // initializeCoordinates

void
IBFEMethod::updateCoordinateMapping(const unsigned int part)
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
    dX_coords.close();
    return;
} // updateCoordinateMapping

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBFEMethod::commonConstructor(const std::string& object_name,
                              Pointer<Database> input_db,
                              const std::vector<libMesh::Mesh*>& meshes,
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
    const bool interp_use_consistent_mass_matrix = true;
    d_default_interp_spec = FEDataManager::InterpSpec(
        "IB_4", QGAUSS, INVALID_ORDER, use_adaptive_quadrature, point_density, interp_use_consistent_mass_matrix);
    d_default_spread_spec =
        FEDataManager::SpreadSpec("IB_4", QGAUSS, INVALID_ORDER, use_adaptive_quadrature, point_density);
    d_ghosts = 0;
    d_split_normal_force = false;
    d_split_tangential_force = false;
    d_use_jump_conditions = false;
    d_modify_vel_interp_jumps = false;
    d_use_higher_order_jump = false;
    d_vel_interp_width = 0.0;
    d_mu = 0.0;
    d_use_consistent_mass_matrix = true;
    d_do_log = false;

    d_fe_family.resize(d_num_parts, INVALID_FE);
    d_fe_order.resize(d_num_parts, INVALID_ORDER);
    d_default_quad_type.resize(d_num_parts, INVALID_Q_RULE);
    d_default_quad_order.resize(d_num_parts, INVALID_ORDER);
    // Initialize function data to NULL.
    d_coordinate_mapping_fcn_data.resize(d_num_parts);
    d_lag_force_fcn_data.resize(d_num_parts);

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
                       << "::IBFEMethod():\n"
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

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_half_time = std::numeric_limits<double>::quiet_NaN();

    // Keep track of the initialization state.
    d_fe_equation_systems_initialized = false;
    d_fe_data_initialized = false;
    d_is_initialized = false;
    return;
} // commonConstructor

void
IBFEMethod::getFromInput(Pointer<Database> db, bool /*is_from_restart*/)
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

    // Force computation settings.
    if (db->isBool("split_normal_force"))
        d_split_normal_force = db->getBool("split_normal_force");
    else if (db->isBool("split_forces"))
        d_split_normal_force = db->getBool("split_forces");
    if (db->isBool("split_tangential_force"))
        d_split_tangential_force = db->getBool("split_tangential_force");
    else if (db->isBool("split_forces"))
        d_split_tangential_force = db->getBool("split_forces");
    if (db->isBool("use_jump_conditions")) d_use_jump_conditions = db->getBool("use_jump_conditions");
    if (db->isDouble("vel_interp_width")) d_vel_interp_width = db->getDouble("vel_interp_width");
    if (db->isBool("modify_vel_interp_jumps")) d_modify_vel_interp_jumps = db->getBool("modify_vel_interp_jumps");
    if (db->isBool("use_higher_order_jump")) d_use_higher_order_jump = db->getBool("use_higher_order_jump");
    if (db->isBool("use_consistent_mass_matrix"))
        d_use_consistent_mass_matrix = db->getBool("use_consistent_mass_matrix");

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
IBFEMethod::getFromRestart()
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
    int ver = db->getInteger("IBFE_METHOD_VERSION");
    if (ver != IBFE_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    db->getIntegerArray("d_ghosts", d_ghosts, NDIM);
    d_split_normal_force = db->getBool("d_split_normal_force");
    d_split_tangential_force = db->getBool("d_split_tangential_force");
    d_use_jump_conditions = db->getBool("d_use_jump_conditions");
    d_modify_vel_interp_jumps = db->getBool("d_modify_interp_jumps");
    d_use_higher_order_jump = db->getBool("d_use_higher_order_jump");
    d_vel_interp_width = db->isDouble("vel_interp_width");
    d_use_consistent_mass_matrix = db->getBool("d_use_consistent_mass_matrix");
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
