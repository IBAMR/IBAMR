// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2021 by the IBAMR developers
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

#include "ibamr/FEMechanicsBase.h"
#include "ibamr/IBFEDirectForcingKinematics.h"
#include "ibamr/IBFEMethod.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/BoxPartitioner.h"
#include "ibtk/CartSideDoubleRT0Refine.h"
#include "ibtk/FEDataInterpolation.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/FEProjector.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/LibMeshSystemIBVectors.h"
#include "ibtk/PartitioningBox.h"
#include "ibtk/QuadratureCache.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/SAMRAIDataCache.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "GriddingAlgorithm.h"
#include "HierarchyCellDataOpsInteger.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "LoadBalanceStrategy.h"
#include "LoadBalancer.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchData.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefineSchedule.h"
#include "RefineTransactionFactory.h"
#include "SideData.h"
#include "SideIndex.h"
#include "StandardTagAndInitialize.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/InputDatabase.h"
#include "tbox/MathUtilities.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include "libmesh/boundary_info.h"
#include "libmesh/compare_types.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/edge.h"
#include "libmesh/elem.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_parallel_type.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/face.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/fem_context.h"
#include "libmesh/id_types.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_config.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameters.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/system.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_tensor.h"
#include "libmesh/type_vector.h"
#include "libmesh/variant_filter_iterator.h"
#include "libmesh/vector_value.h"

#include "petscvec.h"
#include <petsclog.h>
#include <petscsys.h>

#include "ibamr/namespaces.h" // IWYU pragma: keep

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace SAMRAI
{
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefinePatchStrategy;
} // namespace xfer
} // namespace SAMRAI

using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static Timer* t_preprocess_integrate_data;
static Timer* t_postprocess_integrate_data;
static Timer* t_interpolate_velocity;
static Timer* t_compute_lagrangian_force;
static Timer* t_spread_force;
static Timer* t_compute_lagrangian_fluid_source;
static Timer* t_spread_fluid_source;
static Timer* t_add_workload_estimate;
static Timer* t_begin_data_redistribution;
static Timer* t_end_data_redistribution;
static Timer* t_apply_gradient_detector;
// Version of IBFEMethod restart file data.
const int IBFE_METHOD_VERSION = 6;

inline boundary_id_type
get_dirichlet_bdry_ids(const std::vector<boundary_id_type>& bdry_ids)
{
    boundary_id_type dirichlet_bdry_ids = 0;
    for (const auto& bdry_id : bdry_ids)
    {
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
    std::vector<boundary_id_type> bdry_ids;
    boundary_info.boundary_ids(elem, side, bdry_ids);
    bool at_physical_bdry = !elem->neighbor_ptr(side);
    for (const auto& bdry_id : bdry_ids)
    {
        if (dof_map.is_periodic_boundary(bdry_id)) at_physical_bdry = false;
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
    std::vector<boundary_id_type> bdry_ids;
    boundary_info.boundary_ids(elem, side, bdry_ids);
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

const Real PENALTY = 1.e10;

void
assemble_poisson(EquationSystems& es, const std::string& /*system_name*/)
{
    const MeshBase& mesh = es.get_mesh();
    const BoundaryInfo& boundary_info = mesh.get_boundary_info();
    const unsigned int dim = mesh.mesh_dimension();
    auto& system = es.get_system<LinearImplicitSystem>(IBFEMethod::PRESSURE_SYSTEM_NAME);
    const DofMap& dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(0);

    std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
    QGauss qrule(dim, FIFTH);
    fe->attach_quadrature_rule(&qrule);
    std::unique_ptr<FEBase> fe_face(FEBase::build(dim, fe_type));
    QGauss qface(dim - 1, FIFTH);
    fe_face->attach_quadrature_rule(&qface);

    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

    const std::vector<std::vector<Real> >& phi_face = fe_face->get_phi();
    const std::vector<Real>& JxW_face = fe_face->get_JxW();

    DenseMatrix<Number> Ke;
    std::vector<dof_id_type> dof_indices;

    const double epsilon = es.parameters.get<Real>("Phi_epsilon");
    const double epsilon_inv = (std::abs(epsilon) > std::numeric_limits<double>::epsilon() ? 1.0 / epsilon : 0.0);

    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for (; el != end_el; ++el)
    {
        const Elem* elem = *el;
        dof_map.dof_indices(elem, dof_indices);
        fe->reinit(elem);
        auto Ke_size = static_cast<unsigned int>(dof_indices.size());
        Ke.resize(Ke_size, Ke_size);
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            for (unsigned int i = 0; i < phi.size(); i++)
            {
                for (unsigned int j = 0; j < phi.size(); j++)
                {
                    Ke(i, j) += (epsilon_inv * phi[i][qp] * phi[j][qp] + (dphi[i][qp] * dphi[j][qp])) * JxW[qp];
                }
            }
        }
        for (unsigned int side = 0; side < elem->n_sides(); side++)
        {
            if (is_physical_bdry(elem, side, boundary_info, dof_map))
            {
                fe_face->reinit(elem, side);
                for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                {
                    for (unsigned int i = 0; i < phi_face.size(); ++i)
                    {
                        for (unsigned int j = 0; j < phi_face.size(); ++j)
                        {
                            Ke(i, j) += PENALTY * phi_face[i][qp] * phi_face[j][qp] * JxW_face[qp];
                        }
                    }
                }
            }
        }
        dof_map.constrain_element_matrix(Ke, dof_indices);
        system.matrix->add_matrix(Ke, dof_indices);
    }
}

void
build_ib_ghosted_system_data(std::vector<SystemData>& ghosted_system_data,
                             std::vector<std::unique_ptr<NumericVector<double> > >& ib_ghost_system_vecs,
                             const std::vector<SystemData>& unghosted_system_data,
                             FEDataManager* const fe_data_manager)
{
    ghosted_system_data.clear();
    ib_ghost_system_vecs.clear();
    EquationSystems* equation_systems = fe_data_manager->getEquationSystems();
    for (const SystemData& system_data : unghosted_system_data)
    {
        const std::string& system_name = system_data.system_name;
        const System& system = equation_systems->get_system(system_name);
        NumericVector<double>* original_system_vec = system_data.system_vec;
        if (!original_system_vec) original_system_vec = system.current_local_solution.get();
        ib_ghost_system_vecs.emplace_back(fe_data_manager->buildIBGhostedVector(system_name));
        copy_and_synch(*original_system_vec, *ib_ghost_system_vecs.back(), /*close_v_in*/ false);
        ghosted_system_data.emplace_back(
            SystemData(system_name, system_data.vars, system_data.grad_vars, ib_ghost_system_vecs.back().get()));
    }
    return;
}
} // namespace

const std::string IBFEMethod::SOURCE_SYSTEM_NAME = "IB source system";

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEMethod::IBFEMethod(const std::string& object_name,
                       const Pointer<Database>& input_db,
                       MeshBase* mesh,
                       int max_levels,
                       bool register_for_restart,
                       const std::string& restart_read_dirname,
                       unsigned int restart_restore_number)
    : FEMechanicsBase(object_name, input_db, mesh, register_for_restart, restart_read_dirname, restart_restore_number)
{
    commonConstructor(input_db, max_levels);
    return;
} // IBFEMethod

IBFEMethod::IBFEMethod(const std::string& object_name,
                       const Pointer<Database>& input_db,
                       const std::vector<MeshBase*>& meshes,
                       int max_levels,
                       bool register_for_restart,
                       const std::string& restart_read_dirname,
                       unsigned int restart_restore_number)
    : FEMechanicsBase(object_name, input_db, meshes, register_for_restart, restart_read_dirname, restart_restore_number)
{
    commonConstructor(input_db, max_levels);
    return;
} // IBFEMethod

FEDataManager*
IBFEMethod::getFEDataManager(const unsigned int part) const
{
    TBOX_ASSERT(d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_meshes.size());
    return d_primary_fe_data_managers[part];
} // getFEDataManager

void
IBFEMethod::registerStaticPressurePart(PressureProjectionType projection_type,
                                       FEMechanicsBase::VolumetricEnergyDerivativeFcn U_prime,
                                       unsigned int part)
{
    TBOX_ASSERT(d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_meshes.size());
    // The same part can either have a pressure or can use stress normalization, but not both.
    if (d_has_stress_normalization_parts) TBOX_ASSERT(!d_stress_normalization_part[part]);
    FEMechanicsBase::registerStaticPressurePart(projection_type, U_prime, part);
    return;
} // registerStaticPressurePart

void
IBFEMethod::registerStressNormalizationPart(unsigned int part)
{
    TBOX_ASSERT(d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_meshes.size());
    // The same part can either have a pressure or can use stress normalization, but not both.
    if (d_has_static_pressure_parts) TBOX_ASSERT(!d_static_pressure_part[part]);
    if (d_stress_normalization_part[part]) return;
    d_has_stress_normalization_parts = true;
    d_stress_normalization_part[part] = true;
    auto& Phi_system = d_equation_systems[part]->add_system<LinearImplicitSystem>(PRESSURE_SYSTEM_NAME);
    d_equation_systems[part]->parameters.set<Real>("Phi_epsilon") = d_epsilon;
    Phi_system.attach_assemble_function(assemble_poisson);
    // This system has a single variable so we don't need to also specify diagonal coupling
    Phi_system.add_variable("Phi", d_fe_order_pressure[part], d_fe_family_pressure[part]);
    // Setup cached system vectors at restart.
    std::vector<std::string> vector_names;
    switch (d_ib_solver->getTimeSteppingType())
    {
    case MIDPOINT_RULE:
        vector_names = { "current", "half", "tmp", "RHS Vector" };
        break;
    default:
        vector_names = { "current", "half", "new", "tmp", "RHS Vector" };
    }
    IBTK::setup_system_vectors(d_equation_systems[part].get(),
                               { PRESSURE_SYSTEM_NAME },
                               vector_names,
                               RestartManager::getManager()->isFromRestart());
    return;
} // registerStressNormalizationPart

void
IBFEMethod::registerLagBodySourceFunction(const LagBodySourceFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_meshes.size());
    if (d_lag_body_source_part[part]) return;
    d_has_lag_body_source_parts = true;
    d_lag_body_source_part[part] = true;
    d_lag_body_source_fcn_data[part] = data;
    auto& Q_system = d_equation_systems[part]->add_system<ExplicitSystem>(SOURCE_SYSTEM_NAME);
    // This system has a single variable so we don't need to also specify diagonal coupling
    Q_system.add_variable("Q", d_fe_order_pressure[part], d_fe_family_pressure[part]);
    // Setup cached system vectors.
    IBTK::setup_system_vectors(d_equation_systems[part].get(),
                               { SOURCE_SYSTEM_NAME },
                               { "current", "half", "tmp", "RHS Vector" },
                               RestartManager::getManager()->isFromRestart());
    return;
} // registerLagBodySourceFunction

IBFEMethod::LagBodySourceFcnData
IBFEMethod::getLagBodySourceFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_meshes.size());
    return d_lag_body_source_fcn_data[part];
} // getLagBodySourceFunction

void
IBFEMethod::registerDirectForcingKinematics(const Pointer<IBFEDirectForcingKinematics>& data, unsigned int part)
{
    TBOX_ASSERT(part < d_meshes.size());
    d_direct_forcing_kinematics_data[part] = data;
    return;
} // registerDirectForcingKinematics

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
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        const int gcw = d_primary_fe_data_managers[part]->getGhostCellWidth().max();
        // We need to refine cells up to, but not including, those on
        // FEDataManager's finest patch level used for interaction
        const int tag_ln = d_primary_fe_data_managers[part]->getFinestPatchLevelNumber() - 1;
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

double
IBFEMethod::getMaxPointDisplacement() const
{
    double max_displacement = 0.0;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        PetscVector<double>& last = d_X_vecs->get("last_regrid", part);
        PetscVector<double>& diff = d_X_vecs->get("tmp", part);
        diff = d_X_vecs->get("solution", part);
        diff -= last;
        max_displacement = std::max(max_displacement, diff.linfty_norm());
    }
    // That's right: SAMRAI doesn't define BasePatchLevel::getPatch() making
    // this interface always require a cast
    return max_displacement / get_min_patch_dx(dynamic_cast<const PatchLevel<NDIM>&>(
                                  *d_hierarchy->getPatchLevel(getFinestPatchLevelNumber())));
} // getMaxPointDisplacement

void
IBFEMethod::inactivateLagrangianStructure(int structure_number, int /*level_number*/)
{
    // make sure we don't call before things are properly initialized
    TBOX_ASSERT(static_cast<std::size_t>(structure_number) < d_part_is_active.size());
    d_part_is_active[structure_number] = false;
    return;
} // inactivateLagrangianStructures

void
IBFEMethod::activateLagrangianStructure(int structure_number, int /*level_number*/)
{
    // make sure we don't call before things are properly initialized
    TBOX_ASSERT(static_cast<std::size_t>(structure_number) < d_part_is_active.size());
    d_part_is_active[structure_number] = true;
    return;
} // activateLagrangianStructures

bool
IBFEMethod::getLagrangianStructureIsActivated(int structure_number, int /*level_number*/) const
{
    // make sure we don't call before things are properly initialized
    TBOX_ASSERT(static_cast<std::size_t>(structure_number) < d_part_is_active.size());
    return d_part_is_active[structure_number];
} // getLagrangianStructureIsActivated

void
IBFEMethod::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    IBAMR_TIMER_START(t_preprocess_integrate_data);
    FEMechanicsBase::preprocessIntegrateData(current_time, new_time, num_cycles);

    d_started_time_integration = true;

    // Determine if we need to execute reinitElementMappings().
    bool do_reinit_element_mappings = false;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        const double patch_dx = get_min_patch_dx(dynamic_cast<const PatchLevel<NDIM>&>(
            *d_hierarchy->getPatchLevel(d_active_fe_data_managers[part]->getFinestPatchLevelNumber())));
        // To avoid allocating an extra vector we abuse "new" for a moment since
        // it gets overwritten below anyway
        PetscVector<double>& diff = d_X_vecs->get("new", part);
        PetscVector<double>& last = d_X_vecs->get("last_patch_elem_assoc", part);
        VecWAXPY(diff.vec(), -1.0, last.vec(), d_X_vecs->get("solution", part).vec());
        const double max_distance = diff.linfty_norm();

        if (max_distance > d_patch_association_cfl * patch_dx)
        {
            plog << d_object_name << "::preprocessIntegrateData(): "
                 << "Maximum structure node displacement is " << max_distance
                 << " - reinitializing element to patch mappings\n";
            do_reinit_element_mappings = true;
            break;
        }
    }
    if (do_reinit_element_mappings) reinitElementMappings();

    // Update direct forcing data.
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        if (d_direct_forcing_kinematics_data[part])
        {
            d_direct_forcing_kinematics_data[part]->preprocessIntegrateData(current_time, new_time, num_cycles);
        }
    }

    // Initialize variables.
    d_X_vecs->copy("solution", { "current", "new", "half" });
    d_U_vecs->copy("solution", { "current", "new", "half" });
    switch (d_ib_solver->getTimeSteppingType())
    {
    case MIDPOINT_RULE:
        d_F_vecs->copy("solution", { "current", "half" });
        if (d_P_vecs) d_P_vecs->copy("solution", { "current", "half" });
        break;
    default:
        d_F_vecs->copy("solution", { "current", "new", "half" });
        if (d_P_vecs) d_P_vecs->copy("solution", { "current", "new", "half" });
    }
    if (d_Q_vecs) d_Q_vecs->copy("solution", { "current", "half" });

    // Update the mask data.
    getVelocityHierarchyDataOps()->copyData(mask_new_idx, mask_current_idx);
    IBAMR_TIMER_STOP(t_preprocess_integrate_data);
    return;
} // preprocessIntegrateData

void
IBFEMethod::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    IBAMR_TIMER_START(t_postprocess_integrate_data);
    const std::string forcing_data_time_str = get_data_time_str(
        d_ib_solver->getTimeSteppingType() == MIDPOINT_RULE ? d_half_time : d_new_time, d_current_time, d_new_time);
    std::vector<std::vector<PetscVector<double>*> > vecs{ d_X_vecs->get("new"),
                                                          d_U_vecs->get("new"),
                                                          d_F_vecs->get(forcing_data_time_str) };
    if (d_P_vecs) vecs.push_back(d_P_vecs->get(forcing_data_time_str));
    if (d_Q_vecs) vecs.push_back(d_Q_vecs->get("half"));
    batch_vec_ghost_update(vecs, INSERT_VALUES, SCATTER_FORWARD);

    d_X_vecs->copy("new", { "solution", "current" });
    d_U_vecs->copy("new", { "solution", "current" });
    d_F_vecs->copy(forcing_data_time_str, { "solution", "current" });
    if (d_P_vecs) d_P_vecs->copy(forcing_data_time_str, { "solution", "current" });
    if (d_Q_vecs) d_Q_vecs->copy("half", { "solution", "current" });

    // Update direct forcing data.
    for (unsigned part = 0; part < d_meshes.size(); ++part)
    {
        if (d_direct_forcing_kinematics_data[part])
        {
            d_direct_forcing_kinematics_data[part]->postprocessIntegrateData(current_time, new_time, num_cycles);
        }
    }

    FEMechanicsBase::postprocessIntegrateData(current_time, new_time, num_cycles);
    IBAMR_TIMER_STOP(t_postprocess_integrate_data);
    return;
} // postprocessIntegrateData

void
IBFEMethod::interpolateVelocity(const int u_data_idx,
                                const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                const double data_time)
{
    IBAMR_TIMER_START(t_interpolate_velocity);
    const std::string data_time_str = get_data_time_str(data_time, d_current_time, d_new_time);

    // Ensure coarse grid data are consistent with any overlying fine grid data.
    for (int ln = getFinestPatchLevelNumber(); ln > getCoarsestPatchLevelNumber(); --ln)
    {
        if (ln < static_cast<int>(u_synch_scheds.size()) && u_synch_scheds[ln])
        {
            u_synch_scheds[ln]->coarsenData();
        }
    }

    if (d_use_scratch_hierarchy)
    {
        for (int ln = 0; ln <= getFinestPatchLevelNumber(); ++ln)
        {
            d_secondary_hierarchy
                ->transferPrimaryToSecondary(ln, u_data_idx, u_data_idx, data_time, d_ib_solver->getVelocityPhysBdryOp());
        }
    }
    else
    {
        // Communicate ghost data. Note that this only needs to be done when
        // we do *not* use the scratch hierarchy since the schedule will, in
        // that case, fill ghost data.
        for (const auto& u_ghost_fill_sched : u_ghost_fill_scheds)
        {
            if (u_ghost_fill_sched) u_ghost_fill_sched->fillData(data_time);
        }
    }

    std::vector<PetscVector<double>*> U_vecs = d_U_vecs->get(data_time_str);
    std::vector<PetscVector<double>*> X_vecs = d_X_vecs->get(data_time_str);
    std::vector<PetscVector<double>*> U_rhs_vecs = d_U_IB_vecs->getIBGhosted("tmp");
    std::vector<PetscVector<double>*> X_IB_ghost_vecs = d_X_IB_vecs->getIBGhosted("tmp");
    batch_vec_copy(X_vecs, X_IB_ghost_vecs);
    batch_vec_ghost_update(X_IB_ghost_vecs, INSERT_VALUES, SCATTER_FORWARD);

    // Build the right-hand-sides to compute the interpolated data.
    std::vector<Pointer<RefineSchedule<NDIM> > > no_fill(u_ghost_fill_scheds.size());
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        if (d_part_is_active[part])
        {
            d_active_fe_data_managers[part]->interpWeighted(u_data_idx,
                                                            *U_rhs_vecs[part],
                                                            *X_IB_ghost_vecs[part],
                                                            VELOCITY_SYSTEM_NAME,
                                                            no_fill,
                                                            data_time,
                                                            /*close_F*/ false,
                                                            /*close_X*/ false);
        }
    }

    // Note that FEDataManager only reads (and does not modify) Eulerian data
    // during interpolation so nothing needs to be transferred back to
    // d_hierarchy from d_secondary_hierarchy->getSecondaryHierarchy().
    batch_vec_ghost_update(U_rhs_vecs, ADD_VALUES, SCATTER_REVERSE);
    // We need ghost entries for constraint resolution
    batch_vec_ghost_update(U_rhs_vecs, INSERT_VALUES, SCATTER_FORWARD);

    // If the velocity system has constraints then those need to be dealt with
    // here and not during assembly because we do not distribute the
    // constraints corresponding to the IB partitioning.
    std::vector<PetscVector<double>*> vecs_for_second_summation;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        if (d_part_is_active[part])
        {
            EquationSystems& equation_systems = *d_primary_fe_data_managers[part]->getEquationSystems();
            vecs_for_second_summation.push_back(U_rhs_vecs[part]);
            apply_transposed_constraint_matrix(equation_systems.get_system(VELOCITY_SYSTEM_NAME).get_dof_map(),
                                               *U_rhs_vecs[part]);
        }
    }

    batch_vec_ghost_update(vecs_for_second_summation, ADD_VALUES, SCATTER_REVERSE);
    // we don't need ghost values on the RHS for the linear solve so don't scatter forward

    // Solve for the interpolated data.
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        if (d_part_is_active[part])
        {
            // TODO: Commenting out this line changes the solution slightly.
            d_U_vecs->get("solution", part) = *U_vecs[part];
            d_active_fe_data_managers[part]->computeL2Projection(d_U_vecs->get("solution", part),
                                                                 *U_rhs_vecs[part],
                                                                 VELOCITY_SYSTEM_NAME,
                                                                 d_interp_spec[part].use_consistent_mass_matrix,
                                                                 /*close_U*/ false,
                                                                 /*close_F*/ false);
            *U_vecs[part] = d_U_vecs->get("solution", part);
        }
        else
        {
            U_vecs[part]->zero();
        }
    }
    IBAMR_TIMER_STOP(t_interpolate_velocity);
    return;
} // interpolateVelocity

void
IBFEMethod::forwardEulerStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        int ierr = VecWAXPY(d_X_vecs->get("new", part).vec(),
                            dt,
                            d_U_vecs->get("current", part).vec(),
                            d_X_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(d_X_vecs->get("half", part).vec(),
                           0.5,
                           0.5,
                           0.0,
                           d_X_vecs->get("current", part).vec(),
                           d_X_vecs->get("new", part).vec());
        IBTK_CHKERRQ(ierr);
        if (d_direct_forcing_kinematics_data[part])
        {
            d_direct_forcing_kinematics_data[part]->forwardEulerStep(current_time,
                                                                     new_time,
                                                                     d_X_vecs->get("current", part),
                                                                     d_X_vecs->get("half", part),
                                                                     d_X_vecs->get("new", part));
        }
    }
    return;
} // forwardEulerStep

void
IBFEMethod::backwardEulerStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        int ierr = VecWAXPY(d_X_vecs->get("new", part).vec(),
                            dt,
                            d_U_vecs->get("new", part).vec(),
                            d_X_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(d_X_vecs->get("half", part).vec(),
                           0.5,
                           0.5,
                           0.0,
                           d_X_vecs->get("current", part).vec(),
                           d_X_vecs->get("new", part).vec());
        IBTK_CHKERRQ(ierr);
        if (d_direct_forcing_kinematics_data[part])
        {
            TBOX_ERROR("IBFEDirectForcingKinematics does not implement backwardEulerStep().\n");
        }
    }
    return;
} // backwardEulerStep

void
IBFEMethod::midpointStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        int ierr = VecWAXPY(d_X_vecs->get("new", part).vec(),
                            dt,
                            d_U_vecs->get("half", part).vec(),
                            d_X_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(d_X_vecs->get("half", part).vec(),
                           0.5,
                           0.5,
                           0.0,
                           d_X_vecs->get("current", part).vec(),
                           d_X_vecs->get("new", part).vec());
        IBTK_CHKERRQ(ierr);
        if (d_direct_forcing_kinematics_data[part])
        {
            d_direct_forcing_kinematics_data[part]->midpointStep(current_time,
                                                                 new_time,
                                                                 d_X_vecs->get("current", part),
                                                                 d_X_vecs->get("half", part),
                                                                 d_X_vecs->get("new", part));
        }
    }
    return;
} // midpointStep

void
IBFEMethod::trapezoidalStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        int ierr = VecWAXPY(d_X_vecs->get("new", part).vec(),
                            0.5 * dt,
                            d_U_vecs->get("current", part).vec(),
                            d_X_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPY(d_X_vecs->get("new", part).vec(), 0.5 * dt, d_U_vecs->get("new", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(d_X_vecs->get("half", part).vec(),
                           0.5,
                           0.5,
                           0.0,
                           d_X_vecs->get("current", part).vec(),
                           d_X_vecs->get("new", part).vec());
        IBTK_CHKERRQ(ierr);
        if (d_direct_forcing_kinematics_data[part])
        {
            d_direct_forcing_kinematics_data[part]->trapezoidalStep(current_time,
                                                                    new_time,
                                                                    d_X_vecs->get("current", part),
                                                                    d_X_vecs->get("half", part),
                                                                    d_X_vecs->get("new", part));
        }
    }
    return;
} // trapezoidalStep

void
IBFEMethod::computeLagrangianForce(const double data_time)
{
    IBAMR_TIMER_START(t_compute_lagrangian_force);
    const std::string data_time_str = get_data_time_str(data_time, d_current_time, d_new_time);
    batch_vec_ghost_update(d_X_vecs->get(data_time_str), INSERT_VALUES, SCATTER_FORWARD);
    d_F_vecs->zero("RHS Vector");
    d_F_vecs->zero("tmp");
    for (unsigned part = 0; part < d_meshes.size(); ++part)
    {
        if (d_stress_normalization_part[part])
        {
            computeStressNormalization(
                d_P_vecs->get(data_time_str, part), d_X_vecs->get(data_time_str, part), data_time, part);
        }
        if (d_static_pressure_part[part])
        {
            computeStaticPressure(
                d_P_vecs->get(data_time_str, part), d_X_vecs->get(data_time_str, part), data_time, part);
        }
        assembleInteriorForceDensityRHS(d_F_vecs->get("RHS Vector", part),
                                        d_X_vecs->get(data_time_str, part),
                                        d_P_vecs ? &d_P_vecs->get(data_time_str, part) : nullptr,
                                        data_time,
                                        part);
    }
    batch_vec_ghost_update(d_F_vecs->get("RHS Vector"), ADD_VALUES, SCATTER_REVERSE);
    // RHS vectors don't need ghost entries for the solve so we can skip the other scatter
    for (unsigned part = 0; part < d_meshes.size(); ++part)
    {
        d_active_fe_data_managers[part]->computeL2Projection(d_F_vecs->get("solution", part),
                                                             d_F_vecs->get("RHS Vector", part),
                                                             FORCE_SYSTEM_NAME,
                                                             d_use_consistent_mass_matrix,
                                                             /*close_U*/ false,
                                                             /*close_F*/ false);
    }
    d_F_vecs->copy("solution", { data_time_str });
    for (unsigned part = 0; part < d_meshes.size(); ++part)
    {
        if (d_direct_forcing_kinematics_data[part])
        {
            d_direct_forcing_kinematics_data[part]->computeLagrangianForce(d_F_vecs->get("tmp", part),
                                                                           d_X_vecs->get(data_time_str, part),
                                                                           d_U_vecs->get(data_time_str, part),
                                                                           data_time);
            int ierr = VecAXPY(d_F_vecs->get(data_time_str, part).vec(), 1.0, d_F_vecs->get("tmp", part).vec());
            IBTK_CHKERRQ(ierr);
        }
    }
    IBAMR_TIMER_STOP(t_compute_lagrangian_force);
    return;
} // computeLagrangianForce

void
IBFEMethod::spreadForce(const int f_data_idx,
                        RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                        const std::vector<Pointer<RefineSchedule<NDIM> > >& /*f_prolongation_scheds*/,
                        const double data_time)
{
    IBAMR_TIMER_START(t_spread_force);
    const std::string data_time_str = get_data_time_str(data_time, d_current_time, d_new_time);

    // Communicate ghost data.
    std::vector<PetscVector<double>*> X_IB_ghost_vecs = d_X_IB_vecs->getIBGhosted("tmp");
    std::vector<PetscVector<double>*> F_IB_ghost_vecs = d_F_IB_vecs->getIBGhosted("tmp");
    batch_vec_copy({ d_X_vecs->get(data_time_str), d_F_vecs->get(data_time_str) },
                   { X_IB_ghost_vecs, F_IB_ghost_vecs });
    batch_vec_ghost_update({ X_IB_ghost_vecs, F_IB_ghost_vecs }, INSERT_VALUES, SCATTER_FORWARD);

    // set up a new data index for computing forces on the active hierarchy.
    Pointer<PatchHierarchy<NDIM> > hierarchy =
        d_use_scratch_hierarchy ? d_secondary_hierarchy->getSecondaryHierarchy() : d_hierarchy;
    // get the data cache that is associated with whatever hierarchy we are
    // actually performing spreading operations on.
    std::shared_ptr<SAMRAIDataCache> data_cache =
        d_use_scratch_hierarchy ? d_secondary_hierarchy->getSAMRAIDataCache() : d_eulerian_data_cache;
    const int ln = hierarchy->getFinestLevelNumber();
    const auto f_scratch_data_idx = data_cache->getCachedPatchDataIndex(f_data_idx);
    const auto f_prolong_scratch_data_idx = data_cache->getCachedPatchDataIndex(f_data_idx);
    // zero data.
    Pointer<hier::Variable<NDIM> > f_var;
    VariableDatabase<NDIM>::getDatabase()->mapIndexToVariable(f_data_idx, f_var);
    auto f_active_data_ops = HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(f_var, hierarchy, true);
    f_active_data_ops->resetLevels(0, getFinestPatchLevelNumber());
    f_active_data_ops->setToScalar(f_scratch_data_idx,
                                   0.0,
                                   /*interior_only*/ false);
    f_active_data_ops->setToScalar(f_prolong_scratch_data_idx,
                                   0.0,
                                   /*interior_only*/ false);

    // Spread interior force density values.
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        if (!d_part_is_active[part]) continue;
        PetscVector<double>* X_ghost_vec = X_IB_ghost_vecs[part];
        PetscVector<double>* F_ghost_vec = F_IB_ghost_vecs[part];
        d_active_fe_data_managers[part]->spread(f_scratch_data_idx, *F_ghost_vec, *X_ghost_vec, FORCE_SYSTEM_NAME);
    }

    // Handle any transmission conditions.
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        if (!d_part_is_active[part]) continue;
        PetscVector<double>* X_ghost_vec = X_IB_ghost_vecs[part];
        PetscVector<double>* F_ghost_vec = F_IB_ghost_vecs[part];
        if (d_split_normal_force || d_split_tangential_force)
        {
            if (d_use_jump_conditions && d_split_normal_force)
            {
                assertStructureOnFinestLevel();
                imposeJumpConditions(f_scratch_data_idx, *F_ghost_vec, *X_ghost_vec, data_time, part);
            }
            if (!d_use_jump_conditions || d_split_tangential_force)
            {
                assertStructureOnFinestLevel();
                spreadTransmissionForceDensity(f_scratch_data_idx, *X_ghost_vec, data_time, part);
            }
        }
    }

    // Deal with force values spread outside the physical domain. Since these
    // are spread into ghost regions that don't correspond to actual degrees
    // of freedom they are ignored by the accumulation step - we have to
    // handle this before we do that.
    if (f_phys_bdry_op)
    {
        for (int ln = getCoarsestPatchLevelNumber(); ln <= getFinestPatchLevelNumber(); ++ln)
        {
            f_phys_bdry_op->setPatchDataIndex(f_scratch_data_idx);
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                const Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_scratch_data_idx);
                f_phys_bdry_op->accumulateFromPhysicalBoundaryData(*patch, data_time, f_data->getGhostCellWidth());
            }
        }
    }

    // Accumulate forces spread into patch ghost regions.
    {
        if (!d_ghost_data_accumulator)
        {
            // If we have multiple IBMethod objects we may end up with a wider
            // ghost region than the one required by this class. Hence, set the
            // ghost width by just picking whatever the data actually has at the
            // moment.
            const Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            const IntVector<NDIM> gcw =
                level->getPatchDescriptor()->getPatchDataFactory(f_scratch_data_idx)->getGhostCellWidth();

            // TODO - SAMRAIGhostDataAccumulator has not been tested with multiple levels
            d_ghost_data_accumulator.reset(new SAMRAIGhostDataAccumulator(
                hierarchy, f_var, gcw, getCoarsestPatchLevelNumber(), getFinestPatchLevelNumber()));
        }
        d_ghost_data_accumulator->accumulateGhostData(f_scratch_data_idx);
    }

    // Prolong forces spread onto coarser levels onto finer levels.
    for (int coarse_ln = 0; coarse_ln < getFinestPatchLevelNumber(); ++coarse_ln)
    {
        getProlongationSchedule(coarse_ln, f_scratch_data_idx, f_prolong_scratch_data_idx).fillData(data_time);

        // Add the values prolonged from the coarser level onto the scratch index on the finer level.
        auto f_data_ops = HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(f_var, hierarchy, true);
        f_data_ops->resetLevels(coarse_ln + 1, coarse_ln + 1);
        f_data_ops->add(f_scratch_data_idx, f_scratch_data_idx, f_prolong_scratch_data_idx);
    }

    if (d_use_scratch_hierarchy)
    {
        // harder case: we have a scratch index on the scratch hierarchy but
        // no corresponding index on the primary hierarchy.
        //
        // unlike the other data ops, this is always on the primary hierarchy
        auto f_primary_data_ops =
            HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(f_var, d_hierarchy, true);
        for (int ln = 0; ln <= getFinestPatchLevelNumber(); ++ln)
        {
            f_primary_data_ops->resetLevels(ln, ln);
            const auto f_primary_scratch_data_idx = d_eulerian_data_cache->getCachedPatchDataIndex(f_data_idx);
            // we have to zero everything here since the scratch to primary
            // communication does not touch ghost cells, which may have junk
            f_primary_data_ops->setToScalar(f_primary_scratch_data_idx,
                                            0.0,
                                            /*interior_only*/ false);
            d_secondary_hierarchy->transferSecondaryToPrimary(ln, f_primary_scratch_data_idx, f_scratch_data_idx, data_time);
            f_primary_data_ops->add(f_data_idx, f_data_idx, f_primary_scratch_data_idx);
        }
    }
    else
    {
        // easier case: f_scratch_data_idx is on the primary hierarchy so we
        // just need to add its values to those in f_data_idx
        f_active_data_ops->resetLevels(0, getFinestPatchLevelNumber());
        f_active_data_ops->add(f_data_idx, f_data_idx, f_scratch_data_idx);
    }
    IBAMR_TIMER_STOP(t_spread_force);
    return;
} // spreadForce

bool
IBFEMethod::hasFluidSources() const
{
    return d_has_lag_body_source_parts;
}

void
IBFEMethod::computeLagrangianFluidSource(double data_time)
{
    IBAMR_TIMER_START(t_compute_lagrangian_fluid_source);
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        if (!d_lag_body_source_part[part]) continue;

        EquationSystems& equation_systems = *d_primary_fe_data_managers[part]->getEquationSystems();
        const MeshBase& mesh = equation_systems.get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        // Extract the FE systems and DOF maps, and setup the FE object.
        auto& Q_system = equation_systems.get_system<ExplicitSystem>(SOURCE_SYSTEM_NAME);
        const DofMap& Q_dof_map = Q_system.get_dof_map();
        FEDataManager::SystemDofMapCache& Q_dof_map_cache =
            *d_primary_fe_data_managers[part]->getDofMapCache(SOURCE_SYSTEM_NAME);
        FEType Q_fe_type = Q_dof_map.variable_type(0);
        auto& X_system = equation_systems.get_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;

        FEDataInterpolation fe(dim, d_primary_fe_data_managers[part]->getFEData());
        std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, FIFTH);
        qrule->allow_rules_with_negative_weights = d_allow_rules_with_negative_weights;
        fe.attachQuadratureRule(qrule.get());
        fe.evalQuadraturePoints();
        fe.evalQuadratureWeights();
        fe.registerSystem(Q_system);
        NumericVector<double>& X_vec = d_X_vecs->get("half", part);
        const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, vars, vars, &X_vec);
        std::vector<size_t> Q_fcn_system_idxs;
        fe.setupInterpolatedSystemDataIndexes(
            Q_fcn_system_idxs, d_lag_body_source_fcn_data[part].system_data, &equation_systems);
        fe.init();

        const std::vector<libMesh::Point>& q_point = fe.getQuadraturePoints();
        const std::vector<double>& JxW = fe.getQuadratureWeights();
        const std::vector<std::vector<double> >& phi = fe.getPhi(Q_fe_type);

        const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
        const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
            fe.getGradVarInterpolation();

        std::vector<const std::vector<double>*> Q_var_data;
        std::vector<const std::vector<VectorValue<double> >*> Q_grad_var_data;

        // Setup global and elemental right-hand-side vectors.
        NumericVector<double>* Q_rhs_vec = Q_system.rhs;
        Q_rhs_vec->zero();
        DenseVector<double> Q_rhs_e;

        TensorValue<double> FF, FF_inv_trans;
        VectorValue<double> x;
        double Q;
        std::vector<libMesh::dof_id_type> dof_id_scratch;
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            const auto& Q_dof_indices = Q_dof_map_cache.dof_indices(elem);
            Q_rhs_e.resize(static_cast<int>(Q_dof_indices[0].size()));
            fe.reinit(elem);
            fe.collectDataForInterpolation(elem);
            fe.interpolate(elem);
            const unsigned int n_qp = qrule->n_points();
            const size_t n_basis = phi.size();
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const libMesh::Point& X = q_point[qp];
                const std::vector<double>& x_data = fe_interp_var_data[qp][X_sys_idx];
                const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
                get_x_and_FF(x, FF, x_data, grad_x_data);

                fe.setInterpolatedDataPointers(Q_var_data, Q_grad_var_data, Q_fcn_system_idxs, elem, qp);
                d_lag_body_source_fcn_data[part].fcn(
                    Q, FF, x, X, elem, Q_var_data, Q_grad_var_data, data_time, d_lag_body_source_fcn_data[part].ctx);
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    Q_rhs_e(k) += Q * phi[k][qp] * JxW[qp];
                }
            }

            // Apply constraints (e.g., enforce periodic boundary conditions)
            // and add the elemental contributions to the global vector.
            copy_dof_ids_to_vector(0, Q_dof_indices, dof_id_scratch);
            Q_dof_map.constrain_element_vector(Q_rhs_e, dof_id_scratch);
            Q_rhs_vec->add_vector(Q_rhs_e, dof_id_scratch);
        }

        // Solve for Q.
        NumericVector<double>& Q_vec = d_Q_vecs->get("half", part);
        d_primary_fe_data_managers[part]->computeL2Projection(
            Q_vec, *Q_rhs_vec, SOURCE_SYSTEM_NAME, d_use_consistent_mass_matrix);
    }
    IBAMR_TIMER_STOP(t_compute_lagrangian_fluid_source);
    return;
} // computeLagrangianFluidSource

void
IBFEMethod::spreadFluidSource(const int q_data_idx,
                              RobinPhysBdryPatchStrategy* q_phys_bdry_op,
                              const std::vector<Pointer<RefineSchedule<NDIM> > >& /*q_prolongation_scheds*/,
                              const double data_time)
{
    IBAMR_TIMER_START(t_spread_fluid_source);
    std::vector<PetscVector<double>*> X_IB_ghost_vecs = d_X_IB_vecs->getIBGhosted("tmp");
    std::vector<PetscVector<double>*> Q_IB_ghost_vecs = d_Q_IB_vecs->getIBGhosted("tmp");
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
    batch_vec_copy({ d_X_vecs->get("half"), d_Q_vecs->get("half") }, { X_IB_ghost_vecs, Q_IB_ghost_vecs });
    batch_vec_ghost_update({ X_IB_ghost_vecs, Q_IB_ghost_vecs }, INSERT_VALUES, SCATTER_FORWARD);

    if (d_use_scratch_hierarchy)
    {
        assertStructureOnFinestLevel();
        d_secondary_hierarchy
            ->transferPrimaryToSecondary(d_hierarchy->getFinestLevelNumber(), q_data_idx, q_data_idx, data_time);
    }

    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        if (!d_lag_body_source_part[part] || !d_part_is_active[part]) continue;

        d_active_fe_data_managers[part]->spread(q_data_idx,
                                                *Q_IB_ghost_vecs[part],
                                                *X_IB_ghost_vecs[part],
                                                SOURCE_SYSTEM_NAME,
                                                q_phys_bdry_op,
                                                data_time,
                                                /*close_Q*/ false,
                                                /*close_X*/ false);
    }

    if (d_use_scratch_hierarchy)
    {
        assertStructureOnFinestLevel();
        d_secondary_hierarchy
            ->transferSecondaryToPrimary(d_hierarchy->getFinestLevelNumber(), q_data_idx, q_data_idx, data_time);
    }

    IBAMR_TIMER_STOP(t_spread_fluid_source);
    return;
} // spreadFluidSource

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
    TBOX_ASSERT(part < d_meshes.size());
    d_interp_spec[part] = interp_spec;
    return;
}

void
IBFEMethod::setSpreadSpec(const FEDataManager::SpreadSpec& spread_spec, const unsigned int part)
{
    TBOX_ASSERT(!d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_meshes.size());
    d_spread_spec[part] = spread_spec;
    return;
}

void
IBFEMethod::setWorkloadSpec(const FEDataManager::WorkloadSpec& workload_spec, const unsigned int part)
{
    TBOX_ASSERT(!d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_meshes.size());
    d_workload_spec[part] = workload_spec;
    return;
}

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

    d_lagrangian_workload_var = new CellVariable<NDIM, double>(d_object_name + "::lagrangian_workload");
    registerVariable(d_lagrangian_workload_current_idx,
                     d_lagrangian_workload_new_idx,
                     d_lagrangian_workload_scratch_idx,
                     d_lagrangian_workload_var,
                     ghosts,
                     d_lagrangian_workload_coarsen_type,
                     d_lagrangian_workload_refine_type);
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

    // At this point we have not yet regridded, so we have not repartitioned:
    // make the scratch hierarchy a copy of the primary one so that it is not
    // null (which greatly simplifies control flow below)
    if (d_use_scratch_hierarchy)
    {
        d_secondary_hierarchy->reinit(getCoarsestPatchLevelNumber(), getFinestPatchLevelNumber(), d_hierarchy);
    }

    // Initialize the FE data managers.
    reinitElementMappings();

    // The patch hierarchies are finally available so set up the scratch patch data cache:
    d_eulerian_data_cache->setPatchHierarchy(hierarchy);
    d_eulerian_data_cache->resetLevels(0, getFinestPatchLevelNumber());

    d_is_initialized = true;
    return;
} // initializePatchHierarchy

void
IBFEMethod::registerLoadBalancer(Pointer<LoadBalancer<NDIM> > load_balancer, int workload_data_idx)
{
    IBAMR_DEPRECATED_MEMBER_FUNCTION1("IBFEMethod", "registerLoadBalancer");
    TBOX_ASSERT(load_balancer);
    d_load_balancer = load_balancer;
    d_workload_idx = workload_data_idx;
    return;
} // registerLoadBalancer

void
IBFEMethod::addWorkloadEstimate(Pointer<PatchHierarchy<NDIM> > hierarchy, const int workload_data_idx)
{
    IBAMR_TIMER_START(t_add_workload_estimate);
    const bool old_d_do_log = d_do_log;
    if (d_skip_initial_workload_log && !d_started_time_integration)
    {
        d_do_log = false;
        for (unsigned int part = 0; part < d_meshes.size(); ++part)
        {
            if (d_use_scratch_hierarchy) d_scratch_fe_data_managers[part]->setLoggingEnabled(false);
            d_primary_fe_data_managers[part]->setLoggingEnabled(false);
        }
    }

    // TODO: since this function is called before we finish setting everything
    // up the data interdependencies are complicated and its too hard to
    // communicate between the scratch and primary hierarchies. Try to fix
    // this.
    if (d_use_scratch_hierarchy && hierarchy == d_secondary_hierarchy->getSecondaryHierarchy())
    {
        for (unsigned int part = 0; part < d_meshes.size(); ++part)
        {
            d_scratch_fe_data_managers[part]->addWorkloadEstimate(hierarchy, workload_data_idx);
        }
    }
    else
    {
        // If we use the scratch hierarchy then the work on the provided
        // hierarchy is zero, so we have nothing to add
        if (!d_use_scratch_hierarchy)
        {
            for (unsigned int part = 0; part < d_meshes.size(); ++part)
            {
                d_primary_fe_data_managers[part]->addWorkloadEstimate(hierarchy, workload_data_idx);
            }
        }
    }

    if (d_do_log)
    {
        const int n_processes = IBTK_MPI::getNodes();
        const int current_rank = IBTK_MPI::getRank();

        std::vector<double> workload_per_processor(n_processes);
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, 0, hierarchy->getFinestLevelNumber());
        workload_per_processor[current_rank] = hier_cc_data_ops.L1Norm(workload_data_idx, IBTK::invalid_index, true);

        const auto right_padding = std::size_t(std::log10(n_processes)) + 1;

        int ierr = MPI_Allreduce(MPI_IN_PLACE,
                                 workload_per_processor.data(),
                                 workload_per_processor.size(),
                                 MPI_DOUBLE,
                                 MPI_SUM,
                                 IBTK_MPI::getCommunicator());
        TBOX_ASSERT(ierr == 0);
        if (current_rank == 0)
        {
            for (int rank = 0; rank < n_processes; ++rank)
            {
                SAMRAI::tbox::plog << "workload estimate on processor " << std::setw(right_padding) << std::left << rank
                                   << " = " << long(workload_per_processor[rank]) << '\n';
            }
        }

        std::vector<std::size_t> dofs_per_processor(n_processes);
        for (unsigned int part = 0; part < d_meshes.size(); ++part)
        {
            auto& equation_systems = *d_active_fe_data_managers[part]->getEquationSystems();
            for (unsigned int system_n = 0; system_n < equation_systems.n_systems(); ++system_n)
            {
                dofs_per_processor[current_rank] += equation_systems.get_system(system_n).n_local_dofs();
            }
        }

        ierr = MPI_Allreduce(MPI_IN_PLACE,
                             dofs_per_processor.data(),
                             dofs_per_processor.size(),
                             MPI_UNSIGNED_LONG,
                             MPI_SUM,
                             IBTK_MPI::getCommunicator());
        TBOX_ASSERT(ierr == 0);
        if (current_rank == 0)
        {
            for (int rank = 0; rank < n_processes; ++rank)
            {
                SAMRAI::tbox::plog << "local active DoFs on processor " << std::setw(right_padding) << std::left << rank
                                   << " = " << dofs_per_processor[rank] << '\n';
            }
        }
    }

    d_do_log = old_d_do_log;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        if (d_use_scratch_hierarchy) d_scratch_fe_data_managers[part]->setLoggingEnabled(d_do_log);
        d_primary_fe_data_managers[part]->setLoggingEnabled(d_do_log);
    }
    IBAMR_TIMER_STOP(t_add_workload_estimate);
    return;
} // addWorkloadEstimate

void IBFEMethod::beginDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                         Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    IBAMR_TIMER_START(t_begin_data_redistribution);
    // clear some things that contain data specific to the current patch hierarchy
    d_ghost_data_accumulator.reset();
    d_prolongation_schedules.clear();

    // If we are using the secondary hierarchy then communicate workload data
    // back to the primary hierarchy (which we will use in endDataRedistribution
    // to load balance the secondary hierarchy)
    if (d_is_initialized && d_use_scratch_hierarchy)
    {
        {
            for (int ln = 0; ln <= getFinestPatchLevelNumber(); ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = d_secondary_hierarchy->getSecondaryHierarchy()->getPatchLevel(ln);
                if (!level->checkAllocated(d_lagrangian_workload_current_idx))
                    level->allocatePatchData(d_lagrangian_workload_current_idx);
            }

            HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(
                d_secondary_hierarchy->getSecondaryHierarchy(),
                0,
                d_secondary_hierarchy->getSecondaryHierarchy()->getFinestLevelNumber());
            hier_cc_data_ops.setToScalar(d_lagrangian_workload_current_idx, 0.0);
        }
        addWorkloadEstimate(d_secondary_hierarchy->getSecondaryHierarchy(), d_lagrangian_workload_current_idx);

        {
            HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(
                d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
            hier_cc_data_ops.setToScalar(d_lagrangian_workload_current_idx, 0.0);
        }
        for (int ln = 0; ln <= getFinestPatchLevelNumber(); ++ln)
        {
            d_secondary_hierarchy
                ->transferSecondaryToPrimary(
                    ln, d_lagrangian_workload_current_idx, d_lagrangian_workload_current_idx, 0.0);
        }
    }

    IBAMR_TIMER_STOP(t_begin_data_redistribution);
    return;
} // beginDataRedistribution

void IBFEMethod::endDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                       Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    IBAMR_TIMER_START(t_end_data_redistribution);
    // if we are not initialized then there is nothing to do
    if (d_is_initialized)
    {
        if (d_use_scratch_hierarchy)
        {
            if (d_do_log) plog << "IBFEMethod: starting scratch hierarchy regrid" << std::endl;
            d_secondary_hierarchy->reinit(getCoarsestPatchLevelNumber(),
                                          getFinestPatchLevelNumber(),
                                          d_hierarchy,
                                          d_lagrangian_workload_current_idx);
            if (d_do_log) plog << "IBFEMethod: finished scratch hierarchy regrid" << std::endl;
        }

        // At this point in the code SAMRAI has already redistributed the
        // patches (usually by taking into account the number of IB points on
        // each patch). Here is the other half: if requested, we inform
        // libMesh of the updated partitioning so that libMesh Elems and Nodes
        // are on the same processor as the relevant SAMRAI patch.
        if (d_libmesh_partitioner_type == SAMRAI_BOX)
        {
            for (unsigned int part = 0; part < d_meshes.size(); ++part)
            {
                EquationSystems& equation_systems = *d_active_fe_data_managers[part]->getEquationSystems();
                MeshBase& mesh = equation_systems.get_mesh();
                BoxPartitioner partitioner(*d_hierarchy, equation_systems.get_system(COORDS_SYSTEM_NAME));
                partitioner.repartition(mesh);
            }
        }

        // We only need to reinitialize FE data when AMR is enabled (which is
        // not yet implemented)
        if (d_libmesh_use_amr) reinitializeFEData();

        reinitElementMappings();

        // Reinitialize System vector objects that depend on the IB partitioning.
        d_X_IB_vecs->reinit();
        d_U_IB_vecs->reinit();
        d_F_IB_vecs->reinit();
        if (d_Q_IB_vecs) d_Q_IB_vecs->reinit();

        // This condition is complex - we only want to log the workload if we are
        // - logging
        // - using the scratch hierarchy
        //
        // but not before the first timestep should d_skip_initial_workload_log
        // be true.
        if (d_use_scratch_hierarchy && d_do_log &&
            (d_started_time_integration || (!d_started_time_integration && !d_skip_initial_workload_log)))
        {
            plog << "IBFEMethod::scratch hierarchy workload" << std::endl;
            const int n_processes = IBTK_MPI::getNodes();
            const int current_rank = IBTK_MPI::getRank();
            const auto right_padding = std::size_t(std::log10(n_processes)) + 1;
            std::vector<double> workload_per_processor(n_processes);
            HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(
                d_secondary_hierarchy->getSecondaryHierarchy(),
                0,
                d_secondary_hierarchy->getSecondaryHierarchy()->getFinestLevelNumber());
            workload_per_processor[current_rank] =
                hier_cc_data_ops.L1Norm(d_lagrangian_workload_current_idx, IBTK::invalid_index, true);

            int ierr = MPI_Allreduce(MPI_IN_PLACE,
                                     workload_per_processor.data(),
                                     workload_per_processor.size(),
                                     MPI_DOUBLE,
                                     MPI_SUM,
                                     IBTK_MPI::getCommunicator());
            TBOX_ASSERT(ierr == 0);
            if (current_rank == 0)
            {
                for (int rank = 0; rank < n_processes; ++rank)
                {
                    plog << "workload estimate on processor " << std::setw(right_padding) << std::left << rank << " = "
                         << long(workload_per_processor[rank]) << '\n';
                }
            }

            plog << "IBFEMethod:: end scratch hierarchy workload" << std::endl;
        }
    }

    // Store the coordinates at which we performed the regrid for later use in
    // calculating structure displacement.
    d_X_vecs->copy("solution", { "last_regrid" });

    IBAMR_TIMER_STOP(t_end_data_redistribution);
    return;
} // endDataRedistribution

void
IBFEMethod::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
                                int /*level_number*/,
                                double /*init_data_time*/,
                                bool /*can_be_refined*/,
                                bool /*initial_time*/,
                                Pointer<BasePatchLevel<NDIM> > /*old_level*/,
                                bool /*allocate_data*/)
{
} // initializeLevelData

void
IBFEMethod::resetHierarchyConfiguration(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                        int /*coarsest_level*/,
                                        int /*finest_level*/)
{
    if (d_use_scratch_hierarchy && hierarchy == d_secondary_hierarchy->getSecondaryHierarchy()) return;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        d_primary_fe_data_managers[part]->setPatchHierarchy(hierarchy);
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
    IBAMR_TIMER_START(t_apply_gradient_detector);
    const bool old_d_do_log = d_do_log;
    if (d_skip_initial_workload_log && !d_started_time_integration)
    {
        d_do_log = false;
        for (unsigned int part = 0; part < d_meshes.size(); ++part)
        {
            if (d_use_scratch_hierarchy) d_scratch_fe_data_managers[part]->setLoggingEnabled(false);
            d_primary_fe_data_managers[part]->setLoggingEnabled(false);
        }
    }

    Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));

    // If we are working on the scratch hierarchy then this is the first and
    // last place where we will do gradient detection: hence it is our
    // responsibility to zero the array; otherwise SAMRAI tries to refine
    // everything.
    if (d_use_scratch_hierarchy && hierarchy == d_secondary_hierarchy->getSecondaryHierarchy())
    {
        HierarchyCellDataOpsInteger<NDIM> hier_cc_data_ops(hierarchy, level_number, level_number);
        hier_cc_data_ops.setToScalar(tag_index, 0);
    }

    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        // We still apply the gradient detector (i.e., tag cells for
        // refinement that contain Lagrangian points) even if a part is
        // disabled so that it can be reenabled without needing to regrid
        auto* fe_data_manager = hierarchy == d_primary_fe_data_managers[part]->getPatchHierarchy() ?
                                    d_primary_fe_data_managers[part] :
                                    d_scratch_fe_data_managers[part];
        fe_data_manager->applyGradientDetector(
            hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    }

    d_do_log = old_d_do_log;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        if (d_use_scratch_hierarchy) d_scratch_fe_data_managers[part]->setLoggingEnabled(d_do_log);
        d_primary_fe_data_managers[part]->setLoggingEnabled(d_do_log);
    }
    IBAMR_TIMER_STOP(t_apply_gradient_detector);
    return;
} // applyGradientDetector

void
IBFEMethod::putToDatabase(Pointer<Database> db)
{
    FEMechanicsBase::putToDatabase(db);
    db->putInteger("IBFE_METHOD_VERSION", IBFE_METHOD_VERSION);
    db->putIntegerArray("d_ghosts", d_ghosts, NDIM);
    db->putBool("d_split_normal_force", d_split_normal_force);
    db->putBool("d_split_tangential_force", d_split_tangential_force);
    db->putBool("d_use_jump_conditions", d_use_jump_conditions);
    std::unique_ptr<bool[]> part_is_active_arr{ new bool[d_part_is_active.size()] };
    std::copy(d_part_is_active.begin(), d_part_is_active.end(), part_is_active_arr.get());
    db->putInteger("part_is_active_arr_size", d_part_is_active.size());
    db->putBoolArray("part_is_active_arr", part_is_active_arr.get(), d_part_is_active.size());
    return;
} // putToDatabase

SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> >
IBFEMethod::getScratchHierarchy()
{
    return d_secondary_hierarchy->getSecondaryHierarchy();
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBFEMethod::doInitializeFEEquationSystems()
{
    FEMechanicsBase::doInitializeFEEquationSystems();

    const bool from_restart = RestartManager::getManager()->isFromRestart();

    // Setup system vectors.
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        EquationSystems& equation_systems = *d_equation_systems[part];
        IBTK::setup_system_vectors(
            &equation_systems, { VELOCITY_SYSTEM_NAME }, { "current", "half", "new" }, from_restart);
        IBTK::setup_system_vectors(&equation_systems,
                                   { COORDS_SYSTEM_NAME },
                                   { "last_patch_elem_assoc", "last_regrid", "current", "half", "new", "tmp" },
                                   from_restart);
        std::vector<std::string> F_vector_names;
        switch (d_ib_solver->getTimeSteppingType())
        {
        case MIDPOINT_RULE:
            F_vector_names = { "current", "half", "tmp", "RHS Vector" };
            break;
        default:
            F_vector_names = { "current", "half", "new", "tmp", "RHS Vector" };
        }
        IBTK::setup_system_vectors(&equation_systems, { FORCE_SYSTEM_NAME }, F_vector_names, from_restart);
    }

    // Create the FE data managers that manage mappings between the FE mesh
    // parts and the Cartesian grid.
    d_primary_fe_data_managers.resize(d_meshes.size());
    d_scratch_fe_data_managers.resize(d_meshes.size());
    d_active_fe_data_managers.resize(d_meshes.size());
    IntVector<NDIM> min_ghost_width(0);
    if (!d_eulerian_data_cache) d_eulerian_data_cache = std::make_shared<SAMRAIDataCache>();
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        const std::string manager_name = "IBFEMethod FEDataManager::" + std::to_string(part);
        d_primary_fe_data_managers[part] = FEDataManager::getManager(d_fe_data[part],
                                                                     manager_name,
                                                                     d_fe_data_manager_db,
                                                                     d_max_level_number + 1,
                                                                     d_interp_spec[part],
                                                                     d_spread_spec[part],
                                                                     d_workload_spec[part],
                                                                     min_ghost_width,
                                                                     d_eulerian_data_cache);
        if (d_use_scratch_hierarchy)
        {
            d_scratch_fe_data_managers[part] = FEDataManager::getManager(d_fe_data[part],
                                                                         manager_name + "::scratch",
                                                                         d_fe_data_manager_db,
                                                                         d_max_level_number + 1,
                                                                         d_interp_spec[part],
                                                                         d_spread_spec[part],
                                                                         d_workload_spec[part],
                                                                         min_ghost_width,
                                                                         d_secondary_hierarchy->getSAMRAIDataCache());
            d_active_fe_data_managers[part] = d_scratch_fe_data_managers[part];
        }
        else
        {
            d_active_fe_data_managers[part] = d_primary_fe_data_managers[part];
        }

        d_active_fe_data_managers[part]->setLoggingEnabled(d_do_log);
        d_ghosts = IntVector<NDIM>::max(d_ghosts, d_active_fe_data_managers[part]->getGhostCellWidth());

        // Since the scratch and primary FEDataManagers use the same FEData
        // object we only have to do this assignment once
        d_active_fe_data_managers[part]->COORDINATES_SYSTEM_NAME = COORDS_SYSTEM_NAME;
    }

    return;
} // doInitializeFEEquationSystems

void
IBFEMethod::doInitializeFESystemVectors()
{
    FEMechanicsBase::doInitializeFESystemVectors();

    // The choice of FEDataManager set is important here since it determines the
    // IB ghost regions.
    d_X_vecs.reset(new LibMeshSystemIBVectors(d_active_fe_data_managers, COORDS_SYSTEM_NAME));
    d_U_vecs.reset(new LibMeshSystemIBVectors(d_active_fe_data_managers, VELOCITY_SYSTEM_NAME));
    d_F_vecs.reset(new LibMeshSystemIBVectors(d_active_fe_data_managers, FORCE_SYSTEM_NAME));
    d_X_IB_vecs.reset(new LibMeshSystemIBVectors(d_active_fe_data_managers, COORDS_SYSTEM_NAME));
    d_U_IB_vecs.reset(new LibMeshSystemIBVectors(d_active_fe_data_managers, VELOCITY_SYSTEM_NAME));
    d_F_IB_vecs.reset(new LibMeshSystemIBVectors(d_active_fe_data_managers, FORCE_SYSTEM_NAME));
    if (d_has_stress_normalization_parts)
    {
        d_P_vecs.reset(
            new LibMeshSystemIBVectors(d_active_fe_data_managers, d_stress_normalization_part, PRESSURE_SYSTEM_NAME));
    }
    if (d_has_static_pressure_parts)
    {
        d_P_vecs.reset(
            new LibMeshSystemIBVectors(d_active_fe_data_managers, d_static_pressure_part, PRESSURE_SYSTEM_NAME));
    }
    if (d_has_lag_body_source_parts)
    {
        d_Q_vecs.reset(
            new LibMeshSystemIBVectors(d_active_fe_data_managers, d_lag_body_source_part, SOURCE_SYSTEM_NAME));
        d_Q_IB_vecs.reset(
            new LibMeshSystemIBVectors(d_active_fe_data_managers, d_lag_body_source_part, SOURCE_SYSTEM_NAME));
    }
    return;
} // doInitializeFESystemVectors

void
IBFEMethod::doInitializeFEData(const bool use_present_data)
{
    FEMechanicsBase::doInitializeFEData(use_present_data);

    // Initialize FE equation systems.
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        EquationSystems& equation_systems = *d_equation_systems[part];
        if (d_stress_normalization_part[part])
        {
            auto& Phi_system = equation_systems.get_system<LinearImplicitSystem>(PRESSURE_SYSTEM_NAME);
            Phi_system.assemble_before_solve = false;
            Phi_system.assemble();
        }
        if (d_direct_forcing_kinematics_data[part])
        {
            d_direct_forcing_kinematics_data[part]->initializeKinematicsData(!use_present_data);
        }
    }
    return;
} // doInitializeFEData

void
IBFEMethod::computeStressNormalization(PetscVector<double>& Phi_vec,
                                       PetscVector<double>& X_vec,
                                       const double data_time,
                                       const unsigned int part)
{
    // Extract the mesh.
    EquationSystems& equation_systems = *d_primary_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems.get_mesh();
    const BoundaryInfo& boundary_info = mesh.get_boundary_info();
    const unsigned int dim = mesh.mesh_dimension();

    // Setup extra data needed to compute stresses/forces.

    // Extract the FE systems and DOF maps, and setup the FE objects.
    auto& Phi_system = equation_systems.get_system<LinearImplicitSystem>(PRESSURE_SYSTEM_NAME);
    const DofMap& Phi_dof_map = Phi_system.get_dof_map();
    FEDataManager::SystemDofMapCache& Phi_dof_map_cache =
        *d_primary_fe_data_managers[part]->getDofMapCache(PRESSURE_SYSTEM_NAME);
    FEType Phi_fe_type = Phi_dof_map.variable_type(0);
    std::vector<int> Phi_vars(1, 0);

    auto& X_system = equation_systems.get_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
    std::vector<int> X_vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) X_vars[d] = d;

    FEDataInterpolation fe(dim, d_primary_fe_data_managers[part]->getFEData());
    std::unique_ptr<QBase> qrule_face = QBase::build(QGAUSS, dim - 1, FIFTH);
    qrule_face->allow_rules_with_negative_weights = d_allow_rules_with_negative_weights;
    fe.attachQuadratureRuleFace(qrule_face.get());
    fe.evalNormalsFace();
    fe.evalQuadraturePointsFace();
    fe.evalQuadratureWeightsFace();
    fe.registerSystem(Phi_system, Phi_vars,
                      Phi_vars); // compute phi and dphi for the Phi system
    const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, X_vars, X_vars, &X_vec);
    const size_t num_PK1_fcns = d_PK1_stress_fcn_data[part].size();
    std::vector<std::vector<size_t> > PK1_fcn_system_idxs(num_PK1_fcns);
    for (unsigned int k = 0; k < num_PK1_fcns; ++k)
    {
        fe.setupInterpolatedSystemDataIndexes(
            PK1_fcn_system_idxs[k], d_PK1_stress_fcn_data[part][k].system_data, &equation_systems);
    }
    std::vector<size_t> surface_force_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        surface_force_fcn_system_idxs, d_lag_surface_force_fcn_data[part].system_data, &equation_systems);
    std::vector<size_t> surface_pressure_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        surface_pressure_fcn_system_idxs, d_lag_surface_pressure_fcn_data[part].system_data, &equation_systems);
    fe.init();

    const std::vector<libMesh::Point>& q_point_face = fe.getQuadraturePointsFace();
    const std::vector<double>& JxW_face = fe.getQuadratureWeightsFace();
    const std::vector<libMesh::Point>& normal_face = fe.getNormalsFace();
    const std::vector<std::vector<double> >& phi_face = fe.getPhiFace(Phi_fe_type);

    const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
    const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
        fe.getGradVarInterpolation();

    std::vector<std::vector<const std::vector<double>*> > PK1_var_data(num_PK1_fcns);
    std::vector<std::vector<const std::vector<VectorValue<double> >*> > PK1_grad_var_data(num_PK1_fcns);
    std::vector<const std::vector<double>*> surface_force_var_data, surface_pressure_var_data;
    std::vector<const std::vector<VectorValue<double> >*> surface_force_grad_var_data, surface_pressure_grad_var_data;

    // Setup global and elemental right-hand-side vectors.
    NumericVector<double>* Phi_rhs_vec = Phi_system.rhs;
    Phi_rhs_vec->zero();
    DenseVector<double> Phi_rhs_e;

    // Set up boundary conditions for Phi.
    TensorValue<double> PP, FF, FF_trans, FF_inv_trans;
    VectorValue<double> F, F_s, F_qp, n, x;
    std::vector<libMesh::dof_id_type> dof_id_scratch;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        bool reinit_all_data = true;
        for (unsigned int side = 0; side < elem->n_sides(); ++side)
        {
            // Skip non-physical boundaries.
            if (!is_physical_bdry(elem, side, boundary_info, Phi_dof_map)) continue;

            // Determine if we need to integrate surface forces along this
            // part of the physical boundary; if not, skip the present side.
            const bool at_dirichlet_bdry = is_dirichlet_bdry(elem, side, boundary_info, Phi_dof_map);
            if (at_dirichlet_bdry) continue;

            fe.reinit(elem, side);
            const auto& Phi_dof_indices = Phi_dof_map_cache.dof_indices(elem);
            if (reinit_all_data)
            {
                Phi_rhs_e.resize(static_cast<int>(Phi_dof_indices[0].size()));
                fe.collectDataForInterpolation(elem);
                reinit_all_data = false;
            }
            fe.interpolate(elem, side);
            const unsigned int n_qp = qrule_face->n_points();
            const size_t n_basis = phi_face.size();
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                // X:     reference coordinate
                // x:     current coordinate
                // FF:    deformation gradient associated with x = chi(X,t)
                //        (FF = dchi/dX)
                // J:     Jacobian determinant (J = det(FF))
                // N:     unit normal in the reference configuration
                // n:     unit normal in the current configuration
                // dA_da: reference surface area per current surface area (from Nanson's
                //        relation)
                const libMesh::Point& X = q_point_face[qp];
                const std::vector<double>& x_data = fe_interp_var_data[qp][X_sys_idx];
                const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
                get_x_and_FF(x, FF, x_data, grad_x_data);
                const double J = std::abs(FF.det());
                FF_trans = FF.transpose();
                tensor_inverse_transpose(FF_inv_trans, FF, NDIM);
                const libMesh::VectorValue<double>& N = normal_face[qp];
                n = (FF_inv_trans * N).unit();
                const double dA_da = 1.0 / (J * (FF_inv_trans * N) * n);

                // Here we build up the boundary value for Phi.
                double Phi = 0.0;
                for (unsigned int k = 0; k < num_PK1_fcns; ++k)
                {
                    if (d_PK1_stress_fcn_data[part][k].fcn)
                    {
                        // Compute the value of the first Piola-Kirchhoff stress
                        // tensor at the quadrature point and add the corresponding
                        // traction force to the right-hand-side vector.
                        fe.setInterpolatedDataPointers(
                            PK1_var_data[k], PK1_grad_var_data[k], PK1_fcn_system_idxs[k], elem, qp);
                        d_PK1_stress_fcn_data[part][k].fcn(PP,
                                                           FF,
                                                           x,
                                                           X,
                                                           elem,
                                                           PK1_var_data[k],
                                                           PK1_grad_var_data[k],
                                                           data_time,
                                                           d_PK1_stress_fcn_data[part][k].ctx);
                        Phi += n * ((PP * FF_trans) * n) / J;
                    }
                }

                if (d_lag_surface_force_fcn_data[part].fcn)
                {
                    // Compute the value of the surface force at the
                    // quadrature point and add the corresponding force to
                    // the right-hand-side vector.
                    fe.setInterpolatedDataPointers(
                        surface_force_var_data, surface_force_grad_var_data, surface_force_fcn_system_idxs, elem, qp);
                    d_lag_surface_force_fcn_data[part].fcn(F_s,
                                                           n,
                                                           N,
                                                           FF,
                                                           x,
                                                           X,
                                                           elem,
                                                           side,
                                                           surface_force_var_data,
                                                           surface_force_grad_var_data,
                                                           data_time,
                                                           d_lag_surface_force_fcn_data[part].ctx);
                    Phi -= n * F_s * dA_da;
                }

                if (d_lag_surface_pressure_fcn_data[part].fcn)
                {
                    // Compute the value of the pressure at the quadrature
                    // point and add the corresponding force to the
                    // right-hand-side vector.
                    double P = 0.0;
                    fe.setInterpolatedDataPointers(surface_pressure_var_data,
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
                                                              side,
                                                              surface_pressure_var_data,
                                                              surface_pressure_grad_var_data,
                                                              data_time,
                                                              d_lag_surface_pressure_fcn_data[part].ctx);
                    Phi += P;
                }

                // Add the boundary forces to the right-hand-side vector.
                for (unsigned int i = 0; i < n_basis; ++i)
                {
                    Phi_rhs_e(i) += PENALTY * Phi * phi_face[i][qp] * JxW_face[qp];
                }
            }

            // Apply constraints (e.g., enforce periodic boundary conditions)
            // and add the elemental contributions to the global vector.
            copy_dof_ids_to_vector(0, Phi_dof_indices, dof_id_scratch);
            Phi_dof_map.constrain_element_vector(Phi_rhs_e, dof_id_scratch);
            Phi_rhs_vec->add_vector(Phi_rhs_e, dof_id_scratch);
        }
    }

    // Solve for Phi.
    Phi_rhs_vec->close();
    Phi_system.solve();
    Phi_dof_map.enforce_constraints_exactly(Phi_system, Phi_system.solution.get());
    Phi_vec = *Phi_system.solution;
    Phi_vec.close();
    return;
}

void
IBFEMethod::spreadTransmissionForceDensity(const int f_data_idx,
                                           PetscVector<double>& X_ghost_vec,
                                           const double data_time,
                                           const unsigned int part)
{
    if (!d_split_normal_force && !d_split_tangential_force) return;

    // Check to see if we need to integrate the surface forces.
    const bool integrate_normal_force =
        d_split_normal_force && !d_use_jump_conditions && !d_stress_normalization_part[part];
    const bool integrate_tangential_force = d_split_tangential_force;
    if (!integrate_normal_force && !integrate_tangential_force) return;

    assertStructureOnFinestLevel();
    const int level_num = d_primary_fe_data_managers[part]->getFinestPatchLevelNumber();
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);

    // Extract the mesh.
    EquationSystems& equation_systems = *d_primary_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems.get_mesh();
    const BoundaryInfo& boundary_info = mesh.get_boundary_info();
    const unsigned int dim = mesh.mesh_dimension();

    // Extract the FE systems and DOF maps, and setup the FE object.
    auto& G_system = equation_systems.get_system<ExplicitSystem>(FORCE_SYSTEM_NAME);
    const DofMap& G_dof_map = G_system.get_dof_map();
    FEType G_fe_type = G_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(G_dof_map.variable_type(d) == G_fe_type);
    }
    std::vector<std::vector<unsigned int> > G_dof_indices(NDIM);
    auto& X_system = equation_systems.get_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
    std::vector<int> vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;

    FEDataInterpolation fe(dim, d_primary_fe_data_managers[part]->getFEData());
    QuadratureType default_quad_type = (d_spread_spec[part].quad_type != INVALID_Q_RULE) ?
                                           d_spread_spec[part].quad_type :
                                           d_default_quad_type_stress[part];
    Order default_quad_order = (d_spread_spec[part].quad_order != INVALID_ORDER) ? d_spread_spec[part].quad_order :
                                                                                   d_default_quad_order_stress[part];
    std::unique_ptr<QBase> default_qrule_face = QBase::build(default_quad_type, dim - 1, default_quad_order);
    default_qrule_face->allow_rules_with_negative_weights = d_allow_rules_with_negative_weights;
    fe.attachQuadratureRuleFace(default_qrule_face.get());
    fe.evalNormalsFace();
    fe.evalQuadraturePointsFace();
    fe.evalQuadratureWeightsFace();
    const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, vars, vars, &X_ghost_vec);

    const size_t num_PK1_fcns = d_PK1_stress_fcn_data[part].size();
    std::vector<std::vector<size_t> > PK1_fcn_system_idxs(num_PK1_fcns);
    std::vector<std::vector<SystemData> > PK1_stress_fcn_ghosted_system_data(num_PK1_fcns);
    std::vector<std::vector<std::unique_ptr<NumericVector<double> > > > PK1_stress_fcn_ghosted_system_vecs(
        num_PK1_fcns);
    for (unsigned int k = 0; k < num_PK1_fcns; ++k)
    {
        build_ib_ghosted_system_data(PK1_stress_fcn_ghosted_system_data[k],
                                     PK1_stress_fcn_ghosted_system_vecs[k],
                                     d_PK1_stress_fcn_data[part][k].system_data,
                                     d_primary_fe_data_managers[part]);
        fe.setupInterpolatedSystemDataIndexes(
            PK1_fcn_system_idxs[k], PK1_stress_fcn_ghosted_system_data[k], &equation_systems);
    }

    std::vector<size_t> surface_force_fcn_system_idxs;
    std::vector<SystemData> surface_force_fcn_ghosted_system_data;
    std::vector<std::unique_ptr<NumericVector<double> > > surface_force_fcn_ghosted_system_vecs;
    build_ib_ghosted_system_data(surface_force_fcn_ghosted_system_data,
                                 surface_force_fcn_ghosted_system_vecs,
                                 d_lag_surface_force_fcn_data[part].system_data,
                                 d_primary_fe_data_managers[part]);
    fe.setupInterpolatedSystemDataIndexes(
        surface_force_fcn_system_idxs, surface_force_fcn_ghosted_system_data, &equation_systems);

    std::vector<size_t> surface_pressure_fcn_system_idxs;
    std::vector<SystemData> surface_pressure_fcn_ghosted_system_data;
    std::vector<std::unique_ptr<NumericVector<double> > > surface_pressure_fcn_ghosted_system_vecs;
    build_ib_ghosted_system_data(surface_pressure_fcn_ghosted_system_data,
                                 surface_pressure_fcn_ghosted_system_vecs,
                                 d_lag_surface_pressure_fcn_data[part].system_data,
                                 d_primary_fe_data_managers[part]);
    fe.setupInterpolatedSystemDataIndexes(
        surface_pressure_fcn_system_idxs, surface_pressure_fcn_ghosted_system_data, &equation_systems);

    fe.init();

    const std::vector<libMesh::Point>& q_point_face = fe.getQuadraturePointsFace();
    const std::vector<double>& JxW_face = fe.getQuadratureWeightsFace();
    const std::vector<libMesh::Point>& normal_face = fe.getNormalsFace();

    const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
    const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
        fe.getGradVarInterpolation();

    std::vector<std::vector<const std::vector<double>*> > PK1_var_data(num_PK1_fcns);
    std::vector<std::vector<const std::vector<VectorValue<double> >*> > PK1_grad_var_data(num_PK1_fcns);
    std::vector<const std::vector<double>*> surface_force_var_data, surface_pressure_var_data;
    std::vector<const std::vector<VectorValue<double> >*> surface_force_grad_var_data, surface_pressure_grad_var_data;

    // Loop over the patches to spread the transmission elastic force density
    // onto the grid.
    const std::vector<std::vector<Elem*> >& active_patch_element_map =
        d_primary_fe_data_managers[part]->getActivePatchElementMap();
    TensorValue<double> PP, FF, FF_inv_trans;
    VectorValue<double> F, F_s, n, x;
    double P;
    std::vector<double> T_bdry, x_bdry;
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const std::vector<Elem*>& patch_elems = active_patch_element_map[local_patch_num];
        const size_t num_active_patch_elems = patch_elems.size();
        if (num_active_patch_elems == 0) continue;

        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);

        // Loop over the elements and compute the values to be spread and the
        // positions of the quadrature points.
        T_bdry.clear();
        x_bdry.clear();
        int qp_offset = 0;

        IBTK::QuadratureCache side_quad_cache(dim - 1);
        for (size_t e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const bool touches_physical_bdry = has_physical_bdry(elem, boundary_info, G_dof_map);
            if (!touches_physical_bdry) continue;

            fe.reinit(elem);
            fe.collectDataForInterpolation(elem);
            const boost::multi_array<double, 2>& X_node = fe.getElemData(elem, X_sys_idx);

            // Set the length scale for each side (and, therefore, the quadrature
            // rule) to be the same as the length scale for the element.
            using quad_key_type = quadrature_key_type;
            // This duplicates the logic in updateSpreadQuadratureRule so that
            // we can get the key here
            const quad_key_type elem_quad_key = getQuadratureKey(d_spread_spec[part].quad_type,
                                                                 d_spread_spec[part].quad_order,
                                                                 d_spread_spec[part].use_adaptive_quadrature,
                                                                 d_spread_spec[part].point_density,
                                                                 d_spread_spec[part].allow_rules_with_negative_weights,
                                                                 elem,
                                                                 X_node,
                                                                 patch_dx_min);
            // TODO: surely there is a better way to get the type of a side
            // than this!
            std::unique_ptr<Elem> side_ptr = elem->build_side_ptr(0);
            const quad_key_type side_quad_key = std::make_tuple(side_ptr->type(),
                                                                d_spread_spec[part].quad_type,
                                                                std::get<2>(elem_quad_key),
                                                                std::get<3>(elem_quad_key));
            libMesh::QBase& side_quadrature = side_quad_cache[side_quad_key];

            // Loop over the element boundaries.
            for (unsigned short int side = 0; side < elem->n_sides(); ++side)
            {
                // Skip non-physical boundaries.
                if (!is_physical_bdry(elem, side, boundary_info, G_dof_map)) continue;

                // Skip Dirichlet boundaries.
                if (is_dirichlet_bdry(elem, side, boundary_info, G_dof_map)) continue;

                fe.attachQuadratureRuleFace(&side_quadrature);
                fe.reinit(elem, side);
                fe.interpolate(elem, side);
                const unsigned int n_qp = side_quadrature.n_points();
                T_bdry.resize(T_bdry.size() + NDIM * n_qp);
                x_bdry.resize(x_bdry.size() + NDIM * n_qp);
                for (unsigned int qp = 0; qp < n_qp; ++qp, ++qp_offset)
                {
                    const libMesh::Point& X = q_point_face[qp];
                    const std::vector<double>& x_data = fe_interp_var_data[qp][X_sys_idx];
                    const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
                    get_x_and_FF(x, FF, x_data, grad_x_data);
                    const double J = std::abs(FF.det());
                    tensor_inverse_transpose(FF_inv_trans, FF, NDIM);
                    const libMesh::VectorValue<double>& N = normal_face[qp];
                    n = (FF_inv_trans * N).unit();

                    F.zero();

                    for (unsigned int k = 0; k < num_PK1_fcns; ++k)
                    {
                        if (d_PK1_stress_fcn_data[part][k].fcn)
                        {
                            // Compute the value of the first Piola-Kirchhoff stress
                            // tensor at the quadrature point and compute the
                            // corresponding force.
                            fe.setInterpolatedDataPointers(
                                PK1_var_data[k], PK1_grad_var_data[k], PK1_fcn_system_idxs[k], elem, qp);
                            d_PK1_stress_fcn_data[part][k].fcn(PP,
                                                               FF,
                                                               x,
                                                               X,
                                                               elem,
                                                               PK1_var_data[k],
                                                               PK1_grad_var_data[k],
                                                               data_time,
                                                               d_PK1_stress_fcn_data[part][k].ctx);
                            F -= PP * normal_face[qp] * JxW_face[qp];
                        }
                    }

                    if (d_lag_surface_pressure_fcn_data[part].fcn)
                    {
                        // Compute the value of the pressure at the quadrature
                        // point and compute the corresponding force.
                        fe.setInterpolatedDataPointers(surface_pressure_var_data,
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
                                                                  side,
                                                                  surface_pressure_var_data,
                                                                  surface_pressure_grad_var_data,
                                                                  data_time,
                                                                  d_lag_surface_pressure_fcn_data[part].ctx);
                        F -= P * J * FF_inv_trans * normal_face[qp] * JxW_face[qp];
                    }

                    if (d_lag_surface_force_fcn_data[part].fcn)
                    {
                        // Compute the value of the surface force at the
                        // quadrature point and compute the corresponding force.
                        fe.setInterpolatedDataPointers(surface_force_var_data,
                                                       surface_force_grad_var_data,
                                                       surface_force_fcn_system_idxs,
                                                       elem,
                                                       qp);
                        d_lag_surface_force_fcn_data[part].fcn(F_s,
                                                               n,
                                                               N,
                                                               FF,
                                                               x,
                                                               X,
                                                               elem,
                                                               side,
                                                               surface_force_var_data,
                                                               surface_force_grad_var_data,
                                                               data_time,
                                                               d_lag_surface_force_fcn_data[part].ctx);
                        F += F_s * JxW_face[qp];
                    }

                    // Remote the normal component of the boundary force when needed.
                    if (!integrate_normal_force) F -= (F * n) * n;

                    // Remote the tangential component of the boundary force when needed.
                    if (!integrate_tangential_force) F -= (F - (F * n) * n);

                    const int idx = NDIM * qp_offset;
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        T_bdry[idx + i] = F(i);
                    }
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        x_bdry[idx + i] = x(i);
                    }
                }
            }
        }

        if (qp_offset == 0) continue;

        // Spread the boundary forces to the grid.
        const std::string& spread_kernel_fcn = d_spread_spec[part].kernel_fcn;
        const Box<NDIM> spread_box = patch->getBox();
        Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_data_idx);

        FEDataManager::zeroExteriorValues(*patch_geom, x_bdry, T_bdry, NDIM);
        LEInteractor::spread(f_data, T_bdry, NDIM, x_bdry, NDIM, patch, spread_box, spread_kernel_fcn);
    }

    return;
} // spreadTransmissionForceDensity

void
IBFEMethod::imposeJumpConditions(const int f_data_idx,
                                 PetscVector<double>& G_ghost_vec,
                                 PetscVector<double>& X_ghost_vec,
                                 const double data_time,
                                 const unsigned int part)
{
    if (!d_split_normal_force && !d_split_tangential_force) return;

    // Check to see if we need to integrate the normal surface force.
    const bool integrate_normal_force =
        d_split_normal_force && d_use_jump_conditions && !d_stress_normalization_part[part];
    if (!integrate_normal_force) return;

    // Extract the mesh.
    EquationSystems& equation_systems = *d_primary_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems.get_mesh();
    const BoundaryInfo& boundary_info = mesh.get_boundary_info();
    const unsigned int dim = mesh.mesh_dimension();
    TBOX_ASSERT(dim == NDIM);

    // Extract the FE systems and DOF maps, and setup the FE object.
    auto& G_system = equation_systems.get_system<ExplicitSystem>(FORCE_SYSTEM_NAME);
    DofMap& G_dof_map = G_system.get_dof_map();
    FEType G_fe_type = G_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(G_dof_map.variable_type(d) == G_fe_type);
    }
    std::vector<std::vector<unsigned int> > G_dof_indices(NDIM);
    auto& X_system = equation_systems.get_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
    DofMap& X_dof_map = X_system.get_dof_map();
    std::vector<int> vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
    std::vector<int> no_vars;

    FEDataInterpolation fe(dim, d_primary_fe_data_managers[part]->getFEData());
    fe.evalQuadraturePointsFace();
    fe.evalNormalsFace();
    const size_t G_sys_idx = fe.registerInterpolatedSystem(G_system, vars, no_vars, &G_ghost_vec);
    const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, vars, vars, &X_ghost_vec);

    const size_t num_PK1_fcns = d_PK1_stress_fcn_data[part].size();
    std::vector<std::vector<size_t> > PK1_fcn_system_idxs(num_PK1_fcns);
    std::vector<std::vector<SystemData> > PK1_stress_fcn_ghosted_system_data(num_PK1_fcns);
    std::vector<std::vector<std::unique_ptr<NumericVector<double> > > > PK1_stress_fcn_ghosted_system_vecs(
        num_PK1_fcns);
    for (unsigned int k = 0; k < num_PK1_fcns; ++k)
    {
        build_ib_ghosted_system_data(PK1_stress_fcn_ghosted_system_data[k],
                                     PK1_stress_fcn_ghosted_system_vecs[k],
                                     d_PK1_stress_fcn_data[part][k].system_data,
                                     d_primary_fe_data_managers[part]);
        fe.setupInterpolatedSystemDataIndexes(
            PK1_fcn_system_idxs[k], PK1_stress_fcn_ghosted_system_data[k], &equation_systems);
    }

    std::vector<size_t> surface_force_fcn_system_idxs;
    std::vector<SystemData> surface_force_fcn_ghosted_system_data;
    std::vector<std::unique_ptr<NumericVector<double> > > surface_force_fcn_ghosted_system_vecs;
    build_ib_ghosted_system_data(surface_force_fcn_ghosted_system_data,
                                 surface_force_fcn_ghosted_system_vecs,
                                 d_lag_surface_force_fcn_data[part].system_data,
                                 d_primary_fe_data_managers[part]);
    fe.setupInterpolatedSystemDataIndexes(
        surface_force_fcn_system_idxs, surface_force_fcn_ghosted_system_data, &equation_systems);

    std::vector<size_t> surface_pressure_fcn_system_idxs;
    std::vector<SystemData> surface_pressure_fcn_ghosted_system_data;
    std::vector<std::unique_ptr<NumericVector<double> > > surface_pressure_fcn_ghosted_system_vecs;
    build_ib_ghosted_system_data(surface_pressure_fcn_ghosted_system_data,
                                 surface_pressure_fcn_ghosted_system_vecs,
                                 d_lag_surface_pressure_fcn_data[part].system_data,
                                 d_primary_fe_data_managers[part]);
    fe.setupInterpolatedSystemDataIndexes(
        surface_pressure_fcn_system_idxs, surface_pressure_fcn_ghosted_system_data, &equation_systems);

    fe.init();

    const std::vector<libMesh::Point>& q_point_face = fe.getQuadraturePointsFace();
    const std::vector<libMesh::Point>& normal_face = fe.getNormalsFace();

    const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
    const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
        fe.getGradVarInterpolation();

    std::vector<std::vector<const std::vector<double>*> > PK1_var_data(num_PK1_fcns);
    std::vector<std::vector<const std::vector<VectorValue<double> >*> > PK1_grad_var_data(num_PK1_fcns);
    std::vector<const std::vector<double>*> surface_force_var_data, surface_pressure_var_data;
    std::vector<const std::vector<VectorValue<double> >*> surface_force_grad_var_data, surface_pressure_grad_var_data;

    // Loop over the patches to impose jump conditions on the Eulerian grid that
    // are determined from the interior and transmission elastic force
    // densities.
    const std::vector<std::vector<Elem*> >& active_patch_element_map =
        d_primary_fe_data_managers[part]->getActivePatchElementMap();
    assertStructureOnFinestLevel();
    const int level_num = d_primary_fe_data_managers[part]->getFinestPatchLevelNumber();
    TensorValue<double> PP, FF, FF_inv_trans;
    VectorValue<double> G, F, F_s, n;
    std::vector<libMesh::Point> X_node_cache, x_node_cache;
    IBTK::Point x_min, x_max;
    std::vector<std::vector<unsigned int> > side_dof_indices(NDIM);
    std::vector<libMesh::Point> intersection_ref_coords;
    std::vector<SideIndex<NDIM> > intersection_indices;
    std::vector<std::pair<double, libMesh::Point> > intersections;
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
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const x_lower = patch_geom->getXLower();
        const double* const dx = patch_geom->getDx();

        SideData<NDIM, int> num_intersections(patch_box, 1, IntVector<NDIM>(0));
        num_intersections.fillAll(0);

        // Loop over the elements.
        for (size_t e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const bool touches_physical_bdry = has_physical_bdry(elem, boundary_info, G_dof_map);
            if (!touches_physical_bdry) continue;

            fe.reinit(elem);
            fe.collectDataForInterpolation(elem);

            // Loop over the element boundaries.
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                // Skip non-physical boundaries.
                if (!is_physical_bdry(elem, side, boundary_info, G_dof_map)) continue;

                // Skip Dirichlet boundaries.
                if (is_dirichlet_bdry(elem, side, boundary_info, G_dof_map)) continue;

                // Construct a side element.
                std::unique_ptr<Elem> side_elem = elem->build_side_ptr(side, /*proxy*/ false);
                const unsigned int n_node_side = side_elem->n_nodes();
                for (int d = 0; d < NDIM; ++d)
                {
                    X_dof_map.dof_indices(side_elem.get(), side_dof_indices[d], d);
                }

                // Cache the nodal and physical coordinates of the side element,
                // determine the bounding box of the current configuration of
                // the side element, and set the nodal coordinates to correspond
                // to the physical coordinates.
                X_node_cache.resize(n_node_side);
                x_node_cache.resize(n_node_side);
                x_min = IBTK::Point::Constant(std::numeric_limits<double>::max());
                x_max = IBTK::Point::Constant(-std::numeric_limits<double>::max());
                for (unsigned int k = 0; k < n_node_side; ++k)
                {
                    X_node_cache[k] = side_elem->point(k);
                    libMesh::Point& x = x_node_cache[k];
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        x(d) = X_ghost_vec(side_dof_indices[d][k]);
                        x_min[d] = std::min(x_min[d], x(d));
                        x_max[d] = std::max(x_max[d], x(d));
                    }
                    side_elem->point(k) = x;
                }
                Box<NDIM> box(IndexUtilities::getCellIndex(&x_min[0], grid_geom, ratio),
                              IndexUtilities::getCellIndex(&x_max[0], grid_geom, ratio));
                box.grow(IntVector<NDIM>(1));
                box = box * patch_box;

                // Loop over coordinate directions and look for intersections
                // with the background fluid grid.
                intersection_ref_coords.clear();
                intersection_indices.clear();
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    // Setup a unit vector pointing in the coordinate direction
                    // of interest.
                    VectorValue<double> q;
                    q(axis) = 1.0;

                    // Loop over the relevant range of indices.
                    Box<NDIM> axis_box = box;
                    axis_box.lower(axis) = 0;
                    axis_box.upper(axis) = 0;
                    for (BoxIterator<NDIM> b(axis_box); b; b++)
                    {
                        const hier::Index<NDIM>& i_c = b();
                        libMesh::Point r;
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            r(d) =
                                (d == axis ? 0.0 :
                                             x_lower[d] + dx[d] * (static_cast<double>(i_c(d) - patch_lower[d]) + 0.5));
                        }
#if (NDIM == 2)
                        intersect_line_with_edge(intersections, static_cast<Edge*>(side_elem.get()), r, q);
#endif
#if (NDIM == 3)
                        intersect_line_with_face(intersections, static_cast<Face*>(side_elem.get()), r, q);
#endif
                        for (const auto& intersection : intersections)
                        {
                            const libMesh::Point x = r + intersection.first * q;
                            SideIndex<NDIM> i_s(i_c, axis, 0);
                            i_s(axis) = std::floor((x(axis) - x_lower[axis]) / dx[axis] + 0.5) + patch_lower[axis];
                            intersection_ref_coords.push_back(intersection.second);
                            intersection_indices.push_back(i_s);
                            num_intersections(i_s) += 1;
                        }
                    }
                }

                // Restore the element coordinates.
                for (unsigned int k = 0; k < n_node_side; ++k)
                {
                    side_elem->point(k) = X_node_cache[k];
                }

                // If there are no intersection points, then continue to the
                // next side.
                if (intersection_ref_coords.empty()) continue;

                // Evaluate the jump conditions and apply them to the Eulerian
                // grid.
                const bool impose_dp_dn_jumps = false;
                static const double TOL = std::sqrt(std::numeric_limits<double>::epsilon());
                fe.reinit(elem, side, TOL, &intersection_ref_coords);
                fe.interpolate(elem, side);
                const size_t n_qp = intersection_ref_coords.size();
                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    const SideIndex<NDIM>& i_s = intersection_indices[qp];
                    const unsigned int axis = i_s.getAxis();
                    const libMesh::Point& X = q_point_face[qp];
                    const std::vector<double>& x_data = fe_interp_var_data[qp][X_sys_idx];
                    const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
                    libMesh::VectorValue<double> x;
                    get_x_and_FF(x, FF, x_data, grad_x_data);
                    const double J = std::abs(FF.det());
                    tensor_inverse_transpose(FF_inv_trans, FF, NDIM);
                    const libMesh::VectorValue<double>& N = normal_face[qp];
                    n = (FF_inv_trans * N).unit();
                    const double dA_da = 1.0 / (J * (FF_inv_trans * N) * n);
                    const std::vector<double>& G_data = fe_interp_var_data[qp][G_sys_idx];
                    std::copy(G_data.begin(), G_data.end(), &G(0));
#if !defined(NDEBUG)
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            const double x_lower_bound = x_lower[d] +
                                                         (static_cast<double>(i_s(d) - patch_lower[d]) - 0.5) * dx[d] -
                                                         std::sqrt(std::numeric_limits<double>::epsilon());
                            const double x_upper_bound = x_lower[d] +
                                                         (static_cast<double>(i_s(d) - patch_lower[d]) + 0.5) * dx[d] +
                                                         std::sqrt(std::numeric_limits<double>::epsilon());
                            TBOX_ASSERT(x_lower_bound <= x(d) && x(d) <= x_upper_bound);
                        }
                        else
                        {
                            const double x_intersection =
                                x_lower[d] + (static_cast<double>(i_s(d) - patch_lower[d]) + 0.5) * dx[d];
                            const double x_interp = x(d);
                            const double rel_diff =
                                std::abs(x_intersection - x_interp) /
                                std::max(1.0, std::max(std::abs(x_intersection), std::abs(x_interp)));
                            TBOX_ASSERT(rel_diff <= std::sqrt(std::numeric_limits<double>::epsilon()));
                        }
                    }
#endif
                    F.zero();

                    for (unsigned int k = 0; k < num_PK1_fcns; ++k)
                    {
                        if (d_PK1_stress_fcn_data[part][k].fcn)
                        {
                            // Compute the value of the first Piola-Kirchhoff
                            // stress tensor at the quadrature point and compute
                            // the corresponding force.
                            fe.setInterpolatedDataPointers(
                                PK1_var_data[k], PK1_grad_var_data[k], PK1_fcn_system_idxs[k], elem, qp);
                            d_PK1_stress_fcn_data[part][k].fcn(PP,
                                                               FF,
                                                               x,
                                                               X,
                                                               elem,
                                                               PK1_var_data[k],
                                                               PK1_grad_var_data[k],
                                                               data_time,
                                                               d_PK1_stress_fcn_data[part][k].ctx);
                            F -= PP * normal_face[qp];
                        }
                    }

                    if (d_lag_surface_pressure_fcn_data[part].fcn)
                    {
                        // Compute the value of the pressure at the quadrature
                        // point and compute the corresponding force.
                        double P = 0.0;
                        fe.setInterpolatedDataPointers(surface_pressure_var_data,
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
                                                                  side,
                                                                  surface_pressure_var_data,
                                                                  surface_pressure_grad_var_data,
                                                                  data_time,
                                                                  d_lag_surface_pressure_fcn_data[part].ctx);
                        F -= P * J * FF_inv_trans * normal_face[qp];
                    }

                    if (d_lag_surface_force_fcn_data[part].fcn)
                    {
                        // Compute the value of the surface force at the
                        // quadrature point and compute the corresponding force.
                        fe.setInterpolatedDataPointers(surface_force_var_data,
                                                       surface_force_grad_var_data,
                                                       surface_force_fcn_system_idxs,
                                                       elem,
                                                       qp);
                        d_lag_surface_force_fcn_data[part].fcn(F_s,
                                                               n,
                                                               N,
                                                               FF,
                                                               x,
                                                               X,
                                                               elem,
                                                               side,
                                                               surface_force_var_data,
                                                               surface_force_grad_var_data,
                                                               data_time,
                                                               d_lag_surface_force_fcn_data[part].ctx);
                        F += F_s;
                    }

                    F *= dA_da;

                    // Determine the value of the interior force density at the
                    // boundary, and convert it to force per unit volume in the
                    // current configuration.  This value determines the
                    // discontinuity in the normal derivative of the pressure at
                    // the fluid-structure interface.
                    if (impose_dp_dn_jumps)
                    {
                        G /= J;
                    }
                    else
                    {
                        G.zero();
                    }

                    // Impose the jump conditions.
                    const double x_cell_bdry =
                        x_lower[axis] + static_cast<double>(i_s(axis) - patch_lower[axis]) * dx[axis];
                    const double h = x_cell_bdry + (x(axis) > x_cell_bdry ? +0.5 : -0.5) * dx[axis] - x(axis);
                    const double C_p = F * n - h * G(axis);
                    (*f_data)(i_s) += (n(axis) > 0.0 ? +1.0 : -1.0) * (C_p / dx[axis]);
                }
            }
        }
    }
    return;
} // imposeJumpConditions

int
IBFEMethod::getCoarsestPatchLevelNumber() const
{
    int level_number = std::numeric_limits<int>::max();
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        level_number = std::min(d_primary_fe_data_managers[part]->getCoarsestPatchLevelNumber(), level_number);
    }
    TBOX_ASSERT(level_number != std::numeric_limits<int>::max());
    return level_number;
} // getCoarsestPatchLevelNumber

int
IBFEMethod::getFinestPatchLevelNumber() const
{
    int level_number = std::numeric_limits<int>::min();
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        level_number = std::max(d_primary_fe_data_managers[part]->getFinestPatchLevelNumber(), level_number);
    }
    TBOX_ASSERT(level_number != std::numeric_limits<int>::min());
    return level_number;
} // getFinestPatchLevelNumber

SAMRAI::xfer::RefineSchedule<NDIM>&
IBFEMethod::getProlongationSchedule(const int level_number, const int coarse_data_idx, const int fine_data_idx)
{
    const auto key = std::make_pair(level_number, std::make_pair(coarse_data_idx, fine_data_idx));
    if (d_prolongation_schedules.count(key) == 0)
    {
        Pointer<PatchHierarchy<NDIM> > hierarchy =
            d_use_scratch_hierarchy ? d_secondary_hierarchy->getSecondaryHierarchy() : d_hierarchy;
        Pointer<RefineAlgorithm<NDIM> > refine_algorithm = new RefineAlgorithm<NDIM>();
        Pointer<RefineOperator<NDIM> > refine_op;

        Pointer<hier::Variable<NDIM> > f_var;
        VariableDatabase<NDIM>::getDatabase()->mapIndexToVariable(coarse_data_idx, f_var);
        {
            Pointer<hier::Variable<NDIM> > f_var_2;
            VariableDatabase<NDIM>::getDatabase()->mapIndexToVariable(fine_data_idx, f_var_2);
            // These should be the same variable
            TBOX_ASSERT(&*f_var == &*f_var_2);
        }

        Pointer<CellVariable<NDIM, double> > f_cc_var = f_var;
        Pointer<SideVariable<NDIM, double> > f_sc_var = f_var;
        const bool cc_data = f_cc_var;
        const bool sc_data = f_sc_var;
        TBOX_ASSERT(cc_data || sc_data);
        if (cc_data)
        {
            Pointer<CartesianGridGeometry<NDIM> > geometry = hierarchy->getGridGeometry();
            refine_op = geometry->lookupRefineOperator(f_var, "CONSERVATIVE_LINEAR_REFINE");
        }
        else
            refine_op = new CartSideDoubleRT0Refine();
        TBOX_ASSERT(refine_op);
        refine_algorithm->registerRefine(fine_data_idx, coarse_data_idx, fine_data_idx, refine_op);
        Pointer<PatchLevel<NDIM> > fine_level = hierarchy->getPatchLevel(level_number + 1);
        // We can ignore the fifth argument since we don't need to deal with
        // forces outside the physical domain (this was handled previously by
        // f_phys_bdry_op). In particular we don't need it since we don't care
        // about ghost force values (they aren't used in the solver).
        Pointer<RefineSchedule<NDIM> > schedule =
            refine_algorithm->createSchedule(fine_level, nullptr, level_number, hierarchy);
        d_prolongation_schedules[key] = schedule;
    }
    return *d_prolongation_schedules[key];
} // getProlongationSchedule

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBFEMethod::commonConstructor(const Pointer<Database>& input_db, int max_levels)
{
    // Keep track of the maximum possible level number.
    d_max_level_number = max_levels - 1;

    // Indicate that all parts are active.
    const auto n_parts = d_meshes.size();
    d_part_is_active.resize(n_parts, true);

    // Set some default values.
    const bool allow_rules_with_negative_weights = true;
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
                                                      use_nodal_quadrature,
                                                      allow_rules_with_negative_weights);
    d_default_spread_spec = FEDataManager::SpreadSpec("IB_4",
                                                      QGAUSS,
                                                      INVALID_ORDER,
                                                      use_adaptive_quadrature,
                                                      point_density,
                                                      use_nodal_quadrature,
                                                      allow_rules_with_negative_weights);

    // Initialize function data to NULL.
    d_lag_body_source_part.resize(n_parts, false);
    d_lag_body_source_fcn_data.resize(n_parts);
    d_direct_forcing_kinematics_data.resize(n_parts, Pointer<IBFEDirectForcingKinematics>(nullptr));

    // Indicate that all of the parts do NOT use stress normalization by default.
    d_stress_normalization_part.resize(n_parts, false);

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Determine how to compute weak forms.
    d_include_normal_stress_in_weak_form = d_split_normal_force;
    d_include_tangential_stress_in_weak_form = d_split_tangential_force;
    d_include_normal_surface_forces_in_weak_form = !d_split_normal_force;
    d_include_tangential_surface_forces_in_weak_form = !d_split_tangential_force;

    // Set up the interaction spec objects.
    d_interp_spec.resize(n_parts, d_default_interp_spec);
    d_spread_spec.resize(n_parts, d_default_spread_spec);
    d_workload_spec.resize(n_parts, d_default_workload_spec);

    // Setup timers.
    auto set_timer = [&](const char* name) { return TimerManager::getManager()->getTimer(name); };
    IBAMR_DO_ONCE(t_preprocess_integrate_data = set_timer("IBAMR::IBFEMethod::preprocessIntegrateData()");
                  t_postprocess_integrate_data = set_timer("IBAMR::IBFEMethod::postprocessIntegrateData()");
                  t_interpolate_velocity = set_timer("IBAMR::IBFEMethod::interpolateVelocity()");
                  t_compute_lagrangian_force = set_timer("IBAMR::IBFEMethod::computeLagrangianForce()");
                  t_spread_force = set_timer("IBAMR::IBFEMethod::spreadForce()");
                  t_compute_lagrangian_fluid_source = set_timer("IBAMR::IBFEMethod::computeLagrangianFluidSource()");
                  t_spread_fluid_source = set_timer("IBAMR::IBFEMethod::spreadFluidSource()");
                  t_add_workload_estimate = set_timer("IBAMR::IBFEMethod::addWorkloadEstimate()");
                  t_begin_data_redistribution = set_timer("IBAMR::IBFEMethod::beginDataRedistribution()");
                  t_end_data_redistribution = set_timer("IBAMR::IBFEMethod::endDataRedistribution()");
                  t_apply_gradient_detector = set_timer("IBAMR::IBFEMethod::applyGradientDetector()"););

    return;
} // commonConstructor

void
IBFEMethod::getFromInput(const Pointer<Database>& db, bool /*is_from_restart*/)
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

    if (db->isBool("interp_use_nodal_quadrature"))
        d_default_interp_spec.use_nodal_quadrature = db->getBool("interp_use_nodal_quadrature");
    else if (db->isBool("IB_use_nodal_quadrature"))
        d_default_interp_spec.use_nodal_quadrature = db->getBool("IB_use_nodal_quadrature");
    if (d_default_interp_spec.use_nodal_quadrature)
    {
        d_default_interp_spec.quad_type = QTRAP;
        d_default_interp_spec.quad_order = FIRST;
        d_default_interp_spec.use_adaptive_quadrature = false;
        d_default_interp_spec.use_consistent_mass_matrix = false;
        d_default_interp_spec.allow_rules_with_negative_weights = false;
        // we default to using a lumped mass matrix, but allow users to choose to use a consistent mass matrix (even
        // though it seems like a bad idea)
    }

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

    if (db->isBool("interp_allow_rules_with_negative_weights"))
        d_default_interp_spec.allow_rules_with_negative_weights =
            db->getBool("interp_allow_rules_with_negative_weights");
    else if (db->isBool("IB_allow_rules_with_negative_weights"))
        d_default_interp_spec.allow_rules_with_negative_weights = db->getBool("IB_allow_rules_with_negative_weights");

    // Spreading settings.
    if (db->isString("spread_delta_fcn"))
        d_default_spread_spec.kernel_fcn = db->getString("spread_delta_fcn");
    else if (db->isString("IB_delta_fcn"))
        d_default_spread_spec.kernel_fcn = db->getString("IB_delta_fcn");
    else if (db->isString("spread_kernel_fcn"))
        d_default_spread_spec.kernel_fcn = db->getString("spread_kernel_fcn");
    else if (db->isString("IB_kernel_fcn"))
        d_default_spread_spec.kernel_fcn = db->getString("IB_kernel_fcn");

    if (db->isBool("spread_use_nodal_quadrature"))
        d_default_spread_spec.use_nodal_quadrature = db->getBool("spread_use_nodal_quadrature");
    else if (db->isBool("IB_use_nodal_quadrature"))
        d_default_spread_spec.use_nodal_quadrature = db->getBool("IB_use_nodal_quadrature");
    if (d_default_spread_spec.use_nodal_quadrature)
    {
        d_default_spread_spec.quad_type = QTRAP;
        d_default_spread_spec.quad_order = FIRST;
        d_default_spread_spec.use_adaptive_quadrature = false;
        d_default_spread_spec.allow_rules_with_negative_weights = false;
    }

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

    if (db->isBool("spread_allow_rules_with_negative_weights"))
        d_default_spread_spec.allow_rules_with_negative_weights =
            db->getBool("spread_allow_rules_with_negative_weights");
    else if (db->isBool("IB_allow_rules_with_negative_weights"))
        d_default_spread_spec.allow_rules_with_negative_weights = db->getBool("IB_allow_rules_with_negative_weights");

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

    // FEDataManager settings.
    if (db->keyExists("FEDataManager"))
    {
        d_fe_data_manager_db = db->getDatabase("FEDataManager");
    }
    else
    {
        d_fe_data_manager_db = new InputDatabase("FEDataManager");
    }

    // Workaround: support passing FEProjector to FEDataManager as an
    // undocumented feature until FEDataManager loses its copy of FEProjector.
    if (!d_fe_data_manager_db->keyExists("FEProjector"))
    {
        auto db = d_fe_data_manager_db->putDatabase("FEProjector");
        // Technically we should copy every key, but the only key we actually
        // use is num_fischer_vectors, so do that until we remove this
        // workaround
        if (d_fe_projector_db->keyExists("num_fischer_vectors"))
        {
            db->putInteger("num_fischer_vectors", d_fe_projector_db->getInteger("num_fischer_vectors"));
        }
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

    if (db->keyExists("skip_initial_workload_log"))
        d_skip_initial_workload_log = db->getBool("skip_initial_workload_log");

    if (db->isDouble("epsilon")) d_epsilon = db->getDouble("epsilon");

    d_libmesh_partitioner_type =
        string_to_enum<LibmeshPartitionerType>(db->getStringWithDefault("libmesh_partitioner_type", "LIBMESH_DEFAULT"));
    if (db->keyExists("workload_quad_point_weight"))
    {
        d_default_workload_spec.q_point_weight = db->getDouble("workload_quad_point_weight");
    }
    if (db->keyExists("workload_duplicated_node_weight"))
    {
        d_default_workload_spec.duplicated_node_weight = db->getDouble("workload_duplicated_node_weight");
    }

    d_use_scratch_hierarchy = db->getBoolWithDefault("use_scratch_hierarchy", false);
    if (d_use_scratch_hierarchy)
    {
        if (!db->isDatabase("GriddingAlgorithm") || !db->isDatabase("LoadBalancer"))
        {
            TBOX_ERROR(d_object_name << ": if the scratch hierarchy is enabled "
                                        "then the input database should contain entries for the "
                                        "scratch GriddingAlgorithm and LoadBalancer."
                                     << std::endl);
        }

        d_secondary_hierarchy =
            std::unique_ptr<SecondaryHierarchy>(new SecondaryHierarchy(d_object_name + "::scratch_hierarchy",
                                                                       db->getDatabase("GriddingAlgorithm"),
                                                                       db->getDatabase("LoadBalancer")));
    }

    d_patch_association_cfl = db->getDoubleWithDefault("patch_association_cfl", d_patch_association_cfl);
    TBOX_ASSERT(d_patch_association_cfl > 0.0);

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
                                 << " not found in restart file." << std::endl);
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
    const int part_is_active_arr_size = db->getInteger("part_is_active_arr_size");
    d_part_is_active.resize(part_is_active_arr_size);
    std::unique_ptr<bool[]> part_is_active_arr{ new bool[d_part_is_active.size()] };
    db->getBoolArray("part_is_active_arr", part_is_active_arr.get(), d_part_is_active.size());
    std::copy(part_is_active_arr.get(), part_is_active_arr.get() + part_is_active_arr_size, d_part_is_active.begin());
    return;
} // getFromRestart

void
IBFEMethod::assertStructureOnFinestLevel() const
{
    TBOX_ASSERT(getFinestPatchLevelNumber() == d_hierarchy->getFinestLevelNumber() &&
                getCoarsestPatchLevelNumber() == d_hierarchy->getFinestLevelNumber());
}

void
IBFEMethod::reinitElementMappings()
{
    // Store the coordinates at which we performed the last reinit.
    d_X_vecs->copy("solution", { "last_patch_elem_assoc" });
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        d_primary_fe_data_managers[part]->reinitElementMappings();
        if (d_use_scratch_hierarchy)
        {
            d_scratch_fe_data_managers[part]->setPatchHierarchy(d_secondary_hierarchy->getSecondaryHierarchy());
            d_scratch_fe_data_managers[part]->reinitElementMappings();
        }
    }
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
