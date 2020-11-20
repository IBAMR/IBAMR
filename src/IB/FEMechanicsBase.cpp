// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
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
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/FEDataInterpolation.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/FEProjector.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/LibMeshSystemVectors.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"

#include "tbox/RestartManager.h"

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
#include "libmesh/fe_map.h"
#include "libmesh/fe_interface.h"

#include <libmesh/enum_preconditioner_type.h>
#include <libmesh/enum_solver_type.h>



#include <iterator>
#include <utility>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of FEMechanicsBase restart file data.
const int FE_MECHANICS_BASE_VERSION = 1;

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

inline void
get_FF(libMesh::TensorValue<double>& FF,
       const std::vector<VectorValue<double> >& grad_x_data,
       const unsigned int dim = NDIM)
{
    FF.zero();
    for (unsigned int i = 0; i < dim; ++i)
    {
        for (unsigned int j = 0; j < dim; ++j)
        {
            FF(i, j) = grad_x_data[i](j);
        }
    }
    for (unsigned int i = dim; i < LIBMESH_DIM; ++i)
    {
        FF(i, i) = 1.0;
    }
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
}

inline void
get_Grad_U(libMesh::TensorValue<double>& Grad_U,
           const std::vector<VectorValue<double> >& grad_U_data,
           const unsigned int dim = NDIM)
{
    Grad_U.zero();
    for (unsigned int i = 0; i < dim; ++i)
    {
        for (unsigned int j = 0; j < dim; ++j)
        {
            Grad_U(i, j) = grad_U_data[i](j);
        }
    }
}
} // namespace

const std::string FEMechanicsBase::COORDS_SYSTEM_NAME = "IB coordinates system";
const std::string FEMechanicsBase::COORD_MAPPING_SYSTEM_NAME = "IB coordinate mapping system";
const std::string FEMechanicsBase::FORCE_SYSTEM_NAME = "IB force system";
const std::string FEMechanicsBase::PRESSURE_SYSTEM_NAME = "IB pressure system";
const std::string FEMechanicsBase::VELOCITY_SYSTEM_NAME = "IB velocity system";

/////////////////////////////// PUBLIC ///////////////////////////////////////

FEMechanicsBase::FEMechanicsBase(const std::string& object_name,
                                 const Pointer<Database>& input_db,
                                 MeshBase* mesh,
                                 bool register_for_restart,
                                 const std::string& restart_read_dirname,
                                 unsigned int restart_restore_number)
{
    commonConstructor(object_name,
                      input_db,
                      std::vector<MeshBase*>(1, mesh),
                      register_for_restart,
                      restart_read_dirname,
                      restart_restore_number);
}

FEMechanicsBase::FEMechanicsBase(const std::string& object_name,
                                 const Pointer<Database>& input_db,
                                 const std::vector<MeshBase*>& meshes,
                                 bool register_for_restart,
                                 const std::string& restart_read_dirname,
                                 unsigned int restart_restore_number)
{
    commonConstructor(
        object_name, input_db, meshes, register_for_restart, restart_read_dirname, restart_restore_number);
}

FEMechanicsBase::~FEMechanicsBase()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
}

libMesh::EquationSystems*
FEMechanicsBase::getEquationSystems(const unsigned int part) const
{
    TBOX_ASSERT(d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_meshes.size());
    return d_equation_systems[part].get();
}

std::shared_ptr<FEData>
FEMechanicsBase::getFEData(const unsigned int part) const
{
    TBOX_ASSERT(d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_meshes.size());
    return d_fe_data[part];
}

void
FEMechanicsBase::registerInitialCoordinateMappingFunction(const CoordinateMappingFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_meshes.size());
    d_coordinate_mapping_fcn_data[part] = data;
}

FEMechanicsBase::CoordinateMappingFcnData
FEMechanicsBase::getInitialCoordinateMappingFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_meshes.size());
    return d_coordinate_mapping_fcn_data[part];
}

void
FEMechanicsBase::registerInitialVelocityFunction(const InitialVelocityFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_meshes.size());
    d_initial_velocity_fcn_data[part] = data;
}

FEMechanicsBase::InitialVelocityFcnData
FEMechanicsBase::getInitialVelocityFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_meshes.size());
    return d_initial_velocity_fcn_data[part];
}

void
FEMechanicsBase::registerPK1StressFunction(const PK1StressFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_meshes.size());
    d_PK1_stress_fcn_data[part].push_back(data);
    if (data.quad_type == INVALID_Q_RULE)
    {
        d_PK1_stress_fcn_data[part].back().quad_type = d_default_quad_type_stress[part];
    }
    if (data.quad_order == INVALID_ORDER)
    {
        d_PK1_stress_fcn_data[part].back().quad_order = d_default_quad_order_stress[part];
    }
}

std::vector<FEMechanicsBase::PK1StressFcnData>
FEMechanicsBase::getPK1StressFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_meshes.size());
    return d_PK1_stress_fcn_data[part];
}

void
FEMechanicsBase::registerLagBodyForceFunction(const LagBodyForceFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_meshes.size());
    d_lag_body_force_fcn_data[part] = data;
}

FEMechanicsBase::LagBodyForceFcnData
FEMechanicsBase::getLagBodyForceFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_meshes.size());
    return d_lag_body_force_fcn_data[part];
}

void
FEMechanicsBase::registerLagSurfacePressureFunction(const LagSurfacePressureFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_meshes.size());
    d_lag_surface_pressure_fcn_data[part] = data;
}

FEMechanicsBase::LagSurfacePressureFcnData
FEMechanicsBase::getLagSurfacePressureFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_meshes.size());
    return d_lag_surface_pressure_fcn_data[part];
}

void
FEMechanicsBase::registerLagSurfaceForceFunction(const LagSurfaceForceFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_meshes.size());
    d_lag_surface_force_fcn_data[part] = data;
}

FEMechanicsBase::LagSurfaceForceFcnData
FEMechanicsBase::getLagSurfaceForceFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_meshes.size());
    return d_lag_surface_force_fcn_data[part];
}


/*void
FEMechanicsBase::registerStaticPressurePart(PressureProjectionType projection_type,
                                        VolumetricEnergyDerivativeFcn dU_dJ_fcn,
                                        unsigned int part )
{
	registerStaticPressurePart(projection_type, dU_dJ_fcn, part, nullptr);
}*/


/*void
FEMechanicsBase::registerStaticPressurePart(PressureProjectionType projection_type,
                                        VolumetricEnergyDerivativeFcn dU_dJ_fcn,
                                        unsigned int part,
                                        const std::set< subdomain_id_type >& subdomains_set)
{
	std::vector< std::set< subdomain_id_type > * > v(1);
	v[0] = &subdomains_set_ptr;
	registerStaticPressurePart(projection_type, dU_dJ_fcn, part, v);
}
*/


void
FEMechanicsBase::registerStaticPressurePart(PressureProjectionType projection_type,
                                            VolumetricEnergyDerivativeFcn dU_dJ_fcn,
                                            unsigned int part,
											const std::vector< std::set< subdomain_id_type > * >* subdomains_set_ptr_vec )
{
    TBOX_ASSERT(d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_meshes.size());
    TBOX_ASSERT(!d_dynamic_pressure_part[part]);
    if (d_static_pressure_part[part]) return;
    d_has_static_pressure_parts = true;
    d_static_pressure_part[part] = true;
    auto& P_system = d_equation_systems[part]->add_system<LinearImplicitSystem>(PRESSURE_SYSTEM_NAME);

    // setup subdomain restricted pressure in case of different material models
    int n_pressure_subdomains = 0;
    if( subdomains_set_ptr_vec ) 
    {
        n_pressure_subdomains = subdomains_set_ptr_vec->size();
    	std::cout << "registering static pressure on part " << part << " for " << n_pressure_subdomains << " subdomains" << std::endl;
    }
    else
    {
    	std::cout << subdomains_set_ptr_vec << std::endl;
    }
    // This system has a single variable so we don't need to also specify diagonal coupling
    std::string pressure_variable_name = "P";
    if( n_pressure_subdomains == 0 )
    {
        P_system.add_variable(pressure_variable_name, d_fe_order_pressure[part], d_fe_family_pressure[part]);
    }
    else
    {
    	for(auto && subdomain_set_ptr : *subdomains_set_ptr_vec)
    	{
    		std::string subdomain_id = std::to_string( *(subdomain_set_ptr->begin()) );
    		std::string pressure_variable_name = "P_" + subdomain_id;
    		P_system.add_variable(pressure_variable_name, d_fe_order_pressure[part], d_fe_family_pressure[part], subdomain_set_ptr);
    	}
    }

    // Setup cached system vectors at restart.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    IBTK::setup_system_vectors(d_equation_systems[part].get(),
                               { PRESSURE_SYSTEM_NAME },
                               { "current", "half", "new", "tmp", "RHS Vector" },
                               from_restart);
    // Keep track of method parameters.
    d_static_pressure_proj_type[part] = projection_type;
    d_static_pressure_dU_dJ_fcn[part] = dU_dJ_fcn;
}

void
FEMechanicsBase::registerDynamicPressurePart(PressureProjectionType projection_type,
                                             VolumetricEnergyDerivativeFcn d2U_dJ2_fcn,
                                             unsigned int part)
{
    TBOX_ASSERT(d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_meshes.size());
    TBOX_ASSERT(!d_static_pressure_part[part]);
    if (d_dynamic_pressure_part[part]) return;
    d_has_dynamic_pressure_parts = true;
    d_dynamic_pressure_part[part] = true;
    auto& P_system = d_equation_systems[part]->add_system<ExplicitSystem>(PRESSURE_SYSTEM_NAME);
    // This system has a single variable so we don't need to also specify diagonal coupling
    P_system.add_variable("P", d_fe_order_pressure[part], d_fe_family_pressure[part]);
    // Setup cached system vectors at restart.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    IBTK::setup_system_vectors(d_equation_systems[part].get(),
                               { PRESSURE_SYSTEM_NAME },
                               { "current", "half", "new", "tmp", "RHS Vector" },
                               from_restart);
    // Keep track of method parameters.
    d_dynamic_pressure_proj_type[part] = projection_type;
    d_dynamic_pressure_d2U_dJ2_fcn[part] = d2U_dJ2_fcn;
}

void
FEMechanicsBase::preprocessIntegrateData(double current_time, double new_time, int /*num_cycles*/)
{
    // Keep track of the current time step interval.
    d_current_time = current_time;
    d_new_time = new_time;
    d_half_time = current_time + 0.5 * (new_time - current_time);
}

void
FEMechanicsBase::postprocessIntegrateData(double /*current_time*/, double /*new_time*/, int /*num_cycles*/)
{
    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_half_time = std::numeric_limits<double>::quiet_NaN();

    // Update the coordinate mapping dX = X - s.
    for (unsigned part = 0; part < d_meshes.size(); ++part)
    {
        updateCoordinateMapping(part);
    }
}

void
FEMechanicsBase::initializeFEEquationSystems()
{
    if (d_fe_equation_systems_initialized) return;

    // Set up the coupling matrix that will be used by each system.
    d_diagonal_system_coupling.resize(NDIM);
    for (unsigned int i = 0; i < NDIM; ++i)
        for (unsigned int j = 0; j < NDIM; ++j) d_diagonal_system_coupling(i, j) = i == j ? 1 : 0;

    doInitializeFEEquationSystems();
    d_fe_equation_systems_initialized = true;
}

void
FEMechanicsBase::initializeFEData()
{
    if (d_fe_data_initialized) return;

    initializeFEEquationSystems();
    doInitializeFEData(RestartManager::getManager()->isFromRestart());
    d_fe_data_initialized = true;
}

void
FEMechanicsBase::reinitializeFEData()
{
    TBOX_ASSERT(d_fe_data_initialized);
    doInitializeFEData(true);
}

void
FEMechanicsBase::putToDatabase(Pointer<Database> db)
{
    db->putInteger("FE_MECHANICS_BASE_VERSION", FE_MECHANICS_BASE_VERSION);
    db->putBool("d_use_consistent_mass_matrix", d_use_consistent_mass_matrix);
    db->putString("d_libmesh_partitioner_type", enum_to_string<LibmeshPartitionerType>(d_libmesh_partitioner_type));
}

void
FEMechanicsBase::writeFEDataToRestartFile(const std::string& restart_dump_dirname, unsigned int time_step_number)
{
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        const std::string& file_name =
            libmesh_restart_file_name(restart_dump_dirname, time_step_number, part, d_libmesh_restart_file_extension);
        const XdrMODE xdr_mode = (d_libmesh_restart_file_extension == "xdr" ? ENCODE : WRITE);
        const int write_mode = EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA;
        d_equation_systems[part]->write(file_name,
                                        xdr_mode,
                                        write_mode,
                                        /*partition_agnostic*/ true);
    }
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void
FEMechanicsBase::computeStaticPressure(PetscVector<double>& P_vec,
                                       PetscVector<double>& X_vec,
                                       const double data_time,
                                       const unsigned int part)
{
	if (!d_static_pressure_part[part]) return;
	std::cout << "calling computeStaticPressure" << std::endl;

    double tau = d_static_pressure_stab_param;

	bool using_epsilon_map = ( d_static_pressure_stab_param_vector_map[part].size() > 0 ) ? true : false;
	//const PressureProjectionType& proj_type = d_static_pressure_proj_type[part];
    const VolumetricEnergyDerivativeFcn& dU_dJ_fcn = d_static_pressure_dU_dJ_fcn[part];

    int map_size = 0;
    if ( part < d_static_pressure_kappa_vector_map.size() )
    {
        map_size = d_static_pressure_kappa_vector_map[part].size();
    }

    // Extract the mesh.
    EquationSystems& equation_systems = *d_equation_systems[part];
    const MeshBase& mesh = equation_systems.get_mesh();
    const Parallel::Communicator& comm = mesh.comm();
    const unsigned int dim = mesh.mesh_dimension();
    const BoundaryInfo& boundary_info = *mesh.boundary_info;

    // Setup extra data needed to compute stresses/forces.

    // Extract the FE systems and DOF maps, and setup the FE objects.
    auto& P_system = equation_systems.get_system<LinearImplicitSystem>(PRESSURE_SYSTEM_NAME);
    const libMesh::DofMap& P_dof_map = P_system.get_dof_map();
    libMesh::FEType P_fe_type = P_dof_map.variable_type(0);
    auto* P_rhs_vec = static_cast<PetscVector<double>*>(P_system.rhs);
    P_rhs_vec->zero();



    auto& X_system = equation_systems.get_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
    	//const libMesh::DofMap& X_dof_map = X_system.get_dof_map();

    std::vector<int> X_vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) X_vars[d] = d;
    std::vector<int> P_vars(P_system.n_vars());
    for(unsigned int i = 0; i < P_vars.size(); ++i) P_vars[i] = i;
    std::vector<int> no_vars = {};
//
//    // set up RHS
//    {
//
//		std::unique_ptr<libMesh::FEBase> fe (libMesh::FEBase::build(dim, P_fe_type));
//		libMesh::QGauss qrule (dim, FIRST);
//		fe->attach_quadrature_rule (&qrule);
//		const std::vector<libMesh::Real> & JxW = fe->get_JxW();
//		const std::vector<std::vector<libMesh::Real>> & phi = fe->get_phi();
//
//		// Setup global and elemental right-hand-side vectors.
//		auto* P_rhs_vec = static_cast<PetscVector<double>*>(P_system.rhs);
//		P_rhs_vec->zero();
//		DenseVector<double> P_rhs_e;
//
//		// Deformation gradient tensor
//		TensorValue<double> FF;
//		// FF_{ij} = \partial_j x_i
//		//    = [ d1 x_1, d2 x_1, d3 x_1 ]
//		//    = [ d1 x_2, d2 x_2, d3 x_2 ]
//		//    = [ d1 x_3, d2 x_3, d3 x_3 ]
//
//
//		std::vector<libMesh::dof_id_type> pressure_dof_indices;
//		std::vector<libMesh::dof_id_type> X_dof_indices;
//		const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
//		const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
//		for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
//		{
//			Elem* const elem = *el_it;
//			unsigned int blockID = elem->subdomain_id();
//
//			P_dof_map.dof_indices(elem, pressure_dof_indices);
//			const auto n_basis = pressure_dof_indices.size();
//
//			if( n_basis > 0)
//			{
//				// check if there is a pressure associated with this element
//				P_rhs_e.resize(static_cast<int>(n_basis));
//
//				fe->reinit (elem);
//
//				if(map_size > 0) d_static_pressure_kappa = d_static_pressure_kappa_vector_map[part][blockID];
//				const unsigned int n_qp = qrule_pressure_rhs->n_points();
//				for (unsigned int qp = 0; qp < n_qp; ++qp)
//				{
//					FF *= 0.0;
//					for(unsigned int idim = 0; idim < 3; ++idim)
//					{
//						if(idim < dim)
//						{
//							auto dx_dX = X_system.point_gradient (idim, qrule->qp(qp), elem);
//							FF(idim, 0) = dx_dX(0);
//							FF(idim, 1) = dx_dX(1);
//							FF(idim, 2) = dx_dX(2);
//						}
//						else
//						{
//							FF(idim, idim) = 1.0;
//						}
//					}
//					double J = FF.det();
//					const double P = (dU_dJ_fcn ? dU_dJ_fcn(J) : -d_static_pressure_kappa * std::log(J));
//					for (unsigned int k = 0; k < n_basis; ++k)
//					{
//						P_rhs_e(k) += P * phi[k][qp] * JxW[qp];
//					}
//				}
//			}
//			// Apply constraints (e.g., enforce periodic boundary conditions)
//			// and add the elemental contributions to the global vector.
//			//copy_dof_ids_to_vector(/*var_num*/ 0, P_dof_indices, dof_id_scratch);
//			//P_dof_map.constrain_element_vector(P_rhs_e, pressure_dof_indices);
//			P_rhs_vec->add_vector(P_rhs_e, pressure_dof_indices);
//		}
//	}

    // If it is the first time through, create the LHS matrix
    // Build solver components.


	double penalty = 1e8;
    if( !d_pressure_matrix_is_assembled[part] )
    {
    	d_pressure_matrix_is_assembled[part] = true;

    	//d_pressure_solver[part].reset(new PetscLinearSolver<double>(comm));
		//d_pressure_matrix[part].reset(new PetscMatrix<double>(comm));
		//d_pressure_matrix[part]->attach_dof_map(P_dof_map);
		//d_pressure_matrix[part]->init();

		PetscMatrix<double> * M_mat = dynamic_cast<PetscMatrix<double>*>(P_system.matrix);
		MatSetOption(M_mat->mat(), MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
		MatSetOption(M_mat->mat(), MAT_SPD, PETSC_TRUE);
		MatSetOption(M_mat->mat(), MAT_SYMMETRY_ETERNAL, PETSC_TRUE);
		MatSetOption(M_mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

		// Loop over the mesh to construct the system matrix.
		DenseMatrix<double> M_e;
		DenseMatrix<double> K_ee;
		DenseMatrix<double> K_en;
		DenseMatrix<double> K_ne;
		DenseMatrix<double> K_nn;
		DenseVector<double> Pi_phi_e;
		std::vector<libMesh::dof_id_type> dof_indices;
		std::vector<libMesh::dof_id_type> dof_indices_neighbor;


		//
		std::unique_ptr<FEBase> fe_matrix (FEBase::build(dim, P_fe_type));
        // A 5th order Gauss quadrature rule for numerical integration.
        QGauss qrule_matrix (dim, FIFTH);
        // Tell the finite element object to use our quadrature rule.
        fe_matrix->attach_quadrature_rule (&qrule_matrix);

		const std::vector<libMesh::Real> & JxW = fe_matrix->get_JxW();
		const std::vector<std::vector<libMesh::Real>> & phi = fe_matrix->get_phi();


        // Add interface conditions
        std::unique_ptr<FEBase> fe_elem_face(FEBase::build(dim, P_fe_type));
        std::unique_ptr<FEBase> fe_neighbor_face(FEBase::build(dim, P_fe_type));
        libMesh::QGauss qface(dim-1, P_fe_type.default_quadrature_order());
        // Tell the finite element object to use our quadrature rule.
        fe_elem_face->attach_quadrature_rule(&qface);
        fe_neighbor_face->attach_quadrature_rule(&qface);

        // Data for surface integrals on the element boundary
        const std::vector<std::vector<Real>> &  phi_face = fe_elem_face->get_phi();
        //const std::vector<std::vector<RealGradient>> & dphi_face = fe_elem_face->get_dphi();
        const std::vector<Real> & JxW_face = fe_elem_face->get_JxW();
        //const std::vector<libMesh::Point> & qface_normals = fe_elem_face->get_normals();
        //const std::vector<libMesh::Point> & qface_points = fe_elem_face->get_xyz();

        // Data for surface integrals on the neighbor boundary
        const std::vector<std::vector<Real>> &  phi_neighbor_face = fe_neighbor_face->get_phi();
        //const std::vector<std::vector<RealGradient>> & dphi_neighbor_face = fe_neighbor_face->get_dphi();


		const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
		const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();

		for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
		{
			Elem* const elem = *el_it;
			fe_matrix->reinit(elem);
			unsigned int blockID = elem->subdomain_id();
			if(using_epsilon_map) tau = d_static_pressure_stab_param_vector_map[part].at(blockID);

			//const auto& dof_indices = dof_map_cache.dof_indices(elem);
			for (unsigned int var_n = 0; var_n < P_dof_map.n_variables(); ++var_n)
			{
				P_dof_map.dof_indices(elem, dof_indices, var_n);
				const auto n_basis = dof_indices.size();
				if(n_basis <= 0 ) continue;
				//const auto& dof_indices_var = dof_indices[var_n];
				//const auto n_basis = static_cast<unsigned int>(dof_indices_var.size());
				M_e.resize(n_basis, n_basis);
				Pi_phi_e.resize(n_basis);
				const unsigned int n_qp = qrule_matrix.n_points();

				const double vol_e = elem->volume();
				for (unsigned int i = 0; i < n_basis; ++i)
				{
					for (unsigned int qp = 0; qp < n_qp; ++qp)
					{
						Pi_phi_e(i) += phi[i][qp] * JxW[qp] / vol_e;
					}
				}
				for (unsigned int i = 0; i < n_basis; ++i)
				{
					for (unsigned int j = 0; j < n_basis; ++j)
					{
						for (unsigned int qp = 0; qp < n_qp; ++qp)
						{
							M_e(i, j) += ((phi[i][qp] * phi[j][qp]) +
										  tau * (phi[i][qp] - Pi_phi_e(i)) * (phi[j][qp] - Pi_phi_e(j))) *
										 JxW[qp];
						}
					}
				}
				//M_e.print();
				P_system.matrix->add_matrix(M_e, dof_indices, dof_indices);
			}

			for (auto side : elem->side_index_range())
			{
				if (elem->neighbor_ptr(side) == nullptr) continue;
	            const Elem * neighbor = elem->neighbor_ptr(side);
				unsigned int neighbor_blockID = neighbor->subdomain_id();
				// if we are on the same block ID skip the assembly
				if(neighbor_blockID == blockID ) continue;
				for (unsigned int var_n = 0; var_n < P_dof_map.n_variables(); ++var_n)
				{
					P_dof_map.dof_indices(elem, dof_indices, var_n);
					const auto n_basis = dof_indices.size();
					if(n_basis <= 0 ) continue;
					P_dof_map.dof_indices(neighbor, dof_indices_neighbor, var_n);
		            auto n_basis_neighbor = dof_indices_neighbor.size();

		            // if the number of dofs on opposite elements is different
		            // it means we are on an interior interface
		            if( n_basis > 0 && n_basis != n_basis_neighbor )
		            {
		            	  // get the actual dofs for the neighbor element
			              P_dof_map.dof_indices(neighbor, dof_indices_neighbor);
			              n_basis_neighbor = dof_indices_neighbor.size();

		                  K_ne.resize (n_basis_neighbor, n_basis);
		                  K_en.resize (n_basis, n_basis_neighbor);
		                  K_ee.resize (n_basis, n_basis);
		                  K_nn.resize (n_basis_neighbor, n_basis_neighbor);

		                  // Pointer to the element side
		                  std::unique_ptr<const libMesh::Elem> elem_side (elem->build_side_ptr(side));

		                  // The quadrature point locations on the neighbor side
		                  std::vector<libMesh::Point> qface_neighbor_point;
		                  // The quadrature point locations on the element side
		                  std::vector<libMesh::Point > qface_point;
		                  // Reinitialize shape functions on the element side
		                  fe_elem_face->reinit(elem, side);
		                  // Get the physical locations of the element quadrature points
		                  qface_point = fe_elem_face->get_xyz();

		                  libMesh::FEInterface::inverse_map (elem->dim(), P_fe_type, neighbor,
		                		  	  	  	  qface_point,
											  qface_neighbor_point);

		                  // Calculate the neighbor element shape functions at those locations
		                  fe_neighbor_face->reinit(neighbor, &qface_neighbor_point);

		                  for (unsigned int qp=0; qp<qface.n_points(); qp++)
		                  {
		                      for (unsigned int i=0; i<n_basis; i++)
		                          for (unsigned int j=0; j<n_basis; j++)
		                        	  K_ee(i,j) += JxW_face[qp] * penalty * phi_face[j][qp]*phi_face[i][qp];

		                      for (unsigned int i=0; i<n_basis_neighbor; i++)
		                          for (unsigned int j=0; j<n_basis_neighbor; j++)
		                        	  K_nn(i,j) +=
		                                JxW_face[qp] * penalty * phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp];

		                      for (unsigned int i=0; i<n_basis_neighbor; i++)
		                          for (unsigned int j=0; j<n_basis; j++)
		                              K_ne(i,j) -= JxW_face[qp] * penalty * phi_face[j][qp]*phi_neighbor_face[i][qp];

		                      for (unsigned int i=0; i<n_basis; i++)
		                          for (unsigned int j=0; j<n_basis_neighbor; j++)
		                              K_en(i,j) -= JxW_face[qp] * penalty * phi_face[i][qp]*phi_neighbor_face[j][qp];
		                  }
		            }
		            //K_ee.print();
					P_system.matrix->add_matrix(K_ee, dof_indices, dof_indices);
					P_system.matrix->add_matrix(K_nn, dof_indices_neighbor, dof_indices_neighbor);
					P_system.matrix->add_matrix(K_ne, dof_indices_neighbor, dof_indices);
					P_system.matrix->add_matrix(K_en, dof_indices, dof_indices_neighbor);
				}
			}
		}
    	P_system.matrix->close();

        // Setup the solver.
    	PetscLinearSolver<double> * solver = dynamic_cast< PetscLinearSolver<double> * >(P_system.get_linear_solver());
    	solver->reuse_preconditioner(true);
    	solver->set_preconditioner_type(JACOBI_PRECOND);
    	solver->set_solver_type(MINRES);
    	solver->init();

    	//P_system.matrix->print();
    }

	std::cout << "calling computeStaticPressure: matrix done" << std::endl;

	// Assemble RHS
	{
		// handle the stress contributions.  These are handled separately
		// because each stress function may use a different quadrature rule.
		const size_t num_PK1_fcns = d_PK1_stress_fcn_data[part].size();
		for (unsigned int k = 0; k < num_PK1_fcns; ++k)
		{
			if (!d_PK1_stress_fcn_data[part][k].fcn)
				continue;

			FEDataInterpolation fe(dim, d_fe_data[part]);
			std::unique_ptr<QBase> qrule = QBase::build(d_PK1_stress_fcn_data[part][k].quad_type, dim, d_PK1_stress_fcn_data[part][k].quad_order);
			std::unique_ptr<QBase> qrule_face = QBase::build(d_PK1_stress_fcn_data[part][k].quad_type, dim - 1, d_PK1_stress_fcn_data[part][k].quad_order);
			fe.attachQuadratureRule(qrule.get());
			fe.attachQuadratureRuleFace(qrule_face.get());
			fe.evalNormalsFace();
			fe.evalQuadraturePoints();
			fe.evalQuadraturePointsFace();
			fe.evalQuadratureWeights();
			fe.evalQuadratureWeightsFace();
			fe.registerSystem(P_system, P_vars, no_vars);
			const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, no_vars, X_vars, &X_vec);
			std::vector < size_t > PK1_fcn_system_idxs;
			fe.setupInterpolatedSystemDataIndexes(PK1_fcn_system_idxs, d_PK1_stress_fcn_data[part][k].system_data, &equation_systems);
			fe.init();

			//const std::vector<libMesh::Point> &q_point = fe.getQuadraturePoints();
			const std::vector<double> &JxW = fe.getQuadratureWeights();
			const std::vector<std::vector<double> > &phi = fe.getPhi(P_fe_type);
			//const std::vector<std::vector<VectorValue<double> > > &dphi = fe.getDphi(P_fe_type);
			const std::vector < libMesh::Point >&  qface_point = fe.getQuadraturePointsFace();


			const std::vector<libMesh::Point> &q_point_face = fe.getQuadraturePointsFace();
			const std::vector<double> &JxW_face = fe.getQuadratureWeightsFace();
			const std::vector<libMesh::Point> &normal_face = fe.getNormalsFace();
			const std::vector<std::vector<double> > &phi_face = fe.getPhiFace(P_fe_type);

			const std::vector<std::vector<std::vector<double> > > &fe_interp_var_data = fe.getVarInterpolation();
			const std::vector<std::vector<std::vector<VectorValue<double> > > > &fe_interp_grad_var_data = fe.getGradVarInterpolation();

			std::vector<const std::vector<double>*> PK1_var_data;
			std::vector<const std::vector<VectorValue<double> >*> PK1_grad_var_data;

			// Loop over the elements to compute the right-hand side vector.  This
			// is computed via
			//
			//    rhs_k = -int{PP(s,t) grad phi_k(s)}ds + int{PP(s,t) N(s,t)
			//    phi_k(s)}dA(s)
			//
			// This right-hand side vector is used to solve for the nodal values of
			// the interior elastic force density.
			TensorValue<double> PP, FF, FF_inv_trans;
			TensorValue<double> PPe, FFe, FFe_inv_trans;
			TensorValue<double> PPn, FFn, FFn_inv_trans;
			TensorValue<double> sigma_n, sigma_e;
			VectorValue<double> F, F_qp, n, x;
			VectorValue<double> xe, xn, ne, nn;
			DenseVector<double> P_rhs_e;
			DenseVector<double> P_rhs_n;

			std::vector < libMesh::dof_id_type > pressure_dof_indices;
			std::vector < libMesh::dof_id_type > dof_indices_neighbor;

			// neighbor element
			FEDataInterpolation fe_neighbor(dim, d_fe_data[part]);
			fe_neighbor.attachQuadratureRule(qrule.get());
			fe_neighbor.attachQuadratureRuleFace(qrule_face.get());
			fe_neighbor.evalNormalsFace();
			fe_neighbor.evalQuadraturePoints();
			fe_neighbor.evalQuadraturePointsFace();
			fe_neighbor.evalQuadratureWeights();
			fe_neighbor.evalQuadratureWeightsFace();

			fe_neighbor.registerSystem(P_system, P_vars, no_vars);
			const size_t X_sys_idx_neighbor = fe_neighbor.registerInterpolatedSystem(X_system, no_vars, X_vars, &X_vec);
			std::vector < size_t > PK1_fcn_system_idxs_neighbor;
			fe_neighbor.setupInterpolatedSystemDataIndexes(PK1_fcn_system_idxs_neighbor, d_PK1_stress_fcn_data[part][k].system_data, &equation_systems);
			fe_neighbor.init();

			//const std::vector<libMesh::Point> &q_point_face_neighbor = fe_neighbor.getQuadraturePointsFace();
			//const std::vector<double> &JxW_face_neighbor = fe_neighbor.getQuadratureWeightsFace();
			//const std::vector<libMesh::Point> &normal_face_neighbor = fe_neighbor.getNormalsFace();
			const std::vector<std::vector<double> > &phi_face_neighbor = fe_neighbor.getPhiFace(P_fe_type);

			const std::vector<std::vector<std::vector<double> > > &fe_interp_var_data_neighbor = fe_neighbor.getVarInterpolation();
			const std::vector<std::vector<std::vector<VectorValue<double> > > > &fe_interp_grad_var_data_neighbor = fe_neighbor.getGradVarInterpolation();

			std::vector<const std::vector<double>*> PK1_var_data_neighbor;
			std::vector<const std::vector<VectorValue<double> >*> PK1_grad_var_data_neighbor;

			// Loop over elements
		    FEDataManager::SystemDofMapCache& P_dof_map_cache = *d_fe_data[part]->getDofMapCache(PRESSURE_SYSTEM_NAME);
			const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
			const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
			for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
			{
				Elem *const elem = *el_it;
				unsigned int blockID = elem->subdomain_id();


		        const auto& dof_indices_test = P_dof_map_cache.dof_indices(elem);
		        std::cout << " dof indices " << std::endl;
			        for(size_t m = 0; m < dof_indices_test.shape()[0]; ++m)
		        {
			        for(size_t n = 0; n < dof_indices_test.shape()[1]; ++n)
			        {
			        	std::cout << " " << dof_indices_test[m][n] << " " << std::flush;
			        }
			        std::cout << std::endl;
		        }
		        std::cout << std::endl;

				P_dof_map.dof_indices(elem, pressure_dof_indices);
				const auto n_basis = pressure_dof_indices.size();
				if (n_basis <= 0)
					continue;

				P_rhs_e.resize(static_cast<int>(n_basis));

				if (map_size > 0)
					d_static_pressure_kappa = d_static_pressure_kappa_vector_map[part][blockID];

				fe.reinit(elem);
				fe.collectDataForInterpolation(elem);
				fe.interpolate(elem);
				const unsigned int n_qp = qrule->n_points();

				for (unsigned int qp = 0; qp < n_qp; ++qp)
				{
					//const libMesh::Point &X = q_point[qp];
					const std::vector<double> &x_data = fe_interp_var_data[qp][X_sys_idx];
					const std::vector<VectorValue<double> > &grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
					get_x_and_FF(x, FF, x_data, grad_x_data);
					//SR: HACK
					FF *= 0.0;
					for (unsigned int idim = 0; idim < 3; ++idim)
					{
						if (idim < dim)
						{
							x(idim) = X_system.point_value(idim, qrule->qp(qp), elem);
							auto dx_dX = X_system.point_gradient(idim, qrule->qp(qp), elem);
							FF(idim, 0) = dx_dX(0);
							FF(idim, 1) = dx_dX(1);
							FF(idim, 2) = dx_dX(2);
						}
						else
						{
							x(idim) = 0.0;
							FF(idim, idim) = 1.0;
						}
					}

					double J = FF.det();
					const double P = (dU_dJ_fcn ? dU_dJ_fcn(J) : -d_static_pressure_kappa * std::log(J));
					for (unsigned int k = 0; k < n_basis; ++k)
					{
						P_rhs_e(k) += P * phi[k][qp] * JxW[qp];
					}
				}

				P_rhs_vec->add_vector(P_rhs_e, pressure_dof_indices);
				std::cout << "elem: " << elem->id() << " U(J) done " << std::endl;

				// Loop over the element boundaries.
				for (unsigned int side = 0; side < elem->n_sides(); ++side)
				{
					// Skip non-physical boundaries.
					if (elem->neighbor_ptr(side) == nullptr) continue;

					Elem *neighbor = elem->neighbor_ptr(side);
					unsigned int neighbor_blockID = neighbor->subdomain_id();
					std::cout << "side: " << side << ", blockID: " << blockID << ", n blockID: " << neighbor_blockID << std::endl;
					if (neighbor_blockID == blockID)
						continue;

					for (unsigned int var_n = 0; var_n < P_dof_map.n_variables(); ++var_n)
					{
						P_dof_map.dof_indices(elem, pressure_dof_indices, var_n);
						auto ndofs = pressure_dof_indices.size();
						P_dof_map.dof_indices(neighbor, dof_indices_neighbor, var_n);
						auto ndofs_neighbor = dof_indices_neighbor.size();
						if (ndofs > 0 && ndofs != ndofs_neighbor)
						{
							std::cout << "assembling jump in stress  " << std::endl;
							fe.reinit(elem, side);
							std::cout << "interpolate  " << std::endl;
							fe.interpolate(elem, side);

							// get the actual dofs for the neighbor element
							std::cout << "n dofmap  " << std::endl;
							P_dof_map.dof_indices(neighbor, dof_indices_neighbor);
							//auto n_basis_neighbor = dof_indices_neighbor.size();

							// The quadrature point locations on the neighbor side
							std::vector < libMesh::Point > qface_neighbor_point;
							// The quadrature point locations on the element side

							std::cout << "inverse map  " << std::endl;
			                  libMesh::FEInterface::inverse_map (elem->dim(), P_fe_type, neighbor,
			                		  	  	  	  qface_point,
												  qface_neighbor_point);

								std::cout << " n reinit  " << std::endl;
							fe_neighbor.reinit(neighbor, side, libMesh::TOLERANCE, &qface_neighbor_point);
							std::cout << " n interpolate  " << std::endl;
							fe_neighbor.interpolate(neighbor, side);
							std::cout << " n_qp_face  " << std::endl;
							const unsigned int n_qp_face = qrule_face->n_points();
							const size_t n_basis_face = phi_face.size();

							std::cout << "loop over qp  " << std::endl;
							for (unsigned int qp = 0; qp < n_qp_face; ++qp)
							{
								F.zero();

								// Compute the value of the first Piola-Kirchhoff stress
								// tensor at the quadrature point and add the corresponding
								// traction force to the right-hand-side vector.
								if (d_PK1_stress_fcn_data[part][k].fcn)
								{
									PPe *= 0;
									FFe *= 0;
									FFe_inv_trans *= 0;
									xe *= 0;
									const libMesh::Point &Xe = q_point_face[qp];
									const std::vector<double> &xe_data = fe_interp_var_data[qp][X_sys_idx];
									const std::vector<VectorValue<double> > &grad_xe_data = fe_interp_grad_var_data[qp][X_sys_idx];
									get_x_and_FF(xe, FF, xe_data, grad_xe_data);
									tensor_inverse_transpose(FFe_inv_trans, FFe, NDIM);

									fe.setInterpolatedDataPointers(PK1_var_data, PK1_grad_var_data, PK1_fcn_system_idxs, elem, qp);
									d_PK1_stress_fcn_data[part][k].fcn(PP, FF, xe, Xe, elem, PK1_var_data, PK1_grad_var_data, data_time, d_PK1_stress_fcn_data[part][k].ctx);

									PPn *= 0;
									FFn *= 0;
									FFn_inv_trans *= 0;
									xn *= 0;
									const libMesh::Point &Xn = q_point_face[qp];
									const std::vector<double> &xn_data = fe_interp_var_data_neighbor[qp][X_sys_idx_neighbor];
									const std::vector<VectorValue<double> > &grad_xn_data = fe_interp_grad_var_data_neighbor[qp][X_sys_idx_neighbor];
									get_x_and_FF(xn, FFn, xn_data, grad_xn_data);
									tensor_inverse_transpose(FFn_inv_trans, FFn, NDIM);

									fe_neighbor.setInterpolatedDataPointers(PK1_var_data_neighbor, PK1_grad_var_data_neighbor, PK1_fcn_system_idxs, neighbor, qp);
									d_PK1_stress_fcn_data[part][k].fcn(PPn, FFn, xn, Xn, neighbor, PK1_var_data_neighbor, PK1_grad_var_data_neighbor, data_time,
											d_PK1_stress_fcn_data[part][k].ctx);
								} // evaluate PK1

								n = (FFe_inv_trans * normal_face[qp]).unit();

								sigma_e = 1.0 / FFe.det() * PPe * FFe_inv_trans;
								sigma_n = 1.0 / FFn.det() * PPn * FFn_inv_trans;
								double stress_jump = penalty * ((sigma_e - sigma_n) * n) * n;

								P_rhs_e *= 0;
								P_rhs_n *= 0;
								// Add the boundary forces to the right-hand-side vector.
								for (unsigned int i = 0; i < n_basis_face; ++i)
								{
									P_rhs_e(i) += JxW_face[qp] * stress_jump * phi_face[i][qp];
									P_rhs_n(i) += JxW_face[qp] * stress_jump * phi_face_neighbor[i][qp];
								} // assemble local RHS
							} //loop over quadrature points
						} // check if ndofs = ndofs_neighbor
					} // loop over variables in pressure system

					P_rhs_vec->add_vector(P_rhs_e, pressure_dof_indices);
					P_rhs_vec->add_vector(P_rhs_e, dof_indices_neighbor);

				} // loop

			} // loop over elements
		} // loop over PK1 functions
	} // Assemble RHS

	std::cout << "calling computeStaticPressure: rhs done" << std::endl;

	PetscLinearSolver<double> * solver = dynamic_cast< PetscLinearSolver<double> * >(P_system.get_linear_solver());
	PetscMatrix<double> * M_mat = dynamic_cast<PetscMatrix<double>*>(P_system.matrix);
	//P_rhs_vec->print();
    PetscBool rtol_set;
    double runtime_rtol;
    int ierr;
    double tol = 1e-6;
    int max_its = 100;
    ierr = PetscOptionsGetReal(nullptr, "", "-ksp_rtol", &runtime_rtol, &rtol_set);
    IBTK_CHKERRQ(ierr);
    PetscBool max_it_set;
    int runtime_max_it;
    ierr = PetscOptionsGetInt(nullptr, "", "-ksp_max_it", &runtime_max_it, &max_it_set);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetFromOptions( solver->ksp());
    IBTK_CHKERRQ(ierr);
    solver->solve(*M_mat, *M_mat, P_vec, *P_rhs_vec, rtol_set ? runtime_rtol : tol, max_it_set ? runtime_max_it : max_its);
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(solver->ksp(), &reason);
    IBTK_CHKERRQ(ierr);
    //bool converged = reason > 0;



    //if (close_U) U_vec.close();
    //system.get_dof_map().enforce_constraints_exactly(system, &U_vec);


    //std::cout << "static pressure RHS: " << P_rhs_vec->l2_norm() << std::endl;
    //P_rhs_vec->print(std::cout);
    // Solve for P.
    //P_rhs_vec->print();
    //for(auto && v : d_static_pressure_stab_param_vector_map[part]) std::cout << "[ " << v.first << ", " << v.second << std::endl;

    //std::cout << "map_size: " << map_size << std::endl;
//    switch (proj_type)
//    {
//    case CONSISTENT_PROJECTION:
//        d_fe_projectors[part]->computeL2Projection(
//            P_vec, *P_rhs_vec, PRESSURE_SYSTEM_NAME, /*use_consistent_mass_matrix*/ true);
//        break;
//    case LUMPED_PROJECTION:
//        d_fe_projectors[part]->computeL2Projection(
//            P_vec, *P_rhs_vec, PRESSURE_SYSTEM_NAME, /*use_consistent_mass_matrix*/ false);
//        break;
//    case STABILIZED_PROJECTION:
//        //if(map_size > 0)
//        //{
//        //	std::cout << "passing" << std::endl;
//        d_fe_projectors[part]->computeStabilizedL2Projection(
//            P_vec, *P_rhs_vec, PRESSURE_SYSTEM_NAME, d_static_pressure_stab_param, d_static_pressure_stab_param_vector_map[part]);
//        //}
//        //else
//        //{
//        //	std::cout << "not passing" << std::endl;
//        //d_fe_projectors[part]->computeStabilizedL2Projection(
//        //    P_vec, *P_rhs_vec, PRESSURE_SYSTEM_NAME, d_static_pressure_stab_param);
//        //}
//        break;
//    default:
//        TBOX_ERROR("unsupported pressure projection type\n");
//    }
//    //P_vec.print();
}

void
FEMechanicsBase::computeDynamicPressureRateOfChange(PetscVector<double>& dP_dt_vec,
                                                    PetscVector<double>& X_vec,
                                                    PetscVector<double>& U_vec,
                                                    const double /*data_time*/,
                                                    const unsigned int part)
{
    if (!d_dynamic_pressure_part[part]) return;
    const PressureProjectionType& proj_type = d_dynamic_pressure_proj_type[part];
    const VolumetricEnergyDerivativeFcn& d2U_dJ2_fcn = d_dynamic_pressure_d2U_dJ2_fcn[part];

    // Extract the mesh.
    EquationSystems& equation_systems = *d_equation_systems[part];
    const MeshBase& mesh = equation_systems.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // Setup extra data needed to compute stresses/forces.

    // Extract the FE systems and DOF maps, and setup the FE objects.
    auto& P_system = equation_systems.get_system<ExplicitSystem>(PRESSURE_SYSTEM_NAME);
    const DofMap& P_dof_map = P_system.get_dof_map();
    FEDataManager::SystemDofMapCache& P_dof_map_cache = *d_fe_data[part]->getDofMapCache(PRESSURE_SYSTEM_NAME);
    FEType P_fe_type = P_dof_map.variable_type(0);
    std::vector<int> P_vars = { 0 };
    std::vector<int> no_vars = {};
    auto& X_system = equation_systems.get_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
    std::vector<int> X_vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) X_vars[d] = d;
    auto& U_system = equation_systems.get_system<ExplicitSystem>(VELOCITY_SYSTEM_NAME);
    std::vector<int> U_vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) U_vars[d] = d;

    FEDataInterpolation fe(dim, d_fe_data[part]);
    std::unique_ptr<QBase> qrule =
        QBase::build(d_default_quad_type_pressure[part], dim, d_default_quad_order_pressure[part]);
    fe.attachQuadratureRule(qrule.get());
    fe.evalQuadraturePoints();
    fe.evalQuadratureWeights();
    fe.registerSystem(P_system, P_vars, no_vars);
    const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, no_vars, X_vars, &X_vec);
    const size_t U_sys_idx = fe.registerInterpolatedSystem(U_system, no_vars, U_vars, &U_vec);
    fe.init();

    //const std::vector<libMesh::Point>& q_point = fe.getQuadraturePoints();
    const std::vector<double>& JxW = fe.getQuadratureWeights();
    const std::vector<std::vector<double> >& phi = fe.getPhi(P_fe_type);

    //const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
    const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
        fe.getGradVarInterpolation();

    // Setup global and elemental right-hand-side vectors.
    auto* dP_dt_rhs_vec = static_cast<PetscVector<double>*>(P_system.rhs);
    dP_dt_rhs_vec->zero();
    DenseVector<double> dP_dt_rhs_e;

    TensorValue<double> FF, FF_inv_trans, Grad_U;
    std::vector<libMesh::dof_id_type> dof_id_scratch;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        const auto& P_dof_indices = P_dof_map_cache.dof_indices(elem);
        dP_dt_rhs_e.resize(static_cast<int>(P_dof_indices[0].size()));
        fe.reinit(elem);
        fe.collectDataForInterpolation(elem);
        fe.interpolate(elem);
        const unsigned int n_qp = qrule->n_points();
        const size_t n_basis = phi.size();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
            get_FF(FF, grad_x_data);
            FF_inv_trans = tensor_inverse_transpose(FF);
            const std::vector<VectorValue<double> >& grad_U_data = fe_interp_grad_var_data[qp][U_sys_idx];
            get_Grad_U(Grad_U, grad_U_data);
            double J = FF.det();
            const double dP_dt =
                (d2U_dJ2_fcn ? J * d2U_dJ2_fcn(J) : -d_dynamic_pressure_kappa) * FF_inv_trans.contract(Grad_U);
            for (unsigned int k = 0; k < n_basis; ++k)
            {
                dP_dt_rhs_e(k) += dP_dt * phi[k][qp] * JxW[qp];
            }
        }

        // Apply constraints (e.g., enforce periodic boundary conditions)
        // and add the elemental contributions to the global vector.
        copy_dof_ids_to_vector(/*var_num*/ 0, P_dof_indices, dof_id_scratch);
        P_dof_map.constrain_element_vector(dP_dt_rhs_e, dof_id_scratch);
        dP_dt_rhs_vec->add_vector(dP_dt_rhs_e, dof_id_scratch);
    }

    // Solve for P.
    switch (proj_type)
    {
    case CONSISTENT_PROJECTION:
        d_fe_projectors[part]->computeL2Projection(
            dP_dt_vec, *dP_dt_rhs_vec, PRESSURE_SYSTEM_NAME, /*use_consistent_mass_matrix*/ true);
        break;
    case LUMPED_PROJECTION:
        d_fe_projectors[part]->computeL2Projection(
            dP_dt_vec, *dP_dt_rhs_vec, PRESSURE_SYSTEM_NAME, /*use_consistent_mass_matrix*/ false);
        break;
    case STABILIZED_PROJECTION:
        d_fe_projectors[part]->computeStabilizedL2Projection(
            dP_dt_vec, *dP_dt_rhs_vec, PRESSURE_SYSTEM_NAME, d_dynamic_pressure_stab_param, d_dynamic_pressure_stab_param_vector_map[part]);
        break;
    default:
        TBOX_ERROR("unsupported pressure projection type\n");
    }
}

void
FEMechanicsBase::assembleInteriorForceDensityRHS(PetscVector<double>& F_rhs_vec,
                                                 PetscVector<double>& X_vec,
                                                 PetscVector<double>* P_vec,
                                                 const double data_time,
                                                 const unsigned int part)
{
    const bool using_pressure = P_vec != nullptr;

    // Extract the mesh.
    EquationSystems& equation_systems = *d_equation_systems[part];
    const MeshBase& mesh = equation_systems.get_mesh();
    const BoundaryInfo& boundary_info = *mesh.boundary_info;
    const unsigned int dim = mesh.mesh_dimension();
//    // debugging
//    unsigned int n_subdomains = mesh.n_subdomains();
//    std::set< libMesh::subdomain_id_type >  subdomain_ids;
//    mesh.subdomain_ids(subdomain_ids);
//    std::cout << "mesh has " << n_subdomains << std::endl;
//    for(auto && sid : subdomain_ids) std::cout << sid << " ";
//    std::cout << std::endl;
//    //

    // Setup global and elemental right-hand-side vectors.
    auto& F_system = equation_systems.get_system<ExplicitSystem>(FORCE_SYSTEM_NAME);

    // During assembly we sum into ghost regions - this only makes sense if we
    // have a ghosted vector.
    int ierr;
    TBOX_ASSERT(F_rhs_vec.type() == GHOSTED);
    Vec F_rhs_vec_local;
    ierr = VecGhostGetLocalForm(F_rhs_vec.vec(), &F_rhs_vec_local);
    IBTK_CHKERRQ(ierr);
    double* F_rhs_local_soln = nullptr;
    ierr = VecGetArray(F_rhs_vec_local, &F_rhs_local_soln);
    IBTK_CHKERRQ(ierr);
    std::array<DenseVector<double>, NDIM> F_rhs_e;
    std::vector<libMesh::dof_id_type> dof_id_scratch;

//    // debugging
//    std::cout << " assembleInteriorForceDensityRHS test" << std::endl;
//    X_system_tmp.print_info();
//    std::cout << " is F system initialized:" << F_system.is_initialized()  << std::endl;
//    F_system.print_info();
//    unsigned int nvars_F = F_system.n_vars() ;
//    std::cout << " n vars:" << nvars_F << std::endl;
//    for(unsigned int nv = 0; nv < nvars_F; ++nv )
//    {
//    	auto active_on_subdomains = F_system.variable(nv).active_subdomains();
//    	std::cout << "F_" << nv << " is active on " << active_on_subdomains.size() <<  " subdomains: " << std::flush;
//    	for (auto && as : active_on_subdomains) std::cout << as << " " << std::flush;
//    	std::cout << std::endl;
//    }
//    std::cout << " n dofs:" << F_system.n_dofs()  << std::endl;
//    const DofMap& F_dof_map_tmp = F_system.get_dof_map();
//    F_dof_map_tmp.print_info();
//    //
    // First handle the stress contributions.  These are handled separately
    // because each stress function may use a different quadrature rule.
    const size_t num_PK1_fcns = d_PK1_stress_fcn_data[part].size();
    for (unsigned int k = 0; k < num_PK1_fcns; ++k)
    {
        if (!d_PK1_stress_fcn_data[part][k].fcn) continue;

        // Extract the FE systems and DOF maps, and setup the FE object.
        const DofMap& F_dof_map = F_system.get_dof_map();

        FEDataManager::SystemDofMapCache& F_dof_map_cache = *d_fe_data[part]->getDofMapCache(FORCE_SYSTEM_NAME);
        FEType F_fe_type = F_dof_map.variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(F_dof_map.variable_type(d) == F_fe_type);
        }
        auto& X_system = equation_systems.get_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;

        FEDataInterpolation fe(dim, d_fe_data[part]);
        std::unique_ptr<QBase> qrule =
            QBase::build(d_PK1_stress_fcn_data[part][k].quad_type, dim, d_PK1_stress_fcn_data[part][k].quad_order);
        std::unique_ptr<QBase> qrule_face =
            QBase::build(d_PK1_stress_fcn_data[part][k].quad_type, dim - 1, d_PK1_stress_fcn_data[part][k].quad_order);
        fe.attachQuadratureRule(qrule.get());
        fe.attachQuadratureRuleFace(qrule_face.get());
        fe.evalNormalsFace();
        fe.evalQuadraturePoints();
        fe.evalQuadraturePointsFace();
        fe.evalQuadratureWeights();
        fe.evalQuadratureWeightsFace();
        fe.registerSystem(F_system, std::vector<int>(), vars); // compute dphi for the force system
        const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, vars, vars, &X_vec);
        std::vector<size_t> PK1_fcn_system_idxs;
        fe.setupInterpolatedSystemDataIndexes(
            PK1_fcn_system_idxs, d_PK1_stress_fcn_data[part][k].system_data, &equation_systems);
        fe.init();

        const std::vector<libMesh::Point>& q_point = fe.getQuadraturePoints();
        const std::vector<double>& JxW = fe.getQuadratureWeights();
        const std::vector<std::vector<VectorValue<double> > >& dphi = fe.getDphi(F_fe_type);

        const std::vector<libMesh::Point>& q_point_face = fe.getQuadraturePointsFace();
        const std::vector<double>& JxW_face = fe.getQuadratureWeightsFace();
        const std::vector<libMesh::Point>& normal_face = fe.getNormalsFace();
        const std::vector<std::vector<double> >& phi_face = fe.getPhiFace(F_fe_type);

        const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
        const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
            fe.getGradVarInterpolation();

        std::vector<const std::vector<double>*> PK1_var_data;
        std::vector<const std::vector<VectorValue<double> >*> PK1_grad_var_data;

        // Loop over the elements to compute the right-hand side vector.  This
        // is computed via
        //
        //    rhs_k = -int{PP(s,t) grad phi_k(s)}ds + int{PP(s,t) N(s,t)
        //    phi_k(s)}dA(s)
        //
        // This right-hand side vector is used to solve for the nodal values of
        // the interior elastic force density.
        TensorValue<double> PP, FF, FF_inv_trans;
        VectorValue<double> F, F_qp, n, x;
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            const auto& F_dof_indices = F_dof_map_cache.dof_indices(elem);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                F_rhs_e[d].resize(static_cast<int>(F_dof_indices[d].size()));
            }

            fe.reinit(elem);
            fe.collectDataForInterpolation(elem);
            fe.interpolate(elem);
            const unsigned int n_qp = qrule->n_points();
            const size_t n_basis = dphi.size();
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const libMesh::Point& X = q_point[qp];
                const std::vector<double>& x_data = fe_interp_var_data[qp][X_sys_idx];
                const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
                get_x_and_FF(x, FF, x_data, grad_x_data);
                //SR: HACK
                FF *= 0.0;
                for(unsigned int idim = 0; idim < 3; ++idim)
                {
                	if(idim < dim )
					{
                		x(idim) = X_system.point_value(idim, qrule->qp(qp), elem);
						auto dx_dX = X_system.point_gradient (idim, qrule->qp(qp), elem);
						FF(idim, 0) = dx_dX(0);
						FF(idim, 1) = dx_dX(1);
						FF(idim, 2) = dx_dX(2);
					}
                	else
                	{
                		x(idim) = 0.0;
						FF(idim, idim) = 1.0;
                	}
                }
                // Compute the value of the first Piola-Kirchhoff stress tensor
                // at the quadrature point and add the corresponding forces to
                // the right-hand-side vector.
                fe.setInterpolatedDataPointers(PK1_var_data, PK1_grad_var_data, PK1_fcn_system_idxs, elem, qp);
                d_PK1_stress_fcn_data[part][k].fcn(
                    PP, FF, x, X, elem, PK1_var_data, PK1_grad_var_data, data_time, d_PK1_stress_fcn_data[part][k].ctx);
                for (unsigned int basis_n = 0; basis_n < n_basis; ++basis_n)
                {
                    F_qp = -PP * dphi[basis_n][qp] * JxW[qp];
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        F_rhs_e[i](basis_n) += F_qp(i);
                    }
                }
            }

//            for (unsigned int d = 0; d < NDIM; ++d)
//            {
//                F_rhs_e[d].print(std::cout);
//            }


            // Loop over the element boundaries.
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                // Skip non-physical boundaries.
                if (!is_physical_bdry(elem, side, boundary_info, F_dof_map)) continue;

                // Determine if we need to integrate surface forces along this
                // part of the physical boundary; if not, skip the present side.
                const bool at_dirichlet_bdry = is_dirichlet_bdry(elem, side, boundary_info, F_dof_map);
                const bool integrate_normal_stress = (d_include_normal_stress_in_weak_form && !at_dirichlet_bdry) ||
                                                     (!d_include_normal_stress_in_weak_form && at_dirichlet_bdry);
                const bool integrate_tangential_stress =
                    (d_include_tangential_stress_in_weak_form && !at_dirichlet_bdry) ||
                    (!d_include_tangential_stress_in_weak_form && at_dirichlet_bdry);
                if (!integrate_normal_stress && !integrate_tangential_stress) continue;

                fe.reinit(elem, side);
                fe.interpolate(elem, side);
                const unsigned int n_qp_face = qrule_face->n_points();
                const size_t n_basis_face = phi_face.size();
                for (unsigned int qp = 0; qp < n_qp_face; ++qp)
                {
                    const libMesh::Point& X = q_point_face[qp];
                    const std::vector<double>& x_data = fe_interp_var_data[qp][X_sys_idx];
                    const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
                    get_x_and_FF(x, FF, x_data, grad_x_data);
                    tensor_inverse_transpose(FF_inv_trans, FF, NDIM);

                    F.zero();

                    // Compute the value of the first Piola-Kirchhoff stress
                    // tensor at the quadrature point and add the corresponding
                    // traction force to the right-hand-side vector.
                    if (d_PK1_stress_fcn_data[part][k].fcn)
                    {
                        fe.setInterpolatedDataPointers(PK1_var_data, PK1_grad_var_data, PK1_fcn_system_idxs, elem, qp);
                        d_PK1_stress_fcn_data[part][k].fcn(PP,
                                                           FF,
                                                           x,
                                                           X,
                                                           elem,
                                                           PK1_var_data,
                                                           PK1_grad_var_data,
                                                           data_time,
                                                           d_PK1_stress_fcn_data[part][k].ctx);
                        F += PP * normal_face[qp];
                    }

                    n = (FF_inv_trans * normal_face[qp]).unit();

                    if (!integrate_normal_stress)
                    {
                        F -= (F * n) * n; // remove the normal component.
                    }

                    if (!integrate_tangential_stress)
                    {
                        F -= (F - (F * n) * n); // remove the tangential component.
                    }

                    // Add the boundary forces to the right-hand-side vector.
                    for (unsigned int basis_face_n = 0; basis_face_n < n_basis_face; ++basis_face_n)
                    {
                        F_qp = F * phi_face[basis_face_n][qp] * JxW_face[qp];
                        for (unsigned int i = 0; i < NDIM; ++i)
                        {
                            F_rhs_e[i](basis_face_n) += F_qp(i);
                        }
                    }
                }
            }

            // Apply constraints (e.g., enforce periodic boundary conditions)
            // and add the elemental contributions to the global vector.
            for (unsigned int var_n = 0; var_n < NDIM; ++var_n)
            {
                copy_dof_ids_to_vector(var_n, F_dof_indices, dof_id_scratch);
                F_dof_map.constrain_element_vector(F_rhs_e[var_n], dof_id_scratch);
                for (unsigned int j = 0; j < dof_id_scratch.size(); ++j)
                {
                    F_rhs_local_soln[F_rhs_vec.map_global_to_local_index(dof_id_scratch[j])] += F_rhs_e[var_n](j);
                }
            }
        }
    }

    // Now account for any additional force contributions.

    // Extract the FE systems and DOF maps, and setup the FE objects.
    const DofMap& F_dof_map = F_system.get_dof_map();
    FEDataManager::SystemDofMapCache& F_dof_map_cache = *d_fe_data[part]->getDofMapCache(FORCE_SYSTEM_NAME);
    FEType F_fe_type = F_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(F_dof_map.variable_type(d) == F_fe_type);
    }
    auto& X_system = equation_systems.get_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
    System* P_system = using_pressure ? &equation_systems.get_system<ExplicitSystem>(PRESSURE_SYSTEM_NAME) : nullptr;
    std::vector<int> vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
    std::vector<int> P_vars(P_system->n_vars());
    for(unsigned int i = 0; i < P_vars.size(); ++i) P_vars[i] = i;
    std::vector<int> no_vars;

    FEDataInterpolation fe(dim, d_fe_data[part]);
    std::unique_ptr<QBase> qrule = QBase::build(d_default_quad_type_force[part], dim, d_default_quad_order_force[part]);
    std::unique_ptr<QBase> qrule_face =
        QBase::build(d_default_quad_type_force[part], dim - 1, d_default_quad_order_force[part]);
    fe.attachQuadratureRule(qrule.get());
    fe.attachQuadratureRuleFace(qrule_face.get());
    fe.evalNormalsFace();
    fe.evalQuadraturePoints();
    fe.evalQuadraturePointsFace();
    fe.evalQuadratureWeights();
    fe.evalQuadratureWeightsFace();
    fe.registerSystem(F_system, vars, vars); // compute phi and dphi for the force system
    const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, vars, vars, &X_vec);
    //const size_t P_sys_idx = using_pressure ? fe.registerInterpolatedSystem(*P_system, P_vars, no_vars, P_vec) :
    //                                          std::numeric_limits<size_t>::max();
    std::vector<size_t> body_force_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        body_force_fcn_system_idxs, d_lag_body_force_fcn_data[part].system_data, &equation_systems);
    std::vector<size_t> surface_force_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        surface_force_fcn_system_idxs, d_lag_surface_force_fcn_data[part].system_data, &equation_systems);
    std::vector<size_t> surface_pressure_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        surface_pressure_fcn_system_idxs, d_lag_surface_pressure_fcn_data[part].system_data, &equation_systems);
    fe.init();

    const std::vector<libMesh::Point>& q_point = fe.getQuadraturePoints();
    const std::vector<double>& JxW = fe.getQuadratureWeights();
    const std::vector<std::vector<double> >& phi = fe.getPhi(F_fe_type);
    const std::vector<std::vector<VectorValue<double> > >& dphi = fe.getDphi(F_fe_type);

    const std::vector<libMesh::Point>& q_point_face = fe.getQuadraturePointsFace();
    const std::vector<double>& JxW_face = fe.getQuadratureWeightsFace();
    const std::vector<libMesh::Point>& normal_face = fe.getNormalsFace();
    const std::vector<std::vector<double> >& phi_face = fe.getPhiFace(F_fe_type);

    const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
    const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
        fe.getGradVarInterpolation();

    std::vector<const std::vector<double>*> body_force_var_data, surface_force_var_data, surface_pressure_var_data;
    std::vector<const std::vector<VectorValue<double> >*> body_force_grad_var_data, surface_force_grad_var_data,
        surface_pressure_grad_var_data;

    // Loop over the elements to compute the right-hand side vector.
    TensorValue<double> PP, FF, FF_inv_trans;
    VectorValue<double> F, F_b, F_s, F_qp, n, x;
    boost::multi_array<double, 2> X_node;
    boost::multi_array<double, 1> P_node;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        const auto& F_dof_indices = F_dof_map_cache.dof_indices(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            F_rhs_e[d].resize(static_cast<int>(F_dof_indices[d].size()));
        }
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

            // SR: HACK
            FF *= 0.0;
            for(unsigned int idim = 0; idim < 3; ++idim)
            {
            	if(idim < dim )
				{
            		x(idim) = X_system.point_value(idim, qrule->qp(qp), elem);
					auto dx_dX = X_system.point_gradient (idim, qrule->qp(qp), elem);
					FF(idim, 0) = dx_dX(0);
					FF(idim, 1) = dx_dX(1);
					FF(idim, 2) = dx_dX(2);
				}
            	else
            	{
            		x(idim) = 0.0;
					FF(idim, idim) = 1.0;
            	}
            }


            const double J = std::abs(FF.det());
            tensor_inverse_transpose(FF_inv_trans, FF, NDIM);

            if (using_pressure)
            {
                //const double P = fe_interp_var_data[qp][P_sys_idx][0];
            	double P = 0.0;
            	std::vector<libMesh::dof_id_type> dofs;
            	for(unsigned int varn = 0; varn < P_system->n_vars(); varn++)
            	{
                	P_system->get_dof_map().dof_indices(elem, dofs, varn);
                	if(dofs.size() > 0)
            		P = P_system->point_value(varn, qrule->qp(qp), elem);
            	}
                // Compute the value of the first Piola-Kirchhoff stress tensor
                // at the quadrature point and add the corresponding forces to
                // the right-hand-side vector.
                PP = -J * P * FF_inv_trans;
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    F_qp = -PP * dphi[k][qp] * JxW[qp];
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        F_rhs_e[i](k) += F_qp(i);
                    }
                }
            }

            if (d_lag_body_force_fcn_data[part].fcn)
            {
                // Compute the value of the body force at the quadrature
                // point and add the corresponding forces to the
                // right-hand-side vector.
                fe.setInterpolatedDataPointers(
                    body_force_var_data, body_force_grad_var_data, body_force_fcn_system_idxs, elem, qp);
                d_lag_body_force_fcn_data[part].fcn(F_b,
                                                    FF,
                                                    x,
                                                    X,
                                                    elem,
                                                    body_force_var_data,
                                                    body_force_grad_var_data,
                                                    data_time,
                                                    d_lag_body_force_fcn_data[part].ctx);
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    F_qp = F_b * phi[k][qp] * JxW[qp];
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        F_rhs_e[i](k) += F_qp(i);
                    }
                }
            }
        }

        // Loop over the element boundaries.
        for (unsigned int side = 0; side < elem->n_sides(); ++side)
        {
            // Skip non-physical boundaries.
            if (!is_physical_bdry(elem, side, boundary_info, F_dof_map)) continue;

            // Determine if we need to compute surface forces along this
            // part of the physical boundary; if not, skip the present side.
            const bool at_dirichlet_bdry = is_dirichlet_bdry(elem, side, boundary_info, F_dof_map);
            const bool integrate_normal_force = d_include_normal_surface_forces_in_weak_form && !at_dirichlet_bdry;
            const bool integrate_tangential_force =
                d_include_tangential_surface_forces_in_weak_form && !at_dirichlet_bdry;
            if (!integrate_normal_force && !integrate_tangential_force) continue;

            fe.reinit(elem, side);
            fe.interpolate(elem, side);
            const unsigned int n_qp_face = qrule_face->n_points();
            const size_t n_basis_face = phi_face.size();
            for (unsigned int qp = 0; qp < n_qp_face; ++qp)
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

                if (d_lag_surface_pressure_fcn_data[part].fcn)
                {
                    // Compute the value of the pressure at the quadrature
                    // point and add the corresponding force to the
                    // right-hand-side vector.
                    double P = 0;
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
                    F += F_s;
                }

                // Remote the normal component of the boundary force when needed.
                if (!integrate_normal_force) F -= (F * n) * n;

                // Remote the tangential component of the boundary force when needed.
                if (!integrate_tangential_force) F -= (F - (F * n) * n);

                // Add the boundary forces to the right-hand-side vector.
                for (unsigned int k = 0; k < n_basis_face; ++k)
                {
                    F_qp = F * phi_face[k][qp] * JxW_face[qp];
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        F_rhs_e[i](k) += F_qp(i);
                    }
                }
            }
        }

        // Apply constraints (e.g., enforce periodic boundary conditions)
        // and add the elemental contributions to the global vector.
        for (unsigned int var_n = 0; var_n < NDIM; ++var_n)
        {
            copy_dof_ids_to_vector(var_n, F_dof_indices, dof_id_scratch);
            F_dof_map.constrain_element_vector(F_rhs_e[var_n], dof_id_scratch);
            for (unsigned int j = 0; j < dof_id_scratch.size(); ++j)
            {
                F_rhs_local_soln[F_rhs_vec.map_global_to_local_index(dof_id_scratch[j])] += F_rhs_e[var_n](j);
            }
        }
    }

    ierr = VecRestoreArray(F_rhs_vec_local, &F_rhs_local_soln);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostRestoreLocalForm(F_rhs_vec.vec(), &F_rhs_vec_local);
    IBTK_CHKERRQ(ierr);
}

void
FEMechanicsBase::initializeCoordinates(const unsigned int part)
{
    EquationSystems& equation_systems = *d_fe_data[part]->getEquationSystems();
    MeshBase& mesh = equation_systems.get_mesh();
    auto& X_system = equation_systems.get_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
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
    copy_and_synch(X_coords,
                   *X_system.current_local_solution,
                   /*close_v_in*/ false);
}

void
FEMechanicsBase::updateCoordinateMapping(const unsigned int part)
{
    EquationSystems& equation_systems = *d_fe_data[part]->getEquationSystems();
    MeshBase& mesh = equation_systems.get_mesh();
    auto& X_system = equation_systems.get_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    auto& dX_system = equation_systems.get_system<ExplicitSystem>(COORD_MAPPING_SYSTEM_NAME);
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
}

void
FEMechanicsBase::initializeVelocity(const unsigned int part)
{
    EquationSystems& equation_systems = *d_fe_data[part]->getEquationSystems();
    MeshBase& mesh = equation_systems.get_mesh();
    auto& U_system = equation_systems.get_system<ExplicitSystem>(VELOCITY_SYSTEM_NAME);
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
    copy_and_synch(U_vec, *U_system.current_local_solution, /*close_v_in*/ false);
}

void
FEMechanicsBase::setup_system_vectors(EquationSystems* equation_systems,
                                      const std::vector<std::string>& system_names,
                                      const std::vector<std::string>& vector_names)
{
    IBAMR_DEPRECATED_MEMBER_FUNCTION1("FEMechanicsBase", "setup_system_vectors");
    IBTK::setup_system_vectors(
        equation_systems, system_names, vector_names, RestartManager::getManager()->isFromRestart());
}

void
FEMechanicsBase::setup_system_vector(System& system, const std::string& vector_name)
{
    IBAMR_DEPRECATED_MEMBER_FUNCTION1("FEMechanicsBase", "setup_system_vector");
    IBTK::setup_system_vector(system, vector_name, RestartManager::getManager()->isFromRestart());
}

std::string
FEMechanicsBase::libmesh_restart_file_name(const std::string& restart_dump_dirname,
                                           unsigned int time_step_number,
                                           unsigned int part,
                                           const std::string& extension)
{
    std::ostringstream file_name_prefix;
    file_name_prefix << restart_dump_dirname << "/libmesh_data_part_" << part << "." << std::setw(6)
                     << std::setfill('0') << std::right << time_step_number << "." << extension;
    return file_name_prefix.str();
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void
FEMechanicsBase::commonConstructor(const std::string& object_name,
                                   const Pointer<Database>& input_db,
                                   const std::vector<libMesh::MeshBase*>& meshes,
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
    const auto n_parts = d_meshes.size();

    // Set some default values.
    d_fe_order_position.resize(n_parts, INVALID_ORDER);
    d_fe_order_force.resize(n_parts, INVALID_ORDER);
    d_fe_order_pressure.resize(n_parts, INVALID_ORDER);
    d_fe_family_position.resize(n_parts, INVALID_FE);
    d_fe_family_force.resize(n_parts, INVALID_FE);
    d_fe_family_pressure.resize(n_parts, INVALID_FE);
    d_default_quad_type_stress.resize(n_parts, INVALID_Q_RULE);
    d_default_quad_type_force.resize(n_parts, INVALID_Q_RULE);
    d_default_quad_type_pressure.resize(n_parts, INVALID_Q_RULE);
    d_default_quad_order_stress.resize(n_parts, INVALID_ORDER);
    d_default_quad_order_force.resize(n_parts, INVALID_ORDER);
    d_default_quad_order_pressure.resize(n_parts, INVALID_ORDER);

    // Initialize function data to NULL.
    d_coordinate_mapping_fcn_data.resize(n_parts);
    d_initial_velocity_fcn_data.resize(n_parts);
    d_PK1_stress_fcn_data.resize(n_parts);
    d_lag_body_force_fcn_data.resize(n_parts);
    d_lag_surface_pressure_fcn_data.resize(n_parts);
    d_lag_surface_force_fcn_data.resize(n_parts);

    // Indicate that all of the parts do NOT use static pressures by default.
    d_static_pressure_part.resize(n_parts, false);
    d_static_pressure_proj_type.resize(n_parts, UNKNOWN_PRESSURE_TYPE);
    d_static_pressure_dU_dJ_fcn.resize(n_parts, nullptr);

    // Indicate that all of the parts do NOT use dynamic pressures by default.
    d_dynamic_pressure_part.resize(n_parts, false);
    d_dynamic_pressure_proj_type.resize(n_parts, UNKNOWN_PRESSURE_TYPE);
    d_dynamic_pressure_d2U_dJ2_fcn.resize(n_parts, nullptr);

    d_pressure_matrix.resize(n_parts);
    d_pressure_matrix_is_assembled.resize(n_parts, false);
    d_pressure_solver.resize(n_parts);
    // Determine whether we should use first-order or second-order shape
    // functions for each part of the structure.
    for (unsigned int part = 0; part < n_parts; ++part)
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
        mesh_has_first_order_elems = IBTK_MPI::maxReduction(mesh_has_first_order_elems ? 1 : 0);
        mesh_has_second_order_elems = IBTK_MPI::maxReduction(mesh_has_second_order_elems ? 1 : 0);
        if ((mesh_has_first_order_elems && mesh_has_second_order_elems) ||
            (!mesh_has_first_order_elems && !mesh_has_second_order_elems))
        {
            TBOX_ERROR(d_object_name << "::FEMechanicsBase():\n"
                                     << "  each FE mesh part must contain only FIRST "
                                        "order elements or only SECOND order elements"
                                     << std::endl);
        }
        d_fe_family_position[part] = LAGRANGE;
        d_fe_family_force[part] = LAGRANGE;
        d_fe_family_pressure[part] = LAGRANGE;
        d_default_quad_type_stress[part] = QGAUSS;
        d_default_quad_type_force[part] = QGAUSS;
        d_default_quad_type_pressure[part] = QGAUSS;
        if (mesh_has_first_order_elems)
        {
            d_fe_order_position[part] = FIRST;
            d_fe_order_force[part] = FIRST;
            d_fe_order_pressure[part] = FIRST;
            d_default_quad_order_stress[part] = THIRD;
            d_default_quad_order_force[part] = THIRD;
            d_default_quad_order_pressure[part] = THIRD;
        }
        if (mesh_has_second_order_elems)
        {
            d_fe_order_position[part] = SECOND;
            d_fe_order_force[part] = SECOND;
            d_fe_order_pressure[part] = SECOND;
            d_default_quad_order_stress[part] = FIFTH;
            d_default_quad_order_force[part] = FIFTH;
            d_default_quad_order_pressure[part] = FIFTH;
        }
    }

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Report configuration for each part.
    for (unsigned int part = 0; part < n_parts; ++part)
    {
        pout << "\n";
        pout << d_object_name << ": mesh part " << part << " is using ";
        if ((d_fe_family_position[part] == d_fe_family_force[part]) &&
            (d_fe_family_force[part] == d_fe_family_pressure[part]) &&
            (d_fe_order_position[part] == d_fe_order_force[part]) &&
            (d_fe_order_force[part] == d_fe_order_pressure[part]))
        {
            pout << Utility::enum_to_string<Order>(d_fe_order_position[part]) << " order "
                 << Utility::enum_to_string<FEFamily>(d_fe_family_position[part]) << " finite elements.\n";
        }
        else
        {
            pout << "  " << Utility::enum_to_string<Order>(d_fe_order_position[part]) << " order "
                 << Utility::enum_to_string<FEFamily>(d_fe_family_position[part]) << " elements for positions\n"
                 << "  " << Utility::enum_to_string<Order>(d_fe_order_force[part]) << " order "
                 << Utility::enum_to_string<FEFamily>(d_fe_family_pressure[part]) << " elements for forces\n"
                 << "  " << Utility::enum_to_string<Order>(d_fe_order_pressure[part]) << " order "
                 << Utility::enum_to_string<FEFamily>(d_fe_family_pressure[part]) << " elements for pressures\n";
        }
        pout << "\n";
    }
}

void
FEMechanicsBase::getFromInput(const Pointer<Database>& db, bool /*is_from_restart*/)
{
    // libMesh parallelization settings.
    d_libmesh_partitioner_type =
        string_to_enum<LibmeshPartitionerType>(db->getStringWithDefault("libmesh_partitioner_type", "LIBMESH_DEFAULT"));

    // Force computation settings.
    if (db->isBool("use_consistent_mass_matrix"))
        d_use_consistent_mass_matrix = db->getBool("use_consistent_mass_matrix");

    // Pressure settings.
    if (db->isDouble("static_pressure_kappa")) d_static_pressure_kappa = db->getDouble("static_pressure_kappa");
    if (db->isDouble("static_pressure_stab_param"))
        d_static_pressure_stab_param = db->getDouble("static_pressure_stab_param");

    if (db->isDouble("dynamic_pressure_kappa")) d_dynamic_pressure_kappa = db->getDouble("dynamic_pressure_kappa");
    if (db->isDouble("dynamic_pressure_stab_param"))
        d_dynamic_pressure_stab_param = db->getDouble("dynamic_pressure_stab_param");

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
    if (db->keyExists("do_log"))
        d_do_log = db->getBool("do_log");
    else if (db->keyExists("enable_logging"))
        d_do_log = db->getBool("enable_logging");
}

void
FEMechanicsBase::getFromRestart()
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
    int ver = db->getInteger("FE_MECHANICS_BASE_VERSION");
    if (ver != FE_MECHANICS_BASE_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    d_use_consistent_mass_matrix = db->getBool("d_use_consistent_mass_matrix");
    d_libmesh_partitioner_type = string_to_enum<LibmeshPartitionerType>(db->getString("d_libmesh_partitioner_type"));
}


void
FEMechanicsBase::add_static_pressure_kappa_vector_map_entry(double kappa, unsigned int blockID, int part)
{
    int vec_size = d_static_pressure_kappa_vector_map.size();
    // if we have not created the map for this part
    // then create it
    if(vec_size == part )
    {
      d_static_pressure_kappa_vector_map.push_back( std::map<unsigned int, double>() );
    }
    // insert new entry
    if( vec_size >= part )
    {
      d_static_pressure_kappa_vector_map[part].insert(std::pair<unsigned int, double>(blockID, kappa));
    }
    // throw an error
    else
    {
      TBOX_ERROR(d_object_name << " adding kappa entry to part " << part << ", but the vector size is only " << vec_size << std::endl);
    }
}

void
FEMechanicsBase::add_static_pressure_stab_param_vector_map_entry(double tau, unsigned int blockID, int part)
{
    int vec_size = d_static_pressure_stab_param_vector_map.size();
    // if we have not created the map for this part
    // then create it
    if(vec_size == part )
    {
      d_static_pressure_stab_param_vector_map.push_back( std::map<unsigned int, double>() );
    }
    // insert new entry
    if( vec_size >= part )
    {
      d_static_pressure_stab_param_vector_map[part].insert(std::pair<unsigned int, double>(blockID, tau));
    }
    // throw an error
    else
    {
      TBOX_ERROR(d_object_name << " adding kappa entry to part " << part << ", but the vector size is only " << vec_size << std::endl);
    }
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
