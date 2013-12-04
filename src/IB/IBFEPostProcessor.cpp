// Filename: IBFEPostProcessor.cpp
// Created on 4 Dec 2013 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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
#include "Eigen/Dense"
#include "ibamr/IBFEPostProcessor.h"
#include "SAMRAI_config.h"
#include "boost/multi_array.hpp"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/libmesh_utilities.h"
#include "libmesh/boundary_info.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/mesh.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/quadrature.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/string_to_enum.h"

using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEPostProcessor::IBFEPostProcessor(
    MeshBase* mesh,
    FEDataManager* fe_data_manager)
    : d_mesh(mesh),
      d_fe_data_manager(fe_data_manager),
      d_fe_data_initialized(false)
{
    // intentionally blank
    return;
}// IBFEPostProcessor

IBFEPostProcessor::~IBFEPostProcessor()
{
    // intentionally blank
    return;
}// ~IBFEPostProcessor

void
IBFEPostProcessor::registerScalarVariable(
    const std::string& var_name,
    libMeshEnums::FEFamily var_fe_family,
    libMeshEnums::Order var_fe_order,
    ScalarMeshFcnPtr var_fcn,
    std::vector<unsigned int> var_fcn_systems,
    void* var_fcn_ctx)
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    System& system = equation_systems->add_system<System>(var_name + " reconstruction system");
    system.add_variable(var_name, var_fe_order, var_fe_family);
    d_scalar_var_systems.push_back(&system);
    d_scalar_var_fcns.push_back(var_fcn);
    d_scalar_var_fcn_systems.push_back(var_fcn_systems);
    d_scalar_var_fcn_ctxs.push_back(var_fcn_ctx);
    d_var_systems.push_back(&system);
    d_var_fcn_systems.insert(var_fcn_systems.begin(), var_fcn_systems.end());
    return;
}// registerScalarVariable

void
IBFEPostProcessor::registerVectorVariable(
    const std::string& var_name,
    libMeshEnums::FEFamily var_fe_family,
    libMeshEnums::Order var_fe_order,
    VectorMeshFcnPtr var_fcn,
    std::vector<unsigned int> var_fcn_systems,
    void* var_fcn_ctx,
    unsigned int var_dim)
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    System& system = equation_systems->add_system<System>(var_name + " reconstruction system");
    for (unsigned int i = 0; i < var_dim; ++i)
    {
        for (unsigned int j = 0; j < var_dim; ++j)
        {
            std::ostringstream os;
            os << var_name << "_" << i << j;
            system.add_variable(os.str(), var_fe_order, var_fe_family);
        }
    }
    d_vector_var_systems.push_back(&system);
    d_vector_var_fcns.push_back(var_fcn);
    d_vector_var_fcn_systems.push_back(var_fcn_systems);
    d_vector_var_fcn_ctxs.push_back(var_fcn_ctx);
    d_vector_var_dims.push_back(var_dim);
    d_var_systems.push_back(&system);
    d_var_fcn_systems.insert(var_fcn_systems.begin(), var_fcn_systems.end());
    return;
}// registerVectorVariable

void
IBFEPostProcessor::registerTensorVariable(
    const std::string& var_name,
    libMeshEnums::FEFamily var_fe_family,
    libMeshEnums::Order var_fe_order,
    TensorMeshFcnPtr var_fcn,
    std::vector<unsigned int> var_fcn_systems,
    void* var_fcn_ctx,
    unsigned int var_dim)
{
    EquationSystems* equation_systems = d_fe_data_manager->getEquationSystems();
    System& system = equation_systems->add_system<System>(var_name + " reconstruction system");
    for (unsigned int i = 0; i < var_dim; ++i)
    {
        for (unsigned int j = 0; j < var_dim; ++j)
        {
            std::ostringstream os;
            os << var_name << "_" << i << j;
            system.add_variable(os.str(), var_fe_order, var_fe_family);
        }
    }
    d_tensor_var_systems.push_back(&system);
    d_tensor_var_fcns.push_back(var_fcn);
    d_tensor_var_fcn_systems.push_back(var_fcn_systems);
    d_tensor_var_fcn_ctxs.push_back(var_fcn_ctx);
    d_tensor_var_dims.push_back(var_dim);
    d_var_systems.push_back(&system);
    d_var_fcn_systems.insert(var_fcn_systems.begin(), var_fcn_systems.end());
    return;
}// registerTensorVariable

void
IBFEPostProcessor::initializeFEData()
{
    if (d_fe_data_initialized) return;
    for (unsigned int k = 0; k < d_var_systems.size(); ++k)
    {
        System& system = *d_var_systems[k];
        system.assemble_before_solve = false;
        system.assemble();
    }
    d_fe_data_initialized = true;
    return;
}// initializeFEData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
