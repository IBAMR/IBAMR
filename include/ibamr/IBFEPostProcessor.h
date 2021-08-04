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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_IBFEPostProcessor
#define included_IBAMR_IBFEPostProcessor

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#ifdef IBAMR_HAVE_LIBMESH

#include "ibamr/IBFEMethod.h"

#include "ibtk/FEDataManager.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/libmesh_utilities.h"

#include "Variable.h"
#include "VariableContext.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include "libmesh/compare_types.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/mesh.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/point.h"
#include "libmesh/system.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_tensor.h"
#include "libmesh/type_vector.h"
#include "libmesh/vector_value.h"

#include "petscsys.h"

#include <string>
#include <vector>

namespace libMesh
{
class Elem;
class MeshBase;
class System;
} // namespace libMesh

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBFEPostProcessor is a generic interface for specifying the
 * implementation details of a particular post processing algorithm for the
 * IB/FE scheme.
 */
class IBFEPostProcessor
{
public:
    /*!
     * \brief Function for reconstructing the deformation gradient tensor FF =
     * dX/ds.
     */
    static inline void
    FF_fcn(libMesh::TensorValue<double>& FF_out,
           const libMesh::TensorValue<double>& FF_in,
           const libMesh::Point& /*X*/,
           const libMesh::Point& /*s*/,
           libMesh::Elem* /*elem*/,
           const std::vector<const std::vector<double>*>& /*system_var_data*/,
           const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
           double /*data_time*/,
           void* /*ctx*/)
    {
        FF_out = FF_in;
        return;
    } // FF_fcn

    /*!
     * \brief Function for reconstructing the Green-Lagrangian strain tensor EE
     * = 0.5*(CC - II), with CC = FF^T FF and FF = dX/ds.
     */
    static inline void
    EE_fcn(libMesh::TensorValue<double>& EE,
           const libMesh::TensorValue<double>& FF,
           const libMesh::Point& /*X*/,
           const libMesh::Point& /*s*/,
           libMesh::Elem* /*elem*/,
           const std::vector<const std::vector<double>*>& /*system_var_data*/,
           const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
           double /*data_time*/,
           void* /*ctx*/)
    {
        const libMesh::TensorValue<double> CC = FF.transpose() * FF;
        static const libMesh::TensorValue<double> II(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        EE = 0.5 * (CC - II);
        return;
    } // EE_fcn

    /*!
     * \brief Function for reconstructing the Cauchy stress from the PK1 stress,
     * using the PK1 stress function data provided by the ctx argument.
     */
    static inline void cauchy_stress_from_PK1_stress_fcn(
        libMesh::TensorValue<double>& sigma,
        const libMesh::TensorValue<double>& FF,
        const libMesh::Point& X,
        const libMesh::Point& s,
        libMesh::Elem* elem,
        const std::vector<const std::vector<double>*>& system_var_data,
        const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
        double data_time,
        void* ctx)
    {
        TBOX_ASSERT(ctx);
        auto PK1_stress_fcn_data = static_cast<IBFEMethod::PK1StressFcnData*>(ctx);
        TBOX_ASSERT(PK1_stress_fcn_data);
        IBTK::TensorMeshFcnPtr PK1_stress_fcn = PK1_stress_fcn_data->fcn;
        void* PK1_stress_fcn_ctx = PK1_stress_fcn_data->ctx;
        libMesh::TensorValue<double> PP;
        PK1_stress_fcn(PP, FF, X, s, elem, system_var_data, system_grad_var_data, data_time, PK1_stress_fcn_ctx);
        sigma = PP * FF.transpose() / FF.det();
        return;
    } // cauchy_stress_from_PK1_stress_fcn

    /*!
     * \brief Function for reconstructing a deformed material axis.  A pointer
     * to the system number must be passed as the ctx argument.
     */
    static inline void deformed_material_axis_fcn(
        libMesh::VectorValue<double>& f,
        const libMesh::TensorValue<double>& FF,
        const libMesh::Point& /*X*/,
        const libMesh::Point& /*s*/,
        libMesh::Elem* /*elem*/,
        const std::vector<const std::vector<double>*>& system_var_data,
        const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
        double /*data_time*/,
        void* /*ctx*/)
    {
        TBOX_ASSERT(system_var_data.size() == 1);
        libMesh::VectorValue<double> f0;
        for (unsigned int d = 0; d < NDIM; ++d) f0(d) = (*system_var_data[0])[d];
        f = FF * f0;
        return;
    } // deformed_material_axis_fcn

    /*!
     * \brief Function for reconstructing a deformed, normalized material axis.
     * A pointer to the system number must be passed as the ctx argument.
     */
    static inline void deformed_normalized_material_axis_fcn(
        libMesh::VectorValue<double>& f,
        const libMesh::TensorValue<double>& FF,
        const libMesh::Point& /*X*/,
        const libMesh::Point& /*s*/,
        libMesh::Elem* /*elem*/,
        const std::vector<const std::vector<double>*>& system_var_data,
        const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
        double /*data_time*/,
        void* /*ctx*/)
    {
        TBOX_ASSERT(system_var_data.size() == 1);
        libMesh::VectorValue<double> f0;
        for (unsigned int d = 0; d < NDIM; ++d) f0(d) = (*system_var_data[0])[d];
        f = (FF * f0).unit();
        return;
    } // deformed_normalized_material_axis_fcn

    /*!
     * \brief Function for reconstructing the stretch in a material axis.  A
     * pointer to the system number must be passed as the ctx argument.
     */
    static inline void material_axis_stretch_fcn(
        double& lambda,
        const libMesh::TensorValue<double>& FF,
        const libMesh::Point& /*X*/,
        const libMesh::Point& /*s*/,
        libMesh::Elem* /*elem*/,
        const std::vector<const std::vector<double>*>& system_var_data,
        const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
        double /*data_time*/,
        void* /*ctx*/)
    {
        TBOX_ASSERT(system_var_data.size() == 1);
        libMesh::VectorValue<double> f0;
        for (unsigned int d = 0; d < NDIM; ++d) f0(d) = (*system_var_data[0])[d];
        const libMesh::VectorValue<double> f = FF * f0;
        lambda = f.norm() / f0.norm();
        return;
    } // material_axis_stretch_fcn

    /*!
     * Constructor.
     */
    IBFEPostProcessor(std::string name, IBTK::FEDataManager* fe_data_manager);

    /*!
     * Virtual destructor.
     */
    virtual ~IBFEPostProcessor() = default;

    /*!
     * Register a scalar-valued variable for reconstruction.
     */
    virtual void
    registerScalarVariable(const std::string& name,
                           libMesh::FEFamily fe_family,
                           libMesh::Order fe_order,
                           IBTK::ScalarMeshFcnPtr fcn,
                           const std::vector<IBTK::SystemData>& system_data = std::vector<IBTK::SystemData>(),
                           void* var_fcn_ctx = nullptr);

    /*!
     * Register a vector-valued variable for reconstruction.
     */
    virtual void
    registerVectorVariable(const std::string& var_name,
                           libMesh::FEFamily var_fe_family,
                           libMesh::Order var_fe_order,
                           IBTK::VectorMeshFcnPtr var_fcn,
                           const std::vector<IBTK::SystemData>& system_data = std::vector<IBTK::SystemData>(),
                           void* var_fcn_ctx = nullptr,
                           unsigned int var_dim = NDIM);

    /*!
     * Register a tensor-valued variable for reconstruction.
     */
    virtual void
    registerTensorVariable(const std::string& var_name,
                           libMesh::FEFamily var_fe_family,
                           libMesh::Order var_fe_order,
                           IBTK::TensorMeshFcnPtr var_fcn,
                           const std::vector<IBTK::SystemData>& system_data = std::vector<IBTK::SystemData>(),
                           void* var_fcn_ctx = nullptr,
                           unsigned int var_dim = NDIM);

    /*!
     * Register a scalar-valued Eulerian field for reconstruction on the FE
     * mesh.  The variable is interpolated using the default interp spec
     * provided by the associated FEDataManager object.
     */
    virtual void registerInterpolatedScalarEulerianVariable(
        const std::string& var_name,
        libMesh::FEFamily var_fe_family,
        libMesh::Order var_fe_order,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx,
        const IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent& ghost_fill_transaction);

    /*!
     * Register a scalar-valued Eulerian field for reconstruction on the FE
     * mesh.  The variable is interpolated using the specified interp spec.
     */
    virtual void registerInterpolatedScalarEulerianVariable(
        const std::string& var_name,
        libMesh::FEFamily var_fe_family,
        libMesh::Order var_fe_order,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx,
        const IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent& ghost_fill_transaction,
        const IBTK::FEDataManager::InterpSpec& interp_spec);

    /*!
     * Initialize data used by the post processor.
     */
    virtual void initializeFEData();

    /*!
     * Execute all reconstruction and interpolation operations.
     */
    virtual void postProcessData(double data_time);

protected:
    /*!
     * Virtual function to interpolate Eulerian data to the mesh.
     */
    virtual void interpolateVariables(double data_time);

    /*!
     * Pure virtual function to reconstruct the data on the mesh.
     */
    virtual void reconstructVariables(double data_time) = 0;

    /*!
     * Name of the post processor object (used for internal variable context).
     */
    const std::string d_name;

    /*!
     * Mesh data.
     */
    libMesh::MeshBase* d_mesh;
    IBTK::FEDataManager* d_fe_data_manager;
    bool d_fe_data_initialized = false;

    /*!
     * Scalar-valued reconstruction data.
     */
    std::vector<libMesh::System*> d_scalar_var_systems;
    std::vector<IBTK::ScalarMeshFcnPtr> d_scalar_var_fcns;
    std::vector<std::vector<IBTK::SystemData> > d_scalar_var_system_data;
    std::vector<void*> d_scalar_var_fcn_ctxs;

    /*!
     * Vector-valued reconstruction data.
     */
    std::vector<libMesh::System*> d_vector_var_systems;
    std::vector<IBTK::VectorMeshFcnPtr> d_vector_var_fcns;
    std::vector<std::vector<IBTK::SystemData> > d_vector_var_system_data;
    std::vector<void*> d_vector_var_fcn_ctxs;
    std::vector<unsigned int> d_vector_var_dims;

    /*!
     * Tensor-valued reconstruction data.
     */
    std::vector<libMesh::System*> d_tensor_var_systems;
    std::vector<IBTK::TensorMeshFcnPtr> d_tensor_var_fcns;
    std::vector<std::vector<IBTK::SystemData> > d_tensor_var_system_data;
    std::vector<void*> d_tensor_var_fcn_ctxs;
    std::vector<unsigned int> d_tensor_var_dims;

    /*!
     * Eulerian interpolation data.
     */
    std::vector<libMesh::System*> d_scalar_interp_var_systems;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_scalar_interp_vars;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> > d_scalar_interp_ctxs;
    std::vector<int> d_scalar_interp_data_idxs, d_scalar_interp_scratch_idxs;
    std::vector<IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent>
        d_scalar_interp_fill_transactions;
    std::vector<IBTK::FEDataManager::InterpSpec> d_scalar_interp_specs;

    /*!
     * Collection of all systems managed by this object.
     */
    std::vector<libMesh::System*> d_var_systems;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBFEPostProcessor() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBFEPostProcessor(const IBFEPostProcessor& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBFEPostProcessor& operator=(const IBFEPostProcessor& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifdef IBAMR_HAVE_LIBMESH
#endif //#ifndef included_IBAMR_IBFEPostProcessor
