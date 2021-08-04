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

#ifndef included_IBAMR_IBFECentroidPostProcessor
#define included_IBAMR_IBFECentroidPostProcessor

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#ifdef IBAMR_HAVE_LIBMESH

#include "ibamr/IBFEPostProcessor.h"

#include "ibtk/libmesh_utilities.h"

#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"

#include <string>
#include <vector>

namespace IBTK
{
class FEDataManager;
} // namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBFECentroidPostProcessor is a simple post processor that
 * reconstructs piecewise constant field data by direct evaluation at element
 * centroids.
 */
class IBFECentroidPostProcessor : public IBFEPostProcessor
{
public:
    /*!
     * Constructor.
     */
    IBFECentroidPostProcessor(std::string name, IBTK::FEDataManager* fe_data_manager);

    /*!
     * Register a scalar-valued variable for reconstruction.
     *
     * \note This method checks that the requested FE family and order are valid
     * and then calls the corresponding method in the IBFEPostProcessor base
     * class.
     */
    void registerScalarVariable(const std::string& name,
                                libMesh::FEFamily fe_family,
                                libMesh::Order fe_order,
                                IBTK::ScalarMeshFcnPtr fcn,
                                const std::vector<IBTK::SystemData>& system_data = std::vector<IBTK::SystemData>(),
                                void* fcn_ctx = nullptr) override;

    /*!
     * Register a vector-valued variable for reconstruction.
     *
     * \note This method checks that the requested FE family and order are valid
     * and then calls the corresponding method in the IBFEPostProcessor base
     * class.
     */
    void registerVectorVariable(const std::string& var_name,
                                libMesh::FEFamily var_fe_family,
                                libMesh::Order var_fe_order,
                                IBTK::VectorMeshFcnPtr var_fcn,
                                const std::vector<IBTK::SystemData>& system_data = std::vector<IBTK::SystemData>(),
                                void* var_fcn_ctx = nullptr,
                                unsigned int var_dim = NDIM) override;

    /*!
     * Register a tensor-valued variable for reconstruction.
     *
     * \note This method checks that the requested FE family and order are valid
     * and then calls the corresponding method in the IBFEPostProcessor base
     * class.
     */
    void registerTensorVariable(const std::string& var_name,
                                libMesh::FEFamily var_fe_family,
                                libMesh::Order var_fe_order,
                                IBTK::TensorMeshFcnPtr var_fcn,
                                const std::vector<IBTK::SystemData>& system_data = std::vector<IBTK::SystemData>(),
                                void* var_fcn_ctx = nullptr,
                                unsigned int var_dim = NDIM) override;

    /*!
     * Reconstruct the data on the mesh.
     */
    void reconstructVariables(double data_time) override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBFECentroidPostProcessor() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBFECentroidPostProcessor(const IBFECentroidPostProcessor& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBFECentroidPostProcessor& operator=(const IBFECentroidPostProcessor& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifdef IBAMR_HAVE_LIBMESH
#endif //#ifndef included_IBAMR_IBFECentroidPostProcessor
