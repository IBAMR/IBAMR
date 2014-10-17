// Filename: IBFECentroidPostProcessor.h
// Created on 4 Dec 2013 by Boyce Griffith
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

#ifndef included_IBFECentroidPostProcessor
#define included_IBFECentroidPostProcessor

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/IBFEPostProcessor.h"

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
    IBFECentroidPostProcessor(const std::string& name, IBTK::FEDataManager* fe_data_manager);

    /*!
     * Destructor.
     */
    ~IBFECentroidPostProcessor();

    /*!
     * Register a scalar-valued variable for reconstruction.
     *
     * \note This method checks that the requested FE family and order are valid
     * and then calls the corresponding method in the IBFEPostProcessor base
     * class.
     */
    void registerScalarVariable(const std::string& var_name,
                                libMeshEnums::FEFamily var_fe_family,
                                libMeshEnums::Order var_fe_order,
                                IBTK::ScalarMeshFcnPtr var_fcn,
                                std::vector<unsigned int> var_fcn_systems = std::vector<unsigned int>(),
                                void* var_fcn_ctx = NULL);

    /*!
     * Register a vector-valued variable for reconstruction.
     *
     * \note This method checks that the requested FE family and order are valid
     * and then calls the corresponding method in the IBFEPostProcessor base
     * class.
     */
    void registerVectorVariable(const std::string& var_name,
                                libMeshEnums::FEFamily var_fe_family,
                                libMeshEnums::Order var_fe_order,
                                IBTK::VectorMeshFcnPtr var_fcn,
                                std::vector<unsigned int> var_fcn_systems = std::vector<unsigned int>(),
                                void* var_fcn_ctx = NULL,
                                unsigned int var_dim = NDIM);

    /*!
     * Register a tensor-valued variable for reconstruction.
     *
     * \note This method checks that the requested FE family and order are valid
     * and then calls the corresponding method in the IBFEPostProcessor base
     * class.
     */
    void registerTensorVariable(const std::string& var_name,
                                libMeshEnums::FEFamily var_fe_family,
                                libMeshEnums::Order var_fe_order,
                                IBTK::TensorMeshFcnPtr var_fcn,
                                std::vector<unsigned int> var_fcn_systems = std::vector<unsigned int>(),
                                void* var_fcn_ctx = NULL,
                                unsigned int var_dim = NDIM);

    /*!
     * Reconstruct the data on the mesh.
     */
    void reconstructVariables(double data_time);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBFECentroidPostProcessor();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBFECentroidPostProcessor(const IBFECentroidPostProcessor& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBFECentroidPostProcessor& operator=(const IBFECentroidPostProcessor& that);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBFECentroidPostProcessor
