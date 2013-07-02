// Filename: IBFEPatchRecoveryPostProcessor.h
// Created on 2 Jul 2013 by Boyce Griffith
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

#ifndef included_IBFEPatchRecoveryPostProcessor
#define included_IBFEPatchRecoveryPostProcessor

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/FEDataManager.h"
#include "libmesh/mesh.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/point.h"
#include "libmesh/system.h"
#include "libmesh/vector_value.h"
#include "petscsys.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBFEPatchRecoveryPostProcessor uses local least squares recovery
 * techniques to evaluate stresses on the nodes of the FE mesh.
 */
class IBFEPatchRecoveryPostProcessor
{
public:
    /*!
     * Constructor.
     */
    IBFEPatchRecoveryPostProcessor(
        libMesh::MeshBase* mesh,
        IBTK::FEDataManager* fe_data_manager);

    /*!
     * Destructor.
     */
    ~IBFEPatchRecoveryPostProcessor();

    /*!
     * Register the PK1 stress associated with an element and quadrature point.
     */
    void
    registerPK1StressValue(
        const libMesh::Elem* elem,
        const libMesh::QBase* qrule,
        unsigned int qp,
        const libMesh::TensorValue<double>& PP,
        const libMesh::TensorValue<double>& FF);

    /*!
     * Reconstruct the Cauchy stress at the nodes of the mesh.
     */
    void
    reconstructCauchyStress();

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBFEPatchRecoveryPostProcessor();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBFEPatchRecoveryPostProcessor(
        const IBFEPatchRecoveryPostProcessor& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBFEPatchRecoveryPostProcessor&
    operator=(
        const IBFEPatchRecoveryPostProcessor& that);

    /*
     * FE data associated with this object.
     */
    libMesh::MeshBase* d_mesh;
    IBTK::FEDataManager* d_fe_data_manager;

    /*
     * Map from local nodes to element patches.
     */
    typedef std::set<const Elem*> ElemPatch;
    std::map<libMesh::dof_id_type,ElemPatch> d_local_elem_patches;
    std::vector<const Elem*> d_local_elems;

    /*
     * Interpolation point indexing data for each element.
     */
    int d_n_qp_global, d_n_qp_local, d_qp_global_offset;
    std::vector<int> d_elem_n_qp, d_elem_qp_global_offset, d_elem_qp_local_offset;

    /*
     * Stress data at interpolation points.
     */
    typedef std::vector<libMesh::TensorValue<double> > ElemStress;
    std::map<libMesh::dof_id_type,ElemStress> d_elem_sigma;
};
}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBFEPatchRecoveryPostProcessor
