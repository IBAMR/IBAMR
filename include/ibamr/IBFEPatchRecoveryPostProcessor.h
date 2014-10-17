// Filename: IBFEPatchRecoveryPostProcessor.h
// Created on 2 Jul 2013 by Boyce Griffith
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

#ifndef included_IBFEPatchRecoveryPostProcessor
#define included_IBFEPatchRecoveryPostProcessor

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "boost/tuple/tuple.hpp"
#include "ibtk/FEDataManager.h"
#include "libmesh/mesh.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/point.h"
#include "libmesh/system.h"
#include "libmesh/vector_value.h"
#include "libmesh/periodic_boundary.h"
#include "petscsys.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBFEPatchRecoveryPostProcessor uses least-squares
 * reconstructions on element patches to evaluate stresses at the nodes of the
 * FE mesh.
 */
class IBFEPatchRecoveryPostProcessor
{
public:
    /*!
     * Constructor.
     */
    IBFEPatchRecoveryPostProcessor(libMesh::MeshBase* mesh, IBTK::FEDataManager* fe_data_manager);

    /*!
     * Destructor.
     */
    ~IBFEPatchRecoveryPostProcessor();

    /*!
     * Initialize data used by the post processor.
     */
    void initializeFEData(const libMesh::PeriodicBoundaries* periodic_boundaries = NULL);

    /*!
     * Initialize system to store reconstructed Cauchy stress values.
     */
    libMesh::System* initializeCauchyStressSystem();

    /*!
     * Initialize system to store reconstructed pressure values.
     */
    libMesh::System* initializePressureSystem();

    /*!
     * Register the Cauchy stress associated with an element and quadrature
     * point.
     */
    void registerCauchyStressValue(const libMesh::Elem* elem,
                                   const libMesh::QBase* qrule,
                                   unsigned int qp,
                                   const libMesh::TensorValue<double>& sigma);

    /*!
     * Register the pressure associated with an element and quadrature point.
     */
    void registerPressureValue(const libMesh::Elem* elem, const libMesh::QBase* qrule, unsigned int qp, double p);

    /*!
     * Reconstruct the Cauchy stress at the nodes of the mesh.
     */
    void reconstructCauchyStress(libMesh::System& sigma_system);

    /*!
     * Reconstruct the pressure at the nodes of the mesh.
     */
    void reconstructPressure(libMesh::System& p_system);

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
    IBFEPatchRecoveryPostProcessor(const IBFEPatchRecoveryPostProcessor& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBFEPatchRecoveryPostProcessor& operator=(const IBFEPatchRecoveryPostProcessor& that);

    /*
     * FE data associated with this object.
     */
    libMesh::MeshBase* d_mesh;
    IBTK::FEDataManager* d_fe_data_manager;
    const PeriodicBoundaries* d_periodic_boundaries;
    libMeshEnums::Order d_interp_order, d_quad_order;

    /*
     * Map from local nodes to element patches.
     */
    typedef std::vector<SAMRAI::tbox::Pointer<libMesh::PeriodicBoundaryBase> > CompositePeriodicMapping;
    typedef boost::tuple<const Elem*, CompositePeriodicMapping, CompositePeriodicMapping> ElemPatchItem;
    struct ElemPatchItemComp : std::binary_function<const ElemPatchItem&, const ElemPatchItem&, bool>
    {
        inline bool operator()(const ElemPatchItem& lhs, const ElemPatchItem& rhs)
        {
            return lhs.get<0>() < rhs.get<0>();
        }
    };
    typedef std::set<ElemPatchItem, ElemPatchItemComp> ElemPatch;
    std::map<libMesh::dof_id_type, ElemPatch> d_local_elem_patches;

    static inline libMesh::Point apply_composite_periodic_mapping(const CompositePeriodicMapping& mapping,
                                                                  const libMesh::Point& p)
    {
        if (mapping.empty()) return p;
        libMesh::Point periodic_image = p;
        for (unsigned int k = 0; k < mapping.size(); ++k)
        {
            periodic_image = mapping[k]->get_corresponding_pos(periodic_image);
        }
        return periodic_image;
    } // apply_composite_periodic_mapping

    /*
     * Interpolation point indexing data for each element.
     */
    unsigned int d_n_qp_global, d_n_qp_local, d_qp_global_offset;
    std::vector<unsigned int> d_elem_n_qp, d_elem_qp_global_offset, d_elem_qp_local_offset;

    /*
     * Element patch L2 projection matrices.
     */
    std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd> > d_local_patch_proj_solver;

    /*
     * Stress data at interpolation points.
     */
    typedef std::vector<libMesh::TensorValue<double> > ElemStress;
    std::map<libMesh::dof_id_type, ElemStress> d_elem_sigma;
    typedef std::vector<double> ElemPressure;
    std::map<libMesh::dof_id_type, ElemPressure> d_elem_pressure;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBFEPatchRecoveryPostProcessor
