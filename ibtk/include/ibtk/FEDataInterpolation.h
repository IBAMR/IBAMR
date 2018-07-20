// Filename: FEDataInterpolation.h
// Created on 9 Oct 2015 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#ifndef included_FEDataInterpolation
#define included_FEDataInterpolation

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "boost/multi_array.hpp"
#include "ibtk/FEDataManager.h"
#include "ibtk/libmesh_utilities.h"
#include "libmesh/equation_systems.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class FEDataInterpolation manages data requred to evaluate one or more FE field variables at a collection of
 * points, possibly (not not necessarily) corresponding to the points of a quadrature rule.
 */
class FEDataInterpolation
{
public:
    FEDataInterpolation(unsigned int dim, FEDataManager* const fe_data_manager);

    ~FEDataInterpolation();

    inline void attachQuadratureRule(libMesh::QBase* qrule)
    {
        d_qrule = qrule;
        for (size_t k = 0; k < d_fe.size(); ++k)
        {
            SAMRAI::tbox::Pointer<libMesh::FEBase> fe = d_fe[k];
            if (fe)
            {
                fe->attach_quadrature_rule(d_qrule);
            }
        }
        return;
    }

    inline void attachQuadratureRuleFace(libMesh::QBase* qrule_face)
    {
        d_qrule_face = qrule_face;
        for (size_t k = 0; k < d_fe_face.size(); ++k)
        {
            SAMRAI::tbox::Pointer<libMesh::FEBase> fe_face = d_fe_face[k];
            if (fe_face)
            {
                fe_face->attach_quadrature_rule(d_qrule_face);
            }
        }
        return;
    }

    inline void evalQuadraturePoints()
    {
        TBOX_ASSERT(!d_initialized);
        d_eval_q_point = true;
        return;
    }

    inline void evalQuadratureWeights()
    {
        TBOX_ASSERT(!d_initialized);
        d_eval_JxW = true;
        return;
    }

    inline void evalQuadraturePointsFace()
    {
        TBOX_ASSERT(!d_initialized);
        d_eval_q_point_face = true;
        return;
    }

    inline void evalQuadratureWeightsFace()
    {
        TBOX_ASSERT(!d_initialized);
        d_eval_JxW_face = true;
        return;
    }

    inline void evalNormalsFace()
    {
        TBOX_ASSERT(!d_initialized);
        d_eval_normal_face = true;
        return;
    }

    inline const libMesh::QBase* getQrule() const
    {
        return d_qrule;
    }

    inline const libMesh::QBase* getQruleFace() const
    {
        return d_qrule_face;
    }

    inline const std::vector<libMesh::Point>& getQuadraturePoints() const
    {
        TBOX_ASSERT(d_initialized);
        TBOX_ASSERT(d_eval_q_point);
        return *d_q_point;
    }

    inline const std::vector<double>& getQuadratureWeights() const
    {
        TBOX_ASSERT(d_initialized);
        TBOX_ASSERT(d_eval_JxW);
        return *d_JxW;
    }

    inline const std::vector<libMesh::Point>& getQuadraturePointsFace() const
    {
        TBOX_ASSERT(d_initialized);
        TBOX_ASSERT(d_eval_q_point_face);
        return *d_q_point_face;
    }

    inline const std::vector<double>& getQuadratureWeightsFace() const
    {
        TBOX_ASSERT(d_initialized);
        TBOX_ASSERT(d_eval_JxW_face);
        return *d_JxW_face;
    }

    inline const std::vector<libMesh::Point>& getNormalsFace() const
    {
        TBOX_ASSERT(d_initialized);
        TBOX_ASSERT(d_eval_normal_face);
        return *d_normal_face;
    }

    inline const std::vector<std::vector<double> >& getPhi(const libMesh::FEType& fe_type) const
    {
        TBOX_ASSERT(d_initialized);
        const size_t fe_type_idx = getFETypeIndex(fe_type);
        TBOX_ASSERT(fe_type_idx < d_fe_types.size());
        TBOX_ASSERT(d_eval_phi[fe_type_idx]);
        return *d_phi[fe_type_idx];
    }

    inline const std::vector<std::vector<libMesh::VectorValue<double> > >& getDphi(const libMesh::FEType& fe_type) const
    {
        TBOX_ASSERT(d_initialized);
        const size_t fe_type_idx = getFETypeIndex(fe_type);
        TBOX_ASSERT(fe_type_idx < d_fe_types.size());
        TBOX_ASSERT(d_eval_dphi[fe_type_idx]);
        return *d_dphi[fe_type_idx];
    }

    inline const std::vector<std::vector<double> >& getPhiFace(const libMesh::FEType& fe_type) const
    {
        TBOX_ASSERT(d_initialized);
        const size_t fe_type_idx = getFETypeIndex(fe_type);
        TBOX_ASSERT(fe_type_idx < d_fe_types.size());
        TBOX_ASSERT(d_eval_phi[fe_type_idx]);
        return *d_phi_face[fe_type_idx];
    }

    inline const std::vector<std::vector<libMesh::VectorValue<double> > >&
    getDphiFace(const libMesh::FEType& fe_type) const
    {
        TBOX_ASSERT(d_initialized);
        const size_t fe_type_idx = getFETypeIndex(fe_type);
        TBOX_ASSERT(fe_type_idx < d_fe_types.size());
        TBOX_ASSERT(d_eval_dphi[fe_type_idx]);
        return *d_dphi_face[fe_type_idx];
    }

    /*!
     * \brief Configure the class to evaluate the requested shape functions / derivatives for the given system.
     *
     * NOTE: To interpolate variables associated with the system, use registerInterpolatedSystem().
     */
    void registerSystem(const libMesh::System& system,
                        const std::vector<int>& phi_vars = std::vector<int>(1, 0),
                        const std::vector<int>& dphi_vars = std::vector<int>(1, 0));

    /*!
     * \brief Configure the class to interpolate the requested variables / gradients for the given system.
     *
     * If the system data vector is NULL, then this class will use system.current_local_solution.
     *
     * NOTE: The same system can be registered multiple times with different sets of variables/gradients and system data
     * vectors.
     *
     * \return Returns an index that is used to access the interpolated data.
     */
    size_t registerInterpolatedSystem(const libMesh::System& system,
                                      const std::vector<int>& vars = std::vector<int>(1, 0),
                                      const std::vector<int>& grad_vars = std::vector<int>(),
                                      libMesh::NumericVector<double>* system_data = nullptr);

    /*!
     * \brief Get the variable data for all of the systems.
     *
     * NOTE: Data are indexed via [qp][system_idx][var_idx].
     */
    inline const std::vector<std::vector<std::vector<double> > >& getVarInterpolation()
    {
        return d_system_var_data;
    }

    /*!
     * \brief Get the gradient variable data for all of the systems.
     *
     * NOTE: Data are indexed via [qp][system_idx][var_idx].
     */
    inline const std::vector<std::vector<std::vector<libMesh::VectorValue<double> > > >& getGradVarInterpolation()
    {
        return d_system_grad_var_data;
    }

    /*!
     * \brief Register systems to be interpolated and get the corresponding system indices.
     */
    void setupInterpolatedSystemDataIndexes(std::vector<size_t>& system_idxs,
                                            const std::vector<SystemData>& system_data,
                                            const libMesh::EquationSystems* const equation_systems);

    /*!
     * \brief Set up pointers to the interpolated data to be evaluated.
     */
    void setInterpolatedDataPointers(std::vector<const std::vector<double>*>& var_data,
                                     std::vector<const std::vector<libMesh::VectorValue<double> >*>& grad_var_data,
                                     const std::vector<size_t>& system_idxs,
                                     const libMesh::Elem* const elem,
                                     const unsigned int qp);

    /*!
     * \brief Initialize all of the data structures requred to evaluate the FE shape functions, quadrature rules, etc.
     *
     * NOTE: This method must be called before reinitializing data on individual elements.
     */
    void init(bool use_IB_ghosted_vecs);

    /*!
     * \brief Reinitialize the FE shape functions, quadrature rules, etc. for the specified element.
     *
     * NOTE: Nodal values are set by calling collectDataForInterpolation().
     */
    void reinit(const libMesh::Elem* elem,
                const std::vector<libMesh::Point>* const points = nullptr,
                const std::vector<double>* weights = nullptr);

    /*!
     * \brief Reinitialize the FE shape functions, quadrature rules, etc. for the specified side of the specificed
     * element.
     *
     * NOTE: Nodal values are set by calling collectDataForInterpolation().
     */
    void reinit(const libMesh::Elem* elem,
                unsigned int side,
                double tol = libMesh::TOLERANCE,
                const std::vector<libMesh::Point>* points = nullptr,
                const std::vector<double>* weights = nullptr);

    /*!
     * \brief Get the local (element) data to be interpolated from the global vectors.
     */
    void collectDataForInterpolation(const libMesh::Elem* elem);

    /*!
     * \brief Provide the elemental data associated with the given system index and element.
     */
    const boost::multi_array<double, 2>& getElemData(const libMesh::Elem* elem, size_t system_idx);

    /*!
     * \brief Interpolate FE data on the specified element.
     *
     * NOTE: Nodal values are set by calling collectDataForInterpolation().
     */
    void interpolate(const libMesh::Elem* elem);
    /*!
     * \brief Interpolate FE data on the specified side of the specified element.
     *
     * NOTE: Nodal values are set by calling collectDataForInterpolation().
     */
    void interpolate(const libMesh::Elem* elem, unsigned int side);

private:
    FEDataInterpolation() = delete;
    FEDataInterpolation(const FEDataInterpolation&) = delete;
    FEDataInterpolation& operator=(const FEDataInterpolation&) = delete;

    size_t getFETypeIndex(const libMesh::FEType& fe_type) const;

    void
    interpolateCommon(std::vector<std::vector<std::vector<double> > >& system_var_data,
                      std::vector<std::vector<std::vector<libMesh::VectorValue<double> > > >& system_grad_var_data,
                      const std::vector<const std::vector<std::vector<double> >*>& phi_data,
                      const std::vector<const std::vector<std::vector<libMesh::VectorValue<double> > >*>& dphi_data);

    const unsigned int d_dim;
    FEDataManager* const d_fe_data_manager;
    bool d_initialized;
    bool d_eval_q_point, d_eval_JxW, d_eval_q_point_face, d_eval_JxW_face, d_eval_normal_face;
    libMesh::QBase *d_qrule, *d_qrule_face;
    const std::vector<libMesh::Point> *d_q_point, *d_q_point_face;
    const std::vector<double> *d_JxW, *d_JxW_face;
    const std::vector<libMesh::Point>* d_normal_face;

    // Data associated with systems.
    std::vector<const libMesh::System*> d_systems;
    std::vector<FEDataManager::SystemDofMapCache*> d_system_dof_map_caches;
    std::vector<std::vector<int> > d_system_all_vars, d_system_vars, d_system_grad_vars;
    std::vector<std::vector<size_t> > d_system_var_idx, d_system_grad_var_idx;
    std::vector<libMesh::NumericVector<double>*> d_system_data;
    std::vector<std::vector<size_t> > d_system_var_fe_type_idx, d_system_grad_var_fe_type_idx;
    std::vector<std::vector<std::vector<double> > > d_system_var_data;
    std::vector<std::vector<std::vector<libMesh::VectorValue<double> > > > d_system_grad_var_data;
    std::vector<const libMesh::System*> d_noninterp_systems;
    std::vector<std::vector<int> > d_noninterp_system_all_vars, d_noninterp_system_phi_vars,
        d_noninterp_system_dphi_vars;

    // Data associated with FETypes.
    std::vector<libMesh::FEType> d_fe_types;
    std::vector<SAMRAI::tbox::Pointer<libMesh::FEBase> > d_fe, d_fe_face;
    std::vector<bool> d_eval_phi, d_eval_dphi;
    std::vector<const std::vector<std::vector<double> > *> d_phi, d_phi_face;
    std::vector<const std::vector<std::vector<libMesh::VectorValue<double> > > *> d_dphi, d_dphi_face;

    // Data associated with the current element.
    const libMesh::Elem* d_current_elem;
    unsigned int d_current_side;
    std::vector<boost::multi_array<double, 2> > d_system_elem_data;
    unsigned int d_n_qp;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_FEDataInterpolation
