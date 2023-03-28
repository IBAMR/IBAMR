// ---------------------------------------------------------------------
//
// Copyright (c) 2015 - 2021 by the IBAMR developers
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

#ifndef included_FEDataInterpolation
#define included_FEDataInterpolation

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#ifdef IBTK_HAVE_LIBMESH

#include "ibtk/FEDataManager.h"
#include "ibtk/FEValues.h"
#include "ibtk/libmesh_utilities.h"

#include "tbox/Utilities.h"

#include "libmesh/equation_systems.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/fem_context.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/vector_value.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <limits>
#include <memory>
#include <vector>

namespace IBTK
{
struct SystemData;
} // namespace IBTK

namespace libMesh
{
class Elem;
class EquationSystems;
class Point;
class QBase;
class System;
template <typename T>
class NumericVector;
} // namespace libMesh

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class FEDataInterpolation manages data required to evaluate one or more FE field variables at a collection of
 * points, possibly (not not necessarily) corresponding to the points of a quadrature rule.
 */
class FEDataInterpolation
{
public:
    FEDataInterpolation(unsigned int dim, std::shared_ptr<FEData> fe_data);

    ~FEDataInterpolation() = default;

    inline void attachQuadratureRule(libMesh::QBase* qrule)
    {
        d_qrule = qrule;
        return;
    }

    inline void attachQuadratureRuleFace(libMesh::QBase* qrule_face)
    {
        d_qrule_face = qrule_face;
        for (auto& fe_face : d_fe_face)
        {
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
                                      libMesh::NumericVector<double>* system_vec = nullptr);

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
     * \brief Initialize all of the data structures required to evaluate the FE shape functions, quadrature rules, etc.
     *
     * NOTE: This method must be called before reinitializing data on individual elements.
     */
    void init();

    /*!
     * \brief Reinitialize the FE shape functions, quadrature rules, etc. for the specified element.
     *
     * NOTE: Nodal values are set by calling collectDataForInterpolation().
     */
    void reinit(const libMesh::Elem* elem,
                const std::vector<libMesh::Point>* const points = nullptr,
                const std::vector<double>* weights = nullptr);

    /*!
     * \brief Reinitialize the FE shape functions, quadrature rules, etc. for the specified side of the specified
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
    std::shared_ptr<FEData> d_fe_data;
    bool d_initialized = false;
    bool d_eval_q_point = false, d_eval_JxW = false, d_eval_q_point_face = false, d_eval_JxW_face = false,
         d_eval_normal_face = false;
    libMesh::QBase *d_qrule = nullptr, *d_qrule_face = nullptr;
    const std::vector<libMesh::Point>*d_q_point = nullptr, *d_q_point_face = nullptr;
    const std::vector<double>*d_JxW = nullptr, *d_JxW_face = nullptr;
    const std::vector<libMesh::Point>* d_normal_face = nullptr;

    // Data associated with systems.
    std::vector<const libMesh::System*> d_systems;
    std::vector<FEDataManager::SystemDofMapCache*> d_system_dof_map_caches;
    std::vector<std::vector<int> > d_system_all_vars, d_system_vars, d_system_grad_vars;
    std::vector<std::vector<size_t> > d_system_var_idx, d_system_grad_var_idx;
    std::vector<libMesh::NumericVector<double>*> d_system_vecs;
    std::vector<std::vector<size_t> > d_system_var_fe_type_idx, d_system_grad_var_fe_type_idx;
    std::vector<std::vector<std::vector<double> > > d_system_var_data;
    std::vector<std::vector<std::vector<libMesh::VectorValue<double> > > > d_system_grad_var_data;
    std::vector<const libMesh::System*> d_noninterp_systems;
    std::vector<std::vector<int> > d_noninterp_system_all_vars, d_noninterp_system_phi_vars,
        d_noninterp_system_dphi_vars;

    // Data associated with FETypes.
    std::vector<libMesh::FEType> d_fe_types;
    std::vector<std::unique_ptr<IBTK::FEValuesBase> > d_fe;
    std::vector<std::unique_ptr<libMesh::FEBase> > d_fe_face;
    std::vector<bool> d_eval_phi, d_eval_dphi;
    std::vector<const std::vector<std::vector<double> >*> d_phi, d_phi_face;
    std::vector<const std::vector<std::vector<libMesh::VectorValue<double> > >*> d_dphi, d_dphi_face;

    // Data associated with the current element.
    const libMesh::Elem* d_current_elem = nullptr;
    unsigned int d_current_side = std::numeric_limits<unsigned int>::max();
    std::vector<boost::multi_array<double, 2> > d_system_elem_data;
    unsigned int d_n_qp = std::numeric_limits<unsigned int>::max();
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifdef IBTK_HAVE_LIBMESH
#endif //#ifndef included_FEDataInterpolation
