// ---------------------------------------------------------------------
//
// Copyright (c) 2015 - 2022 by the IBAMR developers
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

#include "ibtk/FECache.h"
#include "ibtk/FEDataInterpolation.h"
#include "ibtk/libmesh_utilities.h"

#include "libmesh/compare_types.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/quadrature.h"
#include "libmesh/system.h"
#include "libmesh/type_vector.h"

#include "ibtk/namespaces.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <algorithm>
#include <iterator>
#include <set>

namespace libMesh
{
class Elem;
class Point;
template <typename T>
class NumericVector;
} // namespace libMesh

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

FEDataInterpolation::FEDataInterpolation(const unsigned int dim, std::shared_ptr<FEData> fe_data)
    : d_dim(dim), d_fe_data(std::move(fe_data))
{
    return;
}

void
FEDataInterpolation::registerSystem(const System& system,
                                    const std::vector<int>& phi_vars,
                                    const std::vector<int>& dphi_vars)
{
    TBOX_ASSERT(!d_initialized && (!phi_vars.empty() || !dphi_vars.empty()));
    const unsigned int sys_num = system.number();
    for (auto noninterp_system : d_noninterp_systems)
    {
        if (noninterp_system->number() == sys_num)
        {
            // This system has already been registered.  If so, just merge the collection of variables.
            const size_t system_idx =
                std::distance(d_noninterp_systems.begin(),
                              std::find(d_noninterp_systems.begin(), d_noninterp_systems.end(), noninterp_system));

            std::vector<int>& orig_all_vars = d_noninterp_system_all_vars[system_idx];
            std::set<int> all_vars_set(orig_all_vars.begin(), orig_all_vars.end());
            all_vars_set.insert(phi_vars.begin(), phi_vars.end());
            all_vars_set.insert(dphi_vars.begin(), dphi_vars.end());
            d_noninterp_system_all_vars[system_idx].assign(all_vars_set.begin(), all_vars_set.end());

            std::vector<int>& orig_phi_vars = d_noninterp_system_phi_vars[system_idx];
            std::set<int> phi_vars_set(orig_phi_vars.begin(), orig_phi_vars.end());
            phi_vars_set.insert(phi_vars.begin(), phi_vars.end());
            d_noninterp_system_phi_vars[system_idx].assign(phi_vars_set.begin(), phi_vars_set.end());

            std::vector<int>& orig_dphi_vars = d_noninterp_system_dphi_vars[system_idx];
            std::set<int> dphi_vars_set(orig_dphi_vars.begin(), orig_dphi_vars.end());
            dphi_vars_set.insert(dphi_vars.begin(), dphi_vars.end());
            d_noninterp_system_dphi_vars[system_idx].assign(dphi_vars_set.begin(), dphi_vars_set.end());
            return;
        }
    }

    // We have not previously registered the system, so we need to register it here.
    d_noninterp_systems.push_back(&system);

    std::set<int> all_vars_set;
    all_vars_set.insert(phi_vars.begin(), phi_vars.end());
    all_vars_set.insert(dphi_vars.begin(), dphi_vars.end());
    d_noninterp_system_all_vars.push_back(std::vector<int>(all_vars_set.begin(), all_vars_set.end()));

    std::set<int> phi_vars_set(phi_vars.begin(), phi_vars.end());
    d_noninterp_system_phi_vars.push_back(std::vector<int>(phi_vars_set.begin(), phi_vars_set.end()));

    std::set<int> dphi_vars_set(dphi_vars.begin(), dphi_vars.end());
    d_noninterp_system_dphi_vars.push_back(std::vector<int>(dphi_vars_set.begin(), dphi_vars_set.end()));
    return;
}

size_t
FEDataInterpolation::registerInterpolatedSystem(const System& system,
                                                const std::vector<int>& vars,
                                                const std::vector<int>& grad_vars,
                                                NumericVector<double>* system_vec)
{
    TBOX_ASSERT(!d_initialized && (!vars.empty() || !grad_vars.empty()));
    const unsigned int sys_num = system.number();
    for (auto system : d_systems)
    {
        if (system->number() == sys_num)
        {
            // This system has already been registered.  Check to see if the same collection of variables (etc.)
            // were used in the previous registration action; if so, do not re-register the system.
            //
            // NOTE: The same system may be registered multiple times with different collections of variables,
            // gradient variables, system data, etc.  We want to make sure we check all of the registered systems to
            // see if any of them match the function arguments.
            const size_t system_idx =
                std::distance(d_systems.begin(), std::find(d_systems.begin(), d_systems.end(), system));
            bool same_data = true;
            same_data = same_data && (vars == d_system_vars[system_idx]);
            same_data = same_data && (grad_vars == d_system_grad_vars[system_idx]);
            same_data = same_data && (system_vec == d_system_vecs[system_idx]);
            if (same_data) return system_idx;
        }
    }

    // Either we have not previously registered this system, or we are registering it here with a different
    // collection of variables/data.  In either case, we need to register it here.
    const size_t system_idx = d_systems.size();
    d_systems.push_back(&system);
    d_system_dof_map_caches.push_back(d_fe_data->getDofMapCache(system.name()));
    std::set<int> all_vars_set;
    all_vars_set.insert(vars.begin(), vars.end());
    all_vars_set.insert(grad_vars.begin(), grad_vars.end());
    std::vector<int> all_vars(all_vars_set.begin(), all_vars_set.end());
    d_system_all_vars.push_back(all_vars);
    d_system_vars.push_back(vars);
    std::vector<size_t> var_idx(vars.size());
    for (size_t k = 0; k < vars.size(); ++k)
    {
        var_idx[k] = std::distance(all_vars.begin(), std::find(all_vars.begin(), all_vars.end(), vars[k]));
    }
    d_system_var_idx.push_back(var_idx);
    d_system_grad_vars.push_back(grad_vars);
    std::vector<size_t> grad_var_idx(grad_vars.size());
    for (size_t k = 0; k < grad_vars.size(); ++k)
    {
        grad_var_idx[k] = std::distance(all_vars.begin(), std::find(all_vars.begin(), all_vars.end(), grad_vars[k]));
    }
    d_system_grad_var_idx.push_back(grad_var_idx);
    d_system_vecs.push_back(system_vec);
    return system_idx;
}

void
FEDataInterpolation::setupInterpolatedSystemDataIndexes(std::vector<size_t>& system_idxs,
                                                        const std::vector<SystemData>& system_data,
                                                        const EquationSystems* const equation_systems)
{
    TBOX_ASSERT(!d_initialized);
    const size_t n_systems = system_data.size();
    system_idxs.resize(n_systems);
    for (unsigned int k = 0; k < system_idxs.size(); ++k)
    {
        const System& system = equation_systems->get_system(system_data[k].system_name);
        const std::vector<int>& vars = system_data[k].vars;
        const std::vector<int>& grad_vars = system_data[k].grad_vars;
        NumericVector<double>* const system_vec = system_data[k].system_vec;
        system_idxs[k] = registerInterpolatedSystem(system, vars, grad_vars, system_vec);
    }
    return;
}

void
FEDataInterpolation::setInterpolatedDataPointers(std::vector<const std::vector<double>*>& var_data,
                                                 std::vector<const std::vector<VectorValue<double> >*>& grad_var_data,
                                                 const std::vector<size_t>& system_idxs,
                                                 const Elem* const elem,
                                                 const unsigned int qp)
{
    TBOX_ASSERT(d_initialized);
    TBOX_ASSERT(elem == d_current_elem);
    TBOX_ASSERT(qp < d_n_qp);
    var_data.resize(system_idxs.size());
    grad_var_data.resize(system_idxs.size());
    for (unsigned k = 0; k < system_idxs.size(); ++k)
    {
        var_data[k] = &d_system_var_data[qp][system_idxs[k]];
        grad_var_data[k] = &d_system_grad_var_data[qp][system_idxs[k]];
    }
    return;
}

void
FEDataInterpolation::init()
{
    TBOX_ASSERT(!d_initialized);

    // Collect the distinct FETypes to be used.
    std::set<FEType> fe_type_set;
    const size_t num_systems = d_systems.size();
    for (size_t system_idx = 0; system_idx < num_systems; ++system_idx)
    {
        const System& system = *d_systems[system_idx];
        const DofMap& system_dof_map = system.get_dof_map();
        const std::vector<int>& all_vars = d_system_all_vars[system_idx];
        for (const auto& var : all_vars)
        {
            fe_type_set.insert(system_dof_map.variable_type(var));
        }
    }
    const size_t num_noninterp_systems = d_noninterp_systems.size();
    for (size_t system_idx = 0; system_idx < num_noninterp_systems; ++system_idx)
    {
        const System& system = *d_noninterp_systems[system_idx];
        const DofMap& system_dof_map = system.get_dof_map();
        const std::vector<int>& all_vars = d_noninterp_system_all_vars[system_idx];
        for (const auto& var : all_vars)
        {
            fe_type_set.insert(system_dof_map.variable_type(var));
        }
    }
    d_fe_types.assign(fe_type_set.begin(), fe_type_set.end());

    // Determine which FETypes are used for which variables, and which shape functions and shape function
    // derivatives need to be computed.
    const size_t num_fe_types = d_fe_types.size();
    d_eval_phi.resize(num_fe_types, false);
    d_eval_dphi.resize(num_fe_types, false);
    d_system_var_fe_type_idx.resize(num_systems);
    d_system_grad_var_fe_type_idx.resize(num_systems);
    for (size_t system_idx = 0; system_idx < num_systems; ++system_idx)
    {
        const System& system = *d_systems[system_idx];
        const DofMap& system_dof_map = system.get_dof_map();

        NumericVector<double>*& system_vec = d_system_vecs[system_idx];
        if (!system_vec) system_vec = system.current_local_solution.get();

        const std::vector<int>& vars = d_system_vars[system_idx];
        d_system_var_fe_type_idx[system_idx].resize(vars.size(), std::numeric_limits<std::size_t>::max());
        for (unsigned int k = 0; k < vars.size(); ++k)
        {
            const FEType& fe_type = system_dof_map.variable_type(vars[k]);
            const size_t fe_type_idx = getFETypeIndex(fe_type);
            d_eval_phi[fe_type_idx] = true;
            d_system_var_fe_type_idx[system_idx][k] = fe_type_idx;
        }

        const std::vector<int>& grad_vars = d_system_grad_vars[system_idx];
        d_system_grad_var_fe_type_idx[system_idx].resize(grad_vars.size(), std::numeric_limits<std::size_t>::max());
        for (unsigned int k = 0; k < grad_vars.size(); ++k)
        {
            const FEType& fe_type = system_dof_map.variable_type(grad_vars[k]);
            const size_t fe_type_idx = getFETypeIndex(fe_type);
            d_eval_dphi[fe_type_idx] = true;
            d_system_grad_var_fe_type_idx[system_idx][k] = fe_type_idx;
        }
    }
    for (size_t system_idx = 0; system_idx < num_noninterp_systems; ++system_idx)
    {
        const System& system = *d_noninterp_systems[system_idx];
        const DofMap& system_dof_map = system.get_dof_map();

        const std::vector<int>& phi_vars = d_noninterp_system_phi_vars[system_idx];
        for (const auto& phi_var : phi_vars)
        {
            const FEType& fe_type = system_dof_map.variable_type(phi_var);
            const size_t fe_type_idx = getFETypeIndex(fe_type);
            d_eval_phi[fe_type_idx] = true;
        }

        const std::vector<int>& dphi_vars = d_noninterp_system_dphi_vars[system_idx];
        for (const auto& dphi_var : dphi_vars)
        {
            const FEType& fe_type = system_dof_map.variable_type(dphi_var);
            const size_t fe_type_idx = getFETypeIndex(fe_type);
            d_eval_dphi[fe_type_idx] = true;
        }
    }

    // Set up FE objects and request access to shape functions / gradients, as needed.
    d_fe.resize(num_fe_types);
    d_fe_face.resize(num_fe_types);
    d_phi.resize(num_fe_types, nullptr);
    d_dphi.resize(num_fe_types, nullptr);
    d_phi_face.resize(num_fe_types, nullptr);
    d_dphi_face.resize(num_fe_types, nullptr);
    for (unsigned int fe_type_idx = 0; fe_type_idx < num_fe_types; ++fe_type_idx)
    {
        std::unique_ptr<IBTK::FEValuesBase>& fe = d_fe[fe_type_idx];
        if (!fe)
        {
            FEUpdateFlags update_flags = FEUpdateFlags::update_default;
            if (d_eval_q_point) update_flags |= update_quadrature_points;
            if (d_eval_JxW) update_flags |= update_JxW;
            if (d_eval_phi[fe_type_idx]) update_flags |= update_phi;
            if (d_eval_dphi[fe_type_idx]) update_flags |= update_dphi;

            fe = FEValuesBase::build(d_dim, NDIM, d_qrule, d_fe_types[fe_type_idx], update_flags);

            if (d_eval_q_point && !d_q_point) d_q_point = &fe->getQuadraturePoints();
            if (d_eval_JxW && !d_JxW) d_JxW = &fe->getJxW();
            if (d_eval_phi[fe_type_idx]) d_phi[fe_type_idx] = &fe->getShapeValues();
            if (d_eval_dphi[fe_type_idx]) d_dphi[fe_type_idx] = &fe->getShapeGradients();
        }

        std::unique_ptr<FEBase>& fe_face = d_fe_face[fe_type_idx];
        if (!fe_face)
        {
            const FEType& fe_type = d_fe_types[fe_type_idx];
            fe_face = FEBase::build(d_dim, fe_type);
            if (d_qrule_face) fe_face->attach_quadrature_rule(d_qrule_face);
            if (d_eval_q_point_face && !d_q_point_face) d_q_point_face = &fe_face->get_xyz();
            if (d_eval_JxW_face && !d_JxW_face) d_JxW_face = &fe_face->get_JxW();
            if (d_eval_normal_face && !d_normal_face) d_normal_face = &fe_face->get_normals();
            if (d_eval_phi[fe_type_idx]) d_phi_face[fe_type_idx] = &fe_face->get_phi();
            if (d_eval_dphi[fe_type_idx]) d_dphi_face[fe_type_idx] = &fe_face->get_dphi();
        }
    }

    // Indicate that we have initialized the class.
    d_initialized = true;
    return;
}

void
FEDataInterpolation::reinit(const Elem* elem,
                            const std::vector<libMesh::Point>* const points,
                            const std::vector<double>* weights)
{
    TBOX_ASSERT(d_initialized);
    d_current_elem = elem;
    if (d_qrule || points)
    {
        for (const auto& fe : d_fe)
        {
            TBOX_ASSERT(points == nullptr);
            TBOX_ASSERT(weights == nullptr);
            fe->reinit(elem);
        }
    }
    d_n_qp = static_cast<unsigned int>(points ? points->size() : d_qrule ? d_qrule->n_points() : 0);
    return;
}

void
FEDataInterpolation::reinit(const Elem* const elem,
                            const unsigned int side,
                            const double tol,
                            const std::vector<libMesh::Point>* const points,
                            const std::vector<double>* weights)
{
    TBOX_ASSERT(d_initialized);
    d_current_elem = elem;
    d_current_side = side;
    if (d_qrule_face || points)
    {
        for (const auto& fe_face : d_fe_face)
        {
            fe_face->reinit(elem, side, tol, points, weights);
        }
    }
    d_n_qp = static_cast<unsigned int>(points ? points->size() : d_qrule_face ? d_qrule_face->n_points() : 0);
    return;
}

void
FEDataInterpolation::collectDataForInterpolation(const Elem* const elem)
{
    TBOX_ASSERT(d_initialized);
    TBOX_ASSERT(elem == d_current_elem);

    // Collect local DOF data for the element.
    const size_t num_systems = d_systems.size();
    d_system_elem_data.resize(num_systems);
    for (size_t system_idx = 0; system_idx < num_systems; ++system_idx)
    {
        FEDataManager::SystemDofMapCache* system_dof_map_cache = d_system_dof_map_caches[system_idx];
        // Get the DOF mappings and local data for all variables.
        NumericVector<double>* system_vec = d_system_vecs[system_idx];
        const auto& dof_indices = system_dof_map_cache->dof_indices(d_current_elem);
        boost::multi_array<double, 2>& elem_data = d_system_elem_data[system_idx];
        get_values_for_interpolation(elem_data, *system_vec, dof_indices);
    }
    return;
}

const boost::multi_array<double, 2>&
FEDataInterpolation::getElemData(const Elem* const elem, const size_t system_idx)
{
    TBOX_ASSERT(d_initialized);
    TBOX_ASSERT(elem == d_current_elem);
    TBOX_ASSERT(system_idx < d_systems.size());
    return d_system_elem_data[system_idx];
}

void
FEDataInterpolation::interpolate(const Elem* const elem)
{
    TBOX_ASSERT(d_initialized);
    TBOX_ASSERT(elem == d_current_elem);
    interpolateCommon(d_system_var_data, d_system_grad_var_data, d_phi, d_dphi);
    return;
}

void
FEDataInterpolation::interpolate(const Elem* const elem, const unsigned int side)
{
    TBOX_ASSERT(d_initialized);
    TBOX_ASSERT(elem == d_current_elem);
    TBOX_ASSERT(side == d_current_side);
    interpolateCommon(d_system_var_data, d_system_grad_var_data, d_phi_face, d_dphi_face);
    return;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

size_t
FEDataInterpolation::getFETypeIndex(const FEType& fe_type) const
{
    return std::distance(d_fe_types.begin(), std::find(d_fe_types.begin(), d_fe_types.end(), fe_type));
}

void
FEDataInterpolation::interpolateCommon(
    std::vector<std::vector<std::vector<double> > >& system_var_data,
    std::vector<std::vector<std::vector<VectorValue<double> > > >& system_grad_var_data,
    const std::vector<const std::vector<std::vector<double> >*>& phi_data,
    const std::vector<const std::vector<std::vector<VectorValue<double> > >*>& dphi_data)
{
    // Determine the number of quadrature points where the interpolation will be evaluated.
    const unsigned int n_qp = d_n_qp;
    const size_t n_systems = d_systems.size();
    system_var_data.resize(n_qp, std::vector<std::vector<double> >(n_systems));
    system_grad_var_data.resize(n_qp, std::vector<std::vector<VectorValue<double> > >(n_systems));

    // Interpolate data for each system.
    for (size_t system_idx = 0; system_idx < n_systems; ++system_idx)
    {
        const boost::multi_array<double, 2>& elem_data = d_system_elem_data[system_idx];

        // Interpolate regular variables.
        {
            const std::vector<size_t>& var_idxs = d_system_var_idx[system_idx];
            const auto n_vars = static_cast<unsigned int>(var_idxs.size());
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                system_var_data[qp][system_idx].resize(n_vars);
            }
            const std::vector<size_t>& var_fe_type_idxs = d_system_var_fe_type_idx[system_idx];
            for (unsigned int k = 0; k < n_vars; ++k)
            {
                const size_t var_idx = var_idxs[k];
                const size_t fe_type_idx = var_fe_type_idxs[k];
                const std::vector<std::vector<double> >& phi = *phi_data[fe_type_idx];
                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    system_var_data[qp][system_idx][k] = 0.0;
                    for (unsigned int i = 0; i < phi.size(); ++i)
                    {
                        system_var_data[qp][system_idx][k] += elem_data[i][var_idx] * phi[i][qp];
                    }
                }
            }
        }

        // Interpolate gradient variables.
        {
            const std::vector<size_t>& grad_var_idxs = d_system_grad_var_idx[system_idx];
            const auto n_grad_vars = static_cast<unsigned int>(grad_var_idxs.size());
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                system_grad_var_data[qp][system_idx].resize(n_grad_vars);
            }
            const std::vector<size_t>& grad_var_fe_type_idxs = d_system_grad_var_fe_type_idx[system_idx];
            for (unsigned int k = 0; k < n_grad_vars; ++k)
            {
                const size_t var_idx = grad_var_idxs[k];
                const size_t fe_type_idx = grad_var_fe_type_idxs[k];
                const std::vector<std::vector<VectorValue<double> > >& dphi = *dphi_data[fe_type_idx];
                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    system_grad_var_data[qp][system_idx][k].zero();
                    for (unsigned int i = 0; i < dphi.size(); ++i)
                    {
                        system_grad_var_data[qp][system_idx][k] += elem_data[i][var_idx] * dphi[i][qp];
                    }
                }
            }
        }
    }
    return;
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
