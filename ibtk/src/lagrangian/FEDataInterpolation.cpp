// Filename: FEDataInterpolation.cpp
// Created on 9 Oct 2015 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/FEDataInterpolation.h"
#include "ibtk/namespaces.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

FEDataInterpolation::FEDataInterpolation(const unsigned int dim, FEDataManager* const fe_data_manager)
    : d_dim(dim),
      d_fe_data_manager(fe_data_manager),
      d_initialized(false),
      d_eval_q_point(false),
      d_eval_JxW(false),
      d_eval_q_point_face(false),
      d_eval_JxW_face(false),
      d_eval_normal_face(false),
      d_qrule(NULL),
      d_qrule_face(NULL),
      d_q_point(NULL),
      d_q_point_face(NULL),
      d_JxW(NULL),
      d_JxW_face(NULL),
      d_normal_face(NULL),
      d_current_elem(NULL)
{
    return;
}

FEDataInterpolation::~FEDataInterpolation()
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
    for (std::vector<const System *>::iterator it = d_noninterp_systems.begin(), it_end = d_noninterp_systems.end();
         it != it_end;
         ++it)
    {
        if ((*it)->number() == sys_num)
        {
            // This system has already been registered.  If so, just merge the collection of variables.
            size_t system_idx = std::distance(d_noninterp_systems.begin(), it);

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
                                                NumericVector<double>* system_data)
{
    TBOX_ASSERT(!d_initialized && (!vars.empty() || !grad_vars.empty()));
    const unsigned int sys_num = system.number();
    for (std::vector<const System *>::iterator it = d_systems.begin(), it_end = d_systems.end(); it != it_end; ++it)
    {
        if ((*it)->number() == sys_num)
        {
            // This system has already been registered.  Check to see if the same collection of variables (etc.)
            // were used in the previous registration action; if so, do not re-register the system.
            //
            // NOTE: The same system may be registered multiple times with different collections of variables,
            // gradient variables, system data, etc.  We want to make sure we check all of the registered systems to
            // see if any of them match the function arguments.
            const size_t system_idx = std::distance(d_systems.begin(), it);
            bool same_data = true;
            same_data = same_data && (vars == d_system_vars[system_idx]);
            same_data = same_data && (grad_vars == d_system_grad_vars[system_idx]);
            same_data = same_data && (system_data == d_system_data[system_idx]);
            if (same_data) return system_idx;
        }
    }

    // Either we have not previously registered this system, or we are registering it here with a different
    // collection of variables/data.  In either case, we need to register it here.
    const size_t system_idx = d_systems.size();
    d_systems.push_back(&system);
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
    d_system_data.push_back(system_data);
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
FEDataInterpolation::init(const bool use_IB_ghosted_vecs)
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
        for (unsigned int k = 0; k < all_vars.size(); ++k)
        {
            fe_type_set.insert(system_dof_map.variable_type(all_vars[k]));
        }
    }
    const size_t num_noninterp_systems = d_noninterp_systems.size();
    for (size_t system_idx = 0; system_idx < num_noninterp_systems; ++system_idx)
    {
        const System& system = *d_noninterp_systems[system_idx];
        const DofMap& system_dof_map = system.get_dof_map();
        const std::vector<int>& all_vars = d_noninterp_system_all_vars[system_idx];
        for (unsigned int k = 0; k < all_vars.size(); ++k)
        {
            fe_type_set.insert(system_dof_map.variable_type(all_vars[k]));
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

        NumericVector<double>*& system_data = d_system_data[system_idx];
        if (!system_data) system_data = system.current_local_solution.get();
        if (use_IB_ghosted_vecs)
        {
            NumericVector<double>* ghost_data =
                d_fe_data_manager->buildGhostedSolutionVector(system.name(), /*synch_data*/ false);
            system_data->localize(*ghost_data);
            ghost_data->close();
            system_data = ghost_data;
        }

        const std::vector<int>& vars = d_system_vars[system_idx];
        d_system_var_fe_type_idx[system_idx].resize(vars.size(), -1);
        for (unsigned int k = 0; k < vars.size(); ++k)
        {
            const FEType& fe_type = system_dof_map.variable_type(vars[k]);
            const size_t fe_type_idx = getFETypeIndex(fe_type);
            d_eval_phi[fe_type_idx] = true;
            d_system_var_fe_type_idx[system_idx][k] = fe_type_idx;
        }

        const std::vector<int>& grad_vars = d_system_grad_vars[system_idx];
        d_system_grad_var_fe_type_idx[system_idx].resize(grad_vars.size(), -1);
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
        for (unsigned int k = 0; k < phi_vars.size(); ++k)
        {
            const FEType& fe_type = system_dof_map.variable_type(phi_vars[k]);
            const size_t fe_type_idx = getFETypeIndex(fe_type);
            d_eval_phi[fe_type_idx] = true;
        }

        const std::vector<int>& dphi_vars = d_noninterp_system_dphi_vars[system_idx];
        for (unsigned int k = 0; k < dphi_vars.size(); ++k)
        {
            const FEType& fe_type = system_dof_map.variable_type(dphi_vars[k]);
            const size_t fe_type_idx = getFETypeIndex(fe_type);
            d_eval_dphi[fe_type_idx] = true;
        }
    }

    // Set up FE objects and request access to shape functions / gradients, as needed.
    d_fe.resize(num_fe_types);
    d_fe_face.resize(num_fe_types);
    d_phi.resize(num_fe_types, NULL);
    d_dphi.resize(num_fe_types, NULL);
    d_phi_face.resize(num_fe_types, NULL);
    d_dphi_face.resize(num_fe_types, NULL);
    for (unsigned int fe_type_idx = 0; fe_type_idx < num_fe_types; ++fe_type_idx)
    {
        Pointer<FEBase>& fe = d_fe[fe_type_idx];
        if (!fe)
        {
            const FEType& fe_type = d_fe_types[fe_type_idx];
            fe = FEBase::build(d_dim, fe_type).release();
            if (d_qrule) fe->attach_quadrature_rule(d_qrule);
            if (d_eval_q_point && !d_q_point) d_q_point = &fe->get_xyz();
            if (d_eval_JxW && !d_JxW) d_JxW = &fe->get_JxW();
            if (d_eval_phi[fe_type_idx]) d_phi[fe_type_idx] = &fe->get_phi();
            if (d_eval_dphi[fe_type_idx]) d_dphi[fe_type_idx] = &fe->get_dphi();
        }

        Pointer<FEBase>& fe_face = d_fe_face[fe_type_idx];
        if (!fe_face)
        {
            const FEType& fe_type = d_fe_types[fe_type_idx];
            fe_face = FEBase::build(d_dim, fe_type).release();
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
        const size_t num_fe_types = d_fe_types.size();
        for (unsigned int fe_type_idx = 0; fe_type_idx < num_fe_types; ++fe_type_idx)
        {
            Pointer<FEBase>& fe = d_fe[fe_type_idx];
            fe->reinit(elem, points, weights);
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
        const size_t num_fe_types = d_fe_types.size();
        for (unsigned int fe_type_idx = 0; fe_type_idx < num_fe_types; ++fe_type_idx)
        {
            Pointer<FEBase>& fe_face = d_fe_face[fe_type_idx];
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
        const System& system = *d_systems[system_idx];
        const DofMap& system_dof_map = system.get_dof_map();
        const std::vector<int>& all_vars = d_system_all_vars[system_idx];
        const size_t num_vars = all_vars.size();

        // Get the DOF mappings and local data for all variables.
        std::vector<std::vector<unsigned int> > dof_indices(num_vars);
        NumericVector<double>* system_data = d_system_data[system_idx];
        for (size_t k = 0; k < num_vars; ++k)
        {
            system_dof_map.dof_indices(d_current_elem, dof_indices[k], all_vars[k]);
        }
        boost::multi_array<double, 2>& elem_data = d_system_elem_data[system_idx];
        get_values_for_interpolation(elem_data, *system_data, dof_indices);
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
    TBOX_ASSERT(d_qrule);
    TBOX_ASSERT(elem == d_current_elem);
    interpolateCommon(d_system_var_data, d_system_grad_var_data, d_phi, d_dphi);
    return;
}

void
FEDataInterpolation::interpolate(const Elem* const elem, const unsigned int side)
{
    TBOX_ASSERT(d_initialized);
    TBOX_ASSERT(d_qrule_face);
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
            const unsigned int n_vars = static_cast<unsigned int>(var_idxs.size());
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
            const unsigned int n_grad_vars = static_cast<unsigned int>(grad_var_idxs.size());
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
