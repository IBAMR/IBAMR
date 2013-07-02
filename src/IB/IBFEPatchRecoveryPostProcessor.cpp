// Filename: IBFEPatchRecoveryPostProcessor.cpp
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "IBAMR_config.h"
#include "ibamr/IBFEPatchRecoveryPostProcessor.h"
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
#include "libmesh/string_to_enum.h"

using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEPatchRecoveryPostProcessor::IBFEPatchRecoveryPostProcessor(
    MeshBase* mesh,
    FEDataManager* fe_data_manager)
    : d_mesh(mesh),
      d_fe_data_manager(fe_data_manager)
{
    const int rank = libMesh::processor_id();
    const int size = libMesh::n_processors();

    // Active local elements.
    const MeshBase::const_element_iterator el_begin = d_mesh->active_local_elements_begin();
    const MeshBase::const_element_iterator el_end   = d_mesh->active_local_elements_end();

    // Determine the element patches associated with each node N, which is
    // defined to be the collection of elements that contain node N.
    //
    // Unlike the standard Z-Z patch recovery algorithm, we use "tight" element
    // patches for non-vertex nodes.
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        for (unsigned int n = 0; n < elem->n_nodes(); ++n)
        {
            // Only set up patches for local nodes.
            const Node* const node = elem->get_node(n);
            if (node->processor_id() != rank) continue;

            // Only set up patches once for each node.
            const dof_id_type node_id = node->id();
            if (d_local_elem_patches.find(node_id) != d_local_elem_patches.end()) continue;

            // Find the elements that touch this node.
            std::set<const Elem*>& elem_patch = d_local_elem_patches[node_id];
            elem->find_point_neighbors(*node, elem_patch);
            d_local_elems.insert(d_local_elems.end(), elem_patch.begin(), elem_patch.end());
        }
    }

    // Determine the number of quadrature/interpolation points in each element
    // and setup mappings used to fill global data structures.
    //
    // We use full Gaussian quadrature rules (i.e. third-order Gauss quadrature
    // for first-order elements and fifth-order Gauss quadrature for
    // second-order elements) in all elements to avoid special treatment at
    // boundary nodes.
    d_n_qp_global = 0;
    d_n_qp_local  = 0;
    d_qp_global_offset = 0;
    const unsigned int n_elem = d_mesh->n_elem();
    d_elem_n_qp.resize(n_elem,0);
    d_elem_qp_global_offset.resize(n_elem,0);
    d_elem_qp_local_offset .resize(n_elem,0);
    qrule.reset();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const unsigned int dim = elem->dim();
        TBOX_ASSERT(elem->default_order() == FIRST || elem->default_order() == SECOND);
        const Order order = (elem->default_order() == FIRST ? THIRD : FIFTH);
        bool reinit_qrule = false;
        if (!qrule.get() || qrule->get_dim() != dim || qrule->get_order() != order)
        {
            qrule = QBase::build(QGAUSS, dim, order);
            reinit_qrule = true;
        }
        else if (qrule->get_elem_type() != elem->type() || qrule->get_p_level() != elem->p_level())
        {
            reinit_qrule = true;
        }
        if (reinit_qrule) qrule->init(elem->type(), elem->p_level());
        unsigned int n_qp = qrule->n_points();
        const dof_id_type elem_id = elem->id();
        d_elem_n_qp[elem_id] = n_qp;
        d_elem_qp_local_offset[elem_id] = d_n_qp_local;
        d_n_qp_local += n_qp;
        d_elem_sigma[elem_id].resize(n_qp);
    }
    std::vector<int> n_qp_per_proc(size);
    n_qp_per_proc[rank] = d_n_qp_local;
    Parallel::sum(n_qp_per_proc);
    d_qp_global_offset = std::accumulate(n_qp_per_proc.begin(), n_qp_per_proc.begin()+rank, 0);
    d_n_qp_global = std::accumulate(n_qp_per_proc.begin()+rank, n_qp_per_proc.end(), d_qp_global_offset);
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const dof_id_type elem_id = elem->id();
        d_elem_qp_global_offset[elem_id] = d_elem_qp_local_offset[elem_id]+d_qp_global_offset;
    }
    Parallel::sum(d_elem_qp_global_offset);
    Parallel::sum(d_elem_qp_local_offset);
    return;
}// IBFEPatchRecoveryPostProcessor

IBFEPatchRecoveryPostProcessor::~IBFEPatchRecoveryPostProcessor()
{
    // intentionally blank
    return;
}// IBFEPatchRecoveryPostProcessor

void
IBFEPatchRecoveryPostProcessor::registerPK1StressValue(
    const Elem* const elem,
    const QBase* const qrule,
    const unsigned int qp,
    const TensorValue<double>& PP,
    const TensorValue<double>& FF)
{
    if (elem->processor_id() != libMesh::processor_id() || !elem->active()) return;
    TBOX_ASSERT(qrule->type() == QGAUSS);
    TBOX_ASSERT(elem->default_order() == FIRST || elem->default_order() == SECOND);
    const Order order = (elem->default_order() == FIRST ? THIRD : FIFTH);
    TBOX_ASSERT(order == qrule->get_order());
    TBOX_ASSERT(0 >= qp && qp < qrule->n_points());
    d_elem_sigma[elem->id()][qp] = PP*FF.transpose()/FF.det();
    return;
}// registerPK1StressValue

void
IBFEPatchRecoveryPostProcessor::reconstructCauchyStress()
{
    // Communicate the values of the Cauchy stress.
    static const int NVALS = 2*(NDIM+1);
    std::vector<double> sigma_vals(NVALS*d_n_qp_global,0.0);
    const MeshBase::const_element_iterator el_begin = d_mesh->active_local_elements_begin();
    const MeshBase::const_element_iterator el_end   = d_mesh->active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const dof_id_type elem_id = elem->id();
        const int global_offset = d_elem_qp_global_offset[elem_id];
        for (int qp = 0; qp < d_elem_n_qp[elem_id]; ++qp)
        {
            const TensorValue<double>& sigma = d_elem_sigma[elem_id][qp];
            for (unsigned int i = 0, k = 0; i < NDIM; ++i)
            {
                for (unsigned int j = i; j < NDIM; ++j, ++k)
                {
                    sigma_vals[NVALS*(global_offset+qp)+k] = sigma(i,j);
                }
            }
        }
    }
    Parallel::sum(sigma_vals);


    return;
}// reconstructCauchyStress

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
