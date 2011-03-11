// Filename: FESystemDataCache.C
// Created on 10 Mar 2011 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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

#include "FESystemDataCache.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/FEDataManager.h>
#include <ibtk/namespaces.h>

// LIBMESH INCLUDES
#include <boundary_info.h>
#include <dof_map.h>
#include <fe_base.h>
using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

FESystemDataCache::FESystemDataCache(
    const System& system,
    const MeshBase& mesh)
    : d_system(system),
      d_mesh(mesh),
      d_compute_q_point(false),
      d_compute_JxW(false),
      d_compute_phi(false),
      d_compute_phi_JxW(false),
      d_compute_dphi(false),
      d_compute_dphi_JxW(false),
      d_compute_q_point_face(false),
      d_compute_JxW_face(false),
      d_compute_normal_face(false),
      d_compute_normal_JxW_face(false),
      d_compute_phi_face(false),
      d_compute_phi_JxW_face(false),
      d_compute_dphi_face(false),
      d_compute_dphi_JxW_face(false),
      d_elems(),
      d_dof_indices(),
      d_q_point(),
      d_JxW(),
      d_phi(),
      d_phi_JxW(),
      d_dphi(),
      d_dphi_JxW(),
      d_elem_touches_physical_bdry(),
      d_q_point_face(),
      d_JxW_face(),
      d_normal_face(),
      d_normal_JxW_face(),
      d_phi_face(),
      d_phi_JxW_face(),
      d_dphi_face(),
      d_dphi_JxW_face()
{
    // intentionally blank
    return;
}// FEMeshDataCache

FESystemDataCache::~FESystemDataCache()
{
    // intentionally blank
    return;
}// ~FESystemDataCache

template<typename ElementIterator>
void
FESystemDataCache::computeCachedData(
    const ElementIterator& elem_begin,
    const ElementIterator& elem_end,
    QBase* qrule,
    QBase* qrule_face)
{
    clearCachedData();

    const int dim = d_mesh.mesh_dimension();
    const int n_vars = d_system.n_vars();

    const DofMap& dof_map = d_system.get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(n_vars);
#ifdef DEBUG_CHECK_ASSERTIONS
    for (int d = 0; d < n_vars; ++d)
    {
        TBOX_ASSERT(dof_map.variable_type(0) == dof_map.variable_type(d));
    }
#endif

    const bool compute_q_point = qrule != NULL &&  d_compute_q_point;
    const bool compute_JxW     = qrule != NULL && (d_compute_JxW || d_compute_phi_JxW || d_compute_dphi_JxW);
    const bool compute_phi     = qrule != NULL && (d_compute_phi || d_compute_phi_JxW);
    const bool compute_dphi    = qrule != NULL && (d_compute_dphi || d_compute_dphi_JxW);

    AutoPtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    if (qrule != NULL) fe->attach_quadrature_rule(qrule);
    const std::vector<Point>* const                           q_point = (compute_q_point ? &fe->get_xyz()  : NULL);
    const std::vector<double>* const                              JxW = (compute_JxW     ? &fe->get_JxW()  : NULL);
    const std::vector<std::vector<double> >* const                phi = (compute_phi     ? &fe->get_phi()  : NULL);
    const std::vector<std::vector<VectorValue<double> > >* const dphi = (compute_dphi    ? &fe->get_dphi() : NULL);

    const bool compute_q_point_face = qrule != NULL &&  d_compute_q_point_face;
    const bool compute_JxW_face     = qrule != NULL && (d_compute_JxW_face || d_compute_phi_JxW_face || d_compute_dphi_JxW_face || d_compute_normal_JxW_face);
    const bool compute_normal_face  = qrule != NULL && (d_compute_normal_face || d_compute_normal_JxW_face);
    const bool compute_phi_face     = qrule != NULL && (d_compute_phi_face || d_compute_phi_JxW_face);
    const bool compute_dphi_face    = qrule != NULL && (d_compute_dphi_face || d_compute_dphi_JxW_face);

    AutoPtr<FEBase> fe_face(FEBase::build(dim, dof_map.variable_type(0)));
    if (qrule_face != NULL) fe_face->attach_quadrature_rule(qrule_face);
    const std::vector<Point>* const                           q_point_face = (compute_q_point_face ? &fe_face->get_xyz()     : NULL);
    const std::vector<double>* const                              JxW_face = (compute_JxW_face     ? &fe_face->get_JxW()     : NULL);
    const std::vector<Point>* const                            normal_face = (compute_normal_face  ? &fe_face->get_normals() : NULL);
    const std::vector<std::vector<double> >* const                phi_face = (compute_phi_face     ? &fe_face->get_phi()     : NULL);
    const std::vector<std::vector<VectorValue<double> > >* const dphi_face = (compute_dphi_face    ? &fe_face->get_dphi()    : NULL);

    // Compute cached values.
    const unsigned int num_elems = std::distance(elem_begin,elem_end);

    d_elems.resize(num_elems);
    d_dof_indices.resize(num_elems);

    if (qrule != NULL)
    {
        if (d_compute_q_point)  d_q_point .resize(num_elems);
        if (d_compute_JxW)      d_JxW     .resize(num_elems);
        if (d_compute_phi)      d_phi     .resize(num_elems);
        if (d_compute_phi_JxW)  d_phi_JxW .resize(num_elems);
        if (d_compute_dphi)     d_dphi    .resize(num_elems);
        if (d_compute_dphi_JxW) d_dphi_JxW.resize(num_elems);
    }

    if (qrule_face != NULL)
    {
        d_elem_touches_physical_bdry .resize(num_elems);
        d_elem_side_at_physical_bdry .resize(num_elems);
        d_elem_side_at_dirichlet_bdry.resize(num_elems);
        if (d_compute_q_point_face)    d_q_point_face   .resize(num_elems);
        if (d_compute_JxW_face)        d_JxW_face       .resize(num_elems);
        if (d_compute_normal_face)     d_normal_face    .resize(num_elems);
        if (d_compute_normal_JxW_face) d_normal_JxW_face.resize(num_elems);
        if (d_compute_phi_face)        d_phi_face       .resize(num_elems);
        if (d_compute_phi_JxW_face)    d_phi_JxW_face   .resize(num_elems);
        if (d_compute_dphi_face)       d_dphi_face      .resize(num_elems);
        if (d_compute_dphi_JxW_face)   d_dphi_JxW_face  .resize(num_elems);
    }

    for (ElementIterator elem_it = elem_begin; elem_it != elem_end; ++elem_it)
    {
        // Compute cached element data.
        Elem* const elem = *elem_it;
        const unsigned int e = std::distance(elem_begin,elem_it);
        fe->reinit(elem);
        d_elems(e) = elem;
        d_dof_indices(e).resize(n_vars);
        for (int i = 0; i < n_vars; ++i)
        {
            dof_map.dof_indices(elem, dof_indices[i], i);
            d_dof_indices(e)(i) = dof_indices[i];
        }
        const int num_basis = d_dof_indices(e)(0).size();
        if (qrule != NULL)
        {
            const int num_qp = qrule->n_points();
            if (d_compute_q_point)  d_q_point (e).resize(num_qp);
            if (d_compute_JxW)      d_JxW     (e).resize(num_qp);
            if (d_compute_phi)      d_phi     (e).resize(num_qp,num_basis);
            if (d_compute_phi_JxW)  d_phi_JxW (e).resize(num_qp,num_basis);
            if (d_compute_dphi)     d_dphi    (e).resize(num_qp,num_basis);
            if (d_compute_dphi_JxW) d_dphi_JxW(e).resize(num_qp,num_basis);
            for (int qp = 0; qp < num_qp; ++qp)
            {
                if (d_compute_q_point) d_q_point (e)(qp) = (*q_point)[qp];
                if (d_compute_JxW)     d_JxW     (e)(qp) = (*JxW    )[qp];
                for (int k = 0; d_compute_phi && k < num_basis; ++k)
                {
                    d_phi(e)(qp,k) = (*phi)[k][qp];
                }
                for (int k = 0; d_compute_phi_JxW && k < num_basis; ++k)
                {
                    d_phi_JxW(e)(qp,k) = (*JxW)[qp]*(*phi)[k][qp];
                }
                for (int k = 0; d_compute_dphi && k < num_basis; ++k)
                {
                    d_dphi(e)(qp,k) = (*dphi)[k][qp];
                }
                for (int k = 0; d_compute_dphi_JxW && k < num_basis; ++k)
                {
                    d_dphi_JxW(e)(qp,k) = (*JxW)[qp]*(*dphi)[k][qp];
                }
            }
        }

        if (qrule_face != NULL)
        {
            const int num_sides = elem->n_sides();

            // Loop over the element boundaries to determine if there are any
            // sides that are along the physical boundary.
            d_elem_touches_physical_bdry(e) = false;
            blitz::Array<bool,1> side_at_physical_bdry (num_sides);
            blitz::Array<bool,1> side_at_dirichlet_bdry(num_sides);
            for (int side = 0; side < num_sides; ++side)
            {
                // Determine whether this side touches the physical boundary.
                const std::vector<short int>& bdry_ids = d_mesh.boundary_info->boundary_ids(elem, side);
                bool at_physical_bdry  = false;
                bool at_dirichlet_bdry = false;
                for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end(); ++cit)
                {
                    const short int bdry_id = *cit;
                    at_physical_bdry  = at_physical_bdry  || (elem->neighbor(side) == NULL && !dof_map.is_periodic_boundary(bdry_id));
                    at_dirichlet_bdry = at_dirichlet_bdry || (bdry_id == FEDataManager::DIRICHLET_BDRY_ID);
                }
                side_at_physical_bdry (side) = at_physical_bdry ;
                side_at_dirichlet_bdry(side) = at_dirichlet_bdry;
                if (!at_physical_bdry) continue;
                d_elem_touches_physical_bdry(e) = true;
            }

            if (!d_elem_touches_physical_bdry(e)) continue;

            d_elem_side_at_physical_bdry(e).resize(num_sides);
            d_elem_side_at_physical_bdry(e) = side_at_physical_bdry;

            d_elem_side_at_dirichlet_bdry(e).resize(num_sides);
            d_elem_side_at_dirichlet_bdry(e) = side_at_dirichlet_bdry;

            if (d_compute_q_point_face   )    d_q_point_face(e).resize(num_sides);
            if (d_compute_JxW_face       )        d_JxW_face(e).resize(num_sides);
            if (d_compute_normal_face    )     d_normal_face(e).resize(num_sides);
            if (d_compute_normal_JxW_face) d_normal_JxW_face(e).resize(num_sides);
            if (d_compute_phi_face       )        d_phi_face(e).resize(num_sides);
            if (d_compute_phi_JxW_face   )    d_phi_JxW_face(e).resize(num_sides);
            if (d_compute_dphi_face      )       d_dphi_face(e).resize(num_sides);
            if (d_compute_dphi_JxW_face  )   d_dphi_JxW_face(e).resize(num_sides);

            for (int side = 0; side < num_sides; ++side)
            {
                if (!d_elem_side_at_physical_bdry(e)(side)) continue;

                // Compute cached side data along the physical boundary.
                fe_face->reinit(elem, side);
                const int num_qp_face = qrule_face->n_points();
                if (d_compute_q_point_face   )    d_q_point_face(e)(side).resize(num_qp_face);
                if (d_compute_JxW_face       )        d_JxW_face(e)(side).resize(num_qp_face);
                if (d_compute_normal_face    )     d_normal_face(e)(side).resize(num_qp_face);
                if (d_compute_normal_JxW_face) d_normal_JxW_face(e)(side).resize(num_qp_face);
                if (d_compute_phi_face       )        d_phi_face(e)(side).resize(num_qp_face,num_basis);
                if (d_compute_phi_JxW_face   )    d_phi_JxW_face(e)(side).resize(num_qp_face,num_basis);
                if (d_compute_dphi_face      )       d_dphi_face(e)(side).resize(num_qp_face,num_basis);
                if (d_compute_dphi_JxW_face  )   d_dphi_JxW_face(e)(side).resize(num_qp_face,num_basis);
                for (int qp = 0; qp < num_qp_face; ++qp)
                {
                    if (d_compute_q_point_face   )    d_q_point_face(e)(side)(qp) = (*q_point_face)[qp];
                    if (d_compute_JxW_face       )        d_JxW_face(e)(side)(qp) = (*JxW_face)[qp];
                    if (d_compute_normal_face    )     d_normal_face(e)(side)(qp) = (*normal_face)[qp];
                    if (d_compute_normal_JxW_face) d_normal_JxW_face(e)(side)(qp) = (*JxW_face)[qp]*(*normal_face)[qp];
                    for (int k = 0; d_compute_phi_face && k < num_basis; ++k)
                    {
                        d_phi_face(e)(side)(qp,k) = (*phi_face)[k][qp];
                    }
                    for (int k = 0; d_compute_phi_JxW_face && k < num_basis; ++k)
                    {
                        d_phi_JxW_face(e)(side)(qp,k) = (*JxW_face)[qp]*(*phi_face)[k][qp];
                    }
                    for (int k = 0; d_compute_dphi_face && k < num_basis; ++k)
                    {
                        d_dphi_face(e)(side)(qp,k) = (*dphi_face)[k][qp];
                    }
                    for (int k = 0; d_compute_dphi_JxW_face && k < num_basis; ++k)
                    {
                        d_dphi_JxW_face(e)(side)(qp,k) = (*JxW_face)[qp]*(*dphi_face)[k][qp];
                    }
                }
            }
        }
    }
    return;
}// computeCachedData

void
FESystemDataCache::clearCachedData()
{
    d_elems.free();
    d_dof_indices.free();
    d_q_point.free();
    d_JxW.free();
    d_phi.free();
    d_phi_JxW.free();
    d_dphi.free();
    d_dphi_JxW.free();
    d_elem_touches_physical_bdry.free();
    d_elem_side_at_physical_bdry.free();
    d_elem_side_at_dirichlet_bdry.free();
    d_q_point_face.free();
    d_JxW_face.free();
    d_normal_face.free();
    d_normal_JxW_face.free();
    d_phi_face.free();
    d_phi_JxW_face.free();
    d_dphi_face.free();
    d_dphi_JxW_face.free();
    return;
}// clearCachedData

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

template
void
IBTK::FESystemDataCache::computeCachedData<libMesh::MeshBase::const_element_iterator>(
    const libMesh::MeshBase::const_element_iterator& elem_begin,
    const libMesh::MeshBase::const_element_iterator& elem_end,
    libMesh::QBase* qrule=NULL,
    libMesh::QBase* qrule_face=NULL);

template
void
IBTK::FESystemDataCache::computeCachedData<libMesh::MeshBase::element_iterator>(
    const libMesh::MeshBase::element_iterator& elem_begin,
    const libMesh::MeshBase::element_iterator& elem_end,
    libMesh::QBase* qrule=NULL,
    libMesh::QBase* qrule_face=NULL);

//////////////////////////////////////////////////////////////////////////////
