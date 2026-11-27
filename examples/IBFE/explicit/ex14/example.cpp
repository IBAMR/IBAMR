// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// 3D closed curved TRI6 spherical membrane in simple shear with irreversible
// Marigo continuum surface damage and Peerlings implicit-gradient nonlocal
// energy filtering (`IBFESurfaceDamage.h`).

#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFESurfaceMethod.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <tbox/MemoryDatabase.h>
#include <tbox/RestartManager.h>

#include <libmesh/dense_vector.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/face_tri6.h>
#include <libmesh/fe_base.h>
#include <libmesh/libmesh.h>
#include <libmesh/node.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/point.h>
#include <libmesh/quadrature.h>
#include <libmesh/replicated_mesh.h>

#include <petscsys.h>

#include <boost/multi_array.hpp>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "IBFESurfaceDamage.h"

#include <ibamr/app_namespaces.h>

namespace ModelData
{
using namespace IBAMR::FESurfaceDamage;

std::unique_ptr<libMesh::ReplicatedMesh>
buildSphereTri6Mesh(const libMesh::LibMeshInit& init,
                    const double radius,
                    const libMesh::Point& center,
                    const int n_subdiv)
{
    auto mesh = std::make_unique<libMesh::ReplicatedMesh>(init.comm(), 2);
    mesh->set_spatial_dimension(3);
    mesh->allow_renumbering(false);

    const double phi = 0.5 * (1.0 + std::sqrt(5.0));
    std::vector<std::array<double, 3>> V = {
        { -1, phi, 0 },  { 1, phi, 0 },  { -1, -phi, 0 }, { 1, -phi, 0 }, { 0, -1, phi },  { 0, 1, phi },
        { 0, -1, -phi }, { 0, 1, -phi }, { phi, 0, -1 },  { phi, 0, 1 },  { -phi, 0, -1 }, { -phi, 0, 1 },
    };
    std::vector<std::array<unsigned int, 3>> T = {
        { 0, 11, 5 },  { 0, 5, 1 },  { 0, 1, 7 },  { 0, 7, 10 }, { 0, 10, 11 }, { 1, 5, 9 }, { 5, 11, 4 },
        { 11, 10, 2 }, { 10, 7, 6 }, { 7, 1, 8 },  { 3, 9, 4 },  { 3, 4, 2 },   { 3, 2, 6 }, { 3, 6, 8 },
        { 3, 8, 9 },   { 4, 9, 5 },  { 2, 4, 11 }, { 6, 2, 10 }, { 8, 6, 7 },   { 9, 8, 1 },
    };

    for (int s = 0; s < n_subdiv; ++s)
    {
        std::map<std::pair<unsigned int, unsigned int>, unsigned int> mc;
        auto gm = [&](unsigned int a, unsigned int b) -> unsigned int
        {
            auto key = std::make_pair(std::min(a, b), std::max(a, b));
            auto it = mc.find(key);
            if (it != mc.end()) return it->second;
            const unsigned int id = static_cast<unsigned int>(V.size());
            V.push_back({ 0.5 * (V[a][0] + V[b][0]), 0.5 * (V[a][1] + V[b][1]), 0.5 * (V[a][2] + V[b][2]) });
            mc[key] = id;
            return id;
        };
        std::vector<std::array<unsigned int, 3>> T2;
        T2.reserve(4 * T.size());
        for (const auto& t : T)
        {
            const unsigned int ab = gm(t[0], t[1]);
            const unsigned int bc = gm(t[1], t[2]);
            const unsigned int ca = gm(t[2], t[0]);
            T2.push_back({ t[0], ab, ca });
            T2.push_back({ t[1], bc, ab });
            T2.push_back({ t[2], ca, bc });
            T2.push_back({ ab, bc, ca });
        }
        T.swap(T2);
    }

    auto mapPoint = [&](const std::array<double, 3>& p) -> libMesh::Point
    {
        const double rn = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        return libMesh::Point(
            center(0) + radius * p[0] / rn, center(1) + radius * p[1] / rn, center(2) + radius * p[2] / rn);
    };

    for (std::size_t i = 0; i < V.size(); ++i)
    {
        mesh->add_point(mapPoint(V[i]), static_cast<libMesh::dof_id_type>(i));
    }

    std::map<std::pair<unsigned int, unsigned int>, libMesh::dof_id_type> edge_mid;
    auto getMid = [&](unsigned int a, unsigned int b) -> libMesh::dof_id_type
    {
        auto key = std::make_pair(std::min(a, b), std::max(a, b));
        auto it = edge_mid.find(key);
        if (it != edge_mid.end()) return it->second;
        const libMesh::dof_id_type mid_id = static_cast<libMesh::dof_id_type>(mesh->n_nodes());
        const std::array<double, 3> avg = { 0.5 * (V[a][0] + V[b][0]),
                                            0.5 * (V[a][1] + V[b][1]),
                                            0.5 * (V[a][2] + V[b][2]) };
        mesh->add_point(mapPoint(avg), mid_id);
        edge_mid[key] = mid_id;
        return mid_id;
    };

    for (std::size_t ti = 0; ti < T.size(); ++ti)
    {
        unsigned int a = T[ti][0], b = T[ti][1], c = T[ti][2];
        const libMesh::Point pa = mesh->point(a) - center;
        const libMesh::Point pb = mesh->point(b) - center;
        const libMesh::Point pc = mesh->point(c) - center;
        if ((pb - pa).cross(pc - pa) * ((pa + pb + pc) * (1.0 / 3.0)) < 0.0)
        {
            std::swap(b, c);
        }
        const libMesh::dof_id_type mab = getMid(a, b);
        const libMesh::dof_id_type mbc = getMid(b, c);
        const libMesh::dof_id_type mca = getMid(c, a);
        auto* tri = new libMesh::Tri6();
        tri->set_id() = static_cast<libMesh::dof_id_type>(ti);
        const std::array<libMesh::dof_id_type, 6> conn = { a, b, c, mab, mbc, mca };
        for (unsigned int k = 0; k < 6; ++k) tri->set_node(k) = mesh->node_ptr(conn[k]);
        mesh->add_elem(tri);
    }
    mesh->prepare_for_use();
    return mesh;
}

class DamagedMembraneIBFESurfaceMethod : public IBAMR::IBFESurfaceMethod
{
public:
    DamagedMembraneIBFESurfaceMethod(const std::string& object_name,
                                     SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                     libMesh::MeshBase* mesh,
                                     int max_levels,
                                     bool register_for_restart = true,
                                     const std::string& restart_read_dirname = "",
                                     unsigned int restart_restore_number = 0)
        : IBAMR::IBFESurfaceMethod(object_name,
                                   input_db,
                                   mesh,
                                   max_levels,
                                   register_for_restart,
                                   restart_read_dirname,
                                   restart_restore_number),
          d_damage_trackers(1),
          d_nonlocal_filters(1)
    {
        if (register_for_restart && SAMRAI::tbox::RestartManager::getManager()->isFromRestart())
        {
            getSurfaceDamageFromRestart(object_name);
        }
    }

    void setSurfaceDamageParams(const SurfaceDamageParams& params, unsigned int part = 0)
    {
        TBOX_ASSERT(part < d_damage_trackers.size());
        d_damage_trackers[part].configure(params);
    }

    const SurfaceDamageTracker& getSurfaceDamageTracker(unsigned int part = 0) const
    {
        TBOX_ASSERT(part < d_damage_trackers.size());
        return d_damage_trackers[part];
    }

    SurfaceDamageTracker& getSurfaceDamageTracker(unsigned int part = 0)
    {
        TBOX_ASSERT(part < d_damage_trackers.size());
        return d_damage_trackers[part];
    }

    void pruneDamageSlotsToLocalMesh(unsigned int part = 0)
    {
        TBOX_ASSERT(part < d_damage_trackers.size());
        if (part < d_nonlocal_filters.size()) d_nonlocal_filters[part].destroy();
        if (!d_damage_trackers[part].params().enabled) return;
        if (part >= d_fe_data_managers.size() || !d_fe_data_managers[part]) return;
        const libMesh::EquationSystems* const es = d_fe_data_managers[part]->getEquationSystems();
        if (!es) return;
        const libMesh::MeshBase& mesh = es->get_mesh();
        std::unordered_set<std::uint64_t> local_ids;
        local_ids.reserve(static_cast<std::size_t>(mesh.n_active_local_elem()) * 2 + 1);
        for (auto el_it = mesh.active_local_elements_begin(); el_it != mesh.active_local_elements_end(); ++el_it)
        {
            local_ids.insert(static_cast<std::uint64_t>((*el_it)->id()));
        }
        d_damage_trackers[part].pruneToLocalElements(local_ids);
    }

    void initializePatchHierarchy(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM>> gridding_alg,
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM>>>& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>>>& u_ghost_fill_scheds,
        int integrator_step,
        double init_data_time,
        bool initial_time) override
    {
        IBAMR::IBFESurfaceMethod::initializePatchHierarchy(hierarchy,
                                                           gridding_alg,
                                                           u_data_idx,
                                                           u_synch_scheds,
                                                           u_ghost_fill_scheds,
                                                           integrator_step,
                                                           init_data_time,
                                                           initial_time);
        for (unsigned int part = 0; part < d_damage_trackers.size(); ++part)
        {
            pruneDamageSlotsToLocalMesh(part);
        }
    }

    void endDataRedistribution(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                               SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM>> gridding_alg) override
    {
        IBAMR::IBFESurfaceMethod::endDataRedistribution(hierarchy, gridding_alg);
        for (unsigned int part = 0; part < d_damage_trackers.size(); ++part)
        {
            pruneDamageSlotsToLocalMesh(part);
        }
    }

    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override
    {
        IBAMR::IBFESurfaceMethod::putToDatabase(db);
        if (!db) return;
        for (unsigned int part = 0; part < d_damage_trackers.size(); ++part)
        {
            if (!d_damage_trackers[part].params().enabled) continue;
            auto part_db = db->putDatabase("SurfaceDamagePart_" + std::to_string(part));
            d_damage_trackers[part].putToDatabase(part_db);
        }
    }

    void computeLagrangianForce(double data_time) override
    {
        using namespace libMesh;
        TBOX_ASSERT(IBTK::rel_equal_eps(data_time, d_half_time));
        IBTK::batch_vec_ghost_update(d_X_vecs->get("half"), INSERT_VALUES, SCATTER_FORWARD);

        for (unsigned int part = 0; part < d_num_parts; ++part)
        {
            EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
            MeshBase& mesh = equation_systems->get_mesh();
            const unsigned int dim = mesh.mesh_dimension();
            TBOX_ASSERT(dim == 2);

            PetscVector<double>& F_vec = d_F_vecs->get("half", part);
            PetscVector<double>& F_rhs_vec = d_F_vecs->get("RHS Vector", part);
            F_rhs_vec.zero();
            DenseVector<double> F_rhs_e[NDIM];
            PetscVector<double>& X_vec = d_X_vecs->get("half", part);

            System& F_system = equation_systems->get_system(FORCE_SYSTEM_NAME);
            const DofMap& F_dof_map = F_system.get_dof_map();
            IBTK::FEDataManager::SystemDofMapCache& F_dof_map_cache =
                *d_fe_data_managers[part]->getDofMapCache(FORCE_SYSTEM_NAME);
            IBTK::FEDataManager::SystemDofMapCache& X_dof_map_cache =
                *d_fe_data_managers[part]->getDofMapCache(COORDS_SYSTEM_NAME);
            const FEType F_fe_type = F_dof_map.variable_type(0);

            std::unique_ptr<QBase> qrule = QBase::build(d_default_quad_type[part], dim, d_default_quad_order[part]);
            std::unique_ptr<FEBase> fe = FEBase::build(dim, F_fe_type);
            fe->attach_quadrature_rule(qrule.get());
            const std::vector<double>& JxW = fe->get_JxW();
            const std::vector<RealGradient>& dX_dxi_ref = fe->get_dxyzdxi();
            const std::vector<RealGradient>& dX_deta_ref = fe->get_dxyzdeta();
            const std::vector<std::vector<double>>& dphi_dxi = fe->get_dphidxi();
            const std::vector<std::vector<double>>& dphi_deta = fe->get_dphideta();

            boost::multi_array<double, 2> x_node;
            VectorValue<double> a1, a2;
            SurfaceTangents G_0, G_t, P_tilde;
            std::vector<libMesh::dof_id_type> dof_id_scratch;
            const SurfaceDamageParams& dp = d_damage_trackers[part].params();

            for (auto el_it = mesh.active_local_elements_begin(); el_it != mesh.active_local_elements_end(); ++el_it)
            {
                Elem* const elem = *el_it;
                const auto& F_dof_indices = F_dof_map_cache.dof_indices(elem);
                const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    F_rhs_e[d].resize(static_cast<int>(F_dof_indices[d].size()));
                    F_rhs_e[d].zero();
                }

                fe->reinit(elem);
                IBTK::get_values_for_interpolation(x_node, X_vec, X_dof_indices);

                const unsigned int n_qp = qrule->n_points();
                const std::size_t n_basis = dphi_dxi.size();
                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    IBTK::interpolate(a1, qp, x_node, dphi_dxi);
                    IBTK::interpolate(a2, qp, x_node, dphi_deta);
                    for (unsigned int i = 0; i < 3; ++i)
                    {
                        G_0(i, 0) = dX_dxi_ref[qp](i);
                        G_0(i, 1) = dX_deta_ref[qp](i);
                        G_t(i, 0) = a1(i);
                        G_t(i, 1) = a2(i);
                    }
                    const double mult =
                        d_damage_trackers[part].getElasticMultiplier(static_cast<std::uint64_t>(elem->id()), qp);
                    evalDamagedSurfacePK1(dp, G_0, G_t, P_tilde, mult);

                    const double w_dA0 = JxW[qp];
                    for (std::size_t k = 0; k < n_basis; ++k)
                    {
                        const double g0 = dphi_dxi[k][qp];
                        const double g1 = dphi_deta[k][qp];
                        for (unsigned int i = 0; i < NDIM; ++i)
                        {
                            F_rhs_e[i](static_cast<unsigned int>(k)) -=
                                (P_tilde(i, 0) * g0 + P_tilde(i, 1) * g1) * w_dA0;
                        }
                    }
                }

                for (unsigned int var_n = 0; var_n < NDIM; ++var_n)
                {
                    IBTK::copy_dof_ids_to_vector(var_n, F_dof_indices, dof_id_scratch);
                    F_dof_map.constrain_element_vector(F_rhs_e[var_n], dof_id_scratch);
                    F_rhs_vec.add_vector(F_rhs_e[var_n], dof_id_scratch);
                }
            }

            F_rhs_vec.close();
            d_fe_data_managers[part]->computeL2Projection(
                F_vec, F_rhs_vec, FORCE_SYSTEM_NAME, d_use_consistent_mass_matrix);
        }
    }

    int commitSurfaceDamageStep(double current_time, unsigned int part = 0)
    {
        using namespace libMesh;
        TBOX_ASSERT(part < d_damage_trackers.size());
        if (!d_damage_trackers[part].params().enabled) return 0;

        EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
        MeshBase& mesh = equation_systems->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();
        System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
        IBTK::copy_and_synch(*X_system.solution, *X_system.current_local_solution);
        NumericVector<double>& X_ghost = *X_system.current_local_solution;
        const DofMap& X_dof_map = X_system.get_dof_map();
        const FEType fe_type = X_dof_map.variable_type(0);

        std::unique_ptr<QBase> qrule = QBase::build(d_default_quad_type[part], dim, d_default_quad_order[part]);
        std::unique_ptr<FEBase> fe = FEBase::build(dim, fe_type);
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<RealGradient>& dX_dxi_ref = fe->get_dxyzdxi();
        const std::vector<RealGradient>& dX_deta_ref = fe->get_dxyzdeta();
        const std::vector<std::vector<double>>& dphi_dxi = fe->get_dphidxi();
        const std::vector<std::vector<double>>& dphi_deta = fe->get_dphideta();

        std::unique_ptr<FEBase> fe_nod = FEBase::build(dim, fe_type);
        const std::vector<RealGradient>& dX_dxi_nod = fe_nod->get_dxyzdxi();
        const std::vector<RealGradient>& dX_deta_nod = fe_nod->get_dxyzdeta();
        const std::vector<std::vector<double>>& dphi_dxi_nod = fe_nod->get_dphidxi();
        const std::vector<std::vector<double>>& dphi_deta_nod = fe_nod->get_dphideta();
        const std::vector<libMesh::Point> tri6_ref_pts = {
            libMesh::Point(0.0, 0.0, 0.0), libMesh::Point(1.0, 0.0, 0.0), libMesh::Point(0.0, 1.0, 0.0),
            libMesh::Point(0.5, 0.0, 0.0), libMesh::Point(0.5, 0.5, 0.0), libMesh::Point(0.0, 0.5, 0.0)
        };
        const std::vector<libMesh::Point> tri3_ref_pts = { libMesh::Point(0.0, 0.0, 0.0),
                                                           libMesh::Point(1.0, 0.0, 0.0),
                                                           libMesh::Point(0.0, 1.0, 0.0) };

        std::vector<std::vector<libMesh::dof_id_type>> X_dofs(NDIM);
        boost::multi_array<double, 2> x_node;
        VectorValue<double> a1, a2;
        SurfaceTangents G_0, G_t;

        const double l_d = d_damage_trackers[part].params().nonlocal_length;
        const bool use_nonlocal = (l_d > 0.0);
        if (use_nonlocal)
        {
            d_nonlocal_filters[part].ensureInitialized(
                mesh, fe_type, d_default_quad_type[part], d_default_quad_order[part], l_d);
        }

        std::vector<std::uint64_t> nl_elem_ids;
        std::vector<bool> nl_elem_invalid;
        std::vector<std::vector<double>> nl_y_local_qp;
        if (use_nonlocal)
        {
            const std::size_t n_loc = d_nonlocal_filters[part].numLocalElems();
            nl_elem_ids.reserve(n_loc);
            nl_elem_invalid.reserve(n_loc);
            nl_y_local_qp.reserve(n_loc);
        }

        int loc_folded_elems = 0;
        for (auto el_it = mesh.active_local_elements_begin(); el_it != mesh.active_local_elements_end(); ++el_it)
        {
            Elem* const elem = *el_it;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(elem, X_dofs[d], d);
            }
            IBTK::get_values_for_interpolation(x_node, X_ghost, X_dofs);

            const VectorValue<double> p0(x_node[0][0], x_node[0][1], x_node[0][2]);
            const VectorValue<double> p1(x_node[1][0], x_node[1][1], x_node[1][2]);
            const VectorValue<double> p2(x_node[2][0], x_node[2][1], x_node[2][2]);
            const VectorValue<double> n_sec_raw = (p1 - p0).cross(p2 - p0);
            const double n_sec_norm = n_sec_raw.norm();
            const VectorValue<double> n_sec =
                (n_sec_norm > 1.0e-15) ? (n_sec_raw * (1.0 / n_sec_norm)) : VectorValue<double>(0.0, 0.0, 0.0);
            bool elem_fold_invalid = (n_sec_norm <= 1.0e-15);

            if (!elem_fold_invalid)
            {
                const auto& ref_pts = (elem->n_nodes() >= 6) ? tri6_ref_pts : tri3_ref_pts;
                fe_nod->reinit(elem, &ref_pts);
                for (unsigned int np = 0; np < ref_pts.size(); ++np)
                {
                    IBTK::interpolate(a1, np, x_node, dphi_dxi_nod);
                    IBTK::interpolate(a2, np, x_node, dphi_deta_nod);
                    const VectorValue<double> A1_n(dX_dxi_nod[np](0), dX_dxi_nod[np](1), dX_dxi_nod[np](2));
                    const VectorValue<double> A2_n(dX_deta_nod[np](0), dX_deta_nod[np](1), dX_deta_nod[np](2));
                    const double dA0_n = std::max(1.0e-15, A1_n.cross(A2_n).norm());
                    const double js_nod = (a1.cross(a2) * n_sec) / dA0_n;
                    if (!std::isfinite(js_nod) || js_nod <= 0.0)
                    {
                        elem_fold_invalid = true;
                        break;
                    }
                }
            }
            if (elem_fold_invalid) ++loc_folded_elems;

            fe->reinit(elem);
            const unsigned int n_qp = qrule->n_points();
            const std::uint64_t elem_id = static_cast<std::uint64_t>(elem->id());
            std::vector<double> nl_yq;
            if (use_nonlocal) nl_yq.assign(n_qp, 0.0);
            bool nl_any_invalid = elem_fold_invalid;

            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                IBTK::interpolate(a1, qp, x_node, dphi_dxi);
                IBTK::interpolate(a2, qp, x_node, dphi_deta);
                for (unsigned int i = 0; i < 3; ++i)
                {
                    G_0(i, 0) = dX_dxi_ref[qp](i);
                    G_0(i, 1) = dX_deta_ref[qp](i);
                    G_t(i, 0) = a1(i);
                    G_t(i, 1) = a2(i);
                }
                if (use_nonlocal)
                {
                    double js_raw = 1.0;
                    const double w0_loc =
                        evalUndamagedSurfaceStrainEnergy(d_damage_trackers[part].params(), G_0, G_t, &js_raw);
                    if (elem_fold_invalid || !std::isfinite(w0_loc) || !std::isfinite(js_raw) || js_raw <= 0.0)
                    {
                        nl_any_invalid = true;
                        nl_yq[qp] = 0.0;
                    }
                    else
                    {
                        nl_yq[qp] = w0_loc;
                    }
                }
                else
                {
                    d_damage_trackers[part].updateSlot(elem_id, qp, G_0, G_t, elem_fold_invalid);
                }
            }
            if (use_nonlocal)
            {
                nl_elem_ids.push_back(elem_id);
                nl_elem_invalid.push_back(nl_any_invalid);
                nl_y_local_qp.push_back(std::move(nl_yq));
            }
        }

        if (use_nonlocal)
        {
            std::vector<std::vector<double>> ybar_local_qp;
            d_nonlocal_filters[part].filterQuadratureField(nl_y_local_qp, ybar_local_qp);
            TBOX_ASSERT(nl_elem_ids.size() == ybar_local_qp.size());
            for (std::size_t idx = 0; idx < nl_elem_ids.size(); ++idx)
            {
                const std::uint64_t elem_id = nl_elem_ids[idx];
                const bool inv = nl_elem_invalid[idx];
                const auto& ybar_q = ybar_local_qp[idx];
                for (unsigned int qp = 0; qp < ybar_q.size(); ++qp)
                {
                    d_damage_trackers[part].updateSlotFromEnergy(elem_id, qp, ybar_q[qp], inv);
                }
            }
        }
        d_damage_trackers[part].setCommittedTime(current_time);
        return IBTK::IBTK_MPI::sumReduction(loc_folded_elems);
    }

private:
    void getSurfaceDamageFromRestart(const std::string& object_name)
    {
        auto root_db = SAMRAI::tbox::RestartManager::getManager()->getRootDatabase();
        if (!root_db || !root_db->keyExists(object_name)) return;
        auto db = root_db->getDatabase(object_name);
        for (unsigned int part = 0; part < d_damage_trackers.size(); ++part)
        {
            const std::string key = "SurfaceDamagePart_" + std::to_string(part);
            if (db->keyExists(key))
            {
                d_damage_trackers[part].getFromDatabase(db->getDatabase(key));
            }
        }
    }

    std::vector<SurfaceDamageTracker> d_damage_trackers;
    std::vector<SurfaceNonlocalHelmholtzFilter> d_nonlocal_filters;
};

} // namespace ModelData

using namespace ModelData;

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB3d_ex14.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const std::string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
        const std::string restart_read_dirname = app_initializer->getRestartReadDirectory();
        const int restart_restore_num = app_initializer->getRestartRestoreNumber();

        const double radius = input_db->getDoubleWithDefault("SPHERE_RADIUS", 1.0);
        const int n_subdiv = input_db->getIntegerWithDefault("SPHERE_SUBDIV", 2);
        const double l_box = input_db->getDoubleWithDefault("L", 6.0);
        const libMesh::Point center(0.5 * l_box, 0.5 * l_box, 0.5 * l_box);
        std::unique_ptr<libMesh::ReplicatedMesh> mesh = buildSphereTri6Mesh(init, radius, center, n_subdiv);

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        Pointer<DamagedMembraneIBFESurfaceMethod> ib_method_ops = new DamagedMembraneIBFESurfaceMethod(
            "IBFESurfaceMethod",
            app_initializer->getComponentDatabase("IBFESurfaceMethod"),
            mesh.get(),
            app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
            /*register_for_restart*/ true,
            restart_read_dirname,
            static_cast<unsigned int>(std::max(0, restart_restore_num)));

        SurfaceDamageParams dp;
        dp.enabled = true;
        dp.two_way_coupling = true;
        dp.law = SurfaceDamageLawType::BARTHES_BIESEL_NEO_HOOKEAN;
        dp.G_s = input_db->getDouble("G_S");
        dp.Y_D = input_db->getDouble("Y_D");
        dp.Y_C = input_db->getDouble("Y_C");
        dp.d_cap = input_db->getDoubleWithDefault("D_CAP", 1.0);
        dp.nonlocal_length = input_db->getDouble("NONLOCAL_LENGTH");
        ib_method_ops->setSurfaceDamageParams(dp, 0);

        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);

        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
        if (periodic_shift.min() == 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const std::string bc_name = "u_bc_coefs_" + std::to_string(d);
                const std::string bc_db = "VelocityBcCoefs_" + std::to_string(d);
                u_bc_coefs[d] =
                    new muParserRobinBcCoefs(bc_name, app_initializer->getComponentDatabase(bc_db), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        ib_method_ops->initializeFEEquationSystems();
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        const double loop_time_end = time_integrator->getEndTime();
        int total_folded_elems = 0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();
            const double dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;
            total_folded_elems += ib_method_ops->commitSurfaceDamageStep(loop_time, 0);

            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                ib_method_ops->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
            }
        }

        const SurfaceDamageSummary loc_sum = ib_method_ops->getSurfaceDamageTracker(0).summarizeLocal();
        const double d_eff_max = IBTK_MPI::maxReduction(loc_sum.d_eff_max);
        const int n_damaged = IBTK_MPI::sumReduction(loc_sum.n_damaged_qp);
        const int n_ruptured = IBTK_MPI::sumReduction(loc_sum.n_ruptured_qp);
        const int n_total_qp = IBTK_MPI::sumReduction(loc_sum.n_total_qp);

        Pointer<MemoryDatabase> restart_db = new MemoryDatabase("ex14_restart_check");
        ib_method_ops->putToDatabase(restart_db);
        Pointer<MemoryDatabase> part0_db = restart_db->getDatabase("SurfaceDamagePart_0");
        SurfaceDamageTracker restored;
        restored.configure(dp);
        restored.getFromDatabase(part0_db);
        TBOX_ASSERT(part0_db->getInteger("damage_n_slots") == n_total_qp);

        const SurfaceDamageSummary rest_sum = restored.summarizeLocal();
        const double restart_diff = IBTK_MPI::maxReduction(std::abs(rest_sum.d_eff_max - d_eff_max));
        TBOX_ASSERT(restart_diff <= 1.0e-14);
        TBOX_ASSERT(d_eff_max > 0.0 && n_damaged > 0 && n_ruptured == 0 && total_folded_elems == 0);

        pout << std::fixed << std::setprecision(6) << "EX14_DAMAGE_COMPLETE t=" << loop_time
             << " d_eff_max=" << d_eff_max << " n_damaged=" << n_damaged << " n_total_qp=" << n_total_qp
             << " n_ruptured=" << n_ruptured << " folded_elems=" << total_folded_elems
             << " restart_diff=" << std::scientific << std::setprecision(2) << restart_diff << "\n";

        for (auto* bc : u_bc_coefs) delete bc;
    }
    return 0;
}
