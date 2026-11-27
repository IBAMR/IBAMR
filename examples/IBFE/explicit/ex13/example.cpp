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

// Minimal 2D IBFE explicit viscoelastic Standard Linear Solid (SLS) example
// demonstrating:
//   1. Registration via IBFEMethod::registerPK1StressFunction using fixed
//      non-adaptive QGAUSS FIFTH quadrature;
//   2. Nearest-reference-coordinate slot lookup (findNearestSlotIndex) with
//      centroid fallback for non-Gauss postprocessing or helper calls;
//   3. Step-end committed internal-variable advance after advanceHierarchy(dt).

#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/libmesh_utilities.h>

#include <tbox/MemoryDatabase.h>

#include <libmesh/equation_systems.h>
#include <libmesh/fe_base.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/quadrature.h>

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
#include <memory>
#include <unordered_map>
#include <vector>

#include "IBFEViscoelasticity.h"

#include <ibamr/app_namespaces.h>

namespace ModelData
{
using namespace IBAMR::FEViscoelasticity;

static constexpr double G_INF = 1.0;
static constexpr double G_SLS = 2.0;
static constexpr double TAU_R = 0.05;
static constexpr double GAMMA_0 = 0.12;
static const QuadratureType QP_TYPE = QGAUSS;
static const Order QP_ORDER = FIFTH;

struct ViscoelasticContext
{
    GeneralizedMaxwellHistoryStore store;
    std::unordered_map<std::uint64_t, std::vector<std::array<double, 3>>> elem_qp_coords;
    std::unordered_map<std::uint64_t, double> elem_tol2;
    unsigned int qp_hits = 0U;
    unsigned int qp_fallbacks = 0U;
};

void
coordinate_mapping_function(libMesh::Point& X, const libMesh::Point& s, void* /*ctx*/)
{
    X(0) = s(0) + GAMMA_0 * (s(1) - 0.5);
    X(1) = s(1);
    return;
}

void
compute_equilibrium_pk1(TensorValue<double>& P_eq, const TensorValue<double>& FF, const double G_shear)
{
    SymmetricTensor3 S_iso, C, C_inv;
    computeIsochoric2PK(FF, G_shear, S_iso, C, C_inv);
    pushForward2PKToPK1(FF, S_iso, P_eq);
    return;
}

void
pk1_viscoelastic_stress(TensorValue<double>& PP,
                        const TensorValue<double>& FF,
                        const libMesh::Point& /*X*/,
                        const libMesh::Point& s,
                        Elem* const elem,
                        const std::vector<const std::vector<double>*>& /*var_data*/,
                        const std::vector<const std::vector<VectorValue<double>>*>& /*grad_var_data*/,
                        const double data_time,
                        void* const ctx)
{
    auto* const ve = static_cast<ViscoelasticContext*>(ctx);
    TBOX_ASSERT(ve != nullptr);

    // Baseline hyperelastic equilibrium branch.
    compute_equilibrium_pk1(PP, FF, G_INF);
    if (elem == nullptr)
    {
        ++ve->qp_fallbacks;
        return;
    }

    const std::uint64_t eid = static_cast<std::uint64_t>(elem->id());
    const auto it_coords = ve->elem_qp_coords.find(eid);
    const auto it_tol = ve->elem_tol2.find(eid);
    int qp = -1;
    if (it_coords != ve->elem_qp_coords.end() && it_tol != ve->elem_tol2.end())
    {
        qp = GeneralizedMaxwellHistoryStore::findNearestSlotIndex(it_coords->second, s(0), s(1), s(2), it_tol->second);
    }

    if (qp >= 0)
    {
        ++ve->qp_hits;
        TensorValue<double> P_branch;
        ve->store.evaluatePK1(P_branch, eid, static_cast<unsigned int>(qp), FF, data_time);
        PP += P_branch;
    }
    else
    {
        // Fallback for non-Gauss evaluation coordinates (e.g., centroid postprocessing).
        ++ve->qp_fallbacks;
    }
    return;
}

void
initialize_viscoelastic_slots(ViscoelasticContext& ve,
                              const MeshBase& mesh,
                              EquationSystems* const equation_systems,
                              const std::string& coords_sys_name)
{
    ve.store.configure({ ViscoelasticBranchParams{ G_SLS, TAU_R, InitialBranchState::LOADED } });
    ve.elem_qp_coords.clear();
    ve.elem_tol2.clear();

    System& X_system = equation_systems->get_system<System>(coords_sys_name);
    IBTK::copy_and_synch(*X_system.solution, *X_system.current_local_solution);
    NumericVector<double>* const X_ghost_vec = X_system.current_local_solution.get();
    DofMap& X_dof_map = X_system.get_dof_map();

    std::unique_ptr<FEBase> fe(FEBase::build(NDIM, X_dof_map.variable_type(0)));
    std::unique_ptr<QBase> qrule = QBase::build(QP_TYPE, NDIM, QP_ORDER);
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<libMesh::Point>& q_point = fe->get_xyz();
    const std::vector<std::vector<VectorValue<double>>>& dphi = fe->get_dphi();

    std::vector<std::vector<dof_id_type>> X_dof_indices(NDIM);
    boost::multi_array<double, 2> X_node;
    TensorValue<double> FF;

    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const std::uint64_t eid = static_cast<std::uint64_t>(elem->id());
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[d], d);
        }
        get_values_for_interpolation(X_node, *X_ghost_vec, X_dof_indices);

        const double h_e = elem->hmax();
        ve.elem_tol2[eid] = 1.0e-12 * std::max(1.0e-30, h_e * h_e);
        const unsigned int n_qp = qrule->n_points();
        auto& coords = ve.elem_qp_coords[eid];
        coords.resize(n_qp);
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            coords[qp] = { q_point[qp](0), q_point[qp](1), q_point[qp](2) };
            jacobian(FF, qp, X_node, dphi);
            ve.store.initializeSlot(eid, qp, FF);
        }
    }
    return;
}

void
commit_viscoelastic_step(ViscoelasticContext& ve,
                         const MeshBase& mesh,
                         EquationSystems* const equation_systems,
                         const std::string& coords_sys_name,
                         const double new_time)
{
    System& X_system = equation_systems->get_system<System>(coords_sys_name);
    IBTK::copy_and_synch(*X_system.solution, *X_system.current_local_solution);
    NumericVector<double>* const X_ghost_vec = X_system.current_local_solution.get();
    DofMap& X_dof_map = X_system.get_dof_map();

    std::unique_ptr<FEBase> fe(FEBase::build(NDIM, X_dof_map.variable_type(0)));
    std::unique_ptr<QBase> qrule = QBase::build(QP_TYPE, NDIM, QP_ORDER);
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<std::vector<VectorValue<double>>>& dphi = fe->get_dphi();

    std::vector<std::vector<dof_id_type>> X_dof_indices(NDIM);
    boost::multi_array<double, 2> X_node;
    TensorValue<double> FF;

    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const std::uint64_t eid = static_cast<std::uint64_t>(elem->id());
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[d], d);
        }
        get_values_for_interpolation(X_node, *X_ghost_vec, X_dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            jacobian(FF, qp, X_node, dphi);
            ve.store.commitQuadraturePoint(eid, qp, FF, new_time);
        }
    }
    ve.store.finishCommitStep(new_time);
    return;
}
} // namespace ModelData

using namespace ModelData;

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        Mesh mesh(init.comm(), NDIM);
        MeshTools::Generation::build_square(mesh, 8, 4, 0.25, 0.75, 0.375, 0.625, QUAD9);
        mesh.prepare_for_use();

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           &mesh,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);

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

        ViscoelasticContext ve_ctx;
        ib_method_ops->initializeFEEquationSystems();
        ib_method_ops->registerInitialCoordinateMappingFunction(
            IBFEMethod::CoordinateMappingFcnData(coordinate_mapping_function, nullptr));
        ib_method_ops->registerPK1StressFunction(
            IBFEMethod::PK1StressFcnData(pk1_viscoelastic_stress, {}, &ve_ctx, QP_TYPE, QP_ORDER));

        EquationSystems* equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();
        ib_method_ops->initializeFEData();
        initialize_viscoelastic_slots(ve_ctx, mesh, equation_systems, ib_method_ops->getCurrentCoordinatesSystemName());
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        double loop_time = time_integrator->getIntegratorTime();
        const double loop_time_end = time_integrator->getEndTime();
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            const double dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;
            commit_viscoelastic_step(
                ve_ctx, mesh, equation_systems, ib_method_ops->getCurrentCoordinatesSystemName(), loop_time);
        }

        Pointer<MemoryDatabase> db = new MemoryDatabase("ex13_restart_check");
        ve_ctx.store.putToDatabase(db);
        GeneralizedMaxwellHistoryStore restored;
        restored.getFromDatabase(db);
        const unsigned int local_slots = static_cast<unsigned int>(ve_ctx.store.numSlots());
        const unsigned int total_slots = IBTK_MPI::sumReduction(local_slots);
        TBOX_ASSERT(restored.numSlots() == static_cast<size_t>(total_slots));

        double restart_diff = 0.0;
        const auto el_it0 = mesh.active_local_elements_begin();
        if (el_it0 != mesh.active_local_elements_end())
        {
            const std::uint64_t eid0 = static_cast<std::uint64_t>((*el_it0)->id());
            TensorValue<double> F_probe, P_orig, P_rest;
            F_probe.zero();
            F_probe(0, 0) = 1.0;
            F_probe(1, 1) = 1.0;
            F_probe(2, 2) = 1.0;
            F_probe(0, 1) = GAMMA_0;
            ve_ctx.store.evaluatePK1(P_orig, eid0, 0U, F_probe, loop_time);
            restored.evaluatePK1(P_rest, eid0, 0U, F_probe, loop_time);
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                for (unsigned int j = 0; j < NDIM; ++j)
                {
                    restart_diff = std::max(restart_diff, std::abs(P_orig(i, j) - P_rest(i, j)));
                }
            }
        }
        restart_diff = IBTK_MPI::maxReduction(restart_diff);
        TBOX_ASSERT(restart_diff <= 1.0e-14);

        const unsigned int total_hits = IBTK_MPI::sumReduction(ve_ctx.qp_hits);
        const unsigned int total_fallbacks = IBTK_MPI::sumReduction(ve_ctx.qp_fallbacks);
        pout << std::fixed << std::setprecision(6) << "EX13_SLS_COMPLETE t=" << loop_time << " qp_hits=" << total_hits
             << " qp_fallbacks=" << total_fallbacks << " restart_diff=" << std::scientific << std::setprecision(2)
             << restart_diff << "\n";
        TBOX_ASSERT(total_hits > 0U);
    }
    return 0;
}
