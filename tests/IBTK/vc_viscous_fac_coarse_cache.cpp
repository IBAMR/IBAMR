// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files

#include <SAMRAI_config.h>

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <BoundaryBox.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <GriddingAlgorithm.h>
#include <HierarchySideDataOpsReal.h>
#include <LoadBalancer.h>
#include <RobinBcCoefStrategy.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScKrylovPoissonSolver.h>
#include <ibtk/PoissonFACPreconditioner.h>
#include <ibtk/SCPoissonSolverManager.h>
#include <ibtk/VCSCViscousOpPointRelaxationFACOperator.h>
#include <ibtk/VCSCViscousOperator.h>
#include <ibtk/VCSCViscousPETScLevelSolver.h>
#include <ibtk/ibtk_enums.h>
#include <ibtk/muParserCartGridFunction.h>

// Set up application namespace declarations
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include <ibtk/app_namespaces.h>

namespace
{
class CycleRobinBcCoef : public RobinBcCoefStrategy<NDIM>
{
public:
    CycleRobinBcCoef(const double a, const double b, const double g_scale) : d_a(a), d_b(b), d_g_scale(g_scale)
    {
    }

    void setParameters(const double a, const double b, const double g_scale)
    {
        d_a = a;
        d_b = b;
        d_g_scale = g_scale;
    }

    void setBcCoefs(Pointer<ArrayData<NDIM, double>>& acoef_data,
                    Pointer<ArrayData<NDIM, double>>& bcoef_data,
                    Pointer<ArrayData<NDIM, double>>& gcoef_data,
                    const Pointer<hier::Variable<NDIM>>& /*variable*/,
                    const Patch<NDIM>& patch,
                    const BoundaryBox<NDIM>& bdry_box,
                    double /*fill_time*/ = 0.0) const override
    {
        if (acoef_data)
        {
            acoef_data->fillAll(d_a);
        }
        if (bcoef_data)
        {
            bcoef_data->fillAll(d_b);
        }
        if (gcoef_data)
        {
            const Box<NDIM>& patch_box = patch.getBox();
            const hier::Index<NDIM>& patch_lower = patch_box.lower();
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
            const double* const x_lower = pgeom->getXLower();
            const double* const dx = pgeom->getDx();
            const unsigned int normal_axis = bdry_box.getLocationIndex() / 2;
            for (Box<NDIM>::Iterator it(gcoef_data->getBox()); it; it++)
            {
                const hier::Index<NDIM>& idx = it();
                double coord[NDIM];
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    const double shift = (d == normal_axis) ? 0.0 : 0.5;
                    coord[d] = x_lower[d] + dx[d] * (static_cast<double>(idx(d) - patch_lower(d)) + shift);
                }
                constexpr double pi = 3.14159265358979323846;
                (*gcoef_data)(idx, 0) =
                    d_g_scale * (1.0 + 0.4 * std::cos(2.0 * pi * coord[0]) * std::sin(2.0 * pi * coord[1]));
            }
        }
    }

    IntVector<NDIM> numberOfExtensionsFillable() const override
    {
        return IntVector<NDIM>(2);
    }

private:
    double d_a = 1.0;
    double d_b = 0.0;
    double d_g_scale = 0.0;
};

struct SolverBundle
{
    Pointer<VCSCViscousOperator> op;
    Pointer<VCSCViscousOpPointRelaxationFACOperator> fac_op;
    Pointer<PoissonFACPreconditioner> pc;
    Pointer<PETScKrylovPoissonSolver> solver;
};

SolverBundle
make_solver(const std::string& name,
            const std::string& prefix,
            Pointer<Database> solver_db,
            Pointer<Database> precond_db)
{
    SolverBundle bundle;
    bundle.op = new VCSCViscousOperator(name + "::op", false);
    bundle.fac_op = new VCSCViscousOpPointRelaxationFACOperator(name + "::fac_op", precond_db, prefix + "pc_");
    bundle.pc = new PoissonFACPreconditioner(name + "::pc", bundle.fac_op, precond_db, prefix + "pc_");
    bundle.solver = new PETScKrylovPoissonSolver(name + "::solver", solver_db, prefix);
    bundle.solver->setOperator(bundle.op);
    bundle.solver->setPreconditioner(bundle.pc);
    return bundle;
}
} // namespace

// Exercise coarsest-level solver reuse and invalidation in VCSCViscousOpPointRelaxationFACOperator across
// deallocate/initialize resets while changing inhomogeneous BC values, Robin coefficients, and viscosity interpolation.
int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "vc_viscous_fac_coarse_cache.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        SCPoissonSolverManager::getManager()->registerSolverFactoryFunction(
            "VC_VELOCITY_PETSC_LEVEL_SOLVER", VCSCViscousPETScLevelSolver::allocate_solver);

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        Pointer<SideVariable<NDIM, double>> u_cached_var = new SideVariable<NDIM, double>("u_cached");
        Pointer<SideVariable<NDIM, double>> u_fresh_var = new SideVariable<NDIM, double>("u_fresh");
        Pointer<SideVariable<NDIM, double>> u_diff_var = new SideVariable<NDIM, double>("u_diff");
        Pointer<SideVariable<NDIM, double>> f_var = new SideVariable<NDIM, double>("f");
        Pointer<NodeVariable<NDIM, double>> mu_var = new NodeVariable<NDIM, double>("mu");

        const int u_cached_idx = var_db->registerVariableAndContext(u_cached_var, ctx, IntVector<NDIM>(1));
        const int u_fresh_idx = var_db->registerVariableAndContext(u_fresh_var, ctx, IntVector<NDIM>(1));
        const int u_diff_idx = var_db->registerVariableAndContext(u_diff_var, ctx, IntVector<NDIM>(1));
        const int f_idx = var_db->registerVariableAndContext(f_var, ctx, IntVector<NDIM>(1));
        const int mu_idx = var_db->registerVariableAndContext(mu_var, ctx, IntVector<NDIM>(1));

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        int level_number = 0;
        bool done = false;
        while (!done && gridding_algorithm->levelCanBeRefined(level_number))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, 1);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_cached_idx, 0.0);
            level->allocatePatchData(u_fresh_idx, 0.0);
            level->allocatePatchData(u_diff_idx, 0.0);
            level->allocatePatchData(f_idx, 0.0);
            level->allocatePatchData(mu_idx, 0.0);
        }

        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int h_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> u_cached_vec("u_cached", patch_hierarchy, 0, finest_ln);
        SAMRAIVectorReal<NDIM, double> u_fresh_vec("u_fresh", patch_hierarchy, 0, finest_ln);
        SAMRAIVectorReal<NDIM, double> f_vec("f", patch_hierarchy, 0, finest_ln);
        u_cached_vec.addComponent(u_cached_var, u_cached_idx, h_sc_idx);
        u_fresh_vec.addComponent(u_fresh_var, u_fresh_idx, h_sc_idx);
        f_vec.addComponent(f_var, f_idx, h_sc_idx);

        muParserCartGridFunction f_fcn("f", app_initializer->getComponentDatabase("f"), grid_geometry);
        muParserCartGridFunction mu_fcn("mu", app_initializer->getComponentDatabase("mu"), grid_geometry);
        f_fcn.setDataOnPatchHierarchy(f_idx, f_var, patch_hierarchy, 0.0);
        mu_fcn.setDataOnPatchHierarchy(mu_idx, mu_var, patch_hierarchy, 0.0);

        std::vector<HierarchyGhostCellInterpolation::InterpolationTransactionComponent> mu_tx(1);
        mu_tx[0] = HierarchyGhostCellInterpolation::InterpolationTransactionComponent(
            mu_idx,
            "LINEAR_REFINE",
            false,
            "CONSTANT_COARSEN",
            "LINEAR",
            false,
            nullptr,
            Pointer<VariableFillPattern<NDIM>>(nullptr));
        Pointer<HierarchyGhostCellInterpolation> mu_fill = new HierarchyGhostCellInterpolation();
        mu_fill->initializeOperatorState(mu_tx, patch_hierarchy);
        mu_fill->setHomogeneousBc(false);
        mu_fill->fillData(0.0);

        PoissonSpecifications vc_spec("vc_spec");
        vc_spec.setCConstant(input_db->getDouble("C"));
        vc_spec.setDPatchDataId(mu_idx);

        const double cycle_dt = input_db->getDouble("cycle_dt");
        Pointer<Database> solver_db = input_db->getDatabase("solver_db");
        Pointer<Database> cached_precond_db = input_db->getDatabase("cached_precond_db");
        Pointer<Database> fresh_precond_db = input_db->getDatabase("fresh_precond_db");

        SolverBundle cached = make_solver("cached", "cached_", solver_db, cached_precond_db);
        SolverBundle fresh = make_solver("fresh", "fresh_", solver_db, fresh_precond_db);

        std::vector<CycleRobinBcCoef> bc_storage;
        bc_storage.reserve(NDIM);
        std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs(NDIM, nullptr);
        for (int d = 0; d < NDIM; ++d)
        {
            bc_storage.emplace_back(1.0, 0.2, 1.0);
            bc_coefs[d] = &bc_storage.back();
        }

        struct CycleSpec
        {
            double a;
            double b;
            double g_scale;
            double time;
            IBTK::VCInterpType mu_interp;
        };

        const std::array<CycleSpec, 5> cycles = { {
            { 1.0, 0.2, 1.0, 0.0, IBTK::VC_HARMONIC_INTERP },
            { 1.0, 0.2, 2.5, 0.5, IBTK::VC_HARMONIC_INTERP },
            { 1.0, 1.0e-8, 2.5, 1.0, IBTK::VC_HARMONIC_INTERP },
            { 1.0, 5.0e-7, 2.5, 1.5, IBTK::VC_HARMONIC_INTERP },
            { 1.0, 5.0e-7, 2.5, 2.0, IBTK::VC_AVERAGE_INTERP },
        } };

        const std::array<long, 5> expected_reuse = { 0, 1, 1, 1, 1 };
        const std::array<long, 5> expected_reinit = { 1, 1, 2, 3, 4 };

        HierarchySideDataOpsReal<NDIM, double> sc_ops(patch_hierarchy, 0, finest_ln);
        std::ofstream out;
        if (IBTK_MPI::getRank() == 0)
        {
            out.open("output");
        }

        for (size_t k = 0; k < cycles.size(); ++k)
        {
            const CycleSpec& spec = cycles[k];
            for (int d = 0; d < NDIM; ++d)
            {
                bc_storage[d].setParameters(spec.a, spec.b, spec.g_scale);
            }

            auto configure_and_solve = [&](SolverBundle& bnd, SAMRAIVectorReal<NDIM, double>& u_vec, const int u_idx)
            {
                bnd.solver->deallocateSolverState();
                bnd.op->setDPatchDataInterpolationType(spec.mu_interp);
                bnd.fac_op->setDPatchDataInterpolationType(spec.mu_interp);
                auto* coarse_vc =
                    dynamic_cast<VCSCViscousPETScLevelSolver*>(bnd.fac_op->getCoarseSolver().getPointer());
                if (!coarse_vc)
                {
                    TBOX_ERROR("Expected VCSCViscousPETScLevelSolver coarsest-level solver.\n");
                }
                coarse_vc->setViscosityInterpolationType(spec.mu_interp);

                bnd.pc->setSolutionTime(spec.time);
                bnd.pc->setTimeInterval(spec.time, spec.time + cycle_dt);
                bnd.solver->setSolutionTime(spec.time);
                bnd.solver->setTimeInterval(spec.time, spec.time + cycle_dt);
                bnd.solver->setPoissonSpecifications(vc_spec);
                bnd.solver->setPhysicalBcCoefs(bc_coefs);
                bnd.solver->setHomogeneousBc(false);

                sc_ops.setToScalar(u_idx, 0.0, false);
                bnd.solver->initializeSolverState(u_vec, f_vec);
                bnd.solver->solveSystem(u_vec, f_vec);
            };

            configure_and_solve(cached, u_cached_vec, u_cached_idx);
            configure_and_solve(fresh, u_fresh_vec, u_fresh_idx);

            sc_ops.subtract(u_diff_idx, u_cached_idx, u_fresh_idx);
            const double diff_max = sc_ops.maxNorm(u_diff_idx);
            const double u_l2 = sc_ops.L2Norm(u_cached_idx, h_sc_idx);
            const double u_max = sc_ops.maxNorm(u_cached_idx);
            if (!std::isfinite(diff_max) || !std::isfinite(u_l2) || !std::isfinite(u_max))
            {
                TBOX_ERROR("Nonfinite norm observed at cycle " << k << ".\n");
            }
            if (diff_max > 1.0e-10)
            {
                TBOX_ERROR("Cached and uncached solutions diverged at cycle " << k << ": diff_max=" << diff_max
                                                                              << ".\n");
            }
            const long reuse_count = cached.fac_op->getCoarseSolverReuseCount();
            const long reinit_count = cached.fac_op->getCoarseSolverReinitCount();
            if (reuse_count != expected_reuse[k] || reinit_count != expected_reinit[k])
            {
                TBOX_ERROR("Unexpected coarse solver cache counts at cycle " << k << ": reuse=" << reuse_count
                                                                             << " reinit=" << reinit_count << ".\n");
            }
            if (fresh.fac_op->getCoarseSolverReuseCount() != 0 || fresh.fac_op->getCoarseSolverReinitCount() != 0)
            {
                TBOX_ERROR("Uncached solver reported nonzero reuse/reinit counts at cycle " << k << ".\n");
            }

            if (IBTK_MPI::getRank() == 0)
            {
                out << "cycle_" << k << ": reuse = " << reuse_count << ", reinit = " << reinit_count
                    << ", u_l2 = " << std::setprecision(10) << u_l2 << ", u_max = " << u_max << '\n';
            }
        }
    }
    return 0;
}
