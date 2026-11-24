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
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <IntVector.h>
#include <LoadBalancer.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/BrinkmanPenalizationStrategy.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
void
set_constant_density(int rho_idx,
                     Pointer<hier::Variable<NDIM>> /*rho_var*/,
                     Pointer<HierarchyMathOps> hier_math_ops,
                     const int /*cycle_num*/,
                     const double /*time*/,
                     const double /*current_time*/,
                     const double /*new_time*/,
                     void* ctx)
{
    const double rho_val = *static_cast<const double*>(ctx);
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> rho_data = patch->getPatchData(rho_idx);
            rho_data->fillAll(rho_val);
        }
    }
    return;
}

class BoxBrinkmanPenalization : public BrinkmanPenalizationStrategy
{
public:
    BoxBrinkmanPenalization(std::string object_name,
                            Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                            const double penalty_factor)
        : BrinkmanPenalizationStrategy(std::move(object_name), /*register_for_restart*/ false),
          d_patch_hierarchy(std::move(patch_hierarchy))
    {
        d_penalty_factor = penalty_factor;
    }

    void computeBrinkmanVelocity(int /*u_idx*/, double /*time*/, int /*cycle_num*/) override
    {
        return;
    }

    void demarcateBrinkmanZone(int u_idx, double /*time*/, int /*cycle_num*/) override
    {
        const int coarsest_ln = 0;
        const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const hier::Index<NDIM>& patch_lower = patch_box.lower();
                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                const double* const x_lo = pgeom->getXLower();
                const double* const dx = pgeom->getDx();
                Pointer<SideData<NDIM, double>> L_data = patch->getPatchData(u_idx);
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                    {
                        const SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);
                        double coord[NDIM];
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            const double offset = (d == axis) ? 0.0 : 0.5;
                            coord[d] = x_lo[d] + dx[d] * (static_cast<double>(s_i(d) - patch_lower(d)) + offset);
                        }
                        if (coord[0] >= 0.25 && coord[0] <= 0.75 && coord[1] >= 0.25 && coord[1] <= 0.75)
                        {
                            (*L_data)(s_i) += d_penalty_factor;
                        }
                    }
                }
            }
        }
        return;
    }

private:
    Pointer<PatchHierarchy<NDIM>> d_patch_hierarchy;
};

class PressureCheckIntegrator : public INSVCStaggeredConservativeHierarchyIntegrator
{
public:
    PressureCheckIntegrator(const std::string& object_name,
                            Pointer<Database> input_db,
                            Pointer<PatchHierarchy<NDIM>> hierarchy,
                            const double rel_res_tol)
        : INSVCStaggeredConservativeHierarchyIntegrator(object_name, input_db, /*register_for_restart*/ false),
          d_hier_ref(std::move(hierarchy)),
          d_rel_res_tol(rel_res_tol)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext(object_name + "::PRESSURE_CHECK");
        d_phi_var = new CellVariable<NDIM, double>(object_name + "::phi");
        d_rhs_var = new CellVariable<NDIM, double>(object_name + "::rhs");
        d_res_var = new CellVariable<NDIM, double>(object_name + "::res");
        d_phi_idx = var_db->registerVariableAndContext(d_phi_var, ctx, IntVector<NDIM>(1));
        d_rhs_idx = var_db->registerVariableAndContext(d_rhs_var, ctx, IntVector<NDIM>(1));
        d_res_idx = var_db->registerVariableAndContext(d_res_var, ctx, IntVector<NDIM>(1));
    }

    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override
    {
        const int finest_ln = d_hier_ref->getFinestLevelNumber();
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hier_ref->getPatchLevel(ln);
            if (!level->checkAllocated(d_phi_idx))
            {
                level->allocatePatchData(d_phi_idx, current_time);
            }
            if (!level->checkAllocated(d_rhs_idx))
            {
                level->allocatePatchData(d_rhs_idx, current_time);
            }
            if (!level->checkAllocated(d_res_idx))
            {
                level->allocatePatchData(d_res_idx, current_time);
            }
        }

        const int cw_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
        Pointer<SAMRAIVectorReal<NDIM, double>> phi_vec =
            new SAMRAIVectorReal<NDIM, double>("phi_vec", d_hier_ref, 0, finest_ln);
        phi_vec->addComponent(d_phi_var, d_phi_idx, cw_idx);
        Pointer<SAMRAIVectorReal<NDIM, double>> rhs_vec =
            new SAMRAIVectorReal<NDIM, double>("rhs_vec", d_hier_ref, 0, finest_ln);
        rhs_vec->addComponent(d_rhs_var, d_rhs_idx, cw_idx);
        Pointer<SAMRAIVectorReal<NDIM, double>> res_vec =
            new SAMRAIVectorReal<NDIM, double>("res_vec", d_hier_ref, 0, finest_ln);
        res_vec->addComponent(d_res_var, d_res_idx, cw_idx);

        constexpr double pi = 3.14159265358979323846;
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hier_ref->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> phi_data = patch->getPatchData(d_phi_idx);
                Pointer<CellData<NDIM, double>> rhs_data = patch->getPatchData(d_rhs_idx);
                phi_data->fillAll(0.0);
                rhs_data->fillAll(0.0);
                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                const double* const x_lo = pgeom->getXLower();
                const double* const dx = pgeom->getDx();
                const Box<NDIM>& pbox = patch->getBox();
                const hier::Index<NDIM>& plo = pbox.lower();
                for (Box<NDIM>::Iterator b(pbox); b; b++)
                {
                    const hier::Index<NDIM>& c = b();
                    const double x = x_lo[0] + (static_cast<double>(c(0) - plo(0)) + 0.5) * dx[0];
                    const double y = x_lo[1] + (static_cast<double>(c(1) - plo(1)) + 0.5) * dx[1];
                    (*rhs_data)(c) = 1.0 + 0.5 * std::sin(pi * x) * std::cos(pi * y);
                }
            }
        }

        d_pressure_solver->setHomogeneousBc(true);
        d_pressure_solver->solveSystem(*phi_vec, *rhs_vec);

        PoissonSpecifications current_p_spec("current_p_spec");
        current_p_spec.setCZero();
        current_p_spec.setDPatchDataId(d_pressure_D_idx);
        CCLaplaceOperator laplace_op("current_laplace", nullptr, true);
        laplace_op.setPoissonSpecifications(current_p_spec);
        laplace_op.setPhysicalBcCoef(d_Phi_bc_coef.get());
        laplace_op.setHomogeneousBc(true);
        laplace_op.setSolutionTime(d_pressure_solver->getSolutionTime());
        laplace_op.setTimeInterval(current_time, new_time);
        laplace_op.initializeOperatorState(*phi_vec, *res_vec);
        laplace_op.apply(*phi_vec, *res_vec);
        laplace_op.deallocateOperatorState();

        HierarchyCellDataOpsReal<NDIM, double> cc_ops(d_hier_ref, 0, finest_ln);
        cc_ops.axpy(d_res_idx, -1.0, d_res_idx, d_rhs_idx, false);
        const double norm_r = cc_ops.L2Norm(d_res_idx, cw_idx);
        const double norm_b = cc_ops.L2Norm(d_rhs_idx, cw_idx);
        const double rel_res = norm_r / ((norm_b > 1.0e-30) ? norm_b : 1.0e-30);
        const double phi_l2 = cc_ops.L2Norm(d_phi_idx, cw_idx);
        const double phi_max = cc_ops.maxNorm(d_phi_idx);
        const size_t step_idx = d_step_phi_l2.size();
        if (!std::isfinite(rel_res) || !std::isfinite(phi_l2) || !std::isfinite(phi_max))
        {
            TBOX_ERROR("Nonfinite pressure residual or norm after step " << step_idx << ".\n");
        }
        if (rel_res > d_rel_res_tol)
        {
            TBOX_ERROR("Pressure subdomain solver residual exceeded tolerance after step "
                       << step_idx << ": rel_res = " << rel_res << ", tol = " << d_rel_res_tol << ".\n");
        }
        d_step_phi_l2.push_back(phi_l2);
        d_step_phi_max.push_back(phi_max);

        INSVCStaggeredConservativeHierarchyIntegrator::postprocessIntegrateHierarchy(
            current_time, new_time, skip_synchronize_new_state_data, num_cycles);
        return;
    }

    const std::vector<double>& getStepPhiL2() const
    {
        return d_step_phi_l2;
    }

    const std::vector<double>& getStepPhiMax() const
    {
        return d_step_phi_max;
    }

private:
    Pointer<PatchHierarchy<NDIM>> d_hier_ref;
    double d_rel_res_tol = 1.0e-8;
    Pointer<CellVariable<NDIM, double>> d_phi_var;
    Pointer<CellVariable<NDIM, double>> d_rhs_var;
    Pointer<CellVariable<NDIM, double>> d_res_var;
    int d_phi_idx = -1;
    int d_rhs_idx = -1;
    int d_res_idx = -1;
    std::vector<double> d_step_phi_l2;
    std::vector<double> d_step_phi_max;
};
} // namespace

// Verify that changing the time step size with Brinkman penalization updates the pressure subdomain solver operator.
int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "vc_brinkman_dt_pressure_init.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

        std::vector<std::unique_ptr<RobinBcCoefStrategy<NDIM>>> u_bc_storage;
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            u_bc_storage.push_back(std::make_unique<muParserRobinBcCoefs>(
                "u_bc_" + std::to_string(d),
                app_initializer->getComponentDatabase("VelocityBcCoefs_" + std::to_string(d)),
                grid_geometry));
            u_bc_coefs[d] = u_bc_storage.back().get();
        }
        std::unique_ptr<RobinBcCoefStrategy<NDIM>> rho_bc_coef = std::make_unique<muParserRobinBcCoefs>(
            "rho_bc", app_initializer->getComponentDatabase("DensityBcCoefs"), grid_geometry);

        Pointer<PressureCheckIntegrator> time_integrator = new PressureCheckIntegrator(
            "INSVCStaggeredConservativeHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSVCStaggeredConservativeHierarchyIntegrator"),
            patch_hierarchy,
            input_db->getDouble("rel_res_tol"));
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

        Pointer<hier::Variable<NDIM>> rho_var = new SideVariable<NDIM, double>("rho");
        double rho_val = input_db->getDouble("rho_val");
        time_integrator->registerMassDensityVariable(rho_var);
        time_integrator->registerResetFluidDensityFcn(&set_constant_density, &rho_val);

        Pointer<BoxBrinkmanPenalization> brinkman = new BoxBrinkmanPenalization(
            "BoxBrinkmanPenalization", patch_hierarchy, input_db->getDouble("brinkman_penalty_factor"));
        time_integrator->registerBrinkmanPenalizationStrategy(brinkman);

        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        time_integrator->registerVelocityInitialConditions(u_init);
        time_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        time_integrator->registerMassDensityBoundaryConditions(rho_bc_coef.get());

        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        const std::array<double, 3> step_dts = {
            input_db->getDouble("dt_0"),
            input_db->getDouble("dt_1"),
            input_db->getDouble("dt_2"),
        };
        for (const double step_dt : step_dts)
        {
            time_integrator->advanceHierarchy(step_dt);
        }

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");
            const std::vector<double>& step_phi_l2 = time_integrator->getStepPhiL2();
            const std::vector<double>& step_phi_max = time_integrator->getStepPhiMax();
            for (size_t k = 0; k < step_phi_l2.size(); ++k)
            {
                out << "step_" << k << ": phi_l2 = " << std::setprecision(10) << step_phi_l2[k]
                    << ", phi_max = " << step_phi_max[k] << '\n';
            }
        }
    }
    return 0;
}
