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

#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/INSStaggeredPressureBcCoef.h>
#include <ibamr/INSStaggeredVelocityBcCoef.h>
#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <LoadBalancer.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
// Smooth, nonzero values of the exact solution.
double
exact_value(const hier::Index<NDIM>& i, const int component)
{
    double arg = 0.8 * component;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        arg += (1.7 + 1.2 * d) * i(d);
    }
    return std::sin(arg);
}

// Apply the staggered Stokes operator to an exact solution with inhomogeneous boundary data to obtain a right-hand
// side, solve with the level solver from a zero initial guess, and compare the solution with the exact solution.
void
check_level_solver(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                   Pointer<INSHierarchyIntegrator> ins_integrator,
                   const vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                   Pointer<Database> solver_db)
{
    const double solution_time = 0.5;
    const double mu = ins_integrator->getStokesSpecifications()->getMu();

    // Configure the integrator's boundary condition objects.
    const vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs = ins_integrator->getVelocityBoundaryConditions();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto U_bc_coef = dynamic_cast<INSStaggeredVelocityBcCoef*>(U_bc_coefs[d]);
        U_bc_coef->setStokesSpecifications(ins_integrator->getStokesSpecifications());
        U_bc_coef->setPhysicalBcCoefs(u_bc_coefs);
        U_bc_coef->setSolutionTime(solution_time);
    }
    RobinBcCoefStrategy<NDIM>* P_bc_coef = ins_integrator->getPressureBoundaryConditions();
    auto P_ins_bc_coef = dynamic_cast<INSStaggeredPressureBcCoef*>(P_bc_coef);
    P_ins_bc_coef->setPhysicalBcCoefs(u_bc_coefs);
    P_ins_bc_coef->setSolutionTime(solution_time);

    // The exact solution (u, p), the right-hand side (f, h), and the computed solution (u_sol, p_sol).
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("level_solver_open_boundary");
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
    Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");
    Pointer<SideVariable<NDIM, double>> f_var = new SideVariable<NDIM, double>("f");
    Pointer<CellVariable<NDIM, double>> h_var = new CellVariable<NDIM, double>("h");
    Pointer<SideVariable<NDIM, double>> u_sol_var = new SideVariable<NDIM, double>("u_sol");
    Pointer<CellVariable<NDIM, double>> p_sol_var = new CellVariable<NDIM, double>("p_sol");
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);
    auto register_data = [&](Pointer<Variable<NDIM>> var)
    {
        const int idx = var_db->registerVariableAndContext(var, ctx, IntVector<NDIM>(1));
        level->allocatePatchData(idx, solution_time);
        return idx;
    };
    const int u_idx = register_data(u_var), p_idx = register_data(p_var);
    const int f_idx = register_data(f_var), h_idx = register_data(h_var);
    const int u_sol_idx = register_data(u_sol_var), p_sol_idx = register_data(p_sol_var);
    SAMRAIVectorReal<NDIM, double> x("x", patch_hierarchy, 0, 0), b("b", patch_hierarchy, 0, 0),
        x_sol("x_sol", patch_hierarchy, 0, 0);
    x.addComponent(u_var, u_idx);
    x.addComponent(p_var, p_idx);
    b.addComponent(f_var, f_idx);
    b.addComponent(h_var, h_idx);
    x_sol.addComponent(u_sol_var, u_sol_idx);
    x_sol.addComponent(p_sol_var, p_sol_idx);
    x_sol.setToScalar(0.0);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
        u_data->fillAll(0.0);
        p_data->fillAll(0.0);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); it; it++)
            {
                (*u_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) = exact_value(it(), axis);
            }
        }
        for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
        {
            (*p_data)(it()) = exact_value(it(), NDIM);
        }
    }

    // Prescribe the normal velocity at the physical boundary according to the boundary data.
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
    bc_helper->cacheBcCoefData(u_bc_coefs, solution_time, patch_hierarchy);
    bc_helper->enforceNormalVelocityBoundaryConditions(u_idx, p_idx, U_bc_coefs, solution_time, false);

    // The right-hand side is the operator applied to the exact solution with inhomogeneous boundary conditions.
    PoissonSpecifications spec("spec");
    spec.setCConstant(1.0);
    spec.setDConstant(-mu);
    StaggeredStokesOperator op("StaggeredStokesOperator", /*homogeneous_bc*/ false);
    op.setVelocityPoissonSpecifications(spec);
    op.setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    op.setPhysicalBoundaryHelper(bc_helper);
    op.setSolutionTime(solution_time);
    op.setTimeInterval(solution_time, solution_time);
    op.initializeOperatorState(x, b);
    op.apply(x, b);

    // Solve from a zero initial guess.
    StaggeredStokesPETScLevelSolver solver("StaggeredStokesPETScLevelSolver", solver_db, "level_solver_open_boundary_");
    solver.setVelocityPoissonSpecifications(spec);
    solver.setComponentsHaveNullSpace(false, false);
    solver.setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    solver.setPhysicalBoundaryHelper(bc_helper);
    solver.setSolutionTime(solution_time);
    solver.setTimeInterval(solution_time, solution_time);
    solver.setHomogeneousBc(false);
    solver.setInitialGuessNonzero(false);
    solver.initializeSolverState(x_sol, b);
    const bool converged = solver.solveSystem(x_sol, b);
    solver.deallocateSolverState();

    // Compare. The velocity on a physical boundary at which the normal velocity is prescribed is compared too.
    double f_norm = 0.0, h_norm = 0.0, u_error = 0.0, p_error = 0.0;
    bool solution_is_finite = true;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_idx);
        Pointer<SideData<NDIM, double>> u_sol_data = patch->getPatchData(u_sol_idx);
        Pointer<CellData<NDIM, double>> h_data = patch->getPatchData(h_idx);
        Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
        Pointer<CellData<NDIM, double>> p_sol_data = patch->getPatchData(p_sol_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); it; it++)
            {
                const SideIndex<NDIM> s(it(), axis, SideIndex<NDIM>::Lower);
                const double u_difference = std::abs((*u_sol_data)(s) - (*u_data)(s));
                solution_is_finite = solution_is_finite && std::isfinite(u_difference);
                f_norm = std::max(f_norm, std::abs((*f_data)(s)));
                u_error = std::max(u_error, u_difference);
            }
        }
        for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
        {
            const double p_difference = (*p_sol_data)(it()) - (*p_data)(it());
            solution_is_finite = solution_is_finite && std::isfinite(p_difference);
            h_norm = std::max(h_norm, std::abs((*h_data)(it())));
            p_error = std::max(p_error, std::abs(p_difference));
        }
    }
    if (!converged || !solution_is_finite)
    {
        TBOX_ERROR("the level solver did not converge or its solution is not finite\n");
    }
    if (u_error > 1.0e-9 * f_norm || p_error > 1.0e-9 * f_norm)
    {
        TBOX_ERROR("the level solution differs from the exact solution: the largest velocity error is "
                   << u_error << " and the largest pressure error is " << p_error
                   << "; the largest magnitude of the right-hand side is " << f_norm << "\n");
    }
    plog << std::setprecision(10) << "Largest magnitude of the momentum right-hand side: " << f_norm
         << "\nLargest magnitude of the continuity right-hand side: " << h_norm << "\n";
} // check_level_solver
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");

        Pointer<INSHierarchyIntegrator> ins_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               ins_integrator,
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

        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            u_bc_coefs[d] =
                new muParserRobinBcCoefs("u_bc_coefs_" + std::to_string(d),
                                         app_initializer->getComponentDatabase("VelocityBcCoefs_" + std::to_string(d)),
                                         grid_geometry);
        }
        ins_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);

        ins_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        check_level_solver(
            patch_hierarchy, ins_integrator, u_bc_coefs, app_initializer->getComponentDatabase("level_solver_db"));

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }
    }
} // main
