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
#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/Database.h>
#include <tbox/Logger.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include <ibamr/app_namespaces.h>

// Solve a Stokes problem on a periodic level of several patches along a traction boundary, with a coarse-fine
// boundary opposite. Then change the coefficients, initialize the solver again, and solve again. A solver that is
// initialized once with the changed coefficients must give the same solution, bit for bit.
// Counts the initializations of the solver that compute all the data and the initializations that compute only the data
// that depend on the coefficients.
class CheckedSolver : public StaggeredStokesPETScLevelSolver
{
public:
    using StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver;
    int d_number_of_full_initializations = 0;
    int d_number_of_coefficient_initializations = 0;

    void check(const int full, const int coefficient) const
    {
        if (d_number_of_full_initializations != full || d_number_of_coefficient_initializations != coefficient)
        {
            TBOX_ERROR("The solver ran " << d_number_of_full_initializations << " full initializations and "
                                         << d_number_of_coefficient_initializations
                                         << " coefficient initializations instead of " << full << " and " << coefficient
                                         << ".\n");
        }
    }

protected:
    void initializeSolverStateSpecialized(const SAMRAIVectorReal<NDIM, double>& x,
                                          const SAMRAIVectorReal<NDIM, double>& b) override
    {
        ++d_number_of_full_initializations;
        StaggeredStokesPETScLevelSolver::initializeSolverStateSpecialized(x, b);
    }
    void initializeCoefficientsSpecialized(const SAMRAIVectorReal<NDIM, double>& x,
                                           const SAMRAIVectorReal<NDIM, double>& b) override
    {
        ++d_number_of_coefficient_initializations;
        StaggeredStokesPETScLevelSolver::initializeCoefficientsSpecialized(x, b);
    }
};

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app->getInputDatabase();

    Pointer<INSStaggeredHierarchyIntegrator> integrator = new INSStaggeredHierarchyIntegrator(
        "INSStaggeredHierarchyIntegrator", input_db->getDatabase("INSStaggeredHierarchyIntegrator"), false);
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry =
        new CartesianGridGeometry<NDIM>("CartesianGeometry", input_db->getDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", integrator, input_db->getDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> load_balancer =
        new LoadBalancer<NDIM>("LoadBalancer", input_db->getDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm = new GriddingAlgorithm<NDIM>(
        "GriddingAlgorithm", input_db->getDatabase("GriddingAlgorithm"), error_detector, box_generator, load_balancer);

    // The boundary condition objects of the integrator: Neumann conditions on the velocity at the lower boundary are
    // traction conditions.
    std::vector<std::unique_ptr<LocationIndexRobinBcCoefs<NDIM>>> physical_storage;
    std::vector<RobinBcCoefStrategy<NDIM>*> physical_bcs(NDIM, nullptr);
    for (int axis = 0; axis < NDIM; ++axis)
    {
        physical_storage.push_back(std::make_unique<LocationIndexRobinBcCoefs<NDIM>>("physical_bc", nullptr));
        for (int face = 0; face < 2 * NDIM; ++face)
        {
            physical_storage.back()->setBoundarySlope(face, 0.5);
        }
        physical_bcs[axis] = physical_storage.back().get();
    }
    std::vector<std::unique_ptr<INSStaggeredVelocityBcCoef>> velocity_storage;
    std::vector<RobinBcCoefStrategy<NDIM>*> velocity_bcs(NDIM, nullptr);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        velocity_storage.push_back(
            std::make_unique<INSStaggeredVelocityBcCoef>(axis, integrator, physical_bcs, TRACTION));
        velocity_bcs[axis] = velocity_storage.back().get();
    }
    INSStaggeredPressureBcCoef pressure_bc(integrator, physical_bcs, TRACTION);

    // The integrator initializes the hierarchy after the boundary condition objects that refer to it exist.
    integrator->initializePatchHierarchy(hierarchy, gridding_algorithm);
    const int ln = hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("stokes_petsc_level_solver_reinitialize");
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
    Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");
    Pointer<SideVariable<NDIM, double>> f_var = new SideVariable<NDIM, double>("f");
    Pointer<CellVariable<NDIM, double>> h_var = new CellVariable<NDIM, double>("h");
    const int u_idx = var_db->registerVariableAndContext(u_var, context, IntVector<NDIM>(1));
    const int p_idx = var_db->registerVariableAndContext(p_var, context, IntVector<NDIM>(1));
    const int f_idx = var_db->registerVariableAndContext(f_var, context, IntVector<NDIM>(1));
    const int h_idx = var_db->registerVariableAndContext(h_var, context, IntVector<NDIM>(1));
    for (const int idx : { u_idx, p_idx, f_idx, h_idx })
    {
        level->allocatePatchData(idx);
    }

    Pointer<StaggeredStokesPhysicalBoundaryHelper> helper = new StaggeredStokesPhysicalBoundaryHelper();
    helper->cacheBcCoefData(physical_bcs, 0.0, hierarchy);

    SAMRAIVectorReal<NDIM, double> x("x", hierarchy, ln, ln), b("b", hierarchy, ln, ln);
    x.addComponent(u_var, u_idx);
    x.addComponent(p_var, p_idx);
    b.addComponent(f_var, f_idx);
    b.addComponent(h_var, h_idx);

    // The coefficients of the first and of the second problem.
    const double c[2] = { 1.0, 2.5 };
    const double d[2] = { -1.0, -1.5 };
    PoissonSpecifications coefficients("coefficients");
    auto configure = [&](StaggeredStokesPETScLevelSolver& solver, const int problem)
    {
        coefficients.setCConstant(c[problem]);
        coefficients.setDConstant(d[problem]);
        solver.setVelocityPoissonSpecifications(coefficients);
        solver.setComponentsHaveNullSpace(false, false);
        solver.setPhysicalBcCoefs(velocity_bcs, &pressure_bc);
        solver.setPhysicalBoundaryHelper(helper);
        solver.setHomogeneousBc(false);
    };

    // Solve with the initialized solver from the same vectors, and return the solution values in the interior.
    auto solve = [&](StaggeredStokesPETScLevelSolver& solver, int& iterations)
    {
        int counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
            Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_idx);
            Pointer<CellData<NDIM, double>> h_data = patch->getPatchData(h_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(u_data->getArrayData(axis).getBox()); it; it++)
                {
                    (*u_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) = 1.0 + 0.125 * (counter++ % 13);
                }
                for (Box<NDIM>::Iterator it(f_data->getArrayData(axis).getBox()); it; it++)
                {
                    (*f_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) = 0.5 + 0.25 * (counter++ % 7);
                }
            }
            for (Box<NDIM>::Iterator it(p_data->getGhostBox()); it; it++)
            {
                (*p_data)(CellIndex<NDIM>(it())) = 0.25 * (counter++ % 5);
                (*h_data)(CellIndex<NDIM>(it())) = 0.0;
            }
        }
        solver.solveSystem(x, b);
        iterations = solver.getNumIterations();
        std::vector<double> solution;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); it; it++)
                {
                    solution.push_back((*u_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)));
                }
            }
            for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
            {
                solution.push_back((*p_data)(CellIndex<NDIM>(it())));
            }
        }
        return solution;
    };

    // The solver that is initialized twice.
    CheckedSolver reinitialized_solver(
        "reinitialized_solver", input_db->getDatabase("solver_db"), "reinitialized_solver_");
    int iterations_first, iterations_reinitialized, iterations_fresh;
    configure(reinitialized_solver, 0);
    reinitialized_solver.initializeSolverState(x, b);
    reinitialized_solver.check(1, 1);
    const std::vector<double> first = solve(reinitialized_solver, iterations_first);
    configure(reinitialized_solver, 1);
    reinitialized_solver.initializeSolverState(x, b);
    reinitialized_solver.check(1, 2);
    const std::vector<double> reinitialized = solve(reinitialized_solver, iterations_reinitialized);
    reinitialized_solver.deallocateSolverState();

    // The solver that is initialized once.
    StaggeredStokesPETScLevelSolver fresh_solver("fresh_solver", input_db->getDatabase("solver_db"), "fresh_solver_");
    configure(fresh_solver, 1);
    fresh_solver.initializeSolverState(x, b);
    const std::vector<double> fresh = solve(fresh_solver, iterations_fresh);
    fresh_solver.deallocateSolverState();

    double max_difference = 0.0, max_change = 0.0, solution_min = std::numeric_limits<double>::max(),
           solution_max = std::numeric_limits<double>::lowest();
    for (size_t k = 0; k < fresh.size(); ++k)
    {
        if (!std::isfinite(fresh[k]) || !std::isfinite(reinitialized[k]))
        {
            TBOX_ERROR("A solution value is not finite.\n");
        }
        max_difference = std::max(max_difference, std::abs(reinitialized[k] - fresh[k]));
        max_change = std::max(max_change, std::abs(fresh[k] - first[k]));
        solution_min = std::min(solution_min, fresh[k]);
        solution_max = std::max(solution_max, fresh[k]);
    }
    max_difference = IBTK_MPI::maxReduction(max_difference);
    max_change = IBTK_MPI::maxReduction(max_change);
    solution_min = IBTK_MPI::minReduction(solution_min);
    solution_max = IBTK_MPI::maxReduction(solution_max);
    const int number_of_values = IBTK_MPI::sumReduction(static_cast<int>(fresh.size()));
    if (reinitialized != fresh || iterations_reinitialized != iterations_fresh)
    {
        TBOX_ERROR(
            "The solver that was initialized twice differs from a solver that was initialized once: the "
            "solutions differ by "
            << max_difference << ", and the iteration counts are " << iterations_reinitialized << " and "
            << iterations_fresh << ".\n");
    }
    if (!(max_change > 1.0e-3))
    {
        TBOX_ERROR("The solution did not change with the coefficients: " << max_change << ".\n");
    }

    plog << "number_of_patches = " << level->getNumberOfPatches() << '\n';
    plog << "number_of_values = " << number_of_values << '\n';
    plog << "iterations_first = " << iterations_first << '\n';
    plog << "iterations_reinitialized = " << iterations_reinitialized << '\n';
    plog << "iterations_fresh = " << iterations_fresh << '\n';
    plog << "max_difference_reinitialized_fresh = " << max_difference << '\n';
    plog << "max_change_from_first_problem = " << max_change << '\n';
    plog << "solution_min = " << solution_min << '\n';
    plog << "solution_max = " << solution_max << '\n';
    return 0;
}
