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
#include <ibamr/StaggeredStokesBoxRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesFACPreconditioner.h>
#include <ibamr/StaggeredStokesLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/Database.h>
#include <tbox/Logger.h>

#include <BergerRigoutsos.h>
#include <BoxList.h>
#include <CartesianGridGeometry.h>
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
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

// Apply a Stokes FAC preconditioner on a two-level hierarchy that is periodic in the first direction, then change the
// coefficients and the boundary conditions, initialize the preconditioner again, and apply it again, in three problems:
// slip walls, traction conditions (so that the box relaxation smoother adds the normal stress term at the physical
// boundary), and slip walls. Each result must equal, bit for bit, that of a preconditioner of the same class
// initialized once for that problem.
namespace
{
// Counts the initializations of the operator that compute all the data and the initializations that rebuild only the
// data that depend on the coefficients and the boundary conditions.
template <class Operator>
class CheckedOperator : public Operator
{
public:
    using Operator::Operator;
    int d_number_of_full_initializations = 0;
    int d_number_of_coefficient_initializations = 0;

    void check(const int full, const int coefficient) const
    {
        if (d_number_of_full_initializations != full || d_number_of_coefficient_initializations != coefficient)
        {
            TBOX_ERROR("The operator ran " << d_number_of_full_initializations << " full initializations and "
                                           << d_number_of_coefficient_initializations
                                           << " coefficient initializations instead of " << full << " and "
                                           << coefficient << ".\n");
        }
    }

protected:
    void initializeOperatorStateSpecialized(const SAMRAIVectorReal<NDIM, double>& solution,
                                            const SAMRAIVectorReal<NDIM, double>& rhs,
                                            const int coarsest_reset_ln,
                                            const int finest_reset_ln) override
    {
        ++d_number_of_full_initializations;
        Operator::initializeOperatorStateSpecialized(solution, rhs, coarsest_reset_ln, finest_reset_ln);
    }
    void initializeOperatorCoefficientsSpecialized(const SAMRAIVectorReal<NDIM, double>& solution,
                                                   const SAMRAIVectorReal<NDIM, double>& rhs) override
    {
        ++d_number_of_coefficient_initializations;
        Operator::initializeOperatorCoefficientsSpecialized(solution, rhs);
    }
};

template <class Operator>
void
run_case(Pointer<PatchHierarchy<NDIM>> hierarchy,
         Pointer<Database> input_db,
         const std::vector<RobinBcCoefStrategy<NDIM>*> (&velocity_bcs)[3],
         RobinBcCoefStrategy<NDIM>* const (&pressure_bcs)[3],
         const Pointer<StaggeredStokesPhysicalBoundaryHelper> (&helpers)[3])
{
    const int finest_ln = hierarchy->getFinestLevelNumber();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("stokes_fac_reinitialize");
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
    Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");
    Pointer<SideVariable<NDIM, double>> f_var = new SideVariable<NDIM, double>("f");
    Pointer<CellVariable<NDIM, double>> h_var = new CellVariable<NDIM, double>("h");
    const int u_idx = var_db->registerVariableAndContext(u_var, context, IntVector<NDIM>(1));
    const int p_idx = var_db->registerVariableAndContext(p_var, context, IntVector<NDIM>(1));
    const int f_idx = var_db->registerVariableAndContext(f_var, context, IntVector<NDIM>(1));
    const int h_idx = var_db->registerVariableAndContext(h_var, context, IntVector<NDIM>(1));
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        for (const int idx : { u_idx, p_idx, f_idx, h_idx })
        {
            hierarchy->getPatchLevel(ln)->allocatePatchData(idx);
        }
    }
    SAMRAIVectorReal<NDIM, double> x("x", hierarchy, 0, finest_ln), b("b", hierarchy, 0, finest_ln);
    x.addComponent(u_var, u_idx);
    x.addComponent(p_var, p_idx);
    b.addComponent(f_var, f_idx);
    b.addComponent(h_var, h_idx);

    // The coefficients of the three problems. Slip walls leave the pressure determined up to a constant.
    const double c[3] = { 1.0, 2.5, 1.7 };
    const double d[3] = { -1.0, -1.5, -1.2 };
    const bool has_pressure_nullspace[3] = { true, false, true };
    PoissonSpecifications coefficients("coefficients");
    auto configure = [&](StaggeredStokesSolver& solver, const int problem)
    {
        coefficients.setCConstant(c[problem]);
        coefficients.setDConstant(d[problem]);
        solver.setVelocityPoissonSpecifications(coefficients);
        solver.setComponentsHaveNullSpace(false, has_pressure_nullspace[problem]);
        solver.setPhysicalBcCoefs(velocity_bcs[problem], pressure_bcs[problem]);
        solver.setPhysicalBoundaryHelper(helpers[problem]);
        solver.setSolutionTime(1.0);
        solver.setTimeInterval(0.0, 1.0);
    };

    // Apply the preconditioner to the same right-hand side, and return the values in the interiors of the patches.
    auto apply = [&](StaggeredStokesFACPreconditioner& preconditioner, const int problem)
    {
        int counter = 0;
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
                Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
                Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_idx);
                Pointer<CellData<NDIM, double>> h_data = patch->getPatchData(h_idx);
                u_data->fillAll(0.0);
                p_data->fillAll(0.0);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator it(f_data->getArrayData(axis).getBox()); it; it++)
                    {
                        (*f_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) = 0.5 + 0.25 * (counter++ % 7);
                    }
                }
                for (Box<NDIM>::Iterator it(h_data->getGhostBox()); it; it++)
                {
                    (*h_data)(CellIndex<NDIM>(it())) = 0.125 * (counter++ % 5);
                }
            }
        }
        // The problems with a pressure null space have a right-hand side h of zero mean over the composite grid.
        if (has_pressure_nullspace[problem])
        {
            double sum = 0.0, volume = 0.0;
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
                BoxList<NDIM> covered;
                if (ln < finest_ln)
                {
                    Pointer<PatchLevel<NDIM>> finer_level = hierarchy->getPatchLevel(ln + 1);
                    covered = BoxList<NDIM>(finer_level->getBoxes());
                    covered.coarsen(finer_level->getRatioToCoarserLevel());
                }
                double cell_volume = 1.0;
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    cell_volume /= level->getRatio()(axis);
                }
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    Pointer<CellData<NDIM, double>> h_data = patch->getPatchData(h_idx);
                    BoxList<NDIM> uncovered(patch->getBox());
                    uncovered.removeIntersections(covered);
                    for (BoxList<NDIM>::Iterator boxes(uncovered); boxes; boxes++)
                    {
                        for (Box<NDIM>::Iterator it(boxes()); it; it++)
                        {
                            sum += cell_volume * (*h_data)(CellIndex<NDIM>(it()));
                            volume += cell_volume;
                        }
                    }
                }
            }
            const double mean = IBTK_MPI::sumReduction(sum) / IBTK_MPI::sumReduction(volume);
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<CellData<NDIM, double>> h_data = level->getPatch(p())->getPatchData(h_idx);
                    for (Box<NDIM>::Iterator it(h_data->getGhostBox()); it; it++)
                    {
                        (*h_data)(CellIndex<NDIM>(it())) -= mean;
                    }
                }
            }
        }
        preconditioner.solveSystem(x, b);
        std::vector<double> solution;
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
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
        }
        return solution;
    };

    // The preconditioner that is initialized again for each problem.
    Pointer<Database> operator_db = input_db->getDatabase("operator_db");
    Pointer<Database> preconditioner_db = input_db->getDatabase("preconditioner_db");
    Pointer<CheckedOperator<Operator>> reinitialized_operator =
        new CheckedOperator<Operator>("reinitialized_operator", operator_db, "reinitialized_operator_");
    StaggeredStokesFACPreconditioner reinitialized_preconditioner(
        "reinitialized_preconditioner", reinitialized_operator, preconditioner_db, "reinitialized_preconditioner_");
    std::vector<double> reinitialized[3];
    for (int problem = 0; problem < 3; ++problem)
    {
        configure(reinitialized_preconditioner, problem);
        reinitialized_preconditioner.initializeSolverState(x, b);
        reinitialized_operator->check(1, problem);
        reinitialized[problem] = apply(reinitialized_preconditioner, problem);
    }
    reinitialized_preconditioner.deallocateSolverState();

    // The preconditioners that are initialized once.
    std::vector<double> fresh[3];
    for (int problem = 0; problem < 3; ++problem)
    {
        Pointer<Operator> fresh_operator = new Operator("fresh_operator", operator_db, "fresh_operator_");
        StaggeredStokesFACPreconditioner fresh_preconditioner(
            "fresh_preconditioner", fresh_operator, preconditioner_db, "fresh_preconditioner_");
        configure(fresh_preconditioner, problem);
        fresh_preconditioner.initializeSolverState(x, b);
        fresh[problem] = apply(fresh_preconditioner, problem);
        fresh_preconditioner.deallocateSolverState();
    }

    for (int problem = 0; problem < 3; ++problem)
    {
        double max_difference = 0.0, max_change = 0.0, solution_min = std::numeric_limits<double>::max(),
               solution_max = std::numeric_limits<double>::lowest();
        for (std::size_t k = 0; k < fresh[problem].size(); ++k)
        {
            if (!std::isfinite(fresh[problem][k]) || !std::isfinite(reinitialized[problem][k]))
            {
                TBOX_ERROR("A solution value is not finite.\n");
            }
            max_difference = std::max(max_difference, std::abs(reinitialized[problem][k] - fresh[problem][k]));
            if (problem > 0)
            {
                max_change = std::max(max_change, std::abs(fresh[problem][k] - fresh[problem - 1][k]));
            }
            solution_min = std::min(solution_min, fresh[problem][k]);
            solution_max = std::max(solution_max, fresh[problem][k]);
        }
        max_difference = IBTK_MPI::maxReduction(max_difference);
        max_change = IBTK_MPI::maxReduction(max_change);
        solution_min = IBTK_MPI::minReduction(solution_min);
        solution_max = IBTK_MPI::maxReduction(solution_max);
        const int number_of_values = IBTK_MPI::sumReduction(static_cast<int>(fresh[problem].size()));
        if (reinitialized[problem] != fresh[problem])
        {
            TBOX_ERROR("The preconditioner that was initialized again for problem "
                       << problem << " differs from a preconditioner that was initialized once: the results differ by "
                       << max_difference << ".\n");
        }
        if (problem > 0 && !(max_change > 1.0e-3))
        {
            TBOX_ERROR("The result did not change from problem " << problem - 1 << " to problem " << problem << ": "
                                                                 << max_change << ".\n");
        }

        plog << "problem = " << problem << '\n';
        plog << "number_of_values = " << number_of_values << '\n';
        plog << "max_difference_reinitialized_fresh = " << max_difference << '\n';
        plog << "max_change_from_previous_problem = " << max_change << '\n';
        plog << "solution_min = " << solution_min << '\n';
        plog << "solution_max = " << solution_max << '\n';
    }
    plog << "full_initializations = " << reinitialized_operator->d_number_of_full_initializations << '\n';
    plog << "coefficient_initializations = " << reinitialized_operator->d_number_of_coefficient_initializations << '\n';
}
} // namespace

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

    // Slip walls: Dirichlet conditions on the normal velocity, and Neumann conditions on the tangential velocity.
    std::vector<std::unique_ptr<LocationIndexRobinBcCoefs<NDIM>>> slip_storage;
    std::vector<RobinBcCoefStrategy<NDIM>*> slip_bcs(NDIM, nullptr);
    for (int axis = 0; axis < NDIM; ++axis)
    {
        slip_storage.push_back(std::make_unique<LocationIndexRobinBcCoefs<NDIM>>("slip_bc", nullptr));
        for (int face = 0; face < 2 * NDIM; ++face)
        {
            if (face / 2 == axis)
            {
                slip_storage.back()->setBoundaryValue(face, 0.0);
            }
            else
            {
                slip_storage.back()->setBoundarySlope(face, 0.0);
            }
        }
        slip_bcs[axis] = slip_storage.back().get();
    }

    // The boundary condition objects of the integrator: Neumann conditions on the velocity are traction conditions.
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
    std::vector<std::unique_ptr<INSStaggeredVelocityBcCoef>> traction_storage;
    std::vector<RobinBcCoefStrategy<NDIM>*> traction_bcs(NDIM, nullptr);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        traction_storage.push_back(
            std::make_unique<INSStaggeredVelocityBcCoef>(axis, integrator, physical_bcs, TRACTION));
        traction_bcs[axis] = traction_storage.back().get();
    }
    INSStaggeredPressureBcCoef pressure_bc(integrator, physical_bcs, TRACTION);

    // The integrator initializes the hierarchy after the boundary condition objects that refer to it exist.
    integrator->initializePatchHierarchy(hierarchy, gridding_algorithm);

    Pointer<StaggeredStokesPhysicalBoundaryHelper> slip_helper = new StaggeredStokesPhysicalBoundaryHelper();
    slip_helper->cacheBcCoefData(slip_bcs, 0.0, hierarchy);
    Pointer<StaggeredStokesPhysicalBoundaryHelper> traction_helper = new StaggeredStokesPhysicalBoundaryHelper();
    traction_helper->cacheBcCoefData(physical_bcs, 0.0, hierarchy);

    // The boundary condition objects and the helper of each of the three problems.
    const std::vector<RobinBcCoefStrategy<NDIM>*> velocity_bcs[3] = { slip_bcs, traction_bcs, slip_bcs };
    RobinBcCoefStrategy<NDIM>* const pressure_bcs[3] = { nullptr, &pressure_bc, nullptr };
    const Pointer<StaggeredStokesPhysicalBoundaryHelper> helpers[3] = { slip_helper, traction_helper, slip_helper };

    plog << "number_of_levels = " << hierarchy->getNumberOfLevels() << '\n';
    Pointer<PatchLevel<NDIM>> finest_level = hierarchy->getPatchLevel(hierarchy->getFinestLevelNumber());
    plog << "number_of_patches_on_finest_level = " << finest_level->getNumberOfPatches() << '\n';

    const std::string operator_class = input_db->getDatabase("test")->getString("operator_class");
    if (operator_class == "StaggeredStokesBoxRelaxationFACOperator")
    {
        run_case<StaggeredStokesBoxRelaxationFACOperator>(hierarchy, input_db, velocity_bcs, pressure_bcs, helpers);
    }
    else if (operator_class == "StaggeredStokesLevelRelaxationFACOperator")
    {
        run_case<StaggeredStokesLevelRelaxationFACOperator>(hierarchy, input_db, velocity_bcs, pressure_bcs, helpers);
    }
    else
    {
        TBOX_ERROR("Unknown operator_class \"" << operator_class << "\"\n"
                                               << "  valid values are the Stokes box and level relaxation FAC "
                                                  "operators\n");
    }
    return 0;
}
