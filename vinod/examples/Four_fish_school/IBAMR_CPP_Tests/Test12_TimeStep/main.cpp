// =============================================================================
// TEST 12: Time-step Convergence Study
// =============================================================================
// Validates temporal convergence by running with multiple time steps
// Tests: dt → dt/2 → dt/4 convergence, CFL constraints, temporal accuracy
// Expected: 2nd-order convergence (Crank-Nicolson)
// =============================================================================

#include <SAMRAI_config.h>
#include <petscsys.h>

// SAMRAI headers
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>
#include <CellData.h>
#include <CellVariable.h>

// IBAMR headers
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

// IBTK headers
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Common test utilities
#include "TestUtilities.h"
#include "ErrorCalculator.h"
#include "AnalyticalSolutions.h"

#include <ibamr/app_namespaces.h>
#include <vector>

using namespace TestUtilities;

// Structure to store convergence test results
struct ConvergenceResult
{
    double dt;
    int num_steps;
    double L2_error;
    double Linf_error;
    double CFL_advection;
    double CFL_diffusion;
};

int main(int argc, char* argv[])
{
    // =========================================================================
    // Initialize IBAMR/SAMRAI/PETSc
    // =========================================================================
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    TestUtils::printTestHeader("Time-step Convergence Study", 12);

    Timer total_timer;
    total_timer.start();

    bool test_passed = true;
    ResultLogger logger("test12_results.txt");

    // Time steps to test (dt, dt/2, dt/4)
    std::vector<double> dt_values = {0.02, 0.01, 0.005};
    std::vector<ConvergenceResult> results;

    logger.log("Testing time-step convergence with dt = 0.02, 0.01, 0.005");

    {
        // =====================================================================
        // Get common parameters (read once)
        // =====================================================================
        TestUtils::printProgress("Reading common parameters...");
        Pointer<AppInitializer> temp_app = new AppInitializer(argc, argv, "test12_temp.log");
        Pointer<Database> temp_db = temp_app->getInputDatabase();

        const double kappa = temp_db->getDouble("diffusion_coefficient");
        const double u = temp_db->getDoubleWithDefault("advection_velocity_u", 0.3);
        const double v = temp_db->getDoubleWithDefault("advection_velocity_v", 0.2);
        const double end_time = temp_app->getComponentDatabase("Main")->getDouble("END_TIME");

        // Grid parameters for CFL computation
        Pointer<Database> geom_db = temp_app->getComponentDatabase("CartesianGeometry");
        tbox::Array<int> domain_boxes = geom_db->getIntegerArray("domain_boxes");
        tbox::Array<double> x_lo_arr = geom_db->getDoubleArray("x_lo");
        tbox::Array<double> x_up_arr = geom_db->getDoubleArray("x_up");

        int Nx = domain_boxes[NDIM] - domain_boxes[0] + 1;
        double Lx = x_up_arr[0] - x_lo_arr[0];
        double dx = Lx / Nx;

        logger.logParameter("Diffusion coefficient", kappa);
        logger.logParameter("Advection velocity u", u);
        logger.logParameter("Advection velocity v", v);
        logger.logParameter("Grid spacing dx", dx);
        logger.logParameter("End time", end_time);

        std::cout << "\n  Convergence Study Configuration:\n";
        std::cout << "    Grid spacing: dx = " << dx << "\n";
        std::cout << "    Diffusion coefficient: κ = " << kappa << "\n";
        std::cout << "    Advection velocity: (u,v) = (" << u << ", " << v << ")\n";
        std::cout << "    End time: T = " << end_time << "\n\n";

        // =====================================================================
        // Loop over time steps
        // =====================================================================
        for (size_t dt_idx = 0; dt_idx < dt_values.size(); ++dt_idx)
        {
            const double dt = dt_values[dt_idx];
            const int num_steps = static_cast<int>(std::ceil(end_time / dt));

            TestUtils::printSectionSeparator("Time-step Test: dt = " + std::to_string(dt));
            std::cout << "  Number of steps: " << num_steps << "\n";

            // Compute CFL numbers
            double U_mag = sqrt(u*u + v*v);
            double CFL_adv = TestUtils::computeCFL_advection(U_mag, dt, dx);
            double CFL_diff = TestUtils::computeCFL_diffusion(kappa, dt, dx);

            std::cout << "  CFL (advection): " << CFL_adv << "\n";
            std::cout << "  CFL (diffusion): " << CFL_diff << "\n\n";

            ConvergenceResult conv_result;
            conv_result.dt = dt;
            conv_result.num_steps = num_steps;
            conv_result.CFL_advection = CFL_adv;
            conv_result.CFL_diffusion = CFL_diff;

            // =================================================================
            // Initialize application for this dt
            // =================================================================
            Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test12.log");
            Pointer<Database> input_db = app_initializer->getInputDatabase();

            // Override dt in Main database
            Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");
            main_db->putDouble("DT", dt);

            // =================================================================
            // Create integrators
            // =================================================================
            Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

            Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = new AdvDiffHierarchyIntegrator(
                "AdvDiffHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffHierarchyIntegrator"));

            adv_diff_integrator->setAdvectionVelocity(navier_stokes_integrator->getAdvectionVelocityVariable());
            navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);

            // =================================================================
            // Create grid geometry
            // =================================================================
            Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
                "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
            Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

            // =================================================================
            // Create gridding algorithm
            // =================================================================
            Pointer<StandardTagAndInitialize<NDIM>> error_detector =
                new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                                   navier_stokes_integrator,
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

            // =================================================================
            // Set initial conditions
            // =================================================================
            if (input_db->keyExists("OdorInitialConditions"))
            {
                Pointer<CartGridFunction> C_init = new muParserCartGridFunction(
                    "C_init", app_initializer->getComponentDatabase("OdorInitialConditions"), grid_geometry);
                adv_diff_integrator->setInitialConditions(C_init);
            }

            // =================================================================
            // Set boundary conditions
            // =================================================================
            const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
            RobinBcCoefStrategy<NDIM>* C_bc_coef = nullptr;
            if (!(periodic_shift.min() > 0) && input_db->keyExists("OdorBcCoefs"))
            {
                C_bc_coef = new muParserRobinBcCoefs(
                    "C_bc_coef", app_initializer->getComponentDatabase("OdorBcCoefs"), grid_geometry);
                adv_diff_integrator->setPhysicalBcCoef(C_bc_coef);
            }

            // =================================================================
            // Set diffusion coefficient
            // =================================================================
            adv_diff_integrator->setDiffusionCoefficient(kappa);

            // =================================================================
            // Initialize hierarchy
            // =================================================================
            navier_stokes_integrator->initializeHierarchyIntegrator(gridding_algorithm);
            gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
            int tag_buffer = 1;
            int level_number = 0;
            bool done = false;
            while (!done && (gridding_algorithm->levelCanBeRefined(level_number)))
            {
                gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
                done = !patch_hierarchy->finerLevelExists(level_number);
                ++level_number;
            }

            const int finest_level = patch_hierarchy->getFinestLevelNumber();
            for (int ln = 0; ln <= finest_level; ++ln)
            {
                patch_hierarchy->getPatchLevel(ln)->allocatePatchData(navier_stokes_integrator->getPatchDataIndex(), 0.0);
            }

            navier_stokes_integrator->initializeIntegratorData(0.0);

            const int C_idx = adv_diff_integrator->getTransportedQuantity()->getPatchDataIndex();

            // =================================================================
            // Time integration
            // =================================================================
            std::cout << "  Running time integration...\n";

            double loop_time = navier_stokes_integrator->getIntegratorTime();
            double loop_time_end = navier_stokes_integrator->getEndTime();
            int iteration_num = navier_stokes_integrator->getIntegratorStep();

            while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) &&
                   navier_stokes_integrator->stepsRemaining())
            {
                iteration_num = navier_stokes_integrator->getIntegratorStep();
                loop_time = navier_stokes_integrator->getIntegratorTime();

                if (iteration_num % 20 == 0)
                {
                    std::cout << "    Step " << iteration_num << "/" << num_steps
                              << ", t = " << loop_time << std::endl;
                }

                navier_stokes_integrator->advanceHierarchy(dt);
                loop_time += dt;
            }

            // =================================================================
            // Compute analytical solution
            // =================================================================
            std::cout << "  Computing exact solution...\n";

            Pointer<CellVariable<NDIM, double>> C_exact_var = new CellVariable<NDIM, double>("C_exact");
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            const int C_exact_idx = var_db->registerVariableAndContext(
                C_exact_var, var_db->getContext("EXACT_" + std::to_string(dt_idx)), IntVector<NDIM>(1));

            for (int ln = 0; ln <= finest_level; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
                level->allocatePatchData(C_exact_idx, loop_time);
            }

            const double* x_lo = grid_geometry->getXLower();

            for (int ln = 0; ln <= finest_level; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    const Box<NDIM>& patch_box = patch->getBox();
                    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();

                    Pointer<CellData<NDIM, double>> C_exact_data = patch->getPatchData(C_exact_idx);
                    const double* dx_arr = pgeom->getDx();

                    for (CellIterator<NDIM> ci(patch_box); ci; ci++)
                    {
                        const CellIndex<NDIM>& idx = ci();

                        double X[NDIM];
                        for (int d = 0; d < NDIM; ++d)
                        {
                            X[d] = x_lo[d] + dx_arr[d] * (static_cast<double>(idx(d) - patch_box.lower()(d)) + 0.5);
                        }

                        // Use Gaussian diffusion solution
                        (*C_exact_data)(idx) = AnalyticalSolutions::gaussianDiffusion2D(
                            X[0], X[1], loop_time, kappa, 0.0, 0.0, 1.0);
                    }
                }
            }

            // =================================================================
            // Compute errors
            // =================================================================
            std::cout << "  Computing errors...\n";

            double L2_error = ErrorCalculator::computeL2Error(patch_hierarchy, C_idx, C_exact_idx, loop_time);
            double Linf_error = ErrorCalculator::computeLinfError(patch_hierarchy, C_idx, C_exact_idx);

            conv_result.L2_error = L2_error;
            conv_result.Linf_error = Linf_error;

            std::cout << "  L2 error:   " << TestUtils::formatScientific(L2_error) << "\n";
            std::cout << "  Linf error: " << TestUtils::formatScientific(Linf_error) << "\n\n";

            results.push_back(conv_result);

            // Cleanup
            for (int ln = 0; ln <= finest_level; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
                level->deallocatePatchData(C_exact_idx);
            }

            if (C_bc_coef) delete C_bc_coef;

        } // End dt loop

        // =====================================================================
        // Analyze convergence rates
        // =====================================================================
        TestUtils::printSectionSeparator("Temporal Convergence Analysis");

        std::cout << "\n";
        std::cout << "  dt        Steps   CFL_adv  CFL_diff   L2_error      Linf_error\n";
        std::cout << "  --------  ------  -------  --------  ------------  ------------\n";

        for (const auto& res : results)
        {
            std::cout << "  " << std::setw(8) << std::fixed << std::setprecision(4) << res.dt
                      << "  " << std::setw(6) << res.num_steps
                      << "  " << std::setw(7) << std::setprecision(3) << res.CFL_advection
                      << "  " << std::setw(8) << std::setprecision(3) << res.CFL_diffusion
                      << "  " << std::setw(12) << std::scientific << std::setprecision(4) << res.L2_error
                      << "  " << std::setw(12) << std::scientific << res.Linf_error
                      << "\n";

            logger.logParameter("dt_" + std::to_string(res.dt) + "_L2_error", res.L2_error);
            logger.logParameter("dt_" + std::to_string(res.dt) + "_Linf_error", res.Linf_error);
        }
        std::cout << "\n";

        // Compute convergence rates (L2 norm)
        if (results.size() >= 2)
        {
            std::cout << "  Temporal Convergence Rates (L2):\n";
            for (size_t i = 1; i < results.size(); ++i)
            {
                double rate = ErrorCalculator::computeConvergenceRate(
                    results[i-1].L2_error, results[i].L2_error,
                    results[i-1].dt, results[i].dt);

                std::cout << "    dt " << results[i-1].dt << " → " << results[i].dt
                          << ": rate = " << std::setprecision(3) << rate
                          << " (expected ~2.0 for Crank-Nicolson)\n";

                logger.logParameter("convergence_rate_" + std::to_string(i), rate);
            }
            std::cout << "\n";
        }

        // =====================================================================
        // Test verdict
        // =====================================================================
        TestUtils::printSectionSeparator("Test Verdict");

        std::vector<std::string> checks;
        std::vector<bool> results_pass;

        // Check 1: Errors decrease with finer dt
        checks.push_back("Errors decrease as dt → 0");
        results_pass.push_back(results[2].L2_error < results[0].L2_error);

        // Check 2: Convergence rate is approximately 2nd order
        if (results.size() >= 2)
        {
            double rate = ErrorCalculator::computeConvergenceRate(
                results[0].L2_error, results[2].L2_error,
                results[0].dt, results[2].dt);
            checks.push_back("Convergence rate ~2.0 (1.5 < rate < 2.5)");
            results_pass.push_back(rate > 1.5 && rate < 2.5);
        }

        // Check 3: Finest dt has acceptable error
        checks.push_back("Finest dt error < 0.01");
        results_pass.push_back(results.back().L2_error < 0.01);

        test_passed = true;
        for (size_t i = 0; i < checks.size(); ++i)
        {
            std::string status = results_pass[i] ? "[PASS]" : "[FAIL]";
            std::cout << "  " << status << " " << checks[i] << "\n";
            test_passed = test_passed && results_pass[i];
        }

    } // End scope

    total_timer.stop();

    // =========================================================================
    // Final report
    // =========================================================================
    std::cout << "\n";
    std::cout << "Elapsed time: " << total_timer.getFormattedTime() << "\n";

    logger.logTestResult("Test 12: Time-step Convergence", test_passed);
    logger.log("Total elapsed time: " + total_timer.getFormattedTime());

    TestUtils::printTestFooter("Test 12: Time-step Convergence", test_passed,
                              test_passed ? "Temporal convergence validated" : "Convergence rate not achieved");

    return test_passed ? 0 : 1;
}
