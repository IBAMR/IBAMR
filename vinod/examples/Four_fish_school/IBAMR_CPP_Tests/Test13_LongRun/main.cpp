// =============================================================================
// TEST 13: Long Run Stability Test
// =============================================================================
// Validates solver stability over extended time integration
// Tests: Mass conservation, no drift, no NaN/Inf, solution boundedness
// Critical for production runs with long simulation times
// =============================================================================

#include <SAMRAI_config.h>
#include <petscsys.h>

// SAMRAI headers
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

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
#include <fstream>

using namespace TestUtilities;

int main(int argc, char* argv[])
{
    // =========================================================================
    // Initialize IBAMR/SAMRAI/PETSc
    // =========================================================================
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    TestUtils::printTestHeader("Long Run Stability Test", 13);

    Timer total_timer;
    total_timer.start();

    bool test_passed = true;
    ResultLogger logger("test13_results.txt");

    {
        // =====================================================================
        // Initialize application
        // =====================================================================
        TestUtils::printProgress("Initializing application...");
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test13.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get solver parameters
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");
        const double end_time = main_db->getDouble("END_TIME");
        const double dt = main_db->getDouble("DT");
        const int num_steps = static_cast<int>(std::ceil(end_time / dt));

        const double kappa = input_db->getDouble("diffusion_coefficient");

        logger.logParameter("End time", end_time);
        logger.logParameter("Time step", dt);
        logger.logParameter("Number of steps", num_steps);
        logger.logParameter("Diffusion coefficient", kappa);

        std::cout << "\n  Long Run Test Configuration:\n";
        std::cout << "    Duration: T = " << end_time << " (extended integration)\n";
        std::cout << "    Time step: dt = " << dt << "\n";
        std::cout << "    Total steps: " << num_steps << "\n";
        std::cout << "    Diffusion coefficient: κ = " << kappa << "\n\n";

        // =====================================================================
        // Create integrators
        // =====================================================================
        TestUtils::printProgress("Creating integrators...");

        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = new AdvDiffHierarchyIntegrator(
            "AdvDiffHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffHierarchyIntegrator"));

        adv_diff_integrator->setAdvectionVelocity(navier_stokes_integrator->getAdvectionVelocityVariable());
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);

        // =====================================================================
        // Create grid geometry
        // =====================================================================
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

        // =====================================================================
        // Create gridding algorithm
        // =====================================================================
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

        // =====================================================================
        // Set initial conditions
        // =====================================================================
        TestUtils::printProgress("Setting initial conditions...");

        if (input_db->keyExists("OdorInitialConditions"))
        {
            Pointer<CartGridFunction> C_init = new muParserCartGridFunction(
                "C_init", app_initializer->getComponentDatabase("OdorInitialConditions"), grid_geometry);
            adv_diff_integrator->setInitialConditions(C_init);
        }

        // =====================================================================
        // Set boundary conditions
        // =====================================================================
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        RobinBcCoefStrategy<NDIM>* C_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("OdorBcCoefs"))
        {
            C_bc_coef = new muParserRobinBcCoefs(
                "C_bc_coef", app_initializer->getComponentDatabase("OdorBcCoefs"), grid_geometry);
            adv_diff_integrator->setPhysicalBcCoef(C_bc_coef);
        }

        // =====================================================================
        // Set diffusion coefficient
        // =====================================================================
        adv_diff_integrator->setDiffusionCoefficient(kappa);

        // =====================================================================
        // Set visualization (reduced frequency for long runs)
        // =====================================================================
        const bool uses_visit = app_initializer->dumpVizData() && !app_initializer->getVisItDataWriter().isNull();
        if (uses_visit)
        {
            Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
            navier_stokes_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // =====================================================================
        // Initialize hierarchy
        // =====================================================================
        TestUtils::printProgress("Initializing patch hierarchy...");
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

        // =====================================================================
        // Get scalar variable index
        // =====================================================================
        const int C_idx = adv_diff_integrator->getTransportedQuantity()->getPatchDataIndex();
        logger.log("Scalar variable registered successfully");

        // Get initial statistics
        double initial_mass = ErrorCalculator::computeTotalMass(patch_hierarchy, C_idx);
        double C_min_init, C_max_init;
        ErrorCalculator::getMinMax(patch_hierarchy, C_idx, C_min_init, C_max_init);

        logger.logParameter("Initial mass", initial_mass);
        logger.logParameter("Initial min C", C_min_init);
        logger.logParameter("Initial max C", C_max_init);

        std::cout << "\n  Initial Condition Statistics:\n";
        std::cout << "    Mass: " << TestUtils::formatScientific(initial_mass) << "\n";
        std::cout << "    Min C: " << C_min_init << "\n";
        std::cout << "    Max C: " << C_max_init << "\n\n";

        // =====================================================================
        // Time-series tracking for long-term monitoring
        // =====================================================================
        std::vector<double> time_series;
        std::vector<double> mass_series;
        std::vector<double> min_series;
        std::vector<double> max_series;

        time_series.push_back(0.0);
        mass_series.push_back(initial_mass);
        min_series.push_back(C_min_init);
        max_series.push_back(C_max_init);

        // =====================================================================
        // Main time integration loop
        // =====================================================================
        TestUtils::printSectionSeparator("Long-Term Time Integration");

        double loop_time = navier_stokes_integrator->getIntegratorTime();
        double loop_time_end = navier_stokes_integrator->getEndTime();
        int iteration_num = navier_stokes_integrator->getIntegratorStep();

        bool solver_stable = true;
        int checkpoint_interval = std::max(1, num_steps / 20);  // 20 checkpoints

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) &&
               navier_stokes_integrator->stepsRemaining())
        {
            iteration_num = navier_stokes_integrator->getIntegratorStep();
            loop_time = navier_stokes_integrator->getIntegratorTime();

            if (iteration_num % checkpoint_interval == 0)
            {
                // Checkpoint monitoring
                double current_mass = ErrorCalculator::computeTotalMass(patch_hierarchy, C_idx);
                double mass_drift = TestUtils::computeMassDrift(initial_mass, current_mass);
                double C_min, C_max;
                ErrorCalculator::getMinMax(patch_hierarchy, C_idx, C_min, C_max);

                time_series.push_back(loop_time);
                mass_series.push_back(current_mass);
                min_series.push_back(C_min);
                max_series.push_back(C_max);

                std::cout << "Checkpoint " << iteration_num << "/" << num_steps
                          << ", t = " << std::setprecision(4) << std::fixed << loop_time
                          << ", Mass drift = " << std::scientific << std::setprecision(2) << mass_drift
                          << ", C ∈ [" << C_min << ", " << C_max << "]"
                          << std::endl;

                // Check for catastrophic failure
                bool has_nan = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);
                if (has_nan)
                {
                    TestUtils::printError("NaN/Inf detected at iteration " + std::to_string(iteration_num));
                    solver_stable = false;
                    break;
                }

                if (mass_drift > 0.1)
                {
                    TestUtils::printWarning("Excessive mass drift at t=" + std::to_string(loop_time));
                }
            }

            navier_stokes_integrator->advanceHierarchy(dt);
            loop_time += dt;
        }

        // =====================================================================
        // Write long-term time series
        // =====================================================================
        TestUtils::printProgress("Writing long-term time series...");

        std::ofstream long_run_file("long_run_series.dat");
        long_run_file << "# Long-term stability time series\n";
        long_run_file << "# time  mass  mass_drift  min_C  max_C\n";
        for (size_t i = 0; i < time_series.size(); ++i)
        {
            double drift = (i == 0) ? 0.0 : TestUtils::computeMassDrift(initial_mass, mass_series[i]);
            long_run_file << time_series[i] << "  "
                         << mass_series[i] << "  "
                         << drift << "  "
                         << min_series[i] << "  "
                         << max_series[i] << "\n";
        }
        long_run_file.close();
        logger.log("Long-term series written to long_run_series.dat");

        // =====================================================================
        // Final stability analysis
        // =====================================================================
        TestUtils::printSectionSeparator("Long-Term Stability Analysis");

        double final_mass = mass_series.back();
        double final_mass_drift = TestUtils::computeMassDrift(initial_mass, final_mass);
        double C_min_final = min_series.back();
        double C_max_final = max_series.back();
        bool has_nan_final = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);

        // Compute maximum drift over entire run
        double max_drift = 0.0;
        for (size_t i = 0; i < mass_series.size(); ++i)
        {
            double drift = TestUtils::computeMassDrift(initial_mass, mass_series[i]);
            max_drift = std::max(max_drift, fabs(drift));
        }

        std::cout << "\n  Long-Term Stability Summary:\n";
        std::cout << "    Duration: " << end_time << " time units (" << num_steps << " steps)\n";
        std::cout << "    Solver stable: " << (solver_stable ? "YES" : "NO") << "\n";
        std::cout << "    NaN/Inf detected: " << (has_nan_final ? "YES (FAIL)" : "NO (PASS)") << "\n";
        std::cout << "    Final mass drift: " << TestUtils::formatScientific(final_mass_drift) << "\n";
        std::cout << "    Max drift (entire run): " << TestUtils::formatScientific(max_drift) << "\n";
        std::cout << "    Final C range: [" << C_min_final << ", " << C_max_final << "]\n\n";

        logger.logParameter("Final mass", final_mass);
        logger.logParameter("Final mass drift", final_mass_drift);
        logger.logParameter("Max drift", max_drift);
        logger.logParameter("Final min C", C_min_final);
        logger.logParameter("Final max C", C_max_final);

        // =====================================================================
        // Test verdict
        // =====================================================================
        TestUtils::printSectionSeparator("Test Verdict");

        std::vector<std::string> checks;
        std::vector<bool> results;

        // Check 1: Solver remained stable (no NaN/Inf)
        checks.push_back("Solver stable (no NaN/Inf)");
        results.push_back(solver_stable && !has_nan_final);

        // Check 2: Mass conservation
        checks.push_back("Final mass drift < 0.01");
        results.push_back(final_mass_drift < 0.01);

        // Check 3: Maximum drift over entire run
        checks.push_back("Max drift (entire run) < 0.05");
        results.push_back(max_drift < 0.05);

        // Check 4: Solution remains bounded
        checks.push_back("Solution remains non-negative");
        results.push_back(C_min_final >= -0.01);

        // Check 5: No excessive growth
        checks.push_back("Max C does not grow excessively");
        results.push_back(C_max_final < C_max_init * 1.1);  // Allow 10% growth

        test_passed = true;
        for (size_t i = 0; i < checks.size(); ++i)
        {
            std::string status = results[i] ? "[PASS]" : "[FAIL]";
            std::cout << "  " << status << " " << checks[i] << "\n";
            test_passed = test_passed && results[i];
        }

        if (test_passed)
        {
            std::cout << "\n  Solver is stable for long-term integration!\n";
        }
        else
        {
            std::cout << "\n  WARNING: Solver shows instability in long-term integration\n";
        }

        if (C_bc_coef) delete C_bc_coef;

    } // End scope

    total_timer.stop();

    // =========================================================================
    // Final report
    // =========================================================================
    std::cout << "\n";
    std::cout << "Elapsed time: " << total_timer.getFormattedTime() << "\n";

    logger.logTestResult("Test 13: Long Run Stability", test_passed);
    logger.log("Total elapsed time: " + total_timer.getFormattedTime());

    TestUtils::printTestFooter("Test 13: Long Run", test_passed,
                              test_passed ? "Long-term stability validated" : "Stability issues detected");

    return test_passed ? 0 : 1;
}
