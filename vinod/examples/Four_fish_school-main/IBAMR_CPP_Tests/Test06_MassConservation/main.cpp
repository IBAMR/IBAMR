// =============================================================================
// TEST 06: Global Mass Conservation
// =============================================================================
// Validates that total scalar mass is conserved over time
// Monitors: M(t) = integral(C dV) with sources/sinks accounting
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
    TestUtils::printTestHeader("Global Mass Conservation Test", 6);

    Timer total_timer;
    total_timer.start();

    bool test_passed = true;
    ResultLogger logger("test06_results.txt");

    {
        // =====================================================================
        // Initialize application
        // =====================================================================
        TestUtils::printProgress("Initializing application...");
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test06.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get solver parameters
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");
        const double end_time = main_db->getDouble("END_TIME");
        const double dt = main_db->getDouble("DT");
        const int num_steps = static_cast<int>(std::ceil(end_time / dt));

        logger.logParameter("End time", end_time);
        logger.logParameter("Time step", dt);
        logger.logParameter("Number of steps", num_steps);

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
        if (input_db->keyExists("diffusion_coefficient"))
        {
            double kappa = input_db->getDouble("diffusion_coefficient");
            adv_diff_integrator->setDiffusionCoefficient(kappa);
            logger.logParameter("Diffusion coefficient", kappa);
        }

        // =====================================================================
        // Set visualization
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

        // =====================================================================
        // Get initial mass
        // =====================================================================
        double initial_mass = ErrorCalculator::computeTotalMass(patch_hierarchy, C_idx);
        logger.logParameter("Initial mass", initial_mass);

        TestUtils::printProgress("Initial mass: " + TestUtils::formatScientific(initial_mass));

        // =====================================================================
        // Time-series mass tracking
        // =====================================================================
        std::vector<double> time_series;
        std::vector<double> mass_series;
        std::vector<double> drift_series;

        time_series.push_back(0.0);
        mass_series.push_back(initial_mass);
        drift_series.push_back(0.0);

        // =====================================================================
        // Main time integration loop with mass tracking
        // =====================================================================
        TestUtils::printSectionSeparator("Time Integration - Mass Conservation");

        double loop_time = navier_stokes_integrator->getIntegratorTime();
        double loop_time_end = navier_stokes_integrator->getEndTime();
        int iteration_num = navier_stokes_integrator->getIntegratorStep();

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) &&
               navier_stokes_integrator->stepsRemaining())
        {
            iteration_num = navier_stokes_integrator->getIntegratorStep();
            loop_time = navier_stokes_integrator->getIntegratorTime();

            navier_stokes_integrator->advanceHierarchy(dt);
            loop_time += dt;

            // Compute mass at this timestep
            double current_mass = ErrorCalculator::computeTotalMass(patch_hierarchy, C_idx);
            double mass_drift = TestUtils::computeMassDrift(initial_mass, current_mass);

            // Store time series
            time_series.push_back(loop_time);
            mass_series.push_back(current_mass);
            drift_series.push_back(mass_drift);

            if (iteration_num % 10 == 0)
            {
                std::cout << "Iteration " << iteration_num
                          << ", t = " << loop_time
                          << ", Mass = " << TestUtils::formatScientific(current_mass)
                          << ", Drift = " << TestUtils::formatScientific(mass_drift)
                          << std::endl;
            }

            // Check for excessive drift
            if (mass_drift > 1.0e-4)
            {
                TestUtils::printWarning("Large mass drift detected at t=" +
                                       std::to_string(loop_time) + ": " +
                                       TestUtils::formatScientific(mass_drift));
            }
        }

        // =====================================================================
        // Write mass time series to file
        // =====================================================================
        TestUtils::printProgress("Writing mass time series...");

        std::ofstream mass_file("mass_conservation_series.dat");
        mass_file << "# Time-series mass conservation data\n";
        mass_file << "# time  mass  drift  relative_drift\n";
        for (size_t i = 0; i < time_series.size(); ++i)
        {
            double rel_drift = (initial_mass > 0) ? drift_series[i] / initial_mass : 0.0;
            mass_file << time_series[i] << "  "
                     << mass_series[i] << "  "
                     << drift_series[i] << "  "
                     << rel_drift << "\n";
        }
        mass_file.close();
        logger.log("Mass time series written to mass_conservation_series.dat");

        // =====================================================================
        // Analyze mass conservation
        // =====================================================================
        TestUtils::printSectionSeparator("Mass Conservation Analysis");

        double final_mass = mass_series.back();
        double max_drift = *std::max_element(drift_series.begin(), drift_series.end());
        double mean_drift = 0.0;
        for (double d : drift_series) mean_drift += d;
        mean_drift /= drift_series.size();

        std::cout << "  Initial mass:  " << TestUtils::formatScientific(initial_mass) << "\n";
        std::cout << "  Final mass:    " << TestUtils::formatScientific(final_mass) << "\n";
        std::cout << "  Final drift:   " << TestUtils::formatScientific(drift_series.back()) << "\n";
        std::cout << "  Max drift:     " << TestUtils::formatScientific(max_drift) << "\n";
        std::cout << "  Mean drift:    " << TestUtils::formatScientific(mean_drift) << "\n";

        logger.logParameter("Final mass", final_mass);
        logger.logParameter("Final drift", drift_series.back());
        logger.logParameter("Max drift", max_drift);
        logger.logParameter("Mean drift", mean_drift);

        // =====================================================================
        // Test verdict
        // =====================================================================
        TestUtils::printSectionSeparator("Test Verdict");

        std::vector<std::string> checks;
        std::vector<bool> results;

        // Check 1: Final drift
        checks.push_back("Final drift < 1e-6");
        results.push_back(drift_series.back() < 1.0e-6);

        // Check 2: Max drift
        checks.push_back("Max drift < 1e-5");
        results.push_back(max_drift < 1.0e-5);

        // Check 3: Mean drift
        checks.push_back("Mean drift < 1e-6");
        results.push_back(mean_drift < 1.0e-6);

        // Check 4: No NaN
        bool has_nan = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);
        checks.push_back("No NaN/Inf");
        results.push_back(!has_nan);

        test_passed = true;
        for (size_t i = 0; i < checks.size(); ++i)
        {
            std::string status = results[i] ? "[PASS]" : "[FAIL]";
            std::cout << "  " << status << " " << checks[i] << "\n";
            test_passed = test_passed && results[i];
        }

        if (C_bc_coef) delete C_bc_coef;

    } // End scope

    total_timer.stop();

    // =========================================================================
    // Final report
    // =========================================================================
    std::cout << "\n";
    std::cout << "Elapsed time: " << total_timer.getFormattedTime() << "\n";

    logger.logTestResult("Test 06: Global Mass Conservation", test_passed);
    logger.log("Total elapsed time: " + total_timer.getFormattedTime());

    TestUtils::printTestFooter("Test 06: Mass Conservation", test_passed,
                              test_passed ? "Mass conservation validated" : "Mass drift exceeds tolerance");

    return test_passed ? 0 : 1;
}
