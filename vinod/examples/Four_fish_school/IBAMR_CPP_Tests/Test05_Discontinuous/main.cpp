// =============================================================================
// TEST 05: Discontinuous Initial Condition - Stability Test
// =============================================================================
// Tests solver stability with top-hat (step) function initial condition
// Validates: Monotonicity, no oscillations, no overshoots, robustness
// Critical for ensuring solver handles sharp concentration fronts
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

using namespace TestUtilities;

int main(int argc, char* argv[])
{
    // =========================================================================
    // Initialize IBAMR/SAMRAI/PETSc
    // =========================================================================
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    TestUtils::printTestHeader("Discontinuous Initial Condition - Stability", 5);

    Timer total_timer;
    total_timer.start();

    bool test_passed = true;
    ResultLogger logger("test05_results.txt");

    {
        // =====================================================================
        // Initialize application
        // =====================================================================
        TestUtils::printProgress("Initializing application...");
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test05.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get solver parameters
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");
        const double end_time = main_db->getDouble("END_TIME");
        const double dt = main_db->getDouble("DT");

        // Get physical parameters
        const double kappa = input_db->getDouble("diffusion_coefficient");

        // Top-hat parameters
        const double x0 = input_db->getDoubleWithDefault("tophat_x0", 0.0);
        const double y0 = input_db->getDoubleWithDefault("tophat_y0", 0.0);
        const double radius = input_db->getDoubleWithDefault("tophat_radius", 0.3);
        const double C0 = input_db->getDoubleWithDefault("tophat_amplitude", 1.0);

        logger.logParameter("Diffusion coefficient", kappa);
        logger.logParameter("Top-hat center x0", x0);
        logger.logParameter("Top-hat center y0", y0);
        logger.logParameter("Top-hat radius", radius);
        logger.logParameter("Top-hat amplitude", C0);
        logger.logParameter("End time", end_time);
        logger.logParameter("Time step", dt);

        std::cout << "\n  Discontinuous IC Test Configuration:\n";
        std::cout << "    Initial condition: Top-hat (step function)\n";
        std::cout << "    Center: (" << x0 << ", " << y0 << "), Radius: " << radius << "\n";
        std::cout << "    Amplitude: " << C0 << " inside, 0.0 outside (sharp discontinuity)\n";
        std::cout << "    Diffusion coefficient: " << kappa << "\n\n";

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
        // Set initial conditions (top-hat / step function)
        // =====================================================================
        TestUtils::printProgress("Setting discontinuous top-hat initial condition...");

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

        // Get initial statistics
        double C_min_init, C_max_init;
        ErrorCalculator::getMinMax(patch_hierarchy, C_idx, C_min_init, C_max_init);
        double initial_mass = ErrorCalculator::computeTotalMass(patch_hierarchy, C_idx);

        std::cout << "\n  Initial Condition Statistics:\n";
        std::cout << "    Min C: " << C_min_init << "\n";
        std::cout << "    Max C: " << C_max_init << " (expected: " << C0 << ")\n";
        std::cout << "    Total mass: " << TestUtils::formatScientific(initial_mass) << "\n\n";

        logger.logParameter("Initial min C", C_min_init);
        logger.logParameter("Initial max C", C_max_init);
        logger.logParameter("Initial mass", initial_mass);

        // =====================================================================
        // Main time integration loop
        // =====================================================================
        TestUtils::printSectionSeparator("Time Integration - Discontinuous IC");

        double loop_time = navier_stokes_integrator->getIntegratorTime();
        double loop_time_end = navier_stokes_integrator->getEndTime();
        int iteration_num = navier_stokes_integrator->getIntegratorStep();

        bool solver_stable = true;

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) &&
               navier_stokes_integrator->stepsRemaining())
        {
            iteration_num = navier_stokes_integrator->getIntegratorStep();
            loop_time = navier_stokes_integrator->getIntegratorTime();

            if (iteration_num % 10 == 0)
            {
                std::cout << "Iteration " << iteration_num
                          << ", t = " << loop_time << std::endl;
            }

            navier_stokes_integrator->advanceHierarchy(dt);
            loop_time += dt;

            // Check for instability during integration
            if (iteration_num % 5 == 0)
            {
                bool has_nan = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);
                if (has_nan)
                {
                    TestUtils::printError("NaN/Inf detected at iteration " + std::to_string(iteration_num));
                    solver_stable = false;
                    break;
                }
            }
        }

        // =====================================================================
        // Stability analysis
        // =====================================================================
        TestUtils::printSectionSeparator("Stability Analysis");
        TestUtils::printProgress("Analyzing solution for oscillations and overshoots...");

        double C_min, C_max;
        ErrorCalculator::getMinMax(patch_hierarchy, C_idx, C_min, C_max);
        int neg_count = ErrorCalculator::checkForNegatives(patch_hierarchy, C_idx);
        bool has_nan = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);

        std::cout << "\n  Final Solution Statistics:\n";
        std::cout << "    Min C: " << C_min << " (should be >= 0.0)\n";
        std::cout << "    Max C: " << C_max << " (should be <= " << C0 << ")\n";
        std::cout << "    Negative cells: " << neg_count << "\n";
        std::cout << "    NaN/Inf detected: " << (has_nan ? "YES (FAIL)" : "NO (PASS)") << "\n";
        std::cout << "    Solver stable: " << (solver_stable ? "YES" : "NO") << "\n\n";

        logger.logParameter("Final min C", C_min);
        logger.logParameter("Final max C", C_max);
        logger.logParameter("Negative cells", neg_count);

        // Check for overshoots/undershoots (Gibbs phenomenon)
        double undershoot = std::min(0.0, C_min);
        double overshoot = std::max(0.0, C_max - C0);

        std::cout << "  Oscillation Analysis:\n";
        std::cout << "    Undershoot: " << undershoot << " (should be 0.0)\n";
        std::cout << "    Overshoot: " << overshoot << " (should be 0.0)\n\n";

        logger.logParameter("Undershoot", undershoot);
        logger.logParameter("Overshoot", overshoot);

        // =====================================================================
        // Mass conservation check
        // =====================================================================
        double final_mass = ErrorCalculator::computeTotalMass(patch_hierarchy, C_idx);
        double mass_drift = TestUtils::computeMassDrift(initial_mass, final_mass);

        std::cout << "  Mass Conservation:\n";
        std::cout << "    Initial mass: " << TestUtils::formatScientific(initial_mass) << "\n";
        std::cout << "    Final mass:   " << TestUtils::formatScientific(final_mass) << "\n";
        std::cout << "    Mass drift:   " << TestUtils::formatScientific(mass_drift) << "\n\n";

        logger.logParameter("Final mass", final_mass);
        logger.logParameter("Mass drift", mass_drift);

        // =====================================================================
        // Test verdict
        // =====================================================================
        TestUtils::printSectionSeparator("Test Verdict");

        std::vector<std::string> checks;
        std::vector<bool> results;

        // Check 1: Solver remained stable
        checks.push_back("Solver stable (no NaN/Inf)");
        results.push_back(solver_stable && !has_nan);

        // Check 2: No negative values (monotonicity)
        checks.push_back("No negative values");
        results.push_back(neg_count == 0);

        // Check 3: No significant undershoot
        checks.push_back("Undershoot < 0.01 (numerical oscillations)");
        results.push_back(fabs(undershoot) < 0.01);

        // Check 4: No significant overshoot
        checks.push_back("Overshoot < 0.05 (Gibbs phenomenon)");
        results.push_back(overshoot < 0.05);

        // Check 5: Mass conservation
        checks.push_back("Mass drift < 1e-5");
        results.push_back(mass_drift < 1.0e-5);

        // Check 6: Solution bounded
        checks.push_back("Solution bounded [0, C0]");
        results.push_back(C_min >= -0.01 && C_max <= C0 + 0.05);

        test_passed = true;
        for (size_t i = 0; i < checks.size(); ++i)
        {
            std::string status = results[i] ? "[PASS]" : "[FAIL]";
            std::cout << "  " << status << " " << checks[i] << "\n";
            test_passed = test_passed && results[i];
        }

        if (test_passed)
        {
            std::cout << "\n  Solver handles discontinuous initial conditions robustly!\n";
        }
        else
        {
            std::cout << "\n  WARNING: Solver shows oscillations or instability with discontinuous IC\n";
        }

        if (C_bc_coef) delete C_bc_coef;

    } // End scope

    total_timer.stop();

    // =========================================================================
    // Final report
    // =========================================================================
    std::cout << "\n";
    std::cout << "Elapsed time: " << total_timer.getFormattedTime() << "\n";

    logger.logTestResult("Test 05: Discontinuous Initial Condition", test_passed);
    logger.log("Total elapsed time: " + total_timer.getFormattedTime());

    TestUtils::printTestFooter("Test 05: Discontinuous IC", test_passed,
                              test_passed ? "Solver stable with sharp fronts" : "Oscillations or instability detected");

    return test_passed ? 0 : 1;
}
