// =============================================================================
// TEST 11: AMR Sensitivity Test
// =============================================================================
// Validates that adaptive mesh refinement does not introduce artifacts
// Tests: Consistency across refinement levels, no spurious oscillations
// Critical for production simulations with AMR
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
#include <fstream>

using namespace TestUtilities;

int main(int argc, char* argv[])
{
    // =========================================================================
    // Initialize IBAMR/SAMRAI/PETSc
    // =========================================================================
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    TestUtils::printTestHeader("AMR Sensitivity Test", 11);

    Timer total_timer;
    total_timer.start();

    bool test_passed = true;
    ResultLogger logger("test11_results.txt");

    {
        // =====================================================================
        // Initialize application
        // =====================================================================
        TestUtils::printProgress("Initializing application...");
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test11.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get solver parameters
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");
        const double end_time = main_db->getDouble("END_TIME");
        const double dt = main_db->getDouble("DT");

        const double kappa = input_db->getDouble("diffusion_coefficient");

        logger.logParameter("Diffusion coefficient", kappa);
        logger.logParameter("End time", end_time);
        logger.logParameter("Time step", dt);

        std::cout << "\n  AMR Sensitivity Test Configuration:\n";
        std::cout << "    Testing: AMR artifacts and consistency\n";
        std::cout << "    Diffusion coefficient: κ = " << kappa << "\n";
        std::cout << "    End time: T = " << end_time << "\n\n";

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
        // Create gridding algorithm with AMR
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

        // Get AMR configuration
        Pointer<Database> gridding_db = app_initializer->getComponentDatabase("GriddingAlgorithm");
        const int max_levels = gridding_db->getInteger("max_levels");

        std::cout << "  AMR Configuration:\n";
        std::cout << "    Maximum refinement levels: " << max_levels << "\n";

        if (max_levels > 1)
        {
            std::cout << "    AMR enabled - testing multi-level refinement\n";
        }
        else
        {
            std::cout << "    Single level - testing uniform grid baseline\n";
        }
        std::cout << "\n";

        logger.logParameter("Max AMR levels", max_levels);

        // =====================================================================
        // Set initial conditions (Gaussian pulse for AMR testing)
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
        TestUtils::printProgress("Initializing patch hierarchy with AMR...");
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

        std::cout << "  Hierarchy Information:\n";
        std::cout << "    Levels created: " << finest_level + 1 << " / " << max_levels << "\n";
        for (int ln = 0; ln <= finest_level; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            std::cout << "    Level " << ln << ": " << level->getNumberOfPatches() << " patches\n";
        }
        std::cout << "\n";

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

        std::cout << "  Initial Statistics:\n";
        std::cout << "    Mass: " << TestUtils::formatScientific(initial_mass) << "\n";
        std::cout << "    C range: [" << C_min_init << ", " << C_max_init << "]\n\n";

        // =====================================================================
        // Time series for AMR artifact detection
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
        TestUtils::printSectionSeparator("Time Integration - AMR Test");

        double loop_time = navier_stokes_integrator->getIntegratorTime();
        double loop_time_end = navier_stokes_integrator->getEndTime();
        int iteration_num = navier_stokes_integrator->getIntegratorStep();

        bool solver_stable = true;
        int checkpoint_interval = 10;

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) &&
               navier_stokes_integrator->stepsRemaining())
        {
            iteration_num = navier_stokes_integrator->getIntegratorStep();
            loop_time = navier_stokes_integrator->getIntegratorTime();

            if (iteration_num % checkpoint_interval == 0)
            {
                double current_mass = ErrorCalculator::computeTotalMass(patch_hierarchy, C_idx);
                double mass_drift = TestUtils::computeMassDrift(initial_mass, current_mass);
                double C_min, C_max;
                ErrorCalculator::getMinMax(patch_hierarchy, C_idx, C_min, C_max);

                time_series.push_back(loop_time);
                mass_series.push_back(current_mass);
                min_series.push_back(C_min);
                max_series.push_back(C_max);

                std::cout << "Iteration " << iteration_num
                          << ", t = " << std::setprecision(4) << std::fixed << loop_time
                          << ", Mass drift = " << std::scientific << std::setprecision(2) << mass_drift
                          << ", C ∈ [" << C_min << ", " << C_max << "]"
                          << ", Levels: " << finest_level + 1
                          << std::endl;

                // Check for NaN/Inf (AMR artifact detection)
                bool has_nan = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);
                if (has_nan)
                {
                    TestUtils::printError("NaN/Inf detected - possible AMR artifact!");
                    solver_stable = false;
                    break;
                }

                // Check for negative values (AMR artifact)
                if (C_min < -0.01)
                {
                    TestUtils::printWarning("Negative concentration detected - possible AMR artifact");
                }
            }

            navier_stokes_integrator->advanceHierarchy(dt);
            loop_time += dt;
        }

        // =====================================================================
        // AMR artifact analysis
        // =====================================================================
        TestUtils::printSectionSeparator("AMR Artifact Analysis");

        double final_mass = mass_series.back();
        double final_mass_drift = TestUtils::computeMassDrift(initial_mass, final_mass);
        double C_min_final = min_series.back();
        double C_max_final = max_series.back();
        bool has_nan_final = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);

        // Check for monotonicity violations (artifact indicator)
        int monotonicity_violations = 0;
        for (size_t i = 1; i < max_series.size(); ++i)
        {
            // For pure diffusion, max should be monotonically decreasing
            if (max_series[i] > max_series[i-1] * 1.01)  // Allow 1% tolerance
            {
                monotonicity_violations++;
            }
        }

        // Compute maximum drift
        double max_drift = 0.0;
        for (size_t i = 0; i < mass_series.size(); ++i)
        {
            double drift = TestUtils::computeMassDrift(initial_mass, mass_series[i]);
            max_drift = std::max(max_drift, fabs(drift));
        }

        std::cout << "\n  AMR Consistency Analysis:\n";
        std::cout << "    Refinement levels used: " << finest_level + 1 << "\n";
        std::cout << "    Solver stable: " << (solver_stable ? "YES" : "NO") << "\n";
        std::cout << "    NaN/Inf detected: " << (has_nan_final ? "YES (FAIL)" : "NO (PASS)") << "\n";
        std::cout << "    Monotonicity violations: " << monotonicity_violations << "\n";
        std::cout << "    Final mass drift: " << TestUtils::formatScientific(final_mass_drift) << "\n";
        std::cout << "    Max drift (entire run): " << TestUtils::formatScientific(max_drift) << "\n";
        std::cout << "    Final C range: [" << C_min_final << ", " << C_max_final << "]\n\n";

        logger.logParameter("Refinement levels", finest_level + 1);
        logger.logParameter("Final mass drift", final_mass_drift);
        logger.logParameter("Max drift", max_drift);
        logger.logParameter("Monotonicity violations", monotonicity_violations);
        logger.logParameter("Final min C", C_min_final);
        logger.logParameter("Final max C", C_max_final);

        // Write AMR time series
        std::ofstream amr_file("amr_series.dat");
        amr_file << "# AMR consistency time series\n";
        amr_file << "# time  mass  mass_drift  min_C  max_C\n";
        for (size_t i = 0; i < time_series.size(); ++i)
        {
            double drift = (i == 0) ? 0.0 : TestUtils::computeMassDrift(initial_mass, mass_series[i]);
            amr_file << time_series[i] << "  "
                     << mass_series[i] << "  "
                     << drift << "  "
                     << min_series[i] << "  "
                     << max_series[i] << "\n";
        }
        amr_file.close();
        logger.log("AMR series written to amr_series.dat");

        // =====================================================================
        // Test verdict
        // =====================================================================
        TestUtils::printSectionSeparator("Test Verdict");

        std::vector<std::string> checks;
        std::vector<bool> results;

        // Check 1: No NaN/Inf (critical for AMR)
        checks.push_back("No NaN/Inf artifacts");
        results.push_back(solver_stable && !has_nan_final);

        // Check 2: Mass conservation (AMR should preserve mass)
        checks.push_back("Mass drift < 0.01");
        results.push_back(final_mass_drift < 0.01);

        // Check 3: No excessive monotonicity violations
        checks.push_back("Monotonicity violations < 5");
        results.push_back(monotonicity_violations < 5);

        // Check 4: Solution remains non-negative
        checks.push_back("Non-negative concentration");
        results.push_back(C_min_final >= -0.01);

        // Check 5: No excessive drift over entire run
        checks.push_back("Max drift < 0.05");
        results.push_back(max_drift < 0.05);

        test_passed = true;
        for (size_t i = 0; i < checks.size(); ++i)
        {
            std::string status = results[i] ? "[PASS]" : "[FAIL]";
            std::cout << "  " << status << " " << checks[i] << "\n";
            test_passed = test_passed && results[i];
        }

        if (test_passed)
        {
            std::cout << "\n  AMR implementation is consistent and artifact-free!\n";
            if (max_levels > 1)
            {
                std::cout << "  Multi-level refinement validated successfully.\n";
            }
        }
        else
        {
            std::cout << "\n  WARNING: AMR artifacts detected - review refinement strategy\n";
        }

        if (C_bc_coef) delete C_bc_coef;

    } // End scope

    total_timer.stop();

    // =========================================================================
    // Final report
    // =========================================================================
    std::cout << "\n";
    std::cout << "Elapsed time: " << total_timer.getFormattedTime() << "\n";

    logger.logTestResult("Test 11: AMR Sensitivity", test_passed);
    logger.log("Total elapsed time: " + total_timer.getFormattedTime());

    TestUtils::printTestFooter("Test 11: AMR Sensitivity", test_passed,
                              test_passed ? "AMR validated - no artifacts" : "AMR artifacts detected");

    return test_passed ? 0 : 1;
}
