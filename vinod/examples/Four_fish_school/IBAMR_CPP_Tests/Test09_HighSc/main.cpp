// =============================================================================
// TEST 09: High Schmidt Number Test
// =============================================================================
// Validates solver stability and accuracy for high Schmidt numbers
// Sc = ν/κ = [0.7 (air), 100, 340 (water), 1000]
// Tests implicit diffusion solver at extreme parameter ranges
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

using namespace TestUtilities;

// Structure to store test results for each Sc
struct SchmidtTestResult
{
    double Sc;
    double kappa;
    double nu;
    bool stable;
    bool no_negatives;
    bool no_nan;
    double final_mass_drift;
    double min_C;
    double max_C;
};

int main(int argc, char* argv[])
{
    // =========================================================================
    // Initialize IBAMR/SAMRAI/PETSc
    // =========================================================================
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    TestUtils::printTestHeader("High Schmidt Number Validation", 9);

    Timer total_timer;
    total_timer.start();

    bool test_passed = true;
    ResultLogger logger("test09_results.txt");

    // Schmidt numbers to test (from user's test plan Section 6.1)
    std::vector<double> schmidt_numbers = {0.7, 100.0, 340.0, 1000.0};
    std::vector<SchmidtTestResult> results;

    logger.log("Testing Schmidt numbers: 0.7 (air), 100, 340 (water), 1000");

    {
        // =====================================================================
        // Initialize application (common setup)
        // =====================================================================
        TestUtils::printProgress("Initializing application...");
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test09.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get solver parameters
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");
        const double end_time = main_db->getDouble("END_TIME");
        const double dt = main_db->getDouble("DT");

        // Get base viscosity
        const double nu = input_db->getDoubleWithDefault("kinematic_viscosity", 0.01);

        logger.logParameter("Base kinematic viscosity", nu);
        logger.logParameter("End time", end_time);
        logger.logParameter("Time step", dt);

        // =====================================================================
        // Loop over each Schmidt number
        // =====================================================================
        for (size_t sc_idx = 0; sc_idx < schmidt_numbers.size(); ++sc_idx)
        {
            const double Sc = schmidt_numbers[sc_idx];
            const double kappa = nu / Sc;  // κ = ν / Sc

            TestUtils::printSectionSeparator("Schmidt Number Test: Sc = " + std::to_string(Sc));
            std::cout << "  Kinematic viscosity: " << nu << "\n";
            std::cout << "  Diffusion coefficient: " << kappa << "\n";
            std::cout << "  Schmidt number: " << Sc << "\n\n";

            SchmidtTestResult sc_result;
            sc_result.Sc = Sc;
            sc_result.kappa = kappa;
            sc_result.nu = nu;

            // =================================================================
            // Create integrators for this Sc test
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
            // Set diffusion coefficient for this Sc
            // =================================================================
            adv_diff_integrator->setDiffusionCoefficient(kappa);
            std::cout << "  Set diffusion coefficient: κ = " << kappa
                      << " (Sc = " << Sc << ")\n";

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

            // =================================================================
            // Get scalar variable index
            // =================================================================
            const int C_idx = adv_diff_integrator->getTransportedQuantity()->getPatchDataIndex();

            // Get initial mass
            double initial_mass = ErrorCalculator::computeTotalMass(patch_hierarchy, C_idx);
            std::cout << "  Initial mass: " << initial_mass << "\n\n";

            // =================================================================
            // Time integration loop
            // =================================================================
            std::cout << "  Running time integration...\n";

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
                    std::cout << "    Iteration " << iteration_num
                              << ", t = " << loop_time << std::endl;
                }

                navier_stokes_integrator->advanceHierarchy(dt);
                loop_time += dt;

                // Check for instability
                if (iteration_num % 5 == 0)
                {
                    bool has_nan = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);
                    if (has_nan)
                    {
                        TestUtils::printError("NaN/Inf detected at Sc=" + std::to_string(Sc) +
                                             ", iteration " + std::to_string(iteration_num));
                        solver_stable = false;
                        break;
                    }
                }
            }

            // =================================================================
            // Analyze results for this Sc
            // =================================================================
            std::cout << "\n  Analysis for Sc = " << Sc << ":\n";

            double C_min, C_max;
            ErrorCalculator::getMinMax(patch_hierarchy, C_idx, C_min, C_max);
            double final_mass = ErrorCalculator::computeTotalMass(patch_hierarchy, C_idx);
            double mass_drift = TestUtils::computeMassDrift(initial_mass, final_mass);
            int neg_count = ErrorCalculator::checkForNegatives(patch_hierarchy, C_idx);
            bool has_nan = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);

            sc_result.stable = solver_stable && !has_nan;
            sc_result.no_negatives = (neg_count == 0);
            sc_result.no_nan = !has_nan;
            sc_result.final_mass_drift = mass_drift;
            sc_result.min_C = C_min;
            sc_result.max_C = C_max;

            std::cout << "    Min C: " << C_min << "\n";
            std::cout << "    Max C: " << C_max << "\n";
            std::cout << "    Mass drift: " << TestUtils::formatScientific(mass_drift) << "\n";
            std::cout << "    Negative cells: " << neg_count << "\n";
            std::cout << "    Stable: " << (sc_result.stable ? "YES" : "NO") << "\n";
            std::cout << "    No NaN/Inf: " << (sc_result.no_nan ? "YES" : "NO") << "\n\n";

            results.push_back(sc_result);

            if (C_bc_coef) delete C_bc_coef;

        } // End Sc loop

        // =====================================================================
        // Summary of all Sc tests
        // =====================================================================
        TestUtils::printSectionSeparator("High-Sc Test Summary");

        std::cout << "\n";
        std::cout << "  Sc      κ           Stable  No-Neg  No-NaN  Mass-Drift\n";
        std::cout << "  ----    --------    ------  ------  ------  ----------\n";

        for (const auto& res : results)
        {
            std::cout << "  " << std::setw(6) << res.Sc
                      << "  " << std::setw(10) << std::scientific << std::setprecision(2) << res.kappa
                      << "  " << (res.stable ? "  YES  " : "  NO   ")
                      << "  " << (res.no_negatives ? "  YES " : "  NO  ")
                      << "  " << (res.no_nan ? "  YES " : "  NO  ")
                      << "  " << std::setw(10) << std::scientific << res.final_mass_drift
                      << "\n";

            logger.logParameter("Sc_" + std::to_string((int)res.Sc) + "_stable",
                               res.stable ? 1.0 : 0.0);
            logger.logParameter("Sc_" + std::to_string((int)res.Sc) + "_mass_drift",
                               res.final_mass_drift);
        }
        std::cout << "\n";

        // =====================================================================
        // Test verdict
        // =====================================================================
        TestUtils::printSectionSeparator("Test Verdict");

        test_passed = true;
        for (const auto& res : results)
        {
            bool sc_passed = res.stable && res.no_nan && (res.final_mass_drift < 0.01);

            std::string status = sc_passed ? "[PASS]" : "[FAIL]";
            std::cout << "  " << status << " Sc = " << res.Sc;
            if (!res.stable) std::cout << " (UNSTABLE)";
            if (!res.no_nan) std::cout << " (NaN)";
            if (res.final_mass_drift > 0.01) std::cout << " (MASS DRIFT)";
            std::cout << "\n";

            test_passed = test_passed && sc_passed;
        }

        if (test_passed)
        {
            std::cout << "\n  All Schmidt numbers tested successfully!\n";
            std::cout << "  Solver is stable for Sc = 0.7 (air) to Sc = 1000\n";
        }

    } // End scope

    total_timer.stop();

    // =========================================================================
    // Final report
    // =========================================================================
    std::cout << "\n";
    std::cout << "Elapsed time: " << total_timer.getFormattedTime() << "\n";

    logger.logTestResult("Test 09: High Schmidt Number", test_passed);
    logger.log("Total elapsed time: " + total_timer.getFormattedTime());

    TestUtils::printTestFooter("Test 09: High Schmidt", test_passed,
                              test_passed ? "Solver validated for high-Sc regime" :
                                           "Solver unstable at some Sc values");

    return test_passed ? 0 : 1;
}
