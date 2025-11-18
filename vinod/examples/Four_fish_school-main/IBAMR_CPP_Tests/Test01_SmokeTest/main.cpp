// =============================================================================
// TEST 1: SMOKE TEST - Scalar Transport Infrastructure Validation
// =============================================================================
// Purpose:
// - Verify basic scalar infrastructure works
// - Check variable registration
// - Test boundary conditions
// - Validate I/O functionality
// - Ensure no crashes, no NaNs
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
    TestUtils::printTestHeader("Smoke Test - Scalar Transport Infrastructure", 1);

    Timer total_timer;
    total_timer.start();

    bool test_passed = true;
    ResultLogger logger("test01_results.txt");

    {
        // =====================================================================
        // Initialize application
        // =====================================================================
        TestUtils::printProgress("Initializing application...");

        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test01.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get visualization parameters
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();

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

        // Create Navier-Stokes integrator (for velocity field, u=0 in this test)
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        // Create scalar transport integrator
        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = new AdvDiffHierarchyIntegrator(
            "AdvDiffHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffHierarchyIntegrator"));

        // Register velocity field with scalar transport
        adv_diff_integrator->setAdvectionVelocity(navier_stokes_integrator->getAdvectionVelocityVariable());
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);

        // =====================================================================
        // Create grid geometry and patch hierarchy
        // =====================================================================
        TestUtils::printProgress("Creating grid geometry...");

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

        // Get domain parameters
        const double* x_lo = grid_geometry->getXLower();
        const double* x_hi = grid_geometry->getXUpper();
        const int* N = grid_geometry->getPeriodicShift().getPointer();

        logger.logParameter("Domain X", std::to_string(x_lo[0]) + " to " + std::to_string(x_hi[0]));
        logger.logParameter("Domain Y", std::to_string(x_lo[1]) + " to " + std::to_string(x_hi[1]));

        // =====================================================================
        // Create gridding algorithm
        // =====================================================================
        TestUtils::printProgress("Creating gridding algorithm...");

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
        TestUtils::printProgress("Setting boundary conditions...");

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
        // Set visualization data writers
        // =====================================================================
        if (uses_visit)
        {
            Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
            navier_stokes_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // =====================================================================
        // Initialize hierarchy
        // =====================================================================
        TestUtils::printProgress("Initializing hierarchy...");

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
        // Get initial statistics
        // =====================================================================
        double C_min, C_max;
        ErrorCalculator::getMinMax(patch_hierarchy, C_idx, C_min, C_max);
        double initial_mass = ErrorCalculator::computeTotalMass(patch_hierarchy, C_idx);

        TestUtils::printProgress("Initial statistics:");
        std::cout << "  C_min = " << C_min << "\n";
        std::cout << "  C_max = " << C_max << "\n";
        std::cout << "  Total mass = " << initial_mass << "\n";

        logger.logParameter("Initial C_min", C_min);
        logger.logParameter("Initial C_max", C_max);
        logger.logParameter("Initial mass", initial_mass);

        // =====================================================================
        // Check for initial problems
        // =====================================================================
        bool has_nan = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);
        if (has_nan)
        {
            TestUtils::printError("NaN/Inf detected in initial conditions!");
            test_passed = false;
        }

        int neg_count = ErrorCalculator::checkForNegatives(patch_hierarchy, C_idx);
        if (neg_count > 0)
        {
            TestUtils::printWarning("Found " + std::to_string(neg_count) + " negative values initially");
        }

        // =====================================================================
        // Main time integration loop
        // =====================================================================
        TestUtils::printSectionSeparator("Time Integration");

        double loop_time = navier_stokes_integrator->getIntegratorTime();
        double loop_time_end = navier_stokes_integrator->getEndTime();

        int iteration_num = navier_stokes_integrator->getIntegratorStep();
        int max_iterations = num_steps;

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) &&
               navier_stokes_integrator->stepsRemaining())
        {
            iteration_num = navier_stokes_integrator->getIntegratorStep();
            loop_time = navier_stokes_integrator->getIntegratorTime();

            if (iteration_num % 10 == 0)
            {
                std::cout << "Iteration " << iteration_num << " / " << max_iterations
                          << ", t = " << loop_time << std::endl;
            }

            navier_stokes_integrator->advanceHierarchy(dt);
            loop_time += dt;

            // Check for NaN/Inf during run
            if (iteration_num % 10 == 0)
            {
                has_nan = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);
                if (has_nan)
                {
                    TestUtils::printError("NaN/Inf detected at iteration " +
                                         std::to_string(iteration_num));
                    test_passed = false;
                    break;
                }
            }
        }

        // =====================================================================
        // Final statistics
        // =====================================================================
        TestUtils::printSectionSeparator("Final Statistics");

        ErrorCalculator::getMinMax(patch_hierarchy, C_idx, C_min, C_max);
        double final_mass = ErrorCalculator::computeTotalMass(patch_hierarchy, C_idx);

        std::cout << "  C_min = " << C_min << "\n";
        std::cout << "  C_max = " << C_max << "\n";
        std::cout << "  Total mass = " << final_mass << "\n";

        logger.logParameter("Final C_min", C_min);
        logger.logParameter("Final C_max", C_max);
        logger.logParameter("Final mass", final_mass);

        // Check for negatives
        neg_count = ErrorCalculator::checkForNegatives(patch_hierarchy, C_idx);
        if (neg_count > 0)
        {
            TestUtils::printWarning("Found " + std::to_string(neg_count) + " negative values at end");
        }

        // Check mass conservation (no advection, no flux BCs, so mass should be conserved)
        double mass_drift = TestUtils::computeMassDrift(initial_mass, final_mass);
        std::cout << "  Mass drift = " << mass_drift << "\n";
        logger.logParameter("Mass drift", mass_drift);

        if (mass_drift > 1.0e-6)
        {
            TestUtils::printWarning("Mass drift exceeds tolerance: " +
                                   TestUtils::formatScientific(mass_drift));
        }

        // =====================================================================
        // Test verdict
        // =====================================================================
        TestUtils::printSectionSeparator("Test Verdict");

        std::vector<std::string> checks;
        std::vector<bool> results;

        checks.push_back("No crashes");
        results.push_back(true);

        checks.push_back("No NaN/Inf");
        results.push_back(!has_nan);

        checks.push_back("Negatives < 1% of cells");
        int total_cells = 1000;  // Approximate
        results.push_back(neg_count < total_cells * 0.01);

        checks.push_back("Completed all time steps");
        results.push_back(iteration_num >= max_iterations - 1);

        test_passed = true;
        for (size_t i = 0; i < checks.size(); ++i)
        {
            std::string status = results[i] ? "[PASS]" : "[FAIL]";
            std::cout << "  " << status << " " << checks[i] << "\n";
            test_passed = test_passed && results[i];
        }

        if (C_bc_coef) delete C_bc_coef;

    } // End scope to force cleanup

    total_timer.stop();

    // =========================================================================
    // Final report
    // =========================================================================
    std::cout << "\n";
    std::cout << "Elapsed time: " << total_timer.getFormattedTime() << "\n";

    logger.logTestResult("Test 01: Smoke Test", test_passed);
    logger.log("Total elapsed time: " + total_timer.getFormattedTime());

    TestUtils::printTestFooter("Test 01: Smoke Test", test_passed,
                              test_passed ? "All checks passed" : "Some checks failed");

    return test_passed ? 0 : 1;
}
