// =============================================================================
// TEST 02: Pure Diffusion Analytic
// =============================================================================
// Validates diffusion operator against analytical Gaussian solution
// C(x,y,t) = (C0/(4*pi*kappa*t)) * exp(-r^2/(4*kappa*t))
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

using namespace TestUtilities;

int main(int argc, char* argv[])
{
    // =========================================================================
    // Initialize IBAMR/SAMRAI/PETSc
    // =========================================================================
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    TestUtils::printTestHeader("Pure Diffusion Analytic - Gaussian Solution", 2);

    Timer total_timer;
    total_timer.start();

    bool test_passed = true;
    ResultLogger logger("test02_results.txt");

    {
        // =====================================================================
        // Initialize application
        // =====================================================================
        TestUtils::printProgress("Initializing application...");
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test02.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get visualization parameters
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();

        // Get solver parameters
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");
        const double end_time = main_db->getDouble("END_TIME");
        const double dt = main_db->getDouble("DT");

        // Get diffusion coefficient
        const double kappa = input_db->getDouble("diffusion_coefficient");
        logger.logParameter("Diffusion coefficient (kappa)", kappa);
        logger.logParameter("End time", end_time);
        logger.logParameter("Time step", dt);

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
        TestUtils::printProgress("Creating grid geometry...");
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
        // Set initial conditions (Gaussian at t=0)
        // =====================================================================
        TestUtils::printProgress("Setting initial conditions...");

        if (input_db->keyExists("OdorInitialConditions"))
        {
            Pointer<CartGridFunction> C_init = new muParserCartGridFunction(
                "C_init", app_initializer->getComponentDatabase("OdorInitialConditions"), grid_geometry);
            adv_diff_integrator->setInitialConditions(C_init);
        }

        // =====================================================================
        // Set boundary conditions (homogeneous Dirichlet)
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
        // Main time integration loop
        // =====================================================================
        TestUtils::printSectionSeparator("Time Integration - Pure Diffusion");

        double loop_time = navier_stokes_integrator->getIntegratorTime();
        double loop_time_end = navier_stokes_integrator->getEndTime();
        int iteration_num = navier_stokes_integrator->getIntegratorStep();

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

            // Check for NaN/Inf
            if (iteration_num % 10 == 0)
            {
                bool has_nan = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);
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
        // Compute analytical solution at final time
        // =====================================================================
        TestUtils::printSectionSeparator("Error Analysis");
        TestUtils::printProgress("Computing analytical solution...");

        // Create variable for exact solution
        Pointer<CellVariable<NDIM, double>> C_exact_var = new CellVariable<NDIM, double>("C_exact");
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int C_exact_idx = var_db->registerVariableAndContext(
            C_exact_var, var_db->getContext("EXACT"), IntVector<NDIM>(1));

        // Allocate exact solution variable
        for (int ln = 0; ln <= finest_level; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(C_exact_idx, loop_time);
        }

        // Set exact solution
        const double* x_lo = grid_geometry->getXLower();
        const double* x_up = grid_geometry->getXUpper();
        const double x0 = 0.0;  // Center of Gaussian
        const double y0 = 0.0;
        const double C0 = initial_mass;  // Total mass preserved

        for (int ln = 0; ln <= finest_level; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();

                Pointer<CellData<NDIM, double>> C_exact_data = patch->getPatchData(C_exact_idx);
                const double* dx = pgeom->getDx();

                for (CellIterator<NDIM> ci(patch_box); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();

                    // Get cell center coordinates
                    double X[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        X[d] = x_lo[d] + dx[d] * (static_cast<double>(idx(d) - patch_box.lower()(d)) + 0.5);
                    }

                    // Compute analytical solution
                    (*C_exact_data)(idx) = AnalyticalSolutions::gaussianDiffusion2D(
                        X[0], X[1], loop_time, kappa, x0, y0, C0);
                }
            }
        }

        // =====================================================================
        // Compute errors
        // =====================================================================
        TestUtils::printProgress("Computing errors...");

        double L1_error = ErrorCalculator::computeL1Error(patch_hierarchy, C_idx, C_exact_idx);
        double L2_error = ErrorCalculator::computeL2Error(patch_hierarchy, C_idx, C_exact_idx, loop_time);
        double Linf_error = ErrorCalculator::computeLinfError(patch_hierarchy, C_idx, C_exact_idx);

        ErrorCalculator::printErrorSummary("Test02: Pure Diffusion", L1_error, L2_error, Linf_error);

        logger.logError("L1 error", L1_error);
        logger.logError("L2 error", L2_error);
        logger.logError("Linf error", Linf_error);

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

        // Check mass conservation
        double mass_drift = TestUtils::computeMassDrift(initial_mass, final_mass);
        std::cout << "  Mass drift = " << mass_drift << "\n";
        logger.logParameter("Mass drift", mass_drift);

        if (mass_drift > 1.0e-6)
        {
            TestUtils::printWarning("Mass drift exceeds tolerance: " +
                                   TestUtils::formatScientific(mass_drift));
        }

        // Check for negatives
        int neg_count = ErrorCalculator::checkForNegatives(patch_hierarchy, C_idx);
        if (neg_count > 0)
        {
            TestUtils::printWarning("Found " + std::to_string(neg_count) + " negative values");
        }

        // =====================================================================
        // Test verdict
        // =====================================================================
        TestUtils::printSectionSeparator("Test Verdict");

        std::vector<std::string> checks;
        std::vector<bool> results;

        checks.push_back("No NaN/Inf");
        bool has_nan = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);
        results.push_back(!has_nan);

        checks.push_back("L2 error < 0.01");
        results.push_back(L2_error < 0.01);

        checks.push_back("Mass drift < 1e-6");
        results.push_back(mass_drift < 1.0e-6);

        checks.push_back("No negative values");
        results.push_back(neg_count == 0);

        test_passed = true;
        for (size_t i = 0; i < checks.size(); ++i)
        {
            std::string status = results[i] ? "[PASS]" : "[FAIL]";
            std::cout << "  " << status << " " << checks[i] << "\n";
            test_passed = test_passed && results[i];
        }

        // Cleanup
        for (int ln = 0; ln <= finest_level; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(C_exact_idx);
        }

        if (C_bc_coef) delete C_bc_coef;

    } // End scope

    total_timer.stop();

    // =========================================================================
    // Final report
    // =========================================================================
    std::cout << "\n";
    std::cout << "Elapsed time: " << total_timer.getFormattedTime() << "\n";

    logger.logTestResult("Test 02: Pure Diffusion Analytic", test_passed);
    logger.log("Total elapsed time: " + total_timer.getFormattedTime());

    TestUtils::printTestFooter("Test 02: Pure Diffusion Analytic", test_passed,
                              test_passed ? "Diffusion operator validated" : "Some checks failed");

    return test_passed ? 0 : 1;
}
