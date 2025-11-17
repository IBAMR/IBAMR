// =============================================================================
// TEST 03: Pure Advection - Analytical Validation
// =============================================================================
// Validates advection operator by comparing to analytical moving Gaussian
// C(x,y,t) = C0 * exp(-((x - x0 - u*t)^2 + (y - y0 - v*t)^2) / (2*sigma^2))
// Tests: Conservative advection, shape preservation, no numerical diffusion
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
    TestUtils::printTestHeader("Pure Advection - Analytical Validation", 3);

    Timer total_timer;
    total_timer.start();

    bool test_passed = true;
    ResultLogger logger("test03_results.txt");

    {
        // =====================================================================
        // Initialize application
        // =====================================================================
        TestUtils::printProgress("Initializing application...");
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test03.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get solver parameters
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");
        const double end_time = main_db->getDouble("END_TIME");
        const double dt = main_db->getDouble("DT");

        // Get physical parameters
        const double kappa = input_db->getDoubleWithDefault("diffusion_coefficient", 1e-6);  // Minimal diffusion
        const double u = input_db->getDoubleWithDefault("advection_velocity_u", 0.5);
        const double v = input_db->getDoubleWithDefault("advection_velocity_v", 0.3);

        // Gaussian parameters
        const double x0 = input_db->getDoubleWithDefault("gaussian_x0", -0.3);
        const double y0 = input_db->getDoubleWithDefault("gaussian_y0", -0.2);
        const double sigma = input_db->getDoubleWithDefault("gaussian_sigma", 0.15);
        const double C0 = input_db->getDoubleWithDefault("gaussian_amplitude", 1.0);

        logger.logParameter("Diffusion coefficient", kappa);
        logger.logParameter("Advection velocity u", u);
        logger.logParameter("Advection velocity v", v);
        logger.logParameter("Gaussian center x0", x0);
        logger.logParameter("Gaussian center y0", y0);
        logger.logParameter("Gaussian sigma", sigma);
        logger.logParameter("End time", end_time);
        logger.logParameter("Time step", dt);

        // Compute Peclet number (should be high for pure advection)
        double L = 2.0;  // Domain size
        double Pe = TestUtils::computePecletNumber(sqrt(u*u + v*v), L, kappa);
        logger.logParameter("Peclet number", Pe);

        std::cout << "\n  Pure Advection Test Configuration:\n";
        std::cout << "    Velocity: u = " << u << ", v = " << v << "\n";
        std::cout << "    Gaussian: center = (" << x0 << ", " << y0 << "), sigma = " << sigma << "\n";
        std::cout << "    Peclet number: " << Pe << " (high Pe => advection-dominated)\n\n";

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
        // Set initial conditions (Gaussian profile)
        // =====================================================================
        TestUtils::printProgress("Setting initial Gaussian profile...");

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
        // Set diffusion coefficient (minimal for pure advection)
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

        // Get initial mass for conservation check
        double initial_mass = ErrorCalculator::computeTotalMass(patch_hierarchy, C_idx);
        logger.logParameter("Initial mass", initial_mass);

        // =====================================================================
        // Main time integration loop
        // =====================================================================
        TestUtils::printSectionSeparator("Time Integration - Pure Advection");

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
        }

        // =====================================================================
        // Compute exact analytical solution at final time
        // =====================================================================
        TestUtils::printSectionSeparator("Error Analysis - Advected Gaussian");
        TestUtils::printProgress("Computing exact analytical solution...");

        Pointer<CellVariable<NDIM, double>> C_exact_var = new CellVariable<NDIM, double>("C_exact");
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int C_exact_idx = var_db->registerVariableAndContext(
            C_exact_var, var_db->getContext("EXACT"), IntVector<NDIM>(1));

        for (int ln = 0; ln <= finest_level; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(C_exact_idx, loop_time);
        }

        // Set exact advected Gaussian solution
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
                const double* dx = pgeom->getDx();

                for (CellIterator<NDIM> ci(patch_box); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();

                    double X[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        X[d] = x_lo[d] + dx[d] * (static_cast<double>(idx(d) - patch_box.lower()(d)) + 0.5);
                    }

                    // Analytical advected Gaussian
                    (*C_exact_data)(idx) = AnalyticalSolutions::advectedGaussian2D(
                        X[0], X[1], loop_time, u, v, x0, y0, sigma, C0);
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

        ErrorCalculator::printErrorSummary("Test03: Pure Advection", L1_error, L2_error, Linf_error);

        logger.logError("L1 error", L1_error);
        logger.logError("L2 error", L2_error);
        logger.logError("Linf error", Linf_error);

        // =====================================================================
        // Mass conservation check
        // =====================================================================
        double final_mass = ErrorCalculator::computeTotalMass(patch_hierarchy, C_idx);
        double mass_drift = TestUtils::computeMassDrift(initial_mass, final_mass);

        std::cout << "\n  Mass Conservation:\n";
        std::cout << "    Initial mass: " << TestUtils::formatScientific(initial_mass) << "\n";
        std::cout << "    Final mass:   " << TestUtils::formatScientific(final_mass) << "\n";
        std::cout << "    Mass drift:   " << TestUtils::formatScientific(mass_drift) << "\n\n";

        logger.logParameter("Final mass", final_mass);
        logger.logParameter("Mass drift", mass_drift);

        // =====================================================================
        // Shape preservation check (check for numerical diffusion)
        // =====================================================================
        double C_min, C_max;
        ErrorCalculator::getMinMax(patch_hierarchy, C_idx, C_min, C_max);

        std::cout << "  Shape Preservation:\n";
        std::cout << "    Expected max: " << C0 << "\n";
        std::cout << "    Computed max: " << C_max << "\n";
        std::cout << "    Difference:   " << fabs(C_max - C0) << " (numerical diffusion indicator)\n\n";

        logger.logParameter("Expected max C", C0);
        logger.logParameter("Computed max C", C_max);
        logger.logParameter("Numerical diffusion", fabs(C_max - C0));

        // =====================================================================
        // Test verdict
        // =====================================================================
        TestUtils::printSectionSeparator("Test Verdict");

        std::vector<std::string> checks;
        std::vector<bool> results;

        // Check 1: No NaN
        bool has_nan = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);
        checks.push_back("No NaN/Inf");
        results.push_back(!has_nan);

        // Check 2: L2 error (pure advection should have low error)
        checks.push_back("L2 error < 0.05");
        results.push_back(L2_error < 0.05);

        // Check 3: Linf error
        checks.push_back("Linf error < 0.2");
        results.push_back(Linf_error < 0.2);

        // Check 4: Mass conservation (should be near-perfect)
        checks.push_back("Mass drift < 1e-5");
        results.push_back(mass_drift < 1.0e-5);

        // Check 5: Shape preservation (minimal numerical diffusion)
        checks.push_back("Numerical diffusion < 0.05");
        results.push_back(fabs(C_max - C0) < 0.05);

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

    logger.logTestResult("Test 03: Pure Advection", test_passed);
    logger.log("Total elapsed time: " + total_timer.getFormattedTime());

    TestUtils::printTestFooter("Test 03: Pure Advection", test_passed,
                              test_passed ? "Advection operator validated" : "Some checks failed");

    return test_passed ? 0 : 1;
}
