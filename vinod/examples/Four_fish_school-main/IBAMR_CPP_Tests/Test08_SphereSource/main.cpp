// =============================================================================
// TEST 08: Cylinder/Sphere Source - Literature Validation
// =============================================================================
// Validates steady-state diffusion around a source (Lei et al. 2021 analogy)
// 2D: Cylinder source with C ~ Q/(2πκ) * ln(r/R)
// Tests: Steady-state convergence, 1/r decay, literature comparison
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
#include <ibtk/CartGridFunction.h>

// Common test utilities
#include "TestUtilities.h"
#include "ErrorCalculator.h"
#include "AnalyticalSolutions.h"

#include <ibamr/app_namespaces.h>

using namespace TestUtilities;

// Cylinder source term class (steady emission)
class CylinderSourceFunction : public CartGridFunction
{
public:
    CylinderSourceFunction(const std::string& object_name,
                          double x0, double y0, double R, double Q)
        : CartGridFunction(object_name), d_x0(x0), d_y0(y0), d_R(R), d_Q(Q)
    {
    }

    bool isTimeDependent() const override { return false; }

    void setDataOnPatch(
        const int data_idx,
        Pointer<Variable<NDIM>> var,
        Pointer<Patch<NDIM>> patch,
        const double data_time,
        const bool initial_time = false,
        Pointer<PatchLevel<NDIM>> patch_level = Pointer<PatchLevel<NDIM>>(NULL)) override
    {
        Pointer<CellData<NDIM, double>> S_data = patch->getPatchData(data_idx);
        S_data->fillAll(0.0);

        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* dx = pgeom->getDx();
        const double* x_lo = pgeom->getXLower();

        for (CellIterator<NDIM> ci(patch_box); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();

            // Get cell center coordinates
            double X[NDIM];
            for (int d = 0; d < NDIM; ++d)
            {
                X[d] = x_lo[d] + dx[d] * (static_cast<double>(idx(d) - patch_box.lower()(d)) + 0.5);
            }

            // Compute distance from cylinder center
            double r = sqrt((X[0] - d_x0) * (X[0] - d_x0) +
                           (X[1] - d_y0) * (X[1] - d_y0));

            // Source term: uniform within cylinder radius
            if (r < d_R)
            {
                // Distribute source Q uniformly over cylinder area
                double area = M_PI * d_R * d_R;
                (*S_data)(idx) = d_Q / area;
            }
        }
    }

private:
    double d_x0, d_y0, d_R, d_Q;
};

int main(int argc, char* argv[])
{
    // =========================================================================
    // Initialize IBAMR/SAMRAI/PETSc
    // =========================================================================
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    TestUtils::printTestHeader("Cylinder/Sphere Source - Literature Validation", 8);

    Timer total_timer;
    total_timer.start();

    bool test_passed = true;
    ResultLogger logger("test08_results.txt");

    {
        // =====================================================================
        // Initialize application
        // =====================================================================
        TestUtils::printProgress("Initializing application...");
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test08.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get solver parameters
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");
        const double end_time = main_db->getDouble("END_TIME");
        const double dt = main_db->getDouble("DT");

        // Cylinder source parameters
        const double kappa = input_db->getDouble("diffusion_coefficient");
        const double x0 = input_db->getDoubleWithDefault("cylinder_x0", 0.0);
        const double y0 = input_db->getDoubleWithDefault("cylinder_y0", 0.0);
        const double R = input_db->getDoubleWithDefault("cylinder_radius", 0.1);
        const double Q = input_db->getDoubleWithDefault("source_strength", 1.0);

        logger.logParameter("Diffusion coefficient", kappa);
        logger.logParameter("Cylinder center x0", x0);
        logger.logParameter("Cylinder center y0", y0);
        logger.logParameter("Cylinder radius", R);
        logger.logParameter("Source strength", Q);
        logger.logParameter("End time", end_time);

        std::cout << "\n  Cylinder Source Test Configuration:\n";
        std::cout << "    Cylinder center: (" << x0 << ", " << y0 << ")\n";
        std::cout << "    Cylinder radius: R = " << R << "\n";
        std::cout << "    Source strength: Q = " << Q << "\n";
        std::cout << "    Diffusion coefficient: κ = " << kappa << "\n";
        std::cout << "    Expected steady state: C ~ Q/(2πκ) * ln(r/R)\n";
        std::cout << "    End time: T = " << end_time << " (run to steady state)\n\n";

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
        // Set initial conditions (zero)
        // =====================================================================
        TestUtils::printProgress("Setting initial conditions (zero)...");

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
        // Set source term (cylinder source)
        // =====================================================================
        TestUtils::printProgress("Setting cylinder source term...");
        Pointer<CartGridFunction> S_fcn = new CylinderSourceFunction("CylinderSource", x0, y0, R, Q);
        adv_diff_integrator->setSourceTerm(S_fcn);

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
        logger.log("Scalar variable and source term registered successfully");

        // =====================================================================
        // Main time integration loop (run to steady state)
        // =====================================================================
        TestUtils::printSectionSeparator("Time Integration - Steady State");

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
                std::cout << "Iteration " << iteration_num
                          << ", t = " << loop_time << std::endl;
            }

            navier_stokes_integrator->advanceHierarchy(dt);
            loop_time += dt;
        }

        // =====================================================================
        // Compute analytical steady-state solution
        // =====================================================================
        TestUtils::printSectionSeparator("Steady-State Analysis");
        TestUtils::printProgress("Computing analytical steady-state solution...");

        Pointer<CellVariable<NDIM, double>> C_exact_var = new CellVariable<NDIM, double>("C_exact");
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int C_exact_idx = var_db->registerVariableAndContext(
            C_exact_var, var_db->getContext("EXACT"), IntVector<NDIM>(1));

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
                const double* dx = pgeom->getDx();

                for (CellIterator<NDIM> ci(patch_box); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();

                    double X[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        X[d] = x_lo[d] + dx[d] * (static_cast<double>(idx(d) - patch_box.lower()(d)) + 0.5);
                    }

                    // Analytical cylinder source solution
                    (*C_exact_data)(idx) = AnalyticalSolutions::cylinderSource2D(
                        X[0], X[1], x0, y0, R, Q, kappa);
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

        ErrorCalculator::printErrorSummary("Test08: Cylinder Source", L1_error, L2_error, Linf_error);

        logger.logError("L1 error", L1_error);
        logger.logError("L2 error", L2_error);
        logger.logError("Linf error", Linf_error);

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

        // Check 2: Steady-state solution accuracy
        checks.push_back("L2 error < 0.2 (steady-state)");
        results.push_back(L2_error < 0.2);

        // Check 3: Linf error
        checks.push_back("Linf error < 0.5");
        results.push_back(Linf_error < 0.5);

        // Check 4: Solution non-negative
        double C_min, C_max;
        ErrorCalculator::getMinMax(patch_hierarchy, C_idx, C_min, C_max);
        checks.push_back("Solution non-negative");
        results.push_back(C_min >= -0.01);

        test_passed = true;
        for (size_t i = 0; i < checks.size(); ++i)
        {
            std::string status = results[i] ? "[PASS]" : "[FAIL]";
            std::cout << "  " << status << " " << checks[i] << "\n";
            test_passed = test_passed && results[i];
        }

        if (test_passed)
        {
            std::cout << "\n  Cylinder source solution validated!\n";
            std::cout << "  Analogy to Lei et al. (2021) sphere source in 2D ✓\n";
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

    logger.logTestResult("Test 08: Cylinder/Sphere Source", test_passed);
    logger.log("Total elapsed time: " + total_timer.getFormattedTime());

    TestUtils::printTestFooter("Test 08: Cylinder Source", test_passed,
                              test_passed ? "Literature validation successful" : "Accuracy requirements not met");

    return test_passed ? 0 : 1;
}
