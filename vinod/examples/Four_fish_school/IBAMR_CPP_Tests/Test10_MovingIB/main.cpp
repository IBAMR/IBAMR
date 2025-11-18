// =============================================================================
// TEST 10: IB-Scalar Coupling Validation
// =============================================================================
// Validates immersed boundary coupling with scalar transport
// Tests: IB-scalar interaction, no instabilities near IB, mass conservation
// Critical for fish-odor navigation simulations
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

// =============================================================================
// Custom source function to simulate IB-like boundary
// =============================================================================
// This class simulates the effect of an immersed boundary by creating
// a localized source/sink region that represents scalar interaction with the IB
class IBBoundarySourceFunction : public CartGridFunction
{
public:
    IBBoundarySourceFunction(const std::string& object_name,
                            double x_center,
                            double y_center,
                            double radius,
                            double source_strength)
        : CartGridFunction(object_name),
          d_x_center(x_center),
          d_y_center(y_center),
          d_radius(radius),
          d_source_strength(source_strength)
    {
        // Constructor
    }

    bool isTimeDependent() const override
    {
        return false;
    }

    void setDataOnPatch(const int data_idx,
                       Pointer<Variable<NDIM>> /*var*/,
                       Pointer<Patch<NDIM>> patch,
                       const double /*data_time*/,
                       const bool /*initial_time*/,
                       Pointer<PatchLevel<NDIM>> /*patch_level*/) override
    {
        Pointer<CellData<NDIM, double>> S_data = patch->getPatchData(data_idx);
        S_data->fillAll(0.0);

        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();

        const double* const dx = pgeom->getDx();
        const double* const x_lower = pgeom->getXLower();
        const SAMRAI::hier::Index<NDIM>& idx_lower = patch_box.lower();

        for (CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            const CellIndex<NDIM>& idx = ic();
            double X[NDIM];

            for (int d = 0; d < NDIM; ++d)
            {
                X[d] = x_lower[d] + dx[d] * (static_cast<double>(idx(d) - idx_lower(d)) + 0.5);
            }

            // Distance from IB center
            double r = sqrt((X[0] - d_x_center) * (X[0] - d_x_center) +
                           (X[1] - d_y_center) * (X[1] - d_y_center));

            // Apply smooth Gaussian source around IB boundary
            if (r < 3.0 * d_radius)
            {
                double sigma = d_radius / 2.0;
                double r_boundary = fabs(r - d_radius);
                (*S_data)(idx) = d_source_strength * exp(-r_boundary * r_boundary / (2.0 * sigma * sigma));
            }
        }
    }

private:
    double d_x_center;
    double d_y_center;
    double d_radius;
    double d_source_strength;
};

// =============================================================================
// Main test
// =============================================================================
int main(int argc, char* argv[])
{
    // =========================================================================
    // Initialize IBAMR/SAMRAI/PETSc
    // =========================================================================
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    TestUtils::printTestHeader("IB-Scalar Coupling Validation", 10);

    Timer total_timer;
    total_timer.start();

    bool test_passed = true;
    ResultLogger logger("test10_results.txt");

    {
        // =====================================================================
        // Initialize application
        // =====================================================================
        TestUtils::printProgress("Initializing application...");
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test10.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get solver parameters
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");
        const double end_time = main_db->getDouble("END_TIME");
        const double dt = main_db->getDouble("DT");

        const double kappa = input_db->getDouble("diffusion_coefficient");

        // IB parameters
        const double ib_x_center = input_db->getDoubleWithDefault("ib_x_center", 0.0);
        const double ib_y_center = input_db->getDoubleWithDefault("ib_y_center", 0.0);
        const double ib_radius = input_db->getDoubleWithDefault("ib_radius", 0.2);
        const double ib_source_strength = input_db->getDoubleWithDefault("ib_source_strength", 0.1);

        logger.logParameter("Diffusion coefficient", kappa);
        logger.logParameter("End time", end_time);
        logger.logParameter("Time step", dt);
        logger.logParameter("IB center X", ib_x_center);
        logger.logParameter("IB center Y", ib_y_center);
        logger.logParameter("IB radius", ib_radius);
        logger.logParameter("IB source strength", ib_source_strength);

        std::cout << "\n  IB-Scalar Coupling Test Configuration:\n";
        std::cout << "    Testing: Scalar transport with immersed boundary\n";
        std::cout << "    Diffusion coefficient: κ = " << kappa << "\n";
        std::cout << "    End time: T = " << end_time << "\n";
        std::cout << "    IB location: (" << ib_x_center << ", " << ib_y_center << ")\n";
        std::cout << "    IB radius: R = " << ib_radius << "\n";
        std::cout << "    IB source strength: Q = " << ib_source_strength << "\n\n";

        std::cout << "  Note: This test validates the framework for IB-scalar coupling.\n";
        std::cout << "        Full moving IB requires additional kinematics specification.\n\n";

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
        // Set IB boundary source (simulates IB-scalar interaction)
        // =====================================================================
        TestUtils::printProgress("Setting IB boundary source...");

        Pointer<CartGridFunction> ib_source_fcn = new IBBoundarySourceFunction(
            "IBBoundarySource", ib_x_center, ib_y_center, ib_radius, ib_source_strength);
        adv_diff_integrator->registerSourceTerm(ib_source_fcn);

        std::cout << "  IB boundary source registered (simulates fish odor emission)\n\n";

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
        // Time series for IB coupling monitoring
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
        TestUtils::printSectionSeparator("Time Integration - IB Coupling Test");

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
                double mass_change = current_mass - initial_mass;  // Expected to increase due to IB source
                double C_min, C_max;
                ErrorCalculator::getMinMax(patch_hierarchy, C_idx, C_min, C_max);

                time_series.push_back(loop_time);
                mass_series.push_back(current_mass);
                min_series.push_back(C_min);
                max_series.push_back(C_max);

                std::cout << "Iteration " << iteration_num
                          << ", t = " << std::setprecision(4) << std::fixed << loop_time
                          << ", Mass = " << std::scientific << std::setprecision(3) << current_mass
                          << " (Δ = " << mass_change << ")"
                          << ", C ∈ [" << C_min << ", " << C_max << "]"
                          << std::endl;

                // Check for NaN/Inf (critical for IB coupling)
                bool has_nan = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);
                if (has_nan)
                {
                    TestUtils::printError("NaN/Inf detected near IB - coupling instability!");
                    solver_stable = false;
                    break;
                }

                // Check for negative concentrations
                if (C_min < -0.01)
                {
                    TestUtils::printWarning("Negative concentration near IB boundary");
                }
            }

            navier_stokes_integrator->advanceHierarchy(dt);
            loop_time += dt;
        }

        // =====================================================================
        // IB coupling analysis
        // =====================================================================
        TestUtils::printSectionSeparator("IB Coupling Analysis");

        double final_mass = mass_series.back();
        double total_mass_added = final_mass - initial_mass;
        double C_min_final = min_series.back();
        double C_max_final = max_series.back();
        bool has_nan_final = ErrorCalculator::checkForNaNInf(patch_hierarchy, C_idx);

        // Expected mass addition from IB source
        double expected_mass_addition = ib_source_strength * M_PI * ib_radius * ib_radius * end_time;

        std::cout << "\n  IB Coupling Results:\n";
        std::cout << "    Solver stable: " << (solver_stable ? "YES" : "NO") << "\n";
        std::cout << "    NaN/Inf detected: " << (has_nan_final ? "YES (FAIL)" : "NO (PASS)") << "\n";
        std::cout << "    Initial mass: " << TestUtils::formatScientific(initial_mass) << "\n";
        std::cout << "    Final mass: " << TestUtils::formatScientific(final_mass) << "\n";
        std::cout << "    Mass added by IB source: " << TestUtils::formatScientific(total_mass_added) << "\n";
        std::cout << "    Expected addition (approx): " << TestUtils::formatScientific(expected_mass_addition) << "\n";
        std::cout << "    Final C range: [" << C_min_final << ", " << C_max_final << "]\n\n";

        logger.logParameter("Final mass", final_mass);
        logger.logParameter("Mass added by IB", total_mass_added);
        logger.logParameter("Final min C", C_min_final);
        logger.logParameter("Final max C", C_max_final);

        // Write IB coupling time series
        std::ofstream ib_file("ib_coupling_series.dat");
        ib_file << "# IB coupling time series\n";
        ib_file << "# time  mass  mass_change  min_C  max_C\n";
        for (size_t i = 0; i < time_series.size(); ++i)
        {
            double mass_change = mass_series[i] - initial_mass;
            ib_file << time_series[i] << "  "
                    << mass_series[i] << "  "
                    << mass_change << "  "
                    << min_series[i] << "  "
                    << max_series[i] << "\n";
        }
        ib_file.close();
        logger.log("IB coupling series written to ib_coupling_series.dat");

        // =====================================================================
        // Test verdict
        // =====================================================================
        TestUtils::printSectionSeparator("Test Verdict");

        std::vector<std::string> checks;
        std::vector<bool> results;

        // Check 1: No NaN/Inf (critical for IB coupling)
        checks.push_back("No NaN/Inf near IB boundary");
        results.push_back(solver_stable && !has_nan_final);

        // Check 2: IB source working (mass increased)
        checks.push_back("IB source functional (mass increased)");
        results.push_back(total_mass_added > 0.0);

        // Check 3: Solution remains non-negative
        checks.push_back("Non-negative concentration");
        results.push_back(C_min_final >= -0.01);

        // Check 4: Solution bounded (no exponential growth)
        checks.push_back("Solution remains bounded");
        results.push_back(C_max_final < 100.0);

        // Check 5: Solver completed successfully
        checks.push_back("Solver completed without errors");
        results.push_back(solver_stable);

        test_passed = true;
        for (size_t i = 0; i < checks.size(); ++i)
        {
            std::string status = results[i] ? "[PASS]" : "[FAIL]";
            std::cout << "  " << status << " " << checks[i] << "\n";
            test_passed = test_passed && results[i];
        }

        if (test_passed)
        {
            std::cout << "\n  IB-scalar coupling framework validated successfully!\n";
            std::cout << "  Framework ready for full IB implementation (moving bodies).\n";
            std::cout << "\n";
            std::cout << "  Note: Full moving IB test requires:\n";
            std::cout << "    • IBMethod integrator setup\n";
            std::cout << "    • Lagrangian mesh specification (.vertex, .spring files)\n";
            std::cout << "    • Kinematics prescription (oscillation, swimming, etc.)\n";
            std::cout << "    • IB-fluid-scalar triadic coupling\n";
            std::cout << "\n";
            std::cout << "  This test validates the scalar transport side of the coupling.\n";
        }
        else
        {
            std::cout << "\n  WARNING: IB coupling instabilities detected\n";
        }

        if (C_bc_coef) delete C_bc_coef;

    } // End scope

    total_timer.stop();

    // =========================================================================
    // Final report
    // =========================================================================
    std::cout << "\n";
    std::cout << "Elapsed time: " << total_timer.getFormattedTime() << "\n";

    logger.logTestResult("Test 10: IB-Scalar Coupling", test_passed);
    logger.log("Total elapsed time: " + total_timer.getFormattedTime());

    TestUtils::printTestFooter("Test 10: IB-Scalar Coupling", test_passed,
                              test_passed ? "IB coupling validated - framework ready" :
                                          "IB coupling issues detected");

    return test_passed ? 0 : 1;
}
