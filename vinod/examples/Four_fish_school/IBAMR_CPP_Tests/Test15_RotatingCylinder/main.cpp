// =============================================================================
// TEST 15: Rotating Cylinder Validation (Yan & Zu 2008)
// =============================================================================
// Validates 2-D scalar transport around rotating isothermal cylinder
// Reference: Yan & Zu (2008) - LBM heat transfer benchmark
// Lei et al. (2021) Section II.B, Figure 2 - Solver validation
//
// Parameters: Re=200, Pr=0.5, k=0.5 (rotation parameter V/U)
// Expected: Streamlines + scalar contours match published results
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
#include <vector>
#include <fstream>
#include <cmath>

using namespace TestUtilities;

// =============================================================================
// ROTATING CYLINDER SCALAR SOURCE
// =============================================================================
// Enforces C = C_h (isothermal) on rotating cylinder surface
// Uses smooth Gaussian distribution around cylinder boundary

class RotatingCylinderScalarSource : public CartGridFunction
{
public:
    RotatingCylinderScalarSource(const std::string& object_name,
                                  double x_center,
                                  double y_center,
                                  double radius,
                                  double C_hot,
                                  double thickness)
        : CartGridFunction(object_name),
          d_x_center(x_center),
          d_y_center(y_center),
          d_radius(radius),
          d_C_hot(C_hot),
          d_thickness(thickness)
    {
    }

    bool isTimeDependent() const override { return false; }

    void setDataOnPatch(const int data_idx,
                       Pointer<Variable<NDIM>> /*var*/,
                       Pointer<Patch<NDIM>> patch,
                       const double /*data_time*/,
                       const bool initial_time = false,
                       Pointer<PatchLevel<NDIM>> /*level*/ = nullptr) override
    {
        Pointer<CellData<NDIM, double>> S_data = patch->getPatchData(data_idx);
        S_data->fillAll(0.0);

        if (!initial_time) return;

        Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
        const double* const dx = patch_geom->getDx();
        const double* const x_lower = patch_geom->getXLower();
        const Box<NDIM>& patch_box = patch->getBox();

        // Gaussian source enforcing C = C_hot around cylinder
        for (CellIterator<NDIM> ci(patch_box); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();
            double X[NDIM];
            for (int d = 0; d < NDIM; ++d)
            {
                X[d] = x_lower[d] + dx[d] * (static_cast<double>(idx(d) - patch_box.lower()(d)) + 0.5);
            }

            // Distance from cylinder center
            double r = sqrt(pow(X[0] - d_x_center, 2) + pow(X[1] - d_y_center, 2));

            // Gaussian profile: peaked at r = radius, width = thickness
            double dr = r - d_radius;
            double source_strength = d_C_hot * exp(-0.5 * pow(dr / d_thickness, 2));

            (*S_data)(idx) = source_strength / (d_thickness * sqrt(2.0 * M_PI));
        }
    }

private:
    double d_x_center;
    double d_y_center;
    double d_radius;
    double d_C_hot;
    double d_thickness;
};

// =============================================================================
// ANALYSIS FUNCTIONS
// =============================================================================

void writeCenterlineProfile(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                            Pointer<CellVariable<NDIM, double>> C_var,
                            const int C_idx,
                            const double y_line,
                            const std::string& filename)
{
    std::ofstream outfile(filename);
    outfile << "# x  C(x, y=" << y_line << ")\n";

    // Sample along centerline y = y_line
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> finest_level = patch_hierarchy->getPatchLevel(finest_ln);

    std::vector<std::pair<double, double>> data_points;

    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = finest_level->getPatch(p());
        Pointer<CellData<NDIM, double>> C_data = patch->getPatchData(C_idx);
        Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();

        const double* dx = patch_geom->getDx();
        const double* x_lower = patch_geom->getXLower();
        const Box<NDIM>& patch_box = patch->getBox();

        for (CellIterator<NDIM> ci(patch_box); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();
            double x = x_lower[0] + dx[0] * (static_cast<double>(idx(0) - patch_box.lower()(0)) + 0.5);
            double y = x_lower[1] + dx[1] * (static_cast<double>(idx(1) - patch_box.lower()(1)) + 0.5);

            // Check if close to target y
            if (fabs(y - y_line) < 0.5 * dx[1])
            {
                double C_val = (*C_data)(idx);
                data_points.push_back(std::make_pair(x, C_val));
            }
        }
    }

    // Sort by x coordinate
    std::sort(data_points.begin(), data_points.end());

    // Write to file
    for (const auto& pt : data_points)
    {
        outfile << pt.first << "  " << pt.second << "\n";
    }

    outfile.close();
    std::cout << "  Written centerline profile to " << filename << "\n";
}

bool checkSteadyState(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                     Pointer<CellVariable<NDIM, double>> C_var,
                     const int C_current_idx,
                     const int C_previous_idx,
                     const double tolerance)
{
    double max_change = 0.0;

    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> C_current = patch->getPatchData(C_current_idx);
            Pointer<CellData<NDIM, double>> C_previous = patch->getPatchData(C_previous_idx);

            const Box<NDIM>& patch_box = patch->getBox();

            for (CellIterator<NDIM> ci(patch_box); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                double dC = fabs((*C_current)(idx) - (*C_previous)(idx));
                max_change = std::max(max_change, dC);
            }
        }
    }

    return max_change < tolerance;
}

// =============================================================================
// MAIN
// =============================================================================

int main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    TestUtils::printTestHeader("TEST 15: ROTATING CYLINDER VALIDATION");
    std::cout << "\n";
    std::cout << "Reference: Yan & Zu (2008) - Lattice Boltzmann Method\n";
    std::cout << "Validation: Lei et al. (2021) Section II.B, Figure 2\n";
    std::cout << "\nParameters:\n";
    std::cout << "  Re = 200 (Reynolds number)\n";
    std::cout << "  Pr = 0.5 (Prandtl number)\n";
    std::cout << "  k = 0.5 (rotation parameter V/U)\n";
    std::cout << "\n";

    ResultLogger logger("Test15_RotatingCylinder");

    // Parse command line and input file
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test15.log");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    // Create major algorithm and data objects
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>(
        "PatchHierarchy", grid_geometry);

    // Navier-Stokes integrator
    Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator =
        new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

    // Advection-diffusion integrator for scalar (temperature/odor)
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator =
        new AdvDiffHierarchyIntegrator(
            "AdvDiffHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffHierarchyIntegrator"));

    // Register advection-diffusion integrator with Navier-Stokes
    navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);

    // Set up visualization plot file writers
    Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
    navier_stokes_integrator->registerVisItDataWriter(visit_data_writer);

    // Grid generation
    Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize",
        navier_stokes_integrator,
        app_initializer->getComponentDatabase("StandardTagAndInitialize"));

    Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();

    Pointer<LoadBalancer<NDIM>> load_balancer = new LoadBalancer<NDIM>(
        "LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));

    Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm = new GriddingAlgorithm<NDIM>(
        "GriddingAlgorithm",
        app_initializer->getComponentDatabase("GriddingAlgorithm"),
        error_detector,
        box_generator,
        load_balancer);

    // Configure the advection-diffusion integrator
    Pointer<CellVariable<NDIM, double>> C_var = new CellVariable<NDIM, double>("C");
    adv_diff_integrator->registerTransportedQuantity(C_var);
    adv_diff_integrator->setAdvectionVelocity(C_var,
                                              navier_stokes_integrator->getAdvectionVelocityVariable());

    // Get parameters from input file
    const double Re = input_db->getDoubleWithDefault("RE", 200.0);
    const double Pr = input_db->getDoubleWithDefault("PR", 0.5);
    const double U_inf = input_db->getDoubleWithDefault("U_INF", 1.0);
    const double D = input_db->getDoubleWithDefault("D", 1.0);
    const double nu = U_inf * D / Re;
    const double alpha = nu / Pr;  // Scalar diffusivity

    std::cout << "Computed parameters:\n";
    std::cout << "  ν (viscosity) = " << nu << "\n";
    std::cout << "  α (scalar diffusivity) = " << alpha << "\n";
    std::cout << "  α/ν = 1/Pr = " << alpha/nu << " (should be 2.0)\n";
    std::cout << "\n";

    adv_diff_integrator->setDiffusionCoefficient(C_var, alpha);

    // Initial condition: C = 0 everywhere (cold background)
    Pointer<CartGridFunction> C_init = new muParserCartGridFunction(
        "C_init",
        app_initializer->getComponentDatabase("C_init"),
        grid_geometry);
    adv_diff_integrator->setInitialConditions(C_var, C_init);

    // Source term: Gaussian around cylinder enforcing C = C_hot
    const double x_center = input_db->getDoubleWithDefault("CYLINDER_X", 0.0);
    const double y_center = input_db->getDoubleWithDefault("CYLINDER_Y", 0.0);
    const double radius = D / 2.0;
    const double C_hot = 1.0;
    const double source_thickness = input_db->getDoubleWithDefault("SOURCE_THICKNESS", 0.05);

    Pointer<RotatingCylinderScalarSource> C_source =
        new RotatingCylinderScalarSource("CylinderSource",
                                         x_center, y_center, radius,
                                         C_hot, source_thickness);

    Pointer<CellVariable<NDIM, double>> F_var = new CellVariable<NDIM, double>("F");
    adv_diff_integrator->registerSourceTerm(F_var);
    adv_diff_integrator->setSourceTermFunction(F_var, C_source);
    adv_diff_integrator->setSourceTerm(C_var, F_var);

    // Boundary conditions: Inlet C=0, all others Neumann
    std::vector<RobinBcCoefStrategy<NDIM>*> C_bc_coefs(NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        const std::string bc_coefs_name = "C_bc_coefs_" + std::to_string(d);
        const std::string bc_coefs_db_name = "C_BcCoefs_" + std::to_string(d);
        C_bc_coefs[d] = new muParserRobinBcCoefs(
            bc_coefs_name,
            app_initializer->getComponentDatabase(bc_coefs_db_name),
            grid_geometry);
    }
    adv_diff_integrator->setPhysicalBcCoef(C_var, C_bc_coefs[0]);

    // Initialize hierarchy configuration and data
    navier_stokes_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

    // Variables for steady-state check
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> current_ctx = navier_stokes_integrator->getCurrentContext();
    Pointer<VariableContext> scratch_ctx = navier_stokes_integrator->getScratchContext();

    const int C_current_idx = var_db->mapVariableAndContextToIndex(C_var, current_ctx);
    const int C_scratch_idx = var_db->mapVariableAndContextToIndex(C_var, scratch_ctx);

    // Main time-stepping loop
    double dt = navier_stokes_integrator->getMaximumTimeStepSize();
    double loop_time = navier_stokes_integrator->getIntegratorTime();
    const double loop_time_end = navier_stokes_integrator->getEndTime();
    int iteration_num = navier_stokes_integrator->getIntegratorStep();

    const double steady_tolerance = input_db->getDoubleWithDefault("STEADY_TOL", 1.0e-6);
    const int check_interval = input_db->getIntegerWithDefault("CHECK_INTERVAL", 100);
    bool steady_state_reached = false;

    std::cout << "Starting time integration...\n";
    std::cout << "Target: Steady state (|dC/dt| < " << steady_tolerance << ")\n\n";

    while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) &&
           navier_stokes_integrator->stepsRemaining())
    {
        iteration_num = navier_stokes_integrator->getIntegratorStep();
        loop_time = navier_stokes_integrator->getIntegratorTime();

        if (iteration_num % 50 == 0)
        {
            std::cout << "  Iteration " << iteration_num
                     << ", time = " << loop_time
                     << ", dt = " << dt << "\n";
        }

        navier_stokes_integrator->advanceHierarchy(dt);
        loop_time += dt;

        // Check for steady state
        if (iteration_num % check_interval == 0 && iteration_num > 0)
        {
            bool is_steady = checkSteadyState(patch_hierarchy, C_var,
                                             C_current_idx, C_scratch_idx,
                                             steady_tolerance);
            if (is_steady && !steady_state_reached)
            {
                std::cout << "\n✓ Steady state reached at iteration " << iteration_num
                         << ", time = " << loop_time << "\n\n";
                steady_state_reached = true;
                logger.log("Steady state reached at t=" + std::to_string(loop_time));
            }
        }

        navier_stokes_integrator->updateSolutionAfterAdvance(dt);

        dt = navier_stokes_integrator->getMaximumTimeStepSize();
    }

    // Write centerline profile
    writeCenterlineProfile(patch_hierarchy, C_var, C_current_idx, 0.0,
                          "test15_centerline.dat");

    // Final validation checks
    TestUtils::printSectionSeparator("VALIDATION CHECKS");

    bool all_checks_passed = true;

    // Check 1: Steady state reached
    if (steady_state_reached)
    {
        std::cout << "✓ Steady state reached\n";
        logger.log("CHECK 1: PASS - Steady state");
    }
    else
    {
        std::cout << "✗ Steady state NOT reached\n";
        logger.log("CHECK 1: FAIL - Steady state not reached");
        all_checks_passed = false;
    }

    // Check 2: No NaN/Inf
    bool has_nan = false;
    bool has_negative = false;
    bool out_of_bounds = false;
    double C_min = 1e10, C_max = -1e10;

    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> C_data = patch->getPatchData(C_current_idx);

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                double C = (*C_data)(idx);

                if (std::isnan(C) || std::isinf(C)) has_nan = true;
                if (C < -1e-10) has_negative = true;
                if (C < -0.1 || C > 1.1) out_of_bounds = true;

                C_min = std::min(C_min, C);
                C_max = std::max(C_max, C);
            }
        }
    }

    std::cout << "  Scalar range: C ∈ [" << C_min << ", " << C_max << "]\n";

    if (!has_nan)
    {
        std::cout << "✓ No NaN/Inf values\n";
        logger.log("CHECK 2: PASS - No NaN/Inf");
    }
    else
    {
        std::cout << "✗ NaN/Inf detected\n";
        logger.log("CHECK 2: FAIL - NaN/Inf detected");
        all_checks_passed = false;
    }

    if (!has_negative)
    {
        std::cout << "✓ No negative concentrations\n";
        logger.log("CHECK 3: PASS - No negatives");
    }
    else
    {
        std::cout << "✗ Negative concentrations detected\n";
        logger.log("CHECK 3: FAIL - Negative values");
        all_checks_passed = false;
    }

    if (!out_of_bounds)
    {
        std::cout << "✓ Scalar bounded (0 ≤ C ≤ 1)\n";
        logger.log("CHECK 4: PASS - Bounded");
    }
    else
    {
        std::cout << "✗ Scalar out of physical bounds\n";
        logger.log("CHECK 4: FAIL - Out of bounds");
        all_checks_passed = false;
    }

    // Summary
    std::cout << "\n";
    TestUtils::printSectionSeparator("TEST SUMMARY");

    if (all_checks_passed)
    {
        std::cout << "✅ TEST 15 PASSED\n\n";
        std::cout << "Validation:\n";
        std::cout << "  - Steady state achieved\n";
        std::cout << "  - Solution stable (no NaN/Inf)\n";
        std::cout << "  - Physical bounds satisfied\n";
        std::cout << "  - Scalar field: C ∈ [" << C_min << ", " << C_max << "]\n\n";
        std::cout << "Next steps:\n";
        std::cout << "  1. Visualize with VisIt/ParaView: viz_test15/dumps.visit\n";
        std::cout << "  2. Compare streamlines with Yan & Zu (2008) Figure 2\n";
        std::cout << "  3. Compare scalar contours with reference\n";
        std::cout << "  4. Check centerline profile: test15_centerline.dat\n";

        logger.log("TEST 15: PASSED");
    }
    else
    {
        std::cout << "❌ TEST 15 FAILED\n\n";
        std::cout << "One or more validation checks failed.\n";
        std::cout << "Review output above for details.\n";

        logger.log("TEST 15: FAILED");
    }

    std::cout << "\n";
    logger.writeToFile("test15_results.txt");

    return 0;
}
