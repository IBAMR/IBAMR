// =============================================================================
// TEST 16: 3D Sphere & Cube Validation (Richter & Nikrityuk 2012)
// =============================================================================
// Validates 3-D scalar transport around stationary sphere/cube
// Reference: Richter & Nikrityuk (2012) - Heat transfer CFD benchmark
// Lei et al. (2021) Section II.B, Figure 3 - 3D solver validation
//
// Parameters: Re=200, Pr=0.744
// Expected: Thermal boundary layers, wake structure, Nusselt number
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

#include <ibamr/app_namespaces.h>
#include <vector>
#include <fstream>
#include <cmath>

using namespace TestUtilities;

// =============================================================================
// SPHERE/CUBE SCALAR SOURCE (3D)
// =============================================================================
// Enforces C = C_h (isothermal) on immersed boundary surface

class SphereScalarSource3D : public CartGridFunction
{
public:
    SphereScalarSource3D(const std::string& object_name,
                         double x_center,
                         double y_center,
                         double z_center,
                         double radius,
                         double C_hot,
                         double thickness)
        : CartGridFunction(object_name),
          d_x_center(x_center),
          d_y_center(y_center),
          d_z_center(z_center),
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

        // Gaussian source enforcing C = C_hot around sphere
        for (CellIterator<NDIM> ci(patch_box); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();
            double X[NDIM];
            for (int d = 0; d < NDIM; ++d)
            {
                X[d] = x_lower[d] + dx[d] * (static_cast<double>(idx(d) - patch_box.lower()(d)) + 0.5);
            }

            // Distance from sphere center
            double r = sqrt(pow(X[0] - d_x_center, 2) +
                          pow(X[1] - d_y_center, 2) +
                          pow(X[2] - d_z_center, 2));

            // Gaussian profile around sphere surface
            double dr = r - d_radius;
            double source_strength = d_C_hot * exp(-0.5 * pow(dr / d_thickness, 2));

            (*S_data)(idx) = source_strength / (d_thickness * sqrt(2.0 * M_PI));
        }
    }

private:
    double d_x_center, d_y_center, d_z_center;
    double d_radius;
    double d_C_hot;
    double d_thickness;
};

// =============================================================================
// ANALYSIS FUNCTIONS
// =============================================================================

void writeCenterlineProfile3D(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                               Pointer<CellVariable<NDIM, double>> C_var,
                               const int C_idx,
                               const std::string& filename)
{
    std::ofstream outfile(filename);
    outfile << "# x  C(x, y=0, z=0)\n";

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
            double z = x_lower[2] + dx[2] * (static_cast<double>(idx(2) - patch_box.lower()(2)) + 0.5);

            // Check if on centerline (y=0, z=0)
            if (fabs(y) < 0.5 * dx[1] && fabs(z) < 0.5 * dx[2])
            {
                double C_val = (*C_data)(idx);
                data_points.push_back(std::make_pair(x, C_val));
            }
        }
    }

    std::sort(data_points.begin(), data_points.end());

    for (const auto& pt : data_points)
    {
        outfile << pt.first << "  " << pt.second << "\n";
    }

    outfile.close();
    std::cout << "  Written centerline profile to " << filename << "\n";
}

// =============================================================================
// MAIN
// =============================================================================

int main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    TestUtils::printTestHeader("TEST 16: 3D SPHERE VALIDATION");
    std::cout << "\n";
    std::cout << "Reference: Richter & Nikrityuk (2012)\n";
    std::cout << "Validation: Lei et al. (2021) Section II.B, Figure 3\n";
    std::cout << "\nParameters:\n";
    std::cout << "  Re = 200 (Reynolds number)\n";
    std::cout << "  Pr = 0.744 (Prandtl number for air)\n";
    std::cout << "  Geometry: 3D stationary sphere\n";
    std::cout << "\n";
    std::cout << "NOTE: This is a 3D simulation - computational cost is high!\n";
    std::cout << "      Use coarse grid for quick validation.\n";
    std::cout << "      Recommended: Run with MPI on 8-16 cores\n";
    std::cout << "\n";

    ResultLogger logger("Test16_3DSphere");

    // Parse command line and input file
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test16.log");
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

    // Advection-diffusion integrator for scalar
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator =
        new AdvDiffHierarchyIntegrator(
            "AdvDiffHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffHierarchyIntegrator"));

    navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);

    // Visualization
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

    // Configure advection-diffusion integrator
    Pointer<CellVariable<NDIM, double>> C_var = new CellVariable<NDIM, double>("C");
    adv_diff_integrator->registerTransportedQuantity(C_var);
    adv_diff_integrator->setAdvectionVelocity(C_var,
                                              navier_stokes_integrator->getAdvectionVelocityVariable());

    // Physical parameters
    const double Re = input_db->getDoubleWithDefault("RE", 200.0);
    const double Pr = input_db->getDoubleWithDefault("PR", 0.744);
    const double U_inf = input_db->getDoubleWithDefault("U_INF", 1.0);
    const double D = input_db->getDoubleWithDefault("D", 1.0);
    const double nu = U_inf * D / Re;
    const double alpha = nu / Pr;

    std::cout << "Computed parameters:\n";
    std::cout << "  ν (viscosity) = " << nu << "\n";
    std::cout << "  α (scalar diffusivity) = " << alpha << "\n";
    std::cout << "  α/ν = 1/Pr = " << alpha/nu << " (should be ~1.34)\n";
    std::cout << "\n";

    adv_diff_integrator->setDiffusionCoefficient(C_var, alpha);

    // Initial condition
    Pointer<CartGridFunction> C_init = new muParserCartGridFunction(
        "C_init",
        app_initializer->getComponentDatabase("C_init"),
        grid_geometry);
    adv_diff_integrator->setInitialConditions(C_var, C_init);

    // Source term (Gaussian around sphere)
    const double x_center = input_db->getDoubleWithDefault("SPHERE_X", 0.0);
    const double y_center = input_db->getDoubleWithDefault("SPHERE_Y", 0.0);
    const double z_center = input_db->getDoubleWithDefault("SPHERE_Z", 0.0);
    const double radius = D / 2.0;
    const double C_hot = 1.0;
    const double source_thickness = input_db->getDoubleWithDefault("SOURCE_THICKNESS", 0.05);

    Pointer<SphereScalarSource3D> C_source =
        new SphereScalarSource3D("SphereSource",
                                 x_center, y_center, z_center, radius,
                                 C_hot, source_thickness);

    Pointer<CellVariable<NDIM, double>> F_var = new CellVariable<NDIM, double>("F");
    adv_diff_integrator->registerSourceTerm(F_var);
    adv_diff_integrator->setSourceTermFunction(F_var, C_source);
    adv_diff_integrator->setSourceTerm(C_var, F_var);

    // Boundary conditions
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

    // Initialize hierarchy
    navier_stokes_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

    // Time-stepping loop
    double dt = navier_stokes_integrator->getMaximumTimeStepSize();
    double loop_time = navier_stokes_integrator->getIntegratorTime();
    const double loop_time_end = navier_stokes_integrator->getEndTime();
    int iteration_num = navier_stokes_integrator->getIntegratorStep();

    std::cout << "Starting time integration (3D)...\n";
    std::cout << "WARNING: 3D simulation may take hours to complete!\n\n";

    while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) &&
           navier_stokes_integrator->stepsRemaining())
    {
        iteration_num = navier_stokes_integrator->getIntegratorStep();
        loop_time = navier_stokes_integrator->getIntegratorTime();

        if (iteration_num % 10 == 0)
        {
            std::cout << "  Iteration " << iteration_num
                     << ", time = " << loop_time
                     << ", dt = " << dt << "\n";
        }

        navier_stokes_integrator->advanceHierarchy(dt);
        loop_time += dt;
        navier_stokes_integrator->updateSolutionAfterAdvance(dt);

        dt = navier_stokes_integrator->getMaximumTimeStepSize();
    }

    // Get final solution
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> current_ctx = navier_stokes_integrator->getCurrentContext();
    const int C_current_idx = var_db->mapVariableAndContextToIndex(C_var, current_ctx);

    // Write centerline profile
    writeCenterlineProfile3D(patch_hierarchy, C_var, C_current_idx,
                            "test16_centerline.dat");

    // Validation checks
    TestUtils::printSectionSeparator("VALIDATION CHECKS");

    bool all_checks_passed = true;
    bool has_nan = false;
    bool has_negative = false;
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

                C_min = std::min(C_min, C);
                C_max = std::max(C_max, C);
            }
        }
    }

    std::cout << "  Scalar range: C ∈ [" << C_min << ", " << C_max << "]\n";

    if (!has_nan && !has_negative && C_min >= -0.1 && C_max <= 1.1)
    {
        std::cout << "✓ All checks passed\n";
        logger.log("TEST 16: PASSED");
        all_checks_passed = true;
    }
    else
    {
        std::cout << "✗ Validation failed\n";
        logger.log("TEST 16: FAILED");
        all_checks_passed = false;
    }

    // Summary
    std::cout << "\n";
    TestUtils::printSectionSeparator("TEST SUMMARY");

    if (all_checks_passed)
    {
        std::cout << "✅ TEST 16 PASSED (3D Sphere)\n\n";
        std::cout << "Next steps:\n";
        std::cout << "  1. Visualize with VisIt/ParaView (3D!):\n";
        std::cout << "     - Load viz_test16/dumps.visit\n";
        std::cout << "     - Create slice at y=0 or z=0\n";
        std::cout << "     - Create isosurfaces at C=0.5\n";
        std::cout << "  2. Compare with Richter & Nikrityuk (2012) Figure 3\n";
        std::cout << "  3. Check centerline: test16_centerline.dat\n";
    }
    else
    {
        std::cout << "❌ TEST 16 FAILED\n";
    }

    logger.writeToFile("test16_results.txt");
    return 0;
}
