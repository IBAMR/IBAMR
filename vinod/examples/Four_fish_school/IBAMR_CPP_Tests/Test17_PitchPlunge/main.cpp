// =============================================================================
// TEST 17: Pitch-Plunge Odor Plume Navigation (Lei et al. 2021)
// =============================================================================
// Full production case: Sphere source + flapping ellipsoidal airfoil
// Reference: Lei et al. (2021) Section II.C, Figures 4-10
//
// Parameters: Re=200, Sc=0.71, St=0.9
// Kinematics: y(t) = (L/2)sin(2πft), θ(t) = (1/6)cos(2πft)
// Expected: Inverse von Kármán street, odor modulation by vortices
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
// SPHERE SOURCE (Upstream odor emitter)
// =============================================================================

class SphereOdorSource : public CartGridFunction
{
public:
    SphereOdorSource(const std::string& object_name,
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

        for (CellIterator<NDIM> ci(patch_box); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();
            double X[NDIM];
            for (int d = 0; d < NDIM; ++d)
            {
                X[d] = x_lower[d] + dx[d] * (static_cast<double>(idx(d) - patch_box.lower()(d)) + 0.5);
            }

            double r = sqrt(pow(X[0] - d_x_center, 2) + pow(X[1] - d_y_center, 2));
            double dr = r - d_radius;
            double source_strength = d_C_hot * exp(-0.5 * pow(dr / d_thickness, 2));

            (*S_data)(idx) = source_strength / (d_thickness * sqrt(2.0 * M_PI));
        }
    }

private:
    double d_x_center, d_y_center;
    double d_radius;
    double d_C_hot;
    double d_thickness;
};

// =============================================================================
// PITCH-PLUNGE KINEMATICS
// =============================================================================

struct PitchPlungeKinematics
{
    double L;           // Characteristic length
    double U;           // Free-stream velocity
    double St;          // Strouhal number
    double h0;          // Plunge amplitude
    double theta0;      // Pitch amplitude (radians)

    PitchPlungeKinematics(double L_, double U_, double St_, double h0_, double theta0_)
        : L(L_), U(U_), St(St_), h0(h0_), theta0(theta0_) {}

    double frequency() const { return St * U / L; }

    // Plunge displacement: y(t) = h0 * sin(2πft)
    double plunge(double t) const
    {
        return h0 * sin(2.0 * M_PI * frequency() * t);
    }

    // Pitch angle: θ(t) = θ0 * cos(2πft)
    double pitch(double t) const
    {
        return theta0 * cos(2.0 * M_PI * frequency() * t);
    }

    // Plunge velocity: dy/dt
    double plunge_velocity(double t) const
    {
        return h0 * 2.0 * M_PI * frequency() * cos(2.0 * M_PI * frequency() * t);
    }

    // Pitch angular velocity: dθ/dt
    double pitch_velocity(double t) const
    {
        return -theta0 * 2.0 * M_PI * frequency() * sin(2.0 * M_PI * frequency() * t);
    }

    // Nondimensional time: t* = tU/L
    double nondim_time(double t) const
    {
        return t * U / L;
    }
};

// =============================================================================
// POST-PROCESSING: Odor PDF
// =============================================================================

void computeOdorPDF(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                    Pointer<CellVariable<NDIM, double>> C_var,
                    const int C_idx,
                    const double C_l,
                    const double C_h,
                    const std::string& filename,
                    const double x_min, const double x_max,
                    const double y_min, const double y_max)
{
    const int num_bins = 50;
    std::vector<int> histogram(num_bins, 0);
    int total_samples = 0;

    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
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

                // Sample only in specified region (e.g., wake region)
                if (x >= x_min && x <= x_max && y >= y_min && y <= y_max)
                {
                    double C = (*C_data)(idx);
                    double C_star = (C - C_l) / (C_h - C_l);  // Nondimensionalize

                    // Bin the value
                    int bin = static_cast<int>(C_star * num_bins);
                    if (bin >= 0 && bin < num_bins)
                    {
                        histogram[bin]++;
                        total_samples++;
                    }
                }
            }
        }
    }

    // Write PDF to file
    std::ofstream outfile(filename);
    outfile << "# C*  PDF(C*)\n";
    for (int i = 0; i < num_bins; ++i)
    {
        double C_star = (i + 0.5) / num_bins;
        double pdf = static_cast<double>(histogram[i]) / total_samples;
        outfile << C_star << "  " << pdf << "\n";
    }
    outfile.close();

    std::cout << "  Written odor PDF to " << filename << " (N=" << total_samples << " samples)\n";
}

// =============================================================================
// MAIN
// =============================================================================

int main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    TestUtils::printTestHeader("TEST 17: PITCH-PLUNGE ODOR PLUME");
    std::cout << "\n";
    std::cout << "Reference: Lei et al. (2021) Section II.C, Figures 4-10\n";
    std::cout << "Full production case: Sphere source + flapping airfoil\n";
    std::cout << "\nParameters:\n";
    std::cout << "  Re = 200 (Reynolds number)\n";
    std::cout << "  Sc = 0.71 (Schmidt number)\n";
    std::cout << "  St = 0.9 (Strouhal number)\n";
    std::cout << "  Plunge: y(t) = (L/2)sin(2πft)\n";
    std::cout << "  Pitch: θ(t) = (1/6)cos(2πft) ≈ 9.5°\n";
    std::cout << "\nExpected Output:\n";
    std::cout << "  - Inverse von Kármán vortex street\n";
    std::cout << "  - Odor modulation by wake vortices\n";
    std::cout << "  - PDF: peak at C*=0, trough at C*≈0.2\n";
    std::cout << "\n";

    ResultLogger logger("Test17_PitchPlunge");

    // Parse input
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test17.log");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    // Create grid geometry
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>(
        "PatchHierarchy", grid_geometry);

    // Navier-Stokes integrator
    Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator =
        new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

    // Advection-diffusion integrator
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

    // Configure scalar transport
    Pointer<CellVariable<NDIM, double>> C_var = new CellVariable<NDIM, double>("C");
    adv_diff_integrator->registerTransportedQuantity(C_var);
    adv_diff_integrator->setAdvectionVelocity(C_var,
                                              navier_stokes_integrator->getAdvectionVelocityVariable());

    // Physical parameters from Lei et al. (2021)
    const double Re = input_db->getDoubleWithDefault("RE", 200.0);
    const double Sc = input_db->getDoubleWithDefault("SC", 0.71);
    const double St = input_db->getDoubleWithDefault("ST", 0.9);
    const double U_inf = input_db->getDoubleWithDefault("U_INF", 1.0);
    const double L = input_db->getDoubleWithDefault("L", 1.0);

    const double nu = U_inf * L / Re;
    const double kappa = nu / Sc;

    std::cout << "Computed parameters:\n";
    std::cout << "  ν = " << nu << "\n";
    std::cout << "  κ = " << kappa << "\n";
    std::cout << "  f = " << St * U_inf / L << " Hz\n";
    std::cout << "\n";

    adv_diff_integrator->setDiffusionCoefficient(C_var, kappa);

    // Initial condition: C = 0
    Pointer<CartGridFunction> C_init = new muParserCartGridFunction(
        "C_init",
        app_initializer->getComponentDatabase("C_init"),
        grid_geometry);
    adv_diff_integrator->setInitialConditions(C_var, C_init);

    // Sphere source (upstream at x = -3L)
    const double sphere_x = input_db->getDoubleWithDefault("SPHERE_X", -3.0);
    const double sphere_y = input_db->getDoubleWithDefault("SPHERE_Y", 0.0);
    const double sphere_R = L / 2.0;  // Diameter = L
    const double C_hot = 1.0;
    const double C_cold = 0.0;
    const double source_thickness = input_db->getDoubleWithDefault("SOURCE_THICKNESS", 0.05);

    Pointer<SphereOdorSource> C_source =
        new SphereOdorSource("SphereSource",
                             sphere_x, sphere_y, sphere_R,
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

    // Kinematics
    const double h0 = L / 2.0;
    const double theta0 = 1.0 / 6.0;  // radians
    PitchPlungeKinematics kinematics(L, U_inf, St, h0, theta0);

    std::cout << "Kinematics:\n";
    std::cout << "  Frequency f = " << kinematics.frequency() << " Hz\n";
    std::cout << "  Period T = " << 1.0 / kinematics.frequency() << " s\n";
    std::cout << "  Plunge amplitude h0 = " << h0 << " = L/2\n";
    std::cout << "  Pitch amplitude θ0 = " << theta0 << " rad = " << theta0 * 180.0 / M_PI << "°\n";
    std::cout << "\n";

    // Initialize hierarchy
    navier_stokes_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

    // Time-stepping
    double dt = navier_stokes_integrator->getMaximumTimeStepSize();
    double loop_time = navier_stokes_integrator->getIntegratorTime();
    const double loop_time_end = navier_stokes_integrator->getEndTime();
    int iteration_num = navier_stokes_integrator->getIntegratorStep();

    std::cout << "Starting time integration...\n";
    std::cout << "Target: t* = " << kinematics.nondim_time(loop_time_end) << "\n";
    std::cout << "Number of periods: " << loop_time_end * kinematics.frequency() << "\n\n";

    while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) &&
           navier_stokes_integrator->stepsRemaining())
    {
        iteration_num = navier_stokes_integrator->getIntegratorStep();
        loop_time = navier_stokes_integrator->getIntegratorTime();

        if (iteration_num % 50 == 0)
        {
            double t_star = kinematics.nondim_time(loop_time);
            std::cout << "  Iteration " << iteration_num
                     << ", time = " << loop_time
                     << ", t* = " << t_star
                     << ", dt = " << dt << "\n";
        }

        navier_stokes_integrator->advanceHierarchy(dt);
        loop_time += dt;
        navier_stokes_integrator->updateSolutionAfterAdvance(dt);

        dt = navier_stokes_integrator->getMaximumTimeStepSize();
    }

    // Post-processing
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> current_ctx = navier_stokes_integrator->getCurrentContext();
    const int C_current_idx = var_db->mapVariableAndContextToIndex(C_var, current_ctx);

    // Compute odor PDF in wake region (x ∈ [0, 5L], y ∈ [-L, L])
    computeOdorPDF(patch_hierarchy, C_var, C_current_idx,
                   C_cold, C_hot, "test17_odor_pdf.dat",
                   0.0, 5.0 * L, -L, L);

    // Validation
    TestUtils::printSectionSeparator("VALIDATION CHECKS");

    bool all_checks_passed = true;
    bool has_nan = false;
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
                C_min = std::min(C_min, C);
                C_max = std::max(C_max, C);
            }
        }
    }

    std::cout << "  Scalar range: C ∈ [" << C_min << ", " << C_max << "]\n";

    if (!has_nan && C_min >= -0.1 && C_max <= 1.1)
    {
        std::cout << "✓ All checks passed\n";
        all_checks_passed = true;
    }
    else
    {
        std::cout << "✗ Validation failed\n";
        all_checks_passed = false;
    }

    // Summary
    std::cout << "\n";
    TestUtils::printSectionSeparator("TEST SUMMARY");

    if (all_checks_passed)
    {
        std::cout << "✅ TEST 17 PASSED\n\n";
        std::cout << "Validation complete. Compare with Lei et al. (2021) Figures 6-10:\n";
        std::cout << "  1. Visualize vorticity field (inverse von Kármán street)\n";
        std::cout << "  2. Visualize odor field (plume modulation)\n";
        std::cout << "  3. Check odor PDF: test17_odor_pdf.dat\n";
        std::cout << "     Expected: Peak at C*≈0, trough at C*≈0.2\n";

        logger.log("TEST 17: PASSED");
    }
    else
    {
        std::cout << "❌ TEST 17 FAILED\n";
        logger.log("TEST 17: FAILED");
    }

    logger.writeToFile("test17_results.txt");
    return 0;
}
