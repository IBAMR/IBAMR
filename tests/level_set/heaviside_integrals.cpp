// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/LevelSetUtilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserCartGridFunction.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
void
check_scalar_functions()
{
    const double eps = std::numeric_limits<double>::epsilon();
    const double tiny = std::numeric_limits<double>::denorm_min();
    // (phi, H(phi), delta(phi)) at alpha = 1, computed to 100 decimal digits
    // from (1+phi+sin(pi*phi)/pi)/2 and (1+cos(pi*phi))/2, then rounded to double.
    // Hexadecimal literals and power-of-two scaling preserve the sampled inputs.
    const std::array<double, 3> references[] = {
        { -0x1.0000000000000p+0, 0x0.0p+0, 0x0.0p+0 },
        { -0x1.fffffffffffffp-1, 0x1.a51a6625307d3p-160, 0x1.3bd3cc9be45dep-105 },
        { -0x1.fffffff768fa1p-1, 0x1.04a67286b87bap-90, 0x1.6c1fc5d475b7ep-59 },
        { -0x1.ffffffaa19c47p-1, 0x1.fd151b048acb0p-81, 0x1.1c78d3cbe80d6p-52 },
        { -0x1.fffffca501acbp-1, 0x1.f1269bd788883p-71, 0x1.bc7cca9bcfc09p-46 },
        { -0x1.ffffde7210be9p-1, 0x1.e57fb43e0adf3p-61, 0x1.5b417e5030efdp-39 },
        { -0x1.fffeb074a771dp-1, 0x1.da1eb60367a9bp-51, 0x1.0f4b2aadf9bd2p-32 },
        { -0x1.fff2e48e8a71ep-1, 0x1.cf01fd997a7f2p-41, 0x1.a7e57275fe1f6p-26 },
        { -0x1.ff7ced916872bp-1, 0x1.c427e32dd18e2p-31, 0x1.4b2b2fc02d8d6p-19 },
        { -0x1.ccccccccccccdp-1, 0x1.ad162f0eed3f6p-11, 0x1.90f1ecbbab007p-6 },
        { -0x1.8000000000000p-1, 0x1.984deb20ccef7p-7, 0x1.2bec333018867p-3 },
        { -0x1.0000000000001p-1, 0x1.7419f246c6ef7p-4, 0x1.ffffffffffffdp-2 },
        { -0x1.0000000000000p-1, 0x1.7419f246c6efbp-4, 0x1.0000000000000p-1 },
        { -0x1.fffffffffffffp-2, 0x1.7419f246c6efdp-4, 0x1.0000000000001p-1 },
        { -0x1.999999999999ap-4, 0x1.9a7024b121103p-2, 0x1.f378709a22a80p-1 },
        { 0x0.0p+0, 0x1.0000000000000p-1, 0x1.0000000000000p+0 },
    };
    for (const std::array<double, 3>& reference : references)
    {
        for (const int exponent : { -1000, -333, 0, 333, 1000 })
        {
            const double alpha = std::ldexp(1.0, exponent);
            const double phi = std::ldexp(reference[0], exponent);
            const double ref_delta = std::ldexp(reference[2], -exponent);
            if (!(std::abs(smooth_heaviside(phi, alpha) - reference[1]) <= 32.0 * eps * reference[1] + 4.0 * tiny))
            {
                TBOX_ERROR("Heaviside regression: small phase-fraction accuracy\n");
            }
            if (!(std::abs(smooth_heaviside(-phi, alpha) - (1.0 - reference[1])) <= 32.0 * eps))
            {
                TBOX_ERROR("Heaviside regression: reflected phase-fraction accuracy\n");
            }
            if (!(std::abs(smooth_delta(phi, alpha) - ref_delta) <= 32.0 * eps * ref_delta + 4.0 * tiny))
            {
                TBOX_ERROR("Heaviside regression: delta accuracy including underflow\n");
            }
        }
    }
    for (const double alpha : { 1.0e-300, 1.0e-100, 0.1, 1.0, 10.0, 1.0e100, 1.0e300 })
    {
        std::vector<double> samples{
            -alpha, std::nextafter(-alpha, -2.0 * alpha), std::nextafter(-alpha, 0.0),       0.0,
            alpha,  std::nextafter(alpha, 0.0),           std::nextafter(alpha, 2.0 * alpha)
        };
        for (const double eta : { 0.0, 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 0.1, 0.5, 0.9, 1.0 })
        {
            const double phi = alpha * (-1.0 + eta);
            samples.push_back(phi);
            samples.push_back(-phi);
        }
        for (const double phi : { -0.5 * alpha, 0.5 * alpha })
        {
            samples.push_back(std::nextafter(phi, 0.0));
            samples.push_back(std::nextafter(phi, 2.0 * phi));
        }
        std::sort(samples.begin(), samples.end());
        double previous = 0.0;
        for (const double phi : samples)
        {
            const double h = smooth_heaviside(phi, alpha);
            const double delta = smooth_delta(phi, alpha);
            if (!(std::isfinite(h) && h >= 0.0 && h <= 1.0))
            {
                TBOX_ERROR("Heaviside regression: phase fraction outside [0,1]\n");
            }
            if (!(std::isfinite(delta) && delta >= 0.0))
            {
                TBOX_ERROR("Heaviside regression: invalid delta\n");
            }
            if (!(h + 4.0 * eps >= previous))
            {
                TBOX_ERROR("Heaviside regression: nonmonotone phase fraction\n");
            }
            if (!(std::abs(h + smooth_heaviside(-phi, alpha) - 1.0) <= 2.0 * eps))
            {
                TBOX_ERROR("Heaviside regression: complement symmetry\n");
            }
            if (!(delta == smooth_delta(-phi, alpha)))
            {
                TBOX_ERROR("Heaviside regression: delta symmetry\n");
            }
            previous = h;
            if (phi <= -alpha || phi >= alpha)
            {
                if (!(h == (phi < 0.0 ? 0.0 : 1.0) && delta == 0.0))
                {
                    TBOX_ERROR("Heaviside regression: inexact cutoff or saturation\n");
                }
            }
        }
        if (!(smooth_heaviside(0.0, alpha) == 0.5))
        {
            TBOX_ERROR("Heaviside regression: center value\n");
        }
        previous = 0.0;
        for (int i = 0; i <= 4096; ++i)
        {
            const double phi = alpha * (-1.0 + i / 2048.0);
            const double h = smooth_heaviside(phi, alpha);
            if (!(h + 4.0 * eps >= previous))
            {
                TBOX_ERROR("Heaviside regression: interior monotonicity\n");
            }
            previous = h;
        }
        for (const double r : { -0.8, -0.3, 0.0, 0.3, 0.8 })
        {
            const double step = 1.0e-5;
            const double derivative =
                (smooth_heaviside(alpha * (r + step), alpha) - smooth_heaviside(alpha * (r - step), alpha)) /
                (2.0 * step);
            if (!(std::abs(derivative - alpha * smooth_delta(alpha * r, alpha)) <= 5.0e-9))
            {
                TBOX_ERROR("Heaviside regression: resolved analytic derivative\n");
            }
        }
    }
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    {
        Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
        plog << std::setprecision(16);
        const bool check_scalars = app->getInputDatabase()->getBoolWithDefault("check_scalars", true);
        if (check_scalars)
        {
            check_scalar_functions();
            plog << "H(-1), H(0), H(1): " << smooth_heaviside(-1.0, 1.0) << ' ' << smooth_heaviside(0.0, 1.0) << ' '
                 << smooth_heaviside(1.0, 1.0) << '\n';
            plog << "H(-1+1e-6), delta(-1+1e-6): " << smooth_heaviside(-1.0 + 1.0e-6, 1.0) << ' '
                 << smooth_delta(-1.0 + 1.0e-6, 1.0) << '\n';
        }
        Pointer<AdvDiffHierarchyIntegrator> integrator =
            new AdvDiffSemiImplicitHierarchyIntegrator("AdvDiff", app->getComponentDatabase("AdvDiff"));
        Pointer<CartesianGridGeometry<NDIM>> geometry =
            new CartesianGridGeometry<NDIM>("CartesianGeometry", app->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("Hierarchy", geometry);
        Pointer<StandardTagAndInitialize<NDIM>> tagger = new StandardTagAndInitialize<NDIM>(
            "Tagger", integrator, app->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> boxes = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load = new LoadBalancer<NDIM>("Load", app->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding = new GriddingAlgorithm<NDIM>(
            "Gridding", app->getComponentDatabase("GriddingAlgorithm"), tagger, boxes, load);
        Pointer<FaceVariable<NDIM, double>> velocity = new FaceVariable<NDIM, double>("velocity");
        integrator->registerAdvectionVelocity(velocity);
        integrator->setAdvectionVelocityFunction(
            velocity,
            new muParserCartGridFunction("velocity_function", app->getComponentDatabase("Velocity"), geometry));
        LocationIndexRobinBcCoefs<NDIM> bc;
        for (int face = 0; face < 2 * NDIM; ++face)
        {
            bc.setBoundarySlope(face, 0.0);
        }
        std::vector<Pointer<CellVariable<NDIM, double>>> variables;
        for (int k = 0; k < 2; ++k)
        {
            Pointer<CellVariable<NDIM, double>> var = new CellVariable<NDIM, double>("phi" + std::to_string(k));
            integrator->registerTransportedQuantity(var);
            integrator->setDiffusionCoefficient(var, 0.0);
            integrator->setAdvectionVelocity(var, velocity);
            integrator->setPhysicalBcCoef(var, &bc);
            integrator->setInitialConditions(var,
                                             new muParserCartGridFunction("initial" + std::to_string(k),
                                                                          app->getComponentDatabase("Initial"),
                                                                          geometry));
            variables.push_back(var);
        }
        integrator->initializePatchHierarchy(hierarchy, gridding);
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int phi_idx = var_db->mapVariableAndContextToIndex(variables[0], integrator->getCurrentContext());
        const int psi_idx = var_db->mapVariableAndContextToIndex(variables[1], integrator->getCurrentContext());
        Pointer<HierarchyMathOps> ops = integrator->getHierarchyMathOps();
        const int weight_idx = ops->getCellWeightPatchDescriptorIndex();
        LevelSetUtilities::LevelSetContainer two_phase(integrator, variables[0]);
        LevelSetUtilities::LevelSetContainer three_phase(integrator, variables);
        double domain_volume = 1.0;
        for (int d = 0; d < NDIM; ++d)
        {
            domain_volume *= geometry->getXUpper()[d] - geometry->getXLower()[d];
        }
        int local_cells = 0;
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                local_cells += level->getPatch(p())->getBox().size();
            }
        }
        // Each nonnegative phase sum accumulates at most n local terms. Allow
        // gamma_n = n*eps/(1-n*eps) for each side of an identity, plus
        // reduction and fraction rounding.
        const int n = 2 * IBTK_MPI::maxReduction(local_cells) + 2 * IBTK_MPI::getNodes() + 8;
        const double n_eps = n * std::numeric_limits<double>::epsilon();
        const double tol = n_eps / (1.0 - n_eps) * domain_volume;
        double maximum_sum_error = 0.0;
        for (const double ncells : { 0.5, 1.0, 2.0 })
        {
            two_phase.setInterfaceHalfWidth(ncells);
            three_phase.setInterfaceHalfWidth(ncells);
            for (const double eta : { 0.0, 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 0.5, 1.0 })
            {
                double capacity = 0.0;
                for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
                {
                    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM>> patch = level->getPatch(p());
                        Pointer<CartesianPatchGeometry<NDIM>> geom = patch->getPatchGeometry();
                        double cell_volume = 1.0;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            cell_volume *= geom->getDx()[d];
                        }
                        const double alpha = ncells * std::pow(cell_volume, 1.0 / NDIM);
                        Pointer<CellData<NDIM, double>> phi = patch->getPatchData(phi_idx);
                        Pointer<CellData<NDIM, double>> psi = patch->getPatchData(psi_idx);
                        Pointer<CellData<NDIM, double>> weight = patch->getPatchData(weight_idx);
                        for (Box<NDIM>::Iterator i(patch->getBox()); i; i++)
                        {
                            (*phi)(i()) = alpha * (1.0 - eta);
                            (*psi)(i()) = alpha * (-1.0 + eta);
                            const double h = smooth_heaviside((*psi)(i()), alpha);
                            const double complement = smooth_heaviside(-(*psi)(i()), alpha);
                            const double contribution =
                                smooth_delta((*psi)(i()), alpha) * h * complement * (*weight)(i());
                            if (!(h >= 0.0 && h <= 1.0 && complement >= 0.0 && complement <= 1.0 &&
                                  std::isfinite(contribution) && contribution >= 0.0))
                            {
                                TBOX_ERROR(
                                    "Heaviside regression: negative live phase or redistribution contribution\n");
                            }
                            capacity += h * (*weight)(i());
                        }
                    }
                }
                capacity = IBTK_MPI::sumReduction(capacity);
                const std::vector<double> v2 = LevelSetUtilities::computeHeavisideIntegrals2PhaseFlows(two_phase);
                const std::vector<double> v3 = LevelSetUtilities::computeHeavisideIntegrals3PhaseFlows(three_phase);
                if (!(v2[0] >= 0.0 && v2[1] >= 0.0 && v2[2] >= 0.0 && v3[0] >= 0.0 && v3[1] >= 0.0 && v3[2] >= 0.0 &&
                      v3[3] >= 0.0))
                {
                    TBOX_ERROR("Heaviside regression: negative production integral\n");
                }
                maximum_sum_error = std::max({ maximum_sum_error,
                                               std::abs(v2[0] + v2[1] - domain_volume),
                                               std::abs(v3[0] + v3[1] - capacity),
                                               std::abs(v3[0] + v3[1] + v3[2] - domain_volume) });
                if (!(std::abs(v2[0] + v2[1] - domain_volume) <= tol))
                {
                    TBOX_ERROR("Heaviside regression: two-phase composite identity\n");
                }
                if (!(std::abs(v3[0] + v3[1] - capacity) <= tol))
                {
                    TBOX_ERROR("Heaviside regression: three-phase fluid capacity\n");
                }
                if (!(std::abs(v3[0] + v3[1] + v3[2] - domain_volume) <= tol))
                {
                    TBOX_ERROR("Heaviside regression: three-phase composite identity\n");
                }
                if (ncells == 1.0 && (eta == 0.0 || eta == 1.0e-6 || eta == 0.5))
                {
                    plog << "eta = " << eta << "; normalized gas, liquid, solid volumes: " << v3[0] / domain_volume
                         << ' ' << v3[1] / domain_volume << ' ' << v3[2] / domain_volume << '\n';
                }
                if (eta == 0.0)
                {
                    if (!(capacity == 0.0 && v3[0] == 0.0 && v3[1] == 0.0 && v2[0] == 0.0))
                    {
                        TBOX_ERROR("Heaviside regression: nonzero cutoff volume or capacity\n");
                    }
                }
                else if (eta <= 1.0e-5)
                {
                    if (!(v2[0] > 0.0 && v3[0] > 0.0 && v3[1] > 0.0))
                    {
                        TBOX_ERROR("Heaviside regression: small complementary phase was lost\n");
                    }
                }
            }
        }
        plog << "maximum normalized phase-sum error: " << maximum_sum_error / domain_volume << '\n';
    }
}
