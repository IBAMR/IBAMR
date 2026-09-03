// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/AppInitializer.h>
#include <ibtk/IBKernel.h>
#include <ibtk/IBKernelTensorProduct.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/LEInteractor.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <SideData.h>
#include <SideIterator.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include <ibtk/app_namespaces.h>

namespace
{
void
require(bool valid, const char* message)
{
    if (!valid) TBOX_ERROR(message << '\n');
}

void
check_patch(Pointer<Patch<NDIM>> patch, bool expect_error)
{
    const Box<NDIM>& box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
    const double* const dx = geometry->getDx();
    const double* const x_lower = geometry->getXLower();
    Pointer<SideData<NDIM, double>> field = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(5));
    Pointer<SideData<NDIM, double>> string_spread = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(5));
    Pointer<SideData<NDIM, double>> typed_spread = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(5));

    std::vector<double> X(NDIM), force(NDIM), string_value(NDIM), typed_value(NDIM);
    double cell_volume = 1.0;
    for (int d = 0; d < NDIM; ++d)
    {
        X[d] = x_lower[d] + (3.27 + 0.17 * d) * dx[d];
        force[d] = 1.3 + 0.7 * d;
        cell_volume *= dx[d];
    }
    for (int component = 0; component < NDIM; ++component)
        for (SideIterator<NDIM> i(field->getGhostBox(), component); i; i++)
        {
            double value = 0.4 + component;
            for (int d = 0; d < NDIM; ++d) value += std::sin((0.13 + 0.09 * d) * (i()(d) + 0.3 * component));
            (*field)(i()) = value;
        }

    const IBKernelTensorProduct unsupported({ IBKernel("UNREGISTERED_KERNEL") });
    require(!LEInteractor::isKnownKernel(unsupported), "unsupported typed query must return false");
    if (expect_error)
    {
        std::ofstream("output") << "LEInteractor rejects unsupported kernel descriptions before interpolation\n";
        string_value.assign(NDIM, 17.0);
        LEInteractor::interpolate(string_value, NDIM, X, NDIM, field, patch, box, unsupported);
        TBOX_ERROR("LEInteractor unsupported-kernel test continued unexpectedly\n");
    }

    for (const std::string& name : { "IB_4", "piecewise_constant", "discontinuous_linear", "composite_bspline_23" })
    {
        const IBKernelTensorProduct kernel(name);
        require(LEInteractor::isKnownKernel(name) && LEInteractor::isKnownKernel(kernel),
                "known kernel query rejected");
        require(LEInteractor::getStencilSize(name) == LEInteractor::getStencilSize(kernel),
                "typed stencil query differs");
        LEInteractor::interpolate(string_value, NDIM, X, NDIM, field, patch, box, name);
        LEInteractor::interpolate(typed_value, NDIM, X, NDIM, field, patch, box, kernel);
        require(string_value == typed_value, "typed/string interpolation differs");
        string_spread->fillAll(0.0);
        typed_spread->fillAll(0.0);
        LEInteractor::spread(string_spread, force, NDIM, X, NDIM, patch, box, name);
        LEInteractor::spread(typed_spread, force, NDIM, X, NDIM, patch, box, kernel);
        for (int component = 0; component < NDIM; ++component)
            for (SideIterator<NDIM> i(field->getGhostBox(), component); i; i++)
                require((*string_spread)(i()) == (*typed_spread)(i()), "typed/string spreading differs");
    }

    for (const std::string& name : { "", "BAD NAME", "UNREGISTERED_KERNEL", "COMPOSITE_BSPLINE_12" })
        require(!LEInteractor::isKnownKernel(name), "invalid or unsupported string query must return false");
    require(!LEInteractor::isKnownKernel(
                IBKernelTensorProduct({ IBKernel::IB_4, IBKernel::IB_4, IBKernel::IB_4, IBKernel::IB_4 })),
            "overlong sequence query must return false");

    const IBKernelTensorProduct linear_constant("DISCONTINUOUS_LINEAR");
    const IBKernelTensorProduct explicit_linear_constant({ IBKernel::BSPLINE_2, IBKernel::BSPLINE_1 });
    require(linear_constant == explicit_linear_constant, "discontinuous-linear factor order");
    LEInteractor::interpolate(string_value, NDIM, X, NDIM, field, patch, box, linear_constant);
    typed_spread->fillAll(0.0);
    LEInteractor::spread(typed_spread, force, NDIM, X, NDIM, patch, box, linear_constant);
    for (int component = 0; component < NDIM; ++component)
    {
        double expected_value = 0.0;
        for (SideIterator<NDIM> i(field->getGhostBox(), component); i; i++)
        {
            double weight = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                const double grid_x = x_lower[d] + (i()(d) - box.lower()(d) + (d == component ? 0.0 : 0.5)) * dx[d];
                const double distance = std::abs((X[d] - grid_x) / dx[d]);
                weight *= d == component ? std::max(0.0, 1.0 - distance) : (distance < 0.5 ? 1.0 : 0.0);
            }
            expected_value += weight * (*field)(i());
            require(std::abs((*typed_spread)(i()) * cell_volume - force[component] * weight) < 1.0e-12,
                    "discontinuous-linear spreading orientation");
        }
        require(std::abs(string_value[component] - expected_value) < 1.0e-12,
                "discontinuous-linear interpolation orientation");
    }

    const std::vector<IBKernel> factors
    {
        IBKernel::IB_4, IBKernel::IB_3
#if (NDIM == 3)
            ,
            IBKernel::IB_6
#endif
    };
    const IBKernelTensorProduct ordered(factors), shorthand({ IBKernel::IB_4, IBKernel::IB_3 });
    for (unsigned int component = 0; component < NDIM; ++component)
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const unsigned int slot = axis == component ? 0 : axis < component ? axis + 1 : axis;
            require(LEInteractor::get_kernel_factor(ordered, axis, component) == factors[slot],
                    "full directional mapping");
            require(LEInteractor::get_kernel_factor(shorthand, axis, component) ==
                        IBKernel(component == axis ? IBKernel::IB_4 : IBKernel::IB_3),
                    "two-factor directional mapping");
        }
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const std::string input_file = argc > 1 ? argv[1] : "";
    const bool expect_error = input_file.find("unsupported") != std::string::npos;

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "le_interactor_kernels.log");
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> load_balancer =
        new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
        new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);
    gridding_algorithm->makeCoarsestLevel(hierarchy, 0.0);

    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    PatchLevel<NDIM>::Iterator p(level);
    require(static_cast<bool>(p), "test hierarchy contains no patches");
    check_patch(level->getPatch(p()), expect_error);

    std::ofstream output("output");
    output << "LEInteractor typed/string interpolation and spreading: PASS\n"
           << "Composite alias orientation and directional mapping: PASS\n"
           << "Kernel capability queries: PASS\n";
    return 0;
}
