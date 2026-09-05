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

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
std::size_t user_defined_calls = 0;

double
user_defined_kernel(double r)
{
    ++user_defined_calls;
    r = std::abs(r);
    return r < 1.0 ? 1.0 - r : 0.0;
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

    // This is a valid, extensible tensor-product description, but
    // LEInteractor does not yet provide its numerical implementation.
    const IBKernelTensorProduct unsupported("COMPOSITE_BSPLINE_12");
    TBOX_ASSERT(!LEInteractor::isKnownKernel(unsupported));
    TBOX_ASSERT(LEInteractor::isKnownKernel("IB_4"));
    TBOX_ASSERT(!LEInteractor::isKnownKernel("BAD NAME"));
    TBOX_ASSERT(!LEInteractor::isKnownKernel(static_cast<const char*>(nullptr)));
    if (expect_error)
    {
        Pointer<Logger::Appender> abort_appender = new TestAppender();
        Logger::getInstance()->setAbortAppender(abort_appender);
        PIO::logOnlyNodeZero("output");
        string_value.assign(NDIM, 17.0);
        LEInteractor::interpolate(string_value, NDIM, X, NDIM, field, patch, box, unsupported);
        return;
    }

    for (const std::string name : { "IB_4", "piecewise_constant", "discontinuous_linear", "composite_bspline_23" })
    {
        const IBKernelTensorProduct kernel(name);
        TBOX_ASSERT(LEInteractor::isKnownKernel(name) && LEInteractor::isKnownKernel(kernel));
        TBOX_ASSERT(LEInteractor::getStencilSize(name) == LEInteractor::getStencilSize(kernel));
        LEInteractor::interpolate(string_value, NDIM, X, NDIM, field, patch, box, name);
        LEInteractor::interpolate(typed_value, NDIM, X, NDIM, field, patch, box, kernel);
        TBOX_ASSERT(string_value == typed_value);
        string_spread->fillAll(0.0);
        typed_spread->fillAll(0.0);
        LEInteractor::spread(string_spread, force, NDIM, X, NDIM, patch, box, name);
        LEInteractor::spread(typed_spread, force, NDIM, X, NDIM, patch, box, kernel);
        for (int component = 0; component < NDIM; ++component)
            for (SideIterator<NDIM> i(field->getGhostBox(), component); i; i++)
                TBOX_ASSERT((*string_spread)(i()) == (*typed_spread)(i()));
    }

    for (const char* name : { "", "BAD NAME", "UNREGISTERED_KERNEL", "COMPOSITE_BSPLINE_12" })
        TBOX_ASSERT(!LEInteractor::isKnownKernel(name));
    const IBKernelTensorProduct linear_constant("DISCONTINUOUS_LINEAR");
    const IBKernelTensorProduct explicit_linear_constant({ IBKernel::BSPLINE_2, IBKernel::BSPLINE_1 });
    TBOX_ASSERT(linear_constant == explicit_linear_constant);
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
            TBOX_ASSERT(std::abs((*typed_spread)(i()) * cell_volume - force[component] * weight) < 1.0e-12);
        }
        TBOX_ASSERT(std::abs(string_value[component] - expected_value) < 1.0e-12);
    }

    const auto saved_kernel = LEInteractor::s_kernel_fcn;
    const int saved_stencil_size = LEInteractor::s_kernel_fcn_stencil_size;
    LEInteractor::s_kernel_fcn = &user_defined_kernel;
    LEInteractor::s_kernel_fcn_stencil_size = 2;
    const IBKernelTensorProduct user_defined(IBKernel("USER_DEFINED"));
    TBOX_ASSERT(LEInteractor::isKnownKernel(user_defined) && LEInteractor::isKnownKernel("USER_DEFINED"));
    TBOX_ASSERT(LEInteractor::getStencilSize(user_defined) == 2);

    user_defined_calls = 0;
    LEInteractor::interpolate(string_value, NDIM, X, NDIM, field, patch, box, "USER_DEFINED");
    LEInteractor::interpolate(typed_value, NDIM, X, NDIM, field, patch, box, user_defined);
    TBOX_ASSERT(user_defined_calls > 0 && string_value == typed_value);
    string_spread->fillAll(0.0);
    typed_spread->fillAll(0.0);
    user_defined_calls = 0;
    LEInteractor::spread(string_spread, force, NDIM, X, NDIM, patch, box, "USER_DEFINED");
    LEInteractor::spread(typed_spread, force, NDIM, X, NDIM, patch, box, user_defined);
    TBOX_ASSERT(user_defined_calls > 0);
    for (int component = 0; component < NDIM; ++component)
        for (SideIterator<NDIM> i(field->getGhostBox(), component); i; i++)
            TBOX_ASSERT((*string_spread)(i()) == (*typed_spread)(i()));

    Pointer<SideData<NDIM, double>> mask = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(5));
    Pointer<SideData<NDIM, double>> string_masked_spread = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(5));
    Pointer<SideData<NDIM, double>> typed_masked_spread = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(5));
    mask->fillAll(1.0);
    std::vector<double> string_masked_value(NDIM), typed_masked_value(NDIM);
    user_defined_calls = 0;
    LEInteractor::interpolate(string_masked_value, NDIM, X, NDIM, mask, field, patch, box, "USER_DEFINED");
    LEInteractor::interpolate(typed_masked_value, NDIM, X, NDIM, mask, field, patch, box, user_defined);
    TBOX_ASSERT(user_defined_calls > 0 && string_masked_value == typed_masked_value);
    string_masked_spread->fillAll(0.0);
    typed_masked_spread->fillAll(0.0);
    user_defined_calls = 0;
    LEInteractor::spread(mask, string_masked_spread, force, NDIM, X, NDIM, patch, box, "USER_DEFINED");
    LEInteractor::spread(mask, typed_masked_spread, force, NDIM, X, NDIM, patch, box, user_defined);
    TBOX_ASSERT(user_defined_calls > 0);
    for (int component = 0; component < NDIM; ++component)
        for (SideIterator<NDIM> i(field->getGhostBox(), component); i; i++)
            TBOX_ASSERT((*string_masked_spread)(i()) == (*typed_masked_spread)(i()));
    LEInteractor::s_kernel_fcn = saved_kernel;
    LEInteractor::s_kernel_fcn_stencil_size = saved_stencil_size;
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
    TBOX_ASSERT(static_cast<bool>(p));
    check_patch(level->getPatch(p()), expect_error);

    std::ofstream output("output");
    TBOX_ASSERT(static_cast<bool>(output));
    return 0;
}
