// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2026 by the IBAMR developers
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
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/interpolation_utilities.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include <stdexcept>

#include <ibtk/app_namespaces.h>

namespace
{
void
require_kernel_check(bool valid, const char* message)
{
    if (!valid) throw std::runtime_error(message);
}

template <class F>
void
require_unsupported(F operation)
{
    bool rejected = false;
    try
    {
        operation();
    }
    catch (const std::invalid_argument& error)
    {
        rejected = std::string(error.what()) == "LEInteractor: unsupported kernel description";
    }
    require_kernel_check(rejected, "unsupported description did not produce the intended rejection");
}

void
check_kernel_consumption(Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    int patches = 0;
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        ++patches;
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const auto& box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
        const auto dx = geometry->getDx();
        const auto xlow = geometry->getXLower();
        Pointer<SideData<NDIM, double>> field = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(5));
        Pointer<SideData<NDIM, double>> legacy = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(5));
        Pointer<SideData<NDIM, double>> typed = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(5));
        std::vector<double> X(NDIM), force(NDIM), a(NDIM), b(NDIM);
        double volume = 1.0;
        for (int d = 0; d < NDIM; ++d)
        {
            X[d] = xlow[d] + (3.27 + 0.17 * d) * dx[d];
            force[d] = 1.3 + 0.7 * d;
            volume *= dx[d];
        }
        for (int component = 0; component < NDIM; ++component)
            for (SideIterator<NDIM> i(field->getGhostBox(), component); i; i++)
            {
                double value = 0.4 + component;
                for (int d = 0; d < NDIM; ++d) value += std::sin((0.13 + 0.09 * d) * (i()(d) + 0.3 * component));
                (*field)(i()) = value;
            }
        const char* names[] = { "PIECEWISE_CONSTANT",
                                "PIECEWISE_LINEAR",
                                "PIECEWISE_CUBIC",
                                "IB_3",
                                "IB_4",
                                "IB_4_W8",
                                "IB_5",
                                "IB_6",
                                "BSPLINE_3",
                                "BSPLINE_4",
                                "BSPLINE_5",
                                "BSPLINE_6",
                                "USER_DEFINED",
                                "DISCONTINUOUS_LINEAR",
                                "COMPOSITE_BSPLINE_23",
                                "COMPOSITE_BSPLINE_32",
                                "COMPOSITE_BSPLINE_34",
                                "COMPOSITE_BSPLINE_43",
                                "COMPOSITE_BSPLINE_45",
                                "COMPOSITE_BSPLINE_54",
                                "COMPOSITE_BSPLINE_56",
                                "COMPOSITE_BSPLINE_65" };
        for (const char* name : names)
        {
            const auto kernel = IBKernelTensorProduct::from_name(name);
            require_kernel_check(LEInteractor::isKnownKernel(kernel), "legacy kernel unsupported");
            require_kernel_check(LEInteractor::getStencilSize(kernel) == LEInteractor::getStencilSize(name),
                                 "stencil query changed");
            LEInteractor::interpolate(a, NDIM, X, NDIM, field, patch, box, name);
            LEInteractor::interpolate(b, NDIM, X, NDIM, field, patch, box, kernel);
            require_kernel_check(a == b, "typed/string interpolation differs");
            legacy->fillAll(0.0);
            typed->fillAll(0.0);
            LEInteractor::spread(legacy, force, NDIM, X, NDIM, patch, box, name);
            LEInteractor::spread(typed, force, NDIM, X, NDIM, patch, box, kernel);
            for (int component = 0; component < NDIM; ++component)
                for (SideIterator<NDIM> i(field->getGhostBox(), component); i; i++)
                    require_kernel_check((*legacy)(i()) == (*typed)(i()), "typed/string spreading differs");
        }

        // Independent hat/nearest-grid weights establish the 21 alias orientation.
        const IBKernelTensorProduct linear_constant(IBKernel::BSPLINE_2, IBKernel::BSPLINE_1);
        LEInteractor::interpolate(a, NDIM, X, NDIM, field, patch, box, "DISCONTINUOUS_LINEAR");
        typed->fillAll(0.0);
        LEInteractor::spread(typed, force, NDIM, X, NDIM, patch, box, linear_constant);
        for (int component = 0; component < NDIM; ++component)
        {
            double expected = 0.0;
            for (SideIterator<NDIM> i(field->getGhostBox(), component); i; i++)
            {
                double weight = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    const double grid_x = xlow[d] + (i()(d) - box.lower()(d) + (d == component ? 0.0 : 0.5)) * dx[d];
                    const double distance = std::abs((X[d] - grid_x) / dx[d]);
                    weight *= d == component ? std::max(0.0, 1.0 - distance) : (distance < 0.5 ? 1.0 : 0.0);
                }
                expected += weight * (*field)(i());
                require_kernel_check(std::abs((*typed)(i()) * volume - force[component] * weight) < 1.e-12,
                                     "21 spreading orientation");
            }
            require_kernel_check(std::abs(a[component] - expected) < 1.e-12, "21 interpolation orientation");
        }

        const std::array<IBKernel, NDIM> axes{ { IBKernel::IB_4,
                                                 IBKernel::IB_3
#if (NDIM == 3)
                                                 ,
                                                 IBKernel::IB_6
#endif
        } };
        const IBKernelTensorProduct unsupported[] = { IBKernelTensorProduct(IBKernel("UnregisteredTrialKernel")),
                                                      IBKernelTensorProduct(axes),
                                                      IBKernelTensorProduct(IBKernel::BSPLINE_1, IBKernel::BSPLINE_2) };
        for (const auto& kernel : unsupported)
        {
            a.assign(NDIM, 17.0);
            typed->fillAll(19.0);
            require_unsupported([&] { LEInteractor::interpolate(a, NDIM, X, NDIM, field, patch, box, kernel); });
            require_unsupported([&] { LEInteractor::spread(typed, force, NDIM, X, NDIM, patch, box, kernel); });
            require_kernel_check(a == std::vector<double>(NDIM, 17.0), "rejected interpolation changed output");
            for (int component = 0; component < NDIM; ++component)
                for (SideIterator<NDIM> i(typed->getGhostBox(), component); i; i++)
                    require_kernel_check((*typed)(i()) == 19.0, "rejected spreading changed output");
        }
    }
    require_kernel_check(IBTK_MPI::sumReduction(patches) > 1, "trial needs multiple patches");
}
} // namespace

double
exact_fcn(const VectorNd& x)
{
    double ret = 1.0;
    for (int d = 0; d < NDIM; ++d) ret += x[d] * static_cast<double>(d + 1);
    return ret;
}

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_poisson.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create cell-centered data and extrapolate that data at physical
        // boundaries to obtain ghost cell values.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = var_db->getContext("CONTEXT");
        Pointer<CellVariable<NDIM, double>> cc_var = new CellVariable<NDIM, double>("cc");
        Pointer<SideVariable<NDIM, double>> sc_var = new SideVariable<NDIM, double>("sc");
        Pointer<NodeVariable<NDIM, double>> nc_var = new NodeVariable<NDIM, double>("nc", NDIM);
        const int gcw = 4;
        const int cc_idx = var_db->registerVariableAndContext(cc_var, context, gcw);
        const int sc_idx = var_db->registerVariableAndContext(sc_var, context, gcw);
        const int nc_idx = var_db->registerVariableAndContext(nc_var, context, gcw);

        // Initialize the AMR patch hierarchy.
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

        // Allocate and fill in patch data
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(cc_idx);
            level->allocatePatchData(nc_idx);
            level->allocatePatchData(sc_idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> cc_data = patch->getPatchData(cc_idx);
                Pointer<SideData<NDIM, double>> sc_data = patch->getPatchData(sc_idx);
                Pointer<NodeData<NDIM, double>> nc_data = patch->getPatchData(nc_idx);
                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double* const xlow = pgeom->getXLower();
                const hier::Index<NDIM>& idx_low = patch->getBox().lower();
                for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    VectorNd x;
                    for (int d = 0; d < NDIM; ++d)
                        x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                    (*cc_data)(idx) = exact_fcn(x);
                }

                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                    {
                        const SideIndex<NDIM>& idx = si();
                        VectorNd x;
                        for (int d = 0; d < NDIM; ++d)
                            x[d] =
                                xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
                        (*sc_data)(idx) = exact_fcn(x);
                    }
                }

                for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
                {
                    const NodeIndex<NDIM>& idx = ni();
                    VectorNd x;
                    for (int d = 0; d < NDIM; ++d) x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)));
                    for (int d = 0; d < NDIM; ++d) (*nc_data)(idx, d) = exact_fcn(x);
                }
            }
        }

        // Now fill ghost cells
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_comps{ ITC(cc_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE"),
                                           ITC(sc_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE"),
                                           ITC(nc_idx, "LINEAR_REFINE", false, "NONE") };
        HierarchyGhostCellInterpolation ghost_cell_fill;
        ghost_cell_fill.initializeOperatorState(ghost_cell_comps, patch_hierarchy, coarsest_ln, finest_ln);
        ghost_cell_fill.fillData(0.0);

        if (input_db->getBoolWithDefault("kernel_trial", false))
        {
            check_kernel_consumption(patch_hierarchy);
            pout << "LEInteractor typed/string execution, alias orientation and rejection: PASS\n";
        }

        // Now interpolate to the specified point.
        std::vector<VectorNd> x_pt(2);
        for (int d = 0; d < NDIM; ++d) x_pt[0][d] = 0.7;
        for (int d = 0; d < NDIM; ++d) x_pt[1][d] = 0.2;

        // Cell centered
        pout << "Interpolating cell centered values\n";
        std::vector<double> interped_val = interpolate(x_pt, cc_idx, cc_var, 1, patch_hierarchy, "IB_4");
        const IBKernelTensorProduct resolved(IBKernel::IB_4);
        require_kernel_check(interped_val == interpolate(x_pt, cc_idx, cc_var, 1, patch_hierarchy, resolved),
                             "cell utility typed route");
        require_kernel_check(interpolate(x_pt[0], cc_idx, cc_var, 1, patch_hierarchy, resolved)[0] == interped_val[0],
                             "single-point utility route");
        for (int i = 0; i < 2; ++i)
        {
            bool correct = std::abs(interped_val[i] - exact_fcn(x_pt[i])) < 1.0e-12;
            correct = IBTK_MPI::maxReduction(correct ? 0 : 1) == 0;
            if (!correct)
            {
                plog << "Interpolant number " << i << " was not exact!\n";
                plog << "Expected " << exact_fcn(x_pt[i]) << " and got " << interped_val[i] << "\n";
                plog << "Error: " << interped_val[i] - exact_fcn(x_pt[i]) << "\n";
            }
        }

        // Side centered
        pout << "Interpolating side centered values\n";
        interped_val = interpolate(x_pt, sc_idx, sc_var, 1, patch_hierarchy, "IB_4");
        require_kernel_check(interped_val == interpolate(x_pt, sc_idx, sc_var, 1, patch_hierarchy, resolved),
                             "side utility typed route");
        for (int i = 0; i < 2; ++i)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                bool correct = std::abs(interped_val[i * NDIM + d] - exact_fcn(x_pt[i])) < 1.0e-12;
                correct = IBTK_MPI::maxReduction(correct ? 0 : 1) == 0;
                if (!correct)
                {
                    plog << "Interpolant number " << i << " and depth " << d << " was not exact!\n";
                    plog << "Expected " << exact_fcn(x_pt[i]) << " and got " << interped_val[i] << "\n";
                    plog << "Error: " << interped_val[i] - exact_fcn(x_pt[i]) << "\n";
                }
            }
        }

        // Node centered
        pout << "Interpolating node centered values\n";
        interped_val = interpolate(x_pt, nc_idx, nc_var, NDIM, patch_hierarchy, "IB_4");
        require_kernel_check(interped_val == interpolate(x_pt, nc_idx, nc_var, NDIM, patch_hierarchy, resolved),
                             "node utility typed route");
        for (int i = 0; i < 2; ++i)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                bool correct = std::abs(interped_val[i * NDIM + d] - exact_fcn(x_pt[i])) < 1.0e-12;
                correct = IBTK_MPI::maxReduction(correct ? 0 : 1) == 0;
                if (!correct)
                {
                    plog << "Interpolant number " << i << " and depth " << d << " was not exact!\n";
                    plog << "Expected " << exact_fcn(x_pt[i]) << " and got " << interped_val[i] << "\n";
                    plog << "Error: " << interped_val[i] - exact_fcn(x_pt[i]) << "\n";
                }
            }
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main
