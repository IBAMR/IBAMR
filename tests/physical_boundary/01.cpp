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

// Config files

#include <SAMRAI_config.h>

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <SideGeometry.h>

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <utility>
#include <vector>

#include "../tests.h"

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

double
linear_f(const double* const x)
{
    double val = 0.0;
    for (int d = 0; d < NDIM; ++d)
    {
        val += static_cast<double>(d + 1) * x[d];
    }
    return val;
}

double
quadratic_f(const double* const x)
{
    double val = 0.0;
    for (int d = 0; d < NDIM; ++d)
    {
        val += static_cast<double>(d + 1) * x[d] * x[d];
        for (int e = d + 1; e < NDIM; ++e)
        {
            val += x[d] * x[e];
        }
    }
    return val;
}

// Return the arrays of the components of the data and the boxes of their values in the patch.
std::vector<std::pair<ArrayData<NDIM, double>*, Box<NDIM>>>
component_arrays(CellData<NDIM, double>& data, const Box<NDIM>& patch_box)
{
    return { { &data.getArrayData(), patch_box } };
}

std::vector<std::pair<ArrayData<NDIM, double>*, Box<NDIM>>>
component_arrays(SideData<NDIM, double>& data, const Box<NDIM>& patch_box)
{
    std::vector<std::pair<ArrayData<NDIM, double>*, Box<NDIM>>> arrays;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        arrays.emplace_back(&data.getArrayData(axis), SideGeometry<NDIM>::toSideBox(patch_box, axis));
    }
    return arrays;
}

// Check that the accumulation is the transpose of the fill: (f, E u) = (E^T f, u), in which E is the fill with
// homogeneous conditions, u vanishes outside the patch and at Dirichlet boundary faces, and f is nonzero everywhere.
template <class DataType>
void
check_transpose(Patch<NDIM>& patch, const int u_idx, const int f_idx, RobinPhysBdryPatchStrategy& homogeneous_bdry_op)
{
    Pointer<DataType> u_data = patch.getPatchData(u_idx);
    Pointer<DataType> f_data = patch.getPatchData(f_idx);
    const auto u_arrays = component_arrays(*u_data, patch.getBox());
    const auto f_arrays = component_arrays(*f_data, patch.getBox());
    // Call fcn(u, f, i, d, phase, is_interior) for each value of the data.
    const auto for_each_value = [&](const auto& fcn)
    {
        for (std::size_t c = 0; c < u_arrays.size(); ++c)
        {
            for (Box<NDIM>::Iterator b(u_arrays[c].first->getBox()); b; b++)
            {
                for (int d = 0; d < u_data->getDepth(); ++d)
                {
                    double phase = 0.5 * c + 0.3 * d;
                    for (unsigned int k = 0; k < NDIM; ++k)
                    {
                        phase += (0.7 + 0.6 * k) * b()(k);
                    }
                    fcn(*u_arrays[c].first, *f_arrays[c].first, b(), d, phase, u_arrays[c].second.contains(b()));
                }
            }
        }
    };

    // A first fill zeroes the values at Dirichlet boundary faces. The second fill uses only the interior values.
    for_each_value(
        [](ArrayData<NDIM, double>& u,
           ArrayData<NDIM, double>& f,
           const hier::Index<NDIM>& i,
           const int d,
           const double phase,
           const bool is_interior)
        {
            u(i, d) = is_interior ? std::sin(phase) : 0.0;
            f(i, d) = std::cos(phase);
        });
    homogeneous_bdry_op.setPatchDataIndex(u_idx);
    homogeneous_bdry_op.setPhysicalBoundaryConditions(patch, 0.0, u_data->getGhostCellWidth());
    for_each_value(
        [](ArrayData<NDIM, double>& u,
           ArrayData<NDIM, double>& /*f*/,
           const hier::Index<NDIM>& i,
           const int d,
           const double /*phase*/,
           const bool is_interior)
        {
            if (!is_interior) u(i, d) = 0.0;
        });
    homogeneous_bdry_op.setPhysicalBoundaryConditions(patch, 0.0, u_data->getGhostCellWidth());
    double f_dot_Eu = 0.0;
    for_each_value([&](ArrayData<NDIM, double>& u,
                       ArrayData<NDIM, double>& f,
                       const hier::Index<NDIM>& i,
                       const int d,
                       const double /*phase*/,
                       const bool /*is_interior*/) { f_dot_Eu += f(i, d) * u(i, d); });
    homogeneous_bdry_op.setPatchDataIndex(f_idx);
    homogeneous_bdry_op.accumulateFromPhysicalBoundaryData(patch, 0.0, f_data->getGhostCellWidth());
    double ETf_dot_u = 0.0;
    for_each_value(
        [&](ArrayData<NDIM, double>& u,
            ArrayData<NDIM, double>& f,
            const hier::Index<NDIM>& i,
            const int d,
            const double /*phase*/,
            const bool is_interior)
        {
            if (is_interior) ETf_dot_u += f(i, d) * u(i, d);
        });

    if (!std::isfinite(f_dot_Eu) || !(std::abs(f_dot_Eu - ETf_dot_u) <= 1.0e-12 * std::abs(f_dot_Eu)))
    {
        TBOX_ERROR("(f, E u) = " << f_dot_Eu << " and (E^T f, u) = " << ETf_dot_u << " differ\n");
    }
    pout << std::setprecision(10) << std::fixed << "(f, E u)   = " << f_dot_Eu << '\n'
         << "(E^T f, u) = " << ETf_dot_u << '\n';
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

    // Do not print file names or line numbers to the output file.
    Logger::getInstance()->setAbortAppender(Pointer<Logger::Appender>(new TestAppender()));

    { // cleanup dynamically allocated objects prior to shutdown
        // prevent a warning about timer initialization
        TimerManager::createManager(nullptr);

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

        std::string var_centering = input_db->getString("VAR_CENTERING");
        std::string extrap_type = input_db->getString("EXTRAP_TYPE");
        const int depth = input_db->getIntegerWithDefault("DEPTH", 1);
        // The data are a LINEAR or QUADRATIC function, of the degree of the extrapolation by default.
        const std::string function_type = input_db->getStringWithDefault("FUNCTION_TYPE", extrap_type);
        const int gcw = input_db->getIntegerWithDefault("GHOST_WIDTH", 1);
        const bool check_transpose_of_fill = input_db->getBoolWithDefault("CHECK_TRANSPOSE", false);
        TBOX_ASSERT(var_centering == "CELL" || var_centering == "SIDE");
        TBOX_ASSERT(extrap_type == "LINEAR" || extrap_type == "QUADRATIC");
        TBOX_ASSERT(function_type == "LINEAR" || function_type == "QUADRATIC");

        using fcn_type = double (*)(const double* const);
        std::map<std::string, fcn_type> fcn_map;
        fcn_map["LINEAR"] = &linear_f;
        fcn_map["QUADRATIC"] = &quadratic_f;

        // Create cell-centered data and extrapolate that data at physical
        // boundaries to obtain ghost cell values.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = var_db->getContext("CONTEXT");
        Pointer<CellVariable<NDIM, double>> c_var = new CellVariable<NDIM, double>("c_v");
        Pointer<SideVariable<NDIM, double>> s_var = new SideVariable<NDIM, double>("s_v", depth);
        const int c_idx = var_db->registerVariableAndContext(c_var, context, gcw);
        const int s_idx = var_db->registerVariableAndContext(s_var, context, gcw);
        const int c_f_idx = var_db->registerClonedPatchDataIndex(c_var, c_idx);
        const int s_f_idx = var_db->registerClonedPatchDataIndex(s_var, s_idx);
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(c_idx);
            level->allocatePatchData(s_idx);
            level->allocatePatchData(c_f_idx);
            level->allocatePatchData(s_f_idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const hier::Index<NDIM>& patch_lower = patch_box.lower();

                pout << "checking robin bc handling . . .\n";

                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                const double* const x_lower = pgeom->getXLower();
                const double* const dx = pgeom->getDx();
                if (var_centering == "CELL")
                {
                    Pointer<CellData<NDIM, double>> data = patch->getPatchData(c_idx);
                    for (CellIterator<NDIM> ci(patch_box); ci; ci++)
                    {
                        const CellIndex<NDIM>& i = ci();
                        double X[NDIM];
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            X[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                        }
                        (*data)(i) = fcn_map[function_type](X);
                    }
                }
                else if (var_centering == "SIDE")
                {
                    Pointer<SideData<NDIM, double>> data = patch->getPatchData(s_idx);
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        for (SideIterator<NDIM> si(patch_box, axis); si; si++)
                        {
                            const SideIndex<NDIM>& i = si();
                            double X[NDIM];
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                X[d] = x_lower[d] +
                                       dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + (axis == d ? 0.0 : 0.5));
                            }
                            for (int comp = 0; comp < depth; ++comp)
                            {
                                (*data)(i, comp) = static_cast<double>(comp + 1) * fcn_map[function_type](X);
                            }
                        }
                    }
                }
                else
                {
                    TBOX_ERROR("UNKNOWN DATA CENTERING: " << var_centering << "\n");
                }

                std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs;
                if (var_centering == "CELL")
                {
                    bc_coefs.push_back(new muParserRobinBcCoefs(
                        "bc_coefs", app_initializer->getComponentDatabase("bc_coefs_0"), grid_geometry));
                }
                else if (var_centering == "SIDE")
                {
                    // Each depth is a separate MAC vector field with one coefficient object per
                    // staggered component, ordered as [NDIM * depth_index + axis].
                    for (int d = 0; d < NDIM * depth; ++d)
                    {
                        std::string bc_name = "bc_coefs_" + std::to_string(d);
                        bc_coefs.push_back(new muParserRobinBcCoefs(
                            bc_name, app_initializer->getComponentDatabase(bc_name), grid_geometry));
                    }
                }
                bool warning = false;
                if (var_centering == "CELL")
                {
                    CartCellRobinPhysBdryOp bc_fill_op(c_idx, bc_coefs, false, extrap_type);
                    Pointer<CellData<NDIM, double>> data = patch->getPatchData(c_idx);
                    bc_fill_op.setPhysicalBoundaryConditions(*patch, 0.0, data->getGhostCellWidth());

                    warning = false;
                    for (CellIterator<NDIM> ci(data->getGhostBox()); ci; ci++)
                    {
                        const CellIndex<NDIM>& i = ci();
                        double X[NDIM];
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            X[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                        }
                        double val = fcn_map[function_type](X);
                        if (!IBTK::rel_equal_eps(val, (*data)(i)))
                        {
                            warning = true;
                            pout << "warning: value at location " << i << " is not correct\n";
                            pout << "  expected value = " << val << "   computed value = " << (*data)(i) << "\n";
                        }
                    }

                    if (check_transpose_of_fill)
                    {
                        CartCellRobinPhysBdryOp homogeneous_bc_fill_op(c_idx, bc_coefs, true, extrap_type);
                        check_transpose<CellData<NDIM, double>>(*patch, c_idx, c_f_idx, homogeneous_bc_fill_op);
                    }
                }
                else if (var_centering == "SIDE")
                {
                    CartSideRobinPhysBdryOp bc_fill_op(s_idx, bc_coefs, false, extrap_type);
                    Pointer<SideData<NDIM, double>> data = patch->getPatchData(s_idx);
                    bc_fill_op.setPhysicalBoundaryConditions(*patch, 0.0, data->getGhostCellWidth());

                    warning = false;
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        for (SideIterator<NDIM> si(data->getGhostBox(), axis); si; si++)
                        {
                            const SideIndex<NDIM>& i = si();
                            double X[NDIM];
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                X[d] = x_lower[d] +
                                       dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + (axis == d ? 0.0 : 0.5));
                            }
                            for (int comp = 0; comp < depth; ++comp)
                            {
                                const double val = static_cast<double>(comp + 1) * fcn_map[function_type](X);

                                if (!IBTK::rel_equal_eps(val, (*data)(i, comp)))
                                {
                                    warning = true;
                                    pout << "warning: value at location " << i << " depth " << comp
                                         << " is not correct\n";
                                    pout << "  expected value = " << val << "   computed value = " << (*data)(i, comp)
                                         << "\n";
                                }
                            }
                        }
                    }

                    if (check_transpose_of_fill)
                    {
                        CartSideRobinPhysBdryOp homogeneous_bc_fill_op(s_idx, bc_coefs, true, extrap_type);
                        check_transpose<SideData<NDIM, double>>(*patch, s_idx, s_f_idx, homogeneous_bc_fill_op);
                    }
                }

                if (!warning)
                {
                    pout << "boundary data appears to be correct.\n";
                }
                else
                {
                    pout << "possible errors encountered in extrapolated boundary data.\n";
                }
            }
        }

    } // cleanup dynamically allocated objects prior to shutdown

    return 0;
} // main
