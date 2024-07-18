// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
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
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <tbox/Array.h>

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
    }
    return val;
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
        // prevent a warning about timer initialization
        TimerManager::createManager(nullptr);

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        auto app_initializer = make_samrai_shared<AppInitializer>(argc, argv, "cc_poisson.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        auto grid_geometry = make_samrai_shared<CartesianGridGeometryNd>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        auto patch_hierarchy = make_samrai_shared<PatchHierarchyNd>("PatchHierarchy", grid_geometry);
        auto error_detector = make_samrai_shared<StandardTagAndInitializeNd>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        auto box_generator = make_samrai_shared<BergerRigoutsosNd>();
        auto load_balancer =
            make_samrai_shared<LoadBalancerNd>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        auto gridding_algorithm =
            make_samrai_shared<GriddingAlgorithmNd>("GriddingAlgorithm",
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
        TBOX_ASSERT(var_centering == "CELL" || var_centering == "SIDE");
        TBOX_ASSERT(extrap_type == "LINEAR" || extrap_type == "QUADRATIC");

        using fcn_type = double (*)(const double* const);
        std::map<std::string, fcn_type> fcn_map;
        fcn_map["LINEAR"] = &linear_f;
        fcn_map["QUADRATIC"] = &quadratic_f;

        // Create cell-centered data and extrapolate that data at physical
        // boundaries to obtain ghost cell values.
        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
        SAMRAIPointer<VariableContext> context = var_db->getContext("CONTEXT");
        SAMRAIPointer<CellVariableNd<double> > c_var = make_samrai_shared<CellVariableNd<double> >("c_v");
        SAMRAIPointer<SideVariableNd<double> > s_var = make_samrai_shared<SideVariableNd<double> >("s_v");
        const int gcw = 1;
        const int c_idx = var_db->registerVariableAndContext(c_var, context, gcw);
        const int s_idx = var_db->registerVariableAndContext(s_var, context, gcw);
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(c_idx);
            level->allocatePatchData(s_idx);
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                SAMRAIPointer<PatchNd> patch = level->getPatch(p());
                const BoxNd& patch_box = patch->getBox();
                const hier::IndexNd& patch_lower = patch_box.lower();

                pout << "checking robin bc handling . . .\n";

                SAMRAIPointer<CartesianPatchGeometryNd> pgeom = patch->getPatchGeometry();
                const double* const x_lower = pgeom->getXLower();
                const double* const dx = pgeom->getDx();
                if (var_centering == "CELL")
                {
                    SAMRAIPointer<CellDataNd<double> > data = patch->getPatchData(c_idx);
                    for (CellIteratorNd ci(patch_box); ci; ci++)
                    {
                        const CellIndexNd& i = ci();
                        double X[NDIM];
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            X[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                        }
                        (*data)(i) = fcn_map[extrap_type](X);
                    }
                }
                else if (var_centering == "SIDE")
                {
                    SAMRAIPointer<SideDataNd<double> > data = patch->getPatchData(s_idx);
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        for (SideIteratorNd si(patch_box, axis); si; si++)
                        {
                            const SideIndexNd& i = si();
                            double X[NDIM];
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                X[d] = x_lower[d] +
                                       dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + (axis == d ? 0.0 : 0.5));
                            }
                            (*data)(i) = fcn_map[extrap_type](X);
                        }
                    }
                }
                else
                {
                    TBOX_ERROR("UNKNOWN DATA CENTERING: " << var_centering << "\n");
                }

                std::vector<RobinBcCoefStrategyNd*> bc_coefs;
                if (var_centering == "CELL")
                {
                    bc_coefs.push_back(new muParserRobinBcCoefs(
                        "bc_coefs", app_initializer->getComponentDatabase("bc_coefs_0"), grid_geometry));
                }
                else if (var_centering == "SIDE")
                {
                    for (int d = 0; d < NDIM; ++d)
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
                    SAMRAIPointer<CellDataNd<double> > data = patch->getPatchData(c_idx);
                    bc_fill_op.setPhysicalBoundaryConditions(*patch, 0.0, data->getGhostCellWidth());

                    warning = false;
                    std::vector<BoxNd> ghost_boxes;
                    if (extrap_type == "LINEAR")
                    {
                        ghost_boxes.push_back(data->getGhostBox());
                    }
                    else if (extrap_type == "QUADRATIC")
                    {
                        ghost_boxes.push_back(patch->getBox());
                        const tbox::Array<BoundaryBoxNd> codim1_boxes =
                            PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
                        if (codim1_boxes.size() != 0)
                        {
                            for (int n = 0; n < codim1_boxes.size(); ++n)
                            {
                                const BoundaryBoxNd& bdry_box = codim1_boxes[n];
                                ghost_boxes.push_back(pgeom->getBoundaryFillBox(bdry_box, patch->getBox(), gcw));
                            }
                        }
                    }
                    for (const auto& box : ghost_boxes)
                    {
                        for (CellIteratorNd ci(box); ci; ci++)
                        {
                            const CellIndexNd& i = ci();
                            double X[NDIM];
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                X[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                            }
                            double val = fcn_map[extrap_type](X);
                            if (!IBTK::rel_equal_eps(val, (*data)(i)))
                            {
                                warning = true;
                                pout << "warning: value at location " << i << " is not correct\n";
                                pout << "  expected value = " << val << "   computed value = " << (*data)(i) << "\n";
                            }
                        }
                    }
                }
                else if (var_centering == "SIDE")
                {
                    CartSideRobinPhysBdryOp bc_fill_op(s_idx, bc_coefs, false, extrap_type);
                    SAMRAIPointer<SideDataNd<double> > data = patch->getPatchData(s_idx);
                    bc_fill_op.setPhysicalBoundaryConditions(*patch, 0.0, data->getGhostCellWidth());

                    warning = false;
                    std::vector<BoxNd> ghost_boxes;
                    if (extrap_type == "LINEAR")
                    {
                        ghost_boxes.push_back(data->getGhostBox());
                    }
                    else if (extrap_type == "QUADRATIC")
                    {
                        ghost_boxes.push_back(patch->getBox());
                        const tbox::Array<BoundaryBoxNd> codim1_boxes =
                            PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
                        if (codim1_boxes.size() != 0)
                        {
                            for (int n = 0; n < codim1_boxes.size(); ++n)
                            {
                                const BoundaryBoxNd& bdry_box = codim1_boxes[n];
                                ghost_boxes.push_back(pgeom->getBoundaryFillBox(bdry_box, patch->getBox(), gcw));
                            }
                        }
                    }
                    for (const auto& box : ghost_boxes)
                    {
                        for (unsigned int axis = 0; axis < NDIM; ++axis)
                        {
                            for (SideIteratorNd si(box, axis); si; si++)
                            {
                                const SideIndexNd& i = si();
                                double X[NDIM];
                                for (unsigned int d = 0; d < NDIM; ++d)
                                {
                                    X[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) +
                                                                 (axis == d ? 0.0 : 0.5));
                                }
                                double val = fcn_map[extrap_type](X);

                                if (!IBTK::rel_equal_eps(val, (*data)(i)))
                                {
                                    warning = true;
                                    pout << "warning: value at location " << i << " is not correct\n";
                                    pout << "  expected value = " << val << "   computed value = " << (*data)(i)
                                         << "\n";
                                }
                            }
                        }
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
