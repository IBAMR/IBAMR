// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2022 by the IBAMR developers
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
#include <IBTK_config.h>
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
#include <ibtk/IBTKInit.h>

#include <LocationIndexRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

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
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
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

        // Create cell-centered data and extrapolate that data at physical
        // boundaries to obtain ghost cell values.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = var_db->getContext("CONTEXT");
        Pointer<CellVariable<NDIM, double> > var = new CellVariable<NDIM, double>("v");
        const int gcw = 4;
        const int idx = var_db->registerVariableAndContext(var, context, gcw);
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const hier::Index<NDIM>& patch_lower = patch_box.lower();
                Pointer<CellData<NDIM, double> > data = patch->getPatchData(idx);
                for (Box<NDIM>::Iterator b(patch_box); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    (*data)(i) = 0;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        (*data)(i) += 4 * (d + 1) * (d + 1) * i(d);
                    }
                }

                pout << "level number = " << ln << "\n";
                pout << "patch_box = " << patch_box << "\n";
                pout << "\n";

                plog << "interior data:\n";
                data->print(data->getBox());
                plog << "\n";

                CartExtrapPhysBdryOp constant_fill_op(idx, "CONSTANT");
                constant_fill_op.setPhysicalBoundaryConditions(*patch, 0.0, data->getGhostCellWidth());
                plog << "constant extrapolated ghost data:\n";
                data->print(data->getGhostBox());
                plog << "\n";

                CartExtrapPhysBdryOp linear_fill_op(idx, "LINEAR");
                linear_fill_op.setPhysicalBoundaryConditions(*patch, 0.0, data->getGhostCellWidth());
                plog << "linear extrapolated ghost data:\n";
                data->print(data->getGhostBox());
                plog << "\n";

                bool warning = false;
                for (Box<NDIM>::Iterator b(data->getGhostBox()); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    double val = 0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        val += 4 * (d + 1) * (d + 1) * i(d);
                    }

                    if (!IBTK::rel_equal_eps(val, (*data)(i)))
                    {
                        warning = true;
                        pout << "warning: value at location " << i << " is not correct\n";
                        pout << "  expected value = " << val << "   computed value = " << (*data)(i) << "\n";
                    }
                }

                if (!warning)
                {
                    pout << "linearly extrapolated boundary data appears to be correct.\n";
                }
                else
                {
                    pout << "possible errors encountered in linearly extrapolated boundary data.\n";
                }

                pout << "checking robin bc handling . . .\n";

                Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const x_lower = pgeom->getXLower();
                const double* const x_upper = grid_geometry->getXUpper();
                const double* const dx = pgeom->getDx();
                const double shift = 3.14159;
                for (Box<NDIM>::Iterator b(patch_box); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    double X[NDIM];
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        X[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                    }
                    (*data)(i) = 2.0 * X[NDIM - 1] + shift;
                }

                plog << "interior data:\n";
                data->print(data->getBox());
                plog << "\n";

                LocationIndexRobinBcCoefs<NDIM> dirichlet_bc_coef("dirichlet_bc_coef", NULL);
                for (unsigned int d = 0; d < NDIM - 1; ++d)
                {
                    dirichlet_bc_coef.setBoundarySlope(2 * d, 0.0);
                    dirichlet_bc_coef.setBoundarySlope(2 * d + 1, 0.0);
                }
                dirichlet_bc_coef.setBoundaryValue(2 * (NDIM - 1), shift);
                dirichlet_bc_coef.setBoundaryValue(2 * (NDIM - 1) + 1, 2.0 * x_upper[NDIM - 1] + shift);

                CartCellRobinPhysBdryOp dirichlet_bc_fill_op(idx, &dirichlet_bc_coef);
                dirichlet_bc_fill_op.setPhysicalBoundaryConditions(*patch, 0.0, data->getGhostCellWidth());
                plog << "extrapolated ghost data:\n";
                data->print(data->getGhostBox());
                plog << "\n";

                warning = false;
                for (Box<NDIM>::Iterator b(data->getGhostBox()); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    double X[NDIM];
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        X[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                    }
                    double val = 2.0 * X[NDIM - 1] + shift;

                    if (!IBTK::rel_equal_eps(val, (*data)(i)))
                    {
                        warning = true;
                        pout << "warning: value at location " << i << " is not correct\n";
                        pout << "  expected value = " << val << "   computed value = " << (*data)(i) << "\n";
                    }
                }

                if (!warning)
                {
                    pout << "dirichlet boundary data appears to be correct.\n";
                }
                else
                {
                    pout << "possible errors encountered in extrapolated dirichlet boundary data.\n";
                }

                LocationIndexRobinBcCoefs<NDIM> neumann_bc_coef("neumann_bc_coef", NULL);
                for (unsigned int d = 0; d < NDIM - 1; ++d)
                {
                    neumann_bc_coef.setBoundarySlope(2 * d, 0.0);
                    neumann_bc_coef.setBoundarySlope(2 * d + 1, 0.0);
                }
                neumann_bc_coef.setBoundarySlope(2 * (NDIM - 1), -2.0);
                neumann_bc_coef.setBoundarySlope(2 * (NDIM - 1) + 1, +2.0);

                CartCellRobinPhysBdryOp neumann_bc_fill_op(idx, &neumann_bc_coef);
                neumann_bc_fill_op.setPhysicalBoundaryConditions(*patch, 0.0, data->getGhostCellWidth());
                plog << "extrapolated ghost data:\n";
                data->print(data->getGhostBox());
                plog << "\n";

                warning = false;
                for (Box<NDIM>::Iterator b(data->getGhostBox()); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    double X[NDIM];
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        X[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                    }
                    double val = 2.0 * X[NDIM - 1] + shift;

                    if (!IBTK::rel_equal_eps(val, (*data)(i)))
                    {
                        warning = true;
                        pout << "warning: value at location " << i << " is not correct\n";
                        pout << "  expected value = " << val << "   computed value = " << (*data)(i) << "\n";
                    }
                }

                if (!warning)
                {
                    pout << "neumann boundary data appears to be correct.\n";
                }
                else
                {
                    pout << "possible errors encountered in extrapolated neumann boundary data.\n";
                }
            }
        }

    } // cleanup dynamically allocated objects prior to shutdown
} // main
