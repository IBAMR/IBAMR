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

#include <ibtk/AppInitializer.h>
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
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

#include <ibtk/app_namespaces.h>

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

        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();

        // Create cell-centered data and extrapolate that data at physical
        // boundaries to obtain ghost cell values.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = var_db->getContext("CONTEXT");
        Pointer<CellVariable<NDIM, double> > var = new CellVariable<NDIM, double>("v");
        const int gcw = 4;
        const int data_idx = var_db->registerVariableAndContext(var, context, gcw);
        visit_data_writer->registerPlotQuantity("Q", "SCALAR", data_idx);

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
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(data_idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM, double> > data = patch->getPatchData(data_idx);
                Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double* const xlow = pgeom->getXLower();
                const hier::Index<NDIM>& idx_low = patch->getBox().lower();
                for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    VectorNd x;
                    for (int d = 0; d < NDIM; ++d)
                        x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                    (*data)(idx) = exact_fcn(x);
                }
            }
        }

        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

        // Now fill ghost cells
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        ITC ghost_cell_comp(data_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE");
        HierarchyGhostCellInterpolation ghost_cell_fill;
        ghost_cell_fill.initializeOperatorState(ghost_cell_comp, patch_hierarchy, coarsest_ln, finest_ln);
        ghost_cell_fill.fillData(0.0);

        // Now interpolate to the specified point.
        VectorNd x_pt;
        for (int d = 0; d < NDIM; ++d) x_pt[d] = 0.1;
        double bilinear_interp = Interpolation::interpolate(x_pt, data_idx, var, patch_hierarchy, 0);
        double l2_interp = Interpolation::interpolateL2(x_pt, data_idx, var, patch_hierarchy, 1, 1);
        if (std::abs(bilinear_interp - exact_fcn(x_pt)) > 1.0e-14) plog << "Bilinear interpolant not exact!\n";
        if (std::abs(l2_interp - exact_fcn(x_pt)) > 1.0e-14) plog << "Least squares interpolant not exact!\n";
    } // cleanup dynamically allocated objects prior to shutdown
} // main
