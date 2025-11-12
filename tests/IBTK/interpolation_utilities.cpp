// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2025 by the IBAMR developers
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

        // Create cell-centered data and extrapolate that data at physical
        // boundaries to obtain ghost cell values.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = var_db->getContext("CONTEXT");
        Pointer<CellVariable<NDIM, double> > cc_var = new CellVariable<NDIM, double>("cc");
        Pointer<SideVariable<NDIM, double> > sc_var = new SideVariable<NDIM, double>("sc");
        Pointer<NodeVariable<NDIM, double> > nc_var = new NodeVariable<NDIM, double>("nc", NDIM);
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
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(cc_idx);
            level->allocatePatchData(nc_idx);
            level->allocatePatchData(sc_idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM, double> > cc_data = patch->getPatchData(cc_idx);
                Pointer<SideData<NDIM, double> > sc_data = patch->getPatchData(sc_idx);
                Pointer<NodeData<NDIM, double> > nc_data = patch->getPatchData(nc_idx);
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

        // Now interpolate to the specified point.
        std::vector<VectorNd> x_pt(2);
        for (int d = 0; d < NDIM; ++d) x_pt[0][d] = 0.7;
        for (int d = 0; d < NDIM; ++d) x_pt[1][d] = 0.2;

        // Cell centered
        pout << "Interpolating cell centered values\n";
        std::vector<double> interped_val = interpolate(x_pt, cc_idx, cc_var, 1, patch_hierarchy, "IB_4");
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
