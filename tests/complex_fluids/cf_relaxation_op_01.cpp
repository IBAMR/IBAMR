// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
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
#include "ibtk/muParserCartGridFunction.h"

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include "ibamr/CFGiesekusStrategy.h"
#include "ibamr/CFOldroydBStrategy.h"
#include "ibamr/CFRoliePolyStrategy.h"
#include "ibamr/ibamr_enums.h"

#include <ibtk/AppInitializer.h>
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <tbox/Array.h>

#include "ibamr/namespaces.h"

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
        // prevent a warning about timer initialization
        TimerManager::createManager(nullptr);

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_poisson.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometryNd> grid_geometry = new CartesianGridGeometryNd(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchyNd> patch_hierarchy = new PatchHierarchyNd("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitializeNd> error_detector = new StandardTagAndInitializeNd(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsosNd> box_generator = new BergerRigoutsosNd();
        Pointer<LoadBalancerNd> load_balancer =
            new LoadBalancerNd("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithmNd> gridding_algorithm =
            new GriddingAlgorithmNd("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);
        Pointer<CFStrategy> cf_op;
        std::string relax_op = input_db->getString("RELAX_OP");
        if (relax_op == "OLDROYDB")
        {
            cf_op = new CFOldroydBStrategy("OldroydB", app_initializer->getComponentDatabase("ComplexFluid"));
        }
        else if (relax_op == "GIESEKUS")
        {
            cf_op = new CFGiesekusStrategy("Giesekus", app_initializer->getComponentDatabase("ComplexFluid"));
        }
        else if (relax_op == "ROLIEPOLY")
        {
            cf_op = new CFRoliePolyStrategy("RoliePoly", app_initializer->getComponentDatabase("ComplexFluid"));
        }
        else
        {
            TBOX_ERROR("Unknown type");
        }
        TensorEvolutionType evolve_type =
            IBAMR::string_to_enum<TensorEvolutionType>(input_db->getString("EVOLUTION_TYPE"));
        Pointer<CartGridFunction> exact_fcn = new muParserCartGridFunction(
            "ComplexFluid", app_initializer->getComponentDatabase("ComplexFluid")->getDatabase("FCN"), grid_geometry);

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

        Pointer<CellVariableNd<double> > c_var = new CellVariableNd<double>("C", NDIM * (NDIM + 1) / 2);
        auto var_db = VariableDatabaseNd::getDatabase();
        int c_idx = var_db->registerVariableAndContext(c_var, var_db->getContext("CTX"));
        int r_idx = var_db->registerClonedPatchDataIndex(c_var, c_idx);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(c_idx);
            level->allocatePatchData(r_idx);
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                Pointer<PatchNd> patch = level->getPatch(p());
                Pointer<CellDataNd<double> > c_data = patch->getPatchData(c_idx);
                const BoxNd& box = patch->getBox();

                for (CellIteratorNd ci(box); ci; ci++)
                {
                    const CellIndexNd& idx = ci();
                    switch (evolve_type)
                    {
                    case IBAMR::STANDARD:
                    case SQUARE_ROOT:
#if (NDIM == 2)
                        (*c_data)(idx, 0) = 1.0;
                        (*c_data)(idx, 1) = 1.0;
                        (*c_data)(idx, 2) = 0.0;
#endif
#if (NDIM == 3)
                        (*c_data)(idx, 0) = 1.0;
                        (*c_data)(idx, 1) = 1.0;
                        (*c_data)(idx, 2) = 1.0;
                        (*c_data)(idx, 3) = 0.0;
                        (*c_data)(idx, 4) = 0.0;
                        (*c_data)(idx, 5) = 0.0;
#endif
                        break;
                    case LOGARITHM:
#if (NDIM == 2)
                        (*c_data)(idx, 0) = 0.0;
                        (*c_data)(idx, 1) = 0.0;
                        (*c_data)(idx, 2) = 0.0;
#endif
#if (NDIM == 3)
                        (*c_data)(idx, 0) = 0.0;
                        (*c_data)(idx, 1) = 0.0;
                        (*c_data)(idx, 2) = 0.0;
                        (*c_data)(idx, 3) = 0.0;
                        (*c_data)(idx, 4) = 0.0;
                        (*c_data)(idx, 5) = 0.0;
#endif
                        break;
                    case UNKNOWN_TENSOR_EVOLUTION_TYPE:
                        TBOX_ERROR("Unknown evolution type");
                        break;
                    }
                }
            }
        }

        cf_op->computeRelaxation(r_idx, c_var, c_idx, c_var, evolve_type, patch_hierarchy, 0.0);
        exact_fcn->setDataOnPatchHierarchy(
            c_idx, c_var, patch_hierarchy, 0.0, false, 0, patch_hierarchy->getFinestLevelNumber());

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                Pointer<PatchNd> patch = level->getPatch(p());
                const BoxNd& box = patch->getBox();
                Pointer<CellDataNd<double> > exact_data = patch->getPatchData(c_idx);
                Pointer<CellDataNd<double> > relax_data = patch->getPatchData(r_idx);

                for (CellIteratorNd ci(box); ci; ci++)
                {
                    const CellIndexNd& idx = ci();
                    for (int d = 0; d < (NDIM * (NDIM + 1) / 2); ++d)
                    {
                        double relax_val = (*relax_data)(idx, d);
                        double exact_val = (*exact_data)(idx, d);
                        if (!IBTK::rel_equal_eps(exact_val, relax_val))
                        {
                            pout << "Incorrect value found at patch index " << idx << " at depth " << d << "\n";
                            pout << "Found " << relax_val << ". Expected " << exact_val << "\n";
                        }
                    }
                }
            }
        }

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(c_idx);
            level->deallocatePatchData(r_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown
    return 0;
} // main
