// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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
#include "ibamr/CFGiesekusRelaxation.h"
#include "ibamr/CFOldroydBRelaxation.h"
#include "ibamr/CFRoliePolyRelaxation.h"
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
        Pointer<CFRelaxationOperator> cf_op;
        std::string relax_op = input_db->getString("RELAX_OP");
        if (relax_op == "OLDROYDB")
        {
            cf_op = new CFOldroydBRelaxation("OldroydB", app_initializer->getComponentDatabase("ComplexFluid"));
        }
        else if (relax_op == "GIESEKUS")
        {
            cf_op = new CFGiesekusRelaxation("Giesekus", app_initializer->getComponentDatabase("ComplexFluid"));
        }
        else if (relax_op == "ROLIEPOLY")
        {
            cf_op = new CFRoliePolyRelaxation("RoliePoly", app_initializer->getComponentDatabase("ComplexFluid"));
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

        Pointer<CellVariable<NDIM, double> > c_var = new CellVariable<NDIM, double>("C", NDIM * (NDIM + 1) / 2);
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        int c_idx = var_db->registerVariableAndContext(c_var, var_db->getContext("CTX"));
        int r_idx = var_db->registerClonedPatchDataIndex(c_var, c_idx);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(c_idx);
            level->allocatePatchData(r_idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM, double> > c_data = patch->getPatchData(c_idx);
                const Box<NDIM>& box = patch->getBox();

                for (CellIterator<NDIM> ci(box); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
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

        cf_op->setPatchDataIndex(c_idx);
        cf_op->setDataOnPatchHierarchy(
            r_idx, c_var, patch_hierarchy, 0.0, false, 0, patch_hierarchy->getFinestLevelNumber());
        exact_fcn->setDataOnPatchHierarchy(
            c_idx, c_var, patch_hierarchy, 0.0, false, 0, patch_hierarchy->getFinestLevelNumber());

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& box = patch->getBox();
                Pointer<CellData<NDIM, double> > exact_data = patch->getPatchData(c_idx);
                Pointer<CellData<NDIM, double> > relax_data = patch->getPatchData(r_idx);

                for (CellIterator<NDIM> ci(box); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    for (int d = 0; d < (NDIM * (NDIM + 1) / 2); ++d)
                    {
                        double relax_val = (*relax_data)(idx, d);
                        double exact_val = (*exact_data)(idx, d);
                        if (!MathUtilities<double>::equalEps(exact_val, relax_val))
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
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(c_idx);
            level->deallocatePatchData(r_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown
    return 0;
} // main
