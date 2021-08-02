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
#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/CFGiesekusRelaxation.h"
#include "ibamr/CFINSForcing.h"
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
        Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator =
            new AdvDiffSemiImplicitHierarchyIntegrator(
                "AdvDiffHierarchyIntegrator", app_initializer->getComponentDatabase("AdvDiffHierarchyIntegrator"));
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               adv_diff_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        Pointer<CartGridFunction> u_fcn = nullptr;
        Pointer<CFINSForcing> cf_forcing = new CFINSForcing("ComplexFluid",
                                                            app_initializer->getComponentDatabase("ComplexFluid"),
                                                            u_fcn,
                                                            grid_geometry,
                                                            adv_diff_integrator,
                                                            adv_diff_integrator->getVisItDataWriter());
        Pointer<CartGridFunction> exact_fcn = new muParserCartGridFunction(
            "ComplexFluid", app_initializer->getComponentDatabase("ComplexFluid")->getDatabase("FCN"), grid_geometry);

        // Initialize the AMR patch hierarchy.
        adv_diff_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        Pointer<Variable<NDIM> > c_var;
        std::string var_centering = input_db->getString("VAR_CENTERING");
        if (var_centering == "SIDE")
        {
            c_var = new SideVariable<NDIM, double>("DIV");
        }
        else if (var_centering == "CELL")
        {
            c_var = new CellVariable<NDIM, double>("DIV", NDIM);
        }
        else
        {
            TBOX_ERROR("Incorrect centering.");
        }
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        const int c_idx = var_db->registerVariableAndContext(c_var, var_db->getContext("CTX"));
        const int c_cloned_idx = var_db->registerClonedPatchDataIndex(c_var, c_idx);
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(c_idx);
            level->allocatePatchData(c_cloned_idx);
        }

        cf_forcing->setDataOnPatchHierarchy(
            c_idx, c_var, patch_hierarchy, 0.0, false, 0, patch_hierarchy->getFinestLevelNumber());
        exact_fcn->setDataOnPatchHierarchy(
            c_cloned_idx, c_var, patch_hierarchy, 0.0, false, 0, patch_hierarchy->getFinestLevelNumber());

        HierarchyMathOps hier_math_ops("HierMathOps", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(0, patch_hierarchy->getFinestLevelNumber());
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(
            patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(
            patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        if (var_centering == "CELL")
        {
            hier_cc_data_ops.subtract(c_idx, c_idx, c_cloned_idx);
            plog << "Error norms:\n"
                 << "  L1-norm:  " << std::setprecision(8) << hier_cc_data_ops.L1Norm(c_idx, wgt_cc_idx) << "\n"
                 << "  L2-norm:  " << std::setprecision(8) << hier_cc_data_ops.L2Norm(c_idx, wgt_cc_idx) << "\n"
                 << "  max-norm: " << std::setprecision(8) << hier_cc_data_ops.maxNorm(c_idx, wgt_cc_idx) << "\n";
        }
        else
        {
            hier_sc_data_ops.subtract(c_idx, c_idx, c_cloned_idx);
            plog << "Error norms:\n"
                 << "  L1-norm:  " << std::setprecision(8) << hier_sc_data_ops.L1Norm(c_idx, wgt_sc_idx) << "\n"
                 << "  L2-norm:  " << std::setprecision(8) << hier_sc_data_ops.L2Norm(c_idx, wgt_sc_idx) << "\n"
                 << "  max-norm: " << std::setprecision(8) << hier_sc_data_ops.maxNorm(c_idx, wgt_sc_idx) << "\n";
        }

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(c_idx);
            level->deallocatePatchData(c_cloned_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown
    return 0;
} // main
