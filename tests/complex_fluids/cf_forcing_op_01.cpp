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

#include "ibtk/samrai_compatibility_names.h"

#include <SAMRAI_config.h>

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for main SAMRAI objects
#include "ibtk/muParserCartGridFunction.h"

#include "SAMRAIArray.h"
#include "SAMRAIBergerRigoutsos.h"
#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAICellData.h"
#include "SAMRAICellIndex.h"
#include "SAMRAICellIterator.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIGriddingAlgorithm.h"
#include "SAMRAIHierarchyCellDataOpsReal.h"
#include "SAMRAIHierarchySideDataOpsReal.h"
#include "SAMRAILoadBalancer.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAISideVariable.h"
#include "SAMRAIStandardTagAndInitialize.h"
#include "SAMRAIVariable.h"
#include "SAMRAIVariableDatabase.h"

// Headers for application-specific algorithm/data structure objects
#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/CFGiesekusStrategy.h"
#include "ibamr/CFINSForcing.h"
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
        SAMRAIPointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_poisson.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        SAMRAIPointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator =
            new AdvDiffSemiImplicitHierarchyIntegrator(
                "AdvDiffHierarchyIntegrator", app_initializer->getComponentDatabase("AdvDiffHierarchyIntegrator"));
        SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geometry = new SAMRAICartesianGridGeometry(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy = new SAMRAIPatchHierarchy("PatchHierarchy", grid_geometry);
        SAMRAIPointer<SAMRAIStandardTagAndInitialize> error_detector =
            new SAMRAIStandardTagAndInitialize("StandardTagAndInitialize",
                                               adv_diff_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        SAMRAIPointer<SAMRAIBergerRigoutsos> box_generator = new SAMRAIBergerRigoutsos();
        SAMRAIPointer<SAMRAILoadBalancer> load_balancer =
            new SAMRAILoadBalancer("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        SAMRAIPointer<SAMRAIGriddingAlgorithm> gridding_algorithm =
            new SAMRAIGriddingAlgorithm("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        SAMRAIPointer<CartGridFunction> u_fcn = nullptr;
        SAMRAIPointer<CFINSForcing> cf_forcing = new CFINSForcing("ComplexFluid",
                                                                  app_initializer->getComponentDatabase("ComplexFluid"),
                                                                  u_fcn,
                                                                  grid_geometry,
                                                                  adv_diff_integrator,
                                                                  app_initializer->getVisItDataWriter());
        SAMRAIPointer<CartGridFunction> exact_fcn = new muParserCartGridFunction(
            "ComplexFluid", app_initializer->getComponentDatabase("ComplexFluid")->getDatabase("FCN"), grid_geometry);

        // Initialize the AMR patch hierarchy.
        adv_diff_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        SAMRAIPointer<SAMRAIVariable> c_var;
        std::string var_centering = input_db->getString("VAR_CENTERING");
        if (var_centering == "SIDE")
        {
            c_var = new SAMRAISideVariable<double>("DIV");
        }
        else if (var_centering == "CELL")
        {
            c_var = new SAMRAICellVariable<double>("DIV", NDIM);
        }
        else
        {
            TBOX_ERROR("Incorrect centering.");
        }
        auto var_db = SAMRAIVariableDatabase::getDatabase();
        const int c_idx = var_db->registerVariableAndContext(c_var, var_db->getContext("CTX"));
        const int c_cloned_idx = var_db->registerClonedPatchDataIndex(c_var, c_idx);
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(c_idx);
            level->allocatePatchData(c_cloned_idx);
        }

        bool test_draw = input_db->getBoolWithDefault("TEST_DRAW", false);

        cf_forcing->setDataOnPatchHierarchy(
            c_idx, c_var, patch_hierarchy, 0.0, false, 0, patch_hierarchy->getFinestLevelNumber());
        exact_fcn->setDataOnPatchHierarchy(
            c_cloned_idx, c_var, patch_hierarchy, 0.0, false, 0, patch_hierarchy->getFinestLevelNumber());

        HierarchyMathOps hier_math_ops("HierMathOps", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(0, patch_hierarchy->getFinestLevelNumber());
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        SAMRAIHierarchyCellDataOpsReal<double> hier_cc_data_ops(
            patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIHierarchySideDataOpsReal<double> hier_sc_data_ops(
            patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        if (test_draw)
        {
            bool error = false;
            for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
                for (SAMRAIPatchLevel::Iterator p(level); p; p++)
                {
                    SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                    SAMRAIPointer<SAMRAICellData<double>> draw_data = patch->getPatchData(
                        var_db->getVariable("ComplexFluid::conform_draw"), var_db->getContext("ComplexFluid::CONTEXT"));
                    SAMRAIPointer<SAMRAICellData<double>> C_data = patch->getPatchData(cf_forcing->getVariableIdx());
                    for (SAMRAICellIterator ci(patch->getBox()); ci; ci++)
                    {
                        const SAMRAICellIndex& idx = ci();
#if (NDIM == 2)
                        if (!(abs_equal_eps((*C_data)(idx, 0), (*draw_data)(idx, 0)) &&
                              abs_equal_eps((*C_data)(idx, 2), (*draw_data)(idx, 1)) &&
                              abs_equal_eps((*C_data)(idx, 2), (*draw_data)(idx, 2)) &&
                              abs_equal_eps((*C_data)(idx, 1), (*draw_data)(idx, 3))))
                            error = true;
#endif
#if (NDIM == 3)
                        if (!(abs_equal_eps((*C_data)(idx, 0), (*draw_data)(idx, 0)) &&
                              abs_equal_eps((*C_data)(idx, 5), (*draw_data)(idx, 1)) &&
                              abs_equal_eps((*C_data)(idx, 4), (*draw_data)(idx, 2)) &&
                              abs_equal_eps((*C_data)(idx, 5), (*draw_data)(idx, 3)) &&
                              abs_equal_eps((*C_data)(idx, 1), (*draw_data)(idx, 4)) &&
                              abs_equal_eps((*C_data)(idx, 3), (*draw_data)(idx, 5)) &&
                              abs_equal_eps((*C_data)(idx, 4), (*draw_data)(idx, 6)) &&
                              abs_equal_eps((*C_data)(idx, 3), (*draw_data)(idx, 7)) &&
                              abs_equal_eps((*C_data)(idx, 2), (*draw_data)(idx, 8))))
                            error = true;
#endif
                    }
                }
            }

            if (error) pout << "Drawing data is incorrect\n";
        }
        else if (var_centering == "CELL")
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
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(c_idx);
            level->deallocatePatchData(c_cloned_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown
    return 0;
} // main
