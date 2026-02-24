// ---------------------------------------------------------------------
//
// Copyright (c) 2025 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "ibtk/samrai_compatibility_names.h"
// SAMRAI INCLUDES
#include "ibamr/CFGiesekusStrategy.h"
#include "ibamr/CFOldroydBStrategy.h"
#include "ibamr/CFRoliePolyStrategy.h"
#include "ibamr/ibamr_enums.h"
#include <ibamr/CFUpperConvectiveOperator.h>

#include "ibtk/muParserCartGridFunction.h"
#include <ibtk/AppInitializer.h>
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include "SAMRAIArray.h"
#include "SAMRAIBergerRigoutsos.h"
#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAICellData.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIFaceVariable.h"
#include "SAMRAIGriddingAlgorithm.h"
#include "SAMRAIHierarchyCellDataOpsReal.h"
#include "SAMRAILoadBalancer.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIRobinBcCoefStrategy.h"
#include "SAMRAISAMRAIVectorReal.h"
#include "SAMRAIStandardTagAndInitialize.h"
#include "SAMRAIVariableDatabase.h"

#include <petscsys.h>

#include <SAMRAI_config.h>

#include "ibamr/namespaces.h"

// Dummy CF Strategy class that does nothing
class DummyCFStrategy : public CFStrategy
{
public:
    DummyCFStrategy(std::string object_name) : CFStrategy(std::move(object_name))
    {
    }

    void computeRelaxation(int R_idx,
                           SAMRAIPointer<SAMRAICellVariable<double>> /* R_var */,
                           int /* C_idx */,
                           SAMRAIPointer<SAMRAICellVariable<double>> /* C_var */,
                           TensorEvolutionType /* evolve_type */,
                           SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy,
                           double /* data_time */) override
    {
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = hierarchy->getPatchLevel(ln);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                SAMRAIPointer<SAMRAICellData<double>> sig_data = patch->getPatchData(R_idx);
                sig_data->fillAll(0.0);
            }
        }
    }

    void computeStress(int /* sig_idx */,
                       SAMRAIPointer<SAMRAICellVariable<double>> /* sig_var */,
                       SAMRAIPointer<SAMRAIPatchHierarchy> /* hierarchy */,
                       double /* data_time */) override
    {
        // Do nothing
    }
};

// Dummy Convective Operator
class DummyConvectiveOperator : public ConvectiveOperator
{
public:
    DummyConvectiveOperator(std::string object_name) : ConvectiveOperator(std::move(object_name), ADVECTIVE)
    {
    }

    void applyConvectiveOperator(int /* Q_idx */, int N_idx) override
    {
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = d_hierarchy->getPatchLevel(ln);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                SAMRAIPointer<SAMRAICellData<double>> N_data = patch->getPatchData(N_idx);
                N_data->fillAll(0.0);
            }
        }
    }

    void initializeOperatorState(const SAMRAISAMRAIVectorReal<double>& in,
                                 const SAMRAISAMRAIVectorReal<double>& out) override
    {
        ConvectiveOperator::initializeOperatorState(in, out);
        d_hierarchy = in.getPatchHierarchy();
    }

private:
    SAMRAIPointer<SAMRAIPatchHierarchy> d_hierarchy;
};

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
        SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geometry = new SAMRAICartesianGridGeometry(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy = new SAMRAIPatchHierarchy("PatchHierarchy", grid_geometry);
        SAMRAIPointer<SAMRAIStandardTagAndInitialize> error_detector = new SAMRAIStandardTagAndInitialize(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        SAMRAIPointer<SAMRAIBergerRigoutsos> box_generator = new SAMRAIBergerRigoutsos();
        SAMRAIPointer<SAMRAILoadBalancer> load_balancer =
            new SAMRAILoadBalancer("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        SAMRAIPointer<SAMRAIGriddingAlgorithm> gridding_algorithm =
            new SAMRAIGriddingAlgorithm("GriddingAlgorithm",
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

        // Only setup for periodic boundary conditions
        std::vector<SAMRAIRobinBcCoefStrategy*> c_bc_coefs(NDIM * (NDIM + 1) / 2, nullptr);
        std::vector<SAMRAIRobinBcCoefStrategy*> u_bc_coefs(NDIM, nullptr);

        SAMRAIPointer<SAMRAICellVariable<double>> c_var = new SAMRAICellVariable<double>("C", NDIM * (NDIM + 1) / 2);
        SAMRAIPointer<SAMRAIFaceVariable<double>> u_var = new SAMRAIFaceVariable<double>("U_var");
        auto var_db = SAMRAIVariableDatabase::getDatabase();
        int c_idx = var_db->registerVariableAndContext(c_var, var_db->getContext("CTX"));
        int r_idx = var_db->registerClonedPatchDataIndex(c_var, c_idx);
        int e_idx = var_db->registerClonedPatchDataIndex(c_var, c_idx);
        int u_idx = var_db->registerVariableAndContext(u_var, var_db->getContext("CTX"), 1);

        SAMRAIPointer<CFStrategy> cf_strat = new DummyCFStrategy("CF_STRATEGY");
        SAMRAIPointer<ConvectiveOperator> cf_co = new DummyConvectiveOperator("CF_CO");
        CFUpperConvectiveOperator cf_ucd(
            "CF_UCD", c_var, input_db->getDatabase("CF_UCD"), cf_co, ADVECTIVE, c_bc_coefs, u_bc_coefs);
        cf_ucd.registerCFStrategy(cf_strat);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(c_idx, 0.0);
            level->allocatePatchData(r_idx, 0.0);
            level->allocatePatchData(e_idx, 0.0);
            level->allocatePatchData(u_idx, 0.0);
        }

        SAMRAIPointer<CartGridFunction> exact_fcn =
            new muParserCartGridFunction("ExactFcn", input_db->getDatabase("ExactFcn"), grid_geometry);
        SAMRAIPointer<CartGridFunction> u_fcn =
            new muParserCartGridFunction("UFcn", input_db->getDatabase("UFcn"), grid_geometry);
        SAMRAIPointer<CartGridFunction> stress_fcn =
            new muParserCartGridFunction("StressFcn", input_db->getDatabase("StressFcn"), grid_geometry);

        stress_fcn->setDataOnPatchHierarchy(
            c_idx, c_var, patch_hierarchy, 0.0, false, 0, patch_hierarchy->getFinestLevelNumber());
        exact_fcn->setDataOnPatchHierarchy(
            e_idx, c_var, patch_hierarchy, 0.0, false, 0, patch_hierarchy->getFinestLevelNumber());
        u_fcn->setDataOnPatchHierarchy(
            u_idx, u_var, patch_hierarchy, 0.0, false, 0, patch_hierarchy->getFinestLevelNumber());

        SAMRAISAMRAIVectorReal<double> c_vec("C", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber()),
            r_vec("R", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        c_vec.addComponent(c_var, c_idx);
        r_vec.addComponent(c_var, r_idx);

        cf_ucd.setAdvectionVelocity(u_idx);
        cf_ucd.initializeOperatorState(c_vec, r_vec);
        cf_ucd.applyConvectiveOperator(c_idx, r_idx);

        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAIHierarchyCellDataOpsReal<double> hier_cc_data_ops(patch_hierarchy);
        hier_cc_data_ops.subtract(e_idx, e_idx, r_idx);
        pout << "Errors in UCD:\n";
        pout << "L1 norm:  " << hier_cc_data_ops.L1Norm(e_idx, wgt_cc_idx) << "\n";
        pout << "L2 norm:  " << hier_cc_data_ops.L2Norm(e_idx, wgt_cc_idx) << "\n";
        pout << "max norm: " << hier_cc_data_ops.maxNorm(e_idx, wgt_cc_idx) << "\n";

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(c_idx);
            level->deallocatePatchData(r_idx);
            level->deallocatePatchData(e_idx);
            level->deallocatePatchData(u_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown
    return 0;
} // main
