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

#include <tbox/Array.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include "ibamr/namespaces.h"

// Dummy CF Strategy class that does nothing
class DummyCFStrategy : public CFStrategy
{
public:
    DummyCFStrategy(std::string object_name) : CFStrategy(std::move(object_name))
    {
    }

    void computeRelaxation(int R_idx,
                           Pointer<CellVariable<NDIM, double> > /* R_var */,
                           int /* C_idx */,
                           Pointer<CellVariable<NDIM, double> > /* C_var */,
                           TensorEvolutionType /* evolve_type */,
                           Pointer<PatchHierarchy<NDIM> > hierarchy,
                           double /* data_time */) override
    {
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM, double> > sig_data = patch->getPatchData(R_idx);
                sig_data->fillAll(0.0);
            }
        }
    }

    void computeStress(int /* sig_idx */,
                       Pointer<CellVariable<NDIM, double> > /* sig_var */,
                       Pointer<PatchHierarchy<NDIM> > /* hierarchy */,
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
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM, double> > N_data = patch->getPatchData(N_idx);
                N_data->fillAll(0.0);
            }
        }
    }

    void initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAIVectorReal<NDIM, double>& out) override
    {
        ConvectiveOperator::initializeOperatorState(in, out);
        d_hierarchy = in.getPatchHierarchy();
    }

private:
    Pointer<PatchHierarchy<NDIM> > d_hierarchy;
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
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_poisson.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
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

        // Only setup for periodic boundary conditions
        std::vector<RobinBcCoefStrategy<NDIM>*> c_bc_coefs(NDIM * (NDIM + 1) / 2, nullptr);
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);

        Pointer<CellVariable<NDIM, double> > c_var = new CellVariable<NDIM, double>("C", NDIM * (NDIM + 1) / 2);
        Pointer<FaceVariable<NDIM, double> > u_var = new FaceVariable<NDIM, double>("U_var");
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        int c_idx = var_db->registerVariableAndContext(c_var, var_db->getContext("CTX"));
        int r_idx = var_db->registerClonedPatchDataIndex(c_var, c_idx);
        int e_idx = var_db->registerClonedPatchDataIndex(c_var, c_idx);
        int u_idx = var_db->registerVariableAndContext(u_var, var_db->getContext("CTX"), 1);

        Pointer<CFStrategy> cf_strat = new DummyCFStrategy("CF_STRATEGY");
        Pointer<ConvectiveOperator> cf_co = new DummyConvectiveOperator("CF_CO");
        CFUpperConvectiveOperator cf_ucd(
            "CF_UCD", c_var, input_db->getDatabase("CF_UCD"), cf_co, ADVECTIVE, c_bc_coefs, u_bc_coefs);
        cf_ucd.registerCFStrategy(cf_strat);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(c_idx, 0.0);
            level->allocatePatchData(r_idx, 0.0);
            level->allocatePatchData(e_idx, 0.0);
            level->allocatePatchData(u_idx, 0.0);
        }

        Pointer<CartGridFunction> exact_fcn =
            new muParserCartGridFunction("ExactFcn", input_db->getDatabase("ExactFcn"), grid_geometry);
        Pointer<CartGridFunction> u_fcn =
            new muParserCartGridFunction("UFcn", input_db->getDatabase("UFcn"), grid_geometry);
        Pointer<CartGridFunction> stress_fcn =
            new muParserCartGridFunction("StressFcn", input_db->getDatabase("StressFcn"), grid_geometry);

        stress_fcn->setDataOnPatchHierarchy(
            c_idx, c_var, patch_hierarchy, 0.0, false, 0, patch_hierarchy->getFinestLevelNumber());
        exact_fcn->setDataOnPatchHierarchy(
            e_idx, c_var, patch_hierarchy, 0.0, false, 0, patch_hierarchy->getFinestLevelNumber());
        u_fcn->setDataOnPatchHierarchy(
            u_idx, u_var, patch_hierarchy, 0.0, false, 0, patch_hierarchy->getFinestLevelNumber());

        SAMRAIVectorReal<NDIM, double> c_vec("C", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber()),
            r_vec("R", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        c_vec.addComponent(c_var, c_idx);
        r_vec.addComponent(c_var, r_idx);

        cf_ucd.setAdvectionVelocity(u_idx);
        cf_ucd.initializeOperatorState(c_vec, r_vec);
        cf_ucd.applyConvectiveOperator(c_idx, r_idx);

        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy);
        hier_cc_data_ops.subtract(e_idx, e_idx, r_idx);
        pout << "Errors in UCD:\n";
        pout << "L1 norm:  " << hier_cc_data_ops.L1Norm(e_idx, wgt_cc_idx) << "\n";
        pout << "L2 norm:  " << hier_cc_data_ops.L2Norm(e_idx, wgt_cc_idx) << "\n";
        pout << "max norm: " << hier_cc_data_ops.maxNorm(e_idx, wgt_cc_idx) << "\n";

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(c_idx);
            level->deallocatePatchData(r_idx);
            level->deallocatePatchData(e_idx);
            level->deallocatePatchData(u_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown
    return 0;
} // main
