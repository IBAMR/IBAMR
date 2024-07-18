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

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/SAMRAIDataCache.h>

#include <CellVariable.h>
#include <EdgeVariable.h>
#include <FaceVariable.h>
#include <NodeVariable.h>
#include <OuteredgeVariable.h>
#include <OuterfaceVariable.h>
#include <OuternodeVariable.h>
#include <OutersideVariable.h>
#include <SideVariable.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

void
print_patch_descriptor_data(SAMRAIPointer<PatchDescriptorNd> descriptor)
{
    int max_number_registered_components = descriptor->getMaxNumberRegisteredComponents();
    pout << "patch descriptor configuration:\n";
    pout << "max_number_registered_components = " << max_number_registered_components << "\n";
    pout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    for (int i = 0; i < max_number_registered_components; i++)
    {
        auto factory = descriptor->getPatchDataFactory(i);
        pout << "Patch Data Index = " << i << "\n";
        if (factory)
        {
            pout << "Patch Data Factory Name = " << descriptor->mapIndexToName(i) << "\n";
        }
        else
        {
            pout << "Patch Data Factory is null\n";
        }
    }
    pout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

// Test SAMRAIDataCache for various data types.
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        SAMRAIPointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_poisson.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        SAMRAIPointer<CartesianGridGeometryNd> grid_geometry = new CartesianGridGeometryNd(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        SAMRAIPointer<PatchHierarchyNd> patch_hierarchy = new PatchHierarchyNd("PatchHierarchy", grid_geometry);
        SAMRAIPointer<StandardTagAndInitializeNd> error_detector = new StandardTagAndInitializeNd(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        SAMRAIPointer<BergerRigoutsosNd> box_generator = new BergerRigoutsosNd();
        SAMRAIPointer<LoadBalancerNd> load_balancer =
            new LoadBalancerNd("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        SAMRAIPointer<GriddingAlgorithmNd> gridding_algorithm =
            new GriddingAlgorithmNd("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

        // Create variables and register them with the variable database.
        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
        SAMRAIPointer<VariableContext> ctx = var_db->getContext("context");

        pout << "adding double variables\n";
        SAMRAIPointer<CellVariableNd<double> > cc_double_var = new CellVariableNd<double>("cc_double");
        SAMRAIPointer<EdgeVariableNd<double> > ec_double_var = new EdgeVariableNd<double>("ec_double");
        SAMRAIPointer<FaceVariableNd<double> > fc_double_var = new FaceVariableNd<double>("fc_double");
        SAMRAIPointer<NodeVariableNd<double> > nc_double_var = new NodeVariableNd<double>("nc_double");
        SAMRAIPointer<OuteredgeVariableNd<double> > oec_double_var = new OuteredgeVariableNd<double>("oec_double");
        SAMRAIPointer<OuterfaceVariableNd<double> > ofc_double_var = new OuterfaceVariableNd<double>("ofc_double");
        SAMRAIPointer<OuternodeVariableNd<double> > onc_double_var = new OuternodeVariableNd<double>("onc_double");
        SAMRAIPointer<OutersideVariableNd<double> > osc_double_var = new OutersideVariableNd<double>("osc_double");
        SAMRAIPointer<SideVariableNd<double> > sc_double_var = new SideVariableNd<double>("sc_double");
        std::vector<SAMRAIPointer<VariableNd> > double_vars{ cc_double_var,  ec_double_var,  fc_double_var,
                                                             nc_double_var,  oec_double_var, ofc_double_var,
                                                             onc_double_var, osc_double_var, sc_double_var };
        std::vector<int> double_idxs;
        for (const auto& var : double_vars)
        {
            double_idxs.push_back(var_db->registerVariableAndContext(var, ctx, IntVectorNd(1)));
        }
        print_patch_descriptor_data(patch_hierarchy->getPatchDescriptor());

        pout << "adding float variables\n";
        SAMRAIPointer<CellVariableNd<float> > cc_float_var = new CellVariableNd<float>("cc_float");
        SAMRAIPointer<EdgeVariableNd<float> > ec_float_var = new EdgeVariableNd<float>("ec_float");
        SAMRAIPointer<FaceVariableNd<float> > fc_float_var = new FaceVariableNd<float>("fc_float");
        SAMRAIPointer<NodeVariableNd<float> > nc_float_var = new NodeVariableNd<float>("nc_float");
        SAMRAIPointer<OuteredgeVariableNd<float> > oec_float_var = new OuteredgeVariableNd<float>("oec_float");
        SAMRAIPointer<OuterfaceVariableNd<float> > ofc_float_var = new OuterfaceVariableNd<float>("ofc_float");
        SAMRAIPointer<OuternodeVariableNd<float> > onc_float_var = new OuternodeVariableNd<float>("onc_float");
        SAMRAIPointer<OutersideVariableNd<float> > osc_float_var = new OutersideVariableNd<float>("osc_float");
        SAMRAIPointer<SideVariableNd<float> > sc_float_var = new SideVariableNd<float>("sc_float");

        std::vector<SAMRAIPointer<VariableNd> > float_vars{ cc_float_var,  ec_float_var,  fc_float_var,
                                                            nc_float_var,  oec_float_var, ofc_float_var,
                                                            onc_float_var, osc_float_var, sc_float_var };
        std::vector<int> float_idxs;
        for (const auto& var : float_vars)
        {
            float_idxs.push_back(var_db->registerVariableAndContext(var, ctx, IntVectorNd(1)));
        }
        print_patch_descriptor_data(patch_hierarchy->getPatchDescriptor());

        pout << "adding int variables\n";
        SAMRAIPointer<CellVariableNd<int> > cc_int_var = new CellVariableNd<int>("cc_int");
        SAMRAIPointer<EdgeVariableNd<int> > ec_int_var = new EdgeVariableNd<int>("ec_int");
        SAMRAIPointer<FaceVariableNd<int> > fc_int_var = new FaceVariableNd<int>("fc_int");
        SAMRAIPointer<NodeVariableNd<int> > nc_int_var = new NodeVariableNd<int>("nc_int");
        SAMRAIPointer<OuteredgeVariableNd<int> > oec_int_var = new OuteredgeVariableNd<int>("oec_int");
        SAMRAIPointer<OuterfaceVariableNd<int> > ofc_int_var = new OuterfaceVariableNd<int>("ofc_int");
        SAMRAIPointer<OuternodeVariableNd<int> > onc_int_var = new OuternodeVariableNd<int>("onc_int");
        SAMRAIPointer<OutersideVariableNd<int> > osc_int_var = new OutersideVariableNd<int>("osc_int");
        SAMRAIPointer<SideVariableNd<int> > sc_int_var = new SideVariableNd<int>("sc_int");

        std::vector<SAMRAIPointer<VariableNd> > int_vars{ cc_int_var,  ec_int_var,  fc_int_var,
                                                          nc_int_var,  oec_int_var, ofc_int_var,
                                                          onc_int_var, osc_int_var, sc_int_var };
        std::vector<int> int_idxs;
        for (const auto& var : int_vars)
        {
            int_idxs.push_back(var_db->registerVariableAndContext(var, ctx, IntVectorNd(1)));
        }
        print_patch_descriptor_data(patch_hierarchy->getPatchDescriptor());

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

        {
            SAMRAIDataCache cached_data;
            cached_data.setPatchHierarchy(patch_hierarchy);
            cached_data.resetLevels(0, patch_hierarchy->getFinestLevelNumber());

            std::vector<SAMRAIDataCache::CachedPatchDataIndex> cached_double_idxs, cached_float_idxs, cached_int_idxs;
            for (unsigned int k = 0; k < 2; ++k)
            {
                pout << "cloning pass " << k << "\n";

                for (auto idx : double_idxs) cached_double_idxs.push_back(cached_data.getCachedPatchDataIndex(idx));

                for (auto idx : float_idxs) cached_float_idxs.push_back(cached_data.getCachedPatchDataIndex(idx));

                for (auto idx : int_idxs) cached_int_idxs.push_back(cached_data.getCachedPatchDataIndex(idx));

                print_patch_descriptor_data(patch_hierarchy->getPatchDescriptor());
            }
            pout << "allowing the cache to go out of scope\n";
        }
        print_patch_descriptor_data(patch_hierarchy->getPatchDescriptor());

    } // cleanup dynamically allocated objects prior to shutdown
} // main
