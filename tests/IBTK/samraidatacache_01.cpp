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

// Headers for main SAMRAI objects
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIBergerRigoutsos.h"
#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIEdgeVariable.h"
#include "SAMRAIFaceVariable.h"
#include "SAMRAIGriddingAlgorithm.h"
#include "SAMRAIIntVector.h"
#include "SAMRAILoadBalancer.h"
#include "SAMRAINodeVariable.h"
#include "SAMRAIOuteredgeVariable.h"
#include "SAMRAIOuterfaceVariable.h"
#include "SAMRAIOuternodeVariable.h"
#include "SAMRAIOutersideVariable.h"
#include "SAMRAIPatchDescriptor.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPointer.h"
#include "SAMRAISideVariable.h"
#include "SAMRAIStandardTagAndInitialize.h"
#include "SAMRAIVariable.h"
#include "SAMRAIVariableDatabase.h"

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/SAMRAIDataCache.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

void
print_patch_descriptor_data(SAMRAIPointer<SAMRAIPatchDescriptor> descriptor)
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

        // Create variables and register them with the variable database.
        SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
        SAMRAIPointer<VariableContext> ctx = var_db->getContext("context");

        pout << "adding double variables\n";
        SAMRAIPointer<SAMRAICellVariable<double>> cc_double_var = new SAMRAICellVariable<double>("cc_double");
        SAMRAIPointer<SAMRAIEdgeVariable<double>> ec_double_var = new SAMRAIEdgeVariable<double>("ec_double");
        SAMRAIPointer<SAMRAIFaceVariable<double>> fc_double_var = new SAMRAIFaceVariable<double>("fc_double");
        SAMRAIPointer<SAMRAINodeVariable<double>> nc_double_var = new SAMRAINodeVariable<double>("nc_double");
        SAMRAIPointer<SAMRAIOuteredgeVariable<double>> oec_double_var =
            new SAMRAIOuteredgeVariable<double>("oec_double");
        SAMRAIPointer<SAMRAIOuterfaceVariable<double>> ofc_double_var =
            new SAMRAIOuterfaceVariable<double>("ofc_double");
        SAMRAIPointer<SAMRAIOuternodeVariable<double>> onc_double_var =
            new SAMRAIOuternodeVariable<double>("onc_double");
        SAMRAIPointer<SAMRAIOutersideVariable<double>> osc_double_var =
            new SAMRAIOutersideVariable<double>("osc_double");
        SAMRAIPointer<SAMRAISideVariable<double>> sc_double_var = new SAMRAISideVariable<double>("sc_double");
        std::vector<SAMRAIPointer<SAMRAIVariable>> double_vars{ cc_double_var,  ec_double_var,  fc_double_var,
                                                                nc_double_var,  oec_double_var, ofc_double_var,
                                                                onc_double_var, osc_double_var, sc_double_var };
        std::vector<int> double_idxs;
        for (const auto& var : double_vars)
        {
            double_idxs.push_back(var_db->registerVariableAndContext(var, ctx, SAMRAIIntVector(1)));
        }
        print_patch_descriptor_data(patch_hierarchy->getPatchDescriptor());

        pout << "adding float variables\n";
        SAMRAIPointer<SAMRAICellVariable<float>> cc_float_var = new SAMRAICellVariable<float>("cc_float");
        SAMRAIPointer<SAMRAIEdgeVariable<float>> ec_float_var = new SAMRAIEdgeVariable<float>("ec_float");
        SAMRAIPointer<SAMRAIFaceVariable<float>> fc_float_var = new SAMRAIFaceVariable<float>("fc_float");
        SAMRAIPointer<SAMRAINodeVariable<float>> nc_float_var = new SAMRAINodeVariable<float>("nc_float");
        SAMRAIPointer<SAMRAIOuteredgeVariable<float>> oec_float_var = new SAMRAIOuteredgeVariable<float>("oec_float");
        SAMRAIPointer<SAMRAIOuterfaceVariable<float>> ofc_float_var = new SAMRAIOuterfaceVariable<float>("ofc_float");
        SAMRAIPointer<SAMRAIOuternodeVariable<float>> onc_float_var = new SAMRAIOuternodeVariable<float>("onc_float");
        SAMRAIPointer<SAMRAIOutersideVariable<float>> osc_float_var = new SAMRAIOutersideVariable<float>("osc_float");
        SAMRAIPointer<SAMRAISideVariable<float>> sc_float_var = new SAMRAISideVariable<float>("sc_float");

        std::vector<SAMRAIPointer<SAMRAIVariable>> float_vars{ cc_float_var,  ec_float_var,  fc_float_var,
                                                               nc_float_var,  oec_float_var, ofc_float_var,
                                                               onc_float_var, osc_float_var, sc_float_var };
        std::vector<int> float_idxs;
        for (const auto& var : float_vars)
        {
            float_idxs.push_back(var_db->registerVariableAndContext(var, ctx, SAMRAIIntVector(1)));
        }
        print_patch_descriptor_data(patch_hierarchy->getPatchDescriptor());

        pout << "adding int variables\n";
        SAMRAIPointer<SAMRAICellVariable<int>> cc_int_var = new SAMRAICellVariable<int>("cc_int");
        SAMRAIPointer<SAMRAIEdgeVariable<int>> ec_int_var = new SAMRAIEdgeVariable<int>("ec_int");
        SAMRAIPointer<SAMRAIFaceVariable<int>> fc_int_var = new SAMRAIFaceVariable<int>("fc_int");
        SAMRAIPointer<SAMRAINodeVariable<int>> nc_int_var = new SAMRAINodeVariable<int>("nc_int");
        SAMRAIPointer<SAMRAIOuteredgeVariable<int>> oec_int_var = new SAMRAIOuteredgeVariable<int>("oec_int");
        SAMRAIPointer<SAMRAIOuterfaceVariable<int>> ofc_int_var = new SAMRAIOuterfaceVariable<int>("ofc_int");
        SAMRAIPointer<SAMRAIOuternodeVariable<int>> onc_int_var = new SAMRAIOuternodeVariable<int>("onc_int");
        SAMRAIPointer<SAMRAIOutersideVariable<int>> osc_int_var = new SAMRAIOutersideVariable<int>("osc_int");
        SAMRAIPointer<SAMRAISideVariable<int>> sc_int_var = new SAMRAISideVariable<int>("sc_int");

        std::vector<SAMRAIPointer<SAMRAIVariable>> int_vars{ cc_int_var,  ec_int_var,  fc_int_var,
                                                             nc_int_var,  oec_int_var, ofc_int_var,
                                                             onc_int_var, osc_int_var, sc_int_var };
        std::vector<int> int_idxs;
        for (const auto& var : int_vars)
        {
            int_idxs.push_back(var_db->registerVariableAndContext(var, ctx, SAMRAIIntVector(1)));
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
