// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Config files
#include <IBTK_config.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <CellVariable.h>
#include <EdgeVariable.h>
#include <FaceVariable.h>
#include <NodeVariable.h>
#include <OuteredgeVariable.h>
#include <OuterfaceVariable.h>
#include <OuternodeVariable.h>
#include <OutersideVariable.h>
#include <SideVariable.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/SAMRAIDataCache.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

void
print_patch_descriptor_data(Pointer<PatchDescriptor<NDIM> > descriptor)
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
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
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

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        pout << "adding double variables\n";
        Pointer<CellVariable<NDIM, double> > cc_double_var = new CellVariable<NDIM, double>("cc_double");
        Pointer<EdgeVariable<NDIM, double> > ec_double_var = new EdgeVariable<NDIM, double>("ec_double");
        Pointer<FaceVariable<NDIM, double> > fc_double_var = new FaceVariable<NDIM, double>("fc_double");
        Pointer<NodeVariable<NDIM, double> > nc_double_var = new NodeVariable<NDIM, double>("nc_double");
        Pointer<OuteredgeVariable<NDIM, double> > oec_double_var = new OuteredgeVariable<NDIM, double>("oec_double");
        Pointer<OuterfaceVariable<NDIM, double> > ofc_double_var = new OuterfaceVariable<NDIM, double>("ofc_double");
        Pointer<OuternodeVariable<NDIM, double> > onc_double_var = new OuternodeVariable<NDIM, double>("onc_double");
        Pointer<OutersideVariable<NDIM, double> > osc_double_var = new OutersideVariable<NDIM, double>("osc_double");
        Pointer<SideVariable<NDIM, double> > sc_double_var = new SideVariable<NDIM, double>("sc_double");
        std::vector<Pointer<Variable<NDIM> > > double_vars{ cc_double_var,  ec_double_var,  fc_double_var,
                                                            nc_double_var,  oec_double_var, ofc_double_var,
                                                            onc_double_var, osc_double_var, sc_double_var };
        std::vector<int> double_idxs;
        for (const auto& var : double_vars)
        {
            double_idxs.push_back(var_db->registerVariableAndContext(var, ctx, IntVector<NDIM>(1)));
        }
        print_patch_descriptor_data(patch_hierarchy->getPatchDescriptor());

        pout << "adding float variables\n";
        Pointer<CellVariable<NDIM, float> > cc_float_var = new CellVariable<NDIM, float>("cc_float");
        Pointer<EdgeVariable<NDIM, float> > ec_float_var = new EdgeVariable<NDIM, float>("ec_float");
        Pointer<FaceVariable<NDIM, float> > fc_float_var = new FaceVariable<NDIM, float>("fc_float");
        Pointer<NodeVariable<NDIM, float> > nc_float_var = new NodeVariable<NDIM, float>("nc_float");
        Pointer<OuteredgeVariable<NDIM, float> > oec_float_var = new OuteredgeVariable<NDIM, float>("oec_float");
        Pointer<OuterfaceVariable<NDIM, float> > ofc_float_var = new OuterfaceVariable<NDIM, float>("ofc_float");
        Pointer<OuternodeVariable<NDIM, float> > onc_float_var = new OuternodeVariable<NDIM, float>("onc_float");
        Pointer<OutersideVariable<NDIM, float> > osc_float_var = new OutersideVariable<NDIM, float>("osc_float");
        Pointer<SideVariable<NDIM, float> > sc_float_var = new SideVariable<NDIM, float>("sc_float");

        std::vector<Pointer<Variable<NDIM> > > float_vars{ cc_float_var,  ec_float_var,  fc_float_var,
                                                           nc_float_var,  oec_float_var, ofc_float_var,
                                                           onc_float_var, osc_float_var, sc_float_var };
        std::vector<int> float_idxs;
        for (const auto& var : float_vars)
        {
            float_idxs.push_back(var_db->registerVariableAndContext(var, ctx, IntVector<NDIM>(1)));
        }
        print_patch_descriptor_data(patch_hierarchy->getPatchDescriptor());

        pout << "adding int variables\n";
        Pointer<CellVariable<NDIM, int> > cc_int_var = new CellVariable<NDIM, int>("cc_int");
        Pointer<EdgeVariable<NDIM, int> > ec_int_var = new EdgeVariable<NDIM, int>("ec_int");
        Pointer<FaceVariable<NDIM, int> > fc_int_var = new FaceVariable<NDIM, int>("fc_int");
        Pointer<NodeVariable<NDIM, int> > nc_int_var = new NodeVariable<NDIM, int>("nc_int");
        Pointer<OuteredgeVariable<NDIM, int> > oec_int_var = new OuteredgeVariable<NDIM, int>("oec_int");
        Pointer<OuterfaceVariable<NDIM, int> > ofc_int_var = new OuterfaceVariable<NDIM, int>("ofc_int");
        Pointer<OuternodeVariable<NDIM, int> > onc_int_var = new OuternodeVariable<NDIM, int>("onc_int");
        Pointer<OutersideVariable<NDIM, int> > osc_int_var = new OutersideVariable<NDIM, int>("osc_int");
        Pointer<SideVariable<NDIM, int> > sc_int_var = new SideVariable<NDIM, int>("sc_int");

        std::vector<Pointer<Variable<NDIM> > > int_vars{ cc_int_var,  ec_int_var,  fc_int_var,  nc_int_var, oec_int_var,
                                                         ofc_int_var, onc_int_var, osc_int_var, sc_int_var };
        std::vector<int> int_idxs;
        for (const auto& var : int_vars)
        {
            int_idxs.push_back(var_db->registerVariableAndContext(var, ctx, IntVector<NDIM>(1)));
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

    SAMRAIManager::shutdown();
    PetscFinalize();
} // main
