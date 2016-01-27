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
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <HierarchyDataOpsManager.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/muParserCartGridFunction.h>
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
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "vc_laplace.log");
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

        Pointer<SideVariable<NDIM, double> > u_side_var = new SideVariable<NDIM, double>("u_side");
        Pointer<SideVariable<NDIM, double> > f_side_var = new SideVariable<NDIM, double>("f_side");
        Pointer<SideVariable<NDIM, double> > e_side_var = new SideVariable<NDIM, double>("e_side");

        const int u_side_idx = var_db->registerVariableAndContext(u_side_var, ctx, IntVector<NDIM>(1));
        const int f_side_idx = var_db->registerVariableAndContext(f_side_var, ctx, IntVector<NDIM>(1));
        const int e_side_idx = var_db->registerVariableAndContext(e_side_var, ctx, IntVector<NDIM>(1));

        Pointer<CellVariable<NDIM, double> > u_cell_var = new CellVariable<NDIM, double>("u_cell", NDIM);
        Pointer<CellVariable<NDIM, double> > f_cell_var = new CellVariable<NDIM, double>("f_cell", NDIM);
        Pointer<CellVariable<NDIM, double> > e_cell_var = new CellVariable<NDIM, double>("e_cell", NDIM);

        const int u_cell_idx = var_db->registerVariableAndContext(u_cell_var, ctx, IntVector<NDIM>(0));
        const int f_cell_idx = var_db->registerVariableAndContext(f_cell_var, ctx, IntVector<NDIM>(0));
        const int e_cell_idx = var_db->registerVariableAndContext(e_cell_var, ctx, IntVector<NDIM>(0));

        Pointer<NodeVariable<NDIM, double> > mu_node_var = new NodeVariable<NDIM, double>("mu_node");

        const int mu_node_idx = var_db->registerVariableAndContext(mu_node_var, ctx, IntVector<NDIM>(1));

        // Register variables for plotting.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity(u_cell_var->getName(), "VECTOR", u_cell_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            ostringstream stream;
            stream << d;
            visit_data_writer->registerPlotQuantity(u_cell_var->getName() + stream.str(), "SCALAR", u_cell_idx, d);
        }

        visit_data_writer->registerPlotQuantity(f_cell_var->getName(), "VECTOR", f_cell_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            ostringstream stream;
            stream << d;
            visit_data_writer->registerPlotQuantity(f_cell_var->getName() + stream.str(), "SCALAR", f_cell_idx, d);
        }

        visit_data_writer->registerPlotQuantity(e_cell_var->getName(), "VECTOR", e_cell_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            ostringstream stream;
            stream << d;
            visit_data_writer->registerPlotQuantity(e_cell_var->getName() + stream.str(), "SCALAR", e_cell_idx, d);
        }

        visit_data_writer->registerPlotQuantity(mu_node_var->getName(), "SCALAR", mu_node_idx);

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

        // Set the simulation time to be zero.
        const double data_time = 0.0;

        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_side_idx, data_time);
            level->allocatePatchData(f_side_idx, data_time);
            level->allocatePatchData(e_side_idx, data_time);
            level->allocatePatchData(u_cell_idx, data_time);
            level->allocatePatchData(f_cell_idx, data_time);
            level->allocatePatchData(e_cell_idx, data_time);
            level->allocatePatchData(mu_node_idx, data_time);
        }

        // Setup exact solution data.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        muParserCartGridFunction f_fcn("f", app_initializer->getComponentDatabase("f"), grid_geometry);
        muParserCartGridFunction mu_fcn("mu", app_initializer->getComponentDatabase("mu"), grid_geometry);

        u_fcn.setDataOnPatchHierarchy(u_side_idx, u_side_var, patch_hierarchy, data_time);
        f_fcn.setDataOnPatchHierarchy(e_side_idx, e_side_var, patch_hierarchy, data_time);
        mu_fcn.setDataOnPatchHierarchy(mu_node_idx, mu_node_var, patch_hierarchy, data_time);

        // Create an object to communicate ghost cell data.
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        InterpolationTransactionComponent u_transaction(u_side_idx, "CUBIC_COARSEN", "LINEAR");
        InterpolationTransactionComponent mu_transaction(mu_node_idx, "CONSTANT_COARSEN", "LINEAR");
        vector<InterpolationTransactionComponent> transactions(2);
        transactions[0] = u_transaction;
        transactions[1] = mu_transaction;
        Pointer<HierarchyGhostCellInterpolation> bdry_fill_op = new HierarchyGhostCellInterpolation();
        bdry_fill_op->initializeOperatorState(transactions, patch_hierarchy);

        // Create the math operations object and get the patch data index for
        // the side-centered weighting factor.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int dx_side_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        // Compute (f0,f1) := div mu (grad(u0,u1) + grad(u0,u1)^T).
        hier_math_ops.vc_laplace(f_side_idx,
                                 f_side_var,
                                 1.0,
                                 0.0,
                                 mu_node_idx,
                                 mu_node_var,
                                 u_side_idx,
                                 u_side_var,
                                 bdry_fill_op,
                                 data_time);

        // Compute error and print error norms.
        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_side_data_ops =
            HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(u_side_var, patch_hierarchy, true);
        hier_side_data_ops->subtract(e_side_idx, e_side_idx, f_side_idx); // computes e := e - f
        pout << "|e|_oo = " << hier_side_data_ops->maxNorm(e_side_idx, dx_side_idx) << "\n";
        pout << "|e|_2  = " << hier_side_data_ops->L2Norm(e_side_idx, dx_side_idx) << "\n";
        pout << "|e|_1  = " << hier_side_data_ops->L1Norm(e_side_idx, dx_side_idx) << "\n";

        // Interpolate the side-centered data to cell centers for output.
        static const bool synch_cf_interface = true;
        hier_math_ops.interp(u_cell_idx, u_cell_var, u_side_idx, u_side_var, NULL, data_time, synch_cf_interface);
        hier_math_ops.interp(f_cell_idx, f_cell_var, f_side_idx, f_side_var, NULL, data_time, synch_cf_interface);
        hier_math_ops.interp(e_cell_idx, e_cell_var, e_side_idx, e_side_var, NULL, data_time, synch_cf_interface);

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, data_time);

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
} // main
