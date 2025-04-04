// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2024 by the IBAMR developers
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
#include <BergerRigoutsos.h>
#include <CartesianCellDoubleLinearRefine.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/SAMRAIScopedVectorDuplicate.h"
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/SAMRAIScopedVectorCopy.h>
#include <ibtk/muParserCartGridFunction.h>

// Set up application namespace declarations
#include <random>

#include <ibtk/app_namespaces.h>

// test stuff
#include "../tests.h"

// This file is the main driver for type conversion tests.

int
main(int argc, char** argv)
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "type_conversion.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

        // we don't want to use a conservative refinement scheme
        Pointer<RefineOperator<NDIM> > linear_refine = new CartesianCellDoubleLinearRefine<NDIM>();
        grid_geometry->addSpatialRefineOperator(linear_refine);

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

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        const int n_ghosts = 2;
        Pointer<CellVariable<NDIM, double> > u_cc_f64_var = new CellVariable<NDIM, double>("u_cc_f64", NDIM);
        const int u_cc_f64_idx = var_db->registerVariableAndContext(u_cc_f64_var, ctx, IntVector<NDIM>(n_ghosts));

        // Initialize the AMR patch hierarchy.
        const int tag_buffer = 1;
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_cc_f64_idx, 0.0);
        }

        // Setup exact solutions.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        u_fcn.setDataOnPatchHierarchy(u_cc_f64_idx, u_cc_f64_var, patch_hierarchy, 0.0);

        // Fill ghost data.
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_components(1);
        ghost_cell_components[0] = ITC(u_cc_f64_idx,
                                       "LINEAR_REFINE",
                                       true,
                                       "CONSERVATIVE_COARSEN",
                                       "LINEAR",
                                       false,
                                       {}, // u_bc_coefs
                                       nullptr);
        HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ghost_cell_components, patch_hierarchy);
        ghost_fill_op.fillData(/*time*/ 0.0);

        HierarchyMathOps math_ops("HierarchyMathOps", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        int cc_f64_wgt_idx = math_ops.getCellWeightPatchDescriptorIndex();

        Pointer<SAMRAIVectorReal<NDIM, double> > u_cc_f64_vector =
            new SAMRAIVectorReal<NDIM, double>("u_vector", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        u_cc_f64_vector->addComponent(u_cc_f64_var, u_cc_f64_idx, cc_f64_wgt_idx);

        SAMRAIScopedVectorDuplicate<float> u_cc_f32_duplicated_vector(u_cc_f64_vector, "u_cc_f32_duplicated_vector");
        Pointer<CellVariable<NDIM, float> > u_cc_f32_duplicated_var =
            u_cc_f32_duplicated_vector.getComponentVectors()[0]->getComponentVariable(0);
        TBOX_ASSERT(u_cc_f32_duplicated_var);

        // Convert from f64 to f32:
        SAMRAIScopedVectorCopy<float> u_cc_f32_copied_vector(u_cc_f64_vector, "u_cc_f32_copied_vector");
        Pointer<CellVariable<NDIM, float> > u_cc_f32_copied_var =
            u_cc_f32_copied_vector.getComponentVectors()[0]->getComponentVariable(0);
        TBOX_ASSERT(u_cc_f32_copied_var);
        Pointer<SAMRAIVectorReal<NDIM, float> > u_cc_f32_copied_vector_ptr = u_cc_f32_copied_vector;

        // Convert back from f32 to f64; we expect to get loss of precision (evaluated below):
        SAMRAIScopedVectorCopy<double> u_cc_f64_copied_vector(u_cc_f32_copied_vector_ptr, "u_cc_f64_copied_vector");
        Pointer<CellVariable<NDIM, double> > u_cc_f64_copied_var =
            u_cc_f64_copied_vector.getComponentVectors()[0]->getComponentVariable(0);
        TBOX_ASSERT(u_cc_f64_copied_var);
        int u_cc_f64_copied_idx = u_cc_f64_copied_vector.getComponentVectors()[0]->getComponentDescriptorIndex(0);

        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity(u_cc_f64_var->getName(), "VECTOR", u_cc_f64_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(
                u_cc_f64_var->getName() + std::to_string(d), "SCALAR", u_cc_f64_idx, d);
        }

        visit_data_writer->registerPlotQuantity("u_f64_copy", "VECTOR", u_cc_f64_copied_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("u_f64_copy" + std::to_string(d), "SCALAR", u_cc_f64_copied_idx, d);
        }

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

        // This is the real test:
        double u_norm = u_cc_f64_vector->L2Norm();
        Pointer<SAMRAIVectorReal<NDIM, double> > u_cc_f64_copied_vector_ptr = u_cc_f64_copied_vector;
        u_cc_f64_vector->subtract(u_cc_f64_copied_vector_ptr, u_cc_f64_vector);
        double diff_norm = u_cc_f64_vector->L2Norm();

        std::ofstream out("output");
        out.precision(16);
        out << "absolute difference: " << diff_norm << "\n";
        out << "relative difference: " << diff_norm / u_norm << "\n";

    } // cleanup dynamically allocated objects prior to shutdown
} // main
