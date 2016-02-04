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

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CCPoissonSolverManager.h>
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

        Pointer<CellVariable<NDIM, double> > u_cc_var = new CellVariable<NDIM, double>("u_cc");
        Pointer<CellVariable<NDIM, double> > f_cc_var = new CellVariable<NDIM, double>("f_cc");
        Pointer<CellVariable<NDIM, double> > e_cc_var = new CellVariable<NDIM, double>("e_cc");
        Pointer<CellVariable<NDIM, double> > r_cc_var = new CellVariable<NDIM, double>("r_cc");

        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, IntVector<NDIM>(1));
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, IntVector<NDIM>(1));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVector<NDIM>(1));
        const int r_cc_idx = var_db->registerVariableAndContext(r_cc_var, ctx, IntVector<NDIM>(1));

        // Register variables for plotting.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity(u_cc_var->getName(), "SCALAR", u_cc_idx);
        visit_data_writer->registerPlotQuantity(f_cc_var->getName(), "SCALAR", f_cc_idx);
        visit_data_writer->registerPlotQuantity(e_cc_var->getName(), "SCALAR", e_cc_idx);
        visit_data_writer->registerPlotQuantity(r_cc_var->getName(), "SCALAR", r_cc_idx);

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

        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
            level->allocatePatchData(r_cc_idx, 0.0);
        }

        // Setup vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int h_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> r_vec("r", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        u_vec.addComponent(u_cc_var, u_cc_idx, h_cc_idx);
        f_vec.addComponent(f_cc_var, f_cc_idx, h_cc_idx);
        e_vec.addComponent(e_cc_var, e_cc_idx, h_cc_idx);
        r_vec.addComponent(r_cc_var, r_cc_idx, h_cc_idx);

        u_vec.setToScalar(0.0);
        f_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);
        r_vec.setToScalar(1.0);

        // Setup exact solutions.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        muParserCartGridFunction f_fcn("f", app_initializer->getComponentDatabase("f"), grid_geometry);

        u_fcn.setDataOnPatchHierarchy(e_cc_idx, e_cc_var, patch_hierarchy, 0.0);
        f_fcn.setDataOnPatchHierarchy(f_cc_idx, f_cc_var, patch_hierarchy, 0.0);

        // Ensure that the right-hand-side vector has no components in the
        // nullspace of the operator.
        f_vec.addScalar(Pointer<SAMRAIVectorReal<NDIM, double> >(&f_vec, false),
                        -f_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double> >(&r_vec, false)) /
                            r_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double> >(&r_vec, false)));

        // Setup the Poisson solver.
        PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setCZero();
        poisson_spec.setDConstant(-1.0);
        RobinBcCoefStrategy<NDIM>* bc_coef = NULL;
        CCLaplaceOperator laplace_op("laplace_op");
        laplace_op.setPoissonSpecifications(poisson_spec);
        laplace_op.setPhysicalBcCoef(bc_coef);
        laplace_op.initializeOperatorState(u_vec, f_vec);

        string solver_type = input_db->getString("solver_type");
        Pointer<Database> solver_db = input_db->getDatabase("solver_db");
        string precond_type = input_db->getString("precond_type");
        Pointer<Database> precond_db = input_db->getDatabase("precond_db");
        Pointer<PoissonSolver> poisson_solver = CCPoissonSolverManager::getManager()->allocateSolver(
            solver_type, "poisson_solver", solver_db, "", precond_type, "poisson_precond", precond_db, "");
        poisson_solver->setPoissonSpecifications(poisson_spec);
        poisson_solver->setPhysicalBcCoef(bc_coef);
        poisson_solver->initializeSolverState(u_vec, f_vec);

        // Solve -L*u = f.
        u_vec.setToScalar(0.0);
        poisson_solver->solveSystem(u_vec, f_vec);

        // Compute error and print error norms.
        e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&e_vec, false),
                       Pointer<SAMRAIVectorReal<NDIM, double> >(&u_vec, false));
        pout << "|e|_oo = " << e_vec.maxNorm() << "\n";
        pout << "|e|_2  = " << e_vec.L2Norm() << "\n";
        pout << "|e|_1  = " << e_vec.L1Norm() << "\n";

        // Compute the residual and print residual norms.
        laplace_op.apply(u_vec, r_vec);
        r_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&f_vec, false),
                       Pointer<SAMRAIVectorReal<NDIM, double> >(&r_vec, false));
        pout << "|r|_oo = " << r_vec.maxNorm() << "\n";
        pout << "|r|_2  = " << r_vec.L2Norm() << "\n";
        pout << "|r|_1  = " << r_vec.L1Norm() << "\n";

        // Set invalid values on coarse levels (i.e., coarse-grid values that
        // are covered by finer grid patches) to equal zero.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber() - 1; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            BoxArray<NDIM> refined_region_boxes;
            Pointer<PatchLevel<NDIM> > next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, double> > e_cc_data = patch->getPatchData(e_cc_idx);
                Pointer<CellData<NDIM, double> > r_cc_data = patch->getPatchData(r_cc_idx);
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const Box<NDIM> refined_box = refined_region_boxes[i];
                    const Box<NDIM> intersection = Box<NDIM>::grow(patch_box, 1) * refined_box;
                    if (!intersection.empty())
                    {
                        e_cc_data->fillAll(0.0, intersection);
                        r_cc_data->fillAll(0.0, intersection);
                    }
                }
            }
        }

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
} // main
