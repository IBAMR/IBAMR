// Created on 06 Sep 2017 by Amneet Bhalla
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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
#include <ibtk/PETScKrylovPoissonSolver.h>
#include <ibtk/SCLaplaceOperator.h>
#include <ibtk/SCPoissonSolverManager.h>
#include <ibtk/VCSCViscousOpPointRelaxationFACOperator.h>
#include <ibtk/VCSCViscousOperator.h>
#include <ibtk/VCSCViscousPETScLevelSolver.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

Pointer<PoissonSolver>
allocate_vc_velocity_krylov_solver(const std::string& solver_object_name,
                                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> solver_input_db,
                                   const std::string& solver_default_options_prefix)
{
    Pointer<PETScKrylovPoissonSolver> krylov_solver =
        new PETScKrylovPoissonSolver(solver_object_name, solver_input_db, solver_default_options_prefix);
    krylov_solver->setOperator(new VCSCViscousOperator(solver_object_name + "::vc_viscous_operator"));
    return krylov_solver;
}

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
bool
run_example(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "sc_poisson.log");
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

        Pointer<SideVariable<NDIM, double> > u_sc_var = new SideVariable<NDIM, double>("u_sc");
        Pointer<SideVariable<NDIM, double> > f_sc_var = new SideVariable<NDIM, double>("f_sc");
        Pointer<SideVariable<NDIM, double> > e_sc_var = new SideVariable<NDIM, double>("e_sc");
        Pointer<SideVariable<NDIM, double> > r_sc_var = new SideVariable<NDIM, double>("r_sc");
#if (NDIM == 2)
        Pointer<NodeVariable<NDIM, double> > mu_nc_var = new NodeVariable<NDIM, double>("mu_node");
        const int mu_nc_idx = var_db->registerVariableAndContext(mu_nc_var, ctx, IntVector<NDIM>(1));
#elif (NDIM == 3)
        Pointer<EdgeVariable<NDIM, double> > mu_ec_var = new EdgeVariable<NDIM, double>("mu_edge");
        const int mu_ec_idx = var_db->registerVariableAndContext(mu_ec_var, ctx, IntVector<NDIM>(1));
        Pointer<CellVariable<NDIM, double> > mu_cc_var = new CellVariable<NDIM, double>("mu_cc");
        const int mu_cc_idx = var_db->registerVariableAndContext(mu_cc_var, ctx, IntVector<NDIM>(0));
#endif

        const int u_sc_idx = var_db->registerVariableAndContext(u_sc_var, ctx, IntVector<NDIM>(1));
        const int f_sc_idx = var_db->registerVariableAndContext(f_sc_var, ctx, IntVector<NDIM>(1));
        const int e_sc_idx = var_db->registerVariableAndContext(e_sc_var, ctx, IntVector<NDIM>(1));
        const int r_sc_idx = var_db->registerVariableAndContext(r_sc_var, ctx, IntVector<NDIM>(1));

        Pointer<CellVariable<NDIM, double> > u_cc_var = new CellVariable<NDIM, double>("u_cc", NDIM);
        Pointer<CellVariable<NDIM, double> > f_cc_var = new CellVariable<NDIM, double>("f_cc", NDIM);
        Pointer<CellVariable<NDIM, double> > e_cc_var = new CellVariable<NDIM, double>("e_cc", NDIM);
        Pointer<CellVariable<NDIM, double> > r_cc_var = new CellVariable<NDIM, double>("r_cc", NDIM);

        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, IntVector<NDIM>(0));
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, IntVector<NDIM>(0));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVector<NDIM>(0));
        const int r_cc_idx = var_db->registerVariableAndContext(r_cc_var, ctx, IntVector<NDIM>(0));

        // Register variables for plotting.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity(u_cc_var->getName(), "VECTOR", u_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            ostringstream stream;
            stream << d;
            visit_data_writer->registerPlotQuantity(u_cc_var->getName() + stream.str(), "SCALAR", u_cc_idx, d);
        }

        visit_data_writer->registerPlotQuantity(f_cc_var->getName(), "VECTOR", f_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            ostringstream stream;
            stream << d;
            visit_data_writer->registerPlotQuantity(f_cc_var->getName() + stream.str(), "SCALAR", f_cc_idx, d);
        }

        visit_data_writer->registerPlotQuantity(e_cc_var->getName(), "VECTOR", e_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            ostringstream stream;
            stream << d;
            visit_data_writer->registerPlotQuantity(e_cc_var->getName() + stream.str(), "SCALAR", e_cc_idx, d);
        }

        visit_data_writer->registerPlotQuantity(r_cc_var->getName(), "VECTOR", r_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            ostringstream stream;
            stream << d;
            visit_data_writer->registerPlotQuantity(r_cc_var->getName() + stream.str(), "SCALAR", r_cc_idx, d);
        }

#if (NDIM == 2)
        visit_data_writer->registerPlotQuantity(mu_nc_var->getName(), "SCALAR", mu_nc_idx);
#elif (NDIM == 3)
        visit_data_writer->registerPlotQuantity(mu_cc_var->getName(), "SCALAR", mu_cc_idx);
#endif

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
            level->allocatePatchData(u_sc_idx, 0.0);
            level->allocatePatchData(f_sc_idx, 0.0);
            level->allocatePatchData(e_sc_idx, 0.0);
            level->allocatePatchData(r_sc_idx, 0.0);
#if (NDIM == 2)
            level->allocatePatchData(mu_nc_idx, 0.0);
#elif (NDIM == 3)
            level->allocatePatchData(mu_ec_idx, 0.0);
            level->allocatePatchData(mu_cc_idx, 0.0);
#endif
            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
            level->allocatePatchData(r_cc_idx, 0.0);
        }

        // Setup vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int h_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> r_vec("r", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        u_vec.addComponent(u_sc_var, u_sc_idx, h_sc_idx);
        f_vec.addComponent(f_sc_var, f_sc_idx, h_sc_idx);
        e_vec.addComponent(e_sc_var, e_sc_idx, h_sc_idx);
        r_vec.addComponent(r_sc_var, r_sc_idx, h_sc_idx);

        u_vec.setToScalar(0.0);
        f_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);
        r_vec.setToScalar(0.0);

        // Setup exact solutions.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        muParserCartGridFunction f_fcn("f", app_initializer->getComponentDatabase("f"), grid_geometry);
        muParserCartGridFunction mu_fcn("mu", app_initializer->getComponentDatabase("mu"), grid_geometry);

        u_fcn.setDataOnPatchHierarchy(e_sc_idx, e_sc_var, patch_hierarchy, 0.0);
        f_fcn.setDataOnPatchHierarchy(f_sc_idx, f_sc_var, patch_hierarchy, 0.0);
#if (NDIM == 2)
        mu_fcn.setDataOnPatchHierarchy(mu_nc_idx, mu_nc_var, patch_hierarchy, 0.0);
#elif (NDIM == 3)
        mu_fcn.setDataOnPatchHierarchy(mu_ec_idx, mu_ec_var, patch_hierarchy, 0.0);
#endif
        // Fill ghost cells of viscosity.
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> transaction_comp(1);
#if (NDIM == 2)
        transaction_comp[0] = InterpolationTransactionComponent(mu_nc_idx,
                                                                /*DATA_REFINE_TYPE*/ "LINEAR_REFINE",
                                                                /*USE_CF_INTERPOLATION*/ false,
                                                                /*DATA_COARSEN_TYPE*/ "CONSTANT_COARSEN",
                                                                /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                                /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                                /*mu_bc_coef*/ NULL,
                                                                Pointer<VariableFillPattern<NDIM> >(NULL));
#elif (NDIM == 3)
        transaction_comp[0] = InterpolationTransactionComponent(mu_ec_idx,
                                                                /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                                /*USE_CF_INTERPOLATION*/ false,
                                                                /*DATA_COARSEN_TYPE*/ "CONSERVATIVE_COARSEN",
                                                                /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                                /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                                /*mu_bc_coef*/ NULL,
                                                                Pointer<VariableFillPattern<NDIM> >(NULL));
#endif
        Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_bdry_fill->initializeOperatorState(transaction_comp, patch_hierarchy);
        hier_bdry_fill->setHomogeneousBc(false);
        hier_bdry_fill->fillData(/*time*/ 0.0);

        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > > U_nul_vecs(NDIM);
        const bool has_velocity_nullspace = periodic_shift.min() > 0;
        if (has_velocity_nullspace)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }

            for (unsigned int k = 0; k < NDIM; ++k)
            {
                std::ostringstream stream;
                stream << k;
                U_nul_vecs[k] = f_vec.cloneVector("nul_vec_U_" + stream.str());
                U_nul_vecs[k]->allocateVectorData(/*0.0*/);
                U_nul_vecs[k]->setToScalar(0.0);
                for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
                {
                    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        Pointer<SideData<NDIM, double> > U_nul_data =
                            patch->getPatchData(U_nul_vecs[k]->getComponentDescriptorIndex(0));
                        U_nul_data->getArrayData(k).fillAll(1.0);
                    }
                }
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                Pointer<Database> bc_coefs_db = app_initializer->getComponentDatabase(bc_coefs_db_name);
                u_bc_coefs[d] = new muParserRobinBcCoefs(bc_coefs_name, bc_coefs_db, grid_geometry);
            }
        }

        // Setup the implicit viscous solver.
        PoissonSpecifications vc_vel_spec("vc_vel_spec");
        vc_vel_spec.setCConstant(input_db->getDouble("C"));
#if (NDIM == 2)
        vc_vel_spec.setDPatchDataId(mu_nc_idx);
#elif (NDIM == 3)
        vc_vel_spec.setDPatchDataId(mu_ec_idx);
#endif
        VCSCViscousOperator viscous_op("viscous_op");
        viscous_op.setPoissonSpecifications(vc_vel_spec);
        viscous_op.setPhysicalBcCoefs(u_bc_coefs);
        viscous_op.initializeOperatorState(u_vec, f_vec);
        viscous_op.setSolutionTime(0.0);

        string solver_type = input_db->getString("solver_type");
        Pointer<Database> solver_db = input_db->getDatabase("solver_db");
        string precond_type = input_db->getString("precond_type");
        Pointer<Database> precond_db = input_db->getDatabase("precond_db");

        SCPoissonSolverManager* solver_manager = SCPoissonSolverManager::getManager();
        solver_manager->registerSolverFactoryFunction("VC_VELOCITY_PETSC_KRYLOV_SOLVER",
                                                      allocate_vc_velocity_krylov_solver);
        solver_manager->registerSolverFactoryFunction("VC_VELOCITY_POINT_RELAXATION_FAC_PRECONDITIONER",
                                                      VCSCViscousOpPointRelaxationFACOperator::allocate_solver);
        solver_manager->registerSolverFactoryFunction("VC_VELOCITY_PETSC_LEVEL_SOLVER",
                                                      VCSCViscousPETScLevelSolver::allocate_solver);
        Pointer<PoissonSolver> poisson_solver = solver_manager->allocateSolver(
            solver_type, "vc_velocity_solver", solver_db, "vc_velocity_", precond_type, "vc_velocity_precond",
           precond_db, "vc_velocity_pc_");
        /*Pointer<PoissonSolver> poisson_solver = solver_manager->allocateSolver(solver_type,
                                                                               "vc_velocity_solver",
                                                                               solver_db,
                                                                               "vc_velocity_",
                                                                               precond_type,
                                                                               "vc_velocity_precond",
                                                                               precond_db,
                                                                               "vc_velocity_fac_pc");*/

        poisson_solver->setPoissonSpecifications(vc_vel_spec);
        poisson_solver->setPhysicalBcCoefs(u_bc_coefs);
        poisson_solver->setSolutionTime(0.0);

        // Remove nullspace from RHS
        HierarchySideDataOpsReal<NDIM, double> sc_data_ops(patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        if (has_velocity_nullspace)
        {
            // Ensure that the right-hand-side vector has no components in the
            // nullspace of the operator.
            for (unsigned int k = 0; k < NDIM; ++k)
            {
                r_vec.setToScalar(0.0);
                for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
                {
                    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        Pointer<SideData<NDIM, double> > r_data =
                            patch->getPatchData(r_vec.getComponentDescriptorIndex(0));
                        r_data->getArrayData(k).fillAll(1.0);
                    }
                }

                double fk_mean = f_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double> >(&r_vec, false)) /
                                 r_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double> >(&r_vec, false));
                sc_data_ops.axpy(f_sc_idx, -fk_mean, r_sc_idx, f_sc_idx);
            }

            Pointer<PETScKrylovLinearSolver> krylov_solver = poisson_solver;
            if (krylov_solver)
            {
                krylov_solver->setNullspace(false, U_nul_vecs);
            }
        }

        // Solve L*u = f.
        u_vec.setToScalar(0.0);
        poisson_solver->initializeSolverState(u_vec, f_vec);
        poisson_solver->solveSystem(u_vec, f_vec);

        // Subtract mean velocity from velocity components.
        if (has_velocity_nullspace)
        {
            for (unsigned int k = 0; k < NDIM; ++k)
            {
                r_vec.setToScalar(0.0);
                for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
                {
                    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        Pointer<SideData<NDIM, double> > r_data =
                            patch->getPatchData(r_vec.getComponentDescriptorIndex(0));
                        r_data->getArrayData(k).fillAll(1.0);
                    }
                }

                const double uk_mean = u_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double> >(&r_vec, false)) /
                                       r_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double> >(&r_vec, false));
                sc_data_ops.axpy(u_sc_idx, -uk_mean, r_sc_idx, u_sc_idx);
            }
        }

        // Compute error and print error norms.
        e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&e_vec, false),
                       Pointer<SAMRAIVectorReal<NDIM, double> >(&u_vec, false));
        pout << "|e|_oo = " << e_vec.maxNorm() << "\n";
        pout << "|e|_2  = " << e_vec.L2Norm() << "\n";
        pout << "|e|_1  = " << e_vec.L1Norm() << "\n";

        // Compute the residual and print residual norms.
        viscous_op.setHomogeneousBc(false);
        viscous_op.apply(u_vec, r_vec);
        r_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&f_vec, false),
                       Pointer<SAMRAIVectorReal<NDIM, double> >(&r_vec, false));
        pout << "|r|_oo = " << r_vec.maxNorm() << "\n";
        pout << "|r|_2  = " << r_vec.L2Norm() << "\n";
        pout << "|r|_1  = " << r_vec.L1Norm() << "\n";

        // Interpolate the side-centered data to cell centers for output.
        static const bool synch_cf_interface = true;
        hier_math_ops.interp(u_cc_idx, u_cc_var, u_sc_idx, u_sc_var, NULL, 0.0, synch_cf_interface);
        hier_math_ops.interp(f_cc_idx, f_cc_var, f_sc_idx, f_sc_var, NULL, 0.0, synch_cf_interface);
        hier_math_ops.interp(e_cc_idx, e_cc_var, e_sc_idx, e_sc_var, NULL, 0.0, synch_cf_interface);
        hier_math_ops.interp(r_cc_idx, r_cc_var, r_sc_idx, r_sc_var, NULL, 0.0, synch_cf_interface);

#if (NDIM == 3)
        // Interpolate the edge-centered data to cell centers for output.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& box = patch->getBox();
                Pointer<EdgeData<NDIM, double> > mu_ec_data = patch->getPatchData(mu_ec_idx);
                Pointer<CellData<NDIM, double> > mu_cc_data = patch->getPatchData(mu_cc_idx);
                for (Box<NDIM>::Iterator it(CellGeometry<NDIM>::toCellBox(box)); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    Box<NDIM> edge_box(ci, ci);
                    double avg_mu = 0.0;
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for(EdgeIterator<NDIM> e(edge_box, axis); e; e++)
                        {
                            EdgeIndex<NDIM> ei(e());
                            avg_mu += (*mu_ec_data)(ei);
                        }
                    }
                    (*mu_cc_data)(ci) = avg_mu/12.0;
                }
            }
        }
#endif
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
    return true;
} // run_example
