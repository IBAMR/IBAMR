// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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
#include <petsc.h>

// Headers for basic SAMRAI objects
#include <PatchLevel.h>
#include <VariableDatabase.h>
#include <tbox/Database.h>
#include <tbox/InputDatabase.h>
#include <tbox/InputManager.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>
#include <tbox/RestartManager.h>
#include <tbox/SAMRAIManager.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// Headers for major algorithm/data structure objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <StandardTagAndInitialize.h>
#include <VisItDataWriter.h>

// Headers for application-specific algorithm/data structure objects
#include <CellVariable.h>
#include <HierarchyDataOpsManager.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <Variable.h>
#include <VariableDatabase.h>

#include <ibtk/CCPoissonFACOperator.h>
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CCLaplaceOperator2.h>
#include <ibtk/FACPreconditioner.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/NormOps.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PETScNewtonKrylovSolver.h>

#include "PoissonTester.h"

using namespace SAMRAI;
using namespace IBTK;
using namespace std;

/************************************************************************
 *                                                                      *
 * For each run, the input filename must be given on the command line.  *
 * In all cases, the command line is:                                   *
 *                                                                      *
 *    executable <input file name> <PETSc options>                      *
 *                                                                      *
 ************************************************************************
 */

int
main(
    int argc,
    char *argv[])
{
    /*
     * Initialize PETSc, MPI, and SAMRAI.
     */
    PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
    tbox::SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    tbox::SAMRAIManager::startup();

    {// ensure all smart Pointers are properly deleted
        string input_filename;
        input_filename = argv[1];

        tbox::plog << "input_filename = " << input_filename << endl;

        /*
         * Create input database and parse all data in input file.
         */
        tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
        tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

        /*
         * Retrieve "Main" section of the input database.
         */
        tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

        string log_file_name = "fac_test.log";
        if (main_db->keyExists("log_file_name"))
        {
            log_file_name = main_db->getString("log_file_name");
        }
        bool log_all_nodes = false;
        if (main_db->keyExists("log_all_nodes"))
        {
            log_all_nodes = main_db->getBool("log_all_nodes");
        }
        if (log_all_nodes)
        {
            tbox::PIO::logAllNodes(log_file_name);
        }
        else
        {
            tbox::PIO::logOnlyNodeZero(log_file_name);
        }

        bool viz_dump_enabled = false;
        tbox::Array<string> viz_writer(1);
        viz_writer[0] = "VisIt";
        string viz_dump_filename;
        string visit_dump_dirname;
        bool uses_visit = false;
        int visit_number_procs_per_file = 1;
        if (main_db->keyExists("viz_dump_enabled"))
        {
            viz_dump_enabled = main_db->getBool("viz_dump_enabled");
        }
        if (viz_dump_enabled)
        {
            if (main_db->keyExists("viz_writer"))
            {
                viz_writer = main_db->getStringArray("viz_writer");
            }
            if (main_db->keyExists("viz_dump_filename"))
            {
                viz_dump_filename = main_db->getString("viz_dump_filename");
            }
            string viz_dump_dirname;
            if (main_db->keyExists("viz_dump_dirname"))
            {
                viz_dump_dirname = main_db->getString("viz_dump_dirname");
            }
            for (int i = 0; i < viz_writer.getSize(); ++i)
            {
                if (viz_writer[i] == "VisIt") uses_visit = true;
            }
            if (uses_visit)
            {
                visit_dump_dirname = viz_dump_dirname;
            }
            else
            {
                TBOX_ERROR("main(): "
                           << "\nUnrecognized 'viz_writer' entry..."
                           << "\nOnly valid option is 'VisIt'"
                           << endl);
            }
            if (uses_visit)
            {
                if (viz_dump_dirname.empty())
                {
                    TBOX_ERROR("main(): "
                               << "\nviz_dump_dirname is null ... "
                               << "\nThis must be specified for use with VisIt"
                               << endl);
                }
                if (main_db->keyExists("visit_number_procs_per_file"))
                {
                    visit_number_procs_per_file =
                        main_db->getInteger("visit_number_procs_per_file");
                }
            }
        }

        const bool viz_dump_data = viz_dump_enabled;

        bool timer_enabled = false;
        if (main_db->keyExists("timer_enabled"))
        {
            timer_enabled = main_db->getBool("timer_enabled");
        }

        const bool write_timer_data = timer_enabled;

        if (write_timer_data)
        {
            tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
        }

        /*
         * Create major algorithm and data objects which comprise application.
         * Each object will be initialized either from input data or restart
         * files, or a combination of both.  Refer to each class constructor for
         * details.  For more information on the composition of objects for this
         * application, see comments at top of file.
         */
        tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry =
            new geom::CartesianGridGeometry<NDIM>("CartesianGeometry",
                                                  input_db->getDatabase("CartesianGeometry"));

        tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
            new hier::PatchHierarchy<NDIM>("PatchHierarchy",grid_geometry);

        tbox::Pointer<PoissonTester> poisson_tester = new PoissonTester();

        tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector =
            new mesh::StandardTagAndInitialize<NDIM>(
                "StandardTagAndInitialize",
                poisson_tester,
                input_db->getDatabase("StandardTagAndInitialize"));

        tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator = new mesh::BergerRigoutsos<NDIM>();

        tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer =
            new mesh::LoadBalancer<NDIM>("LoadBalancer",
                                         input_db->getDatabase("LoadBalancer"));

        tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > gridding_algorithm =
            new mesh::GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                              input_db->getDatabase("GriddingAlgorithm"),
                                              error_detector,
                                              box_generator,
                                              load_balancer);

        /*
         * Set up visualization plot file writer.
         */
        tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer;

        if (uses_visit)
        {
            visit_data_writer =
                new appu::VisItDataWriter<NDIM>("VisIt Writer",
                                                visit_dump_dirname,
                                                visit_number_procs_per_file);

            visit_data_writer->registerPlotQuantity("U", "SCALAR", poisson_tester->getUIndex());
            visit_data_writer->registerPlotQuantity("U_true - U", "SCALAR", poisson_tester->getVIndex());
            visit_data_writer->registerPlotQuantity("F", "SCALAR", poisson_tester->getFIndex());
        }

        /*
         * Initialize hierarchy configuration and data on all patches.
         */
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);

        int tag_buffer = 1;
        int level_number = 0;
        bool done = false;
        while (!done &&
               (gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->
                makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);

            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        /*
         * Initialize the initial guess for the solution and the right-hand
         * side.
         */
        const int U_idx = poisson_tester->getUIndex();
        const int V_idx = poisson_tester->getVIndex();
        const int F_idx = poisson_tester->getFIndex();

        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);

        const int h_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        solv::SAMRAIVectorReal<NDIM,double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        solv::SAMRAIVectorReal<NDIM,double> v_vec("v", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        solv::SAMRAIVectorReal<NDIM,double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
        tbox::Pointer<hier::Variable<NDIM> > var;

        var_db->mapIndexToVariable(U_idx, var);
        tbox::Pointer<pdat::CellVariable<NDIM,double> > U_var = var;
        u_vec.addComponent(U_var, U_idx, h_idx);

        var_db->mapIndexToVariable(V_idx, var);
        tbox::Pointer<pdat::CellVariable<NDIM,double> > V_var = var;
        v_vec.addComponent(V_var, V_idx, h_idx);

        var_db->mapIndexToVariable(F_idx, var);
        tbox::Pointer<pdat::CellVariable<NDIM,double> > F_var = var;
        f_vec.addComponent(F_var, F_idx, h_idx);

        math::HierarchyDataOpsManager<NDIM>* hier_ops_manager = math::HierarchyDataOpsManager<NDIM>::getManager();
        tbox::Pointer<math::HierarchyCellDataOpsReal<NDIM,double> > hier_cc_data_ops =
            hier_ops_manager->getOperationsDouble(F_var, patch_hierarchy);
        const double volume = hier_math_ops.getVolumeOfPhysicalDomain();
        const double F_mean = (1.0/volume)*hier_cc_data_ops->integral(F_idx, h_idx);
        hier_cc_data_ops->addScalar(F_idx, F_idx, -F_mean);

        /*
         * Write out data files for plotting.
         */
        if (viz_dump_data)
        {
            if (uses_visit) visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);
        }

        /*
         * Setup the linear solver and preconditioner.
         */
        solv::PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setDConstant(-1.0);
        poisson_spec.setCZero();

        static const bool homogeneous_bc = false;
        solv::LocationIndexRobinBcCoefs<NDIM>* physical_bc_coef =
            new solv::LocationIndexRobinBcCoefs<NDIM>(
                "bc_coef", input_db->getDatabase("LocationIndexRobinBcCoefs"));

        tbox::Pointer<LinearOperator> laplace_op;
        const bool use_fortran_kernels = input_db->getBoolWithDefault("use_fortran_kernels",true);
        if (use_fortran_kernels)
        {
            laplace_op = new CCLaplaceOperator(
                "laplace_op", poisson_spec, physical_bc_coef, homogeneous_bc);
        }
        else
        {
            laplace_op = new CCLaplaceOperator2(
                "laplace_op", poisson_spec, physical_bc_coef, homogeneous_bc);
        }

        tbox::Pointer<PETScKrylovLinearSolver> petsc_linear_solver =
            new PETScKrylovLinearSolver("petsc_linear_solver");
        petsc_linear_solver->setOperator(laplace_op);

        tbox::Pointer<CCPoissonFACOperator> cc_poisson_fac_op =
            new CCPoissonFACOperator(
                "poisson_fac_op", input_db->getDatabase("FACOp"));
        cc_poisson_fac_op->setPoissonSpecifications(poisson_spec);
        cc_poisson_fac_op->setPhysicalBcCoef(physical_bc_coef);

        tbox::Pointer<tbox::Database> fac_db = input_db->getDatabase("FACPreconditioner");
        tbox::Pointer<FACPreconditioner> fac_pc =
            new FACPreconditioner("poisson_fac_pc", *cc_poisson_fac_op, fac_db);
        petsc_linear_solver->setPreconditioner(fac_pc);

        /*
         * Setup the nonlinear solver.
         */
        tbox::Pointer<PETScNewtonKrylovSolver> petsc_nonlinear_solver =
            new PETScNewtonKrylovSolver("petsc_nonlinear_solver");
        petsc_nonlinear_solver->setOperator(laplace_op);

        /*
         * After creating all objects and initializing their state, we print the
         * input database contents to the log file.
         */
        tbox::plog << "\nCheck input data and variables before computation:" << endl;
        tbox::plog << "Input database..." << endl;
        input_db->printClassData(tbox::plog);

        /*
         * Solve the system.
         */
        tbox::TimerManager* timer_manager = tbox::TimerManager::getManager();
        tbox::Pointer<tbox::Timer> t_initialize_solver_state = timer_manager->getTimer("IBTK::main::initializeSolverState()");
        tbox::Pointer<tbox::Timer> t_solve_system = timer_manager->getTimer("IBTK::main::solveSystem()");

        static const int reps = 5;
        static const bool use_linear_solver = true;
        if (use_linear_solver)
        {
            t_initialize_solver_state->start();
            petsc_linear_solver->initializeSolverState(u_vec, f_vec);
            t_initialize_solver_state->stop();

            for (int k = 0; k < reps; ++k)
            {
                u_vec.setToScalar(0.0,false);

                t_solve_system->start();
                petsc_linear_solver->solveSystem(u_vec, f_vec);
                t_solve_system->stop();
            }
        }
        else
        {
            t_initialize_solver_state->start();
            petsc_nonlinear_solver->initializeSolverState(u_vec, f_vec);
            petsc_nonlinear_solver->getLinearSolver()->setPreconditioner(fac_pc);
            laplace_op->modifyRhsForInhomogeneousBc(f_vec);
            t_initialize_solver_state->stop();

            t_solve_system->start();
            for (int k = 0; k < reps; ++k)
            {
                u_vec.setToScalar(0.0,false);
                petsc_nonlinear_solver->solveSystem(u_vec, f_vec);
            }
            t_solve_system->stop();
        }

        /*
         * Compute error, error norms.
         */
        v_vec.subtract(
            tbox::Pointer<solv::SAMRAIVectorReal<NDIM,double> >(&v_vec,false),
            tbox::Pointer<solv::SAMRAIVectorReal<NDIM,double> >(&u_vec,false));

        tbox::pout << "|e|_oo = " << NormOps::maxNorm(&v_vec) << endl;
        tbox::pout << "|e|_2  = " << NormOps::L2Norm(&v_vec) << endl;
        tbox::pout << "|e|_1  = " << NormOps::L1Norm(&v_vec) << endl;

        /*
         * Write out data files for plotting.
         */
        if (viz_dump_data)
        {
            if (uses_visit) visit_data_writer->writePlotData(patch_hierarchy, 1, 1.0);
        }

        /*
         * Write timer data to tbox::plog.
         */
        if (write_timer_data) tbox::TimerManager::getManager()->print(tbox::plog);

        /*
         * Deallocate solver data.
         */
        petsc_nonlinear_solver.setNull();
        petsc_linear_solver.setNull();
        fac_pc.setNull();
        cc_poisson_fac_op.setNull();
        delete physical_bc_coef;
        laplace_op.setNull();
    }// ensure all smart Pointers are properly deleted

    tbox::SAMRAIManager::shutdown();
    PetscFinalize();

    return 0;
}// main
