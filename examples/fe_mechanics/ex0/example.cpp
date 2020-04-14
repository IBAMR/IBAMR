// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
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
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic libMesh objects
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/FEMechanicsExplicitIntegrator.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Elasticity model data.
namespace ModelData
{
// Problem parameters.
static const double mu = 1.0;

// Stress tensor function.
void
PK1_stress_function(TensorValue<double>& PP,
                    const TensorValue<double>& FF,
                    const libMesh::Point& /*x*/,
                    const libMesh::Point& /*X*/,
                    Elem* const /*elem*/,
                    const std::vector<const std::vector<double>*>& /*var_data*/,
                    const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
                    double /*time*/,
                    void* /*ctx*/)
{
    PP = mu * FF;
}
} // namespace ModelData
using namespace ModelData;

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
#ifdef LIBMESH_HAVE_EXODUS_API
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
#else
        const bool uses_exodus = false;
        if (!app_initializer->getExodusIIFilename().empty())
        {
            plog << "WARNING: libMesh was compiled without Exodus support, so no "
                 << "Exodus output will be written in this program.\n";
        }
#endif
        const string exodus_filename = app_initializer->getExodusIIFilename();

        // Create a simple FE mesh.
        Mesh mesh(init.comm(), NDIM);
        const double dx = input_db->getDouble("DX");
        string elem_type = input_db->getString("ELEM_TYPE");
        const int n_x = n_y = nested_meshes ? ceil(1.0 / dx);
        MeshTools::Generation::build_square(
            mesh, n_x, n_y, 0.0, 1.0, 0.0, 1.0, Utility::string_to_enum<ElemType>(elem_type));

        // Configure the FE solver.
        FEMechanicsExplicitIntegrator time_integrator("FEMechanicsExplicitIntegrator",
                           app_initializer->getComponentDatabase("FEMechanicsExplicitIntegrator"),
                           &mesh,
                           /*register_for_restart*/ true,
                           restart_read_dirname,
                           restart_restore_num);

        time_integrator.initializeFEEquationSystems();
        FEDataManager* fe_data_manager = time_integrator.getFEDataManager();
        time_integrator.registerPK1StressFunction(IBFEMethod::PK1StressFcnData(PK1_stress_function));

        // Set up visualization plot file writers.
        std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);

        // Check to see if this is a restarted run to append current exodus files
        if (uses_exodus)
        {
            const bool from_restart = RestartManager::getManager()->isFromRestart();
            exodus_io->append(from_restart);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator.initializeFEData();
        EquationSystems* equation_systems = time_integrator.getFEData()->getEquationSystems();

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = 0; //time_integrator.getIntegratorStep();
        double loop_time = 0.0; //time_integrator.getIntegratorTime();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_exodus)
            {
                if (ib_post_processor) ib_post_processor->postProcessData(loop_time);
                exodus_io->write_timestep(
                    exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        // Main time step loop.
        double loop_time_end = input_db->getDouble("END_TIME");
        double dt = input_db->getDouble("DT");
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end))
        {
            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator.getMaximumTimeStepSize();
            time_integrator.advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator.stepsRemaining();
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator.setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
                if (uses_exodus)
                {
                    if (ib_post_processor) ib_post_processor->postProcessData(loop_time);
                    exodus_io->write_timestep(
                        exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                }
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                time_integrator.writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                pout << "\nWriting state data...\n\n";
                output_data(patch_hierarchy,
                            navier_stokes_integrator,
                            mesh,
                            equation_systems,
                            iteration_num,
                            loop_time,
                            postproc_data_dump_dirname);
            }

            // Compute the volume of the structure.
            double J_integral = 0.0;
            System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
            NumericVector<double>* X_vec = X_system.solution.get();
            NumericVector<double>* X_ghost_vec = X_system.current_local_solution.get();
            copy_and_synch(*X_vec, *X_ghost_vec);
            DofMap& X_dof_map = X_system.get_dof_map();
            std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
            std::unique_ptr<FEBase> fe(FEBase::build(NDIM, X_dof_map.variable_type(0)));
            std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, NDIM, FIFTH);
            fe->attach_quadrature_rule(qrule.get());
            const std::vector<double>& JxW = fe->get_JxW();
            const std::vector<std::vector<VectorValue<double> > >& dphi = fe->get_dphi();
            TensorValue<double> FF;
            boost::multi_array<double, 2> X_node;
            const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
            const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
            for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
            {
                Elem* const elem = *el_it;
                fe->reinit(elem);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X_dof_map.dof_indices(elem, X_dof_indices[d], d);
                }
                const int n_qp = qrule->n_points();
                get_values_for_interpolation(X_node, *X_ghost_vec, X_dof_indices);
                for (int qp = 0; qp < n_qp; ++qp)
                {
                    jacobian(FF, qp, X_node, dphi);
                    J_integral += abs(FF.det()) * JxW[qp];
                }
            }
            J_integral = SAMRAI_MPI::sumReduction(J_integral);
            if (SAMRAI_MPI::getRank() == 0)
            {
                volume_stream.precision(12);
                volume_stream.setf(ios::fixed, ios::floatfield);
                volume_stream << loop_time << " " << J_integral << endl;
            }

            ++iteration_num;
            loop_time += dt;
        }

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
} // main
