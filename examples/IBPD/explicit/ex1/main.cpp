// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBPDForceGen.h>
#include <ibamr/IBPDMethod.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h> 
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                 LDataManager* l_data_manager,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

static const int N = 128;
static const double MFAC = 1.0;
static const double ds = (1.0 / N) * MFAC;
static const double w = 0.2;
static const double xc = 0.4;
static const double yc = 0.4;
static const double zc = 0.0;

double
sheet_inf_fcn(double R0, double /*delta*/)
{
    static const double A = 2;
    static const double C = 24.0 / (2.0 * M_PI * A * A * A);

    double W;

    double r = R0 / ds;
    if (r < 1)
    {
        W = C * (2.0 / 3.0 - r * r + 0.5 * r * r * r);
    }
    else if (r <= 2)
    {
        W = C * std::pow((2.0 - r), 3) / 6.0;
    }
    else
    {
        W = 0.0;
    }

    return W;

} // sheet_inf_fcn

double
sheet_vol_frac_fcn(double R0, double horizon, double delta)
{
    double vol_frac;
    if (R0 <= (horizon - delta))
    {
        vol_frac = 1.0;
    }
    else if (R0 <= (horizon + delta))
    {
        vol_frac = (horizon + delta - R0) / (2.0 * delta);
    }
    else
    {
        vol_frac = 0.0;
    }

    return vol_frac;

} // sheet_vol_frac_fcn

// void
// sheet_residual_PK1_fcn(Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor>& PK1,
//                        const Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> >& FF,
//                        const Eigen::Map<const IBTK::Vector>& X0,
//                        int /*lag_idx*/)
// {
//     const double yy = X0(1) - yc;
//     const double xx = X0(0) - xc;

//     const double dis = sqrt(1.0 - 4.0 * xx + 2.0 * yy + yy * yy);
//     const double s0 = 0.5 * (1.0 + yy - dis);
//     const double s1 = 0.5 * (-1.0 + yy + dis);

//     if ((s0 < 0.0 && abs(s0) > 1e-8) || (s1 < 0.0 && abs(s1) > 1e-8))
//     {
//         std::cout << "s0 = " << s0 << "\ts1 = " << s1 << std::endl;
//         TBOX_ERROR("Negative s' encountered \n");
//     }

//     Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> FF0;
//     FF0 << 1.0 + s1, s0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0;
//     const double J0 = FF0.determinant();

//     PK1 = (1.0 / w) * (1 / J0) * FF * FF0 * FF0.transpose();

//     const double x_ana = s0 + s0 * s1 + xc;
//     const double y_ana = s0 + s1 + yc;
//     if (!MathUtilities<double>::equalEps(x_ana, X0(0)) || !MathUtilities<double>::equalEps(y_ana, X0(1)))
//     {
//         TBOX_ERROR("Mismatch in values.\n"
//                    << " x_ana = " << x_ana << " x_num = " << X0(0) << " y_ana = " << y_ana << " y_num = " << X0(1)
//                    << "\n");
//     }

//     return;

// } // sheet_residual_PK1_fcn

void
coordTransform(Eigen::Map<IBTK::Point>& X,
               Eigen::Map<const IBTK::Point>& X0,
               int /*lag_idx*/,
               int /*level_number*/,
               void* /*ctx*/)
{
    const int Nndim = X.size();
    X(0) = X0(0) + X0(1)*xc;
    X(1) = X0(0)*yc + X0(1);
    if (Nndim == 3)
    {
        X(2) = X0(2) + zc;
    }
    return;

} // coordTransform

void
sheet_PK1_fcn(Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor>& PK1,
              const Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> >& FF,
              const Eigen::Map<const IBTK::Vector>& /*X0*/,
              int /*lag_idx*/)
{
    PK1 = (1.0 / w) * FF;
    return;

} // sheet_PK1_fcn

Eigen::Vector4d
sheet_force_damage_fcn(const double /*horizon*/,
                       const double /*delta*/,
                       const double W,
                       const double vol_frac,
                       double* parameters,
                       const Eigen::Map<const IBTK::Vector>& X0_mastr,
                       const Eigen::Map<const IBTK::Vector>& X0_slave,
                       const Eigen::Map<const IBTK::Vector>& /*X_mastr*/,
                       const Eigen::Map<const IBTK::Vector>& /*X_slave*/,
                       const Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> >& FF_mastr,
                       const Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> >& FF_slave,
                       const Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> >& B_mastr,
                       const Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> >& B_slave,
                       Eigen::Map<IBTK::Vector>& F_mastr,
                       Eigen::Map<IBTK::Vector>& F_slave,
                       const int lag_mastr_node_idx,
                       const int lag_slave_node_idx)
{
    // Bond parameters
    // 0 --> Kappa, 1 --> R0, 2--> user defined
    const double& vol_mastr = parameters[2];
    const double& vol_slave = parameters[3];
    double& fail = parameters[4];

    // PK1 stress tensor
    using vec_type = IBTK::Vector;
    using mat_type = Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor>;
    mat_type PK1_mastr, PK1_slave;
    sheet_PK1_fcn(PK1_mastr, FF_mastr, X0_mastr, lag_mastr_node_idx);
    sheet_PK1_fcn(PK1_slave, FF_slave, X0_slave, lag_slave_node_idx);

    // Compute PD force.
    vec_type trac = W * (PK1_mastr * B_mastr + PK1_slave * B_slave) * (X0_slave - X0_mastr);
    trac(2) = 0.0;
    F_mastr += fail * (vol_frac * vol_slave) * trac * (vol_frac * vol_slave);
    F_slave += -fail * (vol_frac * vol_mastr) * trac * (vol_frac * vol_mastr);

    // Compute damage.
    Eigen::Vector4d D;
    D(0) = vol_slave * vol_frac * fail;
    D(1) = vol_slave * vol_frac;
    D(2) = vol_mastr * vol_frac * fail;
    D(3) = vol_mastr * vol_frac;

    return D;

} // sheet_force_damage_fcn

class SheetIBPDMethod : public IBPDMethod
{
public:
    SheetIBPDMethod(std::string object_name, Pointer<Database> input_db, bool register_for_restart = true)
        : IBPDMethod(std::move(object_name), input_db, register_for_restart)
    {
        return;
    } // SheetIBPDMethod

    ~SheetIBPDMethod() = default;

    void preprocessIntegrateData(double current_time, double new_time, int num_cycles) override
    {
        IBPDMethod::preprocessIntegrateData(current_time, new_time, num_cycles);

        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

            Pointer<LData> U_current_data = d_U_current_data[ln];
            boost::multi_array_ref<double, 2>& U_current_data_array = *U_current_data->getLocalFormVecArray();

            const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
            const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

            for (const auto& node : local_nodes)
            {
                const int local_idx = node->getLocalPETScIndex();

                double* U_current = &U_current_data_array[local_idx][0];
                U_current[2] = 0.0;
            }
        }

        return;
    } // preprocessIntegrateData

    void midpointStep(const double current_time, const double new_time) override
    {
        IBPDMethod::midpointStep(current_time, new_time);

        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

            Pointer<LData> X_0_data = d_l_data_manager->getLData("X0_unshifted", ln);
            Pointer<LData> X_new_data = d_X_new_data[ln];
            Pointer<LData> U_new_data = d_U_new_data[ln];

            boost::multi_array_ref<double, 2>& X_0_data_array = *X_0_data->getLocalFormVecArray();
            boost::multi_array_ref<double, 2>& X_new_data_array = *X_new_data->getLocalFormVecArray();
            boost::multi_array_ref<double, 2>& U_new_data_array = *U_new_data->getLocalFormVecArray();

            const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
            const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
            for (const auto& node : local_nodes)
            {
                const int local_idx = node->getLocalPETScIndex();

                double* U_new = &U_new_data_array[local_idx][0];
                const double* X_0 = &X_0_data_array[local_idx][0];
                double* X_new = &X_new_data_array[local_idx][0];

                U_new[2] = 0.0;
                X_new[2] = X_0[2];
            }

            X_0_data->restoreArrays();
            X_new_data->restoreArrays();
            U_new_data->restoreArrays();
        }

        return;

    } // midpointStep
};

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
main(int argc, char* argv[]) // argc : counts the numer of arguments
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IBPD.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        Pointer<SheetIBPDMethod> ib_method_ops =
            new SheetIBPDMethod("MyIBPDMethod", app_initializer->getComponentDatabase("IBPDMethod"));
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Configure the IB solver.
        Pointer<IBStandardInitializer> ib_initializer = new IBStandardInitializer(
            "IBStandardInitializer", app_initializer->getComponentDatabase("IBStandardInitializer"));
        ib_method_ops->registerLInitStrategy(ib_initializer);
        Pointer<IBPDForceGen> ib_force_fcn = new IBPDForceGen(app_initializer->getComponentDatabase("IBPDForceGen"));
        ib_force_fcn->registerBondForceSpecificationFunction(
            0, &sheet_PK1_fcn, &sheet_force_damage_fcn, &sheet_inf_fcn, &sheet_vol_frac_fcn);
        ib_method_ops->registerIBPDForceGen(ib_force_fcn);
        ib_method_ops->registerInitialCoordinateMappingFunction(&coordTransform, nullptr);

        // Create Eulerian initial condition specification objects.  These
        // objects also are used to specify exact solution values for error
        // analysis.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);

        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = nullptr;
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

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Initialize PD data.
        ib_method_ops->initializePDData();

        // Deallocate initialization objects.
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            silo_data_writer->writePlotData(iteration_num, loop_time);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();

            bool at_regrid = time_integrator->atRegridPoint();
            if (at_regrid) ib_method_ops->setIBPDForceGenNeedsInit();
            time_integrator->advanceHierarchy(dt);
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
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
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
                            ib_method_ops->getLDataManager(),
                            iteration_num,
                            loop_time,
                            postproc_data_dump_dirname);
            }
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    return 0;
} // main

void
output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
            LDataManager* l_data_manager,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    // Write Cartesian data.
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, IBTK_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();

    // Write Lagrangian data.
    const int finest_hier_level = patch_hierarchy->getFinestLevelNumber();
    Pointer<LData> X_data = l_data_manager->getLData("X", finest_hier_level);
    Vec X_petsc_vec = X_data->getVec();
    Vec X_lag_vec;
    VecDuplicate(X_petsc_vec, &X_lag_vec);
    l_data_manager->scatterPETScToLagrangian(X_petsc_vec, X_lag_vec, finest_hier_level);
    file_name = data_dump_dirname + "/" + "X.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
    VecView(X_lag_vec, viewer);
    PetscViewerDestroy(&viewer);
    VecDestroy(&X_lag_vec);
    return;
} // output_data
