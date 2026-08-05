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
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                 LDataManager* l_data_manager,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

static const int nx = 13;
static const int ny = 3;
static const int nz = 3;

static int M_temp;
static double Horizon_size_temp;
static double Poisson_ratio_temp;
static double horizon;
static double P;
static double DX;
static double DX0;
static double k1;
static double k2;
static double K_bulk;
static double current_time;
static double load_time;

// 3d beam
static const double x_begin = -3.0;                                         // left end of the beam
static const double x_end = 3.0;                                            // right end of the beam

static double indi = 1.0;

double
my_inf_fcn(double R0, double /*delta*/)
{
    static const double A = 2.0;
    static const double C = 24.0 / (2.0 * M_PI * A * A * A);

    double W;

    double r = 2.0 *  R0 / horizon;
    if (r < 1.0)
    {
     	W = C * (2.0 / 3.0 - r * r + 0.5 * r * r * r);
    }
    else if (r <= 2.0)
    {
     	W = C * std::pow((2.0 - r), 3) / 6.0;
    }
    else
    {
     	W = 0.0;
    }

    return W;

} // my_inf_fcn

double
my_vol_frac_fcn(double R0, double /*horizon*/, double /*delta*/)
{
    double delta = DX;
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

} // my_vol_frac_fcn

void
my_PK1_fcn(Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor>& PK1,
           const Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> >& FF,
           const Eigen::Map<const IBTK::Vector>& /*X0*/,
           int /*lag_idx*/)
{
    using mat_type = Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor>;

    // stabilized neo-hookean model
    mat_type CC = FF.transpose()*FF;
    mat_type FF_trans = FF.transpose();
    mat_type FF_inv_trans = FF_trans.inverse();
    const double I1 = CC.trace();
    const double J = std::abs(FF.determinant());
    // PK1 = G * pow(J,-2.0/3.0) * (FF - I1 * FF_inv_trans / 3.0) + K_bulk * log(J) * FF_inv_trans;
    
    // modified Mooney-Rivlin
    const double I2 = 0.5 * (I1 * I1 - (CC * CC).trace());
    PK1 = 2* k1 * pow(J,-2.0/3.0) * (FF - I1 * FF_inv_trans / 3.0) + 2 * k2 * pow(J,-4.0/3.0) * (I1 * FF - FF * CC - 2.0 * I2 * FF_inv_trans / 3.0) + K_bulk * log(J) * FF_inv_trans;

    return;
} // my_PK1_fcn

Eigen::Vector4d
my_force_damage_fcn(const double /*horizon*/,
                    const double /*delta*/,
                    const double W,
                    const double vol_frac,
                    double* parameters,
                    const Eigen::Map<const IBTK::Vector>& X0_mastr,
                    const Eigen::Map<const IBTK::Vector>& X0_slave,
                    const Eigen::Map<const IBTK::Vector>& X_mastr,
                    const Eigen::Map<const IBTK::Vector>& X_slave,
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
    const double& R0 = parameters[1];
    const double& vol_mastr = parameters[2];
    const double& vol_slave = parameters[3];
    double& fail = parameters[4];
    const double& critical_stretch = parameters[5];

    // Estimate failure.
    const double R = (X_slave - X_mastr).norm();
    const double stretch = (R - R0) / R0;
    if (!MathUtilities<double>::equalEps(fail, 0.0) && stretch > critical_stretch)
    {
        fail = 0.0;
    }

    // PK1 stress tensor
    using vec_type = IBTK::Vector;
    using mat_type = Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor>;
    mat_type PK1_mastr, PK1_slave;
    my_PK1_fcn(PK1_mastr, FF_mastr, X0_mastr, lag_mastr_node_idx);
    my_PK1_fcn(PK1_slave, FF_slave, X0_slave, lag_slave_node_idx);

    // Compute PD force.
    // const double penalty_fac = 0.0*18.0*K_bulk/(M_PI * pow(horizon,4.0));
    // vec_type pen_trac = penalty_fac * W * ((X_slave - X0_slave) - (X_mastr - X0_mastr));
    // vec_type pen_trac = penalty_fac * stretch * (X0_slave - X0_mastr) / R0;
    vec_type trac = W * (PK1_mastr * B_mastr + PK1_slave * B_slave) * (X0_slave - X0_mastr);

    F_mastr += fail * vol_frac * vol_slave * (trac);
    F_slave += -fail * vol_frac * vol_mastr * (trac);

    // Compute damage.
    Eigen::Vector4d D;
    D(0) = vol_slave * vol_frac * fail;
    D(1) = vol_slave * vol_frac;
    D(2) = vol_mastr * vol_frac * fail;
    D(3) = vol_mastr * vol_frac;

    return D;
    
} // my_force_damage_fcn

void
my_surface_force_func(const Eigen::Map<const IBTK::Vector>& X,
                          const Eigen::Map<const IBTK::Vector>& X_target,
                          const Eigen::Map<const IBTK::Vector>& U,
                          int lag_idx,
                          Eigen::Map<IBTK::Vector>& F)
{
    //X_target is the material variable, X is the spatial variable
    static double kappa = 1.0e8;

    // 3d beam
    if (X_target(0) == x_begin)
    {
        F += kappa * (X_target - X);
    }
    else if (X_target(0) == x_end)
    {   
        double angle = current_time / load_time * M_PI * 2.5;

        if (current_time >= load_time)
        {
            angle = M_PI * 2.5;
        }

        double y = cos(angle) * X_target(1) + sin(angle) * X_target(2);
        double z = - sin(angle) * X_target(1) + cos(angle) * X_target(2);

        F(1) += kappa * (y - X(1));
        F(2) += kappa * (z - X(2));
    }

    return;
} // my_surface_force_func

void
my_target_point_force_fcn(const Eigen::Map<const IBTK::Vector>& X,
                          const Eigen::Map<const IBTK::Vector>& X_target,
                          const Eigen::Map<const IBTK::Vector>& U,
                          double K,
                          double E,
                          int lag_idx,
                          Eigen::Map<IBTK::Vector>& F)
{
    // F += K * (X_target - X) - E * U;
    // F += - E * U;

    my_surface_force_func(X, X_target, U, lag_idx, F);

    return;
} // my_target_point_force_fcn

class MyIBPDMethod : public IBPDMethod
{
private:
    std::fstream d_ss_stream;

public:
    MyIBPDMethod(std::string object_name, Pointer<Database> input_db, bool register_for_restart = true)
        : IBPDMethod(std::move(object_name), input_db, register_for_restart)
    {
        bool from_restart = RestartManager::getManager()->isFromRestart();
        if (!from_restart)
        {
            tbox::Utilities::recursiveMkdir("./data");

            std::ostringstream sstream;
            sstream << "./data/check_ss.txt";
            d_ss_stream.open(sstream.str().c_str(), std::fstream::out);
        }

        return;

    } // MyIBPDMethod

    ~MyIBPDMethod() = default;

    // void preprocessIntegrateData(double current_time, double new_time, int num_cycles) override
    // {
    //     IBPDMethod::preprocessIntegrateData(current_time, new_time, num_cycles);
    //     return;
    // } // preprocessIntegrateData

    void postprocessIntegrateData(double current_time, double new_time, int num_cycles) override
    {
        const int struct_ln = d_hierarchy->getFinestLevelNumber();
        const int step_no = d_ib_solver->getIntegratorStep() + 1;

        if (step_no % 1000 == 0)
        {
            Pointer<LData> D_LData = d_l_data_manager->getLData("damage", struct_ln);
            Vec D_petsc_vec_parallel = D_LData->getVec();
            Vec D_lag_vec_parallel = nullptr;
            Vec D_lag_vec_seq = nullptr;
            VecDuplicate(D_petsc_vec_parallel, &D_lag_vec_parallel);
            d_l_data_manager->scatterPETScToLagrangian(D_petsc_vec_parallel, D_lag_vec_parallel, struct_ln);
            d_l_data_manager->scatterToZero(D_lag_vec_parallel, D_lag_vec_seq);

            Pointer<LData> X0_LData = d_l_data_manager->getLData("X0", struct_ln);
            Vec X0_petsc_vec_parallel = X0_LData->getVec();
            Vec X0_lag_vec_parallel = nullptr;
            Vec X0_lag_vec_seq = nullptr;
            VecDuplicate(X0_petsc_vec_parallel, &X0_lag_vec_parallel);
            d_l_data_manager->scatterPETScToLagrangian(X0_petsc_vec_parallel, X0_lag_vec_parallel, struct_ln);
            d_l_data_manager->scatterToZero(X0_lag_vec_parallel, X0_lag_vec_seq);

            Pointer<LData> X_LData = d_X_new_data[struct_ln];
            Vec X_petsc_vec_parallel = X_LData->getVec();
            Vec X_lag_vec_parallel = nullptr;
            Vec X_lag_vec_seq = nullptr;
            VecDuplicate(X_petsc_vec_parallel, &X_lag_vec_parallel);
            d_l_data_manager->scatterPETScToLagrangian(X_petsc_vec_parallel, X_lag_vec_parallel, struct_ln);
            d_l_data_manager->scatterToZero(X_lag_vec_parallel, X_lag_vec_seq);

            if (IBTK_MPI::getRank() == 0)
            {
                const PetscScalar* D;
                VecGetArrayRead(D_lag_vec_seq, &D);
                int counter_D = -1;

                const PetscScalar* X0;
                VecGetArrayRead(X0_lag_vec_seq, &X0);
                int counter_X0 = -1;

                const PetscScalar* X;
                VecGetArrayRead(X_lag_vec_seq, &X);
                int counter_X = -1;

                std::fstream D_stream;
                std::ostringstream D_sstream;
                D_sstream << "./data/D_" << step_no << "_" << new_time;
                D_stream.open(D_sstream.str().c_str(), std::fstream::out);

                int ib_pts;
                VecGetSize(D_lag_vec_seq, &ib_pts);
                for (int i = 0; i < ib_pts; ++i)
                {
                    const double X0_0 = X0[++counter_X0];
                    const double X0_1 = X0[++counter_X0];
                    const double X0_2 = X0[++counter_X0];

                    const double X_0 = X[++counter_X];
                    const double X_1 = X[++counter_X];
                    const double X_2 = X[++counter_X];

                    const double dmg = D[++counter_D];
                    D_stream << X0_0 << "\t" << X0_1 << "\t" << X0_2 << "\t" << X_0 - X0_0 << "\t" << X_1 - X0_1 << "\t" << X_2 - X0_2 << "\t" << dmg
                             << std::endl;
                }

                VecRestoreArrayRead(D_lag_vec_seq, &D);
                VecRestoreArrayRead(X0_lag_vec_seq, &X0);
                VecRestoreArrayRead(X_lag_vec_seq, &X);
            }
            VecDestroy(&D_lag_vec_parallel);
            VecDestroy(&D_lag_vec_seq);
            VecDestroy(&X0_lag_vec_parallel);
            VecDestroy(&X0_lag_vec_seq);
            VecDestroy(&X_lag_vec_parallel);
            VecDestroy(&X_lag_vec_seq);
        }

        IBPDMethod::postprocessIntegrateData(current_time, new_time, num_cycles);

        return;
    } // postprocessIntegrateData

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
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        M_temp = input_db->getInteger("M");
        Horizon_size_temp = input_db->getDouble("HORIZON_SIZE");
        Poisson_ratio_temp = input_db->getDouble("POISSON_RATIO");
        k1 = input_db->getDouble("K1");
        k2 = input_db->getDouble("K2");
        K_bulk = input_db->getDouble("BULK_MOD");
        DX0 = input_db->getDouble("DX");
        load_time = input_db->getDouble("LOAD_TIME");
        DX = (x_end - x_begin) / double(6.0 * M_temp);
        horizon = Horizon_size_temp * DX;
        P = Poisson_ratio_temp;


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
        Pointer<MyIBPDMethod> ib_method_ops =
            new MyIBPDMethod("MyIBPDMethod", app_initializer->getComponentDatabase("IBPDMethod"));
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
            0, &my_PK1_fcn, &my_force_damage_fcn, &my_inf_fcn, &my_vol_frac_fcn);
        ib_force_fcn->registerTargetPointForceFunction(&my_target_point_force_fcn);
        ib_method_ops->registerIBPDForceGen(ib_force_fcn);

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
            current_time = loop_time;

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
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
