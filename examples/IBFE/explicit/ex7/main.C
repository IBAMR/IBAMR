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
#include <IBAMR_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petsc.h>

// Headers for basic libMesh objects
#include <../base/variable.h>
#include <boundary_info.h>
#include <cell_hex8.h>
#include <exodusII_io.h>
#include <explicit_system.h>
#include <fe_base.h>
#include <mesh.h>
#include <mesh_generation.h>
#include <string_to_enum.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBFEHierarchyIntegrator.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>
#include <BergerRigoutsos.h>
#include <StandardTagAndInitialize.h>

// C++ namespace delcarations
#include <ibamr/namespaces.h>
using namespace IBAMR;
using namespace IBTK;
using namespace libMesh;
using namespace std;

// Elasticity model data.
namespace ModelData
{
// Stress tensor function.
unsigned int fiber_axes_sys_num;
void
PK1_stress_function(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
    const Point& X,
    const Point& s,
    Elem* const elem,
    const int& e,
    NumericVector<double>& X_vec,
    const std::vector<NumericVector<double>*>& system_data,
    const double& time,
    void* ctx)
{
    // Compute the matrix-type contribution to the stress tensor.
    static const double c = 7.64 * 10000.0;
    const TensorValue<double> CC = FF.transpose() * FF;
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);
    const double I3 = CC.det();
    const double p0 = -c;
    const double beta = 1.0e5*c;
    PP = c*FF + (p0 + beta*log(I3))*FF_inv_trans;

    // Determine the local fiber axes in the reference configuration and compute
    // the fiber contribution to the stress tensor.
    NumericVector<double>* E_ghost = system_data[0];

    static const double k1 = 996.6 * 10000.0;
    static const double k2 = 524.6;

    VectorValue<double> a0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        const int dof_index = elem->dof_number(fiber_axes_sys_num,d+0*NDIM,0);
        a0(d) = (*E_ghost)(dof_index);
    }
    const TensorValue<double> A0 = outer_product(a0,a0);
    const double I4 = CC.contract(A0);
    if (I4 >= 1.0)
    {
        PP += 2.0*k1*(I4-1)*exp(k2*(I4-1)*(I4-1))*FF*A0;
    }

    VectorValue<double> b0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        const int dof_index = elem->dof_number(fiber_axes_sys_num,d+1*NDIM,0);
        b0(d) = (*E_ghost)(dof_index);
    }
    const TensorValue<double> B0 = outer_product(b0,b0);
    const double I6 = CC.contract(B0);
    if (I6 >= 1.0)
    {
        PP += 2.0*k1*(I6-1)*exp(k2*(I6-1)*(I6-1))*FF*B0;
    }
    return;
}// PK1_stress_function

// Applied surface force function.
double W = 0.0;
double T = 0.0;
double L = 0.0;
double F_load = 0.0;
double kappa = 1.0e12;
int X_sys_num = -1;
void
surface_force_function(
    VectorValue<double>& F,
    const TensorValue<double>& dX_ds,
    const Point& X,
    const Point& s,
    Elem* const elem,
    const int& e,
    const unsigned short int side,
    NumericVector<double>& X_vec,
    const std::vector<NumericVector<double>*>& system_data,
    const double& time,
    void* ctx)
{
    F.zero();

    // Only apply force at the top/bottom boundaries.
    if (!(side == 0 || side == 5)) return;

    // Apply a normal traction in the z direction.
    F(2) = 0.5*F_load*1.0e5/(W*T);
    if (side == 0) F(2) = -F(2);

    // Penalize displacements in the x/y directions.
    F(0) = kappa*(s(0)-X(0));
    F(1) = kappa*(s(1)-X(1));

    // Penalize deviations from planar.
    AutoPtr<DofObject> side_dofs = elem->side(side);
    Elem* face = dynamic_cast<Elem*>(side_dofs.get());
    double z_mean = 0.0;
    const int n_nodes = face->n_nodes();
    for (int n = 0; n < n_nodes; ++n)
    {
        const int z_dof_index = face->get_node(n)->dof_number(X_sys_num, NDIM-1, 0);
        z_mean += X_vec(z_dof_index);
    }
    z_mean /= double(n_nodes);
    F(2) += kappa*(z_mean - X(2));
    return;
}// surface_force_function

double
compute_displacement(
    Mesh& mesh,
    EquationSystems& equation_systems)
{
    AutoPtr<QBase> qrule = QBase::build(QGAUSS, NDIM, CONSTANT);
    System& X_system = equation_systems.get_system<System>(IBFEHierarchyIntegrator::COORDS_SYSTEM_NAME);
    DofMap& X_dof_map = X_system.get_dof_map();
    blitz::Array<vector<unsigned int>,1> X_dof_indices(NDIM);
    AutoPtr<FEBase> X_fe(FEBase::build(NDIM, X_dof_map.variable_type(0)));
    X_fe->attach_quadrature_rule(qrule.get());
    X_system.update();
    NumericVector<double>& X_vec = *X_system.current_local_solution;

    MeshBase::const_element_iterator       el_it  = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    double X_max = -0.5*numeric_limits<double>::max();
    double X_min =  0.5*numeric_limits<double>::max();
    for ( ; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;

        X_fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(elem, X_dof_indices(d), d);
        }

        for (unsigned k = 0; k < X_dof_indices(2).size(); ++k)
        {
            X_max = max(X_max,X_vec(X_dof_indices(2)[k]));
            X_min = min(X_min,X_vec(X_dof_indices(2)[k]));
        }
    }
    X_max = SAMRAI_MPI::maxReduction(X_max);
    X_min = SAMRAI_MPI::minReduction(X_min);
    return (X_max - X_min) - L;
}
}
using namespace ModelData;

namespace Postprocessing
{
void
dump_hier_data(
    const int iteration_num,
    const double loop_time,
    const double dt,
    const string hier_dump_dirname,
    const int hier_dump_interval,
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
    const ComponentSelector& hier_data,
    Mesh& mesh,
    EquationSystems& equation_systems)
{
    string file_name;

    // Write Cartesian data.
    file_name = hier_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;

    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putDouble("dt", dt);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->putInteger("hier_dump_interval", hier_dump_interval);
    hier_db->close();

    // Write Lagrangian data.
    file_name = hier_dump_dirname + "/" + "fe_mesh.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    file_name += ".xda";
    mesh.write(file_name);

    file_name = hier_dump_dirname + "/" + "fe_equation_systems.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    equation_systems.write(file_name, (EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA));
    return;
}// dump_hier_data
}
using namespace Postprocessing;

int
main(
    int argc,
    char* argv[])
{
    if (argc == 1)
    {
        pout << "USAGE:  " << argv[0] << " <input filename> <restart dir> <restore number> [options]\n"
             << "OPTIONS: PETSc command line options; use -help for more information\n";
        return -1;
    }

    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAIManager::startup();

    // Process command line options and enable logging.
    const string input_filename = argv[1];
    string restart_read_dirname;
    int restore_num = 0;
    bool is_from_restart = false;
    if (argc >= 4)
    {
        // Check whether this appears to be a restarted run.
        FILE* fstream = (SAMRAI_MPI::getRank() == 0 ? fopen(argv[2], "r") : NULL);
        if (SAMRAI_MPI::bcast(fstream != NULL ? 1 : 0, 0) == 1)
        {
            restart_read_dirname = argv[2];
            restore_num = atoi(argv[3]);
            is_from_restart = true;
        }
        if (fstream != NULL)
        {
            fclose(fstream);
        }
    }

    // Create input database and parse all data in input file.
    Pointer<Database> input_db = new InputDatabase("input_db");
    InputManager::getManager()->parseInputFile(input_filename, input_db);

    // Create a FE mesh corresponding to the single-layer fiber-reinforced strip
    // model of Gasser, Holzapfel and Ogden.
    static const double W =  3.0 * 0.1;  // width     (cm)
    static const double T =  0.5 * 0.1;  // thickness (cm)
    static const double L = 10.0 * 0.1;  // length    (cm)
    static const double dx = input_db->getDouble("DX");
    static const double ds = input_db->getDoubleWithDefault("DS",4.0*dx);
    static const int n1 = ceil(W / ds);
    static const int n2 = ceil(T / ds);
    static const int n3 = ceil(L / ds);
    string elem_type = input_db->getStringWithDefault("elem_type", "HEX8");
    Mesh mesh(NDIM);
    MeshTools::Generation::build_cube(mesh,
                                      n1, n2, n3,
                                      -0.5*W, 0.5*W, -0.5*T, 0.5*T, -0.5*L, 0.5*L,
                                      Utility::string_to_enum<ElemType>(elem_type));
    ModelData::W = W;
    ModelData::T = T;
    ModelData::L = L;

    // Create the FE data manager that manages mappings between the FE mesh and
    // the Cartesian grid.
    const string quad_type = input_db->getStringWithDefault("quad_type", "QGAUSS");
    const string quad_order = input_db->getStringWithDefault("quad_order", "FIFTEENTH");
    AutoPtr<QBase> qrule = QBase::build(Utility::string_to_enum<QuadratureType>(quad_type),NDIM,Utility::string_to_enum<Order>(quad_order));
    const string weighting_fcn = input_db->getStringWithDefault("weighting_fcn", "IB_4");
    const bool use_consistent_mass_matrix = input_db->getBoolWithDefault("use_consistent_mass_matrix", true);
    FEDataManager* const fe_data_manager = FEDataManager::getManager("IBFE Manager", weighting_fcn, weighting_fcn, qrule.get(), use_consistent_mass_matrix);
    const int mesh_level_number = input_db->getInteger("MAX_LEVELS")-1;
    EquationSystems equation_systems(mesh);
    fe_data_manager->setEquationSystems(&equation_systems, mesh_level_number);

    // Build an equation system that stores the material axes.
    ExplicitSystem& fiber_axes_system = equation_systems.add_system<ExplicitSystem>("fiber axes");
    fiber_axes_sys_num = fiber_axes_system.number();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        std::ostringstream os;
        os << "A_" << d;
        fiber_axes_system.add_variable(os.str(), CONSTANT, MONOMIAL);
    }

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        std::ostringstream os;
        os << "B_" << d;
        fiber_axes_system.add_variable(os.str(), CONSTANT, MONOMIAL);
    }

    // Process "Main" section of the input database.
    Pointer<Database> main_db = input_db->getDatabase("Main");

    // Configure logging options.
    const string log_file_name = main_db->getStringWithDefault("log_file_name","IB.log");
    const bool log_all_nodes = main_db->getBoolWithDefault("log_all_nodes",false);
    if (log_all_nodes)
    {
        PIO::logAllNodes(log_file_name);
    }
    else
    {
        PIO::logOnlyNodeZero(log_file_name);
    }

    // Configure visualization options.
    const int viz_dump_interval = main_db->getIntegerWithDefault("viz_dump_interval",0);
    const bool viz_dump_data = viz_dump_interval > 0;
    string viz_dump_dirname;
    bool uses_visit = false;
    bool uses_exodus = false;
    int visit_number_procs_per_file = 1;
    string exodus_filename;
    if (viz_dump_data)
    {
        Array<string> viz_writer;
        if (main_db->keyExists("viz_writer"))
        {
            viz_writer = main_db->getStringArray("viz_writer");
        }
        for (int i = 0; i < viz_writer.getSize(); i++)
        {
            if (viz_writer[i] == "VisIt"   ) uses_visit  = true;
            if (viz_writer[i] == "ExodusII") uses_exodus = true;
        }

        if (main_db->keyExists("viz_dump_dirname"))
        {
            viz_dump_dirname = main_db->getString("viz_dump_dirname");
            if (viz_dump_dirname.empty())
            {
                TBOX_ERROR("viz_dump_interval > 0, but `viz_dump_dirname' is empty" << endl);
            }
        }
        else
        {
            TBOX_ERROR("viz_dump_interval > 0, but key `viz_dump_dirname' not specifed in input file" << endl);
        }

        if (uses_visit)
        {
            visit_number_procs_per_file = main_db->getIntegerWithDefault("visit_number_procs_per_file",visit_number_procs_per_file);
        }

        if (uses_exodus)
        {
            exodus_filename = main_db->getStringWithDefault("exodus_filename","output.ex2");
            ostringstream os;
            os << viz_dump_dirname << "/" << exodus_filename;
            exodus_filename = os.str();
        }
    }

    // Configure restart options.
    const int restart_interval = main_db->getIntegerWithDefault("restart_interval",0);
    const bool write_restart = restart_interval > 0;
    string restart_write_dirname;
    if (write_restart)
    {
        if (main_db->keyExists("restart_write_dirname"))
        {
            restart_write_dirname = main_db->getString("restart_write_dirname");
            if (restart_write_dirname.empty())
            {
                TBOX_ERROR("restart_interval > 0, but `restart_write_dirname' is empty" << endl);
            }
        }
        else
        {
            TBOX_ERROR("restart_interval > 0, but key `restart_write_dirname' not specifed in input file" << endl);
        }
    }

    // Configure timing options.
    const int timer_dump_interval = main_db->getIntegerWithDefault("timer_dump_interval",0);
    const bool write_timer_data = timer_dump_interval > 0;
    if (write_timer_data)
    {
        TimerManager::createManager(input_db->getDatabase("TimerManager"));
    }

    // Configure load balancing options.
    const bool use_nonuniform_load_balancer = main_db->getBoolWithDefault("use_nonuniform_load_balancer", false);

    // Configure postprocessing options.
    const int hier_dump_interval = main_db->getIntegerWithDefault("hier_dump_interval",0);
    const bool write_hier_data = hier_dump_interval > 0;
    string hier_dump_dirname;
    if (write_hier_data)
    {
        if (main_db->keyExists("hier_dump_dirname"))
        {
            hier_dump_dirname = main_db->getString("hier_dump_dirname");
            if (hier_dump_dirname.empty())
            {
                TBOX_ERROR("hier_dump_interval > 0, but `hier_dump_dirname' is empty" << endl);
            }
        }
        else
        {
            TBOX_ERROR("hier_dump_interval > 0, but key `hier_dump_dirname' not specifed in input file" << endl);
        }
    }

    // Process restart data if this is a restarted run.
    if (is_from_restart)
    {
        RestartManager::getManager()->openRestartFile(restart_read_dirname, restore_num, SAMRAI_MPI::getNodes());
    }

    // Create major algorithm and data objects which comprise the application.
    // These objects are configured from the input database and, if this is a
    // restarted run, from the restart database.
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", input_db->getDatabase("CartesianGeometry"));

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>(
        "PatchHierarchy", grid_geometry);

    Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
        "INSStaggeredHierarchyIntegrator", input_db->getDatabase("INSStaggeredHierarchyIntegrator"), patch_hierarchy);

    Pointer<IBFEHierarchyIntegrator> time_integrator = new IBFEHierarchyIntegrator(
        "IBFEHierarchyIntegrator", input_db->getDatabase("IBFEHierarchyIntegrator"),
        patch_hierarchy, navier_stokes_integrator, fe_data_manager);

    Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", time_integrator, input_db->getDatabase("StandardTagAndInitialize"));

    Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();

    Pointer<LoadBalancer<NDIM> > load_balancer = new LoadBalancer<NDIM>(
        "LoadBalancer", input_db->getDatabase("LoadBalancer"));

    Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = new GriddingAlgorithm<NDIM>(
        "GriddingAlgorithm", input_db->getDatabase("GriddingAlgorithm"), error_detector, box_generator, load_balancer);

    // Configure the IBFE solver.
    vector<unsigned int> PK1_stress_function_systems;
    PK1_stress_function_systems.push_back(fiber_axes_sys_num);
    time_integrator->registerPK1StressTensorFunction(&PK1_stress_function, PK1_stress_function_systems);
    time_integrator->registerLagSurfaceForceFunction(&surface_force_function);
    if (use_nonuniform_load_balancer)
    {
        time_integrator->registerLoadBalancer(load_balancer);
    }

    // Create initial conditions specification objects.
    if (input_db->isDatabase("VelocityInitialConditions"))
    {
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction("u_init", input_db->getDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);
    }

    if (input_db->isDatabase("PressureInitialConditions"))
    {
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction("p_init", input_db->getDatabase("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);
    }

    // Create boundary condition specification objects (when necessary).
    const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
    const bool has_physical_boundaries = periodic_shift.min() == 0;
    vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        if (periodic_shift(d) == 0)
        {
            ostringstream bc_coefs_name_stream;
            bc_coefs_name_stream << "u_bc_coefs_" << d;
            const string bc_coefs_name = bc_coefs_name_stream.str();

            ostringstream bc_coefs_db_name_stream;
            bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
            const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

            if (input_db->isDatabase(bc_coefs_db_name))
            {
                u_bc_coefs[d] = new muParserRobinBcCoefs(bc_coefs_name, input_db->getDatabase(bc_coefs_db_name), grid_geometry);
            }
            else
            {
                TBOX_ERROR("dimension " << d << "is nonperiodic.\n" <<
                           "boundary condition data must be provided in database named " << bc_coefs_db_name << endl);
            }
        }
    }
    if (has_physical_boundaries) time_integrator->registerVelocityPhysicalBcCoefs(u_bc_coefs);

    // Setup visualization plot file writers.
    Pointer<VisItDataWriter<NDIM> > visit_data_writer;
    if (uses_visit)
    {
        visit_data_writer = new VisItDataWriter<NDIM>(
            "VisIt Writer",
            viz_dump_dirname, visit_number_procs_per_file);
        time_integrator->registerVisItDataWriter(visit_data_writer);
    }
    AutoPtr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);

    // Initialize hierarchy configuration and data on all patches.
    time_integrator->initializeHierarchyIntegrator(gridding_algorithm);
    double dt_now = time_integrator->initializeHierarchy();

    // Setup the loading data.
    ModelData::kappa = input_db->getDoubleWithDefault("kappa",kappa);
    ModelData::X_sys_num = equation_systems.get_system<System>(fe_data_manager->COORDINATES_SYSTEM_NAME).number();
    double t_step = input_db->getDoubleWithDefault("t_step",6.25e-3);
    double F_step = input_db->getDoubleWithDefault("F_step",0.0625);
    ModelData::F_load = F_step;

    // Setup the fiber axes.
    static const double gamma = 49.98 * M_PI / 180.0;
    const Point a0(cos(gamma),0.0,sin(gamma));
    const Point b0(-cos(-gamma),0.0,-sin(-gamma));
    NumericVector<double>& E = *fiber_axes_system.solution;
    MeshBase::element_iterator elements_begin = mesh.elements_begin();
    MeshBase::element_iterator elements_end = mesh.elements_end();
    for (MeshBase::element_iterator it = elements_begin; it != elements_end; ++it)
    {
        const Elem* const elem = *it;

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const int dof_index = elem->dof_number(fiber_axes_sys_num,d+0*NDIM,0);
            E.set(dof_index, a0(d));
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const int dof_index = elem->dof_number(fiber_axes_sys_num,d+1*NDIM,0);
            E.set(dof_index, b0(d));
        }
    }
    E.close();
    E.localize(*fiber_axes_system.current_local_solution);

    // Close the restart manager.
    RestartManager::getManager()->closeRestartFile();

    // Print the input database contents to the log file.
    plog << "Input database:" << endl;
    input_db->printClassData(plog);

    // Indicate the Eulerian data components to be saved for postprocessing.
    hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
    const int U_data_idx = var_db->mapVariableAndContextToIndex(
        navier_stokes_integrator->getVelocityVar(), navier_stokes_integrator->getCurrentContext());
    const int P_data_idx = var_db->mapVariableAndContextToIndex(
        navier_stokes_integrator->getPressureVar(), navier_stokes_integrator->getCurrentContext());
    const int P_extrap_data_idx = var_db->mapVariableAndContextToIndex(
        navier_stokes_integrator->getExtrapolatedPressureVar(), navier_stokes_integrator->getCurrentContext());
    hier::ComponentSelector hier_data_comps;
    hier_data_comps.setFlag(U_data_idx);
    hier_data_comps.setFlag(P_data_idx);
    hier_data_comps.setFlag(P_extrap_data_idx);

    // Write out initial data.
    double loop_time = time_integrator->getIntegratorTime();
    int iteration_num = time_integrator->getIntegratorStep();
    if (viz_dump_data && iteration_num%viz_dump_interval == 0)
    {
        pout << "\nWriting visualization files...\n\n";
        if (uses_visit)
        {
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
        }
        if (uses_exodus)
        {
            exodus_io->write_timestep(exodus_filename, equation_systems, iteration_num/viz_dump_interval+1, loop_time);
        }
    }
    if (write_hier_data && iteration_num%hier_dump_interval == 0)
    {
        pout << "\nWriting data files for postprocessing...\n\n";
        dump_hier_data(iteration_num, loop_time, dt_now, hier_dump_dirname, hier_dump_interval,
                       patch_hierarchy, hier_data_comps, mesh, equation_systems);
    }

    // Main time step loop.
    const double loop_time_end = time_integrator->getEndTime();
    while (!MathUtilities<double>::equalEps(loop_time,loop_time_end) && time_integrator->stepsRemaining())
    {
        iteration_num = time_integrator->getIntegratorStep() + 1;

        pout <<                                                       endl;
        pout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        pout << "At beginning of timestep # " << iteration_num - 1 << endl;
        pout << "Simulation time is " << loop_time                 << endl;

        const double old_loop_time = loop_time;
        double dt_new = time_integrator->advanceHierarchy(dt_now);
        loop_time += dt_now;
        dt_now = dt_new;

        if (fmod(loop_time,t_step) <= fmod(old_loop_time,t_step))
        {
            ModelData::F_load += F_step;
        }

        pout <<                                                       endl;
        pout << "At end       of timestep # " << iteration_num - 1 << endl;
        pout << "Simulation time is " << loop_time                 << endl;
        pout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        pout <<                                                       endl;

        // At specified intervals, write visualization, restart, and
        // postprocessing files, and print out timer data.
        if (viz_dump_data && iteration_num%viz_dump_interval == 0)
        {
            pout << "\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus)
            {
                exodus_io->write_timestep(exodus_filename, equation_systems, iteration_num/viz_dump_interval+1, loop_time);
            }
        }

        if (write_restart && iteration_num%restart_interval == 0)
        {
            pout << "\nWriting restart files...\n\n";
            RestartManager::getManager()->writeRestartFile(restart_write_dirname, iteration_num);
        }

        if (write_hier_data && iteration_num%hier_dump_interval == 0)
        {
            pout << "\nWriting data files for postprocessing...\n\n";
            dump_hier_data(iteration_num, loop_time, dt_now, hier_dump_dirname, hier_dump_interval,
                           patch_hierarchy, hier_data_comps, mesh, equation_systems);
        }

        if (write_timer_data && iteration_num%timer_dump_interval == 0)
        {
            pout << "\nWriting timer data...\n\n";
            TimerManager::getManager()->print(plog);
        }

        // Compute current displacements.
        const int precision = pout.precision();
        pout.setf(ios_base::scientific);
        pout.precision(8);
        pout << "time = " << loop_time;
        pout << " F_load = ";
        pout.precision(8);
        pout << F_load;
        pout << " displacement = ";
        pout.precision(8);
        pout << compute_displacement(mesh,equation_systems) << "\n";
        pout.unsetf(ios_base::scientific);
        pout.precision(precision);
    }

    // Shutdown SAMRAI.
    SAMRAIManager::shutdown();
    return 0;
}// main
