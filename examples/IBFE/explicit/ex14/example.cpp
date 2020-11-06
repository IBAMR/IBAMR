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
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>
#include "HierarchyDataOpsManager.h"

// Headers for basic libMesh objects
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/periodic_boundary.h>
#include <libmesh/boundary_info.h>
#include <libmesh/mesh_function.h>
#include <libmesh/point.h>
#include "libmesh/point_locator_base.h"
#include "libmesh/point_locator_tree.h"
#include "libmesh/fe_interface.h"

// Headers for application-specific algorithm/data structure objects
#include <boost/multi_array.hpp>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFECentroidPostProcessor.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Elasticity model data.
namespace ModelData
{
    // Problem parameters.
	static double kappa_tether;
	static double eta_tether;
	static double eta_tether_volume;

    // surface pressure function parameters
    double P_load;
    double t_load;
    double period;

    static libMesh::TensorValue <double> identity( 1.0, 0.0, 0.0,
                                                   0.0, 1.0, 0.0,
                                                   0.0, 0.0, 1.0 );

    static BoundaryInfo* boundary_info;
    static double C;
    void PK1(libMesh::TensorValue<double> &PP, const libMesh::TensorValue<double> &FF, const libMesh::Point& X, const libMesh::Point& /*s*/,
            Elem *const elem, const vector<const vector<double>*> &var_data, const vector<const vector<libMesh::VectorValue<double> >*>& /*grad_var_data*/,
            double time, void* /*ctx*/)
    {
    	// From: https://royalsocietypublishing.org/doi/full/10.1098/rspa.2015.0641
    	// Verification of cardiac mechanics software: benchmark problems and solutions for testing active and passive material behaviour
    	//
    	// W = C / 2 ( e^Q - 1 ),
    	// Q = bf E11^2 + bt * ( E22^2 + E33^2 + E32^2 + E23^2 ) + bfs ( E12^2 + E21^2 + E13^2 + E31^2 )

    	// SS = dW / dE = C / 2 * e^Q dQ/dE
    	//double C   = 2.0;
    	double E1 = 250.0;
    	double E2 = 80.0;
    	double mu1 = E1 / 3.0;
    	double mu2 = E2 / 3.0;
    	auto blockID = elem->subdomain_id();
    	double mu = mu1;
    	if(blockID == 2) mu = mu2;
    	PP = 2.0 * mu * FF;
    	auto FFinv = FF.inverse();
    	auto FFinvT = FFinv.transpose();
    	// Evaluate the deviatoric part
    	PP -= 1.0 / 3.0 * PP.contract(FF) * FFinvT;
    }


    void pressure_load(double& F,
            const libMesh::VectorValue<double>& n,
            const libMesh::VectorValue<double>& N,
            const libMesh::TensorValue<double>& FF,
            const libMesh::Point& x,
            const libMesh::Point& X,
            libMesh::Elem* elem,
            unsigned short int side,
            const std::vector<const std::vector<double>*>& system_var_data,
            const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
            double data_time,
            void* ctx)
    {
       if (elem->neighbor_ptr(side) == libmesh_nullptr)
       {
        const libMesh::boundary_id_type boundary_id = boundary_info->boundary_id (elem, side);
        if ( boundary_id == 1)
        {
        	F = P_load;
        }
        else
        {
        	F = 0.0;
        }
       }
    }



    void
    solid_body_force_function(VectorValue<double>& F,
                              const TensorValue<double>& /*FF*/,
                              const libMesh::Point& /*x*/,
                              const libMesh::Point& /*X*/,
                              Elem* const elem,
                              const vector<const vector<double>*>& var_data,
                              const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                              double time,
                              void* /*ctx*/)
    {
        VectorValue<double> U;
        std::copy(var_data[0]->begin(), var_data[0]->end(), &U(0));
        F = -eta_tether * U;
    }


    static double kappa;
    void
   PK1_dil_stress_function(TensorValue<double>& PP,
           const TensorValue<double>& FF,
           const libMesh::Point& /*X*/,
           const libMesh::Point& /*s*/,
           Elem* const /*elem*/,
           const std::vector<const std::vector<double>*>& /*var_data*/,
           const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
           double /*time*/,
           void* /*ctx*/)
   {
       const double J = FF.det();
       const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);
       double E1 = 250.0;
       double E2 = 80.0;
       double nu = 0.49995;
		auto blockID = elem->subdomain_id();
		double E = Eu1;
		if(blockID == 2) E = E2;
       kappa = E / 3.0 / (1-2*nu);
       PP = kappa * log(J) * FF_inv_trans;
       return;
   } // PK1_dil_stress_function


    void surface_tether_force_function(VectorValue<double>& F,
                                       const libMesh::VectorValue<double>& n,
                                       const libMesh::VectorValue<double>& N,
                                       const libMesh::TensorValue<double>& FF,
                                       const libMesh::Point& x,
                                       const libMesh::Point& X,
                                       libMesh::Elem* elem,
                                       unsigned short int side,
                                       const std::vector<const std::vector<double>*>& system_var_data,
                                       const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
                                       double time,
                                       void* ctx)
    {
        VectorValue<double> U;
        std::copy(system_var_data[0]->begin(), system_var_data[0]->end(), &U(0));

        if (elem->neighbor_ptr(side) == libmesh_nullptr)
        {
         const libMesh::boundary_id_type boundary_id = boundary_info->boundary_id (elem, side);
         if ( boundary_id == 4)
         {
         	F = kappa_tether * (X - x) + eta_tether * U;;
         }
         else
         {
         	F = 0.0;
         }
        }


    }
}



using namespace ModelData;

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                 Mesh& mesh,
                 EquationSystems* equation_systems,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

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

int main(int argc, char** argv)
{
    std::cout << "Starting ... " << std::endl;
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();
    std::cout << "Start up ... " << std::endl;

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        Pointer<Database> ibfe_db = app_initializer->getComponentDatabase("IBFEMethod");

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        string exodus_filename = app_initializer->getExodusIIFilename();

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

        // for computing velocity errors
        std::vector<double> u_err;
        u_err.resize(3);

        // Create a simple FE mesh
        ModelData::eta_tether = input_db->getDouble("eta_tether");
        ModelData::kappa_tether = input_db->getDouble("kappa_tether");
        ModelData::eta_tether_volume = input_db->getDouble("eta_tether_volume");
        Mesh mesh(init.comm(), NDIM);
        const double dx = input_db->getDouble("DX");
        const double MFAC = input_db->getDouble("MFAC");
        const double ds = MFAC * dx;
        string elem_type = input_db->getString("ELEM_TYPE");

        double xmin = 0.0;
        double ymin = 0.0;
        double zmin = 0.0;
        double xmax =  input_db->getDouble("XMAX");
        double ymax =  input_db->getDouble("YMAX");
        double zmax =  input_db->getDouble("ZMAX");
        ModelData::C =  input_db->getDouble("C");
        ModelData::P_load =  input_db->getDouble("P_load");

        const int n_x = round( (xmax-xmin) / ds );
        const int n_y = round( (ymax-ymin) / ds );
        const int n_z = round( (zmax-zmin) / ds );

        MeshTools::Generation::build_cube(
            mesh, n_x, n_y, n_z,
			xmin, xmax,
			ymin, ymax,
			zmin, zmax,
			Utility::string_to_enum<ElemType>(elem_type));

        // get boundary information and surface pressure info
        boundary_info = &mesh.get_boundary_info();


        // Change sideset ID
        // change sideset ID
        const libMesh::MeshBase::const_element_iterator end_el = mesh.elements_end();
        for (libMesh::MeshBase::const_element_iterator el = mesh.elements_begin(); el != end_el; ++el)
        {
            libMesh::Elem* const elem = *el;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor_ptr(side);
                if (at_mesh_bdry)
                {
                    if (boundary_info->has_boundary_id(elem, side, 4) )
                    {
                        boundary_info->add_side(elem, side, FEDataManager::ZERO_DISPLACEMENT_XYZ_BDRY_ID);
                    }
                }
            }
        }

        std::cout << "INSHierarchyIntegrator" << std::endl;
        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }
        std::cout << "IBFEMethod" << std::endl;
        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           &mesh,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                           true,
                           app_initializer->getRestartReadDirectory(),
                           app_initializer->getRestartRestoreNumber());
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

        // get velocity system data for damping in tether forces
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        vector<SystemData> sys_data(1);
        sys_data[0] = SystemData(IBFEMethod::VELOCITY_SYSTEM_NAME, vars);

        // Configure the IBFE solver.
        std::cout << "IBFEMethod initialize" << std::endl;
        ib_method_ops->initializeFEEquationSystems();
        std::cout << "FEDataMandager" << std::endl;
        FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager();
        //std::cout << "Register Coordinate Mapping" << std::endl;
        //ib_method_ops->registerInitialCoordinateMappingFunction(coordinate_mapping_function);
        std::cout << "Register PK1 stress" << std::endl;
        {
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data(ModelData::PK1);
            PK1_dev_stress_data.quad_order = Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "FIFTH"));
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data);
            //damping
            IBFEMethod::LagBodyForceFcnData solid_body_force_data(ModelData::solid_body_force_function, sys_data);
            ib_method_ops->registerLagBodyForceFunction(solid_body_force_data);
            IBFEMethod::LagSurfaceForceFcnData surface_tether_force_data(ModelData::surface_tether_force_function, sys_data);
            ib_method_ops->registerLagSurfaceForceFunction(surface_tether_force_data);
            // pressure bc
            IBFEMethod::LagSurfacePressureFcnData surface_pressure_data(ModelData::pressure_load);
            ib_method_ops->registerLagSurfacePressureFunction(surface_pressure_data,0);
        }



        ModelData::kappa = input_db->getDoubleWithDefault("kappa", 2e8);
        // setup libmesh things for eliminating pressure jumps
        if (input_db->getBoolWithDefault("USE_PRESSURE_FIELD", false))
        {
            std::cout << "Registering pressure field" << std::endl;
            std::string projection_type = input_db->getStringWithDefault("P_PROJECTION", "CONSISTENT_PROJECTION");
            ib_method_ops->registerStaticPressurePart(IBAMR::string_to_enum<PressureProjectionType>(projection_type));
            double ctau = input_db->getDoubleWithDefault("CTAU", 1.0);
            double tau =  ctau * ModelData::kappa / ModelData::C;
            ib_method_ops->set_static_pressure_stab_param(tau);
            ib_method_ops->set_static_pressure_kappa(ModelData::kappa);
        }
        else
        {
			IBFEMethod::PK1StressFcnData PK1_dil_stress_data(PK1_dil_stress_function);
			PK1_dil_stress_data.quad_order =
					Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DIL_QUAD_ORDER", "FIFTH"));
			ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data);
        }

        // Set up post processor to recover computed stresses.
        Pointer<IBFEPostProcessor> ib_post_processor =
                new IBFECentroidPostProcessor("IBFEPostProcessor", fe_data_manager);
        FEDataManager::InterpSpec p_interp_spec("PIECEWISE_LINEAR",
                QGAUSS,
                FIFTH,
                /*use_adaptive_quadrature*/ false,
                /*point_density*/ 2.0,
                /*use_consistent_mass_matrix*/ true,
                /*use_nodal_quadrature*/ false);

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
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
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        UniquePtr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        EquationSystems* equation_systems = fe_data_manager->getEquationSystems();

	// add system to look at pressure field on FE mesh
	libMesh::System& system1 = equation_systems->add_system<System>("p_f system");
	system1.add_variable("p_f", FIRST, LAGRANGE);

        ib_method_ops->initializeFEData();
        if (ib_post_processor) ib_post_processor->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Setup data used to determine the accuracy of the computed solution.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const Pointer<hier::Variable<NDIM> > u_var = navier_stokes_integrator->getVelocityVariable();
        const Pointer<VariableContext> u_ctx = navier_stokes_integrator->getCurrentContext();
        const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);
        double loop_time = time_integrator->getIntegratorTime();
        const Pointer<hier::Variable<NDIM> > p_var = navier_stokes_integrator->getPressureVariable();
        Pointer<VariableContext> p_current_ctx = navier_stokes_integrator->getCurrentContext();
        const int p_current_idx = var_db->mapVariableAndContextToIndex(p_var, p_current_ctx);
        Pointer<CellVariable<NDIM, double> > p_copy_var = new CellVariable<NDIM, double>("p_copy");
        const IntVector<NDIM> ib_ghosts = ib_method_ops->getMinimumGhostCellWidth();
        const int p_copy_idx =
        var_db->registerVariableAndContext(p_copy_var, time_integrator->getScratchContext(), ib_ghosts);
        visit_data_writer->registerPlotQuantity("P copy", "SCALAR", p_copy_idx);
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(p_copy_idx, loop_time);
        }
        HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_cc_data_ops = hier_data_ops_manager->getOperationsDouble(p_var, patch_hierarchy, true);
        hier_cc_data_ops->copyData(p_copy_idx, p_current_idx, true);

        // fill ghost cells for pressure
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> p_transaction_comp(1);
        p_transaction_comp[0] = InterpolationTransactionComponent(p_copy_idx,
                "CONSERVATIVE_LINEAR_REFINE",
                /*use_cf_bdry_interpolation*/ false,
                "CONSERVATIVE_COARSEN",
                "LINEAR");

        Pointer<HierarchyGhostCellInterpolation> p_hier_bdry_fill = new HierarchyGhostCellInterpolation();
        p_hier_bdry_fill->initializeOperatorState(p_transaction_comp, patch_hierarchy);
        p_hier_bdry_fill->fillData(loop_time);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus)
            {
                if (ib_post_processor) ib_post_processor->postProcessData(loop_time);
                exodus_io->write_timestep(
                    exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        // Open streams to save volume of structure.
        //ofstream volume_stream;
        //ofstream stuff_stream;
        ofstream dX_stream;
        ofstream disp_stream;
        ofstream pressure_field_stream;
        ofstream vel_error_stream;
        disp_stream.precision(12);
        DenseVector<double> dX_center_top_serial; // for looking at displacement of the top point
        DenseVector<double> dX_center_top_parallel; // for looking at displacement of the top point
        dX_center_top_serial.resize(NDIM);
        dX_center_top_parallel.resize(NDIM);

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;

        double umax_old = 0.0;
        double umax = 0.0;
        bool stop = false;

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining() && !stop)
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();


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

            HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
            hier_math_ops.setPatchHierarchy(patch_hierarchy);
            hier_math_ops.resetLevels(coarsest_ln, finest_ln);
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                if (!level->checkAllocated(p_copy_idx)) level->allocatePatchData(p_copy_idx);
            }

            hier_cc_data_ops->copyData(p_copy_idx, p_current_idx, true);
            p_hier_bdry_fill->initializeOperatorState(p_transaction_comp, patch_hierarchy);
            p_hier_bdry_fill->fillData(loop_time);

            // get the displacement system and build a mesh function for evaluating the displacement
            System& dX_system = equation_systems->get_system<System>(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
            NumericVector<double>* dX_vec = dX_system.solution.get();
            NumericVector<double>* dX_ghost_vec = dX_system.current_local_solution.get();
            dX_vec->localize(*dX_ghost_vec);
            std::vector<unsigned int> vars(NDIM);
            for(int dd = 0; dd < NDIM; ++dd) vars[dd] = dd;
            MeshFunction mesh_fcn(*equation_systems,
                                  *dX_system.current_local_solution,
                                  dX_system.get_dof_map(),
                                  vars);
            mesh_fcn.init();
            mesh_fcn.set_point_locator_tolerance(1e-12);
            libMesh::Point center_top(10.0,10.0);
            mesh_fcn(center_top, 0.0, dX_center_top_serial);

            //mesh.sub_point_locator()->set_close_to_point_tol(1e-15);
            dX_center_top_parallel(0) = dX_system.point_value(0, center_top);
            dX_center_top_parallel(1) = dX_system.point_value(1, center_top);

            System& U_system = equation_systems->get_system<System>(IBFEMethod::VELOCITY_SYSTEM_NAME);
            umax = U_system.solution->linfty_norm();

            double u_sum = umax + umax_old;

            if( u_sum < 1e-5 && iteration_num > 1000 ) stop = true;
            umax_old = umax;
            //System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
            PetscVector<double>* X_ghost_vec = dynamic_cast<PetscVector<double>*>(
              fe_data_manager->buildGhostedSolutionVector(IBFEMethod::COORDS_SYSTEM_NAME, true));

            //const double volume = hier_math_ops.getVolumeOfPhysicalDomain();
            const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
            const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

            // get interpolated pressure field on FE mesh
            System& pf_system = equation_systems->get_system<System>("p_f system");
            NumericVector<double>* pf_vec = pf_system.solution.get();
            fe_data_manager->interp(p_copy_idx,
                                    *pf_vec,
                                    *X_ghost_vec,
                                    "p_f system",
                                    p_interp_spec,
                                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                                    loop_time);

            Pointer<CellVariable<NDIM, double> > u_cc_var = u_var;
            if (u_cc_var)
            {
                HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            /*    pout << "Error in u at time " << loop_time << ":\n"
                     << "  L1-norm:  "
                     << std::setprecision(10) << hier_cc_data_ops.L1Norm(u_idx, wgt_cc_idx) << "\n"
                     << "  L2-norm:  " << hier_cc_data_ops.L2Norm(u_idx, wgt_cc_idx) << "\n"
                     << "  max-norm: " << hier_cc_data_ops.maxNorm(u_idx, wgt_cc_idx) << "\n";
            */
                     u_err[0] = hier_cc_data_ops.L1Norm(u_idx, wgt_cc_idx);
                     u_err[1] = hier_cc_data_ops.L2Norm(u_idx, wgt_cc_idx);
                     u_err[2] = hier_cc_data_ops.maxNorm(u_idx, wgt_cc_idx);
            }

            Pointer<SideVariable<NDIM, double> > u_sc_var = u_var;
            if (u_sc_var)
            {
                HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            /*    pout << "Error in u at time " << loop_time << ":\n"
                     << "  L1-norm:  "
                     << std::setprecision(10) << hier_sc_data_ops.L1Norm(u_idx, wgt_sc_idx) << "\n"
                     << "  L2-norm:  " << hier_sc_data_ops.L2Norm(u_idx, wgt_sc_idx) << "\n"
                     << "  max-norm: " << hier_sc_data_ops.maxNorm(u_idx, wgt_sc_idx) << "\n";
            */
                     u_err[0] = hier_sc_data_ops.L1Norm(u_idx, wgt_sc_idx);
                     u_err[1] = hier_sc_data_ops.L2Norm(u_idx, wgt_sc_idx);
                     u_err[2] = hier_sc_data_ops.maxNorm(u_idx, wgt_sc_idx);
            }

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator->setupPlotData();
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
                ib_method_ops->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
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
        }


        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

     } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return 0;
}


void
output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
            Mesh& mesh,
            EquationSystems* equation_systems,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    // Write Cartesian data.
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
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
    file_name = data_dump_dirname + "/" + "fe_mesh.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    file_name += ".xda";
    mesh.write(file_name);
    file_name = data_dump_dirname + "/" + "fe_equation_systems.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    equation_systems->write(file_name, (EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA));
    return;
} // output_data
