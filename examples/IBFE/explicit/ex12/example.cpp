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

void sigma_xx_fcn(double& F,
           const libMesh::TensorValue<double>& FF,
           const libMesh::Point& /*x*/,
           const libMesh::Point& /*X*/,
           libMesh::Elem* /*elem*/,
           const std::vector<const std::vector<double>*>& /*system_var_data*/,
           const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
           double /*data_time*/,
           void* /*ctx*/);

void sigma_xy_fcn(double& F,
           const libMesh::TensorValue<double>& FF,
           const libMesh::Point& /*x*/,
           const libMesh::Point& /*X*/,
           libMesh::Elem* /*elem*/,
           const std::vector<const std::vector<double>*>& /*system_var_data*/,
           const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
           double /*data_time*/,
           void* /*ctx*/);

void sigma_yx_fcn(double& F,
           const libMesh::TensorValue<double>& FF,
           const libMesh::Point& /*x*/,
           const libMesh::Point& /*X*/,
           libMesh::Elem* /*elem*/,
           const std::vector<const std::vector<double>*>& /*system_var_data*/,
           const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
           double /*data_time*/,
           void* /*ctx*/);

void sigma_yy_fcn(double& F,
           const libMesh::TensorValue<double>& FF,
           const libMesh::Point& /*x*/,
           const libMesh::Point& /*X*/,
           libMesh::Elem* /*elem*/,
           const std::vector<const std::vector<double>*>& /*system_var_data*/,
           const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
           double /*data_time*/,
           void* /*ctx*/);

// Elasticity model data.
namespace ModelData
{
    // Problem parameters.
    double mu_e;
    double nu;
    double Lambda;
    double kappa_tether;
    double eta_tether;
    int PK1_dev_flag;
    
    // surface pressure function parameters
    double P_load;
    double t_load;
    double period;
    BoundaryInfo* boundary_info;
    
    void compute_deviatoric_projection(TensorValue<double>& PP,
                                       const TensorValue<double>& FF)
    {
        // compute deviatoric projection
        const TensorValue<double> Sigma = (1.0/FF.det()) * PP * FF.transpose();
        const TensorValue<double> Eye (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0); 
        const double pressure = (1.0/3.0) * (Sigma(0,0) + Sigma(1,1) + Sigma(2,2));
        PP = FF.det() * (Sigma - pressure * Eye) * FF.transpose().inverse();
        
        return;
    }    
    
    double loading_pressure(double time)
    {   
        return P_load * min(time, t_load) / t_load;
    }
        
    // Coordinate mapping function.                                                                                                   
    void
    coordinate_mapping_function(libMesh::Point& X, const libMesh::Point& s, void* /*ctx*/)
    {
        X(0) = s(0) + 5;
        X(1) = s(1) + 10;
#if (NDIM == 3)
        X(2) = s(2) + 0.4;
#endif
        return;
    } // coordinate_mapping_function
    
    // Stress tensor functions.
    void
    PK1_dev_stress_function_1(TensorValue<double>& PP,
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
        PP = mu_e * ( FF - FF_inv_trans );
             
        return;
    } // PK1_dev_stress_function_1
    
    void
    PK1_dev_stress_function_2(TensorValue<double>& PP,
            const TensorValue<double>& FF,
            const libMesh::Point& /*X*/,
            const libMesh::Point& /*s*/,
            Elem* const /*elem*/,
            const std::vector<const std::vector<double>*>& /*var_data*/,
            const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
            double /*time*/,
            void* /*ctx*/)
    {
        PP = mu_e * FF ;
        
        return;
    } // PK1_dev_stress_function_2
    
    // PP = mu_e * d I1_bar / d FF
    void
    PK1_dev_stress_function_3(TensorValue<double>& PP,
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
         double I1 = (FF.transpose() * FF).tr();
         
         PP = mu_e * pow(J, -2.0 / 3.0) * (FF - (1.0 / 3.0) * I1 * tensor_inverse_transpose(FF));        
         
         return;
    } // PK1_dev_stress_function_3
    
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
              
        if (!MathUtilities<double>::equalEps(Lambda, 0.0))
        {
            //PP = 0.5 * Lambda * FF_inv_trans * (J*J - 1);
            PP = Lambda * log(J) * FF_inv_trans;
        }

        return;
    } // PK1_dil_stress_function
    
    // surface pressure functions
    void loading_force_function(double& P,
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
        if(! ( 5.0 - 1.0e-15 < X(0) && X(0) < 15.0 + 1.0e-15 ) )
        {
            P = 0.0; 
            return;
        }
        
        const vector<short int>& bdry_ids = boundary_info->boundary_ids(elem, side);
        std::vector<short int> ids_for_BCs;
        ids_for_BCs.push_back(2);
        if (find_first_of(bdry_ids.begin(), bdry_ids.end(), ids_for_BCs.begin(), ids_for_BCs.end()) != bdry_ids.end())
        {
            P = loading_pressure(data_time);
        }
        else
        {
            P = 0.0;
        }
        return;
    }
    
    void loading_force_function_smooth(double& P,
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
        const double a = 4.0;
        const double b = 16.0;
        
        if(! ( a - 1.0e-15 < X(0) && X(0) < b + 1.0e-15 ) )
        {
            P = 0.0; 
            return;
        }
        
        const vector<short int>& bdry_ids = boundary_info->boundary_ids(elem, side);
        std::vector<short int> ids_for_BCs;
        ids_for_BCs.push_back(2);
        if (find_first_of(bdry_ids.begin(), bdry_ids.end(), ids_for_BCs.begin(), ids_for_BCs.end()) != bdry_ids.end())
        {
            // using a bump function to define the pressure
	  P = loading_pressure(data_time)*exp(1/(pow((2*X(0) - a - b)/(b - a), 2.0) - 1.0))/exp(-1.0);
        }
        else
        {
            P = 0.0;
        }
        return;
    }
    
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
        
        // get velocity
        
        libMesh::Point velocity((*system_var_data[0])[0], (*system_var_data[0])[1]);
               
        F.zero();
        
        libMesh::Point transformed_X;
        
        coordinate_mapping_function(transformed_X, X, NULL);
        
        const vector<short int>& bdry_ids = boundary_info->boundary_ids(elem, side);
        std::vector<short int> ids_for_BCs;

        // boundary conditions for bottom
        ids_for_BCs.push_back(0);
        if (find_first_of(bdry_ids.begin(), bdry_ids.end(), ids_for_BCs.begin(), ids_for_BCs.end()) != bdry_ids.end())
        {
            F = kappa_tether * (transformed_X - x) + eta_tether * velocity;
            F(0) = 0.0; // zero out the horizontal component of the restoring force
        }
        
        // boundary conditions for top
        ids_for_BCs.clear();   
        ids_for_BCs.push_back(2);
        if (find_first_of(bdry_ids.begin(), bdry_ids.end(), ids_for_BCs.begin(), ids_for_BCs.end()) != bdry_ids.end())
        {
            F = kappa_tether * (transformed_X - x) + eta_tether * velocity;;
            F(1) = 0.0; // zero out the vertical component of the restoring force
        }
     
        return;
    
    } 
    
    void trace_sigma_dev_fcn(double& one_third_trace_sigma,
                       const libMesh::TensorValue<double>& FF,
                       const libMesh::Point& x,
                       const libMesh::Point& X,
                       libMesh::Elem* elem,
                       const std::vector<const std::vector<double>*>& system_var_data,
                       const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
                       double data_time,
                       void* ctx)
{
        
        TensorValue<double> PP;
        TensorValue<double> Sigma;
        PP.zero();
        Sigma.zero();
        
        if(PK1_dev_flag == 1)
        {
            PK1_dev_stress_function_1(PP, FF, x, X, elem, system_var_data, system_grad_var_data, data_time, ctx);
        }
        else if(PK1_dev_flag == 2)
        {
            PK1_dev_stress_function_2(PP, FF, x, X, elem, system_var_data, system_grad_var_data, data_time, ctx);
        }
        else if(PK1_dev_flag == 3)
        {
            PK1_dev_stress_function_3(PP, FF, x, X, elem, system_var_data, system_grad_var_data, data_time, ctx);
        }
        else
        {
            PK1_dev_stress_function_1(PP, FF, x, X, elem, system_var_data, system_grad_var_data, data_time, ctx);
        }
        
        Sigma = (1.0/FF.det()) * PP * FF.transpose();
        one_third_trace_sigma = (1.0/3.0) * Sigma.tr();
            
    return;
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
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();
           
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

        if(input_db->getInteger("N") == 8)
        {
            std::string foo("_N_8.ex2");
            exodus_filename.append(foo);
        }
        if(input_db->getInteger("N") == 16)
        {
            std::string foo("_N_16.ex2");
            exodus_filename.append(foo);
        }
        if(input_db->getInteger("N") == 32)
        {
            std::string foo("_N_32.ex2");
            exodus_filename.append(foo);        
        }
        if(input_db->getInteger("N") == 64)
        {
            std::string foo("_N_64.ex2");
            exodus_filename.append(foo);        
        }
        if(input_db->getInteger("N") == 128)
        {
            std::string foo("_N_128.ex2");
            exodus_filename.append(foo);        
        }
        if(input_db->getInteger("N") == 256)
        {
            std::string foo("_N_256.ex2");
            exodus_filename.append(foo);        
        }
              
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
        Mesh mesh(init.comm(), NDIM);
        const double dx = input_db->getDouble("DX");
        const double MFAC = input_db->getDouble("MFAC");
        const double ds = MFAC * dx;
        string elem_type = input_db->getString("ELEM_TYPE");
        const int n_x = round( 20.0 / ds );
        const int n_y = round( 10.0 / ds );
        MeshTools::Generation::build_square(
            mesh, n_x, n_y, 0.0, 20.0, 0.0, 10.0, Utility::string_to_enum<ElemType>(elem_type));

        // get boundary information and surface pressure info
        boundary_info = &mesh.get_boundary_info();
        P_load = input_db->getDouble("P_load");
        t_load = input_db->getDouble("t_load");
        mu_e = input_db->getDouble("mu_e");
        nu = input_db->getDouble("nu");
        Lambda = 2.0*mu_e*(1.0 + nu)/(3*(1.0 - 2.0*nu));
        kappa_tether = input_db->getDouble("kappa_tether");
        eta_tether = input_db->getDouble("eta_tether");
        PK1_dev_flag = input_db->getInteger("PK1_DEV_STRESS_FLAG");
     
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
        ib_method_ops->initializeFEEquationSystems();
        FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager();
        ib_method_ops->registerInitialCoordinateMappingFunction(coordinate_mapping_function);
        if(input_db->getIntegerWithDefault("PK1_DEV_STRESS_FLAG",1) == 1)
        {
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function_1);
            PK1_dev_stress_data.quad_order =
                    Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "FIFTH"));
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data);
        }
        else if(input_db->getIntegerWithDefault("PK1_DEV_STRESS_FLAG",1) == 2)
        {
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function_2);
            PK1_dev_stress_data.quad_order =
                    Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "FIFTH"));
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data);
        }
        else if(input_db->getIntegerWithDefault("PK1_DEV_STRESS_FLAG",1) == 3)
        {
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function_3);
            PK1_dev_stress_data.quad_order =
                    Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "FIFTH"));
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data);
        }
        else
        {
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function_1);
            PK1_dev_stress_data.quad_order =
                    Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "FIFTH"));
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data);
        }
        IBFEMethod::PK1StressFcnData PK1_dil_stress_data(PK1_dil_stress_function);
        PK1_dil_stress_data.quad_order =
                Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DIL_QUAD_ORDER", "FIFTH"));
        ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data);
        
        // register surface pressure function
        if(input_db->getBoolWithDefault("SMOOTH_LOADING_FORCE", true))
        {            
            IBFEMethod::LagSurfacePressureFcnData surface_pressure_data(loading_force_function_smooth, sys_data);
            ib_method_ops->registerLagSurfacePressureFunction(surface_pressure_data, 0);
        }
        else
        {
            IBFEMethod::LagSurfacePressureFcnData surface_pressure_data(loading_force_function, sys_data);
            ib_method_ops->registerLagSurfacePressureFunction(surface_pressure_data, 0);
        }
        
        // register tether force function
        IBFEMethod::LagSurfaceForceFcnData surface_tether_force_data(surface_tether_force_function, sys_data);
        ib_method_ops->registerLagSurfaceForceFunction(surface_tether_force_data, 0);
        
        // setup libmesh things for eliminating pressure jumps
        if (input_db->getBoolWithDefault("ELIMINATE_PRESSURE_JUMPS", false))
        {
            ib_method_ops->registerStressNormalizationPart();
        }
        // setup libmesh things for eliminating pressure jumps
        if (input_db->getBoolWithDefault("DO_PRESSURE_STABILIZATION", false))
        {
            ib_method_ops->registerPressureStabilizationPart();
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
        
        // set up evaluation of the determinant of the deformation gradient
        ib_post_processor->registerScalarVariable("trace_sigma_dev", MONOMIAL, CONSTANT, &trace_sigma_dev_fcn);
        ib_post_processor->registerScalarVariable("sigma_xx", MONOMIAL, CONSTANT, &sigma_xx_fcn);
        ib_post_processor->registerScalarVariable("sigma_xy", MONOMIAL, CONSTANT, &sigma_xy_fcn);  
        ib_post_processor->registerScalarVariable("sigma_yx", MONOMIAL, CONSTANT, &sigma_yx_fcn);  
        ib_post_processor->registerScalarVariable("sigma_yy", MONOMIAL, CONSTANT, &sigma_yy_fcn);  
        
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
        
	if (SAMRAI_MPI::getRank() == 0)
        {
	  if(nu == -1.0)
	    {
	      //         volume_stream.open("volume_nu=m1.curve", ios_base::out | ios_base::trunc);
	      //         stuff_stream.open("stuff_nu=m1.dat", ios_base::out | ios_base::app);
	      //         dX_stream.open("dX_nu=m1.dat", ios_base::out | ios_base::app);
	      disp_stream.open("disp_nu=m1.dat", ios_base::out | ios_base::app);
	      pressure_field_stream.open("pressure_field_nu=m1.dat", ios_base::out | ios_base::app);
	      vel_error_stream.open("vel_error_nu=m1.dat", ios_base::out | ios_base::app);
	    }
	  if(nu == 0.0)
	    {
	      //         volume_stream.open("volume_nu=0.curve", ios_base::out | ios_base::trunc);
	      //         stuff_stream.open("stuff_nu=0.dat", ios_base::out | ios_base::app);
	      //         dX_stream.open("dX_nu=0.dat", ios_base::out | ios_base::app);
	      disp_stream.open("disp_nu=0.dat", ios_base::out | ios_base::app);
	      pressure_field_stream.open("pressure_field_nu=0.dat", ios_base::out | ios_base::app);
	      vel_error_stream.open("vel_error_nu=0.dat", ios_base::out | ios_base::app);
	    }
	  if(nu == 0.4)
	    {
	      //          volume_stream.open("volume_nu=0p4.curve", ios_base::out | ios_base::trunc);
	      //          stuff_stream.open("stuff_nu=0p4.dat", ios_base::out | ios_base::app);
	      //          dX_stream.open("dX_nu=0p4.dat", ios_base::out | ios_base::app);
	      disp_stream.open("disp_nu=0p4.dat", ios_base::out | ios_base::app);
	      pressure_field_stream.open("pressure_field_nu=0p4.dat", ios_base::out | ios_base::app);
	      vel_error_stream.open("vel_error_nu=0p4.dat", ios_base::out | ios_base::app);
	    }
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
                               
            if (SAMRAI_MPI::getRank() == 0)
            {
                dX_stream.precision(12);
                dX_stream.setf(ios::fixed, ios::floatfield);
                dX_stream << loop_time << " " << dX_center_top_serial(0) - 5.0 << " " << dX_center_top_serial(1) - 10.0 << endl;
            }
                     
            System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
            PetscVector<double>* X_ghost_vec = dynamic_cast<PetscVector<double>*>(
              fe_data_manager->buildGhostedSolutionVector(IBFEMethod::COORDS_SYSTEM_NAME, true));
            
            const double volume = hier_math_ops.getVolumeOfPhysicalDomain();
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
            
         
            // Compute the volume of the structure.
            /*double J_integral = 0.0;
            DofMap& X_dof_map = X_system.get_dof_map();
            std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
            UniquePtr<FEBase> fe(FEBase::build(NDIM, X_dof_map.variable_type(0)));
            UniquePtr<QBase> qrule = QBase::build(QGAUSS, NDIM, FIFTH);
            fe->attach_quadrature_rule(qrule.get());
            const std::vector<double>& JxW = fe->get_JxW();
            const std::vector<std::vector<VectorValue<double> > >& dphi = fe->get_dphi();
            TensorValue<double> FF;
            boost::multi_array<double, 2> X_node;
            const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
            const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
            double max_J = -std::numeric_limits<double>::max();
            double min_J = std::numeric_limits<double>::max();
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
                    min_J = std::min(min_J, FF.det());
                    max_J = std::max(max_J, FF.det());
                }
            }
            J_integral = SAMRAI_MPI::sumReduction(J_integral);
	    min_J = SAMRAI_MPI::minReduction(min_J);
	    max_J = SAMRAI_MPI::minReduction(max_J);
            if (SAMRAI_MPI::getRank() == 0)
            {
                volume_stream.precision(12);
                volume_stream.setf(ios::fixed, ios::floatfield);
                volume_stream << loop_time << " " << J_integral << endl;
                
                stuff_stream.precision(12);
                stuff_stream.setf(ios::fixed, ios::floatfield);
                stuff_stream << loop_time << " " << min_J << " " << max_J << endl;
            }
            */
        }
        
        // write values of the pressure field to a file.
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(0);
        const int pressure_idx = p_current_idx;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, double> > P_data = patch->getPatchData(pressure_idx);
                                                                   
                // Read the Cartesian grid values.
                for (Box<NDIM>::Iterator b(patch_box); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    pressure_field_stream << (*P_data)(i) << "\n";
                }
        }
        
        // right out final displacement to file
        disp_stream << dx << " " << dX_center_top_serial(0) - 5.0 << " " << dX_center_top_serial(1) - 10.0;
        disp_stream << " " << dX_center_top_parallel(0) - 5.0 << " " << dX_center_top_parallel(1) - 10.0 << std::endl;
      
        // write out velocity errors
        vel_error_stream << dx << " " << u_err[0] << " " << u_err[1] << " " << u_err[2] << "\n";
        
        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
        //    volume_stream.close();
        //    stuff_stream.close();
        //    dX_stream.close();
            disp_stream.close();
            pressure_field_stream.close();
            vel_error_stream.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        
     } // cleanup dynamically allocated objects prior to shutdown
    
    SAMRAIManager::shutdown();
    return 0;
} 

void sigma_xx_fcn(double& F,
           const libMesh::TensorValue<double>& FF,
           const libMesh::Point& /*x*/,
           const libMesh::Point& /*X*/,
           libMesh::Elem* /*elem*/,
           const std::vector<const std::vector<double>*>& /*system_var_data*/,
           const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
           double /*data_time*/,
           void* /*ctx*/)
{
    const double J = FF.det();
    double I1 = (FF.transpose() * FF).tr();
    TensorValue<double> PP = ModelData::mu_e * pow(J, -2.0 / 3.0) * (FF - (1.0 / 3.0) * I1 * tensor_inverse_transpose(FF));
    if (!MathUtilities<double>::equalEps(ModelData::Lambda, 0.0))
    {
        PP += ModelData::Lambda * log(J) * tensor_inverse_transpose(FF);
    }
    const TensorValue<double> sigma = (PP * FF.transpose()) / J;
    F = sigma(1,1);
    return;
}

void sigma_xy_fcn(double& F,
           const libMesh::TensorValue<double>& FF,
           const libMesh::Point& /*x*/,
           const libMesh::Point& /*X*/,
           libMesh::Elem* /*elem*/,
           const std::vector<const std::vector<double>*>& /*system_var_data*/,
           const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
           double /*data_time*/,
           void* /*ctx*/)
{
    const double J = FF.det();
    double I1 = (FF.transpose() * FF).tr();
    TensorValue<double> PP = ModelData::mu_e * pow(J, -2.0 / 3.0) * (FF - (1.0 / 3.0) * I1 * tensor_inverse_transpose(FF));
    if (!MathUtilities<double>::equalEps(ModelData::Lambda, 0.0))
    {
        PP += ModelData::Lambda * log(J) * tensor_inverse_transpose(FF);
    }
    const TensorValue<double> sigma = (PP * FF.transpose()) / J;
    F = sigma(1,2);
    return;
}

void sigma_yx_fcn(double& F,
           const libMesh::TensorValue<double>& FF,
           const libMesh::Point& /*x*/,
           const libMesh::Point& /*X*/,
           libMesh::Elem* /*elem*/,
           const std::vector<const std::vector<double>*>& /*system_var_data*/,
           const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
           double /*data_time*/,
           void* /*ctx*/)
{
    const double J = FF.det();
    double I1 = (FF.transpose() * FF).tr();
    TensorValue<double> PP = ModelData::mu_e * pow(J, -2.0 / 3.0) * (FF - (1.0 / 3.0) * I1 * tensor_inverse_transpose(FF));
    if (!MathUtilities<double>::equalEps(ModelData::Lambda, 0.0))
    {
        PP += ModelData::Lambda * log(J) * tensor_inverse_transpose(FF);
    }
    const TensorValue<double> sigma = (PP * FF.transpose()) / J;
    F = sigma(2,1);
    return;
}

void sigma_yy_fcn(double& F,
           const libMesh::TensorValue<double>& FF,
           const libMesh::Point& /*x*/,
           const libMesh::Point& /*X*/,
           libMesh::Elem* /*elem*/,
           const std::vector<const std::vector<double>*>& /*system_var_data*/,
           const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
           double /*data_time*/,
           void* /*ctx*/)
{
    const double J = FF.det();
    double I1 = (FF.transpose() * FF).tr();
    TensorValue<double> PP = ModelData::mu_e * pow(J, -2.0 / 3.0) * (FF - (1.0 / 3.0) * I1 * tensor_inverse_transpose(FF));
    if (!MathUtilities<double>::equalEps(ModelData::Lambda, 0.0))
    {
        PP += ModelData::Lambda * log(J) * tensor_inverse_transpose(FF);
    }
    const TensorValue<double> sigma = (PP * FF.transpose()) / J;
    F = sigma(2,2);
    return;
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
