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

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

#include <libmesh/explicit_system.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/parmetis_partitioner.h>
#include <libmesh/hilbert_sfc_partitioner.h>
#include <libmesh/metis_partitioner.h>

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

double structure_volume(EquationSystems *equation_systems)
{
      double J_integral = 0.0;
      System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
      NumericVector<double>* X_vec = X_system.solution.get();
      NumericVector<double>* X_ghost_vec = X_system.current_local_solution.get();
      X_vec->localize(*X_ghost_vec);
      pout << "L2 norm: " << X_vec->l2_norm() << '\n';
      pout << "size: " << X_vec->size() << '\n';
      DofMap& X_dof_map = X_system.get_dof_map();
      std::vector<vector<unsigned int> > X_dof_indices(NDIM);
      libMesh::UniquePtr<FEBase> fe(FEBase::build(NDIM, X_dof_map.variable_type(0)));
      libMesh::UniquePtr<QBase> qrule = QBase::build(QGAUSS, NDIM, FIFTH);
      fe->attach_quadrature_rule(qrule.get());
      const std::vector<double>& JxW = fe->get_JxW();
      const std::vector<vector<VectorValue<double> > >& dphi = fe->get_dphi();
      TensorValue<double> FF;
      boost::multi_array<double, 2> X_node;
      const MeshBase &mesh = equation_systems->get_mesh();
      MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
      MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
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

    return SAMRAI_MPI::sumReduction(J_integral);
}


void create_mesh(const double ds,
                 const std::string elem_type,
                 UnstructuredMesh &mesh)
{
  const double R = 0.2;
  if (NDIM == 2 && (elem_type == "TRI3" || elem_type == "TRI6"))
    {
#ifdef LIBMESH_HAVE_TRIANGLE
      const int num_circum_nodes = ceil(2.0 * M_PI * R / ds);
      for (int k = 0; k < num_circum_nodes; ++k)
        {
          const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_circum_nodes);
          mesh.add_point(libMesh::Point(R * cos(theta), R * sin(theta)));
        }
      TriangleInterface triangle(mesh);
      triangle.triangulation_type() = TriangleInterface::GENERATE_CONVEX_HULL;
      triangle.elem_type() = Utility::string_to_enum<ElemType>(elem_type);
      triangle.desired_area() = 1.5 * sqrt(3.0) / 4.0 * ds * ds;
      triangle.insert_extra_points() = true;
      triangle.smooth_after_generating() = true;
      triangle.triangulate();
#else
      TBOX_ERROR("ERROR: libMesh appears to have been configured without support for Triangle,\n"
                 << "       but Triangle is required for TRI3 or TRI6 elements.\n");
#endif
    }
  else
    {
      // NOTE: number of segments along boundary is 4*2^r.
      const double num_circum_segments = 2.0 * M_PI * R / ds;
      const int r = log2(0.25 * num_circum_segments);
      MeshTools::Generation::build_sphere(mesh, R, r, Utility::string_to_enum<ElemType>(elem_type));
    }

  // Ensure nodes on the surface are on the analytic boundary.
  MeshBase::element_iterator el_end = mesh.elements_end();
  for (MeshBase::element_iterator el = mesh.elements_begin(); el != el_end; ++el)
    {
      Elem* const elem = *el;
      for (unsigned int side = 0; side < elem->n_sides(); ++side)
        {
          const bool at_mesh_bdry = !elem->neighbor(side);
          if (!at_mesh_bdry) continue;
          for (unsigned int k = 0; k < elem->n_nodes(); ++k)
            {
              if (!elem->is_node_on_side(k, side)) continue;
              Node& n = *elem->get_node(k);
              n = R * n.unit();
            }
        }
    }
}



void
output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
            MeshBase& mesh,
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


// Elasticity model data.
namespace ModelData
{
// Coordinate mapping function.
void
coordinate_mapping_function(libMesh::Point& X, const libMesh::Point& s, void* /*ctx*/)
{
    X(0) = s(0) + 0.6;
    X(1) = s(1) + 0.5;
#if (NDIM == 3)
    X(2) = s(2) + 0.5;
#endif
    return;
} // coordinate_mapping_function

// Stress tensor functions.
static double c1_s = 0.05;
static double p0_s = 0.0;
static double beta_s = 0.0;
void
PK1_dev_stress_function(TensorValue<double>& PP,
                        const TensorValue<double>& FF,
                        const libMesh::Point& /*X*/,
                        const libMesh::Point& /*s*/,
                        Elem* const /*elem*/,
                        const vector<const vector<double>*>& /*var_data*/,
                        const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                        double /*time*/,
                        void* /*ctx*/)
{
    PP = 2.0 * c1_s * FF;
    return;
} // PK1_dev_stress_function

void
PK1_dil_stress_function(TensorValue<double>& PP,
                        const TensorValue<double>& FF,
                        const libMesh::Point& /*X*/,
                        const libMesh::Point& /*s*/,
                        Elem* const /*elem*/,
                        const vector<const vector<double>*>& /*var_data*/,
                        const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                        double /*time*/,
                        void* /*ctx*/)
{
    PP = 2.0 * (-p0_s + beta_s * log(FF.det())) * tensor_inverse_transpose(FF, NDIM);
    return;
} // PK1_dil_stress_function
}
using namespace ModelData;

// Function prototypes
class Solver
{
public:
  Solver(int argc, char **argv, const LibMeshInit &init);

  ~Solver();

  void run();

protected:
  AppInitializer app_initializer;
  Database &input_db;

  // Parameters controlling visualization output.
  const bool dump_viz_data;
  const int viz_dump_interval;
  const bool uses_visit;
  const bool uses_exodus;
  const std::string exodus_filename;
  const bool dump_postproc_data;
  const int postproc_data_dump_interval;
  const std::string postproc_data_dump_dirname;

  // Structural (finite element) data.
  ReplicatedMesh mesh;

  // Fundamental SAMRAI data objects.
  Pointer<CartesianGridGeometry<NDIM> > grid_geometry;
  Pointer<PatchHierarchy<NDIM> > patch_hierarchy;
  Pointer<StandardTagAndInitialize<NDIM> > error_detector;
  Pointer<BergerRigoutsos<NDIM> > box_generator;
  Pointer<LoadBalancer<NDIM> > load_balancer;
  Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm;

  // Fluid data.
  Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
  std::vector<RobinBcCoefStrategy<NDIM> *> u_bc_coefs;

  // coupled data.
  Pointer<IBFEMethod> ib_method_ops;
  Pointer<IBHierarchyIntegrator> time_integrator;

  // Input-output objects.
  std::unique_ptr<IBFEPostProcessor> postprocessor;
  std::unique_ptr<ExodusII_IO> exodus_io;

  void setup_eulerian_data();

  void setup_lagrangian_data();

  void setup_coupled_data();

  void setup_output_writers();
};

Solver::Solver(int argc, char **argv, const LibMeshInit &init)
  : app_initializer(argc, argv, "IB.log")
  , input_db(*app_initializer.getInputDatabase())
  , dump_viz_data(app_initializer.dumpVizData())
  , viz_dump_interval(app_initializer.getVizDumpInterval())
  , uses_visit(dump_viz_data && app_initializer.getVisItDataWriter())
#ifdef LIBMESH_HAVE_EXODUS_API
  , uses_exodus(dump_viz_data && !app_initializer.getExodusIIFilename().empty())
#else
  , uses_exodus(false)
#endif
  , exodus_filename(app_initializer.getExodusIIFilename())
  , dump_postproc_data(app_initializer.dumpPostProcessingData())
  , postproc_data_dump_interval(app_initializer.getPostProcessingDataDumpInterval())
  , postproc_data_dump_dirname(app_initializer.getPostProcessingDataDumpDirectory())
  , mesh(init.comm(), NDIM)

{
#ifndef LIBMESH_HAVE_EXODUS_API
  if (!app_initializer.getExodusIIFilename().empty())
    {
      plog << "WARNING: libMesh was compiled without Exodus support, so no "
           << "Exodus output will be written in this program.\n";
    }
#endif

  if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
    {
      Utilities::recursiveMkdir(postproc_data_dump_dirname);
    }

  // Setup static model data. TODO: move these into context objects.
  c1_s = input_db.getDouble("C1_S");
  p0_s = input_db.getDouble("P0_S");
  beta_s = input_db.getDouble("BETA_S");
}


Solver::~Solver()
{
  for (auto ptr : u_bc_coefs)
    delete ptr;
}



void
Solver::setup_eulerian_data()
{
  // convenience lambda
  auto get_comp_db = [&](const std::string &component)
    {
      return app_initializer.getComponentDatabase(component);
    };

  // First part: set up the basic SAMRAI objects:
  {
    grid_geometry = new CartesianGridGeometry<NDIM>("CartesianGeometry",
                                                    get_comp_db("CartesianGeometry"));
    patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    box_generator = new BergerRigoutsos<NDIM>();
    load_balancer = new LoadBalancer<NDIM>("LoadBalancer",
                                           get_comp_db("LoadBalancer"));
  }


  // Second part: set up the relevant IBAMR object:
  {
    const std::string solver_type = get_comp_db("Main")->getString("solver_type");
    if (solver_type == "STAGGERED")
    {
        navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            get_comp_db("INSStaggeredHierarchyIntegrator"));
    }
    else if (solver_type == "COLLOCATED")
    {
        navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
            "INSCollocatedHierarchyIntegrator",
            get_comp_db("INSCollocatedHierarchyIntegrator"));
    }
    else
    {
        TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                               << "Valid options are: COLLOCATED, STAGGERED");
    }

    // Create Eulerian initial condition specification objects.
    if (input_db.keyExists("VelocityInitialConditions"))
      {
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction
          ("u_init", get_comp_db("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);
      }

    if (input_db.keyExists("PressureInitialConditions"))
      {
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction
          ("p_init", get_comp_db("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);
      }

    // Create Eulerian boundary condition specification objects (when necessary).
    const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
    // TODO what does this do? Double check this part.
    if (periodic_shift.min() <= 0)
      {
        for (unsigned int d = 0; d < NDIM; ++d)
          {
            std::ostringstream bc_coefs_name_stream;
            bc_coefs_name_stream << "u_bc_coefs_" << d;
            const std::string bc_coefs_name = bc_coefs_name_stream.str();

            std::ostringstream bc_coefs_db_name_stream;
            bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
            const std::string bc_coefs_db_name = bc_coefs_db_name_stream.str();

            u_bc_coefs.push_back
              (new muParserRobinBcCoefs(bc_coefs_name, get_comp_db(bc_coefs_db_name),
                                        grid_geometry));
          }

        navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
      }
  }
}



void
Solver::setup_lagrangian_data()
{
  // convenience lambda
  auto get_comp_db = [&](const std::string &component)
    {
      return app_initializer.getComponentDatabase(component);
    };

  // Set up the mesh:
  {
    const double dx = input_db.getDouble("DX");
    const double ds = input_db.getDouble("MFAC") * dx;
    const std::string elem_type = input_db.getString("ELEM_TYPE");
    create_mesh(ds, elem_type, mesh);

    // TODO: we should remove this once we get the partitioning code working.
    mesh.prepare_for_use();
  }

  // Note that even though this takes a pointer to the mesh, ib_method_ops
  // does not read any libMesh data until we get to
  // initializeFEEquationSystems.
  {
    ib_method_ops = new IBFEMethod("IBFEMethod",
                                   get_comp_db("IBFEMethod"),
                                   &mesh,
                                   get_comp_db("GriddingAlgorithm")->getInteger("max_levels"));

    ib_method_ops->registerInitialCoordinateMappingFunction(coordinate_mapping_function);
    IBFEMethod::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function);
    IBFEMethod::PK1StressFcnData PK1_dil_stress_data(PK1_dil_stress_function);
    PK1_dev_stress_data.quad_order =
      Utility::string_to_enum<libMesh::Order>(input_db.getStringWithDefault("PK1_DEV_QUAD_ORDER",
                                                                            "THIRD"));
    PK1_dil_stress_data.quad_order =
      Utility::string_to_enum<libMesh::Order>(input_db.getStringWithDefault("PK1_DIL_QUAD_ORDER",
                                                                            "FIRST"));
    ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data);
    ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data);
    if (input_db.getBoolWithDefault("ELIMINATE_PRESSURE_JUMPS", false))
      {
        ib_method_ops->registerStressNormalizationPart();
      }

    // TODO: we will have to move this call to somewhere else (since patch
    // data is not available until we call
    // time_integrator->initializePatchHierarchy()) once we get the SAMRAI
    // partitioner working.
    //
    // Note that initializeFEData immediately calls
    // initializeFEEquationSystems.
    //
    // TODO: we cannot call initializeFEData at this point: we *have* to call
    // that after the postprocessor is set up, since the postprocessor adds
    // more variables to the systems
    ib_method_ops->initializeFEEquationSystems();
  }
}

void
Solver::setup_coupled_data()
{
  // convenience lambda
  auto get_comp_db = [&](const std::string &component)
    {
      return app_initializer.getComponentDatabase(component);
    };

  time_integrator = new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                                      get_comp_db("IBHierarchyIntegrator"),
                                                      ib_method_ops,
                                                      navier_stokes_integrator);

  // Create Eulerian body force function specification objects.
  if (input_db.keyExists("ForcingFunction"))
    {
      Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction
        ("f_fcn", get_comp_db("ForcingFunction"), grid_geometry);
      time_integrator->registerBodyForceFunction(f_fcn);
    }

  error_detector = new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                                      time_integrator,
                                                      get_comp_db("StandardTagAndInitialize"));

  gridding_algorithm = new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                                   get_comp_db("GriddingAlgorithm"),
                                                   error_detector,
                                                   box_generator,
                                                   load_balancer);

  // Now that we have the gridding algorithm we can initialize patches and
  // (TODO) partition the Lagrangian mesh:

  // The hierarchy initialization also registers quantities with the VisIt
  // writer, so we have to register that here before initialization:
  if (uses_visit)
  {
      time_integrator->registerVisItDataWriter(app_initializer.getVisItDataWriter());
  }

  time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
}



void
Solver::setup_output_writers()
{
  // Set up visualization plot file writers.
  if (uses_exodus)
  {
      exodus_io = std::unique_ptr<ExodusII_IO>(new ExodusII_IO(mesh));
  }

  // Configure the postprocessor.
  postprocessor = std::unique_ptr<IBFEPostProcessor>
    (new IBFECentroidPostProcessor("IBFEPostProcessor", ib_method_ops->getFEDataManager()));

  // Add postprocessing actions for stress tensors.
  postprocessor->registerTensorVariable("FF", MONOMIAL, CONSTANT, IBFEPostProcessor::FF_fcn);
  static std::pair<IBTK::TensorMeshFcnPtr, void*> PK1_dev_stress_fcn_data
    (PK1_dev_stress_function, static_cast<void*>(NULL));

  postprocessor->registerTensorVariable("sigma_dev",
                                        MONOMIAL,
                                        CONSTANT,
                                        IBFEPostProcessor::cauchy_stress_from_PK1_stress_fcn,
                                        std::vector<SystemData>(),
                                        &PK1_dev_stress_fcn_data);

  static std::pair<IBTK::TensorMeshFcnPtr, void*> PK1_dil_stress_fcn_data
    (PK1_dil_stress_function, static_cast<void*>(NULL));
  postprocessor->registerTensorVariable("sigma_dil",
                                        MONOMIAL,
                                        CONSTANT,
                                        IBFEPostProcessor::cauchy_stress_from_PK1_stress_fcn,
                                        std::vector<SystemData>(),
                                        &PK1_dil_stress_fcn_data);

  // and an interpolated pressure.
  Pointer<hier::Variable<NDIM> > p_var = navier_stokes_integrator->getPressureVariable();
  Pointer<VariableContext> p_current_ctx = navier_stokes_integrator->getCurrentContext();
  HierarchyGhostCellInterpolation::InterpolationTransactionComponent p_ghostfill
    (/*data_idx*/ -1, "LINEAR_REFINE", /*use_cf_bdry_interpolation*/ false,
     "CONSERVATIVE_COARSEN", "LINEAR");
  FEDataManager::InterpSpec p_interp_spec("PIECEWISE_LINEAR",
                                          QGAUSS,
                                          FIFTH,
                                          /*use_adaptive_quadrature*/ false,
                                          /*point_density*/ 2.0,
                                          /*use_consistent_mass_matrix*/ true,
                                          /*use_nodal_quadrature*/ false);
  postprocessor->registerInterpolatedScalarEulerianVariable
    ("p_f", LAGRANGE, FIRST, p_var, p_current_ctx, p_ghostfill, p_interp_spec);
  // TODO: Clean up the initialization so that this call is not necessary.
  ib_method_ops->initializeFEData();
  // this is needed
  postprocessor->initializeFEData();
}

void
Solver::run()
{
  setup_eulerian_data();

  setup_lagrangian_data();

  // TODO: we have to satisfy the following constraints:
  //
  // - The postprocessor requires that any Eulerian variables (that are to be
  //   postprocessed) exist, but we do not need patches to be allocated.
  //
  // - The postprocessor needs ib_method_ops to have called initializeFEData
  //   (and initializeFEEquationSystems). However, since it adds more system
  //   variables, we have to call ib_method_ops->initializeFEData again before
  //   calling postprocessor->initializeFEData().
  //
  // - time_integrator->initializePatchHierarchy requires ib_method_ops to
  //   have called initializeFEData (which completely sets it up)
  //
  setup_output_writers();

  setup_coupled_data();

  plog << "Input database:\n";
  input_db.printClassData(plog);

  // Write out initial visualization data.
  int iteration_num = time_integrator->getIntegratorStep();
  double loop_time = time_integrator->getIntegratorTime();
  EquationSystems* equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();
  if (dump_viz_data)
    {
      pout << "\n\nWriting visualization files...\n\n";
      if (uses_visit)
        {
          time_integrator->setupPlotData();
          app_initializer.getVisItDataWriter()->writePlotData(patch_hierarchy,
                                                              iteration_num, loop_time);
        }
      if (uses_exodus)
        {
          postprocessor->postProcessData(loop_time);
          exodus_io->write_timestep
            (exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1,
             loop_time);
        }
    }

  // Open streams to save volume of structure.
  std::ofstream volume_stream;
  if (SAMRAI_MPI::getRank() == 0)
    {
      volume_stream.open("volume.curve", ios_base::out | ios_base::trunc);
    }

  // Main time step loop.
  const double loop_time_end = time_integrator->getEndTime();
  double dt = 0.0;
  while (!MathUtilities<double>::equalEps(loop_time, loop_time_end)
         && time_integrator->stepsRemaining())
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

      // At specified intervals, write visualization and restart files, print
      // out timer data, and store hierarchy data for post processing.
      iteration_num += 1;
      const bool last_step = !time_integrator->stepsRemaining();
      if (dump_viz_data && (iteration_num % viz_dump_interval == 0 || last_step))
        {
          pout << "\nWriting visualization files...\n\n";
          if (uses_visit)
            {
              time_integrator->setupPlotData();
              app_initializer.getVisItDataWriter()->writePlotData(patch_hierarchy,
                                                                  iteration_num, loop_time);
            }
          if (uses_exodus)
            {
              postprocessor->postProcessData(loop_time);
              exodus_io->write_timestep
                (exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1,
                 loop_time);
            }
        }
      if (app_initializer.dumpRestartData() && (iteration_num % app_initializer.getRestartDumpInterval() == 0 || last_step))
        {
          pout << "\nWriting restart files...\n\n";
          RestartManager::getManager()->writeRestartFile(app_initializer.getRestartDumpDirectory(), iteration_num);
        }
      if (app_initializer.dumpTimerData() && (iteration_num % app_initializer.getTimerDumpInterval() == 0 || last_step))
        {
          pout << "\nWriting timer data...\n\n";
          TimerManager::getManager()->print(plog);
        }
      if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
        {
          output_data(patch_hierarchy,
                      navier_stokes_integrator,
                      mesh,
                      equation_systems,
                      iteration_num,
                      loop_time,
                      postproc_data_dump_dirname);
        }

      // Compute the volume of the structure.
      const double J_integral = structure_volume(equation_systems);
      if (SAMRAI_MPI::getRank() == 0)
        {
          volume_stream.precision(12);
          volume_stream.setf(ios::fixed, ios::floatfield);
          volume_stream << loop_time << " " << J_integral << endl;
        }
    }
}



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

bool run_example(int argc, char** argv)
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    {
      Solver solver(argc, argv, init);

      solver.run();
    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return true;
} // run_example
