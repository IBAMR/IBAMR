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
// #include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// other samrai stuff
#include <HierarchyDataOpsManager.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

// Headers for application-specific algorithm/data structure objects
#include <boost/multi_array.hpp>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>
#include <ibtk/BoxPartitioner.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// This file is the main driver for velocity interpolation tests (i.e.,
// IBFEmethod::interpolateVelocity).
//
// Due to limitations in how we set up and initialize the Eulerian and
// Lagrangian discretizations, the tests that utilize this executable only
// solve a single interpolation problem. We check the convergence rate by
// running several tests with different grid resolutions (e.g.,
// interpolate_velocity_01.a, interpolate_velocity_01.b, etc.).


// A scalar function parser that behaves similarly to our own
// muParserCartGridFunction.
class ParsedFunction
{
public:
    ParsedFunction(std::vector<std::string> expressions,
                   const unsigned int n_vars = 0)
        :
        string_functions(std::move(expressions)),
        vars(n_vars == 0 ? string_functions.size() : n_vars, 0.0),
        parsers(string_functions.size())
    {
        for (unsigned int var_n = 0; var_n < vars.size(); ++var_n)
            for (mu::Parser &parser : parsers)
                parser.DefineVar("X_" + std::to_string(var_n), &vars[var_n]);

        for (unsigned int p_n = 0; p_n < string_functions.size(); ++p_n)
            parsers[p_n].SetExpr(string_functions[p_n]);
    }

    // Evaluate a single component.
    double
    value(const libMesh::Point &p, const unsigned int component_n) const
    {
        TBOX_ASSERT(component_n < parsers.size());
        for (unsigned int var_n = 0; var_n < vars.size(); ++var_n)
            vars[var_n] = p(var_n);

        try
        {
            return parsers[component_n].Eval();
        }
        catch (mu::ParserError &e)
        {
            std::cerr << "Message:  <" << e.GetMsg() << ">\n";
            std::cerr << "Formula:  <" << e.GetExpr() << ">\n";
            std::cerr << "Token:    <" << e.GetToken() << ">\n";
            std::cerr << "Position: <" << e.GetPos() << ">\n";
            std::cerr << "Errc:     <" << e.GetCode() << ">" << std::endl;
            throw e;
        }
        return 1.0;
    }

    // Evaluate all components. Only makes sense if the number of variables is
    // equal to the number of functions.
    libMesh::Point
    value(const libMesh::Point &p) const
    {
        TBOX_ASSERT(vars.size() == parsers.size());
        for (unsigned int var_n = 0; var_n < vars.size(); ++var_n)
            vars[var_n] = p(var_n);

        libMesh::Point out;
        for (unsigned int component_n = 0; component_n < parsers.size(); ++component_n)
            out(component_n) = parsers[component_n].Eval();
        return out;
    }

private:
    const std::vector<std::string> string_functions;
    mutable std::vector<double> vars;

    std::vector<mu::Parser> parsers;
};



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


int main(int argc, char** argv)
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    PetscOptionsSetValue(nullptr, "-ksp_rtol", "1e-16");
    PetscOptionsSetValue(nullptr, "-ksp_atol", "1e-16");

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");

        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create a simple FE mesh.
        ReplicatedMesh mesh(init.comm(), NDIM);
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        string elem_type = input_db->getString("ELEM_TYPE");
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
        mesh.prepare_for_use();
        plog << "Number of elements: " << mesh.n_active_elem() << std::endl;

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"), false);
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>(
            "PatchHierarchy", grid_geometry, false);
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();

        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"), false);
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"), false);
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
                           false);
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator,
                                              false);

        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer,
                                        false);

        // Configure the IBFE solver.
        ib_method_ops->registerInitialCoordinateMappingFunction(coordinate_mapping_function);
        ib_method_ops->initializeFEEquationSystems();
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);

        // Now for the actual test. The stored velocity field does not contain
        // ghost data, so we set up a new context and velocity field that does:
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const Pointer<SAMRAI::hier::Variable<NDIM> > u_var = time_integrator->getVelocityVariable();
        const Pointer<VariableContext> u_ghost_ctx = var_db->getContext("u_ghost");
        const int u_ghost_idx = var_db->registerVariableAndContext(u_var, u_ghost_ctx, 3);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_ghost_idx);
        }
        u_init->setDataOnPatchHierarchy(u_ghost_idx, u_var, patch_hierarchy, 0.0);

        // synch ghost data
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ghost_cell_components(1);
        ghost_cell_components[0] = InterpolationTransactionComponent(u_ghost_idx,
                                                                     "CONSERVATIVE_LINEAR_REFINE",
                                                                     true,
                                                                     "CONSERVATIVE_COARSEN",
                                                                     "LINEAR",
                                                                     false,
                                                                     {}, // u_bc_coefs
                                                                     nullptr);
        HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ghost_cell_components, patch_hierarchy);
        ghost_fill_op.fillData(/*time*/0.0);

        // TODO what is the best way to set the linear solver tolerance in the
        // L2 projection, besides, e.g., setting -ksp_rtol 1e-12 at the
        // command line?
        const double dt = time_integrator->getMaximumTimeStepSize();
        time_integrator->preprocessIntegrateHierarchy(time_integrator->getIntegratorTime(),
                                                      time_integrator->getIntegratorTime() + dt,
                                                      1/*???*/);
        ib_method_ops->interpolateVelocity(u_ghost_idx, {}, {}, 0.0);

        {
            const std::string& velocity_name = IBFEMethod::VELOCITY_SYSTEM_NAME;
            const std::string& displacement_name = IBFEMethod::COORDS_SYSTEM_NAME;
            EquationSystems* equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();
            System& velocity_system = equation_systems->get_system(velocity_name);
            System& displacement_system = equation_systems->get_system(displacement_name);

            const unsigned int n_vars = velocity_system.n_vars();
            const DofMap& dof_map = velocity_system.get_dof_map();
            FEDataManager::SystemDofMapCache& X_dof_map_cache = *(ib_method_ops->getFEDataManager()
                ->getDofMapCache(displacement_name));
            FEDataManager::SystemDofMapCache& velocity_dof_map_cache = *(ib_method_ops->getFEDataManager()
                ->getDofMapCache(velocity_name));
            FEType fe_type = dof_map.variable_type(0);

            PetscVector<double> &current_velocity =
                *dynamic_cast<PetscVector<double> *>(velocity_system.current_local_solution.get());
            PetscVector<double> &current_X =
                *dynamic_cast<PetscVector<double> *>(displacement_system.current_local_solution.get());

            std::unique_ptr<FEBase> v_fe(FEBase::build(NDIM, fe_type));
            std::unique_ptr<FEBase> X_fe(FEBase::build(NDIM, fe_type));
            QGauss v_qrule(NDIM, FIFTH);
            QGauss X_qrule(NDIM, FIFTH);
            v_fe->attach_quadrature_rule(&v_qrule);
            X_fe->attach_quadrature_rule(&X_qrule);
            const std::vector<std::vector<Real>>& v_phi = v_fe->get_phi();
            const std::vector<std::vector<Real>>& X_phi = X_fe->get_phi();
            std::vector<dof_id_type> v_dof_indices;
            std::array<double, NDIM> max_norm_errors;
            std::fill(max_norm_errors.begin(), max_norm_errors.end(), 0.0);
            TBOX_ASSERT(n_vars == NDIM);

            std::vector<std::string> fs;
            {
                Pointer<Database> v_db = input_db->getDatabase("VelocityInitialConditions");
                for (unsigned int var_n = 0; var_n < n_vars; ++var_n)
                    fs.push_back(v_db->getString("function_" + std::to_string(var_n)));
            }
            ParsedFunction exact_solution(fs);

            std::vector<libMesh::Point> mapped_q_points;
            double max_vertex_distance = 0.0;
            for (auto elem_iter = mesh.active_local_elements_begin();
                 elem_iter != mesh.active_local_elements_end(); ++elem_iter)
            {
                auto elem = *elem_iter;
                max_vertex_distance = std::max(max_vertex_distance, elem->hmax());
                v_fe->reinit(elem);
                X_fe->reinit(elem);
                const auto& v_dof_indices = velocity_dof_map_cache.dof_indices(elem);
                const auto& X_dof_indices = X_dof_map_cache.dof_indices(elem);
                const std::size_t v_n_qp = v_qrule.n_points();
                const std::size_t X_n_qp = X_qrule.n_points();
                TBOX_ASSERT(v_n_qp == X_n_qp);

                // 1. Figure out where the quadrature points really are:
                mapped_q_points.resize(v_n_qp);
                std::fill(mapped_q_points.begin(), mapped_q_points.end(), libMesh::Point());
                for (unsigned int var_n = 0; var_n < n_vars; ++var_n)
                {
                    const std::size_t n_basis = X_dof_indices[var_n].size();
                    for (unsigned int q_point_n = 0; q_point_n < X_n_qp; ++q_point_n)
                    {
                        for (unsigned int i = 0; i < n_basis; ++i)
                        {
                            mapped_q_points[q_point_n](var_n) += current_X(X_dof_indices[var_n][i])
                                * X_phi[i][q_point_n];
                        }
                    }
                }

                // 2. Compute the difference:
                for (unsigned int var_n = 0; var_n < n_vars; ++var_n)
                {
                    const std::size_t n_basis = v_dof_indices[var_n].size();
                    for (unsigned int q_point_n = 0; q_point_n < v_n_qp; ++q_point_n)
                    {
                        const double exact_point_value = exact_solution.value(mapped_q_points[q_point_n], var_n);
                        double fe_point_value = 0.0;
                        for (unsigned int i = 0; i < n_basis; ++i)
                        {
                            fe_point_value += current_velocity(v_dof_indices[var_n][i]) * v_phi[i][q_point_n];
                        }
                        max_norm_errors[var_n] = std::max(max_norm_errors[var_n],
                                                          std::abs(fe_point_value - exact_point_value));
                    }
                }
            }

            plog << "max vertex distance: " << SAMRAI_MPI::maxReduction(max_vertex_distance) << std::endl;
            plog << "max norm errors: ";
            for (unsigned int i = 0; i < n_vars - 1; ++i)
            {
                plog << std::setprecision(20) << max_norm_errors[i] << "   ";
            }
            plog << max_norm_errors[n_vars - 1] << std::endl;
        }
    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
} // main
