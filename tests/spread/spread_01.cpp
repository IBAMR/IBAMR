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

// This file is the main driver for force spreading tests (i.e.,
// IBFEmethod::spreadForce). At the moment it simply prints out the force
// values.

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

        // Now for the actual test. Set up a new variable containing ghost data:
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const Pointer<SAMRAI::hier::Variable<NDIM> > f_var = time_integrator->getBodyForceVariable();
        const Pointer<VariableContext> f_ghost_ctx = var_db->getContext("f_ghost");
        const int f_ghost_idx = var_db->registerVariableAndContext(f_var, f_ghost_ctx, 3);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(f_ghost_idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_ghost_idx);
                f_data->fillAll(0.0);
            }
        }

        // synch ghost data
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ghost_cell_components(1);
        ghost_cell_components[0] = InterpolationTransactionComponent(f_ghost_idx,
                                                                     "CONSERVATIVE_LINEAR_REFINE",
                                                                     true,
                                                                     "CONSERVATIVE_COARSEN",
                                                                     "LINEAR",
                                                                     false,
                                                                     {}, // f_bc_coefs
                                                                     nullptr);
        HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ghost_cell_components, patch_hierarchy);
        ghost_fill_op.fillData(/*time*/0.0);

        const double dt = time_integrator->getMaximumTimeStepSize();
        time_integrator->preprocessIntegrateHierarchy(time_integrator->getIntegratorTime(),
                                                      time_integrator->getIntegratorTime() + dt,
                                                      1/*???*/);
        auto &fe_data_manager = *ib_method_ops->getFEDataManager();
        auto &equation_systems = *fe_data_manager.getEquationSystems();
        auto &force_system = equation_systems.get_system(IBFEMethod::FORCE_SYSTEM_NAME);
        auto &half_f_vector = dynamic_cast<libMesh::PetscVector<double> &>(*force_system.current_local_solution);
        for (unsigned int i = half_f_vector.first_local_index(); i < half_f_vector.last_local_index(); ++i)
        {
            half_f_vector.set(i, i % 10);
        }
        half_f_vector.close();

        ib_method_ops->spreadForce(f_ghost_idx, nullptr, {}, time_integrator->getIntegratorTime() + dt/2);
        {
            const int ln = patch_hierarchy->getFinestLevelNumber();
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                plog << "patch number " << p() << '\n';
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_ghost_idx);
                const Box<NDIM> patch_box = patch->getBox();

                // same as SideData::print, but elides zero values
                plog.precision(16);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    plog << "Array side normal = " << axis << std::endl;
                    for (int d = 0; d < f_data->getDepth(); ++d)
                    {
                        plog << "Array depth = " << d << std::endl;
                        const ArrayData<NDIM, double> &data = f_data->getArrayData(axis);
                        for (SideIterator<NDIM> i(patch_box, axis); i; i++)
                        {
                            const double value = data(i(), d);
                            if (value != 0.0)
                                plog << "array" << i() << " = "
                                     << value << '\n';
                        }
                    }
                }
                // f_data->print(patch_box, plog, 16);
            }
        }
    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
} // main
