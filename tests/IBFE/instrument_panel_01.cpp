// ---------------------------------------------------------------------
//
// Copyright (c) 2022 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <HierarchyDataOpsManager.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/numeric_vector.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEInstrumentPanel.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include "ibtk/LEInteractor.h"
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/StableCentroidPartitioner.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Test IBFEInstrumentPanel by computing fluxes

// A scalar function parser that behaves similarly to our own
// muParserCartGridFunction.
class ParsedFunction
{
public:
    ParsedFunction(std::vector<std::string> expressions, const unsigned int n_vars = 0)
        : string_functions(std::move(expressions)),
          vars(n_vars == 0 ? NDIM : n_vars, 0.0),
          parsers(string_functions.size())
    {
        for (mu::Parser& parser : parsers)
        {
            parser.DefineConst("PI", M_PI);
            for (unsigned int var_n = 0; var_n < NDIM; ++var_n)
            {
                parser.DefineVar("X_" + std::to_string(var_n), &vars[var_n]);
            }
        }

        for (unsigned int p_n = 0; p_n < string_functions.size(); ++p_n) parsers[p_n].SetExpr(string_functions[p_n]);
    }

    // Evaluate a single component.
    double value(const libMesh::Point& p, const unsigned int component_n) const
    {
        TBOX_ASSERT(component_n < parsers.size());
        for (unsigned int var_n = 0; var_n < vars.size(); ++var_n) vars[var_n] = p(var_n);

        try
        {
            return parsers[component_n].Eval();
        }
        catch (mu::ParserError& e)
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
    libMesh::Point value(const libMesh::Point& p) const
    {
        TBOX_ASSERT(vars.size() == parsers.size());
        for (unsigned int var_n = 0; var_n < vars.size(); ++var_n) vars[var_n] = p(var_n);

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

int
main(int argc, char** argv)
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    // set up options for the linear solver to not depend on the parallel partitioning
    PetscOptionsSetValue(nullptr, "-ksp_rtol", "1e-14");
    PetscOptionsSetValue(nullptr, "-ksp_atol", "1e-12");
    PetscOptionsSetValue(nullptr, "-ksp_type", "cg");
    PetscOptionsSetValue(nullptr, "-pc_type", "jacobi");
    PetscOptionsSetValue(nullptr, "-pc_jacobi_type", "diagonal");

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);

    // suppress warnings about using QGRID with too high of an order
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);

    {
        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");

        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create a simple FE mesh.
        ReplicatedMesh mesh(init.comm(), NDIM);
        MeshTools::Generation::build_cube(mesh, 4, 4, 4, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, libMesh::HEX8);
        BoundaryInfo& boundary_info = mesh.get_boundary_info();
        boundary_info.clear();

        const auto end_node_it = mesh.nodes_end();
        const bool skip_some_nodes = input_db->getBoolWithDefault("skip_some_nodes", false);
        if (skip_some_nodes)
        {
            for (auto node_it = mesh.nodes_begin(); node_it != end_node_it; ++node_it)
            {
                Node& node = **node_it;
                // outline the bottom face with id 10: skip some nodes on the top side
                if (node(2) == 0.0 && (node(0) == 0.0 || node(0) == 1.0 || node(1) == 0.0 || node(1) == 1.0))
                {
                    if (node(1) == 1.0)
                    {
                        if (node(0) == 0.0 || node(0) == 1.0) boundary_info.add_node(&node, 10);
                    }
                    else
                    {
                        boundary_info.add_node(&node, 10);
                    }
                }

                // outline the top face with id 100: skip some nodes on the left side
                if (node(2) == 1.0 && (node(0) == 0.0 || node(0) == 1.0 || node(1) == 0.0 || node(1) == 1.0))
                {
                    if (node(0) == 0.0)
                    {
                        if (node(1) == 0.0 || node(1) == 1.0) boundary_info.add_node(&node, 100);
                    }
                    else
                    {
                        boundary_info.add_node(&node, 100);
                    }
                }
            }
        }
        else
        {
            for (auto node_it = mesh.nodes_begin(); node_it != end_node_it; ++node_it)
            {
                Node& node = **node_it;
                // outline the bottom face with id 10
                if (node(2) == 0.0 && (node(0) == 0.0 || node(0) == 1.0 || node(1) == 0.0 || node(1) == 1.0))
                    boundary_info.add_node(&node, 10);

                // outline the top face with id 100
                if (node(2) == 1.0 && (node(0) == 0.0 || node(0) == 1.0 || node(1) == 0.0 || node(1) == 1.0))
                    boundary_info.add_node(&node, 100);
            }
        }

        plog << "Number of elements: " << mesh.n_active_elem() << std::endl;

        // Create major algorithm and data objects that comprise the
        // application. These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"), false);
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy =
            new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry, false);
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();

        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const std::string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"),
                false);
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"),
                false);
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
        ib_method_ops->initializeFEEquationSystems();

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

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
                std::ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const std::string bc_coefs_name = bc_coefs_name_stream.str();

                std::ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const std::string bc_coefs_db_name = bc_coefs_db_name_stream.str();

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

        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        time_integrator->registerVisItDataWriter(visit_data_writer);

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Configure the IBFE solver.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Force the velocity and pressure to agree with the analytic solutions.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            Pointer<hier::Variable<NDIM> > u_var = navier_stokes_integrator->getVelocityVariable();
            Pointer<VariableContext> u_current_ctx = navier_stokes_integrator->getCurrentContext();
            const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, u_current_ctx);
            u_init->setDataOnPatchHierarchy(u_current_idx, u_var, patch_hierarchy, 0.0);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            Pointer<hier::Variable<NDIM> > p_var = navier_stokes_integrator->getPressureVariable();
            Pointer<VariableContext> p_current_ctx = navier_stokes_integrator->getCurrentContext();
            const int p_current_idx = var_db->mapVariableAndContextToIndex(p_var, p_current_ctx);
            p_init->setDataOnPatchHierarchy(p_current_idx, p_var, patch_hierarchy, 0.0);
        }

        //********************************************************************************
        // setting up some objects for measuring fluxes and mean pressures
        //********************************************************************************
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<hier::Variable<NDIM> > p_var = navier_stokes_integrator->getPressureVariable();
        Pointer<VariableContext> p_current_ctx = navier_stokes_integrator->getCurrentContext();
        Pointer<hier::Variable<NDIM> > u_var = navier_stokes_integrator->getVelocityVariable();
        Pointer<VariableContext> u_current_ctx = navier_stokes_integrator->getCurrentContext();
        const int p_current_idx = var_db->mapVariableAndContextToIndex(p_var, p_current_ctx);
        const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, u_current_ctx);
        Pointer<SideVariable<NDIM, double> > u_copy_var = new SideVariable<NDIM, double>("u_copy");
        Pointer<CellVariable<NDIM, double> > p_copy_var = new CellVariable<NDIM, double>("p_copy");
        const IntVector<NDIM> ib_ghosts = ib_method_ops->getMinimumGhostCellWidth();
        const int u_copy_idx =
            var_db->registerVariableAndContext(u_copy_var, time_integrator->getScratchContext(), ib_ghosts);
        const int p_copy_idx =
            var_db->registerVariableAndContext(p_copy_var, time_integrator->getScratchContext(), ib_ghosts);
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_copy_idx, 0.0);
            level->allocatePatchData(p_copy_idx, 0.0);
        }

        //***********************************************
        // get mean pressure and velocity on surface mesh
        //***********************************************
        HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_cc_data_ops =
            hier_data_ops_manager->getOperationsDouble(p_var, patch_hierarchy, true);
        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_sc_data_ops =
            hier_data_ops_manager->getOperationsDouble(u_var, patch_hierarchy, true);
        hier_cc_data_ops->copyData(p_copy_idx, p_current_idx, true);
        hier_sc_data_ops->copyData(u_copy_idx, u_current_idx, true);

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
        p_hier_bdry_fill->fillData(0.0);

        // fill ghost cells for velocity
        std::vector<InterpolationTransactionComponent> u_transaction_comp(1);
        u_transaction_comp[0] = InterpolationTransactionComponent(u_copy_idx,
                                                                  "CONSERVATIVE_LINEAR_REFINE",
                                                                  /*use_cf_bdry_interpolation*/ false,
                                                                  "CONSERVATIVE_COARSEN",
                                                                  "LINEAR");

        Pointer<HierarchyGhostCellInterpolation> u_hier_bdry_fill = new HierarchyGhostCellInterpolation();
        u_hier_bdry_fill->initializeOperatorState(u_transaction_comp, patch_hierarchy);
        u_hier_bdry_fill->fillData(0.0);

        //*******************
        // Do the actual test
        //*******************
        IBFEInstrumentPanel instrument_panel(input_db->getDatabase("IBFEInstrumentPanel"), 0);

        instrument_panel.initializeHierarchyIndependentData(&*ib_method_ops);

        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);
        }

        auto do_instrument_panel = [&](const int data_time) {
            // reallocate in case we regridded
            HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
            hier_math_ops.setPatchHierarchy(patch_hierarchy);
            hier_math_ops.resetLevels(coarsest_ln, finest_ln);
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                if (!level->checkAllocated(p_copy_idx)) level->allocatePatchData(p_copy_idx);
                if (!level->checkAllocated(u_copy_idx)) level->allocatePatchData(u_copy_idx);
            }

            hier_cc_data_ops->copyData(p_copy_idx, p_current_idx, true);
            hier_sc_data_ops->copyData(u_copy_idx, u_current_idx, true);

            p_hier_bdry_fill->initializeOperatorState(p_transaction_comp, patch_hierarchy);
            p_hier_bdry_fill->fillData(double(data_time));

            u_hier_bdry_fill->initializeOperatorState(u_transaction_comp, patch_hierarchy);
            u_hier_bdry_fill->fillData(double(data_time));

            // update the displacement of IBFEMethod
            instrument_panel.readInstrumentData(u_copy_idx, p_copy_idx, patch_hierarchy, ib_method_ops, data_time);

            const int n_meters = instrument_panel.getNumberOfMeterMeshes();

            for (int meter_n = 0; meter_n < n_meters; ++meter_n)
            {
                ExodusII_IO meter_exodus_io(instrument_panel.getMeterMesh(meter_n));
                meter_exodus_io.write_timestep(instrument_panel.getPlotDirectoryName() + "/" +
                                                   instrument_panel.getMeterMeshName(meter_n) + ".e",
                                               instrument_panel.getMeterMeshEquationSystems(meter_n),
                                               data_time + 1,
                                               double(data_time));
            }

            // compute the fluxes through the faces with surface quadrature directly
            std::unique_ptr<QBase> qrule_face = QBase::build(QGAUSS, NDIM - 1, TENTH);
            std::unique_ptr<FEBase> fe = FEBase::build(NDIM, FEType());
            fe->attach_quadrature_rule(qrule_face.get());
            const std::vector<double>& JxW = fe->get_JxW();
            const std::vector<libMesh::Point>& q_points = fe->get_xyz();

            double bottom_analytic_flux = 0.0;
            double top_analytic_flux = 0.0;

            std::vector<std::string> velocity_strings;
            {
                Pointer<Database> v_db = input_db->getDatabase("VelocityInitialConditions");
                for (unsigned int var_n = 0; var_n < NDIM; ++var_n)
                    velocity_strings.push_back(v_db->getString("function_" + std::to_string(var_n)));
            }
            std::vector<std::string> pressure_strings;
            {
                Pointer<Database> v_db = input_db->getDatabase("PressureInitialConditions");
                pressure_strings.push_back(v_db->getString("function"));
            }
            ParsedFunction exact_velocity(velocity_strings);
            ParsedFunction exact_pressure(pressure_strings);

            {
                const auto end_node_it = mesh.nodes_end();
                for (auto node_it = mesh.nodes_begin(); node_it != end_node_it; ++node_it)
                {
                    Node& node = **node_it;
                    // second test: translation
                    if (data_time == 1)
                    {
                        node(0) += 1;
                        node(1) += 1;
                        node(2) += 1;
                    }
                    // third test: dilation
                    else if (data_time == 2)
                    {
                        node(0) *= 0.5625;
                        node(2) += 2;
                    }
                    // fourth test: translation + background velocity
                    else if (data_time == 3)
                    {
                        node(0) += 3;
                        node(1) += 3;
                        node(2) += 3;
                    }
                }
            }

            const auto end_el = mesh.active_local_elements_end();
            for (auto el_it = mesh.active_local_elements_begin(); el_it != end_el; ++el_it)
            {
                Elem& elem = **el_it;
                for (unsigned int side_n = 0; side_n < elem.n_sides(); ++side_n)
                {
                    // at data_time = t the bottom face is at z = t and the top face is at z = 1 + t
                    std::unique_ptr<Elem> side_ptr = elem.side_ptr(side_n);

#if LIBMESH_VERSION_LESS_THAN(1, 7, 0)
                    const double z = side_ptr->centroid()(2);
#else
                    const double z = side_ptr->vertex_average()(2);
#endif
                    // include a simple correction velocity for the 4th test only (note that it needs to be a linear
                    // field)
                    double u_corr = (data_time == 3 ? 1.0 : 0.0);
                    if (std::abs(z - (data_time + 0.0)) < 1.0e-5)
                    {
                        fe->reinit(&elem, side_n);
                        // hard-codes in the flux as the z velocity * +1 (to match sign convention in the meter)
                        for (std::size_t qp = 0; qp < JxW.size(); ++qp)
                        {
                            const auto x = q_points[qp];
                            bottom_analytic_flux += +1.0 * JxW[qp] * (exact_velocity.value(x, 2) - u_corr * 3 * x(2));
                        }
                    }
                    else if (std::abs(z - (data_time + 1.0)) < 1.0e-5)
                    {
                        fe->reinit(&elem, side_n);
                        // hard-codes in the flux as the z velocity * +1 (to match sign convention in the meter)
                        for (std::size_t qp = 0; qp < JxW.size(); ++qp)
                        {
                            const auto x = q_points[qp];
                            top_analytic_flux += +1.0 * JxW[qp] * (exact_velocity.value(x, 2) - u_corr * 3 * x(2));
                        }
                    }
                }
            }

            // and move them back
            {
                const auto end_node_it = mesh.nodes_end();
                for (auto node_it = mesh.nodes_begin(); node_it != end_node_it; ++node_it)
                {
                    Node& node = **node_it;
                    if (data_time == 1)
                    {
                        node(0) -= 1;
                        node(1) -= 1;
                        node(2) -= 1;
                    }
                    else if (data_time == 2)
                    {
                        node(0) /= 0.5625;
                        node(2) -= 2;
                    }
                }
            }

            if (input_db->getBoolWithDefault("test_flux", false))
            {
                bottom_analytic_flux = IBTK_MPI::sumReduction(bottom_analytic_flux);
                top_analytic_flux = IBTK_MPI::sumReduction(top_analytic_flux);

                plog << "bottom analytic flux = " << bottom_analytic_flux << std::endl;
                plog << "top analytic flux = " << top_analytic_flux << std::endl;

                const double bottom_meter_flux = instrument_panel.getMeterFlowRates()[0];
                const double top_meter_flux = instrument_panel.getMeterFlowRates()[1];

                plog << "bottom meter flux = " << bottom_meter_flux << std::endl;
                plog << "top meter flux = " << top_meter_flux << std::endl;

                plog << "bottom flux absolute difference = " << std::abs(bottom_analytic_flux - bottom_meter_flux)
                     << std::endl;
                plog << "top flux absolute difference = " << std::abs(top_analytic_flux - top_meter_flux) << std::endl;
            }

            if (input_db->getBoolWithDefault("test_centroid_pressure", false))
            {
                for (int meter_n = 0; meter_n < instrument_panel.getNumberOfMeterMeshes(); ++meter_n)
                {
                    const MeshBase& mesh = instrument_panel.getMeterMesh(meter_n);
                    // We assume that the centroid is the last node
                    const Node& centroid = mesh.node_ref(mesh.n_nodes() - 1);
                    plog << "mesh " << meter_n << "  centroid = " << centroid << std::endl
                         << "meter " << meter_n
                         << " centroid pressure = " << instrument_panel.getMeterCentroidPressures()[meter_n]
                         << std::endl
                         << "analytic centroid " << meter_n << " pressure = " << exact_pressure.value(centroid, 0)
                         << std::endl;
                }
            }
        };

        // first test: print out values as-is
        {
            plog << "\ntest 0\n";
            do_instrument_panel(0);
        }

        // second test: translate the mesh
        {
            plog << "\ntest 1\n";
            FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager();
            EquationSystems* equation_systems = fe_data_manager->getEquationSystems();
            System& dX_system = equation_systems->get_system(ib_method_ops->getDisplacementSystemName());
            NumericVector<double>& dX_vec = *dX_system.solution;
            dX_vec.add(1.0);
            System& X_system = equation_systems->get_system(ib_method_ops->getCurrentCoordinatesSystemName());
            NumericVector<double>& X_vec = *X_system.solution;
            X_vec += dX_vec;
            do_instrument_panel(1);
            X_vec -= dX_vec;
            dX_vec = 0.0;
        }

        // third test: dilate the mesh
        {
            plog << "\ntest 2\n";
            FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager();
            EquationSystems* equation_systems = fe_data_manager->getEquationSystems();
            System& dX_system = equation_systems->get_system(ib_method_ops->getDisplacementSystemName());
            NumericVector<double>& dX_vec = *dX_system.solution;

            std::vector<dof_id_type> component_dofs;
            for (auto node_it = mesh.nodes_begin(); node_it != end_node_it; ++node_it)
            {
                Node& node = **node_it;

                // shrink x
                IBTK::get_nodal_dof_indices(dX_system.get_dof_map(), &node, 0, component_dofs);
                TBOX_ASSERT(component_dofs.size() == 1);
                const double current_x = node(0);
                const double new_x = node(0) * 0.5625;
                dX_vec.set(component_dofs[0], new_x - current_x);

                // translate z
                IBTK::get_nodal_dof_indices(dX_system.get_dof_map(), &node, 2, component_dofs);
                TBOX_ASSERT(component_dofs.size() == 1);
                dX_vec.set(component_dofs[0], 2.0);
            }
            dX_vec.close();

            System& X_system = equation_systems->get_system(ib_method_ops->getCurrentCoordinatesSystemName());
            NumericVector<double>& X_vec = *X_system.solution;
            X_vec += dX_vec;
            do_instrument_panel(2);
            X_vec -= dX_vec;
            dX_vec = 0.0;
        }

        // fourth test: translate and add a background velocity
        {
            plog << "\ntest 3\n";
            FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager();
            EquationSystems* equation_systems = fe_data_manager->getEquationSystems();
            System& dX_system = equation_systems->get_system(ib_method_ops->getDisplacementSystemName());
            NumericVector<double>& dX_vec = *dX_system.solution;
            dX_vec = 3;
            System& X_system = equation_systems->get_system(ib_method_ops->getCurrentCoordinatesSystemName());
            NumericVector<double>& X_vec = *X_system.solution;
            X_vec += dX_vec;
            System& U_system = equation_systems->get_system(ib_method_ops->getVelocitySystemName());
            NumericVector<double>& U_vec = *U_system.solution;
            std::vector<dof_id_type> component_dofs;
            for (auto node_it = mesh.nodes_begin(); node_it != end_node_it; ++node_it)
            {
                Node& node = **node_it;
                for (int d = 0; d < NDIM; ++d)
                {
                    IBTK::get_nodal_dof_indices(U_system.get_dof_map(), &node, d, component_dofs);
                    TBOX_ASSERT(component_dofs.size() == 1);
                    // NOTE: velocity needs to be evaluated using the shifted coordinates!
                    U_vec.set(component_dofs[0], (d + 1) * (node(d) + 3));
                }
            }
            U_vec.close();
            do_instrument_panel(3);
            X_vec -= dX_vec;
            dX_vec = 0.0;
            U_vec = 0.0;
        }

    } // cleanup dynamically allocated objects prior to shutdown
} // main
