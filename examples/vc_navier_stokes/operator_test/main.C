// Copyright (c) 2002-2010, Boyce Griffith, Thomas Fai
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
#include <IBTK_config.h>
#include <SAMRAI_config.h>
#include <petsc.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <GriddingAlgorithm.h>
#include <HierarchyDataOpsManager.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>
#include <VisItDataWriter.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibamr/INSStaggeredVCStokesOperator.h>

using namespace SAMRAI;
using namespace IBTK;
using namespace IBAMR;
using namespace std;

/************************************************************************
 *                                                                      *
 * For each run, the input filename must be given on the command line.  *
 * In all cases, the command line is:                                   *
 *                                                                      *
 *    executable <input file name> <PETSc options>                      *
 *                                                                      *
 ************************************************************************
 */

int
main(
    int argc,
    char *argv[])
{
    /*
     * Initialize PETSc, MPI, and SAMRAI.
     */
    PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
    tbox::SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    tbox::SAMRAIManager::startup();

    {// ensure all smart Pointers are properly deleted

        string input_filename = argv[1];
        tbox::plog << "input_filename = " << input_filename << "\n";

        /*
         * Create input database and parse all data in input file.
         */
        tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
        tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

        /*
         * Retrieve "Main" section of the input database.
         */
        tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

        string log_file_name = "laplace_test.log";
        if (main_db->keyExists("log_file_name"))
        {
            log_file_name = main_db->getString("log_file_name");
        }
        bool log_all_nodes = false;
        if (main_db->keyExists("log_all_nodes"))
        {
            log_all_nodes = main_db->getBool("log_all_nodes");
        }
        if (log_all_nodes)
        {
            tbox::PIO::logAllNodes(log_file_name);
        }
        else
        {
            tbox::PIO::logOnlyNodeZero(log_file_name);
        }

        string visit_dump_dirname = "viz";
        if (main_db->keyExists("visit_dump_dirname"))
        {
            visit_dump_dirname = main_db->getString("visit_dump_dirname");
        }

        int visit_number_procs_per_file = 1;
        if (main_db->keyExists("visit_number_procs_per_file"))
        {
            visit_number_procs_per_file = main_db->getInteger("visit_number_procs_per_file");
        }

        bool timer_enabled = false;
        if (main_db->keyExists("timer_enabled"))
        {
            timer_enabled = main_db->getBool("timer_enabled");
        }
        if (timer_enabled)
        {
            tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
        }

        /*
         * Create major algorithm and data objects which comprise the
         * application.  Each object will be initialized from input data from
         * the input database.  Refer to each class constructor for details.
         */
        tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry =
            new geom::CartesianGridGeometry<NDIM>("CartesianGeometry",
                                                  input_db->getDatabase("CartesianGeometry"));

        tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
            new hier::PatchHierarchy<NDIM>("PatchHierarchy",
                                           grid_geometry);

        tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector =
            new mesh::StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                                     NULL,
                                                     input_db->getDatabase("StandardTagAndInitialize"));

        tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator =
            new mesh::BergerRigoutsos<NDIM>();

        tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer =
            new mesh::LoadBalancer<NDIM>("LoadBalancer",
                                         input_db->getDatabase("LoadBalancer"));

        tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > gridding_algorithm =
            new mesh::GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                              input_db->getDatabase("GriddingAlgorithm"),
                                              error_detector,
                                              box_generator,
                                              load_balancer);

        /*
         * Create variables and register then with the variable database.
         */
        hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
        tbox::Pointer<hier::VariableContext> ctx = var_db->getContext("context");

        tbox::Pointer<pdat::SideVariable<NDIM,double> > u_side_var = new pdat::SideVariable<NDIM,double>("u_side");
        tbox::Pointer<pdat::SideVariable<NDIM,double> > f_side_var = new pdat::SideVariable<NDIM,double>("f_side");
        tbox::Pointer<pdat::SideVariable<NDIM,double> > e_side_var = new pdat::SideVariable<NDIM,double>("e_side");

        const int u_side_idx = var_db->registerVariableAndContext(u_side_var, ctx, hier::IntVector<NDIM>((USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1)));
        const int f_side_idx = var_db->registerVariableAndContext(f_side_var, ctx, hier::IntVector<NDIM>((USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1)));
        const int e_side_idx = var_db->registerVariableAndContext(e_side_var, ctx, hier::IntVector<NDIM>((USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1)));

        tbox::Pointer<pdat::CellVariable<NDIM,double> > u_cell_var = new pdat::CellVariable<NDIM,double>("u_cell",NDIM);
        tbox::Pointer<pdat::CellVariable<NDIM,double> > f_cell_var = new pdat::CellVariable<NDIM,double>("f_cell",NDIM);
        tbox::Pointer<pdat::CellVariable<NDIM,double> > e_cell_var = new pdat::CellVariable<NDIM,double>("e_cell",NDIM);
        tbox::Pointer<pdat::CellVariable<NDIM,double> > p_cell_var = new pdat::CellVariable<NDIM,double>("p_cell");
        tbox::Pointer<pdat::CellVariable<NDIM,double> > du_cell_var = new pdat::CellVariable<NDIM,double>("du_cell");

        const int u_cell_idx = var_db->registerVariableAndContext(u_cell_var, ctx, hier::IntVector<NDIM>(0));
        const int f_cell_idx = var_db->registerVariableAndContext(f_cell_var, ctx, hier::IntVector<NDIM>(0));
        const int e_cell_idx = var_db->registerVariableAndContext(e_cell_var, ctx, hier::IntVector<NDIM>(0));
        const int p_cell_idx = var_db->registerVariableAndContext(p_cell_var, ctx, hier::IntVector<NDIM>((USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1)));
        const int du_cell_idx = var_db->registerVariableAndContext(du_cell_var, ctx, hier::IntVector<NDIM>(0));

        tbox::Pointer<pdat::NodeVariable<NDIM,double> > mu_node_var = new pdat::NodeVariable<NDIM,double>("mu_node");
        tbox::Pointer<pdat::SideVariable<NDIM,double> > rho_side_var = new pdat::SideVariable<NDIM,double>("rho_side");

        const int mu_node_idx = var_db->registerVariableAndContext(mu_node_var, ctx, hier::IntVector<NDIM>((USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1)));
        const int rho_side_idx = var_db->registerVariableAndContext(rho_side_var, ctx, hier::IntVector<NDIM>((USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1)));

        /*
         * Create visualization plot file writer and register variables for
         * plotting.
         */
        tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer =
            new appu::VisItDataWriter<NDIM>("VisIt Writer",
                                            visit_dump_dirname,
                                            visit_number_procs_per_file);

        visit_data_writer->registerPlotQuantity(u_cell_var->getName(), "VECTOR", u_cell_idx);
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream stream;
            stream << d;
            visit_data_writer->registerPlotQuantity(u_cell_var->getName()+stream.str(), "SCALAR", u_cell_idx, d);
        }

        visit_data_writer->registerPlotQuantity(f_cell_var->getName(), "VECTOR", f_cell_idx);
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream stream;
            stream << d;
            visit_data_writer->registerPlotQuantity(f_cell_var->getName()+stream.str(), "SCALAR", f_cell_idx, d);
        }

        visit_data_writer->registerPlotQuantity(e_cell_var->getName(), "VECTOR", e_cell_idx);
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream stream;
            stream << d;
            visit_data_writer->registerPlotQuantity(e_cell_var->getName()+stream.str(), "SCALAR", e_cell_idx, d);
        }

        visit_data_writer->registerPlotQuantity(mu_node_var->getName(), "SCALAR", mu_node_idx);

        //visit_data_writer->registerPlotQuantity(rho_side_var->getName(), "SCALAR", rho_side_idx);

        //visit_data_writer->registerPlotQuantity(p_cell_var->getName(), "SCALAR", p_cell_idx);

        /*
         * Initialize the AMR patch hierarchy.
         */
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);

        int tag_buffer = 1;
        int level_number = 0;
        bool done = false;
        while (!done && (gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        /*
         * Set the simulation time to be zero.
         */
        const double data_time = 0.0;

        /*
         * Allocate data on each level of the patch hierarchy.
         */
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            tbox::Pointer<hier::PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData( u_side_idx, data_time);
            level->allocatePatchData( f_side_idx, data_time);
            level->allocatePatchData( e_side_idx, data_time);
            level->allocatePatchData( u_cell_idx, data_time);
            level->allocatePatchData( f_cell_idx, data_time);
            level->allocatePatchData( e_cell_idx, data_time);
            level->allocatePatchData(mu_node_idx, data_time);
            level->allocatePatchData(rho_side_idx, data_time);
            level->allocatePatchData(p_cell_idx, data_time);
            level->allocatePatchData(du_cell_idx, data_time);
	}

        /*
         * Setup exact solution data.
         */
        muParserCartGridFunction u_fcn("u", input_db->getDatabase("u"), grid_geometry);
        muParserCartGridFunction f_fcn("f", input_db->getDatabase("f"), grid_geometry);
        muParserCartGridFunction mu_fcn("mu", input_db->getDatabase("mu"), grid_geometry);
        muParserCartGridFunction rho_fcn("rho", input_db->getDatabase("rho"), grid_geometry);
        muParserCartGridFunction p_fcn("p", input_db->getDatabase("p"), grid_geometry);
        muParserCartGridFunction du_fcn("du", input_db->getDatabase("du"), grid_geometry);

        u_fcn.setDataOnPatchHierarchy(u_side_idx, u_side_var, patch_hierarchy, data_time);
        f_fcn.setDataOnPatchHierarchy(e_side_idx, e_side_var, patch_hierarchy, data_time);
        mu_fcn.setDataOnPatchHierarchy(mu_node_idx, mu_node_var, patch_hierarchy, data_time);
        rho_fcn.setDataOnPatchHierarchy(rho_side_idx, rho_side_var, patch_hierarchy, data_time);
        p_fcn.setDataOnPatchHierarchy(p_cell_idx, p_cell_var, patch_hierarchy, data_time);
        du_fcn.setDataOnPatchHierarchy(du_cell_idx, du_cell_var, patch_hierarchy, data_time);

        /*
         * Create the math operations object and get the patch data index for
         * the side-centered weighting factor.
         */
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int dx_side_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
        const int dx_cell_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        /*
         * Compute f = [(rho/dt)*u-0.5*div*(mu*(grad u + (grad u)^T)) + grad p; -div u]
         */
        INSStaggeredVCStokesOperator vc_stokes_op(tbox::Pointer<HierarchyMathOps>(&hier_math_ops,false));  // Create a Pointer to hier_math_ops which does NOT handle memory management for hier_math_ops (i.e., which does NOT delete hier_math_ops when the number of references drops to zero).

	solv::SAMRAIVectorReal<NDIM,double> x("x", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
	solv::SAMRAIVectorReal<NDIM,double> y("y", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        x.addComponent(u_side_var,u_side_idx,dx_side_idx);
	x.addComponent(p_cell_var,p_cell_idx,dx_cell_idx);
	y.addComponent(f_side_var,f_side_idx,dx_side_idx);
        y.addComponent(du_cell_var,du_cell_idx,dx_cell_idx);

        vc_stokes_op.setTimeInterval(0.1,0.2);
        vc_stokes_op.registerViscosityVariable(mu_node_var,mu_node_idx);
        vc_stokes_op.registerDensityVariable(rho_side_var,rho_side_idx);
        vc_stokes_op.initializeOperatorState(x,y);
        vc_stokes_op.apply(x,y);

        /*
         * Compute error and print error norms.
         */
        tbox::Pointer<math::HierarchyDataOpsReal<NDIM,double> > hier_side_data_ops =
            math::HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(
                u_side_var, patch_hierarchy, true);
        hier_side_data_ops->subtract(e_side_idx, e_side_idx, f_side_idx);  // computes e := e - f

        tbox::pout << "|e|_oo = " << hier_side_data_ops->maxNorm(e_side_idx, dx_side_idx) << "\n";
        tbox::pout << "|e|_2  = " << hier_side_data_ops-> L2Norm(e_side_idx, dx_side_idx) << "\n";
        tbox::pout << "|e|_1  = " << hier_side_data_ops-> L1Norm(e_side_idx, dx_side_idx) << "\n";

        /*
         * Interpolate the side-centered data to cell centers (for output
         * purposes only).
         */
        static const bool synch_cf_interface = true;
        hier_math_ops.interp(u_cell_idx, u_cell_var, u_side_idx, u_side_var, NULL, data_time, synch_cf_interface);
        hier_math_ops.interp(f_cell_idx, f_cell_var, f_side_idx, f_side_var, NULL, data_time, synch_cf_interface);
        hier_math_ops.interp(e_cell_idx, e_cell_var, e_side_idx, e_side_var, NULL, data_time, synch_cf_interface);

        /*
         * Output data for plotting.
         */
        visit_data_writer->writePlotData(patch_hierarchy, 0, data_time);

    }// ensure all smart Pointers are properly deleted

    tbox::SAMRAIManager::shutdown();
    PetscFinalize();

    return 0;
}// main
