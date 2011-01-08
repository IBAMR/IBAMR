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
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc objects
#include <petsc.h>

// Headers for basic SAMRAI objects
#include <PatchLevel.h>
#include <VariableDatabase.h>
#include <tbox/Database.h>
#include <tbox/InputDatabase.h>
#include <tbox/InputManager.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>
#include <tbox/RestartManager.h>
#include <tbox/SAMRAIManager.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// Headers for major algorithm/data structure objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <StandardTagAndInitialize.h>
#include <VisItDataWriter.h>

// Headers for application-specific algorithm/data structure objects
#include <CellVariable.h>
#include <HierarchyDataOpsManager.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <ibtk/NormOps.h>
#include <ibtk/SCLaplaceOperator.h>
#include <ibtk/muParserCartGridFunction.h>

using namespace SAMRAI;
using namespace IBTK;
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
        string input_filename;
        input_filename = argv[1];

        tbox::plog << "input_filename = " << input_filename << endl;

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

        bool viz_dump_enabled = false;
        tbox::Array<string> viz_writer(1);
        viz_writer[0] = "VisIt";
        string viz_dump_filename;
        string visit_dump_dirname;
        bool uses_visit = false;
        int visit_number_procs_per_file = 1;
        if (main_db->keyExists("viz_dump_enabled"))
        {
            viz_dump_enabled = main_db->getBool("viz_dump_enabled");
        }
        if (viz_dump_enabled)
        {
            if (main_db->keyExists("viz_writer"))
            {
                viz_writer = main_db->getStringArray("viz_writer");
            }
            if (main_db->keyExists("viz_dump_filename"))
            {
                viz_dump_filename = main_db->getString("viz_dump_filename");
            }
            string viz_dump_dirname;
            if (main_db->keyExists("viz_dump_dirname"))
            {
                viz_dump_dirname = main_db->getString("viz_dump_dirname");
            }
            for (int i = 0; i < viz_writer.getSize(); ++i)
            {
                if (viz_writer[i] == "VisIt") uses_visit = true;
            }
            if (uses_visit)
            {
                visit_dump_dirname = viz_dump_dirname;
            }
            else
            {
                TBOX_ERROR("main(): "
                           << "\nUnrecognized 'viz_writer' entry..."
                           << "\nOnly valid option is 'VisIt'"
                           << endl);
            }
            if (uses_visit)
            {
                if (viz_dump_dirname.empty())
                {
                    TBOX_ERROR("main(): "
                               << "\nviz_dump_dirname is null ... "
                               << "\nThis must be specified for use with VisIt"
                               << endl);
                }
                if (main_db->keyExists("visit_number_procs_per_file"))
                {
                    visit_number_procs_per_file =
                        main_db->getInteger("visit_number_procs_per_file");
                }
            }
        }

        const bool viz_dump_data = viz_dump_enabled;

        bool timer_enabled = false;
        if (main_db->keyExists("timer_enabled"))
        {
            timer_enabled = main_db->getBool("timer_enabled");
        }

        const bool write_timer_data = timer_enabled;

        if (write_timer_data)
        {
            tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
        }

        /*
         * Create major algorithm and data objects which comprise application.
         * Each object will be initialized either from input data or restart
         * files, or a combination of both.  Refer to each class constructor for
         * details.  For more information on the composition of objects for this
         * application, see comments at top of file.
         */
        tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry =
            new geom::CartesianGridGeometry<NDIM>("CartesianGeometry",
                                                  input_db->getDatabase("CartesianGeometry"));

        tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
            new hier::PatchHierarchy<NDIM>("PatchHierarchy",grid_geometry);

        tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector =
            new mesh::StandardTagAndInitialize<NDIM>(
                "StandardTagAndInitialize",
                NULL,
                input_db->getDatabase("StandardTagAndInitialize"));

        tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator = new mesh::BergerRigoutsos<NDIM>();

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
         * Set up variables.
         */
        tbox::Pointer<pdat::SideVariable<NDIM,double> > u_sc_var = new pdat::SideVariable<NDIM,double>("u_sc");
        tbox::Pointer<pdat::SideVariable<NDIM,double> > f_sc_var = new pdat::SideVariable<NDIM,double>("f_sc");
        tbox::Pointer<pdat::SideVariable<NDIM,double> > e_sc_var = new pdat::SideVariable<NDIM,double>("e_sc");

        tbox::Pointer<pdat::CellVariable<NDIM,double> > u_cc_var = new pdat::CellVariable<NDIM,double>("u_cc",NDIM);
        tbox::Pointer<pdat::CellVariable<NDIM,double> > f_cc_var = new pdat::CellVariable<NDIM,double>("f_cc",NDIM);
        tbox::Pointer<pdat::CellVariable<NDIM,double> > e_cc_var = new pdat::CellVariable<NDIM,double>("e_cc",NDIM);

        hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
        tbox::Pointer<hier::VariableContext> ctx = var_db->getContext("context");
        const int u_sc_idx = var_db->registerVariableAndContext(u_sc_var, ctx, hier::IntVector<NDIM>((USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1)));
        const int f_sc_idx = var_db->registerVariableAndContext(f_sc_var, ctx, hier::IntVector<NDIM>((USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1)));
        const int e_sc_idx = var_db->registerVariableAndContext(e_sc_var, ctx, hier::IntVector<NDIM>((USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1)));

        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, hier::IntVector<NDIM>(0));
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, hier::IntVector<NDIM>(0));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, hier::IntVector<NDIM>(0));

        /*
         * Set up visualization plot file writer.
         */
        tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer;

        if (uses_visit)
        {
            visit_data_writer = new appu::VisItDataWriter<NDIM>("VisIt Writer",
                                                                visit_dump_dirname,
                                                                visit_number_procs_per_file);

            visit_data_writer->registerPlotQuantity(u_cc_var->getName(), "VECTOR", u_cc_idx);
            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                visit_data_writer->registerPlotQuantity(u_cc_var->getName()+stream.str(), "SCALAR", u_cc_idx, d);
            }

            visit_data_writer->registerPlotQuantity(f_cc_var->getName(), "VECTOR", f_cc_idx);
            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                visit_data_writer->registerPlotQuantity(f_cc_var->getName()+stream.str(), "SCALAR", f_cc_idx, d);
            }

            visit_data_writer->registerPlotQuantity(e_cc_var->getName(), "VECTOR", e_cc_idx);
            for (int d = 0; d < NDIM; ++d)
            {
                std::ostringstream stream;
                stream << d;
                visit_data_writer->registerPlotQuantity(e_cc_var->getName()+stream.str(), "SCALAR", e_cc_idx, d);
            }
        }

        /*
         * Initialize hierarchy configuration and data on all patches.
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

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            tbox::Pointer<hier::PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_sc_idx, 0.0);
            level->allocatePatchData(f_sc_idx, 0.0);
            level->allocatePatchData(e_sc_idx, 0.0);
            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
        }

        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int h_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        solv::SAMRAIVectorReal<NDIM,double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        solv::SAMRAIVectorReal<NDIM,double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        solv::SAMRAIVectorReal<NDIM,double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        u_vec.addComponent(u_sc_var, u_sc_idx, h_sc_idx);
        f_vec.addComponent(f_sc_var, f_sc_idx, h_sc_idx);
        e_vec.addComponent(e_sc_var, e_sc_idx, h_sc_idx);

        u_vec.setToScalar(0.0);
        f_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);

        /*
         * Setup exact solutions.
         */
        muParserCartGridFunction u_fcn("u", input_db->getDatabase("u"), grid_geometry);
        muParserCartGridFunction f_fcn("f", input_db->getDatabase("f"), grid_geometry);

        u_fcn.setDataOnPatchHierarchy(u_sc_idx, u_sc_var, patch_hierarchy, 0.0);
        f_fcn.setDataOnPatchHierarchy(e_sc_idx, e_sc_var, patch_hierarchy, 0.0);

        /*
         * Compute (I - L)*u = f.
         */
        solv::PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setCConstant( 0.0);
        poisson_spec.setDConstant(-1.0);
        vector<solv::RobinBcCoefStrategy<NDIM>*> bc_coefs(NDIM);
        SCLaplaceOperator laplace_op("laplace op", poisson_spec, bc_coefs);
        laplace_op.initializeOperatorState(u_vec,f_vec);
        laplace_op.apply(u_vec,f_vec);

        /*
         * Compute error, error norms.
         */
#if 0
        tbox::Pointer<math::HierarchyDataOpsReal<NDIM,double> > hier_sc_data_ops = math::HierarchyDataOpsManager<NDIM>::getManager()->
            getOperationsDouble(u_sc_var, patch_hierarchy, true);
        hier_sc_data_ops->resetLevels(0,0);
        hier_sc_data_ops->setToScalar(u_sc_idx, 0.0);
        hier_sc_data_ops->setToScalar(f_sc_idx, 0.0);
        hier_sc_data_ops->setToScalar(e_sc_idx, 0.0);
#endif
        e_vec.subtract(tbox::Pointer<solv::SAMRAIVectorReal<NDIM,double> >(&e_vec,false),
                       tbox::Pointer<solv::SAMRAIVectorReal<NDIM,double> >(&f_vec,false));
        tbox::pout << "|e|_oo = " << NormOps::maxNorm(&e_vec) << endl;
        tbox::pout << "|e|_2  = " << NormOps:: L2Norm(&e_vec) << endl;
        tbox::pout << "|e|_1  = " << NormOps:: L1Norm(&e_vec) << endl;

        /*
         * Write out data files for plotting.
         */
        static const bool synch_cf_interface = true;
        hier_math_ops.interp(u_cc_idx, u_cc_var, u_sc_idx, u_sc_var, NULL, 0.0, synch_cf_interface);
        hier_math_ops.interp(f_cc_idx, f_cc_var, f_sc_idx, f_sc_var, NULL, 0.0, synch_cf_interface);
        hier_math_ops.interp(e_cc_idx, e_cc_var, e_sc_idx, e_sc_var, NULL, 0.0, synch_cf_interface);
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber()-1; ++ln)
        {
            tbox::Pointer<hier::PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            hier::BoxArray<NDIM> refined_region_boxes;
            tbox::Pointer<hier::PatchLevel<NDIM> > next_finer_level = patch_hierarchy->getPatchLevel(ln+1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
            for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
                const hier::Box<NDIM>& patch_box = patch->getBox();
                tbox::Pointer<pdat::CellData<NDIM,double> > e_cc_data = patch->getPatchData(e_cc_idx);
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const hier::Box<NDIM> refined_box = refined_region_boxes[i];
                    const hier::Box<NDIM> intersection = hier::Box<NDIM>::grow(patch_box,1)*refined_box;
                    if (!intersection.empty())
                    {
                        e_cc_data->fillAll(0.0, intersection);
                    }
                }
            }
        }
        if (viz_dump_data)
        {
            if (uses_visit) visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);
        }
    }// ensure all smart Pointers are properly deleted

    tbox::SAMRAIManager::shutdown();
    PetscFinalize();

    return 0;
}// main
