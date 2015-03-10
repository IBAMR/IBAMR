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
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <ChopAndPackLoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/SCLaplaceOperator.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/app_namespaces.h>

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
int main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        auto app_initializer = boost::make_shared<AppInitializer>(argc, argv, "sc_laplace.log");
        auto input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        auto grid_geometry = boost::make_shared<CartesianGridGeometry>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        auto patch_hierarchy = boost::make_shared<PatchHierarchy>("PatchHierarchy", grid_geometry);
        auto error_detector = boost::make_shared<StandardTagAndInitialize>(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        auto box_generator = boost::make_shared<BergerRigoutsos>();
        auto load_balancer = boost::make_shared<ChopAndPackLoadBalancer>("ChopAndPackLoadBalancer", app_initializer->getComponentDatabase("ChopAndPackLoadBalancer"));
        auto gridding_algorithm = boost::make_shared<GriddingAlgorithm>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create variables and register them with the variable database.
        VariableDatabase* var_db = VariableDatabase::getDatabase();
        auto ctx = var_db->getContext("context");

        auto u_sc_var = boost::make_shared<SideVariable<double> >("u_sc");
        auto f_sc_var = boost::make_shared<SideVariable<double> >("f_sc");
        auto e_sc_var = boost::make_shared<SideVariable<double> >("e_sc");

        const int u_sc_idx = var_db->registerVariableAndContext(u_sc_var, ctx, IntVector(1));
        const int f_sc_idx = var_db->registerVariableAndContext(f_sc_var, ctx, IntVector(1));
        const int e_sc_idx = var_db->registerVariableAndContext(e_sc_var, ctx, IntVector(1));

        auto u_cc_var = boost::make_shared<CellVariable<double> >("u_cc", NDIM);
        auto f_cc_var = boost::make_shared<CellVariable<double> >("f_cc", NDIM);
        auto e_cc_var = boost::make_shared<CellVariable<double> >("e_cc", NDIM);

        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, IntVector(0));
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, IntVector(0));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVector(0));

        // Register variables for plotting.
        boost::shared_ptr<VisItDataWriter > visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity(u_cc_var->getName(), "VECTOR", u_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            ostringstream stream;
            stream << d;
            visit_data_writer->registerPlotQuantity(u_cc_var->getName() + stream.str(), "SCALAR", u_cc_idx, d);
        }

        visit_data_writer->registerPlotQuantity(f_cc_var->getName(), "VECTOR", f_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            ostringstream stream;
            stream << d;
            visit_data_writer->registerPlotQuantity(f_cc_var->getName() + stream.str(), "SCALAR", f_cc_idx, d);
        }

        visit_data_writer->registerPlotQuantity(e_cc_var->getName(), "VECTOR", e_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            ostringstream stream;
            stream << d;
            visit_data_writer->registerPlotQuantity(e_cc_var->getName() + stream.str(), "SCALAR", e_cc_idx, d);
        }

        // Initialize the AMR patch hierarchy.
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

        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            auto level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_sc_idx, 0.0);
            level->allocatePatchData(f_sc_idx, 0.0);
            level->allocatePatchData(e_sc_idx, 0.0);
            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
        }

        // Setup vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int h_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        SAMRAIVectorReal<double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        u_vec.addComponent(u_sc_var, u_sc_idx, h_sc_idx);
        f_vec.addComponent(f_sc_var, f_sc_idx, h_sc_idx);
        e_vec.addComponent(e_sc_var, e_sc_idx, h_sc_idx);

        u_vec.setToScalar(0.0);
        f_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);

        // Setup exact solutions.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        muParserCartGridFunction f_fcn("f", app_initializer->getComponentDatabase("f"), grid_geometry);

        u_fcn.setDataOnPatchHierarchy(u_sc_idx, u_sc_var, patch_hierarchy, 0.0);
        f_fcn.setDataOnPatchHierarchy(e_sc_idx, e_sc_var, patch_hierarchy, 0.0);

        // Compute -L*u = f.
        PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setCConstant(0.0);
        poisson_spec.setDConstant(-1.0);
        std::vector<RobinBcCoefStrategy*> bc_coefs(NDIM);
        SCLaplaceOperator laplace_op("laplace op");
        laplace_op.setPoissonSpecifications(poisson_spec);
        laplace_op.setPhysicalBcCoefs(bc_coefs);
        laplace_op.initializeOperatorState(u_vec, f_vec);
        laplace_op.apply(u_vec, f_vec);

        // Compute error and print error norms.
        e_vec.subtract(boost::shared_ptr<SAMRAIVectorReal<double> >(&e_vec, NullDeleter()),
                       boost::shared_ptr<SAMRAIVectorReal<double> >(&f_vec, NullDeleter()));
        pout << "|e|_oo = " << e_vec.maxNorm() << "\n";
        pout << "|e|_2  = " << e_vec.L2Norm() << "\n";
        pout << "|e|_1  = " << e_vec.L1Norm() << "\n";

        // Interpolate the side-centered data to cell centers for output.
        static const bool synch_cf_interface = true;
        hier_math_ops.interp(u_cc_idx, u_cc_var, u_sc_idx, u_sc_var, NULL, 0.0, synch_cf_interface);
        hier_math_ops.interp(f_cc_idx, f_cc_var, f_sc_idx, f_sc_var, NULL, 0.0, synch_cf_interface);
        hier_math_ops.interp(e_cc_idx, e_cc_var, e_sc_idx, e_sc_var, NULL, 0.0, synch_cf_interface);

        // Set invalid values on coarse levels (i.e., coarse-grid values that
        // are covered by finer grid patches) to equal zero.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber() - 1; ++ln)
        {
            auto level = patch_hierarchy->getPatchLevel(ln);
            BoxArray refined_region_boxes;
            auto next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
            for (auto p = level->begin(); p != level->end(); ++p)
            {
                auto patch = *p;
                const Box& patch_box = patch->getBox();
                boost::shared_ptr<CellData<double> > e_cc_data = patch->getPatchData(e_cc_idx);
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const Box refined_box = refined_region_boxes[i];
                    const Box intersection = Box::grow(patch_box, 1) * refined_box;
                    if (!intersection.empty())
                    {
                        e_cc_data->fillAll(0.0, intersection);
                    }
                }
            }
        }

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

    }

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
}
