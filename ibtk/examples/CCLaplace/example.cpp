// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files
#include <IBTK_config.h>

#include <SAMRAI_config.h>

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Since SAMRAI and PETSc both require finalization routines we have to
    // ensure that no SAMRAI or PETSc objects are active at the point where we
    // call SAMRAIManager::shutdown() or PetscFinalize. Hence, to guarantee
    // that all objects are cleaned up by that point, we put everything we use
    // in an inner scope.
    {
        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_laplace.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application. These objects are configured from the input
        // database. Nearly all SAMRAI applications (at least those in IBAMR)
        // start by setting up the same half-dozen objects.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        // We create a variable for every vector we ultimately declare,
        // instead of creating and then cloning vectors. The rationale for
        // this is given below.
        Pointer<CellVariable<NDIM, double> > u_cc_var = new CellVariable<NDIM, double>("u_cc");
        Pointer<CellVariable<NDIM, double> > f_cc_var = new CellVariable<NDIM, double>("f_cc");
        Pointer<CellVariable<NDIM, double> > e_cc_var = new CellVariable<NDIM, double>("e_cc");
        Pointer<CellVariable<NDIM, double> > f_approx_cc_var = new CellVariable<NDIM, double>("f_approx_cc");

        // Internally, SAMRAI keeps track of variables (and their
        // corresponding vectors, data, etc.) by converting them to
        // indices. Here we get the indices after notifying the variable
        // database about them.
        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, IntVector<NDIM>(1));
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, IntVector<NDIM>(1));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVector<NDIM>(1));
        const int f_approx_cc_idx = var_db->registerVariableAndContext(f_approx_cc_var, ctx, IntVector<NDIM>(1));

        // Register variables for plotting.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity(u_cc_var->getName(), "SCALAR", u_cc_idx);
        visit_data_writer->registerPlotQuantity(f_cc_var->getName(), "SCALAR", f_cc_idx);
        visit_data_writer->registerPlotQuantity(f_approx_cc_var->getName(), "SCALAR", f_approx_cc_idx);
        visit_data_writer->registerPlotQuantity(e_cc_var->getName(), "SCALAR", e_cc_idx);

        // Initialize the AMR patch hierarchy. This sets up the coarsest level
        // (level 0) as well as any other levels specified in the input
        // file. We normally use the value tag_buffer to specify the number of
        // cells between a patch on level N and a patch on level N - 2:
        // however, SAMRAI ignores this value when setting up a hierarchy from
        // an input file so we just set it to an invalid value.
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        const int tag_buffer = std::numeric_limits<int>::max();
        int level_number = 0;
        while ((gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            ++level_number;
        }

        const int finest_level = patch_hierarchy->getFinestLevelNumber();

        // Allocate data for each variable on each level of the patch
        // hierarchy.
        for (int ln = 0; ln <= finest_level; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
            level->allocatePatchData(f_approx_cc_idx, 0.0);
        }

        // By default, the norms defined on SAMRAI vectors are vectors in R^n:
        // however, in IBAMR we almost always want to use a norm that
        // corresponds to a numerical quadrature. To do this we have to
        // associate each vector with a set of cell-centered volumes. Rather
        // than set this up manually, we rely on an IBTK utility class that
        // computes this (as well as many other things!). These values are
        // known as `cell weights' in this context, so we get the ID of the
        // associated data by asking for that. Behind the scenes
        // HierarchyMathOps sets up the necessary cell-centered variables and
        // registers them with the usual SAMRAI objects: all we need to do is
        // ask for the ID. Due to the way SAMRAI works these calls must occur
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int cv_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        // SAMRAI patches do not store data as a single contiguous arrays;
        // instead, each hierarchy contains several contiguous arrays. Hence,
        // to do linear algebra, we rely on SAMRAI's own vector class which
        // understands these relationships. We begin by initializing each
        // vector with the patch hierarchy:
        SAMRAIVectorReal<NDIM, double> u_vec("u", patch_hierarchy, 0, finest_level);
        SAMRAIVectorReal<NDIM, double> f_vec("f", patch_hierarchy, 0, finest_level);
        SAMRAIVectorReal<NDIM, double> f_approx_vec("f_approx", patch_hierarchy, 0, finest_level);
        SAMRAIVectorReal<NDIM, double> e_vec("e", patch_hierarchy, 0, finest_level);

        // and then associate them with, in each case, the relevant
        // component. Note that adding the components in this way will
        // register the vector with the visit data writer declared above and
        // we will compute cell integrals over the entire domain with respect
        // to the control volumes defined by cv_cc_idx.
        u_vec.addComponent(u_cc_var, u_cc_idx, cv_cc_idx);
        f_vec.addComponent(f_cc_var, f_cc_idx, cv_cc_idx);
        f_approx_vec.addComponent(f_approx_cc_var, f_approx_cc_idx, cv_cc_idx);
        e_vec.addComponent(e_cc_var, e_cc_idx, cv_cc_idx);

        // By default, in SAMRAI, if we create another vector as
        //
        //     SAMRAIVectorReal<NDIM, double> f_approx_vec("f_approx", patch_hierarchy, 0, finest_level);
        //
        // we will simply get a shallow copy of the u vector: put another way,
        // the two vectors will have different names but refer to the same
        // numerical values. Unfortunately cloning the vector doesn't work
        // either: the following code
        //
        //    tbox::Pointer<SAMRAIVectorReal<NDIM, double> > f_approx = u_vec.cloneVector("f_approx");
        //    SAMRAIVectorReal<NDIM, double> &f_approx_vec = *f_approx;
        //    f_approx_vec.setToScalar(0.0, false);
        //
        // crashes in SAMRAI 2.4.4 with a failed assertion referring to an
        // unknown variable ID. While ambiguous, the error message is not
        // wrong: we have to explicitly allocate data for each variable, so
        // creating a new anonymous variable for a cloned vector does not make
        // much sense.

        // Zero the vectors, including possible ghost data:
        u_vec.setToScalar(0.0, false);
        f_vec.setToScalar(0.0, false);
        f_approx_vec.setToScalar(0.0, false);
        e_vec.setToScalar(0.0, false);

        // Next, we use functions defined with muParser to set up the right
        // hand side and solution. These functions are read from the input
        // database and can be changed without recompiling.
        {
            muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
            muParserCartGridFunction f_fcn("f", app_initializer->getComponentDatabase("f"), grid_geometry);

            u_fcn.setDataOnPatchHierarchy(u_cc_idx, u_cc_var, patch_hierarchy, 0.0);
            f_fcn.setDataOnPatchHierarchy(f_cc_idx, f_cc_var, patch_hierarchy, 0.0);
        }

        // Compute -L*u = f.
        PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setCConstant(0.0);
        poisson_spec.setDConstant(-1.0);
        RobinBcCoefStrategy<NDIM>* bc_coef = NULL;
        CCLaplaceOperator laplace_op("laplace op");
        laplace_op.setPoissonSpecifications(poisson_spec);
        laplace_op.setPhysicalBcCoef(bc_coef);
        laplace_op.initializeOperatorState(u_vec, f_vec);
        laplace_op.apply(u_vec, f_approx_vec);

        // Compute error and print error norms. Here we create temporary smart
        // pointers that will not delete the underlying object since the
        // second argument to the constructor is false.
        e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&f_vec, false),
                       Pointer<SAMRAIVectorReal<NDIM, double> >(&f_approx_vec, false));
        pout << "|e|_oo = " << e_vec.maxNorm() << "\n";
        pout << "|e|_2  = " << e_vec.L2Norm() << "\n";
        pout << "|e|_1  = " << e_vec.L1Norm() << "\n";

        // Finally, we clean up the output by setting error values on patches
        // on coarser levels which are covered by finer levels to zero.
        for (int ln = 0; ln < finest_level; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            Pointer<PatchLevel<NDIM> > next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
            BoxArray<NDIM> refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                const Patch<NDIM>& patch = *level->getPatch(p());
                const Box<NDIM>& patch_box = patch.getBox();
                Pointer<CellData<NDIM, double> > e_cc_data = patch.getPatchData(e_cc_idx);
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const Box<NDIM>& refined_box = refined_region_boxes[i];
                    // Box::operator* returns the intersection of two boxes.
                    const Box<NDIM>& intersection = patch_box * refined_box;
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
} // main
