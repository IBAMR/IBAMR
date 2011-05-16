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
#include <tbox/MathUtilities.h>
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
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <LocationIndexRobinBcCoefs.h>

#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartExtrapPhysBdryOp.h>

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

        string log_file_name = "phys_bdry_test.log";
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

        /*
         * Create major algorithm and data objects which comprise application.
         * Each object will be initialized either from input data or restart
         * files, or a combination of both.  Refer to each class constructor for
         * details.  For more information on the composition of objects for this
         * application, see comments at top of file.
         */
        tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry =
            new geom::CartesianGridGeometry<NDIM>(
                "CartesianGeometry",
                input_db->getDatabase("CartesianGeometry"));

        tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
            new hier::PatchHierarchy<NDIM>(
                "PatchHierarchy",
                grid_geometry);

        tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector =
            new mesh::StandardTagAndInitialize<NDIM>(
                "StandardTagAndInitialize",
                NULL,
                input_db->getDatabase("StandardTagAndInitialize"));

        tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator =
            new mesh::BergerRigoutsos<NDIM>();

        tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer =
            new mesh::LoadBalancer<NDIM>(
                "LoadBalancer",
                input_db->getDatabase("LoadBalancer"));

        tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > gridding_algorithm =
            new mesh::GriddingAlgorithm<NDIM>(
                "GriddingAlgorithm",
                input_db->getDatabase("GriddingAlgorithm"),
                error_detector,
                box_generator,
                load_balancer);

        /*
         * Initialize hierarchy configuration and data on all patches.
         */
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);

        int tag_buffer = 1;
        int level_number = 0;
        bool done = false;
        while (!done &&
               (gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->
                makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);

            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        /*
         * Create cell-centered data and extrapolate that data at physical
         * boundaries to obtain ghost cell values.
         */
        hier::VariableDatabase<NDIM>* var_db =
            hier::VariableDatabase<NDIM>::getDatabase();

        const int gcw = 4;
        tbox::Pointer<hier::VariableContext> context = var_db->getContext("CONTEXT");
        tbox::Pointer<pdat::CellVariable<NDIM,double> > var = new pdat::CellVariable<NDIM,double>("v");
        const int idx = var_db->registerVariableAndContext(var, context, gcw);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            tbox::Pointer<hier::PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(idx);

            for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
                const hier::Box<NDIM>& patch_box = patch->getBox();
                const hier::Index<NDIM>& patch_lower = patch_box.lower();
                tbox::Pointer<pdat::CellData<NDIM,double> > data = patch->getPatchData(idx);
                for (hier::Box<NDIM>::Iterator b(patch_box); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    (*data)(i) = 0;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        (*data)(i) += 4*(d+1)*(d+1)*i(d);
                    }
                }

                tbox::pout << "level number = " << ln << "\n";
                tbox::pout << "patch_box = " << patch_box << "\n";
                tbox::pout << "\n";

                tbox::plog << "interior data:\n";
                data->print(data->getBox());
                tbox::plog << "\n";

                CartExtrapPhysBdryOp constant_fill_op(idx, "CONSTANT");
                constant_fill_op.setPhysicalBoundaryConditions(
                    *patch, 0.0, data->getGhostCellWidth());
                tbox::plog << "constant extrapolated ghost data:\n";
                data->print(data->getGhostBox());
                tbox::plog << "\n";

                CartExtrapPhysBdryOp linear_fill_op(idx, "LINEAR");
                linear_fill_op.setPhysicalBoundaryConditions(
                    *patch, 0.0, data->getGhostCellWidth());
                tbox::plog << "linear extrapolated ghost data:\n";
                data->print(data->getGhostBox());
                tbox::plog << "\n";

                bool warning = false;
                for (hier::Box<NDIM>::Iterator b(data->getGhostBox()); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    double val = 0;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        val += 4*(d+1)*(d+1)*i(d);
                    }

                    if (!tbox::MathUtilities<double>::equalEps(val,(*data)(i)))
                    {
                        warning = true;
                        tbox::pout << "warning: value at location " << i << " is not correct\n";
                        tbox::pout << "  expected value = " << val << "   computed value = " << (*data)(i) << "\n";
                    }
                }

                if (!warning)
                {
                    tbox::pout << "linearly extrapolated boundary data appears to be correct.\n";
                }
                else
                {
                    tbox::pout << "possible errors encountered in linearly extrapolated boundary data.\n";
                }

                tbox::pout << "checking robin bc handling . . .\n";

                tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom =
                    patch->getPatchGeometry();
                const double* const xLower = pgeom->getXLower();
                const double* const ref_xUpper = grid_geometry->getXUpper();
                const double* const dx = pgeom->getDx();
                const double shift = 3.14159;
                for (hier::Box<NDIM>::Iterator b(patch_box); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    double X[NDIM];
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        X[d] = xLower[d] + dx[d]*(double(i(d)-patch_lower(d))+0.5);
                    }
                    (*data)(i) = 2.0*X[NDIM-1] + shift;
                }

                tbox::plog << "interior data:\n";
                data->print(data->getBox());
                tbox::plog << "\n";

                solv::LocationIndexRobinBcCoefs<NDIM> dirichlet_bc_coef("dirichlet_bc_coef", NULL);
                for (unsigned int d = 0; d < NDIM-1; ++d)
                {
                    dirichlet_bc_coef.setBoundarySlope(2*d  ,0.0);
                    dirichlet_bc_coef.setBoundarySlope(2*d+1,0.0);
                }
                dirichlet_bc_coef.setBoundaryValue(2*(NDIM-1)  ,                       shift);
                dirichlet_bc_coef.setBoundaryValue(2*(NDIM-1)+1,2.0*ref_xUpper[NDIM-1]+shift);

                CartCellRobinPhysBdryOp dirichlet_bc_fill_op(idx, &dirichlet_bc_coef);
                dirichlet_bc_fill_op.setPhysicalBoundaryConditions(
                    *patch, 0.0, data->getGhostCellWidth());
                tbox::plog << "extrapolated ghost data:\n";
                data->print(data->getGhostBox());
                tbox::plog << "\n";

                warning = false;
                for (hier::Box<NDIM>::Iterator b(data->getGhostBox()); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    double X[NDIM];
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        X[d] = xLower[d] + dx[d]*(double(i(d)-patch_lower(d))+0.5);
                    }
                    double val = 2.0*X[NDIM-1] + shift;

                    if (!tbox::MathUtilities<double>::equalEps(val,(*data)(i)))
                    {
                        warning = true;
                        tbox::pout << "warning: value at location " << i << " is not correct\n";
                        tbox::pout << "  expected value = " << val << "   computed value = " << (*data)(i) << "\n";
                    }
                }

                if (!warning)
                {
                    tbox::pout << "dirichlet boundary data appears to be correct.\n";
                }
                else
                {
                    tbox::pout << "possible errors encountered in extrapolated dirichlet boundary data.\n";
                }

                solv::LocationIndexRobinBcCoefs<NDIM> neumann_bc_coef("neumann_bc_coef", NULL);
                for (unsigned int d = 0; d < NDIM-1; ++d)
                {
                    neumann_bc_coef.setBoundarySlope(2*d  ,0.0);
                    neumann_bc_coef.setBoundarySlope(2*d+1,0.0);
                }
                neumann_bc_coef.setBoundarySlope(2*(NDIM-1)  ,-2.0);
                neumann_bc_coef.setBoundarySlope(2*(NDIM-1)+1,+2.0);

                CartCellRobinPhysBdryOp neumann_bc_fill_op(idx, &neumann_bc_coef);
                neumann_bc_fill_op.setPhysicalBoundaryConditions(
                    *patch, 0.0, data->getGhostCellWidth());
                tbox::plog << "extrapolated ghost data:\n";
                data->print(data->getGhostBox());
                tbox::plog << "\n";

                warning = false;
                for (hier::Box<NDIM>::Iterator b(data->getGhostBox()); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    double X[NDIM];
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        X[d] = xLower[d] + dx[d]*(double(i(d)-patch_lower(d))+0.5);
                    }
                    double val = 2.0*X[NDIM-1] + shift;

                    if (!tbox::MathUtilities<double>::equalEps(val,(*data)(i)))
                    {
                        warning = true;
                        tbox::pout << "warning: value at location " << i << " is not correct\n";
                        tbox::pout << "  expected value = " << val << "   computed value = " << (*data)(i) << "\n";
                    }
                }

                if (!warning)
                {
                    tbox::pout << "neumann boundary data appears to be correct.\n";
                }
                else
                {
                    tbox::pout << "possible errors encountered in extrapolated neumann boundary data.\n";
                }
            }
        }

    }// ensure all smart Pointers are properly deleted

    tbox::SAMRAIManager::shutdown();
    PetscFinalize();

    return 0;
}// main
