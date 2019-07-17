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
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartExtrapPhysBdryOp.h>

#include <LocationIndexRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

double
const_f(const double* const /*x*/)
{
    return 1.0;
}

double
linear_f(const double* const x)
{
    double val = 0.0;
    for (int d = 0; d < NDIM; ++d)
    {
        val += static_cast<double>(d + 1) * x[d];
    }
    return val;
}

double
quadratic_f(const double* const x)
{
    double val = 0.0;
    for (int d = 0; d < NDIM; ++d)
    {
        val += static_cast<double>(d + 1) * x[d] * x[d];
    }
    return val;
}

using fcn_type = double (*)(const double* const);

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
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    { // cleanup dynamically allocated objects prior to shutdown
        // prevent a warning about timer initializations
        TimerManager::createManager(nullptr);

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_poisson.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
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

        // Get extrapolation type and variable centering
        std::string extrap_type = input_db->getString("EXTRAP_TYPE");
        std::string var_centering = input_db->getString("VAR_CENTERING");
        TBOX_ASSERT(extrap_type == "CONSTANT" || extrap_type == "LINEAR" || extrap_type == "QUADRATIC");
        TBOX_ASSERT(var_centering == "CELL" || var_centering == "SIDE" || var_centering == "FACE" ||
                    var_centering == "NODE");

        // Create cell-centered data and extrapolate that data at physical
        // boundaries to obtain ghost cell values.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = var_db->getContext("CONTEXT");
        Pointer<CellVariable<NDIM, double> > c_var = new CellVariable<NDIM, double>("c_u");
        Pointer<SideVariable<NDIM, double> > s_var = new SideVariable<NDIM, double>("s_u");
        Pointer<FaceVariable<NDIM, double> > f_var = new FaceVariable<NDIM, double>("f_u");
        Pointer<NodeVariable<NDIM, double> > n_var = new NodeVariable<NDIM, double>("n_u");
        const int gcw = 4;
        const int c_idx = var_db->registerVariableAndContext(c_var, context, gcw);
        const int s_idx = var_db->registerVariableAndContext(s_var, context, gcw);
        const int f_idx = var_db->registerVariableAndContext(f_var, context, gcw);
        const int n_idx = var_db->registerVariableAndContext(n_var, context, gcw);
        std::map<std::string, int> typeMap;
        std::map<std::string, fcn_type> fcnMap;
        typeMap["CELL"] = c_idx;
        typeMap["SIDE"] = s_idx;
        typeMap["FACE"] = f_idx;
        typeMap["NODE"] = n_idx;
        fcnMap["CONSTANT"] = &const_f;
        fcnMap["LINEAR"] = &linear_f;
        fcnMap["QUADRATIC"] = &quadratic_f;
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(typeMap[var_centering]);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const hier::Index<NDIM>& patch_lower = patch_box.lower();
                Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
                const double* const x_low = p_geom->getXLower();
                const double* const dx = p_geom->getDx();
                Pointer<CellData<NDIM, double> > c_data = patch->getPatchData(typeMap[var_centering]);
                Pointer<SideData<NDIM, double> > s_data = patch->getPatchData(typeMap[var_centering]);
                Pointer<FaceData<NDIM, double> > f_data = patch->getPatchData(typeMap[var_centering]);
                Pointer<NodeData<NDIM, double> > n_data = patch->getPatchData(typeMap[var_centering]);
                std::vector<double> x(NDIM);
                if (var_centering == "CELL")
                {
                    for (CellIterator<NDIM> ci(patch_box); ci; ci++)
                    {
                        const CellIndex<NDIM>& i = ci();
                        for (int d = 0; d < NDIM; ++d)
                        {
                            x[d] = x_low[d] + dx[d] * (i(d) - patch_lower(d) + 0.5);
                        }
                        auto fcn = fcnMap[extrap_type];
                        (*c_data)(i) = fcn(x.data());
                    }
                }
                else if (var_centering == "SIDE")
                {
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (SideIterator<NDIM> si(patch_box, axis); si; si++)
                        {
                            const SideIndex<NDIM>& i = si();
                            for (int d = 0; d < NDIM; ++d)
                            {
                                x[d] = x_low[d] + dx[d] * (i(d) - patch_lower(d) + (axis == d ? 0.0 : 0.5));
                            }
                            auto fcn = fcnMap[extrap_type];
                            (*s_data)(i) = fcn(x.data());
                        }
                    }
                }
                else if (var_centering == "FACE")
                {
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (FaceIterator<NDIM> fi(patch_box, axis); fi; fi++)
                        {
                            const FaceIndex<NDIM>& i = fi();
                            for (int d = 0; d < NDIM; ++d)
                            {
                                x[d] = x_low[d] + dx[d] * (i(d) - patch_lower(d) + (axis == d ? 0.0 : 0.5));
                            }
                            auto fcn = fcnMap[extrap_type];
                            (*f_data)(i) = fcn(x.data());
                        }
                    }
                }
                else if (var_centering == "NODE")
                {
                    for (NodeIterator<NDIM> ni(patch_box); ni; ni++)
                    {
                        const NodeIndex<NDIM>& i = ni();
                        for (int d = 0; d < NDIM; ++d)
                        {
                            x[d] = x_low[d] + dx[d] * (i(d) - patch_lower(d));
                        }
                        auto fcn = fcnMap[extrap_type];
                        (*n_data)(i) = fcn(x.data());
                    }
                }
                else
                {
                    TBOX_ERROR("UNKNOWN DATA CENTERING: " << var_centering << "\n");
                }
                CartExtrapPhysBdryOp fill_op(typeMap[var_centering], extrap_type);
                fill_op.setPhysicalBoundaryConditions(*patch, 0.0, IntVector<NDIM>(gcw));
                bool warning = false;
                if (var_centering == "CELL")
                {
                    for (CellIterator<NDIM> ci(c_data->getGhostBox()); ci; ci++)
                    {
                        const CellIndex<NDIM>& i = ci();
                        for (int d = 0; d < NDIM; ++d)
                        {
                            x[d] = x_low[d] + dx[d] * (i(d) - patch_lower(d) + 0.5);
                        }
                        auto fcn = fcnMap[extrap_type];
                        double val = fcn(x.data());
                        if (!MathUtilities<double>::equalEps(val, (*c_data)(i)))
                        {
                            warning = true;
                            pout << "warning: value at location " << i << " is not correct\n";
                            pout << "  expected value = " << val << "   computed value = " << (*c_data)(i) << "\n";
                        }
                    }
                }
                else if (var_centering == "SIDE")
                {
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (SideIterator<NDIM> si(s_data->getGhostBox(), axis); si; si++)
                        {
                            const SideIndex<NDIM>& i = si();
                            for (int d = 0; d < NDIM; ++d)
                            {
                                x[d] = x_low[d] + dx[d] * (i(d) - patch_lower(d) + (axis == d ? 0.0 : 0.5));
                            }
                            auto fcn = fcnMap[extrap_type];
                            double val = fcn(x.data());
                            if (!MathUtilities<double>::equalEps(val, (*s_data)(i)))
                            {
                                warning = true;
                                pout << "warning: value at location " << i << " is not correct\n";
                                pout << "  expected value = " << val << "   computed value = " << (*s_data)(i) << "\n";
                            }
                        }
                    }
                }
                else if (var_centering == "FACE")
                {
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (FaceIterator<NDIM> fi(f_data->getGhostBox(), axis); fi; fi++)
                        {
                            const FaceIndex<NDIM>& i = fi();
                            for (int d = 0; d < NDIM; ++d)
                            {
                                x[d] = x_low[d] + dx[d] * (i(d) - patch_lower(d) + (axis == d ? 0.0 : 0.5));
                            }
                            auto fcn = fcnMap[extrap_type];
                            double val = fcn(x.data());
                            if (!MathUtilities<double>::equalEps(val, (*f_data)(i)))
                            {
                                warning = true;
                                pout << "warning: value at location " << i << " is not correct\n";
                                pout << "  expected value = " << val << "   computed value = " << (*f_data)(i) << "\n";
                            }
                        }
                    }
                }
                else if (var_centering == "NODE")
                {
                    for (NodeIterator<NDIM> ni(n_data->getGhostBox()); ni; ni++)
                    {
                        const NodeIndex<NDIM>& i = ni();
                        for (int d = 0; d < NDIM; ++d)
                        {
                            x[d] = x_low[d] + dx[d] * (i(d) - patch_lower(d));
                        }
                        auto fcn = fcnMap[extrap_type];
                        double val = fcn(x.data());
                        if (!MathUtilities<double>::equalEps(val, (*n_data)(i)))
                        {
                            warning = true;
                            pout << "warning: value at location " << i << " is not correct\n";
                            pout << "  expected value = " << val << "   computed value = " << (*n_data)(i) << "\n";
                        }
                    }
                }
                else
                {
                    TBOX_ERROR("UNKNOWN DATA CENTERING: " << var_centering << "\n");
                }

                if (!warning)
                {
                    pout << "extrapolated boundary data appears to be correct.\n";
                }
                else
                {
                    pout << "possible errors encountered in extrapolated boundary data.\n";
                }
            }
            level->deallocatePatchData(typeMap[var_centering]);
        }

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
} // main
