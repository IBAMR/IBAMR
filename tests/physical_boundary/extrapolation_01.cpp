// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
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

#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAI_config.h>

// Headers for basic PETSc objects
#include <petscsys.h>

// SAMRAI INCLUDES
#include <SAMRAIBergerRigoutsos.h>
#include <SAMRAIBox.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAICartesianPatchGeometry.h>
#include <SAMRAICellData.h>
#include <SAMRAICellIndex.h>
#include <SAMRAICellIterator.h>
#include <SAMRAICellVariable.h>
#include <SAMRAIFaceData.h>
#include <SAMRAIFaceIndex.h>
#include <SAMRAIFaceIterator.h>
#include <SAMRAIFaceVariable.h>
#include <SAMRAIGriddingAlgorithm.h>
#include <SAMRAIIndex.h>
#include <SAMRAIIntVector.h>
#include <SAMRAILoadBalancer.h>
#include <SAMRAILocationIndexRobinBcCoefs.h>
#include <SAMRAINodeData.h>
#include <SAMRAINodeIndex.h>
#include <SAMRAINodeIterator.h>
#include <SAMRAINodeVariable.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIPointer.h>
#include <SAMRAISideData.h>
#include <SAMRAISideIndex.h>
#include <SAMRAISideIterator.h>
#include <SAMRAISideVariable.h>
#include <SAMRAIStandardTagAndInitialize.h>
#include <SAMRAIVariableDatabase.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/ibtk_utilities.h>

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
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown
        // prevent a warning about timer initializations
        TimerManager::createManager(nullptr);

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        SAMRAIPointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_poisson.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geometry = new SAMRAICartesianGridGeometry(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy = new SAMRAIPatchHierarchy("PatchHierarchy", grid_geometry);
        SAMRAIPointer<SAMRAIStandardTagAndInitialize> error_detector = new SAMRAIStandardTagAndInitialize(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        SAMRAIPointer<SAMRAIBergerRigoutsos> box_generator = new SAMRAIBergerRigoutsos();
        SAMRAIPointer<SAMRAILoadBalancer> load_balancer =
            new SAMRAILoadBalancer("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        SAMRAIPointer<SAMRAIGriddingAlgorithm> gridding_algorithm =
            new SAMRAIGriddingAlgorithm("GriddingAlgorithm",
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
        SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
        SAMRAIPointer<VariableContext> context = var_db->getContext("CONTEXT");
        SAMRAIPointer<SAMRAICellVariable<double>> c_var = new SAMRAICellVariable<double>("c_u");
        SAMRAIPointer<SAMRAISideVariable<double>> s_var = new SAMRAISideVariable<double>("s_u");
        SAMRAIPointer<SAMRAIFaceVariable<double>> f_var = new SAMRAIFaceVariable<double>("f_u");
        SAMRAIPointer<SAMRAINodeVariable<double>> n_var = new SAMRAINodeVariable<double>("n_u");
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
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(typeMap[var_centering]);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                const SAMRAIBox& patch_box = patch->getBox();
                const SAMRAIIndex& patch_lower = patch_box.lower();
                SAMRAIPointer<SAMRAICartesianPatchGeometry> p_geom = patch->getPatchGeometry();
                const double* const x_low = p_geom->getXLower();
                const double* const dx = p_geom->getDx();
                SAMRAIPointer<SAMRAICellData<double>> c_data = patch->getPatchData(typeMap[var_centering]);
                SAMRAIPointer<SAMRAISideData<double>> s_data = patch->getPatchData(typeMap[var_centering]);
                SAMRAIPointer<SAMRAIFaceData<double>> f_data = patch->getPatchData(typeMap[var_centering]);
                SAMRAIPointer<SAMRAINodeData<double>> n_data = patch->getPatchData(typeMap[var_centering]);
                std::vector<double> x(NDIM);
                if (var_centering == "CELL")
                {
                    for (SAMRAICellIterator ci(patch_box); ci; ci++)
                    {
                        const SAMRAICellIndex& i = ci();
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
                        for (SAMRAISideIterator si(patch_box, axis); si; si++)
                        {
                            const SAMRAISideIndex& i = si();
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
                        for (SAMRAIFaceIterator fi(patch_box, axis); fi; fi++)
                        {
                            const SAMRAIFaceIndex& i = fi();
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
                    for (SAMRAINodeIterator ni(patch_box); ni; ni++)
                    {
                        const SAMRAINodeIndex& i = ni();
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
                fill_op.setPhysicalBoundaryConditions(*patch, 0.0, SAMRAIIntVector(gcw));
                bool warning = false;
                if (var_centering == "CELL")
                {
                    for (SAMRAICellIterator ci(c_data->getGhostBox()); ci; ci++)
                    {
                        const SAMRAICellIndex& i = ci();
                        for (int d = 0; d < NDIM; ++d)
                        {
                            x[d] = x_low[d] + dx[d] * (i(d) - patch_lower(d) + 0.5);
                        }
                        auto fcn = fcnMap[extrap_type];
                        double val = fcn(x.data());
                        if (!IBTK::rel_equal_eps(val, (*c_data)(i)))
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
                        for (SAMRAISideIterator si(s_data->getGhostBox(), axis); si; si++)
                        {
                            const SAMRAISideIndex& i = si();
                            for (int d = 0; d < NDIM; ++d)
                            {
                                x[d] = x_low[d] + dx[d] * (i(d) - patch_lower(d) + (axis == d ? 0.0 : 0.5));
                            }
                            auto fcn = fcnMap[extrap_type];
                            double val = fcn(x.data());
                            if (!IBTK::rel_equal_eps(val, (*s_data)(i)))
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
                        for (SAMRAIFaceIterator fi(f_data->getGhostBox(), axis); fi; fi++)
                        {
                            const SAMRAIFaceIndex& i = fi();
                            for (int d = 0; d < NDIM; ++d)
                            {
                                x[d] = x_low[d] + dx[d] * (i(d) - patch_lower(d) + (axis == d ? 0.0 : 0.5));
                            }
                            auto fcn = fcnMap[extrap_type];
                            double val = fcn(x.data());
                            if (!IBTK::rel_equal_eps(val, (*f_data)(i)))
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
                    for (SAMRAINodeIterator ni(n_data->getGhostBox()); ni; ni++)
                    {
                        const SAMRAINodeIndex& i = ni();
                        for (int d = 0; d < NDIM; ++d)
                        {
                            x[d] = x_low[d] + dx[d] * (i(d) - patch_lower(d));
                        }
                        auto fcn = fcnMap[extrap_type];
                        double val = fcn(x.data());
                        if (!IBTK::rel_equal_eps(val, (*n_data)(i)))
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

    return 0;
} // main
