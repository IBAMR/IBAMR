// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyAveragedDataManager.h>
#include <ibtk/IBTKInit.h>

#include <tbox/InputDatabase.h>
#include <tbox/Logger.h>

#include <CartesianGridGeometry.h>
#include <CellVariable.h>

#include <limits>
#include <set>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

// Construct an IBTK::HierarchyAveragedDataManager with a NaN snapshot time point, which stops with an error.
int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Do not print file names or line numbers to the output file.
    Pointer<Logger::Appender> abort_append(new TestAppender());
    Logger::getInstance()->setAbortAppender(abort_append);

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

    Pointer<Database> avg_manager_db = new InputDatabase("AvgManager");
    avg_manager_db->putBool("output_data", false);
    avg_manager_db->putBool("enable_logging", false);

    // A std::set<double> holds NaN only as its sole element, since NaN compares equivalent to every value.
    const std::set<double> time_pts{ std::numeric_limits<double>::quiet_NaN() };
    HierarchyAveragedDataManager avg_manager("AvgManager",
                                             new CellVariable<NDIM, double>("c_var"),
                                             avg_manager_db,
                                             time_pts,
                                             0.0,
                                             1.0,
                                             0.2,
                                             grid_geometry,
                                             /*register_for_restart*/ false);
    return 0;
} // main
