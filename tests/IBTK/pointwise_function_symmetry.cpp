// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
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
#include <ibtk/CartGridPointwiseFunction.h>
#include <ibtk/IBTKInit.h>

#include <tbox/Logger.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

// Evaluate a pointwise function that returns a nonsymmetric tensor for data
// with symmetric tensor storage. The library should reject the result.

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv);
    Logger::getInstance()->setAbortAppender(new TestAppender());

    Pointer<CartesianGridGeometry<NDIM>> geometry =
        new CartesianGridGeometry<NDIM>("geometry", app->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("hierarchy", geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector =
        new StandardTagAndInitialize<NDIM>("tagging", nullptr, app->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM>> boxes = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> load = new LoadBalancer<NDIM>("load", app->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding = new GriddingAlgorithm<NDIM>(
        "gridding", app->getComponentDatabase("GriddingAlgorithm"), error_detector, boxes, load);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<CellVariable<NDIM, double>> var = new CellVariable<NDIM, double>("tensor", NDIM * (NDIM + 1) / 2);
    const int idx = var_db->registerVariableAndContext(var, var_db->getContext("context"), IntVector<NDIM>(0));
    gridding->makeCoarsestLevel(hierarchy, 0.0);
    hierarchy->getPatchLevel(0)->allocatePatchData(idx);

    Pointer<CartGridFunction> function = make_cart_grid_pointwise_function<MatrixNd>(
        "nonsymmetric",
        var,
        [](const VectorNd&, double, int, int) -> MatrixNd
        {
            MatrixNd q = MatrixNd::Identity();
            q(0, 1) = 1.0;
            return q;
        },
        TensorStorage::SYMMETRIC);
    function->setDataOnPatchHierarchy(idx, var, hierarchy, 0.0);

    // The evaluation should have failed. Return normally so that the
    // expected-error case is rejected.
    return 0;
}
