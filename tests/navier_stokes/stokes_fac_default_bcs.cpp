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

#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>
#include <ibamr/StaggeredStokesSolver.h>
#include <ibamr/StaggeredStokesSolverManager.h>

#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>

#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <cmath>
#include <iomanip>
#include <memory>
#include <vector>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

// Apply one cycle of a Stokes FAC preconditioner on a two-level hierarchy that has no boundary condition objects. A
// hierarchy with physical boundaries has a boundary helper; a fully periodic hierarchy has none.
int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app->getInputDatabase();
    const auto hierarchy_data = setup_hierarchy<NDIM>(app);
    Pointer<PatchHierarchy<NDIM>> hierarchy = std::get<0>(hierarchy_data);
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = hierarchy->getGridGeometry();
    const bool periodic = grid_geometry->getPeriodicShift().min() > 0;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("default_bcs");
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
    Pointer<SideVariable<NDIM, double>> f_u_var = new SideVariable<NDIM, double>("f_u");
    Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");
    Pointer<CellVariable<NDIM, double>> f_p_var = new CellVariable<NDIM, double>("f_p");
    const IntVector<NDIM> gcw(1);
    const int u_idx = var_db->registerVariableAndContext(u_var, ctx, gcw);
    const int f_u_idx = var_db->registerVariableAndContext(f_u_var, ctx, gcw);
    const int p_idx = var_db->registerVariableAndContext(p_var, ctx, gcw);
    const int f_p_idx = var_db->registerVariableAndContext(f_p_var, ctx, gcw);
    for (int ln = 0; ln <= 1; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (const int idx : { u_idx, f_u_idx, p_idx, f_p_idx })
        {
            level->allocatePatchData(idx);
        }
    }
    SAMRAIVectorReal<NDIM, double> x_vec("x", hierarchy, 0, 1), f_vec("f", hierarchy, 0, 1);
    x_vec.addComponent(u_var, u_idx);
    x_vec.addComponent(p_var, p_idx);
    f_vec.addComponent(f_u_var, f_u_idx);
    f_vec.addComponent(f_p_var, f_p_idx);
    x_vec.setToScalar(0.0, /*interior_only*/ false);
    f_vec.setToScalar(0.0, /*interior_only*/ false);
    muParserCartGridFunction f_u_fcn("f_u", input_db->getDatabase("f_u"), grid_geometry);
    f_u_fcn.setDataOnPatchHierarchy(f_u_idx, f_u_var, hierarchy, 0.0);

    Pointer<Database> fac_db = input_db->getDatabase("FACPreconditioner");
    Pointer<StaggeredStokesSolver> fac = StaggeredStokesSolverManager::getManager()->allocateSolver(
        input_db->getString("solver_type"), "fac", fac_db, "fac_");
    PoissonSpecifications poisson_spec("poisson_spec");
    poisson_spec.setCConstant(1.0);
    poisson_spec.setDConstant(-1.0);
    fac->setVelocityPoissonSpecifications(poisson_spec);
    fac->setComponentsHaveNullSpace(false, true);
    fac->setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, nullptr), nullptr);
    if (!periodic)
    {
        // The helper caches the homogeneous Dirichlet conditions that the solver assumes by default.
        std::vector<std::unique_ptr<LocationIndexRobinBcCoefs<NDIM>>> bc_storage;
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
        for (int d = 0; d < NDIM; ++d)
        {
            bc_storage.push_back(std::make_unique<LocationIndexRobinBcCoefs<NDIM>>("u_bc_coefs", nullptr));
            for (int location = 0; location < 2 * NDIM; ++location)
            {
                bc_storage.back()->setBoundaryValue(location, 0.0);
            }
            u_bc_coefs[d] = bc_storage.back().get();
        }
        Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
        bc_helper->cacheBcCoefData(u_bc_coefs, 0.0, hierarchy);
        fac->setPhysicalBoundaryHelper(bc_helper);
    }
    fac->setTimeInterval(0.0, 1.0);
    fac->setSolutionTime(1.0);
    fac->initializeSolverState(x_vec, f_vec);
    fac->solveSystem(x_vec, f_vec);
    fac->deallocateSolverState();

    HierarchySideDataOpsReal<NDIM, double> u_ops(hierarchy, 0, 1);
    HierarchyCellDataOpsReal<NDIM, double> p_ops(hierarchy, 0, 1);
    const double u_l2 = u_ops.L2Norm(u_idx), p_l2 = p_ops.L2Norm(p_idx);
    if (!std::isfinite(u_l2) || !std::isfinite(p_l2))
    {
        TBOX_ERROR("The FAC cycle produced a nonfinite solution.\n");
    }
    plog << std::setprecision(12) << "velocity: max norm = " << u_ops.maxNorm(u_idx) << ", L2 norm = " << u_l2 << '\n'
         << "pressure: max norm = " << p_ops.maxNorm(p_idx) << ", L2 norm = " << p_l2 << '\n';
    return 0;
}
