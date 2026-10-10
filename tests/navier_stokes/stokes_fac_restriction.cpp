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

#include <ibamr/StaggeredStokesFACPreconditioner.h>
#include <ibamr/StaggeredStokesLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>

#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <HierarchySideDataOpsReal.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <iomanip>
#include <memory>
#include <vector>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

// Restrict fine velocity data whose ghost cells hold arbitrary values with the restriction of the Stokes FAC
// preconditioner on a two-level hierarchy.
int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app->getInputDatabase();
    const auto hierarchy_data = setup_hierarchy<NDIM>(app);
    Pointer<PatchHierarchy<NDIM>> hierarchy = std::get<0>(hierarchy_data);
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = hierarchy->getGridGeometry();

    // The velocity satisfies homogeneous Dirichlet conditions on every boundary.
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

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("restriction");
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
    Pointer<SideVariable<NDIM, double>> f_u_var = new SideVariable<NDIM, double>("f_u");
    Pointer<SideVariable<NDIM, double>> r_u_var = new SideVariable<NDIM, double>("r_u");
    Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");
    Pointer<CellVariable<NDIM, double>> f_p_var = new CellVariable<NDIM, double>("f_p");
    Pointer<CellVariable<NDIM, double>> r_p_var = new CellVariable<NDIM, double>("r_p");
    const IntVector<NDIM> gcw(input_db->getIntegerWithDefault("vector_ghost_cell_width", 1));
    const int u_idx = var_db->registerVariableAndContext(u_var, ctx, gcw);
    const int f_u_idx = var_db->registerVariableAndContext(f_u_var, ctx, gcw);
    const int r_u_idx = var_db->registerVariableAndContext(r_u_var, ctx, gcw);
    const int p_idx = var_db->registerVariableAndContext(p_var, ctx, gcw);
    const int f_p_idx = var_db->registerVariableAndContext(f_p_var, ctx, gcw);
    const int r_p_idx = var_db->registerVariableAndContext(r_p_var, ctx, gcw);
    for (int ln = 0; ln <= 1; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (const int idx : { u_idx, f_u_idx, r_u_idx, p_idx, f_p_idx, r_p_idx })
        {
            level->allocatePatchData(idx);
        }
    }
    SAMRAIVectorReal<NDIM, double> x_vec("x", hierarchy, 0, 1), f_vec("f", hierarchy, 0, 1),
        r_vec("r", hierarchy, 0, 1);
    x_vec.addComponent(u_var, u_idx);
    x_vec.addComponent(p_var, p_idx);
    f_vec.addComponent(f_u_var, f_u_idx);
    f_vec.addComponent(f_p_var, f_p_idx);
    r_vec.addComponent(r_u_var, r_u_idx);
    r_vec.addComponent(r_p_var, r_p_idx);
    x_vec.setToScalar(0.0, /*interior_only*/ false);
    f_vec.setToScalar(0.0, /*interior_only*/ false);
    r_vec.setToScalar(0.0, /*interior_only*/ false);

    // The fine ghost cells hold a large value; the coarse data are zero.
    HierarchySideDataOpsReal<NDIM, double> coarse_ops(hierarchy, 0, 0), fine_ops(hierarchy, 1, 1);
    fine_ops.setToScalar(f_u_idx, 1000.0, /*interior_only*/ false);
    muParserCartGridFunction f_u_fcn("f_u", input_db->getDatabase("f_u"), grid_geometry);
    f_u_fcn.setDataOnPatchHierarchy(f_u_idx, f_u_var, hierarchy, 0.0, /*initial_time*/ false, 1, 1);

    Pointer<Database> fac_db = input_db->getDatabase("FACPreconditioner");
    Pointer<StaggeredStokesLevelRelaxationFACOperator> fac_op =
        new StaggeredStokesLevelRelaxationFACOperator("fac_op", fac_db, "fac_");
    StaggeredStokesFACPreconditioner fac("fac", fac_op, fac_db, "fac_");
    PoissonSpecifications poisson_spec("poisson_spec");
    poisson_spec.setCConstant(1.0);
    poisson_spec.setDConstant(-1.0);
    fac.setVelocityPoissonSpecifications(poisson_spec);
    fac.setComponentsHaveNullSpace(false, true);
    fac.setPhysicalBcCoefs(u_bc_coefs, nullptr);
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
    bc_helper->cacheBcCoefData(u_bc_coefs, 0.0, hierarchy);
    fac.setPhysicalBoundaryHelper(bc_helper);
    fac.setTimeInterval(0.0, 1.0);
    fac.setSolutionTime(1.0);
    fac.initializeSolverState(x_vec, f_vec);
    fac_op->restrictResidual(f_vec, r_vec, 0);
    fac.deallocateSolverState();

    plog << std::setprecision(12) << "restricted velocity: max norm = " << coarse_ops.maxNorm(r_u_idx)
         << ", L2 norm = " << coarse_ops.L2Norm(r_u_idx) << '\n'
         << "fine velocity: max norm = " << fine_ops.maxNorm(f_u_idx) << ", L2 norm = " << fine_ops.L2Norm(f_u_idx)
         << '\n';
    return 0;
}
