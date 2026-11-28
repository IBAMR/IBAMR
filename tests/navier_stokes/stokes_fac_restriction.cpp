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

#include <ibamr/StaggeredStokesBoxRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesFACPreconditioner.h>
#include <ibamr/StaggeredStokesLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/SideNoCornersFillPattern.h>
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

#include <iomanip>
#include <memory>
#include <vector>

#include "../coarse_fine_ghost_values.h"
#include "../tests.h"

#include <ibamr/app_namespaces.h>

using namespace coarse_fine_ghost_values;

// Restrict fine velocity data whose ghost cells contain arbitrary values with the restriction of the Stokes FAC
// preconditioner on a two-level hierarchy. With check_coarse_fine_ghost_values, compare instead the coarse-fine ghost
// values that the smoother sets with the values that the ghost cell fill of the residual sets, after a prolongation and
// a post-smoothing sweep and after a pre-smoothing sweep from zero followed by a post-smoothing sweep.
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

    // The fine ghost cells contain a large value. The coarse data are zero.
    HierarchySideDataOpsReal<NDIM, double> coarse_ops(hierarchy, 0, 0), fine_ops(hierarchy, 1, 1);
    fine_ops.setToScalar(f_u_idx, 1000.0, /*interior_only*/ false);
    muParserCartGridFunction f_u_fcn("f_u", input_db->getDatabase("f_u"), grid_geometry);
    f_u_fcn.setDataOnPatchHierarchy(f_u_idx, f_u_var, hierarchy, 0.0, /*initial_time*/ false, 1, 1);

    Pointer<Database> fac_db = input_db->getDatabase("FACPreconditioner");
    Pointer<Database> test_db = input_db->getDatabase("test");
    const bool check_coarse_fine_ghost_values = test_db->getBoolWithDefault("check_coarse_fine_ghost_values", false);
    Pointer<StaggeredStokesFACPreconditionerStrategy> fac_op;
    if (test_db->getStringWithDefault("operator_type", "LEVEL") == "BOX")
    {
        fac_op = new StaggeredStokesBoxRelaxationFACOperator("fac_op", fac_db, "fac_");
    }
    else
    {
        fac_op = new StaggeredStokesLevelRelaxationFACOperator("fac_op", fac_db, "fac_");
    }
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
    if (check_coarse_fine_ghost_values)
    {
        Pointer<SideVariable<NDIM, double>> w_u_var = new SideVariable<NDIM, double>("w_u");
        Pointer<CellVariable<NDIM, double>> w_p_var = new CellVariable<NDIM, double>("w_p");
        const int w_u_idx = var_db->registerVariableAndContext(w_u_var, ctx, gcw);
        const int w_p_idx = var_db->registerVariableAndContext(w_p_var, ctx, gcw);
        for (int ln = 0; ln <= 1; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            level->allocatePatchData(w_u_idx);
            level->allocatePatchData(w_p_idx);
        }
        HierarchySideDataOpsReal<NDIM, double> side_ops(hierarchy, 0, 1);
        HierarchyCellDataOpsReal<NDIM, double> cell_ops(hierarchy, 0, 1);

        // The ghost cell fill of the residual, without coarsening.
        Pointer<VariableFillPattern<NDIM>> sc_fill_pattern =
            new SideNoCornersFillPattern(1, /*overwrite_interior*/ false);
        Pointer<VariableFillPattern<NDIM>> cc_fill_pattern =
            new CellNoCornersFillPattern(1, /*overwrite_interior*/ false);
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> w_comps(2);
        w_comps[0] = InterpolationTransactionComponent(
            w_u_idx, "NONE", /*use_cf_bdry_interpolation*/ true, "NONE", "LINEAR", false, u_bc_coefs, sc_fill_pattern);
        w_comps[1] = InterpolationTransactionComponent(
            w_p_idx, "NONE", /*use_cf_bdry_interpolation*/ true, "NONE", "LINEAR", false, nullptr, cc_fill_pattern);
        HierarchyGhostCellInterpolation w_fill;
        w_fill.initializeOperatorState(w_comps, hierarchy, 0, 1);
        w_fill.setHomogeneousBc(true);
        auto set_f = [&]()
        {
            set_coarse_fine_ghost_test_data(hierarchy, f_u_idx, 1.0);
            set_coarse_fine_ghost_test_data(hierarchy, f_p_idx, 2.0);
        };
        auto copy_x = [&]()
        {
            side_ops.copyData(w_u_idx, u_idx, /*interior_only*/ false);
            cell_ops.copyData(w_p_idx, p_idx, /*interior_only*/ false);
        };
        auto compare = [&](const std::string& name)
        {
            w_fill.fillData(0.0);
            plog << name << ": compares "
                 << compare_coarse_fine_ghost_values(hierarchy, 1, u_idx, w_u_idx, name + " velocity")
                 << " velocity and "
                 << compare_coarse_fine_ghost_values(hierarchy, 1, p_idx, w_p_idx, name + " pressure")
                 << " pressure coarse-fine ghost values\n";
        };

        // A prolongation and a post-smoothing sweep.
        set_coarse_fine_ghost_test_data(hierarchy, u_idx, 0.0);
        set_coarse_fine_ghost_test_data(hierarchy, p_idx, 3.0);
        set_f();
        fac_op->prolongErrorAndCorrect(x_vec, x_vec, 1);
        copy_x();
        fac_op->smoothError(x_vec, f_vec, 1, 1, /*performing_pre_sweeps*/ false, /*performing_post_sweeps*/ true);
        compare("post-smoothing");

        // Pre-smoothing from zero, then a post-smoothing sweep; the coarse correction is zero.
        x_vec.setToScalar(0.0, /*interior_only*/ false);
        set_f();
        fac_op->smoothError(x_vec, f_vec, 1, 1, /*performing_pre_sweeps*/ true, /*performing_post_sweeps*/ false);
        copy_x();
        fac_op->smoothError(x_vec, f_vec, 1, 1, /*performing_pre_sweeps*/ false, /*performing_post_sweeps*/ true);
        compare("zero-start smoothing");
        w_fill.deallocateOperatorState();
        fac.deallocateSolverState();
        return 0;
    }
    fac_op->restrictResidual(f_vec, r_vec, 0);
    fac.deallocateSolverState();

    plog << std::setprecision(12) << "restricted velocity: max norm = " << coarse_ops.maxNorm(r_u_idx)
         << ", L2 norm = " << coarse_ops.L2Norm(r_u_idx) << '\n'
         << "fine velocity: max norm = " << fine_ops.maxNorm(f_u_idx) << ", L2 norm = " << fine_ops.L2Norm(f_u_idx)
         << '\n';
    return 0;
}
